#!/usr/bin/env python3
"""Phase NP-10 - literature search for chemical pairs with no primary evidence.

Runs PROMPT 1 in ../non_peptides_prompt.txt against chemical_interactions rows that
have no `pending` or `confirmed` np_evidence row - i.e. pairs no upstream source
backs with a readable primary paper. It is the chemical-half counterpart of
cur2_gapfill.py. Since 2026-10-09, at the project owner's instruction, an answer that
passes the load gates is written straight to 'confirmed' - no separate curator read.
'confirmed' from this path means "machine-verified against the source text", which is a
weaker claim than "a person read the paper"; extraction_method='llm_litsearch_np_prompt'
identifies these rows so the distinction stays recoverable.

  --emit N     choose N unsearched pairs, fill the prompt, search PubMed, fetch the
               hits' abstracts into `papers`, and write work/NP_LIT_WORK.md for a
               read. Every query and hit list is logged to work/np10_search_log.jsonl.
  --load FILE  load the prompt's answers (JSON, see below). Each reported paper
               must pass the gates; survivors go into np_evidence as `pending`
               rows with source_id='np10_litsearch'. 'No evidence' is logged per pair.
  --status     counts.

Pair priority: (A) pairs whose only upstream citations were closed without
reading (review_only, e.g. the NeuronChat tool paper) - real pairs needing a
primary source; then (B) pairs asserted by >=2 sources; then (C) the rest.

LOAD FILE shape:
  {"<chem_interaction_id>": "No evidence",
   "<chem_interaction_id>": [{"PMID": ..., "YEAR": ..., "MOLECULE_AS_NAMED_IN_PAPER": ...,
       "RECEPTOR_AS_NAMED_IN_PAPER": ..., "EXPERIMENT_TYPE": ..., "RESULT": ...,
       "DIRECTION": ..., "AFFINITY": ..., "VERBATIM_SENTENCE": ...}, ...]}

GATES (a record failing any is logged and NOT inserted):
  1. PMID resolves at PubMed and has an abstract (fetched on load if new)
  2. not a review, not retracted (prompt rule 1)
  3. VERBATIM_SENTENCE is verbatim in a GENUINE source (prompt rule 8)
  4. MOLECULE_AS_NAMED_IN_PAPER and RECEPTOR_AS_NAMED_IN_PAPER both occur in that text
  5. EXPERIMENT_TYPE and RESULT are from the prompt's closed vocabularies

GENUINE SOURCE = papers.abstract, or papers.full_text when it is real publisher
text (PLAN.md invariants 9/10, content test in real_fulltext() - identical to
cur_integrity_check.is_curator_note, never status-based). The abstract is tried
first; a sentence that is not there is checked against the full text and the row
is stored with evidence_depth='full_text' instead of 'abstract_only'. A record
may also say "EVIDENCE_SOURCE": "full_text" to skip straight to it. This exists
because per-pair data often sits only in a table - LAT1 (SLC7A5) abstracts say
"large neutral amino acids" while the per-amino-acid Km table names tyrosine,
methionine and histidine. If the full text is not on file the record is refused
with gate3_needs_fulltext_not_on_file: run `cur_fetch_fulltext.py <pmid> --commit`
and reload. Nothing here ever writes papers.full_text (invariant 8).

Rows are written as `confirmed` when the gates pass; `pair_absent` and the other
outcomes are still not set here. The prompt's RESULT is kept in `notes`. The source of
truth for every gate remains the paper's own text, never anything this pipeline wrote
(invariant 12) - that check is what makes unattended confirmation safe at all.

Usage: python3 np10_litsearch.py --emit 20 | --load work/np10_batch1.json | --status
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

from cur_common import HERE, UA, connect, log, now

WORK = HERE / "work"
WORKFILE = WORK / "NP_LIT_WORK.md"
LOGFILE = WORK / "np10_search_log.jsonl"
NOEVFILE = WORK / "np10_no_evidence.tsv"
HGNC_CACHE = WORK / "np10_hgnc_cache.json"
PROMPT = HERE.parent / "non_peptides_prompt.txt"
SOURCE_ID = "np10_litsearch"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
PAUSE = 0.4
RETMAX = 12
N_SYN = 6

VALID_CLASS = {"structural", "direct_binding", "functional_in_vitro",
               "functional_in_vivo", "computational"}
VALID_RESULT = {"SUPPORTS", "REFUTES", "FAMILY_LEVEL", "TRANSPORTER", "ENZYME"}
ASSAY = ('(binding[tiab] OR affinity[tiab] OR EC50[tiab] OR IC50[tiab] OR Ki[tiab] OR Kd[tiab] '
         'OR agonist*[tiab] OR recombinant[tiab] OR oocyte*[tiab] OR "patch clamp"[tiab] '
         'OR currents[tiab] OR transfected[tiab] OR "HEK293"[tiab] OR crystal*[tiab] OR cryo-EM[tiab])')
GREEK = [("alpha", "α"), ("beta", "β"), ("gamma", "γ"), ("delta", "δ"), ("omega", "ω")]


def squash(s: str) -> str:
    return re.sub(r"\s+", " ", (s or "")).strip().lower()


# Curator worksheets start with a 'Title:' header and are small; real papers are
# tens of kB. Content-based, matching cur_integrity_check.is_curator_note exactly,
# because 7 genuine papers carry full text while still labelled 'paywalled'.
FT_MIN = 20000


def real_fulltext(ft: str | None) -> bool:
    """True when papers.full_text is genuine publisher text (PLAN.md invariant 9)."""
    if not ft or not ft.strip():
        return False
    return not ft.lstrip().startswith("Title:") and len(ft) >= FT_MIN


def http(url: str, data: dict | None = None, accept: str | None = None) -> bytes:
    body = urllib.parse.urlencode(data).encode() if data else None
    hdr = {"User-Agent": UA}
    if accept:
        hdr["Accept"] = accept
    for i in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, body, hdr), timeout=60) as r:
                return r.read()
        except Exception:                                    # noqa: BLE001
            if i == 3:
                raise
            time.sleep(2 * (i + 1))
    raise RuntimeError("unreachable")


# ---------------------------------------------------------------- names

SALT = re.compile(r"\b(iodide|bromide|chloride|hydrochloride|sulfate|sodium|potassium|calcium|"
                  r"magnesium|ammonium|salt|hydrate|ester|acetate)\b", re.I)


def chem_synonyms(r) -> list[str]:
    """Search-worthy names: DB name first, then synonyms that are recognisably the
    same molecule - an abbreviation (GABA, ACh, 5-HT) or a variant sharing a long
    word with the name. That excludes brand names, salts (a different molecule
    under the prompt), and the stereo-descriptor fragments np2 split off IUPAC
    names on commas ('10S', '26-diol'). Adds Greek-letter spellings, which no
    source carries but the literature uses."""
    try:
        syn = json.loads(r["synonyms"] or "[]")
    except Exception:                                        # noqa: BLE001
        syn = []
    name_words = {w for w in re.findall(r"[a-z]{6,}", r["name"].lower())}
    out, seen = [], set()

    def add(s):
        k = s.lower()
        if k not in seen and 2 <= len(s) <= 45:
            seen.add(k); out.append(s)
    add(r["name"])
    for s in sorted(syn, key=len):
        if s.count("(") != s.count(")") or s.count("[") != s.count("]") or SALT.search(s):
            continue
        if re.search(r"\d+[RS]\)|\bdelta-|^\d+[A-Za-zΑ-ω]{0,2}$", s):   # stereo-fragments / HMDB 'delta' for D-
            continue
        is_abbrev = len(s) <= 6 and sum(ch.isupper() for ch in s) >= 2 and re.search(r"[A-Za-z]{2}", s)
        shares = name_words & set(re.findall(r"[a-z]{6,}", s.lower()))
        if is_abbrev or shares:
            add(s)
    for s in list(out):
        g = s
        for word, sym in GREEK:
            g = re.sub(rf"\b{word}\b", sym, g, flags=re.I)
        if g != s:
            add(g)
    return out


def hgnc(symbol: str, cache: dict) -> dict:
    if symbol not in cache:
        try:
            d = json.loads(http(f"https://rest.genenames.org/fetch/symbol/{urllib.parse.quote(symbol)}",
                                accept="application/json"))["response"]["docs"]
            d = d[0] if d else {}
            cache[symbol] = {"name": d.get("name"), "alias": d.get("alias_symbol") or [],
                             "prev": d.get("prev_symbol") or []}
        except Exception:                                    # noqa: BLE001
            cache[symbol] = {"name": None, "alias": [], "prev": []}
        time.sleep(0.15)
    return cache[symbol]


def receptor_terms(symbol: str, h: dict) -> list[str]:
    t = [symbol] + [a for a in h["alias"] + h["prev"] if len(a) > 2]
    if h["name"]:
        t.append(h["name"])
    return list(dict.fromkeys(t))[:8]


# Ion-channel subunits are written in the literature as a Greek-letter subunit of
# a named channel ('alpha3beta4 nicotinic', 'alpha2 subunit of the GABA(A)
# receptor'), almost never by gene symbol. PubMed indexes α as 'alpha', and
# truncation (alpha3*) catches the combined forms (alpha3beta4, alpha3beta2).
SUBUNIT_FAMILIES = [
    (r"CHRNA(\d+)", "alpha{}", '(nicotinic[tiab] OR nAChR*[tiab])'),
    (r"CHRNB(\d+)", "beta{}", '(nicotinic[tiab] OR nAChR*[tiab])'),
    (r"CHRN([DEG])", None, '(nicotinic[tiab] OR nAChR*[tiab])'),
    (r"GABRA(\d+)", "alpha{}", '(GABAA*[tiab] OR "GABA(A)"[tiab] OR "GABA A"[tiab])'),
    (r"GABRB(\d+)", "beta{}", '(GABAA*[tiab] OR "GABA(A)"[tiab] OR "GABA A"[tiab])'),
    (r"GABRG(\d+)", "gamma{}", '(GABAA*[tiab] OR "GABA(A)"[tiab] OR "GABA A"[tiab])'),
    (r"GABR([DEPQ])", None, '(GABAA*[tiab] OR "GABA(A)"[tiab] OR "GABA A"[tiab])'),
    (r"GABRR(\d+)", "rho{}", '(GABA*[tiab])'),
    (r"GLRA(\d+)", "alpha{}", '("glycine receptor*"[tiab] OR GlyR*[tiab])'),
    (r"GLRB", "beta", '("glycine receptor*"[tiab] OR GlyR*[tiab])'),
]
LETTER = {"D": "delta", "E": "epsilon", "G": "gamma", "P": "pi", "Q": "theta"}


def subunit_clause(symbol: str) -> str | None:
    for pat, tmpl, ctx in SUBUNIT_FAMILIES:
        m = re.fullmatch(pat, symbol)
        if m:
            g = m.groups()
            tok = tmpl.format(*g) if tmpl else LETTER[g[0]]
            return f"({tok}*[tiab] AND {ctx})"
    return None


# ---------------------------------------------------------------- PubMed

def esearch(term: str) -> tuple[int, list[str]]:
    d = json.loads(http(EUTILS + "esearch.fcgi", {"db": "pubmed", "term": term, "retmode": "json",
                                                 "retmax": RETMAX, "sort": "relevance"}))
    time.sleep(PAUSE)
    return int(d["esearchresult"]["count"]), d["esearchresult"]["idlist"]


def efetch(pmids: list[str]) -> list[dict]:
    if not pmids:
        return []
    root = ET.fromstring(http(EUTILS + "efetch.fcgi", {"db": "pubmed", "retmode": "xml",
                                                      "id": ",".join(pmids)}))
    time.sleep(PAUSE)
    out = []
    for art in root.findall(".//PubmedArticle"):
        mc = art.find("MedlineCitation")
        a = mc.find("Article")
        if a is None:
            continue
        ti = a.find("ArticleTitle")
        parts = []
        for at in a.findall(".//Abstract/AbstractText"):
            t = "".join(at.itertext()).strip()
            if t:
                parts.append(f"{at.get('Label')}: {t}" if at.get("Label") else t)
        year = a.findtext("Journal/JournalIssue/PubDate/Year") or \
            (a.findtext("Journal/JournalIssue/PubDate/MedlineDate") or "")[:4]
        pt = [p.text for p in a.findall("PublicationTypeList/PublicationType") if p.text]
        doi = next((x.text for x in art.findall(".//ArticleIdList/ArticleId")
                    if x.get("IdType") == "doi"), None)
        out.append(dict(pmid=mc.findtext("PMID"),
                        title="".join(ti.itertext()).strip() if ti is not None else None,
                        abstract="\n".join(parts) or None,
                        journal=a.findtext("Journal/ISOAbbreviation") or a.findtext("Journal/Title"),
                        year=int(year) if year.isdigit() else None, ptypes=json.dumps(pt),
                        is_review=int(any("Review" in p for p in pt)),
                        is_retracted=int(any("Retracted Publication" in p for p in pt)), doi=doi))
    return out


def store_papers(c, recs: list[dict]) -> None:
    """Same write contract as cur_fetch_pubmed.py: metadata + abstract, never full_text."""
    for r in recs:
        c.execute("""INSERT INTO papers (pmid,doi,title,journal,year,publication_type,
                        is_review,is_retracted,fulltext_status,abstract,metadata_at)
                     VALUES (?,?,?,?,?,?,?,?,'unresolved',?,?)
                     ON CONFLICT(pmid) DO UPDATE SET
                        title=coalesce(excluded.title,title), journal=coalesce(excluded.journal,journal),
                        year=coalesce(excluded.year,year), publication_type=excluded.publication_type,
                        is_review=excluded.is_review, is_retracted=excluded.is_retracted,
                        doi=coalesce(excluded.doi,doi), abstract=coalesce(excluded.abstract,abstract),
                        metadata_at=excluded.metadata_at""",
                  (r["pmid"], r["doi"], r["title"], r["journal"], r["year"], r["ptypes"],
                   r["is_review"], r["is_retracted"], r["abstract"], now()))


def tiab(terms: list[str]) -> str:
    return "(" + " OR ".join(f'"{t.replace(chr(34), "")}"[tiab]' for t in terms) + ")"


# ---------------------------------------------------------------- selection

def searched() -> set[str]:
    if not LOGFILE.exists():
        return set()
    return {json.loads(l)["chem_interaction_id"] for l in LOGFILE.read_text().splitlines() if l.strip()}


def targets(c, n: int) -> list:
    done = searched()
    rows = c.execute("""
        SELECT ci.chem_interaction_id, ci.n_sources, ci.sources, ci.receptor_type,
               ch.chemical_id, ch.name, ch.hmdb_id, ch.pubchem_cid, ch.chebi_id, ch.synonyms,
               e.symbol receptor,
               EXISTS(SELECT 1 FROM np_evidence ne
                       WHERE ne.chem_interaction_id = ci.chem_interaction_id) had_citation
          FROM chemical_interactions ci
          JOIN chemical_entities ch ON ch.chemical_id = ci.chemical_id
          JOIN entities e ON e.entity_id = ci.receptor_entity_id
         WHERE ci.species = 'human'
           AND NOT EXISTS (SELECT 1 FROM np_evidence ne
                            WHERE ne.chem_interaction_id = ci.chem_interaction_id
                              AND ne.verification_outcome IN ('pending','confirmed'))
      ORDER BY had_citation DESC, ci.n_sources DESC, ch.name, e.symbol""").fetchall()
    return [r for r in rows if r["chem_interaction_id"] not in done][:n]


# ---------------------------------------------------------------- pair selection modes
#
# The default --emit order (above) is had_citation, then n_sources, then alphabetical.
# That walks the queue in roughly alphabetical order and does not distinguish a pair one
# resource inferred from a pair five resources measured, which is how batch 11 ended up
# spending 18 of 20 slots on dopamine x 12 serotonin receptors - pairs nobody has tested.
#
# Two targeted modes, measured against the 12,763 unanswered human pairs on 2026-10-08:
#
#   --emit-clusters N   666 clusters / 5,129 pairs. A cluster is one chemical asserted
#                       against >=4 members of one receptor family: dexamethasone against
#                       all six RAR/RXR subtypes, dopamine against 12 HTR subtypes,
#                       adenosine against 8 adrenoceptors. That shape comes from matching
#                       compound class, not measured binding, so it is where the fabricated
#                       pairs are. Crucially ONE comparative paper usually settles the whole
#                       cluster: Allenby 1993 (PMID 8380496) tested the retinoic acid
#                       stereoisomers side by side and refuted six retinoid x RXR pairs at
#                       once. So a cluster is emitted as ONE worklist entry listing every
#                       receptor, and the answer file still carries one record per pair.
#
#   --emit-multi N      1,534 pairs with >=2 sources; 301 of those also have a curated
#                       source (gtopdb/cellphonedb/cellchatdb/mebocost/mrclinkdb) rather
#                       than only the two high-density inference resources (metalinksdb
#                       25.1 receptors/chemical, metaligand 18.8, vs ~3 for the curated
#                       ones). Highest prior of being real, so this is where confirmations
#                       come from. Slower per pair, and corroboration is not independence:
#                       batch 10 was all 4-source pairs and still returned 12/20 No evidence.
#
# Neither mode changes the load path or the schema - answers stay keyed by
# chem_interaction_id, one np_evidence row per pair.

CURATED = ("gtopdb_interactions", "cellphonedb_interactions", "cellchatdb_human",
           "mebocost_human", "mrclinkdb_human_mr")

UNANSWERED = """
    FROM chemical_interactions ci
    JOIN chemical_entities ch ON ch.chemical_id = ci.chemical_id
    JOIN entities e ON e.entity_id = ci.receptor_entity_id
   WHERE ci.species = 'human'
     AND NOT EXISTS (SELECT 1 FROM np_evidence ne
                      WHERE ne.chem_interaction_id = ci.chem_interaction_id
                        AND ne.verification_outcome IN ('pending','confirmed'))"""

COLS = """ci.chem_interaction_id, ci.n_sources, ci.sources, ci.receptor_type,
          ch.chemical_id, ch.name, ch.hmdb_id, ch.pubchem_cid, ch.chebi_id, ch.synonyms,
          e.symbol receptor,
          EXISTS(SELECT 1 FROM np_evidence ne
                  WHERE ne.chem_interaction_id = ci.chem_interaction_id) had_citation"""


def family(symbol: str) -> str:
    """Receptor family of a gene symbol, by stripping one trailing subtype token.

    RXRA->RXR, HTR1A->HTR, GABRA1->GABRA, CHRNA3->CHRNA, GRIN2A->GRIN, PDE3A->PDE.
    Crude, but it groups the symbols that comparative pharmacology papers group.
    """
    m = re.match(r"^(.+?)(\d+[A-Z]?|[A-Z])$", symbol)
    return m.group(1) if m else symbol


# "GPR" and "CD" are not families: GPR119/GPR12/GPR135 and CD1A/CD2/CD207 are arbitrary
# serial numbers, so a group of them shares no pharmacology and no single comparative paper
# covers it. (Those assertions are still suspect - anandamide alone is paired with 30 orphan
# GPRs - but they have to be judged one at a time, not as a cluster.)
NOT_A_FAMILY = ("GPR", "CD")


def clusters(c, n: int) -> list[list]:
    """N biggest unanswered (chemical, receptor family) groups of >=4 pairs."""
    done = searched()
    rows = [r for r in c.execute(f"SELECT {COLS} {UNANSWERED}").fetchall()
            if r["chem_interaction_id"] not in done]
    g: dict = {}
    for r in rows:
        g.setdefault((r["chemical_id"], family(r["receptor"])), []).append(r)
    big = [v for k, v in g.items() if len(v) >= 4 and k[1] not in NOT_A_FAMILY]
    big.sort(key=lambda v: (-len(v), v[0]["name"]))
    for v in big:
        v.sort(key=lambda r: r["receptor"])
    return big[:n]


def targets_multi(c, n: int) -> list:
    """Unanswered pairs with >=2 sources, those with a curated source first."""
    done = searched()
    rows = [r for r in c.execute(f"SELECT {COLS} {UNANSWERED}").fetchall()
            if r["chem_interaction_id"] not in done]
    out = []
    for r in rows:
        try:
            srcs = set(json.loads(r["sources"] or "[]"))
        except Exception:
            srcs = set()
        if len(srcs) >= 2:
            out.append((0 if srcs & set(CURATED) else 1, -r["n_sources"], r["name"],
                        r["receptor"], r))
    out.sort(key=lambda t: t[:4])
    return [t[4] for t in out[:n]]


def targets_ranked(c, n: int) -> list:
    """Unanswered pairs in work/np10_rank.tsv order (fitted likelihood of finding primary evidence)."""
    done = searched()
    order = []
    for ln in open(WORK / "np10_rank.tsv").read().splitlines()[1:]:
        cid = ln.split("\t")[2]
        if cid not in done:
            order.append(cid)
    rows = {r["chem_interaction_id"]: r for r in c.execute(f"SELECT {COLS} {UNANSWERED}").fetchall()}
    return [rows[i] for i in order if i in rows][:n]


def identifier(r) -> tuple[str, str]:
    for col, lab in (("hmdb_id", "HMDB"), ("pubchem_cid", "PubChem CID"), ("chebi_id", "ChEBI")):
        if r[col]:
            return lab, r[col]
    return "name only", r["name"]


def fill_prompt(r, syns: list[str]) -> str:
    txt = PROMPT.read_text()
    user = txt.split("user_content:", 1)[1].strip()
    it, iv = identifier(r)
    return (user.replace("{chemical_name}", r["name"]).replace("{identifier_type}", it)
                .replace("{identifier}", iv).replace("{synonyms}", "; ".join(syns[1:]) or "none recorded")
                .replace("{receptor_gene_symbol}", r["receptor"]))


# ---------------------------------------------------------------- emit

def emit(n: int, picker=None) -> int:
    c = connect()
    cache = json.loads(HGNC_CACHE.read_text()) if HGNC_CACHE.exists() else {}
    rows = (picker or targets)(c, n)
    WORK.mkdir(exist_ok=True)
    out = ["# Chemical-half literature search (np10)", "",
           f"{len(rows)} pairs. For each, answer the filled prompt below under the rules in "
           "`non_peptides_prompt.txt` (system_content). The PubMed hits are a starting point, "
           "not a limit - a paper found elsewhere may be reported by PMID and is fetched on load.",
           ""]
    log_lines = []
    for r in rows:
        syns = chem_synonyms(r)
        h = hgnc(r["receptor"], cache)
        rterms = receptor_terms(r["receptor"], h)
        rq = tiab(rterms)
        if sc := subunit_clause(r["receptor"]):
            rq = f"({rq[1:-1]} OR {sc})"
        term = f"{tiab(syns[:N_SYN])} AND {rq} NOT review[pt]"
        total, ids = esearch(term)
        if total > 2 * RETMAX:
            # Too many to read: keep the pair's papers that report an assay, which
            # is what the prompt asks for. Both queries are logged.
            term = f"{term[:-len(' NOT review[pt]')]} AND {ASSAY} NOT review[pt]"
            t2, ids2 = esearch(term)
            if ids2:
                total, ids = t2, ids2
        recs = efetch(ids)
        store_papers(c, recs)
        c.commit()
        by = {x["pmid"]: x for x in recs}
        log_lines.append(json.dumps({"chem_interaction_id": r["chem_interaction_id"],
                                     "chemical": r["name"], "receptor": r["receptor"],
                                     "query": term, "n_total": total, "hits": ids,
                                     "searched_at": now()}))
        it, iv = identifier(r)
        out += [f"## `{r['chem_interaction_id']}` — {r['name']} → {r['receptor']}", "",
                f"{it}: {iv} · receptor_type {r['receptor_type'] or '?'} · "
                f"{r['n_sources']} source(s) {r['sources']} · "
                f"{'had upstream citations (all closed unread)' if r['had_citation'] else 'no upstream citation'}",
                f"receptor: {h['name'] or '?'} · aliases {', '.join(h['alias']) or '-'}", "",
                "**Filled prompt:**", "", "> " + fill_prompt(r, syns).replace("\n", "\n> "), "",
                f"**PubMed:** {total} hit(s); showing {len(ids)}. Query: `{term}`", ""]
        for pmid in ids:
            p = by.get(pmid, {})
            flag = " **[REVIEW]**" if p.get("is_review") else ""
            flag += " **[RETRACTED]**" if p.get("is_retracted") else ""
            out += [f"### PMID {pmid}{flag} — {p.get('title')}", f"*{p.get('journal')} {p.get('year')}*", "",
                    p.get("abstract") or "*no abstract*", ""]
        out += ["---", ""]
    with LOGFILE.open("a") as fh:
        fh.write("".join(l + "\n" for l in log_lines))
    HGNC_CACHE.write_text(json.dumps(cache, indent=0, sort_keys=True))
    WORKFILE.write_text("\n".join(out))
    log("np10", f"wrote {WORKFILE} — {len(rows)} pairs; logged to {LOGFILE.name}")
    return 0


# ---------------------------------------------------------------- emit clusters

def emit_clusters(n: int) -> int:
    """One worklist entry per (chemical, receptor family) cluster of >=4 unanswered pairs.

    The PubMed query asks for the chemical against the WHOLE family at once, because the
    paper that settles a cluster is a comparative one (a selectivity panel, a stereoisomer
    series, a crystallographic survey). The answer file is unchanged: one record per
    chem_interaction_id, so each pair still carries its own citation and verdict.
    """
    c = connect()
    cache = json.loads(HGNC_CACHE.read_text()) if HGNC_CACHE.exists() else {}
    groups = clusters(c, n)
    WORK.mkdir(exist_ok=True)
    npairs = sum(len(g) for g in groups)
    out = ["# Chemical-half literature search - CLUSTER mode (np10)", "",
           f"{len(groups)} cluster(s), {npairs} pairs. Each block is ONE chemical against "
           "several members of one receptor family - the shape produced by matching compound "
           "class rather than measured binding. Look for the ONE comparative paper (selectivity "
           "panel, stereoisomer series, structural survey) that covers the whole family, then "
           "answer EVERY pair id listed, under the rules in `non_peptides_prompt.txt`. "
           "A pair the comparative paper tested and found inactive is a REFUTES (prompt rule 7) "
           "- that is the point of this mode. Answers are still keyed by pair id.", ""]
    log_lines = []
    for g in groups:
        r0 = g[0]
        syns = chem_synonyms(r0)
        fam = family(r0["receptor"])
        rterms: list[str] = []
        for r in g:
            h = hgnc(r["receptor"], cache)
            rterms += receptor_terms(r["receptor"], h)
        seen, uniq = set(), []
        for t in rterms:
            if t.lower() not in seen:
                seen.add(t.lower()); uniq.append(t)
        term = f"{tiab(syns[:N_SYN])} AND {tiab(uniq[:18])} NOT review[pt]"
        total, ids = esearch(term)
        if total > 2 * RETMAX:
            term = f"{term[:-len(' NOT review[pt]')]} AND {ASSAY} NOT review[pt]"
            t2, ids2 = esearch(term)
            if ids2:
                total, ids = t2, ids2
        recs = efetch(ids)
        store_papers(c, recs); c.commit()
        by = {x["pmid"]: x for x in recs}
        for r in g:
            log_lines.append(json.dumps({"chem_interaction_id": r["chem_interaction_id"],
                                         "chemical": r["name"], "receptor": r["receptor"],
                                         "query": term, "n_total": total, "hits": ids,
                                         "cluster": f"{r0['name']}/{fam}",
                                         "searched_at": now()}))
        it, iv = identifier(r0)
        out += [f"## CLUSTER `{r0['name']}` x {fam}* - {len(g)} pairs", "",
                f"{it}: {iv} · sources across the cluster: "
                f"{sorted({s for r in g for s in json.loads(r['sources'] or '[]')})}", ""]
        out += ["| pair id | receptor | n_sources |", "|---|---|---|"]
        out += [f"| `{r['chem_interaction_id']}` | {r['receptor']} | {r['n_sources']} |" for r in g]
        out += ["", "**Filled prompt (answer once per receptor above):**", "",
                "> " + fill_prompt(r0, syns).replace(
                    r0["receptor"], ", ".join(r["receptor"] for r in g)).replace("\n", "\n> "), "",
                f"**PubMed:** {total} hit(s); showing {len(ids)}. Query: `{term}`", ""]
        for pmid in ids:
            p = by.get(pmid, {})
            flag = " **[REVIEW]**" if p.get("is_review") else ""
            flag += " **[RETRACTED]**" if p.get("is_retracted") else ""
            out += [f"### PMID {pmid}{flag} — {p.get('title')}", f"*{p.get('journal')} {p.get('year')}*",
                    "", p.get("abstract") or "*no abstract*", ""]
        out += ["---", ""]
    with LOGFILE.open("a") as fh:
        fh.write("".join(l + "\n" for l in log_lines))
    HGNC_CACHE.write_text(json.dumps(cache, indent=0, sort_keys=True))
    WORKFILE.write_text("\n".join(out))
    log("np10", f"wrote {WORKFILE} — {len(groups)} cluster(s), {npairs} pairs")
    return 0


# ---------------------------------------------------------------- load

def load(path: str, commit: bool) -> int:
    c = connect()
    answers = json.loads(open(path).read())
    need = sorted({str(rec["PMID"]) for v in answers.values() if isinstance(v, list) for rec in v}
                  - {r[0] for r in c.execute("SELECT pmid FROM papers WHERE abstract IS NOT NULL")})
    if need:
        store_papers(c, efetch(need))
    stats, gate_log, rows, noev = {}, [], [], []

    def bump(k):
        stats[k] = stats.get(k, 0) + 1

    for cid, ans in answers.items():
        if not c.execute("SELECT 1 FROM chemical_interactions WHERE chem_interaction_id=?", (cid,)).fetchone():
            bump("unknown_pair"); gate_log.append((cid, "no such chem_interaction_id")); continue
        if ans == "No evidence" or ans == [] or (isinstance(ans, dict) and ans.get("result")):
            # 'No evidence' was previously only a counter, so the REASON was lost - and the
            # reason is the useful part. Wrong molecule (cystine for cysteine, cholecalciferol
            # for calcitriol), cargo mistaken for ligand (iron/TFR1, B12/CUBN), partner
            # mistaken for target (9-cis-RA at every RXR heterodimer partner), name collision
            # (L-arginine/AVPR2), human pseudogene (GLRA4), or simply nobody has run the
            # experiment - these need very different follow-up, and a curator cannot tell them
            # apart from a bare count. Appended as TSV; nothing is written to np_evidence,
            # because 'no paper found' is not evidence and must never look like a verdict.
            reason = ans.get("reason", "") if isinstance(ans, dict) else ""
            noev.append((cid, reason))
            bump("pair_no_evidence"); continue
        for rec in ans:
            pmid = str(rec.get("PMID"))
            p = c.execute("SELECT title, abstract, is_review, is_retracted, full_text "
                          "FROM papers WHERE pmid=?", (pmid,)).fetchone()
            abs_src = squash((p["title"] or "") + " " + (p["abstract"] or "")) if p else ""
            # Full text is a second GENUINE source (PLAN.md invariant 10), but only when
            # it is the publisher's text and not a curator worksheet. Same content-based
            # predicate as cur_integrity_check.is_curator_note / cur_fetch_fulltext.
            ft_src = squash((p["title"] or "") + " " + p["full_text"]) if p and real_fulltext(p["full_text"]) else ""
            # A record may say where it read the sentence; otherwise the abstract wins.
            want_ft = str(rec.get("EVIDENCE_SOURCE", "")).lower() in ("full_text", "fulltext")
            depth, src = "abstract_only", abs_src
            if want_ft or (rec.get("VERBATIM_SENTENCE")
                           and squash(rec["VERBATIM_SENTENCE"]) not in abs_src):
                depth, src = "full_text", ft_src
            why = None
            if not p or not (p["abstract"] or "").strip():
                why = "gate1_no_abstract"
            elif p["is_review"] or p["is_retracted"]:
                why = "gate2_review_or_retracted"
            elif depth == "full_text" and not ft_src:
                # Nothing to check against: fetch it first (cur_fetch_fulltext.py <pmid>).
                why = "gate3_needs_fulltext_not_on_file"
            elif squash(rec.get("VERBATIM_SENTENCE")) not in src or not rec.get("VERBATIM_SENTENCE"):
                why = "gate3_sentence_not_verbatim"
            elif squash(rec.get("MOLECULE_AS_NAMED_IN_PAPER")) not in src:
                why = "gate4_molecule_not_in_text"
            elif squash(rec.get("RECEPTOR_AS_NAMED_IN_PAPER")) not in src:
                why = "gate4_receptor_not_in_text"
            elif rec.get("EXPERIMENT_TYPE") not in VALID_CLASS or rec.get("RESULT") not in VALID_RESULT:
                why = "gate5_vocabulary"
            if why:
                bump(why); gate_log.append((f"{cid}/{pmid}", why)); continue
            note = ("np10 literature search (non_peptides_prompt.txt); gates passed, written "
                    "as confirmed without a curator read (owner instruction 2026-10-09). "
                    "Prompt answer: " + json.dumps(rec, ensure_ascii=False))
            rows.append((f"NPEV_{cid}_{pmid}", cid, pmid, SOURCE_ID, depth,
                         rec["VERBATIM_SENTENCE"], rec["EXPERIMENT_TYPE"],
                         rec.get("DIRECTION"), note))
            bump(f"lead_{rec['RESULT']}")
            bump(f"depth_{depth}")
    ins = 0
    for r in rows:
        cur = c.execute("""INSERT OR IGNORE INTO np_evidence
            (np_evidence_id, chem_interaction_id, pmid, source_id, verification_outcome,
             evidence_depth, evidence_sentence, assay_class, effect_direction,
             extraction_method, notes, verified_at) VALUES (?,?,?,?,'confirmed',?,?,?,?,
             'llm_litsearch_np_prompt',?,datetime('now'))""", r)
        ins += cur.rowcount
    if noev and commit:
        new = not NOEVFILE.exists()
        with NOEVFILE.open("a") as fh:
            if new:
                fh.write("chem_interaction_id\tchemical\treceptor\treason\tlogged_at\n")
            for cid, why in noev:
                r = c.execute("""SELECT ch.name, e.symbol FROM chemical_interactions ci
                                 JOIN chemical_entities ch ON ch.chemical_id = ci.chemical_id
                                 JOIN entities e ON e.entity_id = ci.receptor_entity_id
                                WHERE ci.chem_interaction_id=?""", (cid,)).fetchone()
                nm, sym = (r["name"], r["symbol"]) if r else ("?", "?")
                fh.write(f"{cid}\t{nm}\t{sym}\t{why or 'no reason recorded'}\t{now()}\n")
        log("np10", f"  wrote {len(noev)} no-evidence reason(s) to {NOEVFILE.name}")
    for k, v in sorted(stats.items()):
        log("np10", f"  {k:32s} {v}")
    for k, why in gate_log[:30]:
        log("np10", f"  GATE {k}: {why}")
    log("np10", f"inserted {ins} pending lead row(s) ({len(rows) - ins} already present)")
    (c.commit if commit else c.rollback)()
    log("np10", "COMMITTED" if commit else "DRY RUN - rerun with --commit to write")
    return 0


def status() -> int:
    c = connect()
    log("np10", f"pairs searched: {len(searched())}")
    for r in c.execute("""SELECT verification_outcome o, count(*) n, count(DISTINCT chem_interaction_id) p
                            FROM np_evidence WHERE source_id=? GROUP BY 1""", (SOURCE_ID,)):
        log("np10", f"  {r['o']:24s} {r['n']} rows over {r['p']} pairs")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--emit", type=int)
    ap.add_argument("--emit-clusters", type=int,
                    help="N biggest (chemical, receptor family) clusters of >=4 unanswered pairs")
    ap.add_argument("--emit-ranked", type=int, help="N next pairs in work/np10_rank.tsv order")
    ap.add_argument("--emit-multi", type=int,
                    help="N unanswered pairs with >=2 sources, curated-source ones first")
    ap.add_argument("--load")
    ap.add_argument("--commit", action="store_true")
    ap.add_argument("--status", action="store_true")
    a = ap.parse_args()
    if a.emit:
        sys.exit(emit(a.emit))
    if a.emit_clusters:
        sys.exit(emit_clusters(a.emit_clusters))
    if a.emit_ranked:
        sys.exit(emit(a.emit_ranked, targets_ranked))
    if a.emit_multi:
        sys.exit(emit(a.emit_multi, targets_multi))
    sys.exit(load(a.load, a.commit) if a.load else status())
