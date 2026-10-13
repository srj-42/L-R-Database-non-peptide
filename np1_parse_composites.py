#!/usr/bin/env python3
"""Phase NP-1 - decompose the CellChatDB / CellPhoneDB composite ligand strings.

Both resources bundle the ligand and its enzymes into one token:

    CellChatDB    5-HT-TPH1_SLC18A1          chemical 5-HT, via TPH1, exported by SLC18A1
    CellPhoneDB   5-HT_byTPH1_and_SLC18A1    same claim, different spelling

That is a denormalised (metabolite x synthase x transporter) cross-product. This
phase splits it back apart so the chemical goes to `chemical_entities` and the
enzymes go to `metabolite_enzymes` - which is the whole point of keeping the two
relations separate (PROJECT_REPORT.md section 7).

The v1 prototype (src/phase6_nonpeptide.py) parsed these by matching against a
hardcoded list of known chemical names. That does not scale to CellPhoneDB's 120
tokens, so this inverts the test: gene symbols are the closed vocabulary we can
enumerate, so tokens are consumed from the RIGHT while they are known genes, and
whatever remains on the left is the chemical.

*** THIS SCRIPT WRITES NO DATABASE ROWS. ***
It emits np_composite_map.tsv for human review. Columns `manual_chem_name` and
`manual_hmdb_id` are yours to fill; np2 reads the reviewed file, not this one's
guesses. Rows with needs_review=1 are the ones that need you.

Usage: python3 np1_parse_composites.py
"""
from __future__ import annotations

import csv
import re
import sqlite3
import sys
from pathlib import Path

from cur_common import HERE, RAW, connect, log, read_csv

sys.path.insert(0, str(HERE.parent / "src"))

OUT = HERE / "np_composite_map.tsv"
ROOT = HERE.parent

# Abbreviation -> (name, chebi, pubchem_cid, class), carried over verbatim from
# src/phase6_nonpeptide.py where they were hand-verified against ChEBI/PubChem.
# CellChatDB writes neurotransmitters as abbreviations that no synonym list in
# any of the seven sources resolves, so this table stays hand-maintained.
CHEM_ABBREV = {
    "Glu":         ("L-glutamate",            "CHEBI:16015", "33032",   "neurotransmitter"),
    "GABA":        ("GABA",                   "CHEBI:16865", "119",     "neurotransmitter"),
    "5-HT":        ("serotonin",              "CHEBI:28790", "5202",    "neurotransmitter"),
    "Ach":         ("acetylcholine",          "CHEBI:15355", "187",     "neurotransmitter"),
    "DOPA":        ("dopamine",               "CHEBI:18243", "681",     "neurotransmitter"),
    "NE":          ("norepinephrine",         "CHEBI:18357", "439260",  "neurotransmitter"),
    "ADR":         ("epinephrine",            "CHEBI:28918", "5816",    "neurotransmitter"),
    "Gly":         ("glycine",                "CHEBI:15428", "750",     "neurotransmitter"),
    "Hist":        ("histamine",              "CHEBI:18295", "774",     "neurotransmitter"),
    "ADO":         ("adenosine",              "CHEBI:16335", "60961",   "metabolite"),
    "RA":          ("retinoic acid",          "CHEBI:15367", "444795",  "metabolite"),
    "PGE2":        ("prostaglandin E2",       "CHEBI:15551", "5280360", "lipid"),
    "Cholesterol": ("cholesterol",            "CHEBI:16113", "5997",    "steroid"),
    "E2":          ("estradiol",              "CHEBI:16469", "5757",    "steroid"),
    "Progesterone": ("progesterone",          "CHEBI:17026", "5994",    "steroid"),
    "Testosterone": ("testosterone",          "CHEBI:17347", "6013",    "steroid"),
    "DHEA":        ("dehydroepiandrosterone", "CHEBI:28689", "5881",    "steroid"),
}

# Tokens that look like gene symbols but are substrate names CellChatDB embeds
# in the route description, e.g. `5-HT-DOPA-DDC` = "via DOPA-decarboxylase".
NOT_A_GENE = {"DOPA", "GABA", "ADO", "NE", "ADR", "RA", "AA", "E2"}


def norm_chem(s: str) -> str:
    """Fold to a comparison key: CellPhoneDB strips punctuation and spaces."""
    return re.sub(r"[^a-z0-9]", "", str(s or "").lower())


def gene_vocabulary(c: sqlite3.Connection) -> set[str]:
    """Every symbol we have independent reason to believe is a human gene.

    Breadth matters here: the right-walk only consumes a token it recognises, so
    a missing symbol silently turns a parseable string into a review item. The
    first pass missed transporters like SLC6A12 and SLC29A3 because the peptide
    corpus never saw them, hence the metabolite-side vocabularies below.
    """
    vocab = {r[0].upper() for r in c.execute(
        "SELECT symbol FROM entities WHERE species='human' AND symbol IS NOT NULL") if r[0]}
    # The enzyme and transporter files name genes the peptide corpus never saw.
    for fn, col in (("mrclinkdb_human_enzyme.txt", "GENE_NAME"),
                    ("mrclinkdb_human_transporter.txt", "GENE_NAME")):
        p = RAW / fn
        if p.exists():
            vocab |= {r[col].strip().upper() for r in read_csv(p, "\t") if r.get(col)}
    p = RAW / "metaligand_LR_data_human.csv"
    if p.exists():
        for r in read_csv(p, ","):
            for col in ("Synthetic_genes", "Synthetic_genes_lv2",
                        "Transporter_genes", "Transporter_genes_in", "Receptor_genes"):
                vocab |= {g.strip().upper() for g in (r.get(col) or "").split(",") if g.strip()}
    p = RAW / "metalinksdb.db"
    if p.exists():
        with sqlite3.connect(p) as s:
            vocab |= {g.upper() for (g,) in s.execute(
                "SELECT gene_symbol FROM proteins WHERE gene_symbol IS NOT NULL") if g}
    p = ROOT / "data" / "raw" / "cellphonedb_genes.csv"
    if p.exists():
        vocab |= {r["hgnc_symbol"].strip().upper() for r in read_csv(p, ",")
                  if r.get("hgnc_symbol")}
    p = RAW / "gtopdb_interactions.csv"
    if p.exists():
        vocab |= {r["Target Gene Symbol"].strip().upper() for r in read_csv(p, ",")
                  if r.get("Target Gene Symbol")}
    return {g for g in vocab if g} - {n.upper() for n in NOT_A_GENE}


def chem_index() -> dict[str, tuple[str, str]]:
    """normalised name -> (canonical name, hmdb_id), from the HMDB-native sources."""
    idx: dict[str, tuple[str, str]] = {}

    def add(name: str, hmdb: str | None) -> None:
        k = norm_chem(name)
        if k and k not in idx:
            idx[k] = (name, hmdb or "")

    p = RAW / "mrclinkdb_human_mr.txt"
    if p.exists():
        for r in read_csv(p, "\t"):
            add(r.get("Metabolite name", ""), r.get("HMDB ID"))
    p = RAW / "metaligand_LR_data_human.csv"
    if p.exists():
        for r in read_csv(p, ","):
            hmdb = (r.get("HMDB_id") or "").strip()
            add(r.get("Compounds", ""), hmdb)
            for syn in (r.get("Synonyms") or "").split(","):
                add(syn.strip(), hmdb)
    p = RAW / "mebocost_human_met_sensor.tsv"          # internal-use source, mapping aid only
    if p.exists():
        for r in read_csv(p, "\t"):
            add(r.get("standard_metName", ""), r.get("HMDB_ID"))
    p = RAW / "metalinksdb.db"
    if p.exists():
        with sqlite3.connect(p) as s:
            for name, hmdb in s.execute("SELECT metabolite, hmdb FROM metabolites"):
                add(name or "", hmdb)
    # Abbreviations inherit the HMDB id of their canonical name where the
    # HMDB-keyed sources already carry it - otherwise 'ADR' resolves to
    # 'epinephrine' with an empty id even though epinephrine is indexed.
    #
    # These OVERWRITE whatever a source happened to index under the same key.
    # CHEM_ABBREV is hand-verified; source names are not. Without this, 'DOPA'
    # matched an indexed '(-)-Dopa' (HMDB0000181, L-DOPA, the PRECURSOR) instead
    # of dopamine (HMDB0000073). CellChatDB's DOPA rows have pathway_name
    # 'Dopamine' and receptors DRD1-DRD5, so the ligand is dopamine and L-DOPA is
    # the precursor - a distinction the enzyme table carries, not the ligand.
    for abbr, (name, _chebi, _cid, _cls) in CHEM_ABBREV.items():
        canon = idx.get(norm_chem(name))
        idx[norm_chem(abbr)] = canon if canon else (name, "")
        add(name, canon[1] if canon else "")
    return idx


def fused_pair(chem: str, idx: dict[str, tuple[str, str]]) -> tuple[str, str] | None:
    """Detect two chemical names concatenated into one token.

    CellPhoneDB ships `AldosteroneCorticosterone_byCYP11B2` and
    `SerotoninDopamin_byDDC` - one row asserting two different ligands. They
    cannot be imported as a single chemical and must be split by hand.
    """
    k = norm_chem(chem)
    if len(k) < 8 or k in idx:
        return None
    for i in range(4, len(k) - 3):
        a, b = idx.get(k[:i]), idx.get(k[i:])
        if a and b:
            return a[0], b[0]
    return None


def parse(token: str, vocab: set[str],
          idx: dict[str, tuple[str, str]]) -> tuple[str, list[str], str, str]:
    """`5-HT-TPH1_SLC18A1` -> ('5-HT', ['TPH1','SLC18A1'], 'right-walk', '').

    Consumes known gene symbols from the right. Whatever is left is the chemical,
    with its internal separators restored. The fourth return value is a precursor
    token where CellChatDB embedded one: `5-HT-DOPA-DDC` is "5-HT, made from DOPA
    by DDC", so DOPA is a precursor, not part of the ligand name.
    """
    raw = token.replace("_and_", "_").replace("_by", "_")
    parts = [p for p in re.split(r"[_-]", raw) if p]
    enzymes: list[str] = []
    while len(parts) > 1 and parts[-1].upper() in vocab:
        enzymes.insert(0, parts.pop())
    if not enzymes:
        return token, [], "no-enzyme-found", ""
    # Rebuild the chemical using the separators it originally had.
    head = token
    for e in reversed(enzymes):
        head = re.sub(r"[_-]+(?:and_|by)?" + re.escape(e) + r"$", "", head)
    head = re.sub(r"[_-]+$", "", head)

    # Substrate embedded in the route description: strip it only when BOTH the
    # head and the trailing token are independently recognisable chemicals.
    precursor = ""
    if "-" in head and norm_chem(head) not in idx:
        stem, _, tail = head.rpartition("-")
        if stem and norm_chem(stem) in idx and norm_chem(tail) in idx:
            head, precursor = stem, tail
    return (head or token), enzymes, "right-walk", precursor


def main() -> int:
    c = connect()
    vocab = gene_vocabulary(c)
    idx = chem_index()
    log("np1", f"gene vocabulary: {len(vocab)} symbols | chemical index: {len(idx)} names")

    rows: list[dict] = []

    # --- CellChatDB v2 -----------------------------------------------------
    try:
        from rdata import read_rdata
    except ImportError:
        log("np1", "ERROR: database/src/rdata.py not importable")
        return 1
    rda = ROOT / "data" / "raw" / "CellChatDB.human.rda"
    obj = read_rdata(rda)
    cc = obj.get("CellChatDB.human") or next(iter(obj.values()))
    it = cc["interaction"]
    counts: dict[tuple[str, str], int] = {}
    for i in range(len(it["interaction_name"])):
        if it["annotation"][i] != "Non-protein Signaling":
            continue
        tok = (it["ligand"][i] or "").strip()
        if tok:
            counts[("cellchatdb_human", tok)] = counts.get(("cellchatdb_human", tok), 0) + 1

    # --- CellPhoneDB v5 ----------------------------------------------------
    cp = read_csv(ROOT / "data" / "raw" / "cellphonedb_interactions.csv", ",")
    for r in cp:
        if str(r.get("is_ppi", "")).strip().lower() != "false":
            continue
        tok = (r.get("partner_a") or "").strip()
        if tok:
            counts[("cellphonedb_interactions", tok)] = counts.get(
                ("cellphonedb_interactions", tok), 0) + 1

    n_peptide = 0
    for (src, tok), n in sorted(counts.items()):
        chem, enzymes, method, precursor = parse(tok, vocab, idx)
        hit = idx.get(norm_chem(chem))
        resolved, hmdb = (hit if hit else ("", ""))
        review, reason = 0, ""
        # A "chemical" that is itself a gene symbol is usually a peptide ligand
        # filed under a non-protein annotation upstream. Not always, though:
        # `Cort` is corticosterone here but CORT is also the cortistatin gene, so
        # this flags for review rather than asserting.
        likely_peptide = int(not hit and chem.upper() in vocab)
        fused = fused_pair(chem, idx)
        if chem.endswith("_complex"):
            review, reason = 1, "protein complex, not a metabolite row - exclude"
        elif fused:
            review, reason = 1, f"two chemicals fused upstream: {fused[0]} + {fused[1]} - split into 2 rows"
        elif likely_peptide:
            n_peptide += 1
            review, reason = 1, "token is also a gene symbol - check whether peptide or chemical"
        elif method == "no-enzyme-found":
            review, reason = 1, "no gene symbol recognised in token"
        elif not hit:
            review, reason = 1, "chemical not found in any HMDB-keyed source"
        elif not hmdb:
            review, reason = 1, "matched by name only, no HMDB id"
        rows.append(dict(
            source_id=src, raw_ligand=tok, n_rows=n, chem_token=chem,
            enzyme_genes=";".join(enzymes), precursor_token=precursor,
            parse_method=method, resolved_name=resolved, hmdb_id=hmdb,
            likely_peptide=likely_peptide, needs_review=review, review_reason=reason,
            manual_chem_name="", manual_hmdb_id=""))

    cols = ["source_id", "raw_ligand", "n_rows", "chem_token", "enzyme_genes",
            "precursor_token", "parse_method", "resolved_name", "hmdb_id",
            "likely_peptide", "needs_review", "review_reason",
            "manual_chem_name", "manual_hmdb_id"]
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    for src in ("cellchatdb_human", "cellphonedb_interactions"):
        sub = [r for r in rows if r["source_id"] == src]
        ok = sum(1 for r in sub if not r["needs_review"])
        enz = sum(len(r["enzyme_genes"].split(";")) if r["enzyme_genes"] else 0 for r in sub)
        log("np1", f"{src:26s} tokens={len(sub):4d} resolved={ok:4d} "
                   f"needs_review={len(sub)-ok:4d} enzymes_extracted={enz}")
    n_prec = sum(1 for r in rows if r["precursor_token"])
    log("np1", f"precursors recovered: {n_prec} | likely peptide ligands misfiled upstream: {n_peptide}")
    log("np1", f"wrote {OUT.relative_to(HERE.parent)} ({len(rows)} rows) - NO database rows written")
    log("np1", "review needs_review=1 rows, fill manual_chem_name / manual_hmdb_id, then run np2")
    return 0


if __name__ == "__main__":
    sys.exit(main())
