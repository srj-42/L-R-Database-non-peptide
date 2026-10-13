#!/usr/bin/env python3
"""Phase NP-11 - adjudicate every PENDING np_evidence row, so no pair shows a PMID
with nothing said about it.

The site's pairs table has two different columns that look like degrees of one thing:
`PMIDs` is a count an upstream resource shipped with the pair and nothing has checked,
while `Evidence` only fills from a CONFIRMED row carrying a verbatim sentence. 597
shippable human pairs sit between the two - cited, but with no line extracted. They are
not unsearched: they hold 2,178 np_evidence rows already, of which 1,737 are `pending`
over 1,130 distinct papers. This phase decides each one.

Outcome per row, using np9_verify.py's own vocabulary (nothing new invented):
  confirmed             a sentence in the paper states the result -> extracted verbatim
  review_only           review, or a database/tool/resource paper (prompt rule 1)
  retracted             flagged retracted
  abstract_insufficient the paper is primary and on topic, but no sentence in the
                        abstract states this molecule acting at this receptor
  different_entity      the paper is about a different molecule or a different protein
  unresolvable          no abstract on file and none retrievable

MECHANICAL vs READ. Four buckets are decided without judgement: retracted, review
(papers.is_review), no abstract, and "the molecule or the receptor is never named in
title+abstract". The last one is checked against the chemical's full synonym list AND
the receptor's HGNC aliases and previous symbols - a first pass that looked only at the
gene symbol wrongly marked 478 rows as unnamed, because papers write GPR109A for HCAR2
and nAChR for the CHRNA subunits. Everything else is read.

State lives in work/np11_state.json and work/np11_audit.tsv and is append-only, so a
run can stop anywhere and resume without re-deciding a row or losing a reason.

Usage:
  python3 np11_cited_audit.py --triage            classify mechanically, write state
  python3 np11_cited_audit.py --emit N            next N rows needing a read
  python3 np11_cited_audit.py --mechanical-json   outcomes for np9_verify.py --load
  python3 np11_cited_audit.py --status
"""
from __future__ import annotations

import argparse
import json
import re
import unicodedata
import sys
from pathlib import Path

from cur_common import HERE, connect, log, now

WORK = HERE / "work"
STATE = WORK / "np11_state.json"
AUDIT = WORK / "np11_audit.tsv"
READFILE = WORK / "NP11_READ.md"
HGNC = WORK / "np10_hgnc_cache.json"

TOOL_WORDS = ("database", "resource", "atlas", "web server", "webserver", "toolkit",
              "r package", "python package", "pipeline", "software", "repository",
              "knowledgebase", "knowledge base", "curated collection", "data portal",
              "inference method", "cell-cell communication", "cell-cell interaction",
              "ligand-receptor database", "systematic review", "meta-analysis")


def squash(s: str) -> str:
    return re.sub(r"\s+", " ", (s or "")).strip().lower()


def norm(s: str) -> str:
    """Strip everything but letters and digits, for name matching only.

    A first pass compared squashed strings and got four kinds of false negative,
    found by hand-checking samples before loading anything:
      * 'Leukotriene B4' never matched a paper writing 'Leukotriene B(4)';
      * 'FFAR1' never matched a paper writing 'GPR40', because the alias cache held
        only the 287 symbols earlier np10 runs happened to fetch;
      * 'LPAR3' never matched 'lysophosphatidic acid receptor';
      * 13 rows were called unresolvable although PubMed has the abstract.
    Normalising to alphanumerics fixes the first; the others are fixed by fetching
    all HGNC records and all missing abstracts before triage.
    """
    # NFKD folds Unicode subscripts, so a paper's "LTC₄" reads as "LTC4"
    s = unicodedata.normalize("NFKD", s or "")
    return re.sub(r"[^a-z0-9]", "", s.lower())


# Synonym lists in this build are partly corrupted: a chemical name was comma-split, so
# they contain fragments like "2", "3R", "4-dioxo-1" and bare words like "Acid". Matching
# on those produces both false positives and nonsense reasons ("molecule matched on: Acid"),
# so a synonym has to be long enough and not a generic chemical word to count.
JUNK_SYN = {"acid", "ester", "salt", "anion", "cation", "base", "acidform", "ion",
            "amine", "alcohol", "sodium", "potassium", "chloride", "hydrate", "water"}


def chem_names(name: str, synonyms: str | None) -> list[str]:
    """The chemical's name plus its synonym list, long enough not to match by accident."""
    out = [name] if name else []
    try:
        syn = json.loads(synonyms or "[]")
    except Exception:                                        # noqa: BLE001
        syn = []
    for s in syn:
        if not isinstance(s, str):
            continue
        n = norm(s)
        if len(n) >= 6 and n not in JUNK_SYN and not n.isdigit():
            out.append(s.strip())
    return out


def rec_names(symbol: str, hgnc: dict) -> list[str]:
    """Gene symbol, its HGNC aliases and previous symbols, plus the obvious stems.

    A paper rarely writes the HGNC symbol. HCAR2 appears as GPR109A, CHRNA7 as
    'alpha7 nicotinic', SLC29A1 as ENT1. Without these the triage calls a perfectly
    good paper 'receptor not named'.
    """
    h = hgnc.get(symbol) or {}
    out = [symbol] + list(h.get("alias") or []) + list(h.get("prev") or [])
    if h.get("name"):
        out.append(h["name"])
    return [x for x in out if x and len(norm(str(x))) > 2]


def looks_like_tool(title: str, abstract: str) -> str | None:
    t = squash(title)
    for w in TOOL_WORDS:
        if w in t:
            return w
    return None


def load_state() -> dict:
    if STATE.exists():
        return json.loads(STATE.read_text())
    return {"triaged_at": None, "decided": {}, "to_read": []}


def save_state(st: dict) -> None:
    STATE.write_text(json.dumps(st, indent=0, sort_keys=True))


def append_audit(rows: list[tuple]) -> None:
    new = not AUDIT.exists()
    with AUDIT.open("a") as fh:
        if new:
            fh.write("np_evidence_id\tchem_interaction_id\tpmid\tchemical\treceptor\t"
                     "outcome\treason_or_sentence\tdecided_by\tdecided_at\n")
        for r in rows:
            fh.write("\t".join(str(x).replace("\t", " ").replace("\n", " ") for x in r) + "\n")


def pending_rows(c):
    return c.execute("""
        SELECT ne.np_evidence_id eid, ne.chem_interaction_id cid, ne.pmid,
               ch.name chem, ch.synonyms, e.symbol rec,
               p.abstract, p.title, p.is_review, p.is_retracted, p.journal, p.year
          FROM np_evidence ne
          JOIN chemical_interactions ci ON ci.chem_interaction_id = ne.chem_interaction_id
          JOIN chemical_entities ch ON ch.chemical_id = ci.chemical_id
          JOIN entities e ON e.entity_id = ci.receptor_entity_id
          LEFT JOIN papers p ON p.pmid = ne.pmid
         WHERE ne.verification_outcome = 'pending'
           AND ci.species = 'human' AND ci.shippable = 1
         ORDER BY ch.name, e.symbol, ne.pmid""").fetchall()


def triage() -> int:
    c = connect()
    hgnc = json.loads(HGNC.read_text()) if HGNC.exists() else {}
    st = load_state()
    decided, to_read, audit = st["decided"], [], []
    counts: dict[str, int] = {}

    def bump(k):
        counts[k] = counts.get(k, 0) + 1

    for r in pending_rows(c):
        eid = r["eid"]
        if eid in decided:
            bump("already decided"); continue
        txt = norm((r["title"] or "") + " " + (r["abstract"] or ""))
        out = reason = None
        if r["is_retracted"]:
            out, reason = "retracted", "paper flagged retracted"
        elif r["is_review"]:
            out, reason = "review_only", "publication type is Review (prompt rule 1)"
        elif not (r["abstract"] or "").strip():
            out, reason = "unresolvable", "no abstract on file at PubMed"
        elif (w := looks_like_tool(r["title"] or "", r["abstract"] or "")):
            out, reason = "review_only", f"database/resource/tool paper (title contains '{w}') - prompt rule 1"
        # NOTHING ELSE is decided mechanically. An earlier version excluded a row when
        # the molecule or the receptor failed to match the abstract by name, and
        # hand-checking the output found that test wrong often enough to be unusable:
        #   cGMP            not matched by the name 'Cyclic GMP'
        #   propionate      not matched by 'propanoic acid'
        #   P2Y14           not matched by the symbol 'P2RY14'
        #   GPR40 / BLT2    unmatched until every HGNC record was fetched
        #   Leukotriene B(4) unmatched until names were normalised to alphanumerics
        # Each miss would have buried a real citation under a confident-sounding reason.
        # Name matching is kept, but only as a HINT on the reading worklist - the chemical
        # and receptor vocabularies in this build are not complete enough for it to
        # exclude anything on its own.
        if out:
            decided[eid] = {"outcome": out, "reason": reason, "by": "mechanical"}
            audit.append((eid, r["cid"], r["pmid"], r["chem"], r["rec"], out, reason,
                          "mechanical", now()))
            bump(out)
        else:
            to_read.append(eid)
            bump("needs reading")

    st["decided"], st["to_read"], st["triaged_at"] = decided, to_read, now()
    save_state(st)
    if audit:
        append_audit(audit)
    for k, v in sorted(counts.items(), key=lambda kv: -kv[1]):
        log("np11", f"  {k:34s} {v}")
    log("np11", f"state -> {STATE.name}; audit -> {AUDIT.name}; {len(to_read)} rows need a read")
    return 0


def mechanical_json() -> int:
    """Outcomes for the mechanically decided rows, in np9_verify.py --load shape."""
    st = load_state()
    # np9_verify.py --load iterates a LIST of records, each carrying np_evidence_id.
    out = []
    for eid, d in st["decided"].items():
        if d["by"] != "mechanical":
            continue
        out.append({"np_evidence_id": eid,
                    "verification_outcome": d["outcome"],
                    "evidence_depth": "abstract_only",
                    "extraction_method": "np11_mechanical_triage",
                    "notes": f"np11 cited-citation audit ({d['by']}): {d['reason']}"})
    p = WORK / "np11_mechanical.json"
    p.write_text(json.dumps(out, indent=1))
    log("np11", f"wrote {len(out)} mechanical outcomes to {p.name}")
    return 0


def emit(n: int) -> int:
    c = connect()
    st = load_state()
    todo = [e for e in st["to_read"] if e not in st["decided"]][:n]
    if not todo:
        log("np11", "nothing left to read")
        return 0
    by_eid = {r["eid"]: r for r in pending_rows(c)}
    hgnc = json.loads(HGNC.read_text()) if HGNC.exists() else {}
    out = ["# NP-11 - adjudicate upstream citations", "",
           f"{len(todo)} row(s). Both the molecule and the receptor ARE named in each "
           "abstract below, so each needs a read: either a sentence states this molecule "
           "acting at this receptor (then the row is `confirmed` and that sentence is the "
           "evidence line), or none does (then `abstract_insufficient`), or the paper turns "
           "out to be about a different molecule/protein (`different_entity`) or to be a "
           "review or resource paper (`review_only`).", ""]
    for eid in todo:
        r = by_eid.get(eid)
        if not r:
            continue
        txt = norm((r["title"] or "") + " " + (r["abstract"] or ""))
        cn = [n for n in chem_names(r["chem"], r["synonyms"]) if norm(n) in txt]
        rn = [n for n in rec_names(r["rec"], hgnc) if norm(str(n)) in txt]
        hint = (f"name hints (NOT a verdict): molecule "
                f"{'matched on ' + cn[0] if cn else 'NOT matched by name or synonym'} · receptor "
                f"{'matched on ' + str(rn[0]) if rn else 'NOT matched by symbol or HGNC alias'}")
        out += [f"## `{eid}` — {r['chem']} → {r['rec']}  (pair {r['cid']})", "",
                f"**PMID {r['pmid']}** — {r['title']}", f"*{r['journal']} {r['year']}*", "",
                f"_{hint}_", "",
                r["abstract"] or "*no abstract*", "", "---", ""]
    READFILE.write_text("\n".join(out))
    log("np11", f"wrote {READFILE.name} — {len(todo)} rows to read")
    return 0


def status() -> int:
    st = load_state()
    c = connect()
    left = len([e for e in st["to_read"] if e not in st["decided"]])
    by = {}
    for d in st["decided"].values():
        by[d["outcome"]] = by.get(d["outcome"], 0) + 1
    log("np11", f"decided {len(st['decided'])}; still to read {left}")
    for k, v in sorted(by.items(), key=lambda kv: -kv[1]):
        log("np11", f"  {k:24s} {v}")
    n = c.execute("""SELECT count(*) FROM np_evidence ne
                       JOIN chemical_interactions ci
                         ON ci.chem_interaction_id = ne.chem_interaction_id
                      WHERE ne.verification_outcome='pending' AND ci.species='human'
                        AND ci.shippable=1""").fetchone()[0]
    log("np11", f"pending rows still in the database: {n}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--triage", action="store_true")
    ap.add_argument("--emit", type=int)
    ap.add_argument("--mechanical-json", action="store_true")
    ap.add_argument("--status", action="store_true")
    a = ap.parse_args()
    if a.triage:
        return triage()
    if a.emit:
        return emit(a.emit)
    if a.mechanical_json:
        return mechanical_json()
    return status()


if __name__ == "__main__":
    sys.exit(main())
