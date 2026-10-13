#!/usr/bin/env python3
"""Phase NP-8 - build the evidence layer for the chemical half.

Until now the chemical build carried citations only as opaque strings inside
chem_source_records.raw_reference. This turns them into one row per (pair, paper)
in np_evidence, which is what a verification pass can actually work through.

It applies only the decisions that can be made WITHOUT reading, and applies them
as the peptide build does:

  * a RETRACTED paper closes the row as `retracted`
  * a REVIEW closes it as `review_only` - a review is a lead, not a receipt
    (PLAN.md invariant 2). This is a large fraction here: the non-peptide
    resources cite reviews far more freely than connectomeDB does.
  * everything else is left `pending` for a curator

Nothing is marked confirmed here. That requires reading, and reading happens in
np9_verify.py behind the same gates the peptide build uses.

np_evidence is deliberately a separate table from `evidence`. The two halves are
not comparable and must never be scored together.

Usage: python3 np8_evidence.py [--commit]
"""
from __future__ import annotations

import re
import sys
from collections import Counter, defaultdict

from cur_common import HERE, connect, log, now

NP_SCHEMA = HERE / "np_schema.sql"
PMID_RE = re.compile(r"(?<!\d)(\d{7,8})(?!\d)")
# MEBOCOST may not be redistributed, so its citations are recorded for internal
# use only; np6_export already withholds anything sourced from it.
NON_SHIPPABLE = {"mebocost_human", "mebocost_human_enzymes"}


def main() -> int:
    commit = "--commit" in sys.argv
    c = connect()
    c.executescript(NP_SCHEMA.read_text())
    c.commit()

    meta = {r["pmid"]: r for r in c.execute(
        "SELECT pmid, is_review, is_retracted, abstract, title FROM papers")}

    # (pair, pmid) -> set of sources citing it
    cited: dict[tuple[str, str], set[str]] = defaultdict(set)
    for r in c.execute("""SELECT chem_interaction_id, source_id, raw_reference
                            FROM chem_source_records
                           WHERE raw_reference IS NOT NULL
                             AND chem_interaction_id IS NOT NULL"""):
        for pmid in set(PMID_RE.findall(r["raw_reference"] or "")):
            cited[(r["chem_interaction_id"], pmid)].add(r["source_id"])

    rows, stats = [], Counter()
    for (cid, pmid), srcs in sorted(cited.items()):
        m = meta.get(pmid)
        if m is None:
            outcome, note = "unresolvable", "PMID not resolvable at PubMed."
        elif m["is_retracted"]:
            outcome, note = "retracted", "Supporting paper is retracted (PubMed publication type)."
        elif m["is_review"]:
            outcome, note = "review_only", ("Cited paper is a review (PubMed publication type). "
                                            "A review is a lead, not a receipt - this pair needs "
                                            "the primary paper the review is citing.")
        elif not (m["abstract"] or "").strip():
            outcome, note = "pending", "No abstract on file; nothing to read yet."
        else:
            outcome, note = "pending", None
        stats[outcome] += 1
        rows.append((f"NPEV_{cid}_{pmid}", cid, pmid, ";".join(sorted(srcs)),
                     outcome, note, now() if outcome != "pending" else None))

    c.execute("DELETE FROM np_evidence")
    c.executemany("""INSERT OR IGNORE INTO np_evidence
        (np_evidence_id, chem_interaction_id, pmid, source_id, verification_outcome,
         notes, verified_at) VALUES (?,?,?,?,?,?,?)""", rows)

    if commit:
        c.commit()
    else:
        c.rollback()

    n_pairs = len({r[1] for r in rows})
    log("np8", f"np_evidence rows: {len(rows)} over {n_pairs} pairs and "
               f"{len({r[2] for r in rows})} papers")
    for k, v in stats.most_common():
        log("np8", f"  {k:16s} {v}")
    ship = sum(1 for r in rows if not set(r[3].split(";")) <= NON_SHIPPABLE)
    log("np8", f"shippable rows: {ship} | MEBOCOST-only citations withheld: {len(rows) - ship}")
    log("np8", "COMMITTED" if commit else "DRY RUN - rerun with --commit to write")
    return 0


if __name__ == "__main__":
    sys.exit(main())
