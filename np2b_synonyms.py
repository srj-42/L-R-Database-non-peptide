#!/usr/bin/env python3
"""Phase NP-2b - apply curated chemical synonyms from np_synonym_overrides.tsv.

np2 builds chemical_entities.synonyms from whatever synonym lists the sources
carry. Those lists are incomplete - several chemicals arrive with none at all,
and no source carries the Greek-letter spellings the primary literature actually
uses (7α,25-dihydroxycholesterol). The np9 verification gate matches a curator's
quoted name against this list, so a missing synonym blocks a legitimate
confirmation.

Idempotent: re-run after any np2 rebuild.

Usage: python3 np2b_synonyms.py [--commit]
"""
from __future__ import annotations
import json, sys
from cur_common import HERE, connect, log

SRC = HERE / "np_synonym_overrides.tsv"


def main() -> int:
    commit = "--commit" in sys.argv
    c = connect()
    added = missing = 0
    for line in SRC.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        hmdb, syn = parts[0].strip(), parts[1].strip()
        row = c.execute("SELECT chemical_id, synonyms FROM chemical_entities WHERE hmdb_id=?",
                        (hmdb,)).fetchone()
        if row is None:
            log("np2b", f"  no chemical with hmdb_id {hmdb} - skipped")
            missing += 1
            continue
        syns = json.loads(row["synonyms"] or "[]")
        if syn not in syns:
            syns.append(syn)
            c.execute("UPDATE chemical_entities SET synonyms=? WHERE chemical_id=?",
                      (json.dumps(syns), row["chemical_id"]))
            added += 1
    if commit:
        c.commit()
    else:
        c.rollback()
    log("np2b", f"synonyms added: {added} | overrides with no matching chemical: {missing}")
    log("np2b", "COMMITTED" if commit else "DRY RUN - rerun with --commit to write")
    return 0


if __name__ == "__main__":
    sys.exit(main())
