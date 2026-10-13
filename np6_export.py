#!/usr/bin/env python3
"""Phase NP-6 - export the chemical half to data/release_v2/chemical_*.csv.

Companion files to the peptide release, never merged into it (decision 1). They
share a directory and nothing else: no chemical row appears in
interactions_scored.csv, and no chemical row carries a confidence_score
(decision 9). What ships in its place is provenance - n_sources, sources,
has_citation, n_pmids - so nothing here can be mistaken for a scored peptide row.

Since 2026-10-09 the chemical half ALSO ships chemical_evidence.csv: confirmed
primary citations with their verbatim sentence, produced by np10_litsearch.py.
These rows are machine-verified against the paper's own text (the load gates),
not read by a curator - extraction_method says which. They stay in their own
file and are never merged into evidence.csv, the peptide half's.

LICENCE ENFORCEMENT. MEBOCOST forbids redistributing its data, or any part of
it, modified or not. Two things follow, and both are enforced by assertion
rather than by care:

  * Rows sourced only from MEBOCOST are withheld entirely.
  * MEBOCOST is STRIPPED FROM THE PROVENANCE LISTS of rows that ship for other
    reasons, and n_sources is recomputed accordingly. A list asserting "MEBOCOST
    contains this pair" is arguably a part of the pair list, which is the thing
    the licence protects. Exported n_sources can therefore be lower than the
    internal value, which is the conservative direction.

The build fails rather than shipping if any assertion trips.

Usage: python3 np6_export.py [--out DIR] [--skip-integrity]
"""
from __future__ import annotations

import csv
import json
import re
import sqlite3
import subprocess
import sys
from datetime import date
from pathlib import Path

from cur_common import HERE, connect, log, read_csv

ROOT = HERE.parent
OUT = ROOT / "data" / "release_v2"
MAP = HERE / "np_composite_map.tsv"

# Must match sources.shippable = 0. Asserted against the database below.
NON_SHIPPABLE = {"mebocost_human", "mebocost_human_enzymes"}


def write(path: Path, header: list[str], rows: list[tuple]) -> int:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    log("np6", f"  {path.name:<34} {len(rows):>7} rows")
    return len(rows)


def main() -> int:
    out = OUT
    if "--out" in sys.argv:
        out = Path(sys.argv[sys.argv.index("--out") + 1])
    out.mkdir(parents=True, exist_ok=True)
    c = connect()

    # ---- gate: the standing integrity check ------------------------------
    # PLAN.md: cur_integrity_check.py must exit 0 before any release or export.
    if "--skip-integrity" not in sys.argv:
        r = subprocess.run([sys.executable, str(HERE / "cur_integrity_check.py")],
                           capture_output=True, text=True)
        if r.returncode != 0:
            log("np6", "ABORT: cur_integrity_check.py did not pass; refusing to export")
            print(r.stdout[-2000:])
            return 1
        log("np6", "gate: cur_integrity_check.py passed")

    # ---- assert the non-shippable set matches the database ---------------
    db_blocked = {r["source_id"] for r in
                  c.execute("SELECT source_id FROM sources WHERE shippable=0")}
    if db_blocked != NON_SHIPPABLE:
        log("np6", f"ABORT: sources.shippable=0 is {sorted(db_blocked)} but this script "
                   f"knows {sorted(NON_SHIPPABLE)}; reconcile before exporting")
        return 1

    def ship_sources(js: str) -> list[str]:
        return [s for s in json.loads(js or "[]") if s not in NON_SHIPPABLE]

    # ---- chemicals whose identity is still awaiting sign-off -------------
    provisional: set[str] = set()
    for r in read_csv(MAP, "\t"):
        if r.get("needs_review") == "1":
            hid = (r.get("manual_hmdb_id") or r.get("hmdb_id") or "").strip()
            for h in hid.split("|"):
                if h.strip():
                    provisional.add(h.strip())

    # ---- chemical_entities ------------------------------------------------
    chem_cols = ["chemical_id", "hmdb_id", "chebi_id", "pubchem_cid", "pubchem_sid",
                 "name", "synonyms", "smiles", "formula", "kingdom", "super_class",
                 "class", "chebi_super_class", "chebi_class", "metabolite_subclass",
                 "molecular_class", "is_neurotransmitter"]
    chem_rows, shipped_chem = [], set()
    for r in c.execute(f"SELECT {','.join(chem_cols)}, sources, shippable "
                       "FROM chemical_entities ORDER BY chemical_id"):
        if not r["shippable"]:
            continue
        srcs = ship_sources(r["sources"])
        if not srcs:
            continue
        shipped_chem.add(r["chemical_id"])
        chem_rows.append(tuple(r[k] for k in chem_cols) +
                         (len(srcs), json.dumps(srcs),
                          1 if (r["hmdb_id"] or "") in provisional else 0))
    counts = {}
    counts["chemical_entities"] = write(
        out / "chemical_entities.csv",
        chem_cols + ["n_sources", "sources", "provisional_identity"], chem_rows)

    # ---- chemical_interactions -------------------------------------------
    inter_rows, shipped_inter = [], set()
    for r in c.execute("""SELECT ci.chem_interaction_id, ci.chemical_id, ch.name chem_name,
                                 ch.hmdb_id, ci.receptor_entity_id, e.symbol receptor_symbol,
                                 e.uniprot receptor_uniprot, ci.species, ci.complex_id,
                                 ci.complex_role, ci.receptor_type, ci.sources,
                                 ci.has_citation, ci.n_pmids, ci.verification_status,
                                 ci.shippable
                            FROM chemical_interactions ci
                            JOIN chemical_entities ch ON ch.chemical_id = ci.chemical_id
                            JOIN entities e ON e.entity_id = ci.receptor_entity_id
                           ORDER BY ci.chem_interaction_id"""):
        if not r["shippable"] or r["chemical_id"] not in shipped_chem:
            continue
        srcs = ship_sources(r["sources"])
        if not srcs:
            continue
        shipped_inter.add(r["chem_interaction_id"])
        inter_rows.append((
            r["chem_interaction_id"], r["chemical_id"], r["chem_name"], r["hmdb_id"],
            r["receptor_entity_id"], r["receptor_symbol"], r["receptor_uniprot"],
            r["species"], r["complex_id"], r["complex_role"], r["receptor_type"],
            len(srcs), json.dumps(srcs), r["has_citation"], r["n_pmids"],
            r["verification_status"]))
    counts["chemical_interactions"] = write(
        out / "chemical_interactions.csv",
        ["chem_interaction_id", "chemical_id", "chemical_name", "hmdb_id",
         "receptor_entity_id", "receptor_symbol", "receptor_uniprot", "species",
         "complex_id", "complex_role", "receptor_type", "n_sources", "sources",
         "has_citation", "n_pmids", "verification_status"], inter_rows)

    # ---- metabolite_enzymes ------------------------------------------------
    enz_rows = [tuple(r) for r in c.execute(
        """SELECT me.chemical_id, ch.hmdb_id, ch.name, me.gene_symbol, me.entity_id,
                  me.role, me.direction_known, me.source_id, me.raw_role
             FROM metabolite_enzymes me
             JOIN chemical_entities ch ON ch.chemical_id = me.chemical_id
            WHERE me.source_id NOT IN (%s) AND ch.shippable = 1
            ORDER BY me.chemical_id, me.role, me.gene_symbol"""
        % ",".join("?" * len(NON_SHIPPABLE)), tuple(sorted(NON_SHIPPABLE)))
        if r[0] in shipped_chem]
    counts["metabolite_enzymes"] = write(
        out / "metabolite_enzymes.csv",
        ["chemical_id", "hmdb_id", "chemical_name", "gene_symbol", "entity_id", "role",
         "direction_known", "source_id", "raw_role"], enz_rows)

    # ---- chemical evidence -------------------------------------------------
    # New 2026-10-09. Until now the chemical half shipped provenance only, because it
    # had no verified evidence; np10_litsearch.py has since confirmed primary citations
    # for 132 pairs, so they ship as their own file (never merged into evidence.csv,
    # which is the peptide half - decision 1).
    #
    # Licence: a pair that is not in shipped_inter is already excluded, which covers
    # MEBOCOST-only pairs. `notes` is deliberately NOT exported - it holds the raw
    # prompt answer and internal dedup markers naming source_ids, which is exactly the
    # kind of string the assertion below blocks. The useful fields inside it
    # (RESULT, AFFINITY, and the names as the paper writes them) are parsed into real
    # columns instead.
    ev_cols = ["np_evidence_id", "chem_interaction_id", "chemical_name", "receptor_symbol",
               "pmid", "paper_year", "paper_journal", "verification_outcome", "evidence_depth",
               "result", "assay_class", "effect_direction", "affinity",
               "molecule_as_named_in_paper", "receptor_as_named_in_paper",
               "evidence_sentence", "extraction_method", "verified_at"]
    ev_rows = []
    for r in c.execute("""SELECT ne.np_evidence_id, ne.chem_interaction_id, ch.name chem_name,
                                 e.symbol receptor_symbol, ne.pmid, p.year, p.journal,
                                 ne.verification_outcome, ne.evidence_depth, ne.assay_class,
                                 ne.effect_direction, ne.evidence_sentence,
                                 ne.extraction_method, ne.verified_at, ne.notes
                            FROM np_evidence ne
                            JOIN chemical_interactions ci
                              ON ci.chem_interaction_id = ne.chem_interaction_id
                            JOIN chemical_entities ch ON ch.chemical_id = ci.chemical_id
                            JOIN entities e ON e.entity_id = ci.receptor_entity_id
                            LEFT JOIN papers p ON p.pmid = ne.pmid
                           WHERE ne.verification_outcome IN ('confirmed', 'pair_absent')
                           ORDER BY ch.name, e.symbol, ne.pmid"""):
        if r["chem_interaction_id"] not in shipped_inter:
            continue
        ans = {}
        if r["notes"] and "Prompt answer: " in r["notes"]:
            try:
                ans = json.loads(r["notes"].split("Prompt answer: ", 1)[1])
            except Exception:                                # noqa: BLE001
                ans = {}
        # RESULT comes from the np10 prompt answer when there is one; rows decided by the
        # np11 audit carry it as a [RESULT=...] tag in the note instead; failing both, a
        # pair_absent row is a refutation and a confirmed row is a support.
        res = ans.get("RESULT", "")
        if not res:
            m = re.search(r"\[RESULT=([A-Z_]+)\]", r["notes"] or "")
            res = m.group(1) if m else ("REFUTES" if r["verification_outcome"] == "pair_absent"
                                      else "SUPPORTS")
        ev_rows.append((
            r["np_evidence_id"], r["chem_interaction_id"], r["chem_name"],
            r["receptor_symbol"], r["pmid"], r["year"], r["journal"],
            r["verification_outcome"], r["evidence_depth"], res,
            r["assay_class"], r["effect_direction"], ans.get("AFFINITY", ""),
            ans.get("MOLECULE_AS_NAMED_IN_PAPER", ""),
            ans.get("RECEPTOR_AS_NAMED_IN_PAPER", ""),
            r["evidence_sentence"], r["extraction_method"], r["verified_at"]))
    counts["chemical_evidence"] = write(out / "chemical_evidence.csv", ev_cols, ev_rows)

    # ---- provenance --------------------------------------------------------
    prov_rows = [tuple(r) for r in c.execute(
        """SELECT source_id, chem_interaction_id, chemical_id, raw_ligand,
                  raw_receptor, raw_reference
             FROM chem_source_records
            WHERE source_id NOT IN (%s)
            ORDER BY chemical_id, source_id"""
        % ",".join("?" * len(NON_SHIPPABLE)), tuple(sorted(NON_SHIPPABLE)))
        if r[2] in shipped_chem and (r[1] is None or r[1] in shipped_inter)]
    counts["chemical_source_records"] = write(
        out / "chemical_source_records.csv",
        ["source_id", "chem_interaction_id", "chemical_id", "raw_ligand",
         "raw_receptor", "raw_reference"], prov_rows)

    # ---- assertions: fail the build rather than ship ----------------------
    problems: list[str] = []
    blocked_lower = {s.lower() for s in NON_SHIPPABLE} | {"mebocost"}
    for name in ("chemical_entities.csv", "chemical_interactions.csv",
                 "metabolite_enzymes.csv", "chemical_source_records.csv",
                 "chemical_evidence.csv"):
        text = (out / name).read_text(encoding="utf-8").lower()
        for token in blocked_lower:
            if token in text:
                problems.append(f"{name} contains the string '{token}' - "
                                "a non-shippable source reached the export")
    hdr = (out / "chemical_interactions.csv").read_text(encoding="utf-8").split("\n", 1)[0]
    if "confidence" in hdr.lower() or "score" in hdr.lower():
        problems.append("chemical_interactions.csv has a score column; decision 9 says "
                        "v1 ships provenance only")
    orphans = [r[0] for r in csv.reader((out / "chemical_interactions.csv")
                                        .open(encoding="utf-8"))][1:]
    known = {r[0] for r in csv.reader((out / "chemical_entities.csv")
                                      .open(encoding="utf-8"))}
    bad = [r for r in inter_rows if r[1] not in known]
    if bad:
        problems.append(f"{len(bad)} interactions reference a chemical not in "
                        "chemical_entities.csv")
    ev_orphan = [r[0] for r in ev_rows if r[1] not in shipped_inter]
    if ev_orphan:
        problems.append(f"{len(ev_orphan)} evidence rows reference a pair that is not "
                        "in chemical_interactions.csv")
    if problems:
        for p in problems:
            log("np6", f"ABORT: {p}")
        for name in ("chemical_entities.csv", "chemical_interactions.csv",
                     "metabolite_enzymes.csv", "chemical_source_records.csv",
                     "chemical_evidence.csv"):
            (out / name).unlink(missing_ok=True)
        log("np6", "export files removed; nothing shipped")
        return 1
    log("np6", f"assertions passed: no non-shippable source in {len(orphans)} exported "
               "interactions, no score column, referential integrity holds")

    # ---- manifest ----------------------------------------------------------
    withheld = {
        "chemical_entities": c.execute(
            "SELECT count(*) FROM chemical_entities WHERE shippable=0").fetchone()[0],
        "chemical_interactions": c.execute(
            "SELECT count(*) FROM chemical_interactions WHERE shippable=0").fetchone()[0],
        "metabolite_enzymes": c.execute(
            "SELECT count(*) FROM metabolite_enzymes WHERE source_id IN (%s)"
            % ",".join("?" * len(NON_SHIPPABLE)), tuple(sorted(NON_SHIPPABLE))).fetchone()[0],
        "chemical_evidence": c.execute(
            """SELECT count(*) FROM np_evidence ne
                 JOIN chemical_interactions ci
                   ON ci.chem_interaction_id = ne.chem_interaction_id
                WHERE ne.verification_outcome = 'confirmed'
                  AND ci.shippable = 0""").fetchone()[0],
    }
    by_nsrc: dict[str, int] = {}
    for r in inter_rows:
        by_nsrc[str(r[11])] = by_nsrc.get(str(r[11]), 0) + 1
    manifest = {
        "build": "Phase NP (non-peptide ligands), companion to the v2 peptide release",
        "release_date": date.today().isoformat(),
        "relationship_to_peptide_release": (
            "Companion tables only. No chemical row appears in interactions_scored.csv "
            "and no chemical row carries a confidence_score. The two halves are not "
            "comparable and are deliberately not merged."),
        "species": "human only",
        "row_counts": counts,
        "withheld_under_source_licence": {
            "source": sorted(NON_SHIPPABLE),
            "reason": "MEBOCOST forbids redistribution of its data or any part of it",
            "rows_withheld": withheld,
            "note": ("MEBOCOST is also stripped from the provenance lists of rows that "
                     "ship for other reasons, so exported n_sources can be lower than "
                     "the internal value."),
        },
        "interactions_by_source_count": dict(sorted(by_nsrc.items())),
        "caveats": [
            "NOT a scored resource. n_sources is provenance, not confidence; "
            "cross-database consensus is not evidence (PLAN.md invariant 4).",
            "Two kinds of evidence, do not conflate them. (a) The PMIDs carried "
            "through from the sources are as-is and have NOT been verified; none has "
            "had the publication-type audit the peptide half received. (b) "
            "chemical_evidence.csv is different: every row passed mechanical gates - "
            "PMID resolves, not a review, not retracted, and the quoted sentence is "
            "verbatim in the paper's own abstract or publisher full text, with the "
            "molecule and receptor both named in that text. Those rows were confirmed "
            "WITHOUT a curator reading the paper (owner instruction 2026-10-09), so "
            "'confirmed' there means machine-verified against the source text; "
            "extraction_method='llm_litsearch_np_prompt' marks them.",
            f"{by_nsrc.get('1', 0)} of {len(inter_rows)} interactions rest on a single "
            "source, most of them MetalinksDB (largely STITCH-derived) or MetaLigand, "
            "neither of which carries citations.",
            "chebi_id is empty throughout. No source in the union ships one, and none "
            "was guessed.",
            "molecular_class is PROVISIONAL; the controlled vocabulary is unsettled.",
            "provisional_identity=1 marks chemicals whose identity is still awaiting "
            "sign-off (22-hydroxycholesterol and 12-oxo-LTB4 stereochemistry/tautomer).",
            "Enzyme roles are a union vocabulary mapped from five different source "
            "schemes; raw_role preserves each source's original wording. Sources "
            "disagree in places - e.g. SLC18A1/2 for dopamine is efflux in three "
            "sources and influx in MetalinksDB - and both readings are retained.",
        ],
        "sources": [dict(source_id=r["source_id"], name=r["name"], licence=r["licence"],
                         n_rows=r["n_rows"], sha256=r["sha256"], shippable=r["shippable"])
                    for r in c.execute(
                        "SELECT source_id,name,licence,n_rows,sha256,shippable FROM sources "
                        "WHERE source_id IN (SELECT DISTINCT source_id FROM np_run_contributions) "
                        "   OR source_id LIKE 'm%' OR source_id LIKE 'cell%' OR source_id LIKE 'gtopdb%'"
                        " ORDER BY source_id")],
    }
    (out / "chemical_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    log("np6", f"  chemical_manifest.json           written")
    log("np6", f"withheld under licence: {withheld}")
    log("np6", f"exported to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
