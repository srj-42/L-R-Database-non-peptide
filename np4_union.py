#!/usr/bin/env python3
"""Phase NP-4 - collapse every source into `chemical_interactions`.

One row per (chemical, receptor, species) - the same grain as the peptide
`interactions` table (decision 5), with the contributing sources carried as a
JSON list rather than as separate rows. Every contributing source row is written
verbatim to `chem_source_records`, including the original composite ligand
string, so nothing is lost by the collapse.

Three expansions happen here, and all three multiply rows:

  * FUSED LIGANDS. MetaLigand puts two chemicals in one pair row (324 of 4,267)
    and CellPhoneDB does it twice. Each becomes one row per constituent.
  * COMPLEX RECEPTORS. CellChatDB, CellPhoneDB and MetaLigand all encode them,
    in three different spellings. np3 already exploded them; each subunit gets
    its own row, tied by complex_id (decision 6).
  * MULTI-PMID ROWS. Kept as a count here; the citations themselves live in
    chem_source_records.raw_reference until the verification phase.

No confidence_score column exists (decision 9). What ships in its place is
n_sources / sources / has_citation / n_pmids.

Usage: python3 np4_union.py [--commit]
"""
from __future__ import annotations

import csv
import importlib.util
import json
import sqlite3
import sys
from collections import defaultdict

from cur_common import HERE, RAW, connect, log, pmids_in, read_csv
from np_common import NP_SOURCES, eqi, guard_contributions, in_ci

ROOT = HERE.parent
NP_SCHEMA = HERE / "np_schema.sql"
MAP = HERE / "np_composite_map.tsv"
DROPPED = HERE / "np_dropped_rows.tsv"
SPECIES = "human"
LIGAND_TYPES = ("Metabolite", "Natural product", "Inorganic", "Synthetic organic")
NON_SHIPPABLE = {"mebocost_human"}


def main() -> int:
    commit = "--commit" in sys.argv
    spec = importlib.util.spec_from_file_location("np2", HERE / "np2_chemicals.py")
    np2 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(np2)
    spec1 = importlib.util.spec_from_file_location("np1", HERE / "np1_parse_composites.py")
    np1 = importlib.util.module_from_spec(spec1)
    spec1.loader.exec_module(np1)

    c = connect()
    c.executescript(NP_SCHEMA.read_text())
    c.commit()

    # --- lookups ------------------------------------------------------------
    chem_by_hmdb: dict[str, str] = {}
    chem_by_key: dict[str, str] = {}
    for r in c.execute("SELECT chemical_id, hmdb_id, name_key FROM chemical_entities"):
        if r["hmdb_id"]:
            chem_by_hmdb[r["hmdb_id"]] = r["chemical_id"]
        chem_by_key[r["name_key"]] = r["chemical_id"]

    rec_map: dict[tuple[str, str], list[sqlite3.Row]] = defaultdict(list)
    for r in c.execute("SELECT * FROM np_receptor_map"):
        rec_map[(r["source_id"], r["raw_receptor"])].append(r)

    reviewed, _ = np2.load_reviewed_map()
    excluded = {k for k, v in
                ((k, (r.get("manual_chem_name") or "").strip())
                 for k, r in (((r["source_id"], r["raw_ligand"]), r)
                              for r in read_csv(MAP, "\t")))
                if v == "EXCLUDE"}

    # np1's index is the name-resolution service: it carries every source's
    # names, MetaLigand's synonym lists and the hand-verified abbreviations.
    # Without it GtoPdb's 'GABA', 'PGE<sub>2</sub>' and '(-)-adrenaline' fail to
    # attach to chemicals that are already in the table under another name.
    name_index = np1.chem_index()

    def chems_for(names_ids: list[tuple[str, str]]) -> list[str]:
        out = []
        for name, hmdb in names_ids:
            cid = chem_by_hmdb.get(hmdb) if hmdb else None
            key = np1.norm_chem(np2.clean_name(name))
            if not cid:
                cid = chem_by_key.get(key)
            if not cid:
                hit = name_index.get(key)
                if hit:
                    cid = (chem_by_hmdb.get(hit[1]) if hit[1] else None) \
                          or chem_by_key.get(np1.norm_chem(hit[0]))
            if cid:
                out.append(cid)
        return out

    # MetaLigand fused ligands: expand by HMDB id, exactly as np2 did.
    ml_ligands = {r["Compounds"]: r for r in read_csv(RAW / "metaligand_LR_data_human.csv", ",")}

    def metaligand_chems(token: str) -> list[str]:
        r = ml_ligands.get(token)
        if not r:
            return chems_for([(token, "")])
        ids = [h.strip() for h in (r.get("HMDB_id") or "").split(",") if h.strip()]
        if len(ids) > 1:
            return [chem_by_hmdb[h] for h in ids if h in chem_by_hmdb]
        return chems_for([(token, ids[0] if ids else "")])

    # --- accumulate ---------------------------------------------------------
    # (chemical_id, entity_id) -> dict
    pairs: dict[tuple[str, str], dict] = {}
    src_records: list[tuple] = []
    dropped: list[dict] = []
    stats: dict[str, int] = defaultdict(int)

    def emit(source_id, chem_ids, raw_receptor, raw_ligand, reference, raw_row):
        recs = rec_map.get((source_id, raw_receptor), [])
        if not chem_ids:
            dropped.append(dict(source_id=source_id, raw_ligand=raw_ligand,
                                raw_receptor=raw_receptor, reason="no chemical resolved"))
            return
        if not recs:
            dropped.append(dict(source_id=source_id, raw_ligand=raw_ligand,
                                raw_receptor=raw_receptor, reason="no receptor resolved"))
            return
        pm = pmids_in(reference)
        for cid in chem_ids:
            for rec in recs:
                key = (cid, rec["entity_id"])
                p = pairs.get(key)
                if p is None:
                    p = pairs[key] = dict(
                        chemical_id=cid, entity_id=rec["entity_id"],
                        complex_id=rec["complex_id"], complex_role=rec["complex_role"],
                        receptor_type=rec["receptor_type"], sources=set(), pmids=set())
                p["sources"].add(source_id)
                p["pmids"].update(pm)
                # Carry the pair key so the record can be tied to its
                # chem_interaction_id once the ids are minted below. Without
                # this, provenance is only recoverable per chemical, and a
                # pair's own citations cannot be shown.
                src_records.append((source_id, cid, rec["entity_id"], raw_ligand,
                                    raw_receptor, reference or None, raw_row))
                stats[f"rows:{source_id}"] += 1

    # 1. MRCLinkdb
    for r in read_csv(RAW / "mrclinkdb_human_mr.txt", "\t"):
        cid = chem_by_hmdb.get((r["HMDB ID"] or "").strip())
        emit("mrclinkdb_human_mr", [cid] if cid else [], (r["Receptor_symbol"] or "").strip(),
             r["Metabolite name"], r.get("PMID"), json.dumps(r))

    # 2. MetaLigand
    for r in read_csv(RAW / "metaligand_NPLRdb_human.csv", ","):
        emit("metaligand_pairs_human", metaligand_chems(r["L"]), (r["R"] or "").strip(),
             r["L"], None, json.dumps(r))

    # 3. MetalinksDB (lr edges only)
    with sqlite3.connect(RAW / "metalinksdb.db") as s:
        for hmdb, uni, source, mor in s.execute(
                "SELECT hmdb, uniprot, source, mor FROM edges WHERE type='lr'"):
            cid = chem_by_hmdb.get(hmdb)
            emit("metalinksdb", [cid] if cid else [], uni, hmdb, None,
                 json.dumps(dict(hmdb=hmdb, uniprot=uni, source=source, mor=mor)))

    # 4. MEBOCOST - internal only, shippable=0 downstream
    for r in read_csv(RAW / "mebocost_human_met_sensor.tsv", "\t"):
        cid = chem_by_hmdb.get((r["HMDB_ID"] or "").strip())
        emit("mebocost_human", [cid] if cid else [], (r["Gene_name"] or "").strip(),
             r["standard_metName"], r.get("Evidence"), json.dumps(r))

    # 5. CellChatDB v2
    sys.path.insert(0, str(ROOT / "src"))
    from rdata import read_rdata
    cc = read_rdata(ROOT / "data" / "raw" / "CellChatDB.human.rda")
    it = (cc.get("CellChatDB.human") or next(iter(cc.values())))["interaction"]
    for i in range(len(it["interaction_name"])):
        if not eqi(it["annotation"][i], "Non-protein Signaling"):
            continue
        lig = (it["ligand"][i] or "").strip()
        key = ("cellchatdb_human", lig)
        if key in excluded:
            stats["excluded_by_review"] += 1
            continue
        emit("cellchatdb_human", chems_for(reviewed.get(key, [])),
             (it["receptor"][i] or "").strip(), lig, it["evidence"][i],
             json.dumps(dict(interaction_name=it["interaction_name"][i],
                             pathway_name=it["pathway_name"][i],
                             ligand=lig, receptor=it["receptor"][i],
                             evidence=it["evidence"][i])))

    # 6. CellPhoneDB v5
    for r in read_csv(ROOT / "data" / "raw" / "cellphonedb_interactions.csv", ","):
        if not eqi(r.get("is_ppi"), "false"):
            continue
        lig = (r.get("partner_a") or "").strip()
        key = ("cellphonedb_interactions", lig)
        if key in excluded:
            stats["excluded_by_review"] += 1
            continue
        emit("cellphonedb_interactions", chems_for(reviewed.get(key, [])),
             (r.get("partner_b") or "").strip(), lig, r.get("source"), json.dumps(r))

    # 7. GtoPdb - name-keyed, no HMDB anywhere in the file
    f = (RAW / "gtopdb_interactions.csv").open(encoding="utf-8", errors="replace")
    next(f)
    for r in csv.DictReader(f):
        if not eqi(r["Endogenous"], "true") or not eqi(r["Target Species"], SPECIES):
            continue
        if not in_ci(r["Ligand Type"], LIGAND_TYPES):
            continue
        raw_rec = (r.get("Target Gene Symbol") or "").strip() or (r.get("Target UniProt ID") or "").strip()
        emit("gtopdb_interactions", chems_for([(r["Ligand"], "")]), raw_rec,
             r["Ligand"], r.get("PubMed ID"), json.dumps(r))

    # --- guard --------------------------------------------------------------
    # See np_common: zero contribution, or a collapse relative to the previous
    # run, stops the build. The 2026-09-24 GtoPdb incident produced neither an
    # error nor a dropped row - just an absence.
    per_source = {s: n for s, n in
                  ((k.split(":", 1)[1], v) for k, v in stats.items() if k.startswith("rows:"))}
    if guard_contributions(c, "np4", per_source, NP_SOURCES, commit):
        c.rollback()
        return 1

    # --- write --------------------------------------------------------------
    rows = []
    for i, ((cid, eid), p) in enumerate(sorted(pairs.items()), 1):
        srcs = sorted(p["sources"])
        rows.append((f"CI{i:06d}", cid, eid, SPECIES, p["complex_id"], p["complex_role"],
                     p["receptor_type"], len(srcs), json.dumps(srcs),
                     1 if p["pmids"] else 0, len(p["pmids"]), "unverified",
                     0 if set(srcs) <= NON_SHIPPABLE else 1))

    c.execute("DELETE FROM chemical_interactions")
    c.execute("DELETE FROM chem_source_records")
    c.executemany(
        "INSERT INTO chemical_interactions (chem_interaction_id,chemical_id,receptor_entity_id,"
        "species,complex_id,complex_role,receptor_type,n_sources,sources,has_citation,"
        "n_pmids,verification_status,shippable) VALUES (" + ",".join("?" * 13) + ")", rows)

    idx = {(r[1], r[2]): r[0] for r in rows}
    c.executemany(
        "INSERT INTO chem_source_records (source_id,chem_interaction_id,chemical_id,"
        "raw_ligand,raw_receptor,raw_reference,raw_row) VALUES (?,?,?,?,?,?,?)",
        [(s, idx.get((cid, eid)), cid, rl, rr, ref, raw)
         for (s, cid, eid, rl, rr, ref, raw) in src_records])

    if dropped:
        with DROPPED.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["source_id", "raw_ligand", "raw_receptor",
                                               "reason"], delimiter="\t")
            w.writeheader()
            w.writerows(dropped)

    if commit:
        c.commit()
    else:
        c.rollback()

    # --- report -------------------------------------------------------------
    log("np4", f"chemical_interactions: {len(rows)} pairs "
               f"(chemicals {len({r[1] for r in rows})}, receptors {len({r[2] for r in rows})})")
    log("np4", f"chem_source_records: {len(src_records)}")
    cited = sum(1 for r in rows if r[9])
    log("np4", f"with >=1 PMID: {cited} ({cited/len(rows):.0%}) | total distinct PMIDs: "
               f"{sum(r[10] for r in rows)}")
    log("np4", f"not shippable (MEBOCOST-only): {sum(1 for r in rows if not r[12])}")
    log("np4", f"in a receptor complex: {sum(1 for r in rows if r[4])}")
    by_n: dict[int, int] = defaultdict(int)
    for r in rows:
        by_n[r[7]] += 1
    log("np4", "pairs by source count: " +
        ", ".join(f"{k}src={by_n[k]}" for k in sorted(by_n)))
    if stats.get("excluded_by_review"):
        log("np4", f"  excluded_by_review: {stats['excluded_by_review']}")
    if dropped:
        reasons: dict[str, int] = defaultdict(int)
        for d in dropped:
            reasons[f"{d['source_id']}: {d['reason']}"] += 1
        log("np4", f"dropped source rows: {len(dropped)} -> {DROPPED.name}")
        for k, v in sorted(reasons.items(), key=lambda x: -x[1]):
            log("np4", f"  {k}: {v}")
    log("np4", "COMMITTED" if commit else "DRY RUN - rerun with --commit to write")
    return 0


if __name__ == "__main__":
    sys.exit(main())
