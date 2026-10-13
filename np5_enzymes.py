#!/usr/bin/env python3
"""Phase NP-5 - build `metabolite_enzymes` under the five-role union vocabulary.

This is the relation that makes a chemical ligand measurable in transcriptomics
at all: you cannot measure a metabolite's mRNA, so its biosynthetic and transport
enzymes stand in for it. It is a DIFFERENT claim from ligand-receptor - "cell X
can make this ligand" is not "X binds Y" - and is kept in its own table
(PROJECT_REPORT.md section 7, decision 8).

    synthesis            makes the ligand itself
    precursor_synthesis  makes its precursor
    degradation          consumes or inactivates it
    efflux               exports it, incl. synaptic-vesicle loading
    influx               imports it or its precursor (salvage)

Every source's native labels map in deterministically, and `raw_role` keeps the
original wording so no mapping decision here is irreversible. The five roles were
chosen because between them the sources distinguish exactly these; collapsing to
produce/consume/transport would discard MetaLigand's precursor and salvage steps
and MetalinksDB's transport direction.

A precursor is NOT the ligand. `5-HT-DOPA-DDC` means serotonin made via DOPA
decarboxylase, and CellChatDB's `DOPA` token means dopamine with L-DOPA as its
precursor - the np2 DOPA trap. Precursor relationships live here, never in the
ligand identity.

Usage: python3 np5_enzymes.py [--commit]
"""
from __future__ import annotations

import importlib.util
import re
import sqlite3
import sys
from collections import defaultdict

from cur_common import HERE, RAW, connect, log, read_csv
from np_common import guard_contributions

ROOT = HERE.parent
NP_SCHEMA = HERE / "np_schema.sql"
MAP = HERE / "np_composite_map.tsv"

# Sources that carry enzyme/transport information. GtoPdb has none - it is a
# pharmacology resource, not a metabolic one.
ENZYME_SOURCES = {
    "mrclinkdb_human_enzyme",
    "mrclinkdb_human_transporter",
    "metaligand_ligands_human",
    "metalinksdb",
    "cellchatdb_human",
    "cellphonedb_interactions",
    "mebocost_human_enzymes",
}
NON_SHIPPABLE = {"mebocost_human", "mebocost_human_enzymes"}

# Fallback only, for the composite-string enzymes where MetalinksDB has no type.
TRANSPORTER_RE = re.compile(r"^(SLC|ABC)")


def main() -> int:
    commit = "--commit" in sys.argv
    spec = importlib.util.spec_from_file_location("np4", HERE / "np4_union.py")
    spec1 = importlib.util.spec_from_file_location("np1", HERE / "np1_parse_composites.py")
    np1 = importlib.util.module_from_spec(spec1)
    spec1.loader.exec_module(np1)

    c = connect()
    c.executescript(NP_SCHEMA.read_text())
    c.commit()

    chem_by_hmdb: dict[str, str] = {}
    chem_by_key: dict[str, str] = {}
    for r in c.execute("SELECT chemical_id, hmdb_id, name_key FROM chemical_entities"):
        if r["hmdb_id"]:
            chem_by_hmdb[r["hmdb_id"]] = r["chemical_id"]
        chem_by_key[r["name_key"]] = r["chemical_id"]

    ent_by_symbol = {(r["symbol"] or "").upper(): r["entity_id"] for r in c.execute(
        "SELECT entity_id, symbol FROM entities WHERE species='human' AND kind='protein'")}

    # MetalinksDB protein_type, used to decide whether a composite-string enzyme
    # is a synthase or a transporter without guessing from its name.
    ptype: dict[str, str] = {}
    with sqlite3.connect(RAW / "metalinksdb.db") as s:
        for sym, t in s.execute("SELECT gene_symbol, protein_type FROM proteins "
                                "WHERE gene_symbol IS NOT NULL AND protein_type IS NOT NULL"):
            ptype.setdefault(sym.upper(), t.strip().strip('"'))

    rows: set[tuple] = set()
    counts: dict[str, int] = defaultdict(int)
    unresolved: dict[str, int] = defaultdict(int)

    def add(chem_id, gene, role, source_id, raw_role, direction_known=1):
        gene = (gene or "").strip()
        if not chem_id or not gene:
            return
        rows.add((chem_id, gene.upper(), ent_by_symbol.get(gene.upper()),
                  role, source_id, raw_role, direction_known))
        counts[source_id] += 1

    # --- 1. MRCLinkdb: signed produce/consume -------------------------------
    for r in read_csv(RAW / "mrclinkdb_human_enzyme.txt", "\t"):
        cid = chem_by_hmdb.get((r.get("HMDB_ID") or "").strip())
        if not cid:
            unresolved["mrclinkdb_human_enzyme"] += 1
            continue
        # The website file names this column 'enzyme product/substrate' (the
        # thinner Zenodo deposit calls it enzyme_p_s); values are p / s / s+p.
        flag = (r.get("enzyme product/substrate") or r.get("enzyme_p_s") or "").strip().lower()
        # 's+p' is a genuine dual-capacity enzyme - 70 rows, 3.6%, matching the
        # paper's figure. It earns BOTH roles, not neither.
        roles = {"p": ["synthesis"], "s": ["degradation"],
                 "s+p": ["synthesis", "degradation"]}.get(flag, [])
        for role in roles:
            add(cid, r.get("GENE_NAME"), role, "mrclinkdb_human_enzyme",
                f"enzyme product/substrate={flag}")

    # MRCLinkdb's transporter file states no direction. Recorded as efflux with
    # direction_known=0 rather than silently asserting export.
    for r in read_csv(RAW / "mrclinkdb_human_transporter.txt", "\t"):
        cid = chem_by_hmdb.get((r.get("HMDB_ID") or "").strip())
        if not cid:
            unresolved["mrclinkdb_human_transporter"] += 1
            continue
        add(cid, r.get("GENE_NAME"), "efflux", "mrclinkdb_human_transporter",
            "transporter file, no direction stated", direction_known=0)

    # --- 2. MetaLigand: four typed gene groups ------------------------------
    ML_ROLES = {"Synthetic_genes": "synthesis",
                "Synthetic_genes_lv2": "precursor_synthesis",
                "Transporter_genes": "efflux",
                "Transporter_genes_in": "influx"}
    for r in read_csv(RAW / "metaligand_LR_data_human.csv", ","):
        ids = [h.strip() for h in (r.get("HMDB_id") or "").split(",") if h.strip()]
        cids = [chem_by_hmdb[h] for h in ids if h in chem_by_hmdb]
        if not cids:
            cids = [chem_by_key.get(np1.norm_chem(r.get("Compounds", "")))]
            cids = [x for x in cids if x]
        if not cids:
            unresolved["metaligand_ligands_human"] += 1
            continue
        for col, role in ML_ROLES.items():
            for g in dict.fromkeys(x.strip() for x in (r.get(col) or "").split(",")):
                for cid in cids:
                    add(cid, g, role, "metaligand_ligands_human", col)

    # --- 3. MetalinksDB: pd edges + transport direction ---------------------
    # Where transport_direction is set it WINS over mor: an edge marked 'out' is
    # an export event whatever its mode-of-regulation sign says.
    with sqlite3.connect(RAW / "metalinksdb.db") as s:
        sym_by_uni = {u: g for u, g in s.execute(
            "SELECT uniprot, gene_symbol FROM proteins WHERE gene_symbol IS NOT NULL")}
        for hmdb, uni, mor, tdir in s.execute(
                "SELECT hmdb, uniprot, mor, transport_direction FROM edges WHERE type='pd'"):
            cid = chem_by_hmdb.get(hmdb)
            if not cid:
                unresolved["metalinksdb"] += 1
                continue
            if tdir == "out":
                role, raw = "efflux", "pd/transport_direction=out"
            elif tdir == "in":
                role, raw = "influx", "pd/transport_direction=in"
            elif mor is not None and mor > 0:
                role, raw = "synthesis", "pd/mor=+1"
            elif mor is not None and mor < 0:
                role, raw = "degradation", "pd/mor=-1"
            else:
                unresolved["metalinksdb"] += 1
                continue
            add(cid, sym_by_uni.get(uni), role, "metalinksdb", raw)

    # --- 4. CellChatDB / CellPhoneDB: enzymes parsed out of the ligand string -
    for r in read_csv(MAP, "\t"):
        if (r.get("manual_chem_name") or "").strip() == "EXCLUDE":
            continue
        genes = [g for g in (r.get("enzyme_genes") or "").split(";") if g.strip()]
        genes += [g for g in (r.get("manual_enzyme_genes") or "").split(";") if g.strip()]
        if not genes:
            continue
        names = ((r.get("manual_chem_name") or "").strip() or r["chem_token"]).split("|")
        ids = ((r.get("manual_hmdb_id") or "").strip() or r.get("hmdb_id") or "").split("|")
        cids = []
        for i, nm in enumerate(names):
            cid = chem_by_hmdb.get(ids[i].strip()) if i < len(ids) and ids[i].strip() else None
            cid = cid or chem_by_key.get(np1.norm_chem(nm))
            if cid:
                cids.append(cid)
        if not cids:
            unresolved[r["source_id"]] += 1
            continue
        for g in dict.fromkeys(x.strip() for x in genes):
            t = ptype.get(g.upper(), "")
            if t == "transporter" or (not t and TRANSPORTER_RE.match(g.upper())):
                role, raw = "efflux", "composite string, transporter"
            else:
                role, raw = "synthesis", "composite string, synthase"
            for cid in cids:
                add(cid, g, role, r["source_id"], raw)

    # --- 5. MEBOCOST enzymes - internal cross-check only --------------------
    for r in read_csv(RAW / "mebocost_human_gene_reaction.tsv", "\t"):
        cid = chem_by_hmdb.get((r.get("HMDB_ID") or "").strip())
        if not cid:
            unresolved["mebocost_human_enzymes"] += 1
            continue
        d = (r.get("direction") or "").strip().lower()
        role = {"product": "synthesis", "substrate": "degradation"}.get(d)
        if not role:
            continue
        for tok in (r.get("gene") or "").split(";"):
            g = tok.split("[")[0].strip()
            if g:
                add(cid, g, role, "mebocost_human_enzymes", f"direction={d}")

    # --- guard + write ------------------------------------------------------
    if guard_contributions(c, "np5", dict(counts), ENZYME_SOURCES, commit):
        c.rollback()
        return 1

    c.execute("DELETE FROM metabolite_enzymes")
    c.executemany(
        "INSERT OR IGNORE INTO metabolite_enzymes (chemical_id,gene_symbol,entity_id,role,"
        "source_id,raw_role,direction_known) VALUES (?,?,?,?,?,?,?)", sorted(rows))

    if commit:
        c.commit()
    else:
        c.rollback()

    # --- report -------------------------------------------------------------
    log("np5", f"metabolite_enzymes: {len(rows)} rows over "
               f"{len({r[0] for r in rows})} chemicals and {len({r[1] for r in rows})} genes")
    by_role: dict[str, int] = defaultdict(int)
    for r in rows:
        by_role[r[3]] += 1
    log("np5", "by role: " + ", ".join(f"{k}={v}" for k, v in
                                       sorted(by_role.items(), key=lambda x: -x[1])))
    log("np5", f"direction not stated by the source: {sum(1 for r in rows if not r[6])}"
               " (MRCLinkdb transporters, recorded as efflux but flagged)")
    log("np5", f"resolved to an entities row: {sum(1 for r in rows if r[2])}/{len(rows)}"
               " (metabolic enzymes are largely absent from a peptide corpus; gene_symbol"
               " is the authoritative key here, entity_id is a convenience)")
    # The metric that actually matters: a chemical ligand with no enzymes cannot
    # be detected in transcriptomic data at all, so this is the usable fraction.
    n_chem = c.execute("SELECT count(*) FROM chemical_entities").fetchone()[0]
    with_enz = len({r[0] for r in rows})
    log("np5", f"chemicals with >=1 enzyme: {with_enz}/{n_chem} ({with_enz/n_chem:.0%}) "
               f"- the rest are not inferable from expression data")
    both = len({r[0] for r in rows if r[3] in ("synthesis", "precursor_synthesis")}
               & {r[0] for r in rows if r[3] in ("efflux", "influx")})
    log("np5", f"chemicals with both a synthesis and a transport route: {both}")
    ship = sum(1 for r in rows if r[4] not in NON_SHIPPABLE)
    log("np5", f"shippable rows: {ship} | withheld (MEBOCOST licence): {len(rows) - ship}")
    if unresolved:
        log("np5", "chemical not in chemical_entities, rows skipped: " +
            ", ".join(f"{k}={v}" for k, v in sorted(unresolved.items(), key=lambda x: -x[1])))
    log("np5", "COMMITTED" if commit else "DRY RUN - rerun with --commit to write")
    return 0


if __name__ == "__main__":
    sys.exit(main())
