#!/usr/bin/env python3
"""Phase NP-2 - build `chemical_entities` from the seven sources.

Merge order matters. HMDB-native sources seed the table (decision 4: HMDB is the
primary key), then the name-only sources attach by name -> synonym -> HMDB. A
chemical with no HMDB id anywhere still gets a row, keyed on its normalised name,
rather than being dropped - GtoPdb carries no HMDB at all and 38% of its ligand
names do not resolve to one.

Two things this script refuses to do:

  * It does not invent identifiers. chebi_id is left NULL throughout; the
    prototype resolved 17 of 48 by hand and the rest stayed NULL, and that is
    still the right behaviour.
  * It does not silently prefer one source when two disagree. Conflicts are
    written to np_identity_conflicts.tsv and the first-seen value is kept, so
    every disagreement is inspectable rather than absorbed.

MEBOCOST (sources.shippable = 0) may contribute chemistry to a row that another
source already established, but a chemical it is the ONLY source for is marked
shippable = 0 and will be withheld at export.

Reads the REVIEWED np_composite_map.tsv: honours EXCLUDE, splits fused tokens on
'|', and prefers manual_* over np1's guesses.

Usage: python3 np2_chemicals.py [--commit]
"""
from __future__ import annotations

import csv
import html
import importlib.util
import json
import re
import sqlite3
import sys
from pathlib import Path

from cur_common import HERE, RAW, connect, log, read_csv
from np_common import NP_SOURCES, eqi, guard_contributions, in_ci

ROOT = HERE.parent
NP_SCHEMA = HERE / "np_schema.sql"
MAP = HERE / "np_composite_map.tsv"
CONFLICTS = HERE / "np_identity_conflicts.tsv"

# Sources whose data must never reach release/. Mirrors sources.shippable.
NON_SHIPPABLE = {"mebocost_human"}
LIGAND_TYPES = ("Metabolite", "Natural product", "Inorganic", "Synthetic organic")

# Fields where two sources disagreeing is EXPECTED rather than a defect:
# PubChem SIDs are per-depositor, so several are correct for one compound.
INFORMATIONAL_FIELDS = {"pubchem_sid"}


def clean_name(s: str) -> str:
    """GtoPdb writes names as HTML: '&alpha;-ketoglutaric acid', 'LTB<sub>4</sub>'."""
    s = re.sub(r"<[^>]+>", "", str(s or ""))
    s = html.unescape(s)
    for uni, plain in (("α", "a"), ("β", "b"), ("γ", "g"),
                       ("δ", "d"), ("Δ", "d"), ("±", "")):
        s = s.replace(uni, plain)
    return s.strip()


def derive_class(row: dict) -> str | None:
    """PROVISIONAL molecular_class. The controlled vocabulary is still an open
    item in PHASE_NONPEPTIDE.md - neither source ontology maps cleanly onto the
    four classes section 7 proposes, so this is a best effort, not a decision."""
    if row.get("is_neurotransmitter"):
        return "neurotransmitter"
    blob = " ".join(str(row.get(k) or "") for k in
                    ("super_class", "class", "metabolite_subclass",
                     "chebi_super_class", "chebi_class")).lower()
    if "steroid" in blob:
        return "steroid"
    if "lipid" in blob or "fatty acy" in blob or "prostagland" in blob or "eicosanoid" in blob:
        return "lipid"
    if "homogeneous metal" in blob or "inorganic" in blob or "non-metal" in blob:
        return "ion"
    return "metabolite" if blob.strip() else None


class ChemStore:
    """Accumulates chemical identity across sources, logging disagreements."""

    def __init__(self, norm) -> None:
        self.norm = norm
        self.by_hmdb: dict[str, dict] = {}
        self.by_key: dict[str, dict] = {}
        self.conflicts: list[dict] = []

    def add(self, source_id: str, name: str, hmdb: str | None = None, **fields) -> dict | None:
        name = clean_name(name)
        key = self.norm(name)
        if not key:
            return None
        hmdb = (hmdb or "").strip() or None

        rec = self.by_hmdb.get(hmdb) if hmdb else None
        if rec is None:
            rec = self.by_key.get(key)
            # A name we have seen before, now arriving with a DIFFERENT hmdb id.
            if rec is not None and hmdb and rec.get("hmdb_id") and rec["hmdb_id"] != hmdb:
                self.conflicts.append(dict(
                    kind="name maps to two HMDB ids", name=name, source_id=source_id,
                    field="hmdb_id", kept=rec["hmdb_id"], discarded=hmdb))
                hmdb = rec["hmdb_id"]
        if rec is None:
            rec = dict(name=name, name_key=key, hmdb_id=hmdb, synonyms=set(),
                       sources=set(), is_neurotransmitter=0)
            self.by_key[key] = rec
        if hmdb and not rec.get("hmdb_id"):
            rec["hmdb_id"] = hmdb
        if rec.get("hmdb_id"):
            self.by_hmdb.setdefault(rec["hmdb_id"], rec)
        self.by_key.setdefault(key, rec)

        rec["sources"].add(source_id)
        for f, v in fields.items():
            if f == "synonyms":
                rec["synonyms"].update(x for x in v if x)
                continue
            if v in (None, "", "NA"):
                continue
            if f == "is_neurotransmitter":
                rec["is_neurotransmitter"] = max(rec["is_neurotransmitter"], int(v))
                continue
            if rec.get(f) and str(rec[f]) != str(v):
                self.conflicts.append(dict(
                    kind=("expected: depositor-specific" if f in INFORMATIONAL_FIELDS
                          else "field disagreement"),
                    name=rec["name"], source_id=source_id,
                    field=f, kept=rec[f], discarded=v))
                continue
            rec[f] = v
        return rec

    def records(self) -> list[dict]:
        seen, out = set(), []
        for rec in self.by_key.values():
            if id(rec) in seen:
                continue
            seen.add(id(rec))
            out.append(rec)
        return out


def load_reviewed_map() -> tuple[dict[tuple[str, str], list[tuple[str, str]]], int]:
    """(source_id, raw_ligand) -> [(chem name, hmdb), ...]. EXCLUDE drops the row."""
    out: dict[tuple[str, str], list[tuple[str, str]]] = {}
    n_excluded = 0
    for r in read_csv(MAP, "\t"):
        manual = (r.get("manual_chem_name") or "").strip()
        if manual == "EXCLUDE":
            n_excluded += 1
            continue
        names = (manual or r["chem_token"]).split("|")
        ids = (r.get("manual_hmdb_id") or r.get("hmdb_id") or "").split("|")
        pairs = [(n.strip(), ids[i].strip() if i < len(ids) else "")
                 for i, n in enumerate(names) if n.strip()]
        out[(r["source_id"], r["raw_ligand"])] = pairs
    return out, n_excluded


def main() -> int:
    commit = "--commit" in sys.argv
    spec = importlib.util.spec_from_file_location("np1", HERE / "np1_parse_composites.py")
    np1 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(np1)

    c = connect()
    c.executescript(NP_SCHEMA.read_text())
    # np_schema.sql uses CREATE TABLE IF NOT EXISTS, so a table built by an
    # earlier run does not pick up new columns. Add them here, idempotently.
    have = {r["name"] for r in c.execute("PRAGMA table_info(chemical_entities)")}
    for col in ("chebi_super_class", "chebi_class"):
        if col not in have:
            c.execute(f"ALTER TABLE chemical_entities ADD COLUMN {col} TEXT")
            log("np2", f"migrated: chemical_entities.{col} added")
    c.commit()

    store = ChemStore(np1.norm_chem)

    # --- 1. HMDB-native sources seed the table ----------------------------
    for r in read_csv(RAW / "mrclinkdb_human_mr.txt", "\t"):
        cid = sid = None
        m = re.search(r"CID:(\d+)", r.get("PubChem CID/SID") or "")
        if m:
            cid = m.group(1)
        m = re.search(r"SID:(\d+)", r.get("PubChem CID/SID") or "")
        if m:
            sid = m.group(1)
        store.add("mrclinkdb_human_mr", r["Metabolite name"], r["HMDB ID"],
                  pubchem_cid=cid, pubchem_sid=sid, smiles=r.get("Canonical SMILES"),
                  formula=r.get("Molecular Formula"), kingdom=r.get("Kingdom"),
                  super_class=r.get("Super Class"), **{"class": r.get("Class")})

    # MetaLigand fuses several distinct chemicals into one row - 13 ligands,
    # 324 of its 4,267 pairs. Expansion is driven off the HMDB id list and names
    # are resolved BY ID, never by splitting the name: '3,5-Dihydroxybenzoic
    # acid' and '5a-Pregnane-3,20-dione' contain commas inside one name, so
    # comma-splitting the name is wrong in both directions.
    name_by_hmdb = {rec["hmdb_id"]: rec["name"] for rec in store.records() if rec.get("hmdb_id")}
    # MetalinksDB's 1,128 names are loaded after this point, but they are the
    # broadest HMDB -> name map available and the fused-ligand expansion needs
    # them now, so read that table up front for naming purposes only.
    with sqlite3.connect(RAW / "metalinksdb.db") as s:
        for hmdb, nm in s.execute("SELECT hmdb, metabolite FROM metabolites"):
            if hmdb and nm:
                name_by_hmdb.setdefault(hmdb, nm)
    cls_by_name = {r["Compounds"]: r for r in read_csv(RAW / "metaligand_Ligand_list.csv", ",")}
    ml_ligands = {r["Compounds"]: r for r in read_csv(RAW / "metaligand_LR_data_human.csv", ",")}
    ml_used = {r["L"] for r in read_csv(RAW / "metaligand_NPLRdb_human.csv", ",")}
    n_ml_fused = 0
    for name in sorted(ml_used):
        r = ml_ligands.get(name)
        if not r:
            continue
        extra = cls_by_name.get(name, {})
        syn = [s.strip() for s in (r.get("Synonyms") or "").split(",")][:40]
        ids = [h.strip() for h in (r.get("HMDB_id") or "").split(",") if h.strip()]
        frags = [p.strip() for p in name.split(",")]
        if len(ids) > 1:
            n_ml_fused += 1
        targets: list[tuple[str, str | None]] = []
        if len(ids) <= 1:
            targets = [(name, ids[0] if ids else None)]
        else:
            for k, h in enumerate(ids):
                # Prefer the canonical name already established for that id.
                nm = name_by_hmdb.get(h)
                if not nm and len(frags) == len(ids):
                    nm = frags[k]
                if not nm:
                    store.conflicts.append(dict(
                        kind="fused ligand: cannot name one constituent",
                        name=name, source_id="metaligand_pairs_human",
                        field="hmdb_id", kept=h, discarded=r.get("HMDB_id")))
                    nm = h
                targets.append((nm, h))
        for nm, h in targets:
            store.add("metaligand_pairs_human", nm, h,
                      synonyms=syn if len(targets) == 1 else [],
                      chebi_super_class=extra.get("SuperClass"),
                      chebi_class=extra.get("Class"))

    with sqlite3.connect(RAW / "metalinksdb.db") as s:
        lr = {h for (h,) in s.execute("SELECT DISTINCT hmdb FROM edges WHERE type='lr'")}
        for hmdb, name, pubchem, sub in s.execute(
                "SELECT hmdb, metabolite, pubchem, metabolite_subclass FROM metabolites"):
            if hmdb in lr:
                store.add("metalinksdb", name or hmdb, hmdb,
                          pubchem_cid=pubchem, metabolite_subclass=sub)

    for r in read_csv(RAW / "mebocost_human_met_sensor.tsv", "\t"):
        store.add("mebocost_human", r["standard_metName"], r["HMDB_ID"])

    # --- 2. name-only sources attach ---------------------------------------
    reviewed, n_excluded = load_reviewed_map()
    nt = set()
    try:
        sys.path.insert(0, str(ROOT / "src"))
        from rdata import read_rdata
        cc = read_rdata(ROOT / "data" / "raw" / "CellChatDB.human.rda")
        it = (cc.get("CellChatDB.human") or next(iter(cc.values())))["interaction"]
        for i in range(len(it["interaction_name"])):
            if it["annotation"][i] == "Non-protein Signaling" and it["is_neurotransmitter"][i] == "TRUE":
                nt.add((it["ligand"][i] or "").strip())
    except Exception as e:                                       # noqa: BLE001
        log("np2", f"WARNING: CellChatDB not read ({e}); is_neurotransmitter left 0")

    for (src, raw_ligand), pairs in reviewed.items():
        for name, hmdb in pairs:
            store.add(src, name, hmdb,
                      is_neurotransmitter=1 if raw_ligand in nt else 0)

    f = (RAW / "gtopdb_interactions.csv").open(encoding="utf-8", errors="replace")
    next(f)
    n_gt = 0
    for r in csv.DictReader(f):
        if not eqi(r["Endogenous"], "true") or not eqi(r["Target Species"], "human"):
            continue
        if not in_ci(r["Ligand Type"], LIGAND_TYPES):
            continue
        store.add("gtopdb_interactions", r["Ligand"], None,
                  pubchem_sid=r.get("Ligand PubChem SID") or None)
        n_gt += 1

    # --- 3. write ----------------------------------------------------------
    recs = store.records()
    recs.sort(key=lambda x: (x.get("hmdb_id") or "zzz", x["name_key"]))
    rows = []
    for i, rec in enumerate(recs, 1):
        srcs = sorted(rec["sources"])
        rows.append((
            f"CHEM{i:05d}", rec["name_key"], rec.get("hmdb_id"), None,
            rec.get("pubchem_cid"), rec.get("pubchem_sid"), rec["name"],
            json.dumps(sorted(rec["synonyms"])[:40]) if rec["synonyms"] else None,
            rec.get("smiles"), rec.get("formula"), rec.get("kingdom"),
            rec.get("super_class"), rec.get("class"),
            rec.get("chebi_super_class"), rec.get("chebi_class"),
            rec.get("metabolite_subclass"),
            derive_class(rec), rec["is_neurotransmitter"], len(srcs), json.dumps(srcs),
            0 if set(srcs) <= NON_SHIPPABLE else 1))

    per_source: dict[str, int] = {}
    for rec in recs:
        for sid in rec["sources"]:
            per_source[sid] = per_source.get(sid, 0) + 1
    if guard_contributions(c, "np2", per_source, NP_SOURCES, commit):
        c.rollback()
        return 1

    c.execute("DELETE FROM chemical_entities")
    c.executemany(
        "INSERT INTO chemical_entities (chemical_id,name_key,hmdb_id,chebi_id,pubchem_cid,"
        "pubchem_sid,name,synonyms,smiles,formula,kingdom,super_class,class,"
        "chebi_super_class,chebi_class,metabolite_subclass,molecular_class,"
        "is_neurotransmitter,n_sources,sources,shippable)"
        " VALUES (" + ",".join("?" * 21) + ")", rows)

    if store.conflicts:
        with CONFLICTS.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["kind", "name", "source_id", "field",
                                               "kept", "discarded"], delimiter="\t")
            w.writeheader()
            w.writerows(store.conflicts)

    if commit:
        c.commit()
    else:
        c.rollback()

    # --- 4. report ---------------------------------------------------------
    n = len(rows)
    with_hmdb = sum(1 for r in rows if r[2])
    log("np2", f"chemical_entities: {n} ({with_hmdb} with an HMDB id, "
               f"{n - with_hmdb} name-keyed only)")
    log("np2", f"chebi_id: 0 populated - left NULL rather than guessed (mapping pass is separate)")
    log("np2", f"excluded by review: {n_excluded} composite rows | GtoPdb rows read: {n_gt}")
    log("np2", f"MetaLigand fused ligands expanded by HMDB id: {n_ml_fused}")
    log("np2", f"not shippable (MEBOCOST-only): {sum(1 for r in rows if not r[20])}")
    by_src: dict[str, int] = {}
    for r in rows:
        for s in json.loads(r[19]):
            by_src[s] = by_src.get(s, 0) + 1
    for s, k in sorted(by_src.items(), key=lambda x: -x[1]):
        log("np2", f"  contributed by {s:26s} {k:5d}")
    multi = sum(1 for r in rows if r[18] > 1)
    log("np2", f"chemicals seen by >1 source: {multi} ({multi/n:.0%})")
    cls: dict[str, int] = {}
    for r in rows:
        cls[r[16] or "(none)"] = cls.get(r[16] or "(none)", 0) + 1
    log("np2", "molecular_class (PROVISIONAL): " +
        ", ".join(f"{k}={v}" for k, v in sorted(cls.items(), key=lambda x: -x[1])))
    if store.conflicts:
        log("np2", f"identity conflicts: {len(store.conflicts)} -> {CONFLICTS.name} "
                   "(first value kept, nothing silently merged)")
    log("np2", "COMMITTED" if commit else "DRY RUN - rerun with --commit to write")
    return 0


if __name__ == "__main__":
    sys.exit(main())
