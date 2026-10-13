#!/usr/bin/env python3
"""Phase NP-3 - resolve every non-peptide receptor to an `entities` row.

Receptors reuse the existing `entities` table (decision 6); only the chemical
LIGANDS get their own table. This script builds `np_receptor_map`, the work table
np4 joins against: one row per (source, raw receptor token, resolved entity),
with a complex exploded into one row per subunit under a shared complex_id.

Receptor tokens arrive in three shapes:

    uniprot-native   MRCLinkdb, CellPhoneDB, MetalinksDB, GtoPdb
    symbol-native    MetaLigand, MEBOCOST, CellChatDB
    complexes        CellChatDB (receptor.symbol, comma-separated subunits)
                     CellPhoneDB (partner_b names a row in complex_input)

*** THIS SCRIPT ADDS ROWS TO `entities`, WHICH THE PEPTIDE BUILD SHARES. ***
Many ionotropic receptor subunits are simply absent from a peptide-ligand
corpus, so without this the pairs that use them would be dropped. Rows created
here take an `N#####` entity_id rather than the peptide build's `E######`, so
every one of them is identifiable and reversible with a single DELETE. No
existing entity row is modified.

receptor_type comes from MetalinksDB `proteins.protein_type` (decision 10).
MEBOCOST's `Annotation` is used only to cross-check and report agreement - it is
a non-shippable source and never supplies a written value.

Usage: python3 np3_receptors.py [--commit]
"""
from __future__ import annotations

import csv
import json
import sqlite3
import sys
from collections import defaultdict

from cur_common import HERE, RAW, connect, log, read_csv
from np_common import NP_SOURCES, eqi, guard_contributions, in_ci

ROOT = HERE.parent
NP_SCHEMA = HERE / "np_schema.sql"
MAP = HERE / "np_composite_map.tsv"
UNRESOLVED = HERE / "np_unresolved_receptors.tsv"

SPECIES = "human"
LIGAND_TYPES = ("Metabolite", "Natural product", "Inorganic", "Synthetic organic")


def norm_sym(s) -> str:
    return str(s or "").strip().upper()


def main() -> int:
    commit = "--commit" in sys.argv
    c = connect()
    c.executescript(NP_SCHEMA.read_text())
    c.commit()

    # --- existing entity indexes ------------------------------------------
    by_symbol: dict[str, str] = {}
    by_uniprot: dict[str, str] = {}
    for r in c.execute("SELECT entity_id, symbol, uniprot FROM entities "
                       "WHERE species=? AND kind='protein'", (SPECIES,)):
        if r["symbol"]:
            by_symbol.setdefault(norm_sym(r["symbol"]), r["entity_id"])
        if r["uniprot"]:
            by_uniprot.setdefault(r["uniprot"].strip(), r["entity_id"])
    n_start = len(by_symbol)
    log("np3", f"existing human protein entities: {n_start} by symbol, {len(by_uniprot)} by uniprot")

    # --- annotation lookups ------------------------------------------------
    rtype: dict[str, str] = {}          # gene symbol -> receptor_type
    uni2sym: dict[str, str] = {}
    with sqlite3.connect(RAW / "metalinksdb.db") as s:
        for uni, sym, ptype in s.execute("SELECT uniprot, gene_symbol, protein_type FROM proteins"):
            if sym:
                uni2sym.setdefault(uni, sym)
                if ptype:
                    rtype.setdefault(norm_sym(sym), ptype.strip().strip('"'))
    mebo_type: dict[str, str] = {}
    for r in read_csv(RAW / "mebocost_human_met_sensor.tsv", "\t"):
        if r.get("Gene_name") and r.get("Annotation"):
            mebo_type.setdefault(norm_sym(r["Gene_name"]), r["Annotation"].strip())

    # --- collect raw receptor tokens per source ---------------------------
    # token -> (symbols, uniprots, complex_name or None)
    tokens: dict[tuple[str, str], tuple[list[str], list[str], str | None]] = {}

    for r in read_csv(RAW / "mrclinkdb_human_mr.txt", "\t"):
        sym, uni = r.get("Receptor_symbol"), r.get("Receptor uniprot_ id") or r.get("Receptor uniprot_id")
        if sym:
            tokens[("mrclinkdb_human_mr", sym.strip())] = ([sym.strip()], [uni.strip()] if uni else [], None)

    # MetaLigand encodes receptor COMPLEXES with semicolons - 19 distinct tokens
    # over 29 pair rows, e.g. 'GABRA1;GABRG2;GABRB2' and the five-subunit
    # 'KCTD16;KCTD8;KCTD12;GABBR1;GABBR2'. Taken literally these would be created
    # as entities whose gene symbol is a semicolon-joined string.
    for r in read_csv(RAW / "metaligand_NPLRdb_human.csv", ","):
        raw = (r.get("R") or "").strip()
        if not raw:
            continue
        subs = [x.strip() for x in raw.split(";") if x.strip()]
        tokens[("metaligand_pairs_human", raw)] = (subs, [], raw if len(subs) > 1 else None)

    with sqlite3.connect(RAW / "metalinksdb.db") as s:
        for (uni,) in s.execute("SELECT DISTINCT uniprot FROM edges WHERE type='lr'"):
            if uni:
                tokens[("metalinksdb", uni)] = ([uni2sym.get(uni, "")], [uni], None)

    for r in read_csv(RAW / "mebocost_human_met_sensor.tsv", "\t"):
        if r.get("Gene_name"):
            tokens[("mebocost_human", r["Gene_name"].strip())] = ([r["Gene_name"].strip()], [], None)

    f = (RAW / "gtopdb_interactions.csv").open(encoding="utf-8", errors="replace")
    next(f)
    for r in csv.DictReader(f):
        if not eqi(r["Endogenous"], "true") or not eqi(r["Target Species"], SPECIES):
            continue
        if not in_ci(r["Ligand Type"], LIGAND_TYPES):
            continue
        sym = (r.get("Target Gene Symbol") or "").strip()
        uni = (r.get("Target UniProt ID") or "").strip()
        key = sym or uni
        if key:
            # GtoPdb writes multi-subunit targets as space/pipe separated lists.
            syms = [x for x in sym.replace("|", " ").split() if x]
            unis = [x for x in uni.replace("|", " ").split() if x]
            tokens[("gtopdb_interactions", key)] = (syms, unis,
                                                    key if len(syms) > 1 else None)

    # CellChatDB: `receptor` is the complex name, `receptor.symbol` the subunits.
    sys.path.insert(0, str(ROOT / "src"))
    from rdata import read_rdata
    cc = read_rdata(ROOT / "data" / "raw" / "CellChatDB.human.rda")
    it = (cc.get("CellChatDB.human") or next(iter(cc.values())))["interaction"]
    for i in range(len(it["interaction_name"])):
        if not eqi(it["annotation"][i], "Non-protein Signaling"):
            continue
        rec = (it["receptor"][i] or "").strip()
        syms = [s.strip() for s in (it["receptor.symbol"][i] or "").split(",") if s.strip()]
        if rec:
            tokens[("cellchatdb_human", rec)] = (syms or [rec], [], rec if len(syms) > 1 else None)

    # CellPhoneDB: partner_b is a uniprot, or a complex named in complex_input.
    cpx: dict[str, list[str]] = {}
    for r in read_csv(ROOT / "data" / "raw" / "cellphonedb_complexes.csv", ","):
        subs = [r.get(f"uniprot_{k}") for k in range(1, 6)]
        cpx[r["complex_name"]] = [x.strip() for x in subs if x and x.strip()]
    for r in read_csv(ROOT / "data" / "raw" / "cellphonedb_interactions.csv", ","):
        if not eqi(r.get("is_ppi"), "false"):
            continue
        b = (r.get("partner_b") or "").strip()
        if not b:
            continue
        if b in cpx:
            tokens[("cellphonedb_interactions", b)] = ([], cpx[b], b)
        else:
            tokens[("cellphonedb_interactions", b)] = ([uni2sym.get(b, "")], [b], None)

    # --- resolve -----------------------------------------------------------
    rows: list[tuple] = []
    unresolved: list[dict] = []
    created: dict[str, str] = {}
    n_next = 1
    counts: dict[str, int] = defaultdict(int)

    def new_entity(symbol: str, uniprot: str | None) -> str:
        nonlocal n_next
        while True:
            eid = f"N{n_next:05d}"
            n_next += 1
            if c.execute("SELECT 1 FROM entities WHERE entity_id=?", (eid,)).fetchone() is None:
                break
        c.execute("INSERT INTO entities (entity_id,kind,species,symbol,uniprot,is_receptor)"
                  " VALUES (?,?,?,?,?,1)", (eid, "protein", SPECIES, symbol, uniprot))
        by_symbol[norm_sym(symbol)] = eid
        if uniprot:
            by_uniprot.setdefault(uniprot, eid)
        created[norm_sym(symbol)] = eid
        return eid

    for (src, raw), (syms, unis, complex_name) in sorted(tokens.items()):
        parts = [(s, None) for s in syms if s] or [(None, u) for u in unis if u]
        if len(syms) > 1 and unis and len(unis) == len(syms):
            parts = list(zip(syms, unis))
        elif not syms and unis:
            parts = [(uni2sym.get(u) or None, u) for u in unis]
        if not parts:
            unresolved.append(dict(source_id=src, raw_receptor=raw, reason="no symbol or uniprot"))
            counts["unresolved"] += 1
            continue
        multi = len(parts) > 1 or bool(complex_name)
        for sym, uni in parts:
            eid = res = None
            if sym and norm_sym(sym) in by_symbol:
                eid, res = by_symbol[norm_sym(sym)], "matched_symbol"
            elif uni and uni in by_uniprot:
                eid, res = by_uniprot[uni], "matched_uniprot"
            elif sym:
                eid, res = new_entity(sym, uni), "created"
            elif uni:
                eid, res = new_entity(uni2sym.get(uni) or uni, uni), "created"
            if eid is None:
                unresolved.append(dict(source_id=src, raw_receptor=raw,
                                       reason="unmappable token"))
                counts["unresolved"] += 1
                continue
            counts[res] += 1
            rows.append((src, raw, eid, sym, uni,
                         complex_name if multi else None,
                         "subunit" if multi else None,
                         rtype.get(norm_sym(sym)) if sym else None, res))

    # Guard before writing: a source that contributes nothing here contributes
    # nothing to np4 either, and np3 is where the loss would start.
    per_source: dict[str, int] = defaultdict(int)
    for r in rows:
        per_source[r[0]] += 1
    if guard_contributions(c, "np3", per_source, NP_SOURCES, commit):
        c.rollback()
        return 1

    c.execute("DELETE FROM np_receptor_map")
    c.executemany("INSERT OR IGNORE INTO np_receptor_map (source_id,raw_receptor,entity_id,"
                  "symbol,uniprot,complex_id,complex_role,receptor_type,resolution)"
                  " VALUES (?,?,?,?,?,?,?,?,?)", rows)

    if unresolved:
        with UNRESOLVED.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["source_id", "raw_receptor", "reason"],
                               delimiter="\t")
            w.writeheader()
            w.writerows(unresolved)

    if commit:
        c.commit()
    else:
        c.rollback()

    # --- report -------------------------------------------------------------
    log("np3", f"receptor tokens: {len(tokens)} -> np_receptor_map rows: {len(rows)}")
    for k in ("matched_symbol", "matched_uniprot", "created", "unresolved"):
        log("np3", f"  {k:18s} {counts[k]:5d}")
    log("np3", f"new entities created (N##### ids): {len(created)}")
    ncx = len({r[5] for r in rows if r[5]})
    log("np3", f"complexes exploded: {ncx} -> {sum(1 for r in rows if r[5])} subunit rows")
    typed = sum(1 for r in rows if r[7])
    log("np3", f"receptor_type filled: {typed}/{len(rows)} from MetalinksDB protein_type")

    # MEBOCOST cross-check only; it never supplies a shipped value.
    agree = dis = 0
    # The two taxonomies OVERLAP rather than conflict: a nicotinic receptor is
    # both a receptor and a ligand-gated ion channel, and a guanylate cyclase
    # receptor is both a receptor and an enzyme. Each MetalinksDB type therefore
    # maps to the SET of MEBOCOST labels that are compatible with it, not one.
    fam = {"gpcr": {"Receptor"},
           "catalytic_receptor": {"Receptor", "Enzyme"},
           "nhr": {"Nuclear Receptor", "Receptor"},
           "lgic": {"Channel", "Receptor"},
           "vgic": {"Channel", "Receptor"},
           "other_ic": {"Channel", "Receptor"},
           "transporter": {"Transporter"},
           "enzyme": {"Enzyme"}}
    for r in rows:
        if not (r[7] and r[3]):
            continue
        m = mebo_type.get(norm_sym(r[3]))
        if not m:
            continue
        if any(lab in m for lab in fam.get(r[7], set())):
            agree += 1
        else:
            dis += 1
    if agree + dis:
        log("np3", f"receptor_type vs MEBOCOST (internal cross-check): "
                   f"{agree} agree, {dis} disagree ({agree/(agree+dis):.0%})")
    if unresolved:
        log("np3", f"unresolved: {len(unresolved)} -> {UNRESOLVED.name}")
    log("np3", "COMMITTED" if commit else "DRY RUN - rerun with --commit to write")
    return 0


if __name__ == "__main__":
    sys.exit(main())
