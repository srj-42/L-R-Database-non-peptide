#!/usr/bin/env python3
"""Build the static browse site for the chemical half, from data/release_v2/.

Same construction as cur7_build_site.py: no server, no framework, no CDN, data
emitted as .js files that assign to globals and load via <script src> rather
than fetch(), so the page works by double-clicking index.html as well as behind
`python3 -m http.server`.

It is a SEPARATE site from the peptide one (decision 1), not a tab inside it.
The two halves are not comparable - the chemical rows have no confidence score -
and putting them in one table would invite exactly the comparison the report
argues against. The styling is shared so they read as siblings, and each links
to the other.

Since 2026-10-09 the chemical half DOES carry verified evidence for a minority of
pairs (chemical_evidence.csv, from np10_litsearch.py), so the site shows the
quoted sentence and its paper in the pair drawer, and a "cited" filter. The rest
of the caveat stands: those rows are machine-verified against the paper's own
text, not curator-read, and 13,307 of 13,440 pairs still have no evidence at all.

Provenance is the large table (21k rows), so it is sharded by pair and injected
on demand when a row is opened.

Usage: python3 np7_build_site.py
"""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REL = ROOT / "data" / "release_v2"
SITE = ROOT / "site_chemical"
PEPTIDE_SITE = ROOT / "site"
SHARDS = 16

csv.field_size_limit(10_000_000)


def read(name):
    with (REL / name).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def shard_of(key: str) -> int:
    return sum(ord(c) for c in key) % SHARDS


def main() -> int:
    (SITE / "data" / "prov").mkdir(parents=True, exist_ok=True)
    manifest = json.loads((REL / "chemical_manifest.json").read_text())

    # ---- chemicals -------------------------------------------------------
    chems = {}
    for r in read("chemical_entities.csv"):
        try:
            syn = json.loads(r["synonyms"]) if r["synonyms"] else []
        except Exception:                                    # noqa: BLE001
            syn = []
        chems[r["chemical_id"]] = {
            "n": r["name"], "h": r["hmdb_id"], "pc": r["pubchem_cid"],
            "sm": r["smiles"], "f": r["formula"],
            "k": r["kingdom"], "sc": r["super_class"], "cl": r["class"],
            "cc": r["chebi_class"], "ms": r["metabolite_subclass"],
            "mc": r["molecular_class"], "nt": int(r["is_neurotransmitter"] or 0),
            "src": json.loads(r["sources"] or "[]"),
            "pv": int(r["provisional_identity"] or 0),
            "syn": syn[:12],
        }
    (SITE / "data" / "chemicals.js").write_text(
        "window.CHEMICALS=" + json.dumps(chems, separators=(",", ":")) + ";")
    print(f"chemicals: {len(chems)}")

    # ---- interactions: compact arrays ------------------------------------
    SRC_ORDER = ["gtopdb_interactions", "cellchatdb_human", "mrclinkdb_human_mr",
                 "cellphonedb_interactions", "metaligand_pairs_human", "metalinksdb"]
    sidx = {s: i for i, s in enumerate(SRC_ORDER)}
    rows = []
    for r in read("chemical_interactions.csv"):
        srcs = json.loads(r["sources"] or "[]")
        rows.append([
            r["chem_interaction_id"], r["chemical_id"], r["chemical_name"],
            r["receptor_symbol"] or "", r["receptor_type"] or "",
            int(r["n_sources"] or 0),
            "".join(str(sidx[s]) for s in srcs if s in sidx),
            int(r["n_pmids"] or 0),
            r["complex_id"] or "",
        ])
    (SITE / "data" / "interactions.js").write_text(
        "window.INTERACTIONS=" + json.dumps(rows, separators=(",", ":")) +
        ";window.SRC_ORDER=" + json.dumps(SRC_ORDER) + ";")
    print(f"interactions: {len(rows)}")

    # ---- enzymes, grouped by chemical ------------------------------------
    enz = defaultdict(list)
    for e in read("metabolite_enzymes.csv"):
        enz[e["chemical_id"]].append([e["gene_symbol"], e["role"], e["source_id"],
                                      int(e["direction_known"] or 1)])
    (SITE / "data" / "enzymes.js").write_text(
        "window.ENZYMES=" + json.dumps(enz, separators=(",", ":")) + ";")
    print(f"enzyme rows: {sum(len(v) for v in enz.values())} over {len(enz)} chemicals")

    # ---- confirmed evidence, grouped by pair -----------------------------
    # Small enough (153 rows) to ship in one file rather than shard like provenance.
    ev = defaultdict(list)
    for r in read("chemical_evidence.csv"):
        ev[r["chem_interaction_id"]].append({
            "pmid": r["pmid"], "yr": r["paper_year"], "jr": r["paper_journal"],
            "res": r["result"], "as": r["assay_class"], "dir": r["effect_direction"],
            "aff": r["affinity"], "mol": r["molecule_as_named_in_paper"],
            "rec": r["receptor_as_named_in_paper"], "sent": r["evidence_sentence"],
            "dep": r["evidence_depth"], "xm": r["extraction_method"],
        })
    (SITE / "data" / "evidence.js").write_text(
        "window.EVIDENCE=" + json.dumps(ev, separators=(",", ":")) + ";")
    print(f"evidence rows: {sum(len(v) for v in ev.values())} over {len(ev)} pairs")

    # ---- provenance, sharded by pair -------------------------------------
    buckets = defaultdict(lambda: defaultdict(list))
    n_prov = 0
    for p in read("chemical_source_records.csv"):
        cid = p["chem_interaction_id"]
        if not cid:
            continue
        n_prov += 1
        buckets[shard_of(cid)][cid].append([
            p["source_id"], p["raw_ligand"] or "", p["raw_receptor"] or "",
            p["raw_reference"] or ""])
    for s in range(SHARDS):
        (SITE / "data" / "prov" / f"p{s}.js").write_text(
            f"window.__PROV[{s}]=" + json.dumps(buckets.get(s, {}), separators=(",", ":")) + ";")
    print(f"provenance rows: {n_prov} in {SHARDS} shards")

    # ---- stats -----------------------------------------------------------
    by_role = defaultdict(int)
    for v in enz.values():
        for g in v:
            by_role[g[1]] += 1
    by_class = defaultdict(int)
    by_rtype = defaultdict(int)
    for r in rows:
        by_rtype[r[4] or "(untyped)"] += 1
    for ch in chems.values():
        by_class[ch["mc"] or "(unclassified)"] += 1
    cited = sum(1 for r in rows if r[7])
    by_result = defaultdict(int)
    for v in ev.values():
        for g in v:
            by_result[g["res"] or "(unlabelled)"] += 1
    stats = {
        "n_evidence": sum(len(v) for v in ev.values()),
        "n_pairs_with_evidence": len(ev),
        "n_evidence_papers": len({g["pmid"] for v in ev.values() for g in v}),
        "by_result": dict(sorted(by_result.items(), key=lambda x: -x[1])),
        "manifest": manifest,
        "shards": SHARDS,
        "n_chem": len(chems), "n_inter": len(rows),
        "n_enz": sum(len(v) for v in enz.values()),
        "n_chem_with_enz": len(enz),
        "cited": cited,
        "by_role": dict(sorted(by_role.items(), key=lambda x: -x[1])),
        "by_class": dict(sorted(by_class.items(), key=lambda x: -x[1])),
        "by_rtype": dict(sorted(by_rtype.items(), key=lambda x: -x[1])[:12]),
        "by_nsrc": {str(k): sum(1 for r in rows if r[5] == k) for k in range(1, 8)},
    }
    (SITE / "data" / "stats.js").write_text(
        "window.STATS=" + json.dumps(stats, separators=(",", ":")) + ";")

    # ---- shared stylesheet ------------------------------------------------
    css = (PEPTIDE_SITE / "styles.css").read_text()
    extra = """
/* --- chemical-half additions ------------------------------------------- */
.prov.chem{background:#4c1d95;}
.pill{display:inline-block;padding:1px 7px;border-radius:9px;font-size:11px;
      border:1px solid #d0d7de;margin-right:4px;white-space:nowrap}
.pill.syn{background:#e7f5ec;border-color:#a8d5ba}
.pill.deg{background:#fdecea;border-color:#f3b5ae}
.pill.eff{background:#e8f0fe;border-color:#aac4f5}
.pill.inf{background:#fff4e5;border-color:#f5cfa0}
.pill.pre{background:#f3e8ff;border-color:#d6b4f5}
.warn{background:#fff8e1;border:1px solid #f0d48a;padding:10px 12px;border-radius:6px;
      margin:0 0 14px;font-size:13px}
.smiles{font-family:ui-monospace,Menlo,monospace;font-size:11px;color:#57606a;
        word-break:break-all}
.srcbar{display:flex;gap:2px;align-items:center}
.srcbar i{width:9px;height:9px;border-radius:2px;background:#d0d7de;display:inline-block}
.srcbar i.on{background:#0969da}
.srcbar i.on.q{background:#1a7f37}
"""
    (SITE / "styles.css").write_text(css + extra)

    print(f"site data written to {SITE/'data'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
