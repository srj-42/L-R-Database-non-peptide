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
import json, re
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
    # HMDB context per metabolite (np20): x = {t: tissue, c: cell location,
    # b: biospecimen, d: disease}. Short keys so they do not clash with "cl" (class).
    KEYMAP = {"tissue": "t", "cell_location": "c", "biospecimen": "b", "disease": "d"}
    if (REL / "chemical_context.csv").exists():
        for r in read("chemical_context.csv"):
            ch = chems.get(r["chemical_id"])
            if ch is not None and r["kind"] in KEYMAP:
                ch.setdefault("x", {}).setdefault(KEYMAP[r["kind"]], []).append(r["value"])
    # Pathway families (x.p) and the non-redundant ontology branches (np21):
    # x.s origin, x.r route of exposure, x.l roles by group. Location, health-condition
    # and pathway branches of the ontology repeat data shown elsewhere and are left out.
    if (REL / "chemical_pathways.csv").exists():
        for r in read("chemical_pathways.csv"):
            ch = chems.get(r["chemical_id"])
            if ch is not None:
                ch.setdefault("x", {}).setdefault("p", []).append(r["pathway"])
    if (REL / "chemical_ontology.csv").exists():
        for r in read("chemical_ontology.csv"):
            ch = chems.get(r["chemical_id"])
            if ch is None:
                continue
            path, term, lvl = r["path"], r["term"], int(r["level"])
            x = ch.setdefault("x", {})
            if path.startswith("Disposition > Source > ") and lvl in (3, 4):
                x.setdefault("s", []).append(term)
            elif path.startswith("Disposition > Route of exposure > "):
                x.setdefault("r", []).append(term)
            elif path.startswith("Role > ") and r["node_type"] == "child" and path.count(" > ") >= 2:
                x.setdefault("l", {}).setdefault(path.split(" > ")[1], []).append(term)
    (SITE / "data" / "chemicals.js").write_text(
        "window.CHEMICALS=" + json.dumps(chems, separators=(",", ":")) + ";")
    print(f"chemicals: {len(chems)}")

    # ---- interactions: compact arrays ------------------------------------
    SRC_ORDER = ["gtopdb_interactions", "cellchatdb_human", "mrclinkdb_human_mr",
                 "cellphonedb_interactions", "metaligand_pairs_human", "metalinksdb"]
    sidx = {s: i for i, s in enumerate(SRC_ORDER)}
    # PMIDs per pair = upstream citations + quoted papers + database-annotation PMIDs, so the table
    # column agrees with what the drawer lists (2026-11-03: it used to count upstream citations only).
    pm_all = defaultdict(set)
    for r in read("chemical_source_records.csv"):
        for p in re.findall(r"(?<!\d)\d{7,8}(?!\d)", r.get("raw_reference") or ""):
            pm_all[r["chem_interaction_id"]].add(p)
    for r in read("chemical_evidence.csv"):
        if r.get("pmid"): pm_all[r["chem_interaction_id"]].add(r["pmid"])
    if (REL / "chemical_class_evidence.csv").exists():
        for r in read("chemical_class_evidence.csv"):
            if r.get("pmid"): pm_all[r["chem_interaction_id"]].add(r["pmid"])
    rows = []
    for r in read("chemical_interactions.csv"):
        srcs = json.loads(r["sources"] or "[]")
        rows.append([
            r["chem_interaction_id"], r["chemical_id"], r["chemical_name"],
            r["receptor_symbol"] or "", r["receptor_type"] or "",
            int(r["n_sources"] or 0),
            "".join(str(sidx[s]) for s in srcs if s in sidx),
            len(pm_all[r["chem_interaction_id"]]),
            r["complex_id"] or "",
        ])
    (SITE / "data" / "interactions.js").write_text(
        "window.INTERACTIONS=" + json.dumps(rows, separators=(",", ":")) +
        ";window.SRC_ORDER=" + json.dumps(SRC_ORDER) + ";")
    print(f"interactions: {len(rows)}")

    # ---- enzymes, grouped by chemical ------------------------------------
    # 2026-11-03 (owner): genes are normalised to the current HGNC symbol (sources mix old and new
    # symbols, e.g. AARS and AARS1, which showed as two genes), junk tokens are dropped, and each gene
    # gets a type tag: channel / transporter / enzyme (HGNC EC number) / other. The HGNC complete set
    # is pinned in raw/hgnc/hgnc_complete_set.txt.
    HG = {}; prev = defaultdict(set); alias = defaultdict(set)
    hp = Path(__file__).parent / "raw" / "hgnc" / "hgnc_complete_set.txt"
    if hp.exists():
        for h in csv.DictReader(open(hp, encoding="utf-8"), delimiter="\t"):
            HG[h["symbol"].upper()] = h
            for a in (h["prev_symbol"] or "").strip('"').split("|"):
                if a: prev[a.upper()].add(h["symbol"])
            for a in (h["alias_symbol"] or "").strip('"').split("|"):
                if a: alias[a.upper()].add(h["symbol"])
    JUNK = {"A", "C", "1", "CYP 20.00"}
    def norm_gene(g):
        out = []
        for p in re.split(r"\s+AND\s+", (g or "").replace("\xa0", " ").strip()):
            u = p.upper()
            if not p or u in JUNK: continue
            if u in HG: out.append(HG[u]["symbol"])
            elif len(prev.get(u, ())) == 1: out.append(next(iter(prev[u])))
            elif len(alias.get(u, ())) == 1: out.append(next(iter(alias[u])))
            else: out.append(p)
        return out
    ENZ_WORD = re.compile(r"\b(?!(?:base|release|disease|phase|case|increase|decrease)\b)[a-z0-9-]*ase\b|\benzyme\b")
    def gene_type(sym):
        h = HG.get(sym.upper())
        if not h: return "other"
        grp, nm = (h["gene_group"] or "").lower(), (h["name"] or "").lower()
        if "channel" in grp or " channel" in nm: return "channel"
        if (sym.startswith(("SLC", "ABC")) or any(k in grp for k in ("solute carrier", "atp binding cassette", "transport"))
                or "transporter" in nm or "transporting" in nm): return "transporter"
        if h["enzyme_id"] or ENZ_WORD.search(nm): return "enzyme"     # HGNC's EC field is sparse: also use '-ase' names
        return "other"
    enz = defaultdict(list); ginfo = {}; seen = set(); hubs = set()
    for e in read("metabolite_enzymes.csv"):
        for g in norm_gene(e["gene_symbol"]):
            k = (e["chemical_id"], g, e["role"], e["source_id"])
            if k in seen: continue
            seen.add(k)
            enz[e["chemical_id"]].append([g, e["role"], e["source_id"], int(e["direction_known"] or 1),
                                          e.get("also_listed_by") or "", e.get("origin") or ""])

            if g not in ginfo:
                h = HG.get(g.upper())
                ginfo[g] = [gene_type(g), (h["name"] if h else "")]
    if (REL / "chemical_gene_hubs.csv").exists():
        hubs |= {r["chemical_id"] for r in read("chemical_gene_hubs.csv")}
    # Root-source genes (np22): GENEROOT[cid][gene] = [[role, srcmask, directed, rhea_exp], ...]
    # srcmask letters: G Human-GEM, R Rhea/UniProt, T TCDB, H HMDB, X Reactome, O GO (experimental). Genes also get a GENEINFO
    # type tag so root-only genes can be shown with the same pills.
    groot = defaultdict(lambda: defaultdict(list))
    SRCL = {"human_gem": "G", "rhea_uniprot": "R", "tcdb": "T", "hmdb": "H", "reactome": "X", "go_exp": "O"}
    if (REL / "chemical_gene_roots.csv").exists():
        for r in read("chemical_gene_roots.csv"):
            g = r["gene_symbol"]
            groot[r["chemical_id"]][g].append([r["role"], "".join(SRCL.get(x, "?") for x in r["sources"].split("|")),
                                              int(r["has_directed_reaction"]), int(r["rhea_experimental"])])
            if g not in ginfo:
                h = HG.get(g.upper())
                ginfo[g] = [gene_type(g), (h["name"] if h else "")]
    (SITE / "data" / "enzymes.js").write_text(
        "window.ENZYMES=" + json.dumps(enz, separators=(",", ":")) + ";window.GENEINFO=" +
        json.dumps(ginfo, separators=(",", ":")) + ";window.GENEROOT=" +
        json.dumps({k: dict(v) for k, v in groot.items()}, separators=(",", ":")) +
        ";window.HUBS=" + json.dumps(sorted(hubs)) + ";")
    print(f"root genes: {sum(len(v) for v in groot.values())} chemical-gene links over {len(groot)} chemicals")
    print("gene types:", dict(sorted(defaultdict(int, {t: sum(1 for v in ginfo.values() if v[0] == t) for t in ("enzyme", "transporter", "channel", "other")}).items())))
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

    # ---- non-literature evidence classes (CITATION_METHODS.md Methods 3, 4) --
    cls = defaultdict(list)
    cls_path = REL / "chemical_class_evidence.csv"
    for r in (read("chemical_class_evidence.csv") if cls_path.exists() else []):
        cls[r["chem_interaction_id"]].append({
            "cls": r["evidence_class"], "conf": r["confidence"], "src": r["source"],
            "acc": r["accession"], "det": r["detail"], "pmid": r["pmid"]})
    (SITE / "data" / "classev.js").write_text(
        "window.CLASSEV=" + json.dumps(cls, separators=(",", ":")) + ";")
    print(f"class evidence rows: {sum(len(v) for v in cls.values())} over {len(cls)} pairs")

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
    vset = {r["chem_interaction_id"] for r in read("chemical_evidence.csv")
            if r["verification_outcome"] == "confirmed"} | set(cls)
    ids_in_site = {r[0] for r in rows}
    src_pairs = defaultdict(set)
    for r in read("chemical_source_records.csv"):
        if r["chem_interaction_id"] in ids_in_site: src_pairs[r["source_id"]].add(r["chem_interaction_id"])
    # "own_verified": pairs a resource ALREADY verifies by this project's rules, using only what the resource
    # itself supplies: its own cited paper yielded a confirmed quoted sentence (np_evidence.source_id names the
    # resource), or the resource's record is itself a curated-database record (GtoPdb). Read-only DB query.
    import sqlite3
    db = sqlite3.connect(f"file:{Path(__file__).parent / 'curation.db'}?mode=ro", uri=True)
    own = {}
    for s_ in SRC_ORDER:
        o = {r[0] for r in db.execute("SELECT DISTINCT chem_interaction_id FROM np_evidence WHERE verification_outcome='confirmed' "
                                      "AND source_id LIKE ?", (f"%{s_}%",))}
        if s_ == "gtopdb_interactions":
            o |= {r[0] for r in db.execute("SELECT DISTINCT chem_interaction_id FROM np_class_evidence WHERE source='GtoPdb'")}
        own[s_] = o & src_pairs[s_]
    db.close()
    by_source = {s: {"pairs": len(src_pairs[s]), "own_verified": len(own[s]), "verified": len(src_pairs[s] & vset)} for s in SRC_ORDER}
    db_src = defaultdict(set)
    for cid, lst in cls.items():
        for g in lst: db_src[g["src"]].add(cid)
    (SITE / "data" / "verified.js").write_text("window.VERIFIED=" + json.dumps(sorted(vset & ids_in_site), separators=(",", ":")) + ";")
    stats = {
        "by_source": by_source,
        "by_dbsource": {k: len(v) for k, v in sorted(db_src.items(), key=lambda x: -len(x[1]))},
        "n_lit_verified": len({r["chem_interaction_id"] for r in read("chemical_evidence.csv") if r["verification_outcome"] == "confirmed"}),
        "n_evidence": sum(len(v) for v in ev.values()),
        "n_pairs_with_evidence": len(ev),
        "n_verified": len({r["chem_interaction_id"] for r in read("chemical_evidence.csv")
                           if r["verification_outcome"] == "confirmed"} | set(cls)),
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
.gcount{font-size:9px;margin-left:2px;color:#57606a}
.gtype{font-size:9px;margin-left:4px;padding:0 4px;border-radius:6px;background:#fff;border:1px solid #d0d7de;color:#57606a}
.gt-transporter{border-color:#aac4f5}.gt-channel{border-color:#d6b4f5}.gt-enzyme{border-color:#a8d5ba}
.donuts{display:flex;flex-wrap:wrap;gap:18px 26px;margin:6px 0 4px}
.donut{margin:0;text-align:center;font-size:12.5px;width:120px}
.donut figcaption span{color:var(--muted);font-size:11.5px}
.donut-key{display:flex;gap:18px;font-size:12px;color:var(--muted);margin-top:8px}
.donut-key i{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:5px;vertical-align:-1px}
table.trust{border-collapse:collapse;width:100%;font-size:12.5px}
table.trust th,table.trust td{border-bottom:1px solid var(--line,#e5e5e0);padding:6px 8px;text-align:left;vertical-align:top}
table.trust th{color:var(--muted);font-weight:600}
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
