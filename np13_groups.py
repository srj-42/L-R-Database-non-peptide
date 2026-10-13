#!/usr/bin/env python3
"""np13_groups.py - collapse component rows of a multi-subunit receptor into one (chemical, complex) group.

Owner rule 2026-10-12: where a complex is split into its components in the table, search for evidence once
for the named complex; a paper naming the complex is enough for every component row.
Writes work/np13_groups.json:  {group_key: {label, chemical_id, chemical, pairs:[cid..], symbols:[..], done:[cid..]}}
"""
import json, re, sqlite3, collections
DB = "/home/sreeja/omics/database/curation/curation.db"
FAM = [  # (regex on HGNC symbol, label, regex a paper uses for the complex)
    (r"^GABR[ABGDEPQ]\d*$", "GABA-A receptor", r"GABA\s*\(?A\)?|GABAA|GABA-A"),
    (r"^GABRR\d$", "GABA-C (rho) receptor", r"GABA\s*\(?C\)?|GABAC|rho"),
    (r"^CHRN[ABDEG]\d*$", "nicotinic acetylcholine receptor", r"nicotinic|nAChR"),
    (r"^GRIN\d[A-D]?$", "NMDA receptor", r"NMDA|GluN|NR[123]"),
    (r"^GRIA\d$", "AMPA receptor", r"AMPA|GluA|GluR"),
    (r"^GRIK\d$", "kainate receptor", r"kainate|GluK|GluR"),
    (r"^GLR[AB]\d?$", "glycine receptor", r"glycine receptor|GlyR"),
    (r"^HTR3[A-E]$", "5-HT3 receptor", r"5-HT\(?3\)?|5-HT3"),
    (r"^P2RX\d$", "P2X receptor", r"P2X"),
    (r"^GABBR\d$", "GABA-B receptor", r"GABA\s*\(?B\)?|GABAB"),
    (r"^SCNN1[ABDG]$", "ENaC", r"ENaC|epithelial sodium channel|epithelial Na"),
    (r"^CATSPER", "CatSper channel", r"CatSper"),
    (r"^(KCNJ8|KCNJ11|ABCC8|ABCC9)$", "KATP channel", r"K\(?ATP\)?|KATP|ATP-sensitive"),
    (r"^ITG[AB]\d+$", "integrin", r"integrin"),
    (r"^SCN\d+[AB]$", "voltage-gated sodium channel", r"sodium channel|Nav"),
    (r"^CACN[ABG]", "voltage-gated calcium channel", r"calcium channel|Cav"),
    (r"^RYR\d$", "ryanodine receptor", r"ryanodine"),
    (r"^ITPR\d$", "IP3 receptor", r"IP3|InsP3|inositol 1,4,5"),
]
def family(sym):
    for rx, lab, pap in FAM:
        if re.match(rx, sym or ""): return lab, pap
    return None, None
def build():
    c = sqlite3.connect(DB); c.row_factory = sqlite3.Row
    rows = c.execute("""SELECT ci.chem_interaction_id cid, ci.chemical_id, ch.name, e.symbol, ci.shippable
        FROM chemical_interactions ci JOIN chemical_entities ch USING(chemical_id)
        JOIN entities e ON e.entity_id=ci.receptor_entity_id WHERE ci.species='human'""").fetchall()
    ev = collections.defaultdict(list)
    for r in c.execute("SELECT chem_interaction_id, verification_outcome FROM np_evidence"):
        ev[r[0]].append(r[1])
    done = set()
    for l in open("/home/sreeja/omics/database/curation/work/np10_search_log.jsonl"):
        try: d = json.loads(l); done.add(d.get("chem_interaction_id") or d.get("pair"))
        except Exception: pass
    G = {}
    for r in rows:
        lab, pap = family(r["symbol"])
        if not lab: continue
        k = f"{r['chemical_id']}|{lab}"
        g = G.setdefault(k, {"label": lab, "paper_rx": pap, "chemical_id": r["chemical_id"], "chemical": r["name"],
                             "pairs": [], "symbols": [], "answered": [], "evidenced": []})
        g["pairs"].append(r["cid"]); g["symbols"].append(r["symbol"])
        if r["cid"] in done or ev.get(r["cid"]): g["answered"].append(r["cid"])
        if any(v in ("confirmed", "pair_absent") for v in ev.get(r["cid"], [])): g["evidenced"].append(r["cid"])
    G = {k: v for k, v in G.items() if len(v["pairs"]) > 1}
    json.dump(G, open("/home/sreeja/omics/database/curation/work/np13_groups.json", "w"), indent=0)
    n_pairs = sum(len(v["pairs"]) for v in G.values())
    open_groups = [g for g in G.values() if len(g["answered"]) < len(g["pairs"])]
    print(f"{len(G)} complex groups covering {n_pairs} pairs; {len(open_groups)} groups have unanswered pairs "
          f"({sum(len(g['pairs'])-len(g['answered']) for g in open_groups)} pairs); "
          f"{sum(1 for g in G.values() if g['evidenced'])} groups already hold evidence")
if __name__ == "__main__": build()
