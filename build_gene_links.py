#!/usr/bin/env python3
"""Emit data/genes.js for the four site copies: per-gene identifiers + aliases used by the drawer's
"Gene links" rows (connectomeDB-style card). Link-out only: nothing here enters a score.
Source: entities table, gaps filled from raw/hgnc/hgnc_complete_set.txt (HGNC id, Ensembl, UniProt, MGI
of the mouse orthologue) and work/np19/uniprot_map.json. Key 'h:SYMBOL' / 'm:Symbol'.
Row: [uniprot, hgnc_id, ensembl_id, aliases[<=8], mgi_id_of_orthologue]
Usage: spice/.venv/bin/python build_gene_links.py"""
import csv, json, sqlite3
from pathlib import Path
HERE = Path(__file__).resolve().parent; OMICS = HERE.parent.parent
csv.field_size_limit(10**7)
hg = {}
with open(HERE / "raw/hgnc/hgnc_complete_set.txt", newline="", encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["status"] == "Approved": hg[r["symbol"]] = r
umap = json.loads((HERE / "work/np19/uniprot_map.json").read_text())
c = sqlite3.connect(f"file:{HERE/'curation.db'}?mode=ro", uri=True)
out = {}; fill = 0
for sym, sp, up, hid, ens, al in c.execute(
        "SELECT symbol,species,uniprot,hgnc_id,ensembl_id,aliases FROM entities WHERE kind='protein'"):
    k = ("h:" if sp == "human" else "m:") + sym
    mgi = ""
    if sp == "human":
        h = hg.get(sym)
        if h:
            if not up: up = h["uniprot_ids"].split("|")[0]
            hid = hid or h["hgnc_id"]; ens = ens or h["ensembl_gene_id"]; mgi = h["mgd_id"].split("|")[0]
        if not up and umap.get(sym): up = umap[sym][0]
    a = [x for x in (json.loads(al) if al else []) if x != sym][:8]
    out[k] = [up or "", hid or "", ens or "", a, mgi]
js = "window.GENES=" + json.dumps(out, separators=(",", ":"), ensure_ascii=False) + ";\n"
for d in ["Peptide-site", "database/site", "non-peptide-site", "database/site_chemical"]:
    (OMICS / d / "data/genes.js").write_text(js)
h = [v for k, v in out.items() if k[0] == "h"]
print(len(out), "genes; human with uniprot", sum(1 for v in h if v[0]), "hgnc", sum(1 for v in h if v[1]),
      "ensembl", sum(1 for v in h if v[2]), "mgi", sum(1 for v in h if v[4]), "/", len(h))
