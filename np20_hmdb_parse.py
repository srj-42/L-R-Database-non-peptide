#!/usr/bin/env python3
"""np20_hmdb_parse.py - stream-parse raw/hmdb_metabolites.zip (HMDB 5.0 'All Metabolites')
without unpacking it and cache context annotations for the chemicals we hold.

Output: work/np20/hmdb_context.json  {our_hmdb_id: {tissue:[], cell_location:[],
        biospecimen:[], disease:[], pathways:[[name,smpdb,kegg]], ontology:[[root,path,term,level,type]],
        hmdb_accession:str}}
Matches our chemical_entities.hmdb_id against the primary accession OR any
secondary accession (retired ids), so old ids still resolve.
Run:  python3 -I np20_hmdb_parse.py
"""
import json, sqlite3, zipfile, os, time
import xml.etree.ElementTree as ET

ZIP = "raw/hmdb_metabolites.zip"
OUT = "work/np20/hmdb_context.json"
NS = "{http://www.hmdb.ca}"

c = sqlite3.connect("file:curation.db?mode=ro", uri=True)
ours = {r[0] for r in c.execute("SELECT hmdb_id FROM chemical_entities WHERE hmdb_id IS NOT NULL AND hmdb_id!=''")}
print("our HMDB ids:", len(ours), flush=True)

def texts(el, path):
    out = []
    for x in el.findall(path):
        t = (x.text or "").strip()
        if t:
            out.append(t)
    return sorted(set(out))

res, n, t0 = {}, 0, time.time()
with zipfile.ZipFile(ZIP) as z, z.open(z.namelist()[0]) as f:
    for ev, el in ET.iterparse(f, events=("end",)):
        if el.tag != NS + "metabolite":
            continue
        n += 1
        acc = (el.findtext(NS + "accession") or "").strip()
        sec = texts(el, f"{NS}secondary_accessions/{NS}accession")
        hit = [a for a in [acc] + sec if a in ours]
        if hit:
            bp = el.find(NS + "biological_properties")
            rec = {"hmdb_accession": acc, "tissue": [], "cell_location": [], "biospecimen": [], "disease": []}
            if bp is not None:
                rec["cell_location"] = texts(bp, f"{NS}cellular_locations/{NS}cellular")
                rec["biospecimen"] = texts(bp, f"{NS}biospecimen_locations/{NS}biospecimen")
                rec["tissue"] = texts(bp, f"{NS}tissue_locations/{NS}tissue")
            rec["disease"] = texts(el, f"{NS}diseases/{NS}disease/{NS}name")
            # pathways (SMPDB / KEGG ids) and the full ontology tree, flattened to
            # [root, path, term, level, type] rows (added 2026-11-06)
            rec["pathways"] = []
            for pw in el.findall(f"{NS}biological_properties/{NS}pathways/{NS}pathway"):
                nm = (pw.findtext(NS + "name") or "").strip()
                if nm:
                    rec["pathways"].append([nm, (pw.findtext(NS + "smpdb_id") or "").strip(),
                                            (pw.findtext(NS + "kegg_map_id") or "").strip()])
            rec["ontology"] = []
            def walk(e, path):
                term = (e.findtext(NS + "term") or "").strip()
                if not term:
                    return
                p = path + [term]
                rec["ontology"].append([p[0], " > ".join(p), term,
                                        int(e.findtext(NS + "level") or len(p)),
                                        (e.findtext(NS + "type") or "").strip()])
                d = e.find(NS + "descendants")
                if d is not None:
                    for x in d:
                        walk(x, p)
            ont = el.find(NS + "ontology")
            if ont is not None:
                for root in ont:
                    walk(root, [])
            # ChEBI id and HMDB protein associations (added 2026-11-07 for np22):
            # [gene_name, uniprot_id, protein_type, hmdbp_accession]. Undirected links.
            rec["chebi_id"] = (el.findtext(NS + "chebi_id") or "").strip()
            rec["proteins"] = [[(p.findtext(NS + "gene_name") or "").strip(),
                                (p.findtext(NS + "uniprot_id") or "").strip(),
                                (p.findtext(NS + "protein_type") or "").strip(),
                                (p.findtext(NS + "protein_accession") or "").strip()]
                               for p in el.findall(f"{NS}protein_associations/{NS}protein")]
            for a in hit:
                res[a] = rec
        el.clear()
        if n % 20000 == 0:
            print(f"  {n} metabolites, {len(res)} matched, {time.time()-t0:.0f}s", flush=True)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
json.dump(res, open(OUT, "w"))
print(f"done: {n} HMDB metabolites read, {len(res)} of {len(ours)} of our ids matched -> {OUT}")
