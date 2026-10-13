#!/usr/bin/env python3
"""np13_chembl.py - structured-source pull: ChEMBL activities (with PubMed ids) for our chemicals x our receptors.
Writes work/np13_chembl.jsonl: one line per (chemical_id, receptor symbol, pmid) with assay type/value.
Usage: np13_chembl.py [chemical_id ...]   (no args = all chemicals)
"""
import json, re, sqlite3, sys, time, urllib.request, urllib.parse, os
DB = "/home/sreeja/omics/database/curation/curation.db"
OUT = "/home/sreeja/omics/database/curation/work/np13_chembl.jsonl"
API = "https://www.ebi.ac.uk/chembl/api/data/"
def get(path, **q):
    u = API + path + ".json?" + urllib.parse.urlencode(q)
    for a in range(4):
        try:
            return json.load(urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "np13"}), timeout=60))
        except Exception:
            time.sleep(2 * (a + 1))
    return {}
def norm(s): return re.sub(r"[^a-z0-9]", "", (s or "").lower())
def molecule_id(name, syns):
    names = {norm(name)} | {norm(x) for x in syns if len(norm(x)) > 3}
    d = get("molecule/search", q=name, limit=10)
    for m in d.get("molecules", []):
        cands = {norm(m.get("pref_name"))} | {norm(s.get("molecule_synonym")) for s in m.get("molecule_synonyms", [])}
        if names & cands: return m["molecule_chembl_id"]
    return None
def main():
    c = sqlite3.connect(DB); c.row_factory = sqlite3.Row
    want = sys.argv[1:]
    chems = c.execute("select chemical_id,name,synonyms from chemical_entities").fetchall()
    if want: chems = [r for r in chems if r["chemical_id"] in want]
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT): done.add(json.loads(l).get("chemical_id"))
    up2sym = {}
    for r in c.execute("select symbol, uniprot from entities where uniprot is not null"):
        for u in str(r[1]).replace(";", ",").split(","):
            up2sym.setdefault(u.strip(), set()).add(r[0])
    tcache = {}
    with open(OUT, "a") as fh:
        for ch in chems:
            if ch["chemical_id"] in done: continue
            try: syns = json.loads(ch["synonyms"] or "[]")
            except Exception: syns = []
            mid = molecule_id(ch["name"], syns)
            rec = {"chemical_id": ch["chemical_id"], "name": ch["name"], "chembl_id": mid, "hits": []}
            if mid:
                acts, off = [], 0
                while off < 3000:
                    d = get("activity", molecule_chembl_id=mid, limit=500, offset=off, target_organism="Homo sapiens")
                    acts += d.get("activities", [])
                    if not d.get("page_meta", {}).get("next"): break
                    off += 500
                tids = sorted({a["target_chembl_id"] for a in acts if a.get("target_chembl_id")} - set(tcache))
                for i in range(0, len(tids), 50):
                    d = get("target", target_chembl_id__in=",".join(tids[i:i + 50]), limit=50)
                    for t in d.get("targets", []):
                        tcache[t["target_chembl_id"]] = {x["accession"] for comp in t.get("target_components", []) for x in [comp] if x.get("accession")}
                docs = {}
                dids = sorted({a["document_chembl_id"] for a in acts if a.get("document_chembl_id")})
                for i in range(0, len(dids), 50):
                    d = get("document", document_chembl_id__in=",".join(dids[i:i + 50]), limit=50)
                    for x in d.get("documents", []): docs[x["document_chembl_id"]] = x.get("pubmed_id")
                seen = set()
                for a in acts:
                    syms = set()
                    for u in tcache.get(a.get("target_chembl_id"), ()): syms |= up2sym.get(u, set())
                    pm = docs.get(a.get("document_chembl_id"))
                    for s in syms:
                        k = (s, pm, a.get("standard_type"))
                        if k in seen or not pm: continue
                        seen.add(k)
                        rec["hits"].append({"receptor": s, "pmid": str(pm), "type": a.get("standard_type"), "value": a.get("standard_value"),
                                            "units": a.get("standard_units"), "relation": a.get("standard_relation"),
                                            "assay_type": a.get("assay_type"), "assay": (a.get("assay_description") or "")[:160]})
            fh.write(json.dumps(rec) + "\n"); fh.flush()
            print(ch["chemical_id"], ch["name"], mid, len(rec["hits"]), flush=True)
main()
