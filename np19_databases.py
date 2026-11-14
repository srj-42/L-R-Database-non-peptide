#!/usr/bin/env python3
"""np19_databases.py - Method 3 (structured databases) for ALL chemical pairs, verified or not (owner instruction 2026-10-31).

Owner decision: a curated database record that links this molecule to this receptor counts as DIRECT evidence
(evidence_class='database_annotated' in np_class_evidence), and the pair becomes verified. Sources:
  resolve   receptor symbol -> human UniProt accession (the 740 symbols with no entities.uniprot were skipped by np17)
  chem      molecule keys: PubChem InChIKey / SMILES (RDKit for SMILES-only), ChEBI ids from np17_chebi.json
  chembl    ChEMBL activities (binding B / functional F assays, human single-protein targets) with a PubMed id
  gtopdb    GtoPdb interaction records already in chem_source_records (affinity, action, PMIDs)
  bindingdb BindingDB ligands per receptor UniProt, molecule matched on InChIKey connectivity block
  pdb       RCSB structures that contain the molecule AND a chain with the receptor's UniProt accession
  tcdb      TCDB substrate annotations (transporters) matched on ChEBI
  uniprot   UniProt FUNCTION / ACTIVITY REGULATION / binding-site / Rhea annotations with experimental PubMed or PDB backing
  load      insert hits into np_class_evidence (database_annotated), one row per source+accession+pmid
Every stage is cached under work/np19/ and resumable.  Run with ../.venv_np19/bin/python (needs rdkit, requests).
  np19_databases.py --resolve --chem --chembl --gtopdb --bindingdb --pdb --tcdb --uniprot --load
"""
import argparse, csv, hashlib, json, re, sqlite3, sys, time, urllib.parse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import requests
from rdkit import Chem
from rdkit.Chem.inchi import MolToInchiKey
import np10_litsearch as L
from cur_common import connect, now

W = Path(__file__).parent / "work" / "np19"; W.mkdir(parents=True, exist_ok=True)
UPMAP = W / "uniprot_map.json"; CHEM = W / "chem.json"
S = requests.Session(); S.headers["User-Agent"] = "np19-lr-curation"


def get(url, params=None, tries=4, as_json=True, timeout=90):
    for a in range(tries):
        try:
            r = S.get(url, params=params, timeout=timeout)
            if r.status_code == 404: return None
            if r.status_code in (429, 500, 502, 503, 504): raise RuntimeError(r.status_code)
            r.raise_for_status()
            return r.json() if as_json else r.text
        except Exception:
            time.sleep(2 * (a + 1))
    return None


def pairs(c):
    """All shippable human chemical pairs, with receptor symbol and uniprot (None until resolved)."""
    return c.execute("""SELECT ci.chem_interaction_id cid, ci.chemical_id, ch.name, ch.synonyms, ch.pubchem_cid, ch.smiles,
                               e.symbol receptor, e.uniprot, ci.receptor_type rtype
                          FROM chemical_interactions ci JOIN chemical_entities ch ON ch.chemical_id=ci.chemical_id
                          JOIN entities e ON e.entity_id=ci.receptor_entity_id
                         WHERE ci.species='human' AND ci.shippable=1""").fetchall()


def accs_of(r, upmap):
    out = [u.strip() for u in str(r["uniprot"] or "").replace(";", ",").split(",") if u.strip()]
    if not out: out = upmap.get(r["receptor"], [])
    return out


# ---------------------------------------------------------------- resolve
def stage_resolve():
    c = connect(); m = json.loads(UPMAP.read_text()) if UPMAP.exists() else {}
    syms = sorted({r["receptor"] for r in pairs(c) if not (r["uniprot"] or "").strip()} - set(m))
    def one(s):
        d = get("https://rest.uniprot.org/uniprotkb/search", {"query": f"gene_exact:{s} AND organism_id:9606 AND reviewed:true",
                                                               "fields": "accession", "size": 3, "format": "json"})
        return s, [x["primaryAccession"] for x in (d or {}).get("results", [])]
    with ThreadPoolExecutor(4) as ex:
        for k, (s, acc) in enumerate(ex.map(one, syms)):
            m[s] = acc
            if k % 100 == 0: UPMAP.write_text(json.dumps(m))
    UPMAP.write_text(json.dumps(m))
    print(f"resolved {sum(1 for v in m.values() if v)} / {len(m)} symbols")


# ------------------------------------------------------------------- chem
def stage_chem():
    c = connect(); m = json.loads(CHEM.read_text()) if CHEM.exists() else {}
    chebi = json.loads((L.WORK / "np17_chebi.json").read_text())
    rows = c.execute("SELECT chemical_id,name,pubchem_cid,smiles FROM chemical_entities WHERE shippable=1").fetchall()
    cids = {}
    for r in rows:
        cid = (r["pubchem_cid"] or "").split(".")[0].strip()
        if cid and r["chemical_id"] not in m: cids[cid] = r["chemical_id"]
    ids = list(cids)
    for i in range(0, len(ids), 100):
        d = get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + ",".join(ids[i:i + 100]) +
                "/property/InChIKey,IsomericSMILES/JSON")
        for p in (d or {}).get("PropertyTable", {}).get("Properties", []):
            ch = cids.get(str(p["CID"]))
            if ch: m[ch] = {"inchikey": p.get("InChIKey"), "smiles": p.get("SMILES") or p.get("IsomericSMILES")}
    for r in rows:
        e = m.setdefault(r["chemical_id"], {})
        if not e.get("inchikey") and r["smiles"]:
            try:
                mol = Chem.MolFromSmiles(r["smiles"]); e["inchikey"] = MolToInchiKey(mol) if mol else None; e["smiles"] = r["smiles"]
            except Exception: pass
        e["ik14"] = (e.get("inchikey") or "")[:14] or None
        e["chebi"] = (chebi.get(r["chemical_id"]) or {}).get("chebi", [])
    CHEM.write_text(json.dumps(m))
    print(f"{len(m)} chemicals; {sum(1 for v in m.values() if v.get('inchikey'))} with InChIKey; {sum(1 for v in m.values() if v.get('chebi'))} with ChEBI")


# ------------------------------------------------------------------ chembl
CH = W / "chembl"; API = "https://www.ebi.ac.uk/chembl/api/data/"

def stage_chembl():
    CH.mkdir(exist_ok=True); c = connect(); chem = json.loads(CHEM.read_text())
    todo = [r["chemical_id"] for r in c.execute("SELECT chemical_id FROM chemical_entities WHERE shippable=1")]
    tcache = {}
    def target_uniprots(tid):
        if tid not in tcache:
            d = get(API + f"target/{tid}.json") or {}
            tcache[tid] = [cc.get("accession") for cc in (d.get("target_components") or []) if cc.get("accession")], d.get("target_type")
        return tcache[tid]
    def one(ch):
        out = CH / f"{ch}.json"
        if out.exists(): return
        ik = (chem.get(ch) or {}).get("inchikey")
        mols = []
        if ik:
            d = get(API + "molecule.json", {"molecule_structures__standard_inchi_key": ik, "limit": 5}) or {}
            mols = [m["molecule_chembl_id"] for m in d.get("molecules", [])]
        if not mols: out.write_text(json.dumps({"mol": None, "acts": []})); return
        acts = []; off = 0
        while off < 20000:
            d = get(API + "activity.json", {"molecule_chembl_id": mols[0], "target_organism": "Homo sapiens", "assay_type__in": "B,F",
                                            "limit": 1000, "offset": off}) or {}
            for a in d.get("activities", []):
                acts.append({"t": a.get("target_chembl_id"), "pm": a.get("document_chembl_id") and a.get("pubmed_id"),
                             "ty": a.get("standard_type"), "v": a.get("standard_value"), "u": a.get("standard_units"),
                             "rel": a.get("standard_relation"), "at": a.get("assay_type"), "ad": (a.get("assay_description") or "")[:240]})
            if not (d.get("page_meta") or {}).get("next"): break
            off += 1000
        out.write_text(json.dumps({"mol": mols[0], "acts": acts}))
    with ThreadPoolExecutor(3) as ex:
        list(ex.map(one, todo))
    # target -> uniprot map (single protein targets only)
    tids = sorted({a["t"] for p in CH.glob("CHEM*.json") for a in json.loads(p.read_text())["acts"] if a["t"]})
    tm = json.loads((W / "chembl_targets.json").read_text()) if (W / "chembl_targets.json").exists() else {}
    need = [t for t in tids if t not in tm]
    def tgt(t):
        d = get(API + f"target/{t}.json") or {}
        return t, {"type": d.get("target_type"), "acc": [cc.get("accession") for cc in d.get("target_components") or [] if cc.get("accession")]}
    with ThreadPoolExecutor(4) as ex:
        for k, (t, v) in enumerate(ex.map(tgt, need)):
            tm[t] = v
            if k % 200 == 0: (W / "chembl_targets.json").write_text(json.dumps(tm))
    (W / "chembl_targets.json").write_text(json.dumps(tm))
    print(f"chembl: {len(todo)} chemicals, {len(tids)} targets")


# ----------------------------------------------------------------- bindingdb
BD = W / "bindingdb"

def stage_bindingdb():
    BD.mkdir(exist_ok=True); c = connect(); chem = json.loads(CHEM.read_text()); upmap = json.loads(UPMAP.read_text())
    accs = sorted({a for r in pairs(c) for a in accs_of(r, upmap)})
    def one(a):
        out = BD / f"{a}.json"
        if out.exists(): return
        d = get("https://www.bindingdb.org/rest/getLigandsByUniprots", {"uniprot": a, "cutoff": 10000, "response": "application/json"}, timeout=180)
        hits = {}
        for x in ((d or {}).get("getLindsByUniprotsResponse") or {}).get("affinities", []):
            try:
                mol = Chem.MolFromSmiles(x.get("smile") or "")
                ik = MolToInchiKey(mol)[:14] if mol else None
            except Exception: ik = None
            if ik: hits.setdefault(ik, []).append([x.get("affinity_type"), x.get("affinity"), x.get("pmid"), x.get("doi")])
        out.write_text(json.dumps(hits))
    with ThreadPoolExecutor(3) as ex:
        list(ex.map(one, accs))
    print(f"bindingdb: {len(accs)} accessions")


# ------------------------------------------------------------------------ pdb
PDB = W / "pdb"
RCSB = "https://search.rcsb.org/rcsbsearch/v2/query"; GQL = "https://data.rcsb.org/graphql"

def stage_pdb():
    PDB.mkdir(exist_ok=True); chem = json.loads(CHEM.read_text())
    def one(ch):
        out = PDB / f"{ch}.json"; ik = chem[ch].get("inchikey")
        if out.exists() or not ik: return
        q = {"query": {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_chem_comp_descriptor.InChIKey",
                                                                            "operator": "exact_match", "value": ik}},
             "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": 10000}}}
        try: r = S.post(RCSB, json=q, timeout=90); ids = [x["identifier"] for x in r.json().get("result_set", [])] if r.status_code == 200 else []
        except Exception: ids = []
        res = {}
        for i in range(0, len(ids), 100):
            gq = {"query": "query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id rcsb_primary_citation{pdbx_database_id_PubMed} "
                           "polymer_entities{rcsb_polymer_entity_container_identifiers{uniprot_ids}}}}", "variables": {"ids": ids[i:i + 100]}}
            try: g = S.post(GQL, json=gq, timeout=90).json()
            except Exception: continue
            for e in ((g or {}).get("data") or {}).get("entries") or []:
                pm = (e.get("rcsb_primary_citation") or {}).get("pdbx_database_id_PubMed")
                for pe in e.get("polymer_entities") or []:
                    for u in (pe.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids") or []:
                        res.setdefault(u, []).append([e["rcsb_id"], pm])
        out.write_text(json.dumps(res))
    with ThreadPoolExecutor(3) as ex:
        list(ex.map(one, list(chem)))
    print("pdb done")


# ----------------------------------------------------------------------- tcdb
def stage_tcdb():
    sub = get("https://tcdb.org/cgi-bin/substrates/getSubstrates.py", as_json=False) or ""
    acc = get("https://tcdb.org/cgi-bin/projectv/public/acc2tcid.py", as_json=False) or ""
    (W / "tcdb_substrates.tsv").write_text(sub); (W / "tcdb_acc2tcid.tsv").write_text(acc)
    print(f"tcdb substrates {len(sub.splitlines())} lines; acc2tcid {len(acc.splitlines())} lines")


# --------------------------------------------------------------------- uniprot
UPD = L.WORK / "np17_uniprot"

def stage_uniprot_fetch():
    UPD.mkdir(exist_ok=True); c = connect(); upmap = json.loads(UPMAP.read_text())
    accs = sorted({a for r in pairs(c) for a in accs_of(r, upmap)}); n = 0
    def one(a):
        p = UPD / f"{a}.json"
        if p.exists(): return 0
        d = get(f"https://rest.uniprot.org/uniprotkb/{a}.json"); p.write_text(json.dumps(d or {})); return 1
    with ThreadPoolExecutor(4) as ex:
        n = sum(ex.map(one, accs))
    print(f"uniprot: {len(accs)} accessions; fetched {n} new")


def evid(evs):
    pm = [e["id"] for e in evs or [] if e.get("source") == "PubMed" and e.get("evidenceCode") == "ECO:0000269"]
    pdb = [e["id"] for e in evs or [] if e.get("source") == "PDB"]
    return pm, pdb


# ------------------------------------------------------------------ collect/load
PDB_SKIP = re.compile(r"^(sulfate|glycerol|polyethylene glycol|chloride( ion)?|ethylene glycol|phosphate|acetate|formate|citrate|"
                      r"dimethyl sulfoxide|tris|calcium|sodium|potassium|magnesium|zinc|iron|copper|manganese|cobalt|nickel|water|hydron|"
                      r"carbonate|nitrate|iodide|bromide|cadmium|lithium|mercury)", re.I)
BIND_VERB = re.compile(r"\b(bind|binds|binding|bound|activat\w*|agonist\w*|antagonist\w*|ligand\w*|receptor for|inhibit\w*|"
                       r"gated|transport\w*|substrate\w*|stimulat\w*|sensor|sensing|respon\w*|modulat\w*)\b", re.I)
BAD_CTX = re.compile(r"\b(levels?|release|synthesis|metabolism|production|uptake in|concentration)\b", re.I)


def collect():
    """Return list of hit dicts: cid, source, accession, detail, pmid."""
    c = connect(); chem = json.loads(CHEM.read_text()); upmap = json.loads(UPMAP.read_text()); P = pairs(c)
    hits = []
    # --- ChEMBL
    tm = json.loads((W / "chembl_targets.json").read_text()) if (W / "chembl_targets.json").exists() else {}
    for r in P:
        f = CH / f"{r['chemical_id']}.json"
        if not f.exists(): continue
        d = json.loads(f.read_text()); accs = set(accs_of(r, upmap))
        seen = {}
        for a in d["acts"]:
            t = tm.get(a["t"]) or {}
            if t.get("type") != "SINGLE PROTEIN" or not (accs & set(t.get("acc", []))): continue
            if not a["ty"] or a["v"] in (None, ""): continue
            if a.get("rel") in (">", ">="): continue                       # inactive / lower bound only
            try: val = float(a["v"])
            except Exception: continue
            if (a.get("u") or "nM") == "nM" and val > 10000: continue      # weaker than 10 uM: not counted
            if a["ty"] not in ("Ki", "Kd", "IC50", "EC50", "AC50", "Potency"): continue
            key = (a["pm"] or "")
            seen.setdefault(key, []).append(a)
        for pm, lst in seen.items():
            a = lst[0]
            best = next((x for x in lst if x["ty"] in ("Ki", "Kd", "IC50", "EC50")), a)
            det = f"{d['mol']}: {best['ty']} {best.get('rel') or ''} {best['v']} {best['u'] or ''} ({best['at']}); {best['ad']}"[:600]
            hits.append(dict(cid=r["cid"], source="ChEMBL", accession=d["mol"], detail=det, pmid=pm or None))
    # --- GtoPdb
    for r in c.execute("""SELECT chem_interaction_id cid, raw_reference, raw_ligand, raw_receptor, raw_row FROM chem_source_records
                          WHERE source_id='gtopdb_interactions' AND chem_interaction_id IS NOT NULL"""):
        try: row = json.loads(r["raw_row"] or "{}")
        except Exception: row = {}
        det = "GtoPdb interaction: " + "; ".join(f"{k}={row[k]}" for k in ("Type", "Action", "Endogenous", "Affinity Units", "Affinity Low", "Affinity High")
                                                  if row.get(k))[:500]
        pms = [p for p in re.split(r"[|;,]", (r["raw_reference"] or "") or str(row.get("PubMed ID") or "")) if p.strip().isdigit()] or [None]
        for pm in pms: hits.append(dict(cid=r["cid"], source="GtoPdb", accession=r["raw_receptor"] or "", detail=det or "GtoPdb interaction record", pmid=pm))
    # --- BindingDB
    for r in P:
        ik = (chem.get(r["chemical_id"]) or {}).get("ik14")
        if not ik: continue
        for a in accs_of(r, upmap):
            f = BD / f"{a}.json"
            if not f.exists(): continue
            h = json.loads(f.read_text()).get(ik)
            if not h: continue
            by = {}
            for ty, v, pm, doi in h: by.setdefault(pm or "", []).append((ty, v))
            for pm, lst in by.items():
                ty, v = lst[0]
                hits.append(dict(cid=r["cid"], source="BindingDB", accession=a, detail=f"BindingDB: {ty} {v} nM", pmid=pm or None))
    # --- PDB
    for r in P:
        f = PDB / f"{r['chemical_id']}.json"
        if not f.exists(): continue
        d = json.loads(f.read_text())
        for a in accs_of(r, upmap):
            if PDB_SKIP.search(r["name"] or ""): continue
            for entry, pm in d.get(a, [])[:5]:
                hits.append(dict(cid=r["cid"], source="PDB", accession=entry, detail=f"PDB {entry}: structure contains this molecule and a chain of {a}", pmid=pm))
    # --- TCDB
    if (W / "tcdb_substrates.tsv").exists():
        sub = {}
        for line in (W / "tcdb_substrates.tsv").read_text().splitlines():
            p = line.split("\t")
            if len(p) == 2: sub[p[0]] = {x.split(";")[0] for x in p[1].split("|")}
        a2t = {}
        for line in (W / "tcdb_acc2tcid.tsv").read_text().splitlines():
            p = line.split()
            if len(p) >= 2 and re.match(r"^[A-Z0-9]{6,10}$", p[0]): a2t.setdefault(p[0], set()).add(p[1])
        for r in P:
            ids = set((chem.get(r["chemical_id"]) or {}).get("chebi", []))
            if not ids: continue
            for a in accs_of(r, upmap):
                for tc in a2t.get(a, ()):
                    if ids & sub.get(tc, set()):
                        hits.append(dict(cid=r["cid"], source="TCDB", accession=tc, detail=f"TCDB {tc}: substrate annotation (ChEBI match) for {a}", pmid=None))
    # --- UniProt (widened: all pairs, receptors resolved through the map)
    for r in P:
        ids = set((chem.get(r["chemical_id"]) or {}).get("chebi", []))
        names = [n for n in L.chem_synonyms(r) if len(n) >= 4]
        nrx = re.compile(r"(?<![A-Za-z0-9-])(" + "|".join(re.escape(n) for n in names) + r")(?![A-Za-z0-9]|-[A-Za-z0-9])", re.I) if names else None
        for a in accs_of(r, upmap):
            p = UPD / f"{a}.json"
            if not p.exists(): continue
            d = json.loads(p.read_text() or "{}")
            for f in d.get("features", []):
                if f.get("type") == "Binding site" and (f.get("ligand") or {}).get("id") in ids:
                    pm, pdb = evid(f.get("evidences"))
                    for q in (pm or pdb or []):
                        hits.append(dict(cid=r["cid"], source="UniProt", accession=a, detail=f"binding site for {f['ligand'].get('name')}", pmid=q if q.isdigit() else None))
            for cm in d.get("comments", []):
                t = cm.get("commentType")
                if t == "CATALYTIC ACTIVITY":
                    rx = cm.get("reaction", {}); xs = {x.get("id") for x in rx.get("reactionCrossReferences", [])}
                    if ids & xs:
                        pm, pdb = evid(rx.get("evidences"))
                        for q in (pm or pdb or []): hits.append(dict(cid=r["cid"], source="UniProt (Rhea)", accession=a, detail=str(rx.get("name"))[:500], pmid=q if q.isdigit() else None))
                elif t in ("FUNCTION", "ACTIVITY REGULATION") and nrx:
                    for tx in cm.get("texts", []):
                        for sent in re.split(r"(?<=[.;])\s+", tx.get("value", "")):
                            if nrx.search(sent) and BIND_VERB.search(sent) and not BAD_CTX.search(sent) and not PDB_SKIP.search(r["name"] or ""):
                                pm, pdb = evid(tx.get("evidences"))
                                for q in (pm or pdb or []): hits.append(dict(cid=r["cid"], source="UniProt", accession=a, detail=f"{t.title()}: {sent.strip()}"[:600], pmid=q if q.isdigit() else None))
    return hits


def stage_load():
    c = connect(); hits = collect(); n = 0; by = {}
    for h in hits:
        key = "DB19_" + hashlib.md5(f"{h['cid']}|{h['source']}|{h['accession']}|{h['pmid']}".encode()).hexdigest()[:16]
        cur = c.execute("INSERT OR IGNORE INTO np_class_evidence VALUES (?,?,?,?,?,?,?,?,?)",
                        (key, h["cid"], "database_annotated", "curated", h["source"], h["accession"], h["detail"][:900], h["pmid"], now()))
        n += cur.rowcount; by.setdefault(h["source"], set()).add(h["cid"])
    c.commit()
    print(f"loaded {n} new database_annotated rows")
    for k, v in sorted(by.items()): print(f"  {k:16s} {len(v)} pairs")
    print("pairs with any database evidence from this run:", len({h['cid'] for h in hits}))


ap = argparse.ArgumentParser()
for s in ("resolve", "chem", "chembl", "gtopdb", "bindingdb", "pdb", "tcdb", "uniprot", "load", "count"): ap.add_argument("--" + s, action="store_true")
a = ap.parse_args()
if a.resolve: stage_resolve()
if a.chem: stage_chem()
if a.chembl: stage_chembl()
if a.bindingdb: stage_bindingdb()
if a.pdb: stage_pdb()
if a.tcdb: stage_tcdb()
if a.uniprot: stage_uniprot_fetch()
if a.count:
    h = collect(); by = {}
    for x in h: by.setdefault(x["source"], set()).add(x["cid"])
    print({k: len(v) for k, v in by.items()}, "union:", len({x['cid'] for x in h}))
if a.load: stage_load()
