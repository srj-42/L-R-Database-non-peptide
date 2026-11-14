#!/usr/bin/env python3
"""np17_databases.py - Method 3 of CITATION_METHODS.md: structured protein databases.

Stages (each cached under work/, resumable):
  --chebi     map chemical_entities -> ChEBI ids via PubChem synonyms (work/np17_chebi.json)
  --uniprot   fetch UniProt entries for receptors of unverified pairs (work/np17_uniprot/)
  --match     find annotations naming the pair's molecule -> work/np17_uniprot_hits.tsv
              * Binding site / Cofactor / Catalytic activity (Rhea): matched by ChEBI id
              * FUNCTION / ACTIVITY REGULATION text: matched by molecule name (whole word)
              only annotations backed by an experimental PubMed evidence (ECO:0000269) or a PDB
              structure (ECO:0007744) are kept
  --load      load hits into np_class_evidence as database_annotated (one row per source+pmid)
"""
import hashlib, argparse, csv, json, re, sqlite3, time, urllib.parse, urllib.request
from pathlib import Path
import np10_litsearch as L
from cur_common import connect, now

W = L.WORK; CHEBI = W / "np17_chebi.json"; UPD = W / "np17_uniprot"; HITS = W / "np17_uniprot_hits.tsv"

def get(url, tries=4):
    for a in range(tries):
        try:
            return json.load(urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "np17"}), timeout=60))
        except Exception as e:
            if getattr(e, "code", None) == 404: return None
            time.sleep(2 * (a + 1))
    return None

def unverified(c):
    return c.execute("""
        SELECT ci.chem_interaction_id cid, ch.chemical_id, ch.name, ch.synonyms, ch.pubchem_cid, e.symbol receptor, e.uniprot
          FROM chemical_interactions ci
          JOIN chemical_entities ch ON ch.chemical_id=ci.chemical_id
          JOIN entities e ON e.entity_id=ci.receptor_entity_id
         WHERE ci.species='human' AND ci.shippable=1
           AND NOT EXISTS (SELECT 1 FROM np_evidence n WHERE n.chem_interaction_id=ci.chem_interaction_id
                            AND n.verification_outcome IN ('confirmed','pair_absent'))""").fetchall()

def stage_chebi():
    c = connect(); m = json.loads(CHEBI.read_text()) if CHEBI.exists() else {}
    for r in c.execute("SELECT chemical_id, name, pubchem_cid FROM chemical_entities"):
        if r["chemical_id"] in m: continue
        cid = (r["pubchem_cid"] or "").split(".")[0].strip()
        if not cid:
            d = get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + urllib.parse.quote(r["name"]) + "/cids/JSON")
            cid = str((d or {}).get("IdentifierList", {}).get("CID", [""])[0])
        ids = []
        if cid:
            d = get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON")
            syn = ((d or {}).get("InformationList", {}).get("Information") or [{}])[0].get("Synonym", [])
            ids = sorted({s for s in syn if re.fullmatch(r"CHEBI:\d+", s)})
        m[r["chemical_id"]] = {"pubchem": cid, "chebi": ids}
        time.sleep(0.25)
        if len(m) % 50 == 0: CHEBI.write_text(json.dumps(m))
    CHEBI.write_text(json.dumps(m))
    print(f"{len(m)} chemicals; {sum(1 for v in m.values() if v['chebi'])} with ChEBI")

def stage_uniprot():
    c = connect(); UPD.mkdir(exist_ok=True)
    accs = sorted({u.strip() for r in unverified(c) for u in str(r["uniprot"] or "").replace(";", ",").split(",") if u.strip()})
    n = 0
    for a in accs:
        p = UPD / f"{a}.json"
        if p.exists(): continue
        d = get(f"https://rest.uniprot.org/uniprotkb/{a}.json")
        p.write_text(json.dumps(d or {})); n += 1
    print(f"{len(accs)} accessions; fetched {n} new")

def evid(evs):
    pm = [e["id"] for e in evs or [] if e.get("source") == "PubMed" and e.get("evidenceCode") == "ECO:0000269"]
    pdb = [e["id"] for e in evs or [] if e.get("source") == "PDB"]
    return pm, pdb

def stage_match():
    c = connect(); chebi = json.loads(CHEBI.read_text())
    out = []
    for r in unverified(c):
        ids = set(chebi.get(r["chemical_id"], {}).get("chebi", []))
        names = [n for n in L.chem_synonyms(r) if len(n) >= 4]
        nrx = re.compile(r"(?<![A-Za-z0-9-])(" + "|".join(re.escape(n) for n in names) + r")(?![A-Za-z0-9]|-[A-Za-z0-9])", re.I) if names else None
        for a in [u.strip() for u in str(r["uniprot"] or "").replace(";", ",").split(",") if u.strip()]:
            p = UPD / f"{a}.json"
            if not p.exists(): continue
            d = json.loads(p.read_text() or "{}")
            hits = []
            for f in d.get("features", []):
                if f.get("type") != "Binding site": continue
                lig = f.get("ligand", {}); lid = lig.get("id", "")
                if lid in ids:
                    hits.append(("Binding site", f"binding site for {lig.get('name')} ({lid})", *evid(f.get("evidences"))))
            for cm in d.get("comments", []):
                t = cm.get("commentType")
                if t == "COFACTOR":
                    for cf in cm.get("cofactors", []):
                        x = (cf.get("cofactorCrossReference") or {}).get("id", "")
                        if x in ids: hits.append(("Cofactor", f"cofactor {cf.get('name')} ({x})", *evid(cf.get("evidences"))))
                elif t == "CATALYTIC ACTIVITY":
                    rx = cm.get("reaction", {})
                    xs = {x.get("id") for x in rx.get("reactionCrossReferences", [])}
                    if ids & xs:
                        rhea = next((x for x in xs if x and x.startswith("RHEA:")), "")
                        hits.append(("Catalytic activity (Rhea)", f"{rx.get('name')} [{rhea}]", *evid(rx.get("evidences"))))
                elif t in ("FUNCTION", "ACTIVITY REGULATION") and nrx:
                    for tx in cm.get("texts", []):
                        for sent in re.split(r"(?<=[.;])\s+", tx.get("value", "")):
                            if nrx.search(sent):
                                pm, pdb = evid(tx.get("evidences"))
                                hits.append((t.title(), sent.strip(), pm, pdb))
            for kind, text, pm, pdb in hits:
                if not pm and not pdb: continue      # only experimental / structural backing
                out.append([r["cid"], r["name"], r["receptor"], a, kind, text, ";".join(pm), ";".join(pdb)])
    with open(HITS, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chem_interaction_id", "chemical", "receptor", "uniprot", "annotation", "text", "pmids", "pdb_ids"])
        w.writerows(out)
    print(f"{len(out)} annotations over {len({x[0] for x in out})} pairs -> {HITS.name}")

def stage_load():
    c = connect(); n = 0
    for r in csv.DictReader(open(HITS), delimiter="\t"):
        if r.get("decision", "accept") != "accept": continue
        pms = [p for p in r["pmids"].split(";") if p] or [None]
        for pm in pms:
            src = "UniProt" + (" (Rhea)" if "Rhea" in r["annotation"] else "")
            det = f"{r['annotation']}: {r['text']}" + (f" | PDB {r['pdb_ids']}" if r["pdb_ids"] else "")
            key = f"DB_{r['chem_interaction_id']}_{r['uniprot']}_{pm or 'PDB'}_{hashlib.md5((r["annotation"]+r["text"]).encode()).hexdigest()[:8]}"
            cur = c.execute("INSERT OR IGNORE INTO np_class_evidence VALUES (?,?,?,?,?,?,?,?,?)",
                (key, r["chem_interaction_id"], "database_annotated", "curated", src, r["uniprot"], det[:900], pm, now()))
            n += cur.rowcount
    c.commit(); print(f"loaded {n} database_annotated rows")

ap = argparse.ArgumentParser()
for s in ("chebi", "uniprot", "match", "load"): ap.add_argument("--" + s, action="store_true")
a = ap.parse_args()
if a.chebi: stage_chebi()
if a.uniprot: stage_uniprot()
if a.match: stage_match()
if a.load: stage_load()
