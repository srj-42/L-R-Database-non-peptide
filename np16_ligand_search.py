#!/usr/bin/env python3
"""np16_ligand_search.py - Method 2 of CITATION_METHODS.md: ligand-only literature search.

For each ligand with unverified pairs: PubMed query  (ligand names)[tiab] AND (assay terms)[tiab] NOT review[pt],
top 200 by relevance. Abstracts are stored in `papers`. In each abstract, sentences naming the ligand AND any
receptor of that ligand's unverified pairs become candidates -> work/np16/chunk_NN.json (reviewer format).
(pair, pmid) combinations that already have an np_evidence row are skipped (already adjudicated).
Resumable via work/np16/done_ligands.txt.   Usage: np16_ligand_search.py [N_LIGANDS]
"""
import json, re, sys, time
import np10_litsearch as L
from np11_cited_audit import chem_names, rec_names, norm
from cur_common import connect, now
OUT = L.WORK / "np16"; DONE = OUT / "done_ligands.txt"
ASSAY = '(receptor[tiab] OR binding[tiab] OR agonist[tiab] OR antagonist[tiab] OR EC50[tiab] OR IC50[tiab] OR Ki[tiab] OR Kd[tiab])'
def sents(t): return [s.strip() for s in re.split(r"(?<=[.!?])\s+(?=[A-Z(\[])", t or "") if s.strip()]
def main():
    n_lig = int(sys.argv[1]) if len(sys.argv) > 1 else 10**6
    c = connect(); hg = json.loads(L.HGNC_CACHE.read_text())
    done = set(DONE.read_text().split()) if DONE.exists() else set()
    rows = c.execute("""
        SELECT ci.chem_interaction_id cid, ch.chemical_id, ch.name, ch.synonyms, e.symbol receptor, e.entity_id
          FROM chemical_interactions ci JOIN chemical_entities ch ON ch.chemical_id=ci.chemical_id
          JOIN entities e ON e.entity_id=ci.receptor_entity_id
         WHERE ci.species='human' AND ci.shippable=1
           AND NOT EXISTS (SELECT 1 FROM np_evidence n WHERE n.chem_interaction_id=ci.chem_interaction_id AND n.verification_outcome IN ('confirmed','pair_absent'))
           AND NOT EXISTS (SELECT 1 FROM np_class_evidence k WHERE k.chem_interaction_id=ci.chem_interaction_id AND k.evidence_class='nomenclature_evident')""").fetchall()
    by = {}
    for r in rows: by.setdefault(r["chemical_id"], []).append(r)
    existing = {(r[0], r[1]) for r in c.execute("SELECT chem_interaction_id, pmid FROM np_evidence")}
    order = sorted(by, key=lambda k: len(by[k]))              # small ligands first: cheap, focused
    chunk, ci_n, ndone = [], len(list(OUT.glob("chunk_*.json"))), 0
    for chem in order:
        if chem in done: continue
        if ndone >= n_lig: break
        prs = by[chem]; r0 = prs[0]
        syns = L.chem_synonyms(r0)
        term = f"{L.tiab(syns[:L.N_SYN])} AND {ASSAY} NOT review[pt]"
        try:
            total, ids = None, []
            d = json.loads(L.http(L.EUTILS + "esearch.fcgi", {"db": "pubmed", "term": term, "retmode": "json", "retmax": 200, "sort": "relevance"}))
            ids = d["esearchresult"]["idlist"]; time.sleep(L.PAUSE)
            recs = []
            for i in range(0, len(ids), 100): recs += L.efetch(ids[i:i + 100])
            L.store_papers(c, recs); c.commit()
        except Exception as e:
            print("skip", chem, e); continue
        cn = chem_names(r0["name"], r0["synonyms"]) + syns
        rec_by_pair = {p["cid"]: [str(x) for x in rec_names(p["receptor"], hg)] for p in prs}
        for p in recs:
            if p["is_review"] or p["is_retracted"] or not p["abstract"]: continue
            ss = sents(p["abstract"])
            for cid, rn in rec_by_pair.items():
                if (cid, p["pmid"]) in existing: continue
                hits = [{"i": i, "text": s} for i, s in enumerate(ss)
                        if any(norm(x) in norm(s) for x in cn) and any(norm(x) in norm(s) for x in rn)]
                if hits:
                    chunk.append({"cid": cid, "molecule": r0["name"], "molecule_names": cn[:8],
                                  "receptor": next(q["receptor"] for q in prs if q["cid"] == cid), "receptor_names": rn[:8],
                                  "candidates": [{"pmid": p["pmid"], "year": p["year"], "title": p["title"], "sentences": hits}]})
        open(OUT / "search_log.jsonl", "a").write(json.dumps({"chemical_id": chem, "name": r0["name"], "n_pairs": len(prs), "query": term,
                                                             "n_papers": len(ids), "searched_at": now()}) + "\n")
        done.add(chem); ndone += 1; DONE.write_text("\n".join(sorted(done)))
        if len(chunk) >= 40:
            json.dump(chunk, open(OUT / f"chunk_{ci_n:03d}.json", "w"), indent=0); ci_n += 1; chunk = []
            print(f"chunk {ci_n} written; {len(done)} ligands done", flush=True)
    if chunk: json.dump(chunk, open(OUT / f"chunk_{ci_n:03d}.json", "w"), indent=0)
    print(f"{len(done)} ligands done; chunks: {len(list(OUT.glob('chunk_*.json')))}")
main()
