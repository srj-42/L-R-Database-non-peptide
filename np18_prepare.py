#!/usr/bin/env python3
"""np18_prepare.py - Method 1b of CITATION_METHODS.md: build the full-text work list.
Papers behind unverified pairs (chembl + upstream abstract_insufficient rows), grouped by paper,
ranked by number of pairs. Writes work/np18/tasks_NN.json (browser) and work/np18/oa_pmids.txt."""
import json, re
import np10_litsearch as L
from np11_cited_audit import chem_names, rec_names
from cur_common import connect
PROFILING = re.compile(r"kinome|selectivity|profiling|kinase inhibitor|panel|library|screen|proteom|atlas|database", re.I)
def main():
    c = connect(); hg = json.loads(L.HGNC_CACHE.read_text()); oa = json.load(open(L.WORK / "np18_oa_check.json"))
    rows = c.execute("""
        SELECT ne.pmid, p.doi, p.title, p.year, p.journal, p.full_text IS NOT NULL AND length(p.full_text)>=20000 hasft,
               ci.chem_interaction_id cid, ch.name chem, ch.synonyms, e.symbol rec, ne.source_id
          FROM np_evidence ne JOIN papers p ON p.pmid=ne.pmid
          JOIN chemical_interactions ci ON ci.chem_interaction_id=ne.chem_interaction_id
          JOIN chemical_entities ch ON ch.chemical_id=ci.chemical_id
          JOIN entities e ON e.entity_id=ci.receptor_entity_id
         WHERE ne.verification_outcome='abstract_insufficient' AND ci.species='human' AND ci.shippable=1
           AND ne.source_id NOT IN ('np10_litsearch','np13_complex')
           AND NOT EXISTS (SELECT 1 FROM np_evidence n2 WHERE n2.chem_interaction_id=ne.chem_interaction_id
                            AND n2.verification_outcome IN ('confirmed','pair_absent'))""").fetchall()
    by = {}
    for r in rows:
        if r["hasft"]: continue                      # already audited at full text (np12)
        t = by.setdefault(r["pmid"], {"pmid": r["pmid"], "doi": r["doi"], "title": r["title"], "year": r["year"],
                                      "journal": r["journal"], "pairs": []})
        t["pairs"].append({"cid": r["cid"], "molecule": r["chem"],
                           "molecule_names": chem_names(r["chem"], r["synonyms"])[:8],
                           "receptor": r["rec"], "receptor_names": [str(x) for x in rec_names(r["rec"], hg)[:8]]})
    open_ids = [p for p, v in oa.items() if v.get("oa") == "Y" and v.get("pmcid")]
    open(L.WORK / "np18" / "oa_pmids.txt", "w").write("\n".join(f"{p}\t{oa[p]['pmcid']}" for p in open_ids if p in by))
    todo = [t for p, t in by.items() if p not in open_ids and t["doi"]]
    skipped = [t for t in todo if len(t["pairs"]) > 15 and PROFILING.search(t["title"] or "")]
    todo = [t for t in todo if t not in skipped]
    # publishers that put a bot check in front of the page (ACS, ASH/Blood): not attempted, logged as blocked
    blocked = [t for t in todo if t["doi"].split("/")[0] in ("10.1021", "10.1182")]
    todo = [t for t in todo if t not in blocked]
    json.dump([{"pmid": t["pmid"], "doi": t["doi"], "journal": t["journal"], "n_pairs": len(t["pairs"])} for t in blocked],
              open(L.WORK / "np18" / "blocked_publisher.json", "w"), indent=0)
    todo.sort(key=lambda t: -len(t["pairs"]))
    json.dump([{"pmid": t["pmid"], "title": t["title"], "n_pairs": len(t["pairs"])} for t in skipped],
              open(L.WORK / "np18" / "skipped_profiling.json", "w"), indent=0)
    for i in range(0, len(todo), 25):
        json.dump(todo[i:i + 25], open(L.WORK / "np18" / f"tasks_{i // 25:02d}.json", "w"), indent=0)
    print(f"{len(blocked)} blocked-publisher papers not attempted; {len(by)} papers; {len(open_ids)} open access; {len(skipped)} skipped as profiling datasets; "
          f"{len(todo)} for the browser in {(len(todo)+24)//25} task files; top pair counts {[len(t['pairs']) for t in todo[:8]]}")
main()
