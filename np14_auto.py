#!/usr/bin/env python3
"""np14_auto.py - automated first pass of the np10 literature search.
For each unsearched pair: same PubMed query as np10 --emit, fetch abstracts, keep only
non-review papers, and split abstracts into sentences. A pair with no paper that has a sentence
naming BOTH the molecule and the receptor is logged "No evidence" (automatic, reason recorded);
the rest go to work/np14/auto_NNN.json for reviewers. Complex rule: pairs in the same np13 group
(same chemical, same receptor family) are searched once; siblings take the same outcome.
Resumable: every pair is logged to np10_search_log.jsonl as soon as it is searched.
Usage: np14_auto.py N_PAIRS_PER_RUN [PAIRS_PER_CHUNK]"""
import json, re, sys, time
import np10_litsearch as L
from np11_cited_audit import chem_names, rec_names, norm
from cur_common import connect, now
def sents(t): return [s.strip() for s in re.split(r"(?<=[.!?])\s+(?=[A-Z(\[])", t or "") if s.strip()]
def main():
    n, per = int(sys.argv[1]), int(sys.argv[2]) if len(sys.argv) > 2 else 60
    c = connect(); cache = json.loads(L.HGNC_CACHE.read_text()) if L.HGNC_CACHE.exists() else {}
    groups = json.load(open("work/np13_groups.json")); gof = {}
    for k, g in groups.items() if isinstance(groups, dict) else enumerate(groups):
        for cid in g["pairs"]: gof[cid] = k
    rows = L.targets(c, n * 3)
    seen_g, picked, sib = {}, [], {}
    for r in rows:
        g = gof.get(r["chem_interaction_id"])
        if g is not None and g in seen_g: sib.setdefault(seen_g[g], []).append(r); continue
        if g is not None: seen_g[g] = r["chem_interaction_id"]
        picked.append(r)
        if len(picked) >= n: break
    noev, chunk, log_lines = [], [], []
    stamp = int(time.time())
    for r in picked:
        syns = L.chem_synonyms(r); h = L.hgnc(r["receptor"], cache)
        rq = L.tiab(L.receptor_terms(r["receptor"], h))
        if sc := L.subunit_clause(r["receptor"]): rq = f"({rq[1:-1]} OR {sc})"
        term = f"{L.tiab(syns[:L.N_SYN])} AND {rq} NOT review[pt]"
        try:
            total, ids = L.esearch(term)
            if total > 2 * L.RETMAX:
                t2 = f"{term[:-len(' NOT review[pt]')]} AND {L.ASSAY} NOT review[pt]"
                tt, i2 = L.esearch(t2)
                if i2: term, total, ids = t2, tt, i2
            recs = L.efetch(ids); L.store_papers(c, recs); c.commit()
        except Exception as e:
            print("skip", r["chem_interaction_id"], e); continue
        cn = chem_names(r["name"], r["synonyms"]) + syns; rn = rec_names(r["receptor"], h)
        cands = []
        for p in recs:
            if p["is_review"] or p["is_retracted"] or not p["abstract"]: continue
            ss = sents(p["abstract"])
            both = [{"i": i, "text": s} for i, s in enumerate(ss)
                    if any(norm(x) in norm(s) for x in cn) and any(norm(str(x)) in norm(s) for x in rn)]
            if both: cands.append({"pmid": p["pmid"], "year": p["year"], "title": p["title"], "sentences": both})
        cid = r["chem_interaction_id"]
        log_lines.append(json.dumps({"chem_interaction_id": cid, "chemical": r["name"], "receptor": r["receptor"],
                                     "query": term, "n_total": total, "hits": ids, "searched_at": now(), "by": "np14_auto"}))
        members = [cid] + [s["chem_interaction_id"] for s in sib.get(cid, [])]
        for s in sib.get(cid, []):
            log_lines.append(json.dumps({"chem_interaction_id": s["chem_interaction_id"], "chemical": s["name"],
                  "receptor": s["receptor"], "query": "(complex rule: searched with " + cid + ")", "n_total": total,
                  "hits": ids, "searched_at": now(), "by": "np14_auto"}))
        if cands:
            chunk.append({"cid": cid, "siblings": members[1:], "molecule": r["name"], "molecule_names": cn[:10],
                          "receptor": r["receptor"], "receptor_names": [str(x) for x in rn[:10]], "candidates": cands})
        else:
            reason = ("AUTO first pass: PubMed query returned %d paper(s); none non-review with a sentence naming both the "
                      "molecule and the receptor (np14_auto)" % total)
            for m in members:
                noev.append((m, r["name"], r["receptor"], reason))
        if len(log_lines) >= 10:
            open(L.LOGFILE, "a").write("".join(l + "\n" for l in log_lines)); log_lines = []
    open(L.LOGFILE, "a").write("".join(l + "\n" for l in log_lines))
    with open("work/np10_no_evidence.tsv", "a") as fh:
        for m, ch, rc, why in noev: fh.write(f"{m}\t{ch}\t{rc}\t{why}\t{now()}\n")
    L.HGNC_CACHE.write_text(json.dumps(cache, indent=0, sort_keys=True))
    for i in range(0, len(chunk), per):
        json.dump(chunk[i:i + per], open(f"work/np14/auto_{stamp}_{i // per:02d}.json", "w"), indent=0)
    print(f"searched {len(picked)} groups/pairs; no-evidence {len(noev)}; with candidates {len(chunk)}")
main()
