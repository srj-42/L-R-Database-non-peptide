#!/usr/bin/env python3
"""np14_chunks.py - automated first pass for np11 to_read rows: split each abstract into sentences,
flag the ones naming the molecule and/or the receptor, write chunk files for reviewers.
Output: work/np14/chunk_NN.json (input), reviewers write work/np14/out_NN.json."""
import json, re, sys
from np11_cited_audit import *   # noqa
CH = 115
def sents(t):
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+(?=[A-Z(\[])", t or "") if s.strip()]
def main():
    c = connect(); st = load_state(); hgnc = json.loads(HGNC.read_text())
    todo = {e for e in st["to_read"] if e not in st["decided"]}
    rows = [r for r in pending_rows(c) if r["eid"] in todo]
    notes = dict(c.execute("select np_evidence_id,notes from np_evidence").fetchall())
    out = []
    for r in rows:
        cn = chem_names(r["chem"], r["synonyms"]); rn = rec_names(r["rec"], hgnc)
        ss = sents(r["abstract"]); flag = []
        for s in ss:
            n = norm(s)
            a = any(norm(x) in n for x in cn); b = any(norm(str(x)) in n for x in rn)
            flag.append("BOTH" if a and b else "MOL" if a else "REC" if b else "")
        out.append({"eid": r["eid"], "molecule": r["chem"], "molecule_names": cn[:12],
                    "receptor": r["rec"], "receptor_names": [str(x) for x in rn[:12]],
                    "pmid": r["pmid"], "title": r["title"], "year": r["year"],
                    "chembl": (notes.get(r["eid"]) or "")[:200],
                    "sentences": [{"i": i, "flag": f, "text": s} for i, (s, f) in enumerate(zip(ss, flag))]})
    def tier1(x):
        fl = [q["flag"] for q in x["sentences"]]
        return "BOTH" in fl or any(fl[i] and fl[i + 1] and fl[i] != fl[i + 1] for i in range(len(fl) - 1))
    auto = [{"np_evidence_id": x["eid"], "verification_outcome": "abstract_insufficient", "evidence_depth": "abstract_only",
             "extraction_method": "np14_auto_first_pass",
             "notes": "np14 automated first pass: no sentence (or adjacent sentence pair) in the abstract names both the molecule and the receptor; ChEMBL lists the pair for this paper so the result is probably in tables/full text."}
            for x in out if not tier1(x)]
    json.dump(auto, open("work/np14/auto_insufficient.json", "w"))
    out = [x for x in out if tier1(x)]
    out.sort(key=lambda x: (x["pmid"], x["molecule"]))
    for i in range(0, len(out), CH):
        json.dump(out[i:i + CH], open(f"work/np14/chunk_{i // CH:02d}.json", "w"), indent=0)
    both = sum(any(s["flag"] == "BOTH" for s in x["sentences"]) for x in out)
    print(len(out), "rows;", both, "have a sentence naming both;", (len(out) + CH - 1) // CH, "chunks")
main()
