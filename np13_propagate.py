#!/usr/bin/env python3
"""np13_propagate.py - apply the complex rule retroactively (owner instruction 2026-10-12).

For a (chemical, complex) group where one component row is already confirmed by a paper whose EVIDENCE
SENTENCE names the complex generically, the sibling component rows get the same citation: the paper names
the complex, which is enough. SUPPORTS / TRANSPORTER / ENZYME only - a REFUTES is never propagated.
Usage: python3 np13_propagate.py [--commit]
"""
import json, re, sqlite3, sys, datetime
DB = "/home/sreeja/omics/database/curation/curation.db"
G = json.load(open("/home/sreeja/omics/database/curation/work/np13_groups.json"))
commit = "--commit" in sys.argv
c = sqlite3.connect(DB); c.row_factory = sqlite3.Row
have = {(r[0], str(r[1])) for r in c.execute("select chem_interaction_id,pmid from np_evidence")}
new, seen = [], set()
for k, g in G.items():
    rx = re.compile(g["paper_rx"], re.I)
    sib = [p for p in g["pairs"] if p not in g["evidenced"]]
    if not g["evidenced"] or not sib: continue
    for src in g["evidenced"]:
        for r in c.execute("select * from np_evidence where chem_interaction_id=? and verification_outcome='confirmed'", (src,)):
            if not rx.search(r["evidence_sentence"] or ""): continue
            m = re.search(r"Prompt answer: (\{.*\})", r["notes"] or "")
            ans = json.loads(m.group(1)) if m else {}
            res = ans.get("RESULT") or (re.search(r"RESULT=(\w+)", r["notes"] or "") or [None, "SUPPORTS"])[1]
            if res not in ("SUPPORTS", "TRANSPORTER", "ENZYME"): continue
            if (k, r["pmid"]) in seen: continue
            seen.add((k, r["pmid"]))
            for s in sib:
                if (s, str(r["pmid"])) in have: continue
                a2 = dict(ans); a2["AFFINITY"] = (f"propagated by the complex rule from {src}: the paper names the {g['label']}; "
                                                  + str(ans.get("AFFINITY", "")))[:400]
                new.append((f"NPEV_{s}_{r['pmid']}", s, r["pmid"], r["evidence_depth"], r["evidence_sentence"],
                            r["assay_class"], r["effect_direction"],
                            "np13 complex rule (owner 2026-10-12): the paper names the " + g["label"] + ", which covers every "
                            "component row; propagated from " + src + ". Prompt answer: " + json.dumps(a2, ensure_ascii=False)))
                have.add((s, str(r["pmid"])))
print(len(new), "sibling rows to add from", len(seen), "(group, paper) pairs")
if commit:
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    for e in new:
        c.execute("""insert into np_evidence (np_evidence_id,chem_interaction_id,pmid,source_id,verification_outcome,evidence_depth,
            evidence_sentence,assay_class,effect_direction,extraction_method,notes,verified_at)
            values (?,?,?,?,?,?,?,?,?,?,?,?)""", (e[0], e[1], e[2], "np13_complex", "confirmed", e[3], e[4], e[5], e[6],
                                                   "llm_complex_propagation", e[7], now))
    c.commit(); print("COMMITTED")
