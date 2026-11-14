#!/usr/bin/env python3
"""np18_load.py - load browser-read Results/Discussion records (Method 1b) into np_evidence.

Input: work/np18/out_NN.json  ({pmid: {status, url, records: {cid: [np10 records + SECTION]}}}).
Convention mirrors the peptide half's invariant 13 (dated exception 2026-10-12): evidence_depth
'full_text_browser', extraction_method 'browser_read', notes carry 'source: <url>'. The quote is
NOT verifiable against stored text (the page text is not stored), so these rows are flagged
'spot-check pending'. Mechanical gates that still apply: row exists as abstract_insufficient for
that (pair, pmid); paper not a review/retracted; sentence long enough, contains a recorded name
of the molecule AND the receptor; closed vocabularies; section is Results or Discussion.
  python3 np18_load.py FILE [--commit]
"""
import json, re, sys, unicodedata
from np11_cited_audit import chem_names, rec_names, norm, JUNK_SYN as JUNK
import np10_litsearch as L
from cur_common import connect, now
GREEK = (("\u03b1", "alpha"), ("\u03b2", "beta"), ("\u03b3", "gamma"), ("\u03b4", "delta"), ("\u03ba", "kappa"), ("\u03c9", "omega"))
def fold(x):
    x = x or ""
    for g, w in GREEK: x = x.replace(g, w)
    return norm(x)
GROUPS = json.load(open(L.WORK / "np13_groups.json"))
FAMRX = {cid: re.compile(g["paper_rx"], re.I) for g in GROUPS.values() for cid in g["pairs"]}
RES = {"SUPPORTS", "REFUTES", "FAMILY_LEVEL", "TRANSPORTER", "ENZYME"}
ET = {"structural", "direct_binding", "functional_in_vitro", "functional_in_vivo", "computational"}
def main():
    path, commit = sys.argv[1], "--commit" in sys.argv
    c = connect(); hg = json.loads(L.HGNC_CACHE.read_text()); data = json.load(open(path))
    n = {"ok": 0}; why = []
    for pmid, d in data.items():
        for cid, recs in (d.get("records") or {}).items():
            for r in recs:
                er = c.execute("SELECT np_evidence_id FROM np_evidence WHERE chem_interaction_id=? AND pmid=?", (cid, pmid)).fetchone()
                eid = er[0] if er else f"NPEV_{cid}_{pmid}"
                row = c.execute("""SELECT ne.verification_outcome o, ch.name, ch.synonyms, e.symbol rec, p.is_review, p.is_retracted
                    FROM np_evidence ne JOIN chemical_interactions ci ON ci.chem_interaction_id=ne.chem_interaction_id
                    JOIN chemical_entities ch ON ch.chemical_id=ci.chemical_id JOIN entities e ON e.entity_id=ci.receptor_entity_id
                    LEFT JOIN papers p ON p.pmid=ne.pmid WHERE ne.np_evidence_id=?""", (eid,)).fetchone()
                def bad(m): n[m] = n.get(m, 0) + 1; why.append((eid, m))
                if not row: bad("no_row"); continue
                if row["o"] != "abstract_insufficient": bad("not_abstract_insufficient"); continue
                if row["is_review"] or row["is_retracted"]: bad("review_or_retracted"); continue
                s = (r.get("VERBATIM_SENTENCE") or "").strip()
                if len(s) < 40: bad("sentence_short"); continue
                if r.get("SECTION") not in ("Results", "Discussion"): bad("section"); continue
                if r.get("RESULT") not in RES or r.get("EXPERIMENT_TYPE") not in ET: bad("vocab"); continue
                ns = fold(s)
                allsyn = [row["name"]] + [x for x in json.loads(row["synonyms"] or "[]") if isinstance(x, str) and len(norm(x)) >= 3 and norm(x) not in JUNK]
                if not any(fold(x) in ns for x in chem_names(row["name"], row["synonyms"]) + L.chem_synonyms(row) + allsyn): bad("molecule_not_in_sentence"); continue
                if not (any(fold(str(x)) in ns for x in rec_names(row["rec"], hg)) or (cid in FAMRX and FAMRX[cid].search(unicodedata.normalize('NFKD', s)))): bad("receptor_not_in_sentence"); continue
                notes = (f"source: {d.get('url')}; section: {r['SECTION']}; spot-check pending (browser-read, quote not "
                         f"verifiable against stored text). Prompt answer: {json.dumps({k: v for k, v in r.items() if k != 'SECTION'}, ensure_ascii=False)}")
                if commit:
                    c.execute("""UPDATE np_evidence SET verification_outcome=?, evidence_depth='full_text_browser', evidence_sentence=?,
                        assay_class=?, effect_direction=?, extraction_method='browser_read', extraction_confidence='medium',
                        notes=?, verified_at=? WHERE np_evidence_id=?""",
                        ("pair_absent" if r["RESULT"] == "REFUTES" else "confirmed", s, r["EXPERIMENT_TYPE"],
                         (r.get("DIRECTION") or "").lower() or None, notes, now(), eid))
                n["ok"] += 1
    if commit: c.commit()
    print(n, "COMMITTED" if commit else "DRY RUN")
    for w in why[:15]: print("  ", w)
main()
