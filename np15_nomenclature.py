#!/usr/bin/env python3
"""np15_nomenclature.py - Method 4 of CITATION_METHODS.md: name-evident pairs (exact name only).

A pair is `nomenclature_evident` when the receptor's official HGNC full name contains the
molecule's own name (or a recorded synonym naming the same molecule) as a whole word/phrase.
No class map (owner decision 2026-10-15). Scope: human shippable pairs with no confirmed line.

  --build    write work/np15_nomenclature.tsv (every match, auto accept/reject + reason)
  --load     load accepted rows of that TSV into np_class_evidence (after review)
"""
import argparse, csv, json, re, sqlite3, sys, unicodedata
import np10_litsearch as L
from cur_common import connect, now

TSV = L.WORK / "np15_nomenclature.tsv"
# name matches that are not ligand-receptor statements (CITATION_METHODS.md, Method 4)
REJECT = [
    (r"\bzinc finger\b", "zinc finger protein: structural zinc, not a ligand-receptor name"),
    (r"\bATP binding cassette\b", "ABC transporter family name: ATP is the hydrolysed cofactor"),
    (r"\bvoltage-gated channel\b", "ion named as the permeant species of a voltage-gated channel"),
    (r"\b(calcium|sodium|potassium|chloride) (channel|transporting|release-activated)\b", "ion named as permeant species"),
    (r"\barginine vasopressin\b", "'arginine' is part of the peptide name arginine vasopressin"),
    (r"\bbinding protein\b", "'binding protein' name: check partner is this molecule"),
]

def fold(s):
    s = unicodedata.normalize("NFKD", s or "").lower()
    for g, w in (("α", "alpha"), ("β", "beta"), ("γ", "gamma"), ("δ", "delta")):
        s = s.replace(g, w)
    return re.sub(r"[\s_]+", " ", s)

def build():
    c = connect()
    cache = json.loads(L.HGNC_CACHE.read_text()) if L.HGNC_CACHE.exists() else {}
    rows = c.execute("""
        SELECT ci.chem_interaction_id cid, ch.chemical_id, ch.name, ch.synonyms, e.symbol receptor
          FROM chemical_interactions ci
          JOIN chemical_entities ch ON ch.chemical_id = ci.chemical_id
          JOIN entities e ON e.entity_id = ci.receptor_entity_id
         WHERE ci.species='human' AND ci.shippable=1
           AND NOT EXISTS (SELECT 1 FROM np_evidence n WHERE n.chem_interaction_id=ci.chem_interaction_id
                            AND n.verification_outcome IN ('confirmed','pair_absent'))""").fetchall()
    out = []
    for r in rows:
        h = L.hgnc(r["receptor"], cache)
        full = h.get("name") or ""
        if not full:
            continue
        f = fold(full)
        hit = None
        for n in L.chem_synonyms(r):
            k = fold(n).strip()
            if len(k) < 3:
                continue
            if re.search(r"(?<![a-z0-9-])" + re.escape(k) + r"(?![a-z0-9]|-[a-z0-9])", f):
                hit = n; break
        if not hit:
            continue
        dec, why = "accept", "receptor full name contains the molecule name"
        for pat, reason in REJECT:
            if re.search(pat, full, re.I):
                dec, why = "reject", reason; break
        out.append([r["cid"], r["name"], r["receptor"], full, hit, dec, why])
    L.HGNC_CACHE.write_text(json.dumps(cache, indent=0, sort_keys=True))
    with open(TSV, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chem_interaction_id", "chemical", "receptor", "hgnc_full_name", "matched_name", "decision", "reason"])
        w.writerows(sorted(out, key=lambda x: (x[5], x[1], x[2])))
    print(f"{len(rows)} unverified pairs; {len(out)} name matches; "
          f"{sum(1 for x in out if x[5]=='accept')} auto-accept, {sum(1 for x in out if x[5]=='reject')} auto-reject -> {TSV.name}")

def load():
    c = connect()
    c.execute("""CREATE TABLE IF NOT EXISTS np_class_evidence (
        class_evidence_id   TEXT PRIMARY KEY,
        chem_interaction_id TEXT NOT NULL REFERENCES chemical_interactions,
        evidence_class      TEXT NOT NULL,      -- nomenclature_evident | database_annotated
        confidence          TEXT NOT NULL,
        source              TEXT NOT NULL,      -- 'HGNC approved name', 'UniProt', ...
        accession           TEXT,               -- HGNC symbol / UniProt acc / PDB id
        detail              TEXT,               -- matched name / annotation text
        pmid                TEXT,
        created_at          TEXT NOT NULL,
        UNIQUE (chem_interaction_id, evidence_class, source, accession, pmid))""")
    n = 0
    for r in csv.DictReader(open(TSV), delimiter="\t"):
        if r["decision"] != "accept":
            continue
        cur = c.execute("""INSERT OR IGNORE INTO np_class_evidence VALUES (?,?,?,?,?,?,?,?,?)""",
            (f"NOMEN_{r['chem_interaction_id']}", r["chem_interaction_id"], "nomenclature_evident", "high",
             "HGNC approved name", r["receptor"],
             f"{r['receptor']} = \"{r['hgnc_full_name']}\"; contains molecule name \"{r['matched_name']}\"",
             None, now()))
        n += cur.rowcount
    c.commit()
    print(f"loaded {n} nomenclature_evident rows")

ap = argparse.ArgumentParser(); ap.add_argument("--build", action="store_true"); ap.add_argument("--load", action="store_true")
a = ap.parse_args()
build() if a.build else load() if a.load else ap.print_help()
