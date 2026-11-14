#!/usr/bin/env python3
"""np21_pathways_ontology.py - HMDB pathways and ontology for the chemical half.

Adds two tables (owner-approved additive exceptions, 2026-11-06):
  chemical_pathway(chemical_id, pathway_name, family, smpdb_id, kegg_map_id, source)
  chemical_ontology(chemical_id, root, path, term, level, node_type, source)
Source is HMDB 5.0 (raw/hmdb_metabolites.zip), cached by np20_hmdb_parse.py in
work/np20/hmdb_context.json. Matched on chemical_entities.hmdb_id.

PATHWAYS. HMDB lists every SMPDB/KEGG pathway the metabolite takes part in (as
substrate, product or cofactor), not only pathways it signals through. SMPDB also
holds one auto-generated entry per lipid species ('Cardiolipin Biosynthesis
CL(a-13:0/i-24:0/...)'), which gives ubiquitous cofactors such as hydrogen ion
~48,000 rows. Every raw row is kept; `family` is the name with that trailing
species token removed (811 families overall) and is what the export and site use.

ONTOLOGY. The HMDB tree is flattened to one row per node with its full path.
Roots: Physiological effect, Disposition, Process, Role. Branches that repeat data
already held (Biological location = chemical_context tissue/biospecimen/cell
location, Health condition = disease, Biochemical pathway = chemical_pathway) are
stored for completeness but not shown on the site.

Both tables describe the MOLECULE, not the ligand-receptor pair, and feed no
verification status.  Usage: python3 np21_pathways_ontology.py [--commit]
"""
import json, re, sqlite3, sys

DB, CACHE, TAG = "curation.db", "work/np20/hmdb_context.json", "hmdb"
commit = "--commit" in sys.argv

def family(name):
    m = re.match(r'^(.*?)\s+[A-Za-z0-9\-]{1,12}\(.*\)\s*$', name)
    return (m.group(1) if m else name).strip()

c = sqlite3.connect(DB)
c.execute("""CREATE TABLE IF NOT EXISTS chemical_pathway(
    chemical_id TEXT NOT NULL, pathway_name TEXT NOT NULL, family TEXT NOT NULL,
    smpdb_id TEXT NOT NULL DEFAULT '', kegg_map_id TEXT NOT NULL DEFAULT '', source TEXT NOT NULL,
    PRIMARY KEY (chemical_id, pathway_name, smpdb_id, kegg_map_id, source))""")
c.execute("""CREATE TABLE IF NOT EXISTS chemical_ontology(
    chemical_id TEXT NOT NULL, root TEXT NOT NULL, path TEXT NOT NULL, term TEXT NOT NULL,
    level INTEGER NOT NULL, node_type TEXT NOT NULL DEFAULT '', source TEXT NOT NULL,
    PRIMARY KEY (chemical_id, path, source))""")
c.execute("DELETE FROM chemical_pathway WHERE source=?", (TAG,))
c.execute("DELETE FROM chemical_ontology WHERE source=?", (TAG,))

d = json.load(open(CACHE))
pw, on = [], []
for cid, hid in c.execute("SELECT chemical_id, hmdb_id FROM chemical_entities "
                          "WHERE hmdb_id IS NOT NULL AND hmdb_id!=''").fetchall():
    rec = d.get(hid)
    if not rec:
        continue
    for name, smpdb, kegg in rec.get("pathways", []):
        pw.append((cid, name, family(name), smpdb, kegg, TAG))
    for root, path, term, level, ntype in rec.get("ontology", []):
        on.append((cid, root, path, term, level, ntype, TAG))
c.executemany("INSERT OR IGNORE INTO chemical_pathway VALUES (?,?,?,?,?,?)", pw)
c.executemany("INSERT OR IGNORE INTO chemical_ontology VALUES (?,?,?,?,?,?,?)", on)

q = lambda s: c.execute(s).fetchone()
print("pathway rows:", q("SELECT COUNT(*) FROM chemical_pathway")[0],
      "| chemicals:", q("SELECT COUNT(DISTINCT chemical_id) FROM chemical_pathway")[0],
      "| families:", q("SELECT COUNT(DISTINCT family) FROM chemical_pathway")[0],
      "| (chemical,family) pairs:", q("SELECT COUNT(*) FROM (SELECT DISTINCT chemical_id,family FROM chemical_pathway)")[0])
print("ontology rows:", q("SELECT COUNT(*) FROM chemical_ontology")[0],
      "| chemicals:", q("SELECT COUNT(DISTINCT chemical_id) FROM chemical_ontology")[0])
if commit:
    c.commit(); print("COMMITTED")
else:
    c.rollback(); print("dry run - nothing written (use --commit)")
