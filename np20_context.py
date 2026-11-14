#!/usr/bin/env python3
"""np20_context.py - metabolite-level context annotations for the chemical half.

Adds the table chemical_context (owner-approved additive exception, 2026-11-05):
    chemical_id, kind, value, source
kind in {tissue, cell_location, biospecimen, disease}.

UPDATE 2026-11-06: the preferred source is now HMDB 5.0 itself, parsed by
np20_hmdb_parse.py into work/np20/hmdb_context.json (source tag 'hmdb'). Run with
--hmdb to load it; that REPLACES the 'hmdb_via_metalinksdb' rows. The text below
describes the earlier MetalinksDB-based route.

SOURCE. HMDB records these per metabolite, but hmdb.ca blocks scripted downloads
(Cloudflare challenge), so the values are read from raw/metalinksdb.db, whose
tissue_location / cell_location / biospecimen_location / disease tables were built
from HMDB (MetalinksDB paper, Brief Bioinform 25:bbae347). source is therefore
'hmdb_via_metalinksdb'. A direct HMDB XML parse would be fuller (see report).
Matching is on chemical_entities.hmdb_id.

NOT loaded: MetalinksDB's `pathway` table (partly corrupted, PHASE_NONPEPTIDE.md).

These annotate the METABOLITE, not the ligand-receptor pair. They do not feed any
verification status. Usage:  python3 np20_context.py [--commit]
"""
import sqlite3, sys

DB = "curation.db"
SRC = "raw/metalinksdb.db"
TAG = "hmdb_via_metalinksdb"
KINDS = {"tissue": ("tissue_location", "tissue_location"),
         "cell_location": ("cell_location", "cell_location"),
         "biospecimen": ("biospecimen_location", "biospecimen_location"),
         "disease": ("disease", "disease")}

commit = "--commit" in sys.argv
use_hmdb = "--hmdb" in sys.argv
c = sqlite3.connect(DB)
c.execute(f"attach 'file:{SRC}?mode=ro' as m")
c.execute("""CREATE TABLE IF NOT EXISTS chemical_context(
    chemical_id TEXT NOT NULL, kind TEXT NOT NULL, value TEXT NOT NULL, source TEXT NOT NULL,
    PRIMARY KEY (chemical_id, kind, value, source))""")
if use_hmdb:
    import json
    TAG = "hmdb"
    d = json.load(open("work/np20/hmdb_context.json"))
    old = {(a, b, v) for a, b, v in c.execute("SELECT chemical_id,kind,value FROM chemical_context WHERE source='hmdb_via_metalinksdb'")}
    c.execute("DELETE FROM chemical_context WHERE source IN ('hmdb','hmdb_via_metalinksdb')")
else:
    c.execute("DELETE FROM chemical_context WHERE source=?", (TAG,))  # idempotent reload

rows = []
if use_hmdb:
    for cid, hid in c.execute("SELECT chemical_id, hmdb_id FROM chemical_entities WHERE hmdb_id IS NOT NULL AND hmdb_id!=''").fetchall():
        rec = d.get(hid)
        for kind in KINDS:
            for val in (rec or {}).get(kind, []):
                rows.append((cid, kind, val, TAG))
    new = {(a, b, v) for a, b, v, _ in rows}
    print(f"previous MetalinksDB-derived rows {len(old)}; kept by HMDB {len(old & new)}; "
          f"only in old {len(old - new)}; new from HMDB {len(new - old)}")
for kind, (tbl, col) in ([] if use_hmdb else KINDS.items()):
    q = f"""SELECT DISTINCT e.chemical_id, TRIM(x.{col}) FROM chemical_entities e
            JOIN m.{tbl} x ON x.hmdb = e.hmdb_id
            WHERE x.{col} IS NOT NULL AND TRIM(x.{col}) != ''"""
    for cid, val in c.execute(q).fetchall():
        rows.append((cid, kind, val, TAG))
c.executemany("INSERT OR IGNORE INTO chemical_context VALUES (?,?,?,?)", rows)

print("rows to write:", len(rows))
for kind in KINDS:
    n, k = c.execute("SELECT COUNT(DISTINCT chemical_id), COUNT(*) FROM chemical_context WHERE kind=?",
                     (kind,)).fetchone()
    print(f"  {kind:14s} chemicals {n:4d}  rows {k:5d}")
print("chemicals with any context:",
      c.execute("SELECT COUNT(DISTINCT chemical_id) FROM chemical_context").fetchone()[0],
      "of", c.execute("SELECT COUNT(*) FROM chemical_entities").fetchone()[0])
if commit:
    c.commit(); print("COMMITTED")
else:
    c.rollback(); print("dry run - nothing written (use --commit)")
