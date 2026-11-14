#!/usr/bin/env python3
"""np23_raw_gene_set.py - the ligand gene set the chemical half RELEASES, rebuilt from raw
sources (owner decision 2026-11-09).

Rules (owner-chosen 'Rebuild from raw sources'):
  * Genes come from Human-GEM v2.1.1 and Rhea/UniProt only (ligand_gene_root, np22): these are
    the raw sources that give a role. GO, TCDB, Reactome and HMDB give no reliable role and do
    not contribute genes.
  * Roles synthesis, degradation, efflux, influx from reactions with a stated direction
    (reversible = 0), plus precursor_synthesis (Human-GEM, directed). Links that exist only
    through reversible / direction-less reactions are left out.
  * Cofactor hubs (metabolite_enzyme_support.hub_cofactor, e.g. ATP, phosphate, NAD) get NO genes:
    hundreds of enzymes make or use them, so their genes are not an expression proxy.
  * FALLBACK (owner, 2026-11-09): a non-hub chemical that had genes in the aggregated
    metabolite_enzymes table but gets none from the raw sources keeps its old genes, with
    origin = 'lr_fallback' (not confirmed by a raw source).
  * HMDB-only rows are dropped (owner decision); they only survive inside the fallback.
  * also_listed_by: the LR resources whose aggregated rows name the same gene for the chemical.

metabolite_enzymes (the aggregated table) is NOT modified: it stays as provenance.
Output table: ligand_gene_set.  Usage: python3 np23_raw_gene_set.py [--commit]
"""
import collections, sqlite3, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import np22_gene_roots as np22                                      # noqa: E402

HERE = Path(__file__).resolve().parent
CORE = ("synthesis", "degradation", "efflux", "influx")

def main():
    c = sqlite3.connect(HERE / "curation.db")
    norm = np22.hgnc_map()
    hubs = {x for (x,) in c.execute("SELECT DISTINCT chemical_id FROM metabolite_enzyme_support WHERE hub_cofactor=1")}
    hubs |= {x for (x,) in c.execute(
        f"SELECT chemical_id FROM chemical_entities WHERE name IN ({','.join('?' * len(np22.HUB_NAMES))})",
        sorted(np22.HUB_NAMES))}

    # LR resources naming each (chemical, gene) in the aggregated table
    listed = collections.defaultdict(set)
    agg = collections.defaultdict(list)
    for cid, gene, role, src, sup_status in c.execute(
            """SELECT m.chemical_id, m.gene_symbol, m.role, m.source_id, COALESCE(s.status,'')
               FROM metabolite_enzymes m LEFT JOIN metabolite_enzyme_support s ON s.me_id = m.id"""):
        g = norm(gene or "")
        if g:
            listed[(cid, g)].add(src)
            agg[cid].append((g, role, src, sup_status))

    raw = collections.defaultdict(lambda: {"src": set(), "ref": set()})
    for cid, g, role, src, ref in c.execute(
            """SELECT chemical_id, gene_symbol, role, source, ref_id FROM ligand_gene_root
               WHERE source IN ('human_gem','rhea_uniprot') AND reversible = 0
                 AND (role IN ('synthesis','degradation','efflux','influx')
                      OR (role = 'precursor_synthesis' AND source = 'human_gem'))"""):
        if cid in hubs:
            continue
        k = (cid, g, role)
        raw[k]["src"].add(src)
        raw[k]["ref"].add(ref)

    rows = []
    for (cid, g, role), v in raw.items():
        rows.append((cid, g, role, "raw", "|".join(sorted(v["src"])),
                     "|".join(sorted(v["ref"]))[:500], "|".join(sorted(listed.get((cid, g), ()))), 0))
    raw_chems = {k[0] for k in raw}
    fallback_chems = sorted(cid for cid in agg if cid not in raw_chems and cid not in hubs)
    for cid in fallback_chems:
        seen = set()
        for g, role, src, st in agg[cid]:
            if (g, role, src) in seen:
                continue
            seen.add((g, role, src))
            rows.append((cid, g, role, "lr_fallback", src, "", "|".join(sorted(listed.get((cid, g), ()))), 0))
    for cid in sorted(hubs):
        rows.append((cid, "", "", "hub_cofactor", "", "", "", 1))         # marker row, no gene

    c.execute("DROP TABLE IF EXISTS ligand_gene_set")
    c.execute("""CREATE TABLE ligand_gene_set(
        chemical_id TEXT NOT NULL, gene_symbol TEXT NOT NULL, role TEXT NOT NULL, origin TEXT NOT NULL,
        sources TEXT, ref_ids TEXT, also_listed_by TEXT, hub_cofactor INTEGER NOT NULL DEFAULT 0)""")
    c.executemany("INSERT INTO ligand_gene_set VALUES (?,?,?,?,?,?,?,?)", rows)

    q = lambda s, a=(): c.execute(s, a).fetchall()
    nm = dict(q("SELECT chemical_id, name FROM chemical_entities"))
    print("rows by origin:", q("SELECT origin, COUNT(*), COUNT(DISTINCT chemical_id) FROM ligand_gene_set GROUP BY 1"))
    print("raw rows by role:", q("SELECT role, COUNT(*) FROM ligand_gene_set WHERE origin='raw' GROUP BY 1 ORDER BY 2 DESC"))
    print("raw genes also listed by an LR db:",
          q("SELECT SUM(also_listed_by!=''), COUNT(*) FROM ligand_gene_set WHERE origin='raw'"))
    print("fallback chemicals:", [nm[x] for x in fallback_chems])
    print("hub chemicals (no genes):", [nm[x] for x in sorted(hubs)])
    old = set(agg); new = {r[0] for r in rows if r[3] in ("raw", "lr_fallback")}
    print(f"chemicals with genes: before {len(old)}, now {len(new)}; gained {len(new - old)}; "
          f"lost (all hubs?) {len(old - new)} -> {sorted(nm[x] for x in old - new)}")
    if "--commit" in sys.argv:
        c.commit(); print("COMMITTED")
    else:
        c.rollback(); print("dry run - nothing written (use --commit)")

if __name__ == "__main__":
    main()
