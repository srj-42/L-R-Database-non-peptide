# L-R-Database-non-peptide
| Script | Role |
|---|---|
| `np0_acquire.py` | Pin the non-peptide sources. |
| `np1_parse_composites.py`, `np1b_fill_review.py` | Split composite ligand strings (a ligand plus the enzymes that make it) and fill the manual review queue. |
| `np2_chemicals.py`, `np2b_synonyms.py`, `np2c_add_abbreviations.py`, `np_dedup_chemicals.py` | Build `chemical_entities` keyed on HMDB; apply curated synonyms and abbreviations; merge name-keyed duplicates. |
| `np3_receptors.py`, `np4_union.py`, `np5_enzymes.py` | Resolve receptors into entities; collapse the sources into `chemical_interactions`; build the enzyme relation under a five-role vocabulary. |
| `np6_export.py`, `np7_build_site.py` | Export `chemical_*.csv`; build the chemical browse site. |
| `np8_evidence.py`, `np9_verify.py` | The evidence layer for the chemical half; the verification loop. |
| `np10_litsearch.py`, `np11_cited_audit.py` | Literature search against `non_peptides_prompt.txt` for pairs with no evidence; adjudicate every pending cited row. |
| `np13_chembl.py`, `np13_groups.py`, `np13_propagate.py` | Pull ChEMBL activities; group multi-subunit receptors into one (chemical, complex) unit; apply that complex rule retroactively. |
| `np14_auto.py`, `np14_chunks.py` | Automated first pass. For each pair it runs the np10 PubMed query and splits the abstracts into sentences. A pair with no sentence naming both the molecule and the receptor is logged "No evidence" with the reason. |
| `np_common.py` | Shared guards, written after a silent-failure incident. |
| `np15_nomenclature.py` | Method 4: the receptor's official name contains the molecule's name. |
| `np16_ligand_search.py` | Method 2: ligand-only literature search. |
| `np17_databases.py`, `np19_databases.py` | Method 3: curated databases. `np19` redoes it for all pairs: ChEMBL, BindingDB, GtoPdb, PDB, UniProt/Rhea, TCDB, stored as `database_annotated` class evidence. |
| `np18_prepare.py`, `np18_load.py` | Method 1b: full-text work list and loader for browser-read records. |
| `np20_hmdb_parse.py`, `np20_context.py`, `np21_pathways_ontology.py` | HMDB 5.0 context for each molecule: tissue, biospecimen, cell location, disease, pathways, ontology. Describes the molecule, never evidence for a pair. |
| `np22_gene_roots.py` | Rebuilds each ligand's gene sets from root sources (Human-GEM v2.1.1, Rhea through UniProt, Reactome, GO, TCDB) and marks which aggregated rows they support. |
| `np23_raw_gene_set.py` | The released gene set: Human-GEM and Rhea/UniProt directed roles only; cofactor hubs get no genes; chemicals with no raw genes keep their aggregated genes, marked as fallback. |
| `build_gene_links.py` | Builds the gene-card link data for the sites. |
| `formula_benchmark/` | Test of the availability formula against known producer tissues and cell types (Human Protein Atlas): `truth.py`, `run_bench.py`, `RESULTS.md`, `summary.csv`, `per_ligand.csv`. |
| `np_schema.sql`, `sources_np.json`, `np_composite_map.tsv`, `np_synonym_overrides.tsv`, `np_identity_conflicts.tsv` | Schema, pinned source list, the reviewed composite-ligand map, curated synonyms, identity conflicts found while merging. |
| `non_peptides_prompt.txt`, `candidate_generation_non_peptide.txt`, `CITATION_METHODS.md` | The literature-search prompt, the candidate-pair prompt, and the four-method evidence search. |

