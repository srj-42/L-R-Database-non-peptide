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
| `np14_auto.py`, `np14_chunks.py` | Automated first pass. For each pair it runs the np10 PubMed query and splits the abstracts into sentences. A pair with no sentence naming both the molecule and the receptor is logged "No evidence" with the reason. Everything else goes to reviewers. |
| `np_common.py` | Shared guards, written after a silent-failure incident. |
