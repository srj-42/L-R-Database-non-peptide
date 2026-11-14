# Citation search methods — chemical (non-peptide) pairs

The reference for automating evidence search on `chemical_interactions`. Agreed with
the project owner on 2026-10-15. Every run follows this file. If a rule here changes,
record the change and its date in the change log at the bottom.

Scope for all four methods: **human, shippable pairs that have no confirmed evidence line yet**
(`np_evidence.verification_outcome IN ('confirmed','pair_absent')` absent). A pair that
already has a confirmed line is not re-searched.

Order of execution: **Method 4 → Method 3 → Method 1 → Method 2.**

Global rules (from `../non_peptides_prompt.txt` and PLAN.md invariants 8–12). These apply
to every literature method:

- **Exact molecule only.** No stereoisomers, analogues, salts, conjugates, or class
  statements such as "catecholamines" or "fatty acids".
- **The paper's own experiment.** No background, review, "known ligand", or others'-work
  sentences, and no comparator-only or reference-drug mentions. No indirect inference such
  as "effects abrogated by antagonist" or "pathway required".
- **Both names in the evidence sentence.** The sentence must name the molecule and the
  receptor, or a named subunit or complex that includes it. A rodent or primate ortholog
  is accepted as the same receptor.
- **Verbatim and from genuine text.** The sentence is quoted verbatim from genuine source
  text (abstract, or publisher full text). Curator prose is never written into
  `papers.full_text`.
- **Closed vocabularies.** `RESULT` is one of SUPPORTS, REFUTES, FAMILY_LEVEL, TRANSPORTER
  or ENZYME. `EXPERIMENT_TYPE` is one of structural, direct_binding, functional_in_vitro,
  functional_in_vivo or computational.
- **Licence.** MEBOCOST-only pairs are withheld, and `mebocost*` never appears in exported
  provenance.
- **Complex rule.** Components of one receptor complex or family group (`np13_groups.py`)
  are searched once. Evidence that names the complex covers every component row
  (`np13_propagate.py`, which marks each copied row as propagated).

Review: candidate sentences are proposed by reviewer subagents
(`work/np14/REVIEWER_PROMPT.md`). Claude checks every accepted record before loading.
Loaders re-check the verbatim text and the names mechanically.

Evidence classes shown on the site, highest rank first:

| class | from | confidence | ranked |
|---|---|---|---|
| `nomenclature_evident` | Method 4 | high | 1 (above literature, owner decision) |
| `confirmed` / `pair_absent` (literature) | Methods 1, 2 | sentence-verified | 2 |
| `database_annotated` | Method 3 | curated database with PMID | 3 |

---

## Method 1 — pair search: abstract, then Results + Discussion

### 1a. What has already been done (recorded as it actually ran)

1. **Upstream citations** (np8, np11). PMIDs that the source databases attached to a pair
   were read at abstract level, giving 2,178 rows. Each row was decided `confirmed`,
   `abstract_insufficient`, `review_only`, `different_entity`, `unresolvable` or
   `retracted`.
2. **Full-text audit of upstream citations** (np12). PMC open-access papers among the
   `abstract_insufficient` rows were fetched (359 papers) and read in full, not by section.
   Confirmed rows carry `evidence_depth = full_text`.
3. **Literature search for pairs without a citation** (np10, then np14):
   - The PubMed query is `(molecule synonyms)[tiab] AND (receptor symbol/HGNC aliases/full
     name, or subunit clause)[tiab] NOT review[pt]`. Above 24 hits, the query adds
     assay terms. Up to 12 hits per pair are kept.
   - **np10** (ranked, targeted and cluster batches): Claude read every abstract.
   - **np14** (automated first pass): abstracts are split into sentences. A pair goes to a
     reviewer only if one sentence names both the molecule and the receptor. Otherwise it is
     logged `No evidence` with the reason, in `work/np10_no_evidence.tsv`.
   - All pairs searched are logged in `work/np10_search_log.jsonl`. Every unsearched
     no-citation pair had been searched by 2026-10-13.
4. **ChEMBL check** (np13_chembl). ChEMBL activities with a PMID were matched to our pairs,
   and each (pair, PMID) was read at abstract level. 3,691 rows went in; 22 were confirmed,
   and most of the rest were `abstract_insufficient` because the data is in tables.

Full text was **not** read for papers found by np10/np14. Discussion-only reading was never
done before this plan.

### 1b. Full-text stage (new)

- **Input.** Papers already linked to an unverified pair, from two sets:
  - (i) ChEMBL-listed rows, `source_id='chembl'` and `abstract_insufficient`;
  - (ii) upstream-cited rows, `abstract_insufficient`, from source databases.
  The ~10,000 pairs closed by the np14 automatic filter are **not** included (owner
  decision).
- **Fetch.**
  1. PMC open access first (Europe PMC `fullTextXML`).
  2. Otherwise publisher full text through the connected Chrome browser and the owner's
     institutional access.
  3. Captcha, login wall, or no access: log `fulltext_unavailable`, never retry in a loop.
- **Read.** The **Results and Discussion** sections only. Sections are found from JATS
  `sec-type`/titles, or from page headings in the browser. The Discussion is read with
  extra care about background sentences.
- **Accept.** Same global rules. The row becomes `confirmed` with
  `evidence_depth = full_text`. A sentence from a browser-read page must be quoted
  verbatim, and the source URL is recorded.
- **Store.** Genuine publisher text only goes into `papers.full_text`. Browser-read text is
  kept in `work/` with its URL and is not written to `papers.full_text` unless it is the
  complete publisher text.

## Method 2 — ligand-only search

- **Input.** The 625 unique ligands that have at least one unverified pair.
- **Query.** `(molecule synonyms)[tiab] AND (receptor[tiab] OR binding[tiab] OR
  agonist[tiab] OR antagonist[tiab] OR EC50[tiab] OR IC50[tiab] OR Ki[tiab] OR
  Kd[tiab]) NOT review[pt]`. Take the **top 200 by relevance** per ligand.
- **Read.** For each paper, look in the abstract and, where full text can be fetched
  (same fetch order as Method 1b), the **Results and Discussion**. Look for **any**
  receptor in that ligand's unverified pairs, using symbol, HGNC aliases, previous symbols,
  full name and subunit forms.
- **Accept.** Same global rules. The row is loaded with `source_id='np16_ligand_search'`.
- **Receptors outside our table.** If a paper shows the ligand acting on a receptor not
  paired with it in our table, log it to `work/np16_candidate_new_pairs.tsv` (ligand,
  receptor, PMID, sentence). It is **not** added to the table.

## Method 3 — protein and structured databases

- **Allowed.** UniProt, Rhea, PDB (ligand-bound structures), TCDB and BindingDB.
- **Not allowed.** DrugBank and BRENDA (licence).
- **Already aggregated.** GtoPdb and ChEMBL; no need to repeat them.
- **Match.** Map database proteins by UniProt accession to `entities.uniprot`. Map
  molecules by ChEBI, PubChem or HMDB id first, and name only as a fallback.
- **Accept.** An annotation stating this molecule binds, activates, inhibits or is
  transported or turned over by this protein, **with a PMID**, gives a
  `database_annotated` record. The record holds the database, accession, annotation text
  and PMID. PDB structures give the PDB id and the primary-citation PMID.
- **Then.** Try to extract a verbatim sentence from the cited paper (Method 1 rules). If
  one is found, the pair also gets a normal `confirmed` literature row.
- **Store.** Table `np_class_evidence` with `evidence_class='database_annotated'`. These
  rows never go into `np_evidence`, because they are not sentence-verified.

## Method 4 — name-evident pairs (exact name only)

- **Rule.** The pair is `nomenclature_evident` when the receptor's **official HGNC full
  name** contains the molecule's **own name**, or one of its recorded synonyms that names
  the same molecule, as a whole word or phrase. Examples: histamine with "histamine
  receptor H1"; adenosine with "adenosine A1 receptor"; glycine with "glycine receptor
  alpha 1". There is **no class map**, by owner decision on 2026-10-15. So thyroxine with
  "thyroid hormone receptor beta" does **not** qualify.
- **Exclusions.** A name match that is not a ligand–receptor statement is rejected, and the
  reason is logged. Examples:
  - "zinc finger protein …" (zinc)
  - "ATP binding cassette …" (ATP as a hydrolysed cofactor, not the ligand searched)
  - ion names inside channel names where the ion is only the permeant species, such as
    "calcium voltage-gated channel …" or "sodium voltage-gated channel …"
  - "… binding protein" names that describe a different partner.

  Every match is listed in `work/np15_nomenclature.tsv` with `accept` or `reject` and a
  reason, before loading.
- **Confidence.** High. Ranked **above** literature evidence on the site. No PMID is
  needed. A literature citation found later is shown underneath, and the class stays.
- **Store.** Table `np_class_evidence` with `evidence_class='nomenclature_evident'`,
  `source='HGNC approved name'`, and the matched name string.

---

## Files and scripts

| step | script | log |
|---|---|---|
| upstream audit | `np11_cited_audit.py`, `np9_verify.py` | `work/np11_audit.tsv` |
| full-text audit (PMC) | `np12_build.py` | `work/np12_audit.tsv` |
| pair search | `np10_litsearch.py`, `np14_auto.py`, `np14_chunks.py` | `work/np10_search_log.jsonl`, `work/np10_no_evidence.tsv` |
| ChEMBL | `np13_chembl.py` | `work/np13_chembl.jsonl` |
| complex rule | `np13_groups.py`, `np13_propagate.py` | `work/np13_groups.json` |
| Method 4 | `np15_nomenclature.py` | `work/np15_nomenclature.tsv` |
| Method 3 | `np17_databases.py` (to write) | `work/np17_*.tsv` |
| Method 1b | `np18_fulltext.py` (to write) | `work/np18_audit.tsv` |
| Method 2 | `np16_ligand_search.py` (to write) | `work/np16_*.tsv` |

After every load: `np13_groups.py` → `np13_propagate.py --commit` → `np6_export.py` →
`np7_build_site.py` → `cur_integrity_check.py` (must print `INTEGRITY: CLEAN`).

## Change log
- 2026-10-15: file created. Owner decisions:
  - Method 1: Results + Discussion; browser for paywalled papers; scope is the ChEMBL and
    upstream-cited unverified rows.
  - Method 2: assay terms, top 200 per ligand, unverified pairs only, new pairs logged
    separately.
  - Method 3: UniProt, Rhea, PDB, TCDB and BindingDB, as a separate class.
  - Method 4: exact name only, ranked above literature.
  - Order 4→3→1→2; subagent review plus Claude's spot-check.
