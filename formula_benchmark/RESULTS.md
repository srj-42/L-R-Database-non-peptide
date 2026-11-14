# NPL availability-formula benchmark: run 1 (2026-11-09)

Question: which formula best picks out the KNOWN producer tissues / cell types of each ligand,
with the gene sets held fixed (`ligand_gene_set`, np23)?

- Expression: HPA v25 `rna_tissue_consensus` (50 tissues) and `rna_single_cell_type` (153 cell types), in `dl/`.
- Truth: `truth.py` (37 ligands at tissue level, 44 at cell-type level; written before any scoring).
- Metric: per-ligand AUROC, producers vs all other samples. Paired Wilcoxon against the site formula.
  Null = same formula on random genes with the same role sizes (200 draws).
- Scaling: log2(x+1) and per-gene percentile rank across samples.
- Code: `run_bench.py` (run with `spice/.venv/bin/python -I`). Outputs: `summary.csv`, `per_ligand.csv`.

## Headline (mean AUROC; log scale; tissues / cell types)
| formula | tissues | cell types |
|---|---|---|
| site formula (usage.js) | 0.836 | 0.738 |
| same, no degradation (MetaLigand-like) | 0.878 | 0.858 (p=0.003) |
| degradation as ratio S*S/(S+D) | 0.878 (p=0.03) | 0.838 (p=0.001) |
| synthesis mean only | 0.892 | 0.879 (p=0.001) |
| CellPhoneDB v5 sqrt(syn*eff) | 0.889 | 0.878 |
| MEBOCOST syn-deg | 0.639 (p=0.02) | 0.653 |
| random genes (null) | 0.53 | 0.47 |
With per-gene rank scaling the best variants were "synthesis max only" (0.902 / 0.894) and
"max, no degradation" (0.895 / 0.883); the site formula was 0.735 / 0.686.

## Findings
1. Subtracting mean(degradation) is the main problem. It clamps producers to 0 (8 of 44 ligands at
   cell level: GABA in inhibitory neurons AUROC 0.02, noradrenaline in adrenal medulla 0.49 because
   PNMT/MAO/COMT are high there, cortexolone, deoxycorticosterone, serotonin, DHEA-S, S1P, PGI2).
   Every formula that subtracts (ours, MEBOCOST, MRCLinkdb) ranks last.
2. Extra multiplicative steps (precursor, salvage, efflux) do not help on average; synthesis alone
   is as good or better. The efflux gate is neutral.
3. Gene-set defects (independent of formula), to review in np23:
   - Creatine: synthesis = creatine kinases (CKB/CKM/CKMT*, which make phosphocreatine), while GAMT
     (the real creatine synthase) is listed as DEGRADATION. Looks like a direction error.
   - Estriol: synthesis = hepatic CYP2C8/CYP3A4 only; placental route (16a-OH-DHEA-S + aromatase) absent.
   - GABA precursor_synthesis: ~50 genes incl. cathepsins, GGTs, ALB (protein-degradation reactions);
     estradiol degradation: ~45 hepatic CYPs/UGTs. Means over such lists are mostly noise.
4. The formula beats the random-gene null for only 20/44 (cell) and 22/37 (tissue) ligands at p<0.05.

## Limits
- Tests WHERE a ligand is made (specificity), not HOW MUCH. Truth set leans to classic hormones /
  transmitters, which favours synthesis-centric formulas. Small n (37/44 ligands).
- Projection transmitters (dopamine made in midbrain, released in striatum) cannot be scored by place.
- Next: CCLE metabolomics (225 metabolites x ~900 lines; quantity) needs a browser download from
  DepMap (scripts are bot-checked); Vicari MALDI+Visium (Mendeley, scriptable) for spatial.
