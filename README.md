# Lean-mass genetic estimands — R2 manifest composites

**READY_FOR_R2_RELEASE — reduced public result scope, 2026-11-15.** This branch provides the original R2 analysis code, synthetic fixtures and four verified original aggregate result CSVs. The author-scoped public grant is recorded; numerical closure passed. Local readiness is not a claim that a commit, tag or public release has already been created.

The four R2 result files in data/derived_r2_candidate/ are:

- weight_grid.csv: all 25 regional weight combinations on each side, including paired ALM contrasts and variance/covariance retention.
- genetic_sd.csv: genetic-SD-standardized ALM sensitivity.
- partial_rg.csv: ALM partial genetic correlations conditional on M.
- alm_restrictions.csv: six ALM necessary-restriction diagnostics, not software SEM fit indices.

The complete composite.csv, bone_global.csv and bone_pairwise.csv are intentionally absent and archived privately. Their FN/LS source-license question is **EXCLUDED_PENDING_CLARIFICATION**, not resolved or newly authorized. Names and hashes in r2/EXCLUDED_R2_RESULTS.tsv describe exclusions only. Historical fixed-loading v1.0.0 assets remain unchanged and are not current R2 result files.

R2 defines F=0.5(ArmFFM+LegFFM), M=0.5(ArmFM+LegFM), gamma=Cov_g(F,M)/Var_g(M), U=F−gamma M. Paired Delta inference uses full V and differentiates gamma, SD scaling and denominators. All 25 weight pairs are retained per side. U is a statistical projection, not pure muscle or a causal intervention; right/left are not independent cohort replications.

All recovered numerical code is retained, including optional bone analyses and private-input validators. Public synthetic checks require none of the excluded study CSVs:

```sh
python tests/test_r2_smoke.py
python code/smoke_test.py
python code/release_validate.py
python code/r2/run_r2.py --format synthetic-json --input-dir data/synthetic_r2 --output-dir local_synthetic_output
```

The runner generates outputs from inputs; it does not read the excluded study CSVs. Synthetic outputs may include the full analysis families but contain no study values. Reproduction of full study results requires authorized real frozen RDS and private references; any resulting FN/LS outputs remain subject to the unresolved source conditions. See r2/README.md.

The previously completed V5 audit covered 798 mapped components including overlapping 624 scalar checks. Prior 3,466 frozen-reference and 940 independent R comparisons passed. Those full-audit counts include results excluded from this public selection. Of the original 798 components, 546 map to the four retained CSVs; this is retained audit evidence, not a new numeric audit or a claim that all workbook inputs are public.

File-specific licenses: original R2 code MIT; the four original R2 aggregate CSVs CC BY 4.0; enumerated original documentation, synthetic fixtures and metadata CC BY 4.0. Historical grants remain unchanged. The root MIT text does not license the entire repository. See LICENSE_SCOPE.md and PUBLICATION_RECOMMENDATIONS.tsv.

Raw GWAS/VCF, participants, real full S/V/I/N, controlled RDS, LD panels, private workbooks/maps, approval records, audit archives and credentials are excluded. Historical fixed-loading code/data/artwork is preserved byte-for-byte in its stated scope. main and v1.0.0 remain at 0795a68f5b2b9a4cc914d516fffd89c535370617. See DATA_AVAILABILITY.md, REPRODUCIBILITY_LIMITATIONS.md and RELEASE_NOTES.md.
