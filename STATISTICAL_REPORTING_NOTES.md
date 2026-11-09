# Statistical reporting notes

- The machine-readable model identifier `D` corresponds to the primary analysis. `R75` and `R50` identify 75% and 50% assumed-N numerical sensitivity specifications, not independent replication. These internal keys remain in archived source-derived tables to preserve their identity.
- `MODEL_FIT_ALL_SCENARIOS.tsv` includes descriptive fit indices. The call `sample.nobs=2` means its software chi-square objective is not a conventional full-sample likelihood-ratio test; the removed AIC and likelihood-ratio P values should not be inferred.
- `FIGURE2_DATA.tsv` includes marginal ALM correlations/covariances and paired correlation contrasts. `DELTACOV_VERIFICATION_SUMMARY.tsv` records the separately verified full-V Delta-method covariance contrast, not a public rerun of the source GWAS.
- `BONE_OUTCOME_SHIFTS_DATA.tsv` contains six primary bone-outcome shift rows corresponding to the manuscript's Figure 4. The retained Figure 3 plotting data are in `code/figures/submitted_figure3/FIGURE3_PLOT_SOURCE.tsv`.
- The analysis is model-conditional. Regional FFM is a fixed-loading factor, and the FM-orthogonal residual is not an anatomical muscle measurement.
- Original GWAS/VCF files, full S/V/I/N, model RDS and provenance-restricted evidence are not included.
