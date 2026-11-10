# Lean-mass genetic estimands | Osteosarcopenia Causal Atlas

> **Licensing:** Original analysis code and model syntax are licensed under MIT. Author-owned aggregate result tables and original research artwork are licensed under CC BY 4.0. Different files have different licence scopes; consult `LICENSE_SCOPE.md` and `DATA_LICENSE.md`. Original GWAS records, LD reference data and controlled model inputs are **not** distributed.

This repository accompanies the Original Article **“Lean-Mass Decomposition Reveals Divergent Genetic Correlation and Covariance Rankings”**, prepared for the *Journal of Cachexia, Sarcopenia and Muscle*.

## Scientific scope

A model-defined regional fat-free-mass (FFM) factor and its fat-mass-orthogonal genetic residual were compared against adjusted appendicular lean mass using genetic correlation and unstandardized genetic covariance. The analysis reports opposite rankings of these constructs across the two genetic estimands. The supplied materials document selected model syntax, code, plotting inputs and publication figures; they do not perform a fresh GWAS, LDSC, or primary GenomicSEM analysis.

## What is included and reproducible

| Component | Included | What can be verified |
|---|---|---|
| Submitted Figure 3 | Plotting source table and `draw_figure3.py` in this repository | Rebuild Figure 3 from the supplied derived plotting records; this is not an independent model refit |
| Model syntax | Two original lavaan/GenomicSEM syntax files | Inspect model definitions and derived parameters; original S/V/I/N inputs are not distributed |
| DeltaCov numerical verification | Archived R verification script | Review the computation; execution requires separately controlled original model objects and matrices |
| Derived result tables | Selected aggregate/model-derived summaries | Inspect selected numerical summaries and validate selected relationships |
| Figure 2 SVG | Editable artwork | Inspect the author-created CC BY 4.0 SVG; this is an artwork release, not an independent GWAS or model reconstruction |
| Smoke test | Repository hygiene test | Basic file and path checks only |
| Content validator | Release-content validator | Checks selected frozen table relationships; it does not independently estimate GWAS or genetic models |

## Run the checks

From the repository root:

```bash
python code/smoke_test.py
python code/release_validate.py
```

For Figure 3, using the supplied plot source table:

```bash
cd code/figures/submitted_figure3
python draw_figure3.py
```

Plotting dependencies are described in `REQUIREMENTS.txt`. Script outputs should be checked against the submitted figure. The R verification script is review/verification source code, not a turnkey command for an independent rerun: its input paths refer to author-controlled inputs that are not distributed here.

## Analysis keys and reporting

`D` is a machine-readable key for the primary specification, while `R75` and `R50` identify 75% and 50% assumed-N numerical sensitivity specifications. These labels do not indicate independent cohorts, replications, or successive published models. They are retained in selected machine-readable tables for provenance matching.

## Availability and credits

The data and results originate from published GWAS sources cited in the manuscript and Supplementary Table S6a. Original SNP-level GWAS/VCF records, LD reference panels, complete LDSC S/V/I/N matrices, primary model RDS objects, reviewer-only historical archives and excluded analyses are not included. The right- and left-side models share ALM and bone GWAS inputs, and the FM-orthogonal residual is a statistical construct rather than measured skeletal muscle.

## Licensing

The root `LICENSE` (MIT) applies **only** to the six original analysis-code/model-syntax files enumerated in `LICENSE_SCOPE.md`. The author-owned aggregate result TSVs and original Figure 2/Figure 3 artwork in this repository are offered under CC BY 4.0, subject to the scoped notice in `DATA_LICENSE.md`. Third-party GWAS records, LD panels, and restricted model inputs are not distributed or sublicensed. The repository-level `CITATION.cff` intentionally does not declare one global SPDX licence because different files have different licences.

**Reproducibility boundary:** the public files support selected Figure 3 regeneration and numerical cross-checks, but not independent reruns of the upstream GWAS, LDSC or GenomicSEM fitting.
