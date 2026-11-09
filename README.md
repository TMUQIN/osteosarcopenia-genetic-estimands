# Lean-mass genetic estimands | Osteosarcopenia Causal Atlas

> **PRIVATE STAGING COPY — NOT APPROVED FOR PUBLIC RELEASE.** Do not change a GitHub repository containing these files to Public until the coauthors, institutional rights holders, source-data terms, and software/data licences have been confirmed. This notice is intentionally retained during preparation and must be removed from a separately verified public release.

This repository is prepared to accompany the Original Article **“Lean-Mass Decomposition Reveals Divergent Genetic Correlation and Covariance Rankings”**, intended for the *Journal of Cachexia, Sarcopenia and Muscle*.

## Scientific scope

A model-defined regional fat-free-mass (FFM) factor and its fat-mass-orthogonal genetic residual were compared against adjusted appendicular lean mass using genetic correlation and unstandardized genetic covariance. The published analysis reports opposite rankings of these constructs on the two genetic estimands. The supplied material documents selected model syntax, code and publication figures; it does not perform a fresh GWAS, LDSC, or primary GenomicSEM analysis.

## What is included and reproducible

| Component | Included | What can be verified |
|---|---|---|
| Submitted Figure 3 | Plotting source table and `draw_figure3.py` in the FULL variant | Rebuild Figure 3 from frozen, derived plotting records; not an independent refit |
| Model syntax | Two original lavaan/GenomicSEM syntax files | Inspect model definitions and derived parameters; no original S/V/I/N distributed |
| DeltaCov numerical verification | Archived R verification script | Review computation; running requires separately controlled original model objects and matrices |
| Derived result tables | FULL variant only | Inspect selected numerical summaries and validate a subset of relationships |
| Figure 2 SVG | FULL variant only | Open the editable art; data and figure rights require confirmation |
| Smoke test | Both variants | Basic file and path hygiene only |
| Content validator | FULL variant | Checks selected frozen table relationships; does not independently estimate GWAS or genetic models |

## Run the checks

From the repository root:

```bash
python code/smoke_test.py
python code/release_validate.py  # FULL variant only
```

For Figure 3, when the FULL plotting source has been lawfully included:

```bash
cd code/figures/submitted_figure3
python draw_figure3.py
```

Plotting dependencies are described in `REQUIREMENTS.txt`. Script outputs should be checked against the submitted figure. The R script is **review-only verification source code**, not a turnkey command for an independent rerun: its input paths refer to an author-controlled archive not distributed here.

## Analysis keys and report wording

`D` is a machine-readable key for the **primary specification**, while `R75` and `R50` are 75% and 50% assumed-N numerical sensitivity specifications. These labels remain in certain raw derived TSVs to permit exact provenance matching; they do not indicate independent cohorts, replications, or successive published models. Historical version identifiers in retained source code are not claims of a new analysis.

## Availability and credits

The data and results originate from published GWAS sources cited in the manuscript and Supplementary Table S6a. Original SNP-level GWAS/VCF records, LD reference panels, complete LDSC S/V/I/N matrices, primary model RDS objects, reviewer-only historical ZIPs, and excluded analyses are not included. The right- and left-side models share ALM and bone GWAS inputs, and the FM-orthogonal residual is a statistical construct rather than measured skeletal muscle.

The data/source restrictions in `DATA_AVAILABILITY.md` and the limitations in `REPRODUCIBILITY_LIMITATIONS.md` apply. This staging archive has **no selected LICENSE, public GitHub URL or Zenodo DOI**. For a final public release, obtain author and any required source rights approvals, choose explicit compatible licences, update citation metadata, verify the approved files, and publish a fixed version.
