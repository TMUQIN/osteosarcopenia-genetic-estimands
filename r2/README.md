# Running the recovered R2 calculations

The original mathematical functions are copied verbatim from the real local R2 closure and full-precision Wald scripts. SOURCE_FUNCTION_MANIFEST.tsv records source identity, line spans, and extracted text SHA256. These are not the old fixed-loading scripts relabelled as R2. run_r2.py is a newly written portable I/O adapter; it is not represented as the historical runner. verify_r2_independent.R changes only input/output paths from the retained independent R script. Complete original local runners are inventoried in the private audit and are not copied with their historical paths or reviewer-archive dependencies.

Trait order: ALM, ArmFFM, LegFFM, ArmFM, LegFM, eBMD, FNBMD, LSBMD. S is 8x8 and V is 36x36. V follows column-major lower-triangle vech, equivalently (S11,S21,...,S81,S22,...,S88). The five-trait V block must be selected on both axes with 1-based indices **1,2,3,4,5,9,10,11,12,16,17,18,22,23,27**. Taking the top-left 15x15 is wrong. VECH_36_MAPPING.tsv contains labels/indices only, not source matrix values.

## Synthetic smoke run

```sh
python tests/test_r2_smoke.py
python code/r2/run_r2.py --format synthetic-json --input-dir data/synthetic_r2 --output-dir local_synthetic_output
```

The synthetic fixtures are random positive definite matrices generated with disclosed seeds. They are not calibrated to the study and do not reproduce paper numbers.

## Author-input reproduction

Obtain permission for the two frozen RDS files. Their exact hashes are in FROZEN_INPUT_IDENTITIES.tsv. Each file must contain S8, V36, and the exact S column names. The adapter reads binary doubles through R sprintf('%.17g'); it does not substitute rounded TSV exports. The public folder deliberately contains no real RDS.

```sh
python code/r2/run_r2.py --input-dir /authorized/frozen/rds --output-dir local_author_output
Rscript code/r2/verify_r2_independent.R /authorized/frozen/rds /existing/independent_output
python code/r2/validate_frozen_references.py --project-root /authorized/original/project --recomputed-dir local_author_output --independent-r-dir /existing/independent_output --output-dir local_reference_validation
```

No GWAS, LDSC, GenomicSEM fitting, or lavaan fitting is called by these commands. R jsonlite and numDeriv are required for RDS reading / independent derivatives. Python versions used are in logs/ENVIRONMENT.json. CI definitions are deliberately retained: the original all-outcome Wald helper uses 1.96; the later grid/SD closure uses norm.ppf(0.975). Do not silently alter these during matching.

The explicit frozen-reference validator is newly written verification code. It needs the original author project reference TSVs identified in LOCAL_REFERENCE_IDENTITIES.tsv; these private historical files are not redistributed. Its detailed execution output is preserved privately; the public summary is logs/FROZEN_NUMERIC_VALIDATION.json. The verified public GitHub repository was cloned at main/tag v1.0.0 commit 0795a68f5b2b9a4cc914d516fffd89c535370617. Branch r2-manifest-composites is local and based on that commit. No commit, push, Release or DOI has been created. See BASELINE_IDENTITY.json; the exact staged diff is retained privately.

## Explicit submitted-Excel validation

```sh
python code/r2/validate_r2_excel.py --excel /author/designated_R2_15_tables.xlsx --mapping /reviewed/main_results_cell_map.tsv --derived-dir data/derived_r2_candidate --output excel_validation.json
```

The TSV mapping columns are: table, sheet, cell, csv_file, key_json, field, absolute_tolerance, relative_tolerance, precision_basis. table must cover S2/S3/S3b/S3c/S4/S5; key_json is an exact CSV row selector such as {"side":"right","metric":"ALM_delta_rg"}. Every main-result numeric cell must be mapped and reviewed for complete coverage. Tolerances must be justified by the actual table's printed precision before comparing. The real author-provided map (SHA256 a02d59b3f1d42f1fd0630088161572dd7b589b80a1f95cabe958438ff13c8602) and V5 workbook (SHA256 783861e83479d4feb5a6cf07271c979df010ca6f9edc88f071c56073b17ebcf4) were used locally. The original validator passed 624 mapped scalar cells. A separate private audit checked printed estimate(SE), interval components and identities, with independent Decimal arithmetic passing all 798 components. The supplied script was also rerun. The private map and workbook are not distributed. Without these inputs the validator still exits nonzero and writes NOT_VERIFIED. It does not infer table equivalence from matching numbers or sheet counts.

## Tests actually executed here

9 smoke tests passed. Fresh frozen-input output passed 3,466 keyed numeric comparisons against retained R2 full-precision reference tables, plus 940 comparisons for 470 independent R rows. Existing local reference tolerance is rtol=1e-10, atol=1e-12; the independent Richardson check uses rtol=1e-8, atol=1e-11. V5 workbook agreement was independently closed on 2026-11-13 for 798 mapped result components, including 624 original-validator scalar cells. The run also freshly regenerated all seven CSVs from frozen RDS with identical SHA256 and repeated 3,466/940 reference comparisons. These are checks on mapped downstream results, not all source/input tables. Sanitized verification summaries are included; raw/private execution and audit logs are excluded. Scripts and files are hashed in SHA256SUMS.tsv and PUBLIC_FILE_INVENTORY.tsv.

Private execution logs, prior missing-workbook records and authorization evidence are preserved outside the publication tree. Public summaries retain the numerical scope and environment; they do not distribute the full audits.


The selected release contains four CSVs: weight_grid.csv, genetic_sd.csv, partial_rg.csv and alm_restrictions.csv. composite.csv, bone_global.csv and bone_pairwise.csv are absent and EXCLUDED_PENDING_CLARIFICATION. Prior full numeric-audit counts above include these now-excluded results; they do not imply that all seven study tables are public. The prior mapped component counts for the four selected tables total 546.

The synthetic tests and runner require no excluded study CSVs. The runner writes all families from supplied synthetic or authorized matrices; synthetic outputs are not study data. validate_frozen_references.py is deliberately unchanged: it reads the separately recomputed full private output directory and private reference tables, not this four-file public folder. A full V5 workbook mapping also requires the private full result set. Do not point these full-scope private validation commands at the reduced public data directory and expect full coverage.

Retained originals use the file-specific MIT/CC BY 4.0 terms in LICENSE_SCOPE.md. Current selected release status is READY_FOR_R2_RELEASE. The excluded full FN/LS-derived tables remain outside the grant and release selection pending source clarification. Complete raw/frozen inputs, audits and workbook/maps stay private.
