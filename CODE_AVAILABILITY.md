# Code availability

All existing R2 analysis code is retained unchanged. code/r2/ contains recovered numerical functions plus clearly identified new adapters/checks; tests/test_r2_smoke.py exercises them with synthetic matrices. SOURCE_FUNCTION_MANIFEST.tsv records source paths, spans and exact function hashes. R_ADAPTER_PROVENANCE.tsv records the independent R adaptation.

Original R2 code, the synthetic test script and new inventory-check portions of code/smoke_test.py are scoped MIT under the author's grant. Historical file-specific licenses remain unchanged. Dependencies are separately installed, not vendored or relicensed. No exclusive ownership or approval on behalf of other authors is asserted.

Removing the three study CSVs does not change the mathematical implementation. The synthetic tests read no excluded study tables, and the runner produces all output families from supplied matrices. Optional full frozen-reference validation needs the privately generated full result set and authorized references; it is not a public-only smoke test. The explicit Excel validator needs its private workbook and mapping.

The reduced public result set is READY_FOR_R2_RELEASE. The FN/LS source condition remains EXCLUDED_PENDING_CLARIFICATION for the omitted CSVs and any proposed later redistribution of study FN/LS output. Retaining generic bone-analysis code is not a grant to the underlying inputs/results. r2/README.md gives synthetic and private-input commands; no upstream GWAS, LDSC or SEM fitting is required by those commands.
