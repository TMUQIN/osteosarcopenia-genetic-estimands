# Running the selected analysis materials

Python 3.10+ is recommended. Install only packages needed by the chosen script.

From the repository root, run the structural and content validators:

```bash
python code/smoke_test.py
python code/release_validate.py
```

With permission to use the retained source table, Figure 3 can be rebuilt:

```bash
cd code/figures/submitted_figure3
python draw_figure3.py
```

The R verification script in `code/verification/` needs author-controlled S/V/I/N inputs and the original model objects. These are not included. This repository does not reproduce GWAS, LDSC or primary GenomicSEM end to end.
