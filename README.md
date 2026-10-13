# DPIM probabilistic error estimation

Source code and supporting numerical data for:

**The probabilistic error estimation of Direct Probability Integral Method
for stochastic mechanics**, Guanyao Wang, Guohai Chen, Dixiong Yang and Rui Li.

This is the current manuscript code collection. The earlier repository contents
have been replaced in the current checkout; the Git history is retained.

## Contents

| Directory | Purpose |
| --- | --- |
| `matlab/` | Existing MATLAB interval and numerical calculation source |
| `python/` | Existing periodic probability assignment, plate FE/KL, and coverage coefficient source |
| `plotting/` | The manuscript's MATLAB and Python plotting source, model diagrams, and label revisions |
| `config/` | Source manifest and numerical dependency information |
| `source_archive/` | Current-paper production scripts preserved as source records with machine paths removed |
| `provenance/` | Source hashes and mappings to figures, tables, and numerical stages |
| `docs/` | Input requirements and reproduction notes |
| `data/publication/` | Archived numerical data underlying the manuscript figures and tables |

The archived production scripts have `.txt` appended to their original names.
They document the exact calculation stages, including pool expansion and
coverage recounting, without accidentally starting a large calculation.
These records are not an alternative, superseded scientific version.

## Numerical scope

| Example | Main coverage experiment | Reference |
| --- | --- | --- |
| Euler beam | Independent shifts, M = 30000, n = 768, R = 32 and 128, B = 399; a P = 30000 pool provides the comparison and moment estimates | Analytical smoothed normal density |
| Kirchhoff plate | Fixed pool P = 30000, M = 30000, n = 200, R = 16 and 128, B = 999 | 16 independent scrambled Sobol groups, 131072 responses per group, total 2097152 |
| Nonlinear algebraic model | Fixed pool P = 50000, M = 50000, n = 768, R = 16 and 128, B = 399 | 500000000 independent MC responses |

The internal-sample-size comparisons use n = 192, 384, 768 and M = 30000.
The recommended-R density intervals use newly generated independent outer
samples. The endpoint-table scripts use a first archived trial at seven
locations, with two R values and three methods when those inputs are supplied.

The methods are Student t, raw percentile bootstrap, and bootstrap t.
Probability weights are held fixed during independent common random shifts.
The numerical targets, finite-B ranks, random streams, and reference choices
are retained from the manuscript calculations. No centering is introduced.

## Using the code

Read `docs/` and the provenance mappings before choosing an entry point.
MATLAB calculations require a separately installed MATLAB and Statistics and
Machine Learning Toolbox. The original environment used MATLAB R2026a and
Python 3.12. Python dependencies are listed in `requirements.txt`.

```text
python -m pip install -r requirements.txt
```

The publication dataset includes coverage counts and predictions, reference
density ordinates, interval endpoints, recommended sample sizes, initial
designs, GF discrepancies and finite-B coefficients. See
[`data/publication/README.md`](data/publication/README.md) for the mapping and
[`provenance/data_sources.json`](provenance/data_sources.json) for source hashes
and selections. Plotting and table inputs can be read directly from that folder.
The full raw curve pools, all trial histories and individual reference-response
samples are not included. Running the historical simulation stages still
requires their preceding inputs; this is not a one-command rerun of all trials.

See `docs/figures.md` for input formats, label conventions, and figure mapping.
Large calculation budgets are never launched by installing the dependencies.

## Release checks

This release organizes existing calculations and sources; it does not rerun
the numerical experiments. File hashes, Python syntax, source
completeness, and sensitive-path checks are recorded separately from numerical
validation. The source and first-trial interval verification performed while
preparing the manuscript are distinguished from these packaging checks.

No generated PDF figure, manuscript draft, account
credentials, machine-specific logs, bundled
MathWorks installation, OS fonts, or unrelated higher-order study is included.
Source identifiers contain historical filenames only to make provenance
traceable. Local deployment and automatic manuscript-publishing scripts are
not release entry points.

## License and citation

Project source is released under the [MIT License](LICENSE). The author has
confirmed the right to include the FE/KL source and its Python port. External
dependencies retain their own licenses; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Please cite the accompanying manuscript and this repository when using the
method or code. The manuscript is not assigned an invented DOI or publication
status here.
