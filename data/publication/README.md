# Manuscript numerical data

These are the archived numerical values supporting the current manuscript
figures and tables. No new simulations were run for this release. CSVs retain
full stored precision; manuscript tables round them for display.

| Example | Main coverage experiment | Reference |
| --- | --- | --- |
| Euler beam | n=768, R=32/128, B=399, M=30000 independent trials | Analytical Gaussian convolution, h=0.08 |
| Kirchhoff plate | n=200, R=16/128, B=999, P=M=30000 fixed-pool trials | 16 independent scrambled Sobol groups of 131072 responses, total 2097152 |
| Nonlinear algebraic model | n=768, R=16/128, B=399, P=M=50000 fixed-pool trials | 500000000 independent MC responses, h=0.08 |

Euler's P=M=30000 pool control is separate. The internal-n comparison uses
P=M=30000 for n=192/384/768. Recommended-R intervals use independent outer
samples; the final nonlinear reference update preserves the archived endpoints.
`rqmc` is the historical code for Sobol initial points with fixed probability
weights and independent common shifts.

See [`../../docs/figures.md`](../../docs/figures.md) for all 23 filenames and roles,
[`../../provenance/figure_table_registry.json`](../../provenance/figure_table_registry.json)
for current figure numbers, and
[`../../provenance/data_sources.json`](../../provenance/data_sources.json) for
hashes, archive-relative sources and selections.

- Coverage: `hit_count/M` is empirical coverage, `formula_prediction` is the
  theoretical approximation, and `cp95_lower/upper` are Clopper-Pearson bounds
  for the frequency. They are not density confidence intervals.
- Intervals: `ci_lower/upper` are density endpoints; `sample_mean`, `sample_sd`
  and `reference_density` retain their original roles.
- First trial: the full file has 450 rows over 25 locations; the table file
  selects seven locations and contains 126 rows.
- References: `y` and `reference_density` are matching response coordinates and
  smoothed density ordinates. `reference_se`, where supplied, concerns the
  numerical reference, not coverage uncertainty.
- GF: the initial design files contain the original unit-cube points and fixed
  weights. The summary takes the maximum marginal discrepancy across dimensions.
- `outer_workflow_curves.csv` is headerless: 256 pilot curves by 201 locations
  matching `reference_linear.csv`.

The full raw curve pools, complete interval-trial batches and individual
reference-response samples are not included. No manuscript, rendered figures,
credentials or workstation paths are distributed. The GF summary's local path
column is replaced by public design filenames; numerical values are unchanged.

The repository MIT license applies. Cite the accompanying manuscript and the
repository version when reusing these files.
