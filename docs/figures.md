# Figures and tables

The numerical figure and table inputs are in `data/publication/`. Their origin,
selection and SHA-256 hashes are recorded in `provenance/data_sources.json`.
`provenance/figure_table_registry.json` maps 17 figure filenames and table groups
to sources and inputs. Historical filenames differ from current figure numbers;
use `manuscript_figure_number` for the current order.

The original plotting algorithms, colours, interval viewports, circular windows
and approved labels are retained. This release exports existing arrays without
new simulations or figure regeneration. Exact PDF equivalence of the portable
adapters has not been verified. Rendered figures, manuscripts and fonts are excluded.

## Inputs

| File | Use |
| --- | --- |
| `coverage.csv` | Euler and plate M=30000 primary coverage; nonlinear P=M=50000 with the final 500-million-response reference |
| `euler_pool_control.csv` | Euler P=M=30000 comparison, separate from fresh independent primary trials |
| `coverage_summary.csv` | Means, prediction RMSE and deviations, including labelled Euler control |
| `reference_linear.csv`, `reference_plate.csv` | Production density ordinates at 201 locations |
| `reference_nonlinear_*.csv` | Final nonlinear density, coverage and internal-n reference grids |
| `internal_n.csv` | n=192/384/768 at R=128, B=399, P=M=30000 |
| `intervals_linear_nonlinear.csv`, `intervals_plate.csv` | Saved independently generated intervals at recommended R |
| `recommended_R.csv` | Current method-specific MC and Sobol recommendations |
| `saved_trial_intervals.csv`, `table_interval_values.csv` | 450 original first-trial rows and 126 rows selected for the three tables |
| `finite_B_summary*.csv`, `coverage_coefficients.csv` | Finite-B baseline simulation and saved calibration coefficients |
| `gf_summary.csv`, `gf_by_dimension.csv`, `initial_design_*.csv` | Weighted GF discrepancies and original base coordinates and weights |
| `outer_workflow_curves.csv` | Headerless 256-by-201 pilot slice for the outer workflow illustration |

Coverage shading is MC uncertainty in an empirical coverage frequency. Density
confidence intervals concern a density ordinate. Original method strings are
retained, including `bootstrap-t`; figure labels use manuscript terminology.
Empty fields indicate quantities not applicable to an experiment, not zero.

## Rendering

From the repository root, choose a new output directory for each invocation:

```text
python plotting/historicalproduction.py --kind python --data-dir data/publication --output-dir outputs/python_figures
python plotting/historicalproduction.py --kind tables --data-dir data/publication --output-dir outputs/tables
```

The first command invokes five Python numerical figure routines. The second
formats the three endpoint tables without interval estimation. For six MATLAB
numerical figures, with MATLAB and the required local fonts installed:

```matlab
addpath('plotting/matlab');
replay_matlab('data/publication', 'outputs/matlab_figures');
```

Static modes are `transport`, `plate`, `euler`, `shift`, `framework` and `outer`.
See `python plotting/historicalproduction.py --help`. The outer workflow needs
original Times font paths; the later label overlay is supplied separately as
`plotting/static/outer_overlay.tex`. Check its final export chain before claiming
exact reproduction. Some static sources require local TeX, ReportLab,
Matplotlib, NumPy and original fonts. Importing does not run simulations.

## Checks

`scripts/check_release.py` checks hashes, syntax, approved data inventory and
sensitive-path patterns. `scripts/check_publication_data.py` checks schemas,
counts, summary arithmetic, reference consistency and original GF values.
These checks do not rerun bootstrap trials or finite element solutions.
`plotting/prepare_publication.py` is for maintainers with the original archive;
it is unnecessary for reading these CSVs. Optional PDF copying is off by default.
