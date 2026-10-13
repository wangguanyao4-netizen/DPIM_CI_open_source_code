# Python Numerical API

Copyright (c) 2026 Guanyao Wang. The supplied FE/KL programs and their Python
ports are released under the MIT license. External libraries retain their
own licenses and are not included in this source package.

## Scope and Dependencies

This code-only release contains the numerical Python source used in the
current manuscript calculation. Research data, coefficient tables and figure
files are not distributed. The FE, random-field and coefficient Python
modules are byte-identical copies of their source
snapshots below. The bridge changes only package imports and moves its native
thread limit from import time to the actual FE/KL calls. It does not mutate
environment variables or `sys.path` when imported, and does not create or write
files on import. The namespace `__init__.py` is packaging scaffolding only.

Install NumPy, SciPy, SymPy and threadpoolctl externally. The locally observed
versions are NumPy 2.4.4, SciPy 1.17.1, SymPy 1.14.0 and threadpoolctl 3.6.0;
these are an environment inventory, not a historical experiment lock or a new
compatibility claim. No libraries, PyMuPDF, fonts or proprietary FE software are
bundled. NumPy/SciPy serve the FE/KL and periodic nearest-neighbour code; SymPy
is used by the finite-B coefficient module. MATLAB supplies its own Python
integration. Configure Python to find the release's `python` directory.

## API Contract

The public bridge and coefficient interfaces are:

```python
from dpim_paper import bridge
from dpim_paper import coverage_coefficients

bridge.assigned_probabilities(base, count, seed)  # ndarray, shape (n,)
bridge.initialize_plate()                       # None; lazily caches MODEL/MODES/AUDIT
bridge.solve_plate(theta, validate=False)        # ndarray, shape (n_samples, 1)
bridge.save_plate_model(path)                   # explicit NPZ write; returns str(AUDIT)

coverage_coefficients.finite_b_design(alpha, B)  # FiniteBDesign
coverage_coefficients.direct_coefficients(alpha, B)  # (rows, checks)
coverage_coefficients.coefficient_lookup(rows)  # {(method, term): value}
coverage_coefficients.predict_coverage(method, R, alpha, B, lambda3_sq, lambda4)
```

`count` is the number of auxiliary uniform points; positional calls are
compatible with the historical bridge. `base` is an `n x d` unit-cube design.
`theta` has shape `n_samples x 10`. The default plate uses a 20 by 20
unit-square mesh, 800 KL support points, ten modes and seed 2026092207.
`MODEL`, `MODES` and `AUDIT` are initialized
only on demand. Native FE/KL operations through the bridge are limited to two
threads and restore the previous native thread setting on exit. Direct
low-level calls retain the original behavior and do not impose a process-wide
limit themselves.

`solve_plate` checks element-modulus positivity and the `1e-8` residual
threshold, with optional direct/iterative comparison. It does not clip the
material distribution. The low-level solver uses a direct solve as fallback
when the iterative solve fails its checks.
`save_plate_model` writes to an explicit path using NumPy's overwrite behavior;
callers must supply an appropriate new output path. `test_periodic_assignment`
and `benchmark` are optional diagnostic functions, not invoked automatically.

Coefficient method strings are `t distribution`, `percentile bootstrap` and
`bootstrap-t`. The raw-percentile `A0`, `A4`, `A33` terms and quadrature rules
are retained. Both `alpha` and `B` are explicit original inputs. Coverage
predictions are asymptotic expressions, not guaranteed probabilities; the
module does not clip them to `[0,1]`. `direct_coefficients(alpha, B)` generates
coefficients from the included symbolic expressions and numerical quadrature
code when explicitly called. The coefficient module does not read CSV files
or require an external coefficient-data path.

## Source IDs and SHA256

Source IDs below are relative to the originating research project; they are
provenance identifiers, not runtime paths or external dependencies.

| Release file (relative to `python/dpim_paper/`) | Project source ID | Source SHA256 |
| --- | --- | --- |
| `bridge.py` | `tc_recalculation_20260922/tc_bridge.py` | `edc89196eb62a1cdc8f89974851e3e211fe27463548d898c97a6258133c0c5dd` |
| `dpim_plate/fem.py` | `tc_recalculation_20260922/vendor/dpim_plate/fem.py` | `29b8b8206eadaa90ae301e99609a223cf7546df69b41d206cb742d2afe05cec5` |
| `dpim_plate/random_field.py` | `tc_recalculation_20260922/vendor/dpim_plate/random_field.py` | `56fd79c9a8ec8588b096f3e8b169cfc4134876c80d666bc62d5b255516de0155` |
| `dpim_plate/__init__.py` | `tc_recalculation_20260922/vendor/dpim_plate/__init__.py` | `8b2a9a35bec136632a2c954d95a3471a918edc2b66077bb0912f4b078dc9eb59` |
| `coverage_coefficients.py` | `tc_recalculation_20260922/vendor/coverage_coefficients.py` | `66ef5a9aa8b5841265186701d02d4ca2c5f61e8c6b75d48e397963f9d6829ca8` |

