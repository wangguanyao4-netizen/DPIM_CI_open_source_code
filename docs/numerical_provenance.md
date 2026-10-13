# Numerical production source

`provenance/numerical_sources.json` maps original source identifiers to the
release records, including original and distributed SHA-256 hashes. Local
absolute paths were replaced by `<configure-local-path>` in non-executable
archive records. Model formulas, budgets, seeds, and finite-B definitions were
not deliberately changed by that packaging operation.

The chronological source stages are included because the published numerical
inputs were built through existing pool expansions and additional coverage
trials. Copying a later driver alone would omit its required preceding inputs.
The new release does not replace these stages with a new simulation algorithm.

The Euler main experiment uses independent outer shifts. Its finite-pool
comparison is retained as a distinct dataset. The plate and nonlinear coverage
experiments sample their respective fixed pools. The nonlinear production
endpoints were recounted using the later 500 million-sample reference; the
earlier reference is not the final coverage target approximation.

The finite-B baseline simulation is preserved as the relevant local functions,
not the unrelated experiment driver that originally contained them. GF values
come from the weighted marginal-CDF supremum implementation in the archived
`calculate_gf.m`; GF result tables are not distributed.

The repository contains source code only, without numerical data or generated
figures. A hash verifies file identity, not a
mathematical conclusion or an independent rerun of a numerical experiment.
