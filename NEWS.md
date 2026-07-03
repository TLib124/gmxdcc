# gmxdcc 0.1.0

First packaged release, refactored from the research script
`integrated noskew GARCH-DCC-MIDAS fit.R` (frozen under `dev/legacy/`).

## New interface

* `fit_garch_midas()` replaces `garchmidas_noskew_fit()` and its wrapper:
  the long-run component now accepts an arbitrary named list of
  low-frequency covariates (`X`), each with its own lag order (`K`), plus
  an optional internally generated realized-volatility covariate (`rv`).
* `fit_dcc()` replaces `dccmidas_noskew_longfocus_fit()` and its wrapper,
  with `correlation = "cDCC" / "aDCC" / "DECO" / "DCC-MIDAS"`, an `rc`
  toggle, a `corr_X` covariate list and `K_c`.
* Structural breaks are arguments (`break_dates`, `corr_break_dates`)
  instead of `_sb` model-name suffixes; intercepts and slopes become
  piecewise-constant across the implied regimes.
* `seed` makes the random multi-start reproducible (the legacy script drew
  starting values without a seed, so fits were not reproducible).
* Optimizer bounds/start ranges are collected in documented tables
  (`gm_bounds_default`, `dcc_bounds_default`) and can be overridden via
  `control$bounds`.

## Behavior changes vs. the legacy script

* **Out-of-sample warm start.** Held-out conditional variances and
  correlations are obtained by filtering the full sample with frozen
  in-sample parameters. The legacy script restarted the recursions cold at
  the first out-of-sample day (short-run variance component reset to 1,
  DCC quasi-correlation reset to the identity) and, for cDCC/aDCC/DECO,
  computed the unconditional target from *out-of-sample* residuals
  (a look-ahead). The C++ kernels gained optional `S_in`/`N_in` arguments
  so all unconditional targets now come from in-sample residuals; the
  `N_c_res_for_oos` patch mechanism is gone.
* In-sample results are numerically identical to the legacy script; the
  regression-test suite verifies likelihoods and components at fixed
  parameters to 1e-12 and end-to-end fits from identical starting values.

## Bug fixes folded into the port

* `integrated_garchmidas_noskew_output_function()` referenced the global
  `univ_model` instead of its own `model` argument (worked only by global
  leakage); the new API has no model strings, removing the bug class.
* `stop(cat(...))` calls (which raised empty error messages) replaced by
  proper `stop()` messages via `assert_xts()` and validators.
* `align_sample()` (né `cut_func`) no longer errors when `period`/`macro`
  are `NULL`.
* MIDAS weights are computed internally (`midas_weights()`), numerically
  identical to `rumidas::beta_function()`/`exp_almon()` (verified by unit
  test); `rumidas` is no longer a hard dependency.
