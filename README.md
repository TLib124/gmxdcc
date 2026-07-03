# gmxdcc

GARCH-MIDAS and DCC-MIDAS models with flexible exogenous covariates and
structural breaks, in R with C++ (RcppArmadillo) likelihood kernels.

GARCH-MIDAS 与 DCC-MIDAS 模型的 R 实现:任意个数低频宏观协变量、
结构突变开关、C++ 加速的似然核心。

## Features

- **Univariate GARCH-MIDAS** (`fit_garch_midas()`): the long-run variance
  component admits *any number* of low-frequency covariates — internally
  generated realized volatility and/or user-supplied macro series — each
  with its own MIDAS lag order (Beta or exponential-Almon weights,
  normal or Student-t innovations).
- **Correlation models** (`fit_dcc()`): cDCC, aDCC, DECO and DCC-MIDAS,
  the latter with a realized-correlation term, exogenous covariates and
  structural breaks in the long-run correlation.
- **Structural breaks** as an argument, not a model name: pass
  `break_dates` / `corr_break_dates` and the long-run intercepts and slopes
  become piecewise-constant. `NULL` = no breaks.
- **Robust inference**: Bollerslev–Wooldridge (1992) QMLE sandwich standard
  errors in both estimation steps. (Second-step standard errors do not
  correct for first-stage estimation uncertainty, which is standard
  practice in the two-step DCC literature — see Engle (2002).)
- **Out-of-sample by filtering**: with `out_of_sample`, held-out
  conditional variances/correlations are produced by filtering the full
  sample with frozen in-sample parameters — warm-started from the last
  in-sample state, with in-sample unconditional targets (no look-ahead).
- **Reproducibility**: the random multi-start optimizer accepts a `seed`.

## Installation

```r
# install.packages("devtools")
devtools::install_github("TLib124/gmxdcc")
```

## Quick start (with the bundled sample data)

The package ships four small datasets so every example below runs as-is:
`topix_market` (daily TOPIX market returns), `topix_autos` (daily TOPIX-17
Autos & Transportation Equipment industry returns), `jp_macro_uncertainty`
(monthly Japan macro uncertainty index, column `JPMU`) and `jp_mu_breaks`
(its structural break dates). 包内置四个小型示例数据集,以下示例可直接运行。

```r
library(gmxdcc)
library(xts)

## --- Univariate GARCH-MIDAS: market volatility driven by macro uncertainty,
## --- with structural breaks in the long-run component
fit <- fit_garch_midas(
  topix_market,
  X  = list(JPMU = jp_macro_uncertainty),  # any number of covariates
  K  = 24,                                 # MIDAS lag order
  rv = FALSE,                              # set TRUE to add realized vol too
  break_dates = jp_mu_breaks,              # piecewise m and theta
  distribution = "norm", lag_fun = "Beta", period = "monthly",
  seed = 1
)
summary(fit)

## --- Two-step DCC-MIDAS: industry vs market conditional correlation
panel <- merge.xts(topix_autos, topix_market)
panel <- panel[complete.cases(panel), ]

dcc <- fit_dcc(
  panel,
  univariate = "garch_midas",
  univariate_args = list(X = jp_macro_uncertainty, K = 24, rv = FALSE,
                         break_dates = jp_mu_breaks),
  correlation = "DCC-MIDAS", rc = FALSE,
  corr_X = list(JPMU = jp_macro_uncertainty), K_c = 24,
  corr_break_dates = jp_mu_breaks,
  distribution = "norm", lag_fun = "Beta", period = "monthly",
  seed = 1
)
summary(dcc)

# Conditional correlation path / 条件相关系数路径
plot(xts(dcc$R_t[1, 2, ], order.by = dcc$Days))
```

## Migrating from the legacy research script

| Legacy model string | New arguments |
|---|---|
| `GM_noskew_RV` | `rv = TRUE, X = NULL` |
| `GM_noskew_X` | `rv = FALSE, X = list(x)` |
| `GM_noskew_RV_X` | `rv = TRUE, X = list(x)` |
| `GM_noskew_*_sb` | add `break_dates = ...` |
| `cDCC` / `aDCC` / `DECO` | `correlation = "cDCC" / "aDCC" / "DECO"` |
| `DCCMIDAS_RC` | `correlation = "DCC-MIDAS", rc = TRUE` |
| `DCCMIDAS_X` | `correlation = "DCC-MIDAS", rc = FALSE, corr_X = list(x)` |
| `DCCMIDAS_RC_X` | `correlation = "DCC-MIDAS", rc = TRUE, corr_X = list(x)` |
| `DCCMIDAS_*_sb` | add `corr_break_dates = ...` |

Other renames: `X_stru_break_date_vec` → `break_dates` / `corr_break_dates`,
`type` → `period`, `R` → `n_starts`.

Behavior changes relative to the legacy script (see `NEWS.md`):
out-of-sample components are warm-started instead of cold-restarted, the
out-of-sample unconditional targets use in-sample residuals only, and the
multi-start draw is seedable.

## Parallel use

The legacy project kept a second copy of the script with `sourceCpp()`
lines commented out for use on parallel workers. With the package this is
unnecessary: call `library(gmxdcc)` inside each worker and the compiled
kernels are available automatically.

## License

GPL (>= 3)
