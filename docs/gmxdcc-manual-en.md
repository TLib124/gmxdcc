# gmxdcc Reference Manual

**Package:** gmxdcc 0.1.0 — GARCH-MIDAS and DCC-MIDAS Models with Exogenous
Covariates and Structural Breaks

This document plays the role of a CRAN reference manual: it explains, function
by function, the design rationale and meaning of every argument, with the full
mathematical model and runnable examples (all based on the bundled datasets).

中文版见 [gmxdcc-manual.md](gmxdcc-manual.md).

---

## Contents

1. [Model framework and the two-step workflow](#1-model-framework-and-the-two-step-workflow)
2. [Mathematical model](#2-mathematical-model)
3. [`fit_garch_midas()` — univariate GARCH-MIDAS](#3-fit_garch_midas)
4. [`fit_dcc()` — two-step DCC / DCC-MIDAS](#4-fit_dcc)
5. [Optimizer design: multistart, constraints, warm restarts](#5-optimizer-design)
6. [Out-of-sample design](#6-out-of-sample-design)
7. [Helper functions](#7-helper-functions)
8. [Bundled datasets](#8-bundled-datasets)
9. [Migration table from the legacy research script](#9-migration-table)

---

## 1. Model framework and the two-step workflow

```
Daily return panel r_t (N assets)
        │
Step 1  │ One univariate volatility model per asset
        │   rugarch family: sGARCH / gjrGARCH / eGARCH / iGARCH / csGARCH
        │   or GARCH-MIDAS (this package: arbitrary low-frequency
        │   covariates + structural breaks)
        ▼
Standardized residuals ε_t = (r_t − μ) / σ_t,  diagonal SD cubes D_t
        │
Step 2  │ A correlation model on ε_t
        │   cDCC / aDCC / DECO (C++ kernels)
        │   or DCC-MIDAS (realized-correlation term RC + arbitrary
        │   covariates + structural breaks)
        ▼
Conditional correlations R_t,  covariances H_t = D_t R_t D_t
```

Both steps are constrained QMLE and both report Bollerslev–Wooldridge sandwich
standard errors. Second-step standard errors do not correct for first-stage
estimation uncertainty, in line with standard practice in the two-step DCC
literature (Engle 2002).

**Design philosophy.** The legacy research script enumerated model variants
through name suffixes (`GM_noskew_RV_X_sb` and friends — a 6×2×2 grid of
copy-pasted branches). The package replaces this with **orthogonal switches**:
covariates come as a list, breaks as a date vector, distribution and lag
weighting as single arguments. Any combination works automatically, and every
legacy model is a special case (see the mapping table in Section 9).

---

## 2. Mathematical model

### 2.1 Univariate GARCH-MIDAS

Daily returns decompose into a slow and a fast volatility component:

$$r_t = \mu + \sqrt{\tau_t\, g_t}\;\epsilon_t,\qquad \epsilon_t \sim \mathcal{N}(0,1)\ \text{or}\ t_\nu$$

**Short-run component** (unit-mean GARCH(1,1) recursion):

$$g_t = (1-\alpha-\beta) + \alpha\,\frac{(r_{t-1}-\mu)^2}{\tau_{t-1}} + \beta\, g_{t-1},\qquad g_1 = 1$$

**Long-run component** (MIDAS filter over an arbitrary number $J$ of
low-frequency covariates):

$$\log \tau_t = m(t) + \sum_{j=1}^{J} \theta_j(t) \sum_{k=1}^{K_j} \varphi_k(\omega_j)\, X_{j,\,t-k}$$

where $X_{j,t-k}$ is the value of covariate $j$ in the $k$-th low-frequency
period (month/quarter) preceding the period containing day $t$, and $K_j$ is a
covariate-specific lag order.

**MIDAS weights** $\varphi_k(\omega)$, numerically identical to rumidas:

- Beta (restricted, $w_1=1$, monotonically decaying):
  $$\varphi_k(\omega) = \frac{(1-k/K)^{\omega-1}}{\sum_{i=1}^{K}(1-i/K)^{\omega-1}},\qquad \omega > 1$$
- Exponential Almon ($w_1=0$):
  $$\varphi_k(\omega) = \frac{\exp(\omega k^2)}{\sum_{i=1}^{K}\exp(\omega i^2)},\qquad \omega < 0$$

**Structural breaks.** Given break dates $d_1 < \dots < d_B$, build the step
dummy matrix $D\in\mathbb{R}^{(1+B)\times T}$ (row 1 all ones; row $b{+}1$
equals 1 for $t \ge d_b$). Then

$$m(t) = \mathbf{m}^\top D_{\cdot t},\qquad \theta_j(t) = \boldsymbol{\theta}_j^\top D_{\cdot t}$$

so the intercept and **all** slopes are piecewise constant. With $B=0$ this
collapses to constants — the no-break model is a special case of the break
model, which is exactly how one code path replaces the twelve legacy variants.

### 2.2 Correlation models

Step 1 delivers standardized residuals $\varepsilon_t\in\mathbb{R}^N$.

**cDCC** (Aielli's corrected DCC):

$$Q_t = (1-a-b)\,S + a\, Q_{t-1}^{*1/2}\varepsilon_{t-1}\varepsilon_{t-1}^\top Q_{t-1}^{*1/2} + b\, Q_{t-1},\qquad R_t = Q_t^{*-1/2} Q_t\, Q_t^{*-1/2}$$

with $S$ the unconditional residual covariance (the anchor) and
$Q_t^* = \mathrm{diag}(Q_t)$.

**aDCC** (asymmetric term):

$$Q_t = (1-a-b)S - \gamma\bar N + a\,\varepsilon_{t-1}\varepsilon_{t-1}^\top + b\,Q_{t-1} + \gamma\,\eta_{t-1}\eta_{t-1}^\top,\qquad \eta_t = \min(\varepsilon_t, 0)$$

**DECO** (equicorrelation): compute the cDCC $R_t$, take the cross-sectional
average correlation $\rho_t$, and set $R_t^{deco} = (1-\rho_t)I + \rho_t J$.

**DCC-MIDAS**: the short-run $Q$ recursion is as above, but the anchor becomes
the slowly-moving long-run correlation $\bar R_t$:

$$Q_t = (1-a-b)\,\bar R_t + a\,\varepsilon_{t-1}\varepsilon_{t-1}^\top + b\, Q_{t-1}$$

The long-run correlation takes one of three functional forms (inherited from
the legacy design):

1. **Pure RC, no breaks** (`rc=TRUE, corr_X=NULL, corr_break_dates=NULL`):
   $$\bar R_{ij,t} = \sum_{k=1}^{K_c} \varphi_k(\omega_{RC})\, C_{ij,\,t-k}$$
   where $C_t$ is the rolling realized correlation matrix over a window of
   $N_c$ days (no tanh transform).
2. **With covariates, no breaks** — pair-specific coefficients through a
   Fisher-type transform guaranteeing $|\bar R_{ij}|<1$:
   $$\bar R_{ij,t} = \tanh\!\Big(\delta_{ij} + \textstyle\sum_j \theta_{ij}^{(j)}\,\mathrm{MIDAS}(X_j;\omega_{ij}^{(j)})_t \;[+\; \theta^{RC}_{ij}\, RC_{ij,t}]\Big)$$
3. **With breaks** — levels/slopes shared across pairs but regime-varying
   (driven by the $D$ matrix), while the $\omega$ weights stay pair-specific:
   $$\bar R_{ij,t} = \tanh\!\Big(m(t) + \theta_{RC}(t)\,RC_{ij,t} + \textstyle\sum_j \theta_j(t)\,\mathrm{MIDAS}(X_j;\omega_{ij}^{(j)})_t\Big)$$

### 2.3 Two-step QMLE and robust standard errors

Each step feeds a per-observation log-likelihood vector to `maxLik`
(inequality-constrained BFGS). The robust variance is the sandwich estimator

$$\widehat{\mathrm{Var}}(\hat\vartheta) = H^{-1} \Big(\textstyle\sum_t s_t s_t^\top\Big) H^{-1}$$

with $H$ the numerical Hessian and $s_t$ the per-observation scores (from
`maxLik`'s `gradientObs`). Innovation distributions in GARCH-type models are
almost surely misspecified, hence plain Hessian standard errors are not used
(Bollerslev & Wooldridge, 1992).

---

## 3. `fit_garch_midas()`

```r
fit_garch_midas(returns, X = NULL, K = 12, rv = is.null(X), K_rv = NULL,
                period = c("monthly", "quarterly"),
                distribution = c("norm", "std"),
                lag_fun = c("Beta", "Almon"),
                break_dates = NULL,
                out_of_sample = NULL, vol_proxy = NULL,
                n_starts = 100, seed = NULL, control = list())
```

### Arguments, one by one

| Argument | Type/default | Meaning and design rationale |
|---|---|---|
| `returns` | xts, one column | Daily returns. xts is enforced because the package needs the calendar for MIDAS alignment |
| `X` | NULL / xts / named list | Low-frequency covariates of the long-run component — **any number**. A list (rather than fixed slots) makes "three macro variables" and "one" the same interface; list names become parameter labels (e.g. `theta_JPMU`) |
| `K` | 12 | MIDAS lag order. A scalar is recycled over covariates; a vector sets a covariate-specific $K_j$ |
| `rv` | `is.null(X)` | Prepend internally generated realized volatility (from `realized_measure()`) as the first covariate. Default logic: without `X` there must be at least RV, otherwise the long-run component has no driver |
| `K_rv` | `K[1]` | Separate lag order for the RV covariate |
| `period` | "monthly" | Low frequency (monthly/quarterly); governs RV aggregation and the lag-matrix calendar |
| `distribution` | "norm" | Innovation distribution; `"std"` adds the degrees-of-freedom parameter $\nu$ (constraint $\nu>2.001$) |
| `lag_fun` | "Beta" | Weighting scheme; determines the start value and constraint direction of $\omega$ (Beta: $\omega>1.001$; Almon: $\omega<-0.001$) |
| `break_dates` | NULL | Structural break dates. **Replaces the legacy `_sb` name suffixes**; NULL = no breaks, $B$ dates split $m$ and all $\theta_j$ into $B{+}1$ regimes |
| `out_of_sample` | NULL | Hold out the last n days from estimation; components are filtered over the full sample (warm start, Section 6) |
| `vol_proxy` | NULL | Variance proxy (e.g. realized variance) for MSE/QLIKE evaluation; defaults to squared returns |
| `n_starts` | 100 | Number of random optimizer starts (legacy `R`) |
| `seed` | NULL | RNG seed for the start draw. **The legacy script had no seed and was irreproducible**; this is a new reproducibility guarantee |
| `control` | list() | Optimizer overrides: `bounds` (override the constants table), `start_matrix` (start injection, Section 5), `iterlim`, `method` |

### Canonical parameter ordering (`gm_par_index()`)

All internal functions, `coef()` and `start_matrix` follow one ordering:

$$\underbrace{\alpha,\ \beta}_{\text{short-run}},\ \underbrace{\omega_1..\omega_J}_{\text{weights}},\ \underbrace{\nu}_{\text{std only}},\ \underbrace{m_1..m_{B+1}}_{\text{regime intercepts}},\ \underbrace{\theta_{1,1}..\theta_{J,B+1}}_{\text{regime slopes}},\ \mu$$

Naming: `alpha, beta, w2_<covariate>, nu, m_1.., theta_<covariate>_1.., mu`.

### Return value (class `gm_fit`)

| Field | Content |
|---|---|
| `rob_coef_mat` | Estimate / robust SE / t / p table (rownames = parameter names) |
| `loglik`, `inf_criteria`, `obs`, `period` | In-sample likelihood, AIC/BIC, sample size, span |
| `est_vol_in_s`, `est_lr_in_s`, `est_sr_in_s` | Total conditional vol $\sqrt{h_t}$, long-run $\tau_t$, short-run $g_t$ |
| `loss_in_s` (and `_oos` fields) | MSE(%), QLIKE |
| `maxlik` | Raw maxLik object (used by `vcov()` etc.) |

S3 methods: `summary`, `print`, `coef`, `logLik`, `vcov`.

### Example

```r
library(gmxdcc)

fit <- fit_garch_midas(
  topix_market,
  X = list(JPMU = jp_macro_uncertainty),
  K = 24, rv = FALSE,
  break_dates = jp_mu_breaks,
  distribution = "norm", lag_fun = "Beta", period = "monthly",
  seed = 1
)
summary(fit)          # coefficient table + significance + losses
coef(fit)             # named vector in canonical order
sqrt(vcov(fit))       # sandwich covariance matrix
```

---

## 4. `fit_dcc()`

```r
fit_dcc(returns,
        univariate = c("sGARCH","gjrGARCH","eGARCH","iGARCH","csGARCH","garch_midas"),
        univariate_args = list(),
        correlation = c("cDCC", "aDCC", "DECO", "DCC-MIDAS"),
        rc = TRUE, N_c = NULL, corr_X = NULL, K_c = NULL,
        corr_break_dates = NULL,
        distribution = c("norm", "std"), lag_fun = c("Beta", "Almon"),
        period = c("monthly", "quarterly"),
        out_of_sample = NULL, vol_proxy = NULL,
        n_starts = 100, seed = NULL, control = list())
```

### Arguments, one by one

| Argument | Type/default | Meaning and design rationale |
|---|---|---|
| `returns` | xts, ≥ 2 columns | Multivariate daily return panel |
| `univariate` | "sGARCH" | First-step model. The rugarch family goes through `rugarch::ugarchfit`; `"garch_midas"` fits this package's model per asset |
| `univariate_args` | list() | **First-step-only** argument bundle: `X`, `K`, `rv`, `K_rv`, `break_dates` as in `fit_garch_midas()`; optionally `start_matrix` (per-asset start injection: one matrix used for all assets, or a list with one matrix per asset) |
| `correlation` | "cDCC" | Second-step model |
| `rc` | TRUE | DCC-MIDAS only: include the realized-correlation term in the long-run correlation? |
| `N_c` | NULL | Rolling window (days) of the realized correlation; required when `rc = TRUE` |
| `corr_X` | NULL | Low-frequency covariates of the long-run correlation (list, same design as `X`) |
| `K_c` | NULL | Correlation-side MIDAS lag order (scalar); required for DCC-MIDAS. **Keep it NULL for cDCC/aDCC/DECO** — the C++ kernels apply their own legacy default offset to NULL, and passing a value shifts the recursion start |
| `corr_break_dates` | NULL | Break dates for the long-run correlation (piecewise $m,\theta$; form 3 in 2.2) |
| others | — | Same meaning as in `fit_garch_midas()`; `control$start_matrix` targets **the second step only** (columns in `dcc_par_index()` order) |

### Second-step parameter ordering (`dcc_par_index()`)

With $P = N(N-1)/2$ asset pairs and $R = 1+B$ regimes:

- **No breaks, covariates present:** `a, b, (w2_RC,) delta_1..P, theta_<X>_1..P, w2_<X>_1..P, (theta_RC_1..P)`
- **With breaks:** `a, b, (w2_RC,) w2_<X>_1..P, m_1..R, (theta_RC_1..R,) theta_<X>_1..R`
- **Pure RC, no breaks:** `a, b, w2_RC` only

### Return value (class `dcc_fit`)

| Field | Content |
|---|---|
| `est_univ_model` | Per-asset first-step coefficient tables (list) |
| `corr_coef_mat` | Second-step coefficient table (robust SEs) |
| `eps_t`, `D_t` | Standardized residuals, conditional-SD cubes |
| `H_t`, `R_t`, `R_t_bar` | Conditional covariance / correlation / long-run correlation arrays ($N\times N\times T$) |
| `*_oos` fields | Out-of-sample slices (warm-started filter) |
| `u_llk`, `c_llk`, `Inf_criteria_2step` | Step likelihoods and system AIC/BIC |

### Example

```r
panel <- xts::merge.xts(topix_autos, topix_market)
panel <- panel[complete.cases(panel), ]

dcc <- fit_dcc(
  panel,
  univariate = "garch_midas",
  univariate_args = list(X = jp_macro_uncertainty, K = 24, rv = FALSE,
                         break_dates = jp_mu_breaks),
  correlation = "DCC-MIDAS", rc = FALSE,
  corr_X = list(JPMU = jp_macro_uncertainty), K_c = 24,
  corr_break_dates = jp_mu_breaks,
  seed = 1
)
summary(dcc)
# conditional correlation path
plot(xts::xts(dcc$R_t[1, 2, ], order.by = dcc$Days))
```

---

## 5. Optimizer design

### 5.1 Multistart mechanism

The likelihood surface is multimodal (especially with breaks). Following the
legacy scheme: draw `n_starts` candidate start vectors → evaluate the summed
log-likelihood at each → run **one** constrained BFGS from the best candidate.
Random parameters are drawn uniformly from the ranges below; fixed parameters
(the $\omega$ weights, the mean $\mu$) take fixed start values.

### 5.2 Constants tables (overridable via `control$bounds`)

**Univariate step** (`gm_bounds_default`):

| Item | Start | Constraint |
|---|---|---|
| $\alpha$ | U(0.001, 0.095) | $>0.001$ |
| $\beta$ | U(0.6, 0.8) | $>0.001$; $\alpha+\beta<0.999$ |
| $m,\theta$ (per regime) | U(−1, 1) | none |
| $\omega$ (Beta / Almon) | 1.01 / −0.1 | $>1.001$ / $<-0.001$ |
| $\nu$ | U(2.01, 10) | $>2.001$ |
| $\mu$ | sample mean (fixed) | none |

**Correlation step** (`dcc_bounds_default`; note the $a$ lower bound differs
from step one, inherited from the legacy script):

| Item | Start | Constraint |
|---|---|---|
| $a$ | U(0.001, 0.095) (multistart) or 0.01 (single-start models) | $>0.0001$ |
| $b$ | U(0.6, 0.8) or 0.8 | $>0.001$; $a+b<0.999$ |
| $\gamma$ (aDCC) | 0.01 | $0.001<\gamma<0.15$ |
| $\delta,\theta,m$ | U(−1, 1) | none |
| $\omega$ (Beta / Almon) | 1.2 / −0.05 | $>1.001$ / $<-0.001$ |

cDCC / aDCC / DECO / pure-RC are single-fixed-start models (as in the legacy
script; `n_starts` has no effect on them).

### 5.3 Start injection and warm restarts (key to rolling re-estimation)

- `control$start_matrix` injects **second-step** starts (rows = candidates,
  columns in `dcc_par_index()` order);
- `univariate_args$start_matrix` injects **first-step** starts (per asset,
  columns in `gm_par_index()` order);
- Injected starts that sit on a constraint boundary (common after the
  6-decimal rounding of `coef()`) are automatically blended minimally into
  the strict interior — `constrOptim`'s log-barrier requires a strictly
  interior start and would otherwise fail with "initial value not finite".

Typical rolling re-estimation pattern:

```r
prev_u <- NULL; prev_c <- NULL
for (w in windows) {
  fit <- fit_dcc(data_w, univariate = "garch_midas",
                 univariate_args = list(..., start_matrix = prev_u),
                 control = list(start_matrix = prev_c),
                 n_starts = if (is.null(prev_c)) 300 else 1, ...)
  prev_u <- lapply(fit$est_univ_model, function(m) matrix(m[,"Estimate"], nrow=1))
  prev_c <- matrix(coef(fit), nrow = 1)
}
```

---

## 6. Out-of-sample design

With `out_of_sample = n`:

1. Parameters are estimated on the first $T-n$ days only (clean information
   set);
2. All conditional quantities ($g_t,\tau_t,Q_t,R_t,H_t$) come from
   **filtering the full sample with frozen parameters** and slicing — the
   out-of-sample recursion continues seamlessly from the last in-sample state
   (**warm start**), equivalent to a sequence of one-step-ahead forecasts
   under the *fixed estimation scheme*;
3. Unconditional anchors ($S$, aDCC's $\bar N$) are computed from in-sample
   residuals only (the `S_in`/`N_in` arguments of the C++ kernels) — **no
   look-ahead**.

The legacy script cold-restarted the out-of-sample segment ($g\!\to\!1$,
$Q\!\to\!I$) and computed $S$ from out-of-sample residuals; the new design is
the one consistent with the forecasting literature. See NEWS.md for the full
record of this behavior change.

---

## 7. Helper functions

### `realized_measure(x, period)`
Aggregates daily returns into low-frequency realized variance
$RV_\tau = \sum_{t\in\tau} r_t^2$, indexed at the first day of each period.
```r
rv <- realized_measure(topix_market, "monthly")
```

### `align_sample(r_t, K, K_c, macro, period)`
Trims the daily series so that every retained day has $\max(K,K_c)+1$ complete
low-frequency periods of covariate history (a prerequisite of the MIDAS
filter).

### `midas_lag_matrix(x, mv, K, period)`
Stacks, for every trading day, the current and $K$ previous low-frequency
values of a covariate into a $(K{+}1)\times T$ lag matrix — the direct input
of the long-run builder.

### `midas_weights(k, K, w1, w2, lag_fun)`
The weight formulas of Section 2.1 (numerically identical to rumidas, pinned
by a unit test).
```r
midas_weights(1:12, 12, w1 = 1, w2 = 3, lag_fun = "Beta")
```

### `build_break_dummies(dates, break_dates)`
The step-dummy matrix $D$ of Section 2.1; returns a single row of ones for
`NULL` breaks.

### `loss_avg(vol_est, vol_proxy)`
$$\mathrm{MSE\%} = 100\cdot\overline{(h - h^{proxy})^2},\qquad \mathrm{QLIKE} = \overline{\log h + h^{proxy}/h}$$

### `qmle_se(est)` / `inf_criteria(est)` / `inf_criteria_2step(res_garch, res_dcc)`
Sandwich standard errors (formula in 2.3), single-model AIC/BIC, and two-step
system AIC/BIC ($k$ and the likelihood summed over both steps).

---

## 8. Bundled datasets

| Dataset | Content | Frequency |
|---|---|---|
| `topix_market` | TOPIX market returns (column `TOPIX`) | daily |
| `topix_autos` | TOPIX-17 Autos & Transportation Equipment returns (column `AUTOS`) | daily |
| `jp_macro_uncertainty` | Japan macro uncertainty index (column `JPMU`, standardized) | monthly |
| `jp_mu_breaks` | Four structural break dates of JPMU (`strucchange::breakpoints`) | — |

All are clean Date-indexed xts objects, available right after
`library(gmxdcc)`.

---

## 9. Migration table

| Legacy model name | New call |
|---|---|
| `GM_noskew_RV` | `rv = TRUE, X = NULL` |
| `GM_noskew_X` | `rv = FALSE, X = list(x)` |
| `GM_noskew_RV_X` | `rv = TRUE, X = list(x)` |
| `GM_noskew_*_sb` | add `break_dates = <dates>` |
| `cDCC` / `aDCC` / `DECO` | `correlation = <same name>` |
| `DCCMIDAS_RC` | `correlation = "DCC-MIDAS", rc = TRUE` |
| `DCCMIDAS_X` | `correlation = "DCC-MIDAS", rc = FALSE, corr_X = list(x)` |
| `DCCMIDAS_RC_X` | `correlation = "DCC-MIDAS", rc = TRUE, corr_X = list(x)` |
| `DCCMIDAS_*_sb` | add `corr_break_dates` |

Argument renames: `X_stru_break_date_vec → break_dates / corr_break_dates`;
`type → period`; `R → n_starts`.

**Behavioral differences** (deliberate improvements, see NEWS.md): warm-started
instead of cold-restarted out-of-sample filtering; unconditional anchors from
in-sample residuals only; seedable multistart. Formula-level equivalence with
the legacy script is guaranteed by the test suite (per-observation agreement
at fixed parameters to 1e-12).

---

*This manual tracks package version 0.1.0; for the authoritative signatures
see the help pages `?fit_garch_midas` and `?fit_dcc`.*
