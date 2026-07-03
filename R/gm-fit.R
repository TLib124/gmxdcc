#' Prepare the data bundle for a GARCH-MIDAS fit
#'
#' Normalizes the covariate set (optional internally generated realized
#' volatility plus any number of exogenous low-frequency series), aligns the
#' daily calendar so every retained day has full MIDAS history for every
#' covariate, and builds the lag matrices and break dummies.
#' 归一化协变量集合(可选内生 RV + 任意个外生低频序列),对齐日历,
#' 构建滞后矩阵与突变虚拟矩阵。
#'
#' @inheritParams fit_garch_midas
#' @return List with `returns` (trimmed xts), `lag_mats`, `K` (per-covariate
#'   vector), `cov_names`, `D_mat`, `full_days`, and the trimmed `vol_proxy`.
#' @keywords internal
gm_prepare_data <- function(returns, X = NULL, K = 12, rv = is.null(X),
                            K_rv = NULL, period = "monthly",
                            break_dates = NULL, vol_proxy = NULL) {
  assert_xts(returns)
  assert_xts(vol_proxy, allow_null = TRUE)

  # --- normalize the covariate list --- 归一化协变量列表
  if (!is.null(X) && inherits(X, "xts")) X <- list(X)
  if (!is.null(X)) {
    for (x_i in X) assert_xts(x_i, arg = "X")
    if (is.null(names(X)) || any(names(X) == "")) {
      names(X) <- paste0("X", seq_along(X))
    }
  }
  if (!rv && length(X) == 0L) {
    stop("At least one long-run covariate is required: set rv = TRUE and/or supply X.",
         call. = FALSE)
  }

  covs <- list()
  K_vec <- integer(0)
  if (rv) {
    # Internally generated realized volatility as first covariate
    # 内生 RV 作为第一个协变量,滞后阶数 K_rv(默认取 K 的首元素)
    covs[["RV"]] <- realized_measure(returns, period)
    K_vec <- c(K_vec, as.integer(K_rv %||% K[1]))
  }
  if (length(X) > 0L) {
    K_x <- rep_len(as.integer(K), length(X))   # recycle K over X / K 循环补齐
    for (j in seq_along(X)) {
      covs[[names(X)[j]]] <- X[[j]]
      K_vec <- c(K_vec, K_x[j])
    }
  }

  # --- calendar alignment: sequential trimming intersects the admissible
  # --- days of all covariates. 顺序裁剪等价于各协变量可用日历的交集。
  r_req <- returns
  for (j in seq_along(covs)) {
    r_req <- align_sample(r_req, K = K_vec[j], macro = covs[[j]], period = period)
  }

  # --- lag matrices on the final calendar --- 最终日历上的滞后矩阵
  lag_mats <- lapply(seq_along(covs), function(j) {
    midas_lag_matrix(r_req, covs[[j]], K = K_vec[j], period = period)
  })
  names(lag_mats) <- names(covs)

  # --- vol_proxy coverage check, as in the legacy wrapper ---
  # vol_proxy 需覆盖全部保留日,与旧封装一致
  if (!is.null(vol_proxy)) {
    ref_dates <- zoo::index(r_req)
    if (!all(ref_dates %in% zoo::index(vol_proxy))) {
      stop("vol_proxy data do not cover the aligned return sample.", call. = FALSE)
    }
    vol_proxy <- vol_proxy[ref_dates]
  }

  full_days <- as.Date(zoo::index(r_req))

  list(
    returns   = r_req,
    lag_mats  = lag_mats,
    K         = K_vec,
    cov_names = names(covs),
    D_mat     = build_break_dummies(full_days, break_dates),
    full_days = full_days,
    vol_proxy = vol_proxy
  )
}

#' Fit a GARCH-MIDAS model with flexible covariates and structural breaks
#'
#' Two-component volatility model \eqn{r_t = \mu + \sqrt{\tau_t g_t}\,\epsilon_t}
#' where the short-run component follows a unit-mean GARCH(1,1) recursion and
#' the long-run component is
#' \deqn{\log \tau_t = m(t) + \sum_j \theta_j(t)\,\mathrm{MIDAS}(X_j)_t.}
#' Any number of low-frequency covariates is allowed, each with its own lag
#' order; `break_dates` makes the intercept `m` and all slopes `theta_j`
#' piecewise-constant. Estimation is constrained QMLE with random
#' multi-start; robust Bollerslev–Wooldridge sandwich standard errors are
#' reported (see [qmle_se()]).
#'
#' Old model names map as: `GM_noskew_RV` = `rv = TRUE`; `GM_noskew_X` =
#' `X = list(x)`; `GM_noskew_RV_X` = both; the `_sb` suffix = `break_dates`.
#'
#' When `out_of_sample` is set, parameters are estimated on the first
#' `N - out_of_sample` days and the variance components for the held-out days
#' are obtained by **filtering over the full sample with frozen parameters**,
#' so the out-of-sample recursion is warm-started from the last in-sample
#' state. (The legacy script restarted the recursion cold at the first
#' out-of-sample day.)
#' 样本外部分用冻结参数在全样本滤波后切片,从样本内末期状态热启动,
#' 不再冷启动。
#'
#' @param returns Daily return series (`xts`, one column).
#' @param X `NULL`, an `xts`, or a named list of low-frequency `xts`
#'   covariates entering the long-run component.
#' @param K MIDAS lag order; scalar (recycled over `X`) or vector.
#' @param rv If `TRUE`, prepend internally generated realized volatility
#'   (from [realized_measure()]) as a covariate. Default: `TRUE` when `X` is
#'   `NULL`.
#' @param K_rv Lag order of the realized-volatility covariate (default
#'   `K[1]`).
#' @param period Frequency of the long-run component: `"monthly"` or
#'   `"quarterly"`.
#' @param distribution Innovation distribution, `"norm"` or `"std"`.
#' @param lag_fun MIDAS weighting scheme, `"Beta"` or `"Almon"`.
#' @param break_dates Optional `Date` vector of structural-break dates
#'   segmenting `m` and all `theta_j` (`NULL` = no breaks).
#' @param out_of_sample Optional number of trailing observations to hold out.
#' @param vol_proxy Optional `xts` variance proxy for MSE/QLIKE loss
#'   evaluation (default: squared returns).
#' @param n_starts Number of random optimizer starts (legacy default 100).
#' @param seed Optional RNG seed for reproducible starts.
#' @param control List of optimizer overrides: `bounds`, `start_matrix`,
#'   `iterlim`, `method`.
#'
#' @return An object of class `gm_fit`; see [summary.gm_fit()].
#' @examples
#' \donttest{
#' # Market volatility driven by macro uncertainty, with structural breaks
#' # (bundled sample data; small n_starts to keep the example quick)
#' fit <- fit_garch_midas(
#'   topix_market["2005/2019"],
#'   X = list(JPMU = jp_macro_uncertainty), K = 12, rv = FALSE,
#'   break_dates = jp_mu_breaks[2:3],
#'   distribution = "norm", lag_fun = "Beta", period = "monthly",
#'   n_starts = 3, seed = 1
#' )
#' summary(fit)
#' coef(fit)
#' }
#' @export
fit_garch_midas <- function(returns, X = NULL, K = 12, rv = is.null(X),
                            K_rv = NULL,
                            period = c("monthly", "quarterly"),
                            distribution = c("norm", "std"),
                            lag_fun = c("Beta", "Almon"),
                            break_dates = NULL,
                            out_of_sample = NULL, vol_proxy = NULL,
                            n_starts = 100, seed = NULL, control = list()) {
  period       <- match.arg(period)
  distribution <- match.arg(distribution)
  lag_fun      <- match.arg(lag_fun)

  prep <- gm_prepare_data(returns, X = X, K = K, rv = rv, K_rv = K_rv,
                          period = period, break_dates = break_dates,
                          vol_proxy = vol_proxy)

  N <- length(prep$returns)
  if (!is.null(out_of_sample) && (out_of_sample <= 0 || out_of_sample >= N)) {
    stop("out_of_sample must be strictly between 0 and the aligned sample size.",
         call. = FALSE)
  }
  n_in <- if (is.null(out_of_sample)) N else N - out_of_sample

  # --- data bundles: in-sample for estimation, full for filtering ---
  # 估计用样本内数据束;成分滤波用全样本数据束(OOS 热启动的实现)。
  slice_bundle <- function(cols) {
    list(
      daily_ret = prep$returns[cols],
      lag_mats  = lapply(prep$lag_mats, function(m) m[, cols, drop = FALSE]),
      K         = prep$K,
      D_mat     = prep$D_mat[, cols, drop = FALSE],
      lag_fun   = lag_fun,
      distribution = distribution,
      idx       = gm_par_index(length(prep$lag_mats), nrow(prep$D_mat),
                               distribution, prep$cov_names)
    )
  }
  data_in   <- slice_bundle(seq_len(n_in))
  data_full <- slice_bundle(seq_len(N))

  # --- estimation --- 估计
  spec <- gm_param_spec(
    n_cov = length(prep$lag_mats), n_regimes = nrow(prep$D_mat),
    distribution = distribution, lag_fun = lag_fun,
    cov_names = prep$cov_names,
    mu_start = mean(data_in$daily_ret, na.rm = TRUE),
    bounds = control$bounds %||% list()
  )
  opt <- run_multistart(gm_loglik, data_in, spec,
                        n_starts = n_starts, seed = seed, control = control)
  est      <- opt$est
  est_coef <- stats::coef(est)
  mat_coef <- qmle_coeftable(est)

  # --- components: full-sample filter at frozen parameters, then slice ---
  # 冻结参数全样本滤波,再按样本内/外切片。
  comp   <- gm_components(est_coef, data_full)
  in_ix  <- seq_len(n_in)

  vol_proxy_all <- if (is.null(prep$vol_proxy)) prep$returns^2 else prep$vol_proxy
  vol_est_lf    <- zoo::coredata(comp$h_it_2[in_ix])
  vol_proxy_lf  <- zoo::coredata(vol_proxy_all[in_ix])

  res <- list(
    model = sprintf(
      "GARCH-MIDAS(%s%s, %s, %s)",
      paste(prep$cov_names, collapse = "+"),
      if (nrow(prep$D_mat) > 1) sprintf(", %d breaks", nrow(prep$D_mat) - 1L) else "",
      distribution, lag_fun
    ),
    rob_coef_mat = mat_coef,
    obs          = n_in,
    period       = range(stats::time(prep$returns[in_ix])),
    loglik       = as.numeric(stats::logLik(est)),
    inf_criteria = inf_criteria(est),
    loss_in_s    = loss_avg(vol_est_lf, vol_proxy_lf),
    est_vol_in_s = vol_est_lf^0.5,
    est_lr_in_s  = zoo::coredata(comp$tau_d[in_ix]),
    est_sr_in_s  = zoo::coredata(comp$g_it[in_ix]),
    # metadata / 元数据
    spec = list(cov_names = prep$cov_names, K = prep$K, period = period,
                distribution = distribution, lag_fun = lag_fun,
                break_dates = break_dates, n_starts = n_starts, seed = seed),
    maxlik = est
  )

  if (!is.null(out_of_sample)) {
    oos_ix <- (n_in + 1L):N
    vol_est_oos_lf   <- zoo::coredata(comp$h_it_2[oos_ix])
    vol_proxy_oos_lf <- zoo::coredata(vol_proxy_all[oos_ix])
    res$loss_oos    <- loss_avg(vol_est_oos_lf, vol_proxy_oos_lf)
    res$est_vol_oos <- vol_est_oos_lf^0.5
    res$est_lr_oos  <- zoo::coredata(comp$tau_d[oos_ix])
    res$est_sr_oos  <- zoo::coredata(comp$g_it[oos_ix])
  }

  class(res) <- "gm_fit"
  res
}
