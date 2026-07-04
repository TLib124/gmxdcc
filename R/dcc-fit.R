#' First (univariate) estimation stage
#'
#' Fits one volatility model per asset and returns everything the
#' correlation stage needs: full-sample conditional standard deviations
#' (in-sample fit, out-of-sample continued by each model's own filter),
#' demeaned returns, per-asset coefficient tables and log-likelihoods.
#' 第一步:逐资产估计波动率模型,输出全样本条件标准差(样本外由各自
#' 滤波器延续)、去均值收益及各资产结果。
#'
#' @param r_al Aligned multivariate xts of daily returns.
#' @param univariate Model name (rugarch family or `"garch_midas"`).
#' @param distribution `"norm"` or `"std"`.
#' @param out_of_sample Number of held-out observations (or `NULL`).
#' @param univariate_args List of arguments for the `garch_midas` path:
#'   `X`, `K`, `rv`, `K_rv`, `break_dates` (same semantics as
#'   [fit_garch_midas()]); `rv` covariates are computed from `r_full`.
#'   Optionally `start_matrix`: optimizer starts for the first step — one
#'   matrix used for every asset, or a list with one matrix per asset.
#' @param r_full Original (untrimmed) returns, used to compute realized
#'   volatility with full history.
#' @param period,lag_fun,n_starts,seed,control Passed to [fit_garch_midas()].
#' @param vol_proxy Optional multivariate variance proxy.
#'
#' @return List with `u_est`, `est_details`, `u_loglik`, `u_inf`, `sd_full`
#'   (TT x N matrix), `db_demeaned` (TT x N matrix).
#' @keywords internal
univ_stage <- function(r_al, univariate, distribution, out_of_sample,
                       univariate_args, r_full, period, lag_fun,
                       n_starts, seed, control, vol_proxy = NULL) {

  N_assets <- ncol(r_al)
  TT   <- nrow(r_al)
  n_in <- if (is.null(out_of_sample)) TT else TT - out_of_sample

  u_est <- est_details <- u_loglik <- u_inf <- vector("list", N_assets)
  sd_full <- matrix(NA_real_, TT, N_assets)
  db_dem  <- zoo::coredata(r_al)

  if (univariate == "garch_midas") {

    ua <- univariate_args
    rv_flag <- ua$rv %||% is.null(ua$X)
    K  <- ua$K %||% 12

    # Optional per-asset optimizer starts (e.g. rolling re-estimation warm
    # starts): a single matrix used for every asset, or a list with one
    # matrix per asset, columns in gm_par_index() order.
    # 可选的逐资产优化起点(如滚动重估热启动):单个矩阵作用于全部资产,
    # 或按资产给列表;列序遵循 gm_par_index()。
    sm <- ua$start_matrix
    if (!is.null(sm) && is.list(sm) && length(sm) != N_assets) {
      stop("univariate_args$start_matrix must be a single matrix or a list with one matrix per asset.",
           call. = FALSE)
    }

    for (i in seq_len(N_assets)) {
      # Realized volatility from the FULL column history so that the aligned
      # sample keeps K+1 periods of RV lags (passing it as an explicit
      # covariate makes fit_garch_midas' internal alignment a no-op).
      # RV 用原始整列历史计算并作为显式协变量传入,保证对齐样本的滞后
      # 历史充足,内部再对齐不再裁剪。
      X_i <- ua$X
      if (!is.null(X_i) && inherits(X_i, "xts")) X_i <- list(X = X_i)
      K_i <- rep_len(as.integer(K), length(X_i) %||% 0L)
      if (rv_flag) {
        X_i <- c(list(RV = realized_measure(r_full[, i], period)), X_i)
        K_i <- c(as.integer(ua$K_rv %||% K[1]), K_i)
      }

      ctrl_i <- control
      ctrl_i$start_matrix <- if (is.list(sm)) sm[[i]] else sm

      u_est[[i]] <- fit_garch_midas(
        r_al[, i], X = X_i, K = K_i, rv = FALSE,
        period = period, distribution = distribution, lag_fun = lag_fun,
        break_dates = ua$break_dates,
        out_of_sample = out_of_sample,
        vol_proxy = if (!is.null(vol_proxy)) vol_proxy[, i] else NULL,
        n_starts = n_starts, seed = seed, control = ctrl_i
      )

      est_details[[i]] <- u_est[[i]]$rob_coef_mat
      u_loglik[[i]]    <- u_est[[i]]$loglik
      u_inf[[i]]       <- u_est[[i]]$inf_criteria
      sd_full[, i]     <- c(u_est[[i]]$est_vol_in_s, u_est[[i]]$est_vol_oos)
      # Demean by the estimated mu (legacy behavior for GM models)
      # 按估计的 mu 去均值(与旧脚本 GM 路径一致)
      db_dem[, i] <- db_dem[, i] - u_est[[i]]$rob_coef_mat["mu", "Estimate"]
    }

  } else {

    # rugarch family / rugarch 传统 GARCH 族
    uspec <- rugarch::ugarchspec(
      variance.model = list(model = univariate, garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
      distribution.model = distribution
    )

    for (i in seq_len(N_assets)) {
      u_est[[i]] <- if (is.null(out_of_sample)) {
        rugarch::ugarchfit(spec = uspec, data = r_al[, i])
      } else {
        rugarch::ugarchfit(spec = uspec, data = r_al[, i], out.sample = out_of_sample)
      }
      est_details[[i]] <- u_est[[i]]@fit$robust.matcoef
      u_loglik[[i]]    <- u_est[[i]]@fit$LLH
      u_inf[[i]]       <- inf_criteria(u_est[[i]])
      sd_full[seq_len(n_in), i] <- u_est[[i]]@fit$sigma
      if (!is.null(out_of_sample)) {
        # Rolling 1-step-ahead sigma for the held-out days (legacy behavior;
        # rugarch's rolling forecast is warm-started by construction).
        # 样本外用滚动一步预测的 sigma(rugarch 自身即为热启动)。
        sd_full[(n_in + 1):TT, i] <- as.numeric(rugarch::sigma(
          rugarch::ugarchforecast(u_est[[i]], n.ahead = 1,
                                  n.roll = c(out_of_sample - 1),
                                  out.sample = c(out_of_sample - 1))
        ))
      }
      # Demean by the unconditional mean (legacy behavior for rugarch path)
      # 按无条件均值去均值(与旧脚本 rugarch 路径一致)
      db_dem[, i] <- db_dem[, i] - mean(zoo::coredata(r_al)[, i], na.rm = TRUE)
    }
  }

  list(u_est = u_est, est_details = est_details, u_loglik = u_loglik,
       u_inf = u_inf, sd_full = sd_full, db_demeaned = db_dem)
}

#' Fit a DCC / DCC-MIDAS model with flexible covariates and breaks
#'
#' Two-step QMLE: univariate volatility models per asset (rugarch family or
#' GARCH-MIDAS), then a correlation model on the standardized residuals.
#' The DCC-MIDAS long-run correlation admits a realized-correlation term
#' (`rc`), any number of exogenous covariates (`corr_X`), and structural
#' breaks (`corr_break_dates`); robust Bollerslev–Wooldridge sandwich
#' standard errors are reported in both steps (without first-stage
#' correction in the second step, as is standard for two-step DCC).
#'
#' Old model names map as: `cDCC`/`aDCC`/`DECO` unchanged;
#' `DCCMIDAS_RC` = `correlation = "DCC-MIDAS", rc = TRUE`;
#' `DCCMIDAS_X` = `rc = FALSE, corr_X = list(x)`;
#' `DCCMIDAS_RC_X` = both; `_sb` suffix = `corr_break_dates`.
#'
#' With `out_of_sample`, all conditional matrices for the held-out days are
#' obtained by filtering the full sample with frozen parameters and
#' in-sample unconditional targets — warm-started, no look-ahead (the
#' legacy script restarted cold and used out-of-sample covariance targets).
#' 样本外矩阵由冻结参数全样本滤波获得,无冷启动、无前视。
#'
#' @param returns Multivariate xts of daily returns (>= 2 columns).
#' @param univariate Univariate model: `"sGARCH"`, `"gjrGARCH"`, `"eGARCH"`,
#'   `"iGARCH"`, `"csGARCH"` or `"garch_midas"`.
#' @param univariate_args List for the `garch_midas` path: `X`, `K`, `rv`,
#'   `K_rv`, `break_dates` (see [fit_garch_midas()]).
#' @param correlation Correlation model: `"cDCC"`, `"aDCC"`, `"DECO"` or
#'   `"DCC-MIDAS"`.
#' @param rc DCC-MIDAS: include the realized-correlation MIDAS term?
#' @param N_c Rolling-window length of the realized correlation (required
#'   when `rc = TRUE`).
#' @param corr_X `NULL`, an xts, or a named list of low-frequency covariates
#'   for the long-run correlation.
#' @param K_c MIDAS lag order of the correlation covariates (scalar).
#' @param corr_break_dates Optional break dates for the long-run correlation.
#' @param distribution,lag_fun,period See [fit_garch_midas()].
#' @param out_of_sample Number of trailing observations to hold out.
#' @param vol_proxy Optional multivariate variance proxy (same columns as
#'   `returns`).
#' @param n_starts,seed,control Optimizer settings shared by both steps.
#'   Exception: `control$start_matrix` targets the correlation step only
#'   (columns in `dcc_par_index()` order); first-step starts go through
#'   `univariate_args$start_matrix`. 注:control$start_matrix 只作用于
#'   第二步,第一步起点经 univariate_args$start_matrix 注入。
#'
#' @return An object of class `dcc_fit`; see [summary.dcc_fit()].
#' @examples
#' \donttest{
#' # Industry vs market conditional correlation, long-run correlation
#' # driven by macro uncertainty with breaks (bundled sample data)
#' panel <- xts::merge.xts(topix_autos, topix_market)
#' panel <- panel[stats::complete.cases(panel), ]["2005/2019"]
#'
#' dcc <- fit_dcc(
#'   panel, univariate = "sGARCH",
#'   correlation = "DCC-MIDAS", rc = FALSE,
#'   corr_X = list(JPMU = jp_macro_uncertainty), K_c = 12,
#'   corr_break_dates = jp_mu_breaks[2:3],
#'   n_starts = 3, seed = 1
#' )
#' summary(dcc)
#' }
#' @export
fit_dcc <- function(returns,
                    univariate = c("sGARCH", "gjrGARCH", "eGARCH", "iGARCH",
                                   "csGARCH", "garch_midas"),
                    univariate_args = list(),
                    correlation = c("cDCC", "aDCC", "DECO", "DCC-MIDAS"),
                    rc = TRUE, N_c = NULL, corr_X = NULL, K_c = NULL,
                    corr_break_dates = NULL,
                    distribution = c("norm", "std"),
                    lag_fun = c("Beta", "Almon"),
                    period = c("monthly", "quarterly"),
                    out_of_sample = NULL, vol_proxy = NULL,
                    n_starts = 100, seed = NULL, control = list()) {

  univariate   <- match.arg(univariate)
  correlation  <- match.arg(correlation)
  distribution <- match.arg(distribution)
  lag_fun      <- match.arg(lag_fun)
  period       <- match.arg(period)

  ############################# validation / 输入校验
  assert_xts(returns)
  if (ncol(returns) < 2) stop("returns must contain at least two assets.", call. = FALSE)
  if (nrow(returns) < 50) stop("returns must be long enough (>= 50 observations).", call. = FALSE)
  assert_xts(vol_proxy, allow_null = TRUE)
  if (!is.null(vol_proxy) && ncol(vol_proxy) != ncol(returns)) {
    stop("vol_proxy must have the same number of columns as returns.", call. = FALSE)
  }
  is_midas_corr <- correlation == "DCC-MIDAS"
  if (is_midas_corr && rc && is.null(N_c)) {
    stop("N_c (realized-correlation window) is required when rc = TRUE.", call. = FALSE)
  }
  if (!is.null(corr_X) && inherits(corr_X, "xts")) corr_X <- list(X = corr_X)
  if (!is.null(corr_X)) {
    for (x_i in corr_X) assert_xts(x_i, arg = "corr_X")
    if (is.null(names(corr_X)) || any(names(corr_X) == "")) {
      names(corr_X) <- paste0("X", seq_along(corr_X))
    }
  }
  if (is_midas_corr && !rc && length(corr_X) == 0L) {
    stop("DCC-MIDAS needs rc = TRUE and/or corr_X covariates.", call. = FALSE)
  }
  if (is_midas_corr && is.null(K_c)) {
    stop("K_c (correlation MIDAS lag order) is required for DCC-MIDAS.", call. = FALSE)
  }
  # For cDCC/aDCC/DECO, K_c stays NULL: the C++ kernels then use their own
  # legacy default (start offset 2). Passing 2 explicitly would shift it.
  # 非 MIDAS 模型 K_c 保持 NULL,由 C++ 内部取旧默认;显式传值会偏移。

  ############################# calendar alignment / 日历对齐
  # Trim the panel so every retained day has full covariate history for the
  # univariate stage (RV / X of asset 1: identical calendars across assets)
  # and the correlation stage (each corr_X with K_c lags).
  # 裁剪面板:单变量与相关两步的所有协变量都有完整滞后历史。
  r_al <- returns
  if (univariate == "garch_midas") {
    ua <- univariate_args
    rv_flag <- ua$rv %||% is.null(ua$X)
    K_u <- ua$K %||% 12
    if (rv_flag) {
      r_al <- align_sample(r_al, K = (ua$K_rv %||% K_u[1]),
                           macro = realized_measure(returns[, 1], period),
                           period = period)
    }
    X_u <- ua$X
    if (!is.null(X_u) && inherits(X_u, "xts")) X_u <- list(X_u)
    if (length(X_u) > 0L) {
      K_ux <- rep_len(as.integer(K_u), length(X_u))
      for (j in seq_along(X_u)) {
        r_al <- align_sample(r_al, K = K_ux[j], macro = X_u[[j]], period = period)
      }
    }
  }
  if (is_midas_corr && length(corr_X) > 0L) {
    for (x_i in corr_X) {
      r_al <- align_sample(r_al, K = K_c, macro = x_i, period = period)
    }
  }

  TT   <- nrow(r_al)
  N_a  <- ncol(r_al)
  if (!is.null(out_of_sample) && (out_of_sample <= 0 || out_of_sample >= TT)) {
    stop("out_of_sample must be strictly between 0 and the aligned sample size.",
         call. = FALSE)
  }
  n_in <- if (is.null(out_of_sample)) TT else TT - out_of_sample
  days      <- zoo::index(r_al)
  full_days <- as.Date(days)
  if (!is.null(vol_proxy)) vol_proxy <- vol_proxy[days]

  ############################# first step / 第一步:单变量
  # control$start_matrix targets the SECOND (correlation) step only: its
  # columns follow dcc_par_index() and would not match the univariate
  # specs. Per-asset starts for the first step are injected via
  # univariate_args$start_matrix instead.
  # control$start_matrix 只作用于第二步(相关性);第一步的起点按资产
  # 经 univariate_args$start_matrix 注入,两步参数空间互不混淆。
  control_univ <- control
  control_univ$start_matrix <- NULL

  st1 <- univ_stage(r_al, univariate, distribution, out_of_sample,
                    univariate_args, r_full = returns, period = period,
                    lag_fun = lag_fun, n_starts = n_starts, seed = seed,
                    control = control_univ, vol_proxy = vol_proxy)

  # Standardized residuals and D_t cubes over the full sample
  # 全样本标准化残差与条件标准差对角阵
  eps_full <- st1$db_demeaned / st1$sd_full
  D_full   <- array(0, dim = c(N_a, N_a, TT))
  for (tt in seq_len(TT)) diag(D_full[, , tt]) <- st1$sd_full[tt, ]

  eps_in <- eps_full[seq_len(n_in), , drop = FALSE]

  ############################# second step / 第二步:相关性
  # Correlation covariate lag matrices on the full calendar
  # 相关协变量滞后矩阵(全样本日历)
  X_mats_full <- lapply(corr_X %||% list(), function(x_i) {
    midas_lag_matrix(r_al[, 1], x_i, K = K_c, period = period)
  })
  D_mat_full <- build_break_dummies(full_days, corr_break_dates)

  dcc_bundle <- function(cols, idx = NULL) {
    list(res = eps_full[cols, , drop = FALSE], n_assets = N_a,
         TT = length(cols), K_c = K_c, N_c = N_c,
         rc = isTRUE(rc) && is_midas_corr,
         X_mats = lapply(X_mats_full, function(m) m[, cols, drop = FALSE]),
         D_mat = D_mat_full[, cols, drop = FALSE],
         lag_fun = lag_fun, idx = idx)
  }

  b2 <- utils::modifyList(dcc_bounds_default, control$bounds %||% list())

  if (!is_midas_corr) {
    # cDCC / aDCC / DECO: single fixed start, exactly as the legacy script.
    # 单一固定初值,与旧脚本一致。
    dc_in <- dcc_bundle(seq_len(n_in))
    unit3 <- function(i, n) { v <- numeric(n); v[i] <- 1; v }
    if (correlation == "aDCC") {
      start_val <- c(a = b2$a_fixed, b = b2$b_fixed, g = b2$g_fixed)
      ui <- rbind(unit3(1, 3), unit3(2, 3), unit3(3, 3),
                  -unit3(1, 3) - unit3(2, 3), -unit3(3, 3))
      ci <- c(-b2$a_min, -b2$b_min, -b2$g_min, b2$ab_sum_max, b2$g_max)
      ll_fn <- adcc_loglik_fn
    } else {
      start_val <- c(a = b2$a_fixed, b = b2$b_fixed)
      ui <- rbind(unit3(1, 2), unit3(2, 2), -unit3(1, 2) - unit3(2, 2))
      ci <- c(-b2$a_min, -b2$b_min, b2$ab_sum_max)
      ll_fn <- if (correlation == "cDCC") cdcc_loglik_fn else deco_loglik_fn
    }
    m_est <- suppressWarnings(maxLik::maxLik(
      logLik = ll_fn, start = start_val, data = dc_in,
      constraints = list(ineqA = ui, ineqB = ci),
      iterlim = control$iterlim %||% 1000, method = control$method %||% "BFGS"
    ))
  } else {
    spec <- dcc_param_spec(N_a, rc = isTRUE(rc), n_X = length(X_mats_full),
                           n_regimes = nrow(D_mat_full), lag_fun = lag_fun,
                           X_names = names(corr_X) %||% character(0),
                           bounds = control$bounds %||% list())
    dc_in <- dcc_bundle(seq_len(n_in), idx = spec$idx)
    opt <- run_multistart(dcc_midas_loglik, dc_in, spec,
                          n_starts = if (spec$single_start) 1 else n_starts,
                          seed = seed, control = control)
    m_est <- opt$est
  }

  est_coef <- stats::coef(m_est)
  mat_coef <- qmle_coeftable(m_est)

  ############################# matrices: full-sample filter, then slice
  ############################# 全样本滤波后切片(样本外热启动)
  S_in <- stats::cov(eps_in)
  N_in <- if (correlation == "aDCC") stats::cov(pmin(eps_in, 0)) else NULL

  dc_full <- dcc_bundle(seq_len(TT), idx = if (is_midas_corr) dc_in$idx else NULL)
  mats <- dcc_matrices(est_coef, correlation, dc_full, D_full,
                       S_in = S_in, N_in = N_in)

  in_ix <- seq_len(n_in)
  slice3 <- function(A, ix) if (is.null(A)) NULL else A[, , ix, drop = FALSE]

  mult_label <- if (is_midas_corr) {
    sprintf("DCC-MIDAS(%s%s%s)",
            if (isTRUE(rc)) "RC" else "",
            if (length(corr_X)) paste0(if (isTRUE(rc)) "+" else "",
                                       paste(names(corr_X), collapse = "+")) else "",
            if (nrow(D_mat_full) > 1) sprintf(", %d breaks", nrow(D_mat_full) - 1L) else "")
  } else correlation

  fin_res <- list(
    assets = colnames(r_al),
    model = univariate,
    est_univ_model = st1$est_details,
    est_univ_model_origin = st1$u_est,
    corr_coef_mat = mat_coef,
    mult_model = mult_label,
    obs = n_in,
    period = range(days[in_ix]),
    eps_t = eps_in,
    D_t = slice3(D_full, in_ix),
    H_t = slice3(mats$H_t, in_ix),
    R_t = slice3(mats$R_t, in_ix),
    R_t_bar = slice3(mats$R_t_bar, in_ix),
    Days = days[in_ix],
    u_llk = st1$u_loglik,
    c_llk = stats::logLik(m_est),
    u_Inf_criteria = st1$u_inf,
    c_inf_criteria = inf_criteria(m_est),
    Inf_criteria_2step = inf_criteria_2step(st1$u_est, m_est),
    spec = list(correlation = correlation, rc = isTRUE(rc) && is_midas_corr,
                corr_X = names(corr_X), K_c = K_c, N_c = N_c,
                corr_break_dates = corr_break_dates,
                distribution = distribution, lag_fun = lag_fun,
                period = period, n_starts = n_starts, seed = seed),
    maxlik = m_est
  )

  if (!is.null(out_of_sample)) {
    oos_ix <- (n_in + 1L):TT
    fin_res$eps_t_oos   <- eps_full[oos_ix, , drop = FALSE]
    fin_res$D_t_oos     <- slice3(D_full, oos_ix)
    fin_res$H_t_oos     <- slice3(mats$H_t, oos_ix)
    fin_res$R_t_oos     <- slice3(mats$R_t, oos_ix)
    fin_res$R_t_bar_oos <- slice3(mats$R_t_bar, oos_ix)
    fin_res$Days_oos    <- days[oos_ix]
  }

  class(fin_res) <- "dcc_fit"
  fin_res
}
