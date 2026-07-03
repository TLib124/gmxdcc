#' Information criteria for a fitted model
#'
#' Computes AIC and BIC for a `uGARCHfit` (rugarch), `maxLik`, or `gm_fit`
#' object (see [fit_garch_midas()]). For `gm_fit` objects the criteria
#' stored at fit time are returned directly when available.
#'
#' @param est A fitted object of class `uGARCHfit`, `maxLik`, or `gm_fit`.
#'
#' @return A named numeric vector with elements `AIC` and `BIC`.
#' @examples
#' set.seed(1)
#' y  <- rnorm(300, 1, 2)
#' ll <- function(p) stats::dnorm(y, p[1], exp(p[2]), log = TRUE)
#' fit <- maxLik::maxLik(ll, start = c(mu = 0, logsd = 0), method = "BFGS")
#' inf_criteria(fit)
#' @export
inf_criteria <- function(est) {

  # --- Case 1: rugarch object (e.g. sGARCH, gjrGARCH) ---
  # 情况 1: rugarch 对象,有专门方法提取似然值与数据
  if (inherits(est, "uGARCHfit")) {
    ll <- rugarch::likelihood(est)
    k  <- length(rugarch::coef(est))
    n  <- length(est@model$modeldata$data)

  # --- Case 2: maxLik object (internal optimizer result) ---
  # 情况 2: maxLik 对象,gradientObs 的行数即样本量
  } else if (inherits(est, "maxLik")) {
    ll <- as.numeric(stats::logLik(est))
    k  <- length(stats::coef(est))
    n  <- nrow(est$gradientObs)

  # --- Case 3: gm_fit object (GARCH-MIDAS result list) ---
  # 情况 3: gm_fit 对象,若已存结果直接返回,避免重复计算
  } else if (inherits(est, "gm_fit")) {
    if (!is.null(est$inf_criteria)) {
      return(est$inf_criteria)
    }
    ll <- est$loglik
    k  <- nrow(est$rob_coef_mat)
    n  <- est$obs

  } else {
    stop("inf_criteria: unknown object class provided.", call. = FALSE)
  }

  # --- Unified AIC / BIC --- 统一计算 AIC 与 BIC
  aic <- 2 * k - 2 * ll
  bic <- k * log(n) - 2 * ll

  round(c(AIC = aic, BIC = bic), 6)
}

#' Information criteria for the two-step system
#'
#' Aggregates the univariate first-step fits and the second-step correlation
#' fit into system-wide AIC/BIC. The sample size is taken from the
#' second-step `maxLik` gradient matrix.
#'
#' @param res_garch List of first-step fits (`uGARCHfit` or `gm_fit`).
#' @param res_dcc Second-step `maxLik` object.
#'
#' @return Named numeric vector `AIC`, `BIC`, `k`, `LogLik`.
#' @export
inf_criteria_2step <- function(res_garch, res_dcc) {
  # --- 1. total number of parameters (k_total) ---
  # 兼容 rugarch S4 对象与普通 S3 类
  k_garch <- sum(sapply(res_garch, function(x) {
    if (isS4(x)) {
      return(length(rugarch::coef(x)))
    } else {
      return(length(stats::coef(x)))
    }
  }))

  # res_dcc is a maxLik object / res_dcc 通常是 maxLik 对象
  k_dcc   <- length(stats::coef(res_dcc))
  k_total <- k_garch + k_dcc

  # --- 2. total log-likelihood --- 计算总对数似然
  loglik_garch <- sum(sapply(res_garch, function(x) {
    if (isS4(x)) {
      return(as.numeric(rugarch::likelihood(x)))
    } else {
      # 对自定义类,若 stats::logLik 失败则直接取 $loglik
      ll <- try(stats::logLik(x), silent = TRUE)
      if (inherits(ll, "try-error")) return(as.numeric(x$loglik))
      return(as.numeric(ll))
    }
  }))

  loglik_dcc   <- as.numeric(stats::logLik(res_dcc))
  loglik_total <- loglik_garch + loglik_dcc

  # --- 3. sample size --- gradientObs 仅在收敛后可靠
  n <- nrow(res_dcc$gradientObs)

  # --- 4. criteria ---
  aic <- 2 * k_total - 2 * loglik_total
  bic <- k_total * log(n) - 2 * loglik_total

  round(c(AIC = aic, BIC = bic, k = k_total, LogLik = loglik_total), 6)
}

#' Average volatility loss functions
#'
#' Mean squared error (in percent) and QLIKE loss of a volatility estimate
#' against a volatility proxy.
#'
#' @param vol_est Numeric vector of estimated variances.
#' @param vol_proxy Numeric vector of proxy variances (e.g. squared returns).
#'
#' @return Named numeric vector `MSE(%)` and `QLIKE`.
#' @examples
#' set.seed(1)
#' v_est   <- exp(rnorm(100, 0, 0.3))
#' v_proxy <- exp(rnorm(100, 0, 0.3))
#' loss_avg(v_est, v_proxy)
#' @export
loss_avg <- function(vol_est, vol_proxy) {
  c(
    "MSE(%)" = 100 * mean((vol_est - vol_proxy)^2),
    "QLIKE"  = mean(log(vol_est) + (vol_proxy / vol_est))
  )
}

#' Robust (QMLE sandwich) standard errors
#'
#' Bollerslev–Wooldridge (1992) sandwich standard errors
#' \eqn{H^{-1} (G'G) H^{-1}} computed from the numerical Hessian and the
#' outer product of the per-observation gradients of a `maxLik` fit. If the
#' Hessian is singular, the Moore–Penrose generalized inverse
#' (`MASS::ginv()`) is used instead.
#'
#' Note: in the two-step DCC estimation these standard errors do not account
#' for first-stage parameter uncertainty, which is standard practice in this
#' literature (Engle 2002; Colacito, Engle and Ghysels 2011).
#'
#' @param est A `maxLik` object with a `gradientObs` matrix.
#'
#' @return Numeric vector of robust standard errors (negative-variance
#'   entries caused by numerical error are returned as `NA`).
#' @examples
#' # Tiny normal MLE with per-observation log-likelihood
#' # 小型正态 MLE(逐观测似然,保证 gradientObs 可用)
#' set.seed(1)
#' y  <- rnorm(300, 1, 2)
#' ll <- function(p) stats::dnorm(y, p[1], exp(p[2]), log = TRUE)
#' fit <- maxLik::maxLik(ll, start = c(mu = 0, logsd = 0), method = "BFGS")
#' qmle_se(fit)
#' @export
qmle_se <- function(est) {

  # Hessian of the negative log-likelihood / 取海森矩阵
  H <- -maxLik::hessian(est)

  # Try a standard inverse, fall back to the generalized inverse
  # 标准求逆失败(奇异矩阵)时退回广义逆
  H_hat_inv <- tryCatch({
    solve(H)
  }, error = function(e) {
    message("Warning: Hessian is singular. Switching to generalized inverse (MASS::ginv).")
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Please install 'MASS' package.")
    return(MASS::ginv(H))
  })

  # Outer product of gradients / 逐观测梯度外积
  OPG <- t(est$gradientObs) %*% est$gradientObs

  # Sandwich variance estimator / 三明治方差估计量
  var_cov <- H_hat_inv %*% OPG %*% H_hat_inv

  diag_vals <- diag(var_cov)

  # Numerical error may yield negative variances; report them as NA
  # 数值误差导致的负方差置为 NA
  diag_vals[diag_vals < 0] <- NA

  sqrt(diag_vals)
}

#' Coefficient table with robust standard errors
#'
#' Builds the standard Estimate / Std. Error / t value / p value table from a
#' `maxLik` fit using [qmle_se()] robust standard errors. Used identically by
#' both estimation steps.
#'
#' @param est A `maxLik` object.
#'
#' @return A `data.frame` with columns `Estimate`, `Std. Error`, `t value`,
#'   `Pr(>|t|)`, one row per parameter.
#' @keywords internal
qmle_coeftable <- function(est) {
  est_coef <- stats::coef(est)
  se       <- qmle_se(est)

  mat_coef <- data.frame(
    Estimate     = round(est_coef, 6),
    `Std. Error` = round(se, 6),
    `t value`    = round(est_coef / se, 6),
    `Pr(>|t|)`   = round(2 * (1 - stats::pnorm(abs(est_coef / se))), 6),
    check.names  = FALSE
  )
  rownames(mat_coef) <- names(est_coef)
  mat_coef
}

#' Evaluate an expression, returning NA on error
#'
#' Thin `tryCatch()` wrapper used around optimizer calls so that a failed
#' start value does not abort the whole multi-start search.
#'
#' @param expr Expression to evaluate.
#' @return The value of `expr`, or `NA` if evaluation failed.
#' @keywords internal
safe_maxlik <- function(expr) {
  tryCatch(
    {
      result <- eval(expr)
      return(result)
    },
    error = function(e) {
      return(NA)
    }
  )
}
