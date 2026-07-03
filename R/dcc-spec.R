# Second-step optimizer constants, verbatim from the legacy correlation
# branches (overridable via `control$bounds`). Note the second step bounds
# alpha at 0.0001, not 0.001 as in the univariate step.
# 第二步优化常数表,与旧脚本一致;注意 alpha 下界为 0.0001(与单变量步不同)。
# 整理:旧脚本 w2 初值散落 2 / 1.2 / 1.01 / -0.1 / -0.01,此处统一为
# Beta 1.2 / Almon -0.05(仅初值,不影响约束;golden 测试用注入起点)。
dcc_bounds_default <- list(
  a_start     = c(0.001, 0.095),  # multistart draw / 多起点抽取区间
  b_start     = c(0.6, 0.8),
  coef_start  = c(-1, 1),         # delta, theta, m ~ U(-1, 1)
  a_fixed     = 0.01,             # single-start models / 单起点模型固定初值
  b_fixed     = 0.8,
  g_fixed     = 0.01,             # aDCC asymmetry start
  w2_start_beta  = 1.2,
  w2_start_almon = -0.05,
  a_min       = 0.0001,           # a > 0.0001
  b_min       = 0.001,            # b > 0.001
  ab_sum_max  = 0.999,            # a + b < 0.999
  g_min       = 0.001,            # aDCC: g > 0.001
  g_max       = 0.15,             # aDCC: g < 0.15
  w2_min_beta = 1.001,
  w2_max_almon = -0.001
)

#' Parameter index map for the DCC(-MIDAS) correlation step
#'
#' Canonical ordering, chosen to coincide with the legacy layout for every
#' configuration it supported:
#' * no breaks, covariates present: `a, b, (w2_rc,) delta_1..P,
#'   theta_<X>_1..P, w2_<X>_1..P, (theta_rc_1..P)`
#' * breaks: `a, b, (w2_rc,) w2_<X>_1..P, m_1..R, (theta_rc_1..R,)
#'   theta_<X>_1..R`
#' * pure realized-correlation model without breaks: `a, b, w2_rc` only.
#'
#' `P = N(N-1)/2` asset pairs, `R = 1 + n_breaks` regimes. Without breaks the
#' level/slope coefficients are pair-specific; with breaks they are shared
#' across pairs but regime-varying, while the `w2` weights stay
#' pair-specific — exactly the legacy design.
#' 相关步参数索引:无突变时系数逐资产对,有突变时跨对共享、按区制变化
#' (w2 始终逐对),与旧脚本设计一致。
#'
#' @param n_assets Number of assets `N`.
#' @param rc Does the long-run component include the realized-correlation
#'   term?
#' @param n_X Number of exogenous correlation covariates `J`.
#' @param n_regimes Number of regimes `R`.
#' @param X_names Covariate names for labelling.
#'
#' @return List of index vectors plus `n_par` and `names`.
#' @keywords internal
dcc_par_index <- function(n_assets, rc, n_X, n_regimes,
                          X_names = paste0("X", seq_len(n_X))) {
  P <- n_assets * (n_assets - 1) / 2
  R <- n_regimes; J <- n_X

  idx <- list(a = 1L, b = 2L, P = P)
  pos <- 2L
  nm  <- c("a", "b")

  if (rc) {
    idx$w2_rc <- 3L; pos <- 3L
    nm <- c(nm, "w2_RC")
  }

  pair_lab <- seq_len(P)
  reg_lab  <- seq_len(R)

  if (R == 1L) {
    if (J > 0L) {
      idx$delta <- pos + seq_len(P); pos <- pos + P
      nm <- c(nm, paste0("delta_", pair_lab))
      idx$theta <- lapply(seq_len(J), function(j) pos + (j - 1L) * P + seq_len(P))
      pos <- pos + J * P
      nm <- c(nm, unlist(lapply(X_names, function(x) paste0("theta_", x, "_", pair_lab))))
      idx$w2 <- lapply(seq_len(J), function(j) pos + (j - 1L) * P + seq_len(P))
      pos <- pos + J * P
      nm <- c(nm, unlist(lapply(X_names, function(x) paste0("w2_", x, "_", pair_lab))))
      if (rc) {
        idx$theta_rc <- pos + seq_len(P); pos <- pos + P
        nm <- c(nm, paste0("theta_RC_", pair_lab))
      }
    }
  } else {
    if (J > 0L) {
      idx$w2 <- lapply(seq_len(J), function(j) pos + (j - 1L) * P + seq_len(P))
      pos <- pos + J * P
      nm <- c(nm, unlist(lapply(X_names, function(x) paste0("w2_", x, "_", pair_lab))))
    }
    idx$m <- pos + seq_len(R); pos <- pos + R
    nm <- c(nm, paste0("m_", reg_lab))
    if (rc) {
      idx$theta_rc <- pos + seq_len(R); pos <- pos + R
      nm <- c(nm, paste0("theta_RC_", reg_lab))
    }
    if (J > 0L) {
      idx$theta <- lapply(seq_len(J), function(j) pos + (j - 1L) * R + seq_len(R))
      pos <- pos + J * R
      nm <- c(nm, unlist(lapply(X_names, function(x) paste0("theta_", x, "_", reg_lab))))
    }
  }

  idx$n_par <- pos
  idx$names <- nm
  idx
}

#' Parameter specification for the DCC-MIDAS optimizer
#'
#' Start values and inequality constraints for the second estimation step,
#' assembled from the model configuration (same contract as
#' [gm_param_spec()]). Models whose legacy branches used a single fixed
#' start (`a = 0.01, b = 0.8`, pure RC) are reproduced by fixing all start
#' values; the multistart models draw `a`, `b` and every level/slope
#' coefficient at random.
#' 第二步参数规格:单起点模型全部固定初值,多起点模型按旧脚本抽取。
#'
#' @inheritParams dcc_par_index
#' @param lag_fun `"Beta"` or `"Almon"`.
#' @param bounds Named list overriding entries of the internal
#'   `dcc_bounds_default` constants table.
#' @return Same structure as [gm_param_spec()].
#' @keywords internal
dcc_param_spec <- function(n_assets, rc, n_X, n_regimes, lag_fun,
                           X_names = paste0("X", seq_len(n_X)),
                           bounds = list()) {
  b <- utils::modifyList(dcc_bounds_default, bounds)
  idx <- dcc_par_index(n_assets, rc, n_X, n_regimes, X_names)
  n_par <- idx$n_par
  R <- n_regimes; J <- n_X; P <- idx$P

  # Pure RC without breaks was a fixed-start model in the legacy script.
  # 纯 RC 无突变模型在旧脚本中为固定初值单起点。
  single_start <- (J == 0L && R == 1L)

  w2_fix <- if (lag_fun == "Beta") b$w2_start_beta else b$w2_start_almon
  starts <- data.frame(lo = rep(NA_real_, n_par), hi = rep(NA_real_, n_par),
                       fixed = rep(NA_real_, n_par))
  if (single_start) {
    starts[idx$a, "fixed"] <- b$a_fixed
    starts[idx$b, "fixed"] <- b$b_fixed
  } else {
    starts[idx$a, c("lo", "hi")] <- rbind(b$a_start)
    starts[idx$b, c("lo", "hi")] <- rbind(b$b_start)
  }
  if (rc) starts[idx$w2_rc, "fixed"] <- w2_fix
  if (J > 0L) {
    for (j in seq_len(J)) starts[idx$w2[[j]], "fixed"] <- w2_fix
    if (R == 1L) {
      starts[idx$delta, c("lo", "hi")] <- matrix(b$coef_start, P, 2, byrow = TRUE)
      for (j in seq_len(J)) starts[idx$theta[[j]], c("lo", "hi")] <- matrix(b$coef_start, P, 2, byrow = TRUE)
      if (rc) starts[idx$theta_rc, c("lo", "hi")] <- matrix(b$coef_start, P, 2, byrow = TRUE)
    }
  }
  if (R > 1L) {
    starts[idx$m, c("lo", "hi")] <- matrix(b$coef_start, R, 2, byrow = TRUE)
    if (rc) starts[idx$theta_rc, c("lo", "hi")] <- matrix(b$coef_start, R, 2, byrow = TRUE)
    if (J > 0L) for (j in seq_len(J)) starts[idx$theta[[j]], c("lo", "hi")] <- matrix(b$coef_start, R, 2, byrow = TRUE)
  }

  # ---- constraints ui %*% theta + ci >= 0 ----
  unit <- function(i) { v <- numeric(n_par); v[i] <- 1; v }
  ui <- rbind(unit(idx$a), unit(idx$b), -unit(idx$a) - unit(idx$b))
  ci <- c(-b$a_min, -b$b_min, b$ab_sum_max)

  w2_all <- c(if (rc) idx$w2_rc, if (J > 0L) unlist(idx$w2))
  for (p in w2_all) {
    if (lag_fun == "Beta") {
      ui <- rbind(ui, unit(p));  ci <- c(ci, -b$w2_min_beta)
    } else {
      ui <- rbind(ui, -unit(p)); ci <- c(ci, b$w2_max_almon)
    }
  }

  list(idx = idx, names = idx$names, ui = ui, ci = ci,
       n_par = n_par, starts = starts, single_start = single_start)
}
