# ---------------------------------------------------------------------------
# Central table of optimizer defaults, extracted verbatim from the 21 legacy
# fit branches. Overridable per fit via `control$bounds`.
# 优化器默认常数表:与旧脚本 21 个分支一致,可经 control$bounds 覆盖。
# 唯一的整理:旧脚本 sb-Almon 分支的 w2 初值在 -0.01 / -0.1 间不一致,
# 此处统一为 -0.1。
# ---------------------------------------------------------------------------
gm_bounds_default <- list(
  alpha_start = c(0.001, 0.095),  # alpha ~ U(0.001, 0.095)
  beta_start  = c(0.6, 0.8),      # beta  ~ U(0.6, 0.8)
  coef_start  = c(-1, 1),         # m, theta ~ U(-1, 1)
  nu_start    = c(2.01, 10),      # Student-t df ~ U(2.01, 10)
  w2_start_beta  = 1.01,          # fixed start / 固定初值
  w2_start_almon = -0.1,
  alpha_min   = 0.001,            # alpha >  0.001
  beta_min    = 0.001,            # beta  >  0.001
  ab_sum_max  = 0.999,            # alpha + beta < 0.999
  w2_min_beta = 1.001,            # Beta:  w2 > 1.001
  w2_max_almon = -0.001,          # Almon: w2 < -0.001
  nu_min      = 2.001             # nu > 2.001
)

#' Parameter specification for the univariate GARCH-MIDAS optimizer
#'
#' Assembles, from the model dimensions alone, everything the optimizer
#' needs: canonical parameter names and index map ([gm_par_index()]),
#' random-start ranges, and the `maxLik` inequality constraints
#' `ui %*% theta + ci >= 0`. Replaces the 21 hand-written legacy fit
#' branches. 由模型维度自动组装参数名、初值区间与不等式约束,
#' 取代旧脚本 21 个手写分支。
#'
#' @param n_cov Number of long-run covariates `J`.
#' @param n_regimes Number of regimes `R = 1 + n_breaks`.
#' @param distribution `"norm"` or `"std"`.
#' @param lag_fun `"Beta"` or `"Almon"`.
#' @param cov_names Covariate names for parameter labelling.
#' @param mu_start Fixed starting value for `mu` (legacy: in-sample mean).
#' @param bounds Optional named list overriding entries of the default
#'   bounds table.
#'
#' @return List with `idx`, `names`, `ui`, `ci`, `n_par` and the start-value
#'   description used by [draw_starts()].
#' @keywords internal
gm_param_spec <- function(n_cov, n_regimes, distribution, lag_fun,
                          cov_names = paste0("X", seq_len(n_cov)),
                          mu_start = 0, bounds = list()) {
  b <- utils::modifyList(gm_bounds_default, bounds)
  J <- n_cov; R <- n_regimes
  idx <- gm_par_index(J, R, distribution, cov_names)
  n_par <- idx$n_par

  # ---- start-value description, one row per parameter ----
  # 每个参数一行:随机区间 (lo,hi) 或固定值 (fixed)
  w2_fix <- if (lag_fun == "Beta") b$w2_start_beta else b$w2_start_almon
  starts <- data.frame(
    lo = rep(NA_real_, n_par), hi = rep(NA_real_, n_par),
    fixed = rep(NA_real_, n_par)
  )
  starts[idx$alpha, c("lo", "hi")] <- rbind(b$alpha_start)
  starts[idx$beta,  c("lo", "hi")] <- rbind(b$beta_start)
  starts[idx$w2, "fixed"] <- w2_fix
  if (!is.null(idx$nu)) starts[idx$nu, c("lo", "hi")] <- rbind(b$nu_start)
  starts[idx$m, c("lo", "hi")] <- matrix(b$coef_start, R, 2, byrow = TRUE)
  for (j in seq_len(J)) {
    starts[idx$theta[[j]], c("lo", "hi")] <- matrix(b$coef_start, R, 2, byrow = TRUE)
  }
  starts[idx$mu, "fixed"] <- mu_start

  # ---- inequality constraints ui %*% theta + ci >= 0 ----
  # alpha > alpha_min; beta > beta_min; alpha + beta < ab_sum_max;
  # per covariate: Beta w2 > w2_min_beta, Almon w2 < w2_max_almon;
  # std: nu > nu_min. 其余参数无约束,与旧脚本一致。
  unit <- function(i) { v <- numeric(n_par); v[i] <- 1; v }
  ui <- rbind(
    unit(idx$alpha),
    unit(idx$beta),
    -unit(idx$alpha) - unit(idx$beta)
  )
  ci <- c(-b$alpha_min, -b$beta_min, b$ab_sum_max)
  for (j in seq_len(J)) {
    if (lag_fun == "Beta") {
      ui <- rbind(ui, unit(idx$w2[j]));  ci <- c(ci, -b$w2_min_beta)
    } else {
      ui <- rbind(ui, -unit(idx$w2[j])); ci <- c(ci, b$w2_max_almon)
    }
  }
  if (!is.null(idx$nu)) {
    ui <- rbind(ui, unit(idx$nu)); ci <- c(ci, -b$nu_min)
  }

  list(idx = idx, names = idx$names, ui = ui, ci = ci,
       n_par = n_par, starts = starts)
}

#' Draw the random multi-start matrix
#'
#' Draws `n_starts` candidate parameter vectors column by column: uniform
#' within the spec ranges for random parameters, constant for fixed ones —
#' the same scheme as the legacy `begin_val` matrices.
#' 按参数规格逐列生成候选初值矩阵,与旧脚本 begin_val 同构。
#'
#' @param spec Output of [gm_param_spec()] (or the DCC analogue).
#' @param n_starts Number of candidate start vectors.
#'
#' @return `n_starts x n_par` matrix with parameter names as column names.
#' @keywords internal
draw_starts <- function(spec, n_starts) {
  s <- spec$starts
  begin_val <- matrix(NA_real_, nrow = n_starts, ncol = spec$n_par)
  colnames(begin_val) <- spec$names
  for (p in seq_len(spec$n_par)) {
    begin_val[, p] <- if (!is.na(s$fixed[p])) {
      s$fixed[p]
    } else {
      stats::runif(n_starts, min = s$lo[p], max = s$hi[p])
    }
  }
  begin_val
}
