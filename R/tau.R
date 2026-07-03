#' Long-run (MIDAS) variance component
#'
#' The single builder of the daily long-run component shared by the
#' likelihood, the component extractor and (in its correlation flavour) the
#' DCC-MIDAS code. Implements
#' \deqn{\log \tau_t = m(t) + \sum_j \theta_j(t)\, \mathrm{MIDAS}(X_j; w_{2,j})_t}
#' where \eqn{m(t)} and every \eqn{\theta_j(t)} are piecewise-constant across
#' the regimes encoded in `D_mat`. Each old model variant (RV / X / RV_X,
#' with or without `_sb`) is a special case with 1–2 covariates and 0+ breaks.
#' 唯一的长期成分构造器:log tau = 分段截距 + 各协变量分段斜率×MIDAS 滤波;
#' 旧模型 6 变体均为其 1–2 个协变量、0+ 突变的特例。
#'
#' @param m_vec Numeric vector of regime intercepts, length `nrow(D_mat)`.
#' @param theta List (one element per covariate) of regime-slope vectors,
#'   each of length `nrow(D_mat)`.
#' @param w2 Numeric vector of second weighting parameters, one per covariate.
#' @param lag_mats List of `(K_j + 1) x TT` covariate lag matrices (from
#'   [midas_lag_matrix()]).
#' @param K Integer vector of MIDAS lag orders, one per covariate.
#' @param D_mat Step-dummy matrix from [build_break_dummies()].
#' @param lag_fun `"Beta"` or `"Almon"`.
#'
#' @return Numeric vector `tau` of length `TT` (daily long-run variance).
#' @keywords internal
build_tau <- function(m_vec, theta, w2, lag_mats, K, D_mat, lag_fun) {

  # Piecewise intercept / 分段截距序列
  log_tau <- as.vector(m_vec %*% D_mat)

  for (j in seq_along(lag_mats)) {
    # MIDAS filter of covariate j: weighted rolling sum over the lag matrix;
    # row K_j + 1 holds the completed weighted sum for every day.
    # 第 j 个协变量的 MIDAS 滤波:对滞后矩阵做加权滚动和,第 K_j+1 行
    # 即每个交易日的完整加权和。suppressWarnings 与旧脚本保持一致
    # (TODO: 显式处理 NA 后移除)。
    betas <- midas_filter_weights(K[j], w2[j], lag_fun)
    filt  <- suppressWarnings(roll::roll_sum(lag_mats[[j]], c(K[j] + 1), weights = betas))
    filt  <- filt[(K[j] + 1), ]

    # Piecewise slope for covariate j / 第 j 个协变量的分段斜率
    theta_t <- as.vector(theta[[j]] %*% D_mat)

    log_tau <- log_tau + theta_t * filt
  }

  exp(log_tau)
}
