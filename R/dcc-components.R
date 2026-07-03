#' Conditional correlation and covariance matrices of a fitted DCC model
#'
#' Runs the correlation recursion at fixed parameters and returns the
#' `H_t`, `R_t` (and `R_t_bar` for DCC-MIDAS) arrays, replacing the legacy
#' `dccmidas_noskew_mat_est()` plus the direct `new_*_mat_est` calls.
#'
#' The function is a pure filter: called on the full (in-sample plus
#' out-of-sample) residual series with parameters estimated in-sample, the
#' out-of-sample segment continues seamlessly from the last in-sample state.
#' Pass `S_in` (and `N_in` for aDCC) computed from in-sample residuals so
#' the unconditional targets carry no look-ahead — this replaces both the
#' legacy cold restart and its `N_c_res_for_oos` patch mechanism.
#' 纯滤波器:对全样本残差与冻结参数运行,样本外自然热启动;S_in/N_in
#' 传入样本内目标以消除前视,取代旧的冷启动与 N_c_res_for_oos 补丁。
#'
#' @param param Estimated parameter vector.
#' @param model One of `"cDCC"`, `"aDCC"`, `"DECO"`, `"DCC-MIDAS"`.
#' @param data DCC data bundle (see [build_longrun_corr()]); for the
#'   non-MIDAS models only `res` and `K_c` are used.
#' @param D_t `N x N x TT` cube of conditional-standard-deviation diagonal
#'   matrices.
#' @param S_in Optional unconditional covariance target (in-sample).
#' @param N_in Optional aDCC asymmetric target (in-sample).
#'
#' @return List with `H_t`, `R_t`, and for DCC-MIDAS also `R_t_bar`.
#' @keywords internal
dcc_matrices <- function(param, model, data, D_t, S_in = NULL, N_in = NULL) {

  if (model == "cDCC") {
    return(new_dcc_mat_est(param, data$res, D_t, K_c_in = data$K_c, S_in = S_in))
  }
  if (model == "aDCC") {
    return(new_a_dcc_mat_est(param, data$res, D_t, K_c_in = data$K_c,
                             S_in = S_in, N_in = N_in))
  }
  if (model == "DECO") {
    return(new_deco_mat_est(param, data$res, D_t, K_c_in = data$K_c, S_in = S_in))
  }

  ## --- DCC-MIDAS: R-side recursion, ported verbatim from the legacy loop ---
  ## DCC-MIDAS:R 端递推,逐行移植旧脚本
  N <- data$n_assets; TT <- data$TT; K_c <- data$K_c
  a <- param[data$idx$a]; b <- param[data$idx$b]

  R_t_bar <- build_longrun_corr(param, data)
  for (i in 1:TT) diag(R_t_bar[, , i]) <- 1

  S <- if (is.null(S_in)) stats::cov(data$res) else S_in
  Q_t <- array(diag(rep(1, N)), dim = c(N, N, TT))
  R_t <- array(diag(rep(1, N)), dim = c(N, N, TT))
  H_t <- array(S, dim = c(N, N, TT))

  for (tt in (K_c + 1):TT) {
    Q_t[, , tt] <- (1 - a - b) * R_t_bar[, , tt] +
      a * tcrossprod(data$res[tt - 1, ]) + b * Q_t[, , tt - 1]

    q_sd <- sqrt(diag(Q_t[, , tt]))
    R_t[, , tt] <- Q_t[, , tt] / tcrossprod(q_sd)
    H_t[, , tt] <- D_t[, , tt] %*% R_t[, , tt] %*% D_t[, , tt]
  }

  list("H_t" = H_t, "R_t" = R_t, "R_t_bar" = R_t_bar)
}
