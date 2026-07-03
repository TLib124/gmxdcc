#' Per-observation log-likelihood of the DCC-MIDAS correlation step
#'
#' Builds the long-run correlation target with [build_longrun_corr()] and
#' evaluates the DCC quasi-likelihood recursion in C++ (`midas_ll_rcpp`),
#' replacing the twelve legacy `dccmidas_noskew_loglik` branches. The
#' per-observation return value feeds `maxLik`'s OPG/sandwich machinery.
#' DCC-MIDAS 相关步统一似然:长期目标出自 build_longrun_corr,递推与
#' 似然在 C++ 完成;逐观测返回供三明治标准误使用。
#'
#' @param param Parameter vector in [dcc_par_index()] order.
#' @param data Data bundle as in [build_longrun_corr()].
#'
#' @return Numeric vector of per-observation log-likelihood contributions.
#' @keywords internal
dcc_midas_loglik <- function(param, data) {
  a <- param[data$idx$a]
  b <- param[data$idx$b]
  R_t_bar <- build_longrun_corr(param, data)
  as.numeric(midas_ll_rcpp(data$TT, data$n_assets, data$K_c, a, b, R_t_bar, data$res))
}

#' Per-observation log-likelihoods of the non-MIDAS correlation models
#'
#' Thin wrappers around the C++ kernels for cDCC, aDCC and DECO, sharing the
#' `loglik(param, data)` contract used by [run_multistart()].
#' cDCC/aDCC/DECO 的 C++ 似然薄封装,统一 (param, data) 调用约定。
#'
#' @param param `c(a, b)` for cDCC/DECO, `c(a, b, g)` for aDCC.
#' @param data List with `res` (residual matrix) and `K_c`.
#' @return Numeric vector of per-observation log-likelihood contributions.
#' @keywords internal
cdcc_loglik_fn <- function(param, data) {
  as.numeric(new_dcc_loglik(param, data$res, data$K_c))
}

#' @rdname cdcc_loglik_fn
#' @keywords internal
adcc_loglik_fn <- function(param, data) {
  as.numeric(new_a_dcc_loglik(param, data$res, data$K_c))
}

#' @rdname cdcc_loglik_fn
#' @keywords internal
deco_loglik_fn <- function(param, data) {
  as.numeric(new_deco_loglik(param, data$res, data$K_c))
}
