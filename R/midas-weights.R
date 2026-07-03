#' MIDAS lag weighting functions
#'
#' Beta and exponential-Almon lag weights, numerically identical to
#' `rumidas::beta_function()` and `rumidas::exp_almon()` (reimplemented here
#' to avoid the dependency; equality is checked in the test suite).
#'
#' @param k Integer vector of lag positions, typically `1:K`.
#' @param K Total number of lags.
#' @param w1 First weighting parameter. In the restricted Beta scheme used by
#'   the models in this package `w1 = 1` (monotonically decaying weights);
#'   for the Almon scheme `w1 = 0`.
#' @param w2 Second weighting parameter.
#' @param lag_fun `"Beta"` or `"Almon"`.
#'
#' @return Numeric vector of weights of the same length as `k`, summing to 1
#'   over `1:K`.
#' @examples
#' # Decaying Beta weights over 12 lags / 12 阶递减 Beta 权重
#' midas_weights(1:12, 12, w1 = 1, w2 = 3, lag_fun = "Beta")
#' # Exponential Almon / 指数 Almon 权重
#' midas_weights(1:12, 12, w1 = 0, w2 = -0.05, lag_fun = "Almon")
#' @export
midas_weights <- function(k, K, w1, w2, lag_fun = c("Beta", "Almon")) {
  lag_fun <- match.arg(lag_fun)
  j <- 1:K
  if (lag_fun == "Beta") {
    # Beta lag polynomial / Beta 滞后多项式 (rumidas::beta_function)
    num <- ((k / K)^(w1 - 1)) * (1 - k / K)^(w2 - 1)
    den <- sum(((j / K)^(w1 - 1)) * (1 - j / K)^(w2 - 1))
  } else {
    # Exponential Almon lag / 指数 Almon 滞后 (rumidas::exp_almon)
    num <- exp(w1 * k + w2 * k^2)
    den <- sum(exp(w1 * j + w2 * j^2))
  }
  num / den
}

#' Rolling-sum weight vector for the MIDAS filter
#'
#' Converts the raw lag weights into the vector expected by
#' `roll::roll_sum()` on a `(K+1) x N` lag matrix: weights are reversed
#' (oldest row first), the current period is excluded (trailing zero), and
#' `w1` is fixed at 1 (Beta) or 0 (Almon) exactly as in the legacy code.
#' 将原始权重转为 roll_sum 使用的形式:反转顺序、剔除当期(末位补 0),
#' w1 按 lag_fun 固定为 1(Beta)或 0(Almon),与旧脚本一致。
#'
#' @param K MIDAS lag order.
#' @param w2 Second weighting parameter.
#' @param lag_fun `"Beta"` or `"Almon"`.
#'
#' @return Numeric weight vector of length `K + 1`.
#' @keywords internal
midas_filter_weights <- function(K, w2, lag_fun) {
  w1 <- ifelse(lag_fun == "Beta", 1, 0)
  c(rev(midas_weights(1:(K + 1), (K + 1), w1, w2, lag_fun))[2:(K + 1)], 0)
}
