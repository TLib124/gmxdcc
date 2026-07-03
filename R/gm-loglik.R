#' Parameter index map for the univariate GARCH-MIDAS model
#'
#' Fixes the canonical parameter ordering used by [gm_loglik()],
#' [gm_components()] and the fitting driver:
#' `alpha, beta, w2_1..w2_J, (nu,) m_1..m_R, theta_1_1..theta_J_R, mu`
#' where `J` is the number of long-run covariates and `R = 1 + n_breaks` the
#' number of regimes. 统一的参数排序约定:短期 α/β、各协变量权重 w2、
#' (std 分布的自由度 ν、)分段截距 m、各协变量分段斜率 θ、均值 μ。
#'
#' @param n_cov Number of long-run covariates `J`.
#' @param n_regimes Number of regimes `R` (`1 + n_breaks`).
#' @param distribution `"norm"` or `"std"`.
#' @param cov_names Optional covariate names used to label parameters.
#'
#' @return List with index vectors `alpha`, `beta`, `w2`, `nu` (or `NULL`),
#'   `m`, `theta` (list of length `J`), `mu`, total length `n_par`, and the
#'   canonical `names` vector.
#' @keywords internal
gm_par_index <- function(n_cov, n_regimes, distribution,
                         cov_names = paste0("X", seq_len(n_cov))) {
  J <- n_cov; R <- n_regimes
  has_nu <- distribution == "std"

  idx <- list(alpha = 1L, beta = 2L)
  pos <- 2L

  idx$w2 <- pos + seq_len(J); pos <- pos + J
  idx$nu <- if (has_nu) pos + 1L else NULL
  pos    <- pos + as.integer(has_nu)

  idx$m <- pos + seq_len(R); pos <- pos + R
  idx$theta <- lapply(seq_len(J), function(j) pos + (j - 1L) * R + seq_len(R))
  pos <- pos + J * R

  idx$mu    <- pos + 1L
  idx$n_par <- pos + 1L

  # Human-readable parameter names / 参数命名
  reg_suffix <- if (R == 1L) "" else paste0("_", seq_len(R))
  idx$names <- c(
    "alpha", "beta",
    paste0("w2_", cov_names),
    if (has_nu) "nu",
    paste0("m", reg_suffix),
    unlist(lapply(cov_names, function(nm) paste0("theta_", nm, reg_suffix))),
    "mu"
  )
  idx
}

#' Per-observation log-likelihood of the GARCH-MIDAS model
#'
#' Replaces the twelve model-by-distribution branches of the legacy
#' `garchmidas_noskew_loglik()`: the long-run component comes from the single
#' generalized builder [build_tau()], the short-run recursion and likelihood
#' from the C++ kernels. Returns the per-observation vector required by
#' `maxLik` for OPG/sandwich standard errors.
#' 取代旧脚本 12 个分支的统一似然:长期成分走 build_tau,短期递推与似然
#' 走 C++ 核心;逐观测返回以便 maxLik 计算三明治标准误。
#'
#' @param param Numeric parameter vector ordered as in [gm_par_index()].
#' @param data List bundle with elements `daily_ret` (xts), `lag_mats`,
#'   `K`, `D_mat`, `lag_fun`, `distribution`, `idx`.
#'
#' @return Numeric vector of per-observation log-likelihood contributions.
#' @keywords internal
gm_loglik <- function(param, data) {
  idx <- data$idx

  alpha <- param[idx$alpha]
  beta  <- param[idx$beta]
  w2    <- param[idx$w2]
  m_vec <- param[idx$m]
  theta <- lapply(idx$theta, function(ii) param[ii])
  mu    <- param[idx$mu]

  daily_ret <- data$daily_ret
  epsilon   <- daily_ret - mu

  # Long-run component / 长期成分
  tau_d <- build_tau(
    m_vec = m_vec, theta = theta, w2 = w2,
    lag_mats = data$lag_mats, K = data$K,
    D_mat = data$D_mat, lag_fun = data$lag_fun
  )

  # Short-run recursion + likelihood in C++ / 短期递推与似然在 C++ 中完成
  if (data$distribution == "norm") {
    ll <- fast_garch_ll_cpp(
      alpha = alpha, beta = beta, mu = mu,
      epsilon = as.numeric(epsilon), tau_d = as.numeric(tau_d),
      daily_ret = as.numeric(daily_ret)
    )
  } else {
    ll <- fast_t_garch_ll(
      alpha = alpha, beta = beta, v = param[idx$nu],
      epsilon = as.numeric(epsilon), tau_d = as.numeric(tau_d),
      daily_ret = as.numeric(daily_ret)
    )
  }

  ll
}
