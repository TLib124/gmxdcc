#' Variance components of a fitted GARCH-MIDAS model
#'
#' Computes the daily long-run component `tau_d`, the short-run component
#' `g_it` and the total conditional variance `h_it_2 = g_it * tau_d` at a
#' given parameter vector, replacing the twelve branches of the legacy
#' `garchmidas_noskew_est()`. The recursion is a pure filter: when `data`
#' spans in-sample plus out-of-sample days, the out-of-sample segment is
#' automatically warm-started from the last in-sample state (no cold reset),
#' which is the basis of the package's out-of-sample design.
#' 计算长期/短期成分与总条件方差;作为纯滤波器,传入全样本数据时
#' 样本外段自然从样本内末期状态热启动,无冷重置。
#'
#' @param param Numeric parameter vector ordered as in [gm_par_index()].
#' @param data Data bundle as in [gm_loglik()].
#'
#' @return List of xts series `h_it_2`, `g_it`, `tau_d` indexed like
#'   `data$daily_ret`.
#' @keywords internal
gm_components <- function(param, data) {
  idx <- data$idx

  alpha <- param[idx$alpha]
  beta  <- param[idx$beta]
  w2    <- param[idx$w2]
  m_vec <- param[idx$m]
  theta <- lapply(idx$theta, function(ii) param[ii])
  mu    <- param[idx$mu]

  daily_ret <- data$daily_ret
  epsilon   <- daily_ret - mu
  TT        <- length(daily_ret)

  ##### daily long-run / 日频长期成分
  tau_d <- build_tau(
    m_vec = m_vec, theta = theta, w2 = w2,
    lag_mats = data$lag_mats, K = data$K,
    D_mat = data$D_mat, lag_fun = data$lag_fun
  )

  ####### short-run: g_t = (1-a-b) + a*eps_{t-1}^2/tau_{t-1} + b*g_{t-1},
  ####### vectorized via stats::filter exactly as in the legacy code.
  ####### 短期递推,与旧脚本相同的 filter 向量化写法。
  step_1 <- (1 - alpha - beta) + (alpha) * (epsilon)^2 / tau_d
  step_1[is.na(step_1)] <- 0
  x_step_1 <- c(0, step_1[-TT])
  x_step_1[2] <- x_step_1[2] + beta * 1

  g_it <- as.numeric(stats::filter(x_step_1, filter = beta, method = "recursive"))
  g_it[1] <- 1

  ###### variance / 总条件方差
  h_it_2 <- zoo::coredata(g_it * tau_d)
  g_it   <- zoo::coredata(g_it)
  tau_d  <- zoo::coredata(tau_d)

  h_it_2 <- xts::as.xts(h_it_2, zoo::index(daily_ret))
  g_it   <- xts::as.xts(g_it,   zoo::index(daily_ret))
  tau_d  <- xts::as.xts(tau_d,  zoo::index(daily_ret))

  list(
    "h_it_2" = h_it_2,
    "g_it"   = g_it,
    "tau_d"  = tau_d
  )
}
