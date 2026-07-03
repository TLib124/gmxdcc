# Gate 7a: dcc_midas_loglik / dcc_matrices must reproduce the legacy
# dccmidas_noskew_loglik / dccmidas_noskew_mat_est for all six correlation
# variants x 2 lag functions at fixed parameters, and the S_in extension of
# the cpp kernels must default to legacy behavior.
# 关卡 7a:相关步似然与矩阵滤波在固定参数处复现旧脚本全部变体。

# Build a shared test setup: standardized residuals with dependence, macro
# lag matrix, D_t cube. 共享测试环境:带相关性的标准化残差与协变量。
dcc_setup <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    d <- sim_data()
    set.seed(99)
    K_c <- 6
    x  <- align_sample(d$r1, K = K_c, macro = d$macro, period = "monthly")
    TT <- length(x); N <- 3
    rho <- 0.5
    Sig <- matrix(rho, N, N); diag(Sig) <- 1
    L <- chol(Sig)
    res <- matrix(rnorm(TT * N), TT, N) %*% L
    mv_c <- midas_lag_matrix(x, d$macro, K = K_c, period = "monthly")
    D_t <- array(0, dim = c(N, N, TT))
    for (tt in seq_len(TT)) diag(D_t[, , tt]) <- abs(rnorm(N, 1, 0.1))
    cache <<- list(x = x, TT = TT, N = N, K_c = K_c, N_c = 20,
                   res = res, mv_c = mv_c, D_t = D_t,
                   days = as.Date(zoo::index(x)),
                   breaks = as.Date(c("2014-06-01", "2017-01-15")))
    cache
  }
})

# Legacy positional parameter vector per corr model / 旧脚本位置参数
legacy_corr_param <- function(model, p, P, n_b) {
  switch(model,
    "DCCMIDAS_RC"      = c(p$a, p$b, p$w2_rc),
    "DCCMIDAS_X"       = c(p$a, p$b, p$delta, p$theta, p$w2),
    "DCCMIDAS_RC_X"    = c(p$a, p$b, p$w2_rc, p$delta, p$theta, p$w2, p$theta_rc),
    "DCCMIDAS_RC_sb"   = c(p$a, p$b, p$w2_rc, p$m_vec, p$th_vec),
    "DCCMIDAS_X_sb"    = c(p$a, p$b, p$w2, p$m_vec, p$th_vec),
    "DCCMIDAS_RC_X_sb" = c(p$a, p$b, p$w2_rc, p$w2, p$m_vec, p$thrc_vec, p$th_vec)
  )
}

# Same parameters in the new canonical order (coincides with legacy layouts)
# 新排序与旧布局在各配置下一致,直接复用
new_corr_param <- legacy_corr_param

test_that("dcc_midas_loglik matches legacy for all variants", {
  skip_if(is.null(legacy_env()), "dev/legacy not available")
  lg <- legacy_env(); s <- dcc_setup()
  P <- s$N * (s$N - 1) / 2; n_b <- length(s$breaks)

  for (lag_fun in c("Beta", "Almon")) {
    p <- list(
      a = 0.03, b = 0.9,
      w2_rc = ifelse(lag_fun == "Beta", 2, -0.05),
      w2    = rep(ifelse(lag_fun == "Beta", 1.5, -0.03), P),
      delta = c(0.2, -0.1, 0.15), theta = c(0.3, -0.2, 0.1),
      theta_rc = c(0.4, 0.2, -0.1),
      m_vec = c(0.2, -0.1, 0.15), th_vec = c(0.3, 0.1, -0.2),
      thrc_vec = c(0.4, -0.15, 0.05)
    )

    for (model in c("DCCMIDAS_RC", "DCCMIDAS_X", "DCCMIDAS_RC_X",
                    "DCCMIDAS_RC_sb", "DCCMIDAS_X_sb", "DCCMIDAS_RC_X_sb")) {
      sb <- grepl("_sb$", model)
      has_rc <- grepl("RC", model)
      has_x  <- grepl("X(_sb)?$", model) | grepl("_X", model)
      n_X <- if (grepl("X", sub("_sb$", "", sub("RC", "", model)))) 1L else 0L
      R <- if (sb) 1L + n_b else 1L

      par_vec <- legacy_corr_param(model, p, P, n_b)
      lab <- paste(model, lag_fun)

      ll_leg <- lg$dccmidas_noskew_loglik(
        par_vec, model = model, res = s$res, days_period = s$days,
        mv_c = if (n_X) s$mv_c else NULL, lag_fun = lag_fun,
        N_c = s$N_c, K_c = s$K_c,
        X_stru_break_date_vec = if (sb) s$breaks else NULL
      )

      dc <- list(res = s$res, n_assets = s$N, TT = s$TT, K_c = s$K_c,
                 N_c = s$N_c, rc = has_rc,
                 X_mats = if (n_X) list(s$mv_c) else list(),
                 D_mat = build_break_dummies(s$days, if (sb) s$breaks else NULL),
                 lag_fun = lag_fun,
                 idx = gmxdcc:::dcc_par_index(s$N, rc = has_rc, n_X = n_X,
                                              n_regimes = R))
      ll_new <- gmxdcc:::dcc_midas_loglik(par_vec, dc)
      expect_equal(ll_new, as.numeric(ll_leg), tolerance = 1e-12, label = lab)

      ## matrices at the same parameters / 同参数下的矩阵滤波
      mats_leg <- lg$dccmidas_noskew_mat_est(
        par_vec, model = model, res = s$res, days_period = s$days,
        D_t = s$D_t, mv_c = if (n_X) s$mv_c else NULL, lag_fun = lag_fun,
        N_c = s$N_c, K_c = s$K_c,
        X_stru_break_date_vec = if (sb) s$breaks else NULL
      )
      mats_new <- gmxdcc:::dcc_matrices(par_vec, "DCC-MIDAS", dc, s$D_t)
      expect_equal(mats_new$H_t, mats_leg$H_t, tolerance = 1e-12, label = paste("H", lab))
      expect_equal(mats_new$R_t, mats_leg$R_t, tolerance = 1e-12, label = paste("R", lab))
    }
  }
})

test_that("cpp mat_est with S_in = NULL keeps legacy behavior; S_in overrides", {
  s <- dcc_setup()

  ref <- new_dcc_mat_est(c(0.03, 0.9), s$res, s$D_t, K_c = s$K_c)
  # Default must equal the pre-S_in behavior; arma::cov and stats::cov agree
  # only to machine precision, hence the 1e-12 tolerance.
  # 缺省行为 = 旧行为;arma 与 R 的 cov 仅在机器精度上一致,容差 1e-12。
  via_sin <- new_dcc_mat_est(c(0.03, 0.9), s$res, s$D_t, K_c = s$K_c,
                             S_in = stats::cov(s$res))
  expect_equal(ref, via_sin, tolerance = 1e-12)

  # A different S changes the recursion / 不同的 S 改变递推
  other <- new_dcc_mat_est(c(0.03, 0.9), s$res, s$D_t, K_c = s$K_c,
                           S_in = diag(s$N))
  expect_false(isTRUE(all.equal(ref$R_t, other$R_t)))
})
