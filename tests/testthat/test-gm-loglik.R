# Gate 5: gm_loglik / gm_components must reproduce the legacy
# garchmidas_noskew_loglik / garchmidas_noskew_est for all 6 model variants
# x 2 distributions x 2 lag functions at fixed feasible parameters.
# 关卡 5:统一似然/成分函数在固定参数处复现旧脚本全部 12+ 个分支。

# Map a named test-parameter set to the LEGACY positional vector of a given
# variant. 将命名参数映射为旧脚本各变体的位置参数向量。
legacy_param <- function(model, distribution, p, n_breaks = 2) {
  nb <- n_breaks
  switch(model,
    "GM_noskew_RV"   = ,
    "GM_noskew_X"    = c(p$alpha, p$beta, p$m, p$th1, p$w2_1,
                         if (distribution == "std") p$nu, p$mu),
    "GM_noskew_RV_X" = c(p$alpha, p$beta, p$m, p$th1, p$th2, p$w2_1, p$w2_2,
                         if (distribution == "std") p$nu, p$mu),
    "GM_noskew_RV_sb" = ,
    "GM_noskew_X_sb"  = c(p$alpha, p$beta, p$w2_1,
                          if (distribution == "std") p$nu,
                          p$m_vec, p$th1_vec, p$mu),
    "GM_noskew_RV_X_sb" = c(p$alpha, p$beta, p$w2_1, p$w2_2,
                            if (distribution == "std") p$nu,
                            p$m_vec, p$th1_vec, p$th2_vec, p$mu)
  )
}

# Same parameter set in the NEW canonical order of gm_par_index().
# 同一组参数按新排序排列。
new_param <- function(model, distribution, p) {
  sb  <- grepl("_sb$", model)
  two <- grepl("RV_X", model)
  m     <- if (sb) p$m_vec   else p$m
  th1   <- if (sb) p$th1_vec else p$th1
  th2   <- if (sb) p$th2_vec else p$th2
  c(p$alpha, p$beta,
    p$w2_1, if (two) p$w2_2,
    if (distribution == "std") p$nu,
    m, th1, if (two) th2,
    p$mu)
}

test_that("gm_loglik and gm_components match legacy for all variants", {
  skip_if(is.null(legacy_env()), "dev/legacy not available")
  lg <- legacy_env(); d <- sim_data()

  K <- 12
  RV_macro <- realized_measure(d$r1, "monthly")
  x    <- align_sample(d$r1, K = K, macro = RV_macro, period = "monthly")
  RV_m <- midas_lag_matrix(x, RV_macro, K = K, period = "monthly")
  mv_m <- midas_lag_matrix(x, d$macro,  K = K, period = "monthly")
  full_days <- as.Date(zoo::index(x))
  breaks <- as.Date(c("2014-06-01", "2017-01-15"))

  models <- c("GM_noskew_RV", "GM_noskew_X", "GM_noskew_RV_X",
              "GM_noskew_RV_sb", "GM_noskew_X_sb", "GM_noskew_RV_X_sb")

  for (lag_fun in c("Beta", "Almon")) {
    p <- list(
      alpha = 0.05, beta = 0.9, mu = 0.02, nu = 8,
      m = 0.1, th1 = 0.15, th2 = -0.1,
      w2_1 = ifelse(lag_fun == "Beta", 3, -0.05),
      w2_2 = ifelse(lag_fun == "Beta", 2, -0.02),
      m_vec = c(0.1, -0.05, 0.08),
      th1_vec = c(0.15, 0.05, -0.03),
      th2_vec = c(-0.1, 0.02, 0.06)
    )

    for (model in models) {
      sb  <- grepl("_sb$", model)
      two <- grepl("RV_X", model)
      # Covariate set per variant: RV, X, or both / 各变体的协变量集合
      lag_mats <- switch(sub("_sb$", "", model),
        "GM_noskew_RV"   = list(RV_m),
        "GM_noskew_X"    = list(mv_m),
        "GM_noskew_RV_X" = list(RV_m, mv_m)
      )
      J <- length(lag_mats)
      R <- if (sb) 1L + length(breaks) else 1L
      D_mat <- build_break_dummies(full_days, if (sb) breaks else NULL)

      for (distribution in c("norm", "std")) {
        data <- list(
          daily_ret = x, lag_mats = lag_mats, K = rep(K, J),
          D_mat = D_mat, lag_fun = lag_fun, distribution = distribution,
          idx = gm_par_index(J, R, distribution)
        )
        lab <- paste(model, distribution, lag_fun)

        ll_new <- gm_loglik(new_param(model, distribution, p), data)
        ll_leg <- lg$garchmidas_noskew_loglik(
          legacy_param(model, distribution, p),
          model = model, daily_ret = x,
          RV_m = if (grepl("RV", model)) RV_m else NULL,
          mv_m = if (grepl("X", sub("_sb$", "", model))) mv_m else NULL,
          K = K, distribution = distribution, lag_fun = lag_fun,
          X_stru_break_date_vec = if (sb) breaks else NULL
        )
        expect_equal(as.numeric(ll_new), as.numeric(ll_leg),
                     tolerance = 1e-12, label = paste("ll", lab))

        cmp_new <- gm_components(new_param(model, distribution, p), data)
        cmp_leg <- lg$garchmidas_noskew_est(
          legacy_param(model, distribution, p),
          model = model, daily_ret = x,
          RV_m = if (grepl("RV", model)) RV_m else NULL,
          mv_m = if (grepl("X", sub("_sb$", "", model))) mv_m else NULL,
          K = K, distribution = distribution, lag_fun = lag_fun,
          X_stru_break_date_vec = if (sb) breaks else NULL
        )
        for (comp in c("h_it_2", "g_it", "tau_d")) {
          expect_equal(as.numeric(cmp_new[[comp]]), as.numeric(cmp_leg[[comp]]),
                       tolerance = 1e-12, label = paste(comp, lab))
        }
      }
    }
  }
})
