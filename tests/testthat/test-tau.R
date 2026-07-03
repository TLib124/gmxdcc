# Gate 4: the generalized weight / break / tau builders must reproduce the
# legacy inline computations for every old model variant they replace.
# 关卡 4:泛化的权重/突变/长期成分构造器必须复现旧脚本各变体的内联计算。

test_that("midas_weights is numerically identical to rumidas", {
  skip_if_not_installed("rumidas")

  for (K in c(6, 13, 25)) {
    k <- 1:K
    for (w2 in c(1.5, 3, 7.2)) {
      expect_identical(midas_weights(k, K, 1, w2, "Beta"),
                       rumidas::beta_function(k, K, 1, w2))
    }
    for (w2 in c(-0.05, -0.2)) {
      expect_identical(midas_weights(k, K, 0, w2, "Almon"),
                       rumidas::exp_almon(k, K, 0, w2))
    }
  }
})

test_that("build_break_dummies reproduces the legacy D_mat construction", {
  d <- sim_data()
  full_days <- as.Date(zoo::index(d$r1))
  TT <- length(full_days)
  breaks <- as.Date(c("2014-06-01", "2017-01-15"))

  # Legacy construction, transcribed from lines 542-548 / 旧脚本逐行转写
  n_n <- length(breaks)
  D_legacy <- matrix(0, nrow = 1 + n_n, ncol = TT)
  D_legacy[1, ] <- 1
  for (i in 1:n_n) D_legacy[i + 1, full_days >= breaks[i]] <- 1

  expect_identical(build_break_dummies(full_days, breaks), D_legacy)

  # No breaks -> single constant row / 无突变退化为单行全 1
  expect_identical(build_break_dummies(full_days, NULL),
                   matrix(1, nrow = 1, ncol = TT))

  # Invalid inputs / 非法输入
  expect_error(build_break_dummies(full_days, rev(breaks)), "increasing")
  expect_error(build_break_dummies(full_days, as.Date("2030-01-01")), "inside")
})

test_that("build_tau reproduces the legacy tau for all structural cases", {
  skip_if(is.null(legacy_env()), "dev/legacy not available")
  lg <- legacy_env(); d <- sim_data()

  # Align the daily calendar to the RV series (which starts with the returns)
  # so that K+1 months of RV history exist for every retained day, mirroring
  # the legacy wrapper. 按 RV 日历裁剪日收益,保证每个交易日都有 K+1 个月
  # 的 RV 历史,与旧封装函数一致。
  K <- 12
  RV_macro <- realized_measure(d$r1, "monthly")
  x  <- align_sample(d$r1, K = K, macro = RV_macro, period = "monthly")
  TT <- length(x)
  RV_m <- midas_lag_matrix(x, RV_macro, K = K, period = "monthly")
  mv_m <- midas_lag_matrix(x, d$macro, K = K, period = "monthly")
  full_days <- as.Date(zoo::index(x))

  for (lag_fun in c("Beta", "Almon")) {
    w1 <- ifelse(lag_fun == "Beta", 1, 0)
    weight_fun <- if (lag_fun == "Beta") lg$weight_func_beta else lg$weight_func_exp
    w2_RV <- ifelse(lag_fun == "Beta", 3, -0.05)
    w2_X  <- ifelse(lag_fun == "Beta", 2, -0.02)

    ## Case (a): one covariate, no break — legacy GM_noskew_RV lines 323-325
    ## 情形 (a):单协变量无突变
    m <- 0.1; theta <- 0.2
    betas   <- c(rev(weight_fun(1:(K + 1), (K + 1), w1, w2_RV))[2:(K + 1)], 0)
    tau_leg <- m + theta * suppressWarnings(roll::roll_sum(RV_m, c(K + 1), weights = betas))
    tau_leg <- exp(tau_leg[(K + 1), ])

    tau_new <- build_tau(
      m_vec = m, theta = list(theta), w2 = w2_RV,
      lag_mats = list(RV_m), K = K,
      D_mat = build_break_dummies(full_days, NULL), lag_fun = lag_fun
    )
    expect_equal(tau_new, tau_leg, tolerance = 1e-12, label = paste("case a", lag_fun))

    ## Case (b): two covariates, no break — legacy GM_noskew_RV_X lines 455-465
    ## 情形 (b):双协变量无突变
    th_RV <- 0.15; th_X <- -0.1
    betas_RV <- c(rev(weight_fun(1:(K + 1), (K + 1), w1, w2_RV))[2:(K + 1)], 0)
    f_RV     <- suppressWarnings(roll::roll_sum(RV_m, c(K + 1), weights = betas_RV))[(K + 1), ]
    betas_X  <- c(rev(weight_fun(1:(K + 1), (K + 1), w1, w2_X))[2:(K + 1)], 0)
    f_X      <- suppressWarnings(roll::roll_sum(mv_m, c(K + 1), weights = betas_X))[(K + 1), ]
    tau_leg2 <- exp(m + th_RV * f_RV + th_X * f_X)

    tau_new2 <- build_tau(
      m_vec = m, theta = list(th_RV, th_X), w2 = c(w2_RV, w2_X),
      lag_mats = list(RV_m, mv_m), K = c(K, K),
      D_mat = build_break_dummies(full_days, NULL), lag_fun = lag_fun
    )
    expect_equal(tau_new2, tau_leg2, tolerance = 1e-12, label = paste("case b", lag_fun))

    ## Case (c): two covariates + two breaks — legacy GM_noskew_RV_X_sb
    ## lines 742-779 / 情形 (c):双协变量 + 两个突变
    breaks    <- as.Date(c("2014-06-01", "2017-01-15"))
    n_n       <- length(breaks)
    D_mat     <- matrix(0, nrow = 1 + n_n, ncol = TT)
    D_mat[1, ] <- 1
    for (i in 1:n_n) D_mat[i + 1, full_days >= breaks[i]] <- 1

    m_vec     <- c(0.1, -0.05, 0.08)
    th_RV_vec <- c(0.15, 0.05, -0.03)
    th_X_vec  <- c(-0.1, 0.02, 0.06)

    m_t     <- as.vector(m_vec %*% D_mat)
    thRV_t  <- as.vector(th_RV_vec %*% D_mat)
    thX_t   <- as.vector(th_X_vec %*% D_mat)
    tau_leg3 <- exp(m_t + thRV_t * f_RV + thX_t * f_X)

    tau_new3 <- build_tau(
      m_vec = m_vec, theta = list(th_RV_vec, th_X_vec), w2 = c(w2_RV, w2_X),
      lag_mats = list(RV_m, mv_m), K = c(K, K),
      D_mat = build_break_dummies(full_days, breaks), lag_fun = lag_fun
    )
    expect_equal(tau_new3, tau_leg3, tolerance = 1e-12, label = paste("case c", lag_fun))
  }
})
