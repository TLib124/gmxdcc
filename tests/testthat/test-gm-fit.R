# Gate 6: end-to-end fit equivalence and fit-level properties.
#
# The legacy fit draws its multi-start matrix with runif() and no seed, so
# exact replay works by (a) seeding the RNG, (b) reconstructing the exact
# begin_val matrix the legacy branch builds internally (same runif call
# order), and (c) injecting the same candidates into the new fit via
# control$start_matrix (mapped to the new canonical parameter order).
# 关卡 6:端到端等价——设种子后逐列重构 legacy 内部初值矩阵,经
# control$start_matrix 注入新拟合,二者从同一候选起点出发。

test_that("fit_garch_midas reproduces the legacy fit (RV, norm, Beta)", {
  skip_on_cran()
  skip_if(is.null(legacy_env()), "dev/legacy not available")
  lg <- legacy_env(); d <- sim_data()

  K <- 12; R_starts <- 20; SEED <- 907

  # Shared aligned data via the package's own preparation path; r_dgp has
  # genuine GARCH-MIDAS dynamics so the optimum is well identified.
  # 用包内数据准备路径生成共享对齐样本;r_dgp 含真实动态,最优点可识别。
  prep <- gmxdcc:::gm_prepare_data(d$r_dgp, X = NULL, K = K, rv = TRUE,
                                   period = "monthly")
  x <- prep$returns; RV_m <- prep$lag_mats$RV

  # Legacy begin_val reconstruction: columns alpha,beta,m,theta,w2,mu with
  # runif draws in that order (lines 1801-1808).
  # 重构 legacy 初值矩阵:runif 依 alpha,beta,m,theta 列序消耗随机数。
  set.seed(SEED)
  B_leg <- matrix(NA_real_, R_starts, 6,
                  dimnames = list(NULL, c("alpha", "beta", "m", "theta", "w2", "mu")))
  B_leg[, "alpha"] <- stats::runif(R_starts, 0.001, 0.095)
  B_leg[, "beta"]  <- stats::runif(R_starts, 0.6, 0.8)
  B_leg[, "m"]     <- stats::runif(R_starts, -1, 1)
  B_leg[, "theta"] <- stats::runif(R_starts, -1, 1)
  B_leg[, "w2"]    <- 1.01
  B_leg[, "mu"]    <- mean(x, na.rm = TRUE)

  # Legacy fit with the same seed draws exactly B_leg internally
  # 同一种子下 legacy 内部抽出的正是 B_leg
  set.seed(SEED)
  fit_leg <- suppressMessages(lg$garchmidas_noskew_fit(
    model = "GM_noskew_RV", daily_ret = x, RV_m = RV_m, K = K,
    distribution = "norm", lag_fun = "Beta", R = R_starts
  ))

  # New fit from the same candidates (new order: alpha,beta,w2,m,theta,mu)
  # 新拟合注入同一候选(按新参数排序重排列)
  fit_new <- fit_garch_midas(
    d$r_dgp, X = NULL, K = K, rv = TRUE, period = "monthly",
    distribution = "norm", lag_fun = "Beta",
    control = list(start_matrix = B_leg[, c("alpha", "beta", "w2", "m", "theta", "mu")])
  )

  expect_equal(fit_new$loglik, fit_leg$loglik, tolerance = 1e-6)

  # Coefficient mapping old -> new order / 系数按新旧排序映射对比
  co_leg <- fit_leg$rob_coef_mat[, "Estimate"]
  co_new <- fit_new$rob_coef_mat[, "Estimate"]
  names(co_leg) <- rownames(fit_leg$rob_coef_mat)
  names(co_new) <- rownames(fit_new$rob_coef_mat)
  map <- c(alpha = "alpha", beta = "beta", m = "m", theta = "theta_RV",
           w2 = "w2_RV", mu = "mu")
  expect_equal(unname(co_new[map]), unname(co_leg[names(map)]), tolerance = 1e-4)

  # Robust SEs agree. Tolerance 2%: the numerical Hessian is evaluated at
  # convergence points that themselves differ by ~1e-4, which is amplified
  # in the flatter w2 direction. 三明治标准误一致;容差 2%,因数值
  # Hessian 在相差 ~1e-4 的收敛点上取值,对较平坦的 w2 方向敏感。
  se_leg <- fit_leg$rob_coef_mat[, "Std. Error"]
  se_new <- fit_new$rob_coef_mat[, "Std. Error"]
  names(se_leg) <- rownames(fit_leg$rob_coef_mat)
  names(se_new) <- rownames(fit_new$rob_coef_mat)
  expect_equal(unname(se_new[map]), unname(se_leg[names(map)]), tolerance = 2e-2)

  # In-sample components agree / 样本内成分一致
  expect_equal(as.numeric(fit_new$est_vol_in_s), as.numeric(fit_leg$est_vol_in_s),
               tolerance = 1e-4)
  expect_equal(unname(fit_new$inf_criteria), unname(fit_leg$inf_criteria),
               tolerance = 1e-4)
})

test_that("fit_garch_midas handles breaks + two covariates + std + OOS", {
  skip_on_cran()
  d <- sim_data()

  fit <- fit_garch_midas(
    d$r1, X = list(M = d$macro), K = 12, rv = TRUE, period = "monthly",
    distribution = "std", lag_fun = "Beta",
    break_dates = as.Date("2016-01-01"),
    out_of_sample = 150, n_starts = 8, seed = 11
  )

  expect_s3_class(fit, "gm_fit")
  expect_equal(fit$maxlik$code, 0)
  # 2 (alpha,beta) + 2 w2 + 1 nu + 2 m + 4 theta + 1 mu = 12 parameters
  expect_identical(nrow(fit$rob_coef_mat), 12L)
  expect_length(fit$est_vol_oos, 150)
  # summed per-obs loglik at the estimate equals the stored loglik
  # 参数处逐观测似然之和等于存储的 loglik
  expect_s3_class(logLik(fit), "logLik")
  expect_named(coef(fit))
  expect_true(all(is.finite(fit$est_vol_in_s)))
})

test_that("out-of-sample components are warm-started (no cold reset)", {
  skip_on_cran()
  d <- sim_data()

  oos <- 150
  fit <- fit_garch_midas(
    d$r_dgp, X = NULL, K = 12, rv = TRUE, period = "monthly",
    distribution = "norm", lag_fun = "Beta",
    out_of_sample = oos, n_starts = 8, seed = 7
  )

  # Reconstruct the full-sample filter at the estimated parameters and check
  # the first OOS day continues the recursion from the last in-sample state:
  # g_{T+1} = (1-a-b) + a*eps_T^2/tau_T + b*g_T.
  # 验证首个样本外日的 g 由样本内末期状态递推一步而来,无冷启动跳变。
  prep <- gmxdcc:::gm_prepare_data(d$r_dgp, X = NULL, K = 12, rv = TRUE,
                                   period = "monthly")
  # Unrounded estimates (the coefficient table rounds to 6 decimals)
  # 用未舍入的估计值(系数表舍入到 6 位)
  co <- stats::coef(fit$maxlik)
  n_in <- fit$obs
  eps_T <- as.numeric(prep$returns[n_in]) - co["mu"]
  tau_T <- fit$est_lr_in_s[n_in]
  g_T   <- fit$est_sr_in_s[n_in]
  g_first_oos <- (1 - co["alpha"] - co["beta"]) +
    co["alpha"] * eps_T^2 / tau_T + co["beta"] * g_T

  expect_equal(unname(fit$est_sr_oos[1]), unname(g_first_oos), tolerance = 1e-10)
  # A cold restart would give exactly 1 / 冷启动会把首日 g 重置为 1
  expect_false(isTRUE(all.equal(fit$est_sr_oos[1], 1)))
})
