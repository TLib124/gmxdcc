# Gate 3: ported utility functions must reproduce the legacy versions.
# 关卡 3:工具函数与 legacy 版本在模拟数据上输出一致。

test_that("align_sample matches legacy cut_func", {
  skip_if(is.null(legacy_env()), "dev/legacy not available")
  lg <- legacy_env(); d <- sim_data()

  expect_identical(
    align_sample(d$r1, K = 12, K_c = 0, macro = d$macro, period = "monthly"),
    lg$cut_func(d$r1, K = 12, K_c = 0, macro = d$macro, type = "monthly")
  )
  expect_identical(
    align_sample(d$r1, K = 4, K_c = 2, macro = d$macro_q, period = "quarterly"),
    lg$cut_func(d$r1, K = 4, K_c = 2, macro = d$macro_q, type = "quarterly")
  )
  # NULL macro/period passthrough branch. The legacy cut_func errored here
  # (`if (NULL == "monthly")`); the reordered NULL check is a deliberate fix,
  # so only the new behavior is asserted.
  # macro/period 为 NULL 时原样返回;legacy 在此本身会报错,分支重排是
  # 有意修复,故只断言新函数行为。
  expect_identical(
    align_sample(d$r1, K = 12, macro = NULL, period = NULL),
    d$r1
  )
})

test_that("midas_lag_matrix matches legacy mv_into_mat2", {
  skip_if(is.null(legacy_env()), "dev/legacy not available")
  lg <- legacy_env(); d <- sim_data()

  x_m <- align_sample(d$r1, K = 12, macro = d$macro, period = "monthly")
  expect_identical(
    midas_lag_matrix(x_m, d$macro, K = 12, period = "monthly"),
    lg$mv_into_mat2(x_m, d$macro, K = 12, type = "monthly")
  )

  x_q <- align_sample(d$r1, K = 4, macro = d$macro_q, period = "quarterly")
  expect_identical(
    midas_lag_matrix(x_q, d$macro_q, K = 4, period = "quarterly"),
    lg$mv_into_mat2(x_q, d$macro_q, K = 4, type = "quarterly")
  )
})

test_that("realized_measure matches legacy RV_gen", {
  skip_if(is.null(legacy_env()), "dev/legacy not available")
  lg <- legacy_env(); d <- sim_data()

  for (p in c("weekly", "monthly", "quarterly", "yearly")) {
    expect_identical(realized_measure(d$r1, p), lg$RV_gen(d$r1, p), label = p)
  }
})

test_that("loss_avg matches legacy LF_f", {
  skip_if(is.null(legacy_env()), "dev/legacy not available")
  lg <- legacy_env()

  set.seed(1)
  v_est   <- exp(rnorm(500, 0, 0.3))
  v_proxy <- exp(rnorm(500, 0, 0.3))
  expect_identical(loss_avg(v_est, v_proxy), lg$LF_f(v_est, v_proxy))
})

test_that("qmle_se and inf_criteria match legacy on a maxLik fit", {
  skip_if(is.null(legacy_env()), "dev/legacy not available")
  lg <- legacy_env()

  # Tiny normal MLE with per-observation log-likelihood so that maxLik
  # stores gradientObs. 小型正态 MLE,逐观测似然保证 gradientObs 存在。
  set.seed(2)
  y  <- rnorm(400, 1, 2)
  ll <- function(par) stats::dnorm(y, par[1], exp(par[2]), log = TRUE)
  fit <- maxLik::maxLik(ll, start = c(mu = 0, logsd = 0), method = "BFGS")

  expect_identical(qmle_se(fit), lg$QMLE_sd(fit))
  expect_identical(inf_criteria(fit), lg$Inf_criteria(fit))

  # Coefficient table reproduces the legacy inline construction
  # 系数表与旧脚本内联构造一致
  # Parameter names live in rownames (as in the legacy tables); columns are
  # plain numeric. 参数名在行名上,列为纯数值,与旧脚本一致。
  tab <- qmle_coeftable(fit)
  expect_identical(rownames(tab), names(stats::coef(fit)))
  expect_identical(tab[, "Estimate"], unname(round(stats::coef(fit), 6)))
  expect_identical(tab[, "Std. Error"], unname(round(lg$QMLE_sd(fit), 6)))
})
