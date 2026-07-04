# Gate 7b: end-to-end fit_dcc — golden comparison against the legacy
# two-step fit (deterministic cDCC path) and structural/warm-start checks
# for the DCC-MIDAS path.
# 关卡 7b:fit_dcc 端到端——cDCC 金标准对照 + DCC-MIDAS 结构与热启动检验。

# Correlated 2-asset panel with GARCH dynamics in the first column.
# 两资产相关面板,第一列含真实 GARCH 动态。
dcc_panel <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    d <- sim_data()
    set.seed(505)
    r1 <- d$r_dgp
    z  <- xts::xts(rnorm(length(r1), 0, 0.8), order.by = zoo::index(r1))
    r2 <- 0.6 * r1 + z                     # correlated second asset
    panel <- xts::merge.xts(A1 = r1, A2 = r2)
    colnames(panel) <- c("A1", "A2")
    cache <<- list(panel = panel, macro = d$macro)
    cache
  }
})

test_that("fit_dcc with sGARCH + cDCC reproduces the legacy two-step fit", {
  skip_on_cran()
  skip_if(is.null(legacy_env()), "dev/legacy not available")
  lg <- legacy_env(); pn <- dcc_panel()

  # Both paths are fully deterministic here (rugarch fit + fixed cDCC start),
  # so no seed reconstruction is needed.
  # 该路径完全确定(rugarch + cDCC 固定初值),无需种子重放。
  fit_leg <- suppressWarnings(suppressMessages(lg$dccmidas_noskew_longfocus_fit(
    r_t = pn$panel, univ_model = "sGARCH", distribution = "norm",
    corr_model = "cDCC"
  )))
  fit_new <- suppressWarnings(fit_dcc(
    pn$panel, univariate = "sGARCH", correlation = "cDCC",
    distribution = "norm"
  ))

  expect_equal(as.numeric(fit_new$c_llk), as.numeric(fit_leg$c_llk),
               tolerance = 1e-8)
  expect_equal(fit_new$corr_coef_mat[, "Estimate"],
               fit_leg$corr_coef_mat[, "Estimate"], tolerance = 1e-6)
  expect_equal(fit_new$corr_coef_mat[, "Std. Error"],
               fit_leg$corr_coef_mat[, "Std. Error"], tolerance = 1e-4)
  expect_equal(unname(fit_new$eps_t), unname(zoo::coredata(fit_leg$eps_t)),
               tolerance = 1e-10)
  expect_equal(fit_new$H_t, fit_leg$H_t, tolerance = 1e-8)
  expect_equal(fit_new$R_t, fit_leg$R_t, tolerance = 1e-8)
  expect_equal(unname(fit_new$Inf_criteria_2step),
               unname(fit_leg$Inf_criteria_2step), tolerance = 1e-4)
})

test_that("fit_dcc: garch_midas + DCC-MIDAS RC+X with breaks and OOS", {
  skip_on_cran()
  pn <- dcc_panel()
  oos <- 120

  fit <- suppressWarnings(fit_dcc(
    pn$panel, univariate = "garch_midas",
    univariate_args = list(K = 12, rv = TRUE),
    correlation = "DCC-MIDAS", rc = TRUE, N_c = 60,
    corr_X = list(M = pn$macro), K_c = 6,
    corr_break_dates = as.Date("2016-01-01"),
    distribution = "norm", lag_fun = "Beta", period = "monthly",
    out_of_sample = oos, n_starts = 6, seed = 21
  ))

  expect_s3_class(fit, "dcc_fit")
  N <- 2; n_in <- fit$obs
  # parameters: a, b, w2_RC, w2_M_1 (P=1), m_1..2, theta_RC_1..2, theta_M_1..2
  expect_identical(nrow(fit$corr_coef_mat), 10L)
  expect_identical(dim(fit$H_t), as.integer(c(N, N, n_in)))
  expect_identical(dim(fit$H_t_oos), as.integer(c(N, N, oos)))
  expect_true(all(is.finite(fit$R_t_oos)))
  # Correlations bounded, unit diagonal / 相关阵对角为 1,元素有界
  expect_true(all(abs(fit$R_t_oos) <= 1 + 1e-8))
  expect_equal(fit$R_t_oos[1, 1, ], rep(1, oos), tolerance = 1e-12)

  # Warm start: the first OOS correlation must continue the in-sample
  # recursion — a legacy-style cold restart on the OOS segment alone gives
  # a different (identity-initialized) answer.
  # 热启动:首个 OOS 相关阵延续样本内递推;旧式冷启动结果不同。
  cold <- new_dcc_mat_est(c(coef(fit)["a"], coef(fit)["b"]),
                          fit$eps_t_oos, fit$D_t_oos, K_c = NULL)
  expect_false(isTRUE(all.equal(cold$R_t[, , 1], fit$R_t_oos[, , 1])))

  # S3 methods run / S3 方法可用
  expect_named(coef(fit))
  expect_s3_class(logLik(fit), "logLik")
  expect_output(print(fit), "DCC-MIDAS")
})

test_that("fit_dcc supports aDCC and DECO with std innovations", {
  skip_on_cran()
  pn <- dcc_panel()

  for (cm in c("aDCC", "DECO")) {
    fit <- suppressWarnings(fit_dcc(
      pn$panel, univariate = "sGARCH", correlation = cm,
      distribution = "std"
    ))
    expect_s3_class(fit, "dcc_fit")
    expect_true(all(is.finite(fit$corr_coef_mat[, "Estimate"])))
    expect_identical(nrow(fit$corr_coef_mat), if (cm == "aDCC") 3L else 2L)
  }
})

test_that("start_matrix routing: corr-step injection + per-asset univ starts", {
  skip_on_cran()
  pn <- dcc_panel()
  panel <- pn$panel["2014/2019"]

  # Reference fit with random multistart / 随机多起点的基准拟合
  base <- suppressWarnings(fit_dcc(
    panel, univariate = "garch_midas",
    univariate_args = list(X = pn$macro, K = 6, rv = FALSE),
    correlation = "DCC-MIDAS", rc = FALSE,
    corr_X = list(M = pn$macro), K_c = 6,
    n_starts = 2, seed = 3
  ))

  # Warm restart from the base optimum: per-asset univariate starts via
  # univariate_args, correlation start via control — the rolling
  # re-estimation pattern. Before the routing fix this errored (the corr
  # start matrix leaked into the univariate stage with wrong columns).
  # 从基准最优热重启:第一步逐资产、第二步经 control 注入,即滚动重估
  # 用法;路由修复前相关步起点会漏进第一步导致列数错误。
  univ_starts <- lapply(base$est_univ_model,
                        function(m) matrix(m[, "Estimate"], nrow = 1))
  corr_start  <- matrix(coef(base), nrow = 1)

  warm <- suppressWarnings(fit_dcc(
    panel, univariate = "garch_midas",
    univariate_args = list(X = pn$macro, K = 6, rv = FALSE,
                           start_matrix = univ_starts),
    correlation = "DCC-MIDAS", rc = FALSE,
    corr_X = list(M = pn$macro), K_c = 6,
    control = list(start_matrix = corr_start)
  ))

  # Restarting at the optimum must converge back to (essentially) it
  # 从最优点出发必须收敛回同一处
  expect_equal(coef(warm), coef(base), tolerance = 1e-3)
  expect_equal(as.numeric(warm$c_llk), as.numeric(base$c_llk), tolerance = 1e-4)

  # Wrong-length list is rejected / 列表长度不符时报错
  expect_error(
    fit_dcc(panel, univariate = "garch_midas",
            univariate_args = list(X = pn$macro, K = 6, rv = FALSE,
                                   start_matrix = univ_starts[1]),
            correlation = "cDCC"),
    "one matrix per asset"
  )
})
