# Gate 2: packaged C++ kernels must be numerically identical to the legacy
# sourceCpp() versions frozen under dev/legacy/.
# 关卡 2:包内 C++ 核心与 dev/legacy 冻结版逐元素一致(唯一改动是导出签名
# 的 arma:: 限定,数值路径未动)。
#
# Compiling the legacy files takes ~1 min, so this test is skipped on CRAN
# and when the legacy sources are absent (e.g. installed package).
# 编译 legacy 源码较慢,CRAN 或安装后的包环境下跳过。

legacy_dir <- normalizePath(
  file.path(testthat::test_path(), "..", "..", "dev", "legacy"),
  mustWork = FALSE
)

test_that("packaged C++ kernels match legacy sourceCpp versions", {
  skip_on_cran()
  skip_if(!dir.exists(legacy_dir), "dev/legacy not available")

  legacy <- new.env()
  legacy_files <- c(
    "garch_midas_ll_norm_turbo.cpp", "garch_midas_ll_std_turbo.cpp",
    "midas_ll_rcpp_turbo.cpp", "midas_rc_rcpp_turbo.cpp",
    "new_dcc_loglik_turbo.cpp", "new_dcc_mat_est_turbo.cpp",
    "new_a_dcc_loglik_turbo.cpp", "new_a_dcc_mat_est_turbo.cpp",
    "new_deco_loglik_turbo.cpp", "new_deco_mat_est_turbo.cpp"
  )
  for (f in legacy_files) {
    Rcpp::sourceCpp(file.path(legacy_dir, f), env = legacy)
  }

  set.seed(42)
  TT <- 300; N <- 3; K_c <- 5; N_c <- 20

  # Simulated standardized residuals and univariate inputs
  # 模拟标准化残差与单变量输入
  res       <- matrix(rnorm(TT * N), TT, N)
  daily_ret <- rnorm(TT, 0, 1.2)
  epsilon   <- daily_ret - 0.05
  tau_d     <- exp(0.5 + 0.2 * sin(seq_len(TT) / 20))

  # D_t cube of conditional sd diagonal matrices / 条件标准差对角阵序列
  Dt <- array(0, dim = c(N, N, TT))
  for (tt in seq_len(TT)) diag(Dt[, , tt]) <- abs(rnorm(N, 1, 0.1))

  # R_t_bar cube of long-run correlations / 长期相关阵序列
  R_bar <- array(0, dim = c(N, N, TT))
  for (tt in seq_len(TT)) {
    A <- matrix(rnorm(N * N, 0, 0.3), N, N)
    S <- crossprod(A) + diag(N)
    R_bar[, , tt] <- cov2cor(S)
  }

  ## --- univariate GARCH-MIDAS likelihood kernels ---
  expect_identical(
    fast_garch_ll_cpp(0.05, 0.9, 0.05, epsilon, tau_d, daily_ret),
    legacy$fast_garch_ll_cpp(0.05, 0.9, 0.05, epsilon, tau_d, daily_ret)
  )
  expect_identical(
    fast_t_garch_ll(0.05, 0.9, 8, epsilon, tau_d, daily_ret),
    legacy$fast_t_garch_ll(0.05, 0.9, 8, epsilon, tau_d, daily_ret)
  )

  ## --- DCC-family likelihood / matrix kernels ---
  expect_identical(
    new_dcc_loglik(c(0.03, 0.9), res, K_c),
    legacy$new_dcc_loglik(c(0.03, 0.9), res, K_c)
  )
  expect_identical(
    new_a_dcc_loglik(c(0.03, 0.9, 0.02), res, K_c),
    legacy$new_a_dcc_loglik(c(0.03, 0.9, 0.02), res, K_c)
  )
  expect_identical(
    new_deco_loglik(c(0.03, 0.9), res, K_c),
    legacy$new_deco_loglik(c(0.03, 0.9), res, K_c)
  )
  expect_identical(
    new_dcc_mat_est(c(0.03, 0.9), res, Dt, K_c),
    legacy$new_dcc_mat_est(c(0.03, 0.9), res, Dt, K_c)
  )
  expect_identical(
    new_a_dcc_mat_est(c(0.03, 0.9, 0.02), res, Dt, K_c),
    legacy$new_a_dcc_mat_est(c(0.03, 0.9, 0.02), res, Dt, K_c)
  )
  expect_identical(
    new_deco_mat_est(c(0.03, 0.9), res, Dt, K_c),
    legacy$new_deco_mat_est(c(0.03, 0.9), res, Dt, K_c)
  )

  ## --- DCC-MIDAS kernels ---
  expect_identical(
    midas_ll_rcpp(TT, N, K_c, 0.03, 0.9, R_bar, res),
    legacy$midas_ll_rcpp(TT, N, K_c, 0.03, 0.9, R_bar, res)
  )
  expect_identical(
    midas_rc_rcpp(N_c, N, TT, res),
    legacy$midas_rc_rcpp(N_c, N, TT, res)
  )
})
