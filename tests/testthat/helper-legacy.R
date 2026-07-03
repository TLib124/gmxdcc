# Load selected function definitions from the frozen legacy script without
# executing it (skips sourceCpp calls and any top-level side effects).
# 从冻结的 legacy 脚本中只取函数定义,不执行 sourceCpp 等顶层副作用。

legacy_script_path <- function() {
  normalizePath(
    file.path(testthat::test_path(), "..", "..", "dev", "legacy", "legacy_integrated_fit.R"),
    mustWork = FALSE
  )
}

# Parse the legacy file once and eval only `name <- function(...)` and
# `name <- pkg::obj` assignments into a fresh environment.
# 只 eval 赋值表达式(函数与包对象别名),缓存于会话内。
legacy_env <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    path <- legacy_script_path()
    if (!file.exists(path)) return(NULL)
    # The legacy script assumed library(xts) etc. were attached; replicate
    # that for its unqualified calls (apply.weekly, hessian, roll_sum, ...).
    # legacy 脚本依赖挂载的包,这里补挂以支持其未加限定的调用。
    for (p in c("xts", "zoo", "lubridate", "maxLik", "roll")) {
      suppressMessages(require(p, character.only = TRUE, quietly = TRUE))
    }
    exprs <- parse(path, encoding = "UTF-8")
    env <- new.env(parent = globalenv())
    for (e in exprs) {
      if (is.call(e) && identical(as.character(e[[1]]), "<-")) {
        rhs <- e[[3]]
        is_fun   <- is.call(rhs) && identical(as.character(rhs[[1]]), "function")
        is_alias <- is.call(rhs) && identical(as.character(rhs[[1]]), "::")
        if (is_fun || is_alias) eval(e, envir = env)
      }
    }
    cache <<- env
    env
  }
})

# Shared simulated dataset: ~8 years of daily returns for 2 assets plus a
# monthly macro covariate that starts well before the returns.
# 共享模拟数据:两资产约 8 年日收益 + 提前起始的月度宏观变量。
sim_data <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    set.seed(20260703)
    days   <- seq(as.Date("2012-01-02"), as.Date("2019-12-31"), by = "day")
    days   <- days[!(format(days, "%u") %in% c("6", "7"))]  # weekdays only
    n      <- length(days)
    r1     <- xts::xts(rnorm(n, 0, 1.1), order.by = days)
    r2     <- xts::xts(rnorm(n, 0, 0.9), order.by = days)
    months <- seq(as.Date("2010-01-01"), as.Date("2019-12-01"), by = "month")
    macro  <- xts::xts(as.numeric(scale(cumsum(rnorm(length(months), 0, 0.3)))),
                       order.by = months)
    quarters <- seq(as.Date("2010-01-01"), as.Date("2019-10-01"), by = "quarter")
    macro_q  <- xts::xts(as.numeric(scale(cumsum(rnorm(length(quarters), 0, 0.3)))),
                         order.by = quarters)

    # r_dgp: returns generated from an actual GARCH-MIDAS process so that
    # fit-level tests face a well-identified, curved likelihood surface
    # (white noise leaves theta/m flat and optimizer paths diverge).
    # r_dgp:按真实 GARCH-MIDAS 过程生成的收益,保证似然曲面弯曲、
    # 参数可识别,端到端拟合测试才有意义。
    cal   <- xts::xts(rep(0, n), order.by = days)
    cal   <- align_sample(cal, K = 12, macro = macro, period = "monthly")
    lagm  <- midas_lag_matrix(cal, macro, K = 12, period = "monthly")
    filt  <- suppressWarnings(
      roll::roll_sum(lagm, 13, weights = gmxdcc:::midas_filter_weights(12, 3, "Beta"))
    )[13, ]
    tau   <- exp(0.2 + 0.5 * filt)
    TTd   <- length(tau)
    a <- 0.06; b <- 0.85; mu <- 0.05
    z <- rnorm(TTd)
    g <- eps <- numeric(TTd)
    g[1] <- 1; eps[1] <- sqrt(g[1] * tau[1]) * z[1]
    for (t in 2:TTd) {
      g[t]   <- (1 - a - b) + a * eps[t - 1]^2 / tau[t - 1] + b * g[t - 1]
      eps[t] <- sqrt(g[t] * tau[t]) * z[t]
    }
    r_dgp <- xts::xts(mu + eps, order.by = zoo::index(cal))

    cache <<- list(r1 = r1, r2 = r2, macro = macro, macro_q = macro_q,
                   r_dgp = r_dgp)
    cache
  }
})
