#' Assert that an argument is an xts object
#'
#' Centralized input check replacing the legacy `stop(cat(...))` pattern
#' (which raised an empty error after printing).
#' 统一的 xts 类型校验,取代旧脚本的 stop(cat(...)) 反模式。
#'
#' @param x Object to check.
#' @param arg Argument name used in the error message.
#' @param allow_null Should `NULL` pass the check?
#' @return Invisibly `TRUE`; raises an error otherwise.
#' @keywords internal
assert_xts <- function(x, arg = deparse(substitute(x)), allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) return(invisible(TRUE))
    stop(sprintf("Parameter '%s' must not be NULL.", arg), call. = FALSE)
  }
  if (!inherits(x, "xts")) {
    stop(sprintf("Parameter '%s' must be an xts object. Please provide it in the correct form.", arg),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Align daily returns to the coverage of a low-frequency series
#'
#' Trims a daily return series so that every retained day has at least
#' `max(K, K_c) + 1` complete low-frequency periods of covariate history
#' available before it. Direct port of the legacy `cut_func()`.
#'
#' @param r_t Daily xts series (returns).
#' @param K MIDAS lag order of the volatility covariates (may be `NULL`).
#' @param K_c MIDAS lag order of the correlation covariates (may be `NULL`).
#' @param macro Low-frequency xts covariate whose calendar defines coverage.
#' @param period `"monthly"` or `"quarterly"`.
#'
#' @return The trimmed daily xts series.
#' @examples
#' # Keep only days with 12 full months of macro history available
#' # 只保留具备 12 个月完整宏观历史的交易日
#' x <- align_sample(topix_market, K = 12, macro = jp_macro_uncertainty,
#'                   period = "monthly")
#' range(zoo::index(x))
#' @export
align_sample <- function(r_t, K = NULL, K_c = NULL, macro = NULL, period = NULL) {

  if (is.null(K))   K   <- 0
  if (is.null(K_c)) K_c <- 0

  #### available for monthly/quarterly / 支持月度与季度
  if (is.null(macro) || is.null(period)) {

    required_r_t <- r_t

  } else if (period == "monthly") {

    # Shift the macro calendar forward by the lag window and keep only the
    # daily observations whose year-month falls inside the shifted range.
    # 宏观日历前移滞后窗口,只保留年月落在其中的交易日。
    cut_size      <- max(K, K_c) + 1
    year_month    <- unique(format(zoo::index(macro), "%Y-%m"))
    dates         <- as.Date(paste0(year_month, "-01"))
    shifted_dates <- seq(dates[1], by = "month", length.out = length(dates) + cut_size)[-c(1:cut_size)]
    shifted_dates <- shifted_dates[-((length(shifted_dates) - cut_size + 1):length(shifted_dates))]
    year_month    <- format(shifted_dates, "%Y-%m")
    required_r_t  <- r_t[format(zoo::index(r_t), "%Y-%m") %in% year_month]

  } else if (period == "quarterly") {

    cut_size      <- max(K, K_c) + 1
    year_quarter  <- unique(format(zoo::index(macro), "%Y-%m"))
    dates         <- as.Date(paste0(year_quarter, "-01"))
    shifted_dates <- seq(dates[1], by = "quarter", length.out = length(dates) + cut_size)[-c(1:(cut_size))]
    shifted_dates <- shifted_dates[-((length(shifted_dates) - cut_size):length(shifted_dates))]
    year_quarter_month_series <- seq(from = shifted_dates[1], to = shifted_dates[length(shifted_dates)], by = "month")
    year_quarter_month <- format(year_quarter_month_series, "%Y-%m")

    required_r_t <- r_t[format(zoo::index(r_t), "%Y-%m") %in% year_quarter_month]

  } else {

    required_r_t <- r_t

  }

  required_r_t
}

#' Build the MIDAS lag matrix of a low-frequency covariate
#'
#' For every high-frequency observation in `x`, stacks the current and `K`
#' previous low-frequency values of `mv` into a `(K+1) x N` matrix, where
#' `N = length(x)`. Direct port of the legacy `mv_into_mat2()`.
#'
#' @param x Daily xts series that defines the target calendar.
#' @param mv Low-frequency xts covariate.
#' @param K MIDAS lag order (number of low-frequency lags).
#' @param period One of `"weekly"`, `"monthly"`, `"quarterly"`, `"yearly"`.
#'
#' @return A numeric `(K+1) x N` matrix of stacked covariate lags.
#' @examples
#' x <- align_sample(topix_market, K = 12, macro = jp_macro_uncertainty,
#'                   period = "monthly")
#' M <- midas_lag_matrix(x, jp_macro_uncertainty, K = 12, period = "monthly")
#' dim(M)  # (K+1) x length(x)
#' @export
midas_lag_matrix <- function(x, mv, K, period) {

  N    <- length(x)
  N_mv <- length(mv)

  first_obs_r_t <- xts::first(zoo::index(x))
  first_obs_mv  <- xts::first(zoo::index(mv))
  diff_time     <- as.numeric(difftime(first_obs_r_t, first_obs_mv, units = "day"))
  # Rough days-per-period used only for the history-length check
  # 每期近似天数,仅用于历史长度检查
  num_obs_k <- ifelse(period == "weekly", 5,
                      ifelse(period == "monthly", 30,
                             ifelse(period == "quarterly", 120, 365)))

  w_days_r_t <- lubridate::wday(zoo::index(xts::last(x)))

  last_mv  <- paste(xts::last(lubridate::year(mv)),  xts::last(lubridate::month(mv)),  "20", sep = "-")
  last_r_t <- paste(xts::last(lubridate::year(x)),   xts::last(lubridate::month(x)),   "20", sep = "-")

  ############# checks / 输入检查
  assert_xts(x)
  assert_xts(mv)

  if (diff_time < 0 || ((diff_time - K * num_obs_k) < 0)) {
    stop("mv must start at least 'K+1' periods before x. Please decrease K or provide more observations (in the past) for mv.",
         call. = FALSE)
  }

  if (last_mv != last_r_t && period == "weekly") stop("x and mv must end in the same period", call. = FALSE)
  if (w_days_r_t != 6 && period == "weekly")     stop("the last day of x must be friday", call. = FALSE)

  mv_m <- matrix(c(rep(NA, (K + 1) * N)), ncol = N)

  if (period == "weekly") {
    for (i in 1:N) {
      start_p_mv1 <- which(format(zoo::index(mv), "%Y-%m-%d") ==
                             format(
                               lubridate::floor_date(zoo::index(x), "week", week_start = 1)[i], "%Y-%m-%d"))
      mv_m[, i]   <- mv[(start_p_mv1 - K):(start_p_mv1 + 1)][2:(K + 2)]
    }

  } else if (period == "monthly") {
    for (i in 1:N) {
      start_p_mv1 <- which(format(zoo::index(mv), "%Y-%m") == format(zoo::index(x)[i], "%Y-%m"))
      mv_m[, i]   <- mv[(start_p_mv1 - K):(start_p_mv1)]
    }
  } else if (period == "quarterly") {
    for (i in 1:N) {
      start_p_mv1 <- which(format(zoo::index(mv), "%Y-%m-%d") ==
                             format(
                               lubridate::floor_date(zoo::index(x), "quarter", week_start = 1)[i], "%Y-%m-%d"))
      mv_m[, i]   <- mv[(start_p_mv1 - K):(start_p_mv1 + 1)][2:(K + 2)]
    }

  } else {
    for (i in 1:N) {
      start_p_mv1 <- which(format(zoo::index(mv), "%Y-%m-%d") ==
                             format(
                               lubridate::floor_date(zoo::index(x), "year", week_start = 1)[i], "%Y-%m-%d"))
      mv_m[, i]   <- c(zoo::coredata(mv[(start_p_mv1 - K):(start_p_mv1)][2:(K + 1)]), 0)
    }
  }

  mv_m
}

#' Realized volatility at a low frequency
#'
#' Aggregates squared daily returns into weekly/monthly/quarterly/yearly
#' realized variance, with the period index aligned to the first calendar day
#' of each period. Direct port of the legacy `RV_gen()`.
#'
#' @param x Daily xts return series.
#' @param period One of `"weekly"`, `"monthly"`, `"quarterly"`, `"yearly"`.
#'
#' @return A low-frequency xts series of realized variances.
#' @examples
#' rv <- realized_measure(topix_market, "monthly")
#' head(rv)
#' @export
realized_measure <- function(x, period) {

  ############# checks / 输入检查
  assert_xts(x)

  calc_vol <- function(x) sum(x^2)

  if (period == "weekly") {

    weekly_vol <- xts::apply.weekly(x, FUN = calc_vol)
    # Align each week to its Monday (cut(..., "week") starts on Monday)
    # 将每周对齐到周一(cut 默认以周一为起始)
    zoo::index(weekly_vol) <- as.Date(cut(zoo::index(weekly_vol), "week"))
    colnames(weekly_vol) <- "Weekly_RV"

    vol <- weekly_vol

  } else if (period == "monthly") {

    monthly_vol <- xts::apply.monthly(x, FUN = calc_vol)
    zoo::index(monthly_vol) <- as.Date(format(zoo::index(monthly_vol), "%Y-%m-01"))
    colnames(monthly_vol) <- "Monthly_RV"

    vol <- monthly_vol

  } else if (period == "quarterly") {

    quarterly_vol <- xts::apply.quarterly(x, FUN = calc_vol)
    # zoo::as.Date on a yearqtr returns the first day of the quarter; the
    # zoo:: prefix is required because zoo ships its own as.Date generic.
    # zoo::as.Date 作用于 yearqtr 得到季度首日;必须加 zoo:: 前缀,
    # 因为 zoo 自带同名泛型,不挂载时 base::as.Date 无法分派。
    zoo::index(quarterly_vol) <- zoo::as.Date(zoo::as.yearqtr(zoo::index(quarterly_vol)))
    colnames(quarterly_vol) <- "Quarterly_RV"

    vol <- quarterly_vol

  } else {

    yearly_vol <- xts::apply.yearly(x, FUN = calc_vol)
    # Force "YYYY-01-01" formatting / 强制对齐到每年 1 月 1 日
    zoo::index(yearly_vol) <- as.Date(format(zoo::index(yearly_vol), "%Y-01-01"))
    colnames(yearly_vol) <- "Yearly_RV"

    vol <- yearly_vol
  }

  vol
}
