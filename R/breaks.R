#' Step-dummy matrix for structural breaks
#'
#' Builds the `(1 + n_breaks) x TT` step-dummy matrix that turns the long-run
#' intercept and slopes into piecewise-constant series. Row 1 is the constant
#' regime; row `i + 1` switches on from break date `i` onwards. With
#' `break_dates = NULL` the matrix collapses to a single row of ones, so the
#' no-break model is simply the `n_breaks = 0` special case.
#' 结构突变阶梯虚拟矩阵:第 1 行为常数基准,第 i+1 行从第 i 个突变日起为 1;
#' `break_dates = NULL` 时退化为单行全 1,无突变即 n_breaks = 0 的特例。
#'
#' @param dates Vector of (daily) observation dates, coercible to `Date`.
#' @param break_dates Vector of break dates, or `NULL` for no breaks.
#'
#' @return Numeric matrix with `1 + length(break_dates)` rows and
#'   `length(dates)` columns.
#' @examples
#' days <- seq(as.Date("2019-01-01"), as.Date("2019-01-10"), by = "day")
#' build_break_dummies(days, as.Date("2019-01-06"))
#' build_break_dummies(days, NULL)  # no breaks -> single row of ones
#' @export
build_break_dummies <- function(dates, break_dates = NULL) {
  full_days <- as.Date(dates)
  TT  <- length(full_days)
  n_n <- length(break_dates)

  D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
  D_mat[1, ] <- 1  # constant regime / 常数基准段

  if (n_n > 0) {
    bd <- as.Date(break_dates)
    if (is.unsorted(bd)) stop("break_dates must be in increasing order.", call. = FALSE)
    if (any(bd <= full_days[1]) || any(bd > full_days[TT])) {
      stop("Every break date must lie strictly inside the sample period.", call. = FALSE)
    }
    for (i in 1:n_n) {
      # Days at or after the i-th break date get a 1
      # 第 i 个突变日及之后置 1
      D_mat[i + 1, full_days >= bd[i]] <- 1
    }
  }

  D_mat
}
