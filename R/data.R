#' Daily TOPIX market returns
#'
#' Daily returns of the TOPIX index (Tokyo Stock Exchange), used as the
#' market portfolio in the package examples.
#' 日频 TOPIX 指数收益,示例中作为市场组合。
#'
#' @format An `xts` object with one column `TOPIX`.
#' @source Tokyo Stock Exchange, via Datastream.
#' @seealso [topix_autos], [jp_macro_uncertainty], [jp_mu_breaks]
"topix_market"

#' Daily TOPIX-17 Autos & Transportation Equipment returns
#'
#' Daily returns of the TOPIX-17 "Autos & Transportation Equipment" industry
#' index, used as the example industry asset.
#' 日频 TOPIX-17 汽车与运输设备行业指数收益,示例中作为行业资产。
#'
#' @format An `xts` object with one column `AUTOS`.
#' @source Tokyo Stock Exchange, via Datastream.
#' @seealso [topix_market]
"topix_autos"

#' Japan macro uncertainty index (monthly)
#'
#' Monthly standardized Japanese macroeconomic uncertainty index, the
#' low-frequency covariate used in the package examples.
#' 月频标准化日本宏观不确定性指数,示例中的低频宏观协变量。
#'
#' @format An `xts` object with one column `JPMU`, monthly observations
#'   dated at the first day of each month.
#' @seealso [jp_mu_breaks]
"jp_macro_uncertainty"

#' Structural break dates of the Japan macro uncertainty index
#'
#' Break dates in the mean of [jp_macro_uncertainty], detected with
#' `strucchange::breakpoints(y ~ 1)`. Suitable as `break_dates` /
#' `corr_break_dates` arguments in [fit_garch_midas()] and [fit_dcc()].
#' JPMU 均值的结构突变日期(strucchange 检测),可直接作为
#' break_dates / corr_break_dates 参数使用。
#'
#' @format A `Date` vector of length 4.
"jp_mu_breaks"
