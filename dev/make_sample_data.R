# One-time script: build the packaged sample datasets from test20260704.RData.
# 一次性脚本:从 test20260704.RData 提取并整理随包发布的示例数据。
# The original xts objects carry POSIXct indices and quirky internal
# structure (column subsetting drops the index), so every dataset is
# rebuilt as a fresh, clean xts with a Date index.
# 原始对象索引为 POSIXct 且内部结构不规范(切列会丢索引),
# 故全部用 coredata + Date 索引重新构造。
# Run from the package root: Rscript dev/make_sample_data.R

suppressMessages({ library(xts); library(zoo) })

e <- new.env()
load("test20260704.RData", envir = e)

clean_col <- function(x, j, new_name) {
  out <- xts::xts(matrix(as.numeric(zoo::coredata(x)[, j]), ncol = 1),
                  order.by = as.Date(zoo::index(x)))
  colnames(out) <- new_name
  out
}

# Daily TOPIX market returns / 日频 TOPIX 市场收益
topix_market <- clean_col(e$rate_topix_xts, 1, "TOPIX")

# Daily TOPIX-17 Autos & Transportation Equipment industry returns
# 日频 TOPIX-17 汽车与运输设备行业收益
topix_autos <- clean_col(e$rate_topix17_xts, 1, "AUTOS")

# Monthly Japan macro uncertainty index (single series)
# 月频日本宏观不确定性指数(单列)
jp_macro_uncertainty <- clean_col(
  e$post_standardized_uncertainty_JP_xts,
  which(colnames(e$post_standardized_uncertainty_JP_xts) == "Japan Macro Uncertainty"),
  "JPMU"
)

# Structural break dates detected in the uncertainty series
# (strucchange::breakpoints on the JPMU mean)
# JPMU 均值的结构突变日期(strucchange::breakpoints 检测)
jp_mu_breaks <- as.Date(e$jp_sb_break_time)

stopifnot(
  inherits(topix_market, "xts"), inherits(topix_autos, "xts"),
  inherits(jp_macro_uncertainty, "xts"),
  inherits(zoo::index(topix_autos), "Date"),
  inherits(zoo::index(jp_macro_uncertainty), "Date")
)

dir.create("data", showWarnings = FALSE)
save(topix_market,         file = "data/topix_market.rda",         compress = "xz")
save(topix_autos,          file = "data/topix_autos.rda",          compress = "xz")
save(jp_macro_uncertainty, file = "data/jp_macro_uncertainty.rda", compress = "xz")
save(jp_mu_breaks,         file = "data/jp_mu_breaks.rda",         compress = "xz")

cat("Saved. Sizes (KB):\n")
print(round(file.size(list.files("data", full.names = TRUE)) / 1024, 1))
