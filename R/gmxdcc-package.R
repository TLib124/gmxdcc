#' @keywords internal
#' @aliases gmxdcc-package
"_PACKAGE"

## usethis namespace: start
#' @useDynLib gmxdcc, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats coef logLik pnorm runif cov
#' @importFrom utils tail
## Formal imports from xts/zoo ensure both namespaces load together with
## gmxdcc, so their S3 methods (index.xts etc.) are registered even when
## the user has not attached xts — required for the packaged datasets.
## 正式 import 保证 xts/zoo 命名空间随包加载、S3 方法注册,
## 否则用户未 attach xts 时 zoo::index() 对 xts 对象不会正确分派。
#' @importFrom xts as.xts first last
#' @importFrom zoo index coredata
## usethis namespace: end
NULL
