#' Summary of a GARCH-MIDAS fit
#'
#' Prints the robust coefficient table with significance codes, sample
#' information and in-sample losses, in the same layout as the legacy
#' `summary.gmnoskew()`.
#'
#' @param object A `gm_fit` object.
#' @param ... Unused.
#' @return `object`, invisibly.
#' @export
summary.gm_fit <- function(object, ...) {

  mat_coef <- object$rob_coef_mat
  Obs      <- object$obs
  Period   <- object$period
  loss     <- object$loss_in_s

  Period <- paste(substr(Period[1], 1, 10), "/",
                  substr(Period[2], 1, 10), sep = "")

  p_value <- mat_coef[, 4]
  sig <- ifelse(p_value <= 0.01, "***",
                ifelse(p_value > 0.01 & p_value <= 0.05, "**",
                       ifelse(p_value > 0.05 & p_value <= 0.1, "*", " ")))

  mat_coef <- round(mat_coef, 4)
  mat_coef <- cbind(mat_coef, Sig. = sig)

  cat("\nModel:", object$model, "\n\n")
  cat("Coefficients:\n")
  cat(utils::capture.output(mat_coef), sep = "\n")
  cat("--- \n")
  cat("Signif. codes: 0.01 '***', 0.05 '**', 0.1 '*' \n\n")
  cat("Obs.:", paste(Obs, ".", sep = ""), "Sample Period:", Period, "\n")
  cat("MSE(%):", paste(as.numeric(round(loss[1], 6)), "; ", sep = ""),
      "QLIKE:", as.numeric(round(loss[2], 6)), "\n")
  if (!is.null(object$loss_oos)) {
    cat("Out-of-sample MSE(%):",
        paste(as.numeric(round(object$loss_oos[1], 6)), "; ", sep = ""),
        "QLIKE:", as.numeric(round(object$loss_oos[2], 6)), "\n")
  }
  cat("\n")
  invisible(object)
}

#' @describeIn summary.gm_fit Compact coefficient printout.
#' @param x A `gm_fit` object.
#' @export
print.gm_fit <- function(x, ...) {
  mat_coef <- x$rob_coef_mat

  coef_char <- as.character(round(mat_coef[, 1], 4))
  row_names <- gsub("\\s", " ", format(rownames(mat_coef), width = 9))
  coef_val  <- gsub("\\s", " ", format(coef_char, width = 9))

  cat("\nModel:", x$model, "\n\n")
  cat("Coefficients: \n")
  cat(row_names, sep = " ", "\n")
  cat(coef_val, sep = " ", "\n")
  cat("\n")
  invisible(x)
}

#' @describeIn summary.gm_fit Log-likelihood with df/nobs attributes.
#' @export
logLik.gm_fit <- function(object, ...) {
  val <- object$loglik
  # df = number of estimated parameters; nobs feeds BIC
  # df 为参数个数,nobs 供 BIC 计算
  attr(val, "df")   <- nrow(object$rob_coef_mat)
  attr(val, "nobs") <- object$obs
  class(val) <- "logLik"
  val
}

#' @describeIn summary.gm_fit Estimated coefficients.
#' @method coef gm_fit
#' @export
coef.gm_fit <- function(object, ...) {
  stats::setNames(object$rob_coef_mat[, "Estimate"],
                  rownames(object$rob_coef_mat))
}

#' @describeIn summary.gm_fit Robust (sandwich) variance-covariance matrix.
#' @method vcov gm_fit
#' @export
vcov.gm_fit <- function(object, ...) {
  est <- object$maxlik
  H   <- -maxLik::hessian(est)
  H_inv <- tryCatch(solve(H), error = function(e) {
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Please install 'MASS' package.")
    MASS::ginv(H)
  })
  OPG <- t(est$gradientObs) %*% est$gradientObs
  V   <- H_inv %*% OPG %*% H_inv
  dimnames(V) <- list(names(stats::coef(est)), names(stats::coef(est)))
  V
}
