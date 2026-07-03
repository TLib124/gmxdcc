#' Summary of a DCC(-MIDAS) fit
#'
#' Prints the univariate coefficient tables per asset, the second-step
#' correlation table with significance codes, and the two-step information
#' criteria, following the layout of the legacy `summary.dmnoskew()`.
#'
#' @param object A `dcc_fit` object.
#' @param ... Unused.
#' @return `object`, invisibly.
#' @export
summary.dcc_fit <- function(object, ...) {

  Period <- paste(substr(object$period[1], 1, 10), "/",
                  substr(object$period[2], 1, 10), sep = "")

  cat("\nUnivariate model:", object$model,
      "  |  Correlation model:", object$mult_model, "\n")
  cat("Obs.:", paste(object$obs, ".", sep = ""), "Sample Period:", Period, "\n\n")

  for (i in seq_along(object$est_univ_model)) {
    cat("--- Asset:", object$assets[i], "---\n")
    print(round(as.data.frame(object$est_univ_model[[i]]), 4))
    cat("\n")
  }

  mat_coef <- object$corr_coef_mat
  p_value  <- mat_coef[, 4]
  sig <- ifelse(p_value <= 0.01, "***",
                ifelse(p_value > 0.01 & p_value <= 0.05, "**",
                       ifelse(p_value > 0.05 & p_value <= 0.1, "*", " ")))
  mat_coef <- cbind(round(mat_coef, 4), Sig. = sig)

  cat("--- Correlation step (robust QMLE standard errors) ---\n")
  cat(utils::capture.output(mat_coef), sep = "\n")
  cat("--- \n")
  cat("Signif. codes: 0.01 '***', 0.05 '**', 0.1 '*' \n\n")
  cat("Two-step information criteria:\n")
  print(object$Inf_criteria_2step)
  cat("\n")
  invisible(object)
}

#' @describeIn summary.dcc_fit Compact printout of the correlation step.
#' @param x A `dcc_fit` object.
#' @export
print.dcc_fit <- function(x, ...) {
  cat("\nUnivariate model:", x$model,
      "  |  Correlation model:", x$mult_model, "\n\n")
  cat("Correlation coefficients:\n")
  print(round(x$corr_coef_mat[, "Estimate", drop = FALSE], 4))
  cat("\n")
  invisible(x)
}

#' @describeIn summary.dcc_fit Total two-step log-likelihood.
#' @export
logLik.dcc_fit <- function(object, ...) {
  ic  <- object$Inf_criteria_2step
  val <- as.numeric(ic["LogLik"])
  attr(val, "df")   <- as.numeric(ic["k"])
  attr(val, "nobs") <- object$obs
  class(val) <- "logLik"
  val
}

#' @describeIn summary.dcc_fit Second-step (correlation) coefficients.
#' @method coef dcc_fit
#' @export
coef.dcc_fit <- function(object, ...) {
  stats::setNames(object$corr_coef_mat[, "Estimate"],
                  rownames(object$corr_coef_mat))
}
