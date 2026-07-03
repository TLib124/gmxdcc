`%||%` <- function(a, b) if (is.null(a)) b else a

#' Multi-start constrained maximum-likelihood driver
#'
#' Reproduces the legacy estimation scheme shared by every fit branch:
#' draw a matrix of candidate starts, evaluate the summed log-likelihood at
#' each, and run a single constrained BFGS `maxLik` from the best candidate.
#' New relative to the legacy script: an optional `seed` makes the draw
#' reproducible, and `control$start_matrix` injects a fixed start matrix
#' (used by the golden regression tests).
#' 多起点估计驱动:生成候选初值、按似然择优、单次带约束 BFGS;新增
#' seed(旧脚本无种子不可复现)与 start_matrix 注入口(供回归测试)。
#'
#' @param loglik_fn Per-observation log-likelihood, called as
#'   `loglik_fn(param, data)`.
#' @param data Data bundle forwarded to `loglik_fn`.
#' @param spec Parameter specification from [gm_param_spec()] (or the DCC
#'   analogue): supplies names, start ranges and `ui`/`ci` constraints.
#' @param n_starts Number of random starts (legacy default 100).
#' @param seed Optional RNG seed for the start draw.
#' @param control List; recognized fields: `start_matrix` (fixed candidate
#'   matrix in canonical parameter order), `iterlim` (default 1000),
#'   `method` (default `"BFGS"`).
#'
#' @return List with `est` (the `maxLik` object), `start_matrix`, and
#'   `start_used` (the selected start vector).
#' @keywords internal
run_multistart <- function(loglik_fn, data, spec, n_starts = 100,
                           seed = NULL, control = list()) {

  begin_val <- control$start_matrix
  if (is.null(begin_val)) {
    if (!is.null(seed)) set.seed(seed)
    begin_val <- draw_starts(spec, n_starts)
  } else {
    begin_val <- as.matrix(begin_val)
    if (ncol(begin_val) != spec$n_par) {
      stop(sprintf("control$start_matrix must have %d columns (canonical parameter order).",
                   spec$n_par), call. = FALSE)
    }
    colnames(begin_val) <- spec$names
  }

  # Rank candidates by summed log-likelihood, exactly as the legacy loop.
  # 按总似然为候选初值排序,与旧脚本一致。
  which_row <- rep(NA_real_, nrow(begin_val))
  for (i in seq_len(nrow(begin_val))) {
    which_row[i] <- sum(loglik_fn(begin_val[i, ], data))
  }

  start_val <- begin_val[which.max(which_row), ]
  names(start_val) <- spec$names

  # Single constrained BFGS run from the best start (legacy behavior).
  # 由最优候选出发的单次带约束 BFGS(与旧脚本一致)。
  est <- suppressWarnings(maxLik::maxLik(
    logLik = loglik_fn,
    start  = start_val,
    data   = data,
    constraints = list(ineqA = spec$ui, ineqB = spec$ci),
    iterlim = control$iterlim %||% 1000,
    method  = control$method %||% "BFGS"
  ))

  list(est = est, start_matrix = begin_val, start_used = start_val)
}
