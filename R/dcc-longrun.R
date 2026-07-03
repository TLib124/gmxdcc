#' Symmetrize a vector of pair coefficients into an N x N matrix
#'
#' Fills the lower triangle (column-major pair order, as in the legacy
#' `lower.tri()` construction) and mirrors it; the diagonal stays zero.
#' 将逐对参数向量填充为对称矩阵(对角线为 0),与旧脚本 lower.tri 顺序一致。
#'
#' @param v Numeric vector of length `N (N - 1) / 2`.
#' @param n_assets Matrix dimension `N`.
#' @return Symmetric `N x N` matrix with zero diagonal.
#' @keywords internal
pair_matrix <- function(v, n_assets) {
  M <- matrix(0, n_assets, n_assets)
  M[lower.tri(M)] <- v
  M[upper.tri(M)] <- t(M)[upper.tri(M)]
  M
}

#' Long-run correlation component of the DCC-MIDAS model
#'
#' Single builder of the slowly-moving correlation target `R_t_bar`,
#' replacing the six legacy `dccmidas_noskew_loglik`/`mat_est` branches. The
#' functional forms follow the legacy code exactly:
#' * realized-correlation only, no breaks: `R_bar` is the MIDAS-weighted
#'   rolling average of the realized correlation matrices (no `tanh`);
#' * otherwise `R_bar = tanh(z)` with
#'   `z = level + sum_j slope_j * MIDAS(X_j) (+ slope_rc * RC)`, where the
#'   level/slopes are pair-specific without breaks and regime-varying
#'   (common across pairs) with breaks; `w2` weights are always
#'   pair-specific.
#' 长期相关构造器:纯 RC 无突变为加权滚动平均(无 tanh);其余为
#' tanh(截距 + 各斜率×MIDAS 项),无突变时逐对、有突变时按区制分段。
#'
#' @param param Parameter vector in [dcc_par_index()] order.
#' @param dc Data bundle: `res` (TT x N residual matrix), `n_assets`, `TT`,
#'   `K_c`, `N_c`, `rc`, `X_mats` (list of `(K_c+1) x TT` lag matrices),
#'   `D_mat`, `lag_fun`, `idx`.
#'
#' @return `N x N x TT` array `R_t_bar` (diagonal handling is left to the
#'   caller, as in the legacy code).
#' @keywords internal
build_longrun_corr <- function(param, dc) {
  N  <- dc$n_assets; TT <- dc$TT; K_c <- dc$K_c
  R  <- nrow(dc$D_mat); J <- length(dc$X_mats)
  idx <- dc$idx
  w1 <- ifelse(dc$lag_fun == "Beta", 1, 0)

  # All N^2 element positions, exactly as the legacy matrix_id_2 loop
  # (diagonal included; its NaN/placeholder values are overwritten later).
  # 全部 N^2 元素位置,与旧脚本一致(对角线值随后会被覆盖)。
  matrix_id   <- matrix(1:N^2, ncol = N)
  matrix_id_2 <- which(matrix_id == 1:N^2, arr.ind = TRUE)

  ## --- realized-correlation term --- 已实现相关项
  if (dc$rc) {
    C_t <- midas_rc_rcpp(dc$N_c, N, TT, dc$res)
    rc_betas <- c(rev(midas_weights(1:(K_c + 1), (K_c + 1), w1,
                                    param[idx$w2_rc], dc$lag_fun))[2:(K_c + 1)], 0)
    RC_tt <- array(diag(N), dim = c(N, N, TT))
    for (i in 1:nrow(matrix_id_2)) {
      RC_tt[matrix_id_2[i, 1], matrix_id_2[i, 2], ] <- suppressWarnings(
        roll::roll_sum(C_t[matrix_id_2[i, 1], matrix_id_2[i, 2], ],
                       c(K_c + 1), weights = rc_betas)
      )
    }
    # Pure RC model without breaks: the smoothed RC IS the long-run target.
    # 纯 RC 无突变:平滑后的 RC 即长期目标,无 tanh。
    if (J == 0L && R == 1L) return(RC_tt)
  }

  ## --- MIDAS filter of each exogenous covariate (pair-specific w2) ---
  ## 各外生协变量的 MIDAS 滤波(w2 逐资产对)
  mid_arrays <- vector("list", J)
  if (J > 0L) {
    for (j in seq_len(J)) {
      w2_pair <- pair_matrix(param[idx$w2[[j]]], N)
      bm <- array(0, dim = c(N, N, TT))
      for (i in 1:nrow(matrix_id_2)) {
        betas_ij <- c(rev(midas_weights(1:(K_c + 1), (K_c + 1), w1,
                                        w2_pair[matrix_id_2[i, 1], matrix_id_2[i, 2]],
                                        dc$lag_fun))[2:(K_c + 1)], 0)
        midas_c <- suppressWarnings(
          roll::roll_sum(dc$X_mats[[j]], c(K_c + 1), weights = betas_ij)
        )
        bm[matrix_id_2[i, 1], matrix_id_2[i, 2], ] <- midas_c[(K_c + 1), ]
      }
      mid_arrays[[j]] <- bm
    }
  }

  ## --- assemble z and apply the Fisher-type transform ---
  ## 组装 z 并做 tanh 变换
  if (R == 1L) {
    # Pair-specific level and slopes / 逐对截距与斜率
    z <- array(pair_matrix(param[idx$delta], N), dim = c(N, N, TT))
    for (j in seq_len(J)) {
      z <- z + array(pair_matrix(param[idx$theta[[j]]], N), dim = c(N, N, TT)) * mid_arrays[[j]]
    }
    if (dc$rc) {
      z <- z + array(pair_matrix(param[idx$theta_rc], N), dim = c(N, N, TT)) * RC_tt
    }
  } else {
    # Regime-varying, pair-common level/slopes via the break dummies;
    # only off-diagonal pair positions are filled (legacy behavior).
    # 区制分段、跨对共享;仅填充非对角资产对位置(与旧脚本一致)。
    m_t <- as.vector(param[idx$m] %*% dc$D_mat)
    z <- array(0, dim = c(N, N, TT))
    pairs_idx <- which(lower.tri(matrix(0, N, N)), arr.ind = TRUE)

    term_series <- list()
    if (dc$rc) term_series$rc <- as.vector(param[idx$theta_rc] %*% dc$D_mat)
    for (j in seq_len(J)) {
      term_series[[paste0("x", j)]] <- as.vector(param[idx$theta[[j]]] %*% dc$D_mat)
    }

    for (k in 1:nrow(pairs_idx)) {
      i_idx <- pairs_idx[k, 1]; j_idx <- pairs_idx[k, 2]
      z_k <- m_t
      if (dc$rc) z_k <- z_k + term_series$rc * RC_tt[i_idx, j_idx, ]
      for (j in seq_len(J)) {
        z_k <- z_k + term_series[[paste0("x", j)]] * mid_arrays[[j]][i_idx, j_idx, ]
      }
      z[i_idx, j_idx, ] <- z_k
      z[j_idx, i_idx, ] <- z_k
    }
  }

  tanh(z)
}
