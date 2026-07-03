#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec midas_ll_rcpp(int TT, int Num_assets, int K_c, double a, double b, arma::cube R_t_bar, arma::mat res) {
  
  // --- 1. 初始化变量 ---
  // Q_t: 初始化为单位矩阵 (Identity Matrix)
  arma::cube Q(Num_assets, Num_assets, TT, fill::zeros);
  for(int i = 0; i < TT; i++) {
    Q.slice(i).eye();
  }
  
  // 结果向量
  arma::vec log_det_R_t(TT, fill::zeros);
  arma::vec Eps_R_Eps(TT, fill::zeros);
  arma::vec ll(TT, fill::zeros);
  
  // --- 2. 预处理 R_t_bar (对应 R 中的 diag(R_t_bar[,,i]) <- 1) ---
  // 注意：这将直接修改传入的 R_t_bar 副本，不影响 R 中的原对象(除非使用 clone=false)
  for(int i = 0; i < TT; i++) {
    R_t_bar.slice(i).diag().ones();
  }
  
  // --- 3. 主循环 ---
  // R: for(tt in (K_c+1):TT) 
  // R 的索引是从 1 开始的，C++ 是从 0 开始的。
  // R 的 index K_c+1 对应 C++ 的 index K_c。
  for (int t = K_c; t < TT; t++) {
    
    // A. 准备数据
    arma::vec eps_prev = res.row(t - 1).t(); // res[tt-1, ]
    arma::mat Q_prev = Q.slice(t - 1);       // Q[,,tt-1]
    
    // B. 更新 Q_t
    // Q_t[,,tt] <- (1-a-b)*R_t_bar[,,tt] + a*tcrossprod(res[tt-1, ]) + b*Q_t[,,tt-1]
    Q.slice(t) = (1.0 - a - b) * R_t_bar.slice(t) + a * (eps_prev * eps_prev.t()) + b * Q_prev;
    
    // C. 计算 R_t (标准化)
    // q_sd <- sqrt(diag(Q_t[,,tt]))
    arma::vec q_sd = sqrt(Q.slice(t).diag());
    
    // =========================================================
    // 修改 2: 防止除以 0 产生 NaN
    // 如果 q_sd 中有极小值，替换为一个极小的正数，防止 NaN 污染
    // =========================================================
    for(int k=0; k<Num_assets; ++k){
      if(q_sd(k) < 1e-12) q_sd(k) = 1e-12; 
    }
    
    
    // R_t[,,tt] <- Q_t[,,tt] / tcrossprod(q_sd)
    // 优化：不生成 R_t 的 Cube，直接计算当前 R_t 矩阵以节省内存
    arma::mat R_curr = Q.slice(t);
    // 利用 Armadillo 的 broadcasting 机制进行行/列除法 (极快)
    R_curr.each_col() /= q_sd;
    R_curr.each_row() /= q_sd.t();
    
    // === FIX: 强制对称化 (解决浮点误差导致的 chol 失败) ===
    // 这一步至关重要，防止 Armadillo 抛出 "matrix not symmetric"
    arma::mat R_sym = 0.5 * (R_curr + R_curr.t());
    R_sym.diag().ones(); // 强制对角线为 1
    // ====================================================
    
    // =========================================================
    // 修改 3: 检查 NaN (如果包含 NaN，chol 必报不对称)
    // =========================================================
    if (!R_sym.is_finite()) {
      ll.fill(-1e10); 
      return ll;
    }
    
    // D. Cholesky 分解
    arma::mat U;
    // chol 返回 bool，对应 R 的 tryCatch
    bool success = chol(U, R_sym); // U 是上三角矩阵
    
    if (!success) {
      // 对应 R 的 if (is.null(chol_R)) return rep(-1e10...)
      // 这里为了不中断整个向量返回，我们填入一个极小的数
      ll.fill(-1e10); 
      return ll; 
    }
    
    // E. 计算 Log Det
    // log_det_R_t[tt] <- 2 * sum(log(diag(chol_R)))
    double log_det = 2.0 * sum(log(U.diag()));
    log_det_R_t(t) = log_det;
    
    // F. 计算 quadratic form (backsolve)
    // tmp <- backsolve(chol_R, res[tt, ], transpose = TRUE)
    // 也就是解方程 U'x = res[tt]
    arma::vec eps_curr = res.row(t).t();
    
    // trimatu(U).t() 明确告诉编译器 U 是上三角，转置后是下三角
    arma::vec tmp = solve(trimatu(U).t(), eps_curr);
    
    double quad_form = dot(tmp, tmp); // sum(tmp^2)
    Eps_R_Eps(t) = quad_form;
    
    // G. 最终似然
    ll(t) = -0.5 * (log_det + quad_form);
  }
  
  return ll;
}