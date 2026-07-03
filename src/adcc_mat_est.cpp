#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List new_a_dcc_mat_est(NumericVector est_param, arma::mat res, arma::cube Dt, Nullable<int> K_c_in = R_NilValue, Nullable<NumericMatrix> S_in = R_NilValue, Nullable<NumericMatrix> N_in = R_NilValue) {
  double a = est_param[0];
  double b = est_param[1];
  double g = est_param[2];
  
  int TT = res.n_rows;
  int Num_col = res.n_cols;
  
  // 严格对应 R 逻辑: K_c+1
  int K_c = K_c_in.isNull() ? 2 : as<int>(K_c_in) + 1;
  
  // 1. 在 C++ 中计算 eta (负残差矩阵)
  mat eta = res;
  eta.transform( [](double x) { return (x < 0.0) ? x : 0.0; } );
  
  // 2. 计算无条件协方差 S 和 N
  // S_in/N_in: optional in-sample targets for look-ahead-free full-sample
  // filtering (default keeps legacy behavior).
  // S_in/N_in:可选样本内目标,供无前视的全样本滤波;缺省保持旧行为。
  mat S     = S_in.isNull() ? cov(res) : as<mat>(S_in.get());
  mat N_mat = N_in.isNull() ? cov(eta) : as<mat>(N_in.get());
  
  // 3. 预分配三维 Cube (对应 R 的 array)
  cube R_t(Num_col, Num_col, TT);
  cube H_t(Num_col, Num_col, TT);
  cube Q_seq(Num_col, Num_col, TT);
  
  // 初始化首层为单位阵和 S
  for(int i = 0; i < TT; i++) {
    R_t.slice(i) = eye(Num_col, Num_col);
    H_t.slice(i) = S;
    Q_seq.slice(i) = eye(Num_col, Num_col);
  }
  
  // 4. 递归计算
  for (int t = K_c - 1; t < TT; t++) {
    rowvec r_p = res.row(t - 1);
    rowvec e_p = eta.row(t - 1);
    
    // 更新 Q_t: 注意转置顺序 r_p.t() * r_p 得到 N x N
    Q_seq.slice(t) = (1 - a - b) * S - g * N_mat + 
      a * (r_p.t() * r_p) + 
      b * Q_seq.slice(t - 1) + 
      g * (e_p.t() * e_p);
    
    // 计算相关矩阵 R_t
    vec inv_q_diag = 1.0 / sqrt(Q_seq.slice(t).diag());
    mat inv_Q_star = diagmat(inv_q_diag);
    R_t.slice(t) = inv_Q_star * Q_seq.slice(t) * inv_Q_star;
    
    // 计算协方差矩阵 H_t = Dt * R_t * Dt
    mat D_curr = Dt.slice(t);
    H_t.slice(t) = D_curr * R_t.slice(t) * D_curr;
  }
  
  return List::create(
    Named("H_t") = H_t,
    Named("R_t") = R_t
  );
}