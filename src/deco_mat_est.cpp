#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List new_deco_mat_est(NumericVector est_param, arma::mat res, arma::cube Dt, Nullable<int> K_c_in = R_NilValue, Nullable<NumericMatrix> S_in = R_NilValue) {
  double a = est_param[0];
  double b = est_param[1];
  
  int TT = res.n_rows;
  int Num_col = res.n_cols;
  int K_c = K_c_in.isNull() ? 2 : as<int>(K_c_in) + 1;
  
  // S_in: optional in-sample target for look-ahead-free full-sample
  // filtering (default keeps legacy behavior cov(res)).
  // S_in:可选样本内目标,供无前视的全样本滤波;缺省保持旧行为。
  mat S = S_in.isNull() ? cov(res) : as<mat>(S_in.get());
  mat I_n = eye(Num_col, Num_col);
  mat J_n = ones(Num_col, Num_col);
  
  // 预分配存储空间
  cube Q_seq(Num_col, Num_col, TT);
  cube R_t_deco(Num_col, Num_col, TT);
  cube H_t(Num_col, Num_col, TT);
  
  // 初始化
  for(int i = 0; i < TT; i++) {
    Q_seq.slice(i) = eye(Num_col, Num_col);
    R_t_deco.slice(i) = eye(Num_col, Num_col);
    H_t.slice(i) = S;
  }
  
  double first_elem = 1.0 / (double(Num_col) * (double(Num_col) - 1.0));
  
  for (int t = K_c - 1; t < TT; t++) {
    rowvec r_prev = res.row(t - 1);
    
    // 1. DCC 递归计算 Q_t
    mat Q_prev = Q_seq.slice(t - 1);
    mat Q_star_p = diagmat(sqrt(Q_prev.diag()));
    Q_seq.slice(t) = (1 - a - b) * S + a * (Q_star_p * r_prev.t() * r_prev * Q_star_p) + b * Q_prev;
    
    // 2. 计算标准 DCC 相关矩阵 R_t
    vec inv_sqrt_q = 1.0 / sqrt(Q_seq.slice(t).diag());
    mat inv_Q_star = diagmat(inv_sqrt_q);
    mat R_t_dcc = inv_Q_star * Q_seq.slice(t) * inv_Q_star;
    
    // 3. 计算 DECO rho (利用 accu 替代矩阵乘法 iota' * R * iota)
    double rho_t = first_elem * (accu(R_t_dcc) - double(Num_col));
    
    // 4. 构建 DECO 相关矩阵与协方差矩阵
    R_t_deco.slice(t) = (1.0 - rho_t) * I_n + rho_t * J_n;
    mat D_curr = Dt.slice(t);
    H_t.slice(t) = D_curr * R_t_deco.slice(t) * D_curr;
  }
  
  // 5. 计算 H_t[,,1] 的均值 (对应 R 的 apply(H_t[,,2:TT], c(1,2), mean))
  // 注意：R 的 2:TT 对应 C++ 的索引 1 到 TT-1
  mat H_sum = zeros(Num_col, Num_col);
  for(int i = 1; i < TT; i++) {
    H_sum += H_t.slice(i);
  }
  H_t.slice(0) = H_sum / double(TT - 1);
  
  return List::create(
    Named("H_t") = H_t,
    Named("R_t") = R_t_deco
  );
}