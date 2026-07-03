#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericVector new_a_dcc_loglik(NumericVector param, mat res, Nullable<int> K_c_in = R_NilValue) {
  double a = param[0];
  double b = param[1];
  double g = param[2];
  
  int TT = res.n_rows;
  int Num_col = res.n_cols;
  
  // 严格对应 R 的 K_c 逻辑
  int K_c = K_c_in.isNull() ? 2 : as<int>(K_c_in) + 1;
  
  // 计算 eta (仅保留负残差)
  mat eta = res;
  eta.transform( [](double x) { return (x < 0.0) ? x : 0.0; } );
  
  mat S = cov(res);
  // 注意：你原 R 代码中 N <- stats::cov(res)，但我推测 A-DCC 中 N 应该是 cov(eta)
  // 这里先严格遵循你代码写的 S 和 N (即 N = S)
  mat N_mat = S; 
  
  // 初始化 Q 序列为单位阵 (对应 R 的 array(diag(1), ...))
  field<mat> Q_seq(TT);
  for(int i = 0; i < TT; i++) Q_seq(i) = eye(Num_col, Num_col);
  
  NumericVector ll(TT);
  
  // 循环从 K_c 开始 (C++ 索引比 R 小 1)
  for (int t = K_c - 1; t < TT; t++) {
    rowvec r_prev = res.row(t - 1);
    rowvec e_prev = eta.row(t - 1);
    rowvec r_curr = res.row(t);
    
    // 更新 Q_t: (1-a-b)*S - g*N + a*(r_p' * r_p) + b*Q_p + g*(e_p' * e_p)
    // 在 Armadillo 中，rowvec.t() * rowvec 得到的是 N x N 矩阵
    Q_seq(t) = (1 - a - b) * S - g * N_mat + 
      a * (r_prev.t() * r_prev) + 
      b * Q_seq(t - 1) + 
      g * (e_prev.t() * e_prev);
    
    // 计算相关矩阵 R_t
    vec inv_sqrt_q = 1.0 / sqrt(Q_seq(t).diag());
    mat inv_Q_star = diagmat(inv_sqrt_q);
    mat R_t = inv_Q_star * Q_seq(t) * inv_Q_star;
    
    // 计算似然项
    double val, sign;
    log_det(val, sign, R_t); // 对数行列式
    mat quad_form = r_curr * solve(R_t, r_curr.t(), solve_opts::fast);
    
    ll[t] = -(val + quad_form(0,0));
  }
  
  return ll;
}