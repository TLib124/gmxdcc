#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericVector new_dcc_loglik(NumericVector param, mat res, Nullable<int> K_c_in = R_NilValue) {
  double a = param[0];
  double b = param[1];
  
  int TT = res.n_rows;
  int Num_col = res.n_cols;
  
  // 严格对应 R: if(is.null(K_c)) K_c<-2 else K_c<-K_c+1
  int K_c;
  if (K_c_in.isNull()) {
    K_c = 2; 
  } else {
    K_c = as<int>(K_c_in) + 1;
  }
  
  mat S = cov(res);
  // 预分配 Q 序列，初始化为单位阵（对应 R 的 array(diag(1), ...)）
  field<mat> Q_seq(TT);
  for(int i = 0; i < TT; i++) {
    Q_seq(i) = eye(Num_col, Num_col);
  }
  
  NumericVector ll(TT); // 初始为 0
  
  // R 循环：for(tt in K_c:TT)
  // 对应 C++ 索引：for(int t = K_c - 1; t < TT; t++)
  for (int t = K_c - 1; t < TT; t++) {
    
    // 对应 R: res_row1 <- matrix(res[tt-1,], nrow = 1)
    // C++ 索引比 R 小 1，所以 res[tt-1,] 对应 res.row(t-1)
    rowvec res_prev = res.row(t - 1);
    rowvec res_curr = res.row(t);
    
    // 对应 R: Q_t_star <- sqrt(diag(diag(Q_t[,,tt-1])))
    mat Q_prev = Q_seq(t - 1);
    mat Q_t_star_prev = diagmat(sqrt(Q_prev.diag()));
    
    // 更新当前的 Q_t: Q_seq(t)
    Q_seq(t) = (1 - a - b) * S + a * (Q_t_star_prev * res_prev.t() * res_prev * Q_t_star_prev) + b * Q_prev;
    
    // 计算 R_t
    vec inv_sqrt_diag = 1.0 / sqrt(Q_seq(t).diag());
    mat inv_Q_star_curr = diagmat(inv_sqrt_diag);
    mat R_t = inv_Q_star_curr * Q_seq(t) * inv_Q_star_curr;
    
    // 计算似然
    double val, sign;
    log_det(val, sign, R_t);
    mat eps_R_inv_eps = res_curr * solve(R_t, res_curr.t(), solve_opts::fast);
    
    ll[t] = -(val + eps_R_inv_eps(0,0));
  }
  
  return ll;
}