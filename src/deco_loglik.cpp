#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericVector new_deco_loglik(NumericVector param, arma::mat res, Nullable<int> K_c_in = R_NilValue) {
  double a = param[0];
  double b = param[1];
  
  int TT = res.n_rows;
  int Num_col = res.n_cols;
  int K_c = K_c_in.isNull() ? 2 : as<int>(K_c_in) + 1;
  
  mat S = cov(res);
  mat I_n = eye(Num_col, Num_col);
  mat J_n = ones(Num_col, Num_col);
  
  // 初始化 Q 序列为单位阵
  field<mat> Q_seq(TT);
  for(int i = 0; i < TT; i++) Q_seq(i) = eye(Num_col, Num_col);
  
  NumericVector ll(TT);
  double first_elem = 1.0 / (double(Num_col) * (double(Num_col) - 1.0));
  
  for (int t = K_c - 1; t < TT; t++) {
    rowvec res_prev = res.row(t - 1);
    rowvec res_curr = res.row(t);
    
    // 1. DCC 递归部分 (Q_t)
    mat Q_prev = Q_seq(t - 1);
    mat Q_t_star_prev = diagmat(sqrt(Q_prev.diag()));
    
    Q_seq(t) = (1 - a - b) * S + a * (Q_t_star_prev * res_prev.t() * res_prev * Q_t_star_prev) + b * Q_prev;
    
    // 2. 计算标准 DCC 相关矩阵 R_t
    vec inv_sqrt_q = 1.0 / sqrt(Q_seq(t).diag());
    mat inv_Q_star = diagmat(inv_sqrt_q);
    mat R_t_dcc = inv_Q_star * Q_seq(t) * inv_Q_star;
    
    // 3. 计算 DECO 均值相关系数 rho_t
    // 对应 R: first_elem * (sum(R_t) - Num_col)
    double rho_t = first_elem * (accu(R_t_dcc) - double(Num_col));
    
    // 4. 构建 DECO 等相关矩阵
    mat R_t_deco = (1.0 - rho_t) * I_n + rho_t * J_n;
    
    // 5. 计算似然项
    double val, sign;
    log_det(val, sign, R_t_deco);
    mat quad_form = res_curr * solve(R_t_deco, res_curr.t(), solve_opts::fast);
    
    ll[t] = -(val + quad_form(0,0));
  }
  
  return ll;
}