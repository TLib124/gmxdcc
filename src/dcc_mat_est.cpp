#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List new_dcc_mat_est(NumericVector est_param, arma::mat res, arma::cube Dt, Nullable<int> K_c_in = R_NilValue, Nullable<NumericMatrix> S_in = R_NilValue) {
  double a = est_param[0];
  double b = est_param[1];
  
  int TT = res.n_rows;
  int Num_col = res.n_cols;
  
  // 严格对应 R: if(is.null(K_c)) K_c<-2 else K_c<-K_c+1
  int K_c;
  if (K_c_in.isNull()) {
    K_c = 2; 
  } else {
    K_c = as<int>(K_c_in) + 1;
  }
  
  // S_in: optional unconditional target. Pass the in-sample covariance when
  // filtering the full (in+out-of-sample) series with frozen parameters, so
  // the intercept has no look-ahead. Default keeps legacy behavior cov(res).
  // S_in:可选无条件目标;冻结参数做全样本滤波时传入样本内协方差,
  // 避免前视;缺省保持旧行为 cov(res)。
  mat S;
  if (S_in.isNull()) {
    S = cov(res);
  } else {
    S = as<mat>(S_in.get());
  }

  // 预分配三维数组 (Cube)
  // R 的 array(diag(1), dim=...) 在 C++ 中对应 slice 初始化
  cube R_t(Num_col, Num_col, TT);
  cube H_t(Num_col, Num_col, TT);
  cube Q_seq(Num_col, Num_col, TT);
  
  // 初始化
  for(int i = 0; i < TT; i++) {
    R_t.slice(i) = eye(Num_col, Num_col);
    H_t.slice(i) = S; // 对应 R: H_t<-array(S, ...)
    Q_seq.slice(i) = eye(Num_col, Num_col);
  }
  
  // 循环：从 K_c 开始 (C++ 索引为 K_c - 1)
  for (int t = K_c - 1; t < TT; t++) {
    rowvec res_prev = res.row(t - 1);
    
    // Q_t_star 对应上一期 Q 的对角线根号
    mat Q_prev = Q_seq.slice(t - 1);
    mat Q_t_star_prev = diagmat(sqrt(Q_prev.diag()));
    
    // 更新 Q_t
    Q_seq.slice(t) = (1 - a - b) * S + a * (Q_t_star_prev * res_prev.t() * res_prev * Q_t_star_prev) + b * Q_prev;
    
    // 计算 R_t
    vec inv_sqrt_diag = 1.0 / sqrt(Q_seq.slice(t).diag());
    mat inv_Q_star_curr = diagmat(inv_sqrt_diag);
    R_t.slice(t) = inv_Q_star_curr * Q_seq.slice(t) * inv_Q_star_curr;
    
    // 计算 H_t = Dt %*% R_t %*% Dt
    mat D_curr = Dt.slice(t);
    H_t.slice(t) = D_curr * R_t.slice(t) * D_curr;
  }
  
  return List::create(
    Named("H_t") = H_t,
    Named("R_t") = R_t
  );
}