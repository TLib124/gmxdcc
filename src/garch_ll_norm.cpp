#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fast_garch_ll_cpp(double alpha, double beta, double mu, 
                                NumericVector epsilon, NumericVector tau_d, 
                                NumericVector daily_ret) {
  int TT = daily_ret.size();
  NumericVector x_step_1(TT);
  NumericVector g_it(TT);
  NumericVector ll(TT);
  
  // 1. 预计算 step_1 并构造 x_step_1
  // R: x_step_1 <- c(0, step_1[-TT])
  for(int t = 1; t < TT; t++) {
    double val = (1.0 - alpha - beta) + (alpha * std::pow(epsilon[t-1], 2.0) / tau_d[t-1]);
    if (NumericVector::is_na(val)) val = 0.0;
    x_step_1[t] = val;
  }
  x_step_1[0] = 0.0; // c(0, ...)
  
  // R: x_step_1[2] <- x_step_1[2] + beta * 1
  if (TT > 1) {
    x_step_1[1] += beta * 1.0;
  }
  
  // 2. 模拟 stats::filter(..., method = "recursive")
  // y[t] = x[t] + beta * y[t-1]
  g_it[0] = x_step_1[0]; 
  for(int t = 1; t < TT; t++) {
    g_it[t] = x_step_1[t] + beta * g_it[t-1];
  }
  
  // 3. 覆盖初始值 R: g_it[1] <- 1
  g_it[0] = 1.0;
  
  // 4. 计算对数似然
  double log2pi = std::log(2.0 * M_PI);
  for(int t = 0; t < TT; t++) {
    double var_t = g_it[t] * tau_d[t];
    ll[t] = -0.5 * (log2pi + std::log(var_t) + std::pow(daily_ret[t] - mu, 2.0) / var_t);
  }
  
  return ll;
}