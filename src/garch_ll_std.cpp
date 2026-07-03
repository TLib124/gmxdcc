#include <RcppArmadillo.h>
#include <cmath> // Ensures std::lgamma and std::log are available

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericVector fast_t_garch_ll(double alpha, double beta, double v, 
                              NumericVector epsilon, NumericVector tau_d, 
                              NumericVector daily_ret) {
  int TT = daily_ret.size();
  NumericVector x_step_1(TT);
  NumericVector g_it(TT);
  NumericVector ll(TT);
  
  // 1. Construct x_step_1 (equivalent to R's c(0, step_1[-TT]))
  for(int t = 1; t < TT; t++) {
    double val = (1.0 - alpha - beta) + (alpha * std::pow(epsilon[t-1], 2.0) / tau_d[t-1]);
    if (NumericVector::is_na(val)) val = 0.0;
    x_step_1[t] = val;
  }
  x_step_1[0] = 0.0;
  
  // Correct the initial offset in R logic
  if (TT > 1) {
    x_step_1[1] += beta * 1.0;
  }
  
  // 2. Simulate stats::filter(method = "recursive")
  g_it[0] = x_step_1[0]; 
  for(int t = 1; t < TT; t++) {
    g_it[t] = x_step_1[t] + beta * g_it[t-1];
  }
  
  // 3. Override the initial value: g_it[1] <- 1
  g_it[0] = 1.0;
  
  // 4. Calculate Student-t log-likelihood
  // Use std::lgamma from <cmath> instead of R::lgamma to avoid namespace conflicts
  double term1 = std::lgamma((v + 1.0) / 2.0);
  double term2 = std::lgamma(v / 2.0);
  double term3 = 0.5 * std::log(v - 2.0);
  double const_part = term1 - term2 - term3;
  
  for(int t = 0; t < TT; t++) {
    double h_it_2 = g_it[t] * tau_d[t];
    double log_c = const_part - 0.5 * std::log(h_it_2);
    double log_d = -((v + 1.0) / 2.0) * std::log(1.0 + std::pow(epsilon[t], 2.0) / (h_it_2 * (v - 2.0)));
    
    ll[t] = log_c + log_d;
  }
  
  return ll;
}