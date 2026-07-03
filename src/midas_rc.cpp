#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::cube midas_rc_rcpp(int N_c, int Num_assets, int TT, arma::mat res) {
  
  // 1. 初始化结果三维数组 (Num_assets x Num_assets x TT)
  // fill::zeros 会自动把初始值填为 0
  arma::cube C_t(Num_assets, Num_assets, TT, fill::zeros);
  
  // 2. 循环遍历时间 t
  // 注意：R 的循环是从 (N_c+1) 到 TT (1-based index)
  // 对应 C++ 的循环是从 N_c 到 TT-1 (0-based index)
  for (int t = N_c; t < TT; t++) {
    
    // 3. 提取窗口数据 window_res
    // R: res[(tt-N_c):tt, ] -> 长度是 N_c + 1
    // Armadillo: rows(start_row, end_row) -> 包含两端
    // 这里的 t 对应 R 中的 tt-1，所以索引逻辑是直接平移的
    arma::mat window = res.rows(t - N_c, t);
    
    // 4. 计算 V_t_0.5 (列平方和的开根号)
    // sum(X, 0) 对列求和，返回一个行向量
    arma::rowvec V = sqrt(sum(square(window), 0));
    
    // 5. 计算 Prod_eps_t (交叉乘积 X'X)
    // 这相当于 R 中的 crossprod(window)
    arma::mat P = window.t() * window;
    
    // 6. 计算分母矩阵 (V 的外积)
    // tcrossprod(V) 在 R 中如果是向量，等价于 V %*% t(V)
    // 也就是生成一个矩阵，元素 (i,j) = V[i] * V[j]
    arma::mat Denom = V.t() * V;
    
    // 7. 逐元素相除并存入结果 Cube
    // Armadillo 的 C_t.slice(t) 直接引用内存切片
    // P / Denom 在 C++ 中是矩阵求解(Solve)，所以这里我们手动做逐元素除法
    // 或者使用 P % (1/Denom) 但手动循环对于几百个资产来说更快且更安全
    for(int i = 0; i < Num_assets; i++) {
      for(int j = 0; j < Num_assets; j++) {
        // 如果分母为0，Armadillo默认会处理为Inf，这里直接除即可
        C_t(i, j, t) = P(i, j) / Denom(i, j);
      }
    }
  }
  
  return C_t;
}