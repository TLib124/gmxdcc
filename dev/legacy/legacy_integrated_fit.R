library(dccmidas)
library(quantmod)
library(zoo)
library(rumidas)
library(xts)
library(lubridate)
library(expm)
library(Rcpp)
library(RcppArmadillo)

# 1. 编译 C++ 函数
# 确保你的工作目录在 .cpp 文件所在的文件夹
sourceCpp("midas_rc_rcpp_turbo.cpp")
sourceCpp("midas_ll_rcpp_turbo.cpp")
sourceCpp("garch_midas_ll_norm_turbo.cpp")
sourceCpp("garch_midas_ll_std_turbo.cpp")
sourceCpp("new_dcc_loglik_turbo.cpp") 
sourceCpp("new_dcc_mat_est_turbo.cpp")
sourceCpp("new_a_dcc_loglik_turbo.cpp")
sourceCpp("new_a_dcc_mat_est_turbo.cpp")
sourceCpp("new_deco_loglik_turbo.cpp")
sourceCpp("new_deco_mat_est_turbo.cpp")


##univ_model:sGARCH/gjrGARCH/eGARCH/iGARCH/csGARCH/
##GM_noskew_RV/GM_noskew_X/GM_noskew_RV_X/GM_noskew_RV_sb/GM_noskew_X_sb/GM_noskew_RV_X_sb

##corr_model:cDCC/DECO/aDCC/
##DCCMIDAS_RC/DCCMIDAS_RC_sb/DCCMIDAS_X/DCCMIDAS_X_sb/DCCMIDAS_RC_X/DCCMIDAS_RC_X_sb


##GM_noskew_RV/GM_noskew_X/GM_noskew_RV_X/GM_noskew_RV_sb/GM_noskew_X_sb/GM_noskew_RV_X_sb

roll_sum_func <- roll::roll_sum
weight_func_beta <- rumidas::beta_function
weight_func_exp <- rumidas::exp_almon

Inf_criteria <- function(est) {
  
  # --- 情况 1: 输入是 rugarch 对象 (比如 sGARCH, gjrGARCH) ---
  if (inherits(est, "uGARCHfit")) {
    # rugarch 包有专门的方法提取似然值和数据
    ll <- rugarch::likelihood(est)
    k <- length(rugarch::coef(est))
    n <- length(est@model$modeldata$data) # 获取样本量
    
    # --- 情况 2: 输入是 maxLik 对象 (GARCH-MIDAS 内部优化结果) ---
  } else if (inherits(est, "maxLik")) {
    ll <- as.numeric(stats::logLik(est))
    k <- length(stats::coef(est))
    # maxLik 对象包含梯度矩阵，行数即为样本量
    n <- nrow(est$gradientObs)
    
    # --- 情况 3: 输入是 gmnoskew 对象 (GARCH-MIDAS 最终结果列表) ---
  } else if (inherits(est, "gmnoskew")) {
    # 如果对象里已经存了结果，直接返回 (避免重复计算)
    if (!is.null(est$inf_criteria)) {
      return(est$inf_criteria)
    }
    # 否则手动提取
    ll <- est$loglik
    k <- nrow(est$rob_coef_mat)
    n <- est$obs
    
  } else {
    # 为了防止再次报错，给一个未知的默认处理或报错
    stop("Inf_criteria: Unknown object class provided.")
  }
  
  # --- 统一计算 AIC 和 BIC ---
  aic <- 2 * k - 2 * ll
  bic <- k * log(n) - 2 * ll
  
  inf <- round(c(AIC = aic, BIC = bic), 6)
  return(inf)
}

LF_f<-function(vol_est,vol_proxy){
  
  LF_avg<-c(
    "MSE(%)"=100*mean( (vol_est-vol_proxy)^2),
    "QLIKE"=mean( log(vol_est) + (vol_proxy/vol_est) )
  )
  
  return(LF_avg)
  
}

QMLE_sd <- function(est){
  
  # 获取海森矩阵
  H <- -hessian(est)
  
  # 尝试标准求逆
  H_hat_inv <- tryCatch({
    solve(H)
  }, error = function(e) {
    # 如果报错（奇异矩阵），使用广义逆 (Generalized Inverse)
    message("Warning: Hessian is singular. Switching to generalized inverse (MASS::ginv).")
    # 确保加载 MASS 包
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Please install 'MASS' package.")
    return(MASS::ginv(H))
  })
  
  # 计算 OPG 矩阵
  OPG <- t(est$gradientObs) %*% est$gradientObs
  
  # 三明治方差估计量
  var_cov <- H_hat_inv %*% OPG %*% H_hat_inv
  
  # 提取对角线（方差）
  diag_vals <- diag(var_cov)
  
  # 处理可能出现的负值（数值误差导致），将其设为 NA
  diag_vals[diag_vals < 0] <- NA
  
  return(sqrt(diag_vals))
}

#cut function
cut_func <- function(r_t,K=NULL,K_c=NULL,macro=NULL,type=NULL){
  
  if (is.null(K)){
    K <- 0
  }
  
  if (is.null(K_c)){
    K_c <- 0
  }
  
  ####available for monthly/quarterly
  if(type=="monthly"){
    
    cut_size <- max(K,K_c)+1
    year_month <- unique(format(index(macro), "%Y-%m"))
    dates <- as.Date(paste0(year_month, "-01"))
    shifted_dates <- seq(dates[1], by = "month", length.out = length(dates) + cut_size)[-c(1:cut_size)]
    shifted_dates <- shifted_dates[-((length(shifted_dates) - cut_size+1):length(shifted_dates))]
    year_month <- format(shifted_dates, "%Y-%m")
    required_r_t <- r_t[format(index(r_t), "%Y-%m") %in% year_month]
    
  }else if(type=="quarterly"){
    
    cut_size <- max(K,K_c)+1
    year_quarter <- unique(format(index(macro), "%Y-%m"))
    dates <- as.Date(paste0(year_quarter, "-01"))
    shifted_dates <- seq(dates[1], by = "quarter", length.out = length(dates) + cut_size)[-c(1:(cut_size))]
    shifted_dates <- shifted_dates[-((length(shifted_dates) - cut_size):length(shifted_dates))]
    year_quarter_month_series <- seq(from = shifted_dates[1], to = shifted_dates[length(shifted_dates)], by = "month")
    year_quarter_month <- format(year_quarter_month_series, "%Y-%m")
    
    required_r_t <- r_t[format(index(r_t), "%Y-%m") %in% year_quarter_month]
    
  } else if(is.null(macro)|is.null(type)){
    
    required_r_t <- r_t
    
  }
  
  return(required_r_t)
}

#MIDAS variable matrix transformation
mv_into_mat2<-function(x,mv,K,type){
  
  N<-length(x)
  N_mv<-length(mv)
  
  cond_r_t<- class(x)[1]
  cond_mv<- class(mv)[1]
  
  first_obs_r_t<- (first(index(x)))
  first_obs_mv<- (first(index(mv)))
  diff_time<-as.numeric(difftime(first_obs_r_t,first_obs_mv,units="day"))
  num_obs_k<-ifelse(type=="weekly",5,
                    ifelse(type=="monthly",30,
                           ifelse(type=="quarterly",120,365)))
  
  
  w_days_r_t<- lubridate::wday(index(xts::last(x)))
  
  
  last_mv<-paste(xts::last(lubridate::year(mv)),xts::last(lubridate::month(mv)),"20",sep="-")
  last_r_t<-paste(xts::last(lubridate::year(x)),xts::last(lubridate::month(x)),"20",sep="-")
  
  
  
  ############# checks 
  
  if(cond_r_t != "xts") { stop("x must be an xts object. Please provide it in the correct form")}
  if(cond_mv != "xts") { stop("mv must be an xts object. Please provide it in the correct form")}
  
  
  if(diff_time<0|((diff_time - K*num_obs_k)<0)) { stop("mv must start at least
'K+1' periods before x. Please decrease K or provide more observations (in the past) for mv")}
  
  
  if(last_mv != last_r_t & type=="weekly") { stop("x and mv must end in the same period")}
  if(w_days_r_t != 6 & type=="weekly") { stop("the last day of x must be friday")}
  
  
  mv_m<-matrix(c(rep(NA,(K+1)*N)),ncol=N)
  
  
  if(type=="weekly"){
    for(i in 1:N){
      start_p_mv1     <- which(format(index(mv),"%Y-%m-%d")==
                                 format(
                                   lubridate::floor_date(index(x), "week",week_start = 1)[i],"%Y-%m-%d"))
      mv_m[,i]         <- mv[(start_p_mv1-K):(start_p_mv1+1)][2:(K+2)]
    }
    
  } else if (type=="monthly") {
    for(i in 1:N){
      start_p_mv1     <- which(format(index(mv),"%Y-%m")==format(index(x)[i],"%Y-%m"))
      mv_m[,i]         <- mv[(start_p_mv1-K):(start_p_mv1)]
    }
  } else if (type=="quarterly"){
    for(i in 1:N){
      start_p_mv1     <- which(format(index(mv),"%Y-%m-%d")==
                                 format(
                                   lubridate::floor_date(index(x), "quarter",week_start = 1)[i],"%Y-%m-%d"))
      mv_m[,i]         <- mv[(start_p_mv1-K):(start_p_mv1+1)][2:(K+2)]
      
    } 
    
  } else {
    for(i in 1:N){
      start_p_mv1     <- which(format(index(mv),"%Y-%m-%d")==
                                 format(
                                   lubridate::floor_date(index(x), "year",week_start = 1)[i],"%Y-%m-%d"))
      mv_m[,i]         <- c(zoo::coredata(mv[(start_p_mv1-K):(start_p_mv1)][2:(K+1)]),0)
      
    } 
    
  }
  
  return(mv_m)
}

#RV matrix transformation
RV_gen <- function(x,type){
  
  ############# checks 
  cond_r_t<- class(x)[1]
  if(cond_r_t != "xts") { stop("x must be an xts object. Please provide it in the correct form")}
  
  calc_vol <- function(x) sum(x^2)
  
  
  if(type=="weekly"){
    
    weekly_vol <- apply.weekly(x, FUN = calc_vol)
    # 使用 cut 函数将日期切割并对齐到每周的周一
    # cut(..., "week") 默认以周一为起始
    index(weekly_vol) <- as.Date(cut(index(weekly_vol), "week"))
    colnames(weekly_vol) <- "Weekly_RV"
    
    vol <- weekly_vol
    
  } else if (type=="monthly") {
    
    monthly_vol <- apply.monthly(x, FUN = calc_vol)
    index(monthly_vol) <- as.Date(format(index(monthly_vol), "%Y-%m-01"))
    colnames(monthly_vol) <- "Monthly_RV"
    
    vol <- monthly_vol
    
  } else if (type=="quarterly"){
    
    quarterly_vol <- apply.quarterly(x, FUN = calc_vol)
    # 利用 zoo 包的 as.yearqtr 转换为季度对象，再转回 Date
    # as.Date(as.yearqtr(...)) 默认返回该季度的第一天
    index(quarterly_vol) <- as.Date(as.yearqtr(index(quarterly_vol)))
    colnames(quarterly_vol) <- "Quarterly_RV"
    
    vol <- quarterly_vol
    
  } else {
    
    yearly_vol <- apply.yearly(x, FUN = calc_vol)
    # 强制格式化为 "年-01-01"
    index(yearly_vol) <- as.Date(format(index(yearly_vol), "%Y-01-01"))
    colnames(yearly_vol) <- "Yearly_RV"
    
    vol <- yearly_vol
  }
  
  
  return(vol)
}

#loglik
garchmidas_noskew_loglik <- function(param,model,daily_ret,RV_m=NULL,mv_m=NULL,K,distribution,lag_fun="Beta",X_stru_break_date_vec=NULL){
  
  
  
  
  if (model=="GM_noskew_RV"){
    
    if(distribution=="norm"){
      
      alpha  <- param[1]
      beta  <- param[2]
      m   <- param[3]
      theta		<- param[4]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2			<- param[5]
      
      mu			<- param[6]
      
      TT    <- length(daily_ret)
      tau_d<-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      epsilon <- daily_ret - mu
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
      tau_d <-m+theta*suppressWarnings(roll_sum_func(RV_m, c(K+1),weights = betas)) 
      tau_d<-exp(tau_d[(K+1),])	 
      
      ll <- fast_garch_ll_cpp(alpha=alpha, beta=beta, mu=mu, 
                              epsilon=as.numeric(epsilon), tau_d=as.numeric(tau_d), daily_ret=as.numeric(daily_ret))
      
    } else if(distribution=="std"){
      alpha <- param[1]
      beta  <- param[2]
      m     <- param[3]
      theta		<- param[4]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2			<- param[5]
      v			<- param[6]
      
      mu			<- param[7]
      
      epsilon <- daily_ret - mu
      
      TT     <- length(daily_ret)
      tau_d<-rep(NA,TT)
      g_it  <- rep(1,TT)            # daily conditional variance 
      ll       <- 0 
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
      
      tau_d   <-  m+theta*suppressWarnings(roll_sum_func(RV_m, c(K+1),weights = betas)) 
      tau_d   <-  exp(tau_d[(K+1),])	 
      
      ll <- fast_t_garch_ll(alpha=alpha, beta=beta, v=v, 
                            epsilon=as.numeric(epsilon), tau_d=as.numeric(tau_d), daily_ret=as.numeric(daily_ret))
      
      
    }
    
  } else if (model=="GM_noskew_X"){
    
    if(distribution=="norm"){
      
      alpha  <- param[1]
      beta  <- param[2]
      m   <- param[3]
      theta		<- param[4]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2			<- param[5]
      
      mu			<- param[6]
      
      
      epsilon <- daily_ret - mu
      
      TT    <- length(daily_ret)
      tau_d<-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
      tau_d <-m+theta*suppressWarnings(roll_sum_func(mv_m, c(K+1),weights = betas)) 
      tau_d<-exp(tau_d[(K+1),])	 
      
      ll <- fast_garch_ll_cpp(alpha=alpha, beta=beta, mu=mu, 
                              epsilon=as.numeric(epsilon), tau_d=as.numeric(tau_d), daily_ret=as.numeric(daily_ret))
      
    } else if(distribution=="std"){
      alpha <- param[1]
      beta  <- param[2]
      m     <- param[3]
      theta		<- param[4]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2			<- param[5]
      v			<- param[6]
      mu			<- param[7]
      
      
      epsilon <- daily_ret - mu
      
      TT     <- length(daily_ret)
      tau_d<-rep(NA,TT)
      g_it  <- rep(1,TT)            # daily conditional variance 
      ll       <- 0 
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
      
      tau_d   <-  m+theta*suppressWarnings(roll_sum_func(mv_m, c(K+1),weights = betas)) 
      tau_d   <-  exp(tau_d[(K+1),])	 
      
      ll <- fast_t_garch_ll(alpha=alpha, beta=beta, v=v, 
                            epsilon=as.numeric(epsilon), tau_d=as.numeric(tau_d), daily_ret=as.numeric(daily_ret))
      
      
    }
    
  } else if (model=="GM_noskew_RV_X"){
    
    if(distribution=="norm"){
      
      alpha  <- param[1]
      beta  <- param[2]
      m   <- param[3]
      theta_RV		<- param[4]
      theta_X		<- param[5]
      
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2_RV			<- param[6]
      w2_X			<- param[7]
      
      mu			<- param[8]
      
      
      epsilon <- daily_ret - mu
      
      TT    <- length(daily_ret)
      
      tau_d_RV <-rep(NA,TT)
      tau_d_X <-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      
      betas_RV<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2_RV))[2:(K+1)],0)
      tau_d_RV <-suppressWarnings(roll_sum_func(RV_m, c(K+1),weights = betas_RV)) 
      tau_d_RV<-tau_d_RV[(K+1),]	 
      
      betas_X<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2_X))[2:(K+1)],0)
      tau_d_X <-suppressWarnings(roll_sum_func(mv_m, c(K+1),weights = betas_X)) 
      tau_d_X<-tau_d_X[(K+1),]
      
      tau_d <- exp(m + theta_RV*tau_d_RV + theta_X*tau_d_X)
      
      ll <- fast_garch_ll_cpp(alpha=alpha, beta=beta, mu=mu, 
                              epsilon=as.numeric(epsilon), tau_d=as.numeric(tau_d), daily_ret=as.numeric(daily_ret))
      
    } else if(distribution=="std"){
      alpha  <- param[1]
      beta  <- param[2]
      m   <- param[3]
      theta_RV		<- param[4]
      theta_X		<- param[5]
      
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2_RV			<- param[6]
      w2_X			<- param[7]
      
      v			<-   param[8]
      
      mu			<- param[9]
      
      
      epsilon <- daily_ret - mu
      
      TT    <- length(daily_ret)
      tau_d<-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      
      betas_RV<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2_RV))[2:(K+1)],0)
      tau_d_RV <-suppressWarnings(roll_sum_func(RV_m, c(K+1),weights = betas_RV)) 
      tau_d_RV<-tau_d_RV[(K+1),]	 
      
      betas_X<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2_X))[2:(K+1)],0)
      tau_d_X <-suppressWarnings(roll_sum_func(mv_m, c(K+1),weights = betas_X)) 
      tau_d_X<-tau_d_X[(K+1),]
      
      tau_d <- exp(m + theta_RV*tau_d_RV + theta_X*tau_d_X)
      
      ll <- fast_t_garch_ll(alpha=alpha, beta=beta, v=v, 
                            epsilon=as.numeric(epsilon), tau_d=as.numeric(tau_d), daily_ret=as.numeric(daily_ret))
      
      
    }
    
  } else if (model=="GM_noskew_RV_sb"){
    
    if(distribution=="norm"){
      
      alpha  <- param[1]
      beta  <- param[2]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2			<- param[3]
      
      TT    <- length(daily_ret)
      tau_d<-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ########structural breaks
      # 提取当前资产对的系数向量 (长度 1+n_n)
      
      n_n <- length(X_stru_break_date_vec)
      full_days <- as.Date(stats::time(daily_ret))
      
      curr_m_vec <- param[4:(n_n+4)]  
      curr_theta_vec <- param[(n_n+5):(2*n_n+5)] 
      
      mu			<- param[(2*n_n+6)]
      
      epsilon <- daily_ret - mu
      
      D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
      D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
      
      for (i in 1:n_n) {
        # 找到处于第 i 个节点之后的日期，设为 1
        D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
      }
      
      
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      theta_t_series <- as.vector(curr_theta_vec %*% D_mat)
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
      tau_d <-suppressWarnings(roll_sum_func(RV_m, c(K+1),weights = betas)) 
      tau_d<-exp(m_t_series+theta_t_series*tau_d[(K+1),])	 
      
      ll <- fast_garch_ll_cpp(alpha=alpha, beta=beta, mu=mu, 
                              epsilon=as.numeric(epsilon), tau_d=as.numeric(tau_d), daily_ret=as.numeric(daily_ret))
      
    } else if(distribution=="std"){
      
      alpha  <- param[1]
      beta  <- param[2]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2			<- param[3]
      v   <- param[4]
      
      TT    <- length(daily_ret)
      tau_d<-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ########structural breaks
      # 提取当前资产对的系数向量 (长度 1+n_n)
      
      n_n <- length(X_stru_break_date_vec)
      full_days <- as.Date(stats::time(daily_ret))
      
      curr_m_vec <- param[5:(n_n+5)]  
      curr_theta_vec <- param[(n_n+6):(2*n_n+6)] 
      
      mu			<- param[(2*n_n+7)]
      
      
      epsilon <- daily_ret - mu
      
      D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
      D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
      
      for (i in 1:n_n) {
        # 找到处于第 i 个节点之后的日期，设为 1
        D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
      }
      
      
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      theta_t_series <- as.vector(curr_theta_vec %*% D_mat)
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
      tau_d <-suppressWarnings(roll_sum_func(RV_m, c(K+1),weights = betas)) 
      tau_d<-exp(m_t_series+theta_t_series*tau_d[(K+1),])	
      
      ll <- fast_t_garch_ll(alpha=alpha, beta=beta, v=v, 
                            epsilon=as.numeric(epsilon), tau_d=as.numeric(tau_d), daily_ret=as.numeric(daily_ret))
      
      
    }
    
  } else if (model=="GM_noskew_X_sb"){
    
    if(distribution=="norm"){
      
      alpha  <- param[1]
      beta  <- param[2]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2			<- param[3]
      
      TT    <- length(daily_ret)
      tau_d<-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ########structural breaks
      # 提取当前资产对的系数向量 (长度 1+n_n)
      
      n_n <- length(X_stru_break_date_vec)
      full_days <- as.Date(stats::time(daily_ret))
      
      curr_m_vec <- param[4:(n_n+4)]  
      curr_theta_vec <- param[(n_n+5):(2*n_n+5)] 
      
      mu			<- param[(2*n_n+6)]
      
      
      epsilon <- daily_ret - mu
      
      D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
      D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
      
      for (i in 1:n_n) {
        # 找到处于第 i 个节点之后的日期，设为 1
        D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
      }
      
      
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      theta_t_series <- as.vector(curr_theta_vec %*% D_mat)
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
      tau_d <-suppressWarnings(roll_sum_func(mv_m, c(K+1),weights = betas)) 
      tau_d<-exp(m_t_series+theta_t_series*tau_d[(K+1),])	 
      
      ll <- fast_garch_ll_cpp(alpha=alpha, beta=beta, mu=mu, 
                              epsilon=as.numeric(epsilon), tau_d=as.numeric(tau_d), daily_ret=as.numeric(daily_ret))
      
    } else if(distribution=="std"){
      
      alpha  <- param[1]
      beta  <- param[2]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2			<- param[3]
      v   <- param[4]
      
      TT    <- length(daily_ret)
      tau_d<-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ########structural breaks
      # 提取当前资产对的系数向量 (长度 1+n_n)
      
      n_n <- length(X_stru_break_date_vec)
      full_days <- as.Date(stats::time(daily_ret))
      
      curr_m_vec <- param[5:(n_n+5)]  
      curr_theta_vec <- param[(n_n+6):(2*n_n+6)] 
      
      mu			<- param[(2*n_n+7)]
      
      
      epsilon <- daily_ret - mu
      
      D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
      D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
      
      for (i in 1:n_n) {
        # 找到处于第 i 个节点之后的日期，设为 1
        D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
      }
      
      
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      theta_t_series <- as.vector(curr_theta_vec %*% D_mat)
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
      tau_d <-suppressWarnings(roll_sum_func(mv_m, c(K+1),weights = betas)) 
      tau_d<-exp(m_t_series+theta_t_series*tau_d[(K+1),])	
      
      ll <- fast_t_garch_ll(alpha=alpha, beta=beta, v=v, 
                            epsilon=as.numeric(epsilon), tau_d=as.numeric(tau_d), daily_ret=as.numeric(daily_ret))
      
      
    }
    
    
  } else if (model=="GM_noskew_RV_X_sb"){
    
    if(distribution=="norm"){
      
      alpha  <- param[1]
      beta  <- param[2]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2_RV			<- param[3]
      w2_X			<- param[4]
      
      TT    <- length(daily_ret)
      tau_d_RV <-rep(NA,TT)
      tau_d_X <-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ########structural breaks
      # 提取当前资产对的系数向量 (长度 1+n_n)
      
      n_n <- length(X_stru_break_date_vec)
      full_days <- as.Date(stats::time(daily_ret))
      
      curr_m_vec <- param[5:(n_n+5)]   
      curr_theta_RV_vec <- param[(n_n+6):(2*n_n+6)]
      curr_theta_X_vec <- param[(2*n_n+7):(3*n_n+7)]
      
      mu			<- param[(3*n_n+8)]
      
      
      epsilon <- daily_ret - mu
      
      D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
      D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
      
      for (i in 1:n_n) {
        # 找到处于第 i 个节点之后的日期，设为 1
        D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
      }
      
      
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      theta_RV_t_series <- as.vector(curr_theta_RV_vec %*% D_mat)
      theta_X_t_series <- as.vector(curr_theta_X_vec %*% D_mat)
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      
      betas_RV<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2_RV))[2:(K+1)],0)
      tau_d_RV <-suppressWarnings(roll_sum_func(RV_m, c(K+1),weights = betas_RV)) 
      tau_d_RV<-tau_d_RV[(K+1),]	 
      
      betas_X<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2_X))[2:(K+1)],0)
      tau_d_X <-suppressWarnings(roll_sum_func(mv_m, c(K+1),weights = betas_X)) 
      tau_d_X<-tau_d_X[(K+1),]
      
      tau_d <- exp(m_t_series+ theta_RV_t_series*tau_d_RV + theta_X_t_series*tau_d_X)
      
      ll <- fast_garch_ll_cpp(alpha=alpha, beta=beta, mu=mu, 
                              epsilon=as.numeric(epsilon), tau_d=as.numeric(tau_d), daily_ret=as.numeric(daily_ret))
      
    } else if(distribution=="std"){
      
      alpha  <- param[1]
      beta  <- param[2]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2_RV			<- param[3]
      w2_X			<- param[4]
      v   <- param[5]
      
      TT    <- length(daily_ret)
      tau_d_RV <-rep(NA,TT)
      tau_d_X <-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ########structural breaks
      # 提取当前资产对的系数向量 (长度 1+n_n)
      
      n_n <- length(X_stru_break_date_vec)
      full_days <- as.Date(stats::time(daily_ret))
      
      curr_m_vec <- param[6:(n_n+6)]  
      curr_theta_RV_vec <- param[(n_n+7):(2*n_n+7)] 
      curr_theta_X_vec <- param[(2*n_n+8):(3*n_n+8)]
      
      mu			<- param[(3*n_n+9)]
      
      
      epsilon <- daily_ret - mu
      
      D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
      D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
      
      for (i in 1:n_n) {
        # 找到处于第 i 个节点之后的日期，设为 1
        D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
      }
      
      
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      theta_RV_t_series <- as.vector(curr_theta_RV_vec %*% D_mat)
      theta_X_t_series <- as.vector(curr_theta_X_vec %*% D_mat)
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      
      betas_RV<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2_RV))[2:(K+1)],0)
      tau_d_RV <-suppressWarnings(roll_sum_func(RV_m, c(K+1),weights = betas_RV)) 
      tau_d_RV<-tau_d_RV[(K+1),]	 
      
      betas_X<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2_X))[2:(K+1)],0)
      tau_d_X <-suppressWarnings(roll_sum_func(mv_m, c(K+1),weights = betas_X)) 
      tau_d_X<-tau_d_X[(K+1),]
      
      tau_d <- exp(m_t_series+ theta_RV_t_series*tau_d_RV + theta_X_t_series*tau_d_X)
      
      
      ll <- fast_t_garch_ll(alpha=alpha, beta=beta, v=v, 
                            epsilon=as.numeric(epsilon), tau_d=as.numeric(tau_d), daily_ret=as.numeric(daily_ret))
      
      
      
    }
    
    
  } 
  
  return(ll)
  ##### end
  
}

#est
garchmidas_noskew_est <- function(param,model,daily_ret,RV_m=NULL,mv_m=NULL,K,distribution,lag_fun="Beta",X_stru_break_date_vec=NULL){
  
  
  
  if (model=="GM_noskew_RV"){
    
    if(distribution=="norm"){
      
      alpha  <- param[1]
      beta  <- param[2]
      m   <- param[3]
      theta		<- param[4]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2			<- param[5]
      
      mu			<- param[6]
      
      
      epsilon <- daily_ret - mu
      
      TT    <- length(daily_ret)
      tau_d<-rep(NA,TT)
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
      tau_d <-m+theta*suppressWarnings(roll_sum_func(RV_m, c(K+1),weights = betas)) 
      tau_d<-exp(tau_d[(K+1),])	 
      
      ####### short-run 
      step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      step_1[is.na(step_1)] <- 0
      x_step_1 <- c(0, step_1[-TT])
      x_step_1[2] <- x_step_1[2] + beta * 1
      
      g_it <- as.numeric(stats::filter(x_step_1, filter = beta, method = "recursive"))
      g_it[1] <- 1
      
      #old and slow version
      # step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      # for(i in 2:TT){
      #   g_it[i]  <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
      # }
      
    } else if(distribution=="std"){
      alpha <- param[1]
      beta  <- param[2]
      m     <- param[3]
      theta		<- param[4]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2			<- param[5]
      v			<- param[6]
      
      mu			<- param[7]
      
      
      epsilon <- daily_ret - mu
      
      TT     <- length(daily_ret)
      tau_d<-rep(NA,TT)
      g_it  <- rep(1,TT)            # daily conditional variance 
      ll       <- 0 
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
      
      tau_d   <-  m+theta*suppressWarnings(roll_sum_func(RV_m, c(K+1),weights = betas)) 
      tau_d   <-  exp(tau_d[(K+1),])	 
      
      ####### short-run 
      step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      step_1[is.na(step_1)] <- 0
      x_step_1 <- c(0, step_1[-TT])
      x_step_1[2] <- x_step_1[2] + beta * 1
      
      g_it <- as.numeric(stats::filter(x_step_1, filter = beta, method = "recursive"))
      g_it[1] <- 1
      
      #old and slow version
      # step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      # 
      # for(i in 2:TT){
      #   g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
      # }
      
    }
    
    ###### variance 
    h_it_2 <- zoo::coredata(g_it*tau_d)
    g_it <- zoo::coredata(g_it)
    tau_d <- zoo::coredata(tau_d)
    
    
    h_it_2 <- as.xts(h_it_2,index(daily_ret))
    g_it <- as.xts(g_it,index(daily_ret))
    tau_d <- as.xts(tau_d,index(daily_ret))
    
    
    results<-list(
      "h_it_2" = h_it_2,
      "g_it" = g_it,
      "tau_d" = tau_d
    )
    return(results)
    
  } else if (model=="GM_noskew_X"){
    
    if(distribution=="norm"){
      
      alpha  <- param[1]
      beta  <- param[2]
      m   <- param[3]
      theta		<- param[4]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2			<- param[5]
      
      mu			<- param[6]
      
      
      epsilon <- daily_ret - mu
      
      TT    <- length(daily_ret)
      tau_d<-rep(NA,TT)
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
      tau_d <-m+theta*suppressWarnings(roll_sum_func(mv_m, c(K+1),weights = betas)) 
      tau_d<-exp(tau_d[(K+1),])	 
      
      ####### short-run 
      step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      step_1[is.na(step_1)] <- 0
      x_step_1 <- c(0, step_1[-TT])
      x_step_1[2] <- x_step_1[2] + beta * 1
      
      g_it <- as.numeric(stats::filter(x_step_1, filter = beta, method = "recursive"))
      g_it[1] <- 1
      
      #old and slow version
      # step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      # for(i in 2:TT){
      #   g_it[i]  <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
      # }
      
      
    } else if(distribution=="std"){
      alpha <- param[1]
      beta  <- param[2]
      m     <- param[3]
      theta		<- param[4]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2			<- param[5]
      v			<- param[6]
      
      mu			<- param[7]
      
      
      epsilon <- daily_ret - mu
      
      TT     <- length(daily_ret)
      tau_d<-rep(NA,TT)
      g_it  <- rep(1,TT)            # daily conditional variance 
      ll       <- 0 
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
      
      tau_d   <-  m+theta*suppressWarnings(roll_sum_func(mv_m, c(K+1),weights = betas)) 
      tau_d   <-  exp(tau_d[(K+1),])	 
      
      ####### short-run 
      step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      step_1[is.na(step_1)] <- 0
      x_step_1 <- c(0, step_1[-TT])
      x_step_1[2] <- x_step_1[2] + beta * 1
      
      g_it <- as.numeric(stats::filter(x_step_1, filter = beta, method = "recursive"))
      g_it[1] <- 1
      
      #old and slow version
      # step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      # 
      # for(i in 2:TT){
      #   g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
      # }
      
      
    }
    
    ###### variance 
    h_it_2 <- zoo::coredata(g_it*tau_d)
    g_it <- zoo::coredata(g_it)
    tau_d <- zoo::coredata(tau_d)
    
    h_it_2 <- as.xts(h_it_2,index(daily_ret))
    g_it <- as.xts(g_it,index(daily_ret))
    tau_d <- as.xts(tau_d,index(daily_ret))
    
    results<-list(
      "h_it_2" = h_it_2,
      "g_it" = g_it,
      "tau_d" = tau_d
    )
    return(results)
    
  } else if (model=="GM_noskew_RV_X"){
    
    if(distribution=="norm"){
      
      alpha  <- param[1]
      beta  <- param[2]
      m   <- param[3]
      theta_RV		<- param[4]
      theta_X		<- param[5]
      
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2_RV			<- param[6]
      w2_X			<- param[7]
      
      mu			<- param[8]
      
      
      epsilon <- daily_ret - mu
      
      TT    <- length(daily_ret)
      
      tau_d_RV <-rep(NA,TT)
      tau_d_X <-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      
      betas_RV<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2_RV))[2:(K+1)],0)
      tau_d_RV <-suppressWarnings(roll_sum_func(RV_m, c(K+1),weights = betas_RV)) 
      tau_d_RV<-tau_d_RV[(K+1),]	 
      
      betas_X<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2_X))[2:(K+1)],0)
      tau_d_X <-suppressWarnings(roll_sum_func(mv_m, c(K+1),weights = betas_X)) 
      tau_d_X<-tau_d_X[(K+1),]
      
      tau_d <- exp(m + theta_RV*tau_d_RV + theta_X*tau_d_X)
      
      ####### short-run 
      step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      step_1[is.na(step_1)] <- 0
      x_step_1 <- c(0, step_1[-TT])
      x_step_1[2] <- x_step_1[2] + beta * 1
      
      g_it <- as.numeric(stats::filter(x_step_1, filter = beta, method = "recursive"))
      g_it[1] <- 1
      
      #old and slow version
      # step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      # for(i in 2:TT){
      #   g_it[i]  <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
      # }
      
      
    } else if(distribution=="std"){
      alpha  <- param[1]
      beta  <- param[2]
      m   <- param[3]
      theta_RV		<- param[4]
      theta_X		<- param[5]
      
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2_RV			<- param[6]
      w2_X			<- param[7]
      
      v			<-   param[8]
      
      mu			<- param[9]
      
      
      epsilon <- daily_ret - mu
      
      TT    <- length(daily_ret)
      tau_d<-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      
      betas_RV<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2_RV))[2:(K+1)],0)
      tau_d_RV <-suppressWarnings(roll_sum_func(RV_m, c(K+1),weights = betas_RV)) 
      tau_d_RV<-tau_d_RV[(K+1),]	 
      
      betas_X<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2_X))[2:(K+1)],0)
      tau_d_X <-suppressWarnings(roll_sum_func(mv_m, c(K+1),weights = betas_X)) 
      tau_d_X<-tau_d_X[(K+1),]
      
      tau_d <- exp(m + theta_RV*tau_d_RV + theta_X*tau_d_X)
      
      ####### short-run 
      step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      step_1[is.na(step_1)] <- 0
      x_step_1 <- c(0, step_1[-TT])
      x_step_1[2] <- x_step_1[2] + beta * 1
      
      g_it <- as.numeric(stats::filter(x_step_1, filter = beta, method = "recursive"))
      g_it[1] <- 1
      
      #old and slow version
      # step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      # 
      # for(i in 2:TT){
      #   g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
      # }
      
      
    }
    
    ###### variance 
    h_it_2 <- zoo::coredata(g_it*tau_d)
    g_it <- zoo::coredata(g_it)
    tau_d <- zoo::coredata(tau_d)
    tau_d_RV <-zoo::coredata(tau_d_RV)
    tau_d_X <-zoo::coredata(tau_d_X)
    
    h_it_2 <- as.xts(h_it_2,index(daily_ret))
    g_it <- as.xts(g_it,index(daily_ret))
    tau_d <- as.xts(tau_d,index(daily_ret))
    tau_d_RV <-as.xts(tau_d_RV,index(daily_ret))
    tau_d_X <-as.xts(tau_d_X,index(daily_ret))
    
    results<-list(
      "h_it_2" = h_it_2,
      "g_it" = g_it,
      "tau_d" = tau_d,
      "tau_d_RV" = tau_d_RV,
      "tau_d_X" = tau_d_X
    )
    return(results)
    
  } else if (model=="GM_noskew_RV_sb"){
    
    if(distribution=="norm"){
      
      alpha  <- param[1]
      beta  <- param[2]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2			<- param[3]
      
      TT    <- length(daily_ret)
      tau_d<-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ########structural breaks
      # 提取当前资产对的系数向量 (长度 1+n_n)
      
      n_n <- length(X_stru_break_date_vec)
      full_days <- as.Date(stats::time(daily_ret))
      
      curr_m_vec <- param[4:(n_n+4)]  
      curr_theta_vec <- param[(n_n+5):(2*n_n+5)] 
      
      mu			<- param[(2*n_n+6)]
      
      
      epsilon <- daily_ret - mu
      
      D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
      D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
      
      for (i in 1:n_n) {
        # 找到处于第 i 个节点之后的日期，设为 1
        D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
      }
      
      
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      theta_t_series <- as.vector(curr_theta_vec %*% D_mat)
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
      tau_d <-suppressWarnings(roll_sum_func(RV_m, c(K+1),weights = betas)) 
      tau_d<-exp(m_t_series+theta_t_series*tau_d[(K+1),])	 
      
      ####### short-run 
      step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      step_1[is.na(step_1)] <- 0
      x_step_1 <- c(0, step_1[-TT])
      x_step_1[2] <- x_step_1[2] + beta * 1
      
      g_it <- as.numeric(stats::filter(x_step_1, filter = beta, method = "recursive"))
      g_it[1] <- 1
      
      #old and slow version
      # step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      # for(i in 2:TT){
      #   g_it[i]  <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
      # }
      
      
    } else if(distribution=="std"){
      
      alpha  <- param[1]
      beta  <- param[2]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2			<- param[3]
      v   <- param[4]
      
      TT    <- length(daily_ret)
      tau_d<-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ########structural breaks
      # 提取当前资产对的系数向量 (长度 1+n_n)
      
      n_n <- length(X_stru_break_date_vec)
      full_days <- as.Date(stats::time(daily_ret))
      
      curr_m_vec <- param[5:(n_n+5)]  
      curr_theta_vec <- param[(n_n+6):(2*n_n+6)] 
      
      mu			<- param[(2*n_n+7)]
      
      
      epsilon <- daily_ret - mu
      
      D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
      D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
      
      for (i in 1:n_n) {
        # 找到处于第 i 个节点之后的日期，设为 1
        D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
      }
      
      
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      theta_t_series <- as.vector(curr_theta_vec %*% D_mat)
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
      tau_d <-suppressWarnings(roll_sum_func(RV_m, c(K+1),weights = betas)) 
      tau_d<-exp(m_t_series+theta_t_series*tau_d[(K+1),])	
      
      ####### short-run 
      step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      step_1[is.na(step_1)] <- 0
      x_step_1 <- c(0, step_1[-TT])
      x_step_1[2] <- x_step_1[2] + beta * 1
      
      g_it <- as.numeric(stats::filter(x_step_1, filter = beta, method = "recursive"))
      g_it[1] <- 1
      
      #old and slow version
      # step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      # 
      # for(i in 2:TT){
      #   g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
      # }
      
      
    }
    
    ###### variance 
    h_it_2 <- zoo::coredata(g_it*tau_d)
    g_it <- zoo::coredata(g_it)
    tau_d <- zoo::coredata(tau_d)
    
    h_it_2 <- as.xts(h_it_2,index(daily_ret))
    g_it <- as.xts(g_it,index(daily_ret))
    tau_d <- as.xts(tau_d,index(daily_ret))
    
    
    results<-list(
      "h_it_2" = h_it_2,
      "g_it" = g_it,
      "tau_d" = tau_d
    )
    return(results)
    
  } else if (model=="GM_noskew_X_sb"){
    
    if(distribution=="norm"){
      
      alpha  <- param[1]
      beta  <- param[2]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2			<- param[3]
      
      TT    <- length(daily_ret)
      tau_d<-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ########structural breaks
      # 提取当前资产对的系数向量 (长度 1+n_n)
      
      n_n <- length(X_stru_break_date_vec)
      full_days <- as.Date(stats::time(daily_ret))
      
      curr_m_vec <- param[4:(n_n+4)]  
      curr_theta_vec <- param[(n_n+5):(2*n_n+5)] 
      
      mu			<- param[(2*n_n+6)]
      
      
      epsilon <- daily_ret - mu
      
      D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
      D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
      
      for (i in 1:n_n) {
        # 找到处于第 i 个节点之后的日期，设为 1
        D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
      }
      
      
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      theta_t_series <- as.vector(curr_theta_vec %*% D_mat)
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
      tau_d <-suppressWarnings(roll_sum_func(mv_m, c(K+1),weights = betas)) 
      tau_d<-exp(m_t_series+theta_t_series*tau_d[(K+1),])	 
      
      ####### short-run 
      step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      step_1[is.na(step_1)] <- 0
      x_step_1 <- c(0, step_1[-TT])
      x_step_1[2] <- x_step_1[2] + beta * 1
      
      g_it <- as.numeric(stats::filter(x_step_1, filter = beta, method = "recursive"))
      g_it[1] <- 1
      
      #old and slow version
      # step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      # for(i in 2:TT){
      #   g_it[i]  <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
      # }
      
      
    } else if(distribution=="std"){
      
      alpha  <- param[1]
      beta  <- param[2]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2			<- param[3]
      v   <- param[4]
      
      TT    <- length(daily_ret)
      tau_d<-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ########structural breaks
      # 提取当前资产对的系数向量 (长度 1+n_n)
      
      n_n <- length(X_stru_break_date_vec)
      full_days <- as.Date(stats::time(daily_ret))
      
      curr_m_vec <- param[5:(n_n+5)]  
      curr_theta_vec <- param[(n_n+6):(2*n_n+6)] 
      
      mu			<- param[(2*n_n+7)]
      
      
      epsilon <- daily_ret - mu
      
      D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
      D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
      
      for (i in 1:n_n) {
        # 找到处于第 i 个节点之后的日期，设为 1
        D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
      }
      
      
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      theta_t_series <- as.vector(curr_theta_vec %*% D_mat)
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
      tau_d <-suppressWarnings(roll_sum_func(mv_m, c(K+1),weights = betas)) 
      tau_d<-exp(m_t_series+theta_t_series*tau_d[(K+1),])	
      
      ####### short-run 
      step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      step_1[is.na(step_1)] <- 0
      x_step_1 <- c(0, step_1[-TT])
      x_step_1[2] <- x_step_1[2] + beta * 1
      
      g_it <- as.numeric(stats::filter(x_step_1, filter = beta, method = "recursive"))
      g_it[1] <- 1
      
      #old and slow version
      # step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      # 
      # for(i in 2:TT){
      #   g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
      # }
      
      
      
    }
    
    ###### variance 
    h_it_2 <- zoo::coredata(g_it*tau_d)
    g_it <- zoo::coredata(g_it)
    tau_d <- zoo::coredata(tau_d)
    
    
    h_it_2 <- as.xts(h_it_2,index(daily_ret))
    g_it <- as.xts(g_it,index(daily_ret))
    tau_d <- as.xts(tau_d,index(daily_ret))
    
    
    results<-list(
      "h_it_2" = h_it_2,
      "g_it" = g_it,
      "tau_d" = tau_d
    )
    return(results)
    
    
  } else if (model=="GM_noskew_RV_X_sb"){
    
    if(distribution=="norm"){
      
      alpha  <- param[1]
      beta  <- param[2]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2_RV			<- param[3]
      w2_X			<- param[4]
      
      TT    <- length(daily_ret)
      tau_d_RV <-rep(NA,TT)
      tau_d_X <-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ########structural breaks
      # 提取当前资产对的系数向量 (长度 1+n_n)
      
      n_n <- length(X_stru_break_date_vec)
      full_days <- as.Date(stats::time(daily_ret))
      
      curr_m_vec <- param[5:(n_n+5)]   
      curr_theta_RV_vec <- param[(n_n+6):(2*n_n+6)]
      curr_theta_X_vec <- param[(2*n_n+7):(3*n_n+7)]
      
      mu			<- param[(3*n_n+8)]
      
      
      epsilon <- daily_ret - mu
      
      D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
      D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
      
      for (i in 1:n_n) {
        # 找到处于第 i 个节点之后的日期，设为 1
        D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
      }
      
      
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      theta_RV_t_series <- as.vector(curr_theta_RV_vec %*% D_mat)
      theta_X_t_series <- as.vector(curr_theta_X_vec %*% D_mat)
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      
      betas_RV<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2_RV))[2:(K+1)],0)
      tau_d_RV <-suppressWarnings(roll_sum_func(RV_m, c(K+1),weights = betas_RV)) 
      tau_d_RV<-tau_d_RV[(K+1),]	 
      
      betas_X<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2_X))[2:(K+1)],0)
      tau_d_X <-suppressWarnings(roll_sum_func(mv_m, c(K+1),weights = betas_X)) 
      tau_d_X<-tau_d_X[(K+1),]
      
      tau_d <- exp(m_t_series+ theta_RV_t_series*tau_d_RV + theta_X_t_series*tau_d_X)
      
      ####### short-run 
      step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      step_1[is.na(step_1)] <- 0
      x_step_1 <- c(0, step_1[-TT])
      x_step_1[2] <- x_step_1[2] + beta * 1
      
      g_it <- as.numeric(stats::filter(x_step_1, filter = beta, method = "recursive"))
      g_it[1] <- 1
      
      #old and slow version
      # step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      # for(i in 2:TT){
      #   g_it[i]  <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
      # }
      
      
    } else if(distribution=="std"){
      
      alpha  <- param[1]
      beta  <- param[2]
      w1			<- ifelse(lag_fun=="Beta",1,0)
      w2_RV			<- param[3]
      w2_X			<- param[4]
      v   <- param[5]
      
      TT    <- length(daily_ret)
      tau_d_RV <-rep(NA,TT)
      tau_d_X <-rep(NA,TT)
      
      g_it  <- rep(1,TT)     # daily conditional variance 
      ll <- 0 
      
      ########structural breaks
      # 提取当前资产对的系数向量 (长度 1+n_n)
      
      n_n <- length(X_stru_break_date_vec)
      full_days <- as.Date(stats::time(daily_ret))
      
      curr_m_vec <- param[6:(n_n+6)]  
      curr_theta_RV_vec <- param[(n_n+7):(2*n_n+7)] 
      curr_theta_X_vec <- param[(2*n_n+8):(3*n_n+8)]
      
      mu			<- param[(3*n_n+9)]
      
      
      epsilon <- daily_ret - mu
      
      D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
      D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
      
      for (i in 1:n_n) {
        # 找到处于第 i 个节点之后的日期，设为 1
        D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
      }
      
      
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      theta_RV_t_series <- as.vector(curr_theta_RV_vec %*% D_mat)
      theta_X_t_series <- as.vector(curr_theta_X_vec %*% D_mat)
      
      ##### daily long-run
      
      weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
      
      betas_RV<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2_RV))[2:(K+1)],0)
      tau_d_RV <-suppressWarnings(roll_sum_func(RV_m, c(K+1),weights = betas_RV)) 
      tau_d_RV<-tau_d_RV[(K+1),]	 
      
      betas_X<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2_X))[2:(K+1)],0)
      tau_d_X <-suppressWarnings(roll_sum_func(mv_m, c(K+1),weights = betas_X)) 
      tau_d_X<-tau_d_X[(K+1),]
      
      tau_d <- exp(m_t_series+ theta_RV_t_series*tau_d_RV + theta_X_t_series*tau_d_X)
      
      ####### short-run 
      step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      step_1[is.na(step_1)] <- 0
      x_step_1 <- c(0, step_1[-TT])
      x_step_1[2] <- x_step_1[2] + beta * 1
      
      g_it <- as.numeric(stats::filter(x_step_1, filter = beta, method = "recursive"))
      g_it[1] <- 1
      
      #old and slow version
      # step_1<-(1-alpha-beta)+(alpha)*(epsilon)^2/tau_d
      # 
      # for(i in 2:TT){
      #   g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
      # }
      
    }
    
    ###### variance 
    h_it_2 <- zoo::coredata(g_it*tau_d)
    g_it <- zoo::coredata(g_it)
    tau_d <- zoo::coredata(tau_d)
    tau_d_RV <-zoo::coredata(tau_d_RV)
    tau_d_X <-zoo::coredata(tau_d_X)
    
    h_it_2 <- as.xts(h_it_2,index(daily_ret))
    g_it <- as.xts(g_it,index(daily_ret))
    tau_d <- as.xts(tau_d,index(daily_ret))
    tau_d_RV <-as.xts(tau_d_RV,index(daily_ret))
    tau_d_X <-as.xts(tau_d_X,index(daily_ret))
    
    results<-list(
      "h_it_2" = h_it_2,
      "g_it" = g_it,
      "tau_d" = tau_d,
      "tau_d_RV" = tau_d_RV,
      "tau_d_X" = tau_d_X
    )
    return(results)
    
    
  } 
  
  
  ##### end
  
}

#fit
garchmidas_noskew_fit <- function(model,daily_ret,RV_m=NULL,mv_m=NULL,K,distribution,lag_fun="Beta",X_stru_break_date_vec=NULL,
                                  out_of_sample=NULL,vol_proxy=NULL,R=NULL){
  
  ############################# check on valid choices 
  
  if((model != "GM_noskew_RV")&(model != "GM_noskew_X")&(model != "GM_noskew_RV_X")&(model != "GM_noskew_RV_sb")&(model != "GM_noskew_X_sb")&(model != "GM_noskew_RV_X_sb")) 
  { stop(cat("#Warning:\n Valid choices for the parameter 'model' are 'GM_noskew_RV'/'GM_noskew_X'/'GM_noskew_RV_X'/'GM_noskew_RV_sb'/'GM_noskew_X_sb'/'GM_noskew_RV_X_sb' \n"))}
  if((model == "GM_noskew_RV"|model == "GM_noskew_RV_X"|model == "GM_noskew_RV_sb"|model == "GM_noskew_RV_X_sb")&(is.null(RV_m))) 
  { stop(cat("#Warning:\n If the model chosen includes the 'RV' term, then the 'RV_m' variable has to be provided \n"))}
  if((model == "GM_noskew_X"|model == "GM_noskew_RV_X"|model == "GM_noskew_X_sb"|model == "GM_noskew_RV_X_sb")&(is.null(mv_m))) 
  { stop(cat("#Warning:\n If the model chosen includes the 'X(Macro Varibales)' term, then the 'mv_m' variable has to be provided \n"))}
  if((model == "GM_noskew_RV_sb"|model == "GM_noskew_X_sb"|model == "GM_noskew_RV_X_sb")&(is.null(X_stru_break_date_vec))) 
  { stop(cat("#Warning:\n If the model chosen includes the 'Structural Breaks' term, then the 'X_stru_break_date_vec' variable has to be provided \n"))}
  if((distribution != "norm")&(distribution != "std")) 
  { stop(cat("#Warning:\n Valid choices for the parameter 'distribution' are 'norm' and 'std' \n"))}
  if((lag_fun != "Beta")&(lag_fun!= "Almon")) 
  { stop(cat("#Warning:\n Valid choices for the parameter 'lag_fun' are 'Beta' and 'Almon' \n"))}
  
  if (is.null(R)) {
    
    R <- 100
  }
  
  
  N<-length(daily_ret)
  cond_r_t<- class(daily_ret)[1]
  
  if(cond_r_t != "xts") { stop(
    cat("#Warning:\n Parameter 'daily_ret' must be an xts object. Please provide it in the correct form \n")
  )}
  
  if(!is.null(RV_m) && !is.matrix(RV_m)) { stop(
    cat("#Warning:\n Parameter 'RV_m' must be a matrix. Please provide it in the correct form \n")
  )}
  
  if(!is.null(mv_m) && !is.matrix(mv_m)) { stop(
    cat("#Warning:\n Parameter 'mv_m' must be a matrix. Please provide it in the correct form \n")
  )}
  
  if(!is.null(RV_m) && dim(RV_m)[2] != N) { stop(
    cat("#Warning:\n The columns of the matrix 'RV_m' must be equal to the length of vector 'daily_ret'. Please provide it in the correct form \n")
  )}
  
  if(!is.null(mv_m) && dim(mv_m)[2] != N) { stop(
    cat("#Warning:\n The columns of the matrix 'mv_m' must be equal to the length of vector 'daily_ret'. Please provide it in the correct form \n")
  )}
  
  ############## check if the vol_proxy is provided and if it has the same time span of daily_ret
  
  if(!is.null(vol_proxy)){ 
    if(any(range(time(daily_ret))!=range(time(vol_proxy)))){
      stop(
        cat("#Warning:\n The vector 'vol_proxy' has to be observed during the same time span of 'daily_ret' \n")
      )}}
  
  ############## in sample and out of sample
  
  if (is.null(out_of_sample)){ # out_of_sample parameter is missing
    
    r_t_in_s<-daily_ret
    
    if(!is.null(RV_m)){
      RV_m_in_s<-RV_m
    }
    
    if(!is.null(mv_m)){
      mv_m_in_s<-mv_m
    }
    
    if(is.null(vol_proxy)){
      
      vol_proxy_in_s<-r_t_in_s^2
      
    } else {
      
      vol_proxy_in_s<-vol_proxy
    }
    
    
  } else {					# out_of_sample parameter is present
    
    r_t_in_s<-daily_ret[1:(N-out_of_sample)]
    
    if(!is.null(RV_m)){
      RV_m_in_s<-RV_m[,1:(N-out_of_sample)]
    }
    
    if(!is.null(mv_m)){
      mv_m_in_s<-mv_m[,1:(N-out_of_sample)]
    }
    
    if(is.null(vol_proxy)){
      
      vol_proxy_in_s<-r_t_in_s^2
      
    } else {
      
      vol_proxy_in_s<-vol_proxy[1:(N-out_of_sample)]
    }
    
    
  }
  
  #########estimate for each model############
  ##GM_noskew_RV/GM_noskew_X/GM_noskew_RV_X/GM_noskew_RV_sb/GM_noskew_X_sb/GM_noskew_RV_X_sb
  
  if(model=="GM_noskew_RV"&distribution=="norm"&lag_fun=="Beta"){
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=6)
    colnames(begin_val)<-c("alpha","beta","m","theta","w2","mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-stats::runif(R,min=-1,max=1)
    begin_val[,4]<-stats::runif(R,min=-1,max=1)
    begin_val[,5]<-1.01
    begin_val[,6]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,
                                                 RV_m=RV_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,0,0),       	 			## alpha>0.001
      c(0,1,0,0,0,0),        		 		## beta>0.001
      c(-1,-1,0,0,0,0),	     				## alpha+beta<1
      c(0,0,0,0,1,0))				 	## w2>1.001
    
    ci<-c(-0.001,-0.001,0.999,-1.001)
    
    # numerical_gradient <- function(param,model,daily_ret,RV_m,K,distribution,lag_fun) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      RV_m=RV_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      RV_m=RV_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun)
    
  } else if(model=="GM_noskew_RV"&distribution=="norm"&lag_fun=="Almon"){
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=6)
    colnames(begin_val)<-c("alpha","beta","m","theta","w2","mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-stats::runif(R,min=-1,max=1)
    begin_val[,4]<-stats::runif(R,min=-1,max=1)
    begin_val[,5]<--0.1
    begin_val[,6]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,
                                                 RV_m=RV_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,0,0),       	 			## alpha>0.001
      c(0,1,0,0,0,0),        		 		## beta>0.001
      c(-1,-1,0,0,0,0),	     				## alpha+beta<1
      c(0,0,0,0,-1,0))				 	## w2<0
    
    ci<-c(-0.001,-0.001,0.999,-0.001)
    
    # numerical_gradient <- function(param,model,daily_ret,RV_m,K,distribution,lag_fun) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      RV_m=RV_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      RV_m=RV_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun)
    
  }else if(model=="GM_noskew_RV"&distribution=="std"&lag_fun=="Beta"){
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=7)
    colnames(begin_val)<-c("alpha","beta","m","theta","w2","v","mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-stats::runif(R,min=-1,max=1)
    begin_val[,4]<-stats::runif(R,min=-1,max=1)
    begin_val[,5]<-1.01
    begin_val[,6]<-stats::runif(R,min=2.01,max=10)
    begin_val[,7]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,
                                                 RV_m=RV_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,0,0,0),       	 			## alpha>0.001
      c(0,1,0,0,0,0,0),        		 		## beta>0.001
      c(-1,-1,0,0,0,0,0),	     				## alpha+beta<1
      c(0,0,0,0,1,0,0),              	## w2>1.001
      c(0,0,0,0,0,1,0))				       ## shape>2.001
    
    ci<-c(-0.001,-0.001,0.999,-1.001,-2.001)
    
    # numerical_gradient <- function(param,model,daily_ret,RV_m,K,distribution,lag_fun) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      RV_m=RV_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      RV_m=RV_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun)
    
    
  }else if(model=="GM_noskew_RV"&distribution=="std"&lag_fun=="Almon"){
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=7)
    colnames(begin_val)<-c("alpha","beta","m","theta","w2","v","mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-stats::runif(R,min=-1,max=1)
    begin_val[,4]<-stats::runif(R,min=-1,max=1)
    begin_val[,5]<--0.1
    begin_val[,6]<-stats::runif(R,min=2.01,max=10)
    begin_val[,7]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,
                                                 RV_m=RV_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,0,0,0),       	 			## alpha>0.001
      c(0,1,0,0,0,0,0),        		 		## beta>0.001
      c(-1,-1,0,0,0,0,0),	     				## alpha+beta<1
      c(0,0,0,0,-1,0,0),             ## w2<0
      c(0,0,0,0,0,1,0) )				 	## shape>2.001
    
    ci<-c(-0.001,-0.001,0.999,-0.001,-2.001)
    
    # numerical_gradient <- function(param,model,daily_ret,RV_m,K,distribution,lag_fun) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      RV_m=RV_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      RV_m=RV_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun)
    
  }else if(model=="GM_noskew_X"&distribution=="norm"&lag_fun=="Beta"){
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=6)
    colnames(begin_val)<-c("alpha","beta","m","theta","w2","mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-stats::runif(R,min=-1,max=1)
    begin_val[,4]<-stats::runif(R,min=-1,max=1)
    begin_val[,5]<-1.01
    begin_val[,6]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,
                                                 mv_m=mv_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,0,0),       	 			## alpha>0.001
      c(0,1,0,0,0,0),        		 		## beta>0.001
      c(-1,-1,0,0,0,0),	     				## alpha+beta<1
      c(0,0,0,0,1,0))				 	## w2>1.001
    
    ci<-c(-0.001,-0.001,0.999,-1.001)
    
    # numerical_gradient <- function(param,model,daily_ret,mv_m,K,distribution,lag_fun) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      mv_m=mv_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      mv_m=mv_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun)
    
  } else if(model=="GM_noskew_X"&distribution=="norm"&lag_fun=="Almon"){
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=6)
    colnames(begin_val)<-c("alpha","beta","m","theta","w2","mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-stats::runif(R,min=-1,max=1)
    begin_val[,4]<-stats::runif(R,min=-1,max=1)
    begin_val[,5]<--0.1
    begin_val[,6]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,
                                                 mv_m=mv_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,0,0),       	 			## alpha>0.001
      c(0,1,0,0,0,0),        		 		## beta>0.001
      c(-1,-1,0,0,0,0),	     				## alpha+beta<1
      c(0,0,0,0,-1,0))				 	## w2<0
    
    ci<-c(-0.001,-0.001,0.999,-0.001)
    
    # numerical_gradient <- function(param,model,daily_ret,mv_m,K,distribution,lag_fun) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      mv_m=mv_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      mv_m=mv_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun)
    
  }else if(model=="GM_noskew_X"&distribution=="std"&lag_fun=="Beta"){
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=7)
    colnames(begin_val)<-c("alpha","beta","m","theta","w2","v","mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-stats::runif(R,min=-1,max=1)
    begin_val[,4]<-stats::runif(R,min=-1,max=1)
    begin_val[,5]<-1.01
    begin_val[,6]<-stats::runif(R,min=2.01,max=10)
    begin_val[,7]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,
                                                 mv_m=mv_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,0,0,0),       	 			## alpha>0.001
      c(0,1,0,0,0,0,0),        		 		## beta>0.001
      c(-1,-1,0,0,0,0,0),	     				## alpha+beta<1
      c(0,0,0,0,1,0,0),              	## w2>1.001
      c(0,0,0,0,0,1,0))				       ## shape>2.001
    
    ci<-c(-0.001,-0.001,0.999,-1.001,-2.001)
    
    # numerical_gradient <- function(param,model,daily_ret,mv_m,K,distribution,lag_fun) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      mv_m=mv_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      mv_m=mv_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun)
    
    
  }else if(model=="GM_noskew_X"&distribution=="std"&lag_fun=="Almon"){
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=7)
    colnames(begin_val)<-c("alpha","beta","m","theta","w2","v","mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-stats::runif(R,min=-1,max=1)
    begin_val[,4]<-stats::runif(R,min=-1,max=1)
    begin_val[,5]<--0.1
    begin_val[,6]<-stats::runif(R,min=2.01,max=10)
    begin_val[,7]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,
                                                 mv_m=mv_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,0,0,0),       	 			## alpha>0.001
      c(0,1,0,0,0,0,0),        		 		## beta>0.001
      c(-1,-1,0,0,0,0,0),	     				## alpha+beta<1
      c(0,0,0,0,-1,0,0),             ## w2<0
      c(0,0,0,0,0,1,0) )				 	## shape>2.001
    
    ci<-c(-0.001,-0.001,0.999,-0.001,-2.001)
    
    # numerical_gradient <- function(param,model,daily_ret,mv_m,K,distribution,lag_fun) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      mv_m=mv_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      mv_m=mv_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun)
    
    
  }else if(model=="GM_noskew_RV_X"&distribution=="norm"&lag_fun=="Beta"){
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=8)
    colnames(begin_val)<-c("alpha","beta","m","theta_RV","theta_X","w2_RV","w2_X","mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-stats::runif(R,min=-1,max=1)
    begin_val[,4]<-stats::runif(R,min=-1,max=1)
    begin_val[,5]<-stats::runif(R,min=-1,max=1)
    begin_val[,6]<-1.01
    begin_val[,7]<-1.01
    begin_val[,8]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,RV_m=RV_m_in_s,
                                                 mv_m=mv_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,0,0,0,0),       	 			## alpha>0.001
      c(0,1,0,0,0,0,0,0),        		 		## beta>0.001
      c(-1,-1,0,0,0,0,0,0),	     				## alpha+beta<1
      c(0,0,0,0,0,1,0,0),				 	## w2_RV>1.001
      c(0,0,0,0,0,0,1,0))				 	## w2_X>1.001
    
    ci<-c(-0.001,-0.001,0.999,-1.001,-1.001)
    
    # numerical_gradient <- function(param,model,daily_ret,RV_m,mv_m,K,distribution,lag_fun) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      RV_m=RV_m_in_s,
      mv_m=mv_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      RV_m=RV_m_in_s,mv_m=mv_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun)
    
  } else if(model=="GM_noskew_RV_X"&distribution=="norm"&lag_fun=="Almon"){
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=8)
    colnames(begin_val)<-c("alpha","beta","m","theta_RV","theta_X","w2_RV","w2_X","mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-stats::runif(R,min=-1,max=1)
    begin_val[,4]<-stats::runif(R,min=-1,max=1)
    begin_val[,5]<-stats::runif(R,min=-1,max=1)
    begin_val[,6]<--0.1
    begin_val[,7]<--0.1
    begin_val[,8]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,RV_m=RV_m_in_s,
                                                 mv_m=mv_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,0,0,0,0),       	 			## alpha>0.001
      c(0,1,0,0,0,0,0,0),        		 		## beta>0.001
      c(-1,-1,0,0,0,0,0,0),	     				## alpha+beta<1
      c(0,0,0,0,0,-1,0,0),				 	## w2_RV<0
      c(0,0,0,0,0,0,-1,0))				 	## w2_X<0
    
    ci<-c(-0.001,-0.001,0.999,-0.001,-0.001)
    
    # numerical_gradient <- function(param,model,daily_ret,RV_m,mv_m,K,distribution,lag_fun) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      RV_m=RV_m_in_s,
      mv_m=mv_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      RV_m=RV_m_in_s,mv_m=mv_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun)
    
  }else if(model=="GM_noskew_RV_X"&distribution=="std"&lag_fun=="Beta"){
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=9)
    colnames(begin_val)<-c("alpha","beta","m","theta_RV","theta_X","w2_RV","w2_X","v","mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-stats::runif(R,min=-1,max=1)
    begin_val[,4]<-stats::runif(R,min=-1,max=1)
    begin_val[,5]<-stats::runif(R,min=-1,max=1)
    begin_val[,6]<-1.01
    begin_val[,7]<-1.01
    begin_val[,8]<-stats::runif(R,min=2.01,max=10)
    begin_val[,9]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,RV_m=RV_m_in_s,
                                                 mv_m=mv_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,0,0,0,0,0),       	 			## alpha>0.001
      c(0,1,0,0,0,0,0,0,0),        		 		## beta>0.001
      c(-1,-1,0,0,0,0,0,0,0),	     				## alpha+beta<1
      c(0,0,0,0,0,1,0,0,0),              	## w2_RV>1.001
      c(0,0,0,0,0,0,1,0,0),              	## w2_X>1.001
      c(0,0,0,0,0,0,0,1,0))				       ## shape>2.001
    
    ci<-c(-0.001,-0.001,0.999,-1.001,-1.001,-2.001)
    
    # numerical_gradient <- function(param,model,daily_ret,RV_m,mv_m,K,distribution,lag_fun) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      RV_m=RV_m_in_s,
      mv_m=mv_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      RV_m=RV_m_in_s,mv_m=mv_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun)
    
    
  }else if(model=="GM_noskew_RV_X"&distribution=="std"&lag_fun=="Almon"){
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=9)
    colnames(begin_val)<-c("alpha","beta","m","theta_RV","theta_X","w2_RV","w2_X","v","mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-stats::runif(R,min=-1,max=1)
    begin_val[,4]<-stats::runif(R,min=-1,max=1)
    begin_val[,5]<-stats::runif(R,min=-1,max=1)
    begin_val[,6]<--0.1
    begin_val[,7]<--0.1
    begin_val[,8]<-stats::runif(R,min=2.01,max=10)
    begin_val[,9]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,RV_m=RV_m_in_s,
                                                 mv_m=mv_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,0,0,0,0,0),       	 			## alpha>0.001
      c(0,1,0,0,0,0,0,0,0),        		 		## beta>0.001
      c(-1,-1,0,0,0,0,0,0,0),	     				## alpha+beta<1
      c(0,0,0,0,0,-1,0,0,0),             ## w2<0
      c(0,0,0,0,0,0,-1,0,0),             ## w2<0
      c(0,0,0,0,0,0,0,1,0) )				 	## shape>2.001
    
    ci<-c(-0.001,-0.001,0.999,-0.001,-0.001,-2.001)
    
    # numerical_gradient <- function(param,model,daily_ret,RV_m,mv_m,K,distribution,lag_fun) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      RV_m=RV_m_in_s,
      mv_m=mv_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      RV_m=RV_m_in_s,mv_m=mv_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun)
    
    
  }else if(model=="GM_noskew_RV_sb"&distribution=="norm"&lag_fun=="Beta"){
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(stats::time(r_t_in_s))
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=2*n_n+6)
    colnames(begin_val)<-c("alpha","beta","w2",paste0("D_m_", seq_len(n_n+1)),paste0("D_theta_", seq_len(n_n+1)),"mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-1.01
    for (i in 4:(2*n_n+5)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    begin_val[,(2*n_n+6)]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,
                                                 RV_m=RV_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun,
                                                 X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,rep(0,(2*n_n+3))),       	 			## alpha>0.001
      c(0,1,0,rep(0,(2*n_n+3))),        		 		## beta>0.001
      c(-1,-1,0,rep(0,(2*n_n+3))),	     				## alpha+beta<1
      c(0,0,1,rep(0,(2*n_n+3))))				 	## w2>1.001
    
    ci<-c(-0.001,-0.001,0.999,-1.001)
    
    # numerical_gradient <- function(param,model,daily_ret,RV_m,K,distribution,lag_fun,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      RV_m=RV_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      X_stru_break_date_vec=X_stru_break_date_vec,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      RV_m=RV_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    
  } else if(model=="GM_noskew_RV_sb"&distribution=="norm"&lag_fun=="Almon"){
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(stats::time(r_t_in_s))
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=2*n_n+6)
    colnames(begin_val)<-c("alpha","beta","w2",paste0("D_m_", seq_len(n_n+1)),paste0("D_theta_", seq_len(n_n+1)),"mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<--0.01
    for (i in 4:(2*n_n+5)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    begin_val[,(2*n_n+6)]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,
                                                 RV_m=RV_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun,
                                                 X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,rep(0,(2*n_n+3))),       	 			## alpha>0.001
      c(0,1,0,rep(0,(2*n_n+3))),        		 		## beta>0.001
      c(-1,-1,0,rep(0,(2*n_n+3))),	     				## alpha+beta<1
      c(0,0,-1,rep(0,(2*n_n+3))))				 	## w2<0
    
    ci<-c(-0.001,-0.001,0.999,-0.001)
    
    # numerical_gradient <- function(param,model,daily_ret,RV_m,K,distribution,lag_fun,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      RV_m=RV_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      X_stru_break_date_vec=X_stru_break_date_vec,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      RV_m=RV_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    
  }else if(model=="GM_noskew_RV_sb"&distribution=="std"&lag_fun=="Beta"){
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(stats::time(r_t_in_s))
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=2*n_n+7)
    colnames(begin_val)<-c("alpha","beta","w2","v",paste0("D_m_", seq_len(n_n+1)),paste0("D_theta_", seq_len(n_n+1)),"mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-1.01
    begin_val[,4]<-stats::runif(R,min=2.01,max=10)
    for (i in 5:(2*n_n+6)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    begin_val[,(2*n_n+7)]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,
                                                 RV_m=RV_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,rep(0,(2*n_n+3))),       	 			## alpha>0.001
      c(0,1,0,0,rep(0,(2*n_n+3))),        		 		## beta>0.001
      c(-1,-1,0,0,rep(0,(2*n_n+3))),	     				## alpha+beta<1
      c(0,0,1,0,rep(0,(2*n_n+3))),              	## w2>1.001
      c(0,0,0,1,rep(0,(2*n_n+3))))				       ## shape>2.001
    
    ci<-c(-0.001,-0.001,0.999,-1.001,-2.001)
    
    # numerical_gradient <- function(param,model,daily_ret,RV_m,K,distribution,lag_fun,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      RV_m=RV_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      X_stru_break_date_vec=X_stru_break_date_vec,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      RV_m=RV_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    
    
  }else if(model=="GM_noskew_RV_sb"&distribution=="std"&lag_fun=="Almon"){
    
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(stats::time(r_t_in_s))
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=2*n_n+7)
    colnames(begin_val)<-c("alpha","beta","w2","v",paste0("D_m_", seq_len(n_n+1)),paste0("D_theta_", seq_len(n_n+1)),"mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<- -0.1
    begin_val[,4]<-stats::runif(R,min=2.01,max=10)
    for (i in 5:(2*n_n+6)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    begin_val[,(2*n_n+7)]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,
                                                 RV_m=RV_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,rep(0,(2*n_n+3))),       	 			## alpha>0.001
      c(0,1,0,0,rep(0,(2*n_n+3))),        		 		## beta>0.001
      c(-1,-1,0,0,rep(0,(2*n_n+3))),	     				## alpha+beta<1
      c(0,0,-1,0,rep(0,(2*n_n+3))),             ## w2<0
      c(0,0,0,1,rep(0,(2*n_n+3))) )				 	## shape>2.001
    
    ci<-c(-0.001,-0.001,0.999,-0.001,-2.001)
    
    # numerical_gradient <- function(param,model,daily_ret,RV_m,K,distribution,lag_fun,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      RV_m=RV_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      X_stru_break_date_vec=X_stru_break_date_vec,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      RV_m=RV_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    
  }else if(model=="GM_noskew_X_sb"&distribution=="norm"&lag_fun=="Beta"){
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(stats::time(r_t_in_s))
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=2*n_n+6)
    colnames(begin_val)<-c("alpha","beta","w2",paste0("D_m_", seq_len(n_n+1)),paste0("D_theta_", seq_len(n_n+1)),"mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-1.01
    for (i in 4:(2*n_n+5)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    begin_val[,(2*n_n+6)]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,
                                                 mv_m=mv_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun,
                                                 X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,rep(0,(2*n_n+3))),       	 			## alpha>0.001
      c(0,1,0,rep(0,(2*n_n+3))),        		 		## beta>0.001
      c(-1,-1,0,rep(0,(2*n_n+3))),	     				## alpha+beta<1
      c(0,0,1,rep(0,(2*n_n+3))))				 	## w2>1.001
    
    ci<-c(-0.001,-0.001,0.999,-1.001)
    
    # numerical_gradient <- function(param,model,daily_ret,mv_m,K,distribution,lag_fun,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      mv_m=mv_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      X_stru_break_date_vec=X_stru_break_date_vec,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      mv_m=mv_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    
  } else if(model=="GM_noskew_X_sb"&distribution=="norm"&lag_fun=="Almon"){
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(stats::time(r_t_in_s))
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=2*n_n+6)
    colnames(begin_val)<-c("alpha","beta","w2",paste0("D_m_", seq_len(n_n+1)),paste0("D_theta_", seq_len(n_n+1)),"mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<--0.01
    for (i in 4:(2*n_n+5)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    begin_val[,(2*n_n+6)]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,
                                                 mv_m=mv_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun,
                                                 X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,rep(0,(2*n_n+3))),       	 			## alpha>0.001
      c(0,1,0,rep(0,(2*n_n+3))),        		 		## beta>0.001
      c(-1,-1,0,rep(0,(2*n_n+3))),	     				## alpha+beta<1
      c(0,0,-1,rep(0,(2*n_n+3))))				 	## w2<0
    
    ci<-c(-0.001,-0.001,0.999,-0.001)
    
    # numerical_gradient <- function(param,model,daily_ret,mv_m,K,distribution,lag_fun,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      mv_m=mv_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      X_stru_break_date_vec=X_stru_break_date_vec,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      mv_m=mv_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    
  }else if(model=="GM_noskew_X_sb"&distribution=="std"&lag_fun=="Beta"){
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(stats::time(r_t_in_s))
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=2*n_n+7)
    colnames(begin_val)<-c("alpha","beta","w2","v",paste0("D_m_", seq_len(n_n+1)),paste0("D_theta_", seq_len(n_n+1)),"mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-1.01
    begin_val[,4]<-stats::runif(R,min=2.01,max=10)
    for (i in 5:(2*n_n+6)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    begin_val[,(2*n_n+7)]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,
                                                 mv_m=mv_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,rep(0,(2*n_n+3))),       	 			## alpha>0.001
      c(0,1,0,0,rep(0,(2*n_n+3))),        		 		## beta>0.001
      c(-1,-1,0,0,rep(0,(2*n_n+3))),	     				## alpha+beta<1
      c(0,0,1,0,rep(0,(2*n_n+3))),              	## w2>1.001
      c(0,0,0,1,rep(0,(2*n_n+3))))				       ## shape>2.001
    
    ci<-c(-0.001,-0.001,0.999,-1.001,-2.001)
    
    # numerical_gradient <- function(param,model,daily_ret,mv_m,K,distribution,lag_fun,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      mv_m=mv_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      X_stru_break_date_vec=X_stru_break_date_vec,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      mv_m=mv_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    
    
  }else if(model=="GM_noskew_X_sb"&distribution=="std"&lag_fun=="Almon"){
    
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(stats::time(r_t_in_s))
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=2*n_n+7)
    colnames(begin_val)<-c("alpha","beta","w2","v",paste0("D_m_", seq_len(n_n+1)),paste0("D_theta_", seq_len(n_n+1)),"mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<- -0.1
    begin_val[,4]<-stats::runif(R,min=2.01,max=10)
    for (i in 5:(2*n_n+6)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    begin_val[,(2*n_n+7)]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,
                                                 mv_m=mv_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,rep(0,(2*n_n+3))),       	 			## alpha>0.001
      c(0,1,0,0,rep(0,(2*n_n+3))),        		 		## beta>0.001
      c(-1,-1,0,0,rep(0,(2*n_n+3))),	     				## alpha+beta<1
      c(0,0,-1,0,rep(0,(2*n_n+3))),             ## w2<0
      c(0,0,0,1,rep(0,(2*n_n+3))) )				 	## shape>2.001
    
    ci<-c(-0.001,-0.001,0.999,-0.001,-2.001)
    
    # numerical_gradient <- function(param,model,daily_ret,mv_m,K,distribution,lag_fun,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      mv_m=mv_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      X_stru_break_date_vec=X_stru_break_date_vec,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      mv_m=mv_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    
  }else if(model=="GM_noskew_RV_X_sb"&distribution=="norm"&lag_fun=="Beta"){
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(stats::time(r_t_in_s))
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=3*n_n+8)
    colnames(begin_val)<-c("alpha","beta","w2_RV","w2_X",
                           paste0("D_m_", seq_len(n_n+1)),
                           paste0("D_theta_RV_", seq_len(n_n+1)),paste0("D_theta_X_", seq_len(n_n+1)),"mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-1.01
    begin_val[,4]<-1.01
    for (i in 5:(3*n_n+7)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    begin_val[,(3*n_n+8)]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,RV_m=RV_m_in_s,
                                                 mv_m=mv_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,rep(0,(3*n_n+4))),       	 			## alpha>0.001
      c(0,1,0,0,rep(0,(3*n_n+4))),        		 		## beta>0.001
      c(-1,-1,0,0,rep(0,(3*n_n+4))),	     				## alpha+beta<1
      c(0,0,1,0,rep(0,(3*n_n+4))),				 	## w2_RV>1.001
      c(0,0,0,1,rep(0,(3*n_n+4))))				 	## w2_X>1.001
    
    ci<-c(-0.001,-0.001,0.999,-1.001,-1.001)
    
    # numerical_gradient <- function(param,model,daily_ret,RV_m,mv_m,K,distribution,lag_fun,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      RV_m=RV_m_in_s,
      mv_m=mv_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      X_stru_break_date_vec=X_stru_break_date_vec,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      RV_m=RV_m_in_s,mv_m=mv_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    
  } else if(model=="GM_noskew_RV_X_sb"&distribution=="norm"&lag_fun=="Almon"){
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(stats::time(r_t_in_s))
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=3*n_n+8)
    colnames(begin_val)<-c("alpha","beta","w2_RV","w2_X",
                           paste0("D_m_", seq_len(n_n+1)),
                           paste0("D_theta_RV_", seq_len(n_n+1)),paste0("D_theta_X_", seq_len(n_n+1)),"mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<--0.1
    begin_val[,4]<--0.1
    for (i in 5:(3*n_n+7)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    begin_val[,(3*n_n+8)]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,RV_m=RV_m_in_s,
                                                 mv_m=mv_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,rep(0,(3*n_n+4))),       	 			## alpha>0.001
      c(0,1,0,0,rep(0,(3*n_n+4))),        		 		## beta>0.001
      c(-1,-1,0,0,rep(0,(3*n_n+4))),	     				## alpha+beta<1
      c(0,0,-1,0,rep(0,(3*n_n+4))),				 	## w2_RV<0
      c(0,0,0,-1,rep(0,(3*n_n+4))))				 	## w2_X<0
    
    
    ci<-c(-0.001,-0.001,0.999,-0.001,-0.001)
    
    # numerical_gradient <- function(param,model,daily_ret,RV_m,mv_m,K,distribution,lag_fun,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      RV_m=RV_m_in_s,
      mv_m=mv_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      X_stru_break_date_vec=X_stru_break_date_vec,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      RV_m=RV_m_in_s,mv_m=mv_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    
  }else if(model=="GM_noskew_RV_X_sb"&distribution=="std"&lag_fun=="Beta"){
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(stats::time(r_t_in_s))
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=3*n_n+9)
    colnames(begin_val)<-c("alpha","beta","w2_RV","w2_X","v",
                           paste0("D_m_", seq_len(n_n+1)),
                           paste0("D_theta_RV_", seq_len(n_n+1)),paste0("D_theta_X_", seq_len(n_n+1)),"mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-1.01
    begin_val[,4]<-1.01
    begin_val[,5]<-stats::runif(R,min=2.01,max=10)
    for (i in 6:(3*n_n+8)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    begin_val[,(3*n_n+9)]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,RV_m=RV_m_in_s,
                                                 mv_m=mv_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,0,rep(0,(3*n_n+4))),       	 			## alpha>0.001
      c(0,1,0,0,0,rep(0,(3*n_n+4))),        		 		## beta>0.001
      c(-1,-1,0,0,0,rep(0,(3*n_n+4))),	     				## alpha+beta<1
      c(0,0,1,0,0,rep(0,(3*n_n+4))),				 	## w2_RV>1.001
      c(0,0,0,1,0,rep(0,(3*n_n+4))),          ## w2_X>1.001
      c(0,0,0,0,1,rep(0,(3*n_n+4))))				 	## shape>2.001
    
    
    ci<-c(-0.001,-0.001,0.999,-1.001,-1.001,-2.001)
    
    # numerical_gradient <- function(param,model,daily_ret,RV_m,mv_m,K,distribution,lag_fun,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      RV_m=RV_m_in_s,
      mv_m=mv_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      X_stru_break_date_vec=X_stru_break_date_vec,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      RV_m=RV_m_in_s,mv_m=mv_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    
    
  }else if(model=="GM_noskew_RV_X_sb"&distribution=="std"&lag_fun=="Almon"){
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(stats::time(r_t_in_s))
    
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    begin_val<-matrix(NA,nrow=R,ncol=3*n_n+9)
    colnames(begin_val)<-c("alpha","beta","w2_RV","w2_X","v",
                           paste0("D_m_", seq_len(n_n+1)),
                           paste0("D_theta_RV_", seq_len(n_n+1)),paste0("D_theta_X_", seq_len(n_n+1)),"mu")
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<--0.1
    begin_val[,4]<--0.1
    begin_val[,5]<-stats::runif(R,min=2.01,max=10)
    for (i in 6:(3*n_n+8)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    begin_val[,(3*n_n+9)]<-mean(r_t_in_s, na.rm = TRUE)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(garchmidas_noskew_loglik(begin_val[i,],model=model,
                                                 daily_ret=r_t_in_s,RV_m=RV_m_in_s,
                                                 mv_m=mv_m_in_s,K=K,
                                                 distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    
    ui<-rbind(
      c(1,0,0,0,0,rep(0,(3*n_n+4))),       	 			## alpha>0.001
      c(0,1,0,0,0,rep(0,(3*n_n+4))),        		 		## beta>0.001
      c(-1,-1,0,0,0,rep(0,(3*n_n+4))),	     				## alpha+beta<1
      c(0,0,-1,0,0,rep(0,(3*n_n+4))),				 	## w2_RV<0
      c(0,0,0,-1,0,rep(0,(3*n_n+4))),          ## w2_X<0
      c(0,0,0,0,1,rep(0,(3*n_n+4))))				 	## shape>2.001
    
    
    ci<-c(-0.001,-0.001,0.999,-0.001,-0.001,-2.001)
    
    # numerical_gradient <- function(param,model,daily_ret,RV_m,mv_m,K,distribution,lag_fun,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(daily_ret)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- garchmidas_noskew_loglik(param = param + step, model=model,daily_ret = daily_ret,
    #                                          RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- garchmidas_noskew_loglik(param = param - step, model=model,daily_ret = daily_ret,
    #                                           RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    est<-suppressWarnings(maxLik(
      logLik=garchmidas_noskew_loglik,
      # grad = numerical_gradient,
      start=start_val,
      
      model=model,
      daily_ret=r_t_in_s,
      RV_m=RV_m_in_s,
      mv_m=mv_m_in_s,
      K=K,
      distribution=distribution,
      lag_fun=lag_fun,
      X_stru_break_date_vec=X_stru_break_date_vec,
      constraints=list(ineqA=ui, ineqB=ci),
      iterlim=1000,
      method="BFGS"))
    
    est_coef <- stats::coef(est)
    vol_est <-  garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_in_s,
                                      RV_m=RV_m_in_s,mv_m=mv_m_in_s,K=K,distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec)
    
  }
  
  N_coef<-length(est_coef)
  
  mat_coef<-data.frame(rep(NA,N_coef),rep(NA,N_coef),rep(NA,N_coef),rep(NA,N_coef))
  colnames(mat_coef)<-c("Estimate","Std. Error","t value","Pr(>|t|)")
  
  rownames(mat_coef)<-names(stats::coef(est))
  
  mat_coef[,1]<-round(est_coef,6)
  mat_coef[,2]<-round(QMLE_sd(est),6)
  mat_coef[,3]<-round(est_coef/QMLE_sd(est),6)
  mat_coef[,4]<-round(apply(rbind(est_coef/QMLE_sd(est)),1,function(x) 2*(1-stats::pnorm(abs(x)))),6)
  
  if (is.null(out_of_sample)){
    
    N_est<-length(vol_est[["h_it_2"]])
    
    vol_est_lf<-zoo::coredata(vol_est[["h_it_2"]])
    vol_proxy_lf<-zoo::coredata(vol_proxy_in_s)
    
    res<-list(
      model=model,
      rob_coef_mat=mat_coef,
      obs=N,
      period=range(stats::time(r_t_in_s)),
      loglik=as.numeric(stats::logLik(est)),
      inf_criteria=Inf_criteria(est),
      loss_in_s=LF_f(vol_est_lf,vol_proxy_lf),
      est_vol_in_s=vol_est_lf^0.5,
      est_lr_in_s=zoo::coredata(vol_est[["tau_d"]]),
      est_sr_in_s=zoo::coredata(vol_est[["g_it"]]))
    
  } else {
    
    r_t_oos<-daily_ret[(N-out_of_sample+1):N]
    
    if(!is.null(RV_m)){
      RV_m_oos<-RV_m[,(N-out_of_sample+1):(N)]
    } else {
      
      RV_m_oos <- NULL
    }
    
    if(!is.null(mv_m)){
      mv_m_oos<-mv_m[,(N-out_of_sample+1):(N)]
    } else {
      
      mv_m_oos <- NULL
    }
    
    if(is.null(vol_proxy)){
      vol_proxy_oos<-r_t_oos^2
    } else {
      vol_proxy_oos<-vol_proxy[(N-out_of_sample+1):N]
    }
    
    N_est<-length(length(vol_est[["h_it_2"]]))
    
    vol_est_lf<-zoo::coredata(vol_est[["h_it_2"]])
    vol_proxy_lf<-zoo::coredata(vol_proxy_in_s)
    
    vol_proxy_oos_lf<-zoo::coredata(vol_proxy_oos)
    
    
    vol_est_oos <- garchmidas_noskew_est(est_coef,model=model,daily_ret=r_t_oos,
                                         RV_m=RV_m_oos,mv_m=mv_m_oos,
                                         K=K,distribution=distribution,lag_fun=lag_fun,
                                         X_stru_break_date_vec=X_stru_break_date_vec)
    
    
    vol_est_oos_lf<-zoo::coredata(vol_est_oos[["h_it_2"]])
    
    res<-list(
      model=model,
      rob_coef_mat=mat_coef,
      obs=length(r_t_in_s),
      period=range(time(r_t_in_s)),
      loglik=as.numeric(stats::logLik(est)),
      inf_criteria=Inf_criteria(est),
      loss_in_s=LF_f(vol_est_lf,vol_proxy_lf),
      est_vol_in_s=vol_est_lf^0.5,
      est_lr_in_s=zoo::coredata(vol_est[["tau_d"]]),
      est_sr_in_s=zoo::coredata(vol_est[["g_it"]]),
      loss_oos=LF_f(vol_est_oos_lf,vol_proxy_oos_lf),
      est_vol_oos=vol_est_oos_lf^0.5,
      est_lr_oos=zoo::coredata(vol_est_oos[["tau_d"]]),
      est_sr_oos=zoo::coredata(vol_est_oos[["g_it"]]))
    
  }
  
  class(res)<-c("gmnoskew")
  return(res)
  
  
}

summary.gmnoskew <- function(object, ...) {
  
  model<-object$model
  mat_coef<-object$rob_coef_mat
  Obs<-object$obs
  Period<-object$period
  loss<-object$loss_in_s
  
  
  Period<-paste(substr(Period[1], 1, 10),"/",
                substr(Period[2], 1, 10),sep="")
  
  p_value<-mat_coef[,4]
  
  sig<-ifelse(p_value<=0.01,"***",ifelse(p_value>0.01&p_value<=0.05,"**",
                                         ifelse(p_value>0.05&p_value<=0.1,"*"," ")))
  
  mat_coef<-round(mat_coef,4)
  
  mat_coef<-cbind(mat_coef,Sig.=sig)
  
  cat(
    cat("\n"),
    cat("Coefficients:\n"),
    cat(utils::capture.output(mat_coef),  sep = '\n'),
    cat("--- \n"),
    cat("Signif. codes: 0.01 '***', 0.05 '**', 0.1 '*' \n"),
    cat("\n"),
    cat("Obs.:", paste(Obs, ".",sep=""), "Sample Period:", Period, "\n"),
    cat("MSE(%):", paste(as.numeric(round(loss[1],6)), "; ",sep=""), 
        "QLIKE:", as.numeric(round(loss[2],6)), "\n"),
    cat("\n"))
  
}

print.gmnoskew <- function(x, ...) {
  
  options(scipen = 999)
  
  model<-x[[1]]
  mat_coef<-x[[2]]
  
  coef_char<-as.character(round(mat_coef[,1],4))
  
  row_names<-gsub("\\s", " ", format(rownames(mat_coef), width=9))
  coef_val<-gsub("\\s", " ", format(coef_char,width=9))
  
  cat(
    cat("\n"),
    cat(paste("Model:",model),"\n"),
    cat("\n"),
    cat(paste("Coefficients: \n",sep="\n")),
    cat(row_names, sep=" ", "\n"),
    cat(coef_val, sep=" ", "\n"),
    cat("\n"))
}

logLik.gmnoskew <- function(object, ...) {
  # 1. 提取存储在对象中的 loglik 值
  val <- object$loglik
  
  # 2. 必须设置自由度 (df)，否则 AIC/BIC 无法计算
  # 这里假设 rob_coef_mat 的行数等于估计参数的个数
  attr(val, "df") <- nrow(object$rob_coef_mat)
  
  # 3. 设置观测值数量 (nobs)，这对 BIC 计算很重要
  attr(val, "nobs") <- object$obs
  
  # 4. 指定返回值的类为 "logLik"
  class(val) <- "logLik"
  
  return(val)
}

coef.gmnoskew <- function(object, ...) {
  # 同时也建议定义 coef 方法，方便提取系数
  return(object$rob_coef_mat[, "Estimate"])
}

integrated_garchmidas_noskew_output_function <- function(model,daily_ret,macro=NULL,type,K,distribution="norm",lag_fun="Beta",
                                                         X_stru_break_date_vec=NULL,out_of_sample=NULL,vol_proxy=NULL,R=NULL){
  ####available for monthly/quarterly
  cond_r_t<- class(daily_ret)[1]
  
  K_c <- 0
  
  if(cond_r_t != "xts") { stop(
    cat("#Warning:\n Parameter 'daily_ret' must be an xts object. Please provide it in the correct form \n")
  )}
  
  if(!is.null(macro)&&class(macro)[1] != "xts") { stop(
    cat("#Warning:\n Parameter 'macro' must be an xts object. Please provide it in the correct form \n")
  )}
  
  if(!is.null(vol_proxy)&&class(vol_proxy)[1] != "xts") { stop(
    cat("#Warning:\n Parameter 'vol_proxy' must be an xts object. Please provide it in the correct form \n")
  )}
  
  RV_macro <- RV_gen(x=daily_ret,type=type)
  
  if(type=="monthly"){
    
    RV_required_r_t <- cut_func(r_t=daily_ret,K=K,K_c=K_c,macro=RV_macro,"monthly")
    
  }else if(type=="quarterly"){
    
    RV_required_r_t <- cut_func(r_t=daily_ret,K=K,K_c=K_c,macro=RV_macro,"quarterly")
    
  }
  
  ##GM_noskew_RV/GM_noskew_X/GM_noskew_RV_X/GM_noskew_RV_sb/GM_noskew_X_sb/GM_noskew_RV_X_sb  
  if(univ_model %in% c("GM_noskew_X", "GM_noskew_RV_X", "GM_noskew_X_sb", "GM_noskew_RV_X_sb")){
    
    if(type=="monthly"){
      
      mv_required_r_t <- cut_func(r_t=daily_ret,K=K,macro=macro,type="monthly")
      
    }else if(type=="quarterly"){
      
      mv_required_r_t <- cut_func(r_t=daily_ret,K=K,macro=macro,type="quarterly")
      
    }
    
    if (model=="GM_noskew_X"|model=="GM_noskew_X_sb"){
      
      
      if(!is.null(vol_proxy)){
        
        ref_dates <- index(mv_required_r_t)
        
        if (!all(ref_dates %in% index(vol_proxy))) {
          
          stop(paste0("vol_proxy data do not cover daily_ret"))
          
        } else {
          
          vol_proxy <- vol_proxy[ref_dates]
        }        
      }
      mv_m<- mv_into_mat2(x=mv_required_r_t,mv=macro,K=K,type=type)
      
      result <- garchmidas_noskew_fit(model=model,daily_ret=mv_required_r_t,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,
                                      X_stru_break_date_vec=X_stru_break_date_vec,
                                      out_of_sample=out_of_sample,vol_proxy=vol_proxy,R=R)
    } else if(model=="GM_noskew_RV_X"|model=="GM_noskew_RV_X_sb"){
      
      common_dates <- as.Date(intersect(index(RV_required_r_t), index(mv_required_r_t)))
      required_r_t <- mv_required_r_t[common_dates]
      
      if(!is.null(vol_proxy)){
        
        ref_dates <- common_dates
        
        if (!all(ref_dates %in% index(vol_proxy))) {
          
          stop(paste0("vol_proxy data do not cover daily_ret"))
          
        } else {
          
          vol_proxy <- vol_proxy[ref_dates]
        }        
      }
      
      mv_m<- mv_into_mat2(x=required_r_t,mv=macro,K=K,type=type)
      RV_m<- mv_into_mat2(x=required_r_t,mv=RV_macro,K=K,type=type)
      
      result <- garchmidas_noskew_fit(model=model,daily_ret= required_r_t,RV_m=RV_m,mv_m=mv_m,K=K,distribution=distribution,lag_fun=lag_fun,
                                      X_stru_break_date_vec=X_stru_break_date_vec,
                                      out_of_sample=out_of_sample,vol_proxy=vol_proxy,R=R)
    }
  } else {
    
    RV_m<- mv_into_mat2(x=RV_required_r_t,mv=RV_macro,K=K,type=type)
    
    if(!is.null(vol_proxy)){
      
      ref_dates <- index(RV_required_r_t)
      
      if (!all(ref_dates %in% index(vol_proxy))) {
        
        stop(paste0("vol_proxy data do not cover daily_ret"))
        
      } else {
        
        vol_proxy <- vol_proxy[ref_dates]
      }        
    }
    
    result <- garchmidas_noskew_fit(model=model,daily_ret=RV_required_r_t,RV_m=RV_m,K=K,distribution=distribution,lag_fun=lag_fun,
                                    X_stru_break_date_vec=X_stru_break_date_vec,
                                    out_of_sample=out_of_sample,vol_proxy=vol_proxy,R=R)
    
  }
  
  
  return(result)
  
  
}



##DCCMIDAS_RC/DCCMIDAS_RC_sb/DCCMIDAS_X/DCCMIDAS_X_sb/DCCMIDAS_RC_X/DCCMIDAS_RC_X_sb

Inf_criteria_2step <- function(res_garch, res_dcc) {
  # --- 1. 计算总参数数量 (k_total) ---
  # 使用 length(coef(x)) 是通用的，但要确保显式调用
  k_garch <- sum(sapply(res_garch, function(x) {
    # 兼容 rugarch S4 对象和普通的类
    if (isS4(x)) {
      return(length(rugarch::coef(x)))
    } else {
      return(length(stats::coef(x)))
    }
  }))
  
  # res_dcc 通常是 maxLik 对象
  k_dcc  <- length(stats::coef(res_dcc))
  k_total <- k_garch + k_dcc
  
  # --- 2. 计算总对数似然值 (LL_total) ---
  loglik_garch <- sum(sapply(res_garch, function(x) {
    if (isS4(x)) {
      return(as.numeric(rugarch::likelihood(x)))
    } else {
      # 针对你自定义的 gmnoskew 类，如果 stats::logLik 报错，直接取 $loglik
      ll <- try(stats::logLik(x), silent = TRUE)
      if (inherits(ll, "try-error")) return(as.numeric(x$loglik))
      return(as.numeric(ll))
    }
  }))
  
  loglik_dcc   <- as.numeric(stats::logLik(res_dcc))
  loglik_total <- loglik_garch + loglik_dcc
  
  # --- 3. 确定样本量 (n) ---
  # 注意：maxLik 对象的 gradientObs 只有在收敛后才可靠
  n <- nrow(res_dcc$gradientObs) 
  
  # --- 4. 计算准则值 ---
  aic <- 2 * k_total - 2 * loglik_total
  bic <- k_total * log(n) - 2 * loglik_total
  
  inf <- round(c(AIC = aic, BIC = bic, k = k_total, LogLik = loglik_total), 6)
  return(inf)
}

safe_function_with_warning <- function(expr) {
  tryCatch(
    {
      result <- eval(expr)
      return(result)
    },
    
    error = function(e) {
      return(NA)
    }
  )
}

#dccmidas_noskew_loglik
dccmidas_noskew_loglik <- function(param,model,res,days_period,mv_c=NULL,lag_fun="Beta",N_c=NULL,K_c,X_stru_break_date_vec=NULL){
  
  
  
  if (model=="DCCMIDAS_RC"){
    
    a<-param[1]
    b<-param[2]
    w1<- ifelse(lag_fun=="Beta",1,0)
    w2<-param[3]
    
    ##### index
    
    Num_assets<-ncol(res)		# number of assets
    TT<-length(as.Date(days_period))			# number of daily observations
    
    ##### matrices and vectors
    
    
    
    Q_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    R_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    log_det_R_t<-rep(0,TT)
    Eps_t_R_t_Eps_t<-rep(0,TT)
    
    ##### first cycle
    
    C_t <- midas_rc_rcpp(N_c, Num_assets, TT, res)
    
    #  C_t<-array(0,dim=c(Num_assets,Num_assets,TT))
    # for(tt in (N_c+1):TT){
    #   
    #   window_res <- res[(tt-N_c):tt, ] 
    #   V_t_0.5 <- sqrt(colSums(window_res^2)) 
    #   Prod_eps_t<-crossprod(res[tt:(tt - N_c), ])
    #   C_t[,,tt] <- Prod_eps_t / tcrossprod(V_t_0.5)
    #   
    # }
    
    #### R_t_bar
    
    
    weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
    
    rc_betas<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,w2))[2:(K_c+1)],0)
    R_t_bar<-array(diag(Num_assets),dim=c(Num_assets,Num_assets,TT))
    
    matrix_id<-matrix(1:Num_assets^2,ncol=Num_assets)
    matrix_id_2<-which(matrix_id==1:Num_assets^2, arr.ind=TRUE)
    
    for(i in 1:nrow(matrix_id_2)){
      R_t_bar[matrix_id_2[i,1],matrix_id_2[i,2],]  <- suppressWarnings(
        roll_sum_func(C_t[matrix_id_2[i,1],matrix_id_2[i,2],], c(K_c+1),weights = rc_betas)
      ) 
    }
    
    ll <-  as.numeric(midas_ll_rcpp(TT, Num_assets, K_c, a, b, R_t_bar, res))  
    
    
    # for (i in 1:TT) {
    #   diag(R_t_bar[,,i]) <- 1
    # }
    # 
    # 
    # ################# likelihood
    # 
    # ll<-rep(0,TT)
    # 
    # for(tt in (K_c+1):TT){
    #   Q_t[,,tt]<-(1-a-b)*R_t_bar[,,tt] + a*tcrossprod(res[tt-1, ]) + b*Q_t[,,tt-1]
    #   
    #   q_sd <- sqrt(diag(Q_t[,,tt]))
    #   R_t[,,tt] <- Q_t[,,tt] / tcrossprod(q_sd)
    #   
    #   chol_R <- tryCatch(chol(R_t[,,tt]), error = function(e) return(NULL)) # 进行 Cholesky 分解 (得到上三角矩阵)
    #   # if (is.null(chol_R)) {
    #   #   return(rep(-1e10, TT)) # 如果非正定，返回一个惩罚值代替似然值
    #   # }
    #   
    #   log_det_R_t[tt] <- 2 * sum(log(diag(chol_R))) # 等价于 log_det_R_t[tt]<-log(Det(R_t[,,tt]))
    #   
    #   tmp <- backsolve(chol_R, res[tt, ], transpose = TRUE)
    #   Eps_t_R_t_Eps_t[tt] <- sum(tmp^2)
    #   
    #   
    # }
    # 
    # ll<- - 0.5*(log_det_R_t+Eps_t_R_t_Eps_t)
    
    #sum(ll)
    
    return(ll)
    
  }else if (model=="DCCMIDAS_X"){
    
    
    
    ##### index and param
    Num_assets<-ncol(res)		# number of assets
    TT<-length(as.Date(days_period))			# number of daily observations
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    
    # 预先引用 roll_sum 避免多次查找 定义在fit函数开头里
    #roll_sum_func <- roll::roll_sum
    #weight_func_beta <- rumidas::beta_function
    #weight_func_exp <- rumidas::exp_almon
    
    a<-param[1]
    b<-param[2]
    w1<- ifelse(lag_fun=="Beta",1,0)
    
    param_delta_c <- param[3:(length_dccmidas_param+2)]
    param_theta_c <- param[(length_dccmidas_param+3):(length_dccmidas_param*2+2)] 
    param_w2_c <- param[(length_dccmidas_param*2+3):(length_dccmidas_param*3+2)]
    
    
    midasdcc_param <- list(delta_c=param_delta_c,
                           theta_c =param_theta_c, 
                           w2_c = param_w2_c)
    midasdcc_param_matrix <- list()
    
    ##### matrices and vectors
    for (pname in names(midasdcc_param)) {
      
      off_diag_elements <- midasdcc_param[[pname]]  
      cov_matrix <- matrix(0, nrow = Num_assets, ncol = Num_assets)
      cov_matrix[lower.tri(cov_matrix)] <- off_diag_elements
      cov_matrix[upper.tri(cov_matrix)] <- t(cov_matrix)[upper.tri(cov_matrix)]
      midasdcc_param_matrix[[pname]] <- cov_matrix
      
    }
    
    delta_c_matrix_tt <- array(midasdcc_param_matrix[["delta_c"]], dim = c(Num_assets, Num_assets, TT))
    theta_c_matrix_tt <- array(midasdcc_param_matrix[["theta_c"]], dim = c(Num_assets, Num_assets, TT))
    betamidas_c_matrix_tt <- array(0,dim=c(Num_assets,Num_assets,TT))
    
    Q_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    R_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    log_det_R_t<-rep(0,TT)
    Eps_t_R_t_Eps_t<-rep(0,TT)
    betas <- array(0,dim=c(Num_assets,Num_assets,(K_c+1)))
    
    
    #### R_t_bar
    weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
    
    matrix_id<-matrix(1:Num_assets^2,ncol=Num_assets)
    matrix_id_2<-which(matrix_id==1:Num_assets^2, arr.ind=TRUE)
    
    for(i in 1:nrow(matrix_id_2)){
      
      betas[matrix_id_2[i,1],matrix_id_2[i,2],]<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,midasdcc_param_matrix[["w2_c"]][matrix_id_2[i,1],matrix_id_2[i,2]]))[2:(K_c+1)],0)
      midas_c <- suppressWarnings(roll_sum_func(mv_c, c(K_c+1),weights = betas[matrix_id_2[i,1],matrix_id_2[i,2],]))
      betamidas_c_matrix_tt[matrix_id_2[i,1],matrix_id_2[i,2],] <-midas_c[(K_c+1),]
      
    }
    
    z_c <- delta_c_matrix_tt+theta_c_matrix_tt*betamidas_c_matrix_tt
    R_t_bar<-tanh(z_c)
    
    ll <-  as.numeric(midas_ll_rcpp(TT, Num_assets, K_c, a, b, R_t_bar, res))  
    
    
    # for (i in 1:TT) {
    #   diag(R_t_bar[,,i]) <- 1
    # }
    # 
    # ################# likelihood
    # 
    # ll<-rep(0,TT)
    # 
    # for(tt in (K_c+1):TT){
    #   
    #   Q_t[,,tt]<-(1-a-b)*R_t_bar[,,tt] + a*tcrossprod(res[tt-1, ]) + b*Q_t[,,tt-1]
    #   q_sd <- sqrt(diag(Q_t[,,tt]))
    #   R_t[,,tt] <- Q_t[,,tt] / tcrossprod(q_sd)
    #   
    #   chol_R <- tryCatch(chol(R_t[,,tt]), error = function(e) return(NULL)) # 进行 Cholesky 分解 (得到上三角矩阵)
    #   # if (is.null(chol_R)) {
    #   #   return(rep(-1e10, TT)) # 如果非正定，返回一个惩罚值代替似然值
    #   # }
    #   log_det_R_t[tt] <- 2 * sum(log(diag(chol_R))) # 等价于 log_det_R_t[tt]<-log(Det(R_t[,,tt]))
    #   
    #   tmp <- backsolve(chol_R, res[tt, ], transpose = TRUE)
    #   Eps_t_R_t_Eps_t[tt] <- sum(tmp^2)
    #   
    # }
    # 
    # ll<- - 0.5*(log_det_R_t+Eps_t_R_t_Eps_t)
    
    return(ll)
    
  }else if (model=="DCCMIDAS_RC_X"){
    
    
    ##### index and param
    Num_assets<-ncol(res)		# number of assets
    TT<-length(as.Date(days_period))			# number of daily observations
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    
    a<-param[1]
    b<-param[2]
    w1<- ifelse(lag_fun=="Beta",1,0)
    w2<-param[3]
    
    # 预先引用 roll_sum 避免多次查找 定义在fit函数开头里
    #roll_sum_func <- roll::roll_sum
    #weight_func_beta <- rumidas::beta_function
    #weight_func_exp <- rumidas::exp_almon
    
    ##### matrices and vectors
    ######. X
    param_delta_c <- param[4:(length_dccmidas_param+3)]
    param_theta_c <- param[(length_dccmidas_param+4):(length_dccmidas_param*2+3)] 
    param_w2_c <- param[(length_dccmidas_param*2+4):(length_dccmidas_param*3+3)]
    
    param_theta_RC_c <- param[(length_dccmidas_param*3+4):(length_dccmidas_param*4+3)] 
    
    midasdcc_param <- list(delta_c=param_delta_c,
                           theta_c =param_theta_c, 
                           w2_c = param_w2_c,
                           theta_RC_c = param_theta_RC_c)
    midasdcc_param_matrix <- list()
    
    
    for (pname in names(midasdcc_param)) {
      
      off_diag_elements <- midasdcc_param[[pname]]  
      cov_matrix <- matrix(0, nrow = Num_assets, ncol = Num_assets)
      cov_matrix[lower.tri(cov_matrix)] <- off_diag_elements
      cov_matrix[upper.tri(cov_matrix)] <- t(cov_matrix)[upper.tri(cov_matrix)]
      midasdcc_param_matrix[[pname]] <- cov_matrix
      
    }
    
    delta_c_matrix_tt <- array(midasdcc_param_matrix[["delta_c"]], dim = c(Num_assets, Num_assets, TT))
    theta_c_matrix_tt <- array(midasdcc_param_matrix[["theta_c"]], dim = c(Num_assets, Num_assets, TT))
    betamidas_c_matrix_tt <- array(0,dim=c(Num_assets,Num_assets,TT))
    
    theta_RC_c_matrix_tt <- array(midasdcc_param_matrix[["theta_RC_c"]], dim = c(Num_assets, Num_assets, TT))
    
    
    betas <- array(0,dim=c(Num_assets,Num_assets,(K_c+1)))
    
    ##################### RC
    
    RC_tt<-array(diag(Num_assets),dim=c(Num_assets,Num_assets,TT))
    
    #####################
    
    Q_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    R_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    log_det_R_t<-rep(0,TT)
    Eps_t_R_t_Eps_t<-rep(0,TT)
    
    ##### first cycle
    
    C_t <- midas_rc_rcpp(N_c, Num_assets, TT, res)
    # C_t<-array(0,dim=c(Num_assets,Num_assets,TT))
    # for(tt in (N_c+1):TT){
    #   
    #   window_res <- res[(tt-N_c):tt, ] 
    #   V_t_0.5 <- sqrt(colSums(window_res^2)) 
    #   Prod_eps_t<-crossprod(res[tt:(tt - N_c), ])
    #   C_t[,,tt] <- Prod_eps_t / tcrossprod(V_t_0.5)
    #   
    # }
    ################
    
    weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
    
    #### R_t_bar
    rc_betas<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,w2))[2:(K_c+1)],0)
    matrix_id<-matrix(1:Num_assets^2,ncol=Num_assets)
    matrix_id_2<-which(matrix_id==1:Num_assets^2, arr.ind=TRUE)
    
    for(i in 1:nrow(matrix_id_2)){
      RC_tt[matrix_id_2[i,1],matrix_id_2[i,2],]  <- suppressWarnings(
        roll_sum_func(C_t[matrix_id_2[i,1],matrix_id_2[i,2],], c(K_c+1),weights = rc_betas)) 
      
      betas[matrix_id_2[i,1],matrix_id_2[i,2],]<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,midasdcc_param_matrix[["w2_c"]][matrix_id_2[i,1],matrix_id_2[i,2]]))[2:(K_c+1)],0)
      midas_c <- suppressWarnings(roll_sum_func(mv_c, c(K_c+1),weights = betas[matrix_id_2[i,1],matrix_id_2[i,2],]))
      betamidas_c_matrix_tt[matrix_id_2[i,1],matrix_id_2[i,2],] <-midas_c[(K_c+1),]
      
    }
    
    
    
    z_c <- delta_c_matrix_tt+theta_c_matrix_tt*betamidas_c_matrix_tt+theta_RC_c_matrix_tt*RC_tt
    
    R_t_bar<-tanh(z_c)
    
    ll <-  as.numeric(midas_ll_rcpp(TT, Num_assets, K_c, a, b, R_t_bar, res))  
    
    # for (i in 1:TT) {
    #   diag(R_t_bar[,,i]) <- 1
    # }
    # 
    # ################# likelihood############
    # ll<-rep(0,TT)
    # 
    # for(tt in (K_c+1):TT){
    #   Q_t[,,tt]<-(1-a-b)*R_t_bar[,,tt] + a*tcrossprod(res[tt-1, ]) + b*Q_t[,,tt-1]
    #   
    #   q_sd <- sqrt(diag(Q_t[,,tt]))
    #   R_t[,,tt] <- Q_t[,,tt] / tcrossprod(q_sd)
    #   
    #   chol_R <- tryCatch(chol(R_t[,,tt]), error = function(e) return(NULL)) # 进行 Cholesky 分解 (得到上三角矩阵)
    #   # if (is.null(chol_R)) {
    #   #   return(rep(-1e10, TT)) # 如果非正定，返回一个惩罚值代替似然值
    #   # }
    #   
    #   log_det_R_t[tt] <- 2 * sum(log(diag(chol_R))) # 等价于 log_det_R_t[tt]<-log(Det(R_t[,,tt]))
    #   
    #   tmp <- backsolve(chol_R, res[tt, ], transpose = TRUE)
    #   Eps_t_R_t_Eps_t[tt] <- sum(tmp^2)
    #   
    #   
    # }
    # 
    # ll<- - 0.5*(log_det_R_t+Eps_t_R_t_Eps_t)
    
    return(ll)
    
  }else if (model=="DCCMIDAS_RC_sb"){
    
    a<-param[1]
    b<-param[2]
    w1<- ifelse(lag_fun=="Beta",1,0)
    w2<-param[3]
    
    
    Num_assets<-ncol(res)		# number of assets
    TT<-length(as.Date(days_period))			# number of daily observations
    
    ##### matrices and vectors
    
    
    
    Q_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    R_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    log_det_R_t<-rep(0,TT)
    Eps_t_R_t_Eps_t<-rep(0,TT)
    
    ########structural breaks
    # 提取当前资产对的系数向量 (长度 1+n_n)
    
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(days_period)
    length_dccmidas_param <- 0.5 * Num_assets * (Num_assets - 1)
    curr_m_vec <- param[4:(n_n+4)]  
    curr_theta_vec <- param[(n_n+5):(2*n_n+5)] 
    
    D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
    D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
    
    for (i in 1:n_n) {
      # 找到处于第 i 个节点之后的日期，设为 1
      D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
    }
    
    # --- 步骤 2: 准备参数矩阵 (从优化的向量中提取) ---
    
    # --- 步骤 3: 矩阵乘法计算并生成张量 ---
    # 初始化目标张量 [N, N, T]
    m_tensor <- array(0, dim = c(Num_assets, Num_assets, TT))
    theta_tensor <- array(0, dim = c(Num_assets, Num_assets, TT))
    pairs_idx <- which(lower.tri(matrix(0, Num_assets, Num_assets)), arr.ind = TRUE)
    
    for (k in 1:length_dccmidas_param) {
      i_idx <- pairs_idx[k, 1]
      j_idx <- pairs_idx[k, 2]
      
      # 核心计算：利用矩阵乘法完成所有时间点的分配
      # 第一部分：截距项时间序列 (Base + Increments)
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      
      # 第二部分：斜率项时间序列 (Base + Increments)
      theta_t_series <- as.vector(curr_theta_vec %*% D_mat)
      
      # 填充入张量的对称位置
      m_tensor[i_idx, j_idx, ] <- m_t_series
      m_tensor[j_idx, i_idx, ] <- m_t_series
      theta_tensor[i_idx, j_idx, ] <- theta_t_series
      theta_tensor[j_idx, i_idx, ] <- theta_t_series
    }
    
    ##### first cycle
    
    C_t <- midas_rc_rcpp(N_c, Num_assets, TT, res)
    
    # C_t<-array(0,dim=c(Num_assets,Num_assets,TT))
    # for(tt in (N_c+1):TT){
    #   
    #   window_res <- res[(tt-N_c):tt, ] 
    #   V_t_0.5 <- sqrt(colSums(window_res^2)) 
    #   Prod_eps_t<-crossprod(res[tt:(tt - N_c), ])
    #   C_t[,,tt] <- Prod_eps_t / tcrossprod(V_t_0.5)
    #   
    # }
    
    #### R_t_bar
    
    
    weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
    
    rc_betas<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,w2))[2:(K_c+1)],0)
    RC_tt<-array(diag(Num_assets),dim=c(Num_assets,Num_assets,TT))
    
    matrix_id<-matrix(1:Num_assets^2,ncol=Num_assets)
    matrix_id_2<-which(matrix_id==1:Num_assets^2, arr.ind=TRUE)
    
    for(i in 1:nrow(matrix_id_2)){
      RC_tt[matrix_id_2[i,1],matrix_id_2[i,2],]  <- suppressWarnings(
        roll_sum_func(C_t[matrix_id_2[i,1],matrix_id_2[i,2],], c(K_c+1),weights = rc_betas)
      ) 
    }
    
    
    
    z_c <- m_tensor+theta_tensor*RC_tt
    
    R_t_bar<-tanh(z_c)
    
    ll <-  as.numeric(midas_ll_rcpp(TT, Num_assets, K_c, a, b, R_t_bar, res))  
    
    # for (i in 1:TT) {
    #   diag(R_t_bar[,,i]) <- 1
    # }
    # 
    # ################# likelihood
    # 
    # ll<-rep(0,TT)
    # 
    # for(tt in (K_c+1):TT){
    #   Q_t[,,tt]<-(1-a-b)*R_t_bar[,,tt] + a*tcrossprod(res[tt-1, ]) + b*Q_t[,,tt-1]
    #   
    #   q_sd <- sqrt(diag(Q_t[,,tt]))
    #   R_t[,,tt] <- Q_t[,,tt] / tcrossprod(q_sd)
    #   
    #   chol_R <- tryCatch(chol(R_t[,,tt]), error = function(e) return(NULL)) # 进行 Cholesky 分解 (得到上三角矩阵)
    #   # if (is.null(chol_R)) {
    #   #   return(rep(-1e10, TT)) # 如果非正定，返回一个惩罚值代替似然值
    #   # }
    #   
    #   log_det_R_t[tt] <- 2 * sum(log(diag(chol_R))) # 等价于 log_det_R_t[tt]<-log(Det(R_t[,,tt]))
    #   
    #   tmp <- backsolve(chol_R, res[tt, ], transpose = TRUE)
    #   Eps_t_R_t_Eps_t[tt] <- sum(tmp^2)
    #   
    #   
    # }
    # 
    # ll<- - 0.5*(log_det_R_t+Eps_t_R_t_Eps_t)
    
    #sum(ll)
    
    return(ll)
    
  }else if (model=="DCCMIDAS_X_sb"){
    
    
    
    ##### index and param
    Num_assets<-ncol(res)		# number of assets
    TT<-length(as.Date(days_period))			# number of daily observations
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    
    # 预先引用 roll_sum 避免多次查找 定义在fit函数开头里
    #roll_sum_func <- roll::roll_sum
    #weight_func_beta <- rumidas::beta_function
    #weight_func_exp <- rumidas::exp_almon
    
    a<-param[1]
    b<-param[2]
    w1<- ifelse(lag_fun=="Beta",1,0)
    
    param_w2_c <- param[3:(length_dccmidas_param+2)] 
    
    
    midasdcc_param <- list(w2_c = param_w2_c)
    midasdcc_param_matrix <- list()
    
    ##### matrices and vectors
    for (pname in names(midasdcc_param)) {
      
      off_diag_elements <- midasdcc_param[[pname]]  
      cov_matrix <- matrix(0, nrow = Num_assets, ncol = Num_assets)
      cov_matrix[lower.tri(cov_matrix)] <- off_diag_elements
      cov_matrix[upper.tri(cov_matrix)] <- t(cov_matrix)[upper.tri(cov_matrix)]
      midasdcc_param_matrix[[pname]] <- cov_matrix
      
    }
    
    
    betamidas_c_matrix_tt <- array(0,dim=c(Num_assets,Num_assets,TT))
    
    Q_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    R_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    log_det_R_t<-rep(0,TT)
    Eps_t_R_t_Eps_t<-rep(0,TT)
    betas <- array(0,dim=c(Num_assets,Num_assets,(K_c+1)))
    
    ########structural breaks
    # 提取当前资产对的系数向量 (长度 1+n_n)
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(days_period)
    curr_m_vec <- param[(length_dccmidas_param+3):(length_dccmidas_param+n_n+3)]  
    curr_theta_vec <- param[(length_dccmidas_param+n_n+4):(length_dccmidas_param+2*n_n+4)] 
    
    D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
    D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
    
    for (i in 1:n_n) {
      # 找到处于第 i 个节点之后的日期，设为 1
      D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
    }
    
    # --- 步骤 2: 准备参数矩阵 (从优化的向量中提取) ---
    
    # --- 步骤 3: 矩阵乘法计算并生成张量 ---
    # 初始化目标张量 [N, N, T]
    m_tensor <- array(0, dim = c(Num_assets, Num_assets, TT))
    theta_tensor <- array(0, dim = c(Num_assets, Num_assets, TT))
    pairs_idx <- which(lower.tri(matrix(0, Num_assets, Num_assets)), arr.ind = TRUE)
    
    for (k in 1:length_dccmidas_param) {
      i_idx <- pairs_idx[k, 1]
      j_idx <- pairs_idx[k, 2]
      
      # 核心计算：利用矩阵乘法完成所有时间点的分配
      # 第一部分：截距项时间序列 (Base + Increments)
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      
      # 第二部分：斜率项时间序列 (Base + Increments)
      theta_t_series <- as.vector(curr_theta_vec %*% D_mat)
      
      # 填充入张量的对称位置
      m_tensor[i_idx, j_idx, ] <- m_t_series
      m_tensor[j_idx, i_idx, ] <- m_t_series
      theta_tensor[i_idx, j_idx, ] <- theta_t_series
      theta_tensor[j_idx, i_idx, ] <- theta_t_series
    }
    
    #### R_t_bar
    weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
    
    matrix_id<-matrix(1:Num_assets^2,ncol=Num_assets)
    matrix_id_2<-which(matrix_id==1:Num_assets^2, arr.ind=TRUE)
    
    for(i in 1:nrow(matrix_id_2)){
      
      betas[matrix_id_2[i,1],matrix_id_2[i,2],]<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,midasdcc_param_matrix[["w2_c"]][matrix_id_2[i,1],matrix_id_2[i,2]]))[2:(K_c+1)],0)
      midas_c <- suppressWarnings(roll_sum_func(mv_c, c(K_c+1),weights = betas[matrix_id_2[i,1],matrix_id_2[i,2],]))
      betamidas_c_matrix_tt[matrix_id_2[i,1],matrix_id_2[i,2],] <-midas_c[(K_c+1),]
      
    }
    
    z_c <- m_tensor+theta_tensor*betamidas_c_matrix_tt
    R_t_bar<-tanh(z_c)
    
    ll <-  as.numeric(midas_ll_rcpp(TT, Num_assets, K_c, a, b, R_t_bar, res))  
    
    # for (i in 1:TT) {
    #   diag(R_t_bar[,,i]) <- 1
    # }
    # 
    # ################# likelihood
    # 
    # ll<-rep(0,TT)
    # 
    # for(tt in (K_c+1):TT){
    #   
    #   Q_t[,,tt]<-(1-a-b)*R_t_bar[,,tt] + a*tcrossprod(res[tt-1, ]) + b*Q_t[,,tt-1]
    #   q_sd <- sqrt(diag(Q_t[,,tt]))
    #   R_t[,,tt] <- Q_t[,,tt] / tcrossprod(q_sd)
    #   
    #   chol_R <- tryCatch(chol(R_t[,,tt]), error = function(e) return(NULL)) # 进行 Cholesky 分解 (得到上三角矩阵)
    #   # if (is.null(chol_R)) {
    #   #   return(rep(-1e10, TT)) # 如果非正定，返回一个惩罚值代替似然值
    #   # }
    #   log_det_R_t[tt] <- 2 * sum(log(diag(chol_R))) # 等价于 log_det_R_t[tt]<-log(Det(R_t[,,tt]))
    #   
    #   tmp <- backsolve(chol_R, res[tt, ], transpose = TRUE)
    #   Eps_t_R_t_Eps_t[tt] <- sum(tmp^2)
    #   
    # }
    # 
    # ll<- - 0.5*(log_det_R_t+Eps_t_R_t_Eps_t)
    
    return(ll)
    
  }else if (model=="DCCMIDAS_RC_X_sb"){
    
    ##### index and param
    Num_assets<-ncol(res)		# number of assets
    TT<-length(as.Date(days_period))			# number of daily observations
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    
    a<-param[1]
    b<-param[2]
    w1<- ifelse(lag_fun=="Beta",1,0)
    w2<-param[3]
    
    # 预先引用 roll_sum 避免多次查找 定义在fit函数开头里
    #roll_sum_func <- roll::roll_sum
    #weight_func_beta <- rumidas::beta_function
    #weight_func_exp <- rumidas::exp_almon
    
    ##### matrices and vectors
    ######. X
    param_w2_c <- param[4:(length_dccmidas_param+3)]
    
    midasdcc_param <- list(w2_c = param_w2_c)
    midasdcc_param_matrix <- list()
    
    for (pname in names(midasdcc_param)) {
      
      off_diag_elements <- midasdcc_param[[pname]]  
      cov_matrix <- matrix(0, nrow = Num_assets, ncol = Num_assets)
      cov_matrix[lower.tri(cov_matrix)] <- off_diag_elements
      cov_matrix[upper.tri(cov_matrix)] <- t(cov_matrix)[upper.tri(cov_matrix)]
      midasdcc_param_matrix[[pname]] <- cov_matrix
      
    }
    
    
    betamidas_c_matrix_tt <- array(0,dim=c(Num_assets,Num_assets,TT))
    
    betas <- array(0,dim=c(Num_assets,Num_assets,(K_c+1)))
    
    ##################### RC
    
    RC_tt<-array(diag(Num_assets),dim=c(Num_assets,Num_assets,TT))
    
    #####################
    
    Q_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    R_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    
    log_det_R_t<-rep(0,TT)
    Eps_t_R_t_Eps_t<-rep(0,TT)
    
    S<-stats::cov(res)
    H_t<-array(S,dim=c(Num_assets,Num_assets,TT))
    
    
    ########structural breaks
    # 提取当前资产对的系数向量 (长度 1+n_n)
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(days_period)
    curr_m_vec <- param[(length_dccmidas_param+4):(length_dccmidas_param+n_n+4)]  
    curr_theta_RC_vec <- param[(length_dccmidas_param+n_n+5):(length_dccmidas_param+2*n_n+5)] 
    curr_theta_c_vec <- param[(length_dccmidas_param+2*n_n+6):(length_dccmidas_param+3*n_n+6)]
    
    D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
    D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
    
    for (i in 1:n_n) {
      # 找到处于第 i 个节点之后的日期，设为 1
      D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
    }
    
    # --- 步骤 2: 准备参数矩阵 (从优化的向量中提取) ---
    
    # --- 步骤 3: 矩阵乘法计算并生成张量 ---
    # 初始化目标张量 [N, N, T]
    m_tensor <- array(0, dim = c(Num_assets, Num_assets, TT))
    theta_RC_tensor <- array(0, dim = c(Num_assets, Num_assets, TT))
    theta_c_tensor <- array(0, dim = c(Num_assets, Num_assets, TT))
    
    pairs_idx <- which(lower.tri(matrix(0, Num_assets, Num_assets)), arr.ind = TRUE)
    
    for (k in 1:length_dccmidas_param) {
      i_idx <- pairs_idx[k, 1]
      j_idx <- pairs_idx[k, 2]
      
      # 核心计算：利用矩阵乘法完成所有时间点的分配
      # 第一部分：截距项时间序列 (Base + Increments)
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      
      # 第二部分：斜率项时间序列 (Base + Increments)
      theta_RC_t_series <- as.vector(curr_theta_RC_vec %*% D_mat)
      
      # 第三部分
      theta_c_t_series <- as.vector(curr_theta_c_vec %*% D_mat)
      
      # 填充入张量的对称位置
      m_tensor[i_idx, j_idx, ] <- m_t_series
      m_tensor[j_idx, i_idx, ] <- m_t_series
      theta_RC_tensor[i_idx, j_idx, ] <- theta_RC_t_series 
      theta_RC_tensor[j_idx, i_idx, ] <- theta_RC_t_series 
      theta_c_tensor[i_idx, j_idx, ] <- theta_c_t_series
      theta_c_tensor[j_idx, i_idx, ] <- theta_c_t_series
    }
    
    
    ##### first cycle
    
    C_t <- midas_rc_rcpp(N_c, Num_assets, TT, res)
    # C_t<-array(0,dim=c(Num_assets,Num_assets,TT))
    # for(tt in (N_c+1):TT){
    #   
    #   window_res <- res[(tt-N_c):tt, ] 
    #   V_t_0.5 <- sqrt(colSums(window_res^2)) 
    #   Prod_eps_t<-crossprod(res[tt:(tt - N_c), ])
    #   C_t[,,tt] <- Prod_eps_t / tcrossprod(V_t_0.5)
    #   
    # }
    ################
    
    weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
    
    #### R_t_bar
    rc_betas<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,w2))[2:(K_c+1)],0)
    matrix_id<-matrix(1:Num_assets^2,ncol=Num_assets)
    matrix_id_2<-which(matrix_id==1:Num_assets^2, arr.ind=TRUE)
    
    for(i in 1:nrow(matrix_id_2)){
      RC_tt[matrix_id_2[i,1],matrix_id_2[i,2],]  <- suppressWarnings(
        roll_sum_func(C_t[matrix_id_2[i,1],matrix_id_2[i,2],], c(K_c+1),weights = rc_betas)) 
      
      betas[matrix_id_2[i,1],matrix_id_2[i,2],]<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,midasdcc_param_matrix[["w2_c"]][matrix_id_2[i,1],matrix_id_2[i,2]]))[2:(K_c+1)],0)
      midas_c <- suppressWarnings(roll_sum_func(mv_c, c(K_c+1),weights = betas[matrix_id_2[i,1],matrix_id_2[i,2],]))
      betamidas_c_matrix_tt[matrix_id_2[i,1],matrix_id_2[i,2],] <-midas_c[(K_c+1),]
      
    }
    
    
    z_c <- m_tensor+theta_c_tensor*betamidas_c_matrix_tt+theta_RC_tensor*RC_tt
    
    R_t_bar<-tanh(z_c)
    
    ll <-  as.numeric(midas_ll_rcpp(TT, Num_assets, K_c, a, b, R_t_bar, res))  
    
    # for (i in 1:TT) {
    #   diag(R_t_bar[,,i]) <- 1
    # }
    # 
    # ################# likelihood
    # 
    # ll<-rep(0,TT)
    # 
    # for(tt in (K_c+1):TT){
    #   Q_t[,,tt]<-(1-a-b)*R_t_bar[,,tt] + a*tcrossprod(res[tt-1, ]) + b*Q_t[,,tt-1]
    #   
    #   q_sd <- sqrt(diag(Q_t[,,tt]))
    #   R_t[,,tt] <- Q_t[,,tt] / tcrossprod(q_sd)
    #   
    #   chol_R <- tryCatch(chol(R_t[,,tt]), error = function(e) return(NULL)) # 进行 Cholesky 分解 (得到上三角矩阵)
    #   # if (is.null(chol_R)) {
    #   #   return(rep(-1e10, TT)) # 如果非正定，返回一个惩罚值代替似然值
    #   # }
    #   log_det_R_t[tt] <- 2 * sum(log(diag(chol_R))) # 等价于 log_det_R_t[tt]<-log(Det(R_t[,,tt]))
    #   
    #   tmp <- backsolve(chol_R, res[tt, ], transpose = TRUE)
    #   Eps_t_R_t_Eps_t[tt] <- sum(tmp^2)
    #   
    #   
    #   
    # }
    # 
    # 
    # ll<- - 0.5*(log_det_R_t+Eps_t_R_t_Eps_t)
    
    return(ll)
    
  }
  
  
}

#dccmidas_noskew_mat_est
dccmidas_noskew_mat_est <- function(param,model,res,days_period,D_t,mv_c=NULL,lag_fun="Beta",N_c=NULL,K_c,X_stru_break_date_vec=NULL,N_c_res_for_oos=NULL){
  
  
  
  if (model=="DCCMIDAS_RC"){
    
    a<-param[1]
    b<-param[2]
    w1<- ifelse(lag_fun=="Beta",1,0)
    w2<-param[3]
    
    ##### index
    
    Num_assets<-ncol(res)		# number of assets
    TT<-length(as.Date(days_period))			# number of daily observations
    
    ##### matrices and vectors
    
    C_t<-array(0,dim=c(Num_assets,Num_assets,TT))
    
    Q_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    R_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    
    S<-stats::cov(res)
    H_t<-array(S,dim=c(Num_assets,Num_assets,TT))
    
    ##### first cycle
    if (is.null(N_c_res_for_oos)) {
      
      
      for(tt in (N_c+1):TT){
        
        window_res <- res[(tt-N_c):tt, ] 
        V_t_0.5 <- sqrt(colSums(window_res^2)) 
        Prod_eps_t<-crossprod(res[tt:(tt - N_c), ])
        C_t[,,tt] <- Prod_eps_t / tcrossprod(V_t_0.5)
        
      }
      
    } else {
      
      first_cycle_res <- rbind(N_c_res_for_oos,res)
      
      for(tt in (N_c+1):(TT+N_c)){
        
        window_res <- first_cycle_res[(tt-N_c):tt, ] 
        V_t_0.5 <- sqrt(colSums(window_res^2)) 
        Prod_eps_t<-crossprod(first_cycle_res[tt:(tt - N_c), ])
        C_t[,,tt-N_c] <- Prod_eps_t / tcrossprod(V_t_0.5)
        
      }
    }
    
    
    
    #### R_t_bar
    
    
    weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
    
    rc_betas<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,w2))[2:(K_c+1)],0)
    R_t_bar<-array(diag(Num_assets),dim=c(Num_assets,Num_assets,TT))
    
    matrix_id<-matrix(1:Num_assets^2,ncol=Num_assets)
    matrix_id_2<-which(matrix_id==1:Num_assets^2, arr.ind=TRUE)
    
    for(i in 1:nrow(matrix_id_2)){
      R_t_bar[matrix_id_2[i,1],matrix_id_2[i,2],]  <- suppressWarnings(
        roll_sum_func(C_t[matrix_id_2[i,1],matrix_id_2[i,2],], c(K_c+1),weights = rc_betas)
      ) 
    }
    
    for (i in 1:TT) {
      
      diag(R_t_bar[,,i]) <- 1
    }
    
    
    for(tt in (K_c+1):TT){
      Q_t[,,tt]<-(1-a-b)*R_t_bar[,,tt] + a*tcrossprod(res[tt-1, ]) + b*Q_t[,,tt-1]
      
      q_sd <- sqrt(diag(Q_t[,,tt]))
      R_t[,,tt] <- Q_t[,,tt] / tcrossprod(q_sd)
      H_t[,,tt]<-D_t[,,tt]%*%R_t[,,tt]%*%D_t[,,tt]
    }
    
    
    
  }else if (model=="DCCMIDAS_X"){
    
    ##### index and param
    Num_assets<-ncol(res)		# number of assets
    TT<-length(as.Date(days_period))			# number of daily observations
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    
    # 预先引用 roll_sum 避免多次查找 定义在fit函数开头里
    #roll_sum_func <- roll::roll_sum
    #weight_func_beta <- rumidas::beta_function
    #weight_func_exp <- rumidas::exp_almon
    
    a<-param[1]
    b<-param[2]
    w1<- ifelse(lag_fun=="Beta",1,0)
    
    param_delta_c <- param[3:(length_dccmidas_param+2)]
    param_theta_c <- param[(length_dccmidas_param+3):(length_dccmidas_param*2+2)] 
    param_w2_c <- param[(length_dccmidas_param*2+3):(length_dccmidas_param*3+2)]
    
    
    midasdcc_param <- list(delta_c=param_delta_c,
                           theta_c =param_theta_c, 
                           w2_c = param_w2_c)
    midasdcc_param_matrix <- list()
    
    ##### matrices and vectors
    for (pname in names(midasdcc_param)) {
      
      off_diag_elements <- midasdcc_param[[pname]]  
      cov_matrix <- matrix(0, nrow = Num_assets, ncol = Num_assets)
      cov_matrix[lower.tri(cov_matrix)] <- off_diag_elements
      cov_matrix[upper.tri(cov_matrix)] <- t(cov_matrix)[upper.tri(cov_matrix)]
      midasdcc_param_matrix[[pname]] <- cov_matrix
      
    }
    
    delta_c_matrix_tt <- array(midasdcc_param_matrix[["delta_c"]], dim = c(Num_assets, Num_assets, TT))
    theta_c_matrix_tt <- array(midasdcc_param_matrix[["theta_c"]], dim = c(Num_assets, Num_assets, TT))
    betamidas_c_matrix_tt <- array(0,dim=c(Num_assets,Num_assets,TT))
    
    Q_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    R_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    
    betas <- array(0,dim=c(Num_assets,Num_assets,(K_c+1)))
    S<-stats::cov(res)
    H_t<-array(S,dim=c(Num_assets,Num_assets,TT))
    
    #### R_t_bar
    weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
    
    matrix_id<-matrix(1:Num_assets^2,ncol=Num_assets)
    matrix_id_2<-which(matrix_id==1:Num_assets^2, arr.ind=TRUE)
    
    for(i in 1:nrow(matrix_id_2)){
      
      betas[matrix_id_2[i,1],matrix_id_2[i,2],]<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,midasdcc_param_matrix[["w2_c"]][matrix_id_2[i,1],matrix_id_2[i,2]]))[2:(K_c+1)],0)
      midas_c <- suppressWarnings(roll_sum_func(mv_c, c(K_c+1),weights = betas[matrix_id_2[i,1],matrix_id_2[i,2],]))
      betamidas_c_matrix_tt[matrix_id_2[i,1],matrix_id_2[i,2],] <-midas_c[(K_c+1),]
      
    }
    
    z_c <- delta_c_matrix_tt+theta_c_matrix_tt*betamidas_c_matrix_tt
    R_t_bar<-tanh(z_c)
    
    for (i in 1:TT) {
      diag(R_t_bar[,,i]) <- 1
    }
    
    
    for(tt in (K_c+1):TT){
      
      Q_t[,,tt]<-(1-a-b)*R_t_bar[,,tt] + a*tcrossprod(res[tt-1, ]) + b*Q_t[,,tt-1]
      q_sd <- sqrt(diag(Q_t[,,tt]))
      R_t[,,tt] <- Q_t[,,tt] / tcrossprod(q_sd)
      H_t[,,tt]<-D_t[,,tt]%*%R_t[,,tt]%*%D_t[,,tt]
      
      
    }
    
  }else if (model=="DCCMIDAS_RC_X"){
    
    ##### index and param
    Num_assets<-ncol(res)		# number of assets
    TT<-length(as.Date(days_period))			# number of daily observations
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    
    a<-param[1]
    b<-param[2]
    w1<- ifelse(lag_fun=="Beta",1,0)
    w2<-param[3]
    
    # 预先引用 roll_sum 避免多次查找 定义在fit函数开头里
    #roll_sum_func <- roll::roll_sum
    #weight_func_beta <- rumidas::beta_function
    #weight_func_exp <- rumidas::exp_almon
    
    ##### matrices and vectors
    ######. X
    param_delta_c <- param[4:(length_dccmidas_param+3)]
    param_theta_c <- param[(length_dccmidas_param+4):(length_dccmidas_param*2+3)] 
    param_w2_c <- param[(length_dccmidas_param*2+4):(length_dccmidas_param*3+3)]
    
    param_theta_RC_c <- param[(length_dccmidas_param*3+4):(length_dccmidas_param*4+3)] 
    
    midasdcc_param <- list(delta_c=param_delta_c,
                           theta_c =param_theta_c, 
                           w2_c = param_w2_c,
                           theta_RC_c = param_theta_RC_c)
    midasdcc_param_matrix <- list()
    
    
    for (pname in names(midasdcc_param)) {
      
      off_diag_elements <- midasdcc_param[[pname]]  
      cov_matrix <- matrix(0, nrow = Num_assets, ncol = Num_assets)
      cov_matrix[lower.tri(cov_matrix)] <- off_diag_elements
      cov_matrix[upper.tri(cov_matrix)] <- t(cov_matrix)[upper.tri(cov_matrix)]
      midasdcc_param_matrix[[pname]] <- cov_matrix
      
    }
    
    delta_c_matrix_tt <- array(midasdcc_param_matrix[["delta_c"]], dim = c(Num_assets, Num_assets, TT))
    theta_c_matrix_tt <- array(midasdcc_param_matrix[["theta_c"]], dim = c(Num_assets, Num_assets, TT))
    betamidas_c_matrix_tt <- array(0,dim=c(Num_assets,Num_assets,TT))
    
    theta_RC_c_matrix_tt <- array(midasdcc_param_matrix[["theta_RC_c"]], dim = c(Num_assets, Num_assets, TT))
    
    
    betas <- array(0,dim=c(Num_assets,Num_assets,(K_c+1)))
    S<-stats::cov(res)
    H_t<-array(S,dim=c(Num_assets,Num_assets,TT))
    
    ##################### RC
    C_t<-array(0,dim=c(Num_assets,Num_assets,TT))
    RC_tt<-array(diag(Num_assets),dim=c(Num_assets,Num_assets,TT))
    
    #####################
    
    Q_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    R_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    
    ##### first cycle
    if (is.null(N_c_res_for_oos)) {
      
      
      for(tt in (N_c+1):TT){
        
        window_res <- res[(tt-N_c):tt, ] 
        V_t_0.5 <- sqrt(colSums(window_res^2)) 
        Prod_eps_t<-crossprod(res[tt:(tt - N_c), ])
        C_t[,,tt] <- Prod_eps_t / tcrossprod(V_t_0.5)
        
      }
      
    } else {
      
      first_cycle_res <- rbind(N_c_res_for_oos,res)
      
      for(tt in (N_c+1):(TT+N_c)){
        
        window_res <- first_cycle_res[(tt-N_c):tt, ] 
        V_t_0.5 <- sqrt(colSums(window_res^2)) 
        Prod_eps_t<-crossprod(first_cycle_res[tt:(tt - N_c), ])
        C_t[,,tt-N_c] <- Prod_eps_t / tcrossprod(V_t_0.5)
        
      }
    }
    ################
    
    weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
    
    #### R_t_bar
    rc_betas<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,w2))[2:(K_c+1)],0)
    matrix_id<-matrix(1:Num_assets^2,ncol=Num_assets)
    matrix_id_2<-which(matrix_id==1:Num_assets^2, arr.ind=TRUE)
    
    for(i in 1:nrow(matrix_id_2)){
      RC_tt[matrix_id_2[i,1],matrix_id_2[i,2],]  <- suppressWarnings(
        roll_sum_func(C_t[matrix_id_2[i,1],matrix_id_2[i,2],], c(K_c+1),weights = rc_betas)) 
      
      betas[matrix_id_2[i,1],matrix_id_2[i,2],]<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,midasdcc_param_matrix[["w2_c"]][matrix_id_2[i,1],matrix_id_2[i,2]]))[2:(K_c+1)],0)
      midas_c <- suppressWarnings(roll_sum_func(mv_c, c(K_c+1),weights = betas[matrix_id_2[i,1],matrix_id_2[i,2],]))
      betamidas_c_matrix_tt[matrix_id_2[i,1],matrix_id_2[i,2],] <-midas_c[(K_c+1),]
      
    }
    
    
    
    z_c <- delta_c_matrix_tt+theta_c_matrix_tt*betamidas_c_matrix_tt+theta_RC_c_matrix_tt*RC_tt
    
    R_t_bar<-tanh(z_c)
    
    for (i in 1:TT) {
      diag(R_t_bar[,,i]) <- 1
    }
    
    
    for(tt in (K_c+1):TT){
      Q_t[,,tt]<-(1-a-b)*R_t_bar[,,tt] + a*tcrossprod(res[tt-1, ]) + b*Q_t[,,tt-1]
      
      q_sd <- sqrt(diag(Q_t[,,tt]))
      R_t[,,tt] <- Q_t[,,tt] / tcrossprod(q_sd)
      H_t[,,tt]<-D_t[,,tt]%*%R_t[,,tt]%*%D_t[,,tt]
      
    }
    
    
  }else if (model=="DCCMIDAS_RC_sb"){
    
    a<-param[1]
    b<-param[2]
    w1<- ifelse(lag_fun=="Beta",1,0)
    w2<-param[3]
    
    ##### index
    
    Num_assets<-ncol(res)		# number of assets
    TT<-length(as.Date(days_period))			# number of daily observations
    
    ##### matrices and vectors
    
    C_t<-array(0,dim=c(Num_assets,Num_assets,TT))
    
    Q_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    R_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    S<-stats::cov(res)
    H_t<-array(S,dim=c(Num_assets,Num_assets,TT))
    
    ########structural breaks
    # 提取当前资产对的系数向量 (长度 1+n_n)
    
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(days_period)
    length_dccmidas_param <- 0.5 * Num_assets * (Num_assets - 1)
    curr_m_vec <- param[4:(n_n+4)]  
    curr_theta_vec <- param[(n_n+5):(2*n_n+5)] 
    
    D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
    D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
    
    for (i in 1:n_n) {
      # 找到处于第 i 个节点之后的日期，设为 1
      D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
    }
    
    # --- 步骤 2: 准备参数矩阵 (从优化的向量中提取) ---
    
    # --- 步骤 3: 矩阵乘法计算并生成张量 ---
    # 初始化目标张量 [N, N, T]
    m_tensor <- array(0, dim = c(Num_assets, Num_assets, TT))
    theta_tensor <- array(0, dim = c(Num_assets, Num_assets, TT))
    pairs_idx <- which(lower.tri(matrix(0, Num_assets, Num_assets)), arr.ind = TRUE)
    
    for (k in 1:length_dccmidas_param) {
      i_idx <- pairs_idx[k, 1]
      j_idx <- pairs_idx[k, 2]
      
      # 核心计算：利用矩阵乘法完成所有时间点的分配
      # 第一部分：截距项时间序列 (Base + Increments)
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      
      # 第二部分：斜率项时间序列 (Base + Increments)
      theta_t_series <- as.vector(curr_theta_vec %*% D_mat)
      
      # 填充入张量的对称位置
      m_tensor[i_idx, j_idx, ] <- m_t_series
      m_tensor[j_idx, i_idx, ] <- m_t_series
      theta_tensor[i_idx, j_idx, ] <- theta_t_series
      theta_tensor[j_idx, i_idx, ] <- theta_t_series
    }
    
    
    ##### first cycle
    if (is.null(N_c_res_for_oos)) {
      
      
      for(tt in (N_c+1):TT){
        
        window_res <- res[(tt-N_c):tt, ] 
        V_t_0.5 <- sqrt(colSums(window_res^2)) 
        Prod_eps_t<-crossprod(res[tt:(tt - N_c), ])
        C_t[,,tt] <- Prod_eps_t / tcrossprod(V_t_0.5)
        
      }
      
    } else {
      
      first_cycle_res <- rbind(N_c_res_for_oos,res)
      
      for(tt in (N_c+1):(TT+N_c)){
        
        window_res <- first_cycle_res[(tt-N_c):tt, ] 
        V_t_0.5 <- sqrt(colSums(window_res^2)) 
        Prod_eps_t<-crossprod(first_cycle_res[tt:(tt - N_c), ])
        C_t[,,tt-N_c] <- Prod_eps_t / tcrossprod(V_t_0.5)
        
      }
    }
    
    #### R_t_bar
    
    
    weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
    
    rc_betas<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,w2))[2:(K_c+1)],0)
    RC_tt<-array(diag(Num_assets),dim=c(Num_assets,Num_assets,TT))
    
    matrix_id<-matrix(1:Num_assets^2,ncol=Num_assets)
    matrix_id_2<-which(matrix_id==1:Num_assets^2, arr.ind=TRUE)
    
    for(i in 1:nrow(matrix_id_2)){
      RC_tt[matrix_id_2[i,1],matrix_id_2[i,2],]  <- suppressWarnings(
        roll_sum_func(C_t[matrix_id_2[i,1],matrix_id_2[i,2],], c(K_c+1),weights = rc_betas)
      ) 
    }
    
    
    
    z_c <- m_tensor+theta_tensor*RC_tt
    
    R_t_bar<-tanh(z_c)
    
    for (i in 1:TT) {
      diag(R_t_bar[,,i]) <- 1
    }
    
    for(tt in (K_c+1):TT){
      Q_t[,,tt]<-(1-a-b)*R_t_bar[,,tt] + a*tcrossprod(res[tt-1, ]) + b*Q_t[,,tt-1]
      
      q_sd <- sqrt(diag(Q_t[,,tt]))
      R_t[,,tt] <- Q_t[,,tt] / tcrossprod(q_sd)
      H_t[,,tt]<-D_t[,,tt]%*%R_t[,,tt]%*%D_t[,,tt]
      
      
    }
    
  }else if (model=="DCCMIDAS_X_sb"){
    
    ##### index and param
    Num_assets<-ncol(res)		# number of assets
    TT<-length(as.Date(days_period))			# number of daily observations
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    
    # 预先引用 roll_sum 避免多次查找 定义在fit函数开头里
    #roll_sum_func <- roll::roll_sum
    #weight_func_beta <- rumidas::beta_function
    #weight_func_exp <- rumidas::exp_almon
    
    a<-param[1]
    b<-param[2]
    w1<- ifelse(lag_fun=="Beta",1,0)
    
    param_w2_c <- param[3:(length_dccmidas_param+2)] 
    
    
    midasdcc_param <- list(w2_c = param_w2_c)
    midasdcc_param_matrix <- list()
    
    ##### matrices and vectors
    for (pname in names(midasdcc_param)) {
      
      off_diag_elements <- midasdcc_param[[pname]]  
      cov_matrix <- matrix(0, nrow = Num_assets, ncol = Num_assets)
      cov_matrix[lower.tri(cov_matrix)] <- off_diag_elements
      cov_matrix[upper.tri(cov_matrix)] <- t(cov_matrix)[upper.tri(cov_matrix)]
      midasdcc_param_matrix[[pname]] <- cov_matrix
      
    }
    
    
    betamidas_c_matrix_tt <- array(0,dim=c(Num_assets,Num_assets,TT))
    
    Q_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    R_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    
    S<-stats::cov(res)
    H_t<-array(S,dim=c(Num_assets,Num_assets,TT))
    
    betas <- array(0,dim=c(Num_assets,Num_assets,(K_c+1)))
    
    ########structural breaks
    # 提取当前资产对的系数向量 (长度 1+n_n)
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(days_period)
    curr_m_vec <- param[(length_dccmidas_param+3):(length_dccmidas_param+n_n+3)]  
    curr_theta_vec <- param[(length_dccmidas_param+n_n+4):(length_dccmidas_param+2*n_n+4)] 
    
    D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
    D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
    
    for (i in 1:n_n) {
      # 找到处于第 i 个节点之后的日期，设为 1
      D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
    }
    
    # --- 步骤 2: 准备参数矩阵 (从优化的向量中提取) ---
    
    # --- 步骤 3: 矩阵乘法计算并生成张量 ---
    # 初始化目标张量 [N, N, T]
    m_tensor <- array(0, dim = c(Num_assets, Num_assets, TT))
    theta_tensor <- array(0, dim = c(Num_assets, Num_assets, TT))
    pairs_idx <- which(lower.tri(matrix(0, Num_assets, Num_assets)), arr.ind = TRUE)
    
    for (k in 1:length_dccmidas_param) {
      i_idx <- pairs_idx[k, 1]
      j_idx <- pairs_idx[k, 2]
      
      # 核心计算：利用矩阵乘法完成所有时间点的分配
      # 第一部分：截距项时间序列 (Base + Increments)
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      
      # 第二部分：斜率项时间序列 (Base + Increments)
      theta_t_series <- as.vector(curr_theta_vec %*% D_mat)
      
      # 填充入张量的对称位置
      m_tensor[i_idx, j_idx, ] <- m_t_series
      m_tensor[j_idx, i_idx, ] <- m_t_series
      theta_tensor[i_idx, j_idx, ] <- theta_t_series
      theta_tensor[j_idx, i_idx, ] <- theta_t_series
    }
    
    #### R_t_bar
    weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
    
    matrix_id<-matrix(1:Num_assets^2,ncol=Num_assets)
    matrix_id_2<-which(matrix_id==1:Num_assets^2, arr.ind=TRUE)
    
    for(i in 1:nrow(matrix_id_2)){
      
      betas[matrix_id_2[i,1],matrix_id_2[i,2],]<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,midasdcc_param_matrix[["w2_c"]][matrix_id_2[i,1],matrix_id_2[i,2]]))[2:(K_c+1)],0)
      midas_c <- suppressWarnings(roll_sum_func(mv_c, c(K_c+1),weights = betas[matrix_id_2[i,1],matrix_id_2[i,2],]))
      betamidas_c_matrix_tt[matrix_id_2[i,1],matrix_id_2[i,2],] <-midas_c[(K_c+1),]
      
    }
    
    z_c <- m_tensor+theta_tensor*betamidas_c_matrix_tt
    R_t_bar<-tanh(z_c)
    
    for (i in 1:TT) {
      diag(R_t_bar[,,i]) <- 1
    }
    
    
    
    for(tt in (K_c+1):TT){
      
      Q_t[,,tt]<-(1-a-b)*R_t_bar[,,tt] + a*tcrossprod(res[tt-1, ]) + b*Q_t[,,tt-1]
      q_sd <- sqrt(diag(Q_t[,,tt]))
      R_t[,,tt] <- Q_t[,,tt] / tcrossprod(q_sd)
      H_t[,,tt]<-D_t[,,tt]%*%R_t[,,tt]%*%D_t[,,tt]
    }
    
    
    
    
  }else if (model=="DCCMIDAS_RC_X_sb"){
    
    ##### index and param
    Num_assets<-ncol(res)		# number of assets
    TT<-length(as.Date(days_period))			# number of daily observations
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    
    a<-param[1]
    b<-param[2]
    w1<- ifelse(lag_fun=="Beta",1,0)
    w2<-param[3]
    
    # 预先引用 roll_sum 避免多次查找 定义在fit函数开头里
    #roll_sum_func <- roll::roll_sum
    #weight_func_beta <- rumidas::beta_function
    #weight_func_exp <- rumidas::exp_almon
    
    ##### matrices and vectors
    ######. X
    param_w2_c <- param[4:(length_dccmidas_param+3)]
    
    midasdcc_param <- list(w2_c = param_w2_c)
    midasdcc_param_matrix <- list()
    
    for (pname in names(midasdcc_param)) {
      
      off_diag_elements <- midasdcc_param[[pname]]  
      cov_matrix <- matrix(0, nrow = Num_assets, ncol = Num_assets)
      cov_matrix[lower.tri(cov_matrix)] <- off_diag_elements
      cov_matrix[upper.tri(cov_matrix)] <- t(cov_matrix)[upper.tri(cov_matrix)]
      midasdcc_param_matrix[[pname]] <- cov_matrix
      
    }
    
    
    betamidas_c_matrix_tt <- array(0,dim=c(Num_assets,Num_assets,TT))
    
    betas <- array(0,dim=c(Num_assets,Num_assets,(K_c+1)))
    
    ##################### RC
    C_t<-array(0,dim=c(Num_assets,Num_assets,TT))
    RC_tt<-array(diag(Num_assets),dim=c(Num_assets,Num_assets,TT))
    
    #####################
    
    Q_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    R_t<-array(diag(rep(1,Num_assets)),dim=c(Num_assets,Num_assets,TT))
    S<-stats::cov(res)
    H_t<-array(S,dim=c(Num_assets,Num_assets,TT))
    
    
    ########structural breaks
    # 提取当前资产对的系数向量 (长度 1+n_n)
    
    n_n <- length(X_stru_break_date_vec)
    full_days <- as.Date(days_period)
    curr_m_vec <- param[(length_dccmidas_param+4):(length_dccmidas_param+n_n+4)]  
    curr_theta_RC_vec <- param[(length_dccmidas_param+n_n+5):(length_dccmidas_param+2*n_n+5)] 
    curr_theta_c_vec <- param[(length_dccmidas_param+2*n_n+6):(length_dccmidas_param+3*n_n+6)]
    
    D_mat <- matrix(0, nrow = 1 + n_n, ncol = TT)
    D_mat[1, ] <- 1  # 第一行作为 Constant (对应 mc 和 theta_c)
    
    for (i in 1:n_n) {
      # 找到处于第 i 个节点之后的日期，设为 1
      D_mat[i + 1, full_days >= as.Date(X_stru_break_date_vec)[i]] <- 1
    }
    
    # --- 步骤 2: 准备参数矩阵 (从优化的向量中提取) ---
    
    # --- 步骤 3: 矩阵乘法计算并生成张量 ---
    # 初始化目标张量 [N, N, T]
    m_tensor <- array(0, dim = c(Num_assets, Num_assets, TT))
    theta_RC_tensor <- array(0, dim = c(Num_assets, Num_assets, TT))
    theta_c_tensor <- array(0, dim = c(Num_assets, Num_assets, TT))
    
    pairs_idx <- which(lower.tri(matrix(0, Num_assets, Num_assets)), arr.ind = TRUE)
    
    for (k in 1:length_dccmidas_param) {
      i_idx <- pairs_idx[k, 1]
      j_idx <- pairs_idx[k, 2]
      
      # 核心计算：利用矩阵乘法完成所有时间点的分配
      # 第一部分：截距项时间序列 (Base + Increments)
      m_t_series <- as.vector(curr_m_vec %*% D_mat)
      
      # 第二部分：斜率项时间序列 (Base + Increments)
      theta_RC_t_series <- as.vector(curr_theta_RC_vec %*% D_mat)
      
      # 第三部分
      theta_c_t_series <- as.vector(curr_theta_c_vec %*% D_mat)
      
      # 填充入张量的对称位置
      m_tensor[i_idx, j_idx, ] <- m_t_series
      m_tensor[j_idx, i_idx, ] <- m_t_series
      theta_RC_tensor[i_idx, j_idx, ] <- theta_RC_t_series 
      theta_RC_tensor[j_idx, i_idx, ] <- theta_RC_t_series 
      theta_c_tensor[i_idx, j_idx, ] <- theta_c_t_series
      theta_c_tensor[j_idx, i_idx, ] <- theta_c_t_series
    }
    
    
    ##### first cycle
    if (is.null(N_c_res_for_oos)) {
      
      
      for(tt in (N_c+1):TT){
        
        window_res <- res[(tt-N_c):tt, ] 
        V_t_0.5 <- sqrt(colSums(window_res^2)) 
        Prod_eps_t<-crossprod(res[tt:(tt - N_c), ])
        C_t[,,tt] <- Prod_eps_t / tcrossprod(V_t_0.5)
        
      }
      
    } else {
      
      first_cycle_res <- rbind(N_c_res_for_oos,res)
      
      for(tt in (N_c+1):(TT+N_c)){
        
        window_res <- first_cycle_res[(tt-N_c):tt, ] 
        V_t_0.5 <- sqrt(colSums(window_res^2)) 
        Prod_eps_t<-crossprod(first_cycle_res[tt:(tt - N_c), ])
        C_t[,,tt-N_c] <- Prod_eps_t / tcrossprod(V_t_0.5)
        
      }
    }
    ################
    
    weight_fun<-ifelse(lag_fun=="Beta",weight_func_beta,weight_func_exp)
    
    #### R_t_bar
    rc_betas<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,w2))[2:(K_c+1)],0)
    matrix_id<-matrix(1:Num_assets^2,ncol=Num_assets)
    matrix_id_2<-which(matrix_id==1:Num_assets^2, arr.ind=TRUE)
    
    for(i in 1:nrow(matrix_id_2)){
      RC_tt[matrix_id_2[i,1],matrix_id_2[i,2],]  <- suppressWarnings(
        roll_sum_func(C_t[matrix_id_2[i,1],matrix_id_2[i,2],], c(K_c+1),weights = rc_betas)) 
      
      betas[matrix_id_2[i,1],matrix_id_2[i,2],]<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,midasdcc_param_matrix[["w2_c"]][matrix_id_2[i,1],matrix_id_2[i,2]]))[2:(K_c+1)],0)
      midas_c <- suppressWarnings(roll_sum_func(mv_c, c(K_c+1),weights = betas[matrix_id_2[i,1],matrix_id_2[i,2],]))
      betamidas_c_matrix_tt[matrix_id_2[i,1],matrix_id_2[i,2],] <-midas_c[(K_c+1),]
      
    }
    
    
    
    z_c <- m_tensor+theta_c_tensor*betamidas_c_matrix_tt+theta_RC_tensor*RC_tt
    
    R_t_bar<-tanh(z_c)
    
    for (i in 1:TT) {
      diag(R_t_bar[,,i]) <- 1
    }
    
    
    for(tt in (K_c+1):TT){
      Q_t[,,tt]<-(1-a-b)*R_t_bar[,,tt] + a*tcrossprod(res[tt-1, ]) + b*Q_t[,,tt-1]
      
      q_sd <- sqrt(diag(Q_t[,,tt]))
      R_t[,,tt] <- Q_t[,,tt] / tcrossprod(q_sd)
      H_t[,,tt]<-D_t[,,tt]%*%R_t[,,tt]%*%D_t[,,tt]
      
    }
  }
  
  results<-list(
    "H_t"=H_t,
    "R_t"=R_t,
    "R_t_bar"=R_t_bar
  )
  
  return(results)
  
}

dccmidas_noskew_longfocus_fit<-function(r_t,univ_model="sGARCH",distribution="norm",X_stru_break_date_vec=NULL,
                                        RV_m=NULL,mv_m=NULL,K=NULL,mv_c=NULL,corr_model="cDCC",lag_fun=NULL,
                                        N_c=NULL,K_c=NULL,out_of_sample=NULL,vol_proxy=NULL,R=NULL){
  
  if (is.null(lag_fun)) {
    
    lag_fun <- "Beta"
  }
  
  if(!is.null(out_of_sample)&&!is.null(N_c)&&out_of_sample < K_c) { stop(
    cat("#Warning:\n Parameter 'out_of_sample' must be larger than K_c \n")
  )}
  
  ################################### checks
  
  # 预先引用 roll_sum 避免多次查找
  roll_sum_func <- roll::roll_sum
  weight_func_beta <- rumidas::beta_function
  weight_func_exp <- rumidas::exp_almon
  
  if(!is.null(RV_m)&&class(RV_m) != "list") { stop(
    cat("#Warning:\n Parameter 'RV_m' must be a list object. Please provide it in the correct form \n")
  )}
  if(!is.null(mv_m)&&class(mv_m) != "list") { stop(
    cat("#Warning:\n Parameter 'mv_m' must be a list object. Please provide it in the correct form \n")
  )}
  
  if(class(r_t)[1] != "xts") { stop(
    cat("#Warning:\n Parameter 'r_t' must be a xts object. Please provide it in the correct form \n")
  )}
  
  if(!is.null(vol_proxy)&&class(vol_proxy)[1]!="xts") { stop(
    cat("#Warning:\n Parameter 'vol_proxy' must be a xts object. Please provide it in the correct form \n")
  )}
  
  if ((univ_model %in% c("GM_noskew_RV", "GM_noskew_X", "GM_noskew_RV_X", "GM_noskew_RV_sb", "GM_noskew_X_sb", "GM_noskew_RV_X_sb")) && is.null(K))
  { stop(
    cat("#Warning:\n If you want to estimate a GARCH-MIDAS model, please provide parameters  K in the correct form \n")
  )}
  
  if(
    (corr_model%in% c("DCCMIDAS_RC", "DCCMIDAS_RC_sb"))&&(is.null(N_c)|is.null(K_c))
  ) { stop(
    cat("#Warning:\n If you want to estimate a DCC-MIDAS model, please provide parameters N_c and K_c \n")
  )}
  if(
    (corr_model%in% c("DCCMIDAS_X", "DCCMIDAS_X_sb"))&&(is.null(mv_c)|is.null(K_c))
  ) { stop(
    cat("#Warning:\n If you want to estimate a DCC-MIDAS-X model, please provide parameters mv_c and K_c \n")
  )}
  if(
    (corr_model%in% c("DCCMIDAS_RC_X", "DCCMIDAS_RC_X_sb"))&(is.null(N_c)|is.null(mv_c)|is.null(K_c))
  ) { stop(
    cat("#Warning:\n If you want to estimate a DCC-MIDAS-RC-X model, please provide parameters N_c and mv_c and K_c \n")
  )}
  
  if(
    (univ_model %in% c("GM_noskew_RV_sb", "GM_noskew_X_sb", "GM_noskew_RV_X_sb")|corr_model%in% c("DCCMIDAS_RV_sb", "DCCMIDAS_X_sb", "DCCMIDAS_RV_X_sb"))&&(is.null(X_stru_break_date_vec))
  ) { stop(
    cat("#Warning:\n If you want to introduce structural breaks based on X, please provide parameters X_stru_break_date_vec \n")
  )}
  
  if (is.null(R)) {
    
    R <- 100
    
  }
  
  
  ################################### merge (eventually)
  
  Num_assets<- ncol(r_t)
  
  db <- r_t
  
  ################################### in-sample and out-of-sample
  
  TT<-nrow(db)
  
  if (is.null(out_of_sample)){ # out_of_sample parameter is is.null
    
    sd_m<-matrix(NA,ncol=Num_assets,nrow=TT)
    sd_long_m <- matrix(NA,ncol=Num_assets,nrow=TT)
    sd_short_m <- matrix(NA,ncol=Num_assets,nrow=TT)
    
    mv_c_in_s <- mv_c
    
    days_period<-index(db) 
    estim_period<-range(days_period)
    
    days_period_oos <- NA
    
    est_obs<-nrow(db)
    
  } else {					# out_of_sample parameter is present
    
    sd_m<-matrix(NA,ncol=Num_assets,nrow=TT-out_of_sample)
    sd_long_m <- matrix(NA,ncol=Num_assets,nrow=TT-out_of_sample)
    sd_short_m <- matrix(NA,ncol=Num_assets,nrow=TT-out_of_sample)
    
    db_oos<-db[(TT-out_of_sample+1):TT,]
    sd_m_oos<-matrix(NA,ncol=Num_assets,nrow=out_of_sample)
    sd_long_m_oos<-matrix(NA,ncol=Num_assets,nrow=out_of_sample)
    sd_short_m_oos<-matrix(NA,ncol=Num_assets,nrow=out_of_sample)
    
    ###for dccmidas part
    
    db_in_s <- db[1:(TT-out_of_sample),]
    
    mv_c_in_s <- mv_c[,1:(TT-out_of_sample)]
    mv_c_oos <- mv_c[,(TT-out_of_sample+1):TT]
    
    days_period_oos<-index(db)[(TT-out_of_sample+1):TT]
    
    days_period <- index(db)[1:(TT-out_of_sample)]
    estim_period<-range(days_period)
    est_obs<-nrow(db[1:(TT-out_of_sample),])
    
  }
  
  
  ################################### setting the estimation list
  
  u_est <- list()
  
  est_details<-list()
  
  u_loglik<-list()
  
  u_Inf_criteria <- list()
  
  start<-Sys.time()
  
  ########################################### first step (univariate GARCH)
  
  if(univ_model=="sGARCH"){
    
    uspec<-rugarch::ugarchspec(
      variance.model=list(model="sGARCH", garchOrder=c(1,1)),
      mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
      distribution.model = distribution)
    
    
    if (is.null(out_of_sample)){
      
      for(i in 1:Num_assets){
        u_est[[i]]<-rugarch::ugarchfit(spec=uspec,data=db[,i])
        est_details[[i]]<-u_est[[i]]@fit$robust.matcoef
        u_loglik[[i]]<-u_est[[i]]@fit$LLH
        u_Inf_criteria[[i]]<-Inf_criteria(u_est[[i]])
        sd_m[,i]<-u_est[[i]]@fit$sigma
        db[,i] <- db[,i]-mean(db[,i], na.rm = TRUE)
        
      }
      
    } else {
      
      for(i in 1:Num_assets){
        u_est[[i]]<-rugarch::ugarchfit(spec=uspec,data=db[,i],out.sample=out_of_sample)
        est_details[[i]]<-u_est[[i]]@fit$robust.matcoef
        u_loglik[[i]]<-u_est[[i]]@fit$LLH
        u_Inf_criteria[[i]]<-Inf_criteria(u_est[[i]])
        sd_m[,i]<-u_est[[i]]@fit$sigma
        sd_m_oos[,i]<-as.numeric(rugarch::sigma(rugarch::ugarchforecast(u_est[[i]],  n.ahead = 1, 
                                                                        n.roll = c(out_of_sample-1), 
                                                                        out.sample = c(out_of_sample-1))))
        db_in_s[,i] <- db_in_s[,i]-mean(db[,i], na.rm = TRUE)
        db_oos[,i] <- db_oos[,i]-mean(db[,i], na.rm = TRUE)
      }
      
      
    } 
  } else if (univ_model=="gjrGARCH"){
    
    uspec<-rugarch::ugarchspec(
      variance.model=list(model="gjrGARCH", garchOrder=c(1,1)),
      mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
      distribution.model = distribution)
    
    if (is.null(out_of_sample)){
      
      for(i in 1:Num_assets){
        u_est[[i]]<-rugarch::ugarchfit(spec=uspec,data=db[,i])
        est_details[[i]]<-u_est[[i]]@fit$robust.matcoef
        u_loglik[[i]]<-u_est[[i]]@fit$LLH
        u_Inf_criteria[[i]]<-Inf_criteria(u_est[[i]])
        sd_m[,i]<-u_est[[i]]@fit$sigma
        db[,i] <- db[,i]-mean(db[,i], na.rm = TRUE)
      }
      
    } else {
      
      for(i in 1:Num_assets){
        u_est[[i]]<-rugarch::ugarchfit(spec=uspec,data=db[,i],out.sample=out_of_sample)
        est_details[[i]]<-u_est[[i]]@fit$robust.matcoef
        u_loglik[[i]]<-u_est[[i]]@fit$LLH
        u_Inf_criteria[[i]]<-Inf_criteria(u_est[[i]])
        sd_m[,i]<-u_est[[i]]@fit$sigma
        sd_m_oos[,i]<-as.numeric(rugarch::sigma(rugarch::ugarchforecast(u_est[[i]],  n.ahead = 1, 
                                                                        n.roll = c(out_of_sample-1), 
                                                                        out.sample = c(out_of_sample-1))))
        db_in_s[,i] <- db_in_s[,i]-mean(db[,i], na.rm = TRUE)
        db_oos[,i] <- db_oos[,i]-mean(db[,i], na.rm = TRUE)
      }
      
      
    }
  } else if (univ_model=="eGARCH"){
    
    uspec<-rugarch::ugarchspec(
      variance.model=list(model="eGARCH", garchOrder=c(1,1)),
      mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
      distribution.model = distribution)
    
    if (is.null(out_of_sample)){
      
      for(i in 1:Num_assets){
        u_est[[i]]<-rugarch::ugarchfit(spec=uspec,data=db[,i])
        est_details[[i]]<-u_est[[i]]@fit$robust.matcoef
        u_loglik[[i]]<-u_est[[i]]@fit$LLH
        u_Inf_criteria[[i]]<-Inf_criteria(u_est[[i]])
        sd_m[,i]<-u_est[[i]]@fit$sigma
        db[,i] <- db[,i]-mean(db[,i], na.rm = TRUE)
      }
      
    } else {
      
      for(i in 1:Num_assets){
        u_est[[i]]<-rugarch::ugarchfit(spec=uspec,data=db[,i],out.sample=out_of_sample)
        est_details[[i]]<-u_est[[i]]@fit$robust.matcoef
        u_loglik[[i]]<-u_est[[i]]@fit$LLH
        u_Inf_criteria[[i]]<-Inf_criteria(u_est[[i]])
        sd_m[,i]<-u_est[[i]]@fit$sigma
        sd_m_oos[,i]<-as.numeric(rugarch::sigma(rugarch::ugarchforecast(u_est[[i]],  n.ahead = 1, 
                                                                        n.roll = c(out_of_sample-1), 
                                                                        out.sample = c(out_of_sample-1))))
        db_in_s[,i] <- db_in_s[,i]-mean(db[,i], na.rm = TRUE)
        db_oos[,i] <- db_oos[,i]-mean(db[,i], na.rm = TRUE)
      }
      
      
    }
  } else if (univ_model=="iGARCH"){
    
    uspec<-rugarch::ugarchspec(
      variance.model=list(model="iGARCH", garchOrder=c(1,1)),
      mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
      distribution.model = distribution)
    
    if (is.null(out_of_sample)){
      
      for(i in 1:Num_assets){
        u_est[[i]]<-rugarch::ugarchfit(spec=uspec,data=db[,i])
        est_details[[i]]<-u_est[[i]]@fit$robust.matcoef
        u_loglik[[i]]<-u_est[[i]]@fit$LLH
        u_Inf_criteria[[i]]<-Inf_criteria(u_est[[i]])
        sd_m[,i]<-u_est[[i]]@fit$sigma
        db[,i] <- db[,i]-mean(db[,i], na.rm = TRUE)
      }
      
    } else {
      
      for(i in 1:Num_assets){
        u_est[[i]]<-rugarch::ugarchfit(spec=uspec,data=db[,i],out.sample=out_of_sample)
        est_details[[i]]<-u_est[[i]]@fit$robust.matcoef
        u_loglik[[i]]<-u_est[[i]]@fit$LLH
        u_Inf_criteria[[i]]<-Inf_criteria(u_est[[i]])
        sd_m[,i]<-u_est[[i]]@fit$sigma
        sd_m_oos[,i]<-as.numeric(rugarch::sigma(rugarch::ugarchforecast(u_est[[i]],  n.ahead = 1, 
                                                                        n.roll = c(out_of_sample-1), 
                                                                        out.sample = c(out_of_sample-1))))
        db_in_s[,i] <- db_in_s[,i]-mean(db[,i], na.rm = TRUE)
        db_oos[,i] <- db_oos[,i]-mean(db[,i], na.rm = TRUE)
      }
      
      
    }
  } else if (univ_model=="csGARCH"){
    
    uspec<-rugarch::ugarchspec(
      variance.model=list(model="csGARCH", garchOrder=c(1,1)),
      mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
      distribution.model = distribution)
    
    if (is.null(out_of_sample)){
      
      for(i in 1:Num_assets){
        u_est[[i]]<-rugarch::ugarchfit(spec=uspec,data=db[,i])
        est_details[[i]]<-u_est[[i]]@fit$robust.matcoef
        u_loglik[[i]]<-u_est[[i]]@fit$LLH
        u_Inf_criteria[[i]]<-Inf_criteria(u_est[[i]])
        sd_m[,i]<-u_est[[i]]@fit$sigma
        db[,i] <- db[,i]-mean(db[,i], na.rm = TRUE)
      }
      
    } else {
      
      for(i in 1:Num_assets){
        u_est[[i]]<-rugarch::ugarchfit(spec=uspec,data=db[,i],out.sample=out_of_sample)
        est_details[[i]]<-u_est[[i]]@fit$robust.matcoef
        u_loglik[[i]]<-u_est[[i]]@fit$LLH
        u_Inf_criteria[[i]]<-Inf_criteria(u_est[[i]])
        sd_m[,i]<-u_est[[i]]@fit$sigma
        sd_m_oos[,i]<-as.numeric(rugarch::sigma(rugarch::ugarchforecast(u_est[[i]],  n.ahead = 1, 
                                                                        n.roll = c(out_of_sample-1), 
                                                                        out.sample = c(out_of_sample-1))))
        db_in_s[,i] <- db_in_s[,i]-mean(db[,i], na.rm = TRUE)
        db_oos[,i] <- db_oos[,i]-mean(db[,i], na.rm = TRUE)
      }
      
      
    }
    ##GM_noskew_RV/GM_noskew_X/GM_noskew_RV_X/GM_noskew_RV_sb/GM_noskew_X_sb/GM_noskew_RV_X_sb
  } else if (univ_model %in% c("GM_noskew_RV","GM_noskew_X","GM_noskew_RV_X","GM_noskew_RV_sb","GM_noskew_X_sb","GM_noskew_RV_X_sb")){ 
    
    if (is.null(out_of_sample)){
      
      for(i in 1:Num_assets){
        u_est[[i]]<-garchmidas_noskew_fit(model=univ_model,daily_ret=db[,i],RV_m=RV_m[[i]],mv_m=mv_m[[i]],K=K,
                                          distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec,
                                          out_of_sample=out_of_sample,vol_proxy=vol_proxy[,i],R=R)
        
        est_details[[i]]<-u_est[[i]]$rob_coef_mat
        u_loglik[[i]]<-u_est[[i]]$loglik
        u_Inf_criteria[[i]]<-Inf_criteria(u_est[[i]])
        sd_m[,i]<-u_est[[i]]$est_vol_in_s
        sd_long_m[,i] <- u_est[[i]]$est_lr_in_s
        sd_short_m[,i] <- u_est[[i]]$est_sr_in_s
        db[,i] <- db[,i] - u_est[[i]][["rob_coef_mat"]]["mu", "Estimate"]
        
      } 
      
    } else {
      
      for(i in 1:Num_assets){
        u_est[[i]]<-garchmidas_noskew_fit(model=univ_model,daily_ret=db[,i],RV_m=RV_m[[i]],mv_m=mv_m[[i]],K=K,
                                          distribution=distribution,lag_fun=lag_fun,X_stru_break_date_vec=X_stru_break_date_vec,
                                          out_of_sample=out_of_sample,vol_proxy=vol_proxy[,i],R=R)
        
        est_details[[i]]<-u_est[[i]]$rob_coef_mat
        u_loglik[[i]]<-u_est[[i]]$loglik
        u_Inf_criteria[[i]]<-Inf_criteria(u_est[[i]])
        sd_m[,i]<-u_est[[i]]$est_vol_in_s
        sd_long_m[,i] <- u_est[[i]]$est_lr_in_s
        sd_short_m[,i] <- u_est[[i]]$est_sr_in_s
        
        db_in_s[,i] <- db_in_s[,i]-u_est[[i]][["rob_coef_mat"]]["mu", "Estimate"]
        db_oos[,i] <- db_oos[,i]-u_est[[i]][["rob_coef_mat"]]["mu", "Estimate"]
        
        sd_m_oos[,i]<-zoo::coredata(u_est[[i]]$est_vol_oos)
        sd_long_m_oos[,i] <- u_est[[i]]$est_lr_oos
        sd_short_m_oos[,i] <- u_est[[i]]$est_sr_oos
      } 
    } 
  } 
  cat("First step: completed \n")
  
  ##### standardized residuals
  
  if (is.null(out_of_sample)){
    
    db_no_xts <- zoo::coredata(db)
    
    # 强制转为矩阵 + 逐元素除法
    eps_t <- as.matrix(db_no_xts) / as.matrix(sd_m)
    
    # 3. 如果后续代码必须要求三维数组格式，再进行转换
    # 但通常建议保持为 Matrix 以节省内存和提高后续速度
    #eps_matrix <- db_no_xts / sd_m
    #eps_t <- array(t(eps_matrix), dim = c(Num_assets, 1, TT))
    
    D_t <- array(0, dim = c(Num_assets, Num_assets, TT))
    for(tt in 1:TT) diag(D_t[,,tt]) <- sd_m[tt,]
    
    eps_t_oos <- NA
    D_t_oos <- NA
    
  } else {
    
    ## in-sample
    TT_ins <- TT - out_of_sample
    db_ins <- zoo::coredata(db_in_s)
    sd_ins <- sd_m[1:TT_ins, ]
    
    
    # 强制转为矩阵 + 逐元素除法
    eps_t <- as.matrix(db_ins) / as.matrix(sd_ins)
    
    D_t <- array(0, dim = c(Num_assets, Num_assets, TT_ins))
    for(tt in 1:TT_ins) diag(D_t[,,tt]) <- sd_m[tt,]
    
    ## out-of-sample
    
    
    db_oos_no_xts <- zoo::coredata(db_oos)
    
    eps_t_oos <- as.matrix(db_oos_no_xts) / as.matrix(sd_m_oos)
    
    D_t_oos <- array(0, dim = c(Num_assets, Num_assets, out_of_sample))
    for(tt in 1:out_of_sample) diag(D_t_oos[,,tt]) <- sd_m_oos[tt,]
    
    
  }
  
  
  ########################################### second step estimation (cDCC, aDCC, DECO or DCC-MIDAS)
  ##corr_model:cDCC/DECO/aDCC
  ##DCCMIDAS_RC/DCCMIDAS_RC_sb/DCCMIDAS_X/DCCMIDAS_X_sb/DCCMIDAS_RC_X/DCCMIDAS_RC_X_sb
  
  if (corr_model=="cDCC"){
    
    start_val<-ui<-ci<-NULL
    
    start_val<-c(a=0.01,b=0.8)
    
    ui<-rbind(
      c(1,0),       	 	 	## alpha>0.0001
      c(0,1),        		 	## beta>0.001
      c(-1,-1))	     		## alpha+beta<1
    
    ci<-c(-0.0001,-0.001,0.999)
    
    # numerical_gradient <- function(param,res,K_c=NULL) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(res)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- new_dcc_loglik(param = param + step, res = res,K_c = K_c)
    #     minus_obj <- new_dcc_loglik(param = param - step, res = res,K_c = K_c)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    
    m_est<-maxLik(logLik=new_dcc_loglik,
                  # grad = numerical_gradient,
                  start=start_val,
                  res=eps_t,
                  K_c=K_c,
                  constraints=list(ineqA=ui, ineqB=ci),
                  iterlim=1000,
                  method="BFGS") 
    
  } else if (corr_model=="DECO"){
    
    start_val<-ui<-ci<-NULL
    
    start_val<-c(a=0.01,b=0.8)
    
    ui<-rbind(
      c(1,0),       	 	 	## alpha>0.0001
      c(0,1),        		 	## beta>0.001
      c(-1,-1))	     		## alpha+beta<1
    
    ci<-c(-0.0001,-0.001,0.999)
    
    
    m_est<-maxLik(logLik=new_deco_loglik,
                  start=start_val,
                  res=eps_t,
                  K_c=K_c,
                  constraints=list(ineqA=ui, ineqB=ci),
                  iterlim=1000,
                  method="BFGS") 
    
  } else if (corr_model=="aDCC"){
    
    start_val<-ui<-ci<-NULL
    
    start_val<-c(a=0.01,b=0.8,g=0.01)
    
    eta <- eps_t
    eta[eps_t >= 0] <- 0
    
    sample_S<-stats::cov(eps_t)
    sample_N<-stats::cov(eta)
    
    inv_sample_S<-sample_S%^%(-0.5)
    inv_sample_N<-Inv(sample_N)
    
    delta<-eigen(inv_sample_S%*%inv_sample_N%*%inv_sample_S)$values[1]
    
    ui<-rbind(
      c(1,0,0),       	 	 	## a>0.0001
      c(0,1,0),        		 	## b>0.001
      c(0,0,1),				## g>0.001
      c(-1,-1,-delta),	     	## a+b+delta*g<1
      c(0,0,-1))				## g<0.15
    
    ci<-c(-0.0001,-0.001,-0.001,0.999,0.15)
    
    m_est<-maxLik(logLik=new_a_dcc_loglik,
                  start=start_val,
                  res=eps_t,
                  K_c=K_c,
                  constraints=list(ineqA=ui, ineqB=ci),
                  iterlim=1000,
                  method="BFGS") 
    
    
    
  } ##DCCMIDAS_RC/DCCMIDAS_RC_sb/DCCMIDAS_X/DCCMIDAS_X_sb/DCCMIDAS_RC_X/DCCMIDAS_RC_X_sb
  else if (corr_model=="DCCMIDAS_RC"&lag_fun=="Beta"){
    
    start_val<-ui<-ci<-NULL
    
    start_val<-c(a=0.01,b=0.8,w2=2)
    
    ui<-rbind(
      c(1,0,0),       	 	 	## alpha>0.0001
      c(0,1,0),        		 	## beta>0.001
      c(-1,-1,0),	     		## alpha+beta<1
      c(0,0,1))				## w2>1.001
    
    ci<-c(-0.0001,-0.001,0.999,-1.001)
    
    # numerical_gradient <- function(param,res,lag_fun,N_c,K_c) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(res)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    # 
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- dccmidas_RC_loglik(param = param + step, res = res,lag_fun=lag_fun,N_c=N_c, K_c = K_c)
    #     minus_obj <- dccmidas_RC_loglik(param = param - step, res = res,lag_fun=lag_fun, N_c=N_c,K_c = K_c)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    m_est<-maxLik(logLik=dccmidas_noskew_loglik,
                  # grad = numerical_gradient,
                  start=start_val,
                  model=corr_model,
                  res=eps_t,
                  
                  days_period=days_period,
                  mv_c=mv_c_in_s,
                  lag_fun=lag_fun,
                  N_c=N_c,
                  K_c=K_c,
                  X_stru_break_date_vec=X_stru_break_date_vec,
                  constraints=list(ineqA=ui, ineqB=ci),
                  iterlim=1000,
                  method="BFGS") 
    
  } else if (corr_model=="DCCMIDAS_RC"&lag_fun=="Almon") {
    
    start_val<-ui<-ci<-NULL
    
    start_val<-c(a=0.01,b=0.8,w2=-0.1)
    
    ui<-rbind(
      c(1,0,0),       	 	 	## alpha>0.0001
      c(0,1,0),        		 	## beta>0.001
      c(-1,-1,0),	     		## alpha+beta<1
      c(0,0,-1))					## w2<0
    
    ci<-c(-0.0001,-0.001,0.999,-0.001)
    
    # numerical_gradient <- function(param,res,lag_fun,N_c,K_c) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(res)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- dccmidas_RC_loglik(param = param + step, res = res,lag_fun=lag_fun,N_c=N_c, K_c = K_c)
    #     minus_obj <- dccmidas_RC_loglik(param = param - step, res = res,lag_fun=lag_fun, N_c=N_c,K_c = K_c)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    m_est<-maxLik(logLik=dccmidas_noskew_loglik,
                  # grad = numerical_gradient,
                  start=start_val,
                  model=corr_model,
                  res=eps_t,
                  
                  days_period=days_period,
                  mv_c=mv_c_in_s,
                  lag_fun=lag_fun,
                  N_c=N_c,
                  K_c=K_c,
                  X_stru_break_date_vec=X_stru_break_date_vec,
                  constraints=list(ineqA=ui, ineqB=ci),
                  iterlim=1000,
                  method="BFGS") 
    
  } else if (corr_model=="DCCMIDAS_X"&lag_fun=="Beta") {
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    
    param_delta_c_named <- paste0("delta_c_", seq_len(length_dccmidas_param))
    param_theta_c_named <- paste0("theta_c_", seq_len(length_dccmidas_param))
    param_w2_c_named <- paste0("w2_c_", seq_len(length_dccmidas_param))
    
    param_names <- c("a","b",param_delta_c_named,param_theta_c_named,param_w2_c_named)
    
    
    
    begin_val<-matrix(NA,nrow=R,ncol=(2+3*length_dccmidas_param))
    colnames(begin_val)<-param_names
    
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    for (i in 3:(2+2*length_dccmidas_param)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    begin_val[,(length_dccmidas_param*2+3):(length_dccmidas_param*3+2)]<-1.01
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(dccmidas_noskew_loglik(begin_val[i,],model=corr_model,res=eps_t,days_period=days_period,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c=K_c,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    names(start_val) <- param_names
    
    w2_ui1 <- diag(length_dccmidas_param)
    w2_ui2 <- matrix(0,length_dccmidas_param,(2*length_dccmidas_param+2))
    w2_ui <- cbind(w2_ui2,w2_ui1)
    
    ui<-rbind(
      c(1,0,rep(0,3*length_dccmidas_param)),       	 	 	## alpha>0.0001
      c(0,1,rep(0,3*length_dccmidas_param)),        		 	## beta>0.001
      c(-1,-1,rep(0,3*length_dccmidas_param)),	     		## alpha+beta<1
      w2_ui	)			                             	## ## all w2>1.001
    
    ci<-c(-0.0001,-0.001,0.999,rep(-1.001,length_dccmidas_param))
    
    # numerical_gradient <- function(param,res,mv_c,lag_fun,K_c) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(res)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- dccmidas_X_loglik(param = param + step, res = res,mv_c=mv_c_in_s,lag_fun=lag_fun,K_c = K_c)
    #     minus_obj <- dccmidas_X_loglik(param = param - step, res = res,mv_c=mv_c_in_s,lag_fun=lag_fun,K_c = K_c)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    m_est<-maxLik(logLik=dccmidas_noskew_loglik,
                  # grad = numerical_gradient,
                  start=start_val,
                  model=corr_model,
                  res=eps_t,
                  
                  days_period=days_period,
                  mv_c=mv_c_in_s,
                  lag_fun=lag_fun,
                  N_c=N_c,
                  K_c=K_c,
                  X_stru_break_date_vec=X_stru_break_date_vec,
                  constraints=list(ineqA=ui, ineqB=ci),
                  iterlim=1000,
                  method="BFGS") 
    
  } else if (corr_model=="DCCMIDAS_X"&lag_fun=="Almon") {
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    
    param_delta_c_named <- paste0("delta_c_", seq_len(length_dccmidas_param))
    param_theta_c_named <- paste0("theta_c_", seq_len(length_dccmidas_param))
    param_w2_c_named <- paste0("w2_c_", seq_len(length_dccmidas_param))
    
    param_names <- c("a","b",param_delta_c_named,param_theta_c_named,param_w2_c_named)
    
    
    
    begin_val<-matrix(NA,nrow=R,ncol=(2+3*length_dccmidas_param))
    colnames(begin_val)<-param_names
    
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    for (i in 3:(2+2*length_dccmidas_param)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    begin_val[,(length_dccmidas_param*2+3):(length_dccmidas_param*3+2)]<- -0.01
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(dccmidas_noskew_loglik(begin_val[i,],model=corr_model,res=eps_t,days_period=days_period,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c=K_c,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    names(start_val) <- param_names
    
    w2_ui1 <- diag(length_dccmidas_param)
    w2_ui2 <- matrix(0,length_dccmidas_param,(2*length_dccmidas_param+2))
    w2_ui <- cbind(w2_ui2,w2_ui1)
    
    ui<-rbind(
      c(1,0,rep(0,3*length_dccmidas_param)),       	 	 	## alpha>0.0001
      c(0,1,rep(0,3*length_dccmidas_param)),        		 	## beta>0.001
      c(-1,-1,rep(0,3*length_dccmidas_param)),	     		## alpha+beta<1
      -w2_ui	)			                             	## ## all w2<0
    
    ci<-c(-0.0001,-0.001,0.999,rep(-0.001,length_dccmidas_param))
    
    # numerical_gradient <- function(param,res,mv_c,lag_fun,K_c) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(res)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- dccmidas_X_loglik(param = param + step, res = res,mv_c=mv_c_in_s,lag_fun=lag_fun,K_c = K_c)
    #     minus_obj <- dccmidas_X_loglik(param = param - step, res = res,mv_c=mv_c_in_s,lag_fun=lag_fun,K_c = K_c)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    m_est<-maxLik(logLik=dccmidas_noskew_loglik,
                  # grad = numerical_gradient,
                  start=start_val,
                  model=corr_model,
                  res=eps_t,
                  
                  days_period=days_period,
                  mv_c=mv_c_in_s,
                  lag_fun=lag_fun,
                  N_c=N_c,
                  K_c=K_c,
                  X_stru_break_date_vec=X_stru_break_date_vec,
                  constraints=list(ineqA=ui, ineqB=ci),
                  iterlim=1000,
                  method="BFGS") 
    
  } else if (corr_model=="DCCMIDAS_RC_X"&lag_fun=="Beta") {
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    
    param_delta_c_named <- paste0("delta_c_", seq_len(length_dccmidas_param))
    param_theta_c_named <- paste0("theta_c_", seq_len(length_dccmidas_param))
    param_w2_c_named <- paste0("w2_c_", seq_len(length_dccmidas_param))
    param_theta_RC_c_named <- paste0("theta_RC_c_", seq_len(length_dccmidas_param))
    
    param_names <- c("a","b","RCw2",param_delta_c_named,param_theta_c_named,param_w2_c_named,param_theta_RC_c_named)
    
    
    
    begin_val<-matrix(NA,nrow=R,ncol=(3+4*length_dccmidas_param))
    colnames(begin_val)<-param_names
    
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-1.01
    for (i in 4:(3+2*length_dccmidas_param)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    begin_val[,(length_dccmidas_param*2+4):(length_dccmidas_param*3+3)]<-1.01
    begin_val[,(length_dccmidas_param*3+4):(length_dccmidas_param*4+3)]<-stats::runif(R,min=-1,max=1)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(dccmidas_noskew_loglik(begin_val[i,],model=corr_model,res=eps_t,days_period=days_period,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c=K_c,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    names(start_val) <- param_names
    
    w2_ui1 <- diag(length_dccmidas_param)
    w2_ui2 <- cbind(matrix(0,length_dccmidas_param,(2*length_dccmidas_param+3)),w2_ui1)
    w2_ui <- cbind(w2_ui2,matrix(0,length_dccmidas_param,length_dccmidas_param))
    
    ui<-rbind(
      c(1,0,0,rep(0,4*length_dccmidas_param)),       	 	 	## alpha>0.0001
      c(0,1,0,rep(0,4*length_dccmidas_param)),        		 	## beta>0.001
      c(-1,-1,0,rep(0,4*length_dccmidas_param)),	     		## alpha+beta<1
      c(0,0,1,rep(0,4*length_dccmidas_param)),             ##RCw2>1.001
      w2_ui	)			                             	## ## all w2>1.001
    
    ci<-c(-0.0001,-0.001,0.999,-1.001,rep(-1.001,length_dccmidas_param))
    
    # numerical_gradient <- function(param,res,mv_c,lag_fun,N_c,K_c) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(res)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- dccmidas_RC_X_loglik(param = param + step, res = res,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c = K_c)
    #     minus_obj <- dccmidas_RC_X_loglik(param = param - step, res = res,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c = K_c)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    m_est<-maxLik(logLik=dccmidas_noskew_loglik,
                  # grad = numerical_gradient,
                  start=start_val,
                  model=corr_model,
                  res=eps_t,
                  
                  days_period=days_period,
                  mv_c=mv_c_in_s,
                  lag_fun=lag_fun,
                  N_c=N_c,
                  K_c=K_c,
                  X_stru_break_date_vec=X_stru_break_date_vec,
                  constraints=list(ineqA=ui, ineqB=ci),
                  iterlim=1000,
                  method="BFGS") 
    
  } else if (corr_model=="DCCMIDAS_RC_X"&lag_fun=="Almon") {
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    
    param_delta_c_named <- paste0("delta_c_", seq_len(length_dccmidas_param))
    param_theta_c_named <- paste0("theta_c_", seq_len(length_dccmidas_param))
    param_w2_c_named <- paste0("w2_c_", seq_len(length_dccmidas_param))
    param_theta_RC_c_named <- paste0("theta_RC_c_", seq_len(length_dccmidas_param))
    
    param_names <- c("a","b","RCw2",param_delta_c_named,param_theta_c_named,param_w2_c_named,param_theta_RC_c_named)
    
    
    
    begin_val<-matrix(NA,nrow=R,ncol=(3+4*length_dccmidas_param))
    colnames(begin_val)<-param_names
    
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<--0.01
    for (i in 4:(3+2*length_dccmidas_param)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    begin_val[,(length_dccmidas_param*2+4):(length_dccmidas_param*3+3)]<--0.01
    begin_val[,(length_dccmidas_param*3+4):(length_dccmidas_param*4+3)]<-stats::runif(R,min=-1,max=1)
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(dccmidas_noskew_loglik(begin_val[i,],model=corr_model,res=eps_t,days_period=days_period,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c=K_c,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    names(start_val) <- param_names
    
    w2_ui1 <- diag(length_dccmidas_param)
    w2_ui2 <- cbind(matrix(0,length_dccmidas_param,(2*length_dccmidas_param+3)),w2_ui1)
    w2_ui <- cbind(w2_ui2,matrix(0,length_dccmidas_param,length_dccmidas_param))
    
    ui<-rbind(
      c(1,0,0,rep(0,4*length_dccmidas_param)),       	 	 	## alpha>0.0001
      c(0,1,0,rep(0,4*length_dccmidas_param)),        		 	## beta>0.001
      c(-1,-1,0,rep(0,4*length_dccmidas_param)),	     		## alpha+beta<1
      c(0,0,-1,rep(0,4*length_dccmidas_param)),             ##RCw2<0
      -w2_ui	)			                             	## ## all w2<0
    
    ci<-c(-0.0001,-0.001,0.999,-0.001,rep(-0.001,length_dccmidas_param))
    
    # numerical_gradient <- function(param,res,mv_c,lag_fun,N_c,K_c) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(res)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- dccmidas_RC_X_loglik(param = param + step, res = res,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c = K_c)
    #     minus_obj <- dccmidas_RC_X_loglik(param = param - step, res = res,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c = K_c)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    m_est<-maxLik(logLik=dccmidas_noskew_loglik,
                  # grad = numerical_gradient,
                  start=start_val,
                  model=corr_model,
                  res=eps_t,
                  
                  days_period=days_period,
                  mv_c=mv_c_in_s,
                  lag_fun=lag_fun,
                  N_c=N_c,
                  K_c=K_c,
                  X_stru_break_date_vec=X_stru_break_date_vec,
                  constraints=list(ineqA=ui, ineqB=ci),
                  iterlim=1000,
                  method="BFGS") 
    
  } #####structural breaks version
  else if (corr_model=="DCCMIDAS_RC_sb"&lag_fun=="Beta"){
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    n_n <- length(X_stru_break_date_vec)
    
    param_curr_m_vec_named <- paste0("m_c_", seq_len(n_n+1))
    param_curr_theta_RC_vec_named <- paste0("theta_RC_c_", seq_len(n_n+1))
    
    param_names <- c("a","b","RCw2",param_curr_m_vec_named,param_curr_theta_RC_vec_named)
    
    
    
    begin_val<-matrix(NA,nrow=R,ncol=(5+2*n_n))
    colnames(begin_val)<-param_names
    
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<-1.2
    for (i in 4:(5+2*n_n)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(dccmidas_noskew_loglik(begin_val[i,],model=corr_model,res=eps_t,days_period=days_period,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c=K_c,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    names(start_val) <- param_names
    
    ui<-rbind(
      c(1,0,0,rep(0,2+2*n_n)),       	 	 	## alpha>0.0001
      c(0,1,0,rep(0,2+2*n_n)),        		 	## beta>0.001
      c(-1,-1,0,rep(0,2+2*n_n)),	     		## alpha+beta<1
      c(0,0,1,rep(0,2+2*n_n)))				## w2>1.001
    
    ci<-c(-0.0001,-0.001,0.999,-1.001)
    
    # numerical_gradient <- function(param,res,lag_fun,N_c,K_c,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(res)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- dccmidas_RC_sb_loglik(param = param + step, res = res,lag_fun=lag_fun,N_c=N_c, K_c = K_c,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- dccmidas_RC_sb_loglik(param = param - step, res = res,lag_fun=lag_fun, N_c=N_c,K_c = K_c,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    m_est<-maxLik(logLik=dccmidas_noskew_loglik,
                  # grad = numerical_gradient,
                  start=start_val,
                  model=corr_model,
                  res=eps_t,
                  
                  days_period=days_period,
                  mv_c=mv_c_in_s,
                  lag_fun=lag_fun,
                  N_c=N_c,
                  K_c=K_c,
                  X_stru_break_date_vec=X_stru_break_date_vec,
                  constraints=list(ineqA=ui, ineqB=ci),
                  iterlim=1000,
                  method="BFGS") 
    
  } else if (corr_model=="DCCMIDAS_RC_sb"&lag_fun=="Almon") {
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    n_n <- length(X_stru_break_date_vec)
    
    param_curr_m_vec_named <- paste0("m_c_", seq_len(n_n+1))
    param_curr_theta_RC_vec_named <- paste0("theta_RC_c_", seq_len(n_n+1))
    
    param_names <- c("a","b","RCw2",param_curr_m_vec_named,param_curr_theta_RC_vec_named)
    
    
    
    begin_val<-matrix(NA,nrow=R,ncol=(5+2*n_n))
    colnames(begin_val)<-param_names
    
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<--0.01
    for (i in 4:(5+2*n_n)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(dccmidas_noskew_loglik(begin_val[i,],model=corr_model,res=eps_t,days_period=days_period,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c=K_c,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    names(start_val) <- param_names
    
    ui<-rbind(
      c(1,0,0,rep(0,2+2*n_n)),       	 	 	## alpha>0.0001
      c(0,1,0,rep(0,2+2*n_n)),        		 	## beta>0.001
      c(-1,-1,0,rep(0,2+2*n_n)),	     		## alpha+beta<1
      c(0,0,-1,rep(0,2+2*n_n)))				## w2<0
    
    ci<-c(-0.0001,-0.001,0.999,-0.001)
    
    # numerical_gradient <- function(param,res,lag_fun,N_c,K_c,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(res)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- dccmidas_RC_sb_loglik(param = param + step, res = res,lag_fun=lag_fun,N_c=N_c, K_c = K_c,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- dccmidas_RC_sb_loglik(param = param - step, res = res,lag_fun=lag_fun, N_c=N_c,K_c = K_c,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    m_est<-maxLik(logLik=dccmidas_noskew_loglik,
                  # grad = numerical_gradient,
                  start=start_val,
                  model=corr_model,
                  res=eps_t,
                  
                  days_period=days_period,
                  mv_c=mv_c_in_s,
                  lag_fun=lag_fun,
                  N_c=N_c,
                  K_c=K_c,
                  X_stru_break_date_vec=X_stru_break_date_vec,
                  constraints=list(ineqA=ui, ineqB=ci),
                  iterlim=1000,
                  method="BFGS") 
    
  } else if (corr_model=="DCCMIDAS_X_sb"&lag_fun=="Beta") {
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    n_n <- length(X_stru_break_date_vec)
    
    
    param_w2_c_named <- paste0("w2_c_", seq_len(length_dccmidas_param))
    
    param_curr_m_vec_named <- paste0("m_c_", seq_len(n_n+1))
    param_curr_theta_vec_named <- paste0("theta_c_", seq_len(n_n+1))
    
    param_names <- c("a","b",param_w2_c_named,param_curr_m_vec_named,param_curr_theta_vec_named)
    
    
    
    begin_val<-matrix(NA,nrow=R,ncol=(4+length_dccmidas_param+2*n_n))
    colnames(begin_val)<-param_names
    
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3:(length_dccmidas_param+2)]<- 1.01
    for (i in (length_dccmidas_param+3):(4+length_dccmidas_param+2*n_n)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(dccmidas_noskew_loglik(begin_val[i,],model=corr_model,res=eps_t,days_period=days_period,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c=K_c,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    names(start_val) <- param_names
    
    w2_ui1 <- diag(length_dccmidas_param)
    w2_ui2 <- matrix(0,length_dccmidas_param,2)
    w2_ui3 <- cbind(w2_ui2,w2_ui1)
    w2_ui4 <- matrix(0,length_dccmidas_param,(2+2*n_n))
    w2_ui <- cbind(w2_ui3,w2_ui4)
    
    ui<-rbind(
      c(1,0,rep(0,(2+length_dccmidas_param+2*n_n))),       	 	 	## alpha>0.0001
      c(0,1,rep(0,(2+length_dccmidas_param+2*n_n))),        		 	## beta>0.001
      c(-1,-1,rep(0,(2+length_dccmidas_param+2*n_n))),	     		## alpha+beta<1
      w2_ui	)			                             	## ## all w2>1.001
    
    ci<-c(-0.0001,-0.001,0.999,rep(-1.001,length_dccmidas_param))
    
    # numerical_gradient <- function(param,res,mv_c,lag_fun,K_c,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(res)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- dccmidas_X_sb_loglik(param = param + step, res = res,mv_c=mv_c_in_s,lag_fun=lag_fun,K_c = K_c,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- dccmidas_X_sb_loglik(param = param - step, res = res,mv_c=mv_c_in_s,lag_fun=lag_fun,K_c = K_c,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    m_est<-maxLik(logLik=dccmidas_noskew_loglik,
                  # grad = numerical_gradient,
                  start=start_val,
                  model=corr_model,
                  res=eps_t,
                  
                  days_period=days_period,
                  mv_c=mv_c_in_s,
                  lag_fun=lag_fun,
                  N_c=N_c,
                  K_c=K_c,
                  X_stru_break_date_vec=X_stru_break_date_vec,
                  constraints=list(ineqA=ui, ineqB=ci),
                  iterlim=1000,
                  method="BFGS") 
    
  } else if (corr_model=="DCCMIDAS_X_sb"&lag_fun=="Almon") {
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    n_n <- length(X_stru_break_date_vec)
    
    
    param_w2_c_named <- paste0("w2_c_", seq_len(length_dccmidas_param))
    
    param_curr_m_vec_named <- paste0("m_c_", seq_len(n_n+1))
    param_curr_theta_vec_named <- paste0("theta_c_", seq_len(n_n+1))
    
    param_names <- c("a","b",param_w2_c_named,param_curr_m_vec_named,param_curr_theta_vec_named)
    
    
    
    begin_val<-matrix(NA,nrow=R,ncol=(4+length_dccmidas_param+2*n_n))
    colnames(begin_val)<-param_names
    
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3:(length_dccmidas_param+2)]<- -0.01
    for (i in (length_dccmidas_param+3):(4+length_dccmidas_param+2*n_n)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(dccmidas_noskew_loglik(begin_val[i,],model=corr_model,res=eps_t,days_period=days_period,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c=K_c,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    names(start_val) <- param_names
    
    w2_ui1 <- diag(length_dccmidas_param)
    w2_ui2 <- matrix(0,length_dccmidas_param,2)
    w2_ui3 <- cbind(w2_ui2,w2_ui1)
    w2_ui4 <- matrix(0,length_dccmidas_param,(2+2*n_n))
    w2_ui <- cbind(w2_ui3,w2_ui4)
    
    ui<-rbind(
      c(1,0,rep(0,(2+length_dccmidas_param+2*n_n))),       	 	 	## alpha>0.0001
      c(0,1,rep(0,(2+length_dccmidas_param+2*n_n))),        		 	## beta>0.001
      c(-1,-1,rep(0,(2+length_dccmidas_param+2*n_n))),	     		## alpha+beta<1
      -w2_ui	)			                             	## ## all w2<0
    
    ci<-c(-0.0001,-0.001,0.999,rep(-0.001,length_dccmidas_param))
    
    # numerical_gradient <- function(param,res,mv_c,lag_fun,K_c,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(res)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- dccmidas_X_sb_loglik(param = param + step, res = res,mv_c=mv_c_in_s,lag_fun=lag_fun,K_c = K_c,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- dccmidas_X_sb_loglik(param = param - step, res = res,mv_c=mv_c_in_s,lag_fun=lag_fun,K_c = K_c,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    m_est<-maxLik(logLik=dccmidas_noskew_loglik,
                  # grad = numerical_gradient,
                  start=start_val,
                  model=corr_model,
                  res=eps_t,
                  days_period=days_period,
                  mv_c=mv_c_in_s,
                  lag_fun=lag_fun,
                  N_c=N_c,
                  K_c=K_c,
                  X_stru_break_date_vec=X_stru_break_date_vec,
                  constraints=list(ineqA=ui, ineqB=ci),
                  iterlim=1000,
                  method="BFGS") 
    
  } else if (corr_model=="DCCMIDAS_RC_X_sb"&lag_fun=="Beta") {
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    n_n <- length(X_stru_break_date_vec)
    
    
    param_w2_c_named <- paste0("w2_c_", seq_len(length_dccmidas_param))
    param_m_c_named <- paste0("m_c_", seq_len(n_n+1))
    param_theta_RC_c_named <- paste0("theta_RC_c_", seq_len(n_n+1))
    param_theta_c_named <- paste0("theta_c_", seq_len(n_n+1))
    
    param_names <- c("a","b","RCw2",param_w2_c_named,param_m_c_named,param_theta_RC_c_named,param_theta_c_named)
    
    
    
    begin_val<-matrix(NA,nrow=R,ncol=(6+length_dccmidas_param+3*n_n))
    colnames(begin_val)<-param_names
    
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<- 1.01
    begin_val[,4:(length_dccmidas_param+3)]<- 1.01
    for (i in (length_dccmidas_param+4):(6+length_dccmidas_param+3*n_n)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(dccmidas_noskew_loglik(begin_val[i,],model=corr_model,res=eps_t,days_period=days_period,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c=K_c,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    names(start_val) <- param_names
    
    w2_ui1 <- diag(length_dccmidas_param)
    w2_ui2 <- cbind(matrix(0,length_dccmidas_param,3),w2_ui1)
    w2_ui <- cbind(w2_ui2,matrix(0,length_dccmidas_param,3*n_n+3))
    
    ui<-rbind(
      c(1,0,0,rep(0,3+length_dccmidas_param+3*n_n)),       	 	 	## alpha>0.0001
      c(0,1,0,rep(0,3+length_dccmidas_param+3*n_n)),        		 	## beta>0.001
      c(-1,-1,0,rep(0,3+length_dccmidas_param+3*n_n)),	     		## alpha+beta<1
      c(0,0,1,rep(0,3+length_dccmidas_param+3*n_n)),             ##RCw2>1.001
      w2_ui	)			                             	## ## all w2>1.001
    
    ci<-c(-0.0001,-0.001,0.999,-1.001,rep(-1.001,length_dccmidas_param))
    
    # numerical_gradient <- function(param,res,mv_c,lag_fun,N_c,K_c,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(res)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- dccmidas_noskew_loglik(param = param + step, res = res,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c = K_c,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- dccmidas_noskew_loglik(param = param - step, res = res,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c = K_c,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    m_est<-maxLik(logLik=dccmidas_noskew_loglik,
                  # grad = numerical_gradient,
                  start=start_val,
                  model=corr_model,
                  res=eps_t,
                  
                  days_period=days_period,
                  mv_c=mv_c_in_s,
                  lag_fun=lag_fun,
                  N_c=N_c,
                  K_c=K_c,
                  X_stru_break_date_vec=X_stru_break_date_vec,
                  constraints=list(ineqA=ui, ineqB=ci),
                  iterlim=1000,
                  method="BFGS") 
    
  } else if (corr_model=="DCCMIDAS_RC_X_sb"&lag_fun=="Almon") {
    
    start_val<-begin_val<-ui<-ci<-NULL
    
    length_dccmidas_param <- 0.5*Num_assets*(Num_assets-1)
    n_n <- length(X_stru_break_date_vec)
    
    
    param_w2_c_named <- paste0("w2_c_", seq_len(length_dccmidas_param))
    param_m_c_named <- paste0("m_c_", seq_len(n_n+1))
    param_theta_RC_c_named <- paste0("theta_RC_c_", seq_len(n_n+1))
    param_theta_c_named <- paste0("theta_c_", seq_len(n_n+1))
    
    param_names <- c("a","b","RCw2",param_w2_c_named,param_m_c_named,param_theta_RC_c_named,param_theta_c_named)
    
    
    
    begin_val<-matrix(NA,nrow=R,ncol=(6+length_dccmidas_param+3*n_n))
    colnames(begin_val)<-param_names
    
    begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
    begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
    begin_val[,3]<- -0.01
    begin_val[,4:(length_dccmidas_param+3)]<- -0.01
    for (i in (length_dccmidas_param+4):(6+length_dccmidas_param+3*n_n)) {
      begin_val[,i]<-stats::runif(R,min=-1,max=1)
    }
    
    which_row<-rep(NA,R)
    
    for(i in 1:R){
      which_row[i]<-sum(dccmidas_noskew_loglik(begin_val[i,],model=corr_model,res=eps_t,days_period=days_period,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c=K_c,X_stru_break_date_vec=X_stru_break_date_vec))
    }
    
    start_val<-begin_val[which.max(which_row),]
    names(start_val) <- param_names
    
    w2_ui1 <- diag(length_dccmidas_param)
    w2_ui2 <- cbind(matrix(0,length_dccmidas_param,3),w2_ui1)
    w2_ui <- cbind(w2_ui2,matrix(0,length_dccmidas_param,3*n_n+3))
    
    ui<-rbind(
      c(1,0,0,rep(0,3+length_dccmidas_param+3*n_n)),       	 	 	## alpha>0.0001
      c(0,1,0,rep(0,3+length_dccmidas_param+3*n_n)),        		 	## beta>0.001
      c(-1,-1,0,rep(0,3+length_dccmidas_param+3*n_n)),	     		## alpha+beta<1
      c(0,0,-1,rep(0,3+length_dccmidas_param+3*n_n)),             ##RCw2<0
      -w2_ui	)			                             	## ## all w2<0
    
    ci<-c(-0.0001,-0.001,0.999,-0.001,rep(-0.001,length_dccmidas_param))
    
    # numerical_gradient <- function(param,res,mv_c,lag_fun,N_c,K_c,X_stru_break_date_vec) {
    #   
    #   eps <- 1e-5
    #   K <- length(param)
    #   TT <- nrow(res)
    #   
    #   jac_matrix <- matrix(0, TT, K)
    #   
    #   for(i in 1:K){
    #     
    #     step <- rep(0, K)
    #     step[i] <- eps
    #     
    #     plus_obj <- dccmidas_noskew_loglik(param = param + step, res = res,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c = K_c,X_stru_break_date_vec=X_stru_break_date_vec)
    #     minus_obj <- dccmidas_noskew_loglik(param = param - step, res = res,mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c = K_c,X_stru_break_date_vec=X_stru_break_date_vec)
    #     
    #     jac_matrix[, i] <- (plus_obj - minus_obj) / (2.0 * eps)
    #   }
    #   return(jac_matrix)
    # }
    
    m_est<-maxLik(logLik=dccmidas_noskew_loglik,
                  # grad = numerical_gradient,
                  start=start_val,
                  model=corr_model,
                  res=eps_t,
                  
                  days_period = days_period,
                  mv_c=mv_c_in_s,
                  lag_fun=lag_fun,
                  N_c = N_c,
                  K_c = K_c,
                  X_stru_break_date_vec=X_stru_break_date_vec,
                  constraints=list(ineqA=ui, ineqB=ci),
                  iterlim=1000,
                  method="BFGS") 
    
  }
  ########################################### end second step
  
  ###### matrix of coefficients (second step)
  
  est_coef<-stats::coef(m_est)
  N_coef<-length(est_coef)
  
  mat_coef<-data.frame(rep(NA,N_coef),rep(NA,N_coef),rep(NA,N_coef),rep(NA,N_coef))
  colnames(mat_coef)<-c("Estimate","Std. Error","t value","Pr(>|t|)")
  
  rownames(mat_coef)<-names(est_coef)
  
  mat_coef[,1]<-round(est_coef,6)
  mat_coef[,2]<-round(QMLE_sd(m_est),6)
  mat_coef[,3]<-round(est_coef/QMLE_sd(m_est),6)
  mat_coef[,4]<-round(apply(rbind(est_coef/QMLE_sd(m_est)),1,function(x) 2*(1-stats::pnorm(abs(x)))),6)
  
  if(corr_model %in% c("DCCMIDAS_RC","DCCMIDAS_RC_sb","DCCMIDAS_X","DCCMIDAS_X_sb","DCCMIDAS_RC_X","DCCMIDAS_RC_X_sb")){
    
    dcc_mat_est_fin<-dccmidas_noskew_mat_est(est_coef,model=corr_model,res=eps_t,days_period=days_period,D_t=D_t,
                                             mv_c=mv_c_in_s,lag_fun=lag_fun,N_c=N_c,K_c=K_c,X_stru_break_date_vec=X_stru_break_date_vec)
    
    if (!is.null(out_of_sample)){
      
      if(corr_model %in% c("DCCMIDAS_RC","DCCMIDAS_RC_sb","DCCMIDAS_RC_X","DCCMIDAS_RC_X_sb")){
        
        N_c_res_for_oos <- eps_t[(nrow(eps_t)-N_c+1):nrow(eps_t),]
        
        dcc_mat_est_fin_oos<-dccmidas_noskew_mat_est(est_coef,model=corr_model,res=eps_t_oos,days_period=days_period_oos,D_t=D_t_oos,
                                                     mv_c=mv_c_oos,lag_fun=lag_fun,N_c=N_c,K_c=K_c,X_stru_break_date_vec=X_stru_break_date_vec,N_c_res_for_oos=N_c_res_for_oos)
        
      } else {
        
        dcc_mat_est_fin_oos<-dccmidas_noskew_mat_est(est_coef,model=corr_model,res=eps_t_oos,days_period=days_period_oos,D_t=D_t_oos,
                                                     mv_c=mv_c_oos,lag_fun=lag_fun,N_c=N_c,K_c=K_c,X_stru_break_date_vec=X_stru_break_date_vec)
        
      }
    } else {
      
      dcc_mat_est_fin_oos<-list(NA,NA,NA)
      
    }
    
    
  }  else if (corr_model=="cDCC") {
    
    dcc_mat_est_fin<-new_dcc_mat_est(est_coef,eps_t,D_t,K_c=K_c)
    
    
    if (!is.null(out_of_sample)){
      
      dcc_mat_est_fin_oos<-new_dcc_mat_est(est_coef,
                                           eps_t_oos,D_t_oos,K_c=K_c)
      
    } else {
      
      dcc_mat_est_fin_oos<-list(NA,NA)
      
    }
    
    
  } else if (corr_model=="aDCC") {
    
    dcc_mat_est_fin<-new_a_dcc_mat_est(est_coef,eps_t,D_t,K_c=K_c)
    
    
    if (!is.null(out_of_sample)){
      
      dcc_mat_est_fin_oos<-new_a_dcc_mat_est(est_coef,
                                             eps_t_oos,D_t_oos,K_c=K_c)
      
    } else {
      
      dcc_mat_est_fin_oos<-list(NA,NA)
      
    }
    
    
    
  } else if (corr_model=="DECO") {
    
    dcc_mat_est_fin<-new_deco_mat_est(est_coef,eps_t,D_t,K_c=K_c)
    
    
    if (!is.null(out_of_sample)){
      
      dcc_mat_est_fin_oos<-new_deco_mat_est(est_coef,
                                            eps_t_oos,D_t_oos,K_c=K_c)
      
    } else {
      
      dcc_mat_est_fin_oos<-list(NA,NA)
      
    }
    
  }
  
  ######
  
  end<-Sys.time()-start
  
  
  if(corr_model %in% c("DCCMIDAS_RC","DCCMIDAS_RC_sb","DCCMIDAS_X","DCCMIDAS_X_sb","DCCMIDAS_RC_X","DCCMIDAS_RC_X_sb")){
    fin_res<-list(
      assets=colnames(db),
      model=univ_model,
      est_univ_model=est_details,
      est_univ_model_origin=u_est,
      corr_coef_mat=mat_coef,
      mult_model=corr_model,
      obs=est_obs,
      period=estim_period,
      eps_t=eps_t,
      D_t=D_t,
      eps_t_oos=eps_t_oos,
      D_t_oos=D_t_oos,
      "H_t"=dcc_mat_est_fin[[1]],
      "R_t"=dcc_mat_est_fin[[2]],
      "R_t_bar"=dcc_mat_est_fin[[3]],
      "H_t_oos"=dcc_mat_est_fin_oos[[1]],
      "R_t_oos"=dcc_mat_est_fin_oos[[2]],
      "R_t_bar_oos"=dcc_mat_est_fin_oos[[3]],
      est_time=end,
      Days=days_period,
      Days_oos=days_period_oos,
      u_llk=u_loglik,
      c_llk=stats::logLik(m_est),
      u_Inf_criteria=u_Inf_criteria,
      c_inf_criteria=Inf_criteria(m_est),
      Inf_criteria_2step=Inf_criteria_2step(u_est,m_est)
    )
  } else {
    fin_res<-list(
      assets=colnames(db),
      model=univ_model,
      est_univ_model=est_details,
      est_univ_model_origin=u_est,
      corr_coef_mat=mat_coef,
      mult_model=corr_model,
      obs=est_obs,
      period=estim_period,
      eps_t=eps_t,
      D_t=D_t,
      eps_t_oos=eps_t_oos,
      D_t_oos=D_t_oos,
      "H_t"=dcc_mat_est_fin[[1]],
      "R_t"=dcc_mat_est_fin[[2]],
      "H_t_oos"=dcc_mat_est_fin_oos[[1]],
      "R_t_oos"=dcc_mat_est_fin_oos[[2]],
      est_time=end,
      Days=days_period,
      Days_oos=days_period_oos,
      u_llk=u_loglik,
      c_llk=stats::logLik(m_est),
      u_Inf_criteria=u_Inf_criteria,
      c_inf_criteria=Inf_criteria(m_est),
      Inf_criteria_2step=Inf_criteria_2step(u_est,m_est)
    )
  }
  
  class(fin_res)<-c("dmnoskew")
  return(fin_res)
  
}

summary.dmnoskew <- function(object, ...) {
  
  model<-object$model
  mult_model<-object$mult_model
  
  mat_coef<-object$corr_coef_mat
  Obs<-object$obs
  Period<-object$period
  
  Period<-paste(substr(Period[1], 1, 10),"/",
                substr(Period[2], 1, 10),sep="")
  
  ################################ multivariate matrix
  
  p_value<-mat_coef[,4]
  
  sig<-ifelse(p_value<=0.01,"***",ifelse(p_value>0.01&p_value<=0.05,"**",
                                         ifelse(p_value>0.05&p_value<=0.1,"*"," ")))
  
  mat_coef<-round(mat_coef,4)
  
  mat_coef<-cbind(mat_coef,Sig.=sig)
  
  ################################ univariate matrix
  
  assets<-object$assets
  
  est_univ<-object$est_univ_model
  
  mat_f<-list()
  for(i in 1:length(assets)){
    mat<-est_univ[[i]]
    p_value<-mat[,4]
    sig<-ifelse(p_value<=0.01,"***",ifelse(p_value>0.01&p_value<=0.05,"**",
                                           ifelse(p_value>0.05&p_value<=0.1,"*"," ")))
    mat<-round(mat,4)
    mat<-data.frame(mat,Sig.=sig)
    colnames(mat)[1:4]<-colnames(est_univ[[i]])
    mat_f[[i]]<-mat
  }
  names(mat_f) <- assets
  
  ################################# information criteria
  
  u_llk <- object$u_llk
  c_llk<- object$c_llk
  u_Inf_criteria<- object$u_Inf_criteria
  c_inf_criteria<- object$c_inf_criteria
  Inf_criteria_2step<- object$Inf_criteria_2step
  
  ################################# final print
  
  cat(
    cat("\n"),
    cat(paste("Univariate model:",model),"\n"),
    cat("\n"),
    cat("Est. coefficients of the univariate models:\n"),
    cat("\n"),
    cat(utils::capture.output(mat_f),  sep = '\n'),
    cat("---------------------------------------------------","\n"),
    cat(paste("Correlation model:",mult_model),"\n"),
    cat("\n"),
    cat("Est. coefficients of the correlation model:\n"),
    cat("\n"),
    cat(utils::capture.output(mat_coef),  sep = '\n'),
    cat("--- \n"),
    cat("Signif. codes: 0.01 '***', 0.05 '**', 0.1 '*' \n"),
    cat("\n"),
    cat("Obs.:", paste(Obs, ".",sep=""), "Sample Period:", Period, "\n"),
    cat("\n"),
    cat(utils::capture.output(u_llk),  sep = '\n'),
    cat("c_LogLik:", paste(c_llk, ".",sep=""),  sep = '\n'),
    cat(utils::capture.output(u_Inf_criteria),  sep = '\n'),
    cat("c_inf_criteria:", paste(c_inf_criteria, ".",sep=""),  sep = '\n'),
    cat(utils::capture.output(Inf_criteria_2step),  sep = '\n'),
    cat("\n"))
  
}

print.dmnoskew <- function(x, ...) {
  
  options(scipen = 999)
  
  u_model<-x[[2]]
  c_model<-x[[5]]
  mat_coef<-x[[4]]
  
  coef_char<-as.character(round(mat_coef[,1],4))
  
  row_names<-gsub("\\s", " ", format(rownames(mat_coef), width=9))
  coef_val<-gsub("\\s", " ", format(coef_char,width=9))
  
  cat(
    cat("\n"),
    cat(paste("Model:",u_model,c_model),"\n"),
    cat("\n"),
    cat(paste("c_Coefficients: \n",sep="\n")),
    cat(row_names, sep=" ", "\n"),
    cat(coef_val, sep=" ", "\n"),
    cat("\n"))
}

integrated_dccmidas_noskew_longfocus_output_function <- function(d_m,univ_model="sGARCH",distribution="norm",
                                                                 X_stru_break_date_vec=NULL,
                                                                 macro=NULL,type=NULL,K=NULL,
                                                                 corr_model="cDCC",lag_fun=NULL,
                                                                 N_c=NULL,K_c=NULL,out_of_sample=NULL,vol_proxy=NULL,R=NULL){
  
  ####available for monthly/quarterly
  
  if(!is.null(vol_proxy)&&class(vol_proxy)[1] != "xts") { stop(
    cat("#Warning:\n Parameter 'vol_proxy' must be an xts object. Please provide it in the correct form \n")
  )}
  
  
  if(class(d_m)[1] != "xts") { stop(
    cat("#Warning:\n Parameter 'd_m' must be an xts object. Please provide it in the correct form \n")
  )}
  
  if(nrow(d_m) < 50) { stop(
    cat("#Warning:\n Parameter 'd_m' must be long enough. Please provide longer one \n")
  )}
  
  if(!is.null(N_c)&& nrow(d_m) < N_c) { stop(
    cat("#Warning:\n Parameter 'd_m' must be long enough. Please provide longer one \n")
  )}
  
  if(!is.null(macro)&&class(macro)[1] != "xts") { stop(
    cat("#Warning:\n Parameter 'macro' must be an xts object. Please provide it in the correct form \n")
  )}
  
  if(!is.null(vol_proxy)&&class(vol_proxy)[1] != "xts") { stop(
    cat("#Warning:\n Parameter 'vol_proxy' must be an xts object. Please provide it in the correct form \n")
  )}
  
  if(!is.null(vol_proxy)&& ncol(vol_proxy) != ncol(d_m)) { stop(
    cat("#Warning:\n Parameter 'vol_proxy' must have as same columns as d_m. \n")
  )}
  
  if (univ_model%in% c("GM_noskew_RV","GM_noskew_RV_sb","GM_noskew_RV_X","GM_noskew_RV_X_sb","GM_noskew_X","GM_noskew_X_sb")|corr_model %in% c("DCCMIDAS_X","DCCMIDAS_X_sb","DCCMIDAS_RC_X","DCCMIDAS_RC_X_sb")){
    
    RV_macro <- list()
    for (i in 1:dim(d_m)[2]) {
      RV_macro[[i]] <- RV_gen(x=d_m[,i],type=type)
      
    }
    
    RV_required_d_m <- cut_func(r_t=d_m,K=K,K_c=K_c,macro=RV_macro[[1]],type=type)
    
    if(univ_model=="GM_noskew_RV_X"|univ_model=="GM_noskew_RV_X_sb"){
      
      X_required_d_m <- cut_func(r_t=d_m,K=K,K_c=K_c,macro=macro,type=type)
      common_dates <- as.Date(intersect(index(RV_required_d_m), index(X_required_d_m)))
      required_d_m <- X_required_d_m[common_dates]
      
    } else if ( univ_model=="GM_noskew_RV"|univ_model=="GM_noskew_RV_sb" ){
      
      if (corr_model %in% c("DCCMIDAS_X","DCCMIDAS_X_sb","DCCMIDAS_RC_X","DCCMIDAS_RC_X_sb")) {
        
        X_required_d_m <- cut_func(r_t=d_m,K=K,K_c=K_c,macro=macro,type=type)
        common_dates <- as.Date(intersect(index(RV_required_d_m), index(X_required_d_m)))
        required_d_m <- X_required_d_m[common_dates]
        
      } else {
        
        required_d_m <- RV_required_d_m
      }
      
      
      
    } else if (univ_model=="GM_noskew_X"|univ_model=="GM_noskew_X_sb"){
      
      X_required_d_m <- cut_func(r_t=d_m,K=K,K_c=K_c,macro=macro,type=type)
      required_d_m <- X_required_d_m
      
    } else {
      
      if (corr_model%in% c("DCCMIDAS_X", "DCCMIDAS_X_sb","DCCMIDAS_RV_X","DCCMIDAS_RV_X_sb")){
        
        X_required_d_m <- cut_func(r_t=d_m,K=K,K_c=K_c,macro=macro,type=type)
        required_d_m <- X_required_d_m
        
      } else {
        
        required_d_m <- d_m
        
      }
      
      
      
    }
    
    time_index <- index(required_d_m)
    
    if(!is.null(vol_proxy)){
      
      ref_dates <- time_index
      
      if (!all(ref_dates %in% index(vol_proxy))) {
        
        stop(paste0("vol_proxy data do not cover daily_ret"))
        
      } else {
        
        vol_proxy <- vol_proxy[ref_dates]
      }        
    }
    
    
    
    if (univ_model%in% c("GM_noskew_RV_X","GM_noskew_RV_X_sb","GM_noskew_X","GM_noskew_X_sb")|corr_model %in% c("DCCMIDAS_X","DCCMIDAS_X_sb","DCCMIDAS_RC_X","DCCMIDAS_RC_X_sb")){
      
      if (univ_model%in% c("GM_noskew_RV_X","GM_noskew_RV_X_sb")) {
        
        RV_m <- list()
        mv_m<-list()
        for (i in 1:dim(d_m)[2]) {
          
          mv_m[[i]]<- mv_into_mat2(x=required_d_m[,i],mv=macro,K=K,type=type)
          RV_m[[i]]<- mv_into_mat2(x=required_d_m[,i],mv=RV_macro[[i]],K=K,type=type)
        }
        
      } else if (univ_model%in% c("GM_noskew_X","GM_noskew_X_sb")) {
        
        
        mv_m<-list()
        for (i in 1:dim(d_m)[2]) {
          
          mv_m[[i]]<- mv_into_mat2(x=required_d_m[,i],mv=macro,K=K,type=type)
          
        }
        RV_m <- NULL
        
      } else if (univ_model%in% c("GM_noskew_RV","GM_noskew_RV_sb")) {
        
        
        RV_m <- list()
        for (i in 1:dim(d_m)[2]) {
          
          
          RV_m[[i]]<- mv_into_mat2(x=required_d_m[,i],mv=RV_macro[[i]],K=K,type=type)
        }
        
        mv_m<-NULL
        
      } else {
        
        RV_m <- NULL
        mv_m<-NULL
      }
      
      
      if (corr_model %in% c("DCCMIDAS_X","DCCMIDAS_X_sb","DCCMIDAS_RC_X","DCCMIDAS_RC_X_sb")) {
        
        mv_c <- mv_into_mat2(x=required_d_m[,1],mv=macro,K=K_c,type=type)
        
      } else {
        
        mv_c <- NULL
      }
      
      
      
      
      
    } else {
      
      RV_m <- list()
      for (i in 1:dim(d_m)[2]) {
        
        RV_m[[i]]<- mv_into_mat2(x=required_d_m[,i],mv=RV_macro[[i]],K=K,type=type)
        
      }
      mv_m <- NULL
      mv_c <- NULL
      
    }
  } else {
    
    time_index <- index(d_m)
    
    if(!is.null(vol_proxy)){
      
      ref_dates <- time_index
      
      if (!all(ref_dates %in% index(vol_proxy))) {
        
        stop(paste0("vol_proxy data do not cover daily_ret"))
        
      } else {
        
        vol_proxy <- vol_proxy[ref_dates]
      }        
    }
    
    required_d_m <- d_m
    RV_m <- NULL
    mv_m <- NULL
    mv_c <- NULL
    
  }
  
  dccmidas_noskew_est<-dccmidas_noskew_longfocus_fit(r_t=required_d_m,univ_model=univ_model,distribution=distribution,X_stru_break_date_vec=X_stru_break_date_vec,
                                                     RV_m=RV_m,mv_m=mv_m,K=K,mv_c=mv_c,corr_model=corr_model,lag_fun=lag_fun,
                                                     N_c=N_c,K_c=K_c,out_of_sample=out_of_sample,vol_proxy=vol_proxy,R=R)
  
  
  return(dccmidas_noskew_est)
  
}


