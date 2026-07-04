# gmxdcc 参考手册 / Reference Manual

**Package:** gmxdcc 0.1.0 — GARCH-MIDAS and DCC-MIDAS Models with Exogenous
Covariates and Structural Breaks

本手册对应 CRAN 参考手册的角色:逐函数解说参数的设计方法与含义,附完整数学模型与可运行示例(全部基于包内置数据)。

English version: [gmxdcc-manual-en.md](gmxdcc-manual-en.md).

---

## 目录

1. [模型体系与两步估计流程](#1-模型体系与两步估计流程)
2. [数学模型](#2-数学模型)
3. [`fit_garch_midas()` — 单变量 GARCH-MIDAS](#3-fit_garch_midas)
4. [`fit_dcc()` — 两步 DCC / DCC-MIDAS](#4-fit_dcc)
5. [优化器设计:多起点、约束与热重启](#5-优化器设计)
6. [样本外(OOS)设计](#6-样本外设计)
7. [辅助函数](#7-辅助函数)
8. [内置数据集](#8-内置数据集)
9. [旧脚本迁移对照表](#9-旧脚本迁移对照表)

---

## 1. 模型体系与两步估计流程

```
日频收益面板 r_t (N 资产)
        │
第一步  │ 逐资产估计单变量波动率模型
        │   rugarch 族: sGARCH / gjrGARCH / eGARCH / iGARCH / csGARCH
        │   或 GARCH-MIDAS (本包实现, 任意低频协变量 + 结构突变)
        ▼
标准化残差 ε_t = (r_t − μ) / σ_t,  条件标准差对角阵 D_t
        │
第二步  │ 在 ε_t 上估计相关性模型
        │   cDCC / aDCC / DECO (C++ 内核)
        │   或 DCC-MIDAS (已实现相关 RC + 任意协变量 + 结构突变)
        ▼
条件相关阵 R_t,  条件协方差阵 H_t = D_t R_t D_t
```

两步均为约束 QMLE(准极大似然),两步都报告 Bollerslev–Wooldridge 三明治稳健标准误。第二步标准误未修正第一步参数不确定性,与该文献的通行做法一致(Engle 2002)。

**设计哲学**:旧研究脚本用模型名后缀枚举变体(`GM_noskew_RV_X_sb` 等 6×2×2 组合)。本包改为**正交的参数开关**——协变量给列表、突变给日期向量、分布和权重函数各给一个参数——任何组合自动成立,旧模型都是特例(见第 9 节映射表)。

---

## 2. 数学模型

### 2.1 单变量 GARCH-MIDAS

日收益分解为长短两个波动率成分:

$$r_t = \mu + \sqrt{\tau_t\, g_t}\;\epsilon_t,\qquad \epsilon_t \sim \mathcal{N}(0,1)\ \text{或}\ t_\nu$$

**短期成分**(单位均值 GARCH(1,1) 递推):

$$g_t = (1-\alpha-\beta) + \alpha\,\frac{(r_{t-1}-\mu)^2}{\tau_{t-1}} + \beta\, g_{t-1},\qquad g_1 = 1$$

**长期成分**(MIDAS 滤波,任意 $J$ 个低频协变量):

$$\log \tau_t = m(t) + \sum_{j=1}^{J} \theta_j(t) \sum_{k=1}^{K_j} \varphi_k(\omega_j)\, X_{j,\,t-k}$$

其中 $X_{j,t-k}$ 是第 $j$ 个协变量在交易日 $t$ 所属低频期(月/季)之前第 $k$ 期的取值;$K_j$ 为各协变量独立的滞后阶数。

**MIDAS 权重** $\varphi_k(\omega)$ 两种方案(与 rumidas 数值一致):

- Beta(受限,$w_1=1$,单调递减):
  $$\varphi_k(\omega) = \frac{(1-k/K)^{\omega-1}}{\sum_{i=1}^{K}(1-i/K)^{\omega-1}},\qquad \omega > 1$$
- 指数 Almon($w_1=0$):
  $$\varphi_k(\omega) = \frac{\exp(\omega k^2)}{\sum_{i=1}^{K}\exp(\omega i^2)},\qquad \omega < 0$$

**结构突变**:给定突变日期 $d_1 < \dots < d_B$,构造阶梯虚拟矩阵 $D\in\mathbb{R}^{(1+B)\times T}$(第 1 行全 1,第 $b{+}1$ 行在 $t \ge d_b$ 时为 1),则

$$m(t) = \mathbf{m}^\top D_{\cdot t},\qquad \theta_j(t) = \boldsymbol{\theta}_j^\top D_{\cdot t}$$

即截距与全部斜率都是分段常数;$B=0$ 时退化为常数——**无突变模型是突变模型的特例**,这正是包里用同一套代码统一新旧 12 个变体的关键。

### 2.2 相关性模型

第一步产出标准化残差 $\varepsilon_t\in\mathbb{R}^N$。

**cDCC**(Aielli 修正的 DCC):

$$Q_t = (1-a-b)\,S + a\, Q_{t-1}^{*1/2}\varepsilon_{t-1}\varepsilon_{t-1}^\top Q_{t-1}^{*1/2} + b\, Q_{t-1},\qquad R_t = Q_t^{*-1/2} Q_t\, Q_t^{*-1/2}$$

其中 $S$ 为残差无条件协方差(锚点),$Q_t^* = \mathrm{diag}(Q_t)$。

**aDCC**(不对称项):

$$Q_t = (1-a-b)S - \gamma\bar N + a\,\varepsilon_{t-1}\varepsilon_{t-1}^\top + b\,Q_{t-1} + \gamma\,\eta_{t-1}\eta_{t-1}^\top,\qquad \eta_t = \min(\varepsilon_t, 0)$$

**DECO**(等相关):先算 cDCC 的 $R_t$,再取横截面平均相关 $\rho_t$,令 $R_t^{deco} = (1-\rho_t)I + \rho_t J$。

**DCC-MIDAS**:短期 Q 递推同上,但锚点换成慢变的长期相关 $\bar R_t$:

$$Q_t = (1-a-b)\,\bar R_t + a\,\varepsilon_{t-1}\varepsilon_{t-1}^\top + b\, Q_{t-1}$$

长期相关 $\bar R_t$ 有三种函数形式(沿袭旧脚本设计):

1. **纯 RC、无突变**(`rc=TRUE, corr_X=NULL, corr_break_dates=NULL`):
   $$\bar R_{ij,t} = \sum_{k=1}^{K_c} \varphi_k(\omega_{RC})\, C_{ij,\,t-k}$$
   其中 $C_t$ 是窗宽 $N_c$ 的滚动已实现相关阵(无 tanh 变换)。
2. **含协变量、无突变**:逐资产对参数化,经 Fisher 型变换保证 $|\bar R_{ij}|<1$:
   $$\bar R_{ij,t} = \tanh\!\Big(\delta_{ij} + \textstyle\sum_j \theta_{ij}^{(j)}\,\mathrm{MIDAS}(X_j;\omega_{ij}^{(j)})_t \;[+\; \theta^{RC}_{ij}\, RC_{ij,t}]\Big)$$
3. **含突变**:水平/斜率系数跨资产对共享、按区制分段($m(t), \theta(t)$ 由 $D$ 矩阵驱动),权重 $\omega$ 仍逐对:
   $$\bar R_{ij,t} = \tanh\!\Big(m(t) + \theta_{RC}(t)\,RC_{ij,t} + \textstyle\sum_j \theta_j(t)\,\mathrm{MIDAS}(X_j;\omega_{ij}^{(j)})_t\Big)$$

### 2.3 两步 QMLE 与稳健标准误

每步的逐观测对数似然向量喂给 `maxLik`(带不等式约束的 BFGS)。稳健方差为三明治估计:

$$\widehat{\mathrm{Var}}(\hat\vartheta) = H^{-1} \Big(\textstyle\sum_t s_t s_t^\top\Big) H^{-1}$$

$H$ 为数值 Hessian,$s_t$ 为逐观测得分(由 `maxLik` 的 `gradientObs` 提供)。GARCH 类模型创新分布几乎必然误设,故不使用普通 Hessian 标准误(Bollerslev & Wooldridge, 1992)。

---

## 3. `fit_garch_midas()`

```r
fit_garch_midas(returns, X = NULL, K = 12, rv = is.null(X), K_rv = NULL,
                period = c("monthly", "quarterly"),
                distribution = c("norm", "std"),
                lag_fun = c("Beta", "Almon"),
                break_dates = NULL,
                out_of_sample = NULL, vol_proxy = NULL,
                n_starts = 100, seed = NULL, control = list())
```

### 参数逐一解说

| 参数 | 类型/默认 | 意义与设计动机 |
|---|---|---|
| `returns` | xts,单列 | 日频收益。强制 xts 是为了让包能读取日历做 MIDAS 对齐 |
| `X` | NULL / xts / 命名列表 | 长期成分的低频协变量,**任意个数**。设计为列表而非固定槽位,使"三个宏观变量"与"一个"用同一接口;列表名字会进入参数名(如 `theta_JPMU`) |
| `K` | 12 | MIDAS 滞后阶数。标量自动循环补齐到每个协变量,或给向量为每个协变量单独设 $K_j$ |
| `rv` | `is.null(X)` | 是否把内生已实现波动率(由 `realized_measure()` 从 `returns` 自动算出)作为第一个协变量。默认逻辑:不给 X 时至少要有 RV,否则长期成分无驱动 |
| `K_rv` | `K[1]` | RV 协变量单独的滞后阶数 |
| `period` | "monthly" | 低频频率(月/季),决定 RV 聚合与滞后矩阵的日历规则 |
| `distribution` | "norm" | 创新分布;`"std"` 增加自由度参数 $\nu$(约束 $\nu>2.001$) |
| `lag_fun` | "Beta" | 权重函数;决定 $\omega$ 的初值与约束方向(Beta: $\omega>1.001$;Almon: $\omega<-0.001$) |
| `break_dates` | NULL | 结构突变日期向量。**取代旧脚本的 `_sb` 模型名后缀**;NULL=无突变,$B$ 个日期把 $m$ 与全部 $\theta_j$ 分成 $B{+}1$ 段 |
| `out_of_sample` | NULL | 留出末尾 n 天不参与估计;成分经全样本滤波热启动(见第 6 节) |
| `vol_proxy` | NULL | 方差代理(如已实现方差)用于 MSE/QLIKE 损失评估;缺省用日收益平方 |
| `n_starts` | 100 | 随机起点个数(对应旧脚本的 `R`) |
| `seed` | NULL | 随机起点种子。**旧脚本无种子、结果不可复现**;设种子是新增的可复现性保障 |
| `control` | list() | 优化覆盖项:`bounds`(覆盖初值区间/约束常数表)、`start_matrix`(注入起点,见第 5 节)、`iterlim`、`method` |

### 参数向量的规范排序(`gm_par_index()`)

所有内部函数、`coef()`、`start_matrix` 都遵循同一排序:

$$\underbrace{\alpha,\ \beta}_{\text{短期}},\ \underbrace{\omega_1..\omega_J}_{\text{各协变量权重}},\ \underbrace{\nu}_{\text{仅 std}},\ \underbrace{m_1..m_{B+1}}_{\text{分段截距}},\ \underbrace{\theta_{1,1}..\theta_{J,B+1}}_{\text{各协变量分段斜率}},\ \mu$$

命名规则:`alpha, beta, w2_<协变量名>, nu, m_1.., theta_<协变量名>_1.., mu`。

### 返回对象(class `gm_fit`)

| 字段 | 内容 |
|---|---|
| `rob_coef_mat` | Estimate / 稳健SE / t / p 四列表(行名=参数名) |
| `loglik`, `inf_criteria`, `obs`, `period` | 样本内似然、AIC/BIC、样本量、区间 |
| `est_vol_in_s`, `est_lr_in_s`, `est_sr_in_s` | 总条件波动率 $\sqrt{h_t}$、长期 $\tau_t$、短期 $g_t$ |
| `loss_in_s`(及 `_oos` 系列) | MSE(%)、QLIKE |
| `maxlik` | 原始 maxLik 对象(供 `vcov()` 等使用) |

S3 方法:`summary` `print` `coef` `logLik` `vcov`。

### 示例

```r
library(gmxdcc)

fit <- fit_garch_midas(
  topix_market,
  X = list(JPMU = jp_macro_uncertainty),
  K = 24, rv = FALSE,
  break_dates = jp_mu_breaks,
  distribution = "norm", lag_fun = "Beta", period = "monthly",
  seed = 1
)
summary(fit)          # 系数表 + 显著性 + 损失
coef(fit)             # 命名系数向量(规范排序)
sqrt(vcov(fit))       # 三明治协方差阵
```

---

## 4. `fit_dcc()`

```r
fit_dcc(returns,
        univariate = c("sGARCH","gjrGARCH","eGARCH","iGARCH","csGARCH","garch_midas"),
        univariate_args = list(),
        correlation = c("cDCC", "aDCC", "DECO", "DCC-MIDAS"),
        rc = TRUE, N_c = NULL, corr_X = NULL, K_c = NULL,
        corr_break_dates = NULL,
        distribution = c("norm", "std"), lag_fun = c("Beta", "Almon"),
        period = c("monthly", "quarterly"),
        out_of_sample = NULL, vol_proxy = NULL,
        n_starts = 100, seed = NULL, control = list())
```

### 参数逐一解说

| 参数 | 类型/默认 | 意义与设计动机 |
|---|---|---|
| `returns` | xts,≥2 列 | 多资产日收益面板 |
| `univariate` | "sGARCH" | 第一步模型。rugarch 族走 `rugarch::ugarchfit`;`"garch_midas"` 走本包 `fit_garch_midas()` 逐资产估计 |
| `univariate_args` | list() | **第一步专属**参数包:`X`、`K`、`rv`、`K_rv`、`break_dates` 同 `fit_garch_midas()`;另可含 `start_matrix`(逐资产起点注入:单个矩阵作用于全部资产,或按资产给列表) |
| `correlation` | "cDCC" | 第二步模型 |
| `rc` | TRUE | DCC-MIDAS 专用:长期相关是否含已实现相关项 |
| `N_c` | NULL | 已实现相关的滚动窗宽(天);`rc=TRUE` 时必填 |
| `corr_X` | NULL | 长期相关的低频协变量(列表,同 `X` 的设计) |
| `K_c` | NULL | 相关侧 MIDAS 滞后阶数(标量);DCC-MIDAS 必填。**cDCC/aDCC/DECO 请保持 NULL**——C++ 内核对 NULL 有自己的旧默认偏移,显式传值会移动递推起点 |
| `corr_break_dates` | NULL | 长期相关的突变日期(分段 $m,\theta$;见 2.2 第 3 式) |
| 其余 | — | 与 `fit_garch_midas()` 同义;`control$start_matrix` **只作用于第二步**(列序按 `dcc_par_index()`) |

### 第二步参数排序(`dcc_par_index()`,与旧脚本各配置逐一吻合)

设 $P = N(N-1)/2$ 个资产对、$R = 1+B$ 个区制:

- **无突变、含协变量**:`a, b, (w2_RC,) delta_1..P, theta_<X>_1..P, w2_<X>_1..P, (theta_RC_1..P)`
- **有突变**:`a, b, (w2_RC,) w2_<X>_1..P, m_1..R, (theta_RC_1..R,) theta_<X>_1..R`
- **纯 RC 无突变**:仅 `a, b, w2_RC`

### 返回对象(class `dcc_fit`)

| 字段 | 内容 |
|---|---|
| `est_univ_model` | 各资产第一步系数表(列表) |
| `corr_coef_mat` | 第二步系数表(稳健SE) |
| `eps_t`, `D_t` | 标准化残差、条件SD对角阵序列 |
| `H_t`, `R_t`, `R_t_bar` | 条件协方差/相关/长期相关阵($N\times N\times T$) |
| `*_oos` 系列 | 样本外切片(热启动滤波所得) |
| `u_llk`, `c_llk`, `Inf_criteria_2step` | 两步似然与系统 AIC/BIC |

### 示例

```r
panel <- xts::merge.xts(topix_autos, topix_market)
panel <- panel[complete.cases(panel), ]

dcc <- fit_dcc(
  panel,
  univariate = "garch_midas",
  univariate_args = list(X = jp_macro_uncertainty, K = 24, rv = FALSE,
                         break_dates = jp_mu_breaks),
  correlation = "DCC-MIDAS", rc = FALSE,
  corr_X = list(JPMU = jp_macro_uncertainty), K_c = 24,
  corr_break_dates = jp_mu_breaks,
  seed = 1
)
summary(dcc)
# 条件相关路径 / conditional correlation path
plot(xts::xts(dcc$R_t[1, 2, ], order.by = dcc$Days))
```

---

## 5. 优化器设计

### 5.1 多起点机制

似然面多峰(尤其含突变时)。沿袭旧脚本方案:抽 `n_starts` 个候选起点 → 逐一算总似然 → 从最优者出发做**一次**带约束 BFGS。随机参数按下表区间均匀抽取,固定参数(权重 $\omega$、均值 $\mu$)取固定初值。

### 5.2 常数表(可经 `control$bounds` 覆盖)

**单变量步**(`gm_bounds_default`):

| 项 | 初值 | 约束 |
|---|---|---|
| $\alpha$ | U(0.001, 0.095) | $>0.001$ |
| $\beta$ | U(0.6, 0.8) | $>0.001$;$\alpha+\beta<0.999$ |
| $m,\theta$(每区制) | U(−1, 1) | 无 |
| $\omega$(Beta / Almon) | 1.01 / −0.1 | $>1.001$ / $<-0.001$ |
| $\nu$ | U(2.01, 10) | $>2.001$ |
| $\mu$ | 样本均值(固定) | 无 |

**相关步**(`dcc_bounds_default`;注意 $a$ 下界与第一步不同,沿袭旧脚本):

| 项 | 初值 | 约束 |
|---|---|---|
| $a$ | U(0.001, 0.095)(多起点)或 0.01(单起点模型) | $>0.0001$ |
| $b$ | U(0.6, 0.8) 或 0.8 | $>0.001$;$a+b<0.999$ |
| $\gamma$(aDCC) | 0.01 | $0.001<\gamma<0.15$ |
| $\delta,\theta,m$ | U(−1, 1) | 无 |
| $\omega$(Beta / Almon) | 1.2 / −0.05 | $>1.001$ / $<-0.001$ |

cDCC / aDCC / DECO / 纯 RC 为单固定起点模型(与旧脚本一致,`n_starts` 对它们无效)。

### 5.3 起点注入与热重启(滚动重估的关键)

- `control$start_matrix`:注入**第二步**起点(行=候选,列序按 `dcc_par_index()`);
- `univariate_args$start_matrix`:注入**第一步**起点(逐资产,列序按 `gm_par_index()`);
- 注入起点若贴着约束边界(舍入后的 `coef()` 常见),会被自动向可行域内部做最小混合修正——`constrOptim` 的对数障碍法要求严格内点,否则报 "initial value not finite"。

典型滚动重估用法:

```r
prev_u <- NULL; prev_c <- NULL
for (w in windows) {
  fit <- fit_dcc(data_w, univariate = "garch_midas",
                 univariate_args = list(..., start_matrix = prev_u),
                 control = list(start_matrix = prev_c),
                 n_starts = if (is.null(prev_c)) 300 else 1, ...)
  prev_u <- lapply(fit$est_univ_model, function(m) matrix(m[,"Estimate"], nrow=1))
  prev_c <- matrix(coef(fit), nrow = 1)
}
```

---

## 6. 样本外设计

`out_of_sample = n` 时:

1. 参数只用前 $T-n$ 天估计(信息集干净);
2. 全部条件量($g_t,\tau_t,Q_t,R_t,H_t$)由**冻结参数在完整样本上滤波**后切片——样本外递推自然衔接样本内末期状态(**热启动**),等价于逐日一步向前预测的"固定方案"(fixed estimation scheme);
3. 无条件锚点($S$、aDCC 的 $\bar N$)只用样本内残差计算(C++ 内核的 `S_in`/`N_in` 参数),**无前视偏差**。

旧脚本在样本外段冷启动($g\!\to\!1$、$Q\!\to\!I$)且 $S$ 用样本外残差计算;新设计才符合文献惯例,NEWS.md 有完整记录。

---

## 7. 辅助函数

### `realized_measure(x, period)`
把日收益聚合为低频已实现方差 $RV_\tau = \sum_{t\in\tau} r_t^2$,索引对齐到期首日。
```r
rv <- realized_measure(topix_market, "monthly")
```

### `align_sample(r_t, K, K_c, macro, period)`
裁剪日收益,使每个保留日都有 $\max(K,K_c)+1$ 个完整低频期的协变量历史(MIDAS 滤波的前提)。

### `midas_lag_matrix(x, mv, K, period)`
为每个交易日堆叠协变量当期与前 $K$ 期取值,得 $(K{+}1)\times T$ 滞后矩阵——`build_tau` 的直接输入。

### `midas_weights(k, K, w1, w2, lag_fun)`
第 2.1 节的权重公式本体(与 rumidas 数值一致,有单元测试钉住)。
```r
midas_weights(1:12, 12, w1 = 1, w2 = 3, lag_fun = "Beta")
```

### `build_break_dummies(dates, break_dates)`
第 2.1 节的阶梯虚拟矩阵 $D$;`NULL` 突变返回单行全 1。

### `loss_avg(vol_est, vol_proxy)`
$$\mathrm{MSE\%} = 100\cdot\overline{(h - h^{proxy})^2},\qquad \mathrm{QLIKE} = \overline{\log h + h^{proxy}/h}$$

### `qmle_se(est)` / `inf_criteria(est)` / `inf_criteria_2step(res_garch, res_dcc)`
三明治标准误(2.3 节公式)、单模型 AIC/BIC、两步系统 AIC/BIC($k$ 与似然按两步加总)。

---

## 8. 内置数据集

| 数据集 | 内容 | 频率 |
|---|---|---|
| `topix_market` | TOPIX 市场收益(列 `TOPIX`) | 日 |
| `topix_autos` | TOPIX-17 汽车与运输设备行业收益(列 `AUTOS`) | 日 |
| `jp_macro_uncertainty` | 日本宏观不确定性指数(列 `JPMU`,标准化) | 月 |
| `jp_mu_breaks` | JPMU 均值的 4 个结构突变日期(`strucchange::breakpoints`) | — |

全部为 Date 索引的干净 xts,`library(gmxdcc)` 后直接可用。

---

## 9. 旧脚本迁移对照表

| 旧模型名 | 新调用 |
|---|---|
| `GM_noskew_RV` | `rv = TRUE, X = NULL` |
| `GM_noskew_X` | `rv = FALSE, X = list(x)` |
| `GM_noskew_RV_X` | `rv = TRUE, X = list(x)` |
| `GM_noskew_*_sb` | 加 `break_dates = 日期向量` |
| `cDCC` / `aDCC` / `DECO` | `correlation = 同名` |
| `DCCMIDAS_RC` | `correlation = "DCC-MIDAS", rc = TRUE` |
| `DCCMIDAS_X` | `correlation = "DCC-MIDAS", rc = FALSE, corr_X = list(x)` |
| `DCCMIDAS_RC_X` | `correlation = "DCC-MIDAS", rc = TRUE, corr_X = list(x)` |
| `DCCMIDAS_*_sb` | 加 `corr_break_dates` |

参数改名:`X_stru_break_date_vec → break_dates / corr_break_dates`;`type → period`;`R → n_starts`。

**行为差异**(有意改进,详见 NEWS.md):样本外热启动取代冷启动;无条件锚点仅用样本内残差;多起点可设种子复现。旧代码的公式层等价性由测试套件保证(固定参数处逐观测对齐至 1e-12)。

---

*本手册与包版本 0.1.0 同步;函数签名以 `?fit_garch_midas`、`?fit_dcc` 的最新帮助页为准。*
