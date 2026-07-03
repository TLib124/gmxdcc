## Submission

This is a new release (first submission of gmxdcc).

gmxdcc provides two-step quasi-maximum-likelihood estimation of
GARCH-MIDAS volatility models and cDCC/aDCC/DECO/DCC-MIDAS correlation
models, with flexible low-frequency covariates and structural breaks in
the long-run components.

## Test environments

* local: Windows 11, R 4.5.2 (x86_64-w64-mingw32)
* win-builder: R-devel  <!-- 跑完 check_win_devel() 后把结果邮件里的 R 版本号补在这里 -->

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes for the reviewer

* The `\donttest{}` examples of `fit_garch_midas()` and `fit_dcc()` run a
  full (constrained multi-start QML) model fit and take roughly 10-60
  seconds each; all other examples run in well under 5 seconds.
* The package bundles ~170 KB of compressed sample data (daily Japanese
  equity index returns and a monthly macro uncertainty index) used by the
  examples and the vignette.
