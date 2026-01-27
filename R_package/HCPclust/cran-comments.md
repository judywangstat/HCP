## Resubmission
- This is a resubmission. In response to the CRAN incoming pretest NOTE on Debian
  (“Examples with CPU time > 2.5 times elapsed time”), the examples for
  `fit_missingness_propensity()` were revised to reduce computational burden and
  to enforce single-thread settings for GRF and xgboost in the optional (donttest)
  section.

## Test environments
- local macOS (R 4.4.2)
- win-builder (R devel)

## R CMD check results
0 errors | 0 warnings | 1 note

## Notes
- The remaining NOTE on win-builder is the standard “New submission” message from
  the CRAN incoming feasibility check.

## Downstream dependencies
There are no downstream dependencies.
