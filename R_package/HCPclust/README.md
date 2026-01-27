# HCPclust

**HCPclust** is an R package for **Hierarchical Conformal Prediction (HCP)**, designed
for clustered data where each subject has repeated measurements with
within-subject dependence, and outcomes can be missing under a
covariate-dependent **Missing At Random (MAR)** mechanism.


## Installation

You can install the development version of **HCPclust** directly from GitHub using:

```r
install.packages("remotes")
remotes::install_github("judywangstat/HCP")
library(HCPclust)
```
## Method overview

The proposed method consists of four steps.

**Step (1) Model fitting**

- Fit a pooled conditional density model $\hat{\pi}(y \mid x)$ using `fit_cond_density_quantile()`.
- Fit a marginal missingness propensity model $\hat{p}(x) = \mathbb{P}(\delta = 1 \mid x)$ using `fit_missingness_propensity()`.

Both models are estimated on a subject-level training split to capture the
distribution of $Y \mid X$ and correct for covariate-dependent missingness.

**Step (2) Subsampled calibration**

- Construct a calibration set at the subject level.
- Repeatedly draw one observation per subject to form calibration samples.

This preserves the clustered structure while reducing the impact of
within-subject dependence on conformal calibration.

**Step (3) Weighted conformal scoring**

- Compute the nonconformity score $R(x, y) = - \hat{\pi}(y \mid x)$.
- Apply inverse-propensity weighting $w(x) = 1 / \hat{p}(x)$.

**Step (4) Aggregation**

- Aggregate p-values across subsamples (B) and subject-level splits (S),
  using either the Cauchy combination test (CCT / ACAT) or a simple mean.

The final conformal prediction region contains all values $y$ such that $p_{\text{final}}(y) > \alpha.$

## What’s included

### 1) Data generator
- **`generate_clustered_mar()`**
  - Simulates clustered (longitudinal) data under a random-effects model with
    covariate-dependent Missing At Random (MAR) outcomes.
  - Supports:
    - subject-specific random intercepts (clustered structure),
    - optional heteroskedastic errors controlled by `hetero_gamma`,
    - optional within-cluster AR(1) dependence via `rho`,
    - optional calibration to a target marginal missing rate via `target_missing`.

### 2) Conditional density estimator
- **`fit_cond_density_quantile()`**
  - Fits a conditional quantile process $Q(\tau \mid x)$ using pooled observed data.
  - Converts estimated quantiles into a conditional density via a quotient estimator
    along the quantile curve.
  - Supported engines:
    - `method = "rq"`: linear quantile regression (`quantreg`)
    - `method = "qrf"`: quantile random forest (`quantregForest`)
  - Optional monotone (isotonic) adjustment and tail decay extrapolation
    for numerical stability.

### 3) Missingness propensity model
- **`fit_missingness_propensity()`**
  - Estimates the missingness propensity $p(x) = \mathbb{P}(\delta = 1 \mid X)$ from pooled data.
  - Supported methods:
    - `logistic`: logistic regression (`glm`)
    - `grf`: generalized random forests (`grf::probability_forest`)
    - `boosting`: gradient boosting (`xgboost`)
  - Predicted propensities are clipped to `[eps, 1 - eps]` to stabilize
    inverse-propensity weighting.

### 4) HCP conformal prediction
- **`hcp_conformal_region()`**
  - Core HCP procedure for constructing a conformal prediction region
    for a single new target point.
  - Implements the four main steps:
    Step (1). Model fitting,
    Step (2). Subsampled calibration,
    Step (3). Weighted conformal scoring,
    Step (4). Aggregation.
  - Returns the prediction region as a subset of `y_grid`, together with
    the interval endpoints `[lo, hi]`.

- **`hcp_predict_targets()`**
  - A flexible wrapper built on `hcp_conformal_region()` that automatically adapts
    to different test-data structures, including:
    - a single patient at a single target point,
    - a single patient at multiple target points,
    - multiple patients with one target point each,
    - multiple patients with multiple target points each.
  - Supports optional **per-patient Bonferroni correction**, using
    `alpha_eff = alpha / M_p` when a patient has `M_p > 1` target points.

### 5) Plotting
- **`plot_hcp_intervals()`**
  - `mode = "band"`: for one subject with repeated measurements (e.g., over time),
    plots an interval band as a function of a covariate.
  - `mode = "pid"`: for multiple subjects with one target point each,
    displays one prediction interval per subject, optionally sorted by a covariate.

## Quick examples

This section shows (i) how to call `hcp_predict_targets()` under the four common test-data scenarios,
and (ii) how to generate and display two plots from `plot_hcp_intervals()`.

### A) Four scenarios for testing data (via `hcp_predict_targets()`)
####  Generate training data (fixed across all cases)
```r
set.seed(1)
dat_train <- generate_clustered_mar(
  n = 200, m = 4, d = 1,
  x_dist = "uniform", x_params = list(min = 0, max = 10),
  target_missing = 0.30,
  seed = 1
)
y_grid <- seq(-6, 6, length.out = 201)
```
####  Case 1: one patient, one target
```r
test_11 <- data.frame(pid = 1, X1 = 2.5)
out_11 <- hcp_predict_targets(
  dat = dat_train, test = test_11,
  x_cols = "X1", y_grid = y_grid,
  alpha = 0.1,
  seed = 1
)
out_11$pred
```

####  Case 2: one patient, multiple targets
```r
test_1M <- data.frame(pid = 1, X1 = c(1, 3, 7, 9))
out_1M <- hcp_predict_targets(
  dat = dat_train, test = test_1M,
  x_cols = "X1", y_grid = y_grid,
  alpha = 0.1,
  seed = 1
)
out_1M$pred
```

####  Case 3: multiple patients, one target each
```r
test_P1 <- data.frame(pid = 1:4, X1 = c(2, 4, 6, 8))
out_P1 <- hcp_predict_targets(
  dat = dat_train, test = test_P1,
  x_cols = "X1", y_grid = y_grid,
  alpha = 0.1,
  seed = 1
)
out_P1$pred
```

####  Case 4: multiple patients, multiple targets per patient
```r
test_PM <- data.frame(
  pid = c(1,1, 2,2,2, 3,3),
  X1  = c(1,6,  2,5,9,  3,8)
)
out_PM <- hcp_predict_targets(
  dat = dat_train, test = test_PM,
  x_cols = "X1", y_grid = y_grid,
  alpha = 0.1,
  seed = 1
)
out_PM$pred
```

### B)  Two plots
####  Generate training data and test data
```r
dat_train <- generate_clustered_mar(
  n = 200, m = 20, d = 1,
  x_dist = "uniform", x_params = list(min = 0, max = 10),
  hetero_gamma = 2.5,
  target_missing = 0.30,
  seed = 1
)
y_grid <- seq(-6, 10, length.out = 201)

# test data with latent truth
dat_test <- generate_clustered_mar(
  n = 100, m = 20, d = 1,
  x_dist = "uniform", x_params = list(min = 0, max = 10),
  hetero_gamma = 2.5,
  seed = 999
)
```
####  Plot A: one patient, multiple targets (band vs time)
```r
pid <- dat_test$id[1]
idx <- which(dat_test$id == pid)
idx <- idx[order(dat_test$X1[idx])][1:10]
test_1M <- data.frame(pid = pid, X1 = dat_test$X1[idx], y_true = dat_test$Y_full[idx])

out_1M <- hcp_predict_targets(
  dat = dat_train, test = test_1M,
  x_cols = "X1", y_grid = y_grid,
  alpha = 0.1,
  S = 3, B = 3,
  seed = 1
)
plot_hcp_intervals(
  out_1M$pred, mode = "band", x_col = "X1",
  y_true_col = "y_true", show_true = TRUE,
  main = "Case A: one patient, multiple time points (band vs time)"
)
```
<p align="center">
  <img src="figures/band_example.png" width="600">
</p>

#### Plot B: multiple patients, one target each (intervals by patient) 
```r
pids <- unique(dat_test$id)[1:20]
test_P1 <- subset(dat_test, id %in% pids & j == 1, select = c(id, X1, Y_full))
names(test_P1) <- c("pid", "X1", "y_true")

out_P1 <- hcp_predict_targets(
  dat = dat_train, test = test_P1,
  x_cols = "X1", y_grid = y_grid,
  alpha = 0.1,
  S = 3, B = 3,
  seed = 1
)
plot_hcp_intervals(
  out_P1$pred, mode = "pid", pid_col = "pid", x_sort_col = "X1",
  y_true_col = "y_true", show_true = TRUE,
  main = "Case B: multiple patients, one time point (by patient)"
)
```
<p align="center">
  <img src="figures/pid_example.png" width="700">
</p>
