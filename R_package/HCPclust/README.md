# HCPclust

**HCPclust** is an R package for **Hierarchical Conformal Prediction (HCP)**. It targets **clustered / longitudinal** data where each subject has repeated measurements with **within-subject dependence**, and the outcome can be missing under a covariate-dependent **Missing At Random (MAR)** mechanism.


## Installation

You can install the development version of **HCPclust** directly from GitHub using:

```r
install.packages("remotes")
remotes::install_github("judywangstat/HCP", subdir = "R_package/HCPclust")
library(HCPclust)
```
## Method overview

The HCP procedure has four main steps.

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
- Apply inverse-propensity weighting $w(x) = 1 / \hat{p}(x)$ under MAR.

**Step (4) Aggregation**

- Aggregate p-values across subsamples (B) and subject-level splits (S),
  using either the Cauchy combination test (CCT / ACAT) or a simple mean.

The final conformal prediction region contains all values $y$ such that $p_{\text{final}}(y) > \alpha.$

## What’s included in 'HCPclust'

### 1) Data generator
- **`generate_clustered_mar()`**
  - Simulates clustered (longitudinal) data under a random-effects model with
    covariate-dependent MAR outcomes.
  - Supports:
    - subject-specific random intercepts (clustered structure),
    - optional heteroskedastic errors controlled by `hetero_gamma`,
    - optional within-cluster AR(1) dependence via `rho`,
    - optional calibration to a target marginal missing rate via `target_missing`.

### 2) Conditional density estimator
- **`fit_cond_density_quantile()`**
  - Fits a conditional quantile process $Q(\tau \mid x)$ using pooled observed data.
  - Converts estimated quantiles into a conditional density via a quotient estimator
    along the quantile process.
  - Supported engines:
    - `method = "rq"`: linear quantile regression (`quantreg`)
    - `method = "qrf"`: quantile random forest (`quantregForest`)
  - Optionally applies monotone (isotonic) adjustment and tail-decay extrapolation for numerical stability.

### 3) Missingness propensity estimator
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
  - Implements the four steps: model fitting, subsampled calibration, weighted scoring, and aggregation.
  - Returns the prediction region as a subset of `y_grid`, together with
    the interval endpoints `[lo, hi]`.

- **`hcp_predict_targets()`**
  - A flexible wrapper built on `hcp_conformal_region()` that automatically adapts
    to different test-data structures, including:
    - one patient, one measurement;
    - one patient, multiple measurements;
    - multiple patients, one measurement;
    - multiple patients, multiple measurements.
  - Supports optional **per-patient Bonferroni correction**, using
    `alpha_eff = alpha / M_p` when a patient has `M_p > 1` measurement points.

### 5) Plotting
- **`plot_hcp_intervals()`**
  - `mode = "band"`: for one subject with repeated measurements (e.g., over time),
    plots an interval band as a function of a covariate.
  - `mode = "pid"`: for multiple subjects with one measurement point each,
    displays one prediction interval per subject, optionally sorted by a covariate.

## Quick examples

This section shows (i) how to call `hcp_predict_targets()` under the four common test-data scenarios,
and (ii) how to generate two plots from `plot_hcp_intervals()`.

### A) Four scenarios for testing data (via `hcp_predict_targets()`)
####  Generate training data (shared across all cases)
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
####  Case 1: test one patient, one measurement
```r
test_11 <- data.frame(pid = 1, X1 = 2.5)
out_11 <- hcp_predict_targets(
  dat = dat_train, test = test_11,
  x_cols = "X1", y_grid = y_grid,
  alpha = 0.1,
  S = 2, B = 2,
  seed = 1
)
out_11$pred
```

####  Case 2: test one patient, multiple measurements
```r
test_1M <- data.frame(pid = 1, X1 = c(1, 3, 7, 9))
out_1M <- hcp_predict_targets(
  dat = dat_train, test = test_1M,
  x_cols = "X1", y_grid = y_grid,
  alpha = 0.1,
  S = 2, B = 2,
  seed = 1
)
out_1M$pred
```

####  Case 3: test multiple patients, one measurement
```r
test_P1 <- data.frame(pid = 1:4, X1 = c(2, 4, 6, 8))
out_P1 <- hcp_predict_targets(
  dat = dat_train, test = test_P1,
  x_cols = "X1", y_grid = y_grid,
  alpha = 0.1,
  S = 2, B = 2,
  seed = 1
)
out_P1$pred
```

####  Case 4: test multiple patients, multiple measurements per patient
```r
test_PM <- data.frame(
  pid = c(1,1, 2,2,2, 3,3),
  X1  = c(1,6,  2,5,9,  3,8)
)
out_PM <- hcp_predict_targets(
  dat = dat_train, test = test_PM,
  x_cols = "X1", y_grid = y_grid,
  alpha = 0.1,
  S = 2, B = 2,
  seed = 1
)
out_PM$pred
```

### B)  Plots for two cases
####  Load the CD4 dataset as an example (lqmix::cd4)
```r
if (!requireNamespace("lqmix", quietly = TRUE)) {
  stop("Package 'lqmix' is required for this example.")
}
data(cd4, package = "lqmix")

dat0 <- data.frame(
  id = cd4$sbj.id,
  X1 = cd4$time,
  X2 = cd4$age,
  X3 = cd4$packs,
  X4 = cd4$partners,
  X5 = cd4$drugs,
  X6 = cd4$cesd,
  Y_full = cd4$count
)
dat0 <- dat0[order(dat0$id, dat0$X1), ]
dat0$delta <- 1L; dat0$Y <- dat0$Y_full
x_cols <- sort(grep("^X\\d+$", names(dat0), value = TRUE))

yr <- range(dat0$Y_full, finite = TRUE)
pad <- 0.10 * diff(yr); if (!is.finite(pad) || pad <= 0) pad <- 50
y_grid <- seq(max(0, yr[1] - pad), yr[2] + pad, length.out = 201)
```

#### Quick check: observed trajectories
```r
tr0 <- setNames(dat0[, c("id","X1","Y_full")], c("pid","X1","Y"))
plot(range(tr0$X1), range(tr0$Y),
     type = "n", xlab = "Time (X1)", ylab = "CD4 count",
     main = sprintf("Training trajectories (%d subjects)", length(unique(tr0$pid))))
for (pp in unique(tr0$pid)) {
  ii <- which(tr0$pid == pp)
  if (length(ii) >= 2) lines(tr0$X1[ii], tr0$Y[ii], col = "grey70", lwd = 1)
}
legend("topright", bty = "n", legend = "Training trajectories", col = "grey70", lwd = 1)
```
<p align="center">
  <img src="../figures/observed.png" width="500">
</p>

#### Case A: plot band vs time for subjects
```r
for (pid in c(54, 92, 186)) {
  if (!any(dat0$id == pid)) next
  ## leave one out
  dat_test  <- dat0[dat0$id == pid, ]
  dat_train <- dat0[dat0$id != pid, ]
  x_test <- data.matrix(dat_test[, x_cols, drop = FALSE])

  res <- hcp_conformal_region(
    dat = dat_train,
    id_col = "id",
    y_col = "Y",
    delta_col = "delta",
    x_cols = x_cols,
    x_test = x_test,
    y_grid = y_grid,
    alpha = 0.1, S = 1, B = 1,
    dens_method = "rq",
    dens_taus = seq(0.05, 0.95, by = 0.01),
    seed = 1
  )
  pred <- data.frame(
    X1 = dat_test$X1,
    lo = pmax(as.numeric(res$lo_hi[, 1]), 0),
    hi = pmax(as.numeric(res$lo_hi[, 2]), 0),
    y_true = pmax(as.numeric(dat_test$Y_full), 0)
  )
  plot_hcp_intervals(
    pred, mode = "band", x_col = "X1",
    y_true_col = "y_true", show_true = TRUE,
main = sprintf("Case A: prediction band vs time: pid=%s (M=%d)", pid, nrow(pred))
  )
}
```
<p align="center">
  <img src="../figures/pid54.png" width="500">
</p>
<p align="center">
  <img src="../figures/pid92.png" width="500">
</p>
<p align="center">
  <img src="../figures/pid186.png" width="500">
</p>


#### Case B: plot one time point per subject for randomly selected subjects
```r
set.seed(2)
idsB <- sample(unique(dat0$id), 20)
datB <- dat0[dat0$id %in% idsB, ]
datB <- datB[!duplicated(datB$id), ]  # earliest (dat0 already ordered)
x_testB <- data.matrix(datB[, x_cols, drop = FALSE])

resB <- hcp_conformal_region(
  dat = dat0[!(dat0$id %in% idsB), ],
  id_col = "id", y_col = "Y", delta_col = "delta",
  x_cols = x_cols, x_test = x_testB,
  y_grid = y_grid, alpha = 0.1, S = 1, B = 1,
  dens_method = "rq", dens_taus = seq(0.05, 0.95, by = 0.01), seed = 2
)
predB <- data.frame(pid = datB$id, X1 = datB$X1,
                    lo = pmax(as.numeric(resB$lo_hi[,1]), 0),
                    hi = pmax(as.numeric(resB$lo_hi[,2]), 0),
                    y_true = pmax(as.numeric(datB$Y_full), 0))
plot_hcp_intervals(predB, mode = "pid", pid_col = "pid", x_sort_col = "X1",
                   y_true_col = "y_true", show_true = TRUE,
                   main = "Case B: random subjects, one time point")
```
<p align="center">
  <img src="../figures/subjects.png" width="500">
</p>
