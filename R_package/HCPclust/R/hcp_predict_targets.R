#' HCP prediction wrapper for multiple measurements with optional per-patient Bonferroni
#'
#' @description
#' Wraps \code{\link{hcp_conformal_region}} to produce conformal prediction regions for
#' a collection of measurements, possibly including multiple measurements per individual.
#'
#' Based on the structure of the test dataset, the prediction mode is determined
#' automatically as follows, where \eqn{P} denotes the number of patients (clusters)
#' and \eqn{M} denotes the number of measurements per patient:
#' \itemize{
#'   \item \strong{\eqn{P=1,\, M=1}}: Predict a single patient with a single measurement.
#'   \item \strong{\eqn{P=1,\, M>1}}: Predict a single patient with multiple measurements
#'   (e.g., repeated or longitudinal measurements for the same patient). If per-patient
#'   simultaneous prediction is desired, optional per-patient Bonferroni calibration
#'   can be applied.
#'   \item \strong{\eqn{P>1,\, M=1}}: Predict multiple patients, each with a single
#'   measurement. Predictions are performed independently at the nominal level
#'   \eqn{\alpha}, without Bonferroni calibration.
#'   \item \strong{\eqn{P>1,\, M>1}}: Predict multiple patients, each with multiple
#'   measurements. When per-patient simultaneous coverage is desired, a Bonferroni
#'   correction can be applied by using an effective level \eqn{\alpha / M_p} for each
#'   measurement, yielding Bonferroni-adjusted marginal prediction regions for patient
#'   \eqn{p}.
#' }
#'
#' @note
#' When per-patient Bonferroni calibration is enabled and a patient has a large number
#' of measurements (e.g., \eqn{M_p > 10}), the effective level \eqn{\alpha / M_p} may be
#' very small, which can lead to extremely wide prediction regions (potentially spanning
#' the entire \code{y_grid}). This behavior is an inherent consequence of Bonferroni
#' adjustment and not a numerical issue.
#'
#' In longitudinal or panel studies, a cluster corresponds to a single individual
#' (subject), and within-cluster points correspond to multiple time points or repeated
#' measurements on the same individual. In this setting, the time variable \code{time}
#' can be treated as a generic covariate. In the examples below, time is represented by
#' \code{X1}.
#'
#' @param dat Training/calibration data passed to \code{\link{hcp_conformal_region}}.
#' @param test A data.frame of test measurements, where each row corresponds to a single
#'   measurement. The test data must follow one of the four clustered settings
#'   \eqn{P=1, M=1}, \eqn{P=1, M>1}, \eqn{P>1, M=1}, or \eqn{P>1, M>1}, where \eqn{P} is the
#'   number of patients (clusters) and \eqn{M} is the number of measurements per patient.
#'
#'   The data.frame must include a patient identifier specified by \code{pid_col} and
#'   all covariate columns listed in \code{x_cols}. Repeated values of \code{pid_col}
#'   indicate multiple measurements (e.g., repeated or longitudinal measurements) for
#'   the same patient.
#' @param pid_col Column in \code{test} giving the patient (cluster/subject) identifier.
#'   Default \code{"pid"}.
#' @param x_cols Covariate column names (e.g., \code{c("X1")}).
#' @param y_grid Candidate y-grid passed to \code{\link{hcp_conformal_region}}.
#' @param alpha Nominal miscoverage level in (0,1) passed to
#'   \code{\link{hcp_conformal_region}}.
#' @param bonferroni Logical; if TRUE, apply per-patient Bonferroni only when a patient
#'   has multiple test measurements (i.e., \eqn{M_p>1}). If FALSE, always use level
#'   \eqn{\alpha}.
#' @param return_region Logical; if TRUE, return the full region (subset of
#'   \code{y_grid}) for each row.
#' @param id_col,y_col,delta_col Column names in \code{dat} for patient ID, outcome, and
#'   missingness indicator.
#' @param ... Additional arguments forwarded to
#'   \code{\link{hcp_conformal_region}} (e.g., \code{S}, \code{B}, \code{combine_B},
#'   \code{combine_S}, \code{dens_method}, \code{prop_method}, \code{seed}).
#'
#' @return A list with:
#' \describe{
#'   \item{pred}{A data.frame in the same row order as \code{test}. It contains all
#'     columns of \code{test} plus the effective level \code{alpha_eff} and the
#'     prediction-band endpoints \code{lo} and \code{hi} for each measurement.}
#'   \item{region}{If \code{return_region=TRUE}, a list of length
#'     \code{nrow(test)} where each element is the subset of \code{y_grid} retained
#'     in the prediction region for the corresponding test row; otherwise
#'     \code{NULL}.}
#'   \item{meta}{A list with summary information, including the number of patients
#'     \code{P}, the per-patient measurement counts \code{M_by_pid}, and the
#'     settings \code{alpha} and \code{bonferroni}.}
#' }
#'
#' @examples
#' ## ------------------------------------------------------------
#' ## Examples illustrating the four test-data settings:
#' ## (P=1, M=1), (P=1, M>1), (P>1, M=1), and (P>1, M>1)
#' ## ------------------------------------------------------------
#' set.seed(1)
#'
#' ## training data (fixed across all cases)
#' dat_train <- generate_clustered_mar(
#'   n = 200, m = 4, d = 1,
#'   x_dist = "uniform", x_params = list(min = 0, max = 10),
#'   target_missing = 0.30,
#'   seed = 1
#' )
#'
#' y_grid <- seq(-6, 6, length.out = 201)
#'
#' ## Case 1: P=1, M=1  (one patient, one measurement)
#' test_11 <- data.frame(
#'   pid = 1,
#'   X1  = 2.5
#' )
#' out_11 <- hcp_predict_targets(
#'   dat = dat_train,
#'   test = test_11,
#'   x_cols = "X1",
#'   y_grid = y_grid,
#'   alpha = 0.1,
#'   S = 2, B = 2,
#'   seed = 1
#' )
#' out_11$pred
#'
#' ## Case 2: P=1, M>1  (one patient, multiple measurements)
#' test_1M <- data.frame(
#'   pid = 1,
#'   X1  = c(1, 3, 7, 9)
#' )
#' out_1M <- hcp_predict_targets(
#'   dat = dat_train,
#'   test = test_1M,
#'   x_cols = "X1",
#'   y_grid = y_grid,
#'   alpha = 0.1,
#'   S = 2, B = 2,
#'   seed = 1
#' )
#' out_1M$pred
#'
#' ## Case 3: P>1, M=1  (multiple patients, one measurement each)
#' test_P1 <- data.frame(
#'   pid = 1:4,
#'   X1  = c(2, 4, 6, 8)
#' )
#' out_P1 <- hcp_predict_targets(
#'   dat = dat_train,
#'   test = test_P1,
#'   x_cols = "X1",
#'   y_grid = y_grid,
#'   alpha = 0.1,
#'   S = 2, B = 2,
#'   seed = 1
#' )
#' out_P1$pred
#'
#' ## Case 4: P>1, M>1  (multiple patients, multiple measurements per patient)
#' test_PM <- data.frame(
#'   pid = c(1,1, 2,2,2, 3,3),
#'   X1  = c(1,6,  2,5,9,  3,8)
#' )
#' out_PM <- hcp_predict_targets(
#'   dat = dat_train,
#'   test = test_PM,
#'   x_cols = "X1",
#'   y_grid = y_grid,
#'   alpha = 0.1,
#'   S = 2, B = 2,
#'   seed = 1
#' )
#' out_PM$pred
#' @export
hcp_predict_targets <- function(
    dat,
    test,
    pid_col = "pid",
    x_cols,
    y_grid,
    alpha = 0.1,
    bonferroni = FALSE,
    return_region = FALSE,
    id_col = "id",
    y_col = "Y",
    delta_col = "delta",
    ...
) {
  stopifnot(is.data.frame(dat))
  stopifnot(is.data.frame(test))
  stopifnot(is.character(x_cols), length(x_cols) >= 1)
  stopifnot(all(x_cols %in% names(test)))
  stopifnot(pid_col %in% names(test))

  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0,1).")
  }
  if (!is.logical(bonferroni) || length(bonferroni) != 1L) stop("bonferroni must be TRUE/FALSE.")
  if (!is.logical(return_region) || length(return_region) != 1L) stop("return_region must be TRUE/FALSE.")

  stopifnot(is.character(id_col), length(id_col) == 1L, id_col %in% names(dat))
  stopifnot(is.character(y_col), length(y_col) == 1L, y_col %in% names(dat))
  stopifnot(is.character(delta_col), length(delta_col) == 1L, delta_col %in% names(dat))

  y_grid <- as.numeric(y_grid)
  if (length(y_grid) < 2 || anyNA(y_grid) || any(!is.finite(y_grid))) {
    stop("y_grid must be finite numeric values with length >= 2.")
  }
  if (!isTRUE(all(diff(y_grid) >= 0))) y_grid <- sort(unique(y_grid))
  if (length(y_grid) < 2) {
    stop("y_grid must contain at least two distinct finite values.")
  }
  K <- nrow(test)
  if (K < 1) stop("test must have at least one row.")

  # ---- infer P and per-patient counts M_p from test ----
  pid_key <- as.character(test[[pid_col]])
  if (anyNA(pid_key)) stop("pid_col contains NA; please provide a valid patient/cluster id for each test row.")
  tab <- table(pid_key)
  P <- length(tab)
  M_vec <- as.integer(tab[pid_key])    # per-row M_p
  any_multi <- any(tab > 1L)

  # ---- per-row alpha_eff following your rules ----
  alpha_eff <- rep(alpha, K)
  if (isTRUE(bonferroni) && isTRUE(any_multi)) {
    alpha_eff <- alpha / M_vec
  }

  if (isTRUE(bonferroni)) {
    large_mp <- unique(M_vec[M_vec > 10L])
    if (length(large_mp) > 0L) {
      message(
        "Bonferroni calibration enabled: some patients have M_p > 10 measurements.\n",
        "The effective level alpha/M_p may be very small, which can lead to very wide\n",
        "prediction regions (possibly spanning the entire y_grid). ",
        "This is an inherent consequence of Bonferroni adjustment."
      )
    }
  }
  # ---- group rows by alpha_eff to minimize repeated calls ----
  key <- sprintf("%.16g", alpha_eff)
  groups <- split(seq_len(K), key)

  lo <- rep(NA_real_, K)
  hi <- rep(NA_real_, K)
  regions <- if (isTRUE(return_region)) vector("list", K) else NULL

  for (kk in names(groups)) {
    idx <- groups[[kk]]
    a_eff <- alpha_eff[idx[1]]

    x_df <- test[idx, x_cols, drop = FALSE]
    x_test <- suppressWarnings(apply(x_df, 2, function(v) as.numeric(as.character(v))))
    x_test <- as.matrix(x_test)
    colnames(x_test) <- x_cols
    if (anyNA(x_test) || any(!is.finite(x_test))) {
      stop("test contains NA/non-finite covariate values after numeric coercion.")
    }

    res <- hcp_conformal_region(
      dat = dat,
      id_col = id_col,
      y_col = y_col,
      delta_col = delta_col,
      x_cols = x_cols,
      x_test = x_test,
      y_grid = y_grid,
      alpha = a_eff,
      ...
    )

    lo[idx] <- res$lo_hi[, "lo"]
    hi[idx] <- res$lo_hi[, "hi"]
    if (isTRUE(return_region)) regions[idx] <- res$region
  }

  out_tab <- test
  out_tab$alpha_eff <- alpha_eff
  out_tab$lo <- lo
  out_tab$hi <- hi

  list(
    pred = out_tab,
    region = regions,
    meta = list(
      P = P,
      M_by_pid = tab,
      any_multi = any_multi,
      alpha = alpha,
      bonferroni = bonferroni
    )
  )
}
