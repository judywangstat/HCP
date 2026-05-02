#' HCP conformal prediction region with repeated subsampling and repeated data splitting
#'
#' @description
#' Constructs a marginal conformal prediction region for a new covariate value
#' \eqn{x_{n+1}} under clustered data with missing outcomes, following the HCP framework:
#' \itemize{
#'   \item \strong{(1) Model fitting.}
#'   Fit a pooled conditional density model \eqn{\widehat\pi(y\mid x)} using
#'   \code{\link{fit_cond_density_quantile}}, together with a marginal missingness
#'   propensity model \eqn{\widehat p(x)=\mathbb{P}(\delta=1\mid x)} using
#'   \code{\link{fit_missingness_propensity}}, both estimated on a subject-level
#'   training split.
#'
#'   \item \strong{(2) Subsampled calibration.}
#'   Repeatedly construct calibration sets by randomly drawing one observation
#'   per subject from the calibration split.
#'
#'   \item \strong{(3) Weighted conformal scoring.}
#'   Compute weighted conformal \eqn{p}-values over a candidate grid using the
#'   nonconformity score \eqn{R(x,y)=-\widehat\pi(y\mid x)} and inverse-propensity
#'   weights \eqn{w(x)=1/\widehat p(x)} under a MAR assumption.
#'
#'   \item \strong{(4) Aggregation.}
#'   Aggregate dependent \eqn{p}-values across subsamples (B) and data splits (S)
#'   using either the Cauchy combination test (CCT/ACAT) or the arithmetic mean.
#' }
#'
#' The prediction region is returned as a subset of the supplied grid:
#' \deqn{\widehat C(x_{n+1};\alpha)=\{y\in\mathcal Y:\ p_{\text{final}}(y)>\alpha\}.}
#'
#' @param dat A data.frame containing clustered observations. Must include \code{id_col}, \code{y_col},
#'   \code{delta_col}, and all columns in \code{x_cols}.
#' @param id_col Subject/cluster identifier column name.
#' @param y_col Outcome column name.
#' @param delta_col Missingness indicator column name (1 observed, 0 missing).
#' @param x_cols Covariate column names used for both density estimation and missingness propensity.
#' @param x_test New covariate value(s). A numeric vector (treated as one row),
#'   or a numeric matrix/data.frame with \code{nrow(x_test)=K} test points and
#'   \code{ncol(x_test)=length(x_cols)} covariates.
#' @param y_grid Numeric vector of candidate \eqn{y} values at which to evaluate conformal \eqn{p}-values.
#' @param alpha Miscoverage level in (0,1). Region keeps \eqn{y} with \eqn{p(y)>\alpha}.
#' @param train_frac Fraction of subjects assigned to training in each split.
#' @param S Number of independent subject-level splits.
#' @param B Number of subsamples per split (one observation per subject per subsample).
#' @param combine_B Combine \eqn{p}-values across B subsamples: \code{"cct"} (default) or \code{"mean"}.
#' @param combine_S Combine \eqn{p}-values across S splits: \code{"cct"} (default) or \code{"mean"}.
#' @param seed Optional seed for reproducibility.
#' @param return_details Logical; if TRUE, also return split-level p-values and split metadata.
#'
#' @param dens_method Density/quantile engine for \code{\link{fit_cond_density_quantile}}: \code{"rq"} or \code{"qrf"}.
#' @param dens_taus Quantile grid passed to \code{\link{fit_cond_density_quantile}}.
#' @param dens_h Bandwidth(s) passed to \code{\link{fit_cond_density_quantile}}.
#' @param enforce_monotone Passed to \code{\link{fit_cond_density_quantile}}.
#' @param tail_decay Passed to \code{\link{fit_cond_density_quantile}}.
#'
#' @param prop_method Missingness propensity method for \code{\link{fit_missingness_propensity}}:
#'   \code{"logistic"}, \code{"grf"}, or \code{"boosting"}.
#' @param prop_eps Clipping level for propensity predictions used by \code{\link{fit_missingness_propensity}}.
#' @param ... Extra arguments passed to \code{\link{fit_missingness_propensity}}.
#'
#' @return If \code{return_details=FALSE} (default), a list with:
#' \describe{
#'   \item{\code{region}}{Length-\code{K} list; \code{region[[k]]} is the subset of \code{y_grid} with \code{p_final[k, ] > alpha}.}
#'   \item{\code{lo_hi}}{\code{K x 2} matrix with columns \code{c("lo","hi")} giving \code{min/max} of \code{region[[k]]} (NA if empty).}
#'   \item{\code{p_final}}{\code{K x length(y_grid)} matrix of final p-values on \code{y_grid}.}
#'   \item{\code{y_grid}}{The candidate grid used.}
#' }
#' If \code{return_details=TRUE}, also includes:
#' \describe{
#'   \item{\code{p_split}}{An array with dimensions \code{c(S, K, length(y_grid))} of split-level p-values.}
#'   \item{\code{split_meta}}{Train subject IDs for each split.}
#' }
#'
#' @examples
#' dat <- generate_clustered_mar(n = 200, m = 4, d = 2, target_missing = 0.30, seed = 1)
#' y_grid <- seq(-4, 4, length.out = 200)
#' x_test <- matrix(c(0.2, -1.0), nrow = 1); colnames(x_test) <- c("X1", "X2")
#'
#' res <- hcp_conformal_region(
#'   dat, id_col = "id",
#'   y_col = "Y", delta_col = "delta",
#'   x_cols = c("X1", "X2"),
#'   x_test = x_test,
#'   y_grid = y_grid,
#'   alpha = 0.1,
#'   S = 2, B = 2,
#'   seed = 1
#' )
#'
#' ## interval endpoints on the y-grid (outer envelope)
#' c(lo = min(res$region[[1]]), hi = max(res$region[[1]]))
#'
#' @export
hcp_conformal_region <- function(
    dat,
    id_col,
    y_col = "Y",
    delta_col = "delta",
    x_cols,
    x_test,
    y_grid,
    alpha = 0.1,
    train_frac = 0.5,
    S = 5,
    B = 5,
    combine_B = c("cct", "mean"),
    combine_S = c("cct", "mean"),
    seed = NULL,
    return_details = FALSE,
    dens_method = c("rq", "qrf"),
    dens_taus = seq(0.05, 0.95, by = 0.02),
    dens_h = NULL,
    enforce_monotone = TRUE,
    tail_decay = TRUE,
    prop_method = c("logistic", "grf", "boosting"),
    prop_eps = 1e-6,
    ...
) {
  stopifnot(is.data.frame(dat))
  stopifnot(id_col %in% names(dat), y_col %in% names(dat), delta_col %in% names(dat))
  stopifnot(all(x_cols %in% names(dat)))

  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0,1).")
  }
  if (!is.numeric(train_frac) || length(train_frac) != 1L || !is.finite(train_frac) ||
      train_frac <= 0 || train_frac >= 1) {
    stop("train_frac must be a single number in (0,1).")
  }
  if (!is.numeric(S) || length(S)!=1L || !is.finite(S) || S < 1 || S != floor(S)) stop("S must be a positive integer.")
  if (!is.numeric(B) || length(B)!=1L || !is.finite(B) || B < 1 || B != floor(B)) stop("B must be a positive integer.")
  S <- as.integer(S); B <- as.integer(B)

  combine_B <- match.arg(combine_B)
  combine_S <- match.arg(combine_S)
  dens_method <- match.arg(dens_method)
  prop_method <- match.arg(prop_method)

  y_grid <- as.numeric(y_grid)
  if (anyNA(y_grid) || any(!is.finite(y_grid))) stop("y_grid must be finite numeric values.")
  if (!isTRUE(all(diff(y_grid) >= 0))) y_grid <- sort(unique(y_grid))
  Ny <- length(y_grid)

  # x_test: allow K x p
  if (is.vector(x_test) && !is.list(x_test)) {
    x_test <- matrix(as.numeric(x_test), nrow = 1)
  } else if (is.data.frame(x_test)) {
    x_test <- as.matrix(x_test)
  } else {
    x_test <- as.matrix(x_test)
  }
  if (ncol(x_test) != length(x_cols)) stop("x_test must have length(x_cols) columns.")
  storage.mode(x_test) <- "numeric"
  colnames(x_test) <- x_cols
  if (anyNA(x_test) || any(!is.finite(x_test))) stop("x_test contains NA/Inf after coercion.")
  K <- nrow(x_test)

  if (!is.logical(return_details) || length(return_details) != 1L) {
    stop("return_details must be TRUE/FALSE.")
  }

  if (!is.null(seed)) set.seed(seed)

  # elementwise combine across first dim of an array: dim = (M, K, Ny)
  combine_array_firstdim <- function(p_arr, method) {
    if (method == "mean") {
      return(apply(p_arr, c(2, 3), function(v) {
        m <- mean(v, na.rm = TRUE)
        if (!is.finite(m)) return(1)
        pmin(pmax(m, 0), 1)
      }))
    }
    # CCT: 0.5 - atan(mean(tan((0.5 - p)*pi)))/pi
    p_arr <- pmin(pmax(p_arr, 1e-15), 1 - 1e-15)
    t_arr <- tan((0.5 - p_arr) * pi)
    tbar <- apply(t_arr, c(2, 3), mean, na.rm = TRUE)
    bad <- !is.finite(tbar)
    if (any(bad)) {
      tbar[is.nan(tbar)] <- 0
      inf <- is.infinite(tbar)
      tbar[inf] <- sign(tbar[inf]) * .Machine$double.xmax
    }
    p <- 0.5 - atan(tbar) / pi
    pmin(pmax(p, 0), 1)
  }

  # ---- split/subsample helpers ----
  split_by_subject <- function(dat0, seed0 = NULL) {
    if (!is.null(seed0)) set.seed(seed0)
    ids <- unique(dat0[[id_col]])
    n_id <- length(ids)
    if (n_id < 2) stop("Need at least 2 unique subjects for splitting.")
    n_tr <- floor(train_frac * n_id)
    n_tr <- max(1L, min(n_tr, n_id - 1L))
    tr_ids <- sample(ids, size = n_tr, replace = FALSE)
    idx_tr <- dat0[[id_col]] %in% tr_ids
    list(train = dat0[idx_tr, , drop = FALSE],
         calib = dat0[!idx_tr, , drop = FALSE],
         train_ids = tr_ids)
  }

  subsample_one_per_subject <- function(dat_cal, seed0 = NULL) {
    if (!is.null(seed0)) set.seed(seed0)
    ids <- unique(dat_cal[[id_col]])
    pick <- lapply(ids, function(id) {
      rows <- which(dat_cal[[id_col]] == id)
      sample(rows, size = 1L)
    })
    dat_cal[unlist(pick), , drop = FALSE]
  }

  # helper: robust delta==1 detection
  is_delta1 <- function(d_raw) {
    d_chr <- toupper(trimws(as.character(d_raw)))
    d_num <- suppressWarnings(as.numeric(d_chr))
    (!is.na(d_num) & d_num == 1) | (d_chr %in% c("1","TRUE","T"))
  }
  is_delta0 <- function(d_raw) {
    d_chr <- toupper(trimws(as.character(d_raw)))
    d_num <- suppressWarnings(as.numeric(d_chr))
    (!is.na(d_num) & d_num == 0) | (d_chr %in% c("0","FALSE","F"))
  }

  # ---- main storage ----
  p_split_arr <- array(NA_real_, dim = c(S, K, Ny))
  split_meta <- if (isTRUE(return_details)) vector("list", S) else NULL

  for (s in seq_len(S)) {
    seed_s <- if (!is.null(seed)) seed + 1000L * s else NULL
    sp <- split_by_subject(dat, seed0 = seed_s)
    dat_tr <- sp$train
    dat_ca <- sp$calib

    # Step (1): fit density once per split
    dens_fit <- fit_cond_density_quantile(
      dat_tr,
      y_col = y_col, delta_col = delta_col, x_cols = x_cols,
      taus = dens_taus, h = dens_h,
      method = dens_method,
      enforce_monotone = enforce_monotone,
      tail_decay = tail_decay
    )

    # AUTO: if fully observed in training split, skip propensity fitting and set weights=1
    dtr <- dat_tr[[delta_col]]
    fully_observed_split <- all(is_delta1(dtr)) && !any(is_delta0(dtr))

    prop_fit <- NULL
    if (!fully_observed_split) {
      prop_fit <- fit_missingness_propensity(
        dat_tr,
        delta_col = delta_col,
        x_cols = x_cols,
        method = prop_method,
        eps = prop_eps,
        ...
      )
    }

    # ---- test-side weights ----
    if (fully_observed_split) {
      w_test_vec <- rep(1, K)
    } else {
      x_test_df <- as.data.frame(x_test); colnames(x_test_df) <- x_cols
      pt <- as.numeric(prop_fit$predict(x_test_df))
      if (any(!is.finite(pt))) stop("predict returned non-finite propensity values.")
      p_test <- pmax(pmin(pt, 1 - prop_eps), prop_eps)
      w_test_vec <- 1 / p_test
    }

    # ---- test-side nonconformity scores R_test_mat ----
    X_big <- x_test[rep(seq_len(K), each = Ny), , drop = FALSE]
    y_big <- rep(y_grid, times = K)
    X_big_df <- as.data.frame(X_big); colnames(X_big_df) <- x_cols
    f_big <- dens_fit$predict_density(X_big_df, y_big)
    f_big <- as.numeric(f_big)
    if (any(!is.finite(f_big))) stop("predict_density returned non-finite values.")
    R_test_mat <- matrix(-f_big, nrow = K, ncol = Ny, byrow = TRUE)

    # ---- calibration subsamples cache ----
    sub_cache <- vector("list", B)
    for (b in seq_len(B)) {
      seed_b <- if (!is.null(seed_s)) seed_s + 10000L + b else NULL
      dat_sub <- subsample_one_per_subject(dat_ca, seed0 = seed_b)

      d_raw <- dat_sub[[delta_col]]
      delta1 <- is_delta1(d_raw)
      d_obs <- dat_sub[delta1 & !is.na(dat_sub[[y_col]]), , drop = FALSE]
      if (nrow(d_obs) == 0) {
        sub_cache[[b]] <- list(empty = TRUE)
        next
      }

      Xc_df <- as.data.frame(as.matrix(d_obs[, x_cols, drop = FALSE])); colnames(Xc_df) <- x_cols
      Yc <- as.numeric(d_obs[[y_col]])

      R_c <- -as.numeric(dens_fit$predict_density(Xc_df, Yc))

      if (fully_observed_split) {
        w_c <- rep(1, nrow(Xc_df))
      } else {
        pc <- as.numeric(prop_fit$predict(Xc_df))
        if (any(!is.finite(pc))) stop("predict returned non-finite propensity values in calibration.")
        p_c <- pmax(pmin(pc, 1 - prop_eps), prop_eps)
        w_c <- 1 / p_c
      }

      ord <- order(R_c)
      R_sorted <- R_c[ord]
      w_sorted <- w_c[ord]
      W_ge <- rev(cumsum(rev(w_sorted)))
      sum_w <- sum(w_c)

      sub_cache[[b]] <- list(
        empty = FALSE,
        R_sorted = R_sorted,
        W_ge = W_ge,
        sum_w = sum_w
      )
    }

    # ---- compute p-values ----
    p_bky <- array(NA_real_, dim = c(B, K, Ny))
    for (b in seq_len(B)) {
      cb <- sub_cache[[b]]
      if (isTRUE(cb$empty)) {
        p_bky[b, , ] <- NA_real_
        next
      }
      R_sorted <- cb$R_sorted
      W_ge <- cb$W_ge
      sum_w <- cb$sum_w

      for (k in seq_len(K)) {
        denom <- sum_w + w_test_vec[k]
        rt <- R_test_mat[k, ]
        j_ge <- findInterval(rt, R_sorted, left.open = TRUE) + 1L
        w_ge <- ifelse(j_ge <= length(W_ge), W_ge[j_ge], 0)
        p_bky[b, k, ] <- (w_ge + w_test_vec[k]) / denom
      }
    }

    p_split_arr[s, , ] <- combine_array_firstdim(p_bky, method = combine_B)

    if (isTRUE(return_details)) split_meta[[s]] <- list(train_ids = sp$train_ids,
                                                        fully_observed = fully_observed_split)
  }

  # ---- combine over S ----
  p_final_kny <- combine_array_firstdim(p_split_arr, method = combine_S)

  regions <- vector("list", K)
  lohi <- matrix(NA_real_, nrow = K, ncol = 2)
  colnames(lohi) <- c("lo", "hi")
  for (k in seq_len(K)) {
    reg <- y_grid[p_final_kny[k, ] > alpha]
    regions[[k]] <- reg
    if (length(reg) == 0) lohi[k, ] <- c(NA_real_, NA_real_) else lohi[k, ] <- c(min(reg), max(reg))
  }

  out <- list(
    region = regions,
    lo_hi = lohi,
    p_final = p_final_kny,
    y_grid = y_grid
  )
  if (isTRUE(return_details)) {
    out$p_split <- p_split_arr
    out$split_meta <- split_meta
  }
  out
}
