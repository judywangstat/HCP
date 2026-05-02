#' Estimate conditional density pi(y|x) via quantile process + quotient estimator
#'
#' @description
#' Fits a conditional quantile function \eqn{\widehat Q_Y(\tau\mid x)} using pooled observed data
#' (working-independence), and estimates the conditional density through the quotient estimator
#' along the quantile curve:
#' \deqn{
#' \widehat\pi\{\widehat Q(\tau\mid x)\mid x\}
#' = \frac{2h(\tau)}{\widehat Q(\tau+h(\tau)\mid x)-\widehat Q(\tau-h(\tau)\mid x)}.
#' }
#' For numerical stability, the quantile curve can be monotone-adjusted (isotonic regression),
#' and tail decay extrapolation can be used before interpolation to \eqn{\pi(y\mid x)}.
#'
#' @param dat data.frame in long format, containing outcome, missingness indicator, and covariates.
#' @param y_col name of outcome column (observed Y, may contain NA).
#' @param delta_col name of missingness indicator (1 observed, 0 missing).
#' @param x_cols character vector of covariate column names (include time if desired).
#' @param taus grid of quantile levels in (0,1) at which the quantile process is evaluated.
#' @param h Bandwidth(s) for quotient. Either a scalar or a numeric vector of length \code{length(taus)}.
#'   If \code{NULL}, a tau-specific bandwidth vector \eqn{h(\tau)} is computed via
#'   \code{quantreg::bandwidth.rq}, and automatically shrunk near the
#'   boundaries to ensure \eqn{\tau\pm h(\tau)\in(0,1)}.
#' @param method quantile engine: \code{"rq"} (linear quantile regression) or \code{"qrf"} (quantile random forest).
#' @param enforce_monotone logical; if TRUE, apply isotonic regression to the predicted quantile curve in \eqn{\tau}
#'   for each \eqn{x} to reduce quantile crossing.
#' @param tail_decay logical; if TRUE, add extra tail points with geometric decay before interpolation.
#' @param num_extra_points number of extra tail points on each side when \code{tail_decay=TRUE}.
#' @param decay_factor decay factor in (0,1) for tail densities when \code{tail_decay=TRUE}.
#' @param dens_floor lower bound for density to avoid numerical issues.
#' @param eps small stabilizer for denominator \code{pmax(Qplus-Qminus, eps)}.
#' @param gap_min minimum spacing for tail extrapolation points.
#' @param seed optional seed.
#' @param ... extra arguments passed to the underlying quantile engine:
#'   \describe{
#'     \item{\code{rq}}{passed to \code{quantreg::rq.fit}, e.g. \code{rq_method="br"}.}
#'     \item{\code{qrf}}{passed to \code{quantregForest::quantregForest}, e.g. \code{ntree=500}.}
#'   }
#'
#' @return A list containing fitted objects and prediction functions:
#' \describe{
#'   \item{\code{predict_Q(x_new, taus_use)}}{
#'     Returns the estimated conditional quantiles
#'     \deqn{
#'         \widehat Q_Y(\tau \mid x)
#'     }
#'     for \eqn{\tau \in (0,1)} specified by \code{taus_use},
#'     evaluated at new covariate values \code{x_new}.
#'     The output is a numeric matrix with one row per covariate vector \eqn{x}
#'     and one column per quantile level \eqn{\tau}.
#'   }
#'
#'   \item{\code{predict_density(x_new, y_new)}}{
#'     Returns the estimated conditional density
#'     \deqn{
#'       \widehat \pi(y \mid x),
#'     }
#'     evaluated at specified (x,y) pairs.
#'     The inputs \code{x_new} and \code{y_new} are paired row-wise, so that the
#'     r-th row of \code{x_new} is evaluated at \code{y_new[r]}.
#'   }
#' }
#'
#' @examples
#' ## ------------------------------------------------------------
#' ## Case A: Conditional density evaluated at a single point (x, y)
#' ## ------------------------------------------------------------
#' ## This illustrates the most basic usage: estimating pi(y | x)
#' ## at one covariate value x and one response value y.
#'
#' dat <- generate_clustered_mar(
#'   n = 200, m = 4, d = 2,
#'   target_missing = 0.3, seed = 1
#' )
#' fit <- fit_cond_density_quantile(
#'   dat,
#'   y_col = "Y", delta_col = "delta",
#'   x_cols = c("X1", "X2"),
#'   taus = seq(0.05, 0.95, by = 0.02),
#'   method = "rq",
#'   seed = 1
#' )
#' ## a single covariate value x
#' x1 <- matrix(c(0.2, -1.0), nrow = 1)
#' colnames(x1) <- c("X1", "X2")
#' ## estimate pi(y | x) at y = 0.5
#' fit$predict_density(x1, y_new = 0.5)
#'
#'
#' ## ------------------------------------------------------------
#' ## Case B: Conditional density as a function of y (density curve)
#' ## ------------------------------------------------------------
#' ## Here we fix x and evaluate pi(y | x) over a grid of y values,
#' ## which produces an estimated conditional density curve.
#'
#' y_grid <- seq(-3, 3, length.out = 201)
#' ## reuse the same x by repeating it to match the y-grid
#' x_rep <- x1[rep(1, length(y_grid)), , drop = FALSE]
#' f_grid <- fit$predict_density(x_rep, y_grid)
#'
#' ## ------------------------------------------------------------
#' ## True conditional density under the data generator
#' ## ------------------------------------------------------------
#' ## Data are generated as:
#' ##   Y = X^T beta + b + eps,
#' ##   b ~ N(0, sigma_b^2),  eps ~ N(0, sigma_eps^2)
#' ## Hence the marginal conditional density is:
#' ##   Y | X = x ~ N(x^T beta, sigma_b^2 + sigma_eps^2)
#'
#' beta_true <- c(0.5, 0.6)
#' sigma_b_true <- 0.7
#' sigma_eps_true <- 1.0
#' mu_true <- drop(x1 %*% beta_true)
#' sd_true <- sqrt(sigma_b_true^2 + sigma_eps_true^2)
#' f_true <- stats::dnorm(y_grid, mean = mu_true, sd = sd_true)
#'
#'
#' ## ------------------------------------------------------------
#' ## Visualization: estimated vs true conditional density
#' ## (use smooth.spline on log-density for a smoother display)
#' ## ------------------------------------------------------------
#'
#' ## smooth the estimated curve for visualization
#' ok <- is.finite(f_grid) & (f_grid > 0)
#' sp <- stats::smooth.spline(y_grid[ok], log(f_grid[ok]), spar = 0.85)
#' f_smooth <- exp(stats::predict(sp, y_grid)$y)
#'
#' ymax <- max(c(f_smooth, f_true), na.rm = TRUE)
#' plot(
#'   y_grid, f_smooth,
#'   type = "l", lwd = 2,
#'   xlab = "y",
#'   ylab = expression(hat(pi)(y ~ "|" ~ x)),
#'   ylim = c(0, 1.2 * ymax),
#'   main = "Conditional density at a fixed x: estimated vs true"
#' )
#' grid(col = "gray85", lty = 1)
#' lines(y_grid, f_true, lwd = 2, lty = 2)
#' legend(
#'   "topright",
#'   legend = c("Estimated (smoothed)", "True (generator)"),
#'   lty = c(1, 2), lwd = c(2, 2), bty = "n"
#' )
#'
#' @export
fit_cond_density_quantile <- function(
    dat,
    y_col = "Y",
    delta_col = "delta",
    x_cols,
    taus = seq(0.05, 0.95, by = 0.01),
    h = NULL,
    method = c("rq", "qrf"),
    enforce_monotone = TRUE,
    tail_decay = TRUE,
    num_extra_points = 10L,
    decay_factor = 0.8,
    dens_floor = 1e-10,
    eps = 1e-8,
    gap_min = 1e-2,
    seed = NULL,
    ...
) {
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)

  stopifnot(is.data.frame(dat))
  stopifnot(all(c(y_col, delta_col, x_cols) %in% names(dat)))

  taus <- sort(unique(as.numeric(taus)))
  if (anyNA(taus)) stop("taus contains NA (non-numeric values).")
  if (!is.numeric(taus) || any(taus <= 0) || any(taus >= 1)) stop("taus must be numeric values in (0,1).")
  if (length(taus) < 5) stop("taus must contain at least 5 distinct values in (0,1).")

  stopifnot(is.logical(enforce_monotone), length(enforce_monotone) == 1)
  stopifnot(is.logical(tail_decay), length(tail_decay) == 1)
  stopifnot(is.numeric(decay_factor), decay_factor > 0, decay_factor < 1)
  stopifnot(is.numeric(dens_floor), dens_floor > 0)
  stopifnot(is.numeric(eps), eps > 0)
  stopifnot(is.numeric(gap_min), gap_min > 0)
  num_extra_points <- as.integer(num_extra_points)
  if (num_extra_points < 0) stop("num_extra_points must be >= 0.")

  # --- observed Y only (MAR => OK) ---
  delta_raw <- dat[[delta_col]]
  delta_chr <- as.character(delta_raw)
  delta_num <- suppressWarnings(as.numeric(delta_chr))
  delta_ok  <- !is.na(delta_raw) & (delta_num == 1 | delta_chr %in% c("1","TRUE","T"))
  dobs <- dat[delta_ok & !is.na(dat[[y_col]]), , drop = FALSE]
  n_obs <- nrow(dobs)
  if (n_obs < 20) stop("Too few observed outcomes after filtering delta==1.")
  if (n_obs < 50) warning("Very small number of observed outcomes; density estimate may be unstable.")

  x_df <- dobs[, x_cols, drop = FALSE]
  if (!all(vapply(x_df, is.numeric, logical(1)))) {
    bad <- names(x_df)[!vapply(x_df, is.numeric, logical(1))]
    stop("Non-numeric covariates in x_cols: ", paste(bad, collapse = ", "),
         ". Please convert to numeric or use model.matrix().")
  }
  X_obs <- as.matrix(x_df)
  y_obs <- as.numeric(dobs[[y_col]])
  if (ncol(X_obs) != length(x_cols)) stop("x_cols selection error.")
  if (any(!is.finite(y_obs))) stop("Non-finite y values found among observed outcomes.")

  # --- choose tau-specific bandwidth h(tau) and ensure tau +/- h in (0,1) ---
  if (is.null(h)) {
    if (!requireNamespace("quantreg", quietly = TRUE)) {
      stop("When h=NULL, package 'quantreg' is required to compute h(tau) via quantreg::bandwidth.rq().")
    }
    h_vec <- tryCatch(
      quantreg::bandwidth.rq(taus, n_obs, hs = TRUE, alpha = 0.1),
      error = function(e) rep(0.25 * n_obs^(-1/5), length(taus))
    )
  } else if (length(h) == 1L) {
    h_vec <- rep(as.numeric(h), length(taus))
  } else {
    h_vec <- as.numeric(h)
    if (length(h_vec) != length(taus)) stop("If h is a vector, it must have length equal to length(taus).")
  }

  h_vec <- as.numeric(h_vec)
  # safety fallback if bandwidth.rq returns bad values
  h_vec[!is.finite(h_vec) | h_vec <= 0] <- 0.5 * n_obs^(-1/5)

  taus_hi <- taus + h_vec
  taus_lo <- taus - h_vec

  # shrink near boundaries
  min_h <- 1e-6
  for (i in seq_along(taus)) {
    iter <- 0L
    while ((taus_lo[i] <= 0 || taus_hi[i] >= 1) && iter < 60L) {
      h_vec[i] <- h_vec[i] / 2
      taus_hi[i] <- taus[i] + h_vec[i]
      taus_lo[i] <- taus[i] - h_vec[i]
      iter <- iter + 1L
      if (h_vec[i] < min_h) break
    }
  }

  ok_tau <- (taus_lo > 0) & (taus_hi < 1) & is.finite(h_vec) & (h_vec > min_h)
  if (sum(ok_tau) < 5) stop("Too few taus remain after enforcing tau +/- h in (0,1). Use interior taus or provide smaller h.")

  taus0 <- taus[ok_tau]
  h0 <- h_vec[ok_tau]
  taus_minus <- taus_lo[ok_tau]
  taus_plus  <- taus_hi[ok_tau]
  taus_pm <- sort(unique(c(taus_minus, taus0, taus_plus)))

  # --------- engine: rq ----------
  dots <- list(...)
  rq_method <- dots$rq_method
  if (is.null(rq_method)) rq_method <- "br"

  fit_rq_engine <- function() {
    if (!requireNamespace("quantreg", quietly = TRUE)) {
      stop("Package 'quantreg' is required for method='rq'.")
    }
    X_int <- cbind(1, X_obs)
    colnames(X_int)[1] <- "(Intercept)"

    fit_one <- function(tau) {
      tryCatch(
        quantreg::rq.fit(x = X_int, y = y_obs, tau = tau, method = rq_method),
        error = function(e) NULL
      )
    }
    fits <- lapply(taus_pm, fit_one)
    bad <- which(vapply(fits, is.null, logical(1)))
    if (length(bad) > 0) {
      stop("rq.fit failed for taus: ", paste(sprintf("%.6f", taus_pm[bad]), collapse=", "),
           ". Try different rq_method or fewer taus.")
    }
    coefs <- do.call(cbind, lapply(fits, function(f) f$coefficients))
    rownames(coefs) <- colnames(X_int)
    colnames(coefs) <- sprintf("tau=%.6f", taus_pm)

    predict_Q <- function(x_new, taus_use = taus_pm) {
      x_new <- assert_numeric_x(x_new, x_cols, name = "x_new")
      if (ncol(x_new) != length(x_cols)) stop("x_new must have length(x_cols) columns.")
      Xn <- cbind(1, x_new)
      idx <- match(taus_use, taus_pm)
      if (anyNA(idx)) stop("taus_use must be a subset of fitted taus.")
      Xn %*% coefs[, idx, drop = FALSE]
    }

    list(engine = "rq", object = coefs, predict_Q = predict_Q)
  }

  # --------- engine: qrf ----------
  fit_qrf_engine <- function() {
    if (!requireNamespace("quantregForest", quietly = TRUE)) {
      stop("Package 'quantregForest' is required for method='qrf'.")
    }
    xdf <- as.data.frame(X_obs)
    colnames(xdf) <- x_cols

    dots_qrf <- dots
    dots_qrf$rq_method <- NULL  # IMPORTANT: rq-only arg, avoid unused argument() in qrf

    qrf_fit <- do.call(
      quantregForest::quantregForest,
      c(list(x = xdf, y = y_obs), dots_qrf)
    )

    predict_Q <- function(x_new, taus_use = taus_pm) {
      x_new <- assert_numeric_x(x_new, x_cols, name = "x_new")
      if (ncol(x_new) != length(x_cols)) stop("x_new must have length(x_cols) columns.")
      xdf_new <- as.data.frame(x_new)
      colnames(xdf_new) <- x_cols
      Q <- stats::predict(qrf_fit, newdata = xdf_new, what = taus_use)
      Q <- as.matrix(Q)
      colnames(Q) <- sprintf("tau=%.6f", taus_use)
      Q
    }

    list(engine = "qrf", object = qrf_fit, predict_Q = predict_Q)
  }

  engine <- switch(
    method,
    rq  = fit_rq_engine(),
    qrf = fit_qrf_engine()
  )

  # ---- helpers ----
  assert_numeric_x <- function(x, x_cols, name = "x_new") {
    if (is.data.frame(x)) {
      if (ncol(x) != length(x_cols)) stop(name, " must have length(x_cols) columns.")
      bad <- names(x)[!vapply(x, is.numeric, logical(1))]
      if (length(bad) > 0) {
        stop(name, " has non-numeric columns: ", paste(bad, collapse = ", "),
             ". Please convert to numeric or use model.matrix().")
      }
      xm <- as.matrix(x)
    } else {
      xm <- as.matrix(x)
      if (ncol(xm) != length(x_cols)) stop(name, " must have length(x_cols) columns.")
      if (!is.numeric(xm)) {
        suppressWarnings(xm2 <- apply(xm, 2, as.numeric))
        xm2 <- as.matrix(xm2)
        if (anyNA(xm2)) stop(name, " cannot be safely coerced to numeric (NA produced).")
        xm <- xm2
      }
    }
    colnames(xm) <- x_cols
    xm
  }

  iso_adjust_vec <- function(v) as.numeric(stats::isoreg(v)$yf)

  # precompute index maps for fast extraction (avoid repeated match calls)
  idx_m <- match(taus_minus, taus_pm)
  idx_0 <- match(taus0,      taus_pm)
  idx_p <- match(taus_plus,  taus_pm)
  if (anyNA(idx_m) || anyNA(idx_0) || anyNA(idx_p)) stop("Internal error: tau grid mismatch.")

  grid_at_x <- function(x_row) {
    x_row <- matrix(as.numeric(x_row), nrow = 1)
    if (anyNA(x_row)) stop("x_row contains NA after numeric coercion.")
    if (ncol(x_row) != length(x_cols)) stop("x_row has wrong length.")

    # Predict the whole quantile curve on taus_pm, then (optionally) monotone-adjust ONCE
    Qall <- as.numeric(engine$predict_Q(x_row, taus_pm)[1, ])
    if (isTRUE(enforce_monotone)) {
      Qall <- iso_adjust_vec(Qall)
    }

    Qm <- Qall[idx_m]
    Q0 <- Qall[idx_0]
    Qp <- Qall[idx_p]

    denom <- pmax(Qp - Qm, eps)
    fhat <- (2 * h0) / denom
    fhat[!is.finite(fhat)] <- NA_real_
    fhat[fhat < dens_floor] <- dens_floor

    list(q = Q0, f = fhat)
  }

  dens_at_xy_one <- function(x_row, y_val_one, rule = 2) {
    g <- grid_at_x(x_row)
    ok <- is.finite(g$q) & is.finite(g$f)
    if (sum(ok) < 5) return(NA_real_)

    xq <- g$q[ok]
    fq <- g$f[ok]

    ord <- order(xq)
    xq <- xq[ord]; fq <- fq[ord]
    keep <- !duplicated(xq)
    xq <- xq[keep]; fq <- fq[keep]
    if (length(xq) < 2) {
      xq <- c(xq[1], xq[1] + 0.1)
      fq <- c(fq[1], fq[1])
    }

    x_all <- xq
    f_all <- fq

    if (isTRUE(tail_decay) && num_extra_points > 0L) {
      dx <- diff(xq)
      gap <- stats::median(dx[is.finite(dx) & dx > 0], na.rm = TRUE)
      if (!is.finite(gap) || gap <= 0) gap <- gap_min
      gap <- max(gap, gap_min)

      extra_left  <- seq(min(xq) - num_extra_points * gap, min(xq) - gap, by = gap)
      extra_right <- seq(max(xq) + gap, max(xq) + num_extra_points * gap, by = gap)

      extra_f_left  <- fq[1] * decay_factor ^ seq(num_extra_points, 1)
      extra_f_right <- fq[length(fq)] * decay_factor ^ seq(1, num_extra_points)

      x_all <- c(extra_left, xq, extra_right)
      f_all <- c(extra_f_left, fq, extra_f_right)
      f_all[f_all < dens_floor] <- dens_floor
    }

    ord2 <- order(x_all)
    x_all <- x_all[ord2]; f_all <- f_all[ord2]
    keep2 <- !duplicated(x_all)
    x_all <- x_all[keep2]; f_all <- f_all[keep2]
    if (length(x_all) < 2) return(NA_real_)

    stats::approx(
      x = x_all, y = f_all,
      xout = as.numeric(y_val_one),
      method = "linear",
      rule = rule, ties = "ordered"
    )$y
  }

  predict_density <- function(x_new, y_new, rule = 2) {
    x_new <- assert_numeric_x(x_new, x_cols, name = "x_new")
    y_new <- as.numeric(y_new)
    if (ncol(x_new) != length(x_cols)) stop("x_new must have length(x_cols) columns.")
    if (nrow(x_new) != length(y_new)) stop("x_new and y_new must have the same number of rows/elements.")

    out <- numeric(length(y_new))
    for (r in seq_len(nrow(x_new))) {
      out[r] <- dens_at_xy_one(x_new[r, ], y_new[r], rule = rule)
    }
    out
  }

  predict_Q_user <- function(x_new, taus_use = taus0) {
    x_new <- assert_numeric_x(x_new, x_cols, name = "x_new")
    engine$predict_Q(x_new, taus_use)
  }

  list(
    method = method,
    x_cols = x_cols,
    taus = taus0,
    h = h0,
    engine = engine$engine,
    fit_object = engine$object,
    predict_Q = predict_Q_user,
    predict_density = predict_density
  )
}
