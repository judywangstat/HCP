#' Simulate clustered continuous outcomes with covariate-dependent MAR missingness
#'
#' @description
#' Simulates clustered data \eqn{\{(X_{i,j},Y_{i,j},\delta_{i,j})\}} under a hierarchical
#' subject-level model with covariate-dependent Missing at Random (MAR) missingness:
#' \eqn{\delta \perp Y \mid X}. Covariates \eqn{X_{i,j}} are fully observed, while outcomes
#' \eqn{Y_{i,j}} may be missing.
#'
#' Data are generated according to the following mechanisms:
#' \itemize{
#'   \item \strong{Between-subject level}: subject random intercepts \eqn{b_i\sim N(0,\sigma_b^2)}
#'   induce within-cluster dependence, corresponding to latent subject-specific laws \eqn{P_i}.
#'
#'   \item \strong{Outcomes}: for each measurement \eqn{j=1,\ldots,m_i},
#'   \deqn{
#'   Y_{i,j} = X_{i,j}^\top \beta + b_i + \varepsilon_{i,j},
#'   }
#'   where, for each subject i, the within-cluster errors
#'   \eqn{\{\varepsilon_{i,j}\}_{j=1}^{m_i}} are mutually independent with
#'   \eqn{\varepsilon_{i,j}\sim N(0,\sigma_\varepsilon^2)} when \code{rho = 0}.
#'   When \code{rho != 0}, they follow a stationary first-order autoregressive process
#'   (AR(1)) within the cluster:
#'   \deqn{
#'   \varepsilon_{i,j} = \rho\,\varepsilon_{i,j-1} + \eta_{i,j}, \quad
#'   \eta_{i,j}\sim N\!\left(0,\sigma_\varepsilon^2(1-\rho^2)\right),
#'   }
#'   which implies \eqn{\mathrm{Var}(\varepsilon_{i,j})=\sigma_\varepsilon^2} and
#'   \eqn{\mathrm{Cov}(\varepsilon_{i,j},\varepsilon_{i,j+k})
#'   = \sigma_\varepsilon^2\rho^{|k|}} for all k.
#'
#'   \item \strong{MAR missingness}: outcomes are observed with probability
#'   \deqn{
#'   \Pr(\delta_{i,j}=1\mid X_{i,j}) = \mathrm{logit}^{-1}(\alpha_0+\alpha^\top X_{i,j}),
#'   }
#'   which depends only on covariates, ensuring \eqn{\delta \perp Y \mid X}.
#'   If \code{target_missing} is provided, the intercept \eqn{\alpha_0} is automatically
#'   calibrated (via a deterministic root-finding procedure on the expected missing proportion)
#'   so that the \emph{marginal missing proportion} is close to \code{target_missing}.
#' }
#'
#' @param n Number of clusters (subjects).
#' @param m Cluster size. Either a single positive integer (common \eqn{m_i=m}) or
#'   an integer vector of length \code{n} specifying \eqn{m_i} for each subject.
#' @param d Covariate dimension.
#' @param beta Population regression coefficients for \eqn{Y\mid X} (length \code{d}).
#'   If \code{NULL}, defaults to \code{seq(0.5, 0.5 + 0.1*(d-1), by=0.1)}.
#' @param sigma_b SD of subject random intercept \eqn{b_i}.
#' @param sigma_eps Marginal SD of within-subject errors \eqn{\varepsilon_{i,j}}.
#' @param rho AR(1) correlation parameter within cluster for \eqn{\varepsilon_{i,j}}.
#' @param hetero_gamma Optional heteroskedasticity parameter; a value of 0 yields the
#'   standard homoskedastic model, while nonzero values induce covariate-dependent
#'   error variance through the first covariate \eqn{X_1}.
#' @param x_dist Distribution for covariates: \code{"normal"}, \code{"bernoulli"}, or \code{"uniform"}.
#' @param x_params Optional list of distribution parameters for \code{x_dist}.
#' @param alpha0 Missingness intercept \eqn{\alpha_0}. If \code{target_missing} is not \code{NULL},
#'   the effective intercept becomes \eqn{\alpha_0 + s}, where \eqn{s} is a calibrated shift.
#' @param alpha Missingness slopes (length \code{d}). If \code{NULL}, defaults to zeros.
#' @param target_missing Target \emph{marginal missing proportion} defined as the empirical
#'   average of the fitted missing probabilities \eqn{1-\pi(X_{i,j})} over all observations,
#'   where \eqn{\pi(x)=\Pr(\delta=1\mid X=x)}.
#'   If \code{NULL}, no calibration.
#' @param seed Optional RNG seed.
#'
#' @return A \code{data.frame} in long format with one row per measurement:
#' \describe{
#'   \item{id}{Cluster index.}
#'   \item{j}{Within-cluster index.}
#'   \item{Y}{Observed outcome; \code{NA} if missing.}
#'   \item{Y_full}{Latent complete outcome.}
#'   \item{delta}{Observation indicator (1 observed, 0 missing).}
#'   \item{X1..Xd}{Covariates.}
#' }
#' Attributes:
#' \describe{
#'   \item{\code{m_i}}{Integer vector of cluster sizes \eqn{(m_1,\ldots,m_n)}.}
#'   \item{\code{target_missing}}{Target marginal missing proportion used for calibration,
#'     defined as the empirical average of missing probabilities over all observations.}
#'   \item{\code{alpha_shift}}{Calibrated global intercept shift \eqn{s} added to the missingness linear predictor
#'     \eqn{\alpha_0 + s + \alpha^\top X_{i,j}} (present only when \code{target_missing} is provided).}
#'   \item{\code{missing_rate}}{Sample missing rate \eqn{N^{-1}\sum I(\delta_{i,j}=0)}.
#'     This may deviate from \code{target_missing} due to Bernoulli sampling variability.}
#' }
#'
#' @examples
#' dat <- generate_clustered_mar(
#'   n = 200, m = 5, d = 2,
#'   alpha0 = -0.2, alpha = c(-1.0, 0.0),
#'   target_missing = 0.30,
#'   seed = 1
#' )
#' mean(dat$delta == 0)      # ~0.30
#' attr(dat, "alpha_shift")  # calibrated shift
#'
#' @export
generate_clustered_mar <- function(
    n,
    m = 4L,
    d = 2L,
    beta = NULL,
    sigma_b = 0.7,
    sigma_eps = 1.0,
    rho = 0,
    hetero_gamma = 0,
    x_dist = c("normal", "bernoulli", "uniform"),
    x_params = NULL,
    alpha0 = -0.2,
    alpha = NULL,
    target_missing = NULL,
    seed = NULL
) {
  x_dist <- match.arg(x_dist)
  if (!is.null(seed)) set.seed(seed)

  stopifnot(is.numeric(n), length(n) == 1, n > 0, n == as.integer(n))
  stopifnot(is.numeric(d), length(d) == 1, d >= 1, d == as.integer(d))
  stopifnot(is.numeric(alpha0), length(alpha0) == 1)
  stopifnot(is.numeric(sigma_b), length(sigma_b) == 1, sigma_b >= 0)
  stopifnot(is.numeric(sigma_eps), length(sigma_eps) == 1, sigma_eps > 0)
  stopifnot(is.numeric(rho), length(rho) == 1, abs(rho) < 1)
  stopifnot(is.numeric(hetero_gamma), length(hetero_gamma) == 1, is.finite(hetero_gamma))

  if (!is.null(target_missing)) {
    stopifnot(is.numeric(target_missing), length(target_missing) == 1,
              target_missing > 0, target_missing < 1)
  }

  # cluster sizes
  if (length(m) == 1) {
    stopifnot(is.numeric(m), length(m) == 1, m >= 1, m == as.integer(m))
    m_i <- rep(as.integer(m), n)
  } else {
    stopifnot(is.numeric(m), length(m) == n, all(m >= 1), all(m == as.integer(m)))
    m_i <- as.integer(m)
  }
  N <- sum(m_i)

  # defaults
  if (is.null(beta)) beta <- seq(0.5, 0.5 + 0.1 * (d - 1), by = 0.1)
  if (length(beta) != d) stop("beta must have length d.")
  beta <- as.numeric(beta)
  stopifnot(is.numeric(beta), all(is.finite(beta)))

  if (is.null(alpha)) alpha <- rep(0, d)
  if (length(alpha) != d) stop("alpha must have length d.")
  alpha <- as.numeric(alpha)
  stopifnot(is.numeric(alpha), all(is.finite(alpha)))

  if (is.null(x_params)) x_params <- list()
  stopifnot(is.list(x_params))

  # subject random intercept
  b <- if (sigma_b > 0) stats::rnorm(n, 0, sigma_b) else rep(0, n)

  # long indices
  id <- rep(seq_len(n), times = m_i)
  j  <- unlist(lapply(m_i, function(mi) seq_len(mi)), use.names = FALSE)

  # covariates X
  gen_X <- function(nn, dd) {
    if (x_dist == "normal") {
      mu <- if (is.null(x_params$mean)) 0 else x_params$mean
      sd <- if (is.null(x_params$sd)) 1 else x_params$sd
      stopifnot(is.numeric(mu), length(mu) == 1)
      stopifnot(is.numeric(sd), length(sd) == 1, sd > 0)
      matrix(stats::rnorm(nn * dd, mean = mu, sd = sd), nrow = nn, ncol = dd)
    } else if (x_dist == "bernoulli") {
      pr <- if (is.null(x_params$prob)) 0.5 else x_params$prob
      stopifnot(is.numeric(pr), length(pr) == 1, pr >= 0, pr <= 1)
      matrix(stats::rbinom(nn * dd, size = 1, prob = pr), nrow = nn, ncol = dd)
    } else {
      x_min <- if (is.null(x_params$min)) -1 else x_params$min
      x_max <- if (is.null(x_params$max))  1 else x_params$max
      stopifnot(is.numeric(x_min), length(x_min) == 1)
      stopifnot(is.numeric(x_max), length(x_max) == 1)
      stopifnot(x_min < x_max)
      matrix(stats::runif(nn * dd, min = x_min, max = x_max), nrow = nn, ncol = dd)
    }
  }
  X <- gen_X(N, d)

  # heteroskedastic scale based on X1 (defaults to 1 when hetero_gamma=0)
  x1 <- X[, 1]
  x01 <- (x1 - min(x1)) / (max(x1) - min(x1) + 1e-12)  # in [0,1]
  sigma_scale <- exp(hetero_gamma * (x01 - mean(x01))) # mean scale ~ 1

  # within-cluster errors (indep or AR1)
  eps <- numeric(N)
  idx_start <- cumsum(c(1L, m_i))[seq_len(n)]
  for (i in seq_len(n)) {
    mi <- m_i[i]
    pos <- idx_start[i]:(idx_start[i] + mi - 1L)
    if (abs(rho) < 1e-12) {
      eps[pos] <- stats::rnorm(mi, 0, sigma_eps * sigma_scale[pos])
    } else {
      e <- numeric(mi)
      e[1] <- stats::rnorm(1, 0, sigma_eps * sigma_scale[pos[1]])
      if (mi >= 2) {
        for (t in 2:mi) {
          sd_eta_t <- sigma_eps * sigma_scale[pos[t]] * sqrt(1 - rho^2)
          e[t] <- rho * e[t - 1] + stats::rnorm(1, 0, sd_eta_t)
        }
      }
      eps[pos] <- e
    }
  }

  # outcome
  Y_full <- as.vector(X %*% beta) + b[id] + eps

  # missingness linear predictor (WITHOUT extra calibrated shift yet)
  lin_miss <- alpha0 + as.vector(X %*% alpha)

  # fix uniforms once for deterministic calibration
  expected_miss_given_shift <- function(s) {
    eta <- lin_miss + s
    eta <- pmin(pmax(eta, -30), 30)
    p_obs <- 1 / (1 + exp(-eta))
    1 - mean(p_obs)   # continuous / smooth target
  }

  alpha_shift <- 0
  if (!is.null(target_missing)) {
    obj <- function(s) expected_miss_given_shift(s) - target_missing
    lo <- -30
    hi <-  30
    f_lo <- obj(lo)
    f_hi <- obj(hi)

    if (f_lo * f_hi > 0) {
      alpha_shift <- if (abs(f_lo) < abs(f_hi)) lo else hi
      warning(sprintf(
        "Failed to bracket root for missingness intercept shift. Using shift=%.3f; achieved missing=%.3f (target=%.3f).",
        alpha_shift, expected_miss_given_shift(alpha_shift), target_missing
      ))
    } else {
      alpha_shift <- stats::uniroot(obj, lower = lo, upper = hi, tol = 1e-4)$root
    }
  }

  eta_final <- lin_miss + alpha_shift
  eta_final <- pmin(pmax(eta_final, -30), 30)
  p_obs <- 1 / (1 + exp(-eta_final))
  U_miss <- stats::runif(N)
  delta <- as.integer(U_miss <= p_obs)

  Y <- Y_full
  Y[delta == 0L] <- NA_real_

  out <- data.frame(
    id = id,
    j = j,
    Y = Y,
    Y_full = Y_full,
    delta = as.integer(delta)
  )
  for (k in seq_len(d)) out[[paste0("X", k)]] <- X[, k]

  attr(out, "m_i") <- m_i

  if (!is.null(target_missing)) {
    attr(out, "target_missing") <- target_missing
    attr(out, "alpha_shift") <- alpha_shift
    attr(out, "missing_rate") <- mean(delta == 0L)
  }

  out
}
