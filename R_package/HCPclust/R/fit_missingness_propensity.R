#' Fit missingness propensity model P(delta=1 | X) from pooled data
#'
#' @description
#' Fits the missingness propensity
#' \eqn{\pi(x)=\mathbb{P}(\delta=1\mid x)}
#' under a marginal missingness model using pooled observations.
#' Estimation can be carried out using logistic regression,
#' Generalized Random Forests (GRF), or
#' gradient boosting (xgboost).
#' Both continuous and discrete covariates are supported; categorical variables
#' are automatically expanded into dummy variables via \code{model.matrix()}.
#'
#' @param dat A \code{data.frame} containing \code{delta_col} and \code{x_cols}.
#'   Can be any user-supplied dataset; \code{generate_clustered_mar()} is used only in examples.
#' @param delta_col Name of missingness indicator column (1 observed, 0 missing).
#' @param x_cols Character vector of covariate column names used to predict missingness.
#' @param method One of \code{"logistic"}, \code{"grf"}, \code{"boosting"}.
#' @param eps Clipping level applied to the estimated missingness propensity
#'   \eqn{\hat\pi(x)}, truncating predictions to \eqn{[\epsilon,1-\epsilon]}.
#' @param ... Extra arguments passed to the learner:
#' \describe{
#'   \item{\code{logistic}}{passed to \code{stats::glm}.}
#'   \item{\code{grf}}{passed to \code{grf::probability_forest}.}
#'   \item{\code{boosting}}{passed to \code{xgboost::xgb.train} via \code{params=} and \code{nrounds=}.}
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{method}}{The estimation method used.}
#'   \item{\code{fit}}{The fitted missingness propensity model.}
#'   \item{\code{predict}}{
#'     A function \code{predict(x_new)} that returns the estimated missingness
#'     propensity \eqn{\hat\pi(x)=\mathbb{P}(\delta=1\mid x)} evaluated at new
#'     covariate values \code{x_new}, with predictions clipped to
#'     \eqn{[\epsilon,1-\epsilon]}.
#'   }
#' }
#'
#' @examples
#' dat <- generate_clustered_mar(
#'   n = 80, m = 4, d = 2,
#'   alpha0 = -0.4, alpha = c(-1.0, 0.8),
#'   target_missing = 0.30,
#'   seed = 1
#' )
#' x_cols <- c("X1", "X2")
#'
#' ## Logistic regression
#' fit_log <- fit_missingness_propensity(dat, "delta", x_cols, method = "logistic")
#' p_log <- fit_log$predict(dat[, x_cols, drop = FALSE])
#' head(p_log)
#' \donttest{
#' ## Compare with other methods
#' ## True propensity under the generator
#' s <- attr(dat, "alpha_shift")
#' eta <- (-0.4) + (-1.0) * dat$X1 + 0.8 * dat$X2
#' pi_true <- 1 / (1 + exp(-pmin(pmax(eta, -30), 30)))
#'
#' fit_grf <- fit_missingness_propensity(
#'   dat, "delta", x_cols,
#'   method = "grf", num.trees = 800, num.threads = 1
#' )
#' fit_xgb <- fit_missingness_propensity(
#'   dat, "delta", x_cols,
#'   method = "boosting",
#'   nrounds = 300,
#'   params = list(max_depth = 3, eta = 0.05, subsample = 0.8, colsample_bytree = 0.8),
#'   nthread = 1
#' )
#'
#' p_grf <- fit_grf$predict(dat[, x_cols, drop = FALSE])
#' p_xgb <- fit_xgb$predict(dat[, x_cols, drop = FALSE])
#'
#' op <- par(mfrow = c(1, 3))
#' plot(pi_true, p_log, pch = 16, cex = 0.5,
#'      xlab = "True pi(x)", ylab = "Estimated pi-hat(x)", main = "Logistic"); abline(0, 1, lwd = 2)
#' plot(pi_true, p_grf, pch = 16, cex = 0.5,
#'      xlab = "True pi(x)", ylab = "Estimated pi-hat(x)", main = "GRF"); abline(0, 1, lwd = 2)
#' plot(pi_true, p_xgb, pch = 16, cex = 0.5,
#'      xlab = "True pi(x)", ylab = "Estimated pi-hat(x)", main = "Boosting"); abline(0, 1, lwd = 2)
#' par(op)
#' }
#' @export
fit_missingness_propensity <- function(
    dat,
    delta_col = "delta",
    x_cols,
    method = c("logistic", "grf", "boosting"),
    eps = 1e-6,
    ...
) {
  method <- match.arg(method)
  stopifnot(is.data.frame(dat))
  stopifnot(delta_col %in% names(dat))
  stopifnot(all(x_cols %in% names(dat)))

  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0 || eps >= 0.5) {
    stop("eps must be a single finite number in (0, 0.5).")
  }

  `%||%` <- function(x, y) if (is.null(x)) y else x

  # ---- robust delta parsing to {0,1} ----
  delta_raw <- dat[[delta_col]]
  delta_chr <- as.character(delta_raw)
  delta_num <- suppressWarnings(as.numeric(delta_chr))

  delta <- rep(NA_integer_, length(delta_raw))
  idx_num <- !is.na(delta_num)
  delta[idx_num] <- as.integer(delta_num[idx_num])

  idx_chr <- is.na(delta) & !is.na(delta_chr)
  if (any(idx_chr)) {
    dch <- delta_chr[idx_chr]
    delta[idx_chr] <- ifelse(
      dch %in% c("TRUE", "T", "1"),
      1L,
      ifelse(dch %in% c("FALSE", "F", "0"), 0L, NA_integer_)
    )
  }

  if (anyNA(delta)) {
    stop("delta contains NA (unparseable or missing) after parsing. Please clean delta_col.")
  }
  if (any(!delta %in% c(0L, 1L))) stop("delta must be binary (0/1) after parsing.")
  delta <- as.integer(delta)

  # ---- build design matrix (supports discrete+continuous) ----
  x_df <- dat[, x_cols, drop = FALSE]
  for (nm in names(x_df)) {
    if (is.character(x_df[[nm]])) x_df[[nm]] <- as.factor(x_df[[nm]])
  }
  x_levels <- lapply(x_df, function(v) if (is.factor(v)) levels(v) else NULL)

  terms_obj <- stats::terms(stats::as.formula("~ ."), data = x_df)
  X <- stats::model.matrix(terms_obj, data = x_df)
  if (nrow(X) != length(delta)) stop("Internal error: X and delta dimension mismatch.")

  # Make column names safe for formula-based learners (glm)
  colnames(X) <- make.names(colnames(X), unique = TRUE)

  dots <- list(...)

  fit <- switch(
    method,
    logistic = {
      df_glm <- as.data.frame(X)
      df_glm$.delta <- delta
      if (is.null(dots$family)) dots$family <- stats::binomial()
      do.call(stats::glm, c(list(formula = .delta ~ . , data = df_glm), dots))
    },
    grf = {
      if (!requireNamespace("grf", quietly = TRUE)) {
        stop("Package 'grf' is required for method='grf'. Please install it.")
      }
      y_fac <- factor(delta, levels = c(0, 1))
      do.call(grf::probability_forest, c(list(X = X, Y = y_fac), dots))
    },
    boosting = {
      if (!requireNamespace("xgboost", quietly = TRUE)) {
        stop("Package 'xgboost' is required for method='boosting'. Please install it.")
      }
      dtrain <- xgboost::xgb.DMatrix(data = X, label = delta)
      params <- dots$params
      if (is.null(params)) params <- list()
      params$objective <- params$objective %||% "binary:logistic"
      params$eval_metric <- params$eval_metric %||% "logloss"
      nrounds <- dots$nrounds %||% 400L
      dots$params <- NULL
      dots$nrounds <- NULL
      do.call(
        xgboost::xgb.train,
        c(list(params = params, data = dtrain, nrounds = nrounds, verbose = 0), dots)
      )
    }
  )

  build_X_new <- function(x_new) {
    if (is.matrix(x_new)) x_new <- as.data.frame(x_new)
    if (!is.data.frame(x_new)) stop("x_new must be a data.frame or matrix.")
    if (!all(x_cols %in% names(x_new))) stop("x_new must contain columns: ", paste(x_cols, collapse = ", "))
    x_new <- x_new[, x_cols, drop = FALSE]

    for (nm in names(x_new)) {
      if (is.character(x_new[[nm]])) x_new[[nm]] <- as.factor(x_new[[nm]])
      if (is.factor(x_new[[nm]]) && !is.null(x_levels[[nm]])) {
        x_new[[nm]] <- factor(as.character(x_new[[nm]]), levels = x_levels[[nm]])
      }
    }

    Xn <- stats::model.matrix(terms_obj, data = x_new)

    # Match training colnames (safe names) and fill missing columns with zeros
    colnames(Xn) <- make.names(colnames(Xn), unique = TRUE)

    miss <- setdiff(colnames(X), colnames(Xn))
    if (length(miss) > 0) {
      add <- matrix(0, nrow = nrow(Xn), ncol = length(miss))
      colnames(add) <- miss
      Xn <- cbind(Xn, add)
    }
    Xn <- Xn[, colnames(X), drop = FALSE]
    Xn
  }

  predict_fun <- function(x_new) {
    Xn <- build_X_new(x_new)
    p <- switch(
      method,
      logistic = as.numeric(stats::predict(fit, newdata = as.data.frame(Xn), type = "response")),
      grf = {
        pred_obj <- stats::predict(fit, newdata = Xn)
        pr <- pred_obj$predictions
        if (is.matrix(pr)) {
          if (ncol(pr) == 1) as.numeric(pr[, 1]) else as.numeric(pr[, 2])
        } else as.numeric(pr)
      },
      boosting = {
        dtest <- xgboost::xgb.DMatrix(data = Xn)
        as.numeric(stats::predict(fit, newdata = dtest))
      }
    )
    pmax(pmin(p, 1 - eps), eps)
  }

  list(
    method = method,
    fit = fit,
    terms = terms_obj,
    x_levels = x_levels,
    x_cols = x_cols,
    predict = predict_fun
  )
}
