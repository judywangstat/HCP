#' Plot HCP prediction intervals (band vs covariate or intervals by patient)
#'
#' @description
#' Unified plotting function for two common visualizations of HCP prediction intervals:
#' \itemize{
#'   \item \code{mode="band"}: plot an interval band (lo/hi) versus a 1D covariate (e.g., time X1).
#'   \item \code{mode="pid"}: plot one interval per patient on the x-axis (patients optionally sorted by a covariate).
#' }
#'
#' @param df A data.frame containing prediction results. It must include the interval
#'   endpoints specified by \code{lo_col} and \code{hi_col}, and the covariate columns
#'   required by the chosen plotting mode.
#'
#' @param mode Plotting mode. Use \code{"band"} to visualize an interval band as a
#'   function of a continuous covariate, or \code{"pid"} to visualize one prediction
#'   interval per patient on the x-axis.
#'
#' @param x_col (mode = "band") Name of the covariate column used as the x-axis in the
#'   interval band plot (e.g., time or another continuous predictor).
#'
#' @param pid_col (mode = "pid") Name of the column identifying patients (or clusters).
#'   Each patient must appear exactly once in \code{df}. Default is \code{"pid"}.
#'
#' @param x_sort_col (mode = "pid") Optional covariate column used to order patients along
#'   the x-axis (e.g., \code{"X1"}). If \code{NULL}, patients are ordered by their IDs.
#'
#' @param lo_col Name of the column containing the lower endpoint of the prediction
#'   interval. Default is \code{"lo"}.
#'
#' @param hi_col Name of the column containing the upper endpoint of the prediction
#'   interval. Default is \code{"hi"}.
#'
#' @param y_true_col Optional name of a column in \code{df} containing the true outcome
#'   values. Used for overlaying truth points when \code{show_true = TRUE}.
#'
#' @param y_true Optional numeric vector of true outcome values with length equal to
#'   \code{nrow(df)}. If provided, this overrides \code{y_true_col}.
#'
#' @param show_center Logical; if \code{TRUE}, draw the midpoint of each interval
#'   (as a dashed line in \code{mode = "band"} or as points in \code{mode = "pid"}).
#'
#' @param show_true Logical; if \code{TRUE}, overlay true outcome values when available.
#'
#' @param max_patients (mode = "pid") Optional maximum number of patients to display.
#'   If specified, only the first \code{max_patients} patients after sorting are plotted.
#'
#' @param ... Additional graphical parameters passed to \code{\link{plot}}, such as
#'   \code{main}, \code{xlab}, \code{ylab}, \code{xlim}, or \code{ylim}.
#'
#' @return Invisibly returns the data.frame used for plotting:
#' \itemize{
#'   \item For \code{mode = "band"}, the input \code{df} sorted by \code{x_col}.
#'   \item For \code{mode = "pid"}, the input \code{df} sorted by \code{pid_col} or
#'   \code{x_sort_col}, if provided.
#' }
#'
#' @examples
#' ## ------------------------------------------------------------
#' ## Two common plots:
#' ## (A) one patient, multiple measurements  -> interval band vs X1
#' ## (B) multiple patients, one measurement -> intervals by patient (sorted by X1)
#' ## ------------------------------------------------------------
#' dat_train <- generate_clustered_mar(
#'   n = 200, m = 20, d = 1,
#'   x_dist = "uniform", x_params = list(min = 0, max = 10),
#'   hetero_gamma = 2.5,
#'   target_missing = 0.30,
#'   seed = 1
#' )
#' y_grid <- seq(-6, 10, length.out = 201)
#'
#' ## test data with latent truth
#' dat_test <- generate_clustered_mar(
#'   n = 100, m = 20, d = 1,
#'   x_dist = "uniform", x_params = list(min = 0, max = 10),
#'   hetero_gamma = 2.5,
#'   seed = 999
#' )
#'
#' ## ---------- Case A: P=1, M>1 (one patient, multiple measurements) ----------
#' pid <- dat_test$id[1]
#' idx <- which(dat_test$id == pid)
#' idx <- idx[order(dat_test$X1[idx])][1:10]
#' test_1M <- data.frame(pid = pid, X1 = dat_test$X1[idx], y_true = dat_test$Y_full[idx])
#'
#' out_1M <- hcp_predict_targets(
#'   dat = dat_train, test = test_1M,
#'   x_cols = "X1", y_grid = y_grid,
#'   alpha = 0.1,
#'   S = 2, B = 2,
#'   seed = 1
#' )
#' plot_hcp_intervals(
#'   out_1M$pred, mode = "band", x_col = "X1",
#'   y_true_col = "y_true", show_true = TRUE,
#'   main = "Case A: one patient, multiple time points (band vs time)"
#' )
#'
#' ## ---------- Case B: P>1, M=1 (multiple patients, one measurement each) ----------
#' ## take one measurement per patient: j==1 for the first 20 patients
#' pids <- unique(dat_test$id)[1:20]
#' test_P1 <- subset(dat_test, id %in% pids & j == 1,
#'                   select = c(id, X1, Y_full))
#' names(test_P1) <- c("pid", "X1", "y_true")
#'
#' out_P1 <- hcp_predict_targets(
#'   dat = dat_train, test = test_P1,
#'   x_cols = "X1", y_grid = y_grid,
#'   alpha = 0.1,
#'   S = 2, B = 2,
#'   seed = 1
#' )
#' plot_hcp_intervals(
#'   out_P1$pred, mode = "pid", pid_col = "pid", x_sort_col = "X1",
#'   y_true_col = "y_true", show_true = TRUE,
#'   main = "Case B: multiple patients, one time point (by patient)"
#' )
#'
#' @export
plot_hcp_intervals <- function(
    df,
    mode = c("band", "pid"),
    # shared
    lo_col = "lo",
    hi_col = "hi",
    y_true_col = NULL,
    y_true = NULL,
    show_center = TRUE,
    show_true = TRUE,
    # band mode
    x_col = NULL,
    # pid mode
    pid_col = "pid",
    x_sort_col = NULL,
    max_patients = NULL,
    ...
) {
  stopifnot(is.data.frame(df))
  mode <- match.arg(mode)

  stopifnot(lo_col %in% names(df), hi_col %in% names(df))

  lo <- as.numeric(df[[lo_col]])
  hi <- as.numeric(df[[hi_col]])

  # truth (vector overrides column)
  has_true <- FALSE
  if (!is.null(y_true)) {
    y_true <- as.numeric(y_true)
    if (length(y_true) != nrow(df)) stop("y_true length must equal nrow(df).")
    has_true <- any(is.finite(y_true))
  } else if (!is.null(y_true_col) && y_true_col %in% names(df)) {
    y_true <- as.numeric(df[[y_true_col]])
    has_true <- any(is.finite(y_true))
  } else {
    y_true <- rep(NA_real_, nrow(df))
  }

  # style
  band_col   <- grDevices::adjustcolor("steelblue3", alpha.f = 0.15)
  bound_col  <- grDevices::adjustcolor("steelblue4", alpha.f = 0.95)
  center_col <- grDevices::adjustcolor("steelblue4", alpha.f = 0.55)
  grid_col   <- grDevices::adjustcolor("gray85", alpha.f = 0.75)
  truth_col  <- grDevices::adjustcolor("firebrick3", alpha.f = 0.80)

  dots <- list(...)
  has_ylim <- "ylim" %in% names(dots)

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(4.2, 4.2, 2.3, 1.2), bty = "l")

  if (mode == "band") {
    if (is.null(x_col) || !(x_col %in% names(df))) stop("For mode='band', x_col must be provided and exist in df.")
    x <- as.numeric(df[[x_col]])
    if (anyNA(x) || any(!is.finite(x))) stop("x_col contains NA/non-finite values.")

    ord <- order(x)
    df2 <- df[ord, , drop = FALSE]
    x <- x[ord]; lo2 <- lo[ord]; hi2 <- hi[ord]
    y2 <- y_true[ord]

    ok_band <- is.finite(lo2) & is.finite(hi2)

    y_all <- c(lo2[is.finite(lo2)], hi2[is.finite(hi2)])
    if (isTRUE(show_true) && has_true) y_all <- c(y_all, y2[is.finite(y2)])
    if (length(y_all) == 0) stop("No finite values to plot.")
    yl <- range(y_all)
    pad <- 0.10 * diff(yl); if (!is.finite(pad) || pad <= 0) pad <- 1
    ylim_auto <- c(yl[1] - pad, yl[2] + pad)

    if (!("xlab" %in% names(dots))) dots$xlab <- x_col
    if (!("ylab" %in% names(dots))) dots$ylab <- "Y"

    plot_args <- c(
      list(x = x, y = hi2, type = "n"),
      if (!has_ylim) list(ylim = ylim_auto) else list(),
      dots
    )
    do.call(graphics::plot, plot_args)

    graphics::abline(h = pretty(graphics::par("usr")[3:4]), v = pretty(graphics::par("usr")[1:2]),
                     col = grid_col, lty = 3)

    if (any(ok_band)) {
      xb <- x[ok_band]; lob <- lo2[ok_band]; hib <- hi2[ok_band]
      graphics::polygon(c(xb, rev(xb)), c(lob, rev(hib)), col = band_col, border = NA)
    }
    graphics::lines(x, lo2, lwd = 1.8, col = bound_col)
    graphics::lines(x, hi2, lwd = 1.8, col = bound_col)

    if (isTRUE(show_center) && any(ok_band)) {
      center <- (lo2 + hi2) / 2
      graphics::lines(x, center, lwd = 1.2, lty = 2, col = center_col)
    }

    if (isTRUE(show_true) && has_true) {
      ok_t <- is.finite(y2)
      if (any(ok_t)) graphics::points(x[ok_t], y2[ok_t], pch = 16, cex = 0.8, col = truth_col)
    }

    # --- legend: 2) lines/points only (no fill) ---
    leg2 <- c("Bounds")
    col2 <- c(bound_col)
    lty2 <- c(1)
    lwd2 <- c(2)
    pch2 <- c(NA)

    if (isTRUE(show_center)) {
      leg2 <- c(leg2, "Center")
      col2 <- c(col2, center_col)
      lty2 <- c(lty2, 2)
      lwd2 <- c(lwd2, 1.5)
      pch2 <- c(pch2, NA)
    }

    if (isTRUE(show_true) && has_true) {
      leg2 <- c(leg2, "True y")
      col2 <- c(col2, truth_col)
      lty2 <- c(lty2, NA)
      lwd2 <- c(lwd2, NA)
      pch2 <- c(pch2, 16)
    }

    graphics::legend("topleft",
                     inset = c(0, 0.12),   # 往下挪一点，避免重叠
                     legend = leg2,
                     bty = "n",
                     col = col2,
                     lty = lty2,
                     lwd = lwd2,
                     pch = pch2,
                     pt.cex = 0.8)

    return(invisible(df2))
  }

  # mode == "pid"
  if (!(pid_col %in% names(df))) stop("For mode='pid', pid_col must exist in df.")
  pid <- as.character(df[[pid_col]])
  if (anyNA(pid)) stop("pid_col contains NA.")

  # NEW: require one row per patient
  if (anyDuplicated(pid) > 0L) {
    stop("mode='pid' requires one row per patient (unique pid). ",
         "If you have multiple rows per patient, aggregate first or use mode='band'.")
  }

  # optional sorting covariate
  if (!is.null(x_sort_col)) {
    if (!(x_sort_col %in% names(df))) stop("x_sort_col not found in df.")
    x_sort <- suppressWarnings(as.numeric(as.character(df[[x_sort_col]])))
  } else {
    x_sort <- rep(NA_real_, nrow(df))
  }

  # directly use df (since pid is unique)
  dfp <- df
  pidp <- pid
  lop <- as.numeric(dfp[[lo_col]])
  hip <- as.numeric(dfp[[hi_col]])
  tp  <- y_true
  xp  <- x_sort

  # sort patients
  ord <- if (!is.null(x_sort_col)) order(xp, pidp) else order(pidp)
  dfp <- dfp[ord, , drop = FALSE]
  pidp <- pidp[ord]; lop <- lop[ord]; hip <- hip[ord]; tp <- tp[ord]; xp <- xp[ord]

  # cap
  if (!is.null(max_patients)) {
    max_patients <- as.integer(max_patients)
    if (max_patients > 0L && length(pidp) > max_patients) {
      keep <- seq_len(max_patients)
      dfp <- dfp[keep, , drop = FALSE]
      pidp <- pidp[keep]; lop <- lop[keep]; hip <- hip[keep]; tp <- tp[keep]; xp <- xp[keep]
    }
  }

  nP <- length(pidp)
  x <- seq_len(nP)

  ok_band <- is.finite(lop) & is.finite(hip)
  y_all <- c(lop[is.finite(lop)], hip[is.finite(hip)])
  if (isTRUE(show_true) && has_true) y_all <- c(y_all, tp[is.finite(tp)])
  if (length(y_all) == 0) stop("No finite values to plot.")
  yl <- range(y_all)
  pad <- 0.10 * diff(yl); if (!is.finite(pad) || pad <= 0) pad <- 1
  ylim_auto <- c(yl[1] - pad, yl[2] + pad)

  if (!("xlab" %in% names(dots))) {
    dots$xlab <- if (!is.null(x_sort_col)) paste0("Patients (sorted by ", x_sort_col, ")") else "Patients"
  }
  if (!("ylab" %in% names(dots))) dots$ylab <- "Y"

  plot_args <- c(
    list(x = x, y = hip, type = "n"),
    if (!has_ylim) list(ylim = ylim_auto) else list(),
    dots
  )
  do.call(graphics::plot, plot_args)

  graphics::abline(h = pretty(graphics::par("usr")[3:4]), col = grid_col, lty = 3)

  # intervals
  if (any(ok_band)) {
    for (i in which(ok_band)) {
      graphics::rect(x[i] - 0.35, lop[i], x[i] + 0.35, hip[i], col = band_col, border = NA)
    }
    graphics::segments(x[ok_band], lop[ok_band], x[ok_band], hip[ok_band], col = bound_col, lwd = 2)
  }

  if (isTRUE(show_center) && any(ok_band)) {
    center <- (lop + hip) / 2
    graphics::points(x[ok_band], center[ok_band], pch = 16, cex = 0.65, col = center_col)
  }

  if (isTRUE(show_true) && has_true) {
    ok_t <- is.finite(tp)
    if (any(ok_t)) graphics::points(x[ok_t], tp[ok_t], pch = 16, cex = 0.75, col = truth_col)
  }

  # ---- legend (pid mode): interval line + center + truth (no fill) ----
  leg <- c("Interval")
  col <- c(bound_col)
  lty <- c(1)
  lwd <- c(2)
  pch <- c(NA)
  pt_cex <- c(NA)

  if (isTRUE(show_center)) {
    leg <- c(leg, "Center")
    col <- c(col, center_col)
    lty <- c(lty, NA)
    lwd <- c(lwd, NA)
    pch <- c(pch, 16)
    pt_cex <- c(pt_cex, 0.7)
  }

  if (isTRUE(show_true) && has_true) {
    leg <- c(leg, "True y")
    col <- c(col, truth_col)
    lty <- c(lty, NA)
    lwd <- c(lwd, NA)
    pch <- c(pch, 16)
    pt_cex <- c(pt_cex, 0.8)
  }

  graphics::legend("topleft",
                   legend = leg,
                   bty = "n",
                   col = col,
                   lty = lty,
                   lwd = lwd,
                   pch = pch,
                   pt.cex = pt_cex)

  invisible(dfp)
}
