#' Plot HCP prediction intervals (band vs covariate or intervals by patient)
#'
#' @description
#' Unified plotting function for two common visualizations of HCP prediction intervals:
#' \itemize{
#'   \item \code{mode="band"}: plot an interval band (lo/hi) versus a 1D covariate (e.g., time X1).
#'     Optionally overlay a subset of training trajectories ("spaghetti") before drawing the band.
#'   \item \code{mode="pid"}: plot one interval per patient on the x-axis (patients optionally sorted by a covariate).
#' }
#'
#' In \code{mode="band"}, if \code{show_train=TRUE}, this function looks for a
#' data.frame named \code{dat_train} in the calling environment and overlays up to
#' \code{max_train_patients} training trajectories using columns \code{id}, \code{x_col},
#' and \code{Y} (preferred) or \code{Y_full}. The y-axis limits are determined from the
#' prediction band (and optional truth), so training trajectories may be clipped.
#'
#' @param df A data.frame containing prediction results. It must include the interval
#'   endpoints specified by \code{lo_col} and \code{hi_col}, and the covariate columns
#'   required by the chosen plotting mode.
#'
#' @param mode Plotting mode. Use \code{"band"} to visualize an interval band as a
#'   function of a continuous covariate, or \code{"pid"} to visualize one prediction
#'   interval per patient on the x-axis.
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
#' @param x_col (mode = "band") Name of the covariate column used as the x-axis in the
#'   interval band plot (e.g., time or another continuous predictor).
#'
#' @param show_train (mode = "band") Logical; if \code{TRUE}, overlay a subset of training
#'   trajectories ("spaghetti") from a data.frame named \code{dat_train} in the calling
#'   environment. Training trajectories are drawn before the prediction band and may be
#'   clipped by the plot y-limits.
#'
#' @param max_train_patients (mode = "band") Maximum number of training patients/trajectories
#'   to overlay when \code{show_train=TRUE}. Default is 30.
#'
#' @param pid_col (mode = "pid") Name of the column identifying patients (or clusters).
#'   Each patient must appear exactly once in \code{df}. Default is \code{"pid"}.
#'
#' @param x_sort_col (mode = "pid") Optional covariate column used to order patients along
#'   the x-axis (e.g., \code{"X1"}). If \code{NULL}, patients are ordered by their IDs.
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
#' ## CD4 example
#' ## ------------------------------------------------------------
#' ## Load MACS CD4 data (lqmix::cd4)
#' if (!requireNamespace("lqmix", quietly = TRUE)) {
#'   stop("Package 'lqmix' is required for this example.")
#' }
#' data(cd4, package = "lqmix")
#'
#' dat0 <- data.frame(
#'   id = cd4$sbj.id,
#'   X1 = cd4$time,
#'   X2 = cd4$age,
#'   X3 = cd4$packs,
#'   X4 = cd4$partners,
#'   X5 = cd4$drugs,
#'   X6 = cd4$cesd,
#'   Y_full = cd4$count
#' )
#' dat0 <- dat0[order(dat0$id, dat0$X1), ]
#' dat0$delta <- 1L; dat0$Y <- dat0$Y_full
#' x_cols <- sort(grep("^X\\d+$", names(dat0), value = TRUE))
#'
#' yr <- range(dat0$Y_full, finite = TRUE)
#' pad <- 0.10 * diff(yr); if (!is.finite(pad) || pad <= 0) pad <- 50
#' y_grid <- seq(max(0, yr[1] - pad), yr[2] + pad, length.out = 201)
#'
#' ## Quick check: training trajectories
#' tr0 <- setNames(dat0[, c("id","X1","Y_full")], c("pid","X1","Y"))
#' plot(range(tr0$X1), range(tr0$Y),
#'      type = "n", xlab = "Time (X1)", ylab = "CD4 count",
#'      main = sprintf("Training trajectories (%d subjects)", length(unique(tr0$pid))))
#' for (pp in unique(tr0$pid)) {
#'   ii <- which(tr0$pid == pp)
#'   if (length(ii) >= 2) lines(tr0$X1[ii], tr0$Y[ii], col = "grey70", lwd = 1)
#' }
#' legend("topright", bty = "n", legend = "Training trajectories", col = "grey70", lwd = 1)
#'
#' ## Case A: P=1, M>1 (band vs time for subjects)
#' for (pid in c(54, 92, 186)) {
#'   if (!any(dat0$id == pid)) next
#'   ## leave one out
#'   dat_test  <- dat0[dat0$id == pid, ]
#'   dat_train <- dat0[dat0$id != pid, ]
#'   x_test <- data.matrix(dat_test[, x_cols, drop = FALSE])
#'
#'   res <- hcp_conformal_region(
#'     dat = dat_train,
#'     id_col = "id",
#'     y_col = "Y",
#'     delta_col = "delta",
#'     x_cols = x_cols,
#'     x_test = x_test,
#'     y_grid = y_grid,
#'     alpha = 0.1, S = 1, B = 1,
#'     dens_method = "rq",
#'     dens_taus = seq(0.05, 0.95, by = 0.01),
#'     seed = 1
#'   )
#'   pred <- data.frame(
#'     X1 = dat_test$X1,
#'     lo = pmax(as.numeric(res$lo_hi[, 1]), 0),
#'     hi = pmax(as.numeric(res$lo_hi[, 2]), 0),
#'     y_true = pmax(as.numeric(dat_test$Y_full), 0)
#'   )
#'   plot_hcp_intervals(
#'     pred, mode = "band", x_col = "X1",
#'     y_true_col = "y_true", show_true = TRUE,
#' main = sprintf("Case A: prediction band vs time: pid=%s (M=%d)", pid, nrow(pred))
#'   )
#' }
#'
#' ## Case B: random subjects; one time point per subject
#' set.seed(2)
#' idsB <- sample(unique(dat0$id), 20)
#' datB <- dat0[dat0$id %in% idsB, ]
#' datB <- datB[!duplicated(datB$id), ]  # earliest (dat0 already ordered)
#' x_testB <- data.matrix(datB[, x_cols, drop = FALSE])
#'
#' resB <- hcp_conformal_region(
#'   dat = dat0[!(dat0$id %in% idsB), ],
#'   id_col = "id", y_col = "Y", delta_col = "delta",
#'   x_cols = x_cols, x_test = x_testB,
#'   y_grid = y_grid, alpha = 0.1, S = 1, B = 1,
#'   dens_method = "rq", dens_taus = seq(0.05, 0.95, by = 0.01), seed = 2
#' )
#' predB <- data.frame(pid = datB$id, X1 = datB$X1,
#'                     lo = pmax(as.numeric(resB$lo_hi[,1]), 0),
#'                     hi = pmax(as.numeric(resB$lo_hi[,2]), 0),
#'                     y_true = pmax(as.numeric(datB$Y_full), 0))
#' plot_hcp_intervals(predB, mode = "pid", pid_col = "pid", x_sort_col = "X1",
#'                    y_true_col = "y_true", show_true = TRUE,
#'                    main = "Case B: random subjects, one time point")
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
    # NEW (minimal): overlay training trajectories (spaghetti)
    show_train = FALSE,
    max_train_patients = 30,
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

  # fixed training style (no extra params)
  train_col <- "grey70"
  train_lwd <- 1
  train_pch <- 16
  train_cex <- 0.5

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

    ## ---- OPTIONAL: training spaghetti (from 'dat_train' in caller env) ----
    tr_sp <- NULL
    has_train <- isTRUE(show_train)
    if (has_train) {
      if (!exists("dat_train", envir = parent.frame(), inherits = TRUE)) {
        stop("show_train=TRUE requires an object named 'dat_train' in the calling environment.")
      }
      dat_train <- get("dat_train", envir = parent.frame(), inherits = TRUE)
      if (!is.data.frame(dat_train)) stop("'dat_train' must be a data.frame.")

      if (!("id" %in% names(dat_train))) stop("dat_train must contain column 'id' (patient id).")
      if (!(x_col %in% names(dat_train))) stop("dat_train must contain x_col column: ", x_col)

      ycol_tr <- if ("Y" %in% names(dat_train)) "Y" else if ("Y_full" %in% names(dat_train)) "Y_full" else NULL
      if (is.null(ycol_tr)) stop("dat_train must contain either 'Y' or 'Y_full' for training trajectories.")

      tr_sp <- dat_train[, c("id", x_col, ycol_tr), drop = FALSE]
      names(tr_sp) <- c("pid_tr", "x_tr", "y_tr")
      tr_sp$pid_tr <- as.character(tr_sp$pid_tr)
      tr_sp$x_tr <- as.numeric(tr_sp$x_tr)
      tr_sp$y_tr <- as.numeric(tr_sp$y_tr)
      tr_sp <- tr_sp[is.finite(tr_sp$x_tr) & is.finite(tr_sp$y_tr) & !is.na(tr_sp$pid_tr), , drop = FALSE]

      if (nrow(tr_sp) == 0) {
        tr_sp <- NULL
        has_train <- FALSE
      } else {
        max_train_patients <- as.integer(max_train_patients)
        if (!is.finite(max_train_patients) || max_train_patients <= 0L) {
          tr_sp <- NULL
          has_train <- FALSE
        } else {
          pids <- unique(tr_sp$pid_tr)
          pids <- pids[seq_len(min(max_train_patients, length(pids)))]
          tr_sp <- tr_sp[tr_sp$pid_tr %in% pids, , drop = FALSE]
          tr_sp <- tr_sp[order(tr_sp$pid_tr, tr_sp$x_tr), , drop = FALSE]
        }
      }
    }

    # y-limits: determined by band (+ optional truth) ONLY.
    # Training trajectories may be clipped (by design).
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

    graphics::abline(h = pretty(graphics::par("usr")[3:4]),
                     v = pretty(graphics::par("usr")[1:2]),
                     col = grid_col, lty = 3)

    # training trajectories FIRST (may be clipped)
    if (has_train && !is.null(tr_sp)) {
      for (pp in unique(tr_sp$pid_tr)) {
        ii <- which(tr_sp$pid_tr == pp)
        if (length(ii) >= 2) {
          graphics::lines(tr_sp$x_tr[ii], tr_sp$y_tr[ii], col = train_col, lwd = train_lwd)
        } else if (length(ii) == 1) {
          graphics::points(tr_sp$x_tr[ii], tr_sp$y_tr[ii], col = train_col,
                           pch = train_pch, cex = train_cex)
        }
      }
    }

    # band fill ON TOP of spaghetti
    if (any(ok_band)) {
      xb <- x[ok_band]; lob <- lo2[ok_band]; hib <- hi2[ok_band]
      graphics::polygon(c(xb, rev(xb)), c(lob, rev(hib)), col = band_col, border = NA)
    }

    # bounds ON TOP
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

    # legend: lines/points only (no fill)
    leg2 <- c("Bounds")
    col2 <- c(bound_col)
    lty2 <- c(1)
    lwd2 <- c(2)
    pch2 <- c(NA)
    pt_cex2 <- c(NA)

    if (has_train && !is.null(tr_sp)) {
      leg2 <- c(leg2, "Training trajectories")
      col2 <- c(col2, train_col)
      lty2 <- c(lty2, 1)
      lwd2 <- c(lwd2, train_lwd)
      pch2 <- c(pch2, NA)
      pt_cex2 <- c(pt_cex2, NA)
    }

    if (isTRUE(show_center)) {
      leg2 <- c(leg2, "Center")
      col2 <- c(col2, center_col)
      lty2 <- c(lty2, 2)
      lwd2 <- c(lwd2, 1.5)
      pch2 <- c(pch2, NA)
      pt_cex2 <- c(pt_cex2, NA)
    }

    if (isTRUE(show_true) && has_true) {
      leg2 <- c(leg2, "True y")
      col2 <- c(col2, truth_col)
      lty2 <- c(lty2, NA)
      lwd2 <- c(lwd2, NA)
      pch2 <- c(pch2, 16)
      pt_cex2 <- c(pt_cex2, 0.8)
    }

    graphics::legend("topright",
                     inset = c(0, 0.12),
                     legend = leg2,
                     bty = "n",
                     col = col2,
                     lty = lty2,
                     lwd = lwd2,
                     pch = pch2,
                     pt.cex = pt_cex2)

    return(invisible(df2))
  }

  # -------------------- mode == "pid" (unchanged) --------------------
  if (!(pid_col %in% names(df))) stop("For mode='pid', pid_col must exist in df.")
  pid <- as.character(df[[pid_col]])
  if (anyNA(pid)) stop("pid_col contains NA.")

  if (anyDuplicated(pid) > 0L) {
    stop("mode='pid' requires one row per patient (unique pid). ",
         "If you have multiple rows per patient, aggregate first or use mode='band'.")
  }

  if (!is.null(x_sort_col)) {
    if (!(x_sort_col %in% names(df))) stop("x_sort_col not found in df.")
    x_sort <- suppressWarnings(as.numeric(as.character(df[[x_sort_col]])))
  } else {
    x_sort <- rep(NA_real_, nrow(df))
  }

  dfp <- df
  pidp <- pid
  lop <- as.numeric(dfp[[lo_col]])
  hip <- as.numeric(dfp[[hi_col]])
  tp  <- y_true
  xp  <- x_sort

  ord <- if (!is.null(x_sort_col)) order(xp, pidp) else order(pidp)
  dfp <- dfp[ord, , drop = FALSE]
  pidp <- pidp[ord]; lop <- lop[ord]; hip <- hip[ord]; tp <- tp[ord]; xp <- xp[ord]

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

  graphics::legend("topright",
                   legend = leg,
                   bty = "n",
                   col = col,
                   lty = lty,
                   lwd = lwd,
                   pch = pch,
                   pt.cex = pt_cex)

  invisible(dfp)
}

