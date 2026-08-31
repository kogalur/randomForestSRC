## mortality-stratified Brier curves used by plot.survival
.surv.brier.strata <- function(mort, brier.matx, brier.all) {
  ntime <- length(brier.all)
  brier.strata <- matrix(NA_real_, nrow = ntime, ncol = 4L)
  mort.ok <- is.finite(mort)
  if (any(mort.ok)) {
    mort.perc <- c(min(mort[mort.ok]) - 1e-5,
                   quantile(mort[mort.ok], (1:4) / 4, na.rm = TRUE))
    for (k in 1:4) {
      mort.pt <- mort.ok & mort > mort.perc[k] &
        mort <= mort.perc[k + 1L]
      if (any(mort.pt)) {
        brier.strata[, k] <- colMeans(
          brier.matx[mort.pt, , drop = FALSE], na.rm = TRUE)
      }
    }
  }
  brier.strata[!is.finite(brier.strata)] <- NA_real_
  out <- data.frame(brier.strata, bs.all = brier.all)
  colnames(out) <- c("bs.q25", "bs.q50", "bs.q75", "bs.q100", "bs.all")
  out
}
## cumulative integrated Brier curves used by plot.survival
.surv.crps.curves <- function(time, brier.score) {
  crps <- matrix(NA_real_, nrow = length(time), ncol = ncol(brier.score))
  for (k in seq_len(ncol(brier.score))) {
    for (j in seq_along(time)) {
      if (j > 1L) {
        denom <- diff(range(time[seq_len(j)]))
        if (is.finite(denom) && denom > 0) {
          crps[j, k] <- trapz(time[seq_len(j)],
                              brier.score[seq_len(j), k]) / denom
        }
      }
    }
  }
  crps <- data.frame(crps)
  colnames(crps) <- c("crps.q25", "crps.q50", "crps.q75",
                      "crps.q100", "crps.all")
  crps
}
## draw individual and mean survival curves
.plot.survival.curves <- function(time, surv.ensb, surv.mean.ensb,
                                  surv.aalen, subset, collapse, ylab, ...) {
  curve <- surv.ensb
  if (!collapse && length(subset) > 500L) {
    r.pt <- sample(seq_along(subset), 500L, replace = FALSE)
    curve <- surv.ensb[, r.pt, drop = FALSE]
  }
  if (!any(is.finite(curve))) {
    plot.new()
    title(xlab = "Time", ylab = ylab)
    return(invisible(NULL))
  }
  matplot(time,
          curve,
          xlab = "Time",
          ylab = ylab,
          type = "l",
          col = 1,
          lty = 3, ...)
  if (!is.null(surv.aalen)) {
    lines(time, surv.aalen, lty = 1, col = 3, lwd = 3)
  }
  lines(time, surv.mean.ensb, lty = 1, col = 2, lwd = 3)
  rug(time, ticksize = -0.03)
  invisible(NULL)
}
plot.survival.rfsrc <- function(x,
                                 show.plots = TRUE,
                                 subset, collapse = FALSE,
                                 cens.model = c("km", "rfsrc"),
                                 ...)
{
  ## incoming parameter checks
  if (is.null(x)) {
    stop("object x is empty!")
  }
  if ((!inherits(x, "rfsrc") || !inherits(x, "grow")) &&
      (!inherits(x, "rfsrc") || !inherits(x, "predict"))) {
    stop("This function only works for objects of class `(rfsrc, grow)' or '(rfsrc, predict)'")
  }
  if (x$family != "surv") {
    stop("this function only supports right-censored survival settings")
  }
  cens.model <- match.arg(cens.model)
  subset.arg <- if (missing(subset)) NULL else subset
  ## Extract predictions before requesting a performance score. This permits
  ## prediction objects without outcomes to display survival curves.
  plot.data <- .get.survival.performance.data(x,
                                               subset = subset.arg,
                                               times = NULL)
  has.outcomes <- !is.null(plot.data$eval.yvar) &&
    any(is.finite(plot.data$eval.event.info$time) &
        is.finite(plot.data$eval.event.info$cens))
  ## A prediction object without outcomes has no observable target for Brier
  ## score, CRPS, or mortality-versus-time diagnostics.
  if (!has.outcomes) {
    surv.ensb <- t(plot.data$survival)
    surv.mean.ensb <- rowMeans(surv.ensb, na.rm = TRUE)
    surv.plot <- if (collapse) surv.mean.ensb else surv.ensb
    if (show.plots) {
      old.par <- par(no.readonly = TRUE)
      on.exit(par(old.par), add = TRUE)
      par(mfrow = c(1, 1), cex = 1.0)
      .plot.survival.curves(plot.data$time,
                            surv.plot,
                            surv.mean.ensb,
                            surv.aalen = NULL,
                            subset = plot.data$subset,
                            collapse = collapse,
                            ylab = "Survival", ...)
    }
    return(invisible(list(time = plot.data$time,
                          survival = surv.ensb,
                          survival.mean = surv.mean.ensb,
                          subset = plot.data$subset)))
  }
  ## Evaluation outcomes are aligned with the predictions, while the
  ## censoring distribution is always estimated from the full grow data.
  brier.obj <- get.brier.survival(x,
                                   subset = subset.arg,
                                   cens.model = cens.model,
                                   conf.int = FALSE,
                                   keep.matrix = TRUE)
  brier.matx <- brier.obj$brier.matx
  brier.all <- brier.obj$brier.score$brier.score
  mort <- brier.obj$mort
  surv.ensb <- brier.obj$surv.ensb
  surv.aalen <- brier.obj$surv.aalen
  event.info <- brier.obj$eval.event.info
  subset <- brier.obj$subset
  time <- brier.obj$time
  ## Brier and cumulative integrated Brier processing
  brier.score <- .surv.brier.strata(mort, brier.matx, brier.all)
  crps <- .surv.crps.curves(time, brier.score)
  ## labels and titles
  if (identical(brier.obj$prediction, "oob")) {
    ylab.1 <- "OOB Survival"
    ylab.2 <- "OOB Brier"
    ylab.3 <- "OOB CRPS"
    ylab.4 <- "OOB Mortality vs Time"
  }
  else {
    ylab.1 <- "Survival"
    ylab.2 <- "Brier Score"
    ylab.3 <- "CRPS"
    ylab.4 <- "Mortality vs Time"
  }
  ## mean ensemble survival
  surv.mean.ensb <- rowMeans(surv.ensb, na.rm = TRUE)
  surv.plot <- if (collapse) surv.mean.ensb else surv.ensb
  ## ------------------------------------------------------------
  ##
  ## plots
  ##
  ## ------------------------------------------------------------
  if (show.plots) {
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par), add = TRUE)
    par(mfrow = c(2, 2), cex = 1.0)
    ## ----survival plot----
    .plot.survival.curves(time,
                          surv.plot,
                          surv.mean.ensb,
                          surv.aalen,
                          subset,
                          collapse,
                          ylab.1, ...)
    ## ----brier plot----
    brier.values <- unlist(brier.score)
    point.x <- max(1L, round(3 * length(time) / 4))
    if (any(is.finite(brier.values))) {
      matplot(time, brier.score,
              xlab = "Time",
              ylab = ylab.2,
              type = "l",
              lwd = c(rep(1, 4), 2),
              col = c(rep(1, 4), 2),
              lty = c(1:4, 1), ...)
      rng <- range(brier.values[is.finite(brier.values)])
      abline(h = seq(rng[1], rng[2], length = 20),
             col = gray(.6), lty = 3, lwd = .85)
      text(time[point.x], brier.score[point.x, 1], "0-25", col = 4)
      text(time[point.x], brier.score[point.x, 2], "25-50", col = 4)
      text(time[point.x], brier.score[point.x, 3], "50-75", col = 4)
      text(time[point.x], brier.score[point.x, 4], "75-100", col = 4)
      rug(time, ticksize = 0.03)
    }
    else {
      plot.new()
      title(xlab = "Time", ylab = ylab.2)
    }
    ## ----CRPS plot----
    crps.values <- unlist(crps)
    if (any(is.finite(crps.values))) {
      matplot(time, crps,
              xlab = "Time",
              ylab = ylab.3,
              type = "l",
              lwd = c(rep(1, 4), 2),
              col = c(rep(1, 4), 2),
              lty = c(1:4, 1), ...)
      rng <- range(crps.values[is.finite(crps.values)])
      abline(h = seq(rng[1], rng[2], length = 20),
             col = gray(.6), lty = 3, lwd = .85)
      text(time[point.x], crps[point.x, 1], "0-25", col = 4)
      text(time[point.x], crps[point.x, 2], "25-50", col = 4)
      text(time[point.x], crps[point.x, 3], "50-75", col = 4)
      text(time[point.x], crps[point.x, 4], "75-100", col = 4)
      rug(time, ticksize = 0.03)
    }
    else {
      plot.new()
      title(xlab = "Time", ylab = ylab.3)
    }
    ## ----mortality plot----
    ytime <- event.info$time
    status <- event.info$cens
    plot.pt <- is.finite(ytime) & is.finite(mort)
    if (any(plot.pt)) {
      plot(ytime[plot.pt], mort[plot.pt],
           xlab = "Time", ylab = ylab.4, type = "n", ...)
      cex <- if (length(mort) > 500L) 0.5 else 0.75
      event.pt <- plot.pt & is.finite(status) & status != 0
      cens.pt <- plot.pt & is.finite(status) & status == 0
      points(ytime[event.pt], mort[event.pt], pch = 16,
             col = 4, cex = cex)
      points(ytime[cens.pt], mort[cens.pt], pch = 16, cex = cex)
      if (sum(event.pt) > 1L) {
        ord <- order(ytime[event.pt])
        lines(supsmu(ytime[event.pt][ord], mort[event.pt][ord]),
              lty = 3, col = 4)
      }
      if (sum(cens.pt) > 1L) {
        ord <- order(ytime[cens.pt])
        lines(supsmu(ytime[cens.pt][ord], mort[cens.pt][ord]),
              lty = 3)
      }
      if (any(event.pt)) {
        rug(ytime[event.pt], ticksize = -0.03)
      }
    }
    else {
      plot.new()
      title(xlab = "Time", ylab = ylab.4)
    }
  }
  ## Invisibly return the Brier and cumulative integrated Brier curves.
  invisible(data.frame(time = time, brier.score, crps))
}
plot.survival <- plot.survival.rfsrc
## draw a confidence polygon over contiguous finite portions of a curve
.surv.plot.band <- function(x, lower, upper, col, border = NA) {
  ok <- is.finite(x) & is.finite(lower) & is.finite(upper)
  if (!any(ok)) {
    return(invisible(NULL))
  }
  run <- rle(ok)
  run.end <- cumsum(run$lengths)
  run.start <- run.end - run$lengths + 1L
  for (k in which(run$values)) {
    idx <- run.start[k]:run.end[k]
    if (length(idx) >= 2L) {
      graphics::polygon(c(x[idx], rev(x[idx])),
                        c(lower[idx], rev(upper[idx])),
                        col = col, border = border)
    }
  }
  invisible(NULL)
}
## draw one Brier or AUCT performance curve using base R graphics
.surv.plot.performance.curve <- function(score, estimate,
                                         ylab, main, reference = NULL,
                                         dots = list()) {
  time <- score$time
  value <- score[[estimate]]
  has.band <- all(c("lower", "upper") %in% names(score))
  yy <- value
  if (has.band) {
    yy <- c(yy, score$lower, score$upper)
  }
  if (!is.null(reference)) {
    yy <- c(yy, reference)
  }
  yy <- yy[is.finite(yy)]
  ylim <- if (length(yy) > 0L) range(yy) else c(0, 1)
  if (diff(ylim) == 0) {
    ylim <- ylim + c(-0.04, 0.04) * max(1, abs(ylim[1L]))
  }
  line.args <- list(col = 2, lwd = 2, lty = 1, type = "l")
  for (nm in intersect(c("col", "lwd", "lty", "type"), names(dots))) {
    line.args[[nm]] <- dots[[nm]]
    dots[[nm]] <- NULL
  }
  band.col <- if (!is.null(dots$band.col)) {
    dots$band.col
  }
  else {
    grDevices::adjustcolor(gray(.5), alpha.f = .25)
  }
  band.border <- if (!is.null(dots$band.border)) dots$band.border else NA
  dots$band.col <- NULL
  dots$band.border <- NULL
  plot.args <- list(x = time, y = value, type = "n",
                    xlab = "Time", ylab = ylab,
                    main = main, ylim = ylim)
  dot.names <- names(dots)
  if (!is.null(dot.names)) {
    named <- nzchar(dot.names)
    for (nm in dot.names[named]) {
      plot.args[[nm]] <- dots[[nm]]
    }
    if (any(!named)) {
      plot.args <- c(plot.args, dots[!named])
    }
  }
  do.call(plot, plot.args)
  if (has.band) {
    .surv.plot.band(time, score$lower, score$upper,
                    col = band.col, border = band.border)
  }
  if (!is.null(reference)) {
    abline(h = reference, lty = 3, col = gray(.5))
  }
  do.call(lines, c(list(x = time, y = value), line.args))
  rug(time, ticksize = -0.025)
}
## base-R plot helper for Brier score and time-dependent AUC
plotBrierAUC <- function(x,
                           subset = NULL,
                           cens.model = c("km", "rfsrc"),
                           papply = lapply,
                           times = NULL,
                           conf.int = TRUE,
                           plots = c("brier", "auct"),
                           show.plots = TRUE,
                           ...) {
  plots <- as.character(plots)
  if (length(plots) == 0L) {
    stop("'plots' must include 'brier', 'auct', or both.")
  }
  plots[plots == "auc"] <- "auct"
  plots <- unique(match.arg(plots, c("brier", "auct"), several.ok = TRUE))
  ans <- .get.survival.performance(x,
                                   subset = subset,
                                   cens.model = cens.model,
                                   papply = papply,
                                   times = times,
                                   conf.int = conf.int,
                                   keep.matrix = FALSE,
                                   metrics = plots)
  if (is.null(ans)) {
    return(invisible(NULL))
  }
  if (show.plots) {
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par), add = TRUE)
    par(mfrow = c(length(plots), 1L), cex = 1.0)
    dots <- list(...)
    for (what in plots) {
      if (what == "brier") {
        .surv.plot.performance.curve(
          ans$brier$brier.score,
          estimate = "brier.score",
          ylab = "Brier Score",
          main = "Brier Score Over Time",
          dots = dots
        )
      }
      else {
        .surv.plot.performance.curve(
          ans$auct$auct.score,
          estimate = "auct",
          ylab = "AUC",
          main = "Time-Dependent AUC",
          reference = 0.5,
          dots = dots
        )
      }
    }
  }
  invisible(ans)
}
