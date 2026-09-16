## ------------------------------------------------------------------
## Plotting helpers. Argument overrides are merged once by name; the
## fitted values and native prediction settings are never replaced by ... .
## ------------------------------------------------------------------
.plot.variable.merge <- function(defaults, overrides) {
  if (length(overrides)) {
    nm <- names(overrides)
    if (is.null(nm) || anyNA(nm) || any(!nzchar(nm))) {
      stop("graphical arguments in '...' must be named")
    }
    for (i in seq_along(overrides)) defaults[nm[i]] <- overrides[i]
  }
  defaults
}
.plot.variable.flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("'", name, "' must be TRUE or FALSE")
  }
  value
}
.plot.variable.count <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop("'", name, "' must be a finite number")
  }
  max(1, round(value))
}
## One summary definition for numeric partial predictions in every family.
## The sample variance and its divisor use exactly the same finite rows.
.plot.variable.summary <- function(pred) {
  if (is.null(dim(pred))) stop("partial predictions must be an observation-by-grid matrix")
  pred[!is.finite(pred)] <- NA_real_
  n.used <- colSums(!is.na(pred))
  center <- colMeans(pred, na.rm = TRUE)
  center[n.used == 0L] <- NA_real_
  spread <- vapply(seq_len(ncol(pred)), function(j) {
    if (n.used[j] > 1L) sd(pred[, j], na.rm = TRUE) else NA_real_
  }, numeric(1L))
  se <- spread / sqrt(n.used)
  se[n.used < 2L] <- NA_real_
  list(mean = center, sd = spread, se = se, n.used = n.used, pred = pred)
}
.plot.variable.range <- function(x, fraction = .005) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(0, 1))
  z <- range(x)
  pad <- fraction * diff(z)
  if (!is.finite(pad) || pad == 0) pad <- fraction * max(1, abs(z))
  z + c(-pad, pad)
}
## Keep box-specific settings off ordinary plot()/points()/lines() calls.
.plot.variable.box.names <- function() {
  c("notch", "outline", "boxfill", "boxcol", "boxlty", "boxlwd",
    "whiskcol", "whisklty", "whisklwd", "staplecol", "staplelty", "staplelwd",
    "medcol", "medlty", "medlwd", "medpch", "medcex", "medbg", "outcol", "outpch",
    "outcex", "outbg", "outlty", "outlwd", "boxwex", "staplewex", "outwex",
    "notch.frac", "varwidth", "width", "border", "horizontal", "at", "show.names", "names", "range")
}
.plot.variable.panel <- function(x, y, dots, xlab, ylab) {
  args <- .plot.variable.merge(list(x = x, y = y, type = "n", xlab = xlab, ylab = ylab),
                               dots[setdiff(names(dots), .plot.variable.box.names())])
  ## Scaffolding is separate from the styled fitted curve and observations.
  args$x <- x
  args$y <- y
  args$type <- "n"
  if (!any(is.finite(x)) && is.null(args$xlim)) args$xlim <- c(0, 1)
  if (!any(is.finite(y)) && is.null(args$ylim)) args$ylim <- c(0, 1)
  do.call("plot", args)
}
## Smooth only finite contiguous pieces, keeping unavailable grid points as gaps.
.plot.variable.line <- function(x, y, smooth, args = list()) {
  good <- is.finite(x) & is.finite(y)
  runs <- rle(good)
  finish <- cumsum(runs$lengths)
  start <- finish - runs$lengths + 1L
  for (j in which(runs$values)) {
    idx <- seq.int(start[j], finish[j])
    if (length(idx) < 2L) next
    xy <- if (smooth && length(unique(x[idx])) > 1L) lowess(x[idx], y[idx]) else
      list(x = x[idx], y = y[idx])
    do.call("lines", .plot.variable.merge(xy, args))
  }
  invisible(NULL)
}
.plot.variable.boxes <- function(groups, labels, dots, xlab, ylab, partial = FALSE) {
  groups <- lapply(groups, function(z) z[is.finite(z)])
  names(groups) <- labels
  if (!any(lengths(groups))) {
    .plot.variable.panel(seq_along(groups), rep(NA_real_, length(groups)), dots, xlab, ylab)
    text(mean(par("usr")[1:2]), mean(par("usr")[3:4]), "No finite predictions")
    return(invisible(NULL))
  }
  box.range <- if (is.null(dots$range)) 1.5 else dots$range
  bp <- boxplot(groups, plot = FALSE, range = box.range)
  labels <- if (is.null(dots$names)) labels else dots$names
  if (length(labels) != length(groups)) stop("'names' must have one label per displayed group")
  bp$names <- as.character(labels)
  horizontal <- if (is.null(dots$horizontal)) FALSE else .plot.variable.flag(dots$horizontal, "horizontal")
  axes <- if (is.null(dots$axes)) TRUE else .plot.variable.flag(dots$axes, "axes")
  show.names <- if (is.null(dots$show.names)) TRUE else .plot.variable.flag(dots$show.names, "show.names")
  xaxt <- if (is.null(dots$xaxt)) par("xaxt") else dots$xaxt
  yaxt <- if (is.null(dots$yaxt)) par("yaxt") else dots$yaxt
  defaults <- list(z = bp, outline = FALSE, notch = TRUE, boxfill = "bisque",
                   horizontal = horizontal,
                   xlab = if (horizontal) ylab else xlab,
                   ylab = if (horizontal) xlab else ylab)
  if (partial) {
    limits <- bp$stats
    if (isTRUE(dots$outline)) limits <- c(limits, bp$out)
    if (is.null(dots$notch) || isTRUE(dots$notch)) limits <- c(limits, bp$conf)
    defaults[[if (horizontal) "xlim" else "ylim"]] <- .plot.variable.range(limits)
  }
  args <- .plot.variable.merge(defaults, dots[setdiff(names(dots), c("type", "range", "names"))])
  args$z <- bp
  args$axes <- FALSE
  args$show.names <- FALSE
  if (is.null(args$frame.plot)) args$frame.plot <- axes
  if (is.null(args$at)) args$at <- seq_along(labels)
  if (!is.numeric(args$at) || length(args$at) != length(labels) || any(!is.finite(args$at))) {
    stop("'at' must give one finite position per displayed group")
  }
  ## bxp reverses its limit arguments for horizontal boxes. Keep public
  ## xlim and ylim attached to the displayed horizontal and vertical axes.
  if (horizontal) {
    for (pair in list(c("xlim", "ylim"), c("xaxs", "yaxs"))) {
      first <- args[[pair[1L]]]
      second <- args[[pair[2L]]]
      args[[pair[1L]]] <- second
      args[[pair[2L]]] <- first
    }
  }
  do.call("bxp", args)
  axis.names <- c("las", "cex.axis", "col.axis", "font.axis", "line", "lwd",
                  "lwd.ticks", "col.ticks", "hadj", "padj", "gap.axis", "mgp",
                  "tck", "tcl", "xpd", "fg")
  axis.args <- dots[intersect(names(dots), axis.names)]
  value.side <- if (horizontal) 1L else 2L
  label.side <- if (horizontal) 2L else 1L
  if (axes && !identical(if (horizontal) xaxt else yaxt, "n")) {
    do.call("axis", .plot.variable.merge(list(side = value.side), axis.args))
  }
  if (axes && show.names && !identical(if (horizontal) yaxt else xaxt, "n")) {
    idx <- unique(round(pretty(seq_along(labels), min(30L, length(labels)))))
    idx <- idx[idx >= 1L & idx <= length(labels)]
    do.call("axis", .plot.variable.merge(list(side = label.side, at = args$at[idx],
                        labels = labels[idx], tick = TRUE), axis.args))
  }
  invisible(bp)
}
plot.variable.rfsrc <- function(
  x,
  xvar.names,
  target,
  m.target = NULL,
  time,
  surv.type = c("mort", "rel.freq", "surv", "years.lost", "cif", "chf"),
  class.type = c("prob", "bayes"),
  partial = FALSE,
  oob = TRUE,
  show.plots = TRUE,
  plots.per.page = 4,
  granule = 5,
  sorted = TRUE,
  nvar,
  npts = 25,
  smooth.lines = FALSE,
  subset,
  ...)
{
  ## --------------------------------------------------------------
  ## Object and display controls. Cached plots need no saved forest.
  ## --------------------------------------------------------------
  if (inherits(x, "rfsrc") && inherits(x, "synthetic")) x <- x$rfSyn
  object <- x
  if (!inherits(object, "rfsrc") ||
      !(inherits(object, "grow") || inherits(object, "plot.variable"))) {
    stop("this function requires an rfsrc grow, synthetic, or plot.variable object")
  }
  cached <- inherits(object, "plot.variable")
  if (identical(object$family, "unsupv")) stop("this function does not apply to unsupervised forests")
  if (inherits(object, "anonymous")) stop("this function does not work with anonymous forests")
  partial <- .plot.variable.flag(partial, "partial")
  oob <- .plot.variable.flag(oob, "oob")
  show.plots <- .plot.variable.flag(show.plots, "show.plots")
  sorted <- .plot.variable.flag(sorted, "sorted")
  smooth.lines <- .plot.variable.flag(smooth.lines, "smooth.lines")
  plots.per.page <- .plot.variable.count(plots.per.page, "plots.per.page")
  granule <- .plot.variable.count(granule, "granule")
  npts <- .plot.variable.count(npts, "npts")
  if (!missing(nvar)) nvar <- .plot.variable.count(nvar, "nvar")
  if (partial && !cached && (is.null(object$forest) || !isTRUE(object$forest$forest))) {
    stop("forest is empty: re-run rfsrc (grow call) with forest=TRUE")
  }
  dots <- .plot.variable.merge(list(), list(...))
  if (!is.null(dots$pars)) {
    if (!is.list(dots$pars)) stop("'pars' must be a named list of graphical arguments")
    pars <- dots$pars
    dots$pars <- NULL
    dots <- .plot.variable.merge(pars, dots)
  }
  if (any(c("x", "y", "z") %in% names(dots))) stop("plotting arguments must not replace stored plot data")
  ## --------------------------------------------------------------
  ## New calculation: retain the training row, factor, and event maps.
  ## --------------------------------------------------------------
  if (!cached) {
    xvar <- object$xvar
    if (is.null(xvar) || is.null(dim(xvar))) stop("stored training predictors are unavailable")
    if (!identical(class(xvar), "data.frame")) xvar <- as.data.frame(xvar, optional = TRUE)
    if (!is.null(object$imputed.indv)) {
      xvar[object$imputed.indv, ] <- object$imputed.data[, object$xvar.names, drop = FALSE]
    }
    n.full <- nrow(xvar)
    subset <- .plot.variable.subset(if (missing(subset)) NULL else subset, n.full)
    xvar <- xvar[subset, , drop = FALSE]
    n <- nrow(xvar)
    family <- object$family
    m.target <- get.univariate.target(object, m.target)
    target.arg <- if (missing(target)) NULL else target
    event.info <- NULL
    VIMP <- NULL
    if (family %in% c("surv", "surv-CR")) {
      event.info <- get.event.info(object, subset)
      ## The full forest map controls positions, even if a subset omits an event.
      full.events <- object$forest$event.info$event.type
      if (is.null(full.events)) full.events <- get.event.info(object)$event.type
      event.info$event.type <- full.events
      allowed <- if (family == "surv") c("mort", "rel.freq", "surv", "chf") else
        c("years.lost", "cif", "chf")
      if (missing(surv.type)) surv.type <- allowed[1L]
      if (length(surv.type) != 1L || is.na(surv.type)) stop("'surv.type' must select one survival summary")
      pred.type <- match.arg(surv.type, allowed)
      if (pred.type %in% c("surv", "cif", "chf")) {
        if (missing(time)) time <- median(event.info$time.interest, na.rm = TRUE)
        .partial.time.index(time, event.info$time.interest, "time", single = TRUE)
      }
      else time <- NULL
      if (family == "surv-CR") {
        if (!is.null(target.arg) && !is.numeric(target.arg)) stop("competing-risk 'target' must be an event index")
        target <- .plot.variable.target(target.arg, seq_along(full.events))
        event.label <- full.events[target]
        if (!is.null(object$importance)) VIMP <- object$importance[, target]
        ylabel <- switch(pred.type,
          "years.lost" = paste("Years lost for event", event.label),
          "cif" = paste0("CIF for event ", event.label, " (time=", time, ")"),
          "chf" = paste0("CHF for event ", event.label, " (time=", time, ")"))
      }
      else {
        target <- 1L
        VIMP <- object$importance
        ylabel <- switch(pred.type, "rel.freq" = "standardized mortality",
          "mort" = "mortality", "chf" = paste0("CHF (time=", time, ")"),
          "surv" = paste0("predicted survival (time=", time, ")"))
      }
    }
    else {
      if (grepl("surv", family)) stop("this survival family is not supported by variable plotting")
      time <- NULL
      univariate <- coerce.multivariate(object, m.target)
      if (is.factor(univariate$yvar)) {
        levels.y <- levels(univariate$yvar)
        pred.type <- match.arg(class.type, c("prob", "bayes"))
        target <- .plot.variable.target(target.arg, levels.y)
        if (pred.type == "prob") {
          if (!is.null(univariate$importance)) VIMP <- univariate$importance[, 1L + target]
          ylabel <- paste("probability", levels.y[target])
        }
        else {
          if (!is.null(univariate$importance)) VIMP <- univariate$importance[, 1L]
          ylabel <- "bayes class code"
        }
      }
      else {
        pred.type <- "y"
        target <- NULL
        VIMP <- univariate$importance
        ylabel <- expression(hat(y))
      }
    }
    ## Keep finite importance values ahead of unavailable ones.
    if (missing(xvar.names)) xvar.names <- object$xvar.names else
      xvar.names <- intersect(xvar.names, object$xvar.names)
    if (!length(xvar.names)) stop("none of the supplied predictors match available variables")
    if (sorted && !is.null(VIMP)) {
      score <- VIMP[xvar.names]
      score[!is.finite(score)] <- NA_real_
      xvar.names <- xvar.names[order(score, decreasing = TRUE, na.last = TRUE)]
    }
    if (!missing(nvar)) xvar.names <- head(xvar.names, nvar)
    nvar <- length(xvar.names)
    plots.per.page <- min(plots.per.page, nvar)
    ## Ordinary prediction dispatch is retained; extraction selects stored data.
    if (!partial) {
      yhat <- extract.pred(predict(object, m.target = m.target), pred.type,
                            subset, time, m.target, target, oob = oob)
    }
    else {
      prtl <- lapply(xvar.names, function(name) {
        xx <- xvar[[name]]
        factor.x <- is.factor(xx)
        if (factor.x) {
          labels <- .partial.factor.levels(object$forest, name)
          if (!length(labels)) labels <- levels(xx)
          codes <- which(labels %in% as.character(xx[!is.na(xx)]))
          grid <- as.numeric(codes)
          display <- factor(labels[codes], levels = labels[codes], ordered = is.ordered(xx))
          xx <- xx[!is.na(xx)]
        }
        else {
          xx <- xx[is.finite(xx)]
          grid <- sort(unique(xx))
          if (length(grid) > npts) {
            grid <- grid[unique(as.integer(seq(1, length(grid), length.out = npts)))]
          }
          display <- grid
        }
        if (!length(grid)) stop("predictor '", name, "' has no available values in the plotting subset")
        n.x <- length(grid)
        discrete <- factor.x || n.x <= granule
        raw <- extract.partial.pred(partial.rfsrc(object$forest,
                     m.target = m.target, partial.type = pred.type,
                     partial.xvar = name, partial.values = grid,
                     partial.time = time, oob = oob),
                   pred.type, subset, m.target, target)
        summary <- .plot.variable.summary(raw)
        if (discrete) {
          ## Retain the boxplot convention on the same SE scale as curves.
          divisor <- sqrt(summary$n.used)
          divisor[divisor == 0] <- NA_real_
          adjusted <- sweep(sweep(summary$pred, 2L, summary$mean, "-"), 2L, divisor, "/")
          adjusted <- sweep(adjusted, 2L, summary$mean, "+")
          yhat <- as.vector(adjusted)
        }
        else yhat <- summary$mean
        list(xvar.names = name, yhat = yhat, yhat.se = summary$se,
             yhat.mean = summary$mean, yhat.sd = summary$sd, n.used = summary$n.used,
             n.x = n.x, x.uniq = display, partial.values = grid, x = xx,
             factor.x = discrete)
      })
      names(prtl) <- xvar.names
    }
    plot.variable.obj <- list(family = family, partial = partial, event.info = event.info,
                 target = target, ylabel = ylabel, n = n, xvar.names = xvar.names,
                 nvar = nvar, plots.per.page = plots.per.page, granule = granule,
                 smooth.lines = smooth.lines)
    if (partial) {
      plot.variable.obj$pData <- prtl
      plot.variable.obj$plotthis <- lapply(prtl, function(pp) {
        data.frame(x = if (pp$factor.x) rep(pp$x.uniq, each = n) else pp$x.uniq,
                   yhat = pp$yhat)
      })
      plot.variable.obj$partial.summary <- lapply(prtl, function(pp) {
        data.frame(x = pp$x.uniq, mean = pp$yhat.mean, sd = pp$yhat.sd,
                   se = pp$yhat.se, n.used = pp$n.used)
      })
    }
    else {
      plot.variable.obj$yhat <- yhat
      plot.variable.obj$xvar <- xvar
    }
    class(plot.variable.obj) <- c("rfsrc", "plot.variable", family)
  }
  else {
    ## Replay uses the saved numerical values and layout settings.
    plot.variable.obj <- object
    family <- object$family
    partial <- object$partial
    event.info <- object$event.info
    target <- object$target
    ylabel <- object$ylabel
    n <- object$n
    if (missing(xvar.names)) xvar.names <- object$xvar.names else
      xvar.names <- intersect(xvar.names, object$xvar.names)
    if (!length(xvar.names)) stop("none of the supplied predictors match the saved plotting object")
    nvar <- length(xvar.names)
    plots.per.page <- min(.plot.variable.count(object$plots.per.page, "plots.per.page"), nvar)
    granule <- object$granule
    smooth.lines <- object$smooth.lines
    if (partial) prtl <- object$pData[xvar.names] else {
      yhat <- object$yhat
      xvar <- object$xvar[, xvar.names, drop = FALSE]
    }
  }
  ## --------------------------------------------------------------
  ## Rendering only. Restore device settings even when a panel errors.
  ## --------------------------------------------------------------
  if (show.plots) {
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par), add = TRUE)
    par(mfrow = c(min(plots.per.page, ceiling(nvar / plots.per.page)), plots.per.page))
    if (!partial) {
      yplot <- if (is.factor(yhat)) as.numeric(yhat) else yhat
      for (name in xvar.names) {
        xx <- xvar[[name]]
        xlabel <- name
        yl <- ylabel
        distinct <- unique(xx[!is.na(xx)])
        if (!is.factor(xx) && length(distinct) > granule) {
          .plot.variable.panel(xx, yplot, dots, xlabel, yl)
          point.args <- .plot.variable.merge(list(pch = 16, cex = if (n > 500) .5 else .75),
                             dots[intersect(names(dots), c("pch", "col", "cex", "bg"))])
          if (grepl("surv", family)) {
            status <- event.info$cens
            event <- if (family == "surv-CR") event.info$event.type[target] else NULL
            use.event <- is.finite(status) & (if (is.null(event)) status != 0 else status == event)
            for (pt in list(which(use.event), which(is.finite(status) & status == 0))) {
              if (!length(pt)) next
              style <- point.args
              if (is.null(dots$col)) style$col <- if (all(status[pt] == 0)) 1 else 4
              do.call("points", .plot.variable.merge(list(x = xx[pt], y = yplot[pt]), style))
            }
          }
          ok <- is.finite(xx) & is.finite(yplot)
          args <- .plot.variable.merge(list(col = 2, lwd = 3, type = "l"),
                     dots[intersect(names(dots), c("col", "lty", "lwd", "type", "pch", "cex", "bg"))])
          if (sum(ok) > 1L && length(unique(xx[ok])) > 1L) {
            do.call("lines", .plot.variable.merge(lowess(xx[ok], yplot[ok]), args))
          }
          else if (any(ok)) do.call("points", .plot.variable.merge(list(x = xx[ok], y = yplot[ok]), point.args))
        }
        else {
          labels <- if (is.factor(xx)) levels(xx)[levels(xx) %in% as.character(distinct)] else sort(distinct)
          groups <- lapply(labels, function(label) yplot[which(xx == label)])
          .plot.variable.boxes(groups, format(labels, trim = TRUE, digits = 4), dots, xlabel, yl)
        }
      }
    }
    else {
      for (pp in prtl) {
        discrete <- if (!is.null(pp$factor.x)) pp$factor.x else is.factor(pp$x) || pp$n.x <= granule
        xlabel <- pp$xvar.names
        yl <- ylabel
        if (discrete) {
          values <- matrix(pp$yhat, nrow = n, ncol = pp$n.x)
          groups <- lapply(seq_len(pp$n.x), function(j) values[, j])
          .plot.variable.boxes(groups, format(pp$x.uniq, trim = TRUE, digits = 4), dots, xlabel, yl, partial = TRUE)
        }
        else {
          grid <- pp$x.uniq
          center <- pp$yhat
          se <- pp$yhat.se
          upper <- center + 2 * se
          lower <- center - 2 * se
          xx <- c(min(pp$x), grid, max(pp$x), grid, grid)
          yy <- c(NA_real_, center, NA_real_, upper, lower)
          .plot.variable.panel(xx, yy, dots, xlabel, yl)
          p.args <- .plot.variable.merge(list(pch = 16, cex = if (pp$n.x > 25) .5 else .75, col = 2),
                       dots[intersect(names(dots), c("pch", "cex", "col", "bg"))])
          style <- if (is.null(dots$type)) "b" else dots$type
          if (style %in% c("p", "b", "o")) do.call("points", .plot.variable.merge(list(x = grid, y = center), p.args))
          if (any(is.finite(se) & se > 0)) {
            band.args <- .plot.variable.merge(list(lty = 3, col = 2),
                           dots[intersect(names(dots), c("col", "lwd"))])
            .plot.variable.line(grid, upper, smooth.lines, band.args)
            .plot.variable.line(grid, lower, smooth.lines, band.args)
          }
          if (!(style %in% c("p", "n"))) {
            line.args <- .plot.variable.merge(list(lty = 2, lwd = 2, type = if (style == "b") "l" else style),
                          dots[intersect(names(dots), c("col", "lty", "lwd"))])
            .plot.variable.line(grid, center, smooth.lines, line.args)
          }
          rug(pp$x[is.finite(pp$x)], ticksize = .03)
        }
      }
    }
  }
  invisible(plot.variable.obj)
}
plot.variable <- plot.variable.rfsrc
