plot.subsample.rfsrc <- function(x, alpha = .01, xvar.names,
                          standardize = TRUE, normal = TRUE, jknife = FALSE,
                          target, m.target = NULL, pmax = 75, main = "",
                          sorted = TRUE, show.plots = TRUE,
                          ...)
{
  ## Validate selectors before extracting or drawing stored intervals.
  .subsample.check.flag(standardize, "standardize")
  .subsample.check.flag(normal, "normal")
  .subsample.check.flag(jknife, "jknife")
  .subsample.check.flag(sorted, "sorted")
  .subsample.check.flag(show.plots, "show.plots")
  if (!is.numeric(pmax) || length(pmax) != 1L || is.na(pmax) ||
      pmax < 1 || (is.finite(pmax) && pmax != floor(pmax))) {
    stop("'pmax' must be a positive integer or Inf.", call. = FALSE)
  }
  ##--------------------------------------------------------------
  ##
  ## was subsampling or double-bootstrap used?
  ##
  ##--------------------------------------------------------------
  if (sum(c(grepl("rfsrc", class(x))), grepl("subsample", class(x))) == 2) {
    subsample <- TRUE
  }
  else if (sum(c(grepl("rfsrc", class(x))), grepl("bootsample", class(x))) == 2) {
    subsample <- FALSE
  }
  else {
    stop("object must be obtained from call to 'subsample' function")
  }
  ##--------------------------------------------------------------
  ##
  ## coerce the (potentially) multivariate rf object
  ##
  ##--------------------------------------------------------------
  m.target <- get.univariate.target(x$rf, m.target)
  x$rf <- coerce.multivariate(x$rf, m.target)
  ##--------------------------------------------------------------
  ##
  ## family specific details
  ## - set the target if not supplied
  ## - assign pretty labels for the horizontal axis
  ##
  ##--------------------------------------------------------------
  fmly <- x$rf$family
  ## labels
  if (standardize) {
    if (fmly == "regr") {
      xlab <- "standardized vimp"
    }
    else {
      xlab <- "vimp"
    }
  }
  else {
    xlab <- "vimp"
  }
  ## over-ride x label when only generalization error is requested
  if (x$performance.only) {
    xlab <- "generalization error"
  }
  ## The selected statistic column supplies both values and its label,
  ## including error-only classification and competing-risk displays.
  if (is.null(m.target)) {
    vmp <- x$vmp[[1]]
  }
  else {
    vmp <- x$vmp[[m.target]]
  }
  if (is.null(vmp) || NROW(vmp) == 0L || NCOL(vmp) == 0L) {
    stop("No reference statistics are available for the selected outcome.",
         call. = FALSE)
  }
  vmp <- as.matrix(vmp)
  vmp.col.names <- colnames(vmp)
  if (is.null(vmp.col.names) && ncol(vmp) == 1L) {
    vmp.col.names <- if (is.null(m.target)) names(x$vmp)[1L] else m.target
  }
  if (fmly == "regr" || fmly == "surv") {
    target <- 0
    vmp <- c(vmp[, 1L])
    xlab <- paste(xlab, " (", vmp.col.names, ")", sep = "")
  }
  else if (fmly == "class") {
    if (missing(target)) {
      target <- 0
    }
    else {
      yvar.levels <- levels(x$rf$yvar)
      if (is.character(target)) {
        target <- match(match.arg(target, yvar.levels), yvar.levels)
      }
      else {
        if (!is.numeric(target) || length(target) != 1L || !is.finite(target) ||
            target != floor(target) || target < 0 || target > length(yvar.levels)) {
          stop("'target' is specified incorrectly")
        }
      }
    }
    vmp <- c(vmp[, 1 + target])
    xlab <- paste(xlab, " (", vmp.col.names[1 + target], ")", sep = "")
  }
  else if (fmly == "surv-CR") {    
    if (missing(target)) {
      target <- 0
    }
    else {
      n.event <- length(get.event.info(x$rf)$event.type)
      if (!is.numeric(target) || length(target) != 1L || !is.finite(target) ||
          target != floor(target) || target < 1 || target > n.event) {
        stop("'target' is specified incorrectly")
      }
      target <- target - 1
    }
    vmp <- c(vmp[, 1 + target])
    xlab <- paste(xlab, " (", vmp.col.names[1 + target], ")", sep = "")
  }
  ##--------------------------------------------------------------
  ##
  ## extract necessary objects
  ##
  ##--------------------------------------------------------------
  if (subsample) {
    oo <- extract.subsample(x, alpha = alpha, target = target, m.target = m.target,
                            standardize = standardize, raw = TRUE)
    boxplot.dta <- oo$boxplot.dta
  }
  else {
    oo <- extract.bootsample(x, alpha = alpha, target = target, m.target = m.target,
                             standardize = standardize, raw = TRUE)
    boxplot.dta <- oo$ci
  }
  if (normal) {
    if (subsample && jknife) {
      ci <- oo$ci.jk.Z
    }
    else {
      ci <- oo$ci.Z
    }
  }
  else {
    ci <- oo$ci
  }
  ##--------------------------------------------------------------
  ##
  ## final processing
  ##
  ##--------------------------------------------------------------
  ## (a) user has provided a list of target variables
  if (!missing(xvar.names)) {
    trim.pt <- colnames(boxplot.dta) %in% xvar.names
    if (sum(trim.pt) > 0) {
      boxplot.dta <- boxplot.dta[, trim.pt, drop = FALSE]
      ci <- ci[, trim.pt, drop = FALSE]
      vmp <- vmp[trim.pt]
    }
    else {
      stop("user supplied variable names do not match any of the variables from analysis")
    }
  }
  ## (b) sort by vimp?
  p <- ncol(boxplot.dta)
  if (sorted || p > pmax) {
    o.pt <- order(vmp, decreasing = TRUE)
  }
  else {
    o.pt <- 1:p
  }
  ## (c) trim if too many variables
  if (p > pmax) {
    o.pt <- o.pt[1:pmax]
  }
  ## redefine the data to match the processing
  boxplot.dta <- boxplot.dta[, o.pt, drop = FALSE]
  ci <- ci[, o.pt, drop = FALSE]
  ##--------------------------------------------------------------
  ##
  ## skeleton boxplot 
  ##
  ##--------------------------------------------------------------
  bp <- boxplot(boxplot.dta,
                yaxt="n",
                outline = FALSE,
                horizontal = TRUE,
                plot = FALSE)
  bp$stats <- ci
  col.pt <- bp$stats[1, ] > 0 
  col.pt[is.na(col.pt)] <- FALSE
  colr <- c("blue", "red")[1 + col.pt]
  ##--------------------------------------------------------------
  ## Resolve graphical settings once, then reuse them for both drawing
  ## passes. User arguments override defaults instead of being duplicated.
  ##--------------------------------------------------------------
  merge.args <- function(defaults, overrides) {
    if (length(overrides)) {
      nms <- names(overrides)
      if (is.null(nms) || anyNA(nms) || any(!nzchar(nms))) {
        stop("Graphical arguments must be named.", call. = FALSE)
      }
      overrides <- overrides[!duplicated(nms, fromLast = TRUE)]
      defaults[names(overrides)] <- overrides
    }
    defaults
  }
  dots <- list(...)
  pars <- dots$pars
  dots$pars <- NULL
  if (!is.null(pars) && !is.list(pars)) {
    stop("'pars' must be a named list of graphical arguments.", call. = FALSE)
  }
  dots <- merge.args(merge.args(list(), pars), dots)
  if ("z" %in% names(dots)) {
    stop("The interval summaries supplied to bxp cannot be replaced by 'z'.",
         call. = FALSE)
  }
  ## Custom row names follow the selected pre-sort rows, as in the original
  ## interface. A scalar ylab is an axis title, not a one-row label override.
  if (!is.null(dots$ylab) && length(dots$ylab) > 1L && length(dots$ylab) == p) {
    bp$names <- dots$ylab[o.pt]
    dots$ylab <- NULL
  }
  if (!is.null(dots$names)) {
    if (length(dots$names) != p) {
      stop("'names' must have one label per selected row before sorting and trimming.",
           call. = FALSE)
    }
    bp$names <- dots$names[o.pt]
    dots$names <- NULL
  }
  if (x$performance.only) bp$names <- rep("", length(bp$names))
  horizontal <- if (is.null(dots$horizontal)) TRUE else dots$horizontal
  .subsample.check.flag(horizontal, "horizontal")
  axes <- if (is.null(dots$axes)) TRUE else dots$axes
  .subsample.check.flag(axes, "axes")
  show.names <- if (is.null(dots$show.names)) TRUE else dots$show.names
  .subsample.check.flag(show.names, "show.names")
  fill <- if (is.null(dots$col)) colr else dots$col
  defaults <- list(z = bp, main = main,
                   boxfill = fill, whisklty = 1, whisklwd = 2,
                   outline = FALSE, horizontal = horizontal,
                   frame.plot = axes)
  defaults[[if (horizontal) "xlab" else "ylab"]] <- xlab
  args <- merge.args(defaults, dots)
  ## Public xlim/ylim refer to the axes as displayed. bxp internally reverses
  ## the limit arguments for horizontal boxes, so translate them here.
  if (horizontal) {
    xy <- args[c("xlim", "ylim")]
    args[c("xlim", "ylim")] <- NULL
    if (!is.null(xy$ylim)) args$xlim <- xy$ylim
    if (!is.null(xy$xlim)) args$ylim <- xy$xlim
  }
  ## Manual axes allow separate suppression and the original row-label
  ## orientation. Both bxp passes have axes disabled to avoid drawing twice.
  xaxt <- if (is.null(dots$xaxt)) "s" else dots$xaxt
  yaxt <- if (is.null(dots$yaxt)) "s" else dots$yaxt
  args$axes <- FALSE
  args$show.names <- FALSE
  if (is.null(args$at)) args$at <- seq_along(bp$names)
  if (!is.numeric(args$at) || length(args$at) != length(bp$names) ||
      any(!is.finite(args$at))) {
    stop("'at' must give one finite position per displayed interval.", call. = FALSE)
  }
  ## Common axis styling is forwarded without replacing the tick positions
  ## or statistic labels. bxp-specific controls remain with bxp.
  axis.names <- c("las", "cex.axis", "col.axis", "font.axis", "line",
                  "lwd", "lwd.ticks", "col", "col.ticks", "hadj", "padj",
                  "gap.axis", "mgp", "tck", "tcl", "xpd", "fg")
  axis.args <- dots[intersect(names(dots), axis.names)]
  if (show.plots) {
    if (is.null(dots$xaxt)) xaxt <- par("xaxt")
    if (is.null(dots$yaxt)) yaxt <- par("yaxt")
    do.call("bxp", args)
    value.side <- if (horizontal) 1L else 2L
    row.side <- if (horizontal) 2L else 1L
    draw.value <- axes && !identical(if (horizontal) xaxt else yaxt, "n")
    draw.rows <- axes && show.names &&
      !identical(if (horizontal) yaxt else xaxt, "n")
    if (draw.value) {
      value.axis <- merge.args(list(side = value.side), axis.args)
      do.call("axis", value.axis)
    }
    if (draw.rows) {
      row.axis <- merge.args(list(side = row.side, at = args$at,
                                  labels = bp$names, las = 2, tick = FALSE),
                             axis.args)
      do.call("axis", row.axis)
    }
    if (horizontal) {
      abline(h = args$at, col = gray(.9), lty = 1)
      if (!par("xlog")) abline(v = 0, lty = 1, lwd = 1.5, col = gray(.8))
    }
    else {
      abline(v = args$at, col = gray(.9), lty = 1)
      if (!par("ylog")) abline(h = 0, lty = 1, lwd = 1.5, col = gray(.8))
    }
    ## Redraw the intervals above the guides with exactly the resolved styles.
    ## Do not repeat titles or frames on the overlay.
    args$add <- TRUE
    args$ann <- FALSE
    args$frame.plot <- FALSE
    do.call("bxp", args)
  }
  invisible(bp)
}
plot.subsample <- plot.subsample.rfsrc
