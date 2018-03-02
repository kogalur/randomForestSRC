plot.rfsrc <- function (x, m.target = NULL, plots.one.page = TRUE, sorted = TRUE, verbose = TRUE, ...)
{
  sf.flag <- FALSE
  ## is this a synthetic forest?  Printing is different in that case
  if (sum(inherits(x, c("rfsrc", "synthetic"), TRUE) == c(1, 2)) == 2) {
    if (sum(inherits(x, c("rfsrc", "synthetic", "oob"), TRUE) == c(1, 2, 3)) != 3) {
      sf.flag <- TRUE
      sf.message <- "OOB was not used for synthetic forests, error rates/VIMP will be unreliable"
    }
    x <- x$rfSyn
  }
  ## check that object is interpretable
  if (sum(inherits(x, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rfsrc", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(rfsrc, grow)' or '(rfsrc, predict)'")
  }
   
  ## coerce the (potentially) multivariate object if necessary.
  m.target <- get.univariate.target(x, m.target)
  x <- coerce.multivariate(x, m.target)
  ## grow objects under non-standard bootstrapping are devoid of
  ## performance values
  if (is.null(x$err.rate)) {
    stop("object is devoid of performance values")
  }
  ## set importance to NA if it is NULL
  if (is.null(x$importance)) {
    x$importance <- NA
  }
  ## return when everything is NA
  if (all(is.na(x$err.rate)) & all(is.na(x$importance))) {
    stop("performance values are all NA")
  }
  ## Check that the error rate vector is assigned for all trees.  If it is not,
  ## fill in the slots for 1:(ntree-1).  This will result in the plot being a flat line.
  if (x$tree.err == FALSE) {
    colnames.err.rate <- colnames(x$err.rate)
    x$err.rate <- cbind(x$err.rate)[x$ntree, ]
    x$err.rate <- matrix(x$err.rate, x$ntree, length(x$err.rate), byrow = TRUE)
    colnames(x$err.rate) <- colnames.err.rate
  }
  ## convert drc two-classifier error rate to one column
  if (!is.null(x$forest$perf.type) && (x$forest$perf.type == "g.mean" || x$forest$perf.type == "g.mean.rfq")) {
    x$err.rate <- x$err.rate[, 1, drop = FALSE]
  }
  ## for surv/CR we do not plot the yvar.names: we set them to blank
  if (x$family == "surv-CR" | x$family == "surv") {
    plot.yvar.names <- ""
  }
  else {
    plot.yvar.names <- paste("(", x$yvar.names, ")", sep = "")
  }
  ## save par for later restoration
  old.par <- par(no.readonly = TRUE)
  cex <- par("cex")
  on.exit(par(old.par))
  ## decide what plots to generate
  if (all(is.na(x$importance))) {
    if (x$ntree > 1 && !all(is.na(x$err.rate))) {
      err <- cbind(x$err.rate)      
      par(cex = cex, mfrow = c(1,1))
      plot.err(err, plot.yvar.names)
    }
  }
    else {
      ## convert err/vimp to matrix format
      err <- cbind(x$err.rate)
      imp <- cbind(x$importance)
      ## convert drc two-classifier vimp to one column
      if (!is.null(x$forest$perf.type) && (x$forest$perf.type == "g.mean" || x$forest$perf.type == "g.mean.rfq")) {
        imp <- imp[, 1, drop = FALSE]
      }
      x.var.names <- rownames(imp)
      n.pred <- nrow(imp)
      if (sorted) pred.order <- order(imp[, 1]) else pred.order <- n.pred:1
      if (ncol(imp) == 1) max.pred <- 100 else max.pred <- 80/ncol(imp)
      if (n.pred > max.pred) {
        dotchart.labels <- rep("",n.pred)
        pretty.pt <- pretty(1:n.pred, n = max.pred)
        dotchart.labels[pretty.pt] <- x.var.names[pred.order][pretty.pt]
      }
        else {
          dotchart.labels <- x.var.names[pred.order]
        }
      if (x$ntree > 1 & !all(is.na(x$err.rate)) & plots.one.page) {
        par(cex = cex, mfrow = c(1,2))
      }
        else {
          par(cex = cex, mfrow = c(1,1))
        }
      if (x$ntree > 1 & !all(is.na(x$err.rate))) {
        plot.err(err, plot.yvar.names)
      }
      ## CR/classification scenarios
      if (ncol(imp) > 1) {
        imp.out <- imp[rev(pred.order),, drop = FALSE]
        dotChart(imp[pred.order,, drop = FALSE], plot.yvar.names, dotchart.labels, cex = cex)
      }
      ## other scenarios
      if (ncol(imp) == 1) {
        dotChart(imp[pred.order, ], plot.yvar.names, dotchart.labels, cex = cex)
        if (!is.null(x$xvar.wt) & length(unique(x$xvar.wt)) > 1 ) {
          if (length(unique(x$xvar.wt)) == 1) x$xvar.wt <- 1
          imp.out <- as.data.frame(cbind(imp, imp/max(abs(imp), na.rm = TRUE), x$xvar.wt),
                                   row.names = x.var.names)[rev(pred.order),]
          if (nrow(imp.out) == 1) imp.out[1 , 2] <- 1
          colnames(imp.out) <- c("Importance","Relative Imp","xvar weight")
        }
          else {
            imp.out=as.data.frame(cbind(imp, imp/max(abs(imp), na.rm = TRUE)),
              row.names=x.var.names)[rev(pred.order),]
            if (nrow(imp.out) == 1) imp.out[1 , 2] <- 1
            colnames(imp.out) <- c("Importance","Relative Imp")
          }
      }
      cat("\n")
      if (verbose) {
        print(round(imp.out[1:min(n.pred, max.pred),, drop = FALSE],4), justify="right", print.gap=3)
      }    
    }
#################################################################################
  ##
  ## synthetic forest flag
  ##
################################################################################# 
  if (sf.flag) {
    message(sf.message)
  }
}
## error rate plot
plot.err <- function(err, yname = NULL, ...) {
  opar <- par("cex")
  on.exit(par(opar))
  matplot(1:nrow(err), err,
          xlab = "Number of Trees",
          ylab = paste("Error rate", yname),
          type = c("p", "l")[1 + 1 * (nrow(err) > 1)], pch = 16, lty = 1, lwd = 3)
  if (ncol(err) > 1) {
    legend("topright",
           legend = colnames(err), col = 1:ncol(err), lty = 1, lwd = 3)
  }
}
## pretty dotchart
dotChart <- function(x, yname = NULL, labels = NULL, cex = cex) {
  if (!is.null(dim(x))) {
    ncol  <- ncol(x)
    x.dot <- NULL
    for (k in ncol(x):1) {x.dot <- c(x.dot, x[, k])}
    gcolor <- 1:ncol
  }
    else {
      x.dot <- x
      gcolor <- par("fg")
    }
  y.dot <- dot.chart.main(x, labels = labels, xlab = paste("Variable Importance", yname),
                          cex = cex, pch="", lwd = 2, lcolor = "white", gcolor = gcolor)
  segments(rep(max(0, min(x.dot, na.rm = TRUE)) - 1e-6, length(y.dot)),
           y.dot, x.dot, y.dot, col=c(2,4)[1 + 1 * (x.dot > 0)], lwd = 4)
  if (min(x.dot, na.rm = TRUE) < 0) abline(v=0, lwd = 2, lty = 2, col = 1)
}
## workhorse for dotchart
dot.chart.main <- function (x, labels = NULL, groups = NULL, gdata = NULL, cex = NULL,
                            pch = 21, gpch = 21, bg = par("bg"), color = par("fg"), gcolor = par("fg"), 
                            lcolor = "gray", xlim = range(x[is.finite(x)]), main = NULL, 
                            xlab = NULL, ylab = NULL, ...) 
{
  opar <- par("mai", "mar", "cex", "yaxs")
  on.exit(par(opar))
  par(yaxs = "i")
  if (!is.numeric(x)) 
    stop("'x' must be a numeric vector or matrix")
  n <- length(x)
  if (is.matrix(x)) {
    if (is.null(labels)) 
      labels <- rownames(x)
    if (is.null(labels)) 
      labels <- as.character(1:nrow(x))
    labels <- rep(labels, length.out = n)
    if (is.null(groups)) 
      groups <- col(x, as.factor = TRUE)
    glabels <- levels(groups)
  }
    else {
      if (is.null(labels)) 
        labels <- names(x)
      glabels <- if (!is.null(groups)) 
                   levels(groups)
    }
  plot.new()
  linch <- if (!is.null(labels)) 
             max(strwidth(labels, "inch"), na.rm = TRUE)
             else 0
  if (is.null(glabels)) {
    ginch <- 0
    goffset <- 0
  }
    else {
      ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
      goffset <- 0.4
    }
  if (!(is.null(labels) && is.null(glabels))) {
    nmai <- par("mai")
    nmai[2] <- nmai[4] + max(linch + goffset, ginch) + 0.1
    par(mai = nmai)
  }
  if (is.null(groups)) {
    o <- 1:n
    y <- o
    ylim <- c(0, n + 1)
  }
    else {
      o <- sort.list(as.numeric(groups), decreasing = TRUE)
      x <- x[o]
      groups <- groups[o]
      color <- rep(color, length.out = length(groups))[o]
      lcolor <- rep(lcolor, length.out = length(groups))[o]
      offset <- cumsum(c(0, diff(as.numeric(groups)) != 0))
      y <- 1:n + 2 * offset
      ylim <- range(0, y + 2)
    }
  plot.window(xlim = xlim, ylim = ylim, log = "")
  lheight <- par("csi")
  if (!is.null(labels)) {
    linch <- max(strwidth(labels, "inch"), na.rm = TRUE)
    loffset <- (linch + 0.1)/lheight
    labs <- labels[o]
    mtext(labs, side = 2, line = loffset, at = y, adj = 0, 
          col = color, las = 2, cex = cex, ...)
  }
  abline(h = y, lty = "dotted", col = lcolor)
  points(x, y, pch = pch, col = color, bg = bg)
  if (!is.null(groups)) {
    gpos <- rev(cumsum(rev(tapply(groups, groups, length)) + 
                         2) - 1)
    ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
    goffset <- (max(linch + 0.2, ginch, na.rm = TRUE) + 0.1)/lheight
    mtext(glabels, side = 2, line = goffset, at = gpos, adj = 0, 
          col = gcolor, las = 2, cex = cex, ...)
    if (!is.null(gdata)) {
      abline(h = gpos, lty = "dotted")
      points(gdata, gpos, pch = gpch, col = gcolor, bg = bg, 
             ...)
    }
  }
  axis(1)
  box()
  title(main = main, xlab = xlab, ylab = ylab, ...)
  invisible(y)
}
