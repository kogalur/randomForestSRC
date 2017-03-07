##  **********************************************************************
##  **********************************************************************
##  
##    RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
##  
##    This program is free software; you can redistribute it and/or
##    modify it under the terms of the GNU General Public License
##    as published by the Free Software Foundation; either version 3
##    of the License, or (at your option) any later version.
##  
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##  
##    You should have received a copy of the GNU General Public
##    License along with this program; if not, write to the Free
##    Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
##    Boston, MA  02110-1301, USA.
##  
##    ----------------------------------------------------------------
##    Project Partially Funded By: 
##    ----------------------------------------------------------------
##    Dr. Ishwaran's work was funded in part by DMS grant 1148991 from the
##    National Science Foundation and grant R01 CA163739 from the National
##    Cancer Institute.
##  
##    Dr. Kogalur's work was funded in part by grant R01 CA163739 from the 
##    National Cancer Institute.
##    ----------------------------------------------------------------
##    Written by:
##    ----------------------------------------------------------------
##      Hemant Ishwaran, Ph.D.
##      Director of Statistical Methodology
##      Professor, Division of Biostatistics
##      Clinical Research Building, Room 1058
##      1120 NW 14th Street
##      University of Miami, Miami FL 33136
##  
##      email:  hemant.ishwaran@gmail.com
##      URL:    http://web.ccs.miami.edu/~hishwaran
##      --------------------------------------------------------------
##      Udaya B. Kogalur, Ph.D.
##      Adjunct Staff
##      Department of Quantitative Health Sciences
##      Cleveland Clinic Foundation
##      
##      Kogalur & Company, Inc.
##      5425 Nestleway Drive, Suite L1
##      Clemmons, NC 27012
##  
##      email:  ubk@kogalur.com
##      URL:    https://github.com/kogalur/randomForestSRC
##      --------------------------------------------------------------
##  
##  **********************************************************************
##  **********************************************************************


plot.rfsrc <- function (x, outcome.target = NULL, plots.one.page = TRUE, sorted = TRUE, verbose = TRUE, ...)
{
  sf.flag <- FALSE
  if (sum(inherits(x, c("rfsrc", "synthetic"), TRUE) == c(1, 2)) == 2) {
    if (sum(inherits(x, c("rfsrc", "synthetic", "oob"), TRUE) == c(1, 2, 3)) != 3) {
      sf.flag <- TRUE
      sf.message <- "OOB was not used for synthetic forests, error rates/VIMP will be unreliable"
    }
    x <- x$rfSyn
  }
  if (sum(inherits(x, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rfsrc", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(rfsrc, grow)' or '(rfsrc, predict)'")
  }
  outcome.target <- get.univariate.target(x, outcome.target)
  x <- coerce.multivariate(x, outcome.target)
  if (is.null(x$err.rate)) {
    stop("object is devoid of performance values")
  }
  if (is.null(x$importance)) {
    x$importance <- NA
  }
  if (all(is.na(x$err.rate)) & all(is.na(x$importance))) {
    stop("performance values are all NA")
  }
  if (x$tree.err == FALSE) {
    colnames.err.rate <- colnames(x$err.rate)
    x$err.rate <- cbind(x$err.rate)[x$ntree, ]
    x$err.rate <- matrix(x$err.rate, x$ntree, length(x$err.rate), byrow = TRUE)
    colnames(x$err.rate) <- colnames.err.rate
  }
  if (x$family == "surv-CR" | x$family == "surv-CR") {
    x$yvar.names <- ""
  }
  old.par <- par(no.readonly = TRUE)
  cex <- par("cex")
  on.exit(par(old.par))
  if (all(is.na(x$importance))) {
    if (x$ntree > 1 && !all(is.na(x$err.rate))) {
      err <- cbind(x$err.rate)      
      par(cex = cex, mfrow = c(1,1))
      plot.err(err, x$yvar.names)    
    }
  }
    else {
      err <- cbind(x$err.rate)
      imp <- cbind(x$importance)
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
        plot.err(err, x$yvar.names)
      }
      if (ncol(imp) > 1) {
        imp.out <- imp[rev(pred.order),, drop = FALSE]
        dotChart(imp[pred.order,, drop = FALSE], x$yvar.names, dotchart.labels, cex = cex)
      }
      if (ncol(imp) == 1) {
        dotChart(imp[pred.order, ], x$yvar.names, dotchart.labels, cex = cex)
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
  if (sf.flag) {
    message(sf.message)
  }
}
plot.err <- function(err, yname = NULL, ...) {
  opar <- par("cex")
  on.exit(par(opar))
  matplot(1:nrow(err), err,
          xlab = "Number of Trees",
          ylab = paste("Error Rate:", yname),
          type = c("p", "l")[1 + 1 * (nrow(err) > 1)], pch = 16, lty = 1, lwd = 3)
  if (ncol(err) > 1) {
    legend("topright",
           legend = colnames(err), col = 1:ncol(err), lty = 1, lwd = 3)
  }
}
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
  y.dot <- dot.chart.main(x, labels = labels, xlab = paste("Variable Importance:", yname),
                          cex = cex, pch="", lwd = 2, lcolor = "white", gcolor = gcolor)
  segments(rep(max(0, min(x.dot, na.rm = TRUE)) - 1e-6, length(y.dot)),
           y.dot, x.dot, y.dot, col=c(2,4)[1 + 1 * (x.dot > 0)], lwd = 4)
  if (min(x.dot, na.rm = TRUE) < 0) abline(v=0, lwd = 2, lty = 2, col = 1)
}
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
