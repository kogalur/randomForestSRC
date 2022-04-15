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
  ##   
  ##   coherence checks
  ##
  ## --------------------------------------------------------------
  ## is this a synthetic forest?
  if (sum(inherits(x, c("rfsrc", "synthetic"), TRUE) == c(1, 2)) == 2) {
    x <- x$rfSyn
  }
  ## check that object is interpretable
  ## first rename x to object to avoid confusion with x matrix
  object <- x
  remove(x)
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(object, c("rfsrc", "plot.variable"), TRUE) == c(1,2)) != 2) {
    stop("this function only works for objects of class `(rfsrc, grow)', '(rfsrc, synthetic)' or '(rfsrc, plot.variable)'")
  }
  if (object$family == "unsupv") {
    stop("this function does not apply to unsupervised forests")
  }
  if (inherits(object, "anonymous")) {
      stop("this function does work with anonymous forests")
  }    
  ## terminate if a partial plot is requested but no forest is found
  if (partial && (is.null(object$forest) || !object$forest$forest)) {
    stop("forest is empty:  re-run rfsrc (grow call) with forest=TRUE")
  }
  ## --------------------------------------------------------------
  ##   
  ##   preliminary processing
  ##
  ## --------------------------------------------------------------
  ## if missing data was imputed then overlay the missing data for x
  ## bug reported by John Ehrlinger
  xvar <- object$xvar
  if (!is.null(object$imputed.indv)) {
    xvar[object$imputed.indv, ] <- object$imputed.data[, object$xvar.names]
  }
  n <- nrow(xvar)
  ## obtain hidden options - set matched options to NULL
  dots <- list(...)
  plot.names <- names(formals(plot.default))
  bxp.names <- unique(c(names(formals(bxp)),
                 "xaxt", "yaxt", "las", "cex.axis", 
                 "col.axis", "main", "cex.main",
                 "col.main", "sub", "cex.sub", "col.sub", 
                 "xlab", "ylab", "cex.lab", "col.lab", "boxfill","varwidth"))
  axis.names <- unique(c(names(formals(axis)), "las", "cex.axis"))
  ## by default xlab and ylab are extracted from the data
  ## over-ridden only by hidden options
  xlab.flag <- ylab.flag <- TRUE
  if (!is.null(dots$xlab)) {
    xlab.flag <- FALSE
  }
  if (!is.null(dots$ylab)) {
    ylab.flag <- FALSE
  }
  ## default notch, boxfill
  if (is.null(dots$notch)) {
    dots$notch <- TRUE
  }
  if (is.null(dots$boxfill)) {
    dots$boxfill <- "bisque"
  }
  ## --------------------------------------------------------------------------------------------
  ## do the following if the class is NOT "plot.variable" 
  ## (i.e. the object has not already been processed by the wraper)
  if (!inherits(object, "plot.variable")) {
    ## process the subsetted index 
    ## assumes the entire data set is to be used if not specified
    if (missing(subset)) {
      subset <- 1:n
    }
      else {
        ## convert the user specified subset into a usable form
        if (is.logical(subset)) subset <- which(subset)
        subset <- unique(subset[subset >= 1 & subset <= n])
        if (length(subset) == 0) {
          stop("'subset' not set properly")
        }
      }
    ## subset the x-data
    xvar <- xvar[subset,, drop = FALSE]
    n <- nrow(xvar)
    ## process the object depending on the underlying family
    family <- object$family
    ## Ensure the coherency of multivariate target.
    m.target <- get.univariate.target(object, m.target)
    ## Survival and competing risk families.
    if (grepl("surv", family)) {
      ##extract event information
      event.info <- get.event.info(object, subset)
      cens <- event.info$cens
      event.type <- event.info$event.type
      ## assign time if missing
      if (missing(time)) {
        time <- median(event.info$time.interest, na.rm = TRUE)
      }
      ## Check for single time point.
      if (length(time) > 1) {
        stop("time must be a single value:  ", time)
      }
      ## CR analysis.
      if (family == "surv-CR") {
        if (missing(target)) {
          target <- 1
        }
          else {
            if (target < 1 || target > max(event.type, na.rm = TRUE)) {
              stop("'target' is specified incorrectly")
            }
          }
        VIMP <- object$importance[, target]
        ## set the surv.type
        surv.type <- setdiff(surv.type, c("mort", "rel.freq", "surv"))[1]
        pred.type <- match.arg(surv.type, c("years.lost", "cif", "chf"))
        ylabel <- switch(pred.type,
                         "years.lost" = paste("Years lost for event ", target),
                         "cif"        = paste("CIF for event ", target, " (time=", time, ")", sep = ""),
                         "chf"        = paste("CHF for event ", target, " (time=", time, ")", sep = ""))
      }
      ## Right-censoring analysis.
        else {
          target <- 1
          VIMP <- object$importance
          ## set the surv.type
          surv.type <- setdiff(surv.type, c("years.lost", "cif", "chf"))[1]
          pred.type <- match.arg(surv.type, c("rel.freq", "mort", "chf", "surv"))
          ylabel <- switch(pred.type,
                           "rel.freq"   = "standardized mortality",
                           "mort"       = "mortality",
                           "chf"        = paste("CHF (time=", time, ")", sep = ""),
                           "surv"       = paste("predicted survival (time=", time, ")", sep = ""))
        }
    }
    ## Univariate and multivariate families.
      else {
        ## assign a null time value
        event.info <- time <- NULL
        ## Factor outcome
        if(is.factor(coerce.multivariate(object, m.target)$yvar)) {
          object.yvar.levels <- levels(coerce.multivariate(object, m.target)$yvar)
          pred.type <- match.arg(class.type, c("prob", "bayes"))
          if (missing(target)) {
            target <- object.yvar.levels[1]
          }
          if (is.character(target)) {
            target <- match(match.arg(target, object.yvar.levels), object.yvar.levels)
          }
            else {
              if ((target > length(object.yvar.levels)) | (target < 1)) {
                stop("target is specified incorrectly:", target)
              }
            }
          if (pred.type == "prob") {
            VIMP <- coerce.multivariate(object, m.target)$importance[, 1 + target]
            ylabel <- paste("probability", object.yvar.levels[target])
            remove(object.yvar.levels)
          }
            else {
              VIMP <- coerce.multivariate(object, m.target)$importance[, 1]
              ylabel <- paste("bayes")
              remove(object.yvar.levels)
            }
        }
        ## Regression outcome.
          else {
            pred.type <- "y"
            target <- NULL
            VIMP <- coerce.multivariate(object, m.target)$importance
            ylabel <- expression(hat(y))
          }
      }
    ## extract x-var names to be plotted
    ## should x-var be sorted by importance?
    if (missing(xvar.names)) {
      xvar.names <- object$xvar.names
    }
      else {
        xvar.names <- intersect(xvar.names, object$xvar.names)
        if (length(xvar.names) ==  0){
          stop("none of the x-variable supplied match available ones:\n", object$xvar.names)
        }
      }
    if (sorted & !is.null(VIMP)) {
      xvar.names <- xvar.names[rev(order(VIMP[xvar.names]))]
    }
    if (!missing(nvar)) {
      nvar <- max(round(nvar), 1)
      xvar.names <- xvar.names[1:min(length(xvar.names), nvar)]
    }
    nvar <- length(xvar.names)
    if (!partial) {
      ## Marginal plot setup only.
      yhat <- extract.pred(predict(object, m.target = m.target),
                           pred.type,
                           subset,
                           time,
                           m.target,
                           target,
                           oob = oob)
    }
      else {
        ## Partial plot set up only.
        if (npts < 1) npts <- 1 else npts <- round(npts)
        ## Loop over all x-variables.
        prtl <- lapply(1:nvar, function(k) {
          ## We now allow subsetting of the x-values in the partial mode call.
          ## x <- na.omit(object$xvar[, object$xvar.names == xvar.names[k]])
          ## Bug reported by Amol Pande 16/11/2017
          x <- na.omit(xvar[, object$xvar.names == xvar.names[k]])
          if (is.factor(x)) x <- factor(x, exclude = NULL)          
          n.x <- length(unique(x))
          if (!is.factor(x) & n.x > npts) {
            x.uniq <- sort(unique(x))[unique(as.integer(seq(1, n.x, length = min(npts, n.x))))]
          }
          else {
            if (is.factor(x)) {
              x.uniq <- 1:length(unique(x))##partial.rfsrc requires factors as integers
            }
            else {
              x.uniq <- sort(unique(x))
            }
          }
          n.x <- length(x.uniq)
          yhat <- yhat.se <- NULL
          factor.x <- is.factor(x) | (n.x <= granule)
          pred.temp <- extract.partial.pred(partial.rfsrc(object$forest,
                                                          m.target = m.target,
                                                          partial.type = pred.type,
                                                          partial.xvar = xvar.names[k],
                                                          partial.values = x.uniq,
                                                          partial.time = time,
                                                          oob = oob),
                                            pred.type,
                                            ## 1:n,
                                            subset,##we now allow subsetting
                                            m.target,
                                            target)
          ## Results in the mean along an x-value over n.
          if (!is.null(dim(pred.temp))) {
            mean.temp <- colMeans(pred.temp, na.rm = TRUE)
          }
          else {
            mean.temp <- mean(pred.temp, na.rm = TRUE)
          }
          if (!factor.x) {
            yhat <- mean.temp
            if (coerce.multivariate(object, m.target)$family == "class") {
              yhat.se <- mean.temp * (1 - mean.temp) / sqrt(n)
            }
            else {
              sd.temp <- apply(pred.temp, 2, sd, na.rm = TRUE)
              yhat.se <- sd.temp / sqrt(n)                
            }
          }
          else {
            if (!is.null(dim(pred.temp))) {
              #pred.temp <- t(mean.temp + (t(pred.temp) - mean.temp) / sqrt(n))
              pred.temp <- t(apply(pred.temp, 1, function(yhat) {
                se <- (yhat - mean.temp) / sqrt(n)
                mean.temp + sign(se) * abs(se)
              }))
            }
            else {
              pred.temp <- mean.temp + (pred.temp - mean.temp) / sqrt(n)
            }
            yhat <- c(yhat, pred.temp)
            x.uniq <- sort(unique(x))##map factor back to original labels
          }
          list(xvar.names = xvar.names[k], yhat = yhat, yhat.se = yhat.se, n.x = n.x, x.uniq = x.uniq, x = x)
        })
    }
    plots.per.page <- max(round(min(plots.per.page,nvar)), 1)
    granule <- max(round(granule), 1)
    ## save the plot.variable object
    plot.variable.obj <- list(family = family,
                              partial = partial,
                              event.info = event.info,
                              target = target,
                              ylabel = ylabel,
                              n = n,
                              xvar.names = xvar.names,
                              nvar = nvar, 
                              plots.per.page = plots.per.page,
                              granule = granule,
                              smooth.lines = smooth.lines)
    if (partial) {
      plot.variable.obj$pData <- prtl
    }
      else {
        plot.variable.obj$yhat <- yhat
        plot.variable.obj$xvar <- xvar
      }
    ## assign the class
    class(plot.variable.obj) <- c("rfsrc", "plot.variable", family)
  }
  ## ---------------------------------------------------------------------------------
  ## TBD2 Currently, the plot.variable family is not implemented!
  ## Just pull the precomputed variables ...
    else {
      plot.variable.obj <- object
      remove(object)
      family <- plot.variable.obj$family
      partial <- plot.variable.obj$partial
      event.info <- plot.variable.obj$event.info
      target <- plot.variable.obj$target
      ylabel <- plot.variable.obj$ylabel
      n <- plot.variable.obj$n
      xvar.names <- plot.variable.obj$xvar.names
      nvar <- plot.variable.obj$nvar 
      plots.per.page <- plot.variable.obj$plots.per.page
      granule <- plot.variable.obj$granule
      smooth.lines <- plot.variable.obj$smooth.lines
      if (partial) {
        prtl <- plot.variable.obj$pData
      }
        else {
          yhat <- plot.variable.obj$yhat
          xvar <- plot.variable.obj$xvar
        }
      if (!is.null(event.info)){
        cens <- event.info$cens
        event.type <- event.info$event.type
      }
    }
  ## save par settings
  if (show.plots) {
    old.par <- par(no.readonly = TRUE)
  }
  ## --------------------------------------------------------------------------------
  ##
  ## marginal plots
  ##
  ## --------------------------------------------------------------------------------
  if (!partial && show.plots) {
    ## graphical set up
    par(mfrow = c(min(plots.per.page, ceiling(nvar / plots.per.page)), plots.per.page))
    if (n > 500) cex.pt <- 0.5 else cex.pt <- 0.75
    ## loop over variables
    for (k in 1:nvar) {
      ##x-stuff
      x <- xvar[, which(colnames(xvar) == xvar.names[k])]
      x.uniq <- unique(x)
      n.x <- length(x.uniq)
      ## hidden options
      if (xlab.flag) {
        dots$xlab <- xvar.names[k]
      }
      if (ylab.flag) {
        dots$ylab <- ylabel
      }
      ## x-continuous-plot
      if (!is.factor(x) & n.x > granule) {
        do.call("plot", c(list(x = x, y = yhat, type = "n"), dots[names(dots) %in% plot.names]))
        if (grepl("surv", family)) {
          points(x[cens == target], yhat[cens == target], pch = 16, col = 4, cex = cex.pt)
          points(x[cens == 0], yhat[cens == 0], pch = 16, cex = cex.pt)
        }
        lines(lowess(x[!is.na(x)], yhat[!is.na(x)]), col = 2, lwd=3)
      }
      ## x-factor-plots
      else {
        if (is.factor(x)) x <- factor(x, exclude = NULL)
        bp <- boxplot(yhat ~ x, na.action = "na.omit", names = rep("", n.x), plot = FALSE)
        do.call("bxp", c(list(z = bp, outline = FALSE, xaxt = "n"), dots[names(dots) %in% bxp.names]))
        at.pretty <- unique(round(pretty(1:n.x, min(30, n.x))))
        at.pretty <- at.pretty[at.pretty >= 1 & at.pretty <= n.x]
        do.call("axis", c(list(side = 1, at = at.pretty,
             labels = format(sort(x.uniq)[at.pretty], trim = TRUE, digits = 4),
             tick = TRUE), dots[names(dots) %in% axis.names]))
        }
    }
  }
  ## --------------------------------------------------------------------------------
  ##
  ## partial plots
  ##
  ## --------------------------------------------------------------------------------
  if (partial && show.plots) {
    ## graphical setup
    plots.per.page <- max(round(min(plots.per.page,nvar)), 1)
    granule <- max(round(granule),1)
    par(mfrow = c(min(plots.per.page, ceiling(nvar/plots.per.page)), plots.per.page))
    ## loop over variables
    for (k in 1:nvar) {
      ## x-stuff
      x <- prtl[[k]]$x
      if (is.factor(x)) x <- factor(x, exclude = NULL)          
      x.uniq <- prtl[[k]]$x.uniq
      n.x <- prtl[[k]]$n.x
      if (n.x > 25) cex.pt <- 0.5 else cex.pt <- 0.75
      ## y-stuff
      yhat <- prtl[[k]]$yhat
      yhat.se <- prtl[[k]]$yhat.se
      factor.x <- !(!is.factor(x) & (n.x > granule))
      ## hidden options
      if (xlab.flag) {
        dots$xlab <- prtl[[k]]$xvar.names
      }
      if (ylab.flag) {
        dots$ylab <- ylabel
      }
      ## x-continuous-plot
      if (!factor.x) {
        do.call("plot", c(list(x = c(min(x), x.uniq, max(x), x.uniq, x.uniq),
                               y = c(NA, yhat, NA, yhat + 2 * yhat.se, yhat - 2 * yhat.se),
                               type = "n"), dots[names(dots) %in% plot.names]))
        points(x.uniq, yhat, pch = 16, cex = cex.pt, col = 2)
        if (any(yhat.se > 0)) {
            if (smooth.lines && !any(is.na(yhat.se))) {
                lines(lowess(x.uniq, yhat + 2 * yhat.se), lty = 3, col = 2)
                lines(lowess(x.uniq, yhat - 2 * yhat.se), lty = 3, col = 2)
            }
            else {
                lines(x.uniq, yhat + 2 * yhat.se, lty = 3, col = 2)
                lines(x.uniq, yhat - 2 * yhat.se, lty = 3, col = 2)
            }
        }
        if (smooth.lines) {
          lines(lowess(x.uniq, yhat), lty = 2, lwd=2)
        }
          else {
            lines(x.uniq, yhat, lty = 2, lwd=2)
          }
        rug(x, ticksize=0.03)
      }
      ## x-factor-plots
        else {
          y.se <- 0.005
          bp <- boxplot(yhat ~ rep(x.uniq, rep(n, n.x)), names = rep("", n.x), plot = FALSE)
          do.call("bxp", c(list(z = bp, outline = FALSE, range = 2,
                                ylim = c(min(bp$stats[1,], na.rm = TRUE) * ( 1 - y.se ),
                                         max(bp$stats[5,], na.rm = TRUE) * ( 1 + y.se )),
                                xaxt = "n"), dots[names(dots) %in% bxp.names]))
          at.pretty <- unique(round(pretty(1:n.x, min(30,n.x))))
          at.pretty <- at.pretty[at.pretty >= 1 & at.pretty <= n.x]
          do.call("axis", c(list(side = 1, at = at.pretty,
             labels = format(sort(x.uniq)[at.pretty], trim = TRUE, digits = 4),
             tick = TRUE), dots[names(dots) %in% axis.names]))
      }
    }
  }
  ## Restore par settings
  if (show.plots) {
    par(old.par)
  }
  ## Return the plot.variable object for reuse
  invisible(plot.variable.obj)
}
plot.variable <- plot.variable.rfsrc
