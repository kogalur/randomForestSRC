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
  ## is this a synthetic forest?
  if (sum(inherits(x, c("rfsrc", "synthetic"), TRUE) == c(1, 2)) == 2) {
    x <- x$rfSyn
  }
  ## check that object is interpretable
  ## first rename x to object to avoid confusion with x matrix
  object <- x
  remove(x)
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(object, c("rfsrc", "predict"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(object, c("rfsrc", "plot.variable"), TRUE) == c(1,2)) != 2) {
    stop("this function only works for objects of class `(rfsrc, grow)', '(rfsrc, predict)' or '(rfsrc, plot.variable)'")
  }
  if (object$family == "unsupv") {
    stop("this function does not apply to unsupervised forests")
  }
  ## terminate if a partial plot is requested but no forest is found
  if (partial && is.null(object$forest)) {
    stop("forest is empty:  re-run rfsrc (grow) call with forest=TRUE")
  }
  ## if missing data was imputed then overlay the missing data for x
  ## bug reported by John Ehrlinger
  xvar <- object$xvar
  if (!is.null(object$imputed.indv)) {
    xvar[object$imputed.indv, ] <- object$imputed.data[, object$xvar.names]
  }
  n <- nrow(xvar)
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
        ## Check for a wierd NULL surv.type.  This just makes life easier later.
      if (is.null(surv.type)) {
        stop("surv.type is specified incorrectly:  ", surv.type)
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
        ## Pick the first occurrence.  This allows the omission of a value for surv.type.
        pred.type <- setdiff(surv.type, c("rel.freq", "mort", "chf", "surv"))[1]
        pred.type <- match.arg(pred.type, c("years.lost", "cif", "chf"))
        ylabel <- switch(pred.type,
                         "years.lost" = paste("Years lost for event ", target),
                         "cif"        = paste("CIF for event ", target, " (time=", time, ")", sep = ""),
                         "chf"        = paste("CHF for event ", target, " (time=", time, ")", sep = ""))
      }
      ## Right-censoring analysis.
        else {
          target <- 1
          VIMP <- object$importance
          ## Pick the first occurrence.  This allows the omission of a value for surv.type.
          pred.type <- setdiff(surv.type, c("years.lost", "cif", "chf"))[1]
          pred.type <- match.arg(pred.type, c("rel.freq", "mort", "chf", "surv"))
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
        ## Factor outcome.
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
      yhat <- extract.pred(predict.rfsrc(object, m.target = m.target, importance = "none"),
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
              x.uniq <- sort(unique(x))
            }
          n.x <- length(x.uniq)
          yhat <- yhat.se <- NULL
          factor.x <- !(!is.factor(x) & (n.x > granule))
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
            mean.temp <- apply(pred.temp, 2, mean, na.rm = TRUE)
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
              pred.temp <- mean.temp + (pred.temp - mean.temp) / sqrt(n)
              yhat <- c(yhat, pred.temp)
            }
          list(xvar.name = xvar.names[k], yhat = yhat, yhat.se = yhat.se, n.x = n.x, x.uniq = x.uniq, x = x)
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
  ## (TBD TBD) Currently, the plot.variable family is not implemented!
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
      ## x-continuous-plot
      if (!is.factor(x) & n.x > granule) {
        plot(x,
             yhat,
             xlab = xvar.names[k],
             ylab = ylabel,
             type = "n", ...) 
        if (grepl("surv", family)) {
          points(x[cens == target], yhat[cens == target], pch = 16, col = 4, cex = cex.pt)
          points(x[cens == 0], yhat[cens == 0], pch = 16, cex = cex.pt)
        }
        lines(lowess(x[!is.na(x)], yhat[!is.na(x)]), col = 2, lwd=3)
      }
      ## x-factor-plots
        else {
          if (is.factor(x)) x <- factor(x, exclude = NULL)          
          boxplot(yhat ~ x, na.action = "na.omit",
                  xlab = xvar.names[k],
                  ylab = ylabel,
                  notch = TRUE,
                  outline = FALSE,
                  col = "bisque",
                  names = rep("", n.x),
                  xaxt = "n", ...)
          at.pretty <- unique(round(pretty(1:n.x, min(30, n.x))))
          at.pretty <- at.pretty[at.pretty >= 1 & at.pretty <= n.x]
          axis(1,
               at = at.pretty,
               labels = format(sort(x.uniq)[at.pretty], trim = TRUE, digits = 4),
               tick = TRUE)
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
      ## x-continuous-plot
      if (!factor.x) {
        plot(c(min(x), x.uniq, max(x), x.uniq, x.uniq),
             c(NA, yhat, NA, yhat + 2 * yhat.se, yhat - 2 * yhat.se),
             xlab = prtl[[k]]$xvar.name,
             ylab = ylabel,
             type = "n", ...)
        points(x.uniq, yhat, pch = 16, cex = cex.pt, col = 2)
        if (!is.na(yhat.se) && any(yhat.se > 0)) {
          if (smooth.lines) {
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
          bxp.call <- boxplot(yhat ~ rep(x.uniq, rep(n, n.x)), range = 2, plot = FALSE)
          boxplot(yhat ~ rep(x.uniq, rep(n, n.x)),
                  xlab = xvar.names[k],
                  ylab = ylabel,
                  notch = TRUE,
                  outline = FALSE,
                  range = 2,
                  ylim = c(min(bxp.call$stats[1,], na.rm=TRUE) * ( 1 - y.se ),
                    max(bxp.call$stats[5,], na.rm=TRUE) * ( 1 + y.se )),
                  col = "bisque",
                  names = rep("",n.x),
                  xaxt = "n", ...)
          at.pretty <- unique(round(pretty(1:n.x, min(30,n.x))))
          at.pretty <- at.pretty[at.pretty >= 1 & at.pretty <= n.x]
          axis(1,
               at = at.pretty,
               labels = format(sort(x.uniq)[at.pretty], trim = TRUE, digits = 4),
               tick = TRUE)
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
## extraction function when partial=FALSE
extract.pred <- function(obj, type, subset, time, m.target, target, oob = oob) {
    ## Coerce the (potentially) multivariate object if necessary.
    obj <- coerce.multivariate(obj, m.target)
    ## decide if OOB or in-bag values are requested
    if (oob == FALSE) {
        pred <- obj$predicted
        surv <- obj$survival
        chf <- obj$chf
        cif <- obj$cif
    }
    else {
        pred <- obj$predicted.oob
        surv <- obj$survival.oob
        chf <- obj$chf.oob
        cif <- obj$cif.oob
    }
    ## Survival and competing risk families.
    if (grepl("surv", obj$family)) {
      ## Competing risks:
      if (obj$family == "surv-CR") {
        n <- dim(pred)[1]
        if (missing(subset)) subset <- 1:n
        ## Default is first event type.
        if (missing(target)) target <- 1
        type <- match.arg(type, c("years.lost", "cif", "chf"))
        ## Get the index of the closest point of the grow time interest vector.
        time.idx <-  max(which(obj$time.interest <= time))
        return(switch(type,
                      "years.lost" = pred[subset, target],
                      "cif"        = cif[subset, time.idx, target],
                      "chf"        = chf[subset, time.idx, target]
                      ))
      }
      ## Right-censored:
      else {
        n <- length(pred)
        if (missing(subset)) subset <- 1:n
        type <- match.arg(type, c("rel.freq", "mort", "chf", "surv"))
        ## Get the index of the closest point of the grow time interest vector.
        time.idx <-  max(which(obj$time.interest <= time))
        return(switch(type,
                      "rel.freq" = pred[subset]/max(n, na.omit(pred)),
                      "mort"     = pred[subset],
                      "chf"      = 100 * chf[subset, time.idx],
                      "surv"     = 100 * surv[subset, time.idx]
                      ))
      }
    }
    ## The object has been coerced.  It will be univariate.
      else {
        if (obj$family == "class") {
          n <- dim(pred)[1]
          if (missing(subset)) subset <- 1:n
          type <- match.arg(type, c("prob", "bayes"))
          ## The default is first class label.
          if (missing(target)) target <- 1
          prob <- pred[subset,, drop = FALSE]
          return(switch(type,
                        "prob" = prob[, target],
                        "bayes" =  bayes.rule(prob)))
        }
          else {
            n <- length(pred)
            if (missing(subset)) subset <- 1:n
            return(pred[subset])
          }
      }
  }
## extraction function when partial=TRUE
## Note that currently only one (1) time point is requested via the
## plot.variable() R-wrapper.  Only one (1) is allowed.
## Additional time points can be specified
## using the partial.rfsrc() wrapper, and then extracted without any
## additional calls to the native code.  The extract call
## is limited to only one (1) time point on each call.
## The same protocol holds for m.target.
extract.partial.pred <- function(obj, type, subset, m.target, target) {
  ## Survival families
  if (grepl("surv", obj$family)) {
    n <- dim(obj$survOutput)[1]
    if (missing(subset)) {
      subset <- 1:n
    }
    time.idx <-  1
    ## Competing risks:
    if (obj$family == "surv-CR") {
      ## Default is first event type.
      if (missing(target)) target <- 1
      type <- match.arg(type, c("years.lost", "cif", "chf"))
      return(switch(type,
                    "years.lost" = obj$survOutput[subset, target, ],
                    "cif" = obj$survOutput[subset, time.idx, target, ],
                    "chf" = obj$survOutput[subset, time.idx, target, ]
                    ))
    }
    ## Right censored:
      else {
        if (type == "rel.freq") {
          sz <- apply(obj$survOutput[subset, ], 2, function(x) {length(na.omit(x))})
          rs <- t(apply(obj$survOutput[subset, ], 1, function(x) {x / sz}))
        }
        return(switch(type,
                      "rel.freq" = rs,
                      "mort"     = obj$survOutput[subset, ],
                      "chf"      =  obj$survOutput[subset, time.idx, ],
                      "surv"     =  obj$survOutput[subset, time.idx, ]
                      ))
      }
  }
    else {
      ## Univariate and multivariate families:
      ## The incoming partial object is always presented in a mulitvariate list, regardless of whether
      ## the grow object was univariate or multivariate.  There is no coersion.  In a multivariate case,
      ## m.target is coherent.  But in a univariate case, it will be NULL.  As a result, we
      ## accomodate for the NULL case by grabbing the first (and actually only) valid output below, in
      ## order to parse the multivariate list correctly.
      ## We are either classification or regression but we don't know which yet.
      regrClassTarget <- NULL
      ## It doesn't matter whether we process regression or classification first.
      ## But one and only one family will be processed.
      if (is.null(regrClassTarget)) {
        if (!is.null(obj$classOutput)) {
          if (is.null(m.target)) {
            ## Grab the first outcome!  It will exist, as this is the univariate case.
            m.target <- names(obj$classOutput)[1]
          }
          regrClassTarget <- which (names(obj$classOutput) == m.target)
          ## Classification.
          if (length(regrClassTarget) > 0) {
            if (length(regrClassTarget) > 1) {
              ## Safety check in case m.target was incorrect coming into this routine.
              stop("Invalid number of target outcomes specified in partial plot extraction.")
            }
            type <- match.arg(type, c("prob", "bayes"))
            ## Default is first class label.
            if (missing(target)) target <- 1
            n <- dim(obj$classOutput[[regrClassTarget]])[1]
            if (missing(subset)) subset <- 1:n
            return(switch(type,
                          "prob" = obj$classOutput[[regrClassTarget]][subset, 1 + target, ],
                          "bayes" =  obj$classOutput[[regrClassTarget]][subset, 1, ]))
          }
        }
      }
      ## It doesn't matter whether we process regression or classification first.
      ## But one and only one family will be processed.
      if (is.null(regrClassTarget)) {
        if (!is.null(obj$regrOutput)) {
          if (is.null(m.target)) {
            ## Grab the first outcome!  It will exist, as this is the univariate case.
            m.target <- names(obj$regrOutput)[1]
          }
          regrClassTarget <- which (names(obj$regrOutput) == m.target)
          ## Regression.
          if (length(regrClassTarget) > 0) {
            if (length(regrClassTarget) > 1) {            
              ## Safety check in case m.target was incorrect coming into this routine.
              stop("Invalid number of target outcomes specified in partial plot extraction.")
            }
            n <- dim(obj$classOutput[[regrClassTarget]])[1]
            if (missing(subset)) subset <- 1:n
            return(obj$regrOutput[[regrClassTarget]][subset, ])
          }
        }
      }
      ## If after parsing both classification and regression we have not found the target, we error.
      if (is.null(regrClassTarget)) {
        stop("Invalid target specified in partial plot extraction:  ", m.target)
      }
    }
}
plot.variable <- plot.variable.rfsrc
