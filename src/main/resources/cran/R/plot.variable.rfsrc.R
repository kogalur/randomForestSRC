plot.variable.rfsrc <- function(
  x,
  xvar.names,
  which.class,
  outcome.target=NULL,
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
  if (sum(inherits(x, c("rfsrc", "synthetic"), TRUE) == c(1, 2)) == 2) {
    x <- x$rfSyn
  }
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
  if (partial && is.null(object$forest)) {
    stop("forest is empty:  re-run rfsrc (grow) call with forest=TRUE")
  }
  xvar <- object$xvar
  if (!is.null(object$imputed.indv)) {
    xvar[object$imputed.indv, ] <- object$imputed.data[, object$xvar.names]
  }
  n <- nrow(xvar)
  if (!inherits(object, "plot.variable")) {
    if (missing(subset)) {
      subset <- 1:n
    }
      else {
        if (is.logical(subset)) subset <- which(subset)
        subset <- unique(subset[subset >= 1 & subset <= n])
        if (length(subset) == 0) {
          stop("'subset' not set properly")
        }
      }
    xvar <- xvar[subset,, drop = FALSE]
    n <- nrow(xvar)
    family <- object$family
    outcome.target <- get.univariate.target(object, outcome.target)
    if (grepl("surv", family)) {
      event.info <- get.event.info(object, subset)
      cens <- event.info$cens
      event.type <- event.info$event.type
      if (missing(time)) {
        time <- median(event.info$time.interest, na.rm = TRUE)
      }
      if (length(time) > 1) {
        stop("time must be a single value:  ", time)
      }
      if (is.null(surv.type)) {
        stop("surv.type is specified incorrectly:  ", surv.type)
      }
      if (family == "surv-CR") {
        if (missing(which.class)) {
          which.class <- 1
        }
          else {
            if (which.class < 1 || which.class > max(event.type, na.rm = TRUE)) {
              stop("'which.class' is specified incorrectly")
            }
          }
        VIMP <- object$importance[, which.class]
        pred.type <- setdiff(surv.type, c("rel.freq", "mort", "chf", "surv"))[1]
        pred.type <- match.arg(pred.type, c("years.lost", "cif", "chf"))
        ylabel <- switch(pred.type,
                         "years.lost" = paste("Years lost for event ", which.class),
                         "cif"        = paste("CIF for event ", which.class, " (time=", time, ")", sep = ""),
                         "chf"        = paste("CHF for event ", which.class, " (time=", time, ")", sep = ""))
      }
        else {
          which.class <- 1
          VIMP <- object$importance
          pred.type <- setdiff(surv.type, c("years.lost", "cif", "chf"))[1]
          pred.type <- match.arg(pred.type, c("rel.freq", "mort", "chf", "surv"))
          ylabel <- switch(pred.type,
                           "rel.freq"   = "standardized mortality",
                           "mort"       = "mortality",
                           "chf"        = paste("CHF (time=", time, ")", sep = ""),
                           "surv"       = paste("predicted survival (time=", time, ")", sep = ""))
        }
    }
      else {
        event.info <- time <- NULL
        if(is.factor(coerce.multivariate(object, outcome.target)$yvar)) {
          object.yvar.levels <- levels(coerce.multivariate(object, outcome.target)$yvar)
          pred.type <- match.arg(class.type, c("prob", "bayes"))
          if (missing(which.class)) {
            which.class <- object.yvar.levels[1]
          }
          if (is.character(which.class)) {
            which.class <- match(match.arg(which.class, object.yvar.levels), object.yvar.levels)
          }
            else {
              if ((which.class > length(object.yvar.levels)) | (which.class < 1)) {
                stop("which.class is specified incorrectly:", which.class)
              }
            }
          if (pred.type == "prob") {
            VIMP <- coerce.multivariate(object, outcome.target)$importance[, 1 + which.class]
            ylabel <- paste("probability", object.yvar.levels[which.class])
            remove(object.yvar.levels)
          }
            else {
              VIMP <- coerce.multivariate(object, outcome.target)$importance[, 1]
              ylabel <- paste("bayes")
              remove(object.yvar.levels)
            }
        }
          else {
            pred.type <- "y"
            which.class <- NULL
            VIMP <- coerce.multivariate(object, outcome.target)$importance
            ylabel <- expression(hat(y))
          }
      }
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
      yhat <- extract.pred(predict.rfsrc(object,
                                         outcome.target = outcome.target,
                                         importance = "none"),
                           pred.type,
                           subset,
                           time,
                           outcome.target,
                           which.class,
                           oob = oob)
    }
      else {
        if (npts < 1) npts <- 1 else npts <- round(npts)
        prtl <- lapply(1:nvar, function(k) {
          x <- na.omit(object$xvar[, object$xvar.names == xvar.names[k]])
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
                                                          outcome.target = outcome.target,
                                                          partial.type = pred.type,
                                                          partial.xvar = xvar.names[k],
                                                          partial.values = x.uniq,
                                                          partial.time = time,
                                                          oob = oob),
                                            pred.type,
                                            1:n,
                                            outcome.target,
                                            which.class)
          if (!is.null(dim(pred.temp))) {
            mean.temp <- apply(pred.temp, 2, mean, na.rm = TRUE)
          }
            else {
              mean.temp <- mean(pred.temp, na.rm = TRUE)
            }
          if (!factor.x) {
            yhat <- mean.temp
            if (coerce.multivariate(object, outcome.target)$family == "class") {
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
    plot.variable.obj <- list(family = family,
                              partial = partial,
                              event.info = event.info,
                              which.class = which.class,
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
    class(plot.variable.obj) <- c("rfsrc", "plot.variable", family)
  }
    else {
      plot.variable.obj <- object
      remove(object)
      family <- plot.variable.obj$family
      partial <- plot.variable.obj$partial
      event.info <- plot.variable.obj$event.info
      which.class <- plot.variable.obj$which.class
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
  if (show.plots) {
    old.par <- par(no.readonly = TRUE)
  }
  if (!partial && show.plots) {
    par(mfrow = c(min(plots.per.page, ceiling(nvar / plots.per.page)), plots.per.page))
    if (n > 500) cex.pt <- 0.5 else cex.pt <- 0.75
    for (k in 1:nvar) {
      x <- xvar[, which(colnames(xvar) == xvar.names[k])]
      x.uniq <- unique(x)
      n.x <- length(x.uniq)
      if (!is.factor(x) & n.x > granule) {
        plot(x,
             yhat,
             xlab = xvar.names[k],
             ylab = ylabel,
             type = "n", ...) 
        if (grepl("surv", family)) {
          points(x[cens == which.class], yhat[cens == which.class], pch = 16, col = 4, cex = cex.pt)
          points(x[cens == 0], yhat[cens == 0], pch = 16, cex = cex.pt)
        }
        lines(lowess(x[!is.na(x)], yhat[!is.na(x)]), col = 2, lwd=3)
      }
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
  if (partial && show.plots) {
    plots.per.page <- max(round(min(plots.per.page,nvar)), 1)
    granule <- max(round(granule),1)
    par(mfrow = c(min(plots.per.page, ceiling(nvar/plots.per.page)), plots.per.page))
    for (k in 1:nvar) {
      x <- prtl[[k]]$x
      if (is.factor(x)) x <- factor(x, exclude = NULL)          
      x.uniq <- prtl[[k]]$x.uniq
      n.x <- prtl[[k]]$n.x
      if (n.x > 25) cex.pt <- 0.5 else cex.pt <- 0.75
      yhat <- prtl[[k]]$yhat
      yhat.se <- prtl[[k]]$yhat.se
      factor.x <- !(!is.factor(x) & (n.x > granule))
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
  if (show.plots) {
    par(old.par)
  }
  invisible(plot.variable.obj)
}
extract.pred <- function(obj, type, subset, time, outcome.target, which.class, oob = oob) {
    obj <- coerce.multivariate(obj, outcome.target)
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
    if (grepl("surv", obj$family)) {
      if (obj$family == "surv-CR") {
        n <- dim(pred)[1]
        if (missing(subset)) subset <- 1:n
        if (missing(which.class)) which.class <- 1
        type <- match.arg(type, c("years.lost", "cif", "chf"))
        time.idx <-  max(which(obj$time.interest <= time))
        return(switch(type,
                      "years.lost" = pred[subset, which.class],
                      "cif"        = cif[subset, time.idx, which.class],
                      "chf"        = chf[subset, time.idx, which.class]
                      ))
      }
      else {
        n <- length(pred)
        if (missing(subset)) subset <- 1:n
        type <- match.arg(type, c("rel.freq", "mort", "chf", "surv"))
        time.idx <-  max(which(obj$time.interest <= time))
        return(switch(type,
                      "rel.freq" = pred[subset]/max(n, na.omit(pred)),
                      "mort"     = pred[subset],
                      "chf"      = 100 * chf[subset, time.idx],
                      "surv"     = 100 * surv[subset, time.idx]
                      ))
      }
    }
      else {
        if (obj$family == "class") {
          n <- dim(pred)[1]
          if (missing(subset)) subset <- 1:n
          type <- match.arg(type, c("prob", "bayes"))
          if (missing(which.class)) which.class <- 1
          prob <- pred[subset,, drop = FALSE]
          return(switch(type,
                        "prob" = prob[, which.class],
                        "bayes" =  bayes.rule(prob)))
        }
          else {
            n <- length(pred)
            if (missing(subset)) subset <- 1:n
            return(pred[subset])
          }
      }
  }
extract.partial.pred <- function(obj, type, subset, outcome.target, which.class) {
  if (grepl("surv", obj$family)) {
    n <- dim(obj$survOutput)[1]
    if (missing(subset)) subset <- 1:n
    time.idx <-  1
    if (obj$family == "surv-CR") {
      if (missing(which.class)) which.class <- 1
      type <- match.arg(type, c("years.lost", "cif", "chf"))
      return(switch(type,
                    "years.lost" = obj$survOutput[subset, which.class, ],
                    "cif" = obj$survOutput[subset, time.idx, which.class, ],
                    "chf" = obj$survOutput[subset, time.idx, which.class, ]
                    ))
    }
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
      target <- NULL
      if (is.null(target)) {
        if (!is.null(obj$classOutput)) {
          if (is.null(outcome.target)) {
            outcome.target <- names(obj$classOutput)[1]
          }
          target <- which (names(obj$classOutput) == outcome.target)
          if (length(target) > 0) {
            if (length(target) > 1) {
              stop("Invalid number of target outcomes specified in partial plot extraction.")
            }
            type <- match.arg(type, c("prob", "bayes"))
            if (missing(which.class)) which.class <- 1
            n <- dim(obj$classOutput[[target]])[1]
            if (missing(subset)) subset <- 1:n
            return(switch(type,
                          "prob" = obj$classOutput[[target]][subset, 1 + which.class, ],
                          "bayes" =  obj$classOutput[[target]][subset, 1, ]))
          }
        }
      }
      if (is.null(target)) {
        if (!is.null(obj$regrOutput)) {
          if (is.null(outcome.target)) {
            outcome.target <- names(obj$regrOutput)[1]
          }
          target <- which (names(obj$regrOutput) == outcome.target)
          if (length(target) > 0) {
            if (length(target) > 1) {            
              stop("Invalid number of target outcomes specified in partial plot extraction.")
            }
            n <- dim(obj$classOutput[[target]])[1]
            if (missing(subset)) subset <- 1:n
            return(obj$regrOutput[[target]][subset, ])
          }
        }
      }
      if (is.null(target)) {
        stop("Invalid target specified in partial plot extraction:  ", outcome.target)
      }
    }
}
plot.variable <- plot.variable.rfsrc
