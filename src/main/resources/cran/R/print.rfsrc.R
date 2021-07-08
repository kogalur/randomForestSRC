print.rfsrc <- function(x, outcome.target = NULL, ...) {
  ## default printing
  if (sum(inherits(x, c("rfsrc", "forest"), TRUE) == c(1, 2)) == 2) {
    print.default(x)
    return()
  }
  if(sum(inherits(x, c("rfsrc", "plot.variable"), TRUE) == c(1, 2)) == 2) {
    print.default(x)
    return()
  }
  if (sum(inherits(x, c("rfsrc", "partial"), TRUE) == c(1, 2)) == 2) {
    print.default(x)
    return()
  }
  if (sum(inherits(x, c("rfsrc", "sidClustering"), TRUE) == c(1, 2)) == 2) {
    print.default(x)
    return()
  }
  ## processing for synthetic forests        
  sf.flag <- FALSE
  if (sum(inherits(x, c("rfsrc", "synthetic"), TRUE) == c(1, 2)) == 2) {
    if (sum(inherits(x, c("rfsrc", "synthetic", "oob"), TRUE) == c(1, 2, 3)) != 3) {
      sf.flag <- TRUE
      sf.message <- "OOB was not used for synthetic forests, error rates/VIMP will be unreliable"
    }
    x <- x$rfSyn
  }
  ## processing for subsampled object
  if (sum(inherits(x, c("rfsrc", "subsample"), TRUE) == c(1, 4)) == 2) {
    print.subsample(x, ...)
    return()
  }
  ## processing for subsampled-bootstrap object?
  if (sum(inherits(x, c("rfsrc", "bootsample"), TRUE) == c(1, 4)) == 2) {
    print.bootsample(x, ...)
    return()
  }
  ## check that the object is interpretable
  if (sum(inherits(x, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rfsrc", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("This function only works for objects of class `(rfsrc, grow)' or '(rfsrc, predict)'.")
  }
  ## are we in grow mode?
  if (sum(inherits(x, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
    grow.mode <- TRUE
  }
  else {
    grow.mode <- FALSE
  }
  ## specify the type of sampling
  if (x$forest$bootstrap == "by.root") {
    sampUsed <- x$forest$samptype
  }
  else {
    sampUsed <- x$forest$bootstrap
  }
  ## x will be processed if it's multivariate - therefore save some values from it
  familyPretty <- family.pretty(x)
  familyOrg <- x$family
  yvar.dim <- ncol(x$yvar)
  ## coerce the (potentially) multivariate object if necessary.
  outcome.target <- get.univariate.target(x, outcome.target)
  x <- coerce.multivariate(x, outcome.target)
  ## survival: event frequencies
  if (grepl("surv", x$family)) {
    event <- x$event.info$event
    n.event <- 1
    if (!is.null(event)) {
      n.event <- length(unique(event))
      if (length(event) > 0) {
        event.freq <- paste(tapply(event, event, length), collapse = ", ")
      }
        else {
          event.freq <- 0
        }
    }
  }
  ## classification: outcome frequency/confusion matrix/brier/auc/gmean
  if (x$family == "class") {
    if (!is.null(x$yvar)) {
      event.freq <- paste(tapply(x$yvar, x$yvar, length), collapse = ", ")
    }
    else {
      event.freq <- NA
    }
    if (!is.null(x$err.rate)) {
      conf.matx <- get.confusion(x$yvar,
             if(!is.null(x$class.oob) && !all(is.na(x$class.oob))) x$class.oob else x$class)
      names(dimnames(conf.matx)) <- c("  observed", "predicted")
      brier.err <- get.brier.error(x$yvar,
             if(!is.null(x$predicted.oob) && !all(is.na(x$predicted.oob))) x$predicted.oob else x$predicted)
      auc.error <- get.auc(x$yvar,
             if(!is.null(x$predicted.oob) && !all(is.na(x$predicted.oob))) x$predicted.oob else x$predicted)
      ## special processing needed to handle class imbalanced rfq classifier
      if (grow.mode) {
        rfqO <- list(rfq = x$forest$rfq, perf.type = x$forest$perf.type)
      }
      else {
        rfqO <- list(rfq = x$forest$rfq, perf.type = x$perf.type)
      }
      if (!is.null(rfqO$perf.type) && (rfqO$perf.type == "g.mean")) {
        gmean.error <- 1 - x$err.rate[nrow(x$err.rate), 1]
        pi.hat <- table(x$yvar) / length(x$yvar)
        iratio <- max(pi.hat, na.rm  = TRUE) / min(pi.hat, na.rm  = TRUE)
      }
      else {
        gmean.error <- NULL
      }
    }
    else {
      conf.matx <- brier.err <- auc.error <- gmean.error <- NULL
    }
  }
  ## error rates 
  if (!is.null(x$err.rate)) {
    err.rate <- cbind(x$err.rate)    
    if (grepl("surv", x$family)) {
      err.rate <- paste(round(100 * err.rate[nrow(err.rate), ], 8), "%", collapse=", ", sep = "")
    }
    else if (x$family == "class") {
      brier.err <- 100 * brier.err
      auc.error <- 100 * auc.error
      overall.err.rate <- paste(round(100 * err.rate[nrow(err.rate), 1], 6), "%", sep = "")
      ## rfq related adjustments
      if (!is.null(gmean.error)) {
        err.rate <- err.rate[nrow(err.rate), 1]
      }
      else {
        err.rate <- paste(round(err.rate[nrow(err.rate), ], 8), collapse=", ", sep = "")
      }
    }
    else if (x$family == "regr") {
      r.squared <- 1 - err.rate[nrow(err.rate), ] / var(x$yvar, na.rm = TRUE)
      err.rate <- err.rate[nrow(err.rate), ]
    }
    else {
      err.rate <- NULL
    }
  }
  else {
    err.rate <- NULL
  }
  ## ensure backward compatibility for nsplit
  if (is.null(x$nsplit)) {
    x$nsplit <- 0
  }
  #################################################################################
  ##
  ##
  ## GROW MODE
  ##
  ##
  #################################################################################
  if (grow.mode) {
    cat("                         Sample size: ", x$n, "\n", sep="")
    if (grepl("surv", x$family)) {
      if (n.event > 1) {
        cat("                    Number of events: ", event.freq, "\n", sep="")
      }
        else {
          cat("                    Number of deaths: ", x$ndead, "\n", sep="")
        }
    }
    if (x$family == "class") {
      cat("           Frequency of class labels: ", event.freq, "\n", sep="")
    }
    if (!is.null(x$imputed.indv)) {
      cat("                    Was data imputed: ", "yes", "\n", sep="")
    }
    cat("                     Number of trees: ", x$ntree,                "\n",sep="")
    cat("           Forest terminal node size: ", x$nodesize,             "\n", sep="")
    cat("       Average no. of terminal nodes: ", mean(x$leaf.count),     "\n", sep="")
    cat("No. of variables tried at each split: ", x$mtry,                 "\n", sep="")
    cat("              Total no. of variables: ", length(x$xvar.names),   "\n", sep="")
    if (!x$univariate) { 
      cat("              Total no. of responses: ", yvar.dim, "\n", sep="")
      cat("         User has requested response: ", outcome.target, "\n", sep="")
    }
    cat("       Resampling used to grow trees: ", sampUsed,                     "\n",sep="")
    cat("    Resample size used to grow trees: ", round(x$forest$sampsize(x$n)),"\n",sep="")
    cat("                            Analysis: ", familyPretty,                 "\n", sep="")
    cat("                              Family: ", familyOrg,                    "\n", sep="")
    if (x$nsplit > 0 & x$splitrule != "random") {
      cat("                      Splitting rule: ", paste(x$splitrule,"*random*"),"\n", sep="")
      cat("       Number of random split points: ", x$nsplit                   ,  "\n", sep="")
    }
      else {
        cat("                      Splitting rule: ", x$splitrule, "\n", sep="")
      } 
    if (!is.null(err.rate)) {
      if (x$family == "regr") {
        cat("                     (OOB) R squared: ", r.squared, "\n", sep="")
      }
      if (x$family == "class" && !is.null(gmean.error)) {
        cat("                    Imbalanced ratio: ", iratio, "\n", sep="")
        cat("                        (OOB) G-mean: ", gmean.error, "\n", sep="")
      }
      if (x$family == "class" && !is.null(brier.err)) {
        cat("        (OOB) Normalized Brier score:", brier.err, "\n")
      }
      if (x$family == "class" && !is.null(auc.error)) {
        cat("                           (OOB) AUC:", auc.error, "\n")
      }
      cat("                    (OOB) Error rate: ", err.rate, "\n\n", sep="")
    }
    if (x$family == "class" && !is.null(conf.matx)) {
      if (!is.null(x$predicted.oob) && any(is.na(x$predicted.oob))) {
        cat("Confusion matrix (cases with missing OOB predicted values have been removed):\n\n")
      }
      else {
        cat("Confusion matrix:\n\n")
      }
      print(conf.matx)
      cat("\n\tOverall (OOB) error rate:", overall.err.rate, "\n")
    }
     
  }
  #################################################################################
  ##
  ##
  ## PREDICT MODE
  ##
  ##
  #################################################################################
  else {
    cat("  Sample size of test (predict) data: ", x$n, "\n", sep="")
    if (grepl(x$family, "surv") && !is.null(event)) {
      if (n.event > 1) {
        cat("       Number of events in test data: ", event.freq, "\n", sep="")
      }
      else {
        cat("       Number of deaths in test data: ", unlist(event.freq), "\n", sep="")
      }
    }
    if (!is.null(x$imputed.data)) {
      cat("               Was test data imputed: ", "yes", "\n", sep="")
    }
    cat("                Number of grow trees: ", x$ntree,              "\n",sep="")
    cat("  Average no. of grow terminal nodes: ", mean(x$leaf.count),   "\n", sep="")
    cat("         Total no. of grow variables: ", length(x$xvar.names), "\n", sep="")  
    if (!x$univariate) { 
      cat("         Total no. of grow responses: ", yvar.dim, "\n", sep="")
      cat("         User has requested response: ", outcome.target, "\n", sep="")
    }
    cat("       Resampling used to grow trees: ", sampUsed,                     "\n",sep="")
    cat("    Resample size used to grow trees: ", round(x$forest$sampsize(x$n)),"\n",sep="")
    cat("                            Analysis: ", familyPretty,                 "\n", sep="")
    cat("                              Family: ", familyOrg,                    "\n", sep="")
    if (!is.null(err.rate)) {
      if (x$family == "regr") {
        cat("                           R squared: ", r.squared, "\n", sep="")
      }
      if (x$family == "class" && !is.null(gmean.error)) {
        cat("                    Imbalanced ratio: ", iratio, "\n", sep="")
        cat("                              G-mean: ", gmean.error, "\n", sep="")
      }
      if (x$family == "class" && !is.null(brier.err)) {
        cat("              Normalized Brier score:", brier.err, "\n")
      }
      if (x$family == "class" && !is.null(auc.error)) {
        cat("                                 AUC:", auc.error, "\n")
      }
      cat("                          error rate: ", err.rate, "\n\n", sep="")
    }
    if (x$family == "class" && !is.null(conf.matx)) {
      if (!is.null(x$predicted.oob) && any(is.na(x$predicted.oob))) {
        cat("Confusion matrix (cases with missing OOB predicted values have been removed):\n\n")
      }
      else {
        cat("Confusion matrix:\n\n")
      }
      print(conf.matx)
      cat("\n\tOverall error rate:", overall.err.rate, "\n")
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
