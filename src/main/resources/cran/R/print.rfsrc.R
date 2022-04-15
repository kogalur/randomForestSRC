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
  ## classification: outcome frequency/confusion matrix/brier/auc/gmean
  rfq.flag <- FALSE
  if (x$family == "class") {
    if (!is.null(x$yvar)) {
      event.freq <- paste(tapply(x$yvar, x$yvar, length), collapse = ", ")
    }
    else {
      event.freq <- NA
    }
    if (!is.null(x$err.rate) && !is.null(x$yvar)) {
      conf.matx <- get.confusion(x$yvar,
                         if(!is.null(x$class.oob) && !all(is.na(x$class.oob))) x$class.oob else x$class)
      miss.err.rate <- 1 - sum(diag(conf.matx[, -ncol(conf.matx), drop = FALSE])) / sum(conf.matx[, -ncol(conf.matx), drop = FALSE])
      names(dimnames(conf.matx)) <- c("  observed", "predicted")
      ## J > 2 class problems
      if (length(levels(x$yvar)) > 2) {
        brier.err <- get.brier.error(x$yvar,
             if(!is.null(x$predicted.oob) && !all(is.na(x$predicted.oob))) x$predicted.oob else x$predicted,
             normalized = FALSE)
        brier.norm.err <- get.brier.error(x$yvar,
             if(!is.null(x$predicted.oob) && !all(is.na(x$predicted.oob))) x$predicted.oob else x$predicted)
        auc.err <- get.auc(x$yvar,
             if(!is.null(x$predicted.oob) && !all(is.na(x$predicted.oob))) x$predicted.oob else x$predicted)
        iratio <- pr.auc.err <- gmean.err <- NULL
      }
      ## performance for two-class setting (returns NULL otherwise)
      else {
        ## use rfq thresholding when gmean performance has been requested
        if (x$forest$rfq) {
          threshold <- get.rfq.threshold(x$forest$yvar)
        }
        else {
          threshold <- 0.5
        }
        perO <- get.imbalanced.performance(x$yvar,
                  if(!is.null(x$predicted.oob) && !all(is.na(x$predicted.oob))) x$predicted.oob else x$predicted,
                  threshold = threshold, confusion = TRUE)
        iratio <- perO$iratio
        brier.err <- perO$brier
        brier.norm.err <- perO$brier.norm
        auc.err <- perO$auc
        pr.auc.err <- perO$pr.auc
        gmean.err <- perO$gmean
        #if (x$forest$rfq) {
        #  conf.matx <- perO$confusion
        #  names(dimnames(conf.matx)) <- c("  observed", "predicted")
        #}
      }
    }
    else {
      conf.matx <- miss.err.rate <- iratio <- brier.err <- brier.norm.err <- auc.err <- pr.auc.err <- gmean.err <- NULL
    }
  }
  ## survival, CR: CRPS and event frequencies
  if (grepl("surv", x$family)) {
    if (!is.null(x$err.rate) &x$family == "surv") {
      crps.err  <- get.brier.survival(x)$crps.std
    }
    else {
      crps.err <- NULL
    }
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
  ## error rates 
  if (!is.null(x$err.rate)) {
    err.rate <- cbind(x$err.rate)
    allcol.na <- apply(err.rate, 2, function(x) {all(is.na(x))})
    err.rate <- na.omit(err.rate[, !allcol.na, drop = FALSE])
    attributes(err.rate)$na.action <- NULL
    if (grepl("surv", x$family)) {
      err.rate <- digits.pretty(err.rate[nrow(err.rate), ], 8)
    }
    else if (x$family == "class") {
      ## rfq related adjustments
      if ((grow.mode && x$forest$perf.type == "gmean") || (!grow.mode && x$perf.type == "gmean")) {
        err.rate <- digits.pretty(err.rate[nrow(err.rate), 1], 8)
      }
      else {
        err.rate <- digits.pretty(err.rate[nrow(err.rate), ], 8)
      }
    }
    else if (x$family == "regr") {
      if (!is.null(x$yvar)) {
        r.squared <- 1 - err.rate[nrow(err.rate), ] / var(x$yvar, na.rm = TRUE)
      }
      else {
        r.squared <- NULL
      }
      err.rate <- digits.pretty(err.rate[nrow(err.rate), ], 8)
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
        cat("                    Number of deaths: ", x$ndead,    "\n", sep="")
      }
    }
    if (x$family == "class") {
      cat("           Frequency of class labels: ", event.freq, "\n", sep="")
    }
    if (!is.null(x$imputed.indv)) {
      cat("                    Was data imputed: ", "yes",      "\n", sep="")
    }
    cat("                     Number of trees: ", x$ntree,                              "\n", sep="")
    cat("           Forest terminal node size: ", x$nodesize,                           "\n", sep="")
    cat("       Average no. of terminal nodes: ", digits.pretty(mean(x$leaf.count), 4), "\n", sep="")
    cat("No. of variables tried at each split: ", x$mtry,                               "\n", sep="")
    cat("              Total no. of variables: ", length(x$xvar.names),                 "\n", sep="")
    if (!x$univariate) { 
      cat("              Total no. of responses: ", yvar.dim,                           "\n", sep="")
      cat("         User has requested response: ", outcome.target,                     "\n", sep="")
    }
    cat("       Resampling used to grow trees: ", sampUsed,                             "\n", sep="")
    cat("    Resample size used to grow trees: ", round(x$forest$sampsize(x$n)),        "\n", sep="")
    cat("                            Analysis: ", familyPretty,                         "\n", sep="")
    cat("                              Family: ", familyOrg,                            "\n", sep="")
    if (x$nsplit > 0 & x$splitrule != "random") {
      cat("                      Splitting rule: ", paste(x$splitrule,"*random*"),      "\n", sep="")
      cat("       Number of random split points: ", x$nsplit,                           "\n", sep="")
    }
    else {
      cat("                      Splitting rule: ", x$splitrule,                        "\n", sep="")
    } 
    if (!is.null(err.rate)) {
      if (x$family == "regr" && !is.null(r.squared)) {
        cat("                     (OOB) R squared: ", digits.pretty(r.squared, 8),      "\n", sep="")
      }
      if (x$family == "class" && !is.null(iratio)) {
        cat("                    Imbalanced ratio: ", digits.pretty(iratio, 4),         "\n", sep="")
      }
      if (x$family == "class" && !is.null(brier.err)) {
        cat("                   (OOB) Brier score: ", digits.pretty(brier.err, 8),      "\n", sep="")
      }
      if (x$family == "class" && !is.null(brier.norm.err)) {
        cat("        (OOB) Normalized Brier score: ", digits.pretty(brier.norm.err, 8), "\n", sep="")
      }
      if (x$family == "class" && !is.null(auc.err)) {
        cat("                           (OOB) AUC: ", digits.pretty(auc.err, 8),        "\n", sep="")
      }
      if (x$family == "class" && !is.null(pr.auc.err)) {
        cat("                        (OOB) PR-AUC: ", digits.pretty(pr.auc.err, 8),     "\n", sep="")
      }
      if (x$family == "class" && !is.null(gmean.err)) {
        cat("                        (OOB) G-mean: ", digits.pretty(gmean.err, 8),      "\n", sep="")
      }
      if (x$family == "surv" && !is.null(crps.err)) {
        cat("                          (OOB) CRPS: ", digits.pretty(crps.err, 8),      "\n", sep="")
      }
           cat("   (OOB) Requested performance error: ", err.rate,                      "\n\n", sep="")
    }
    if (x$family == "class" && !is.null(conf.matx)) {
      if (!is.null(x$predicted.oob) && any(is.na(x$predicted.oob))) {
        cat("Confusion matrix (cases with missing OOB predicted values have been removed):\n\n")
      }
      else {
        cat("Confusion matrix:\n\n")
      }
      print(conf.matx)
        cat("\n      (OOB) Misclassification rate: ", miss.err.rate,   "\n", sep="")
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
        cat("       Number of events in test data: ", event.freq,         "\n", sep="")
      }
      else {
        cat("       Number of deaths in test data: ", unlist(event.freq), "\n", sep="")
      }
    }
    if (!is.null(x$imputed.data)) {
      cat("               Was test data imputed: ", "yes", "\n", sep="")
    }
    cat("                Number of grow trees: ", x$ntree,                                "\n",sep="")
    cat("  Average no. of grow terminal nodes: ", digits.pretty(mean(x$leaf.count), 4),   "\n", sep="")
    cat("         Total no. of grow variables: ", length(x$xvar.names),                   "\n", sep="")  
    if (!x$univariate) { 
      cat("         Total no. of grow responses: ", yvar.dim,                             "\n", sep="")
      cat("         User has requested response: ", outcome.target,                       "\n", sep="")
    }
    cat("       Resampling used to grow trees: ", sampUsed,                               "\n",sep="")
    cat("    Resample size used to grow trees: ", round(x$forest$sampsize(x$forest$n)),   "\n",sep="")
    cat("                            Analysis: ", familyPretty,                           "\n", sep="")
    cat("                              Family: ", familyOrg,                              "\n", sep="")
    if (!is.null(err.rate)) {
      if (x$family == "regr" && !is.null(r.squared)) {
        cat("                           R squared: ", digits.pretty(r.squared, 8),        "\n", sep="")
      }
      if (x$family == "class" && !is.null(iratio)) {
        cat("                    Imbalanced ratio: ", digits.pretty(iratio, 4),           "\n", sep="")
      }
      if (x$family == "class" && !is.null(brier.err)) {
        cat("                         Brier score: ", digits.pretty(brier.err, 8),        "\n", sep="")
      }
      if (x$family == "class" && !is.null(brier.norm.err)) {
        cat("              Normalized Brier score: ", digits.pretty(brier.norm.err, 8),   "\n", sep="")
      }
      if (x$family == "class" && !is.null(auc.err)) {
        cat("                                 AUC: ", digits.pretty(auc.err, 8),          "\n", sep="")
      }
      if (x$family == "class" && !is.null(pr.auc.err)) {
        cat("                              PR-AUC: ", digits.pretty(pr.auc.err, 8),       "\n", sep="")
      }
      if (x$family == "class" && !is.null(gmean.err)) {
        cat("                              G-mean: ", digits.pretty(gmean.err, 8),        "\n", sep="")
      }
      if (x$family == "surv" && !is.null(crps.err)) {
        cat("                                CRPS: ", digits.pretty(crps.err, 8),      "\n", sep="")
      }
        cat("         Requested performance error: ", err.rate,                         "\n\n", sep="")
    }
    if (x$family == "class" && !is.null(conf.matx)) {
      if (!is.null(x$predicted.oob) && any(is.na(x$predicted.oob))) {
        cat("Confusion matrix (cases with missing OOB predicted values have been removed):\n\n")
      }
      else {
        cat("Confusion matrix:\n\n")
      }
      print(conf.matx)
        cat("\n           Misclassification error: ", miss.err.rate,  "\n", sep="")
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
