## --------------------------------------------------------------------------------------------
##
## extraction function when partial=FALSE
##
## --------------------------------------------------------------------------------------------
extract.pred <- function(obj, type, subset, time, m.target, target, oob = oob) {
  ## coerce the (potentially) multivariate object if necessary.
  obj <- coerce.multivariate(obj, m.target)
  ## extract the sample size, assign the default subset
  if (missing(subset)) {
    subset <- 1:obj$n
  }
  n <- length(subset)
  ## decide if OOB or in-bag values are requested
  if (oob == FALSE) {
    if (is.null(obj$predicted)) {
      stop("inbag requested but inbag values are not available in the object")
    }
    pred <- obj$predicted
    surv <- obj$survival
    chf <- obj$chf
    cif <- obj$cif
  }
  else {
    if (is.null(obj$predicted.oob)) {
      stop("oob requested but oob values are not available in the object")
    }
    pred <- obj$predicted.oob
    surv <- obj$survival.oob
    chf <- obj$chf.oob
    cif <- obj$cif.oob
  }
  ## Survival and competing risk families.
  if (grepl("surv", obj$family)) {
    ## Competing risks:
    if (obj$family == "surv-CR") {
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
      type <- match.arg(type, c("bayes", "rfq", "prob"))
      ## The default is first class label.
      if (missing(target)) target <- 1
      ## extract the subsetted values
      prob <- pred[subset,, drop = FALSE]
      y <- obj$yvar[subset]
      ## rfq only permitted in two class problems
      if (type == "rfq" && nlevels(y) != 2) {
        type <- "bayes"
      }
      ## return the requested performance
      return(switch(type,
                    "bayes" =  get.bayes.rule(prob),
                    "rfq"   =  get.gmean.rule(y, prob),
                    "prob"  =  prob[, target]))
    }
    else {
      return(pred[subset])
    }
  }
}
## --------------------------------------------------------------------------------------------
##
## extraction function when partial=TRUE
## Note that currently only one (1) time point is requested via the
## plot.variable() R-wrapper.  Only one (1) is allowed.
## Additional time points can be specified
## using the partial.rfsrc() wrapper, and then extracted without any
## additional calls to the native code.  The extract call
## is limited to only one (1) time point on each call.
## The same protocol holds for m.target.
##
## --------------------------------------------------------------------------------------------
extract.partial.pred <- function(obj, type, subset, m.target, target) {
  ## Survival families
  if (grepl("surv", obj$family)) {
    ## extract the sample size, assign the default subset
    if (missing(subset)) {
      subset <- 1:obj$n
    }
    n <- length(subset)
    ## array dimension
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
          regrClassTarget <- which(names(obj$classOutput) == m.target)
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
