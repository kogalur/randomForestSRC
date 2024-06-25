subsample.rfsrc <- function(obj,                
                            B = 100,
                            block.size = 1,
                            importance,
                            subratio = NULL,
                            stratify = TRUE,
                            performance = FALSE,
                            performance.only = FALSE,
                            joint = FALSE,
                            xvar.names = NULL,
                            bootstrap = FALSE,
                            verbose = TRUE) 
{
  ##--------------------------------------------------------------
  ##
  ## coherence checks
  ##
  ##--------------------------------------------------------------
  ## incoming object must be a grow forest
  if (sum(inherits(obj, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2) {
    stop("This function only works for objects of class `(rfsrc, grow)'")
  }
  ## grow forests must have true forest information
  if (sum(inherits(obj, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.forest.missing(obj)) {
      stop("Forest information for prediction is missing.  Re-run rfsrc (grow call) with forest=TRUE")
    }
  }
  if (inherits(obj, "anonymous")) {
    stop("this function does work with anonymous forests")
  }
  ##--------------------------------------------------------------
  ##
  ## is this a multivariate family?
  ##
  ##--------------------------------------------------------------
  fmly <- obj$family
  mv.fmly <- fmly == "regr+" || fmly == "class+" || fmly == "mix+" 
  ##--------------------------------------------------------------
  ##
  ## set the sample size
  ## set subratio (if not supplied)
  ## confirm subratio is appropriately set
  ##
  ##--------------------------------------------------------------
  n <- obj$n
  if (is.null(subratio)) {
    subratio <- get.subsample.subratio(n)
  }
  if (!bootstrap && (subratio < 0 || subratio > 1)) {
    stop("subratio must be between 0 and 1")
  }
  ##--------------------------------------------------------------
  ##
  ## extract random forest parmaters
  ##
  ##--------------------------------------------------------------
  ## list of parameters to be obtained from the grow object
  grow.prms <- c("call",
                 "ntree",
                 "mtry",
                 "nodesize",
                 "nodedepth",
                 "splitrule",
                 "nsplit",
                 "xvar.wt",
                 "split.wt",
                 "cause.wt",
                 "block.size")
   
  ## list of parameters to be obtained from the forest object
  forest.prms <- c("forest",
                   "bootstrap",
                   "sampsize",
                   "samptype",
                   "case.wt",
                   "perf.type",
                   "vimp.threshold",
                   "rfq")
  ## pull the parameters (grow first, then forest object)
  rf.prms <- c(obj[grow.prms], obj$forest[forest.prms])
  ## rename parameters to their correct random forest names
  rf.prms$call <- formula(rf.prms$call)
  names(rf.prms)[names(rf.prms) == "call"] <- "formula"
  names(rf.prms)[names(rf.prms) == "cause.wt"] <- "cause"
  ## subsamping cannot be applied when case weights are non-standard
  if (!all(diff(rf.prms$case.wt) == 0)) {
    stop("subsampling is not permitted on forests grown under non-standard case weights")
  }
  ## everthing is OK - nuke the case weights
  rf.prms$case.wt <- NULL
  ## need time.interest if this is a survival object
  if (fmly == "surv" || fmly == "surv-CR") {
    rf.prms$ntime <- obj$time.interest
  }
  ## nuke xvar and split weighting if they are standard
  if (all(diff(rf.prms$xvar.wt) == 0)) {
    rf.prms$xvar.wt <- NULL
  }
  if (all(diff(rf.prms$split.wt) == 0)) {
    rf.prms$split.wt <- NULL
  }
  ## nuke nodedepth if it's standard
  if (rf.prms$nodedepth <= 0) {
    rf.prms$nodedepth <- NULL
  }
  ## subsampling is not permitted if bootstrapping is not by root
  if (!bootstrap && rf.prms$bootstrap != "by.root") {
    stop("subsampling is not permitted on grow objects where sampling is not 'by.root'")
  }
  ## for subsampling we need to determine the subsampling tree sample size
  ## this is now handled by forest$sampsize which is a function
  ## for double bootstrapping we do NOT allow custom sampling for grow object
  if (bootstrap && rf.prms$bootstrap != "by.root") {
    stop("double bootstrapping is only permitted for breiman (bootstrap) forests")
  }
  else {
    ##everything is OK - set bootstrap parameter to NULL for later custom bootstrap
    ## leave $samptype and $sampsize in place to be used for making double bootstrap
    rf.prms$bootstrap <- NULL
  }
  ##--------------------------------------------------------------
  ##
  ## performance only?  (user is requesting generalization error only)
  ##
  ##--------------------------------------------------------------
  if (performance.only) {
    performance <- TRUE
    importance <- "none"
    vmp <- NULL
  }
  ##--------------------------------------------------------------
  ##
  ## call the bootstrap subroutine if double bootstrap was requested
  ##
  ##--------------------------------------------------------------
  if (bootstrap) {
    if (missing(importance) && !performance.only) {
      importance <- TRUE
    }
    bootO <- bootsample(obj, rf.prms, B = B, block.size = block.size,
                        joint = joint, xvar.names = xvar.names,
                        importance = importance,
                        performance = performance, performance.only = performance.only,
                        verbose = verbose)
    rO <- list(rf = bootO$rf, vmp = bootO$vmp, vmpB = bootO$vmpB,
               subratio = NULL, performance.only = performance.only)
    class(rO) <- c(class(obj), "bootsample")
    return(rO)
  }
  ##--------------------------------------------------------------
  ##
  ## extract data from the previously grown forest
  ## set the dimension, pull the feature names
  ##
  ##--------------------------------------------------------------
  dta <- data.frame(obj$yvar, obj$xvar)
  colnames(dta)[1:length(obj$yvar.names)] <- obj$yvar.names
  ##--------------------------------------------------------------
  ##
  ## call vimp to extract importance if not available in original object
  ##
  ##--------------------------------------------------------------
  if (!performance.only) {
    vmp <- get.mv.vimp(obj, FALSE, FALSE)
    if (is.null(vmp)) {
      if (verbose) cat("no importance found: calculating it now ...\n")
      if (missing(importance)) {
        importance <- TRUE
      }
      obj <- vimp(obj, perf.type = rf.prms$perf.type, block.size = block.size, importance = importance)
      vmp <- get.mv.vimp(obj, FALSE, FALSE)
      rf.prms$block.size <- block.size
      if (verbose) cat("done\n")
    }
    else {
      importance <- obj$forest$importance
    }
  }
  ##--------------------------------------------------------------
  ##
  ## call joint vimp - reference vimp value for pure noise settings
  ##
  ##--------------------------------------------------------------
  if (joint && !performance.only) {
    vmp.joint <- get.mv.vimp(get.subsample.joint.vimp(obj, rf.prms, xvar.names), FALSE, FALSE)
    vmp <- vimp.subsample.combine(vmp, vmp.joint, "joint")
  }
  ##--------------------------------------------------------------
  ##
  ## obtain performance value - used for error rate confidence intervals
  ##
  ##--------------------------------------------------------------
  if (performance) {
    err <- get.mv.error(obj, FALSE, FALSE)
    vmp <- vimp.subsample.combine(vmp, err, "err")
  }
  ##----------------------------------------------------------
  ##
  ## subsampling loop for calculating VIMP confidence regions
  ##
  ##----------------------------------------------------------
  vmpS <- lapply(1:B, function(b) {
    ## progress bar
    if (verbose && B > 1) {
      custom.progress.bar(b, B)
    }
    ## draw the subsample
    ## use stratified sampling for classification/CR if requested (default)
    ## stratified sampling not allowed for mv-families
    if (!mv.fmly && stratify) {
      if (fmly == "class") {
        pt <- make.strat.sample(obj$yvar, subratio)
      }
      else if (fmly == "surv" || fmly == "surv-CR") {
        pt <- make.strat.sample(obj$yvar[, 2], subratio)
      }
      else {
        pt <- sample(1:n, size = (n * subratio), replace = FALSE)
      }
    }
    else {
      pt <- sample(1:n, size = (n * subratio), replace = FALSE)
    }
    ## calculate VMP on the subsampled data
    rf.b <-  do.call("rfsrc",
                     c(list(data = dta[pt,, drop = FALSE], importance = importance), rf.prms))
    vmp.b <- get.mv.vimp(rf.b, FALSE, FALSE)
    if (joint && !performance.only) {
      vmp.b.joint <- get.mv.vimp(get.subsample.joint.vimp(rf.b, rf.prms, xvar.names), FALSE, FALSE)
      vmp.b <- vimp.subsample.combine(vmp.b, vmp.b.joint, "joint")
    }
    if (performance) {
      err.b <- get.mv.error(rf.b, FALSE, FALSE)
      vmp.b <- vimp.subsample.combine(vmp.b, err.b, "err")
    }
    ## combine
    rO.b <- lapply(1:length(vmp), function(j) {
      vmp.j <- vmp[[j]]
      vmp.b.j <- vmp.b[[j]]
      vmp.mtx <- data.frame(matrix(NA, nrow(vmp.j), ncol(vmp.j)))
      rownames(vmp.mtx) <- rownames(vmp.j)
      colnames(vmp.mtx) <- colnames(vmp.j)
      if (ncol(vmp.b.j) == ncol(vmp.j)) {
        vmp.mtx[rownames(vmp.b.j), ] <- vmp.b.j
        vmp.b.j <- vmp.mtx
      }
      vmp.b.j
    })
    names(rO.b) <- names(vmp)
    rO.b
  })
  ##------------------------------------------------------------------
  ##
  ## return the goodies
  ##
  ##------------------------------------------------------------------
  rO <- list(rf = obj, vmp = vmp, vmpS = vmpS, subratio = subratio, performance.only = performance.only)
  class(rO) <- c(class(obj), "subsample")
  rO
}
subsample <- subsample.rfsrc
