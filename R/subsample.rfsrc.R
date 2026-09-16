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
  ## Validate the saved grow object before replaying its forest settings.
  ##
  ##--------------------------------------------------------------
  ## Validate the outer replicate count without changing the treatment of
  ## a one-replicate bank (its sample-SD summaries remain unavailable).
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 1 || B != floor(B) || B > .Machine$integer.max) {
    stop("'B' must be a positive integer.", call. = FALSE)
  }
  B <- as.integer(B)
  flags <- list(stratify = stratify, performance = performance,
                performance.only = performance.only, joint = joint,
                bootstrap = bootstrap, verbose = verbose)
  for (nm in names(flags)) {
    .subsample.check.flag(flags[[nm]], nm)
  }
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
    stop("this function does not work with anonymous forests")
  }
  ##--------------------------------------------------------------
  ##
  ## Multivariate families use unstratified outer subsampling.
  ##
  ##--------------------------------------------------------------
  fmly <- obj$family
  mv.fmly <- fmly == "regr+" || fmly == "class+" || fmly == "mix+" 
  ##--------------------------------------------------------------
  ##
  ## Set the nominal outer subsample fraction.
  ## This is distinct from the per-tree sampsize function used within
  ## each newly grown forest. The realized row count is not saved here.
  ##
  ##--------------------------------------------------------------
  n <- obj$n
  if (is.null(subratio)) {
    subratio <- get.subsample.subratio(n)
  }
  if (!bootstrap && (!is.numeric(subratio) || length(subratio) != 1L ||
                     !is.finite(subratio) || subratio <= 0 || subratio > 1)) {
    stop("'subratio' must be a finite number greater than zero and at most one.",
         call. = FALSE)
  }
  ##--------------------------------------------------------------
  ##
  ## Collect the settings explicitly replayed in each subsample forest.
  ## This list does not reconstruct every possible option in the grow call.
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
  ## The outer sampling scheme requires equal case weights.
  if (!all(diff(rf.prms$case.wt) == 0)) {
    stop("subsampling is not permitted on forests grown under non-standard case weights")
  }
  ## Equal case weights can be omitted from the new grow calls.
  rf.prms$case.wt <- NULL
  ## need time.interest if this is a survival object
  if (fmly == "surv" || fmly == "surv-CR") {
    rf.prms$ntime <- obj$time.interest
  }
  ## Omit uniform predictor/split weights; retain nonuniform weights.
  if (all(diff(rf.prms$xvar.wt) == 0)) {
    rf.prms$xvar.wt <- NULL
  }
  if (all(diff(rf.prms$split.wt) == 0)) {
    rf.prms$split.wt <- NULL
  }
  ## Omit the nonpositive sentinel for unrestricted tree depth.
  if (rf.prms$nodedepth <= 0) {
    rf.prms$nodedepth <- NULL
  }
  ## subsampling is not permitted if bootstrapping is not by root
  if (!bootstrap && rf.prms$bootstrap != "by.root") {
    stop("subsampling is not permitted on grow objects where sampling is not 'by.root'")
  }
  ## The saved sampsize function is evaluated on each smaller data set.
  ## It specifies the within-forest tree sample size, not the outer sample.
  ## for double bootstrapping we do NOT allow custom sampling for grow object
  if (bootstrap && rf.prms$bootstrap != "by.root") {
    stop("double bootstrapping is only permitted for breiman (bootstrap) forests")
  }
  else {
    ## Remove the explicit bootstrap setting. Ordinary refits use the
    ## grow default; bootsample() constructs its own sampling protocol.
    ## Retain samptype and sampsize for both paths.
    rf.prms$bootstrap <- NULL
  }
  ## Reuse the choice of Uno weighting, not the original row-wise weights.
  rf.prms$use.uno <- !is.null(obj$forest$uno.weights)
  ##--------------------------------------------------------------
  ##
  ## Performance-only requests suppress VIMP and append just the error row.
  ##
  ##--------------------------------------------------------------
  if (performance.only) {
    performance <- TRUE
    importance <- "none"
    vmp <- NULL
  }
  ##--------------------------------------------------------------
  ##
  ## Delegate double-bootstrap generation to bootsample(). Its extraction
  ## method receives the returned vmpB bank through extract.subsample().
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
  ## Reconstruct the row-level data used in outer subsampling.
  ## This uses the saved response/predictor tables before any VIMP restore.
  ##
  ##--------------------------------------------------------------
  dta <- data.frame(obj$yvar, obj$xvar, check.names = FALSE)
  colnames(dta)[1:length(obj$yvar.names)] <- obj$yvar.names
  ##--------------------------------------------------------------
  ##
  ## Establish full-data reference VIMP before generating replicates.
  ## Existing VIMP fixes the importance type and block size. Otherwise
  ## restore VIMP on the saved trees with the requested settings.
  ##
  ##--------------------------------------------------------------
  if (!performance.only) {
    vmp <- get.mv.vimp(obj, FALSE, FALSE)
    if (is.null(vmp)) {
      if (verbose) cat("no importance found: calculating it now ...\n")
      if (missing(importance)) {
        importance <- TRUE
      }
      obj <- vimp(obj,
                  perf.type = rf.prms$perf.type,
                  use.uno = rf.prms$use.uno,
                  block.size = block.size,
                  importance = importance)
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
  ## Optionally append a joint-VIMP row for the selected predictor group.
  ## This does not replace the individual-VIMP rows.
  ##
  ##--------------------------------------------------------------
  if (joint && !performance.only) {
    vmp.joint <- get.mv.vimp(get.subsample.joint.vimp(obj, rf.prms, xvar.names), FALSE, FALSE)
    vmp <- vimp.subsample.combine(vmp, vmp.joint, "joint")
  }
  ##--------------------------------------------------------------
  ##
  ## Optionally append the full-data performance statistic as row "err".
  ## get.mv.error(..., pretty = FALSE) preserves outcome-specific columns.
  ##
  ##--------------------------------------------------------------
  if (performance) {
    err <- get.mv.error(obj, FALSE, FALSE)
    vmp <- vimp.subsample.combine(vmp, err, "err")
  }
  ##----------------------------------------------------------
  ##
  ## Grow a new forest on each outer subsample and retain its statistics.
  ## Resampling is at the observation level, not over the saved trees.
  ##
  ##----------------------------------------------------------
  .subsample.require.reference(vmp)
  vmpS <- lapply(seq_len(B), function(b) {
    ## progress bar
    if (verbose && B > 1) {
      custom.progress.bar(b, B)
    }
    ## Stratify univariate classification by label and survival/competing
    ## risks by status (including censoring). Other paths sample all rows.
    ## Integer allocation within strata can omit very small strata.
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
    if (length(pt) == 0L) {
      stop("The subsample is empty. Increase 'subratio'.", call. = FALSE)
    }
    sample.data <- dta[pt, , drop = FALSE]
    ## Survival refits require at least one observed event after the existing
    ## omission step. This check does not change the sampled rows, the quotas,
    ## or missing-data handling, and is not applied to classification.
    if (fmly == "surv" || fmly == "surv-CR") {
      .subsample.check.events(sample.data, obj$yvar.names, b)
    }
    ## calculate VIMP on the subsampled data using the existing grow defaults
    rf.b <- do.call("rfsrc",
                    c(list(data = sample.data, importance = importance), rf.prms))
    vmp.b <- get.mv.vimp(rf.b, FALSE, FALSE)
    if (joint && !performance.only) {
      vmp.b.joint <- get.mv.vimp(get.subsample.joint.vimp(rf.b, rf.prms, xvar.names), FALSE, FALSE)
      vmp.b <- vimp.subsample.combine(vmp.b, vmp.b.joint, "joint")
    }
    if (performance) {
      err.b <- get.mv.error(rf.b, FALSE, FALSE)
      vmp.b <- vimp.subsample.combine(vmp.b, err.b, "err")
    }
    ## Match outcome, statistic-row, and class/event-column identities to
    ## the full-data reference. Extra replicate rows cannot expand it.
    .subsample.align.vimp(vmp, vmp.b)
  })
  ##------------------------------------------------------------------
  ##
  ## Return the full-data reference and B replicate statistic tables.
  ## Interval construction is deferred to extraction/printing/plotting.
  ##
  ##------------------------------------------------------------------
  rO <- list(rf = obj, vmp = vmp, vmpS = vmpS, subratio = subratio, performance.only = performance.only)
  class(rO) <- c(class(obj), "subsample")
  rO
}
subsample <- subsample.rfsrc
