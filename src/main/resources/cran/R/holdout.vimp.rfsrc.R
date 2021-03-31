holdout.vimp.rfsrc <- function(formula, data,
                               ntree = function(p, vtry){1000 * p / vtry},
                               nsplit = 10,
                               ntime = 50,
                               sampsize = function(x){x * .632},
                               samptype = "swor",
                               block.size = 10,
                               vtry = 1,
                               ...)
{
  ## --------------------------------------------------------------
  ##   
  ##   preliminary processing
  ##
  ## --------------------------------------------------------------
  ## pull unnamed parameters to pass to rfsrc
  ## lock down options that cannot be used
  dots <- list(...)
  dots$formula <- dots$data <- dots$ntree <- dots$nsplit <- dots$ntime <-
    dots$samp <- dots$sampsize <- dots$samptype <- dots$block.size <-
      dots$importance <- NULL
  dots$terminal.qualts <- FALSE
  ## get dimension p - make a fast single stump call for accurate dimension
  p <- length(rfsrc(formula, data, nodedepth = 1, ntime = 1,
              ntree = 1, splitrule = "random")$xvar.names)
  ## process ntree - is it a function or number?
  ## targeted hold out vimp must be handled carefully - checks are put in place
  if (!is.function(ntree) && !is.numeric(ntree)) {
    stop("ntree must be a function or number specifying requested number of grow trees")
  }
  if (is.function(ntree)) {
    if (!is.list(vtry)) {
      ntree <- ntree(p, vtry)
    }
    ## special case when vtry is a list and user wants targeted hold out vimp
    else {
      xvar <- unlist(vtry[1])
      if (any(xvar < 1) | any(xvar > p)) {
        stop("x-variables specified in vtry should be integer values between 1 and p=number of features")
      }      
      joint <- unlist(vtry[2])
      if (!is.logical(joint)) {
        stop("joint vimp in vtry must be specified as being either TRUE or FALSE")
      }
      if (!joint) {
        ntree <- ntree(length(xvar), 1)
      }
      else {
        if (length(xvar) == p) {
          stop("number of hold out variables for joint analysis must be less than number of features")
        }
        ntree <- ntree(1, 1)
      }
    }
  }
  ## initialize the holdout specs object
  holdout.specs <- get.holdout.specs(vtry, p, ntree, block.size, mode = "baseline")
  ## --------------------------------------------------------------
  ##   
  ##   There are two native-code calls to affect holdout vimp:
  ##
  ##   The first call grows a baseline forest.  This is called
  ##   "baseline" mode, and is specificaly recognized by the native
  ##   code via the incoming holdout specs.  In this mode, an x-var IS
  ##   NOT HELD OUT during growth, but performance is calculated
  ##   over the relevant blocks as if the x-vars WAS held out.
  ##
  ##   The second call grows a holdout forest.  This is called
  ##   "holdout" mode, and is specificaly recognized by the native
  ##   code via the incoming parameters.  In this mode, an x-var IS
  ##   HELD OUT during growth, and performance is calculated over
  ##   the relevant blocks.
  ##
  ##   Note that the trees over which a particular x-var is held
  ##   are determined by the holdout array.  Theses trees are
  ##   divided in blocks of block.size.  Incomplete blocks are
  ##   ignored in the resulting performance calculation.
  ##   Blocks are non-contiguous and specific to each x-var.
  ##
  ## --------------------------------------------------------------
  ## Note that we use the same seed for both forests, in an attempt
  ## to reduce variance. This will only result in identical forests
  ## for both modes if the holdout array is the zero matrix, there
  ## is no missing data, and we do deterministic splitting.  It's a
  ## little complicated, and we should probably test the for
  ## variance without identical seeds.  In fact, in the presence of
  ## missing data, it might be the choice of identical seeds is
  ## irrelevant.
  seed <- dots$seed
  dots$seed <- NULL
  ## initialize the seed
  seed <- get.seed(seed)
  ## baseline (harness) forest call:
  base.obj <- do.call(rfsrc, c(list(formula = formula,
                                    data = data,
                                    ntree = ntree,
                                    nsplit = nsplit,
                                    ntime = ntime,
                                    sampsize = sampsize,
                                    samptype = samptype,
                                    holdout.specs = holdout.specs,
                                    importance = FALSE,
                                    forest = FALSE,
                                    terminal.qualts = FALSE,
                                    seed = seed),
                               dots))
  ## change the specifiction object to holdout mode. 
  holdout.specs <- switch.holdout.specs(holdout.specs, "holdout")
  ## holdout forest call:    
  hold.obj <- do.call(rfsrc, c(list(formula = formula,
                                    data = data,
                                    ntree = ntree,
                                    nsplit = nsplit,
                                    ntime = ntime,
                                    sampsize = sampsize,
                                    samptype = samptype,
                                    holdout.specs = holdout.specs,
                                    importance = FALSE,
                                    forest = FALSE,
                                    terminal.qualts = FALSE,
                                    seed = seed),
                               dots))
  ## We are now in a position to analyze the resulting blocked vimp.
  ## The delta between the two call will give you the holdout vimp.
  ## The objects are multivariate compliant and return the performance in
  ## the response specific lists in the multivariate case.  They follow
  ## the simplified list format in the univaratiate case.
  ## Note: obj$holdout.blk is a vector of length [p] indicating how many blocks
  ## each x-var contains.  This is variable and might be zero if the x-var is
  ## never held out.
  ## --------------------------------------------------------------    
  ## The outputs are as follows:
  ## SURV:  obj$holdout.vimp
  ## To the R code:
  ##   -> of dim [[p]] x [RF_eventTypeSize] x [holdout.blk[.]]
  ##              ^^^
  ##   list element can be null
  ##   if no blocks for this x-var  
  ## CLAS:  obj$holdout.vimp
  ## To the R code:
  ##   -> of dim [[p]] x [1 + levels.count[.]] x [holdout.blk[.]]
  ##              ^^^
  ##   list element can be null
  ##   if no blocks for this x-var  
  ## REGR:  obj$holdout.vimp
  ## To the R code:
  ##   -> of dim [[p]] x [holdout.blk[.]]
  ##              ^^^
  ##   list element can be null
  ##   if no blocks for this x-var  
  ##
  ## --------------------------------------------------------------
  ##cat("\nBlocks in each x-var:  \n")
  ##print(base.obj$holdout.blk)
  mv.flag <- FALSE
  ##--------------------------------
  ##
  ## multivariate object
  ##
  ##--------------------------------
  if ((class(base.obj)[3] == "regr+") | (class(base.obj)[3] == "class+") | (class(base.obj)[3] == "mix+")) { 
    ## set the mv flag
    mv.flag <- TRUE
    ## create the holdout vimp object
    hvimp <- list()
    ## regr outcomes
    if (length(base.obj$regrOutput) > 0) {
      hvimp$regrOutput <- lapply(1:length(base.obj$regrOutput), function(i) {
        hv <- sapply(1:length(base.obj$regrOutput[[i]]$holdout.vimp), function(j) {
          mean(hold.obj$regrOutput[[i]]$holdout.vimp[[j]] - base.obj$regrOutput[[i]]$holdout.vimp[[j]], na.rm = TRUE)
        })
        names(hv) <- names(base.obj$regrOutput[[i]]$holdout.vimp)
        hv
      })
      names(hvimp$regrOutput) <- names(base.obj$regrOutput)
    }
    ## class outcomes
    if (length(base.obj$classOutput) > 0) {
      hvimp$classOutput <- lapply(1:length(base.obj$classOutput), function(i) {
        yval.i <- colnames(base.obj$classOutput[[i]]$predicted)
        hv <- do.call(rbind, lapply(1:length(base.obj$classOutput[[i]]$holdout.vimp), function(j) {
          if (!is.null(dim(base.obj$classOutput[[i]]$holdout.vimp[[j]]))) {
            rowMeans(hold.obj$classOutput[[i]]$holdout.vimp[[j]] - base.obj$classOutput[[i]]$holdout.vimp[[j]], na.rm = TRUE)
          }
          else {
            rep(NA, length(yval.i) + 1)
          }
        }))
        colnames(hv) <- c("all", yval.i)
        rownames(hv) <- names(base.obj$classOutput[[i]]$holdout.vimp)
        hv
      })
      names(hvimp$classOutput) <- names(base.obj$classOutput)
    }
  }
  ##--------------------------------
  ##
  ## regression object
  ##
  ##--------------------------------
  else if (class(base.obj)[3] == "regr") {  
    hvimp <- sapply(1:length(base.obj$holdout.vimp), function(j) {
      if (!is.null(base.obj$holdout.vimp[[j]])) {
        mean(hold.obj$holdout.vimp[[j]] - base.obj$holdout.vimp[[j]], na.rm = TRUE)
      }
      else {
        NA
      }
    })
    names(hvimp) <- names(base.obj$holdout.vimp)
  }
  ##--------------------------------
  ##
  ## classification object
  ##
  ##--------------------------------
  else if (class(base.obj)[3] == "class") {  
    yval <- levels(base.obj$yvar)
    hvimp <- do.call(rbind, lapply(1:length(base.obj$holdout.vimp), function(j) {
      if (!is.null(dim(base.obj$holdout.vimp[[j]]))) {
        rowMeans(hold.obj$holdout.vimp[[j]] - base.obj$holdout.vimp[[j]], na.rm = TRUE)
      }
      else {
        rep(NA, length(yval) + 1)
      }
    }))
    colnames(hvimp) <- c("all", yval)
    rownames(hvimp) <- names(base.obj$holdout.vimp)
  }
  ##--------------------------------
  ##
  ## survival, CR objects
  ##
  ##--------------------------------
  else if (class(base.obj)[3] == "surv") {
    hvimp <- sapply(1:length(base.obj$holdout.vimp), function(j) {
      if (!is.null(base.obj$holdout.vimp[[j]])) {
        mean(hold.obj$holdout.vimp[[j]] - base.obj$holdout.vimp[[j]], na.rm = TRUE)
      }
      else {
        rep(NA)
      }
    })
    names(hvimp) <- names(base.obj$holdout.vimp)
  }
  else if (class(base.obj)[3] == "surv-CR") {
    hvimp <- do.call(rbind, lapply(1:length(base.obj$holdout.vimp), function(j) {
      if (!is.null(dim(base.obj$holdout.vimp[[j]]))) {
        rowMeans(hold.obj$holdout.vimp[[j]] - base.obj$holdout.vimp[[j]], na.rm = TRUE)
      }
      else {
        rep(NA, length(base.obj$holdout.vimp))
      }
    }))
    colnames(hvimp) <- paste("event.", 1:length(get.event.info(base.obj)$event.type), sep = "")
    rownames(hvimp) <- names(base.obj$holdout.vimp)
  }
  ## return the promised holdout vimp
  ## processing is different for multivariate families
  rO <- list()
  if (mv.flag) {
    rO$baseline <- rO$holdout <- rO$importance <- list()
    if (length(base.obj$regrOutput) > 0) {
      rO$baseline$regrOutput <- lapply(base.obj$regrOutput, function(o) {o$holdout.vimp})
      rO$holdout$regrOutput <- lapply(hold.obj$regrOutput, function(o) {o$holdout.vimp})
      rO$importance$regrOutput <- hvimp$regrOutput
    }
    if (length(base.obj$classOutput) > 0) {
      rO$baseline$classOutput <- lapply(base.obj$classOutput, function(o) {o$holdout.vimp})
      rO$holdout$classOutput <- lapply(hold.obj$classOutput, function(o) {o$holdout.vimp})
      rO$importance$classOutput <- hvimp$classOutput
    }
  }
  else {
    rO$baseline <- base.obj$holdout.vimp
    rO$holdout <- hold.obj$holdout.vimp
    rO$importance <- hvimp
  }
  return(invisible(rO))
}
holdout.vimp <- holdout.vimp.rfsrc
## --------------------------------------------------------------
##
## function to create holdout native code parameters
## - no checks on block.size are made
##
## --------------------------------------------------------------
get.holdout.specs <- function(vtry, p, ntree, block.size, mode) {
  ## make holdout array
  ## check that vtry is not a list --> default case and default processing
  if (!is.list(vtry)) {
    if (vtry > 0) {
      holdout <- make.holdout.array(vtry, p, ntree, NULL)
    }
    else {
      stop("vtry must be positive")
    }
  }
  ## vtry is a list --> user wants targeted hold out vimp
  else {
    ## parse vtry
    xvar <- unlist(vtry[1])
    joint <- unlist(vtry[2])
    ## set vtry to integer value
    if (!joint) {
      vtry <- 1
    }
    else {
      vtry <- length(xvar)
    }
    if (!joint) {
      holdout <- do.call(cbind, lapply(1:ntree, function(b) {
        ho <- rep(0, p)
        ho[resample(xvar, size = 1, replace = FALSE)] <- 1
        ho
      }))
    }
    else {
      holdout <- matrix(0, p, ntree)
      holdout[xvar, ] <- 1
    }
  }
  ## assemble the return object and return it
  obj <- list(as.integer(vtry), as.integer(holdout), as.integer(block.size), as.integer(0))
  names(obj) = c("vtry", "holdout", "block.size", "mode")
  return(switch.holdout.specs(obj, mode))
}
## --------------------------------------------------------------
##
## function to switch holdout mode for native code parameters
## mode changes between "baseline" <-> "holdout"
##
## --------------------------------------------------------------
switch.holdout.specs<- function(obj, mode) {
  if (mode == "baseline") {
    new.mode = 1
  }
  else if (mode == "holdout") {
    new.mode = 2
  }
  else if (mode == "natural") {
    new.mode = 3
  }
  else {
    stop("vtry mode invalid")
  }
  obj$mode = as.integer(new.mode)
  return(obj)
}
