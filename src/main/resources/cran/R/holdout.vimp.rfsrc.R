holdout.vimp.rfsrc <- function(formula, data,
                               ntree = function(p, vtry){1000 * p / vtry},
                               ntree.max = 2000,
                               ntree.allvars = NULL,
                               nsplit = 10,
                               ntime = 50,
                               mtry = NULL,
                               vtry = 1,
                               fast = FALSE,
                               verbose = TRUE, 
                               ...)
{
  ## --------------------------------------------------------------
  ##   
  ##   preliminary processing
  ##
  ## --------------------------------------------------------------
  ## pull unnamed parameters to pass to rfsrc
  ## some options cannot be used
  dots <- list(...)
  dots$formula <- dots$data <- dots$ntree <- dots$nsplit <-
    dots$ntime <- dots$mtry <- dots$importance <- dots$fast <- NULL
  ## hidden options to speed calculations - missing data not handled for speed
  #dots$terminal.quants <- dots$terminal.qualts <- TRUE
  #dots$na.action <- "na.omit"
  ## add vtry to the list
  dots$vtry <- vtry
  ## determine the grow interface - rfsrc or rfsrc.fast?
  if (!fast) {
    rfsrc.grow <- "rfsrc"
  }
  else {
    rfsrc.grow <- "rfsrc.fast"
  }
  ## get dimension p - make a fast single stump call for accurate dimension
  p <- length(rfsrc(formula, data, nodedepth = 1, ntree = 1, splitrule = "random")$xvar.names)
  ## process ntree - is it a function or number?
  if (!is.function(ntree) && !is.numeric(ntree)) {
    stop("ntree must be a function or number specifying requested number of grow trees")
  }
  if (is.function(ntree)) {
    ntree <- ntree(p, vtry)
  }
  ## specify the hold out array
  if (!is.null(ntree.allvars) && ntree.allvars <= 0) {
    ntree.allvars <- NULL
  }
  dots$holdout.array <- make.holdout.array(vtry, p, ntree, ntree.allvars)
  ntree <- ncol(dots$holdout.array)
  ## set block.size used for get.tree 
  block.size <- dots$block.size
  dots$block.size <- NULL
  ## TBD TBD TBD currently not in effect TBD TBD TBD
  block.size <- NULL
  ## --------------------------------------------------------------
  ##   
  ##   grow call - massive ntree
  ##
  ## --------------------------------------------------------------
  ## grow step
  ## set block.size to ntree to reduce computations
  ## the user specified block.size is only relevant for get.tree in predict
  o <- do.call(rfsrc.grow, c(list(formula = formula, data = data, ntree = ntree,
                 importance = "none", block.size = ntree, nsplit = nsplit, mtry = mtry), dots))
  ##--------------------------------------------------------------
  ##
  ## holdout step - get.tree calls - can be slow 
  ##
  ##--------------------------------------------------------------
  ## extract the holdout array
  ho <- o$holdout.array
  ## used when ntree.alvars is not NULL
  ## determine the columns where all variables were used
  ## pull the corresponding trees
  if (!is.null(ntree.allvars)) {
    v.in <- which(colSums(ho, na.rm = TRUE) == 0)
    pe.in  <- predict(o, get.tree = v.in, block.size = block.size)
  }
  ## pull the yvar names 
  ynms <- o$yvar.names
  if (o$family == "surv" || o$family == "surv-CR") {
    ynms <- ynms[1]
  }
  ## verbose setup
  if (verbose) pb <- txtProgressBar(min = 0, max = length(o$xvar.names), style = 3)
  ## -------------------------------------------------------
  ## 
  ## loop over features
  ##
  ## -------------------------------------------------------
  hvimp <- rbind(sapply(1:length(o$xvar.names), function(j) {  
    ## progress bar
    if (verbose && length(o$xvar.names) > 1) {
      setTxtProgressBar(pb, j)
    }
    if (verbose && j == length(o$xvar.names)) {
      cat("\n")
    }
    ## determine hold out trees for v - this may fail in extreme cases
    pt.out <- ho[j, ] == 1
    ## calculate PE for v held out and compare to all variables
    if (sum(pt.out) > 0 && (!is.null(ntree.allvars) || sum(!pt.out) > 0)) {
      ## define the out trees - restrict size of v.out using ntree.max
      ## extract the ensembles for hold out trees
      v.out <- which(pt.out)
      v.out <- sample(v.out)[1:min(ntree.max, length(v.out))]
      pe.out <- predict(o, get.tree = v.out, block.size = block.size)
      ## same as above for the "in" trees
      if (is.null(ntree.allvars)) {
        v.in <- which(!pt.out)
        v.in <- sample(v.in)[1:min(ntree.max, length(v.in))]
        pe.in <- predict(o, get.tree = v.in, block.size = block.size)
      }
      ## loop over outcomes - in case this is a multivariate object
      do.call(cbind, lapply(ynms, function(nn) {
        ## pull the hold out error rate
        o.out <- coerce.multivariate(pe.out, nn)
        err.out <- get.holdout.err(o.out, length(v.out))
        ## pull the error rate when v "is in"
        o.in <- coerce.multivariate(pe.in, nn)
        err.in <- get.holdout.err(o.in, length(v.in))
        ## here's the vimp
        vmp <- err.out - err.in
        dm <- length(get.holdout.err(coerce.multivariate(pe.out, nn), 1))
        if (length(vmp) == 0) {
          rep(NA, dm)
        }
        else {
          cbind(vmp)
        }
      }))
    }
    else {## no holdout v or hold in v - return NA
      do.call(cbind, lapply(ynms, function(nn) {
        dm <- length(get.holdout.err(coerce.multivariate(o, nn), 1))
        rep(NA, dm)
      }))
    }
  }))
  ## pretty up for return
  if (nrow(hvimp) == 1) {
    hvimp <- c(hvimp)
    names(hvimp) <- o$xvar.names
  }
  else {
    if (o$family == "surv-CR") {
      rownames(hvimp) <- colnames(o$err.rate)
    }
    else {
      rownames(hvimp) <- ynms
    }
    colnames(hvimp) <- o$xvar.names
  }
  ## return the goodies
  hvimp
}
holdout.vimp <- holdout.vimp.rfsrc
## --------------------------------------------------------------
##
## function to extract holdout error rate
##
## --------------------------------------------------------------
get.holdout.err <- function(obj, j) {
  err <- obj$err.rate
  if (obj$family == "class") {
    mean(cbind(err)[1:j, 1], na.rm = TRUE)
  }
  else if (obj$family == "surv-CR") {
    colMeans(err[1:j,, drop = FALSE], na.rm = TRUE) 
  }
  else {
    mean(err[1:j], na.rm = TRUE)
  }
}
