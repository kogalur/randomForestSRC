subsample.rfsrc <- function(obj,                
                            B = 100,
                            block.size = 1,
                            subratio = NULL,
                            stratify = TRUE,
                            performance = FALSE,
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
  if (!bootstrap && is.null(subratio)) {
    subratio <- 1 / sqrt(n)
  }
  if (!bootstrap && (subratio < 0 || subratio > 1)) {
    stop("subratio must be between 0 and 1:", subratio)
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
  ## call the bootstrap subroutine if double bootstrap was requested
  ##
  ##--------------------------------------------------------------
  if (bootstrap) {
    bootO <- bootsample(obj, rf.prms, B = B, block.size = block.size,
                        joint = joint, xvar.names = xvar.names,
                        performance = performance, verbose = verbose)
    rO <- list(rf = bootO$rf, vmp = bootO$vmp, vmpB = bootO$vmpB, subratio = NULL)
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
  vmp <- get.mv.vimp(obj, FALSE, FALSE)
  if (is.null(vmp)) {
    if (verbose) cat("no importance found: calculating it now ...\n")
    obj <- vimp(obj, perf.type = rf.prms$perf.type, block.size = block.size)
    vmp <- get.mv.vimp(obj, FALSE, FALSE)
    rf.prms$block.size <- block.size
    if (verbose) cat("done\n")
  }
  ##--------------------------------------------------------------
  ##
  ## call joint vimp - reference vimp value for pure noise settings
  ##
  ##--------------------------------------------------------------
  if (joint) {
    vmp.joint <- get.mv.vimp(get.joint.vimp(obj, rf.prms, xvar.names), FALSE, FALSE)
    vmp <- vimp.combine(vmp, vmp.joint, "joint")
  }
  ##--------------------------------------------------------------
  ##
  ## obtain performance value - used for error rate confidence intervals
  ##
  ##--------------------------------------------------------------
  if (performance) {
    err <- get.mv.error(obj, FALSE, FALSE)
    vmp <- vimp.combine(vmp, err, "err")
  }
  ##----------------------------------------------------------
  ##
  ## subsampling loop for calculating VIMP confidence regions
  ##
  ##----------------------------------------------------------
  if (verbose) pb <- txtProgressBar(min = 0, max = B, style = 3)
  vmpS <- lapply(1:B, function(b) {
    ## progress bar
    if (verbose && B > 1) {
      setTxtProgressBar(pb, b)
    }
    if (verbose && b == B) {
      cat("\n")
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
    ## calculate VMP on the stratified data
    rf.b <-  do.call("rfsrc",
                     c(list(data = dta[pt,, drop = FALSE], importance = TRUE), rf.prms))
    vmp.b <- get.mv.vimp(rf.b, FALSE, FALSE)
    if (joint) {
      vmp.b.joint <- get.mv.vimp(get.joint.vimp(rf.b, rf.prms, xvar.names), FALSE, FALSE)
      vmp.b <- vimp.combine(vmp.b, vmp.b.joint, "joint")
    }
    if (performance) {
      err.b <- get.mv.error(rf.b, FALSE, FALSE)
      vmp.b <- vimp.combine(vmp.b, err.b, "err")
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
  rO <- list(rf = obj, vmp = vmp, vmpS = vmpS, subratio = subratio)
  class(rO) <- c(class(obj), "subsample")
  rO
}
subsample <- subsample.rfsrc
###################################################################
##
## extract confidence interval + other goodies from subsampled object
##
###################################################################
extract.subsample <- function(obj, alpha = .05, target = 0, m.target = NULL, standardize = TRUE)
{
  ## in case a user applies this to a bootsample object
  if (sum(c(grepl("rfsrc", class(obj))), grepl("subsample", class(obj))) == 2) {
    subsample <- TRUE
  }
  else if (sum(c(grepl("rfsrc", class(obj))), grepl("bootsample", class(obj))) == 2) {
    subsample <- FALSE
  }
  else {
    stop("this must be an object obtained by calling the 'subsample' function")
  }
  if (!subsample) {
    return(extract.bootsample(obj, alpha = alpha,
              target = target, m.target = m.target, standardize = standardize))    
  }
  ## coerce the (potentially) multivariate rf object
  m.target <- get.univariate.target(obj$rf, m.target)
  obj$rf <- coerce.multivariate(obj$rf, m.target)
  ## coerce the (potentially) multivariate original grow vimp
  if (is.null(m.target)) {
    vmp <- obj$vmp[[1]][, 1 + target, drop = FALSE]
  }
  else {
    vmp <- obj$vmp[[m.target]][, 1 + target, drop = FALSE]
  }
  ## coerce the (potentially) multivariate subsampled vimp 
  if (is.null(m.target)) {
    obj$vmpS <- lapply(obj$vmpS, data.frame)
  }
  else {
    obj$vmpS <- lapply(obj$vmpS, function(oo) {
      oo[[m.target]]
    })
  }  
  ## extract necessary objects
  theta.hat <- get.standardize.vimp(obj$rf, vmp, standardize)[, 1]
  theta.star <- get.standardize.vimp(obj$rf, rbind(sapply(obj$vmpS, "[[", 1 + target)), standardize)
  rownames(theta.star) <- xvar.names <- rownames(vmp)
  ## sample sizes
  n <- obj$rf$n
  m <- obj$subratio * n
  ## boxplot subsampling data
  scl <- sqrt(obj$subratio)
  theta.star.scl <- theta.hat - scl * rowMeans(theta.star, na.rm = TRUE)
  boxplot.dta <- data.frame(t(scl * theta.star + theta.star.scl))
  ## subsampled z-statistic 
  z.star <- sqrt(m) * (theta.star - theta.hat) 
  ## asymptotic standard error
  se.Z <- sqrt(m / n) * apply(theta.star, 1, sd, na.rm = TRUE)
  ## asymptotic JK standard error
  se.jk.Z <- sqrt((m / n) * rowMeans((theta.star - theta.hat)^2, na.rm = TRUE))
  ## nonparametric (default) subsampling confidence interval  
  ci <- do.call(cbind, lapply(1:nrow(z.star), function(rowj) {
    c(theta.hat[rowj] - quantile(z.star[rowj, ], 1 - alpha/2, na.rm = TRUE) / sqrt(n),
      theta.hat[rowj] - quantile(z.star[rowj, ], .75, na.rm = TRUE) / sqrt(n),
      theta.hat[rowj] - quantile(z.star[rowj, ], .50, na.rm = TRUE) / sqrt(n),
      theta.hat[rowj] - quantile(z.star[rowj, ], .25, na.rm = TRUE) / sqrt(n),
      theta.hat[rowj] - quantile(z.star[rowj, ], alpha/2, na.rm = TRUE) / sqrt(n))
  }))
  rownames(ci) <- paste(round(100 * c(alpha/2,.25,.50,.75,1-alpha/2), 1), "%", sep = "")
  colnames(ci) <- xvar.names
  ## normal-Z subsampling confidence interval
  ci.Z <- do.call(cbind, lapply(1:nrow(z.star), function(rowj) {
  c(theta.hat[rowj] - qnorm(1 - alpha/2) * se.Z[rowj],
    theta.hat[rowj] - qnorm(.75) * se.Z[rowj],
    theta.hat[rowj],
    theta.hat[rowj] + qnorm(.75) * se.Z[rowj],
    theta.hat[rowj] + qnorm(1 - alpha/2) * se.Z[rowj])
  }))
  rownames(ci.Z) <- paste(round(100 * c(alpha/2,.25,.50,.75,1-alpha/2), 1), "%", sep = "")
  colnames(ci.Z) <- xvar.names
  ## normal-Z subsampling confidence interval with JK
  ci.jk.Z <- do.call(cbind, lapply(1:nrow(z.star), function(rowj) {
  c(theta.hat[rowj] - qnorm(1 - alpha/2) * se.jk.Z[rowj],
    theta.hat[rowj] - qnorm(.75) * se.jk.Z[rowj],
    theta.hat[rowj],
    theta.hat[rowj] + qnorm(.75) * se.jk.Z[rowj],
    theta.hat[rowj] + qnorm(1 - alpha/2) * se.jk.Z[rowj])
  }))
  rownames(ci.jk.Z) <- paste(round(100 * c(alpha/2,.25,.50,.75,1-alpha/2), 1), "%", sep = "")
  colnames(ci.jk.Z) <- xvar.names
  ## nonparametric variable selection object
  var.sel <- data.frame(lower  =  ci[1, ],
                        median =  ci[3, ],
                        upper  =  ci[5, ],
                        pvalue = rowMeans((theta.star - sqrt(n / m) * theta.hat) > 0, na.rm = TRUE),
                        signif =  ci[1, ] > 0)
  rownames(var.sel) <- xvar.names
  ## normal-Z variable selection object
  var.sel.Z <- data.frame(lower  =  ci.Z[1, ],
                          mean   =  ci.Z[3, ],
                          upper  =  ci.Z[5, ],
                          pvalue = pnorm(theta.hat, 0, se.Z, FALSE),
                          signif =  ci.Z[1, ] > 0)
  rownames(var.sel.Z) <- xvar.names
  ## normal-Z variable selection object with JK
  var.jk.sel.Z <- data.frame(lower  =  ci.jk.Z[1, ],
                          mean   =  ci.jk.Z[3, ],
                          upper  =  ci.jk.Z[5, ],
                          pvalue = pnorm(theta.hat, 0, se.jk.Z, FALSE),
                          signif =  ci.jk.Z[1, ] > 0)
  rownames(var.jk.sel.Z) <- xvar.names
  list(boxplot.dta = boxplot.dta,
       ci = ci,
       ci.Z = ci.Z,
       ci.jk.Z = ci.jk.Z,
       vmp = theta.hat,
       vmpS = theta.star,
       se.Z = se.Z,
       se.jk.Z = se.jk.Z,
       var.sel = var.sel,
       var.sel.Z = var.sel.Z,
       var.jk.sel.Z = var.jk.sel.Z)
}
###################################################################
##
## get joint vimp from a forest object
##
###################################################################
get.joint.vimp <- function(o, rf.prms = NULL, xvar.names) {
  class(o)[2] <- "grow"  ##fake the class to trick vimp()
  if (is.null(rf.prms)) {
    return(vimp(o, joint = TRUE, xvar.names = xvar.names))
  }
  else {
    return(vimp(o, joint = TRUE, xvar.names = xvar.names,
                perf.type = rf.prms$perf.type, block.size = rf.prms$block.size))
  }
}
###################################################################
##
## get standarized vimp from a matrix (p x B) of subsampled vimp
##
###################################################################
get.standardize.vimp <- function(obj, vmp, standardize = FALSE) {
  ## nothing to do - standardization not rquested
  if (!standardize) {
    return(vmp)
  }
  else {
    do.call(cbind, lapply(1:ncol(vmp), function(j) {
      v <- vmp[, j]
      if (obj$family != "regr") {
        100 * v 
      }
      else {
        100 * v / var(obj$yvar, na.rm = TRUE)
      }
    }))
  }
}
###################################################################
##
## double bootstrap sampling
##
###################################################################
make.double.boot.sample <- function(ntree, n, smp, size = function(x) {x}, replace = TRUE)
{
  dbl.smp <- do.call(cbind, mclapply(1:ntree, function(bb) {
    inb <- rep(0, n)
    smp.2 <- sample(smp, size = size(n), replace = replace)
    frq <- tapply(smp.2, smp.2, length)
    idx <- as.numeric(names(frq))
    inb[idx] <- frq
    inb
  }))
  if (length(setdiff(1:n, smp)) > 0) {
    dbl.smp[-setdiff(1:n, smp),, drop = FALSE]
  }
  else {
    dbl.smp
  }
}
###################################################################
##
## stratified sampling without replacement 
##
###################################################################
make.strat.sample <- function(y, subratio) {
  frq <- table(y)
  class.labels <- names(frq)
  as.numeric(unlist(sapply(class.labels, function(cl) {
    pt <- which(y == cl)
    sample(pt, size = length(pt) * subratio, replace = FALSE)    
  })))
}
###################################################################
##
## print function for subsampled ci
##
###################################################################
print.subsample.rfsrc <- function(x, alpha = .05, standardize = TRUE) {
  vmp <- x$vmp
  nullO <- lapply(1:length(vmp), function(j) {
    m.target <- names(vmp)[j]
    p.vmp <- ncol(vmp[[j]])
    vmp.col.names <- colnames(vmp[[j]])
    if (length(vmp) > 1) {
      cat("processing a multivariate family, outcome:", m.target, "\n")
    }
    lapply(1:p.vmp, function(k) {
      oo <- extract.subsample(x, alpha = alpha, target = k - 1,
                  m.target = m.target, standardize = standardize)
      cat("===== VIMP confidence regions for", vmp.col.names[k], " =====\n")
      cat("nonparametric:\n")
      print(round(oo$ci, 3))
      cat("parametric:\n")
      print(round(oo$ci.Z, 3))
      cat("parametric (jackknife):\n")
      print(round(oo$ci.jk.Z, 3))
      NULL
    })
  })
}
print.subsample <- print.subsample.rfsrc
###################################################################
##
## combine two lists of vimp: o1 is master, o2 is appended
##
###################################################################
vimp.combine <- function(o1, o2, o2.name = NULL) {
  rO <- lapply(1:length(o1), function(j) {
    rO.j <- rbind(o1[[j]], o2[[j]])
    if (!is.null(o2.name)) {
      rownames(rO.j)[nrow(rO.j)] <- o2.name
    }
    rO.j
  })
  names(rO) <- names(o1)
  rO
}
###################################################################
##
## core bootstrap function
##
###################################################################
bootsample <- function(obj, rf.prms, B, block.size, joint, xvar.names, performance, verbose) {
  ##--------------------------------------------------------------
  ##
  ## extract data from the previously grown forest
  ## set the dimensions, pull the feature names
  ##
  ##--------------------------------------------------------------
  dta <- data.frame(obj$yvar, obj$xvar)
  colnames(dta)[1:length(obj$yvar.names)] <- obj$yvar.names
  n <- nrow(obj$xvar)
  ##--------------------------------------------------------------
  ##
  ## call vimp to extract importance if not available in original object
  ##
  ##--------------------------------------------------------------
  vmp <- get.mv.vimp(obj, FALSE, FALSE)
  if (is.null(vmp)) {
    if (verbose) cat("no importance found: calculating it now ...\n")
    obj <- vimp(obj, perf.type = rf.prms$perf.type, block.size = block.size)
    vmp <- get.mv.vimp(obj, FALSE, FALSE)
    rf.prms$block.size <- block.size
    if (verbose) cat("done\n")
  }
  ##--------------------------------------------------------------
  ##
  ## call joint vimp - reference vimp value for pure noise settings
  ##
  ##--------------------------------------------------------------
  if (joint) {
    vmp.joint <- get.mv.vimp(get.joint.vimp(obj, rf.prms, xvar.names), FALSE, FALSE)
    vmp <- vimp.combine(vmp, vmp.joint)
  }
  ##--------------------------------------------------------------
  ##
  ## obtain performance value - used for error rate confidence intervals
  ##
  ##--------------------------------------------------------------
  if (performance) {
    err <- get.mv.error(obj, FALSE, FALSE)
    vmp <- vimp.combine(vmp, err, "err")
  }
  ##----------------------------------------------------------
  ##
  ## bootstrap loop for calculating VIMP confidence regions
  ##
  ##----------------------------------------------------------
  if (verbose) pb <- txtProgressBar(min = 0, max = B, style = 3)
  vmpB <- lapply(1:B, function(b) {
    ## progress bar
    if (verbose && B > 1) {
      setTxtProgressBar(pb, b)
    }
    if (verbose && b == B) {
      cat("\n")
    }
    ##---------------------------------------------------
    ##
    ## double bootstrap 
    ##
    ##---------------------------------------------------
    ## double bootstrap sample
    smp <- sample(1:n, size = n, replace = TRUE)
    smp.unq <- sort(unique(smp))
    dbl.smp <- make.double.boot.sample(rf.prms$ntree, n, smp,
                             size = rf.prms$sampsize, replace = rf.prms$samptype == "swr")
    ## pull VIMP from the double boostrap forest
    rf.b <- do.call("rfsrc",
     c(list(data = dta[smp.unq,, drop = FALSE], importance = TRUE, bootstrap = "by.user", samp = dbl.smp), rf.prms))
    vmp.b <- get.mv.vimp(rf.b, FALSE, FALSE)
    if (joint) {
      vmp.b.joint <- get.mv.vimp(get.joint.vimp(rf.b, rf.prms, xvar.names), FALSE, FALSE)
      vmp.b <- vimp.combine(vmp.b, vmp.b.joint)
    }
    if (performance) {
      err.b <- get.mv.error(rf.b, FALSE, FALSE)
      vmp.b <- vimp.combine(vmp.b, err.b, "err")
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
  ## return the potentially updated forest object
  ## return the vimp 
  list(rf = obj, vmp = vmp, vmpB = vmpB)
}
###################################################################
##
## extract confidence interval + other goodies from bootstrap object
##
###################################################################
extract.bootsample <- function(obj, alpha = .05, target = 0, m.target = NULL, standardize = TRUE)
{
  ## confirm this is the right kind of object
  if (!sum(c(grepl("rfsrc", class(obj))), grepl("bootsample", class(obj))) == 2) {
    stop("this must be a double-bootstrap object obtained by calling the 'subsample' function")
  }
  ## coerce the (potentially) multivariate rf object
  m.target <- get.univariate.target(obj$rf, m.target)
  obj$rf <- coerce.multivariate(obj$rf, m.target)
  ## coerce the (potentially) multivariate vimp object
  if (is.null(m.target)) {
    obj$vmpB <- lapply(obj$vmpB, data.frame)
  }
  else {
    obj$vmpB <- lapply(obj$vmpB, function(oo) {
      oo[[m.target]]
    })
  }  
  ## extract vimp
  theta.boot <- get.standardize.vimp(obj$rf, rbind(sapply(obj$vmpB, "[[", 1 + target)), standardize)
  rownames(theta.boot) <- rownames(obj$vmp[[1]][, 1 + target, drop = FALSE])
  ## bootstrap VIMP
  vmp.boot <- rowMeans(theta.boot, na.rm = TRUE)
  ## bootstrap standard error
  se.boot <- apply(theta.boot, 1, sd, na.rm = TRUE)
  ## nonparametric bootstrap confidence interval
  ci.boot <- do.call(cbind, lapply(1:length(vmp.boot), function(rowj) {
    c(quantile(theta.boot[rowj,], alpha/2, na.rm = TRUE),
      quantile(theta.boot[rowj,], .25, na.rm = TRUE),
      vmp.boot[rowj],
      quantile(theta.boot[rowj,], .75, na.rm = TRUE),
      quantile(theta.boot[rowj,], 1 - alpha/2, na.rm = TRUE))
  }))
  rownames(ci.boot) <- paste(round(100 * c(alpha/2,.25,.50,.75,1-alpha/2), 1), "%", sep = "")
  colnames(ci.boot) <- names(vmp.boot)
  ## parametric bootstrap confidence interval
  ci.boot.Z <- do.call(cbind, lapply(1:length(vmp.boot), function(rowj) {
    c(vmp.boot[rowj] - qnorm(1 - alpha/2) * se.boot[rowj],
      vmp.boot[rowj] - qnorm(.75) * se.boot[rowj],
      vmp.boot[rowj],
      vmp.boot[rowj] + qnorm(.75) * se.boot[rowj],
      vmp.boot[rowj] + qnorm(1 - alpha/2) * se.boot[rowj])
  }))
  rownames(ci.boot.Z) <- paste(round(100 * c(alpha/2,.25,.50,.75,1-alpha/2), 1), "%", sep = "")
  colnames(ci.boot.Z) <- names(vmp.boot)
  ## bootstrap variable selection object
  var.sel.boot <- data.frame(lower  =  ci.boot[1, ],
                             mean   =  ci.boot[3, ],
                             upper  =  ci.boot[5, ],
                             signif =  ci.boot[1, ] > 0)
  rownames(var.sel.boot) <- names(vmp.boot)
  var.sel.boot.Z <- data.frame(lower  =  ci.boot.Z[1, ],
                               mean   =  ci.boot.Z[3, ],
                               upper  =  ci.boot.Z[5, ],
                               pvalue = pnorm(vmp.boot, 0, se.boot, FALSE),
                               signif =  ci.boot.Z[1, ] > 0)
  rownames(var.sel.boot.Z) <- names(vmp.boot)
  list(ci = ci.boot,
       ci.Z = ci.boot.Z,
       vmp <- vmp.boot,
       vmpS = theta.boot,
       se = se.boot,
       var.sel = var.sel.boot,
       var.sel.Z = var.sel.boot.Z)
}
###################################################################
##
## print function from bootstrap ci
##
###################################################################
print.bootsample.rfsrc <- function(x, alpha = .05, standardize = TRUE) {
  vmp <- x$vmp
  nullO <- lapply(1:length(vmp), function(j) {
    m.target <- names(vmp)[j]
    p.vmp <- ncol(vmp[[j]])
    vmp.col.names <- colnames(vmp[[j]])
    if (length(vmp) > 1) {
      cat("processing a multivariate family, outcome:", m.target, "\n")
    }
    lapply(1:p.vmp, function(k) {
      oo <- extract.bootsample(x, alpha = alpha, target = k - 1,
                   m.target = m.target, standardize = standardize)
      cat("===== VIMP confidence regions for", vmp.col.names[k], " =====\n")
      cat("nonparametric:\n")
      print(round(oo$ci, 3))
      cat("parametric:\n")
      print(round(oo$ci.Z, 3))
      NULL
    })
  })
}
print.bootsample <- print.bootsample.rfsrc
