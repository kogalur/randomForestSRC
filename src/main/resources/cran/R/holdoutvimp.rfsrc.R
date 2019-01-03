holdoutvimp.rfsrc <- function(formula, data,
                              ntree = 1000 * ncol(data) / vtry,
                              ntree.max = 2000,
                              nsplit = 10,
                              ntime = 50,
                              mtry = NULL,
                              vtry = 1,
                              fast = FALSE,
                              verbose = TRUE,
                              ...)
{
  ## pull unnamed parameters to pass to rfsrc
  dots <- list(...)
  dots$formula <- dots$data <- dots$ntree <- dots$nsplit <-
    dots$ntime <- dots$mtry <- dots$importance <- dots$block.size <- NULL
  ## add vtry to the list
  dots$vtry <- vtry
  ## make block size large to reduce computations
  dots$block.size <- ntree
  ## determine the grow interface - rfsrc or rfsrcFast?
  if (!fast) {
    rfsrc.grow <- "rfsrc"
  }
  else {
    rfsrc.grow <- "rfsrcFast"
  }
  ## grow a random forest
  o <- do.call(rfsrc.grow, c(list(formula = formula, data = data, ntree = ntree,
                        nsplit = nsplit, mtry = mtry), dots))
  ## extract the holdout array
  ho <- o$holdout.array
  ## pull the yvar names 
  ynms <- o$yvar.names
  if (o$family == "surv" || o$family == "surv-CR") {
    ynms <- ynms[1]
  }
  ## verbose setup
  if (verbose) pb <- txtProgressBar(min = 0, max = length(o$xvar.names), style = 3)
  ## loop over features
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
    ## calculate PE for in and out v
    if (sum(pt.out) > 0 && sum(!pt.out) > 0) {
      v.out <- which(pt.out)
      v.in  <- which(!pt.out)
      ## restrict size of v.in and v.out using ntree.max
      v.out <- sample(v.out)[1:min(ntree.max, length(v.out))]
      v.in  <- sample(v.in)[1:min(ntree.max, length(v.in))]
      ## extract the ensembles for hold out and hold in trees
      ## block size set for safety but generic does this for us
      pe.out <- predict(o, get.tree = v.out, block.size = length(v.out))
      pe.in  <- predict(o, get.tree = v.in,  block.size = length(v.in))
      ## loop over outcomes - in case this is a multivariate object
      do.call(cbind, lapply(ynms, function(nn) {
        ## pull the hold out error rate
        o.out <- coerce.multivariate(pe.out, nn)
        err.out <- get.holdout.err(o.out, length(v.out))
        ## pull the hold in error rate
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
holdoutvimp <- holdoutvimp.rfsrc
###################################################################
##
## extract the holdout error rate
##
###################################################################
get.holdout.err <- function(obj, j) {
  err <- obj$err.rate
  if (obj$family == "class") {
    cbind(err)[j, 1]
  }
  else if (obj$family == "surv-CR") {
    err[j, ]
  }
  else {
    err[j]
  }
}
