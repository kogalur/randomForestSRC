quantileReg.rfsrc <- function(obj, oob = TRUE, prob = (1:10) / 10, newdata = NULL) {
  ## forest object must contain regression outcomes
  ## this includes multivariate regression and mixed multivariate regression
  if (!(obj$family == "regr" | !is.null(obj$regrOutput))) {
    stop("this function only applies to regression settings\n")
  }
  ## pull the target outcome names
  if (obj$family == "regr") {
    ynames <- obj$yvar.names
  }
  else {
    ynames <- names(obj$regrOutput)
  }
  ## training:
  ## forest weights
  ## test whether forest weights are available
  ## cost saving measure for advanced users
  ## otherwise pull the forest weights by calling predict
  ##
  ## testing:
  ## pull forest weights for test data set
    if (!is.null(obj$forest.wt) && is.null(newdata)) {
    fwt <- obj$forest.wt
  }
  else {
    if (is.null(newdata)) {
      if (oob) {
        fwt <- predict(obj, forest.wt = "oob")$forest.wt
      }
      else {
        fwt <- predict(obj, forest.wt = TRUE)$forest.wt
      }
    }
    else {##test data is available
      fwt <- predict(obj, newdata, forest.wt = TRUE)$forest.wt
    }
  } 
  ## calculate the quantiles
  rO <- lapply(ynames, function(yn) {
    ## extract y
    if (ncol(cbind(obj$yvar)) > 1) {
      y <- obj$yvar[, yn]
    }
    else {
      y <- obj$yvar
    }
    ## cdf calculations
    yunq <- sort(unique(y))
    ind.matx <- do.call(rbind, mclapply(yunq, function(yy) {y <= yy}))
    cdf <- t(apply(fwt, 1, function(wt) {ind.matx %*% wt}))
    ## step function interpolation
    sIndex <- function(x, y) {sapply(1:length(y), function(j) {sum(x <= y[j])})}
    ## quantiles
    quant <- t(apply(cdf, 1, function(pr) {
      c(min(yunq, na.rm = TRUE), yunq)[1 + sIndex(pr, prob)]
    }))
    ## return the object
    list(quantiles = quant,
         cdf = cdf,
         prob = prob,
         density = t(apply(cbind(0,cdf), 1, diff)),
         yunq = yunq)
  })
  ## return the goodies
  if (obj$family == "regr") {
    rO[[1]]
  }
  else {
    names(rO) <- ynames
    rO
  }
}
quantileReg <- quantileReg.rfsrc
