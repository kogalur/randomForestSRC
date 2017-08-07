quantileReg.rfsrc <- function(obj, oob = TRUE, prob = (1:10) / 10, newdata = NULL) {
  if (!(obj$family == "regr" | !is.null(obj$regrOutput))) {
    stop("this function only applies to regression settings\n")
  }
  if (obj$family == "regr") {
    ynames <- obj$yvar.names
  }
  else {
    ynames <- names(obj$regrOutput)
  }
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
    else {
      fwt <- predict(obj, newdata, forest.wt = TRUE)$forest.wt
    }
  } 
  rO <- lapply(ynames, function(yn) {
    if (ncol(cbind(obj$yvar)) > 1) {
      y <- obj$yvar[, yn]
    }
    else {
      y <- obj$yvar
    }
    yunq <- sort(unique(y))
    ind.matx <- do.call(rbind, mclapply(yunq, function(yy) {y <= yy}))
    cdf <- t(apply(fwt, 1, function(wt) {ind.matx %*% wt}))
    sIndex <- function(x, y) {sapply(1:length(y), function(j) {sum(x <= y[j])})}
    quant <- t(apply(cdf, 1, function(pr) {
      c(min(yunq, na.rm = TRUE), yunq)[1 + sIndex(pr, prob)]
    }))
    list(quantiles = quant,
         prob = prob,
         density = t(apply(cbind(0,cdf), 1, diff)),
         yunq = yunq)
  })
  if (obj$family == "regr") {
    rO[[1]]
  }
  else {
    names(rO) <- ynames
    rO
  }
}
quantileReg <- quantileReg.rfsrc
