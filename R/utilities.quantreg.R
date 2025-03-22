## pulls the quantile object and converts to a list
extract.quantile <- function(o) {
  ## confirm this is a quantile regression object
  if (sum(grepl("quantreg", class(o))) == 0) { 
    stop("object must be a quantreg object")
  }
  ## used to determine univariate vs multivariate objects
  mv.y.names <- intersect(names(o$quantreg), o$yvar.names)
  ## we have a univariate quantile regression object
  if (length(mv.y.names) == 0) {
    q <- list(o$quantreg)
    names(q) <- o$yvar.names
  }
  ## the quantile regression object is multivariate
  else {
    q <- lapply(mv.y.names, function(m.target) {
      o$quantreg[[m.target]]
    })
    names(q) <- mv.y.names
  }
  ## return the processed list
  q
}
## extract target quantiles
get.quantile <- function(o, target.prob = NULL, pretty = TRUE) {
  ## extract the quantile object
  qo <- extract.quantile(o)
  ## process the target probabilities
  if (!is.null(target.prob)) {
    target.prob <- sort(unique(target.prob))
  }
  else {## default is to use existing values
    target.prob <- qo[[1]]$prob
  }
  ## pull the target quantiles
  rO <- lapply(qo, function(q) {    
    q.dat <- do.call(cbind, lapply(target.prob, function(pr) {
      q$quant[, which.min(abs(pr - q$prob))]
    }))
    colnames(q.dat) <-  paste("q.", 100 * target.prob, sep = "")
    q.dat
  })
  if (pretty && length(rO) == 1) {
    rO <- rO[[1]]
  }
  rO
}
## extract crps
get.quantile.crps <- function(o, pretty = TRUE, subset = NULL, standardize = TRUE) {
  ## extract the quantile object
  qO <- extract.quantile(o)
  ## does not apply to predict objects without y
  if (sum(grepl("predict", class(o))) > 0 && is.null(o$yvar)) {
    stop("no yvar present in quantreg predict object")
  }
  ## subset assignment
  if (is.null(subset)) {
    subset <- 1:o$n
  }
  else {
    if (is.logical(subset)) {
      subset <- which(subset)
    }
    subset <- subset[subset >=1 & subset <= o$n]
  }
  if (length(subset) == 0) {
    stop("requested subset is empty")
  }
  ## pull the target stats
  rO <- lapply(1:length(qO), function(j) {
    q <- qO[[j]]
    if (is.vector(o$yvar)) {
      y <- cbind(o$yvar)
    }
    else {
      y <- o$yvar[, names(qO)[j]]
    }
    n <- length(y)
    ## brier score
    brS <- colMeans(do.call(rbind, mclapply(subset, function(i) {
      (1 * (y[i] <= q$yunq) - q$cdf[i, ]) ^ 2 
    })), na.rm = TRUE)
    ## crps
    crps <- unlist(lapply(1:length(q$yunq), function(j) {
      if (standardize) {
        trapz(q$yunq[1:j], brS[1:j]) / diff(range(q$yunq[1:j]))
      }
      else {
        trapz(q$yunq[1:j], brS[1:j])
      }
    }))
    data.frame(y = q$yunq, crps = crps)
  })
  names(rO) <- names(qO)
  if (pretty && length(rO) == 1) {
    rO <- rO[[1]]
  }
  rO
}
## extract target stats
get.quantile.stat <- function(o, pretty = TRUE) {
  ## extract the quantile object
  qO <- extract.quantile(o)
  ## conditional median
  mdn <- get.quantile(o, .5, FALSE)
  ## pull the target stats
  rO <- lapply(1:length(qO), function(j) {
    q <- qO[[j]]
    ## conditional mean
    mn <- q$density %*% q$yunq
    ## conditional standard deviation
    std <- sqrt(q$density %*% q$yunq^2 - mn ^ 2)
    data.frame(mean = mn, median = c(mdn[[j]]), std = std)
  })
  names(rO) <- names(qO)
  if (pretty && length(rO) == 1) {
    rO <- rO[[1]]
  }
  rO
}
