get.varselect.imp <- function(f.o, target.dim) {
  if (!is.null(f.o$importance)) {
    c(cbind(f.o$importance)[, target.dim])
  }
    else {
      rep(NA, length(f.o$xvar.names))
    }
}
get.varselect.imp.all <- function(f.o) {
  if (!is.null(f.o$importance)) {
    imp.all <- cbind(f.o$importance)
    if (ncol(imp.all) == 1) {
      colnames(imp.all) <- "vimp"
    }
      else {
        colnames(imp.all) <- paste("vimp.", colnames(imp.all), sep = "")
      }
    imp.all
  }
    else {
      rep(NA, length(f.o$xvar.names))
    }
}
get.varselect.err <- function(f.o) {
  if (!is.null(f.o$err.rate)) {
    if (grepl("surv", f.o$family)) {
      err <- 100 * cbind(f.o$err.rate)[f.o$ntree, ]
    }
      else {
        err <- cbind(f.o$err.rate)[f.o$ntree, ]
      }
  }
    else {
      err = NA
    }
  err
}
get.varselect.length <- function(x, y) {
  (length(x) > 0 & length(y) > 0)
}
get.varselect.mtry <- function(x, y) {
  mtry <- round((length(x) - length(y))/3)
  if (mtry == 0) {
    round(length(x)/3)
  }
    else {
      mtry
    }
}
get.varselect.sd <- function(x) {
  if (all(is.na(x))) {
    NA
  }
    else {
      sd(x, na.rm = TRUE)
    }
}
permute.rows <-function(x) {
  n <- nrow(x)
  p <- ncol(x)
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}
balanced.folds <- function(y, nfolds = min(min(table(y)), 10)) {
  y[is.na(y)] <- resample(y[!is.na(y)], size = sum(is.na(y)), replace = TRUE)
  totals <- table(y)
  if (length(totals) < 2) {
    return(cv.folds(length(y), nfolds))
  }
    else {
      fmax <- max(totals)
      nfolds <- min(nfolds, fmax)     
      nfolds <- max(nfolds, 2)
      folds <- as.list(seq(nfolds))
      yids <- split(seq(y), y) 
      bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
      for(i in seq(totals)) {
        if(length(yids[[i]])>1){bigmat[seq(totals[i]), i] <- sample(yids[[i]])}
        if(length(yids[[i]])==1){bigmat[seq(totals[i]), i] <- yids[[i]]}
      }
      smallmat <- matrix(bigmat, nrow = nfolds)
      smallmat <- permute.rows(t(smallmat)) 
      res <- vector("list", nfolds)
      for(j in 1:nfolds) {
        jj <- !is.na(smallmat[, j])
        res[[j]] <- smallmat[jj, j]
      }
      return(res)
    }
}
