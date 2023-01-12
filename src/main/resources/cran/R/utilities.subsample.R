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
