###################################################################
##
## core bootstrap function
##
###################################################################
bootsample <- function(obj, rf.prms, B, block.size, joint, xvar.names, importance, performance, verbose) {
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
    obj <- vimp(obj, perf.type = rf.prms$perf.type, block.size = block.size, importance = importance)
    vmp <- get.mv.vimp(obj, FALSE, FALSE)
    rf.prms$block.size <- block.size
    if (verbose) cat("done\n")
  }
  else {
    importance <- obj$forest$importance
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
    ## calculate VIMP from the double boostrap forest
    rf.b <- do.call("rfsrc",
                    c(list(data = dta[smp.unq,, drop = FALSE], importance = importance,
                           bootstrap = "by.user", samp = dbl.smp), rf.prms))
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
