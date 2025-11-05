impute.rfsrc <- function(formula, data,
                         ntree = 100, nodesize = 1, nsplit = 10,
                         nimpute = 2, fast = FALSE, blocks,                         
                         mf.q, max.iter = 10, eps = 0.01,
                         ytry = NULL, always.use = NULL, verbose = TRUE,
                         full.sweep = FALSE,      ## NEW: optional final sweep
                         ...)
{
  ##--------------------------------------------------------------
  ##
  ## preliminary checks to make sure impute makes sense
  ##
  ##--------------------------------------------------------------
  ## terminate if there is no data
  if (missing(data)) {
    stop("data is missing")
  }
  ## identify the missing data
  ## if none/all: return the data
  which.na <- is.na(data)
  if (!any(which.na) || all(which.na)) {
    return(invisible(data))
  }
  ##--------------------------------------------------------------
  ##
  ## extract additional options specified by user
  ## we lock this down - here are the options we allow
  ##
  ##--------------------------------------------------------------
  ## list of allowed parameters
  rfnames <- c("mtry",
               "splitrule",
               "bootstrap",
               "sampsize",
               "samptype")
  ## get all user-specified dots
  dots.all <- list(...)
  ## get the permissible hidden options for main impute
  dots <- dots.all[names(dots.all) %in% rfnames]
  ## full sweep options (own list under 'full.sweep.options')
  fs.opts <- dots.all[["full.sweep.options"]]
  if (is.null(fs.opts)) fs.opts <- list()
  fs.dots <- fs.opts[names(fs.opts) %in% rfnames]
  fs.ntree    <- if (!is.null(fs.opts$ntree))    fs.opts$ntree    else 500
  fs.nodesize <- if (!is.null(fs.opts$nodesize)) fs.opts$nodesize else NULL
  fs.nsplit   <- if (!is.null(fs.opts$nsplit))   fs.opts$nsplit   else 10
  ##--------------------------------------------------------------
  ##
  ## process the data
  ##
  ##--------------------------------------------------------------
  p <- ncol(data)
  n <- nrow(data)
  all.r.na <- rowSums(which.na) == p
  all.c.na <- colSums(which.na) == n
  data <- data[!all.r.na, !all.c.na, drop = FALSE]
  which.na <- which.na[!all.r.na, !all.c.na, drop = FALSE]
  if (!any(which.na)) {
    return(data)
  }
  p <- ncol(data)
  n <- nrow(data)
  all.var.names <- colnames(data)
  ## mforest details
  if (missing(mf.q)) {
    mforest <- FALSE
  }
  else {
    mforest <- TRUE
  }
  ## set the number of blocks used to subdivide the data
  if (!missing(blocks)) {
    blocks <- cv.folds(nrow(data), max(1, blocks))
  }
  else {
    blocks <- list(1:nrow(data))
  }
  ##--------------------------------------------------------------
  ## METHOD 1: default impute call
  ##--------------------------------------------------------------
  if (!mforest) {
    if (missing(formula)) {
      if (is.null(ytry)) {
        ytry <- min(p - 1, max(25, ceiling(sqrt(p))))
      }
      dots$formula <- as.formula(paste("Unsupervised(", ytry, ") ~ ."))
      dots$splitrule <- NULL
    }
    else {
      dots$formula <- formula
    }
    nullBlocks <- lapply(blocks, function(blk) {
      dta <- data[blk,, drop = FALSE]
      retO <- tryCatch({do.call("generic.impute.rfsrc",
                      c(list(data = dta,
                             ntree = ntree,
                             nodesize = nodesize,
                             nsplit = nsplit,
                             nimpute = nimpute,
                             fast = fast), dots))}, error = function(e) {NULL})
      if (!is.null(retO)) {
        if (!is.null(retO$missing$row)) {
          blk <- blk[-retO$missing$row]
        }
        if (!is.null(retO$missing$col)) {
          ynames <- all.var.names[-retO$missing$col]
        }
          else {
            ynames <- all.var.names
          }
        data[blk, ynames] <<- retO$data[, ynames, drop = FALSE]
      }
      NULL
    })
    rm(nullBlocks)
  }
  ##--------------------------------------------------------------
  ## METHOD 2: mforest
  ##--------------------------------------------------------------
  if (mforest) {
    x.na <- lapply(1:p, function(k) {
      if (sum(which.na[, k]) > 0) {
        as.numeric(which(which.na[, k]))
      }
        else {
          NULL
        }
    })
    which.x.na <- which(sapply(x.na, length) > 0)
    names(x.na) <- all.var.names <- colnames(data)
    var.names <- all.var.names[which.x.na] 
    if (!is.null(always.use)) {
      always.use <- is.element(all.var.names, always.use)
      if (sum(always.use) > 0) {
        always.use <- which(always.use)
      }
    }
    p0 <- length(which.x.na)
    mfOriginal <- FALSE
    if (mf.q == 0) {
      stop("mf.q must be greater than zero")
    }
    if (mf.q == 1 && is.null(always.use)) {
      mfOriginal <- TRUE
    }
    ## Change D: guard when only one variable has missingness
    if (p0 == 1) {
      K <- 1
    } else {
      if (mf.q >= 1) {
        mf.q <- min(p0 - 1, mf.q) / p0
      }
      K <- max(1, round(max(1 / mf.q, 2)))
    }
    dots.rough <- dots
    dots.rough$mtry <- NULL
    dots.rough$splitrule <- "random"
    dots.rough$formula <- as.formula(zZ999Zz ~.)
    data <- do.call("generic.impute.rfsrc",
                    c(list(data = data.frame(zZ999Zz = rnorm(nrow(data)), data),
                           ntree = 10,
                           nodesize = nodesize,
                           nsplit = nsplit,
                           nimpute = 1,
                           fast = fast), dots.rough))$data
    data$zZ999Zz <- NULL
    diff.err <- Inf
    check <- TRUE
    var.grp <- cv.folds(p0, K)
    K <- length(var.grp)
    if (verbose) {
      if (mfOriginal) {
        cat("missForest parameters:", paste0("(#max.iter, #vars)=(", max.iter, ",", p0, ")"), "\n")
      }
      else {
        cat("multivariate missForest parameters:",
            paste0("(#iter, #vars, #blks)=(", max.iter, ",", p0, ",", K, ")"), "\n")
      }
    }
    nullWhile <- lapply(1:max.iter, function(m) {
      if (!check) {
        return(NULL)
      }
      var.grp <- cv.folds(p0, K)
      if (verbose) {
        if (max.iter > 1) {
          cat("\t", paste0("iteration:", m, "\n"))
        }
      }
      data.old <- data
      nullBlocks <- lapply(blocks, function(blk) {
        nullObj <- lapply(var.grp, function(grp) {
          if (verbose) {
            cat(".")
          }
          if (!mfOriginal) {
            ynames <- unique(c(var.names[grp], all.var.names[always.use]))
            dots$formula <- as.formula(paste("Multivar(", paste(ynames, collapse = ","), paste(") ~ ."), sep = ""))
            dta <- data[blk,, drop = FALSE]
            dta[, ynames] <- lapply(ynames, function(nn) {
              xk <- data[, nn]
              xk[unlist(x.na[nn])] <- NA
              xk[blk]
            })
            mvimpute <- tryCatch({do.call("generic.impute.rfsrc",
                                      c(list(data = dta,
                                             ntree = ntree,
                                             nodesize = nodesize,
                                             nsplit = nsplit,
                                             nimpute = 1,
                                             fast = fast), dots))}, error = function(e) {NULL})
            if (!is.null(mvimpute)) {
              if (!is.null(mvimpute$missing$row)) {
                blk <- blk[-mvimpute$missing$row]
              }
              if (!is.null(mvimpute$missing$col)) {
                ynames <- ynames[-mvimpute$missing$col]
              }
              data[blk, ynames] <<- mvimpute$data[, ynames, drop = FALSE]
              rm(dta)
            }
          }
          if (mfOriginal) {
            yname <- var.names[grp]
            trn <- setdiff(blk,  unlist(x.na[yname]))
            tst <- setdiff(blk,  trn)
            if (length(trn) > 0 && length(tst) > 0) {
              dots$formula <- as.formula(paste0(yname, "~."))
              grow <- tryCatch({do.call("rfsrc",
                                        c(list(data = data[trn,, drop = FALSE],
                                               ntree = ntree,
                                               nodesize = nodesize,
                                               nsplit = nsplit,
                                               perf.type = "none",
                                               fast = fast), dots))}, error = function(e) {NULL})
              if (!is.null(grow)) {
                pred <- predict(grow, data[tst,, drop = FALSE])
                if (grow$family == "regr") {
                  data[tst, yname] <<- pred$predicted
                }
                else {
                  data[tst, yname] <<- pred$class
                }
              }
            }
          }
          NULL
        })
        NULL
      })
      diff.new.err <- mean(sapply(var.names, function(nn) {
        xo <- data.old[unlist(x.na[nn]), nn]
        xn <- data[unlist(x.na[nn]), nn]
        if (!is.numeric(xo)) {
          sum(xn != xo, na.rm = TRUE) / (.001 + length(xn))
        }
          else {
            var.xo <- var(xo, na.rm = TRUE)
            if (is.na(var.xo)) {
              var.xo <- 0
            }
            sqrt(mean((xn - xo)^2, na.rm = TRUE) / (.001 + var.xo))
          }
      }), na.rm = TRUE)
      if (verbose) {
        cat("\n")
        err <- paste("err = " , format(diff.new.err, digits = 3), sep = "")
        drp <- paste("drop = ", format(diff.err - diff.new.err, digits = 3), sep = "")
        cat("         >> ", err, ", ", drp, "\n")
      }
      check <<- ((diff.err - diff.new.err) >= eps)
      diff.err <<- diff.new.err
      rm(data.old)
      NULL
    })
  }
  ##--------------------------------------------------------------
  ##
  ## optional final sweep over originally-missing variables
  ## (applies to both methods) + progress output
  ##
  ##--------------------------------------------------------------
  if (isTRUE(full.sweep)) {
    ## identify row indices of original missing values for each variable
    x.na.sweep <- lapply(1:ncol(which.na), function(k) {
      idx <- which(which.na[, k])
      if (length(idx) > 0) as.numeric(idx) else NULL
    })
    names(x.na.sweep) <- colnames(which.na)
    ## variables with any original missingness
    sweep.vars <- names(which(sapply(x.na.sweep, length) > 0))
    if (length(sweep.vars) > 0) {
      nvars <- length(sweep.vars)
      nullSweep <- lapply(seq_along(sweep.vars), function(j) {
        yname <- sweep.vars[[j]]
        ## progress output
        if (verbose) {
          pct <- round(100 * j / nvars)
          cat("--> full.sweep:", paste0("[", j, "/", nvars, "]"), yname,
              paste0(" (", pct, "%)"), "\n")
        }
        tst <- x.na.sweep[[yname]]
        if (length(tst) == 0) return(NULL)
        trn <- setdiff(seq_len(nrow(data)), tst)
        if (length(trn) == 0) return(NULL)
        ## fit on observed y rows using final data
        grow <- tryCatch({
          do.call("rfsrc",
                  c(list(formula   = as.formula(paste0(yname, "~.")),
                         data      = data[trn,, drop = FALSE],
                         ntree     = fs.ntree,
                         nodesize  = fs.nodesize,
                         nsplit    = fs.nsplit,
                         perf.type = "none",
                         fast      = fast), fs.dots))
        }, error = function(e) { NULL })
        if (!is.null(grow)) {
          pred <- tryCatch(predict(grow, data[tst,, drop = FALSE]),
                           error = function(e) { NULL })
          if (!is.null(pred)) {
            if (grow$family == "regr") {
              data[tst, yname] <<- pred$predicted
            }
            else {
              data[tst, yname] <<- pred$class
            }
          }
        }
        NULL
      })
      rm(nullSweep)
      if (verbose) {
        cat("full.sweep: completed (", length(sweep.vars), " variables)\n", sep = "")
      }
    }
  }
  ##--------------------------------------------------------------
  ##
  ## return the imputed data
  ##
  ##--------------------------------------------------------------
  invisible(data)
}
impute <- impute.rfsrc
