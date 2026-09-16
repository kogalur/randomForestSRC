## Complete a training table using on-the-fly or missForest-style imputation.
##
## This function returns a table, not a reusable learner bank. impute.learn()
## calls it for the initial completed table and then grows its own saved forests.
## which.na is the retained table's original missingness mask; x.na records
## those same cells by variable in the missForest branch.
impute.rfsrc <- function(formula, data,
                         ntree = 100, nodesize = 1, nsplit = 10,
                         nimpute = 2, fast = FALSE, blocks,                         
                         mf.q, max.iter = 10, eps = 0.01,
                         ytry = NULL, always.use = NULL, verbose = TRUE,
                         full.sweep = FALSE,      ## optional final sweep
                         restore.integer = TRUE,
                         ...)
{
  ##--------------------------------------------------------------
  ##
  ## validate the input container and detect whether work is needed
  ##
  ##--------------------------------------------------------------
  ## terminate if there is no data
  if (missing(data)) {
    stop("data is missing")
  }
  data <- .impute.data.frame(data)
  .impute.scalar(max.iter, "max.iter", lower = 1, whole = TRUE,
                 upper = .Machine$integer.max)
  .impute.scalar(eps, "eps", lower = 0)
  if (!is.null(ytry)) .impute.scalar(ytry, "ytry", lower = 1, whole = TRUE,
                                    upper = .Machine$integer.max)
  if (!is.null(always.use) && (!is.character(always.use) || anyNA(always.use))) {
    stop("'always.use' must be NULL or a character vector of variable names.",
         call. = FALSE)
  }
  .impute.forest.options(list(ntree = ntree, nodesize = nodesize,
                              nsplit = nsplit, nimpute = nimpute))
  if (!missing(mf.q)) .impute.scalar(mf.q, "mf.q", lower = 0, open.lower = TRUE)
  if (!missing(blocks)) .impute.scalar(blocks, "blocks", lower = 1, whole = TRUE,
                                      upper = .Machine$integer.max)
  for (nm in c("fast", "verbose", "full.sweep", "restore.integer")) {
    .impute.flag(get(nm), nm)
  }
  ##--------------------------------------------------------------
  ##
  ## extract additional options specified by user
  ## only the allow-listed forest options are forwarded
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
  fs.opts <- .impute.named.options(fs.opts, "full.sweep.options")
  .impute.forest.options(fs.opts)
  fs.dots <- fs.opts[names(fs.opts) %in% rfnames]
  fs.ntree    <- if (!is.null(fs.opts$ntree))    fs.opts$ntree    else 500
  fs.nodesize <- if (!is.null(fs.opts$nodesize)) fs.opts$nodesize else NULL
  fs.nsplit   <- if (!is.null(fs.opts$nsplit))   fs.opts$nsplit   else 10
  ## identify the missing data
  ## if none/all: return the data
  which.na <- is.na(data)
  if (!any(which.na) || all(which.na)) {
    return(invisible(data))
  }
  ##--------------------------------------------------------------
  ##
  ## process the data
  ##
  ##--------------------------------------------------------------
  p <- ncol(data)
  n <- nrow(data)
  ## Drop entirely missing rows/columns together and subset the immutable
  ## missingness mask by the same positions.
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
  ## Learn integer support before filling. Restoration uses which.na to round
  ## generated cells only; observed numeric values retain their original precision.
  ## Integer storage is restored only when every retained value is representable.
  integer.support <- .integer.support.map(data)
  ## Use the same branch specification as impute.learn() metadata.
  method <- .impute.method.spec(
    mf.q = if (missing(mf.q)) NULL else mf.q,
    always.use = always.use,
    variable.names = all.var.names
  )
  mforest <- method$mforest
  ## set the number of blocks used to subdivide the data
  if (!missing(blocks)) {
    blocks <- cv.folds(nrow(data), max(1, blocks))
  }
  else {
    blocks <- list(seq_len(nrow(data)))
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
    ## Process row blocks independently. Successful block results are overlaid
    ## onto their retained row/column positions in the enclosing data table.
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
        if (length(retO$missing$row) > 0L) {
          blk <- blk[-retO$missing$row]
        }
        ## The generic helper can return columns in response/predictor order.
        ## Select retained input columns by name, not by output position.
        ynames <- intersect(all.var.names, colnames(retO$data))
        if (length(blk) > 0L && length(ynames) > 0L) {
          data[blk, ynames] <<- retO$data[, ynames, drop = FALSE]
        }
      }
      NULL
    })
    rm(nullBlocks)
  }
  ##--------------------------------------------------------------
  ## METHOD 2: mforest
  ##--------------------------------------------------------------
  if (mforest) {
    x.na <- lapply(seq_len(p), function(k) {
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
    always.use <- method$always.use
    p0 <- length(which.x.na)
    mfOriginal <- method$univariate
    ## A single missing target requires one nonempty response group.
    if (p0 == 1) {
      K <- 1
    } else {
      if (mf.q >= 1) {
        mf.q <- min(p0 - 1, mf.q) / p0
      }
      K <- min(p0, max(1, round(max(1 / mf.q, 2))))
    }
    ## Initialize all missing cells using a random-split forest with a
    ## temporary numeric response before response-wise iterative updates.
    dots.rough <- dots
    dots.rough$mtry <- NULL
    dots.rough$splitrule <- "random"
    response.name <- .impute.response.name(names(data))
    dots.rough$formula <- .impute.formula(response.name)
    rough.data <- data
    rough.data[[response.name]] <- rnorm(nrow(data))
    rough.data <- rough.data[, c(response.name, names(data)), drop = FALSE]
    data <- do.call("generic.impute.rfsrc",
                    c(list(data = rough.data,
                           ntree = 10,
                           nodesize = nodesize,
                           nsplit = nsplit,
                           nimpute = 1,
                           fast = fast), dots.rough))$data
    data[[response.name]] <- NULL
    rm(rough.data)
    data <- .restore.integer.data(data, integer.support,
                                  restore.integer = restore.integer,
                                  generated = which.na)
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
    ## Regroup response variables at every pass. Within a pass, later groups
    ## use the current completed values written by earlier groups.
    nullWhile <- lapply(seq_len(max.iter), function(m) {
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
            lhs <- as.call(c(list(as.name("Multivar")), lapply(ynames, as.name)))
            dots$formula <- .impute.formula(lhs)
            dta <- data[blk,, drop = FALSE]
            ## Reset response cells to their original missingness for this
            ## group; keep current imputations available on the predictor side.
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
              if (length(mvimpute$missing$row) > 0L) {
                blk <- blk[-mvimpute$missing$row]
              }
              ## Deletion indices refer to the full input table, not ynames.
              ## Retain only requested responses actually present in the result.
              ynames <- intersect(ynames, colnames(mvimpute$data))
              if (length(blk) > 0L && length(ynames) > 0L) {
                data[blk, ynames] <<- mvimpute$data[, ynames, drop = FALSE]
              }
              rm(dta)
            }
          }
          if (mfOriginal) {
            yname <- var.names[grp]
            ## Univariate missForest learns from originally observed target
            ## rows and predicts only its originally missing rows in this block.
            trn <- setdiff(blk,  unlist(x.na[yname]))
            tst <- setdiff(blk,  trn)
            if (length(trn) > 0 && length(tst) > 0) {
              dots$formula <- .impute.formula(yname)
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
      ## Summarize changes only at originally missing entries. Numerics use
      ## the previous imputed values' variance; factors use disagreement rates.
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
      ## Continue while the change statistic decreases by at least eps.
      ## The current pass is retained, including the first non-improving pass.
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
    data <- .restore.integer.data(data, integer.support,
                                  restore.integer = restore.integer,
                                  generated = which.na)
    ## identify row indices of original missing values for each variable
    x.na.sweep <- lapply(seq_len(ncol(which.na)), function(k) {
      idx <- which(which.na[, k])
      if (length(idx) > 0) as.numeric(idx) else NULL
    })
    names(x.na.sweep) <- colnames(which.na)
    ## variables with any original missingness
    sweep.vars <- names(which(sapply(x.na.sweep, length) > 0))
    if (length(sweep.vars) > 0) {
      nvars <- length(sweep.vars)
      ## Unlike the saved bank in impute.learn(), this optional sweep updates
      ## the table immediately, so later targets see earlier sweep predictions.
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
                  c(list(formula   = .impute.formula(yname),
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
    data <- .restore.integer.data(data, integer.support,
                                  restore.integer = restore.integer,
                                  generated = which.na)
  }
  ##--------------------------------------------------------------
  ##
  ## return the imputed data
  ##
  ##--------------------------------------------------------------
  data <- .restore.integer.data(data, integer.support,
                                restore.integer = restore.integer,
                                generated = which.na)
  invisible(data)
}
impute <- impute.rfsrc
