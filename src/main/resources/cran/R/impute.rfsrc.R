impute.rfsrc <- function(formula,
                         data,
                         ntree = 500,
                         mtry = NULL,
                         xvar.wt = NULL,
                         nodesize = 1,
                         splitrule = NULL,
                         nsplit = 10,
                         
                         na.action = c("na.impute"),
                         
                          
                         nimpute = 2,
                         mf.q, blocks,
                         always.use = NULL, 
                         max.iter = 10,
                         eps = 0.01,
                         verbose = TRUE, 
                         do.trace = FALSE,
                         ...)
{
  ## get any hidden options
  user.option <- list(...)
  ytry <- is.hidden.ytry(user.option)
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
  ## acquire various dimensions/information
  p <- ncol(data)
  n <- nrow(data)
  ## we must dispose of all rows, columns with all missigness
  ## check after processing there still exists missing data
  all.r.na <- rowSums(which.na) == p
  all.c.na <- colSums(which.na) == n
  data <- data[!all.r.na, !all.c.na, drop = FALSE]
  which.na <- which.na[!all.r.na, !all.c.na, drop = FALSE]
  if (!any(which.na)) {
    return(data)
  }
  p <- ncol(data)
  n <- nrow(data)
  ## extract the variable names to retain consistency of their order
  all.var.names <- colnames(data)
  ## mforest details
  if (missing(mf.q)) {
    mforest <- FALSE
  }
    else {
      mforest <- TRUE
    }
  ## set the number of blocks used to subdivide the data
  ## hence subdivide the problem into more manageable pieces
  if (!missing(blocks)) {
    blocks <- cv.folds(nrow(data), max(1, blocks))
  }
    else {
      blocks <- list(1:nrow(data))
    }
  ## generic impute call
  if (!mforest) {
    ## if the formula is missing we are implementing unsupervised forests
    if (missing(formula)) {
      if (is.null(ytry)) {
        ytry <- min(p - 1, max(25, ceiling(sqrt(p))))
      }
      formula <- as.formula(paste("Unsupervised(", ytry, ") ~ ."))
    }
    ## loop over data blocks
    nullBlocks <- lapply(blocks, function(blk) {
      ## assign the blocked data
      dta <- data[blk,, drop = FALSE]
      ## make the generic impute call
      retO <- tryCatch({generic.impute.rfsrc(formula = formula,
                                             data = dta,
                                             ntree = ntree,
                                             nimpute = nimpute,
                                             mtry = mtry,
                                             nodesize = nodesize,
                                             splitrule = splitrule,
                                             nsplit = nsplit,
                                             na.action = na.action,
                                             xvar.wt = xvar.wt,
                                             do.trace = do.trace)}, error = function(e) {NULL})
      ## if the impute object is non-null proceed ahead
      if (!is.null(retO)) {
        ## overlay the imputed data
        ## we retain the same column order
        if (!is.null(retO$missing$row)) {
          blk <- blk[-retO$missing$row]
        }
        if (!is.null(retO$missing$col)) {
          ynames <- all.var.names[-retO$missing$col]
        }
          else {
            ynames <- all.var.names
          }
        ## data reassignment
        data[blk, ynames] <<- retO$data[, ynames, drop = FALSE]
      }
      ## memory management
      NULL
    })
    rm(nullBlocks)
  }
  ## mforest call
  if (mforest) {
    ## identify which variables have missing data
    ## store the information as a convenient list
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
    ## always.use? if so - convert to index
    if (!is.null(always.use)) {
      always.use <- is.element(all.var.names, always.use)
      if (sum(always.use) > 0) {
        always.use <- which(always.use)
      }
    }
    ## set the multivariate response dimension
    ## convert mf.q to a fraction
    p0 <- length(which.x.na)
    if (mf.q == 0) {
      stop("mf.q must be greater than zero")
    }
    if (mf.q >= 1) {
      mf.q <- min(p0 - 1, mf.q) / p0
    }
    ## convert mf.q to K-fold selection
    K <- max(1 / mf.q, 2)
    ## quick and rough impute
    ## uncomment the following line for a better initial estimate
    data <- generic.impute.rfsrc(data = data,
                                 nimpute = 3,
                                 ntree = 250,
                                 mtry = mtry,
                                 nodesize = nodesize,
                                 nsplit = nsplit)$data
    ###############################################################
    ## main loop: data blocks/groups of variables
    ## we use lapply to avoid for-looping
    ###############################################################
    ## set flags
    diff.err <- Inf
    check <- TRUE
    nullWhile <- lapply(1:max.iter, function(m) {
      ## break
      if (!check) {
        return(NULL)
      }
      ## verbose
      if (verbose && max.iter > 1) {
        cat("\t iteration", m, "\n")
      }
      ## save the current state of data for assessing convergence
      data.old <- data
      ##--------------------------------------------------------
      ## loop over data blocks
      ##--------------------------------------------------------
      nullBlocks <- lapply(blocks, function(blk) {
        ## determine the grouping of the multivariate response
        var.grp <- cv.folds(p0, K)
        ##--------------------------------------------------------
        ## loop over the multivariate response groupings
        ##--------------------------------------------------------
        nullObj <- lapply(var.grp, function(grp) {
          ## multivariate formula
          ynames <- unique(c(var.names[grp], all.var.names[always.use]))
          f <- as.formula(paste("Multivar(", paste(ynames, collapse = ","), paste(") ~ ."), sep = ""))
          ## reset the chosen y missing data back to NA
          dta <- data[blk,, drop = FALSE]
          dta[, ynames] <- lapply(ynames, function(nn) {
            xk <- data[, nn]
            xk[unlist(x.na[nn])] <- NA
            xk[blk]
          })
          ## multivariate missForest call
          ## first column of returned object contains the imputed "y-values" which is what we want
          ## we do not use the other imputed data: this is very important
          retO <- tryCatch({generic.impute.rfsrc(f,
                                                 dta,
                                                 ntree = ntree,
                                                 nimpute = 1,
                                                 na.action =  na.action,
                                                 mtry = mtry,
                                                 nodesize = nodesize,
                                                 splitrule = splitrule,
                                                 nsplit = nsplit)}, error = function(e) {NULL})
          ## confirm that the impute object is non-null in order to proceed
          if (!is.null(retO)) {
            ## overlay the imputed data
            ## we retain the same column order
            if (!is.null(retO$missing$row)) {
              blk <- blk[-retO$missing$row]
            }
            if (!is.null(retO$missing$col)) {
              ynames <- ynames[-retO$missing$col]
            }
            data[blk, ynames] <<- retO$data[, ynames, drop = FALSE]
            rm(dta)
          }
          ## return NULL (memory management)
          NULL
        })
        ## return NULL (memory management)
        NULL
      })
      ## check whether algorithm has converged
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
        cat("         >> ", diff.new.err, diff.err - diff.new.err, "\n")
      }
      check <<- ((diff.err - diff.new.err) >= eps)
      diff.err <<- diff.new.err
      rm(data.old)
      ## return NULL (memory management)
      NULL
    })
  }
  ## #############################################################
  ## return the imputed data
  ## #############################################################
  invisible(data)
}
impute <- impute.rfsrc
