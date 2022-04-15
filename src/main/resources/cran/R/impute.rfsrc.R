impute.rfsrc <- function(formula, data,
                         ntree = 100, nodesize = 1, nsplit = 10,
                         nimpute = 2, fast = FALSE, blocks,                         
                         mf.q, max.iter = 10, eps = 0.01,
                         ytry = NULL, always.use = NULL, verbose = TRUE, 
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
  ## get the permissible hidden options
  dots <- list(...)
  dots <- dots[names(dots) %in% rfnames]
  ##--------------------------------------------------------------
  ##
  ## process the data
  ##
  ##--------------------------------------------------------------
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
  ##--------------------------------------------------------------
  ##
  ## set values needed for the impute call
  ##
  ##--------------------------------------------------------------
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
  ##--------------------------------------------------------------
  ##
  ## METHOD 1: default impute call
  ##
  ##--------------------------------------------------------------
  if (!mforest) {
    ## if the formula is missing we are implementing unsupervised forests
    if (missing(formula)) {
      if (is.null(ytry)) {
        ytry <- min(p - 1, max(25, ceiling(sqrt(p))))
      }
      dots$formula <- as.formula(paste("Unsupervised(", ytry, ") ~ ."))
      ## TBD2 splitrule not permissible with unsupervised forests?
      dots$splitrule <- NULL
    }
    else {
      dots$formula <- formula
    }
    ## loop over data blocks
    nullBlocks <- lapply(blocks, function(blk) {
      ## assign the blocked data
      dta <- data[blk,, drop = FALSE]
      ## make the generic impute call
      retO <- tryCatch({do.call("generic.impute.rfsrc",
                      c(list(data = dta,
                             ntree = ntree,
                             nodesize = nodesize,
                             nsplit = nsplit,
                             nimpute = nimpute,
                             fast = fast), dots))}, error = function(e) {NULL})
      ## proceed if the impute object is non-null
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
  ##--------------------------------------------------------------
  ##
  ## METHOD 2: mforest
  ##
  ##--------------------------------------------------------------
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
    ## will we do original mforest?
    p0 <- length(which.x.na)
    mfOriginal <- FALSE
    if (mf.q == 0) {
      stop("mf.q must be greater than zero")
    }
    if (mf.q == 1 && is.null(always.use)) {
      mfOriginal <- TRUE
    }
    if (mf.q >= 1) {
      mf.q <- min(p0 - 1, mf.q) / p0
    }
    ## convert mf.q to K-fold selection
    K <- max(1, round(max(1 / mf.q, 2)))
    ## quick and *rough* impute
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
    ###############################################################
    ##
    ## main loop: data blocks/groups of variables
    ## we use lapply to avoid for-looping
    ##
    ###############################################################
    ## set flags
    diff.err <- Inf
    check <- TRUE
    ## determine the grouping of the multivariate response
    var.grp <- cv.folds(p0, K)
    K <- length(var.grp)
    ## verbose
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
      ## break
      if (!check) {
        return(NULL)
      }
      ## determine the grouping of the multivariate response
      var.grp <- cv.folds(p0, K)
      ## verbose
      if (verbose) {
        if (max.iter > 1) {
          cat("\t", paste0("iteration:", m, "\n"))
        }
      }
      ## save the current state of data for assessing convergence
      data.old <- data
      ##--------------------------------------------------------
      ## loop over data blocks
      ##--------------------------------------------------------
      nullBlocks <- lapply(blocks, function(blk) {
        ##--------------------------------------------------------
        ## loop over the multivariate response groupings
        ##--------------------------------------------------------
        nullObj <- lapply(var.grp, function(grp) {
          ## missforest is slow, so provide verbose output for each regression/iteration
          if (verbose) {
            cat(".")
          }
          ## ------------------------------------------------------
          ##
          ## multivariate missForest call
          ## first column of returned object contains the imputed "y-values" which is what we want
          ## we do not use the other imputed data: this is very important
          ##
          ## -----------------------------------------------------
          if (!mfOriginal) {
            ## multivariate formula
            ynames <- unique(c(var.names[grp], all.var.names[always.use]))
            dots$formula <- as.formula(paste("Multivar(", paste(ynames, collapse = ","), paste(") ~ ."), sep = ""))
            ## reset the chosen y missing data back to NA
            dta <- data[blk,, drop = FALSE]
            dta[, ynames] <- lapply(ynames, function(nn) {
              xk <- data[, nn]
              xk[unlist(x.na[nn])] <- NA
              xk[blk]
            })
            ## multivariate generic impute call
            mvimpute <- tryCatch({do.call("generic.impute.rfsrc",
                                      c(list(data = dta,
                                             ntree = ntree,
                                             nodesize = nodesize,
                                             nsplit = nsplit,
                                             nimpute = 1,
                                             fast = fast), dots))}, error = function(e) {NULL})
            ## confirm that the impute object is non-null in order to proceed
            if (!is.null(mvimpute)) {
              ## overlay the imputed data
              ## we retain the same column order
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
          ## ------------------------------------------------------
          ##
          ## original mforest: grow (train) followed by prediction (test)
          ##
          ## -----------------------------------------------------
          if (mfOriginal) {
            ## make train, test data
            yname <- var.names[grp]
            trn <- setdiff(blk,  unlist(x.na[yname]))
            tst <- setdiff(blk,  trn)
            ## we must have train (no missing y) and test (missing y)
            if (length(trn) > 0 && length(tst) > 0) {
              dots$formula <- as.formula(paste0(yname, "~."))
              grow <- tryCatch({do.call("rfsrc",
                                        c(list(data = data[trn,, drop = FALSE],
                                               ntree = ntree,
                                               nodesize = nodesize,
                                               nsplit = nsplit,
                                               fast = fast), dots))}, error = function(e) {NULL})
              ## confirm grow object is non-null in order to proceed
              if (!is.null(grow)) {
                ## impute the data using prediction
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
          ##-------------------------------------
          ##
          ## imputation finished for this block
          ##
          ##-------------------------------------
          ## return NULL (memory management)
          NULL
        })
        ## return NULL (memory management)
        NULL
      })
      ##-------------------------------------
      ##     
      ## check whether algorithm has converged
      ##
      ##-------------------------------------
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
      ## return NULL (memory management)
      NULL
    })
  }
  ##--------------------------------------------------------------
  ##
  ## return the imputed data
  ##
  ##--------------------------------------------------------------
  invisible(data)
}
impute <- impute.rfsrc
