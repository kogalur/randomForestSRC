##  **********************************************************************
##  **********************************************************************
##  
##    RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
##  
##    This program is free software; you can redistribute it and/or
##    modify it under the terms of the GNU General Public License
##    as published by the Free Software Foundation; either version 3
##    of the License, or (at your option) any later version.
##  
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##  
##    You should have received a copy of the GNU General Public
##    License along with this program; if not, write to the Free
##    Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
##    Boston, MA  02110-1301, USA.
##  
##    ----------------------------------------------------------------
##    Project Partially Funded By: 
##    ----------------------------------------------------------------
##    Dr. Ishwaran's work was funded in part by DMS grant 1148991 from the
##    National Science Foundation and grant R01 CA163739 from the National
##    Cancer Institute.
##  
##    Dr. Kogalur's work was funded in part by grant R01 CA163739 from the 
##    National Cancer Institute.
##    ----------------------------------------------------------------
##    Written by:
##    ----------------------------------------------------------------
##      Hemant Ishwaran, Ph.D.
##      Director of Statistical Methodology
##      Professor, Division of Biostatistics
##      Clinical Research Building, Room 1058
##      1120 NW 14th Street
##      University of Miami, Miami FL 33136
##  
##      email:  hemant.ishwaran@gmail.com
##      URL:    http://web.ccs.miami.edu/~hishwaran
##      --------------------------------------------------------------
##      Udaya B. Kogalur, Ph.D.
##      Adjunct Staff
##      Department of Quantitative Health Sciences
##      Cleveland Clinic Foundation
##      
##      Kogalur & Company, Inc.
##      5425 Nestleway Drive, Suite L1
##      Clemmons, NC 27012
##  
##      email:  ubk@kogalur.com
##      URL:    https://github.com/kogalur/randomForestSRC
##      --------------------------------------------------------------
##  
##  **********************************************************************
##  **********************************************************************


impute.rfsrc <- function(formula,
                         data,
                         ntree = 500,
                         mtry = NULL,
                         xvar.wt = NULL,
                         nodesize = 1,
                         splitrule = NULL,
                         nsplit = 1,
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
  if (missing(data)) {
    stop("data is missing")
  }
  which.na <- is.na(data)
  if (!any(which.na) || all(which.na)) {
    return(invisible(data))
  }
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
  if (missing(mf.q)) {
    mforest <- FALSE
  }
    else {
      mforest <- TRUE
    }
  if (!missing(blocks)) {
    blocks <- cv.folds(nrow(data), max(1, blocks))
  }
    else {
      blocks <- list(1:nrow(data))
    }
  if (!mforest) {
    if (missing(formula)) {
      ytry <- min(p - 1, max(25, ceiling(sqrt(p))))
      formula <- as.formula(paste("Unsupervised(", ytry, ") ~ ."))
    }
    nullBlocks <- lapply(blocks, function(blk) {
      dta <- data[blk,, drop = FALSE]
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
    if (mf.q == 0) {
      stop("mf.q must be greater than zero")
    }
    if (mf.q >= 1) {
      mf.q <- min(p0 - 1, mf.q) / p0
    }
    K <- max(1 / mf.q, 2)
    data <- generic.impute.rfsrc(data = data,
                                 nimpute = 3,
                                 ntree = 250,
                                 mtry = mtry,
                                 nodesize = nodesize,
                                 nsplit = nsplit)$data
    diff.err <- Inf
    check <- TRUE
    nullWhile <- lapply(1:max.iter, function(m) {
      if (!check) {
        return(NULL)
      }
      if (verbose && max.iter > 1) {
        cat("\t iteration", m, "\n")
      }
      data.old <- data
      nullBlocks <- lapply(blocks, function(blk) {
        var.grp <- cv.folds(p0, K)
        nullObj <- lapply(var.grp, function(grp) {
          ynames <- unique(c(var.names[grp], all.var.names[always.use]))
          f <- as.formula(paste("Multivar(", paste(ynames, collapse = ","), paste(") ~ ."), sep = ""))
          dta <- data[blk,, drop = FALSE]
          dta[, ynames] <- lapply(ynames, function(nn) {
            xk <- data[, nn]
            xk[unlist(x.na[nn])] <- NA
            xk[blk]
          })
          retO <- tryCatch({generic.impute.rfsrc(f,
                                                 dta,
                                                 ntree = ntree,
                                                 nimpute = 1,
                                                 na.action =  na.action,
                                                 mtry = mtry,
                                                 nodesize = nodesize,
                                                 nsplit = nsplit)}, error = function(e) {NULL})
          if (!is.null(retO)) {
            if (!is.null(retO$missing$row)) {
              blk <- blk[-retO$missing$row]
            }
            if (!is.null(retO$missing$col)) {
              ynames <- ynames[-retO$missing$col]
            }
            data[blk, ynames] <<- retO$data[, ynames, drop = FALSE]
            rm(dta)
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
        cat("         >> ", diff.new.err, diff.err - diff.new.err, "\n")
      }
      check <<- ((diff.err - diff.new.err) >= eps)
      diff.err <<- diff.new.err
      rm(data.old)
      NULL
    })
  }
  invisible(data)
}
impute <- impute.rfsrc
