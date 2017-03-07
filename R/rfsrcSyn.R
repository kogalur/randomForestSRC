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


rfsrcSyn.rfsrc <-
  function(formula,          
           data,
           object,
           newdata,
           ntree = 1000,
           mtry = NULL,
           mtrySeq = NULL,
           nodesize = 5,
           nodesizeSeq = c(1:10,20,30,50,100),
           nsplit = 0,
           min.node = 3,
           use.org.features = TRUE,
           na.action = c("na.omit", "na.impute"),
           oob = TRUE,
           verbose = TRUE,
           ...
           )
{
  na.action <- match.arg(na.action, c("na.omit", "na.impute"))
  if (!missing(object)) {
    if (sum(inherits(object, c("rfsrc", "synthetic"), TRUE) == c(1, 2)) != 2) {
      stop("this function only works for objects of class `(rfsrc, synthetic)'")
    }
    M <- length(object$rfMachines)
    fmly <- object$rfMachines[[1]]$family
    xvar.names <- object$rfMachines[[1]]$xvar.names
    yvar.names <- object$rfMachines[[1]]$yvar.names
    rfMachines <- object$rfMachines
    rfSyn <- object$rfSyn
    synthetic <- object$synthetic
    opt.machine <- object$opt.machine
    list.names <- lapply(synthetic, function(ss) {colnames(ss)})
  }
  else {
    if (missing(formula) || missing(data)) {
      stop("need to specify 'formula' and 'data' or provide a grow forest object")
    }
    f <- as.formula(formula)
    if (na.action == "na.impute" && any(is.na(data))) {
        if (verbose) {
        cat("\t imputing the data\n")
      }
      data <- impute.rfsrc(data = data, ntree = ntree, nodesize = nodesize, nsplit = nsplit)
    }
    preObj <- rfsrc(f, data, ntree = 1, importance = "none",
                    nodesize = nrow(data), splitrule = "random")
    fmly <- preObj$family
    if (!(fmly == "regr" | fmly == "regr+" | fmly == "class" | fmly == "class+" | fmly == "mix+")) {
      stop("this function only applies to regression/classification based families")
    }
    xvar.names <- preObj$xvar.names
    yvar.names <- preObj$yvar.names
    preObj$yvar <- data.frame(preObj$yvar)
    colnames(preObj$yvar) <- yvar.names
    p <- length(xvar.names)
    if (is.null(mtrySeq)) {
      mtrySeq <- ceiling(p/3)
    }
    else {
      mtrySeq <- unique(ceiling(mtrySeq))
      mtrySeq <- mtrySeq[mtrySeq>=1 & mtrySeq <= p]
      if (length(mtrySeq) == 0) {
        stop("invalid choice for mtrySeq:", mtrySeq)
      }
    }
    nodesizeSeq <- sort(nodesizeSeq)
  }
  na.action <- match.arg(na.action, c("na.omit", "na.impute"))
  if (missing(object)) {
    if (oob) {
      samp.size <- nrow(rfsrc(f, data, ntree = 1, nodesize = nrow(data), splitrule = "random")$xvar)
      samp <- make.sample(ntree, samp.size)
    }
    rfMachines <- lapply(nodesizeSeq, function(nn) {
      lapply(mtrySeq, function(mm) {
        if (verbose) {
          cat("\t RF nodesize:", nn, "mtry:", mm, "\r")
        }
        if (oob) {
          rfsrc(f, data, ntree = ntree, mtry = mm, nodesize = nn,
                bootstrap = "by.user", samp = samp, 
                nsplit = nsplit, importance = "none")
        }
        else {
          rfsrc(f, data, ntree = ntree, mtry = mm, nodesize = nn, 
                nsplit = nsplit, importance = "none")
        }
      })
    })
    rfMachines <- unlist(rfMachines, recursive = FALSE)
    list.names <- paste(rep(nodesizeSeq, each = length(mtrySeq)), mtrySeq, sep = ".")
    M <- length(rfMachines)                         
    if (is.numeric(min.node) && min.node > 0) {
      good.machines <- which(sapply(1:M, function(m) {
        mean(rfMachines[[m]]$leaf.count, na.rm = TRUE) > min.node}))
      if (length(good.machines) == 0) {
        good.machines <- 1
      }
      list.names <- list.names[good.machines]
      rfMachines <- lapply(good.machines, function(m) {rfMachines[[m]]})
      M <- length(rfMachines)
    }
    names(rfMachines) <- paste("x.s.", list.names, sep = "")
    opt.machine <- rf.opt(rfMachines)
    if (verbose) {
      cat("\t making the synthetic features\n")
    }
    synthetic <- lapply(1:M, function(m) {
      do.call(cbind, lapply(rfMachines[[m]]$yvar.names, function(nn) {
        o.coerced <- coerce.multivariate(rfMachines[[m]], nn)
        yhat <- cbind(o.coerced$predicted.oob)
        J <- ncol(yhat)
        if (J > 1) {
          yhat <- yhat[ ,1:(J-1), drop = FALSE]
        }
        if (J > 2) {
          if (o.coerced$univariate) {
            colnames(yhat) <- paste(1:(J-1), list.names[m], sep = ".")
          }
            else {
              colnames(yhat) <- paste(nn, 1:(J-1), list.names[m], sep = ".")
            }
        }
          else {
            if (o.coerced$univariate) {
              colnames(yhat) <- paste(list.names[m], sep = ".")
            }
              else {
                colnames(yhat) <- paste(nn, list.names[m], sep = ".")
              }
          }
        yhat
      }))
    })
    x.s <- do.call("cbind", synthetic)
    list.names <- lapply(synthetic, function(ss) {colnames(ss)})
    names(synthetic) <- names(rfMachines)
    if (verbose) {
      cat("\t making the synthetic forest\n")
    }
    if (use.org.features) {
      data <- data.frame(preObj$yvar, preObj$xvar, x.s = x.s)
    }
      else {
        data <- data.frame(preObj$yvar, x.s = x.s)
      }
    rfSyn.f <- as.formula(paste("Multivar(", paste(yvar.names, collapse = ","), paste(") ~ ."), sep = ""))
    if (oob) {
      rfSyn <- rfsrc(rfSyn.f, data, ntree = ntree, mtry = mtry, nodesize = nodesize,
                     bootstrap = "by.user", samp = samp,                      
                     nsplit = nsplit, ... )
    }
    else {
      rfSyn <- rfsrc(rfSyn.f, data, ntree = ntree, mtry = mtry, nodesize = nodesize,
                     nsplit = nsplit, ... )
    }
  }
  if (!missing(newdata)) {
    if (na.action == "na.impute" && any(is.na(newdata))) {
      if (verbose) {
        cat("\t imputing the test data\n")
      }
      newdata <- impute.rfsrc(data = newdata, ntree = ntree, nodesize = nodesize, nsplit = nsplit)
    }
    if (verbose) {
      cat("\t making the test set synthetic features\n")
    }
    xtest <- newdata[, xvar.names, drop = FALSE]
    synthetic <- lapply(1:M, function(m) {
      predO <- predict(rfMachines[[m]], xtest, importance = "none")
      syn.o <- do.call(cbind, lapply(rfMachines[[m]]$yvar.names, function(nn) {
        o.coerced <- coerce.multivariate(predO, nn)
        yhat <- cbind(o.coerced$predicted)
        J <- ncol(yhat)
        if (J > 1) {
          yhat <- yhat[ ,1:(J-1), drop = FALSE]
        }
        yhat
      }))
      colnames(syn.o) <- list.names[[m]]
      syn.o
    })
    xtest.s <- do.call("cbind", synthetic)
    if (length(intersect(colnames(newdata), yvar.names) > 0) &&
        setequal(intersect(colnames(newdata), yvar.names), yvar.names)) {
      data.test <- data.frame(newdata[, yvar.names, drop = FALSE], x.s = xtest.s)
    }
    else {
      data.test <- data.frame(x.s = xtest.s)
    }
    if (use.org.features) {
      data.test <- data.frame(data.test, xtest)
    }
    rfSynPred <- predict(rfSyn, data.test, ...)
  }
  else {
    rfSynPred <- NULL
  }
  retObj <- list(rfMachines = rfMachines,
                 rfSyn = rfSyn,
                 rfSynPred = rfSynPred,
                 synthetic = synthetic,
                 opt.machine = opt.machine)
  if (oob) {
    class(retObj) <- c("rfsrc", "synthetic", "oob")
  }
  else {
    class(retObj) <- c("rfsrc", "synthetic", "inb")
  }
  retObj
}
rf.opt <- function(obj)
{
  which.min(sapply(1:length(obj), function(m) {
    mean(sapply(obj[[m]]$yvar.names, function(nn) {
      o.coerced <- coerce.multivariate(obj[[m]], nn)
      yhat <- o.coerced$predicted.oob
      if (o.coerced$family == "class") {
        brier(o.coerced$yvar, yhat)
      }
        else {
          yvar <- as.numeric(o.coerced$yvar) 
          mean((yvar - yhat)^2, na.rm = TRUE) / var(yvar, na.rm = TRUE) 
        }
    }), na.rm = TRUE)
  }))[1]
}
rfsrcSyn <- rfsrcSyn.rfsrc
