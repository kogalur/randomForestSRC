adrop3d.last <- function(x, d, keepColNames = FALSE) {
  ## this function is for arrays only
  if (!is.array(x)) {
    x
  }
    else {
      if (d > 1) {
        x[,,1:d, drop = FALSE]
      }
        else {
          if (dim(x)[1] == 1) {
            rbind(x[,,1, drop = TRUE])
          }
            else {
              if (dim(x)[2] == 1) {
                if (keepColNames) {
                  ##needed for J=1 classification pathology
                  xnew <- cbind(x[,,1, drop = TRUE])
                  colnames(xnew) <- colnames(x)
                  xnew
                }
                  else {
                    cbind(x[,,1, drop = TRUE])
                  }
              }
                else {
                  x[,,1, drop = TRUE]
                }
            }
        }
    }
}
adrop2d.first <- function(x, d, keepColNames = FALSE) {
  ## this function is for arrays only
  if (!is.array(x)) {
    x
  }
    else {
      if (d > 1) {
        x[1:d,, drop = FALSE]
      }
        else {
          x[1, , drop = TRUE]
        }
    }
}
adrop2d.last <- function(x, d, keepColNames = FALSE) {
  ## this function is for arrays only
  if (!is.array(x)) {
    x
  }
    else {
      if (d > 1) {
        x[,1:d, drop = FALSE]
      }
        else {
          x[,1, drop = TRUE]
        }
    }
}
amatrix <- function(x, d, names) {
  x <- matrix(x, d, dimnames = names)
  if (ncol(x) > 1) {
    x
  }
    else {
      c(x)
    }
}
amatrix.remove.names <- function(x) {
  if (!is.null(dim(x)) && ncol(x) == 1) {
    unlist(c(x), use.names = FALSE)
  }
    else {
      x
    }
}
atmatrix <- function(x, d, names, keep.names = FALSE) {
  x <- t(matrix(x, ncol = d, dimnames = names))
  if (ncol(x) > 1) {
    x
  }
    else {
      if (keep.names == FALSE) {
        c(x)
      }
        else {
          x.names <- rownames(x)
          x <- c(x)
          names(x) <- x.names
          x
        }
    }
}
avector <- function(x, name = FALSE) {
  if (!is.null(dim(x)) && nrow(x) > 1 && ncol(x) == 1) {
    x.names <- rownames(x)
    x <- unlist(c(x))
    if (name) names(x) <- x.names else names(x) <- NULL
    x
  }
    else if (!is.null(dim(x)) && nrow(x) == 1 && ncol(x) > 1) {
      x.names <- colnames(x)
      x <- unlist(c(x))
      if (name) names(x) <- x.names else names(x) <- NULL
      x
    }
      else if (!is.null(dim(x)) && nrow(x) == 1 && ncol(x) == 1) {
        unlist(c(x))
      }
        else {
          x
        }
}
available <- function (package, lib.loc = NULL, quietly = TRUE)
{
  package <- as.character(substitute(package))
  installed <- package %in% installed.packages()
  if (installed) {
    require(package, quietly = TRUE, character.only = TRUE)
  }
    else {
      return(invisible(FALSE))
    }
}
cv.folds <- function (n, folds = 10) {
  split(resample(1:n), rep(1:folds, length = n))
}
data.matrix <- function(x) {
  as.data.frame(lapply(x, function(xi) {
    if (is.integer(xi) || is.numeric(xi)) {
      xi
    }
      else if (is.logical(xi) || is.factor(xi)) {
        as.integer(xi)
      }
        else {
          as.numeric(xi)
        }
  }))
}
family.pretty <- function(x) {
  fmly <- x$family
  if (!is.null(x$forest$rfq) && x$forest$rfq) {
    fmly <- "rfq"
  }
  switch(fmly,
         "surv"     = "RSF",
         "surv-CR"  = "RSF",
         "regr"     = "RF-R",
         "class"    = "RF-C",
         "unsupv"   = "RF-U",
         "regr+"    = "mRF-R",
         "class+"   = "mRF-C",
         "mix+"     = "mRF-RC",
         "rfq"      = "RFQ"
         )
}
finalizeFormula <- function(formula.obj, data) {
  ## parse the formula object
  yvar.names <- formula.obj$yvar.names
  all.names  <- formula.obj$all.names
  index      <- length(yvar.names)
  fmly       <- formula.obj$family
  ytry       <- formula.obj$ytry
  ## total number of variables should exceed number of yvars
  if (length(all.names) <= index) {
    stop("formula is misspecified: total number of variables does not exceed total number of y-variables")
  }
  ## extract the xvar names
  if (all.names[index + 1] == ".") {
    if(index == 0) {
      xvar.names <- names(data)
    }
      else {
        xvar.names <- names(data)[!is.element(names(data), all.names[1:index])]
      }
  }
    else {
      if(index == 0) {
        xvar.names <- all.names
      }
        else {
          xvar.names <- all.names[-c(1:index)]
        }
      not.specified <- !is.element(xvar.names, names(data))
      if (sum(not.specified) > 0) {
        stop("formula is misspecified, object ", xvar.names[not.specified], " not found")
      }
    }
  ## return the goodies
  return (list(family=fmly, yvar.names=yvar.names, xvar.names=xvar.names, ytry=ytry))
}
finalizeData <- function(fnames, data, na.action, miss.flag = TRUE) {
  ## restrict data to the target variables
  data <- data[ , is.element(names(data), fnames), drop = FALSE]
  ## data conversion to numeric mode for factors
  ## no need for conversion if factors are not present
  factor.names <- unlist(lapply(data, is.factor))
  if (sum(factor.names) > 0) {
    data[, factor.names] <- data.matrix(data[, factor.names, drop = FALSE])
  }
  ## omit the na's
  if (miss.flag == TRUE && na.action == "na.omit") {
    ## no need to make the call if there is no missing data
    if (any(is.na(data))) {
      data <- na.omit(data)
    }
  }
  ## canvert nan to na's (github reported bug 06/27/2018)
  if (miss.flag == TRUE && na.action != "na.omit") {
    nan.names <- sapply(data, function(x){any(is.nan(x))})
    if (sum(nan.names) > 0) {
      data[, nan.names] <- data.frame(lapply(which(nan.names), function(j) {
        x <- data[, j]
        x[is.nan(x)] <- NA
        x
      }))
    }
  }
  ## is anything left?
  if (nrow(data) == 0) {
    stop("no records in the NA-processed data: consider using 'na.action=na.impute'")
  }
  ## convert logicals to 0/1 real (bug reported by John Ehrlinger)
  logical.names <- unlist(lapply(data, is.logical))
  if (sum(logical.names) > 0) {
    data[, logical.names] <- 1 * data[, logical.names, drop = FALSE]
  }
  ## characters not allowed (bug reported by John Ehrlinger)
  character.names <- unlist(lapply(data, is.character))
  if (sum(character.names) > 0) {
    stop("data types cannot be character: please convert all characters to factors")
  }
  return (data)
}
### AUC workhorse
get.auc.workhorse <- function(roc.data) {
  x <- roc.data[, 1][roc.data[, 2] == 1]
  y <- roc.data[, 1][roc.data[, 2] == 0]
  if (length(x) > 1 & length(y) > 1) {
    AUC  <- tryCatch({wilcox.test(x, y, exact=F)$stat/(length(x)*length(y))}, error=function(ex){NA})
  }
  else {
    AUC <- NA
  }
  AUC
}
### Multiclass AUC -- Hand & Till (2001) definition
get.auc <- function(y, prob) {
  if (is.factor(y)) {
    y.uniq <- levels(y)
  }
  else {
    y.uniq <- sort(unique(y))
  }
  nclass <- length(y.uniq)
  AUC <- NULL
  for (i in 1:(nclass - 1)) {
    for (j in (i + 1):nclass) {
      pt.ij <- (y == y.uniq[i] | y == y.uniq[j])
      if (sum(pt.ij) > 1) {
        y.ij <- y[pt.ij]
        pij <- prob[pt.ij, j]
        pji <- prob[pt.ij, i]
        Aij <-  get.auc.workhorse(cbind(pij, 1 * (y.ij == y.uniq[j])))
        Aji <-  get.auc.workhorse(cbind(pji, 1 * (y.ij == y.uniq[i])))
        AUC <- c(AUC, (Aij + Aji)/2)
      }
    } 
  }
  if (is.null(AUC)) {
    NA
  }
  else {
    mean(AUC, na.rm = TRUE)
  }
}                          
get.bayes.rule <- function(prob, pi.hat = NULL) {
  class.labels <- colnames(prob)
  if (is.null(pi.hat)) {
    factor(class.labels[apply(prob, 1, function(x) {
      if (!all(is.na(x))) {
        resample(which(x == max(x, na.rm = TRUE)), 1)
      }
      else {
        NA
      }
    })], levels = class.labels)
  }
  ## added to handle the rfq classifier
  else {
    minority <- which.min(pi.hat)
    majority <- setdiff(1:2, minority)      
    rfq.rule <- rep(majority, nrow(prob))
    rfq.rule[prob[, minority] >= min(pi.hat, na.rm = TRUE)] <- minority
    factor(class.labels[rfq.rule], levels = class.labels)
  }
}
## normalized brier (normalized to one for strawman coin toss)
get.brier.error <- function(y, prob) {
  if (is.null(colnames(prob))) {
    colnames(prob) <- levels(y)
  }
  cl <- colnames(prob)
  J <- length(cl)
  bs <- rep(NA, J)
  nullO <- sapply(1:J, function(j) {
    bs[j] <<- mean((1 * (y == cl[j]) - prob[, j]) ^ 2, na.rm = TRUE)
    NULL
  })
  norm.const <- (J / (J - 1))
  sum(bs * norm.const, na.rm = TRUE)
}
## misclassification error
get.misclass.error <- function(y, yhat) {
  cl <- sort(unique(y))
  err <- rep(NA, length(cl))
  for (k in 1:length(cl)) {
    cl.pt  <- (y == cl[k])
    if (sum(cl.pt) > 0) {
        err[k] <- mean(y[cl.pt] != yhat[cl.pt])
    }
  }
  err
}
## get confusion matrix
get.confusion <- function(y, class.or.prob) {
  ## response or probability?
  if (is.factor(class.or.prob)) {
    confusion <- table(y, class.or.prob)
  }
  else {
    if (is.null(colnames(class.or.prob))) {
      colnames(class.or.prob) <- levels(y)
    }
    confusion <- table(y, get.bayes.rule(class.or.prob))
  }
  class.error <- 1 - diag(confusion) / rowSums(confusion, na.rm = TRUE)
  cbind(confusion, class.error = round(class.error, 4))
}
## cindex
get.cindex <- function (time, censoring, predicted, do.trace = FALSE) {
  size <- length(time)
  if (size != length(time) |
      size != length(censoring) |
      size != length(predicted)) {
    stop("time, censoring, and predicted must have the same length")
  }
  miss <- is.na(time) | is.na(censoring) | is.na(predicted)
  nmiss <- sum(miss)
  if (nmiss == size) {
    stop("no valid pairs found, too much missing data")
  }
  ## Flag missing members so we can exclude them in the pairs.
  denom <- sapply(miss, function(x) if (x) 0 else 1)
  nativeOutput <- .Call("rfsrcCIndex",
                        as.integer(do.trace),
                        as.integer(size),
                        as.double(time),
                        as.double(censoring),
                        as.double(predicted),
                        as.integer(denom))
  ## check for error return condition in the native code
  if (is.null(nativeOutput)) {
    stop("An error has occurred in rfsrcCIndex.  Please turn trace on for further analysis.")
  }
  return (nativeOutput$err)
}
get.coerced.survival.fmly <- function(fmly, event.type, splitrule = NULL) {
  if (grepl("surv", fmly)) {
    ## assume no coercion
    coerced.fmly <- "surv"
    if (!is.null(splitrule)) {
      ## either competing risks or right censoring competing risks
      ## is coerced to right censoring in some settings
        if ((length(event.type) > 1) &&
            (splitrule != "l2.impute") &&
            (splitrule != "logrankscore")) {
        coerced.fmly <- "surv-CR"
      }
    }
      else {
        if (length(event.type) > 1) {
          coerced.fmly <- "surv-CR"
        }
      }
  }
    else {
      stop("attempt to coerce a non-survival family")
    }
  coerced.fmly
}
get.event.info <- function(obj, subset = NULL) {
  ## survival case
  if (grepl("surv", obj$family)) {
    if (!is.null(obj$yvar)) {
      if (is.null(subset)) {
        subset <- (1:nrow(cbind(obj$yvar)))
      }
      r.dim <- 2
      time <- obj$yvar[subset, 1]
      cens <- obj$yvar[subset, 2]
      ## censoring must be coded coherently
      if (!all(floor(cens) == abs(cens), na.rm = TRUE)) {
        stop("for survival families censoring variable must be coded as a non-negative integer")
      }
      ## Extract the unique event types.
      event <- na.omit(cens)[na.omit(cens) > 0]
      event.type <- sort(unique(event))
    }
    ##everything else
      else {
        r.dim <- 0
        event <- event.type <- cens <- cens <- time <- NULL
      }
    ## Set grid of time points.
    time.interest <- obj$time.interest
  }
    else {
      ## NULL for other families
      if ((obj$family == "regr+") | (obj$family == "class+")) {
        r.dim <- dim(obj$yvar)[2]
      }
        else {
          r.dim <- 1
        }
      event <- event.type <- cens <- time.interest <- cens <- time <- NULL
    }
  return(list(event = event, event.type = event.type, cens = cens,
              time.interest = time.interest, time = time, r.dim = r.dim))
}
## gmean for imbalanced classification
get.gmean <- function(y, prob, rfq = FALSE, robust = FALSE) {
  ## determine frequencies: exit if this is not a two-class problem
  frq <- table(y)  
  if (length(frq) > 2) {
    return(NULL)
  }
  ## determine threshold
  if (rfq) {
    threshold <- min(frq, na.rm = TRUE) / sum(frq, na.rm = TRUE)
  }
  else {
    threshold <- 0.5
  }
  ## convert yhat to a class label: 0 = majority, 1 = minority
  minority <- which.min(frq)
  majority <- setdiff(1:2, minority)
  yhat <- factor(1 * (prob[, minority] >= threshold), levels = c(0,1))
  ## compute the confusion matrix and extract TN,TP etc.
  confusion.matrix <- table(y, yhat)
  if (nrow(confusion.matrix) > 1) {##check that dimension is correct
    ## calculate the various rates
    TN <- confusion.matrix[minority, 2]
    FP <- confusion.matrix[minority, 1]
    FN <- confusion.matrix[majority, 2]
    TP <- confusion.matrix[majority, 1]
    ## assemble the sensitivity/specificity values
    if (robust) {##modified with 1 added to diagonals
      sensitivity <- (1 + TP) / (1 + TP + FN)
      specificity <- (1 + TN) / (1 + TN + FP)
    }
    else {
      sensitivity <- TP / (TP + FN)
      specificity <- TN / (TN + FP)
    }
    ## return the g mean
    sqrt(sensitivity * specificity)
  }
  else {
    NA
  }
}
get.grow.event.info <- function(yvar, fmly, need.deaths = TRUE, ntime) {
  if (grepl("surv", fmly)) {
    r.dim <- 2
    time <- yvar[, 1]
    cens <- yvar[, 2]
    ## censoring must be coded coherently
    if (!all(floor(cens) == abs(cens), na.rm = TRUE)) {
      stop("for survival families censoring variable must be coded as a non-negative integer (perhaps the formula is set incorrectly?)")
    }
    ## check if deaths are available (if user specified)
    if (need.deaths && (all(na.omit(cens) == 0))) {
      stop("no deaths in data!")
    }
    ## Check for event time consistency.
    ## we over-ride this now to allow for negative time (see Stute)
    ##if (!all(na.omit(time) >= 0)) {
    ##  stop("time must be  positive")
    ##}
    ## Extract the unique event types.
    event.type <- unique(na.omit(cens))
    ## Ensure they are all greater than or equal to zero.
    if (sum(event.type >= 0) != length(event.type)) {
      stop("censoring variable must be coded as NA, 0, or greater than 0.")
    }
    ## Discard the censored state, if it exists.
    event <- na.omit(cens)[na.omit(cens) > 0]
    event.type <- unique(event)
    ## Set grid of time points.
    nonMissingOutcome <- which(!is.na(cens) & !is.na(time))
    nonMissingDeathFlag <- (cens[nonMissingOutcome] != 0)
    time.interest <- sort(unique(time[nonMissingOutcome[nonMissingDeathFlag]]))
    ## trim the time points if the user has requested it
    ## we also allow the user to pass requested time points
    if (!missing(ntime)) {
      if (length(ntime) == 1 && length(time.interest) > ntime) {
        time.interest <- time.interest[
                                       unique(round(seq.int(1, length(time.interest), length.out = ntime)))]
      }
      if (length(ntime) > 1) {
        time.interest <- unique(sapply(ntime, function(tt) {
          time.interest[max(1, sum(tt >= time.interest, na.rm = TRUE))]
        }))
      }
    }
  }
  ## handle other families
    else {
      if ((fmly == "regr+") | (fmly == "class+") | (fmly == "mix+")) {
        r.dim <- dim(yvar)[2]
      }
        else {
          if (fmly == "unsupv") {
            r.dim <- 0
          }
            else {
              r.dim <- 1
            }
        }
      event <- event.type <- cens <- time.interest <- cens <- time <- NULL
    }
  return(list(event = event, event.type = event.type, cens = cens,
              time.interest = time.interest,
              time = time, r.dim = r.dim))
}
get.grow.mtry <- function (mtry = NULL, n.xvar, fmly) {
  if (!is.null(mtry)) {
    mtry <- round(mtry)
    if (mtry < 1 | mtry > n.xvar) mtry <- max(1, min(mtry, n.xvar))
  }
    else {
      if (grepl("regr", fmly)) {
        mtry <- max(ceiling(n.xvar/3), 1)
      }
        else {
          mtry <- max(ceiling(sqrt(n.xvar)), 1)
        }
    }
  return (mtry)
}
get.grow.nodesize <- function(fmly, nodesize) {
  ## Default node size for right-censored survival
  if (fmly == "surv"){
    if (is.null(nodesize)) {
      nodesize <- 15
    }
  }
  ## Default node size for competing risks
    else if (fmly == "surv-CR"){
      if (is.null(nodesize)) {
        nodesize <- 15
      }
    }
  ## Default node size for classification
      else if (fmly == "class" | fmly == "class+") {
        if (is.null(nodesize)) {
          nodesize <- 1
        }
      }
  ## Default node size for regression
        else if (fmly == "regr" | fmly == "regr+") {
          if (is.null(nodesize)) {
            nodesize <- 5
          }
        }
  ## Default node size for mixed outcomes
          else if (fmly == "mix+") {
            if (is.null(nodesize)) {
              nodesize <- 3
            }
          }
  ## Default node size for unsupervised splitting
            else if (fmly == "unsupv") {
              if (is.null(nodesize)) {
                nodesize <- 3
              }
            }
  ## The family is misspecified
              else if (is.null(nodesize)) {
                stop("family is misspecified")
              }
  ## Nodesize should be rounded if non-integer
  nodesize <- round(nodesize)
}
get.grow.splitinfo <- function (formula.detail, splitrule, hdim, nsplit, event.type) {
  ## CAUTION:  HARD CODED ON NATIVE SIDE
  splitrule.names <- c("logrank",              ##  1
                       "logrankscore",         ##  2
                       "logrankCR",            ##  3
                       "random",               ##  4
                       "mse",                  ##  5
                       "gini",                 ##  6
                       "unsupv",               ##  7
                       "mv.mse",               ##  8 --  reg/class/mix
                       "mv.gini",              ##  9 --  reg/class/mix
                       "mv.mix",               ## 10 --  reg/class/mix
                       "custom",               ## 11
                       "quantile.regr",        ## 12
                       "la.quantile.regr",     ## 13
                       "bs.gradient",          ## 14
                       "auc",                  ## 15
                       "entropy",              ## 16
                       "sg.regr",              ## 17
                       "sg.class",             ## 18
                       "sg.surv")              ## 19
fmly <- formula.detail$family
    ## Preliminary check for consistency.
    if (hdim > 0) {
        if(!is.null(nsplit)) {
            nsplit <- round(nsplit)    
            if (nsplit < 0) {
                stop("Invalid nsplit value.  Set nsplit >= 0.")
            }
        }
        else {
            nsplit = 1
        }
    }
    else {
        if(!is.null(nsplit)) {
            nsplit <- round(nsplit)    
            if (nsplit < 0) {
                stop("Invalid nsplit value.  Set nsplit >= 0.")
            }
        }
        else {
            nsplit = 0
        }
    }
  cust.idx <- NULL
  splitpass <- FALSE
  if (!is.null(splitrule)) {
    if(grepl("custom", splitrule)) {
      splitrule.idx <- which(splitrule.names == "custom")
      cust.idx <- as.integer(sub("custom", "", splitrule))
      if (is.na(cust.idx)) cust.idx <- 1
      splitpass <- TRUE
    }
      else if (splitrule == "random") {
        splitrule.idx <- which(splitrule.names == "random")
        ## Override the nsplit value in the case of pure random
        ## splitting.  It is set to the defalut value of one (1). In the native code,
        ## values greater than one are ignored in the case of pure random splitting!
        nsplit <- 1
        splitpass <- TRUE
      }
  }
  if (!splitpass) {
    if (grepl("surv", fmly)) {
      if (is.null(splitrule)) {
        ## No split rule specified, use default.
        if (length(event.type) ==  1) {
          splitrule.idx <- which(splitrule.names == "logrank")
        }
          else {
            splitrule.idx <- which(splitrule.names == "logrankCR")
          }
        splitrule <- splitrule.names[splitrule.idx]
      }
        else {
          ## User split rule specified.
          splitrule.idx <- which(splitrule.names == splitrule)
          if (length(splitrule.idx) != 1) {
            stop("Invalid split rule specified:  ", splitrule)
          }
          if ((length(event.type) ==  1) & (splitrule.idx == which(splitrule.names == "logrankCR"))) {
            stop("Cannot specify logrankCR splitting for right-censored data")
          }
          if ((length(event.type) >   1) & (splitrule.idx == which(splitrule.names == "logrank"))) {
            ## Override the splitrule to access the CR split rule.
            splitrule.idx <- which(splitrule.names == "logrankCR")
          }
        }
    }
    if (fmly == "class") {
      if (is.null(splitrule)) {
        ## No split rule specified, use default.
        splitrule.idx <- which(splitrule.names == "gini")
        splitrule <- splitrule.names[splitrule.idx]
      }
        else {
            ## User specified split rule.
            if ((splitrule != "rps") &
                (splitrule != "auc") &
                (splitrule != "entropy") &
                (splitrule != "gini") &
                (splitrule != "sg.class")) {
                stop("Invalid split rule specified:  ", splitrule)
            }
            splitrule.idx <- which(splitrule.names == splitrule)
        }
    }
    if (fmly == "regr") {
      if (is.null(splitrule)) {
        ## No split rule specified, use default.
        splitrule.idx <- which(splitrule.names == "mse")
        splitrule <- splitrule.names[splitrule.idx]
      }
        else {
          ## User specified split rule.
            if ((splitrule != "mse") &
                (splitrule != "sg.regr") &
                (splitrule != "la.quantile.regr") &
              (splitrule != "quantile.regr")) {
            stop("Invalid split rule specified:  ", splitrule)
          }
          splitrule.idx <- which(splitrule.names == splitrule)
        }
    }
    if (fmly == "regr+") {
      if (is.null(splitrule)) {
        ## No split rule specified, use default.
        splitrule.idx <- which(splitrule.names == "mv.mse")
        splitrule <- splitrule.names[splitrule.idx]
      }
        else {
          ## User specified split rule.
          if ((splitrule != "mv.mse")) {
            stop("Invalid split rule specified:  ", splitrule)
          }
          splitrule.idx <- which(splitrule.names == splitrule)
        }
    }
    if (fmly == "class+") {
      if (is.null(splitrule)) {
        ## No split rule specified, use default.
        splitrule.idx <- which(splitrule.names == "mv.gini")
        splitrule <- splitrule.names[splitrule.idx]
      }
        else {
          ## User specified split rule.
          if ((splitrule != "mv.gini")) {
            stop("Invalid split rule specified:  ", splitrule)
          }
          splitrule.idx <- which(splitrule.names == splitrule)
        }
    }
    if (fmly == "mix+") {
      if (is.null(splitrule)) {
        ## No split rule specified, use default.
        splitrule.idx <- which(splitrule.names == "mv.mse")
        splitrule <- "mv.mix"
      }
        else {
          ## User specified split rule.
          if ((splitrule != "mv.mix")) {
            stop("Invalid split rule specified:  ", splitrule)
          }
          splitrule.idx <- which(splitrule.names == splitrule)
        }
    }
    if (fmly == "unsupv") {
      if (is.null(splitrule)) {
        ## No split rule specified, use default.
        splitrule.idx <- which(splitrule.names == "unsupv")
        splitrule <- splitrule.names[splitrule.idx]
      }
        else {
          ## User specified split rule.
          if ((splitrule != "unsupv")) {
            stop("Invalid split rule specified:  ", splitrule)
          }
          splitrule.idx <- which(splitrule.names == splitrule)
        }
    }
  }
  splitinfo <- list(name = splitrule, index = splitrule.idx, cust = cust.idx, nsplit = nsplit)
  return (splitinfo)
}
get.importance.xvar <- function(importance.xvar, importance, object) {
  ## Check that importance has been requested
  if (!is.null(importance)) {
    ## Map vimp names to columns of GROW x-matrix
    ## Ensure names are coherent
    if (missing(importance.xvar) || is.null(importance.xvar)) {
      importance.xvar <- object$xvar.names
    }
      else {
        importance.xvar <- unique(importance.xvar)
        importance.xvar <- intersect(importance.xvar, object$xvar.names)
      }
    if (length(importance.xvar) == 0) {
      stop("xvar names do not match object xvar matrix")
    }
  }
    else {
      ## This was previously zero.  We set this so as to get length zero to the native code.
      importance.xvar <- NULL
    }
  return (importance.xvar)
}
 
get.mv.error <- function(obj, standardize = FALSE, pretty = TRUE, block = FALSE) {
  ## acquire yvar names - don't want "censoring" for surv and surv-CR
  nms <- NULL
  ynms <- obj$yvar.names
  if (obj$family == "surv" ||  obj$family == "surv-CR") {
    ynms <- ynms[1]
  }
  ## pretty option is not allowed for blocked error
  if (block) {
    pretty <- FALSE
  }
  ## loop over the yvar names acquring error if it is present
  err <- lapply(ynms, function(nn) {
    o.coerce <- coerce.multivariate(obj, nn)
    if (!block) {
      er <- o.coerce$err.rate
    }
    else {
      er <- o.coerce$err.block.rate
    }
    if (!is.null(er)) {
      if (o.coerce$family != "regr" || !standardize) {
        if (pretty && o.coerce$family == "class") {##pretty can only pull first vimp column "all"
            er <- utils::tail(cbind(er)[, 1], 1)
        }
        ## now returns survival error rates more coherently
        else {
          if (!block) {
            er <- utils::tail(er, 1)
            rownames(er) <- NULL
          }
        }
      }
      else {##standardized VIMP only applies to regression
        if (!block) {
          er <- utils::tail(er, 1) / var(o.coerce$yvar, na.rm = TRUE)
        }
        else {
          er <- er / var(o.coerce$yvar, na.rm = TRUE)
        }
      }
    }
    if (is.null(dim(er))) {
      nms <<- c(nms, paste(nn))
    }
    else {
      nms <<- c(nms, colnames(er))
    }
    er
  })
  ## if the first entry is NULL make the entire list NULL
  if (is.null(err[[1]])) {
    err <- NULL
  }  
  ## return as a vector for convenient interpretation?
  ## vector reduces information for classification to "all"
  if (!is.null(err)) {
    if (pretty) {
      err <- unlist(err)
      names(err) <- nms
    }
    else {
      names(err) <- ynms
    }
  }
  ### return the goodies
  err
}
get.mv.error.block <- function(obj, standardize = FALSE) {
  ## use get.mv.error but request blocked error rate
  get.mv.error(obj, standardize = standardize, block = TRUE)
}
get.mv.formula <- function(ynames) {
  as.formula(paste("Multivar(", paste(ynames, collapse = ","),paste(") ~ ."), sep = ""))
}
get.mv.predicted <- function(obj, oob = TRUE) {
  ## loop over the yvar names acquring the vimp if it is present
  ## acquire yvar names - don't want "censoring" for surv and surv-CR
  nms <- NULL
  ynms <- obj$yvar.names
  if (obj$family == "surv" ||  obj$family == "surv-CR") {
    ynms <- ynms[1]
  }
  pred <- do.call(cbind, lapply(ynms, function(nn) {
    o.coerce <- coerce.multivariate(obj, nn)
    if (o.coerce$family == "class" || o.coerce$family == "surv-CR") {
      ## ensembles can now be either oob, inbag, or all 05/06/2018
      if (!is.null(o.coerce$predicted.oob)) {
        nms <<- c(nms, paste(nn, ".", colnames(o.coerce$predicted.oob), sep = ""))
      }
      else {
        nms <<- c(nms, paste(nn, ".", colnames(o.coerce$predicted), sep = ""))
      }
    }
    else {
      nms <<- c(nms, paste(nn))
    }
    ## user may request OOB when it doesn't exist: noted 04/04/2018
    if (oob && !is.null(o.coerce$predicted.oob)) {
      o.coerce$predicted.oob
    }
    else {
      o.coerce$predicted
    }
  }))
  ## pretty names
  colnames(pred) <- nms
  pred
}
get.mv.vimp <- function(obj, standardize = FALSE, pretty = TRUE) {
  ## acquire yvar names - don't want "censoring" for surv and surv-CR
  ynms <- obj$yvar.names
  if (obj$family == "surv" ||  obj$family == "surv-CR") {
    ynms <- ynms[1]
  }
  ## loop over the yvar names acquring the vimp if it is present
  vmp <- lapply(ynms, function(nn) {
    o.coerce <- coerce.multivariate(obj, nn)
    v <- o.coerce$importance
    if (!is.null(v)) {
      if (o.coerce$family != "regr" || !standardize) {
        if (pretty && o.coerce$family == "class") {##pretty can only pull first vimp column "all"
          v <- v[, 1]
        }
      }
      else {##standardized VIMP only applies to regression
        v <- v / var(o.coerce$yvar, na.rm = TRUE)
      }
      if (is.null(dim(v))) {
        v <- cbind(v)
        colnames(v) <- nn
      }
    }
    v
  })
  ## if the first entry is NULL make the entire list NULL
  if (is.null(vmp[[1]])) {
    vmp <- NULL
  }  
  ## return as a matrix for convenient interpretation?
  ## matrix reduces information for classification to "all"
  if (!is.null(vmp)) {
    if (pretty) {
      vmp <- do.call(cbind, vmp)
    }
    else {
      names(vmp) <- ynms
    }
  }
  ### return the goodies
  vmp 
}
get.nmiss <- function(xvar, yvar = NULL) {
  if (!is.null(yvar)) {
    sum(apply(yvar, 1, function(x){any(is.na(x))}) | apply(xvar, 1, function(x){any(is.na(x))}))
  }
    else {
      sum(apply(xvar, 1, function(x){any(is.na(x))}))
    }
}
get.outcome.target <- function(family, yvar.names, outcome.target) {
  if (family == "regr" | family == "regr+" | family == "class" | family == "class+" | family == "mix+") {
    if (is.null(outcome.target)) {
      outcome.target <- yvar.names
    }
    ## Map target names to outcome names and ensure coherency.
    outcome.target <- unique(outcome.target)
    outcome.target <- intersect(outcome.target, yvar.names)
    if (length(outcome.target) == 0) {
      stop("yvar target names do not match object yvar names")
    }
    outcome.target <- match(outcome.target, yvar.names)
  }
    else {
      ## This is surv or surv-CR
      outcome.target <- 0
    }
}
get.univariate.target <- function(x, outcome.target = NULL) {
  ## This function takes a grow, grow-equivalent, or predict object and returns a single coherent target.
  ## That is, if no target has been specified, the first regression outcome with statistics is chosen.
  ## If no regression outcome exists, the first classification outcome with statistics is chosen.
  ## If the target is specified, the object is verified to contain the target outcome statistics
  ## for that y-var.  If none exist, the function will error.
  if (x$family == "regr+" | x$family == "class+" | x$family == "mix+") {
    if (is.null(outcome.target)) {
      ## Check the y-vars against regression and then classification.
      ## We choose the "first" variable, favoring regression, then
      ## classification.
      target <- match(c("regrOutput", "classOutput"), names(x))
      target <- target[!is.na(target)]
      if(length(target) > 0) {
        do.break <- FALSE
        for (i in target) {
          for (j in 1:length(x[[i]])) {
            if (length(x[[i]][[j]]) > 0) {
              ## This is a non-null output.
              outcome.target <- names(x[[i]][j])
              ## Exit the loop.
              do.break <- TRUE
              break
            }
          }
          if (do.break == TRUE) {
            break
          }
        }
      }
      else {
        ## Something would have to be seriously wrong for this to happen.
        stop("No outcomes found in object.  Please contact technical support.")
      }
    }
    else {
      ## Check that one and only one target has been specified.
      if (sum(is.element(outcome.target, x$yvar.names)) != 1) {
        stop("Specified target was not found or too many target outcomes were supplied (only one is allowed).")
      }
      ## A target outcome has been specified.  Verify that it contains outcome statistics.
      target <- match(c("regrOutput", "classOutput"), names(x))
      target <- target[!is.na(target)]
      found = FALSE
      if(length(target) > 0) {
        do.break <- FALSE
        for (i in target) {
          for (j in 1:length(x[[i]])) {
            if (length(x[[i]][[j]]) > 0) {
              ## This is a non-null output.
              if (outcome.target == names(x[[i]][j])) {
                found = TRUE
                ## Exit the loop.
                do.break <- TRUE
                break
              }
            }
          }
          if (do.break == TRUE) {
            break
          }
        }     
      }
      if (!found) {
        stop("Target outcome has been correctly specified but it did not contain outcome statistics.")
      }
    }
  }
  ## This function will return NULL if the function is not
  ## multivariate.  Otherwise, the outcome and its associated statistics is
  ## guaranteed to exist in the object.
  outcome.target
}
get.weight <- function(weight, n) {
  ## set the default weight
  if (!is.null(weight)) {
    if (any(weight < 0)      ||
        all(weight == 0)     ||
        length(weight) != n  ||
        any(is.na(weight))) {
      stop("Invalid weight vector specified.")
    }
  }
    else {
      weight <- rep(1, n)
    }
  return (weight)
}
get.ytry <- function(f) {
}
get.xvar.type <- function(generic.types, xvar.names, coerce.factor = NULL) {
  xvar.type <- generic.types
  if (!is.null(coerce.factor$xvar.names)) {
    xvar.type[is.element(xvar.names, coerce.factor$xvar.names)] <- "C"
  }
  xvar.type
}
get.xvar.nlevels <- function(nlevels, xvar.names, xvar, coerce.factor = NULL) {
  xvar.nlevels <- nlevels
  if (!is.null(coerce.factor$xvar.names)) {
    pt <- is.element(xvar.names, coerce.factor$xvar.names)
    xvar.nlevels[pt] <- sapply(coerce.factor$xvar.names, function(nn) {max(xvar[, nn])})
  }
  xvar.nlevels
}
get.yvar.type <- function(fmly, generic.types, yvar.names, coerce.factor = NULL) {
  if (fmly == "unsupv") {
    yvar.type <- NULL
  }
    else {
      if (grepl("surv", fmly)) {
        yvar.type <- c("T", "S")
      }
        else {
          yvar.type <- generic.types
          ## is coerce.factor at play for the y-outcomes?
          if (!is.null(coerce.factor$yvar.names)) {
            yvar.type[is.element(yvar.names, coerce.factor$yvar.names)] <- "C"
          }
        }
    }
  yvar.type
}
get.yvar.nlevels <- function(fmly, nlevels, yvar.names, yvar, coerce.factor = NULL) {
  if (fmly == "unsupv") {
    yvar.nlevels <- NULL
  }
  else {
    yvar.nlevels <- nlevels
    if (!is.null(coerce.factor$yvar.names)) {
      pt <- is.element(yvar.names, coerce.factor$yvar.names)
      yvar.nlevels[pt] <- sapply(coerce.factor$yvar.names, function(nn) {max(yvar[, nn])})
    }
  }
    yvar.nlevels
}
global.prob.assign <- function(prob, prob.epsilon, gk.quantile, quantile.regr, splitrule, n) {
  ## default values
  prob.default <- (1:99) / 100
  prob.epsilon.default <- .005
  prob.bs.default <- .9
  ## is quantile.regr in effect or is a quantile regression splitting rule in place?
  if (quantile.regr  || grepl("quantile.regr", splitrule) || gk.quantile) {
    ## is gk requested?
    if (gk.quantile) {
      ## gk quantile requires prob and prob.epsilon
      ## if both are missing, set to optimal value -- SLOW SLOW !!!!
      if (is.null(prob) && is.null(prob.epsilon)) {
        n <- max(n, 1)
        prob <- (1 : (2 * n)) / (2 * n + 1)
        prob.epsilon <- diff(prob)[1] * .99
      }
      ## if prob missing, set to default value
      else if (is.null(prob) && !is.null(prob.epsilon)) {
        prob <- prob.default
      }
      ## if prob.epsilon missing, set to default value
      else if (!is.null(prob) && is.null(prob.epsilon)) {
        prob.epsilon <- mean(diff(sort(prob)), na.rm = TRUE) * .99
      }
      ## neither are missing
      else {
        ##nothing        
      }
    }
    ## quantile splitting is in effect but no gk estimation
    else {
      prob.epsilon <- prob.epsilon.default
      ## here's where we set the default prob if its missing
      if (is.null(prob)) {
        prob <-  prob.default
      }
    }
    ## make sure prob values are coherent
    if (sum(prob > 0 & prob < 1) == 0) {
      stop("parameter 'prob' is not between (0,1)", prob)
    }
    prob <- sort(prob[prob>0 & prob<1])
  }
  ## survival has quantile splitting rules, here's where we deal with those
  if (splitrule ==  "bs.gradient") {
    if (is.null(prob)) {
      prob <- prob.bs.default
    }
    ## make sure prob values are coherent
    if (sum(prob > 0 & prob < 1) == 0) {
      stop("parameter 'prob' is not between (0,1)", prob)
    }
    prob <- prob[prob>0 & prob<1][1]
  }
  ## for all other families leave the values at default NULL
  #if (is.null(prob)) {
  #  prob <-  prob.default
  #}
  #if (is.null(prob.epsilon)) {
  #  prob.epsilon <-  prob.epsilon.default
  #}
  return(list(prob = prob, prob.epsilon = prob.epsilon))
}
make.samplesize.function <- function(fraction = 1) {
  f <- paste("x * ", paste(eval(fraction)))
  expr <- parse(text = f)
  function(x) eval(expr, list(x = x))
}
make.holdout.array <- function(vtry = 0, p, ntree, ntree.allvars = NULL) {
  ##default is triggered by vtry = 0
  if (vtry == 0) {
    return(NULL)
  }
  ## check that 1<= vtry <= p-1
  if (vtry < 0 || vtry >= p) {
    stop("vtry must be positive and less than the feature dimension")
  }
  ## everything is OK, go ahead and make the p x ntree array
  holdout <- do.call(cbind, lapply(1:ntree, function(b) {
    ho <- rep(0, p)
    ho[sample(1:p, size = vtry, replace = FALSE)] <- 1
    ho
  }))
  if (is.null(ntree.allvars)) {
    holdout
  }
  ## otherwise pad the array with 0's so that no variables are held out
  ## put them at the front for split-optimization  
  else {
    cbind(matrix(0, p, ntree.allvars), holdout)
  }
}
## imbalanced sampling - uses SWOR
make.imbalanced.sample <- function(ntree, ratio = 0.5, y, replace = FALSE) {
  ## ratio must be between 0 and 1
  if (ratio < 0 | ratio > 1) {
    stop("undersampling ratio must be between 0 and 1:", ratio)
  }
  ## determine minority/majority frequencies
  frq <- table(y)
  class.labels <- names(frq)
  minority <- which.min(frq)
  majority <- setdiff(1:2, minority)
  n.minority <- min(frq, na.rm = TRUE)
  n.majority <- max(frq, na.rm = TRUE)
  samp.size <- length(y)
  ## ratio must be larger than minority ratio
  ratio <- max((2 + n.minority) / length(y), ratio)
  ## use rbind to ensure a matrix
  rbind(sapply(1:ntree, function(bb){
    inb <- rep(0, samp.size)
    smp.min <- sample(which(y == class.labels[minority]), size = n.minority, replace = TRUE)
    smp.maj <- sample(which(y == class.labels[majority]), size = ratio * n.majority, replace = replace)
    smp <- c(smp.min, smp.maj)
    inb.frq <- tapply(smp, smp, length)
    idx <- as.numeric(names(inb.frq))
    inb[idx] <- inb.frq
    inb
  }))
}
## modified to default to "SWR"
## make inbag bootstrap samples (for use with bootstrap = "by.user")
## allows for under/over sampling
make.sample <- function(ntree, samp.size, boot.size = NULL, replace = FALSE) {
  ## samp.size cannot be negative
  if (samp.size < 0) {
    stop("samp.size cannot be negative:", samp.size)
  }
  ## default setting
  if (is.null(boot.size)) {
    if (replace == TRUE) {
      boot.size <- samp.size
    }
    else {
      boot.size <- .632 * samp.size
    }
  }
  ## use rbind at the end to ensure a matrix
  rbind(sapply(1:ntree, function(bb){
    inb <- rep(0, samp.size)
    smp <- sample(1:samp.size, size = boot.size, replace = replace)
    frq <- tapply(smp, smp, length)
    idx <- as.numeric(names(frq))
    inb[idx] <- frq
    inb
  }))
}
 
## make stratified inbag bootstrap samples for imbalanced two class problem
## allows for arbitrary under-sampling ratio of the majority class
## sample size for subsampling
make.size <- function(y) {
  ## extract the relative frequencies
  frq <- table(y)
  min(length(y), min(frq, na.rm = TRUE) * length(frq))
}
## imbalanced weights for subsampling
make.wt <- function(y) {
  ## extract the relative frequencies
  frq <- table(y)
  class.labels <- names(frq)
  wt <- rep(1, length(y))
  ## weights for class j equal product of all non-j class frequencies
  nullO <- sapply(1:length(frq), function(j) {
    wt[y == class.labels[j]] <<- prod(frq[-j], na.rm = TRUE)
    NULL
  })
  ## weights can become large, so divide by the maximum value
  ## this ensures 0 < wt <= 1
  wt / max(wt, na.rm = TRUE)
}
parseFormula <- function(f, data, ytry = NULL, coerce.factor = NULL) {
  ## confirm coherency of the formula
  if (!inherits(f, "formula")) {
    stop("'formula' is not a formula object.")
  }
  if (is.null(data)) {
    stop("'data' is missing.")
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  ## pull the family and y-variable names
  fmly <- all.names(f, max.names = 1e7)[2]
  all.names <- all.vars(f, max.names = 1e7)
  yvar.names <- all.vars(formula(paste(as.character(f)[2], "~ .")), max.names = 1e7)
  yvar.names <- yvar.names[-length(yvar.names)]
  ## is coerce.factor at play for the y-outcomes?
  coerce.factor.org <- coerce.factor
  coerce.factor <- vector("list", 2)
  names(coerce.factor) <- c("xvar.names", "yvar.names")
  if (!is.null(coerce.factor.org)) {
    coerce.factor$yvar.names <- intersect(yvar.names, coerce.factor.org)
    if (length(coerce.factor$yvar.names) == 0) {
      coerce.factor$yvar.names <- NULL
    }
    coerce.factor$xvar.names <- intersect(setdiff(colnames(data), yvar.names), coerce.factor.org)
  }
  ## survival forests
  if ((fmly == "Surv")) {
    if (sum(is.element(yvar.names, names(data))) != 2) {
      stop("Survival formula incorrectly specified.")
    }
    family <- "surv"
    ytry <- 2
  }
  ## multivariate forests
    else if ((fmly == "Multivar" || fmly == "cbind")  && length(yvar.names) > 1) {
      if (sum(is.element(yvar.names, names(data))) < length(yvar.names)) {
        stop("Multivariate formula incorrectly specified: y's listed in formula are not in data.")
      }
      ## determine the family: now handles mixed outcomes
      Y <- data[, yvar.names, drop = FALSE]
      ## Convert to 0/1 real (bug reported by John Ehrlinger)
      logical.names <- unlist(lapply(Y, is.logical))
      if (sum(logical.names) > 0) {
        Y[, logical.names] <- 1 * Y[, logical.names, drop = FALSE]
      }
      ## are all the responses factors?
      ## caution: ordered factors are factors!
      if ((sum(unlist(lapply(Y, is.factor))) + 
          length(coerce.factor$yvar.names)) == length(yvar.names)) {
        family <- "class+"
      }
      ## are all the responses continuous?
      ## caution: ordered factors are factors!
      else if ((sum(unlist(lapply(Y, is.factor))) + 
          length(coerce.factor$yvar.names)) == 0) {
        family <- "regr+"
      }
      ## are the responses a combination of factors and continuous?
      ## caution: ordered factors are factors!
      else if (((sum(unlist(lapply(Y, is.factor))) +
                 length(coerce.factor$yvar.names)) > 0) && 
               ((sum(unlist(lapply(Y, is.factor))) +
                 length(coerce.factor$yvar.names)) < length(yvar.names))) {
        family <- "mix+"
      }
      ## failure
        else {
          stop("y-outcomes must be either real or factors in multivariate forests.")
        }
      if (!is.null(ytry)) {
        ## Check that incoming ytry is consistent.
        if ((ytry < 1) || (ytry > length(yvar.names))) {
          stop("invalid value for ytry:  ", ytry)
        }
      }
        else {
          ytry <- length(yvar.names)
        }
    }
  ## unsupervised forests
      else if (fmly == "Unsupervised") {
        ## unsupervised forests
        if (length(yvar.names) != 0) {
          stop("Unsupervised forests require no y-responses")
        }
        family <- "unsupv"
        yvar.names <- NULL
        ## Strip away the family from the formula, leaving ytry.
        temp <- gsub(fmly, "", as.character(f)[2])
        temp <- gsub("\\(|\\)", "", temp)
        ytry <- as.integer(temp)
        if (is.na(ytry)) {
          ytry <- 1
        }
          else {
            if (ytry <= 0) {
              stop("Unsupervised forests require positive ytry value")
            }
          }
      }
  ## univariate forests (regression or classification)
        else {
          ## must be a (univariate) regresssion or classification
          if (sum(is.element(yvar.names, names(data))) != 1) {
            stop("formula is incorrectly specified.")
          }
          Y <- data[, yvar.names]
          ## logicals are treated as 0/1 real (bug reported by John Ehrlinger)
          if (is.logical(Y)) {
            Y <- as.numeric(Y)
          }
          ## check whether we have a factor or a continuous variable
          if (!(is.factor(Y) | is.numeric(Y))) {
            stop("the y-outcome must be either real or a factor.")
          }
          if (is.factor(Y) || length(coerce.factor$yvar.names) == 1) {
            family <- "class"
          }
            else {
              family <- "regr"
            }
          ytry <- 1
        }
  ## done: return the goodies
  return (list(all.names=all.names, family=family, yvar.names=yvar.names, ytry=ytry,
               coerce.factor = coerce.factor))
}
is.all.na <- function(x) {all(is.na(x))}
parseMissingData <- function(formula.obj, data) {
  ## parse the formula object
  yvar.names <- formula.obj$yvar.names
  if (length(yvar.names) > 0) {
    resp <- data[, yvar.names, drop = FALSE]
    ## determine whether all the outcomes are missing
    ## works for any dimension
    col.resp.na <- unlist(lapply(data[, yvar.names, drop = FALSE], is.all.na))
    if (any(col.resp.na)) {
      stop("All records are missing for one (or more) yvar(s)")
    }
  }
  ## remove all x variables with missing values for all records
  colPt <- unlist(lapply(data, is.all.na))
  ## terminate if all columns have all missing data
  if (sum(colPt) > 0 && sum(colPt) >= (ncol(data) - length(yvar.names))) {
    stop("All x-variables have all missing data:  analysis not meaningful.")
  }
  data <- data[, !colPt, drop = FALSE]
  ## remove all records with missing values for all outcomes(s) and xvar
  rowPt <- apply(data, 1, is.all.na)
  if (sum(rowPt) == nrow(data)) {
    stop("Rows of the data have all missing data:  analysis not meaningful.")
  }
  data <- data[!rowPt,, drop = FALSE]
  ## return the NA processed data
  return(data)
}
## robust resampling
resample <- function(x, size, ...) {
  if (length(x) <= 1) {
    if (!missing(size) && size == 0) x[FALSE] else x
  }
    else {
      sample(x, size, ...)
    }
}
## determine rows and columns missing from dat
row.col.deleted <- function(dat, r.n, c.n)
{
  which.r <- setdiff(r.n, rownames(dat))
  if (length(which.r) > 0) {
    which.r <- match(which.r, r.n)
  }
    else {
      which.r <- NULL
    }
  which.c <- setdiff(c.n, colnames(dat))
  if (length(which.c) > 0) {
    which.c <- match(which.c, c.n)
  }
    else {
      which.c <- NULL
    }
  return(list(row = which.r, col = which.c))
}
#weighted gini/entropy performance metric
#mode is equal to either 'gini' or 'entropy'
sid.perf.metric <- function(truth,cluster,mode=c("entropy", "gini")){
  ## verify mode option
  mode <- match.arg(mode, c("entropy", "gini"))
  ## confusion matrix
  k=length(unique(truth))
  tab=table(truth,cluster)
  clustersizes=colSums(tab)
  clustersizesnorm=clustersizes/sum(clustersizes)
  tabprop=tab
  lapply(1:(dim(tab)[2]),function(i){
    tabprop[,i]<<-tabprop[,i]/clustersizes[i]
    NULL
  })
  if(mode=="entropy"){
    measure=0
    maxmeasure=0
    lapply(1:(dim(tabprop)[2]),function(i){
      clustermeasure=0
      maxclustermeasure=0
      lapply(1:(dim(tabprop)[1]),function(j){
        if(tabprop[j,i]!=0){
          clustermeasure<<-clustermeasure+-tabprop[j,i]*log2(tabprop[j,i])
        }
        maxclustermeasure<<-maxclustermeasure+-1/k*log2(1/k)
        NULL
      })
      measure<<-measure+clustersizesnorm[i]*clustermeasure
      maxmeasure<<-maxmeasure+clustersizesnorm[i]*maxclustermeasure
    })
  }
  if(mode=="gini"){
    measure=0
    maxmeasure=0
    lapply(1:(dim(tabprop)[2]),function(i){
      clustermeasure=1
      maxclustermeasure=1
      lapply(1:(dim(tabprop)[1]),function(j){
        if(tabprop[j,i]!=0){
          clustermeasure<<-clustermeasure+(-tabprop[j,i]^2)
        }
        maxclustermeasure<<-maxclustermeasure+(-1/k^2)
        NULL
      })
      measure<<-measure+clustersizesnorm[i]*clustermeasure
      maxmeasure<<-maxmeasure+clustersizesnorm[i]*maxclustermeasure
      NULL
    })
  }
  list(result=measure,measure=mode,normalized_measure=measure/maxmeasure)
}
