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
assign.impute.mean <- function(data, impute.mean) {
  d <- data.frame(mclapply(colnames(data), function(xnms) {
    x <- data[, xnms]
    is.na.x <- is.na(x)
    if (any(is.na.x)) {
      x[is.na.x] <- impute.mean[[xnms]]
    }
    x
  }), stringsAsFactors = TRUE)
  colnames(d) <- colnames(data)
  d
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
digits.pretty <- function(x, digits = 8) {
  paste(round(x, digits), collapse=", ", sep = "")
}
family.pretty <- function(x) {
  fmly <- x$family
  if (!is.null(x$forest$rfq) && x$forest$rfq) {
    fmly <- "rfq"
  }
  switch(fmly,
         "surv"     = "RSF",
         "surv-CR"  = "RSF",
         "surv-TDC" = "RHF",
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
    subj.names <- formula.obj$subj.names
    all.names  <- formula.obj$all.names
    fmly       <- formula.obj$family
    ytry       <- formula.obj$ytry
    index <- length(yvar.names)
    ## Adjust the index for the presence of subject names.
    if (fmly == "surv") {
        if (!is.null(subj.names)) {
            index <- index + 1
        }
    }
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
    return (list(family=fmly, subj.names=subj.names, yvar.names=yvar.names, xvar.names=xvar.names, ytry=ytry))
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
get.coerced.survival.fmly <- function(fmly, subj, event.type, splitrule = NULL) {
    if (grepl("surv", fmly)) {
        ## assume no coercion
        coerced.fmly <- "surv"
        if (!is.null(splitrule)) {
            if (is.null(subj)) {             
                ## either competing risks or right censoring competing risks
                ## is coerced to right censoring in some settings
                if ((length(event.type) > 1) &&
                    (splitrule != "l2.impute") &&
                    (splitrule != "logrankscore")) {
                    coerced.fmly <- "surv-CR"
                }
            }
            else {
                coerced.fmly <- "surv-TDC"
            }
        }
        else {
            if (is.null(subj)) { 
                if (length(event.type) > 1) {
                    coerced.fmly <- "surv-CR"
                }
            }
            else {
                coerced.fmly <- "surv-TDC"
            }
        }
    }
    else {
        stop("attempt to coerce a non-survival family")
    }
    coerced.fmly
}
## get grow mtry 
get.grow.mtry <- function (mtry = NULL, n.xvar, fmly, splitrule = NULL) {
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
  if (!is.null(splitrule) && splitrule == "random") {
    mtry <- 1
  }
  mtry
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
  ## Default node size for time dependent covariates
    else if (fmly == "surv-TDC"){
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
get.grow.splitinfo <- function (formula.detail, splitrule, hdim, nsplit, event.info) {
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
                       "sg.surv",              ## 19
                       "tdc.gradient",         ## 20
                       "mahalanobis",          ## 21
                       "logrankCRGeneral")     ## 22
  ## set the family
  fmly <- formula.detail$family
  ## initialization
  cust.idx <- NULL
  splitpass <- FALSE
  ## set splitpass
  if (!is.null(splitrule)) {
    if(grepl("custom", splitrule)) {
      splitrule.idx <- which(splitrule.names == "custom")
      cust.idx <- as.integer(sub("custom", "", splitrule))
      if (is.na(cust.idx)) cust.idx <- 1
      splitpass <- TRUE
    }
    else if (splitrule == "random") {
      splitrule.idx <- which(splitrule.names == "random")
      splitpass <- TRUE
    }
  }
  ## do this unless pass given for split 
  if (!splitpass) {
    ## if splitrule is present, confirm it is correctly set
    if (!is.null(splitrule)) {
      splitrule <- match.arg(splitrule, splitrule.names)
    }
    ## survival
    if (grepl("surv", fmly)) {
      if (is.null(splitrule)) {
        ## No split rule specified, use default.
        if (event.info$r.dim == 2) {
          ## Survival or Competing Risk:
          if (length(event.info$event.type) ==  1) {
            splitrule.idx <- which(splitrule.names == "logrank")
          }
          else {
            splitrule.idx <- which(splitrule.names == "logrankCR")
          }
        }
        else if (event.info$r.dim == 3) {
          ## Time Dependent Covariates:                    
          splitrule.idx <- which(splitrule.names == "tdc.gradient")
        }
        else {
          stop("Invalid r.dim encountered in split rule:  ", event.info$r.dim)
        }
        splitrule <- splitrule.names[splitrule.idx]
      }
      else {
        ## User split rule specified.
        splitrule.idx <- which(splitrule.names == splitrule)
        if (length(splitrule.idx) != 1) {
          stop("Invalid split rule specified for survival:  ", splitrule)
        }
        if (event.info$r.dim == 2) {
          if ((length(event.info$event.type) ==  1) & (splitrule.idx == which(splitrule.names == "logrankCR"))) {
            stop("Cannot specify logrankCR splitting for right-censored data")
          }
          if ((length(event.info$event.type) >   1) & (splitrule.idx == which(splitrule.names == "logrank"))) {
              ## Override the splitrule to access the generalized CR split rule (see 3.3.1 in the paper).
              ## The default splitrule is 3.3.2 or Gray's test.
            splitrule.idx <- which(splitrule.names == "logrankCRGeneral")
          }
        }
        else if (event.info$r.dim == 3) {
          if (splitrule.idx != which(splitrule.names == "tdc.gradient")) {
            stop("Must specify tdc.gradient for time dependent covariates")                        
          }
        }
        else {
          stop("Invalid r.dim encountered in split rule:  ", event.info$r.dim)
        }
      }
    }
    ## class
    if (fmly == "class") {
      if (is.null(splitrule)) {
        ## No split rule specified, use default.
        splitrule.idx <- which(splitrule.names == "gini")
        splitrule <- splitrule.names[splitrule.idx]
      }
      else {
        ## User specified split rule.
        if ((splitrule != "auc") &
            (splitrule != "entropy") &
            (splitrule != "gini") &
            (splitrule != "sg.class")) {
          stop("Invalid split rule specified for classification:  ", splitrule)
        }
        splitrule.idx <- which(splitrule.names == splitrule)
      }
    }
    ## regression
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
          stop("Invalid split rule specified for regression:  ", splitrule)
        }
        splitrule.idx <- which(splitrule.names == splitrule)
      }
    }
    ## multivariate
    if (fmly == "regr+") {
      if (is.null(splitrule)) {
        ## No split rule specified, use default
        splitrule.idx <- which(splitrule.names == "mv.mse")
        splitrule <- splitrule.names[splitrule.idx]
      }
      else {
        ## User specified split rule
        if ((splitrule != "mv.mse") & (splitrule != "mahalanobis")) {
          stop("Invalid split rule specified for multivariate regression:  ", splitrule)
        }
        splitrule.idx <- which(splitrule.names == splitrule)
      }
    }
    if (fmly == "class+") {
      if (is.null(splitrule)) {
        ## No split rule specified, use default
        splitrule.idx <- which(splitrule.names == "mv.gini")
        splitrule <- splitrule.names[splitrule.idx]
      }
      else {
        ## User specified split rule.
        if ((splitrule != "mv.gini")) {
          stop("Invalid split rule specified for multivariate classsification:  ", splitrule)
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
          stop("Invalid split rule specified for mixed multivariate regression:  ", splitrule)
        }
        splitrule.idx <- which(splitrule.names == splitrule)
      }
    }
    ## unsupervised
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
  }## completes !splitpass
  ## set nsplit based on splitrule
  if(!is.null(nsplit)) {
    nsplit <- round(nsplit)    
    if (nsplit < 0) {
      stop("Invalid nsplit value: set nsplit >= 0")
    }
  }
  else {
    ## pure random
    if (splitrule.idx == 4) {
      nsplit <- 1
    }
    ## fast families
    else if (splitrule.idx %in% c(5, 6)) {
      nsplit <- 0
    }
    ## sg families
    else if (splitrule.idx %in% c(17, 18, 19)) {
      nsplit <- 10 ## can change this as needed later
    }
    else {
      nsplit <- 10
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
get.impute.mean <- function(data) {
  imean <- mclapply(data, function(x) {
    if (all(is.na(x))) {
      NA
    }
    else {
      if (is.factor(x)) {
        x.table <- table(x)
        names(x.table)[which.max(x.table)]
      }
      else {
        mean(x, na.rm = TRUE)
      }
    }
  })
  names(imean) <- colnames(data)
  imean
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
get.rfnames <- function(hidden = TRUE, stealth = FALSE) {
  rfnames <- names(formals(rfsrc))
  if (hidden) {
    rfnames <- c(rfnames,              
               "impute.only",
               "presort.xvar",
               "experimental",
               "rfq",
               "perf.type",
               "gk.quantile",
               "prob",
               "prob.epsilon",
               "vtry",
               "holdout.array")
  }
   
  rfnames
}
get.weight <- function(weight, n) {
  ## set the default weight
  if (!is.null(weight)) {
    if (any(weight < 0)      ||
##        all(weight == 0)     ||
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
get.yvar.type <- function(fmly, generic.types, yvar.names) {
  if (fmly == "unsupv") {
    yvar.type <- NULL
  }
    else {
      if (grepl("surv", fmly)) {
          if (length(yvar.names) == 2) {
              yvar.type <- c("T", "S")
          }
          else {
              yvar.type <- c("t", "T", "S")
          }
      }
        else {
          yvar.type <- generic.types
        }
    }
  yvar.type
}
get.yvar.nlevels <- function(fmly, nlevels, yvar.names, yvar) {
  if (fmly == "unsupv") {
    NULL
  }
  else {
    nlevels
  }
}
get.numeric.levels <- function(fmly, nlevels, gvar) {
    gvar.numeric.levels  <- lapply(1:length(nlevels),
                                   function(nn) {if(nlevels[nn] > 0) unique(sort(gvar[, nn])) else NULL})
    ## Remove null elements in the list
    gvar.numeric.levels <- gvar.numeric.levels[!sapply(gvar.numeric.levels,is.null)]
    ## We are uncomfortable in sending a ist of length zero into the C-code, so we add an additional check.
    if (length(gvar.numeric.levels) == 0) gvar.numeric.levels = NULL
    gvar.numeric.levels
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
    ho[resample(1:p, size = vtry, replace = FALSE)] <- 1
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
    smp.min <- resample(which(y == class.labels[minority]), size = n.minority, replace = TRUE)
    smp.maj <- resample(which(y == class.labels[majority]), size = ratio * n.majority, replace = replace)
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
    smp <- resample(1:samp.size, size = boot.size, replace = replace)
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
  ## Default scenario, no subject information when family is not
  ## time dependent covariates.  Can be overridden later.
  subj.names <- NULL
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
  if (fmly == "Surv") {
      ## Survival and competing risk will have 2 slots, namely time and censoring.
      ## Time dependent covariates will have 4 slots, namely id, start, stop, and event.
      ## If TDC is in effect, we remove the id from the yvars, and tag is an the subject identifier.
      if ((sum(is.element(yvar.names, names(data))) != 2) &&
          (sum(is.element(yvar.names, names(data))) != 4)) {
          stop("Survival formula incorrectly specified.")
      }
      else {
          if (sum(is.element(yvar.names, names(data))) == 4) {
              ## Time dependent covariates is in effect.
              subj.names <- yvar.names[1]
              yvar.names <- yvar.names[-1]
          }
      }
      family <- "surv"
      ytry <- 0
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
  return (list(all.names=all.names, family=family, subj.names=subj.names, yvar.names=yvar.names, ytry=ytry,
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
