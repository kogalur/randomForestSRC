coerce.multivariate <- function(x, outcome.target) {
  ## Warning:  This functon assumes that get.univariate.target has been called first, to
  ## verify the coherency of the target.  This means that the target exists in the forest object, and that
  ## it contains outcome statistics.
  ## If this is a multivarate family, we coerce the object, based on outcome.target
  ## into a univaritate regression or classification object.
  x$univariate <- TRUE
  if (x$family == "regr+" | x$family == "class+" | x$family == "mix+") {
    ## coerce the mulitvariate object into a univariate object
    x.coerced <- unlist(list(x$classOutput, x$regrOutput), recursive = FALSE)[[outcome.target]]
    x$univariate <- FALSE
    x$yvar <- x$yvar[, outcome.target]
    ## test for factors - ordered factors are treated as factors!
    if (is.factor(x$yvar) || is.ordered(x$yvar)) {
      x$family <- "class"
    }
    else {
      x$family <- "regr"
    }
    ## make various assignments to the coerced object.
    x$predicted <- x.coerced$predicted
    x$predicted.oob <- x.coerced$predicted.oob
    x$class <- x.coerced$class
    x$class.oob <- x.coerced$class.oob
    x$err.rate <- x.coerced$err.rate
    x$err.block.rate <- x.coerced$err.block.rate
    x$importance <- x.coerced$importance
    x$yvar.names <- outcome.target
    x$cse.num <- x.coerced$cse.num
    x$cse.den <- x.coerced$cse.den
    x$csv.num <- x.coerced$csv.num
    x$csv.den <- x.coerced$csv.den
  }
  x$outcome.target <- outcome.target
  x
}
get.mv.cserror <- function(obj, standardize = FALSE) {
  ## acquire y variable names
  ynms <- obj$yvar.names
  ## does not apply to survival or competing risk
  if (obj$family == "surv" ||  obj$family == "surv-CR") {
    return(NULL)
  }
  ## loop over the yvar names acquring the csv if it is present
  err <- lapply(ynms, function(nn) {
    ## coerce to a univariate object for uniformity
    o.coerce <- coerce.multivariate(obj, nn)
    e.num <- o.coerce$cse.num
    e.den <- o.coerce$cse.den
    err <- NULL
    if (!is.null(e.num)) {
      ## classification
      if (!(o.coerce$family == "regr" || o.coerce$family == "regr+")) {
        err <- e.num / e.den
      }
      ## regression: standardization can be applied
      else {
        yvar <- 1
        if (standardize) {
          yvar <- var(o.coerce$yvar, na.rm = TRUE)
          if (is.na(yvar) || yvar == 0) {
            yvar <- 1
          }
        }
        err <- e.num / (e.den * yvar)
      }
    }
    ## return
    err
  })
  ## if the first entry is NULL make the entire list NULL
  if (is.null(err[[1]])) {
    err <- NULL
  }  
  ## add ynames - for multivariate objects only
  if (!is.null(err) && length(err)>1) {
    names(err) <- ynms
  }
  if (!is.null(err) && length(err)==1) {
    err <- err[[1]]
  }
  ### return the goodies
  err 
}
## get case-specific error - this is not the same as
## Loss(predicted.oob,y)
get.mv.cserror <- function(obj, standardize = FALSE) {
  ## acquire y variable names
  ynms <- obj$yvar.names
  ## does not apply to survival or competing risk
  if (obj$family == "surv" ||  obj$family == "surv-CR") {
    return(NULL)
  }
  ## loop over the yvar names acquring the csv if it is present
  err <- lapply(ynms, function(nn) {
    ## coerce to a univariate object for uniformity
    o.coerce <- coerce.multivariate(obj, nn)
    e.num <- o.coerce$cse.num
    e.den <- o.coerce$cse.den
    err <- NULL
    if (!is.null(e.num)) {
      ## classification
      if (!(o.coerce$family == "regr" || o.coerce$family == "regr+")) {
        err <- e.num / e.den
      }
      ## regression: standardization can be applied
      else {
        yvar <- 1
        if (standardize) {
          yvar <- var(o.coerce$yvar, na.rm = TRUE)
          if (is.na(yvar) || yvar == 0) {
            yvar <- 1
          }
        }
        err <- e.num / (e.den * yvar)
      }
    }
    ## return
    err
  })
  ## if the first entry is NULL make the entire list NULL
  if (is.null(err[[1]])) {
    err <- NULL
  }  
  ## add ynames - for multivariate objects only
  if (!is.null(err) && length(err)>1) {
    names(err) <- ynms
  }
  if (!is.null(err) && length(err)==1) {
    err <- err[[1]]
  }
  ### return the goodies
  err 
}
## get case-specific VIMP
get.mv.csvimp <- function(obj, standardize = FALSE) {
  ## acquire y variable names
  ynms <- obj$yvar.names
  ## does not apply to survival or competing risk
  if (obj$family == "surv" ||  obj$family == "surv-CR") {
    return(NULL)
  }
  ## loop over the yvar names acquring the csv if it is present
  vmp <- lapply(ynms, function(nn) {
    ## coerce to a univariate object for uniformity
    o.coerce <- coerce.multivariate(obj, nn)
    v.num <- o.coerce$csv.num
    v.den <- o.coerce$csv.den
    v <- NULL
    if (!is.null(v.num)) {
      ## classification
      if (o.coerce$family != "regr") {
        v <- v.num / v.den
      }
      ## regression: standardization can be applied
      else {
        yvar <- 1
        if (standardize) {
          yvar <- var(o.coerce$yvar, na.rm = TRUE)
          if (is.na(yvar) || yvar == 0) {
            yvar <- 1
          }
        }
        v <- v.num / (v.den * yvar)
      }
      ## pretty the vimp up
      if (is.null(dim(v))) {
        v <- cbind(v)
      }
      ## make allowance for joint vimp
      if (length(o.coerce$importance) == 1) {
        v <- v[, 1, drop = FALSE]
      }
      colnames(v) <- names(o.coerce$importance)
    }
    ## return
    v
  })
  ## if the first entry is NULL make the entire list NULL
  if (is.null(vmp[[1]])) {
    vmp <- NULL
  }  
  ## add ynames - for multivariate objects only
  if (!is.null(vmp) && length(vmp)>1) {
    names(vmp) <- ynms
  }
  if (!is.null(vmp) && length(vmp)==1) {
    vmp <- vmp[[1]]
  }
  ### return the goodies
  vmp 
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
