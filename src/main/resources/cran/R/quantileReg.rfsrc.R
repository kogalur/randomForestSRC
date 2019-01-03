quantileReg.rfsrc <- function(formula, data, object, newdata, method = "forest",
                              prob = NULL, prob.epsilon = NULL,
                              oob = TRUE, fast = FALSE, maxn = 1e3, ...)
{
  ## we intialize the grow primary operation as false
  grow <- FALSE
  ## which method will be used?
  method <- match.arg(method, c("forest", "gk", "GK", "G-K", "g-k"))
  ## we only allow quantile regression splitting 
  ## otherwise prob and prob.epsilon will not be processed by global.prob.assign
  splitrule <- "quantile.regr"
  ## --------------------------------------------------------------------------------
  ##
  ## if no object is present, then the user is requesting a quantile regression forest
  ##
  ## --------------------------------------------------------------------------------
  if (missing(object)) {
    ## the primary operation is grow 
    grow <- TRUE
    ## pull unnamed parameters to pass to rfsrc
    dots <- list(...)
    dots$formula <- dots$data <- NULL
    ## if prob and prob.epsilon are NULL, leave them alone
    dots$prob <- prob
    dots$prob.epsilon <- prob.epsilon
    ## TBD TBD quantile splitrule fails in the case of multivariate regression
    fmly <- parseFormula(formula, data)$family
    if (fmly == "regr+" || fmly == "class+" || fmly == "mix+") {
      splitrule <- NULL
    }
    ## determine the grow interface - rfsrc or rfsrcFast?
    if (!fast) {
      rfsrc.grow <- "rfsrc"
    }
    else {
      rfsrc.grow <- "rfsrcFast"
    }
    ## if the user wants method=forest, set GK to false and request forest weights
    if (method == "forest") {
      dots$gk.quantile <- FALSE
      dots$forest.wt <- if (oob) "oob" else TRUE
    }
    else {
      dots$gk.quantile <- TRUE
    }
    ## grow a quantile regression forest
    object <- do.call(rfsrc.grow, c(list(formula = formula, data = data, splitrule = splitrule), dots))
  }
  ## user has passed in an object - is it a quantile object?
  else {
    if (sum(grepl("quantileReg", class(object))) == 0) { 
      stop("object must be a quantileReg object")
    }
  }
  ## -----------------------------------------------------------------------
  ##
  ## forest object must contain regression outcomes
  ## this includes multivariate regression and mixed multivariate regression
  ##
  ## -----------------------------------------------------------------------
  if (!(object$family == "regr" | !is.null(object$regrOutput))) {
    stop("this function only applies to regression settings\n")
  }
  ## ---------------------------------------
  ##
  ##  save grow prob for later calculations
  ##  can be over-written in predict mode
  ##
  ## -----------------------------------------
  prob.grow <- object$forest$prob
  ## ---------------------------------------
  ##
  ##  we are in prediction mode
  ##
  ## -----------------------------------------
  if (!grow || !missing(newdata)) {
    ## the user can over-ride grow options in predict mode
    ## if prob and prob.epsilon are non-NULL, over-ride grow values
    ## make sure to overwrite prob.grow 
    dots <- list(...)
    if (!is.null(prob)) {
      prob.grow <- dots$prob <- prob
    }
    else {
      dots$prob <- prob.grow
    }
    dots$prob.epsilon <- prob.epsilon
    ## forest method -> set GK to false and request forest weights
    if (method == "forest") {
      dots$gk.quantile <- FALSE
      dots$forest.wt <- TRUE
      if (missing(newdata) && oob) {
        dots$forest.wt <- "oob"
      }
    }
    else {
      dots$gk.quantile <- TRUE
    }
    ## predict call on grow data 
    if (missing(newdata)) {
      object <- do.call(predict, c(list(object = object), dots))
    }
    ## predict call with new data
    else {
      ## save grow yvar before predict over-writes it
      yvar <- object$yvar
      object <- do.call(predict, c(list(object = object, newdata = newdata), dots))
      ## overlay the grow yvar on the object
      object$yvar <- yvar
    }
  }
  ## -----------------------------------
  ##
  ## grow/test done: final value of prob
  ##
  ## -----------------------------------
  prob <- prob.grow
  ## -----------------------------------
  ##
  ## pull the target outcome names
  ##
  ## -----------------------------------
  if (object$family == "regr") {
    ynames <- object$yvar.names
  }
  else {
    ynames <- names(object$regrOutput)
  }
  ## -----------------------------------
  ##
  ## finally - extract quantile information
  ##
  ## -----------------------------------
  ## step function interpolation
  sIndex <- function(x1, x2) {sapply(1:length(x2), function(j) {sum(x1 <= x2[j])})}
  rO <- lapply(ynames, function(yn) {
    ## GK method - also clean up the object   
    if (method != "forest") {
      if (ncol(cbind(object$yvar)) > 1) {
        y <- object$yvar[, yn]
        if (oob && !is.null(object$regrOutput[[yn]]$quantile.oob)) {
        quant <- object$regrOutput[[yn]]$quantile.oob
        }
        else {
          quant <- object$regrOutput[[yn]]$quantile
        }
        object$regrOutput[[yn]]$quantile.oob <<- object$regrOutput[[yn]]$quantile <<- NULL
      }
      else {
        y <- object$yvar
        if (oob && !is.null(object$quantile.oob)) {
          quant <- object$quantile.oob
        }
        else {
          quant <- object$quantile
        }
        object$quantile.oob <<- object$quantile <<- NULL
      }
      ## append a minimum and maximum y value corresponding to prob = 0, 1
      quant.temp <- cbind(min(y, na.rm = TRUE) - 1, quant, max(y, na.rm = TRUE))
      prob.temp <- c(0, prob, 1)
      ## cdf - pull the atoms = unique y-values 
      yunq <- sort(unique(y))
      if (length(yunq) > maxn) {##trim atoms: useful for big data otherwise O(n^2)
        yunq <- yunq[unique(round(seq(1, length(yunq), length.out = maxn)))]
      }
      cdf <- t(apply(quant.temp, 1, function(q) {prob.temp[sIndex(q, yunq)]}))
    }
    ## forest weight method
    else {
      if (ncol(cbind(object$yvar)) > 1) {
        y <- object$yvar[, yn]
      }
      else {
        y <- object$yvar
      }
      ## cdf 
      yunq <- sort(unique(y))
      ind.matx <- do.call(rbind, mclapply(yunq, function(yy) {y <= yy}))
      cdf <- t(apply(object$forest.wt, 1, function(wt) {ind.matx %*% wt}))
      ## quantiles
      quant <- t(apply(cdf, 1, function(pr) {
        c(min(yunq, na.rm = TRUE), yunq)[1 + sIndex(pr, prob)]
      }))
    }
    ## density
    density <- t(apply(cbind(0, cdf), 1, diff))
    ## yhat
    yhat <- density %*% yunq
    ## return the object
    list(quantiles = quant,
         prob = prob,
         cdf = cdf,
         density = density,
         yhat = yhat,
         yunq = yunq)
  })
  ## return the goodies
  if (object$family == "regr") {
    rO <- rO[[1]]
  }
  else {
    names(rO) <- ynames
  }
  object$quantileReg <- rO
  class(object)[4] <- "quantileReg"
  object
}
quantileReg <- quantileReg.rfsrc
