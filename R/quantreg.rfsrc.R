quantreg.rfsrc <- function(formula, data, object, newdata,
                           method = "local", splitrule = NULL,
                           prob = NULL, prob.epsilon = NULL,
                           oob = TRUE, fast = FALSE, maxn = 1e3, ...)
{
  ## we intialize the grow primary operation as false
  grow <- FALSE
  ## which method will be used?
  method <- match.arg(method, c("forest", "local", "gk", "GK", "G-K", "g-k"))
  ## we now allow other regression splitting rules to be used 
  ## splitrule <- "quantile.regr"
  ## --------------------------------------------------------------------------------
  ##
  ## grow mode
  ##
  ## if no object is present, then the user is requesting a quantile regression forest
  ##
  ## --------------------------------------------------------------------------------
  if (missing(object)) {
    ## the primary operation is grow 
    grow <- TRUE
    ## list of forest parameters
    rfnames <- get.rfnames(hidden = TRUE)
    ## restrict to allowed values
    rfnames <- rfnames[rfnames != "formula"            &
                       rfnames != "data"               &
                       rfnames != "splitrule"          ]
    ## get the permissible hidden options
    dots <- list(...)
    dots <- dots[names(dots) %in% rfnames]
    ## manually set key hidden options 
    dots$prob <- prob
    dots$prob.epsilon <- prob.epsilon
    ## set the default split rule
    if (is.null(splitrule)) {
      splitrule <- "la.quantile.regr"
    }
    ## quantile splitrule cannot be used for multivariate regression
    fmly <- parseFormula(formula, data)$family
    if (fmly == "regr+" || fmly == "class+" || fmly == "mix+") {
      if (splitrule == "la.quantile.regr" || splitrule == "la.quantile.regr") {
        splitrule <- NULL
      }
    }
    ## determine the grow interface - rfsrc or rfsrc.fast?
    if (!fast) {
      rfsrc.grow <- "rfsrc"
    }
    else {
      rfsrc.grow <- "rfsrc.fast"
    }
    ## if the user wants method=forest, set GK to false and request forest weights
    if (method == "forest") {
      dots$gk.quantile <- FALSE
      dots$forest.wt <- if (oob) "oob" else TRUE
    }
    ## if user wants local method, set GK to false
    else if (method == "local") {
      dots$gk.quantile <- FALSE
    }
    ## user is requesting GK method
    else {
      dots$gk.quantile <- TRUE
    }
    ## set the hidden quantile.regr flag as TRUE
    dots$quantile.regr <- TRUE
    ## grow a quantile regression forest
    object <- do.call(rfsrc.grow, c(list(formula = formula, data = data, splitrule = splitrule), dots))
    ## save the grow yvar - this will be crucial downstream for prediction
    object$yvar.grow <- object$yvar
    ## ----------------------------------------------------------------------
    ##
    ## save the grow residual - this will be crucial downstream for prediction
    ##
    ## ----------------------------------------------------------------------
    ## pull the target outcome names
    if (object$family == "regr") {
      ynames <- object$yvar.names
    }
    else {
      ynames <- names(object$regrOutput)
    }
    ## cycle over ynames, extract residual
    if (ncol(cbind(object$yvar.grow)) > 1) {
      object$res.grow <- do.call(cbind, lapply(ynames, function(yn) {
        if (oob && !is.null(object$regrOutput[[yn]]$predicted.oob)) {
          object$yvar[, yn] - object$regrOutput[[yn]]$predicted.oob
        }
        else {
          object$yvar[, yn] - object$regrOutput[[yn]]$predicted
        }
      }))
      colnames(object$res.grow) <- ynames
    }
    else {##only one yvar - vectorize it
      if (oob && !is.null(object$predicted.oob)) {
        object$res.grow <- object$yvar - object$predicted.oob
      }
        else {
        object$res.grow <- object$yvar - object$predicted
      }
    }
  }
  ## user has passed in an object - is it a quantile object?
  else {
    if (sum(grepl("quantreg", class(object))) == 0) { 
      stop("object must be a quantreg object")
    }
  }
  ## -----------------------------------------------------------------------
  ##
  ## forest object must contain regression outcomes
  ## this includes multivariate regression and mixed multivariate regression
  ##
  ## -----------------------------------------------------------------------
  if (!(object$family == "regr" | !is.null(object$regrOutput))) {
    stop("this function requires at least one outcome that is continuous\n")
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
  ##  prediction mode
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
    ## grow mode is in effect, but new data is present --> oob cannot be TRUE
    if (!missing(newdata)) {
      oob <- FALSE
    }
    ## forest method -> set GK to false and request forest weights
    if (method == "forest") {
      dots$gk.quantile <- FALSE
      dots$forest.wt <- TRUE
      if (missing(newdata) && oob) {
        dots$forest.wt <- "oob"
      }
    }
    ## local method -> set GK to false and request forest weights
    else if (method == "local") {
      dots$gk.quantile <- FALSE
    }
    ## user is requesting GK method
    else {
      dots$gk.quantile <- TRUE
    }
    ## save grow information before restore/predict over-writes it
    yvar.grow <- object$yvar.grow
    res.grow <- object$res.grow
    ## restore call on grow data 
    if (missing(newdata)) {
      object <- do.call(predict, c(list(object = object), dots))
    }
    ## predict call with new data
    else {
      object <- do.call(predict, c(list(object = object, newdata = newdata), dots))
    }
    ## restore grow information
    object$yvar.grow <- yvar.grow
    object$res.grow <- res.grow
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
  rO <- lapply(ynames, function(yn) {
    ## ------------------------------------
    ##
    ## GK method - also clean up the object   
    ##
    ## ------------------------------------
    if (!(method == "forest" || method == "local")) {
      if (ncol(cbind(object$yvar.grow)) > 1) {
        y <- object$yvar.grow[, yn]
        if (oob && !is.null(object$regrOutput[[yn]]$quantile.oob)) {
        quant <- object$regrOutput[[yn]]$quantile.oob
        }
        else {
          quant <- object$regrOutput[[yn]]$quantile
        }
        object$regrOutput[[yn]]$quantile.oob <<- object$regrOutput[[yn]]$quantile <<- NULL
      }
      else {
        y <- object$yvar.grow
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
      ## unique y values
      yunq <- sort(unique(y))
      ## cdf calculation
      ## trim atoms: useful for big data otherwise O(n^2)
      if (length(yunq) > maxn) {
        yunq <- yunq[unique(round(seq(1, length(yunq), length.out = maxn)))]
      }
      cdf <- t(apply(quant.temp, 1, function(q) {prob.temp[sIndex(q, yunq)]}))
      ## density
      density <- t(apply(cbind(0, cdf), 1, diff))
    }
    ## ------------------------------------
    ##
    ## forest weight or local method
    ##
    ## ------------------------------------
    else {
      ## pull grow y, grow residual, and yhat (either grow/predicted)
      if (ncol(cbind(object$yvar.grow)) > 1) {
        y <- object$yvar.grow[, yn]
        res <- object$res.grow[, yn]
        if (oob && !is.null(object$regrOutput[[yn]]$predicted.oob)) {
          yhat <- object$regrOutput[[yn]]$predicted.oob
        }
        else {
          yhat <- object$regrOutput[[yn]]$predicted
        }
      }
      else {
        y <- object$yvar.grow
        res <- object$res.grow
         if (oob && !is.null(object$predicted.oob)) {
          yhat <- object$predicted.oob
        }
        else {
          yhat <- object$predicted
        }
      }
      ## extract unique y values
      yunq <- sort(unique(y))
      ## forest weight method, locally adjusts cdf using forest weights
      if (method == "forest") {
        cdf <- do.call(rbind, mclapply(1:nrow(object$forest.wt), function(i) {
          ind.matx <- do.call(rbind, lapply(yunq, function(yy) {res <= (yy - yhat[i])}))
          c(ind.matx %*% object$forest.wt[i, ])
        }))
      }
      else {
        ## locally adjusted cdf 
        cdf <- do.call(rbind, mclapply(1:length(yhat), function(i) {
          sapply(yunq, function(yy) {mean(res <= (yy - yhat[i]))})
        }))
      }
      ## quantiles
      quant <- t(apply(cdf, 1, function(pr) {
        c(min(yunq, na.rm = TRUE), yunq)[1 + sIndex(pr, prob)]
      }))
      ## density
      density <- t(apply(cbind(0, cdf), 1, diff))
    }
    ## ------------------------------------
    ##
    ## ends loop for specific yvar
    ## return various quantities as list
    ##
    ## ------------------------------------
    list(quantiles = quant,
         prob = prob,
         cdf = cdf,
         density = density,
         yunq = yunq)
  })
  ## ------------------------------------
  ##
  ## finished -- return final goodies
  ##
  ## ------------------------------------
  if (object$family == "regr") {
    rO <- rO[[1]]
  }
  else {
    names(rO) <- ynames
  }
  object$quantreg <- rO
  class(object)[4] <- "quantreg"
  object
}
quantreg <- quantreg.rfsrc
