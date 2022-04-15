imbalanced.rfsrc <- function(formula, data, ntree = 3000,
                 method = c("rfq", "brf", "standard"),
                 block.size = NULL, perf.type = NULL, fast = FALSE,
                 ratio = NULL,  ...)
{
  ## preliminary checks: all are fatal
  ## parse the formula to ensure this is a two-class problem
  formulaPrelim <- parseFormula(formula, data)
  if (formulaPrelim$family != "class") {
    stop("this function only applies to classification problems")
  }
  yvar <- data[, formulaPrelim$yvar.names]
  if (length(levels(yvar)) != 2) {
    stop("this function only applies to two-class problems")
  }
  ## check that method is set correctly
  method <- match.arg(method, c("rfq", "brf", "standard"))
  rfq.flag <- NULL
  if (method == "rfq") {
    rfq.flag <- TRUE
  }
  ## set default performance
  if (is.null(perf.type)) {
    if (method == "brf" || method == "rfq") {
      perf.type <- "gmean"
    }
    else {
      perf.type <- "default"##equivalent to misclass
    }
  }
  ## check performance type is properly set
  perf.type <- match.arg(perf.type,
        c("none", "default", "standard", "misclass", "brier", "gmean", "g.mean"))
  if (perf.type == "g.mean") {##legacy
    perf.type <- "gmean"
  }
  ##-----------------------------------------
  ##
  ## rfsrc call - depends on what's requested
  ##
  ##-----------------------------------------
  ##----------------------------------------
  ## brf method
  ##----------------------------------------
  if (method == "brf") {
    ## TBD2 currently cannot handle missing values
    data <- na.omit(data)
    yvar <- data[, formulaPrelim$yvar.names]
    ## for legacy reasons we maintain swr here
    ## TBD2 consider swor
    o <- rfsrc(formula, data, ntree = ntree, 
               perf.type = perf.type,
               block.size = block.size,
               case.wt = make.wt(yvar),
               sampsize = make.size(yvar),
               samptype = "swr", ...)
  }
  ##----------------------------------------
  ## standard and rfq method
  ##----------------------------------------
  else {
    ##-------------------------------------------------
    ## determine the grow interface - rfsrc or rfsrc.fast?
    ##-------------------------------------------------
    if (!fast) {
      rfsrc.grow <- "rfsrc"
    }
    else {
      rfsrc.grow <- "rfsrc.fast"
    }
    ##-------------------------------------------------
    ## acquire the user specified additional options
    ##-------------------------------------------------
    dots <- list(...)
    ##-------------------------------------------------
    ## undersampling of the majority class if requested
    ##-------------------------------------------------
    if (!is.null(ratio)) {
      ## TBD2 currently cannot handle missing values
      data <- na.omit(data)
      yvar <- data[, formulaPrelim$yvar.names]
      samp <- make.imbalanced.sample(ntree = ntree, ratio = ratio, y = yvar)
      dots$bootstrap <- dots$samp <- NULL
      o <- do.call(rfsrc.grow, c(list(formula = formula, data = data, ntree = ntree,
                                 rfq = rfq.flag, perf.type = perf.type, block.size = block.size,
                                 samp = samp, bootstrap = "by.user"), dots))
    }
    ##-----------------------------------------------------------
    ## proceed to standard and rfq analysis without undersampling
    ## this is the default scenario
    ## allow fast random forests if requested
    ##-----------------------------------------------------------
    else {
      o <- do.call(rfsrc.grow, c(list(formula = formula, data = data, ntree = ntree,
                                 rfq = rfq.flag, perf.type = perf.type, block.size = block.size), dots))    
    }
  }
  ## return the object
  o
}
imbalanced <- imbalanced.rfsrc
