rfsrc.anonymous <- function(formula, data, forest = TRUE, ...)
{
  ## --------------------------------------------------------------
  ##   
  ##   preliminary processing
  ##
  ## --------------------------------------------------------------
  #if (any(is.na(data))) {
    #stop("missing values not allowed in anonymous mode")
  #}
  ##--------------------------------------------------------------
  ##
  ## extract additional options specified by user
  ## we lock this down to allowed types
  ##
  ##--------------------------------------------------------------
  ## list of forest parameters
  rfnames <- get.rfnames(hidden = TRUE)
   
  ## restrict to allowed values
  rfnames <- rfnames[rfnames != "data" & rfnames != "forest"]
  ## get the permissible hidden options
  ## add formula if present
  dots <- list(...)
  dots <- dots[names(dots) %in% rfnames]
  if (!missing(formula)) {
    dots$formula <- formula
  }
  ## manually over-ride key hidden options
  dots$terminal.qualts <- TRUE
  dots$terminal.quants <- TRUE
  ##--------------------------------------------------------------
  ##
  ## make the grow call
  ##
  ##--------------------------------------------------------------
  retO <- do.call("rfsrc", c(list(data = data, forest = forest, na.action = "na.omit"), dots))
  ##--------------------------------------------------------------
  ##
  ## save impute mean 
  ##
  ##--------------------------------------------------------------
  retO$forest$impute.mean <- get.impute.mean(data)
  ##--------------------------------------------------------------
  ##
  ## strip out the training data
  ##
  ##--------------------------------------------------------------
  retO$xvar <- retO$forest$xvar <- NULL
  ##--------------------------------------------------------------
  ##
  ## add special class distinction --> for both the grow object AND forest
  ##
  ##--------------------------------------------------------------
  class(retO) <- c(class(retO), "anonymous")
  class(retO$forest) <- c(class(retO$forest), "anonymous")
  ##--------------------------------------------------------------
  ##
  ## return the anonymized object
  ##
  ##--------------------------------------------------------------
  retO
}
