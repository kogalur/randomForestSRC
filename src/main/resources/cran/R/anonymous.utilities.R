rfsrc.anonymous <- function(formula, data, forest = TRUE, ...)
{
  ## --------------------------------------------------------------
  ##   
  ##   preliminary processing
  ##
  ## --------------------------------------------------------------
  ## none yet
  ##--------------------------------------------------------------
  ##
  ## extract additional options specified by user
  ## we lock this down to allowed types
  ##
  ##--------------------------------------------------------------
  ## list of forest parameters
  rfnames <- names(formals(rfsrc))
  ## add key hidden parameters
  rfnames <- c(rfnames, "rfq", "perf.type", "gk.quantile", "prob", "prob.epsilon", "vtry")
   
  ## restrict to allowed values
  rfnames <- rfnames[rfnames != "forest"]
  ## get the permissible hidden options
  dots <- list(...)
  ## add formula if present
  dots <- dots[names(dots) %in% rfnames]
  if (!missing(formula)) {
    dots$formula <- formula
  }
  ## manually over-ride these key hidden options
  dots$terminal.qualts <- TRUE
  ## UBK UBK uncomment this when you are ready UBK UBK
  #dots$terminal.quants <- TRUE
  ##--------------------------------------------------------------
  ##
  ## make the grow call
  ##
  ##--------------------------------------------------------------
  retO <- do.call("rfsrc", c(list(data = data, forest = forest), dots))
  ##--------------------------------------------------------------
  ##
  ## strip out the training data
  ##
  ##--------------------------------------------------------------
  retO$xvar <- retO$xvar.names <- retO$forest$xvar <- NULL
  ## add special class distinction
  class(retO) <- c(class(retO), "anonymous")
  retO
}
