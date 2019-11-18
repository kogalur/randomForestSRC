rfsrc.fast <- function(formula, data,
                       ntree = 500,
                       nsplit = 10,
                       bootstrap = "by.root",
                       ensemble = "oob",
                       sampsize = function(x){min(x * .632, max(150, x ^ (3/4)))},
                       samptype = "swor",
                       samp = NULL,
                       ntime = 50,
                       forest = FALSE,
                       terminal.qualts = FALSE,
                       ...)
{
  ## --------------------------------------------------------------
  ##   
  ##   preliminary processing
  ##
  ## --------------------------------------------------------------
  ## verify key options
  if (!is.function(sampsize) && !is.numeric(sampsize)) {
    stop("sampsize must be a function or number specifying size of subsampled data")
  }
  ##if (is.function(sampsize)) {
  ##  sampsize <- sampsize(nrow(data))
  ##}
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
  rfnames <- rfnames[rfnames != "ntree"              &
                     rfnames != "nsplit"             &
                     rfnames != "bootstrap"          &
                     rfnames != "ensemble"           &
                     rfnames != "sampsize"           &
                     rfnames != "samptype"           &
                     rfnames != "ntime"              &
                     rfnames != "forest"             &
                     rfnames != "terminal.qualts"    ]
  ## get the permissible hidden options
  ## add formula if present
  dots <- list(...)
  dots <- dots[names(dots) %in% rfnames]
  if (!missing(formula)) {
    dots$formula <- formula
  }
  ## set bootstrap accordingly if the user has provided their own sampling
  ## ntree and sampsize are handled in rfsrc
  if (!is.null(samp)) {
    bootstrap <- "by.user"
  }
  ##--------------------------------------------------------------
  ##
  ## make the grow call and return the object
  ##
  ##--------------------------------------------------------------
  return(do.call("rfsrc",
                 c(list(data = data,
                 ntree = ntree,
                 nsplit = nsplit,
                 bootstrap = bootstrap,
                 ensemble = ensemble,
                 sampsize = sampsize,
                 samptype = samptype,
                 samp = samp,
                 ntime = ntime,
                 terminal.qualts = terminal.qualts,
                 forest = forest), dots)))
}
