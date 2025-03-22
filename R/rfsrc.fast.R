rfsrc.fast <- function(formula, data,
                       ntree = 500,
                       nsplit = 10,
                       bootstrap = "by.root",
                       sampsize = function(x){min(x * .632, max(150, x ^ (3/4)))},
                       samptype = "swor",
                       samp = NULL,
                       ntime = 50,
                       forest = FALSE,
                       save.memory = TRUE,
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
  ##--------------------------------------------------------------
  ##
  ## extract additional options specified by user
  ## we lock this down to allowed types
  ##
  ##--------------------------------------------------------------
  ## list of forest parameters
  rfnames <- get.rfnames(hidden = TRUE)
   
  ## restrict to allowed values
  rfnames <- rfnames[rfnames != "data"               &
                     rfnames != "ntree"              &
                     rfnames != "nsplit"             &
                     rfnames != "bootstrap"          &
                     rfnames != "sampsize"           &
                     rfnames != "samptype"           &
                     rfnames != "ntime"              &
                     rfnames != "forest"             &
                     rfnames != "save.memory"        ]
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
  ## manually set key hidden options if forest is not requested
  if (!forest) {
    dots$terminal.qualts <- FALSE
    dots$terminal.quants <- FALSE
    dots$save.memory <- TRUE
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
                 sampsize = sampsize,
                 samptype = samptype,
                 samp = samp,
                 ntime = ntime,
                 forest = forest), dots)))
}
