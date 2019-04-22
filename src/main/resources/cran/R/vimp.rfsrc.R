vimp.rfsrc <- function(object,
                       xvar.names,
                       m.target = NULL,
                       importance = c("permute", "random", "anti"),
                       block.size = 10,
                       joint = FALSE,
                       subset,
                       seed = NULL,
                       do.trace = FALSE,
                       ...)
{
  ## incoming parameter checks - all are fatal
  if (missing(object)) {
    stop("object is missing")
  }
  if (object$family == "unsupv") {
    stop("vimp does not apply to unsupervised forests: consider using max.subtree and var.select")
  }
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2) {
    stop("This function only works for objects of class `(rfsrc, grow)' or '(rfsrc, forest)'")
  }
  ## process the importance specification
  if (!is.logical(joint)) {
    stop("joint must be a logical value")
  }
  importance <- importance[1]
  if (joint & importance != "none") {
    i.str <- unlist(strsplit(importance, "\\."))
    if (length(i.str) == 1) {
      importance <- paste(i.str[1], ".joint", sep = "")
    }
      else if (length(i.str) == 2) {
        importance <- paste(i.str[1], ".joint.", i.str[2], sep = "")
      }
  }
  importance <- match.arg(as.character(importance),
                          c("permute", "random", "anti",
                            "permute.joint", "random.joint", "anti.joint"))
  ## grow objects under non-standard bootstrapping are devoid of performance values
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.null(object$forest)) {
      stop("The forest is empty.  Re-run rfsrc (grow) call with forest=TRUE")
    }
      else {
        bootstrap <- object$forest$bootstrap
      }
  }
    else {
      bootstrap <- object$bootstrap
    }
  if (bootstrap == "none" || bootstrap == "by.node") {
    stop("grow objects under non-standard bootstrapping are devoid of performance values")
  }
  ## process the subsetted index 
  ## assumes the entire data set is to be used if not specified
  if (missing(subset)) {
    subset <- NULL
  }
  else {
    ## convert the user specified subset into a usable form
    if (is.logical(subset)) {
        subset <- which(subset)
      }
    subset <- unique(subset[subset >= 1 & subset <= nrow(object$xvar)])
    if (length(subset) == 0) {
        stop("'subset' not set properly")
    }
  }
  ## make the call to generic predict
  result <- generic.predict.rfsrc(object,
                                  m.target = m.target,
                                  importance = importance,
                                  block.size = block.size,
                                  importance.xvar = xvar.names,
                                  seed = seed,
                                  do.trace = do.trace,
                                  membership = FALSE,
                                  subset = subset,
                                  ...)
  return(result)
}
vimp <- vimp.rfsrc
