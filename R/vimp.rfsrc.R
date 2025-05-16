vimp.rfsrc <- function(object,
                       xvar.names,
                       importance = c("anti", "permute", "random"),
                       block.size = 10,
                       joint = FALSE,
                       seed = NULL,
                       do.trace = FALSE,
                       ...)
{
  ## incoming parameter checks - all are fatal
  if (missing(object)) {
    stop("object is missing")
  }
  if (object$family == "unsupv") {
    stop("vimp does not apply to unsupervised forests: consider using max.subtree or varpro")
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
      c(TRUE, "anti", "permute", "random", "anti.joint", "permute.joint", "random.joint"))
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
  ## legacy m.target
  dots <- list(...)
  m.target <- dots$m.target
  dots$m.target <- NULL
  ## generic predict call
  args <- c(list(
    object = object,
    m.target = m.target,
    importance = importance,
    block.size = block.size,
    seed = seed,
    do.trace = do.trace,
    membership = FALSE
    ), dots)
  ## if xvar.names is not missing  
  if (!missing(xvar.names)) {
    args$importance.xvar <- xvar.names
  }
  return(do.call("generic.predict.rfsrc", args))
}
vimp <- vimp.rfsrc
