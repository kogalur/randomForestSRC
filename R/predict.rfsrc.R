predict.rfsrc <-  function(object,
   newdata,
   importance = c(FALSE, TRUE, "none", "anti", "permute", "random"),
   get.tree = NULL,
   block.size = if (any(is.element(as.character(importance), c("none", "FALSE")))) NULL else 10,
   na.action = c("na.omit", "na.impute", "na.random"),
   outcome = c("train", "test"),
   perf.type = NULL,
   proximity = FALSE,
   forest.wt = FALSE,
   ptn.count = 0,
   distance = FALSE,
   var.used = c(FALSE, "all.trees", "by.tree"),
   split.depth = c(FALSE, "all.trees", "by.tree"),
   case.depth = FALSE,
   seed = NULL,
   do.trace = FALSE,
   membership = FALSE,
   marginal.xvar = NULL,
   ...)
{
  dots <- list(...)
  m.target <- dots$m.target
  dots$m.target <- NULL
  importance.xvar <- dots$importance.xvar
  dots$importance.xvar <- NULL
  args <- c(list(
    object = object,
    m.target = m.target,
    importance = importance,
    get.tree = get.tree,
    block.size = block.size,
    importance.xvar = importance.xvar,
    na.action = na.action,
    outcome = outcome,
    perf.type = perf.type,
    proximity = proximity,
    forest.wt = forest.wt,
    ptn.count = ptn.count,
    distance = distance,
    var.used = var.used,
    split.depth = split.depth,
    case.depth = case.depth,
    seed = seed,
    do.trace = do.trace,
    membership = membership,
    marginal.xvar = marginal.xvar
  ), dots)
  ## if newdata is not missing
  if (!missing(newdata)) {
    args$newdata <- newdata
  }
  return(do.call("generic.predict.rfsrc", args))
}
