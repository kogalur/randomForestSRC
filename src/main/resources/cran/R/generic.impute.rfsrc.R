generic.impute.rfsrc <- function(formula,
                                 data,
                                 ntree = 250,
                                 nimpute = 1,
                                 bootstrap = c("by.root", "by.node", "none"),
                                 mtry = NULL,
                                 nodesize = NULL,
                                 splitrule = NULL,
                                 nsplit = 1,
                                 na.action = c("na.impute"),
                                 xvar.wt = NULL,
                                 seed = NULL,
                                 do.trace = FALSE)
{
  bootstrap <- match.arg(bootstrap, c("by.root", "by.node", "none"))
  importance <- "none"
  na.action <- match.arg(na.action, c("na.impute"))
  forest <- FALSE
  proximity <- FALSE
  var.used <- FALSE
  split.depth <- FALSE
  impute.only <- TRUE
  membership <- FALSE
  c.names <- colnames(data)
  r.names <- rownames(data)
  object <- rfsrc(formula = formula,
                  data = data,
                  ntree = ntree,
                  bootstrap = bootstrap,
                  mtry = mtry,
                  nodesize = nodesize,
                  splitrule = splitrule,
                  nsplit = nsplit,
                  nimpute = nimpute,
                  xvar.wt = xvar.wt,
                  seed = seed,
                  do.trace = do.trace,
                  importance = importance,
                  na.action = na.action,
                  forest = forest,
                  proximity = proximity,
                  var.used = var.used,
                  split.depth = split.depth,
                  membership = membership,
                  impute.only = impute.only)
  if (is.null(object)) {
    return(NULL)
  }
  rm(data)
  if (is.data.frame(object)) {
    return(invisible(list(data = object, missing = row.col.deleted(object, r.names, c.names)))) 
  }
  if(is.null(object$yvar.names)) {
    imputed.result <- object$xvar
  }
    else {
      imputed.result <- cbind(object$yvar, object$xvar)
    }
  colnames(imputed.result) <- c(object$yvar.names, object$xvar.names)
  if (nimpute == 1) {
    imputed.result[object$imputed.indv, ] <- object$imputed.data
  }
  rm(object)
  invisible(list(data = imputed.result, missing = row.col.deleted(imputed.result, r.names, c.names)))
}
