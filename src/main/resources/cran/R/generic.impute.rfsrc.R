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
  ## set parameters accordingly
  bootstrap <- match.arg(bootstrap, c("by.root", "by.node", "none"))
  importance <- "none"
  
  na.action <- match.arg(na.action, c("na.impute"))
  
   
  forest <- FALSE
  proximity <- FALSE
     
  var.used <- FALSE
  split.depth <- FALSE
  impute.only <- TRUE
  membership <- FALSE
  ## save the row and column names: later we will check if any rows or columns
  ## were deleted as part of the missing data preprocessing
  c.names <- colnames(data)
  r.names <- rownames(data)
  ## rfsrc grow call
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
  ## confirm that no error has occured
  if (is.null(object)) {
    return(NULL)
  }
  ## the data is no longer needed
  rm(data)
  ## if the return object is a data frame then imputation was not
  ## performed: for example, there was no missing data either before
  ## or after processing
  if (is.data.frame(object)) {
    return(invisible(list(data = object, missing = row.col.deleted(object, r.names, c.names)))) 
  }
  ## preliminary results of imputation
  if(is.null(object$yvar.names)) {
    imputed.result <- object$xvar
  }
    else {
      imputed.result <- cbind(object$yvar, object$xvar)
    }
  colnames(imputed.result) <- c(object$yvar.names, object$xvar.names)
  ## overlay the data (only necessary when nimpute = 1)
  if (nimpute == 1) {
    imputed.result[object$imputed.indv, ] <- object$imputed.data
  }
  ## the object is no longer required
  rm(object)
  ## return the goodies
  invisible(list(data = imputed.result, missing = row.col.deleted(imputed.result, r.names, c.names)))
}
