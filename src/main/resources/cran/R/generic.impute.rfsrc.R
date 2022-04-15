generic.impute.rfsrc <- function(data,
                                 ntree = 250,
                                 nodesize = NULL,
                                 nsplit = 1,
                                 nimpute = 1,
                                 fast = FALSE,
                                 ...)
{
  ## save the row and column names: later we will check if any rows or columns
  ## were deleted as part of the missing data preprocessing
  c.names <- colnames(data)
  r.names <- rownames(data)
  ## acquire the permissible hidden options
  dots <- list(...)
  dots$na.action <- dots$impute.only <- dots$forest <- NULL
  ## rfsrc grow call 
  if (!fast) {
    object <- do.call("rfsrc",
                    c(list(data = data,
                           ntree = ntree,
                           nodesize = nodesize,
                           nsplit = nsplit,
                           nimpute = nimpute,
                           na.action = "na.impute",
                           impute.only = TRUE), dots))
  }
  else {## user has requested the fast forest interface
    object <- do.call("rfsrc.fast",
                    c(list(data = data,
                           ntree = ntree,
                           nodesize = nodesize,
                           nsplit = nsplit,
                           nimpute = nimpute,
                           na.action = "na.impute",
                           impute.only = TRUE), dots))
  }
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
