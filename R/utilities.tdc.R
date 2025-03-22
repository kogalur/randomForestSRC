####################################################################
##
## Time Dependent Covariates Related Functions
##
####################################################################
get.duplicated <- function(x) {
  c(apply(x, 2, function(xx) {
    1 * !all(xx == xx[1])
  }))
}
get.tdc.cov <- function(dta) {
  ## extract covariates
  x <- dta[, !(colnames(dta) %in% c("id", "start", "stop", "event")), drop = FALSE]
  ## for looping across the sorted id
  id <- dta$id
  id.unq <- sort(unique(id))
  ## iterate by id, check if i's data is duplicated - should work for real and character
  tdcm <- do.call(rbind, mclapply(id.unq, function(ii) {
    get.duplicated(x[id == ii,, drop = FALSE])
  }))
  ## return 0 for a covariate if x is duplicated for each id value
  ## value of 1 represents a tdc covariate
  c(apply(tdcm, 2, function(xx) {
    1 * !all(xx == 0)
  }))  
}
get.tdc.subj.time <- function(dta) {
  ## first check if there are time dependent covariates
  if (sum(get.tdc.cov(dta)) == 0) {
    return(rep(0, nrow(dta)))
  }
  ## extract covariates
  x <- dta[, !(colnames(dta) %in% c("id", "start", "stop", "event")), drop = FALSE]
  ## for looping across the id (id should *not* be sorted)
  id <- dta$id
  id.unq <- unique(id)
  ## iterate by id, check if i's data is duplicated - should work for real and character
  tdcm <- do.call(rbind, mclapply(id.unq, function(ii) {
    get.duplicated(x[id == ii,, drop = FALSE])
  }))
  ## return 1 if subject i has any tdc covariate
  c(apply(tdcm, 1, function(xx) {
    1 * any(xx != 0)
  }))
}
