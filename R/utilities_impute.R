assign.impute.mean <- function(data, impute.mean) {
  d <- data.frame(lapply(colnames(data), function(xnms) {
    x <- data[, xnms]
    is.na.x <- is.na(x)
    if (any(is.na.x)) {
      x[is.na.x] <- impute.mean[[xnms]]
    }
    x
  }), stringsAsFactors = TRUE)
  colnames(d) <- colnames(data)
  d
}
get.impute.mean <- function(data) {
  imean <- mclapply(data, function(x) {
    if (all(is.na(x))) {
      NA
    }
    else {
      if (is.factor(x)) {
        x.table <- table(x)
        names(x.table)[which.max(x.table)]
      }
      else {
        mean(x, na.rm = TRUE)
      }
    }
  })
  names(imean) <- colnames(data)
  imean
}
get.na.roughfix <- function(data) {
  assign.impute.mean(data, get.impute.mean(data))
}
