cindex <- function (time,
                    censoring,
                    predicted,
                    do.trace = FALSE)
{
  size <- length(time)
  if (size != length(time) |
      size != length(censoring) |
      size != length(predicted)) {
    stop("time, censoring, and predicted must have the same length")
  }
  miss <- is.na(time) | is.na(censoring) | is.na(predicted)
  nmiss <- sum(miss)
  if (nmiss == size) {
    stop("no valid pairs found, too much missing data")
  }
  ## Flag missing members so we can exclude them in the pairs.
  denom <- sapply(miss, function(x) if (x) 0 else 1)
  nativeOutput <- .Call("rfsrcCIndex",
                        as.integer(do.trace),
                        as.integer(size),
                        as.double(time),
                        as.double(censoring),
                        as.double(predicted),
                        as.integer(denom))
  ## check for error return condition in the native code
  if (is.null(nativeOutput)) {
    stop("An error has occurred in rfsrcCIndex.  Please turn trace on for further analysis.")
  }
  return (nativeOutput$err)
}
