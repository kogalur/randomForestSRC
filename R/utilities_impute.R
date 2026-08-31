.integer.support.info <- function(x,
                                  tol = sqrt(.Machine$double.eps),
                                  frac.cut = 1) {
  if (inherits(x, c("Date", "POSIXt", "difftime")) || is.factor(x)) {
    return(list(
      restore = FALSE,
      storage.integer = FALSE,
      integer.frac = NA_real_
    ))
  }
  if (is.integer(x)) {
    return(list(
      restore = TRUE,
      storage.integer = TRUE,
      integer.frac = 1
    ))
  }
  if (!is.numeric(x)) {
    return(list(
      restore = FALSE,
      storage.integer = FALSE,
      integer.frac = NA_real_
    ))
  }
  z <- as.numeric(x)
  z <- z[!is.na(z) & is.finite(z)]
  if (length(z) == 0L) {
    return(list(
      restore = FALSE,
      storage.integer = FALSE,
      integer.frac = NA_real_
    ))
  }
  scale <- pmax(1, abs(z))
  frac <- mean(abs(z - round(z)) <= tol * scale)
  list(
    restore = isTRUE(frac >= frac.cut),
    storage.integer = FALSE,
    integer.frac = frac
  )
}
.integer.support.map <- function(data,
                                 tol = sqrt(.Machine$double.eps),
                                 frac.cut = 1) {
  out <- lapply(data, .integer.support.info, tol = tol, frac.cut = frac.cut)
  names(out) <- names(data)
  out
}
.restore.integer.vector <- function(x, info, restore.integer = TRUE) {
  if (!isTRUE(restore.integer) || !isTRUE(info$restore)) {
    return(x)
  }
  z <- as.numeric(x)
  ok <- !is.na(z) & is.finite(z)
  z[ok] <- round(z[ok])
  if (isTRUE(info$storage.integer)) {
    out <- rep(NA_integer_, length(z))
    out[ok] <- as.integer(z[ok])
    return(out)
  }
  z
}
.restore.integer.data <- function(data, integer.support, restore.integer = TRUE) {
  if (!isTRUE(restore.integer) || is.null(integer.support) || length(integer.support) == 0L) {
    return(data)
  }
  for (nm in intersect(names(integer.support), names(data))) {
    info <- integer.support[[nm]]
    if (isTRUE(info$restore)) {
      data[[nm]] <- .restore.integer.vector(data[[nm]], info, restore.integer = TRUE)
    }
  }
  data
}
assign.impute.mean <- function(data, impute.mean) {
  cn <- colnames(data)
  p  <- length(cn)
  ## Build columns in a preallocated list (faster than lapply over names)
  out <- vector("list", p)
  names(out) <- cn
  for (j in seq_len(p)) {
    nm <- cn[j]
    x  <- data[[nm]]
    na_idx <- is.na(x)
    if (any(na_idx)) {
      x[na_idx] <- impute.mean[[nm]]
    }
    out[[j]] <- x
  }
  ## Preserve original behavior: coerce character columns to factor
  ## (original code used stringsAsFactors = TRUE)
  data.frame(out, stringsAsFactors = TRUE, check.names = FALSE)
}
get.impute.mean <- function(data) {
  cn <- colnames(data)
  p  <- length(data)
  ## Preallocate list output for speed
  imean <- vector("list", p)
  names(imean) <- cn
  if (p == 0L) return(imean)
  is_fac <- vapply(data, is.factor, logical(1))
  ## Fast path for numeric/integer/logical columns: use a matrix and colMeans
  is_numlike <- vapply(data, function(x) {
    (is.numeric(x) || is.integer(x) || is.logical(x)) && !is.factor(x)
  }, logical(1))
  num_idx <- which(is_numlike)
  fac_idx <- which(is_fac)
  if (length(num_idx)) {
    xm <- data.matrix(data[, num_idx, drop = FALSE])
    ##storage.mode(xm) <- "double" : 
    cnt <- colSums(!is.na(xm))
    mu  <- colMeans(xm, na.rm = TRUE)
    mu[cnt == 0] <- NA_real_
    for (k in seq_along(num_idx)) {
      imean[[num_idx[k]]] <- mu[k]
    }
  }
  ## Factor columns: compute the modal level (fast via tabulate)
  if (length(fac_idx)) {
    for (k in fac_idx) {
      x <- data[[k]]
      if (all(is.na(x))) {
        imean[[k]] <- NA
      } else {
        tab <- tabulate(x, nbins = length(levels(x)))
        imean[[k]] <- levels(x)[which.max(tab)]
      }
    }
  }
  ## Any remaining columns (e.g., character): match original logic
  ## (mean() -> NA with warning). We suppress warnings to avoid slowdown/noise.
  other_idx <- setdiff(seq_len(p), c(num_idx, fac_idx))
  if (length(other_idx)) {
    for (k in other_idx) {
      x <- data[[k]]
      if (all(is.na(x))) {
        imean[[k]] <- NA
      } else {
        imean[[k]] <- suppressWarnings(mean(x, na.rm = TRUE))
      }
    }
  }
  imean
}
get.na.roughfix <- function(data) {
  if (!is.data.frame(data)) data <- data.frame(data)
  assign.impute.mean(data, get.impute.mean(data))
}
