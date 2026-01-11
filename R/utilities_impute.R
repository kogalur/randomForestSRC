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
  assign.impute.mean(data, get.impute.mean(data))
}
