## Shared input contracts for the imputation entry points. Containers are
## normalized before column lookup; each column represents one variable.
.impute.data.frame <- function(data, arg = "data", check.finite = TRUE) {
  if (is.atomic(data) && is.null(dim(data))) {
    stop("'", arg, "' must be a tabular object, not a bare vector.",
         call. = FALSE)
  }
  if (!identical(class(data), "data.frame")) {
    data <- as.data.frame(data, optional = TRUE, stringsAsFactors = FALSE)
  }
  nm <- names(data)
  if (length(nm) != ncol(data) || anyNA(nm) ||
      any(!nzchar(trimws(nm))) || anyDuplicated(nm)) {
    stop("'", arg, "' must have unique, nonempty column names.", call. = FALSE)
  }
  bad <- vapply(data, function(x) !is.atomic(x) || !is.null(dim(x)), logical(1))
  if (any(bad)) {
    stop("Each column of '", arg, "' must be a vector. Problem columns: ",
         paste(nm[bad], collapse = ", "), call. = FALSE)
  }
  if (isTRUE(check.finite)) {
    bad <- vapply(data, function(x) {
      is.numeric(x) && any(!is.na(x) & !is.finite(x))
    }, logical(1))
    if (any(bad)) {
      stop("Numeric values in '", arg, "' must be finite or missing. Problem columns: ",
           paste(nm[bad], collapse = ", "), call. = FALSE)
    }
  }
  data
}
.impute.scalar <- function(x, arg, lower = 0, whole = FALSE,
                            upper = Inf, open.lower = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || (isTRUE(open.lower) && x == lower) || x > upper ||
      (isTRUE(whole) && x != floor(x))) {
    stop("'", arg, "' must be a finite ", if (whole) "integer" else "number",
         if (open.lower) " greater than " else " at least ", lower,
         if (is.finite(upper)) paste0(" and at most ", upper) else "",
         ".", call. = FALSE)
  }
  invisible(x)
}
.impute.flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("'", arg, "' must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(x)
}
.impute.named.options <- function(x, arg) {
  if (is.null(x)) return(list())
  if (!is.list(x) || is.data.frame(x) ||
      (length(x) > 0L && (is.null(names(x)) ||
       anyNA(names(x)) || any(!nzchar(trimws(names(x))))))) {
    stop("'", arg, "' must be a named list.", call. = FALSE)
  }
  if (anyDuplicated(names(x))) {
    warning("Duplicate ", arg, " entries; keeping the last occurrence of each name.",
            call. = FALSE)
    x <- x[!duplicated(names(x), fromLast = TRUE)]
  }
  x
}
## Validate structural forest options once, before target-wise error handling.
.impute.forest.options <- function(x) {
  for (nm in intersect(c("ntree", "nodesize", "nsplit", "nimpute"), names(x))) {
    if (is.null(x[[nm]])) next
    .impute.scalar(x[[nm]], nm, lower = if (nm == "nsplit") 0 else 1,
                   whole = TRUE, upper = .Machine$integer.max)
  }
  invisible(x)
}
## Resolve the imputation branch once for both the engine and its metadata.
## Keep a non-NULL always.use distinct from NULL, even when no names match:
## supplying it selects multivariate missForest under the existing contract.
.impute.method.spec <- function(mf.q = NULL, always.use = NULL,
                                variable.names = character(0)) {
  if (!is.null(mf.q)) {
    .impute.scalar(mf.q, "mf.q", lower = 0, open.lower = TRUE)
  }
  if (!is.null(always.use) && (!is.character(always.use) || anyNA(always.use))) {
    stop("'always.use' must be NULL or a character vector of variable names.",
         call. = FALSE)
  }
  mforest <- !is.null(mf.q)
  always.index <- if (is.null(always.use)) NULL else {
    which(variable.names %in% always.use)
  }
  univariate <- mforest && mf.q == 1 && is.null(always.index)
  list(
    mforest = mforest,
    univariate = univariate,
    always.use = always.index,
    name = if (!mforest) "on-the-fly imputation" else if (univariate) {
      "missForest"
    } else "multivariate missForest"
  )
}
## Build formulas from symbols, so names containing spaces, punctuation or
## backticks are retained without string interpolation or name repair.
.impute.formula <- function(response, predictors = NULL, env = parent.frame()) {
  lhs <- if (is.character(response)) {
    if (length(response) == 1L) as.name(response) else {
      as.call(c(list(as.name("Multivar")), lapply(response, as.name)))
    }
  } else response
  rhs <- if (is.null(predictors)) as.name(".") else {
    if (length(predictors) == 0L) 1 else {
      Reduce(function(a, b) call("+", a, b), lapply(predictors, as.name))
    }
  }
  stats::as.formula(call("~", lhs, rhs), env = env)
}
.impute.response.name <- function(existing.names) {
  nm <- ".impute.learn.response."
  while (nm %in% existing.names) nm <- paste0(nm, ".")
  nm
}
## Shared numeric-support restoration and mean/mode rough-fix utilities.
## Integer support is learned before filling, separately from R integer storage.
## Record whether a column is stored as integers or passes the numeric
## integer-distance criterion. The tolerance is an absolute grid distance.
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
  ## Use an absolute distance from the integer grid. A magnitude-relative
  ## tolerance can otherwise accept every fractional value in large columns.
  .impute.scalar(tol, "tol", lower = 0, upper = 0.25)
  .impute.scalar(frac.cut, "frac.cut", lower = 0, upper = 1)
  frac <- mean(abs(z - round(z)) <= tol)
  list(
    restore = isTRUE(frac >= frac.cut),
    storage.integer = FALSE,
    integer.frac = frac
  )
}
## Compute original support metadata independently for each column.
.integer.support.map <- function(data,
                                 tol = sqrt(.Machine$double.eps),
                                 frac.cut = 1) {
  out <- lapply(data, .integer.support.info, tol = tol, frac.cut = frac.cut)
  names(out) <- names(data)
  out
}
## Data-table restoration needs the mask recorded before values were generated.
## Without a mask, existing table entries are treated as observed. Column names
## align the mask after a forest returns responses before predictors.
.impute.generated.mask <- function(data, generated = NULL) {
  if (is.null(generated)) {
    return(matrix(FALSE, nrow = nrow(data), ncol = ncol(data),
                   dimnames = list(NULL, names(data))))
  }
  generated <- as.matrix(generated)
  if (!is.logical(generated) || length(dim(generated)) != 2L ||
      nrow(generated) != nrow(data) || anyNA(generated) ||
      is.null(colnames(generated)) || anyDuplicated(colnames(generated)) ||
      !all(names(data) %in% colnames(generated))) {
    stop("The generated-value mask must be logical, have one row per observation, ",
         "and contain each data column by name.", call. = FALSE)
  }
  generated[, names(data), drop = FALSE]
}
## Round only generated finite values. Prediction vectors contain generated
## values only and therefore use the default; table callers supply a cell mask.
## Observed fractions require double storage even under an R-integer schema.
.restore.integer.vector <- function(x, info, restore.integer = TRUE,
                                      generated = NULL) {
  if (!isTRUE(restore.integer) || !isTRUE(info$restore)) {
    return(x)
  }
  z <- as.numeric(x)
  if (is.null(generated)) {
    generated <- rep(TRUE, length(z))
  }
  if (!is.logical(generated) || !is.null(dim(generated)) ||
      length(generated) != length(z) || anyNA(generated)) {
    stop("The generated-value mask must be a logical vector with one entry per value.",
         call. = FALSE)
  }
  round.idx <- generated & is.finite(z)
  z[round.idx] <- round(z[round.idx])
  ## Integer storage is safe only when conversion preserves every nonmissing
  ## value exactly, including observations outside the generated-value mask.
  ## R's lowest integer code is reserved for NA.
  if (isTRUE(info$storage.integer) &&
      all(is.na(z) | (is.finite(z) & z == trunc(z) &
                     abs(z) <= .Machine$integer.max))) {
    return(as.integer(z))
  }
  z
}
## Restore generated cells by column name using the original missingness mask.
## The mask is unchanged across initialization, later sweeps and final output.
.restore.integer.data <- function(data, integer.support, restore.integer = TRUE,
                                   generated = NULL) {
  if (!isTRUE(restore.integer) || is.null(integer.support) || length(integer.support) == 0L) {
    return(data)
  }
  generated <- .impute.generated.mask(data, generated)
  for (nm in intersect(names(integer.support), names(data))) {
    info <- integer.support[[nm]]
    if (isTRUE(info$restore)) {
      data[[nm]] <- .restore.integer.vector(
        data[[nm]], info, restore.integer = TRUE, generated = generated[, nm]
      )
    }
  }
  data
}
## Fill NA cells using supplied means/modes, retaining row and column names.
## Character columns are converted to factors in the returned data frame.
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
  data.frame(out, row.names = rownames(data),
             stringsAsFactors = TRUE, check.names = FALSE)
}
## Compute numeric/logical means and factor modes; all-missing columns get NA.
## Factor modes use declared factor-level order to break frequency ties.
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
  ## For other types use their mean() method, suppressing conversion warnings.
  ## Character inputs with observed values therefore receive NA initialization.
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
## Compute and apply one mean/mode initialization without growing a forest.
get.na.roughfix <- function(data) {
  data <- .impute.data.frame(data)
  assign.impute.mean(data, get.impute.mean(data))
}
