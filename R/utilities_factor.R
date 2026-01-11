check.factor <- function(test, gfactor, ignore = TRUE) {
  if (is.null(gfactor)) return(test)
  cn <- colnames(test)
  ## ------------------------------------------------------------------
  ## unordered factors
  ## ------------------------------------------------------------------
  if (length(gfactor$factor) > 0) {
    idx <- match(gfactor$factor, cn)
    new_cols <- lapply(seq_along(idx), function(k) {
      j <- idx[k]
      ## pull vector fast (avoid colnames(test) == ...)
      fk.test <- as.character(test[[j]])
      levs    <- gfactor$levels[[k]]
      ## if test introduces levels not seen in training, optionally coerce
      if (!all(is.element(unique(stats::na.omit(fk.test)), levs))) {
        if (!ignore) {
          stop("levels of factors in test data do not match those in training data\n")
        }
        fake.level.k <- paste0(tail(levs, 1), "+")
        fk.test[!is.element(fk.test, levs)] <- fake.level.k
        factor(fk.test,
               levels  = c(levs, fake.level.k),
               exclude = NULL)
      } else {
        factor(fk.test,
               levels  = levs,
               exclude = NULL)
      }
    })
    ## assign back (column-wise) without constructing a data.frame
    test[idx] <- new_cols
  }
  ## ------------------------------------------------------------------
  ## ordered factors
  ## ------------------------------------------------------------------
  if (length(gfactor$order) > 0) {
    idx <- match(gfactor$order, cn)
    new_cols <- lapply(seq_along(idx), function(k) {
      j <- idx[k]
      fk.test <- as.character(test[[j]])
      levs    <- gfactor$order.levels[[k]]
      if (!all(is.element(unique(stats::na.omit(fk.test)), levs))) {
        if (!ignore) {
          stop("levels of ordered factors in test data do not match those in training data\n")
        }
        fake.level.k <- paste0(tail(levs, 1), "+")
        fk.test[!is.element(fk.test, levs)] <- fake.level.k
        factor(fk.test,
               levels  = c(levs, fake.level.k),
               ordered = TRUE)
      } else {
        factor(fk.test,
               levels  = levs,
               ordered = TRUE)
      }
    })
    test[idx] <- new_cols
  }
  test
}
is.factor.not.ordered <- function(x) { is.factor(x) && !is.ordered(x) }
extract.factor <- function(dat, generic.names = NULL) {
  generic.types <- gfactor <- gfactor.order <- gfactor.levels <- gfactor.order.levels <- NULL
  target.names <- if (is.null(generic.names)) names(dat) else generic.names
  nlevels <- rep(0, length(target.names))
  ## identify factor and ordered factor columns
  gfactor       <- names(dat)[vapply(dat, is.factor.not.ordered, logical(1))]
  gfactor.order <- names(dat)[vapply(dat, is.ordered,          logical(1))]
  ## restrict to generic names if requested
  if (!is.null(generic.names)) {
    gfactor       <- intersect(gfactor, generic.names)
    gfactor.order <- intersect(gfactor.order, generic.names)
  }
  ## collect levels and level counts
  if (length(gfactor) > 0) {
    gfactor.levels <- lapply(gfactor, function(nm) levels(dat[[nm]]))
    nlevels[match(gfactor, target.names)] <- lengths(gfactor.levels)
  }
  if (length(gfactor.order) > 0) {
    gfactor.order.levels <- lapply(gfactor.order, function(nm) levels(dat[[nm]]))
    nlevels[match(gfactor.order, target.names)] <- lengths(gfactor.order.levels)
  }
  ## generic types bookkeeping
  if (!is.null(generic.names)) {
    generic.types <- rep("R", length(generic.names))
    if (length(gfactor) > 0) {
      generic.types[match(gfactor, generic.names)] <- "C"
    }
    if (length(gfactor.order) > 0) {
      generic.types[match(gfactor.order, generic.names)] <- "I"
    }
  } else {
    generic.types <- rep("R", ncol(dat))
    if (length(gfactor) > 0) {
      generic.types[match(gfactor, names(dat))] <- "C"
    }
    if (length(gfactor.order) > 0) {
      generic.types[match(gfactor.order, names(dat))] <- "I"
    }
  }
  ## preserve legacy behaviour: return NULL only if factor lists are NULL
  ## (in practice these are usually character(0), but keep this consistent)
  if (is.null(gfactor) & is.null(gfactor.order)) {
    return(NULL)
  } else {
    return(list(
      factor        = gfactor,
      order         = gfactor.order,
      levels        = gfactor.levels,
      order.levels  = gfactor.order.levels,
      nlevels       = nlevels,
      generic.types = generic.types
    ))
  }
}
map.factor <- function(gvar, gfactor) {
  if (is.null(gfactor)) return(gvar)
  cn <- colnames(gvar)
  ## ------------------------------------------------------------------
  ## unordered factors
  ## ------------------------------------------------------------------
  if (length(gfactor$factor) > 0) {
    idx <- match(gfactor$factor, cn)
    new_cols <- lapply(seq_along(idx), function(k) {
      j <- idx[k]
      levs <- gfactor$levels[[k]]
      codes <- gvar[[j]]
      ## mapping from integer codes back to factor labels
      factor.k <- levs[codes]
      labels.k <- levs[sort(unique(codes))]
      gk <- factor(factor.k,
                   labels  = labels.k,
                   levels  = labels.k,
                   exclude = NULL)
      if (length(setdiff(levs, labels.k)) > 0) {
        gk <- factor(as.character(gk), levels = levs)
      }
      gk
    })
    gvar[idx] <- new_cols
  }
  ## ------------------------------------------------------------------
  ## ordered factors
  ## ------------------------------------------------------------------
  if (length(gfactor$order) > 0) {
    idx <- match(gfactor$order, cn)
    new_cols <- lapply(seq_along(idx), function(k) {
      j <- idx[k]
      levs <- gfactor$order.levels[[k]]
      codes <- gvar[[j]]
      factor.k <- levs[codes]
      labels.k <- levs[sort(unique(codes))]
      gk <- factor(factor.k,
                   labels  = labels.k,
                   levels  = labels.k,
                   exclude = NULL,
                   ordered = TRUE)
      if (length(setdiff(levs, labels.k)) > 0) {
        gk <- factor(as.character(gk), levels = levs, ordered = TRUE)
      }
      gk
    })
    gvar[idx] <- new_cols
  }
  gvar
}
rm.na.levels <- function(dat, xvar.names = NULL) {
  factor.names <- names(dat)[vapply(dat, is.factor, logical(1))]
  if (!is.null(xvar.names)) factor.names <- intersect(factor.names, xvar.names)
  if (length(factor.names) > 0) {
    ## identify factors with a literal "NA" level (string)
    levels.na.pt <- vapply(factor.names, function(nm) {
      any(levels(dat[[nm]]) == "NA", na.rm = TRUE)
    }, logical(1))
    if (any(levels.na.pt)) {
      fn <- factor.names[levels.na.pt]
      idx <- match(fn, names(dat))
      new_cols <- lapply(fn, function(nm) {
        x <- dat[[nm]]
        levels(x)[levels(x) == "NA"] <- NA
        x
      })
      dat[idx] <- new_cols
    }
  }
  dat
}
