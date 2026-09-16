## ------------------------------------------------------------------
## Shared R-side validation for partial and variable-plot extraction.
## These helpers do not select native response slots or change forests.
## ------------------------------------------------------------------
.plot.variable.subset <- function(subset, n) {
  if (is.null(subset)) return(seq_len(n))
  if (is.logical(subset)) {
    if (length(subset) != n) {
      stop("a logical 'subset' must have one entry per stored observation")
    }
    subset <- which(subset)
  }
  else {
    if (!is.numeric(subset) || !is.null(dim(subset))) {
      stop("'subset' must be a numeric or logical vector")
    }
    subset <- subset[is.finite(subset)]
    if (any(subset != floor(subset))) {
      stop("numeric 'subset' entries must be integer row indices")
    }
    subset <- unique(subset[subset >= 1 & subset <= n])
  }
  if (!length(subset)) stop("'subset' does not select any stored observations")
  as.integer(subset)
}
## Resolve a class label or an array position, retaining the original order.
.plot.variable.target <- function(target, labels, what = "target") {
  if (is.null(target)) return(1L)
  if (length(target) != 1L || is.na(target)) {
    stop("'", what, "' must select one value")
  }
  if (is.character(target)) {
    target <- match(match.arg(target, as.character(labels)), as.character(labels))
  }
  if (!is.numeric(target) || !is.finite(target) ||
      target != floor(target) || target < 1L || target > length(labels)) {
    stop("'", what, "' must be a valid label or an index from 1 to ", length(labels))
  }
  as.integer(target)
}
## A time-dependent output uses the last stored time at or before a request.
## Do not call this for mortality or years lost, which have no time axis.
.partial.time.index <- function(times, grid, name = "partial.time", single = FALSE) {
  if (!is.numeric(grid) || !length(grid) || any(!is.finite(grid)) ||
      is.unsorted(grid)) {
    stop("the forest does not contain a valid ordered time.interest grid")
  }
  if (!is.numeric(times) || !is.null(dim(times)) || !length(times) ||
      any(!is.finite(times)) || (single && length(times) != 1L)) {
    stop("'", name, "' must contain ", if (single) "one finite time" else "finite evaluation times")
  }
  if (any(times < grid[1L])) {
    stop("'", name, "' must not precede the first stored time (", grid[1L],
         "). Select a time from the forest's time.interest grid.")
  }
  findInterval(times, grid)
}
## Drop only selected class/event/time axes, never the row or grid axes.
.partial.matrix <- function(x, n, k) {
  if (length(x) != n * k) {
    stop("partial predictions do not match the observation and grid dimensions")
  }
  matrix(x, nrow = n, ncol = k)
}
## A partial object always keeps named lists, including univariate fits.
.partial.response <- function(obj, m.target = NULL) {
  values <- c(obj$classOutput, obj$regrOutput)
  available <- names(values)[!vapply(values, is.null, logical(1L))]
  if (!length(available)) stop("no response predictions are available in the partial object")
  if (is.null(m.target)) {
    candidates <- intersect(obj$yvar.names, available)
    m.target <- if (length(candidates)) candidates[1L] else available[1L]
  }
  else {
    if (!is.character(m.target) || length(m.target) != 1L || is.na(m.target)) {
      stop("'m.target' must select one response name")
    }
    m.target <- match.arg(m.target, available)
  }
  pdta <- values[[m.target]]
  if (is.null(dim(pdta)) || !length(dim(pdta))) {
    stop("partial response predictions do not have array dimensions")
  }
  list(name = m.target, values = pdta,
       classification = m.target %in% names(obj$classOutput))
}
## ------------------------------------------------------------------
## Extract ordinary predictions for marginal variable plots.
## OOB selection remains strict: absent OOB output is not replaced by
## full-ensemble predictions. Preserve the existing display scales.
## ------------------------------------------------------------------
extract.pred <- function(obj, type, subset, time, m.target, target, oob = TRUE) {
  if (missing(m.target)) m.target <- NULL
  m.target <- get.univariate.target(obj, m.target)
  obj <- coerce.multivariate(obj, m.target)
  subset <- .plot.variable.subset(if (missing(subset)) NULL else subset, obj$n)
  n <- length(subset)
  if (!is.logical(oob) || length(oob) != 1L || is.na(oob)) {
    stop("'oob' must be TRUE or FALSE")
  }
  suffix <- if (oob) ".oob" else ""
  pred <- obj[[paste0("predicted", suffix)]]
  surv <- obj[[paste0("survival", suffix)]]
  chf <- obj[[paste0("chf", suffix)]]
  cif <- obj[[paste0("cif", suffix)]]
  if (is.null(pred)) {
    stop(if (oob) "oob requested but oob values are not available in the object" else
           "inbag requested but inbag values are not available in the object")
  }
  target.arg <- if (missing(target)) NULL else target
  if (obj$family == "surv-CR") {
    type <- match.arg(type, c("years.lost", "cif", "chf"))
    target <- .plot.variable.target(target.arg, seq_len(ncol(pred)))
    if (type == "years.lost") return(pred[subset, target])
    time.idx <- .partial.time.index(if (missing(time)) NULL else time,
                                     obj$time.interest, "time", single = TRUE)
    values <- if (type == "cif") cif else chf
    if (is.null(values)) stop("requested survival predictions are unavailable")
    return(as.vector(values[subset, time.idx, target, drop = FALSE]))
  }
  if (obj$family == "surv") {
    type <- match.arg(type, c("rel.freq", "mort", "chf", "surv"))
    if (type == "mort") return(pred[subset])
    if (type == "rel.freq") return(pred[subset] / max(n, na.omit(pred)))
    time.idx <- .partial.time.index(if (missing(time)) NULL else time,
                                     obj$time.interest, "time", single = TRUE)
    values <- if (type == "surv") surv else chf
    if (is.null(values)) stop("requested survival predictions are unavailable")
    ## Preserve the existing marginal survival/CHF percentage-scale output.
    return(100 * as.vector(values[subset, time.idx, drop = FALSE]))
  }
  if (grepl("surv", obj$family)) stop("this survival family is not supported by variable plotting")
  if (obj$family == "class") {
    type <- match.arg(type, c("bayes", "rfq", "prob"))
    prob <- pred[subset, , drop = FALSE]
    y <- obj$yvar[subset]
    if (type == "rfq" && nlevels(y) != 2L) type <- "bayes"
    if (type == "bayes") return(get.bayes.rule(prob))
    if (type == "rfq") return(get.gmean.rule(y, prob))
    target <- .plot.variable.target(target.arg, colnames(prob))
    return(prob[, target])
  }
  pred[subset]
}
## ------------------------------------------------------------------
## Extract one response and one time from an already calculated partial
## object. Every return is selected-observations by primary-grid-values.
## Use get.partial.plot.data() for several stored times at once.
## ------------------------------------------------------------------
extract.partial.pred <- function(obj, type, subset, m.target, target) {
  target.arg <- if (missing(target)) NULL else target
  if (obj$family %in% c("surv", "surv-CR")) {
    values <- obj$survOutput
    if (is.null(values) || is.null(dim(values))) stop("partial survival predictions are unavailable")
    n <- dim(values)[1L]
    k <- tail(dim(values), 1L)
    subset <- .plot.variable.subset(if (missing(subset)) NULL else subset, n)
    if (obj$family == "surv-CR") {
      type <- match.arg(type, c("years.lost", "cif", "chf"))
      j <- if (type == "years.lost") dim(values)[2L] else dim(values)[3L]
      target <- .plot.variable.target(target.arg, seq_len(j))
      if (type == "years.lost") {
        return(.partial.matrix(values[subset, target, , drop = FALSE], length(subset), k))
      }
      if (dim(values)[2L] != 1L) stop("extract.partial.pred() requires exactly one stored evaluation time")
      return(.partial.matrix(values[subset, 1L, target, , drop = FALSE], length(subset), k))
    }
    type <- match.arg(type, c("rel.freq", "mort", "chf", "surv"))
    if (type %in% c("mort", "rel.freq")) {
      out <- .partial.matrix(values[subset, , drop = FALSE], length(subset), k)
      if (type == "rel.freq") {
        ## Preserve column-count scaling; unsupported columns remain NA.
        out[!is.finite(out)] <- NA_real_
        sz <- colSums(!is.na(out))
        sz[sz == 0L] <- NA_integer_
        out <- sweep(out, 2L, sz, "/")
      }
      return(out)
    }
    if (dim(values)[2L] != 1L) stop("extract.partial.pred() requires exactly one stored evaluation time")
    return(.partial.matrix(values[subset, 1L, , drop = FALSE], length(subset), k))
  }
  if (grepl("surv", obj$family)) stop("this survival family is not supported by partial extraction")
  response <- .partial.response(obj, if (missing(m.target)) NULL else m.target)
  values <- response$values
  n <- dim(values)[1L]
  k <- tail(dim(values), 1L)
  subset <- .plot.variable.subset(if (missing(subset)) NULL else subset, n)
  if (response$classification) {
    type <- match.arg(type, c("prob", "bayes"))
    col <- 1L
    if (type == "prob") {
      labels <- dimnames(values)[[2L]][-1L]
      col <- 1L + .plot.variable.target(target.arg, labels)
    }
    return(.partial.matrix(values[subset, col, , drop = FALSE], length(subset), k))
  }
  .partial.matrix(values[subset, , drop = FALSE], length(subset), k)
}
