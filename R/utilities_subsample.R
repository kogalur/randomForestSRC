###################################################################
##
## Select an outcome/statistic and summarize its subsample distribution.
## All three interval summaries are returned; plotting chooses among them.
##
###################################################################
extract.subsample <- function(obj, alpha = .05, target = 0, m.target = NULL, standardize = TRUE, raw = FALSE)
{
  ## Subsample and double-bootstrap objects use separate extraction rules.
  if (sum(c(grepl("rfsrc", class(obj))), grepl("subsample", class(obj))) == 2) {
    subsample <- TRUE
  }
  else if (sum(c(grepl("rfsrc", class(obj))), grepl("bootsample", class(obj))) == 2) {
    subsample <- FALSE
  }
  else {
    stop("this must be an object obtained by calling the 'subsample' function")
  }
  if (!subsample) {
    return(extract.bootsample(obj, alpha = alpha,
              target = target, m.target = m.target, standardize = standardize, raw = raw))    
  }
  .subsample.check.extract(alpha, target, standardize, raw)
  ## Resolve one outcome, then use its response variance/family metadata.
  ## This selects saved output and does not request a new prediction.
  m.target <- get.univariate.target(obj$rf, m.target)
  obj$rf <- coerce.multivariate(obj$rf, m.target)
  ## Extract one fixed, named layout for the selected outcome and column.
  ## This also handles saved banks with differently ordered replicate tables.
  selected <- .subsample.extract.bank(obj$vmp, obj$vmpS, m.target, target)
  vmp <- selected$reference
  theta.hat <- get.subsample.standardize.vimp(obj$rf, vmp, standardize)[, 1]
  theta.star <- get.subsample.standardize.vimp(obj$rf, selected$replicates, standardize)
  rownames(theta.star) <- xvar.names <- rownames(vmp)
  ## Use nominal m = subratio * n, as stored by the current implementation.
  ## The actual integer sample sizes are not available in this object.
  n <- obj$rf$n
  m <- obj$subratio * n
  ## Rescale replicate spread by sqrt(m/n), then center its mean at
  ## theta.hat. These display values are distinct from the centered roots
  ## and the reverse-quantile interval endpoints calculated below.
  scl <- sqrt(obj$subratio)
  theta.star.scl <- theta.hat - scl * rowMeans(theta.star, na.rm = TRUE)
  boxplot.dta <- data.frame(t(scl * theta.star + theta.star.scl), check.names = FALSE)
  ## Centered subsampling roots used by the nonparametric interval.
  z.star <- sqrt(m) * (theta.star - theta.hat) 
  ## Center-at-replicate-mean SE: sqrt(m/n) times the replicate sample SD.
  se.Z <- sqrt(m / n) * apply(theta.star, 1, sd, na.rm = TRUE)
  ## Center-at-full-estimate SE (the existing jackknife convention).
  ## Includes squared mean displacement as well as replicate spread.
  se.jk.Z <- sqrt((m / n) * rowMeans((theta.star - theta.hat)^2, na.rm = TRUE))
  ## Reverse the empirical root quantiles to form interval endpoints.
  ## Rows store lower endpoint, lower quartile, median, upper quartile,
  ## upper endpoint. quantile() uses its existing default interpolation.
  ci <- do.call(cbind, lapply(1:nrow(z.star), function(rowj) {
    c(theta.hat[rowj] - quantile(z.star[rowj, ], 1 - alpha/2, na.rm = TRUE) / sqrt(n),
      theta.hat[rowj] - quantile(z.star[rowj, ], .75, na.rm = TRUE) / sqrt(n),
      theta.hat[rowj] - quantile(z.star[rowj, ], .50, na.rm = TRUE) / sqrt(n),
      theta.hat[rowj] - quantile(z.star[rowj, ], .25, na.rm = TRUE) / sqrt(n),
      theta.hat[rowj] - quantile(z.star[rowj, ], alpha/2, na.rm = TRUE) / sqrt(n))
  }))
  rownames(ci) <- paste(round(100 * c(alpha/2,.25,.50,.75,1-alpha/2), 1), "%", sep = "")
  colnames(ci) <- xvar.names
  ## Normal intervals centered at theta.hat using the subsampling SE.
  ci.Z <- do.call(cbind, lapply(1:nrow(z.star), function(rowj) {
  c(theta.hat[rowj] - qnorm(1 - alpha/2) * se.Z[rowj],
    theta.hat[rowj] - qnorm(.75) * se.Z[rowj],
    theta.hat[rowj],
    theta.hat[rowj] + qnorm(.75) * se.Z[rowj],
    theta.hat[rowj] + qnorm(1 - alpha/2) * se.Z[rowj])
  }))
  rownames(ci.Z) <- paste(round(100 * c(alpha/2,.25,.50,.75,1-alpha/2), 1), "%", sep = "")
  colnames(ci.Z) <- xvar.names
  ## Normal intervals centered at theta.hat using the jackknife SE.
  ci.jk.Z <- do.call(cbind, lapply(1:nrow(z.star), function(rowj) {
  c(theta.hat[rowj] - qnorm(1 - alpha/2) * se.jk.Z[rowj],
    theta.hat[rowj] - qnorm(.75) * se.jk.Z[rowj],
    theta.hat[rowj],
    theta.hat[rowj] + qnorm(.75) * se.jk.Z[rowj],
    theta.hat[rowj] + qnorm(1 - alpha/2) * se.jk.Z[rowj])
  }))
  rownames(ci.jk.Z) <- paste(round(100 * c(alpha/2,.25,.50,.75,1-alpha/2), 1), "%", sep = "")
  colnames(ci.jk.Z) <- xvar.names
  ## Use the same centered subsampling roots for the one-sided test of
  ## zero importance and the nonparametric interval. Keep strict upper-tail
  ## comparisons and the existing quantile interpolation/tie conventions.
  var.sel <- data.frame(lower  =  ci[1, ],
                        median =  ci[3, ],
                        upper  =  ci[5, ],
                        pvalue = rowMeans((z.star - sqrt(n) * theta.hat) > 0, na.rm = TRUE),
                        signif =  ci[1, ] > 0)
  rownames(var.sel) <- xvar.names
  ## Normal one-sided upper-tail p-values for zero VIMP; signif uses
  ## the lower endpoint of the two-sided (1-alpha) interval.
  var.sel.Z <- data.frame(lower  =  ci.Z[1, ],
                          mean   =  ci.Z[3, ],
                          upper  =  ci.Z[5, ],
                          pvalue = pnorm(theta.hat, 0, se.Z, FALSE),
                          signif =  ci.Z[1, ] > 0)
  rownames(var.sel.Z) <- xvar.names
  ## Corresponding upper-tail summary using the jackknife SE.
  var.jk.sel.Z <- data.frame(lower  =  ci.jk.Z[1, ],
                          mean   =  ci.jk.Z[3, ],
                          upper  =  ci.jk.Z[5, ],
                          pvalue = pnorm(theta.hat, 0, se.jk.Z, FALSE),
                          signif =  ci.jk.Z[1, ] > 0)
  rownames(var.jk.sel.Z) <- xvar.names
  if (raw) {
    list(boxplot.dta = boxplot.dta,
         ci = ci,
         ci.Z = ci.Z,
         ci.jk.Z = ci.jk.Z,
         vmp = theta.hat,
         vmpS = theta.star,
         se.Z = se.Z,
         se.jk.Z = se.jk.Z,
         var.sel = var.sel,
         var.sel.Z = var.sel.Z,
         var.jk.sel.Z = var.jk.sel.Z)
  }
  else {
    list(se.Z = se.Z,
         se.jk.Z = se.jk.Z,
         var.sel = var.sel,
         var.sel.Z = var.sel.Z,
         var.jk.sel.Z = var.jk.sel.Z)
  }
}
###################################################################
##
## Default nominal fraction: exp(-1) through n=1000, then proportional
## to n^(-1/2). The same rule is used independently of response family.
##
###################################################################
get.subsample.subratio <- function(n, base.n = 1000, base.r = exp(-1)) {
  if (n < base.n) {
    subratio <- base.r
  }
  else {
    subratio <- base.r * sqrt(base.n) / sqrt(n)
  }
}
###################################################################
##
## Restore joint VIMP for a selected predictor group. No importance
## override is supplied, so vimp() currently selects anti.joint.
##
###################################################################
get.subsample.joint.vimp <- function(o, rf.prms = NULL, xvar.names) {
  class(o)[2] <- "grow"  ## use the grow interface expected by vimp()
  if (is.null(rf.prms)) {
    return(vimp(o, joint = TRUE, xvar.names = xvar.names))
  }
  else {
    return(vimp(o, joint = TRUE, xvar.names = xvar.names,
                perf.type = rf.prms$perf.type, block.size = rf.prms$block.size))
  }
}
###################################################################
##
## Apply one common full-data response-variance scale to regression
## statistics. Other outcome families retain their supplied error scale.
##
###################################################################
get.subsample.standardize.vimp <- function(obj, vmp, standardize = FALSE, sc = 1) {
  if (standardize && !is.null(obj)) {
    vmp <- do.call(cbind, lapply(1:ncol(vmp), function(j) {
      v <- vmp[, j]
      if (obj$family != "regr") {
        sc * v 
      }
      else {
        sc * v / var(obj$yvar, na.rm = TRUE)
      }
    }))
  }
  vmp
}
###################################################################
##
## Construct per-tree counts conditional on a first-stage bootstrap
## multiset. Output rows follow sort(unique(smp)); the caller must use
## that same row order. This helper does not generate the first stage.
##
###################################################################
make.double.boot.sample <- function(ntree, n, smp,
                                    size = function(x) { x },
                                    replace = TRUE)
{
  ntree <- as.integer(ntree)[1]
  n     <- as.integer(n)[1]
  m     <- as.integer(size(n))[1]
  if (is.na(ntree) || ntree < 1L) stop("ntree must be >= 1")
  if (is.na(n)     || n     < 1L) stop("n must be >= 1")
  if (is.na(m)     || m     < 0L) stop("size(n) must be >= 0")
  # Keep only rows that correspond to cases present in the 1st-stage bootstrap sample.
  idx <- sort(unique(smp))  # this is smp.unq, but we need it locally
  k   <- length(idx)
  ## Map first-stage positions to their unique-observation rows.
  ## The map has length(smp) entries, which equals n for an n-draw bootstrap.
  map <- match(smp, idx)  # values in 1..k
  if (replace) {
    ## Probabilities proportional to first-stage frequencies.
    ## sum(freq) equals length(smp).
    freq <- tabulate(map, nbins = k)
    prob <- freq / sum(freq)
    # k x ntree matrix of counts (rows correspond to idx in ascending order)
    dbl.smp <- rmultinom(ntree, size = m, prob = prob)
  } else {
    if (m > length(map)) stop("size(n) cannot exceed length(smp) when replace=FALSE")
    ## Sample positions without replacement from the first-stage multiset.
    ## Different positions can refer to the same original observation.
    dbl.smp <- vapply(seq_len(ntree), function(bb) {
      tabulate(map[sample.int(length(map), size = m, replace = FALSE)], nbins = k)
    }, integer(k))
  }
  ## vapply simplifies when the reduced support has one row.
  dim(dbl.smp) <- c(k, ntree)
  dbl.smp
}
###################################################################
##
## Allocate floor(stratum size * subratio) observations per stratum.
## Missing stratum labels are not represented in table(y).
##
###################################################################
make.strat.sample <- function(y, subratio) {
  class.labels <- names(table(y))
  as.numeric(unlist(lapply(class.labels, function(cl) {
    pt <- which(y == cl)
    size <- as.integer(length(pt) * subratio)
    if (size == 0L) return(integer(0L))
    ## Sample positions so a singleton at row j means row j, not 1:j.
    pt[sample.int(length(pt), size = size, replace = FALSE)]
  }), use.names = FALSE))
}
###################################################################
##
## Print all outcomes/statistic columns and all three interval families.
## The print default alpha=.05 differs from plot.subsample's documented .01.
##
###################################################################
print.subsample.rfsrc <- function(x, alpha = .05, standardize = TRUE, ...) {
  vmp <- x$vmp
  nullO <- lapply(1:length(vmp), function(j) {
    m.target <- names(vmp)[j]
    p.vmp <- ncol(vmp[[j]])
    vmp.col.names <- colnames(vmp[[j]])
    if (length(vmp) > 1) {
      cat("processing a multivariate family, outcome:", m.target, "\n")
    }
    lapply(1:p.vmp, function(k) {
      oo <- extract.subsample(x, alpha = alpha, target = k - 1,
                  m.target = m.target, standardize = standardize, raw = TRUE)
      cat("===== VIMP/error confidence regions for", vmp.col.names[k], " =====\n")
      cat("nonparametric:\n")
      print(round(oo$ci, 3))
      cat("parametric:\n")
      print(round(oo$ci.Z, 3))
      cat("parametric (jackknife):\n")
      print(round(oo$ci.jk.Z, 3))
      NULL
    })
  })
}
print.subsample <- print.subsample.rfsrc
###################################################################
##
## Append statistic rows outcome by outcome, matching names to o1.
## An unavailable append must never relabel an existing predictor.
##
###################################################################
vimp.subsample.combine <- function(o1, o2, o2.name = NULL) {
  if (!length(o1) && !length(o2)) return(NULL)
  master <- if (length(o1)) o1 else o2
  out <- lapply(seq_along(master), function(j) {
    nm <- names(master)[j]
    first <- .subsample.outcome(o1, nm, j)
    added <- .subsample.outcome(o2, nm, j)
    if (!is.null(first)) first <- as.matrix(first)
    if (!is.null(added)) {
      added <- if (is.null(dim(added))) rbind(added) else as.matrix(added)
    }
    if (is.null(first)) {
      if (is.null(added)) return(NULL)
      if (!is.null(o2.name) && nrow(added) > 0L) {
        rownames(added)[nrow(added)] <- o2.name
      }
      return(added)
    }
    if (is.null(added) || nrow(added) == 0L) {
      if (is.null(o2.name)) return(first)
      ## Reserve the requested row as unavailable. Never rename an existing
      ## predictor when no joint/error statistic was returned.
      added <- matrix(NA_real_, 1L, ncol(first),
                      dimnames = list(o2.name, colnames(first)))
    }
    else {
      ii <- .subsample.match.axis(colnames(first), colnames(added),
                                   ncol(first), ncol(added), "statistic columns")
      aligned <- matrix(NA_real_, nrow(added), ncol(first),
                        dimnames = list(rownames(added), colnames(first)))
      ok <- !is.na(ii)
      if (any(ok)) aligned[, ok] <- added[, ii[ok], drop = FALSE]
      added <- aligned
      if (!is.null(o2.name)) rownames(added)[nrow(added)] <- o2.name
    }
    rbind(first, added)
  })
  names(out) <- names(master)
  out
}
###################################################################
##
## Console progress display shared by resampling callers.
##
###################################################################
custom_progress_bar <- function(iteration, total, width = 50) {
  # Calculate the percentage of completion
  percent <- iteration / total
  # Calculate the number of characters to display
  progress <- round(width * percent)
  # ANSI color codes for green
  green_start <- "\033[32m"
  color_end <- "\033[0m"
  # Unicode block character (escaped)
  block_char <- "\u2588"
  # Create the progress bar string with color and blocks
  bar <- paste0("[", green_start, paste(rep(block_char, progress), collapse = ""), 
                color_end, paste(rep(" ", width - progress), collapse = ""), "]")
  # Create the percentage string
  percent_string <- sprintf(" %3d%%", round(100 * percent))
  # Print the progress bar to the console
  cat("\r", bar, percent_string, sep = "")
  # Flush the console output to ensure it is displayed immediately
  flush.console()
  # Clear the progress bar when complete
  if (iteration == total) {
    cat("\r", paste(rep(" ", width + 10), collapse = ""), "\r", sep = "")
  }
}
custom.progress.bar  <- custom_progress_bar
###################################################################
## Small boundary and name-alignment helpers shared by both procedures.
###################################################################
.subsample.check.flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(NULL)
}
.subsample.check.extract <- function(alpha, target, standardize, raw) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a finite number strictly between zero and one.",
         call. = FALSE)
  }
  if (!is.numeric(target) || length(target) != 1L || !is.finite(target) ||
      target < 0 || target != floor(target)) {
    stop("'target' must be a nonnegative integer statistic-column offset.",
         call. = FALSE)
  }
  .subsample.check.flag(standardize, "standardize")
  .subsample.check.flag(raw, "raw")
  invisible(NULL)
}
.subsample.require.reference <- function(x) {
  available <- vapply(x, function(z) {
    !is.null(z) && NROW(z) > 0L && NCOL(z) > 0L
  }, logical(1L))
  if (!length(available) || !any(available)) {
    stop("No reference VIMP or prediction-error statistics are available. ",
         "Check the original forest's importance and performance settings.",
         call. = FALSE)
  }
  invisible(NULL)
}
## Outcome names take precedence; unnamed single-outcome lists retain their
## univariate meaning. Do not pair differently named multivariate outcomes.
.subsample.outcome <- function(x, name, index) {
  if (!length(x)) return(NULL)
  if (!is.null(names(x)) && (anyNA(names(x)) || anyDuplicated(names(x)))) {
    stop("Cannot align ambiguous outcome names in the resampling statistics.",
         call. = FALSE)
  }
  if (length(name) && !is.na(name) && nzchar(name) && !is.null(names(x))) {
    return(x[[name]])
  }
  if (index <= length(x)) x[[index]] else NULL
}
## Match a reference axis by identity. A single unnamed column or two fully
## unnamed axes of equal size have an unambiguous positional interpretation.
.subsample.match.axis <- function(reference, current, nref, ncur, label) {
  if (is.null(reference) || is.null(current)) {
    if (nref == ncur && (nref == 1L ||
                         (is.null(reference) && is.null(current)))) {
      return(seq_len(nref))
    }
    stop("Cannot align unnamed ", label, " with the reference statistics.",
         call. = FALSE)
  }
  if (anyNA(reference) || anyNA(current) ||
      anyDuplicated(reference) || anyDuplicated(current)) {
    stop("Cannot align ambiguous ", label, " in the resampling statistics.",
         call. = FALSE)
  }
  match(reference, current)
}
.subsample.align.table <- function(reference, current) {
  reference <- as.matrix(reference)
  ans <- matrix(NA_real_, nrow(reference), ncol(reference),
                dimnames = dimnames(reference))
  if (!is.null(current) && NROW(current) > 0L && NCOL(current) > 0L) {
    current <- as.matrix(current)
    ir <- .subsample.match.axis(rownames(reference), rownames(current),
                                 nrow(reference), nrow(current), "statistic rows")
    ic <- .subsample.match.axis(colnames(reference), colnames(current),
                                 ncol(reference), ncol(current), "statistic columns")
    rr <- !is.na(ir)
    cc <- !is.na(ic)
    if (any(rr) && any(cc)) ans[rr, cc] <- current[ir[rr], ic[cc], drop = FALSE]
  }
  ans
}
.subsample.align.vimp <- function(reference, current) {
  ans <- lapply(seq_along(reference), function(j) {
    if (is.null(reference[[j]])) return(NULL)
    value <- .subsample.outcome(current, names(reference)[j], j)
    out <- as.data.frame(.subsample.align.table(reference[[j]], value),
                         optional = TRUE)
    ## Retain the established table layout, including unnamed singleton columns.
    colnames(out) <- colnames(reference[[j]])
    out
  })
  names(ans) <- names(reference)
  ans
}
.subsample.extract.bank <- function(reference, bank, m.target, target) {
  .subsample.require.reference(reference)
  j <- if (is.null(m.target)) 1L else match(m.target, names(reference))
  if (length(j) != 1L || is.na(j) || is.null(reference[[j]])) {
    stop("No reference statistics are available for the selected outcome.",
         call. = FALSE)
  }
  ref <- as.matrix(reference[[j]])
  if (target >= ncol(ref)) {
    stop("'target' exceeds the available statistic columns.", call. = FALSE)
  }
  if (!length(bank)) stop("The replicate bank is empty.", call. = FALSE)
  ans <- matrix(NA_real_, nrow(ref), length(bank),
                dimnames = list(rownames(ref), NULL))
  for (b in seq_along(bank)) {
    value <- .subsample.outcome(bank[[b]], names(reference)[j], j)
    ans[, b] <- .subsample.align.table(ref, value)[, 1L + target]
  }
  list(reference = ref[, 1L + target, drop = FALSE], replicates = ans)
}
## Check outcomes after the established omission policy, without imputing or
## changing any rows. Entirely missing predictor columns are removed by grow
## preprocessing; outcome columns must remain observed for an event to count.
.subsample.check.events <- function(data, yvar.names, replicate) {
  keep <- colSums(!is.na(data)) > 0L | names(data) %in% yvar.names
  complete <- stats::complete.cases(data[, keep, drop = FALSE])
  status <- data[[yvar.names[2L]]]
  if (!any(complete & !is.na(status) & status > 0)) {
    stop("Survival subsample ", replicate,
         " has no observed events after omission of missing records. ",
         "Increase 'subratio' so that events are retained, and check missing ",
         "values in the event records.", call. = FALSE)
  }
  invisible(NULL)
}
