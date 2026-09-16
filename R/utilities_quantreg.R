## Internal validation shared by quantile construction and extraction.
.quantreg.check.prob <- function(prob, name = "prob", endpoints = FALSE) {
  if (!is.numeric(prob) || !is.null(dim(prob)) || length(prob) == 0L ||
      any(!is.finite(prob)) ||
      any(if (endpoints) prob < 0 | prob > 1 else prob <= 0 | prob >= 1)) {
    stop("'", name, "' must be a nonempty vector of finite probabilities in ",
         if (endpoints) "[0, 1]." else "(0, 1).", call. = FALSE)
  }
  invisible(TRUE)
}
.quantreg.check.flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}
## Select an outcome without dropping or confusing response dimensions.
.quantreg.response <- function(y, target) {
  if (is.data.frame(y) || is.matrix(y)) {
    if (!(target %in% colnames(y))) {
      stop("Outcome '", target, "' is missing from the saved response data.",
           call. = FALSE)
    }
    y <- y[, target, drop = TRUE]
  }
  if (!is.numeric(y)) {
    stop("Outcome '", target, "' must be numeric.", call. = FALSE)
  }
  as.numeric(y)
}
.quantreg.vector <- function(x, n, label) {
  if (!is.numeric(x) || length(x) != n ||
      (!is.null(dim(x)) && NCOL(x) != 1L)) {
    stop(label, " must contain one numeric value per observation.",
         call. = FALSE)
  }
  as.numeric(x)
}
## A NULL OOB component permits the established full-ensemble fallback.
## Missing entries within an existing OOB component remain unavailable.
.quantreg.prediction <- function(component, oob, quantile = FALSE) {
  field <- if (quantile) "quantile" else "predicted"
  use.oob <- oob && !is.null(component[[paste0(field, ".oob")]])
  value <- component[[if (use.oob) paste0(field, ".oob") else field]]
  if (is.null(value)) {
    stop("The forest did not return the required '", field,
         "' output.", call. = FALSE)
  }
  list(value = value, oob = use.oob)
}
## Training data already belong to the saved forest. On a grow object its
## top-level yvar is also training data; on a predict object it is not.
## The final fallback reads older quantreg objects without writing legacy fields.
.quantreg.training.yvar <- function(o) {
  forest <- o[["forest"]]
  y <- if (is.list(forest)) forest[["yvar"]] else NULL
  if (is.null(y) && inherits(o, "grow")) y <- o[["yvar"]]
  if (is.null(y)) y <- o[["yvar.grow"]]
  y
}
## Read the selected training residuals, never residuals of the evaluated rows.
## A .oob suffix identifies numerical OOB residuals. Selection is kept separately
## in forest$quantreg$residual.source and in each returned quantile component.
.quantreg.training.residual <- function(o, target, source = NULL) {
  forest <- o[["forest"]]
  if (is.null(source) && is.list(forest)) {
    source <- forest[["quantreg"]][["residual.source"]][[target]]
  }
  if (!is.null(source)) {
    if (!is.character(source) || length(source) != 1L || is.na(source) ||
        !(source %in% c("oob", "ensemble"))) {
      stop("Invalid saved residual source for outcome '", target, "'.",
           call. = FALSE)
    }
    field <- if (source == "oob") "residual.oob" else "residual"
    bank <- if (is.list(forest)) forest[[field]] else NULL
    if (is.null(bank) && inherits(o, "grow")) {
      component <- o[["regrOutput"]][[target]]
      if (is.null(component)) component <- o
      bank <- component[[field]]
    }
  } else {
    ## Older objects saved only the selected bank. Keep their scoring usable;
    ## migration of a grow object checks its source before assigning new names.
    bank <- o[["res.grow"]]
  }
  if (is.null(bank)) {
    stop("Saved training residuals are unavailable for outcome '", target,
         "'. Retain the original quantreg grow object and its forest.",
         call. = FALSE)
  }
  residual <- .quantreg.response(bank, target)
  y <- .quantreg.training.yvar(o)
  if (!is.null(y) && length(residual) != length(.quantreg.response(y, target))) {
    stop("Training responses and residuals are not aligned for outcome '",
         target, "'.", call. = FALSE)
  }
  residual
}
## Store both numerical training-residual quantities beside the grow predictions
## and in the saved forest. The initial oob setting selects the bank to reuse;
## changing the current prediction setting later does not change that selection.
## No extra forest fit or prediction call is needed to construct these residuals.
.quantreg.save.residuals <- function(object, oob, initialize = FALSE) {
  forest <- object[["forest"]]
  if (!is.list(forest)) {
    stop("The quantreg grow object must retain its forest settings.", call. = FALSE)
  }
  univariate <- identical(object$family, "regr")
  ynames <- if (univariate) object$yvar.names else names(object$regrOutput)
  if (!length(ynames)) {
    stop("The forest did not return continuous-outcome predictions.", call. = FALSE)
  }
  settings <- forest[["quantreg"]]
  if (is.null(settings)) settings <- list()
  if (!is.list(settings)) {
    stop("Saved quantreg settings must be a list.", call. = FALSE)
  }
  source <- settings[["residual.source"]]
  if (!is.null(source)) {
    ## Already initialized: preserve the original residual source and values.
    if (!is.character(source) || is.null(names(source)) || anyDuplicated(names(source)) ||
        !all(ynames %in% names(source))) {
      stop("Saved residual sources do not identify the continuous outcomes.",
           call. = FALSE)
    }
    for (yn in ynames) .quantreg.training.residual(object, yn, source[[yn]])
  } else {
    n <- object$n
    residual <- residual.oob <- setNames(vector("list", length(ynames)), ynames)
    source <- setNames(character(length(ynames)), ynames)
    legacy <- object[["res.grow"]]
    legacy.flag <- object[["res.grow.oob"]]
    if (is.null(legacy) && !initialize) {
      stop("The quantreg grow object is missing its original residual selection. ",
           "Retain the original grow object or regrow with quantreg().",
           call. = FALSE)
    }
    for (yn in ynames) {
      component <- if (univariate) object else object$regrOutput[[yn]]
      y <- .quantreg.response(object[["yvar"]], yn)
      if (length(y) != n) {
        stop("Training responses and forest rows are not aligned.", call. = FALSE)
      }
      pred <- component[["predicted"]]
      pred.oob <- component[["predicted.oob"]]
      r <- if (is.null(pred)) NULL else {
        y - .quantreg.vector(pred, n, "Full-ensemble training predictions")
      }
      r.oob <- if (is.null(pred.oob)) NULL else {
        y - .quantreg.vector(pred.oob, n, "OOB training predictions")
      }
      if (is.null(legacy)) {
        ## Preserve the existing component-level fallback, including all-NA OOB
        ## components: unavailable individual predictions are not filled in.
        selected <- .quantreg.prediction(component, oob)
        source[[yn]] <- if (selected$oob) "oob" else "ensemble"
      } else {
        ## Migrate an older grow object without changing its selected bank.
        ## Use recorded provenance when available, otherwise identify the bank
        ## from the retained original predictions. Never guess an OOB label.
        saved <- .quantreg.response(legacy, yn)
        if (length(saved) != n) {
          stop("Saved training residuals and forest rows are not aligned.",
               call. = FALSE)
        }
        same <- function(a, b) {
          !is.null(b) && isTRUE(all.equal(a, b, tolerance = 0,
                                         check.attributes = FALSE))
        }
        flag <- if (!is.null(names(legacy.flag)) && yn %in% names(legacy.flag)) {
          legacy.flag[[yn]]
        } else if (length(ynames) == 1L && length(legacy.flag) == 1L) {
          legacy.flag[[1L]]
        } else NULL
        if (is.logical(flag) && length(flag) == 1L && !is.na(flag)) {
          source[[yn]] <- if (flag) "oob" else "ensemble"
          expected <- if (flag) r.oob else r
          if (!is.null(expected) && !same(saved, expected)) {
            stop("Saved residuals disagree with their recorded source for outcome '",
                 yn, "'. Use the unmodified original quantreg grow object.",
                 call. = FALSE)
          }
        } else {
          matches.oob <- same(saved, r.oob)
          matches.ensemble <- same(saved, r)
          if (matches.oob == matches.ensemble) {
            stop("The original residual source cannot be established for outcome '",
                 yn, "'. Regrow the quantreg object to initialize named residuals.",
                 call. = FALSE)
          }
          source[[yn]] <- if (matches.oob) "oob" else "ensemble"
        }
        if (source[[yn]] == "oob") r.oob <- saved else r <- saved
      }
      residual[yn] <- list(r)
      residual.oob[yn] <- list(r.oob)
      if (univariate) {
        object["residual"] <- list(r)
        object["residual.oob"] <- list(r.oob)
      } else {
        component["residual"] <- list(r)
        component["residual.oob"] <- list(r.oob)
        object$regrOutput[[yn]] <- component
      }
    }
    ## Saved multivariate banks use one column per continuous response, with
    ## training rows unchanged. An unavailable component has an NA column;
    ## a wholly unavailable bank is NULL. Grow outputs keep response nesting.
    pack <- function(bank) {
      if (all(vapply(bank, is.null, logical(1)))) return(NULL)
      if (univariate) return(bank[[1L]])
      ans <- do.call(cbind, lapply(bank, function(z) {
        if (is.null(z)) rep(NA_real_, n) else z
      }))
      colnames(ans) <- ynames
      ans
    }
    ## Preserve named NULL slots so $residual cannot partially match
    ## residual.oob when full-ensemble predictions were not retained.
    forest["residual"] <- list(pack(residual))
    forest["residual.oob"] <- list(pack(residual.oob))
    settings[["residual.source"]] <- source
    forest[["quantreg"]] <- settings
    object[["forest"]] <- forest
  }
  ## Legacy fields are read only for migration, not returned as parallel aliases.
  object[["yvar.grow"]] <- NULL
  object[["res.grow"]] <- NULL
  object[["res.grow.oob"]] <- NULL
  object
}
## The reporting grid uses finite training responses, separately for each outcome.
.quantreg.grid <- function(y, maxn) {
  grid <- sort(unique(y[is.finite(y)]))
  if (!length(grid)) {
    stop("No finite training responses are available for the reporting grid.",
         call. = FALSE)
  }
  if (length(grid) > maxn) {
    grid <- grid[unique(round(seq(1, length(grid), length.out = maxn)))]
  }
  grid
}
## Increments are interval masses, not a density divided by grid spacing.
## Explicit matrix indexing preserves both singleton and zero-row dimensions.
.quantreg.increments <- function(cdf) {
  density <- cdf
  k <- ncol(cdf)
  if (k > 1L) {
    density[, 2:k] <- cdf[, 2:k, drop = FALSE] -
      cdf[, seq_len(k - 1L), drop = FALSE]
  }
  density
}
## Lower generalized inverse: the first cumulative weight reaching each level.
## Weights of zero have already been removed. Tied response values are allowed.
.quantreg.inverse <- function(value, cumulative, prob) {
  cumulative[length(cumulative)] <- 1
  index <- 1L + findInterval(prob, cumulative, left.open = TRUE)
  value[pmin(length(value), index)]
}
## Probability-distribution moments (no sample-variance degrees-of-freedom factor).
## Center before computing spread; scaling avoids overflow when squaring deviations.
.quantreg.moments <- function(value, weight = NULL) {
  if (is.null(weight)) weight <- rep(1 / length(value), length(value))
  center <- value[ceiling(length(value) / 2)]
  delta <- value - center
  if (all(is.finite(delta))) {
    delta.mean <- sum(weight * delta)
    mu <- center + delta.mean
    deviation <- delta - delta.mean
  } else {
    scale <- max(abs(value))
    mu.scaled <- sum(weight * (value / scale))
    mu <- scale * mu.scaled
    deviation <- (value / scale) - mu.scaled
    spread <- max(abs(deviation))
    sd <- if (spread == 0) 0 else {
      scale * (spread * sqrt(sum(weight * (deviation / spread)^2)))
    }
    return(c(mean = mu, std = sd))
  }
  spread <- max(abs(deviation))
  sd <- if (spread == 0) 0 else {
    spread * sqrt(sum(weight * (deviation / spread)^2))
  }
  c(mean = mu, std = sd)
}
## Summarize the saved residual distribution and shift by each current prediction.
## NULL weights give the common exchangeable empirical residual distribution;
## forest weights give a row-specific residual distribution.
## Sorting is shared across rows. The shifted observation-by-donor matrix is never
## constructed. An empty grid requests quantiles/moments only (used by printing).
.quantreg.residual <- function(residual, prediction, prob, grid,
                               weight = NULL) {
  n <- length(prediction)
  k <- length(grid)
  if (!is.numeric(residual) || !is.numeric(prediction)) {
    stop("Residuals and predictions must be numeric.", call. = FALSE)
  }
  if (!is.null(weight) &&
      (!is.matrix(weight) || !is.numeric(weight) ||
       nrow(weight) != n || ncol(weight) != length(residual))) {
    stop("Forest weights must have one row per prediction and one column ",
         "per saved training residual.", call. = FALSE)
  }
  out <- list(quantiles = matrix(NA_real_, n, length(prob)),
              cdf = matrix(NA_real_, n, k),
              mean = rep(NA_real_, n), std = rep(NA_real_, n))
  donor <- which(is.finite(residual))
  donor <- donor[order(residual[donor])]
  residual <- residual[donor]
  nr <- length(residual)
  if (nr == 0L || n == 0L) return(out)
  if (is.null(weight)) {
    rq <- .quantreg.inverse(residual, seq_len(nr) / nr, prob)
    moments <- .quantreg.moments(residual)
    for (i in which(is.finite(prediction))) {
      out$quantiles[i, ] <- prediction[i] + rq
      out$mean[i] <- prediction[i] + moments[["mean"]]
      out$std[i] <- moments[["std"]]
      if (k) {
        out$cdf[i, ] <- findInterval(grid - prediction[i], residual) / nr
      }
    }
  } else {
    for (i in seq_len(n)) {
      ## Remove unavailable residual donors and their corresponding columns
      ## together. Missing retained weights invalidate the row, not the donor.
      w <- weight[i, donor]
      if (any(w < 0, na.rm = TRUE)) {
        stop("Forest weights must be nonnegative.", call. = FALSE)
      }
      if (!is.finite(prediction[i]) || any(!is.finite(w))) next
      positive <- w > 0
      if (!any(positive)) next
      r <- residual[positive]
      w <- w[positive]
      w <- w / max(w)
      w <- w / sum(w)
      cumulative <- pmin(1, cumsum(w))
      cumulative[length(cumulative)] <- 1
      out$quantiles[i, ] <- prediction[i] +
        .quantreg.inverse(r, cumulative, prob)
      moments <- .quantreg.moments(r, w)
      out$mean[i] <- prediction[i] + moments[["mean"]]
      out$std[i] <- moments[["std"]]
      if (k) {
        position <- findInterval(grid - prediction[i], r)
        out$cdf[i, ] <- c(0, cumulative)[position + 1L]
      }
    }
  }
  out
}
## Native GK estimates are retained. Reconstruct a step CDF by taking the largest
## requested probability whose quantile is at or below the evaluation point.
## The upper training endpoint completes the reporting distribution at one.
.quantreg.gk.cdf <- function(quantiles, prob, grid, upper) {
  cdf <- matrix(NA_real_, nrow(quantiles), length(grid))
  for (i in seq_len(nrow(quantiles))) {
    q <- quantiles[i, ]
    if (any(!is.finite(q))) next
    if (any(diff(q) < 0)) {
      stop("Native GK quantiles are not nondecreasing in probability.",
           call. = FALSE)
    }
    cdf[i, ] <- c(0, prob)[findInterval(grid, q) + 1L]
    cdf[i, grid >= upper] <- 1
  }
  cdf
}
## Extract the component layout structurally: a response may itself be named
## 'prob', 'cdf', 'quantiles', or any other ordinary component name.
extract.quantile <- function(o) {
  if (!inherits(o, "quantreg") || !is.list(o$quantreg)) {
    stop("object must be a quantreg object", call. = FALSE)
  }
  is.component <- function(q) {
    is.list(q) && is.numeric(q[["prob"]]) &&
      is.matrix(q[["quantiles"]]) && is.numeric(q[["quantiles"]]) &&
      is.matrix(q[["cdf"]]) && is.numeric(q[["cdf"]]) &&
      is.matrix(q[["density"]]) && is.numeric(q[["density"]]) &&
      is.numeric(q[["yunq"]])
  }
  if (is.component(o$quantreg)) {
    if (length(o$yvar.names) != 1L) {
      stop("The univariate quantile component must identify one response.",
           call. = FALSE)
    }
    q <- setNames(list(o$quantreg), o$yvar.names)
  } else {
    q <- o$quantreg
    q <- q[!vapply(q, is.null, logical(1))]
    if (!length(q) || is.null(names(q)) || anyNA(names(q)) ||
        any(!nzchar(names(q))) || anyDuplicated(names(q)) ||
        any(!names(q) %in% o$yvar.names) ||
        !all(vapply(q, is.component, logical(1)))) {
      stop("The object does not contain valid named quantile components.",
           call. = FALSE)
    }
  }
  for (component in q) {
    .quantreg.check.prob(component$prob)
    n <- nrow(component$quantiles)
    if (ncol(component$quantiles) != length(component$prob) ||
        !identical(dim(component$cdf), dim(component$density)) ||
        nrow(component$cdf) != n ||
        ncol(component$cdf) != length(component$yunq) ||
        !length(component$yunq) || any(!is.finite(component$yunq)) ||
        any(diff(component$yunq) <= 0)) {
      stop("Quantile, probability, and reporting-grid dimensions are inconsistent.",
           call. = FALSE)
    }
  }
  q
}
## Nearest stored probability lookup is retained. It does not recompute quantiles.
get.quantile <- function(o, target.prob = NULL, pretty = TRUE) {
  .quantreg.check.flag(pretty, "pretty")
  qo <- extract.quantile(o)
  if (!is.null(target.prob)) {
    .quantreg.check.prob(target.prob, "target.prob", endpoints = TRUE)
    target.prob <- sort(unique(target.prob))
  }
  rO <- lapply(qo, function(q) {
    pr <- if (is.null(target.prob)) q$prob else target.prob
    index <- vapply(pr, function(p) which.min(abs(p - q$prob)), integer(1))
    q.dat <- q[["quantiles"]][, index, drop = FALSE]
    colnames(q.dat) <- paste0("q.", 100 * pr)
    q.dat
  })
  if (pretty && length(rO) == 1L) rO <- rO[[1L]]
  rO
}
## Logical selectors must cover every evaluated row. Numeric selectors retain
## duplicates (repeated inclusion) but must contain valid positive integer indices.
.quantreg.subset <- function(subset, n) {
  if (is.null(subset)) subset <- seq_len(n)
  if (is.logical(subset)) {
    if (length(subset) != n || anyNA(subset)) {
      stop("A logical 'subset' must have one nonmissing entry per evaluated row.",
           call. = FALSE)
    }
    subset <- which(subset)
  } else if (!is.numeric(subset) || !is.null(dim(subset)) ||
             any(!is.finite(subset)) || any(subset != floor(subset)) ||
             any(subset < 1 | subset > n)) {
    stop("'subset' must contain positive integer indices of evaluated rows.",
         call. = FALSE)
  }
  if (!length(subset)) stop("requested subset is empty", call. = FALSE)
  subset
}
## Keep the established finite-grid score curve. Average one threshold at a time
## and accumulate each trapezoid once, avoiding another row-by-grid matrix and
## repeated integration of all earlier prefixes.
get.quantile.crps <- function(o, pretty = TRUE, subset = NULL, standardize = TRUE) {
  .quantreg.check.flag(pretty, "pretty")
  .quantreg.check.flag(standardize, "standardize")
  qO <- extract.quantile(o)
  if (is.null(o$yvar)) {
    stop("Observed responses are required for quantile CRPS.", call. = FALSE)
  }
  n <- nrow(qO[[1L]]$cdf)
  subset <- .quantreg.subset(subset, n)
  rO <- lapply(names(qO), function(target) {
    q <- qO[[target]]
    y <- .quantreg.response(o$yvar, target)
    if (length(y) != n || nrow(q$cdf) != n) {
      stop("Observed responses and CDF rows are not aligned.", call. = FALSE)
    }
    y <- y[subset]
    brS <- vapply(seq_along(q$yunq), function(k) {
      p <- q$cdf[subset, k]
      ok <- is.finite(y) & is.finite(p)
      if (!any(ok)) return(NA_real_)
      mean((as.numeric(y[ok] <= q$yunq[k]) - p[ok])^2)
    }, numeric(1))
    integral <- c(0, cumsum(diff(q$yunq) *
                           (utils::head(brS, -1L) + utils::tail(brS, -1L)) / 2))
    if (standardize) {
      width <- q$yunq - q$yunq[1L]
      value <- rep(NA_real_, length(width))
      positive <- width > 0
      value[positive] <- integral[positive] / width[positive]
    } else {
      value <- integral
    }
    data.frame(y = q$yunq, crps = value)
  })
  names(rO) <- names(qO)
  if (pretty && length(rO) == 1L) rO <- rO[[1L]]
  rO
}
## Grid moment fallback for GK and older objects. This retains the original
## unnormalized grid-moment definition, including incomplete represented mass.
## The centered identity is E_grid[(Y - mu)^2] + (1 - mass) * mu^2.
.quantreg.grid.moments <- function(q) {
  n <- nrow(q$density)
  out <- data.frame(mean = rep(NA_real_, n), std = rep(NA_real_, n))
  for (i in seq_len(n)) {
    w <- q$density[i, ]
    if (any(!is.finite(w)) || any(w < 0)) next
    mass <- sum(w)
    if (mass <= 0 || mass > 1 + 64 * .Machine$double.eps) next
    if (abs(mass - 1) <= 64 * .Machine$double.eps) mass <- 1
    w <- w / sum(w)
    moments <- .quantreg.moments(q$yunq, w)
    mu <- mass * moments[["mean"]]
    ## Algebraically the previous raw-moment result, evaluated without
    ## subtracting two nearly equal large squared response values.
    a <- sqrt(mass) * moments[["std"]]
    b <- sqrt(mass * (1 - mass)) * abs(moments[["mean"]])
    scale <- max(a, b)
    out$mean[i] <- mu
    out$std[i] <- if (scale == 0) 0 else {
      scale * sqrt((a / scale)^2 + (b / scale)^2)
    }
  }
  out
}
## Residual-method moments are computed at construction, independently of yunq.
## The median retains the public nearest-stored-probability extraction policy.
get.quantile.stat <- function(o, pretty = TRUE) {
  .quantreg.check.flag(pretty, "pretty")
  qO <- extract.quantile(o)
  mdn <- get.quantile(o, .5, FALSE)
  rO <- lapply(names(qO), function(target) {
    q <- qO[[target]]
    n <- nrow(q$quantiles)
    if (!is.null(q[["mean"]]) && !is.null(q[["std"]])) {
      mn <- .quantreg.vector(q[["mean"]], n, "Quantile-distribution means")
      sd <- .quantreg.vector(q[["std"]], n, "Quantile-distribution standard deviations")
    } else {
      moments <- .quantreg.grid.moments(q)
      mn <- moments$mean
      sd <- moments$std
    }
    data.frame(mean = mn, median = as.numeric(mdn[[target]]), std = sd)
  })
  names(rO) <- names(qO)
  if (pretty && length(rO) == 1L) rO <- rO[[1L]]
  rO
}
## Resolve reporting levels without changing the session's options or the fit.
## Explicit levels take precedence; the former option is a compatibility fallback.
.quantreg.resolve.tau <- function(tau = NULL, name = "quantreg.tau") {
  if (is.null(tau)) {
    tau <- getOption("quantreg.tau")
    name <- "quantreg.tau"
    if (is.null(tau)) {
      tau <- getOption("rfsrc.pinball.taus")
      if (!is.null(tau)) name <- "rfsrc.pinball.taus"
    }
    if (is.null(tau)) tau <- c(.1, .5, .9)
  }
  .quantreg.check.prob(tau, name)
  as.numeric(tau)
}
## Mean pinball loss at each requested level, using the object's current quantiles
## or its saved predictive distribution. Each response is evaluated separately.
## This is extraction/scoring only: no trees are grown and no prediction is called.
get.pinball.error <- function(o, tau = NULL, subset = NULL,
                              m.target = NULL, pretty = TRUE) {
  .quantreg.check.flag(pretty, "pretty")
  tau <- .quantreg.resolve.tau(tau, name = "tau")
  qO <- extract.quantile(o)
  if (is.null(o[["yvar"]])) {
    stop("Observed responses are required for pinball loss.", call. = FALSE)
  }
  if (!is.null(m.target)) {
    if (!is.character(m.target) || length(m.target) != 1L ||
        is.na(m.target) || !(m.target %in% names(qO))) {
      stop("'m.target' must name one available continuous response: ",
           paste(names(qO), collapse = ", "), ".", call. = FALSE)
    }
    qO <- qO[m.target]
  }
  out <- lapply(names(qO), function(target) {
    ## A response-specific view keeps the training banks in the saved forest,
    ## and leaves evaluated responses and forest weights in their original order.
    view <- o
    view$quantreg <- qO[[target]]
    view$family <- "regr"
    view$yvar.names <- target
    view$yvar <- .quantreg.response(o[["yvar"]], target)
    .quantreg.pinball(view, taus = tau, subset = subset)
  })
  names(out) <- names(qO)
  if (pretty && length(out) == 1L) out <- out[[1L]]
  out
}
## Quantile scoring uses the residual law for exchangeable/forest objects.
## The former stored method name "local" follows the same residual path.
## Native GK and legacy objects retain CDF interpolation for reporting pinball
## losses; a one-point CDF is a constant interpolation rather than an approx error.
.quantreg.pinball <- function(o, taus = c(.1, .5, .9), subset = NULL) {
  .quantreg.check.prob(taus, "tau")
  qO <- extract.quantile(o)
  if (length(qO) != 1L) {
    stop("Select one continuous outcome for printed pinball losses.",
         call. = FALSE)
  }
  q <- qO[[1L]]
  y <- .quantreg.response(o$yvar, names(qO)[1L])
  n <- nrow(q$quantiles)
  if (length(y) != n) stop("Responses and quantiles are not aligned.", call. = FALSE)
  index <- .quantreg.subset(subset, n)
  if (!is.null(q$prediction) &&
      q$method %in% c("exchangeable", "forest", "local")) {
    stored <- match(taus, q$prob)
    if (!anyNA(stored)) {
      ## Default reporting levels are often already stored: avoid repeating
      ## the weighted-residual calculation when exact requested levels exist.
      qhat <- q$quantiles[index, stored, drop = FALSE]
    } else {
      residual <- .quantreg.training.residual(o, names(qO)[1L],
                                              source = q[["residual.source"]])
      if (q$method == "forest" && is.null(o$forest.wt)) {
        stop("Forest weights are missing for residual pinball scoring.", call. = FALSE)
      }
      weight <- if (q$method == "forest") o$forest.wt[index, , drop = FALSE] else NULL
      qhat <- .quantreg.residual(residual, q$prediction[index], taus,
                                numeric(0), weight)$quantiles
    }
  } else {
    qhat <- matrix(NA_real_, length(index), length(taus))
    for (i in seq_along(index)) {
      p <- q$cdf[index[i], ]
      ok <- is.finite(p) & is.finite(q$yunq)
      if (sum(ok) == 1L) {
        qhat[i, ] <- q$yunq[ok]
      } else if (sum(ok) > 1L) {
        qhat[i, ] <- stats::approx(p[ok], q$yunq[ok], xout = taus,
                                  ties = "ordered", rule = 2)$y
      }
    }
  }
  y <- y[index]
  out <- vapply(seq_along(taus), function(j) {
    ok <- is.finite(y) & is.finite(qhat[, j])
    if (!any(ok)) return(NA_real_)
    error <- y[ok] - qhat[ok, j]
    mean((taus[j] - (error < 0)) * error)
  }, numeric(1))
  names(out) <- paste0("tau=", taus)
  out
}
