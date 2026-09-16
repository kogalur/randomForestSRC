###################################################################
##
## performance measures (with random-baseline as an attribute)
##
###################################################################
.safe_div <- function(num, den) {
  num <- as.numeric(num); den <- as.numeric(den)
  ifelse(den == 0, NA_real_, num / den)
}
## Convert the supported probability inputs to a two-column matrix ordered
## as (negative, positive). For imbalance summaries, positive is the minority
## class determined from the evaluation response, before excluding scores.
.imbalanced.probability.matrix <- function(prob, n, class.labels, positive,
                                           clamp = TRUE) {
  if (is.data.frame(prob)) {
    prob <- as.data.frame(prob, optional = TRUE)
    numeric.columns <- vapply(prob, function(z) {
      is.numeric(z) || is.logical(z)
    }, logical(1L))
    if (!all(numeric.columns)) {
      stop("'prob' must contain numeric probabilities.", call. = FALSE)
    }
    prob <- as.matrix(prob)
  }
  d <- dim(prob)
  if (is.null(d)) {
    if (!(is.numeric(prob) || is.logical(prob)) || length(prob) != n) {
      stop("'prob' must contain one numeric probability per response.",
           call. = FALSE)
    }
    prob <- matrix(as.numeric(prob), nrow = n, ncol = 1L)
  }
  else if (length(d) != 2L || d[1L] != n || !(d[2L] %in% c(1L, 2L)) ||
           !(is.numeric(prob) || is.logical(prob))) {
    stop("'prob' must be a numeric vector or a one- or two-column ",
         "matrix/data frame with one row per response.", call. = FALSE)
  }
  storage.mode(prob) <- "double"
  prob[!is.finite(prob)] <- NA_real_
  if (clamp) {
    use <- is.finite(prob)
    prob[use] <- pmin(pmax(prob[use], 0), 1)
  }
  if (ncol(prob) == 1L) {
    ## A vector or one column always represents the positive class.
    ## Its name is optional and does not change this convention.
    p <- prob[, 1L]
    prob <- cbind(1 - p, p)
  }
  else {
    labs <- colnames(prob)
    unnamed <- is.null(labs) ||
      (all(!is.na(labs)) && all(!nzchar(labs)))
    if (!unnamed) {
      if (anyNA(labs) || anyDuplicated(labs) ||
          !setequal(labs, class.labels)) {
        stop("The two probability column names must match the response ",
             "factor levels.", call. = FALSE)
      }
      prob <- prob[, match(class.labels, labs), drop = FALSE]
    }
    ## Unnamed columns retain response-factor level order.
    prob <- prob[, c(setdiff(1:2, positive), positive), drop = FALSE]
    if (clamp) {
      complete <- rowSums(is.finite(prob)) == 2L
      if (any(abs(rowSums(prob[complete, , drop = FALSE]) - 1) >
              sqrt(.Machine$double.eps))) {
        stop("Each finite row of two-column 'prob' must sum to one.",
             call. = FALSE)
      }
    }
  }
  colnames(prob) <- c("0", "1")
  prob
}
## Resolve class frequencies from the observed evaluation responses. The
## threshold and random-score reference retain this population even when
## some responses do not have usable probability predictions.
.imbalanced.prepare <- function(yvar, prob) {
  if (!is.factor(yvar) || length(levels(yvar)) != 2L) return(NULL)
  frq <- table(yvar)
  minority <- which.min(frq)
  majority <- setdiff(1:2, minority)
  n <- sum(frq)
  pihat <- if (n > 0) as.numeric(frq[minority]) / n else NA_real_
  y <- as.integer(yvar == levels(yvar)[minority])
  prob.matx <- .imbalanced.probability.matrix(
    prob, length(yvar), levels(yvar), minority)
  ## One common scoring sample; unknown outcomes retain no artificial label.
  use <- !is.na(y) & rowSums(is.finite(prob.matx)) == 2L
  prob.matx <- prob.matx[use, , drop = FALSE]
  list(y = factor(y[use], levels = c(0, 1)),
       prob = prob.matx[, 2L], prob.matx = prob.matx,
       frq = frq, minority = minority, majority = majority,
       n = n, pihat = pihat)
}
## Validate one threshold while retaining the existing finite clamping rule.
.imbalanced.threshold <- function(threshold, default) {
  if (is.null(threshold)) return(default)
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold)) {
    stop("'threshold' must be a single finite number.", call. = FALSE)
  }
  max(0, min(1, threshold))
}
## Harmonic means are zero when all component rates exist and any is zero.
## An undefined component rate remains unavailable, rather than becoming zero.
.imbalanced.harmonic <- function(...) {
  rates <- cbind(...)
  out <- rep(NA_real_, nrow(rates))
  supported <- rowSums(!is.finite(rates) | rates < 0) == 0L
  zero <- rowSums(rates == 0, na.rm = TRUE) > 0L
  out[supported & zero] <- 0
  use <- supported & !zero
  if (any(use)) {
    out[use] <- ncol(rates) / rowSums(1 / rates[use, , drop = FALSE])
  }
  out
}
get.imbalanced.performance <- function(obj,
                                       prob = NULL,
                                       threshold = NULL,
                                       confusion = FALSE,
                                       robust = FALSE) {
  ## determine if a forest object is provided or two vectors (yvar, prob)
  if (is.null(prob)) {
    if (class(obj)[1] != "rfsrc") {
      stop("obj must be a forest object")
    } else {
      yvar <- obj$yvar
      if (!is.null(obj$predicted.oob) && !all(is.na(obj$predicted.oob))) {
        prob <- obj$predicted.oob
      } else {
        prob <- obj$predicted
      }
    }
  } else {
    ## user has supplied (yvar, prob)
    yvar <- obj
  }
  ## if this is not a two-class problem, return NULL
  if (!is.factor(yvar) || length(levels(yvar)) != 2) {
    return(NULL)
  }
  ## Named rules select a cutoff; the workhorse retains the evaluation-based
  ## minority definition and prevalence calculation used by threshold = NULL.
  if (is.character(threshold)) {
    if (length(threshold) != 1L || is.na(threshold) ||
        !(threshold %in% c("auto", "bayes"))) {
      stop("'threshold' must be NULL, \"auto\", \"bayes\", or a single finite number.",
           call. = FALSE)
    }
    rfq <- FALSE
    if (threshold == "auto" && inherits(obj, "rfsrc")) {
      forest <- obj[["forest"]]
      rfq <- is.list(forest) && isTRUE(forest[["rfq"]])
    }
    threshold <- if (rfq) NULL else 0.5
  }
  ## compute
  perf.o <- get.imbalanced.performance.workhorse(
    yvar, prob,
    threshold = threshold,
    confusion = confusion,
    robust = robust
  )
  ## carry the baseline attribute through compact printing
  rand.attr <- attr(perf.o, "rand")
  if (!confusion) {
    v <- unlist(perf.o)        # default behavior preserved
    attr(v, "rand") <- rand.attr
    class(v) <- c("imbalanced.performance", class(v))
    v
  } else {
    attr(perf.o, "rand") <- rand.attr
    class(perf.o) <- c("imbalanced.performance", class(perf.o))
    perf.o
  }
}
get.imbalanced.performance.workhorse <- function (yvar, prob,
                                                  threshold = NULL,
                                                  confusion = FALSE,
                                                  robust = FALSE) {
  dat <- .imbalanced.prepare(yvar, prob)
  if (is.null(dat)) return(NULL)
  y.frq <- dat$frq
  minority <- dat$minority
  majority <- dat$majority
  n.eff <- dat$n
  pihat <- dat$pihat
  iratio <- if (n.eff > 0) max(y.frq) / min(y.frq) else NA_real_
  y <- dat$y
  prob <- dat$prob
  prob.matx <- dat$prob.matx
  threshold <- .imbalanced.threshold(threshold, pihat)
  ## confusion at threshold
  yhat <- factor(1 * (prob >= threshold), levels = c(0, 1))
  confusion.matx <- table(y, yhat)
  N <- sum(confusion.matx)
  if (N > 0) {
    TN <- confusion.matx[1, 1]; FP <- confusion.matx[1, 2]
    FN <- confusion.matx[2, 1]; TP <- confusion.matx[2, 2]
    if (robust) {
      sens <- .safe_div(1 + TP, 1 + FN + TP)
      spec <- .safe_div(1 + TN, 1 + TN + FP)
      prec <- .safe_div(1 + TP, 1 + TP + FP)
      npv  <- .safe_div(1 + TN, 1 + TN + FN)
    } else {
      sens <- .safe_div(TP, FN + TP)
      spec <- .safe_div(TN, TN + FP)
      prec <- .safe_div(TP, TP + FP)
      npv  <- .safe_div(TN, TN + FN)
    }
    misclass   <- (FP + FN) / N
    F1         <- .imbalanced.harmonic(prec, sens)
    F1mod      <- .imbalanced.harmonic(sens, spec, prec, npv)
    gmean      <- sqrt(sens * spec)
    F1gmean    <- (F1 + gmean) / 2
    F1modgmean <- (F1mod + gmean) / 2
  } else {
    sens <- spec <- prec <- npv <- misclass <-
      F1 <- F1mod <- gmean <- F1gmean <- F1modgmean <- NA_real_
  }
  ## Probability scores use the same observed-outcome/finite-prediction rows.
  ## Keep the shared log-loss policy; it is not the rate-smoothing option.
  if (N > 0) {
    brier <- get.brier.error(y, prob.matx, normalized = FALSE)
    brier.norm <- get.brier.error(y, prob.matx)
    auc <- get.auc(y, prob.matx)
    logloss <- get.logloss(y, prob.matx)
    prO <- get.pr.auc(y, prob)
    pr.auc <- prO[1L]
    pr.auc.rand <- prO[2L]
  }
  else {
    brier <- brier.norm <- auc <- logloss <- pr.auc <- pr.auc.rand <- NA_real_
  }
  ## ------------------------------
  ## Analytical random baseline
  ## ------------------------------
  t  <- threshold
  pi <- pihat
  TP.e <- n.eff * pi * (1 - t)
  FN.e <- n.eff * pi * t
  TN.e <- n.eff * (1 - pi) * t
  FP.e <- n.eff * (1 - pi) * (1 - t)
  if (robust) {
    sens.rand <- (1 + TP.e) / (1 + FN.e + TP.e)
    spec.rand <- (1 + TN.e) / (1 + TN.e + FP.e)
    prec.rand <- (1 + TP.e) / (1 + TP.e + FP.e)
    npv.rand  <- (1 + TN.e) / (1 + TN.e + FN.e)
  } else {
    sens.rand <- 1 - t
    spec.rand <- t
    prec.rand <- pi
    npv.rand  <- 1 - pi
  }
  misclass.rand   <- (FP.e + FN.e) / n.eff                     # == pi*t + (1-pi)*(1-t)
  F1.rand         <- ifelse((prec.rand + sens.rand) > 0,
                            2 * (prec.rand * sens.rand) / (prec.rand + sens.rand), NA_real_)
  F1mod.rand      <- ifelse(all(is.finite(c(sens.rand, spec.rand, prec.rand, npv.rand))) &&
                            all(c(sens.rand, spec.rand, prec.rand, npv.rand) > 0),
                            4 / (1/sens.rand + 1/spec.rand + 1/prec.rand + 1/npv.rand), NA_real_)
  gmean.rand      <- sqrt(sens.rand * spec.rand)
  F1gmean.rand    <- (F1.rand + gmean.rand) / 2
  F1modgmean.rand <- (F1mod.rand + gmean.rand) / 2
  auc.rand        <- 0.5
  ## Brier baselines aligned to your get.brier.error():
  ##   - unnormalized (J=2): E[(p - Y)^2] = 1/3
  ##   - normalized: factor is 4x the unnormalized in binary case
  brier.rand      <- 1/3
  brier.norm.rand <- 4 * brier.rand
  logloss.rand    <- 1.0
  ## assemble primary output (unchanged fields)
  rO <- list(n.majority = as.numeric(y.frq[majority]),
             n.minority = as.numeric(y.frq[minority]),
             iratio = iratio,
             threshold = threshold,
             sens = sens,
             spec = spec,
             prec = prec,
             npv = npv,
             misclass = misclass,
             brier = brier,
             brier.norm = brier.norm,
             auc = auc,
             logloss = logloss,
             F1 = F1,
             F1mod = F1mod,
             pr.auc.rand = pr.auc.rand,
             pr.auc = pr.auc,
             F1gmean = F1gmean,
             F1modgmean = F1modgmean,
             gmean = gmean)
  if (confusion) {
    class.error <- 1 - diag(confusion.matx) / rowSums(confusion.matx, na.rm = TRUE)
    rO$confusion <- cbind(confusion.matx, class.error = round(class.error, 4))
  }
  ## attach random-baseline as an attribute so default print stays compact
  rand.metrics <- c(
    sens.rand = sens.rand,
    spec.rand = spec.rand,
    prec.rand = prec.rand,
    npv.rand  = npv.rand,
    misclass.rand = misclass.rand,
    F1.rand = F1.rand,
    F1mod.rand = F1mod.rand,
    F1gmean.rand = F1gmean.rand,
    F1modgmean.rand = F1modgmean.rand,
    gmean.rand = gmean.rand,
    auc.rand = auc.rand,
    brier.rand = brier.rand,
    brier.norm.rand = brier.norm.rand,
    logloss.rand = logloss.rand
  )
  confusion.rand.expected <- matrix(c(TN.e, FP.e, FN.e, TP.e), nrow = 2, byrow = TRUE)
  dimnames(confusion.rand.expected) <- list(Truth = levels(y), Pred = levels(y))
  attr(rO, "rand") <- list(
    metrics = rand.metrics,
    confusion.expected = confusion.rand.expected
  )
  rO
}
# pretty printer for imbalanced performance (model vs random baseline)
print.imbalanced.performance <- function(x, digits = 4, show.confusion = TRUE, ...) {
  # helpers to extract by name from list or named numeric
  .has <- function(nm) nm %in% names(x)
  .get <- function(nm) {
    if (.has(nm)) {
      if (is.list(x)) x[[nm]] else as.numeric(x[nm])
    } else NA_real_
  }
  .fmt <- function(z) formatC(z, digits = digits, format = "f")
  .fmt_pct <- function(z) ifelse(is.na(z), NA, paste0(.fmt(z), " %"))
  rand <- attr(x, "rand")
  rand.metrics <- if (!is.null(rand)) rand$metrics else NULL
  ## Observed prevalence uses response counts, not smoothed reference precision.
  pi.hat <- .safe_div(.get("n.minority"),
                      .get("n.majority") + .get("n.minority"))
  # header ----------------------------------------------------------
  cat("Imbalanced classification performance (model vs random baseline)\n")
  cat(sprintf(" n.majority = %s   n.minority = %s   iratio = %s   threshold = %s   prevalence  = %s\n\n",
              .get("n.majority"),
              .get("n.minority"),
              .fmt(.get("iratio")),
              .fmt(.get("threshold")),
              .fmt(pi.hat)))
  ## Metric direction is independent of the table in which a row appears.
  ## +1 means larger is better; -1 means smaller is better.
  direction <- c(sens = 1, spec = 1, prec = 1, npv = 1,
                 F1 = 1, F1mod = 1, gmean = 1,
                 F1gmean = 1, F1modgmean = 1,
                 auc = 1, pr.auc = 1,
                 misclass = -1, brier = -1, brier.norm = -1, logloss = -1)
  all_metrics <- names(direction)
  ## Display G-mean in both tables, but calculate its comparison only once.
  rate_metrics <- c("sens", "spec", "prec", "npv", "F1", "F1mod",
                    "gmean", "F1gmean", "F1modgmean")
  overall_metrics <- c("gmean", "auc", "pr.auc", "misclass",
                       "brier", "brier.norm", "logloss")
  # mapping to baseline names
  base_map <- stats::setNames(paste0(all_metrics, ".rand"), all_metrics)
  base_map["pr.auc"] <- "pr.auc.rand"  # special case lives in main list
  # build table -----------------------------------------------------
  model_vals <- sapply(all_metrics, .get)
  base_vals  <- sapply(all_metrics, function(m) {
    nm <- base_map[[m]]
    if (m == "pr.auc") {
      .get(nm)                       # "pr.auc.rand" is in the main object
    } else if (!is.null(rand.metrics) && nm %in% names(rand.metrics)) {
      as.numeric(rand.metrics[[nm]]) # most baselines live in the attribute
    } else {
      NA_real_
    }
  })
  ## Positive Delta and Gain always mean improvement over the baseline.
  ## Preserve the existing handling of zero and near-zero baselines.
  delta <- rep(NA_real_, length(all_metrics)); names(delta) <- all_metrics
  gain  <- delta
  for (m in all_metrics) {
    mv <- model_vals[[m]]; bv <- base_vals[[m]]
    if (is.na(mv) || is.na(bv)) next
    if (direction[[m]] < 0) {
      delta[m] <- bv - mv
      gain[m]  <- if (bv != 0) (bv - mv) / bv * 100 else NA_real_
    } else {
      delta[m] <- mv - bv
      gain[m]  <- if (abs(bv) > .Machine$double.eps) (mv - bv) / abs(bv) * 100 else NA_real_
    }
  }
  df <- data.frame(
    Metric   = all_metrics,
    Model    = .fmt(model_vals[all_metrics]),
    Baseline = .fmt(base_vals[all_metrics]),
    Delta    = .fmt(delta[all_metrics]),
    Gain     = .fmt_pct(gain[all_metrics]),
    row.names = NULL,
    check.names = FALSE
  )
  ## Reuse the same formatted rows so G-mean agrees exactly across tables.
  cat("Classification rates and combined measures\n")
  print(df[match(rate_metrics, df$Metric), , drop = FALSE],
        row.names = FALSE, right = TRUE)
  cat("\nOverall performance\n")
  print(df[match(overall_metrics, df$Metric), , drop = FALSE],
        row.names = FALSE, right = TRUE)
  cat("\nPositive Delta and Gain indicate improvement over the random baseline.\n")
  # confusion matrices ---------------------------------------------
  if (show.confusion && is.list(x) && !is.null(x$confusion)) {
    cat("\nConfusion matrix (model):\n")
    print(x$confusion)
    ce <- if (!is.null(rand) && !is.null(rand$confusion.expected)) rand$confusion.expected else NULL
    if (!is.null(ce)) {
      cat("\nExpected confusion (random baseline):\n")
      ce_err <- c(
        ifelse(sum(ce[1,]) > 0, 1 - ce[1,1]/sum(ce[1,]), NA_real_),
        ifelse(sum(ce[2,]) > 0, 1 - ce[2,2]/sum(ce[2,]), NA_real_)
      )
      ce_out <- cbind(round(ce, 2), class.error = round(ce_err, 4))
      print(ce_out)
    }
  }
  invisible(x)
}
###################################################################
##
##
## rfq threshold
##
##
###################################################################
get.rfq.threshold <- function(y) {
  frq <- table(y)
  if (length(frq) != 2) {
    return(NULL)
  }
  as.numeric(min(frq, na.rm = TRUE) / sum(frq, na.rm = TRUE))
}
###################################################################
##
##
## optimize the gmean (or balanced measure) by the threshold
##
##
###################################################################
# internal: fast, vectorized sweep over thresholds (avoids mclapply)
.get.imbalanced.optimize.sweep <- function(yvar, prob, thresholds) {
  dat <- .imbalanced.prepare(yvar, prob)
  if (is.null(dat)) return(NULL)
  y <- as.integer(as.character(dat$y))
  p <- dat$prob
  if (length(p) == 0L) {
    return(data.frame(
      threshold = thresholds,
      gmean = rep(NA_real_, length(thresholds)),
      F1 = rep(NA_real_, length(thresholds)),
      F1mod = rep(NA_real_, length(thresholds)),
      F1modgmean = rep(NA_real_, length(thresholds))
    ))
  }
  ## sort by probability (ascending)
  o <- order(p)
  p <- p[o]
  y <- y[o]
  ## prefix sums (length n+1 to handle k=0)
  n <- length(p)
  cum_pos <- c(0L, cumsum(y))
  cum_neg <- c(0L, cumsum(1L - y))
  total_pos <- cum_pos[n + 1L]
  total_neg <- cum_neg[n + 1L]
  ## for each threshold t: k = number of probs < t  (since prediction is prob >= t)
  k <- findInterval(thresholds, p, left.open = TRUE)
  idx <- k + 1L
  FN <- cum_pos[idx]
  TN <- cum_neg[idx]
  TP <- total_pos - FN
  FP <- total_neg - TN
  sens <- .safe_div(TP, FN + TP)
  spec <- .safe_div(TN, TN + FP)
  prec <- .safe_div(TP, TP + FP)
  npv  <- .safe_div(TN, TN + FN)
  F1 <- .imbalanced.harmonic(prec, sens)
  F1mod <- .imbalanced.harmonic(sens, spec, prec, npv)
  gmean <- sqrt(sens * spec)
  F1modgmean <- (F1mod + gmean) / 2
  data.frame(
    threshold = thresholds,
    gmean = gmean,
    F1 = F1,
    F1mod = F1mod,
    F1modgmean = F1modgmean
  )
}
get.imbalanced.optimize <- function(obj,
                                    prob = NULL,
                                    newdata = NULL,
                                    measure = c("gmean", "F1", "F1mod", "F1modgmean"),
                                    ngrid = 1000,
                                    plot.it = TRUE) {
  ## determine if a forest object is provided or two vectors (yvar, prob)
  if (is.null(prob)) {
    if (class(obj)[1] != "rfsrc") {
      stop("obj must be a forest object")
    }
    ## object is a forest ---> parse for (yvar, prob)
    else {
      ## if newdata is provided, we swap the grow object with the predict object
      if (!is.null(newdata)) {
        obj <- predict(obj, newdata)
      }
      yvar <- obj$yvar
      if (!is.null(obj$predicted.oob) && !all(is.na(obj$predicted.oob))) {
        prob <- obj$predicted.oob
      }
      else {
        prob <- obj$predicted
      }
    }
  }
  ## user has supplied (yvar, prob)
  else {
    yvar <- obj
  }
  ## check that measure requested is coherent
  measure <- match.arg(measure, c("gmean", "F1", "F1mod", "F1modgmean"))
  if (!is.numeric(ngrid) || length(ngrid) != 1L || !is.finite(ngrid) ||
      ngrid < 2 || ngrid != floor(ngrid) || ngrid > .Machine$integer.max) {
    stop("'ngrid' must be an integer of at least two.", call. = FALSE)
  }
  thresholds <- seq(0, 1, length.out = as.integer(ngrid))
  ## fast, vectorized threshold sweep (avoids mclapply and is usually faster)
  sweep <- .get.imbalanced.optimize.sweep(yvar, prob, thresholds)
  if (is.null(sweep)) {
    return(NULL)
  }
  available <- which(is.finite(sweep[[measure]]))
  if (length(available) == 0L) {
    stop("No finite '", measure, "' values are available for threshold ",
         "selection. Check the observed responses and probability predictions, ",
         "or select another measure.", call. = FALSE)
  }
  ## Preserve the first maximum, including a valid optimum at zero.
  best <- available[which.max(sweep[[measure]][available])]
  best.threshold <- sweep$threshold[best]
  if (plot.it) {
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par), add = TRUE)
    par(mfrow = c(2, 2))
    pt <- sweep$threshold < 3 * best.threshold
    if (sum(pt) < 2L) pt <- rep(TRUE, nrow(sweep))
    for (metric in c("gmean", "F1", "F1mod", "F1modgmean")) {
      value <- sweep[[metric]][pt]
      finite <- is.finite(value)
      limits <- if (any(finite)) range(value[finite]) else c(0, 1)
      plot(sweep$threshold[pt], value,
           xlab = "threshold", ylab = metric, ylim = limits)
      abline(v = best.threshold, lty = 2, col = "blue")
    }
  }
  ## return full performance metrics at the optimal threshold (same as previous behavior)
  perf <- get.imbalanced.performance(yvar, prob, threshold = best.threshold)
  if (is.null(perf)) {
    return(NULL)
  }
  out <- data.frame(as.list(perf))
  rownames(out) <- best
  out
}
###################################################################
##
##
## precision recall utilities (these are internal functions)
## assumes y=0--> majority, y=1--> minority
##
##
###################################################################
## Direct PR callers may supply numeric score vectors, not just probabilities.
## Preserve their scale and the positive=1 convention for already coded truth.
## For other two-level factors, resolve the minority before recoding labels.
.imbalanced.pr.input <- function(truth, yhat) {
  if (is.factor(truth)) {
    if (length(levels(truth)) != 2L) return(NULL)
    labs <- levels(truth)
    positive <- if (setequal(labs, c("0", "1"))) {
      match("1", labs)
    }
    else {
      which.min(table(truth))
    }
    y <- as.integer(truth == labs[positive])
  }
  else {
    if (!is.atomic(truth) || !is.null(dim(truth)) ||
        !all(is.na(truth) | truth %in% c(0, 1))) return(NULL)
    labs <- c("0", "1")
    positive <- 2L
    y <- as.integer(truth)
  }
  scores <- .imbalanced.probability.matrix(
    yhat, length(truth), labs, positive, clamp = FALSE)[, 2L]
  use <- !is.na(y) & is.finite(scores)
  list(positive = scores[use & y == 1L],
       negative = scores[use & y == 0L])
}
get.pr.auc <- function(truth, yhat) {
  dat <- .imbalanced.pr.input(truth, yhat)
  if (is.null(dat)) return(c(NA_real_, NA_real_))
  x <- dat$positive
  y <- dat$negative
  rO <- c(NA, NA)
  if (length(x) > 0 && length(y) > 0) {
    ## precision recall with random baseline
    pr.o <- tryCatch({get.pr.workhorse(x, y, rand.compute = TRUE)}, error = function(ex) {NULL})
    if (!is.null(pr.o)) {
      rO <- c(pr.o$auc.integral, pr.o$auc.integral.rand)
    }
  }
  rO
}
get.pr.curve <- function(truth, yhat) {
  dat <- .imbalanced.pr.input(truth, yhat)
  if (is.null(dat)) return(NULL)
  x <- dat$positive
  y <- dat$negative
  rO <- NULL
  if (length(x) > 0 && length(y) > 0) {
    ## precision recall with random baseline
    pr.o <- tryCatch({get.pr.workhorse(x, y, curve = TRUE)}, error = function(ex) {NULL})
    if (!is.null(pr.o)) {
      rO <- pr.o$curve
    }
  }
  rO
}
## copied from PRROC library - we gratefully acknowledge the developer for this function
get.pr.workhorse <- function(scores.class0, scores.class1 = scores.class0, weights.class0 = NULL,
                     weights.class1 = {if (is.null(weights.class0)) {NULL} else {1 - weights.class0}},
                     sorted = FALSE,
                     minStepSize = min(1, ifelse(is.null(weights.class0), 1, sum(weights.class0)/100)),
                     rand.compute = FALSE, curve = FALSE) 
{
  ## preliminary processing
  if (!sorted) {
    o0 <- order(scores.class0)
    scores.class0 <- scores.class0[o0]
    if (!is.null(weights.class0)) {
      weights.class0 <- weights.class0[o0]
    }
    o1 <- order(scores.class1)
    scores.class1 <- scores.class1[o1]
    if (!is.null(weights.class1)) {
      weights.class1 <- weights.class1[o1]
    }
  }
  sorted.scores.class0 <- scores.class0
  sorted.scores.class1 <- scores.class1
  ## main workhorse
  weights.class0 <- c(rep(1, length(sorted.scores.class0)), 
                      rep(0, length(sorted.scores.class1)))
  sorted.scores.class0 <- c(sorted.scores.class0, sorted.scores.class1)
  o0 <- order(sorted.scores.class0)
  sorted.scores.class0 <- sorted.scores.class0[o0]
  weights.class0 <- weights.class0[o0]
  weights.class1 <- 1 - weights.class0
  sorted.scores.class1 <- sorted.scores.class0
  all.scores <- sorted.scores.class0
  all.weights.pos <- weights.class0
  all.weights.neg <- weights.class1
  o <- order(all.scores, decreasing = T)
  all.scores <- all.scores[o]
  all.weights.pos <- all.weights.pos[o]
  all.weights.neg <- all.weights.neg[o]
  cum.weights.pos <- cumsum(all.weights.pos)
  cum.weights.neg <- cumsum(all.weights.neg)
  cum.use <- c(all.scores[-length(all.scores)] != all.scores[-1], TRUE)
  all.scores <- all.scores[cum.use]
  cum.weights.pos <- cum.weights.pos[cum.use]
  cum.weights.neg <- cum.weights.neg[cum.use]
  r.fg <- sum(all.weights.pos)
  tp <- cum.weights.pos
  fp <- cum.weights.neg
  tp.prev <- c(0, cum.weights.pos[-length(cum.weights.pos)])
  fp.prev <- c(0, cum.weights.neg[-length(cum.weights.neg)])
  h <- (fp - fp.prev)/(tp - tp.prev)
  a <- 1 + h
  b <- (fp.prev - h * tp.prev)/r.fg
  h[tp == tp.prev] <- 1
    a[tp == tp.prev] <- 1
  b[tp == tp.prev] <- 0
  v <- (tp/r.fg - tp.prev/r.fg - b/a * (log(a * tp/r.fg + b) - log(a * tp.prev/r.fg + b)))/a
  v2 <- (tp/r.fg - tp.prev/r.fg)/a
  v[b == 0] <- v2[b == 0]
  vals <- v
  auc.integral <- sum(vals)
  res <- list(auc.integral = auc.integral)
  rand.auc <- pr.curve <- NULL
  ## random auc
  if (rand.compute) {
    if (is.null(weights.class0)) {
      rand.auc <- length(sorted.scores.class0)/(length(sorted.scores.class0) + 
                                                length(sorted.scores.class1))
    }
    else {
      rand.auc <- sum(weights.class0)/sum(weights.class0 + weights.class1)
    }
  }
  ## pr curve
  if (curve) {
    minStepSize.2 <- minStepSize/r.fg
    min.curve <- cbind(tp/r.fg, tp/(tp + fp), all.scores)
    idxs <- which((tp - tp.prev)/r.fg > minStepSize.2 & tp/(tp + fp) != tp.prev/(tp.prev + fp.prev))
    idxs <- idxs[idxs > 1]
    if (length(idxs) > 0) {
      m <- sapply(idxs, function(i) {
        x <- seq(0, min.curve[i, 1] - min.curve[i - 1, 1], by = minStepSize.2)
        sns <- min.curve[i - 1, 1] + x
        prcs <- (min.curve[i - 1, 1] + x)/(min.curve[i - 1, 1] + x + fp[i - 1]/r.fg + h[i] * x)
        temp <- rbind(sns, prcs, rep(all.scores[i], length(x)))
        temp
      })
      m <- matrix(unlist(m), ncol = 3, byrow = T)
      m <- rbind(min.curve, m)
    }
    else {
      m <- min.curve
    }
    m <- m[order(m[, 1], -m[, 3], decreasing = T), , drop = FALSE]
    m <- rbind(c(1, m[1, 2:3]), m, c(0, m[nrow(m), 2:3]))
    dimnames(m) <- c(NULL, NULL)
    colnames(m) <- c("recall", "precision", "threshold")
    pr.curve <- m
  }
  c(res, list(auc.integral.rand = rand.auc, curve = pr.curve))
}
###################################################################
##
##
## one versus one  OVO
##
##
###################################################################
# Function to train one-vs-one classifiers using randomForestSRC's imbalanced
train_ovo <- function(formula, data, verbose = TRUE) {
  label_col <- all.vars(formula)[1]
  class_pairs <- combn(as.character(unique(data[[label_col]])), 2, simplify = FALSE)
  classifiers <- list()
  cat("training one-vs-one\n")
  for (pair in class_pairs) {
    # Verbose output
    if (verbose) print(pair)
    # Filter data for the current pair of classes
    binary_data <- data[data[[label_col]] %in% pair, ]
    binary_data[[label_col]] <- droplevels(binary_data[[label_col]])
    # No training if two classes are not present
    if (length(levels(binary_data[, label_col])) < 2) {
      model <- NULL
    }
    # Train random forest classifier with imbalanced function
    else {
      model <- tryCatch({suppressWarnings(imbalanced(formula, data = binary_data))},
                      error = function(ex) {NULL})
    }
    # Assign the classifier
    classifiers[[paste(pair, collapse = "_")]] <- model
  }
  classifiers <- classifiers[lengths(classifiers) != 0]
  attr(classifiers, "formula") <- formula
  attr(classifiers, "data") <- data
  return(classifiers)
}
train.ovo <- train_ovo
# Function to predict using one-vs-one classifiers
predict_ovo <- function(classifiers, data = NULL, verbose = TRUE) {
  ## true testing?
  test <- TRUE
  if (is.null(data)) {
    data <- attr(classifiers, "data")
    test <- FALSE
  }
  ## set up the vote matrix
  formula <- attr(classifiers, "formula")
  label_col <- all.vars(formula)[1]
  votes <- matrix(0, nrow = nrow(data), ncol = length(unique(data[[label_col]])))
  colnames(votes) <- levels(droplevels(data[[label_col]]))
  cat("prediction one-vs-one\n")
  for (pair in names(classifiers)) {
    # Split the pair name to get class labels
    classes <- unlist(strsplit(pair, "_"))
    # Verbose output
    if (verbose) print(classes)
    # Get the relevant model
    model <- classifiers[[pair]]
    # Get predictions for the entire data
    full_predictions <- predict(model, data[, model$xvar.names, drop = FALSE])$class
    # Use OOB predictions for the training data (this is not a true train/test scenario)
    if (!test) {
      binary_data_indices <- which(data[[label_col]] %in% classes)
      full_predictions[binary_data_indices] <- model$class.oob
    }
    # Assign votes based on predictions
    for (i in 1:nrow(data)) {
      cl <- which(colnames(votes) == full_predictions[i])
      if (length(cl) > 0) {
        votes[i, cl] <- votes[i, cl] + 1
      }
    }
  }
  # Determine final class by majority vote
  final_predictions <- apply(votes, 1, function(vote) colnames(votes)[which.max(vote)])
  return(final_predictions)
}
## predict.ovo <- predict_ovo
# Function to evaluate the model and calculate AUC and PR-AUC
one.versus.one <- function(formula, data, verbose = TRUE) {
  # Clean the data
  o <- rfsrc.cart(formula, data, nodedepth = 0, perf.type = "none")
  if (o$family != "class") {
    stop("this function only works for classification families\n")
  }
  data <- droplevels(data.frame(o$yvar, o$xvar))
  colnames(data) <- c(o$yvar.names, o$xvar.names)
  # Train classifiers
  classifiers <- train_ovo(formula, data, verbose)
  attr(classifiers, "formula") <- formula
  # Predict on the entire data
  predictions <- predict_ovo(classifiers, verbose = verbose)
  # Evaluate the predictions
  confusionMatrixCustom(predictions, o$yvar)
}
ovo <- one.versus.one
confusionMatrixCustom <- function(predictions, actuals) {
  # Ensure that predictions and actuals are factors
  predictions <- factor(predictions, levels = levels(actuals))
  # Create the confusion matrix
  confusion <- table(Predicted = predictions, Actual = actuals)
  # Calculate overall accuracy
  accuracy <- sum(diag(confusion)) / sum(confusion)
  n <- sum(confusion)
  # Calculate the 95% confidence interval for accuracy
  se <- sqrt(accuracy * (1 - accuracy) / n)
  ci <- c(accuracy - qnorm(0.975) * se, accuracy + qnorm(0.975) * se)
  # Calculate the no information rate (NIR)
  nir <- max(table(actuals) / n)
  # Calculate the p-value for accuracy greater than NIR
  p_value <- binom.test(sum(diag(confusion)), n, p = nir, alternative = "greater")$p.value
  # Initialize metrics
  sensitivity <- numeric(ncol(confusion))
  specificity <- numeric(ncol(confusion))
  pos_pred_value <- numeric(ncol(confusion))
  neg_pred_value <- numeric(ncol(confusion))
  prevalence <- numeric(ncol(confusion))
  detection_rate <- numeric(ncol(confusion))
  detection_prevalence <- numeric(ncol(confusion))
  balanced_accuracy <- numeric(ncol(confusion))
  for (i in 1:ncol(confusion)) {
    sensitivity[i] <- confusion[i, i] / sum(confusion[, i])
    specificity[i] <- sum(confusion[-i, -i]) / sum(confusion[, -i])
    pos_pred_value[i] <- confusion[i, i] / sum(confusion[i, ])
    neg_pred_value[i] <- sum(confusion[-i, -i]) / sum(confusion[-i, ])
    prevalence[i] <- sum(confusion[, i]) / n
    detection_rate[i] <- confusion[i, i] / n
    detection_prevalence[i] <- sum(confusion[i, ]) / n
    balanced_accuracy[i] <- (sensitivity[i] + specificity[i]) / 2
  }
  # Compile results into a list
  results <- list(
    confusion = confusion,
    overall = c(
      Accuracy = accuracy,
      "95% CI" = paste0("(", round(ci[1], 4), ", ", round(ci[2], 4), ")"),
      "No Information Rate" = nir,
      "P-Value [Acc > NIR]" = format_p_value(p_value)
    ),
    byClass = data.frame(rbind(
      Sensitivity = sensitivity,
      Specificity = specificity,
      Pos_Pred_Value = pos_pred_value,
      Neg_Pred_Value = neg_pred_value,
      Prevalence = prevalence,
      Detection_Rate = detection_rate,
      Detection_Prevalence = detection_prevalence,
      Balanced_Accuracy = balanced_accuracy
    ))
  )
  colnames(results$byClass) <- levels(actuals)
  # Print results
  cat("Overall Statistics\n")
  cat("\n")
  cat(sprintf("               Accuracy : %.4s\n", results$overall["Accuracy"]))
  cat(sprintf("                 95%% CI : %s\n", results$overall["95% CI"]))
  cat(sprintf("    No Information Rate : %.4s\n", results$overall["No Information Rate"]))
  cat(sprintf("    P-Value [Acc > NIR] : %s\n", results$overall["P-Value [Acc > NIR]"]))
  cat("\n")
  cat("Confusion Matrix\n")
  print(results$confusion)
  cat("\n")
  cat("Class Statistics\n")
  print(results$byClass)
  # Return
  invisible(results)
}
# Helper function to format p-values
format_p_value <- function(p_value) {
  if (p_value < 2.2e-16) {
    return("< 2.2e-16")
  } else {
    return(sprintf("%.4g", p_value))
  }
}
