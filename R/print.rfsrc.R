print.rfsrc <- function(x, outcome.target = NULL, ...) {
  has.all.classes <- function(obj, cls) all(cls %in% class(obj))
  ## simple pass-throughs for certain composite objects
  if (has.all.classes(x, c("rfsrc","forest")) ||
      has.all.classes(x, c("rfsrc","plot.variable")) ||
      has.all.classes(x, c("rfsrc","partial")) ||
      has.all.classes(x, c("rfsrc","sidClustering"))) {
    print.default(x)
    return()
  }
  ## mean pinball losses from a quantreg object
  .get.pinball.losses <- function(qobj, taus = c(0.1, 0.5, 0.9), subset = NULL) {
    qL <- extract.quantile(qobj)
    j <- 1L
    q <- qL[[j]]
    y <- if (is.vector(qobj$yvar)) c(qobj$yvar) else qobj$yvar[, names(qL)[j]]
    if (is.null(subset)) subset <- seq_along(y)
    y    <- y[subset]
    cdf  <- q$cdf[subset, , drop = FALSE]
    yunq <- q$yunq
    pinball <- function(y, q, tau) (tau - (y < q)) * (y - q)
    out <- sapply(taus, function(tt) {
      ## invert CDF to get q_tau by monotone interpolation
      qhat <- vapply(seq_len(nrow(cdf)), function(i) {
        suppressWarnings(stats::approx(x = cdf[i, ], y = yunq, xout = tt,
                                ties = "ordered", rule = 2)$y)
      }, numeric(1))
      mean(pinball(y, qhat, tt), na.rm = TRUE)
    })
    names(out) <- paste0("tau=", taus)
    out
  }
  ## subsample/bootsample
  if (has.all.classes(x, c("rfsrc","subsample"))) {
    print.subsample(x, ...)
    return()
  }
  if (has.all.classes(x, c("rfsrc","bootsample"))) {
    print.bootsample(x, ...)
    return()
  }
  ## accepted classes: (rfsrc,grow) or (rfsrc,predict)
  if (!has.all.classes(x, c("rfsrc","grow")) && !has.all.classes(x, c("rfsrc","predict"))) {
    stop("This function only works for objects of class '(rfsrc, grow)' or '(rfsrc, predict)'.")
  }
  grow.mode <- has.all.classes(x, c("rfsrc","grow"))
  ## sampling type label
  if (!is.null(x$forest) && !is.null(x$forest$bootstrap) && x$forest$bootstrap == "by.root") {
    samp.used <- x$forest$samptype
  } else {
    samp.used <- if (!is.null(x$forest)) x$forest$bootstrap else NA
  }
  ## save before coercion
  family.pretty <- family.pretty(x)
  family.org <- x$family
  yvar.dim <- if (!is.null(x$yvar)) ncol(as.data.frame(x$yvar)) else NA_integer_
  ## coerce multivariate -> univariate if requested/needed
  ## save standardized mv error rates if available
  x$univariate <- coerce.multivariate(x, 1)$univariate
  if (is.null(outcome.target) && !x$univariate) {
    mv.err.rate <- as.numeric(get.mv.error(x, FALSE))
    mv.err.rate <- c(mean(mv.err.rate, na.rm=TRUE), mv.err.rate)
  }
  if (!is.null(outcome.target) && !x$univariate) {
    outcome.target <- get.univariate.target(x, outcome.target)
    x <- coerce.multivariate(x, outcome.target)
  }
  ## ---- Metrics containers
  conf.mat <- miss.err <- iratio <- brier.err <- brier.norm.err <- auc.err <- pr.auc.err <- gmean.err <- NULL
  logloss.err <- NULL
  crps.err <- crps.std.err <- NULL
  event.freq.txt <- NULL
  k.class <- NULL
  brier.rand <- brier.norm.rand <- logloss.rand <- NULL
  q.crps <- q.crps.std <- NULL
  q.pinball <- NULL
  ## ---- Classification block
  if (x$family == "class") {
    if (!is.null(x$yvar)) {
      tab <- table(x$yvar)
      event.freq.txt <- paste(paste(names(tab), as.integer(tab), sep = "="), collapse = ", ")
      k.class <- length(tab)
    }
    ## choose predicted probs / classes (prefer OOB)
    prob.mat <- if (!is.null(x$predicted.oob) && !all(is.na(x$predicted.oob))) x$predicted.oob else x$predicted
    class.pred <- if (!is.null(x$class.oob) && !all(is.na(x$class.oob))) x$class.oob else x$class
    if (!is.null(x$err.rate) && !is.null(x$yvar) && !is.null(prob.mat)) {
      conf.mat <- get.confusion(x$yvar, class.pred)
      names(dimnames(conf.mat)) <- c("  observed", "predicted")
      ## robust misclassification: use first K columns
      if (!is.null(k.class) && k.class >= 2) {
        core <- conf.mat[, seq_len(min(k.class, ncol(conf.mat))), drop = FALSE]
        miss.err <- 1 - sum(diag(core)) / sum(core)
      }
      if (k.class > 2) {
        brier.err <- get.brier.error(x$yvar, prob.mat, normalized = FALSE)
        brier.norm.err <- get.brier.error(x$yvar, prob.mat, normalized = TRUE)
        auc.err <- get.auc(x$yvar, prob.mat)
        logloss.err <- get.logloss(x$yvar, prob.mat)
        iratio <- pr.auc.err <- gmean.err <- NULL
      } else if (!is.null(k.class) && k.class == 2) {
        thr <- if (!is.null(x$forest$rfq) && isTRUE(x$forest$rfq)) get.rfq.threshold(x$forest$yvar) else 0.5
        perO <- get.imbalanced.performance(x$yvar, prob.mat, threshold = thr, confusion = TRUE)
        iratio <- perO$iratio
        brier.err <- perO$brier
        brier.norm.err <- perO$brier.norm
        auc.err <- perO$auc
        logloss.err <- perO$logloss
        pr.auc.err <- perO$pr.auc
        gmean.err <- perO$gmean
        ## if you want RFQ confusion, you could set conf.mat <- perO$confusion here
      }
      ## --- Random-classifier strawman (uniform probabilities) ---
      if (!is.null(k.class) && k.class >= 2) {
        lev <- colnames(prob.mat)
        if (is.null(lev) && !is.null(levels(x$yvar))) lev <- levels(x$yvar)
        P0 <- matrix(1 / k.class, nrow = length(x$yvar), ncol = k.class,
                     dimnames = list(NULL, lev))
        brier.rand <- get.brier.error(x$yvar, P0, normalized = FALSE)
        brier.norm.rand <- get.brier.error(x$yvar, P0, normalized = TRUE)
        logloss.rand <- get.logloss(x$yvar, P0)
      }
    }
  }
  ## ---- Survival / CR block
  if (grepl("surv", x$family)) {
    if (!is.null(x$err.rate) && x$family == "surv") {
      bs <- get.brier.survival(x)
      crps.err <- bs$crps
      crps.std.err <- bs$crps.std
    }
    event <- x$event.info$event
    if (!is.null(event)) {
      tab <- table(event)
      event.freq.txt <- if (length(tab)) paste(paste(names(tab), as.integer(tab), sep = "="), collapse = ", ") else "0"
    }
  }
  ## ---- Quantile regression (pinball/CRPS) block
  quantile.flag <- has.all.classes(x, c("rfsrc","quantreg")) && identical(x$family, "regr")
  if (quantile.flag) {
    if (!is.null(x$yvar)) {
      ## CRPS (standardized and raw) from existing internal helper
      qcrps.std.df <- tryCatch(get.quantile.crps(x, pretty = TRUE, standardize = TRUE),
                               error = function(e) NULL)
      qcrps.raw.df <- tryCatch(get.quantile.crps(x, pretty = TRUE, standardize = FALSE),
                               error = function(e) NULL)
      if (is.data.frame(qcrps.std.df) && nrow(qcrps.std.df))
        q.crps.std <- tail(qcrps.std.df$crps, 1L)
      if (is.data.frame(qcrps.raw.df) && nrow(qcrps.raw.df))
        q.crps      <- tail(qcrps.raw.df$crps, 1L)
      ## Pinball loss at user-selectable taus (defaults 0.1, 0.5, 0.9)
      taus <- getOption("rfsrc.pinball.taus", c(0.1, 0.5, 0.9))
      q.pinball <- tryCatch(.get.pinball.losses(x, taus = taus), error = function(e) NULL)
    }
  }
  ## ---- Requested performance error (last row)
  err.rate <- NULL
  r.sq <- NULL
  if (!is.null(x$err.rate)) {
    er <- cbind(x$err.rate)
    if (grepl("surv", x$family)) {
      err.rate <- digits.pretty(er[nrow(er), , drop = TRUE], 8)
    } else if (x$family == "class") {
      if ((grow.mode && x$forest$perf.type == "gmean") || (!grow.mode && x$perf.type == "gmean")) {
        err.rate <- digits.pretty(er[nrow(er), 1], 8)
      } else {
        err.rate <- digits.pretty(er[nrow(er), , drop = TRUE], 8)
      }
    } else if (x$family == "regr") {
      if (!is.null(x$yvar)) {
        r.sq <- 1 - er[nrow(er), ] / var(x$yvar, na.rm = TRUE)
      }
      err.rate <- digits.pretty(er[nrow(er), , drop = TRUE], 8)
    }
  }
  ## adjustment for multivariate families - swap the standardized mv terror
  if (is.null(outcome.target) && !x$univariate) {
    err.rate <- digits.pretty(mv.err.rate, 3)
  }
  ## ensure nsplit exists
  if (is.null(x$nsplit)) x$nsplit <- 0
  ## helper for sampsize printing
  sampsize.pretty <- function(sf, nval) {
    out <- tryCatch({
      if (is.function(sf)) round(sf(nval)) else as.numeric(sf)
    }, error = function(e) NA_real_)
    out
  }
  ## ---------------- GROW MODE ----------------
  if (grow.mode) {
    cat("                         Sample size: ", x$n, "\n", sep = "")
    if (grepl("surv", x$family)) {
      if (!is.null(event.freq.txt) && nzchar(event.freq.txt) && grepl(",", event.freq.txt)) {
        cat("                    Number of events: ", event.freq.txt, "\n", sep = "")
      } else {
        cat("                    Number of deaths: ", x$ndead, "\n", sep = "")
      }
    }
    if (x$family == "class" && !is.null(event.freq.txt))
      cat("           Frequency of class labels: ", event.freq.txt, "\n", sep = "")
    if (!is.null(x$imputed.indv))
      cat("                    Was data imputed: ", "yes", "\n", sep = "")
    cat("                     Number of trees: ", x$ntree, "\n", sep = "")
    cat("           Forest terminal node size: ", x$nodesize, "\n", sep = "")
    cat("       Average no. of terminal nodes: ", digits.pretty(mean(x$leaf.count), 4), "\n", sep = "")
    cat("No. of variables tried at each split: ", x$mtry, "\n", sep = "")
    cat("              Total no. of variables: ", length(x$xvar.names), "\n", sep = "")
    if (!is.null(outcome.target) && !x$univariate) {
      cat("              Total no. of responses: ", yvar.dim, "\n", sep = "")
      cat("         User has requested response: ", outcome.target, "\n", sep = "")
    }
    cat("       Resampling used to grow trees: ", samp.used, "\n", sep = "")
    cat("    Resample size used to grow trees: ", sampsize.pretty(x$forest$sampsize, x$n), "\n", sep = "")
    cat("                            Analysis: ", family.pretty, "\n", sep = "")
    cat("                              Family: ", family.org, "\n", sep = "")
    if (x$nsplit > 0 && x$splitrule != "random") {
      cat("                      Splitting rule: ", paste(x$splitrule, "*random*"), "\n", sep = "")
      cat("       Number of random split points: ", x$nsplit, "\n", sep = "")
    } else {
      cat("                      Splitting rule: ", x$splitrule, "\n", sep = "")
    }
    if (!is.null(err.rate)) {
      if (x$family == "regr" && !is.null(r.sq))
        cat("                     (OOB) R squared: ", digits.pretty(r.sq, 8), "\n", sep = "")
      if (x$family == "class" && !is.null(iratio))
        cat("                    Imbalanced ratio: ", digits.pretty(iratio, 4), "\n", sep = "")
      if (x$family == "class" && !is.null(brier.err))
        cat("                   (OOB) Brier score: ", digits.pretty(brier.err, 8), "\n", sep = "")
      if (x$family == "class" && !is.null(brier.norm.err))
        cat("        (OOB) Normalized Brier score: ", digits.pretty(brier.norm.err, 8), "\n", sep = "")
      if (x$family == "class" && !is.null(auc.err))
        cat("                           (OOB) AUC: ", digits.pretty(auc.err, 8), "\n", sep = "")
      if (x$family == "class" && !is.null(logloss.err))
        cat("                      (OOB) Log-loss: ", digits.pretty(logloss.err, 8), "\n", sep = "")
      if (x$family == "class" && !is.null(pr.auc.err))
        cat("                        (OOB) PR-AUC: ", digits.pretty(pr.auc.err, 8), "\n", sep = "")
      if (x$family == "class" && !is.null(gmean.err))
        cat("                        (OOB) G-mean: ", digits.pretty(gmean.err, 8), "\n", sep = "")
      if (x$family == "surv" && !is.null(crps.err))
        cat("                          (OOB) CRPS: ", digits.pretty(crps.err, 8), "\n", sep = "")
      if (x$family == "surv" && !is.null(crps.std.err))
        cat("             (OOB) standardized CRPS: ", digits.pretty(crps.std.err, 8), "\n", sep = "")
      if (quantile.flag && !is.null(q.crps))
        cat("                 (OOB) Quantile-CRPS: ", digits.pretty(q.crps, 8), "\n", sep = "")
      if (quantile.flag && !is.null(q.crps.std))
        cat("           (OOB) Quantile-CRPS (std): ", digits.pretty(q.crps.std, 8), "\n", sep = "")
      if (quantile.flag && length(q.pinball))
        cat("           (OOB) Pinball loss by tau: ", digits.pretty(q.pinball, 8),
            " (tau = ", digits.pretty(taus), ")", "\n", sep = "")
      cat("   (OOB) Requested performance error: ", err.rate, "\n\n", sep = "")
    }
    if (x$family == "class" && !is.null(conf.mat)) {
      if (!is.null(x$predicted.oob) && any(is.na(x$predicted.oob))) {
        cat("Confusion matrix (cases with missing OOB predicted values have been removed):\n\n")
      } else {
        cat("Confusion matrix:\n\n")
      }
      print(conf.mat)
      if (!is.null(miss.err)) cat("\n      (OOB) Misclassification rate: ", miss.err, "\n", sep = "")
      ## Strawman baselines
      if (!is.null(brier.rand) && !is.null(brier.norm.rand) && !is.null(logloss.rand)) {
        cat("\nRandom-classifier baselines (uniform):\n",
            "   Brier: ", digits.pretty(brier.rand, 8),
            "   Normalized Brier: ", digits.pretty(brier.norm.rand, 8),
            "   Log-loss: ", digits.pretty(logloss.rand, 8), "\n", sep = "")
      }
    }
  }
  ## ---------------- PREDICT MODE ----------------
  else {
    cat("  Sample size of test (predict) data: ", x$n, "\n", sep = "")
    if (grepl("surv", x$family) && !is.null(event.freq.txt)) {
      if (grepl(",", event.freq.txt)) {
        cat("       Number of events in test data: ", event.freq.txt, "\n", sep = "")
      } else {
        cat("       Number of deaths in test data: ", event.freq.txt, "\n", sep = "")
      }
    }
    if (!is.null(x$imputed.data))
      cat("               Was test data imputed: ", "yes", "\n", sep = "")
    cat("                Number of grow trees: ", x$ntree, "\n", sep = "")
    cat("  Average no. of grow terminal nodes: ", digits.pretty(mean(x$leaf.count), 4), "\n", sep = "")
    cat("         Total no. of grow variables: ", length(x$xvar.names), "\n", sep = "")
    if (!x$univariate) {
      cat("         Total no. of grow responses: ", yvar.dim, "\n", sep = "")
      cat("         User has requested response: ", outcome.target, "\n", sep = "")
    }
    cat("       Resampling used to grow trees: ", samp.used, "\n", sep = "")
    if (!is.null(outcome.target) && !is.null(x$forest$n))
      cat("    Resample size used to grow trees: ", sampsize.pretty(x$forest$sampsize, x$forest$n), "\n", sep = "")
    cat("                            Analysis: ", family.pretty, "\n", sep = "")
    cat("                              Family: ", family.org, "\n", sep = "")
    if (!is.null(err.rate)) {
      if (x$family == "regr" && !is.null(r.sq))
        cat("                           R squared: ", digits.pretty(r.sq, 8), "\n", sep = "")
      if (x$family == "class" && !is.null(iratio))
        cat("                    Imbalanced ratio: ", digits.pretty(iratio, 4), "\n", sep = "")
      if (x$family == "class" && !is.null(brier.err))
        cat("                         Brier score: ", digits.pretty(brier.err, 8), "\n", sep = "")
      if (x$family == "class" && !is.null(brier.norm.err))
        cat("              Normalized Brier score: ", digits.pretty(brier.norm.err, 8), "\n", sep = "")
      if (x$family == "class" && !is.null(auc.err))
        cat("                                 AUC: ", digits.pretty(auc.err, 8), "\n", sep = "")
      if (x$family == "class" && !is.null(logloss.err))
        cat("                            Log-loss: ", digits.pretty(logloss.err, 8), "\n", sep = "")
      if (x$family == "class" && !is.null(pr.auc.err))
        cat("                              PR-AUC: ", digits.pretty(pr.auc.err, 8), "\n", sep = "")
      if (x$family == "class" && !is.null(gmean.err))
        cat("                              G-mean: ", digits.pretty(gmean.err, 8), "\n", sep = "")
      if (x$family == "surv" && !is.null(crps.err))
        cat("                                CRPS: ", digits.pretty(crps.err, 8), "\n", sep = "")
      if (x$family == "surv" && !is.null(crps.std.err))
        cat("                     standardized CRPS: ", digits.pretty(crps.std.err, 8), "\n", sep = "")
      if (quantile.flag && !is.null(q.crps))
        cat("                       Quantile-CRPS: ", digits.pretty(q.crps, 8), "\n", sep = "")
      if (quantile.flag && !is.null(q.crps.std))
        cat("                 Quantile-CRPS (std): ", digits.pretty(q.crps.std, 8), "\n", sep = "")
            if (quantile.flag && length(q.pinball))
        cat("                 Pinball loss by tau: ", digits.pretty(q.pinball, 8),
            " (tau = ", digits.pretty(taus), ")", "\n", sep = "")
      cat("         Requested performance error: ", err.rate, "\n\n", sep = "")
    }
    if (x$family == "class" && !is.null(conf.mat)) {
      if (!is.null(x$predicted.oob) && any(is.na(x$predicted.oob))) {
        cat("Confusion matrix (cases with missing OOB predicted values have been removed):\n\n")
      } else {
        cat("Confusion matrix:\n\n")
      }
      print(conf.mat)
      if (!is.null(miss.err)) cat("\n           Misclassification error: ", miss.err, "\n", sep = "")
      ## Strawman baselines
      if (!is.null(brier.rand) && !is.null(brier.norm.rand) && !is.null(logloss.rand)) {
        cat("\nRandom-classifier baselines (uniform):\n",
            "   Brier: ", digits.pretty(brier.rand, 8),
            "   Normalized Brier: ", digits.pretty(brier.norm.rand, 8),
            "   Log-loss: ", digits.pretty(logloss.rand, 8), "\n", sep = "")
      }
    }
  }
}
