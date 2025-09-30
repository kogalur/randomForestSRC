###################################################################
##
## performance measures (with random-baseline as an attribute)
##
###################################################################
.safe_div <- function(num, den) {
  num <- as.numeric(num); den <- as.numeric(den)
  ifelse(den == 0, NA_real_, num / den)
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
  ## not two-class -> NULL
  if (!is.factor(yvar) || length(levels(yvar)) != 2) return(NULL)
  ## recode to {0,1}: 0=majority, 1=minority
  y.frq <- table(yvar)
  if (length(y.frq) != 2 || any(is.na(y.frq))) return(NULL)
  class.labels <- names(y.frq)
  minority <- which.min(y.frq)
  majority <- setdiff(1:2, minority)
  n.eff    <- sum(y.frq)
  pihat    <- as.numeric(y.frq[minority]) / n.eff
  iratio   <- max(y.frq, na.rm = TRUE) / min(y.frq, na.rm = TRUE)
  y <- rep(0, length(yvar))
  y[yvar == class.labels[minority]] <- 1
  y <- factor(y, levels = c(0,1))
  ## map probabilities --> (majority, minority)
  if (!is.null(ncol(prob)) && ncol(prob) == 2) {
    prob.matx <- prob[, c(majority, minority), drop = FALSE]
    ## PRESERVE 1-column MATRIX for minority
    prob <- prob[, minority, drop = FALSE]
    ## clamp to [0,1]
    prob[] <- pmin(pmax(prob[], 0), 1)
  } else {
    ## user supplied a vector => assume it is minority prob
    prob <- as.numeric(prob)
    prob <- pmin(pmax(prob, 0), 1)
    prob.matx <- cbind(1 - prob, prob)
  }
  colnames(prob.matx) <- levels(y)
  ## threshold
  if (is.null(threshold)) threshold <- as.numeric(pihat) else threshold <- as.numeric(threshold)
  threshold <- max(0, min(1, threshold))
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
    F1         <- ifelse((prec + sens) > 0, 2 * (prec * sens) / (prec + sens), NA_real_)
    F1mod      <- ifelse(all(is.finite(c(sens, spec, prec, npv))) &&
                         all(c(sens, spec, prec, npv) > 0),
                         4 / (1/sens + 1/spec + 1/prec + 1/npv), NA_real_)
    gmean      <- sqrt(sens * spec)
    F1gmean    <- (F1 + gmean) / 2
    F1modgmean <- (F1mod + gmean) / 2
  } else {
    sens <- spec <- prec <- npv <- misclass <-
      F1 <- F1mod <- gmean <- F1gmean <- F1modgmean <- NA_real_
  }
  ## probability-based metrics
  brier       <- get.brier.error(y, prob.matx, normalized = FALSE)
  brier.norm  <- get.brier.error(y, prob.matx)   # normalized per your function
  auc         <- get.auc(y, prob.matx)
  logloss     <- get.logloss(y, prob.matx)
  prO         <- get.pr.auc(y, prob)             # vector or 1-col matrix is OK
  pr.auc      <- prO[1]
  pr.auc.rand <- prO[2]
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
  # prevalence pi (prefer prec.rand, else pr.auc.rand)
  pi.hat <- if (!is.null(rand.metrics) && "prec.rand" %in% names(rand.metrics)) {
    as.numeric(rand.metrics["prec.rand"])
  } else {
    .get("pr.auc.rand")
  }
  # header ----------------------------------------------------------
  cat("Imbalanced classification performance (model vs random baseline)\n")
  cat(sprintf(" n.majority = %s   n.minority = %s   iratio = %s   threshold = %s   prevalence  = %s\n\n",
              .get("n.majority"),
              .get("n.minority"),
              .fmt(.get("iratio")),
              .fmt(.get("threshold")),
              .fmt(pi.hat)))
  # metrics to display (grouped for readability)
  rate_metrics   <- c("sens","spec","prec","npv","F1","F1mod","gmean","F1gmean","F1modgmean")
  auc_metrics    <- c("auc","pr.auc")
  error_metrics  <- c("misclass","brier","brier.norm","logloss")
  all_metrics    <- c(rate_metrics, auc_metrics, error_metrics)
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
  # gains: higher-is-better vs lower-is-better
  hib <- c(rate_metrics, auc_metrics)                   # higher is better
  lib <- error_metrics                                  # lower is better
  delta <- rep(NA_real_, length(all_metrics)); names(delta) <- all_metrics
  gain  <- delta
  for (m in all_metrics) {
    mv <- model_vals[[m]]; bv <- base_vals[[m]]
    if (is.na(mv) || is.na(bv)) next
    if (m %in% lib) {
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
  # visual grouping by inserting blank lines (printed as separate blocks)
  print(df[df$Metric %in% rate_metrics,   , drop = FALSE], row.names = FALSE, right = TRUE)
  cat("\n")
  print(df[df$Metric %in% auc_metrics,    , drop = FALSE], row.names = FALSE, right = TRUE)
  cat("\n")
  print(df[df$Metric %in% error_metrics,  , drop = FALSE], row.names = FALSE, right = TRUE)
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
  x <- data.frame(do.call(rbind,
                          mclapply(seq(0,1,length=ngrid),function(th){
                            get.imbalanced.performance(yvar, prob, threshold = th)
                          })))
  if (measure == "gmean") {
    best <- which.max(x$gmean)
  }
  if (measure == "F1") {
    best <- which.max(x$F1)
  }
  if (measure == "F1mod") {
    best <- which.max(x$F1mod)
  }
  if (measure == "F1modgmean") {
    best <- which.max(x$F1modgmean)
  }
  if (plot.it) {
    par(mfrow = c(2,2))
    pt <- x$threshold < 3 * x$threshold[best]
    plot(x$threshold[pt], x$gmean[pt], xlab = "threshold", ylab = "gmean")
    abline(v = x$threshold[best], lty = 2, col = "blue")
    plot(x$threshold[pt], x$F1[pt], xlab = "threshold", ylab = "F1")
    abline(v = x$threshold[best], lty = 2, col = "blue")
    plot(x$threshold[pt], x$F1mod[pt], xlab = "threshold", ylab = "F1mod")
    abline(v = x$threshold[best], lty = 2, col = "blue")
    plot(x$threshold[pt], x$F1modgmean[pt], xlab = "threshold", ylab = "F1modgmean")
    abline(v = x$threshold[best], lty = 2, col = "blue")
  }
  x[best,, drop = FALSE]
}
###################################################################
##
##
## precision recall utilities (these are internal functions)
## assumes y=0--> majority, y=1--> minority
##
##
###################################################################
get.pr.auc <- function(truth, yhat) {
  ## failsafe (in case directly called)
  if (is.factor(truth) && !(sum(levels(truth) == c(0, 1)) == 2)) {
    y.frq <- table(truth)
    class.labels <- names(y.frq)
    minority <- which.min(y.frq)
    majority <- setdiff(1:2, minority)
    y <- rep(NA, length(truth))
    y[truth==class.labels[minority]] <- 1
    y[truth==class.labels[majority]] <- 0
    truth <- y
  }
  ## if yhat is a matrix, extract probabilities for minority class
  if (!is.null(ncol(yhat)) && ncol(yhat) == 2) {  
    yhat <- yhat[, which.min(table(truth))]  
  }
  ## create pr (x,y) data
  x <- yhat[truth == 1]
  y <- yhat[truth == 0]
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
  ## failsafe (in case directly called)
  if (is.factor(truth) && !(sum(levels(truth) == c(0, 1)) == 2)) {
    y.frq <- table(truth)
    class.labels <- names(y.frq)
    minority <- which.min(y.frq)
    majority <- setdiff(1:2, minority)
    y <- rep(NA, length(truth))
    y[truth==class.labels[minority]] <- 1
    y[truth==class.labels[majority]] <- 0
    truth <- y
  }
  ## if yhat is a matrix, extract probabilities for minority class
  if (!is.null(ncol(yhat)) && ncol(yhat) == 2) {  
    yhat <- yhat[, which.min(table(truth))]  
  }
  ## create pr (x,y) data
  x <- yhat[truth == 1]
  y <- yhat[truth == 0]
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
    m <- m[order(m[, 1], -m[, 3], decreasing = T), ]
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
