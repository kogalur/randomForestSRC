get.imbalanced.performance <- function(o, threshold = NULL, robust = FALSE) {
  ## if this is not a two-class forest object, return with NULL
  if (!is.factor(o$yvar) || length(levels(o$yvar)) != 2) {
    return(NULL)
  }
  ## get brier and auc using built in functions
  brier <- get.brier.error(o$yvar,
     if(!is.null(o$predicted.oob) && !all(is.na(o$predicted.oob))) o$predicted.oob else o$predicted)
  auc <- get.auc(o$yvar,
     if(!is.null(o$predicted.oob) && !all(is.na(o$predicted.oob))) o$predicted.oob else o$predicted)
  ## acquire all other metrics from built in performance evaluator
  perf.o <- get.imbalanced.performance.workhorse(o, threshold = threshold, robust = robust) 
  ## assemble output as a list, with gmean last
  rO <- c(perf.o, list(brier = brier, auc = auc))
  unlist(rO)[c(1:7, 13:14, 8:12)]
}
## extract performance values
## gmean, sensitivity, specificity, recall, F1, etc.
get.imbalanced.performance.workhorse <- function (o, threshold = NULL, robust = FALSE) {
  ## if this is not a two-class forest object, return with NULL
  if (!is.factor(o$yvar) || length(levels(o$yvar)) != 2) {
    return(NULL)
  }
  ## process y - this implicitly assumes we are dealing with 0/1 class data
  ## y=0 --> majority
  ## y=1 --> minority
  y.frq <- table(o$yvar)
  class.labels <- names(y.frq)
  minority <- which.min(y.frq)
  majority <- setdiff(1:2, minority)
  pihat <- y.frq[minority]/o$n
  iratio <- max(y.frq, na.rm = TRUE)/min(y.frq, na.rm = TRUE)
  y <- rep(NA, o$n)
  y[o$yvar==class.labels[minority]] <- 1
  y[o$yvar==class.labels[majority]] <- 0
  y <- factor(y)
  ## extract probabilities --> map them so that 0=majority, 1=minority
  if (!is.null(o$predicted.oob) && !all(is.na(o$predicted.oob))) {
    prob <- o$predicted.oob[, c(majority, minority)]
  }
  else {
    prob <- o$predicted[, c(majority, minority)]
  }
  ## get rfq threshold (default) if one is not supplied    
  if (is.null(threshold)) {
    threshold <- as.numeric(pihat)
  }
  else {
    threshold <- as.numeric(threshold)
  }
  ##-------------------------------------------
  ## extract performance values
  ##
  ##             predicted
  ## T          0        1
  ## R     0    TN       FP
  ## U     1    FN       TP
  ## E
  ##
  ## sens=TP/(FN+TP)      =tpr=recall
  ## spec=TN/(TN+FP)      =tnr=1-fpr
  ## prec=TP/(TP+FP)      =ppv
  ## ------------------------------------------
  yhat <- factor(1 * (prob[, 2] >= threshold), levels = c(0, 1))
  confusion.matrix <- table(y, yhat)
  if (nrow(confusion.matrix) > 1) {
    TN <- confusion.matrix[1, 1]
    FP <- confusion.matrix[1, 2]
    FN <- confusion.matrix[2, 1]
    TP <- confusion.matrix[2, 2]
    if (robust) {
      sens <- (1 + TP)/(1 + FN + TP)
      spec <- (1 + TN)/(1 + TN + FP)
      prec <- (1 + TP)/(1 + TP + FP)
    }
    else {
      sens <- TP/(FN + TP)
      spec <- TN/(TN + FP)
      prec <- TP/(TP + FP)
    }
    ## performance values
    gmean <- sqrt(sens * spec)
    F1 <- 2 * (prec * sens) / (prec + sens) 
    balanced<- (sens + spec) / 2
    prO <- get.pr.auc(y, prob[, 2])
    pr.auc <- prO[1]
    pr.auc.rand <- prO[2]
  }
  else {
    sens <- spec <- prec <- gmean <- F1 <- balanced <- pr.auc <- pr.auc.rand <- NA
  }
  list(n.majority = as.numeric(y.frq[majority]),
       n.minority = as.numeric(y.frq[minority]),
       iratio = iratio,
       threshold = threshold,
       sens = sens,
       spec = spec,
       prec = prec,
       F1 = F1,
       balanced = balanced,
       pr.auc.rand = pr.auc.rand,
       pr.auc = pr.auc,
       gmean = gmean)
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
  if (length(frq) > 2) {
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
get.imbalanced.optimize <- function(o,
                                    newdata = NULL,
                                    measure = c("gmean", "balanced", "F1"),
                                    ngrid = 1000,
                                    plot.it = TRUE) {
  measure <- match.arg(measure, c("gmean", "balanced", "F1"))
  ## if new data present, replace object with predicted object
  if (!is.null(newdata)) {
    o <- predict(o, newdata)
  }
  x=data.frame(do.call(rbind,
      mclapply(seq(0,1,length=ngrid),function(th){get.imbalanced.performance(o,threshold=th)})))
  if (measure == "gmean") {
    best <- which.max(x$gmean)
  }
  if (measure == "balanced") {
    best <- which.max(x$balanced)
  }
  if (measure == "F1") {
    best <- which.max(x$F1)
  }
  if (plot.it) {
    par(mfrow = c(2,2))
    pt <- x$threshold < 3 * x$threshold[best]
    plot(x$threshold[pt], x$gmean[pt], xlab = "threshold", ylab = "gmean")
    abline(v = x$threshold[best], lty = 2, col = "blue")
    plot(x$threshold[pt], x$balanced[pt], xlab = "threshold", ylab = "balanced")
    abline(v = x$threshold[best], lty = 2, col = "blue")
    plot(x$threshold[pt], x$F1[pt], xlab = "threshold", ylab = "F1")
    abline(v = x$threshold[best], lty = 2, col = "blue")
  }
  x[best,, drop = FALSE]
}
###################################################################
##
##
## precision recall auc function with random calibration
## assumes y=0--> majority, y=1--> minority
##
##
###################################################################
get.pr.auc <- function(truth, yhat) {
  x <- yhat[truth == 1]
  y <- yhat[truth == 0]
  rO <- c(NA, NA)
  if (length(x) > 0 && length(y) > 0) {
    ## precision recall with random baseline
    pr.o <- tryCatch({get.pr.curve(x, y, rand.compute = TRUE)}, error = function(ex) {NULL})
    if (!is.null(pr.o)) {
      rO <- c(pr.o$auc.integral, pr.o$auc.integral.rand)
    }
  }
  rO
}
## copied from PRROC library - we gratefully acknowledge the developer for this function
get.pr.curve <- function(scores.class0, scores.class1 = scores.class0, weights.class0 = NULL,
                     weights.class1 = {if (is.null(weights.class0)) {NULL} else {1 - weights.class0}},
                     sorted = FALSE,
                     minStepSize = min(1, ifelse(is.null(weights.class0), 1, sum(weights.class0)/100)),
                     rand.compute = TRUE) 
{
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
    get.pr.curve.workhorse(scores.class0, scores.class1, weights.class0, 
                              weights.class1, minStepSize, rand.compute)
}
get.pr.curve.workhorse <- function(sorted.scores.class0, sorted.scores.class1,
                  weights.class0, weights.class1, minStepSize, rand.compute) 
{
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
  rand.auc <- NULL
  if (rand.compute) {
    if (is.null(weights.class0)) {
      rand.auc <- length(sorted.scores.class0)/(length(sorted.scores.class0) + 
                                                length(sorted.scores.class1))
    }
    else {
      rand.auc <- sum(weights.class0)/sum(weights.class0 + weights.class1)
    }
  }
  c(res, list(auc.integral.rand = rand.auc))
}
