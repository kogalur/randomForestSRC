###################################################################
##
##
## performance measures
##
##
###################################################################
get.imbalanced.performance <- function(obj,
                                       prob = NULL,
                                       threshold = NULL,
                                       confusion = FALSE,
                                       robust = FALSE)
{
  ## determine if a forest object is provided or two vectors (yvar, prob)
  if (is.null(prob)) {
    if (class(obj)[1] != "rfsrc") {
      stop("obj must be a forest object")
    }
    else {
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
  ## if this is not a two-class prob, return with NULL
  if (!is.factor(yvar) || length(levels(yvar)) != 2) {
    NULL
  }
  ## acquire performance values from workhorse
  perf.o <- get.imbalanced.performance.workhorse(yvar, prob, threshold = threshold,
                       confusion = confusion, robust = robust)
  ## default is to return as a named vector
  if (!confusion) {
    unlist(perf.o)
  }
  ## when confusion is requested, results are returned as a list
  else {
    perf.o
  }
}
## extract performance values
## gmean, sensitivity, specificity, recall, F1, etc.
get.imbalanced.performance.workhorse <- function (yvar, prob,
                                                  threshold = NULL,
                                                  confusion = FALSE,
                                                  robust = FALSE)
{
  ## if this is not a two-class forest object, return with NULL
  if (!is.factor(yvar) || length(levels(yvar)) != 2) {
    return(NULL)
  }
  ## transform y so that class labels are {0,1}
  ## y=0 --> majority
  ## y=1 --> minority
  y.frq <- table(yvar)
  class.labels <- names(y.frq)
  minority <- which.min(y.frq)
  majority <- setdiff(1:2, minority)
  pihat <- y.frq[minority] / length(yvar)
  iratio <- max(y.frq, na.rm = TRUE) / min(y.frq, na.rm = TRUE)
  y <- rep(0, length(yvar))
  y[yvar==class.labels[minority]] <- 1
  y <- factor(y, levels = c(0,1))
  ## map probabilities --> 0=majority, 1=minority
  ## save probability as a vector and matrix
  if (!is.null(ncol(prob)) && ncol(prob) == 2) {
    prob.matx <- prob[, c(majority, minority), drop = FALSE]
    prob <- prob[, minority, drop = FALSE]
  }
  else {## user supplied a vector, we always assume this is for minority class
    prob.matx <- cbind(1 - prob, prob)
  }
  colnames(prob.matx) <- levels(y)
  ## get rfq threshold (default) if one is not supplied    
  if (is.null(threshold)) {
    threshold <- as.numeric(pihat)
  }
  else {
    threshold <- as.numeric(threshold)
  }
  ##-------------------------------------------
  ## extract table rates
  ##
  ##             predicted
  ## T          0        1
  ## R     0    TN       FP
  ## U     1    FN       TP
  ## E
  ##
  ## sens=TP/(FN+TP)      =tpr=recall
  ## spec=TN/(TN+FP)      =tnr=1-fpr
  ## ppv=TP/(TP+FP)       =prec
  ## npv=TN/(TN+FN)  
  ## ------------------------------------------
  yhat <- factor(1 * (prob >= threshold), levels = c(0, 1))
  confusion.matx <- table(y, yhat)
  if (nrow(confusion.matx) > 1) {
    TN <- confusion.matx[1, 1]
    FP <- confusion.matx[1, 2]
    FN <- confusion.matx[2, 1]
    TP <- confusion.matx[2, 2]
    if (robust) {
      sens <- (1 + TP) / (1 + FN + TP)
      spec <- (1 + TN) / (1 + TN + FP)
      prec  <- (1 + TP) / (1 + TP + FP)
      npv <- (1 + TN) / (1 + TN + FN)
    }
    else {
      sens <- TP / (FN + TP)
      spec <- TN / (TN + FP)
      prec <- TP / (TP + FP)
      npv <- TN / (TN + FN)
    }
    ## performance values based on the confusion matrix
    misclass <- (FP + FN) / length(y)
    F1 <- 2 * (prec * sens) / (prec + sens)
    F1mod <- 4 / (1 / sens + 1 / spec + 1 / prec + 1 / npv)
    gmean <- sqrt(sens * spec)
    F1gmean <- (F1 + gmean) / 2
    F1modgmean <- (F1mod + gmean) / 2
  }
  else {
    sens <- spec <- prec <- npv <- misclass <- 
      F1 <- F1mod <- gmean <- F1gmean <- F1modgmean <- NA
  }
  ## performance values based on probability
  brier <- get.brier.error(y, prob.matx, normalized = FALSE)
  brier.norm <- get.brier.error(y, prob.matx)
  auc <- get.auc(y, prob.matx)
  logloss <- get.logloss(y, prob.matx)
  prO <- get.pr.auc(y, prob)
  pr.auc <- prO[1]
  pr.auc.rand <- prO[2]
  ## assemble output as a list
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
  ## add confusion?
  if (confusion) {
    class.error <- 1 - diag(confusion.matx) / rowSums(confusion.matx, na.rm = TRUE)
    rO$confusion <- cbind(confusion.matx, class.error = round(class.error, 4))
  }
  rO
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
