### AUC workhorse
get.auc.workhorse <- function(roc.data) {
  x <- roc.data[, 1][roc.data[, 2] == 1]
  y <- roc.data[, 1][roc.data[, 2] == 0]
  if (length(x) > 1 & length(y) > 1) {
    AUC  <- tryCatch({wilcox.test(x, y, exact=F)$stat/(length(x)*length(y))}, error=function(ex){NA})
  }
  else {
    AUC <- NA
  }
  AUC
}
### Multiclass AUC -- Hand & Till (2001) definition
get.auc <- function(y, prob) {
  if (is.factor(y)) {
    y.uniq <- levels(y)
  }
  else {
    y.uniq <- sort(unique(y))
  }
  nclass <- length(y.uniq)
  AUC <- NULL
  for (i in 1:(nclass - 1)) {
    for (j in (i + 1):nclass) {
      pt.ij <- (y == y.uniq[i] | y == y.uniq[j])
      if (sum(pt.ij) > 1) {
        y.ij <- y[pt.ij]
        pij <- prob[pt.ij, j]
        pji <- prob[pt.ij, i]
        Aij <-  get.auc.workhorse(cbind(pij, 1 * (y.ij == y.uniq[j])))
        Aji <-  get.auc.workhorse(cbind(pji, 1 * (y.ij == y.uniq[i])))
        AUC <- c(AUC, (Aij + Aji)/2)
      }
    } 
  }
  if (is.null(AUC)) {
    NA
  }
  else {
    mean(AUC, na.rm = TRUE)
  }
}                          
## bayes rule 
get.bayes.rule <- function(prob, class.relfrq = NULL) {
  class.labels <- colnames(prob)
  if (is.null(class.relfrq)) {
    factor(class.labels[apply(prob, 1, function(x) {
      if (!all(is.na(x))) {
        resample(which(x == max(x, na.rm = TRUE)), 1)
      }
      else {
        NA
      }
    })], levels = class.labels)
  }
  ## added to handle the rfq classifier
  else {
    minority <- which.min(class.relfrq)
    majority <- setdiff(1:2, minority)      
    rfq.rule <- rep(majority, nrow(prob))
    rfq.rule[prob[, minority] >= min(class.relfrq, na.rm = TRUE)] <- minority
    factor(class.labels[rfq.rule], levels = class.labels)
  }
}
## normalized brier (normalized to one for strawman coin toss)
get.brier.error <- function(y, prob, normalized = TRUE, vector = FALSE) {
  if (is.null(colnames(prob))) {
    colnames(prob) <- levels(y)
  }
  cl <- colnames(prob)
  J <- length(cl)
  bs <- rep(NA, J)
  nullO <- sapply(1:J, function(j) {
    bs[j] <<- mean((1 * (y == cl[j]) - prob[, j]) ^ 2, na.rm = TRUE)
    NULL
  })
  if (normalized) {##normalized to 100
    norm.const <- (J / (J - 1))
  }
  else {##standard, normalized to 0.25
    norm.const <- (1 / J)
  }
  ## return the summary value 
  if (!vector) {
    sum(bs * norm.const, na.rm = TRUE)
  }
  ## return the vector of brier values
  else {
    bs * norm.const
  }
}
## get confusion matrix
get.confusion <- function(y, class.or.prob) {
  ## response or probability?
  if (is.factor(class.or.prob)) {
    confusion <- table(y, class.or.prob)
  }
  else {
    if (is.null(colnames(class.or.prob))) {
      colnames(class.or.prob) <- levels(y)
    }
    confusion <- table(y, get.bayes.rule(class.or.prob))
  }
  class.error <- 1 - diag(confusion) / rowSums(confusion, na.rm = TRUE)
  cbind(confusion, class.error = round(class.error, 4))
}
## cindex
get.cindex <- function (time, censoring, predicted, do.trace = FALSE) {
  size <- length(time)
  if (size != length(time) |
      size != length(censoring) |
      size != length(predicted)) {
    stop("time, censoring, and predicted must have the same length")
  }
  miss <- is.na(time) | is.na(censoring) | is.na(predicted)
  nmiss <- sum(miss)
  if (nmiss == size) {
    stop("no valid pairs found, too much missing data")
  }
  ## Flag missing members so we can exclude them in the pairs.
  denom <- sapply(miss, function(x) if (x) 0 else 1)
  nativeOutput <- .Call("rfsrcCIndex",
                        as.integer(do.trace),
                        as.integer(size),
                        as.double(time),
                        as.double(censoring),
                        as.double(predicted),
                        as.double(denom))
  ## check for error return condition in the native code
  if (is.null(nativeOutput)) {
    stop("An error has occurred in rfsrcCIndex.  Please turn trace on for further analysis.")
  }
  return (nativeOutput$err)
}
## gmean rule
get.gmean.rule <- function(y, prob) {
  ## determine frequencies: exit if this is not a two-class problem
  frq <- table(y)  
  if (length(frq) > 2) {
    return(NULL)
  }
  ## call bayes rule using class relative frequencies 
  get.bayes.rule(prob, frq / sum(frq, na.rm = TRUE))
}
## gmean for imbalanced classification
get.gmean <- function(y, prob, rfq = FALSE, robust = FALSE) {
  ## determine frequencies: exit if this is not a two-class problem
  frq <- table(y)  
  if (length(frq) > 2) {
    return(NULL)
  }
  ## determine threshold
  if (rfq) {
    threshold <- min(frq, na.rm = TRUE) / sum(frq, na.rm = TRUE)
  }
  else {
    threshold <- 0.5
  }
  ## convert yhat to a class label: 0 = majority, 1 = minority
  minority <- which.min(frq)
  majority <- setdiff(1:2, minority)
  yhat <- factor(1 * (prob[, minority] >= threshold), levels = c(0,1))
  ## compute the confusion matrix and extract TN,TP etc.
  confusion.matrix <- table(y, yhat)
  if (nrow(confusion.matrix) > 1) {##check that dimension is correct
    ## calculate the various rates
    TN <- confusion.matrix[minority, 2]
    FP <- confusion.matrix[minority, 1]
    FN <- confusion.matrix[majority, 2]
    TP <- confusion.matrix[majority, 1]
    ## assemble the sensitivity/specificity values
    if (robust) {##modified with 1 added to diagonals
      sensitivity <- (1 + TP) / (1 + TP + FN)
      specificity <- (1 + TN) / (1 + TN + FP)
    }
    else {
      sensitivity <- TP / (TP + FN)
      specificity <- TN / (TN + FP)
    }
    ## return the g mean
    sqrt(sensitivity * specificity)
  }
  else {
    NA
  }
}
## misclassification error
get.misclass.error <- function(y, yhat) {
  cl <- sort(unique(y))
  err <- rep(NA, length(cl))
  for (k in 1:length(cl)) {
    cl.pt  <- (y == cl[k])
    if (sum(cl.pt) > 0) {
        err[k] <- mean(y[cl.pt] != yhat[cl.pt])
    }
  }
  err
}
