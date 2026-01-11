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
#' @name get.auc
#' @export
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
## cindex - extended to CR + uno/fenwick
get.cindex <- function(time, censoring, predicted, weight, fast, do.trace = FALSE) {
  size <- length(time)
  if (size != length(censoring)) {
    stop("time, censoring, and predicted must have the same length")
  }
  ## determine competing risks (CR) vs right-censoring
  isCR <- any(censoring > 1, na.rm = TRUE)
  ## weight is optional
  if (missing(weight)) {
    weight <- NULL
  } else {
    if (length(weight) != size) {
      stop("weight must have the same length as time")
    }
  }
  ## fast switch is dynamic - now set by the native library
  if (missing(fast)) {
    fast <- -1
  }
  ## -----------------------------
  ## Right-censoring
  ## -----------------------------
  if (!isCR) {
    ## accept predicted as a vector or a 1-column matrix/data.frame
    if (is.matrix(predicted) || is.data.frame(predicted)) {
      if (NROW(predicted) != size) {
        stop("predicted must have the same length as time")
      }
      if (NCOL(predicted) != 1) {
        stop("predicted must be a vector for right-censoring")
      }
      predicted <- predicted[, 1]
    }
    if (length(predicted) != size) {
      stop("time, censoring, and predicted must have the same length")
    }
    miss <- is.na(time) | is.na(censoring) | is.na(predicted)
    if (!is.null(weight)) {
      miss <- miss | is.na(weight)
    }
    if (sum(miss) == size) {
      stop("no valid pairs found, too much missing data")
    }
    ## Flag missing members so we can exclude them in the pairs.
    denom <- as.double(!miss)
    nativeOutput <- .Call("rfsrcCIndex",
                          as.integer(do.trace),
                          as.integer(fast),
                          as.integer(size),
                          as.double(time),
                          as.double(censoring),
                          as.double(predicted),
                          as.double(denom),
                          if (is.null(weight)) NULL else as.double(weight))
    if (is.null(nativeOutput)) {
      stop("An error has occurred in rfsrcCIndex.  Please turn trace on for further analysis.")
    }
    return(nativeOutput$err)
  }
  ## -----------------------------
  ## Competing Risks (CR path)
  ## -----------------------------
  ## Convert predicted into an n x J matrix (one column per event type).
  ## Allow: matrix/data.frame, or list of J vectors.
  baseMiss <- is.na(time) | is.na(censoring)
  if (!is.null(weight)) {
    baseMiss <- baseMiss | is.na(weight)
  }
  if (all(baseMiss)) {
    stop("no valid pairs found, too much missing data")
  }
  J <- max(censoring[!baseMiss], na.rm = TRUE)
  if (is.matrix(predicted) || is.data.frame(predicted)) {
    predMat <- as.matrix(predicted)
  } else if (is.list(predicted) && is.null(dim(predicted))) {
    if (length(predicted) < J) {
      stop("for competing risks, predicted must have at least J event-specific prediction vectors")
    }
    predMat <- do.call(cbind, lapply(predicted[1:J], as.double))
  } else {
    stop("for competing risks, predicted must be a matrix/data.frame with one column per event type, or a list of event-specific prediction vectors")
  }
  if (NROW(predMat) != size) {
    stop("for competing risks, predicted must have nrow == length(time)")
  }
  if (NCOL(predMat) < J) {
    stop("for competing risks, predicted must have at least one column per event type")
  }
  predMat <- predMat[, seq_len(J), drop = FALSE]
  err <- rep(NA_real_, J)
  ## (A) CR + weights: IPCW Fenwick, event-by-event
  if (!is.null(weight)) {
    for (j in seq_len(J)) {
      miss.j <- baseMiss | is.na(predMat[, j])
      denom.j <- as.double(!miss.j)
      nativeOutput <- .Call("rfsrcCIndexFenwick",
                            as.integer(do.trace),
                            as.integer(j),
                            as.integer(size),
                            as.double(time),
                            as.double(censoring),
                            as.double(predMat[, j]),
                            as.double(denom.j),
                            as.double(weight))
      if (is.null(nativeOutput)) {
        stop("An error has occurred in rfsrcCIndexFenwick.  Please turn trace on for further analysis.")
      }
      err[j] <- nativeOutput$err
    }
    return(err)
  }
  ## (B) CR + no weights: legacy conditional approach
  ##     subset to {0, j} and call right-censored concordance
  for (j in seq_len(J)) {
    keep <- !baseMiss & !is.na(predMat[, j])
    idx  <- which(keep & (censoring == 0 | censoring == j))
    if (length(idx) < 2) {
      err[j] <- NA_real_
      next
    }
    denom.sub <- rep(1, length(idx))
    nativeOutput <- .Call("rfsrcCIndex",
                          as.integer(do.trace),
                          as.integer(fast),  
                          as.integer(length(idx)),
                          as.double(time[idx]),
                          as.double(censoring[idx]),      # values are 0 or j (j>0 => event)
                          as.double(predMat[idx, j]),
                          as.double(denom.sub),
                          NULL)
    if (is.null(nativeOutput)) {
      stop("An error has occurred in rfsrcCIndex.  Please turn trace on for further analysis.")
    }
    err[j] <- nativeOutput$err
  }
  return(err)
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
## cross-entropy, log-loss
get.logloss <- function(y, prob, robust = TRUE) {
  lgl <- lapply(levels(y), function(yn) {
    phat <- prob[, yn]
    pt <- y == yn
    if (sum(pt) > 0) {
      -log(phat[pt])
    }
    else {
      0
    }
  })
  if (robust) {
    lgl <- unlist(lgl[lengths(lgl) != 0])
    mean(lgl[!is.infinite(lgl)], na.rm = TRUE)
  }
  else {
    mean(unlist(lgl[lengths(lgl) != 0]), na.rm = TRUE)
  }
}                          
