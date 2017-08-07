adrop3d.last <- function(x, d, keepColNames = FALSE) {
  if (!is.array(x)) {
    x
  }
    else {
      if (d > 1) {
        x[,,1:d, drop = FALSE]
      }
        else {
          if (dim(x)[1] == 1) {
            rbind(x[,,1, drop = TRUE])
          }
            else {
              if (dim(x)[2] == 1) {
                if (keepColNames) {
                  xnew <- cbind(x[,,1, drop = TRUE])
                  colnames(xnew) <- colnames(x)
                  xnew
                }
                  else {
                    cbind(x[,,1, drop = TRUE])
                  }
              }
                else {
                  x[,,1, drop = TRUE]
                }
            }
        }
    }
}
adrop2d.first <- function(x, d, keepColNames = FALSE) {
  if (!is.array(x)) {
    x
  }
    else {
      if (d > 1) {
        x[1:d,, drop = FALSE]
      }
        else {
          x[1, , drop = TRUE]
        }
    }
}
adrop2d.last <- function(x, d, keepColNames = FALSE) {
  if (!is.array(x)) {
    x
  }
    else {
      if (d > 1) {
        x[,1:d, drop = FALSE]
      }
        else {
          x[,1, drop = TRUE]
        }
    }
}
amatrix <- function(x, d, names) {
  x <- matrix(x, d, dimnames = names)
  if (ncol(x) > 1) {
    x
  }
    else {
      c(x)
    }
}
amatrix.remove.names <- function(x) {
  if (!is.null(dim(x)) && ncol(x) == 1) {
    unlist(c(x), use.names = FALSE)
  }
    else {
      x
    }
}
atmatrix <- function(x, d, names, keep.names = FALSE) {
  x <- t(matrix(x, ncol = d, dimnames = names))
  if (ncol(x) > 1) {
    x
  }
    else {
      if (keep.names == FALSE) {
        c(x)
      }
        else {
          x.names <- rownames(x)
          x <- c(x)
          names(x) <- x.names
          x
        }
    }
}
avector <- function(x, name = FALSE) {
  if (!is.null(dim(x)) && nrow(x) > 1 && ncol(x) == 1) {
    x.names <- rownames(x)
    x <- unlist(c(x))
    if (name) names(x) <- x.names else names(x) <- NULL
    x
  }
    else if (!is.null(dim(x)) && nrow(x) == 1 && ncol(x) > 1) {
      x.names <- colnames(x)
      x <- unlist(c(x))
      if (name) names(x) <- x.names else names(x) <- NULL
      x
    }
      else if (!is.null(dim(x)) && nrow(x) == 1 && ncol(x) == 1) {
        unlist(c(x))
      }
        else {
          x
        }
}
available <- function (package, lib.loc = NULL, quietly = TRUE)
{
  package <- as.character(substitute(package))
  installed <- package %in% installed.packages()
  if (installed) {
    require(package, quietly = TRUE, character.only = TRUE)
  }
    else {
      return(invisible(FALSE))
    }
}
bayes.rule <- function(prob, pi.hat = NULL) {
  class.labels <- colnames(prob)
  if (is.null(pi.hat)) {
    factor(class.labels[apply(prob, 1, function(x) {
      if (!all(is.na(x))) {
        resample(which(x == max(x, na.rm = TRUE)), 1)
      }
      else {
        NA
      }
    })], levels = class.labels)
  }
  else {
    minority <- which.min(pi.hat)
    majority <- setdiff(1:2, minority)      
    drc.rule <- rep(majority, nrow(prob))
    drc.rule[prob[, minority] >= min(pi.hat, na.rm = TRUE)] <- minority
    class.labels[drc.rule]
  }
}
brier <- function(ytest, prob) {
  cl <- colnames(prob)
  J <- length(cl)
  bs <- rep(NA, J)
  nullO <- sapply(1:J, function(j) {
    bs[j] <<- mean((1 * (ytest == cl[j]) - prob[, j]) ^ 2, na.rm = TRUE)
    NULL
  })
  norm.const <- (J / (J - 1))
  sum(bs * norm.const, na.rm = TRUE)
}
class.error <- function(y, yhat) {
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
cv.folds <- function (n, folds = 10) {
  split(resample(1:n), rep(1:folds, length = n))
}
data.matrix <- function(x) {
  as.data.frame(lapply(x, function(xi) {
    if (is.integer(xi) || is.numeric(xi)) {
      xi
    }
      else if (is.logical(xi) || is.factor(xi)) {
        as.integer(xi)
      }
        else {
          as.numeric(xi)
        }
  }))
}
family.pretty <- function(fmly) {
  switch(fmly,
         "surv"     = "RSF",
         "surv-CR"  = "RSF",
         "regr"     = "RF-R",
         "class"    = "RF-C",
         "unsupv"   = "RF-U",
         "regr+"    = "mRF-R",
         "class+"   = "mRF-C",
         "mix+"     = "mRF-RC"
         )
}
finalizeFormula <- function(formula.obj, data) {
  yvar.names <- formula.obj$yvar.names
  all.names  <- formula.obj$all.names
  index      <- length(yvar.names)
  fmly       <- formula.obj$family
  ytry       <- formula.obj$ytry
  if (length(all.names) <= index) {
    stop("formula is misspecified: total number of variables does not exceed total number of y-variables")
  }
  if (all.names[index + 1] == ".") {
    if(index == 0) {
      xvar.names <- names(data)
    }
      else {
        xvar.names <- names(data)[!is.element(names(data), all.names[1:index])]
      }
  }
    else {
      if(index == 0) {
        xvar.names <- all.names
      }
        else {
          xvar.names <- all.names[-c(1:index)]
        }
      not.specified <- !is.element(xvar.names, names(data))
      if (sum(not.specified) > 0) {
        stop("formula is misspecified, object ", xvar.names[not.specified], " not found")
      }
    }
  return (list(family=fmly, yvar.names=yvar.names, xvar.names=xvar.names, ytry=ytry))
}
finalizeData <- function(fnames, data, na.action, miss.flag = TRUE) {
  data <- data[ , is.element(names(data), fnames), drop = FALSE]
  factor.names <- unlist(lapply(data, is.factor))
  if (sum(factor.names) > 0) {
    data[, factor.names] <- data.matrix(data[, factor.names, drop = FALSE])
  }
  if (miss.flag == TRUE && na.action == "na.omit") {
    if (any(is.na(data))) {
      data <- na.omit(data)
    }
  }
  if (nrow(data) == 0) {
    stop("no records in the NA-processed data: consider using 'na.action=na.impute'")
  }
  logical.names <- unlist(lapply(data, is.logical))
  if (sum(logical.names) > 0) {
    data[, logical.names] <- 1 * data[, logical.names, drop = FALSE]
  }
  character.names <- unlist(lapply(data, is.character))
  if (sum(character.names) > 0) {
    stop("data types cannot be character: please convert all characters to factors")
  }
  return (data)
}
get.importance.xvar <- function(importance.xvar, importance, object) {
  if (!is.null(importance)) {
    if (missing(importance.xvar) || is.null(importance.xvar)) {
      importance.xvar <- object$xvar.names
    }
      else {
        importance.xvar <- unique(importance.xvar)
        importance.xvar <- intersect(importance.xvar, object$xvar.names)
      }
    if (length(importance.xvar) == 0) {
      stop("xvar names do not match object xvar matrix")
    }
  }
    else {
      importance.xvar <- NULL
    }
  return (importance.xvar)
}
get.nmiss <- function(xvar, yvar = NULL) {
  if (!is.null(yvar)) {
    sum(apply(yvar, 1, function(x){any(is.na(x))}) | apply(xvar, 1, function(x){any(is.na(x))}))
  }
    else {
      sum(apply(xvar, 1, function(x){any(is.na(x))}))
    }
}
get.outcome.target <- function(family, yvar.names, outcome.target) {
  if (family == "regr" | family == "regr+" | family == "class" | family == "class+" | family == "mix+") {
    if (is.null(outcome.target)) {
      outcome.target <- yvar.names
    }
    outcome.target <- unique(outcome.target)
    outcome.target <- intersect(outcome.target, yvar.names)
    if (length(outcome.target) == 0) {
      stop("yvar target names do not match object yvar names")
    }
    outcome.target <- match(outcome.target, yvar.names)
  }
    else {
      outcome.target <- 0
    }
}
get.grow.nodesize <- function(fmly, nodesize) {
  if (fmly == "surv"){
    if (is.null(nodesize)) {
      nodesize <- 3
    }
  }
    else if (fmly == "surv-CR"){
      if (is.null(nodesize)) {
        nodesize <- 6
      }
    }
      else if (fmly == "class" | fmly == "class+") {
        if (is.null(nodesize)) {
          nodesize <- 1
        }
      }
        else if (fmly == "regr" | fmly == "regr+") {
          if (is.null(nodesize)) {
            nodesize <- 5
          }
        }
          else if (fmly == "mix+") {
            if (is.null(nodesize)) {
              nodesize <- 3
            }
          }
            else if (fmly == "unsupv") {
              if (is.null(nodesize)) {
                nodesize <- 3
              }
            }
              else if (is.null(nodesize)) {
                stop("family is misspecified")
              }
  nodesize <- round(nodesize)
}
get.coerced.survival.fmly <- function(fmly, event.type, splitrule = NULL) {
  if (grepl("surv", fmly)) {
    coerced.fmly <- "surv"
    if (!is.null(splitrule)) {
        if ((length(event.type) > 1) &&
            (splitrule != "l2.impute") &&
            (splitrule != "logrankscore")) {
        coerced.fmly <- "surv-CR"
      }
    }
      else {
        if (length(event.type) > 1) {
          coerced.fmly <- "surv-CR"
        }
      }
  }
    else {
      stop("attempt to coerce a non-survival family")
    }
  coerced.fmly
}
get.event.info <- function(obj, subset = NULL) {
  if (grepl("surv", obj$family)) {
    if (!is.null(obj$yvar)) {
      if (is.null(subset)) {
        subset <- (1:nrow(cbind(obj$yvar)))
      }
      r.dim <- 2
      time <- obj$yvar[subset, 1]
      cens <- obj$yvar[subset, 2]
      if (!all(floor(cens) == abs(cens), na.rm = TRUE)) {
        stop("for survival families censoring variable must be coded as a non-negative integer")
      }
      event <- na.omit(cens)[na.omit(cens) > 0]
      event.type <- sort(unique(event))
    }
      else {
        r.dim <- 0
        event <- event.type <- cens <- cens <- time <- NULL
      }
    time.interest <- obj$time.interest
  }
    else {
      if ((obj$family == "regr+") | (obj$family == "class+")) {
        r.dim <- dim(obj$yvar)[2]
      }
        else {
          r.dim <- 1
        }
      event <- event.type <- cens <- time.interest <- cens <- time <- NULL
    }
  return(list(event = event, event.type = event.type, cens = cens,
              time.interest = time.interest, time = time, r.dim = r.dim))
}
get.grow.event.info <- function(yvar, fmly, need.deaths = TRUE, ntime) {
  if (grepl("surv", fmly)) {
    r.dim <- 2
    time <- yvar[, 1]
    cens <- yvar[, 2]
    if (!all(floor(cens) == abs(cens), na.rm = TRUE)) {
      stop("for survival families censoring variable must be coded as a non-negative integer (perhaps the formula is set incorrectly?)")
    }
    if (need.deaths && (all(na.omit(cens) == 0))) {
      stop("no deaths in data!")
    }
    if (!all(na.omit(time) >= 0)) {
      stop("time must be  positive")
    }
    event.type <- unique(na.omit(cens))
    if (sum(event.type >= 0) != length(event.type)) {
      stop("censoring variable must be coded as NA, 0, or greater than 0.")
    }
    event <- na.omit(cens)[na.omit(cens) > 0]
    event.type <- unique(event)
    nonMissingOutcome <- which(!is.na(cens) & !is.na(time))
    nonMissingDeathFlag <- (cens[nonMissingOutcome] != 0)
    time.interest <- sort(unique(time[nonMissingOutcome[nonMissingDeathFlag]]))
    if (!missing(ntime)) {
      if (length(ntime) == 1 && length(time.interest) > ntime) {
        time.interest <- time.interest[
                                       unique(round(seq.int(1, length(time.interest), length.out = ntime)))]
      }
      if (length(ntime) > 1) {
        time.interest <- unique(sapply(ntime, function(tt) {
          time.interest[max(1, sum(tt >= time.interest, na.rm = TRUE))]
        }))
      }
    }
  }
    else {
      if ((fmly == "regr+") | (fmly == "class+") | (fmly == "mix+")) {
        r.dim <- dim(yvar)[2]
      }
        else {
          if (fmly == "unsupv") {
            r.dim <- 0
          }
            else {
              r.dim <- 1
            }
        }
      event <- event.type <- cens <- time.interest <- cens <- time <- NULL
    }
  return(list(event = event, event.type = event.type, cens = cens,
              time.interest = time.interest,
              time = time, r.dim = r.dim))
}
get.grow.splitinfo <- function (formula.detail, splitrule, nsplit, event.type) {
  splitrule.names <- c("logrank",              
                       "logrankscore",         
                       "logrankCR",            
                       "logrankACR",           
                       "random",               
                       "mse",                  
                       "mse.unwt",             
                       "mse.hvwt",             
                       "gini",                 
                       "gini.unwt",            
                       "gini.hvwt",            
                       "unsupv",               
                       "mv.mse",               
                       "mv.gini",              
                       "mv.mix",               
                       "custom",               
                       "l2.impute")            
  fmly <- formula.detail$family
  nsplit <- round(nsplit)
  if (nsplit < 0) {
    stop("Invalid nsplit value specified.")
  }
  cust.idx <- NULL
  splitpass <- FALSE
  if (!is.null(splitrule)) {
    if(grepl("custom", splitrule)) {
      splitrule.idx <- which(splitrule.names == "custom")
      cust.idx <- as.integer(sub("custom", "", splitrule))
      if (is.na(cust.idx)) cust.idx <- 1
      splitpass <- TRUE
    }
      else if (splitrule == "random") {
        splitrule.idx <- which(splitrule.names == "random")
        nsplit <- 1
        splitpass <- TRUE
      }
  }
  if (!splitpass) {
    if (grepl("surv", fmly)) {
      if (is.null(splitrule)) {
        if (length(event.type) ==  1) {
          splitrule.idx <- which(splitrule.names == "logrank")
        }
          else {
            splitrule.idx <- which(splitrule.names == "logrankCR")
          }
        splitrule <- splitrule.names[splitrule.idx]
      }
        else {
          splitrule.idx <- which(splitrule.names == splitrule)
          if (length(splitrule.idx) != 1) {
            stop("Invalid split rule specified:  ", splitrule)
          }
          if ((length(event.type) ==  1) & (splitrule.idx == which(splitrule.names == "logrankCR"))) {
            stop("Cannot specify logrankCR splitting for right-censored data")
          }
          if ((length(event.type) >   1) & (splitrule.idx == which(splitrule.names == "logrank"))) {
            splitrule.idx <- which(splitrule.names == "logrankACR")
          }
        }
    }
    if (fmly == "class") {
      if (is.null(splitrule)) {
        splitrule.idx <- which(splitrule.names == "gini")
        splitrule <- splitrule.names[splitrule.idx]
      }
        else {
          if ((splitrule != "gini") &
              (splitrule != "gini.unwt") &
              (splitrule != "gini.hvwt")) {
            stop("Invalid split rule specified:  ", splitrule)
          }
          splitrule.idx <- which(splitrule.names == splitrule)
        }
    }
    if (fmly == "regr") {
      if (is.null(splitrule)) {
        splitrule.idx <- which(splitrule.names == "mse")
        splitrule <- splitrule.names[splitrule.idx]
      }
        else {
          if ((splitrule != "mse") &
              (splitrule != "mse.unwt") &
              (splitrule != "mse.hvwt")) {
            stop("Invalid split rule specified:  ", splitrule)
          }
          splitrule.idx <- which(splitrule.names == splitrule)
        }
    }
    if (fmly == "regr+") {
      if (is.null(splitrule)) {
        splitrule.idx <- which(splitrule.names == "mv.mse")
        splitrule <- splitrule.names[splitrule.idx]
      }
        else {
          if ((splitrule != "mv.mse")) {
            stop("Invalid split rule specified:  ", splitrule)
          }
          splitrule.idx <- which(splitrule.names == splitrule)
        }
    }
    if (fmly == "class+") {
      if (is.null(splitrule)) {
        splitrule.idx <- which(splitrule.names == "mv.gini")
        splitrule <- splitrule.names[splitrule.idx]
      }
        else {
          if ((splitrule != "mv.gini")) {
            stop("Invalid split rule specified:  ", splitrule)
          }
          splitrule.idx <- which(splitrule.names == splitrule)
        }
    }
    if (fmly == "mix+") {
      if (is.null(splitrule)) {
        splitrule.idx <- which(splitrule.names == "mv.mse")
        splitrule <- "mv.mix"
      }
        else {
          if ((splitrule != "mv.mix")) {
            stop("Invalid split rule specified:  ", splitrule)
          }
          splitrule.idx <- which(splitrule.names == splitrule)
        }
    }
    if (fmly == "unsupv") {
      if (is.null(splitrule)) {
        splitrule.idx <- which(splitrule.names == "unsupv")
        splitrule <- splitrule.names[splitrule.idx]
      }
        else {
          if ((splitrule != "unsupv")) {
            stop("Invalid split rule specified:  ", splitrule)
          }
          splitrule.idx <- which(splitrule.names == splitrule)
        }
    }
  }
  splitinfo <- list(name = splitrule, index = splitrule.idx, cust = cust.idx, nsplit = nsplit)
  return (splitinfo)
}
get.weight <- function(weight, n) {
  if (!is.null(weight)) {
    if (any(weight < 0)      ||
        all(weight == 0)     ||
        length(weight) != n  ||
        any(is.na(weight))) {
      stop("Invalid weight vector specified.")
    }
  }
    else {
      weight <- rep(1, n)
    }
  return (weight)
}
get.grow.mtry <- function (mtry = NULL, n.xvar, fmly) {
  if (!is.null(mtry)) {
    mtry <- round(mtry)
    if (mtry < 1 | mtry > n.xvar) mtry <- max(1, min(mtry, n.xvar))
  }
    else {
      if (grepl("regr", fmly)) {
        mtry <- max(ceiling(n.xvar/3), 1)
      }
        else {
          mtry <- max(ceiling(sqrt(n.xvar)), 1)
        }
    }
  return (mtry)
}
get.ytry <- function(f) {
}
get.xvar.type <- function(generic.types, xvar.names, coerce.factor = NULL) {
  xvar.type <- generic.types
  if (!is.null(coerce.factor$xvar.names)) {
    xvar.type[is.element(xvar.names, coerce.factor$xvar.names)] <- "C"
  }
  xvar.type
}
get.mv.error <- function(obj, std = FALSE) {
  c(sapply(obj$yvar.names, function(nn) {
    o.coerce <- coerce.multivariate(obj, nn)
    if (o.coerce$family == "class") {
      utils::tail(o.coerce$err.rate[, 1], 1)
    }
    else {
      if (std) {
        utils::tail(o.coerce$err.rate, 1) / var(o.coerce$yvar, na.rm = TRUE)
      }
      else {
        utils::tail(o.coerce$err.rate, 1)
      }
    }
  }))
}
get.mv.predicted <- function(obj, oob = FALSE) {
  nms <- NULL
  pred <- do.call(cbind, lapply(obj$yvar.names, function(nn) {
    o.coerce <- coerce.multivariate(obj, nn)
    if (o.coerce$family == "class") {
      nms <<- c(nms, paste(nn, ".", colnames(o.coerce$predicted), sep = ""))
    }
    else {
      nms <<- c(nms, paste(nn))
    }
    if (oob) {
      o.coerce$predicted.oob
    }
    else {
      o.coerce$predicted
    }
  }))
  colnames(pred) <- nms
  pred
}
get.mv.vimp <- function(obj, std = FALSE) {
  vimp <- do.call(cbind, lapply(obj$yvar.names, function(nn) {
    o.coerce <- coerce.multivariate(obj, nn)
    if (o.coerce$family == "class") {
      o.coerce$importance[, 1]
    }
    else {
      if (std) {
        o.coerce$importance / var(o.coerce$yvar, na.rm = TRUE)
      }
      else {
        o.coerce$importance
      }
    }
  }))
  if (!is.null(vimp)) {
    colnames(vimp) <- obj$yvar.names
    return(vimp)
  }
  else {
    NULL
  }
}
get.xvar.nlevels <- function(nlevels, xvar.names, xvar, coerce.factor = NULL) {
  xvar.nlevels <- nlevels
  if (!is.null(coerce.factor$xvar.names)) {
    pt <- is.element(xvar.names, coerce.factor$xvar.names)
    xvar.nlevels[pt] <- sapply(coerce.factor$xvar.names, function(nn) {max(xvar[, nn])})
  }
  xvar.nlevels
}
get.yvar.type <- function(fmly, generic.types, yvar.names, coerce.factor = NULL) {
  if (fmly == "unsupv") {
    yvar.type <- NULL
  }
    else {
      if (grepl("surv", fmly)) {
        yvar.type <- c("T", "S")
      }
        else {
          yvar.type <- generic.types
          if (!is.null(coerce.factor$yvar.names)) {
            yvar.type[is.element(yvar.names, coerce.factor$yvar.names)] <- "C"
          }
        }
    }
  yvar.type
}
get.yvar.nlevels <- function(fmly, nlevels, yvar.names, yvar, coerce.factor = NULL) {
  if (fmly == "unsupv") {
    yvar.nlevels <- NULL
  }
  else {
    yvar.nlevels <- nlevels
    if (!is.null(coerce.factor$yvar.names)) {
      pt <- is.element(yvar.names, coerce.factor$yvar.names)
      yvar.nlevels[pt] <- sapply(coerce.factor$yvar.names, function(nn) {max(yvar[, nn])})
    }
  }
    yvar.nlevels
}
gmean <- function(y, prob, perf.type = NULL, robust = FALSE) {
  frq <- table(y)  
  if (length(frq) > 2) {
    return(NULL)
  }
  if (!is.null(perf.type) && perf.type == "g.mean.drc") {
    threshold <- min(frq, na.rm = TRUE) / sum(frq, na.rm = TRUE)
  }
  else {
    threshold <- 0.5
  }
  minority <- which.min(frq)
  majority <- setdiff(1:2, minority)
  yhat <- factor(1 * (prob[, minority] >= threshold), levels = c(0,1))
  confusion.matrix <- table(y, yhat)
  if (nrow(confusion.matrix) > 1) {
    TN <- confusion.matrix[minority, 2]
    FP <- confusion.matrix[minority, 1]
    FN <- confusion.matrix[majority, 2]
    TP <- confusion.matrix[majority, 1]
    if (robust) {
      sensitivity <- (1 + TP) / (1 + TP + FN)
      specificity <- (1 + TN) / (1 + TN + FP)
    }
    else {
      sensitivity <- TP / (TP + FN)
      specificity <- TN / (TN + FP)
    }
    sqrt(sensitivity * specificity)
  }
  else {
    NA
  }
}
parseFormula <- function(f, data, ytry = NULL, coerce.factor = NULL) {
  if (!inherits(f, "formula")) {
    stop("'formula' is not a formula object.")
  }
  if (is.null(data)) {
    stop("'data' is missing.")
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  fmly <- all.names(f, max.names = 1e7)[2]
  all.names <- all.vars(f, max.names = 1e7)
  yvar.names <- all.vars(formula(paste(as.character(f)[2], "~ .")), max.names = 1e7)
  yvar.names <- yvar.names[-length(yvar.names)]
  coerce.factor.org <- coerce.factor
  coerce.factor <- vector("list", 2)
  names(coerce.factor) <- c("xvar.names", "yvar.names")
  if (!is.null(coerce.factor.org)) {
    coerce.factor$yvar.names <- intersect(yvar.names, coerce.factor.org)
    if (length(coerce.factor$yvar.names) == 0) {
      coerce.factor$yvar.names <- NULL
    }
    coerce.factor$xvar.names <- intersect(setdiff(colnames(data), yvar.names), coerce.factor.org)
  }
  if ((fmly == "Surv")) {
    if (sum(is.element(yvar.names, names(data))) != 2) {
      stop("Survival formula incorrectly specified.")
    }
    family <- "surv"
    ytry <- 2
  }
    else if ((fmly == "Multivar" || fmly == "cbind")  && length(yvar.names) > 1) {
      if (sum(is.element(yvar.names, names(data))) < length(yvar.names)) {
        stop("Multivariate formula incorrectly specified: y's listed in formula are not in data.")
      }
      Y <- data[, yvar.names, drop = FALSE]
      logical.names <- unlist(lapply(Y, is.logical))
      if (sum(logical.names) > 0) {
        Y[, logical.names] <- 1 * Y[, logical.names, drop = FALSE]
      }
      if ((sum(unlist(lapply(Y, is.factor))) + 
          length(coerce.factor$yvar.names)) == length(yvar.names)) {
        family <- "class+"
      }
      else if ((sum(unlist(lapply(Y, is.factor))) + 
          length(coerce.factor$yvar.names)) == 0) {
        family <- "regr+"
      }
      else if (((sum(unlist(lapply(Y, is.factor))) +
                 length(coerce.factor$yvar.names)) > 0) && 
               ((sum(unlist(lapply(Y, is.factor))) +
                 length(coerce.factor$yvar.names)) < length(yvar.names))) {
        family <- "mix+"
      }
        else {
          stop("y-outcomes must be either real or factors in multivariate forests.")
        }
      if (!is.null(ytry)) {
        if ((ytry < 1) || (ytry > length(yvar.names))) {
          stop("invalid value for ytry:  ", ytry)
        }
      }
        else {
          ytry <- length(yvar.names)
        }
    }
      else if (fmly == "Unsupervised") {
        if (length(yvar.names) != 0) {
          stop("Unsupervised forests require no y-responses")
        }
        family <- "unsupv"
        yvar.names <- NULL
        temp <- gsub(fmly, "", as.character(f)[2])
        temp <- gsub("\\(|\\)", "", temp)
        ytry <- as.integer(temp)
        if (is.na(ytry)) {
          ytry <- 1
        }
          else {
            if (ytry <= 0) {
              stop("Unsupervised forests require positive ytry value")
            }
          }
      }
        else {
          if (sum(is.element(yvar.names, names(data))) != 1) {
            stop("formula is incorrectly specified.")
          }
          Y <- data[, yvar.names]
          if (is.logical(Y)) {
            Y <- as.numeric(Y)
          }
          if (!(is.factor(Y) | is.numeric(Y))) {
            stop("the y-outcome must be either real or a factor.")
          }
          if (is.factor(Y) || length(coerce.factor$yvar.names) == 1) {
            family <- "class"
          }
            else {
              family <- "regr"
            }
          ytry <- 1
        }
  return (list(all.names=all.names, family=family, yvar.names=yvar.names, ytry=ytry,
               coerce.factor = coerce.factor))
}
is.all.na <- function(x) {all(is.na(x))}
parseMissingData <- function(formula.obj, data) {
  yvar.names <- formula.obj$yvar.names
  if (length(yvar.names) > 0) {
    resp <- data[, yvar.names, drop = FALSE]
    col.resp.na <- unlist(lapply(data[, yvar.names, drop = FALSE], is.all.na))
    if (any(col.resp.na)) {
      stop("All records are missing for the yvar(s)")
    }
  }
  colPt <- unlist(lapply(data, is.all.na))
  if (sum(colPt) > 0 && sum(colPt) >= (ncol(data) - length(yvar.names))) {
    stop("All x-variables have all missing data:  analysis not meaningful.")
  }
  data <- data[, !colPt, drop = FALSE]
  rowPt <- apply(data, 1, is.all.na)
  if (sum(rowPt) == nrow(data)) {
    stop("Rows of the data have all missing data:  analysis not meaningful.")
  }
  data <- data[!rowPt,, drop = FALSE]
  return(data)
}
make.sample <- function(ntree, samp.size, boot.size = NULL) {
  if (samp.size < 0) {
    stop("samp.size cannot be negative:", samp.size)
  }
  if (is.null(boot.size)) {
    boot.size <- samp.size
  }
  rbind(sapply(1:ntree, function(bb){
    inb <- rep(0, samp.size)
    smp <- sample(1:samp.size, size = boot.size, replace = TRUE)
    frq <- tapply(smp, smp, length)
    idx <- as.numeric(names(frq))
    inb[idx] <- frq
    inb
  }))
}
make.imbalanced.sample <- function(ntree, ratio = 0.5, y) {
  if (ratio < 0 | ratio > 1) {
    stop("undersampling ratio must be between 0 and 1:", ratio)
  }
  frq <- table(y)
  class.labels <- names(frq)
  minority <- which.min(frq)
  majority <- setdiff(1:2, minority)
  n.minority <- min(frq, na.rm = TRUE)
  n.majority <- max(frq, na.rm = TRUE)
  samp.size <- length(y)
  rbind(sapply(1:ntree, function(bb){
    inb <- rep(0, samp.size)
    smp.min <- sample(which(y == class.labels[minority]), size = n.minority, replace = TRUE)
    smp.maj <- sample(which(y == class.labels[majority]), size = ratio * n.majority, replace = TRUE)
    smp <- c(smp.min, smp.maj)
    inb.frq <- tapply(smp, smp, length)
    idx <- as.numeric(names(inb.frq))
    inb[idx] <- inb.frq
    inb
  }))
}
make.size <- function(y) {
  frq <- table(y)
  min(length(y), min(frq, na.rm = TRUE) * length(frq))
}
make.wt <- function(y) {
  frq <- table(y)
  class.labels <- names(frq)
  wt <- rep(1, length(y))
  nullO <- sapply(1:length(frq), function(j) {
    wt[y == class.labels[j]] <<- prod(frq[-j], na.rm = TRUE)
    NULL
  })
  wt / max(wt, na.rm = TRUE)
}
resample <- function(x, size, ...) {
  if (length(x) <= 1) {
    if (!missing(size) && size == 0) x[FALSE] else x
  }
    else {
      sample(x, size, ...)
    }
}
row.col.deleted <- function(dat, r.n, c.n)
{
  which.r <- setdiff(r.n, rownames(dat))
  if (length(which.r) > 0) {
    which.r <- match(which.r, r.n)
  }
    else {
      which.r <- NULL
    }
  which.c <- setdiff(c.n, colnames(dat))
  if (length(which.c) > 0) {
    which.c <- match(which.c, c.n)
  }
    else {
      which.c <- NULL
    }
  return(list(row = which.r, col = which.c))
}
