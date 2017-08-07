partial.rfsrc <- function(
  object,
  outcome.target = NULL,
  partial.type = NULL,
  partial.xvar = NULL,
  partial.values = NULL,
  partial.xvar2 = NULL,
  partial.values2 = NULL,
  partial.time = NULL,
  oob = TRUE,
  seed = NULL,
  do.trace = FALSE,
  ...)
{
  user.option <- list(...)
  terminal.qualts <- is.hidden.terminal.qualts(user.option)
  terminal.quants <- is.hidden.terminal.quants(user.option)
  if (missing(object)) {
    stop("object is missing!")
  }
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(rfsrc, grow)' or '(rfsrc, forest)'")
  }
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.null(object$forest)) {
      stop("The forest is empty.  Re-run rfsrc (grow) call with forest=TRUE")
    }
    object <- object$forest
  }
    else {
    }
  family <- object$family
  splitrule <- object$splitrule
  xvar.names <- object$xvar.names
  yvar.names <- object$yvar.names
  if (length(which(xvar.names == partial.xvar)) != 1) {
    stop("x-variable specified incorrectly:  ", partial.xvar)
  }
  if (!is.null(partial.xvar2)) {   
      if (length(partial.xvar2) != length(partial.values2)) {
          stop("second order x-variable and value vectors not of same length:  ", length(partial.xvar2), "vs", length(partial.values2))
      }
      for (i in 1:length(partial.xvar2)) {
          if (length(which(xvar.names == partial.xvar2[i])) != 1) {
              stop("second order x-variable element", i, "specified incorrectly:  ", partial.xvar2[i])
          }
      }
  }
  object$yvar <- as.data.frame(object$yvar)
  colnames(object$yvar) <- yvar.names
  yfactor <- extract.factor(object$yvar)
  outcome.target.idx <- get.outcome.target(family, yvar.names, outcome.target)
  yvar.types <- get.yvar.type(family, yfactor$generic.types, yvar.names, object$coerce.factor)
  yvar.nlevels <- get.yvar.nlevels(family, yfactor$nlevels, yvar.names, object$yvar, object$coerce.factor)
  event.info <- get.event.info(object)
  cr.bits <- get.cr.bits(family)
  xfactor <- extract.factor(object$xvar)
  xvar.types <- get.xvar.type(xfactor$generic.types, xvar.names, object$coerce.factor)
  xvar.nlevels <- get.xvar.nlevels(xfactor$nlevels, xvar.names, object$xvar, object$coerce.factor)
  ntree <- object$ntree
  na.action = object$na.action
  xvar <- as.matrix(data.matrix(object$xvar))
  yvar <- as.matrix(data.matrix(object$yvar))
  r.dim <- ncol(cbind(yvar))
  rownames(xvar) <- colnames(xvar) <- NULL
  n.xvar <- ncol(xvar)
  n <- nrow(xvar)
  outcome = "train"
  oob.bits <- get.oob(oob)
  bootstrap.bits <- get.bootstrap(object$bootstrap)
  na.action.bits <- get.na.action(na.action)
  samptype.bits <- get.samptype(object$samptype)
  partial.bits <- get.partial(length(partial.values))
  terminal.qualts.bits <- get.terminal.qualts(terminal.qualts, object$terminal.qualts)
  terminal.quants.bits <- get.terminal.quants(terminal.quants, object$terminal.quants)
  seed <- get.seed(seed)
  do.trace <- get.trace(do.trace)
  nativeOutput <- tryCatch({.Call("rfsrcPredict",
                                  as.integer(do.trace),
                                  as.integer(seed),
                                  as.integer(oob.bits + bootstrap.bits + cr.bits), 
                                  as.integer(
                                    samptype.bits +
                                      na.action.bits +
                                        terminal.qualts.bits +
                                          terminal.quants.bits +
                                            partial.bits),
                                  as.integer(ntree),
                                  as.integer(n),
                                  as.integer(r.dim),
                                  as.character(yvar.types),
                                  as.integer(outcome.target.idx),
                                  as.integer(length(outcome.target.idx)),
                                  as.integer(yvar.nlevels),
                                  as.double(as.vector(yvar)),
                                  as.integer(ncol(xvar)),
                                  as.character(xvar.types),
                                  as.integer(xvar.nlevels),
                                  as.double(xvar),
                                  as.integer(object$sampsize),
                                  as.integer(object$samp),
                                  as.double(object$case.wt),
                                  as.integer(length(event.info$time.interest)),
                                  as.double(event.info$time.interest),
                                  as.integer((object$nativeArray)$treeID),
                                  as.integer((object$nativeArray)$nodeID),
                                  as.integer((object$nativeArray)$parmID),
                                  as.double((object$nativeArray)$contPT),
                                  as.integer((object$nativeArray)$mwcpSZ),
                                  as.integer(object$nativeFactorArray),
                                  as.integer(object$nativeArrayTNDS$tnRMBR),
                                  as.integer(object$nativeArrayTNDS$tnAMBR),
                                  as.integer(object$nativeArrayTNDS$tnRCNT),
                                  as.integer(object$nativeArrayTNDS$tnACNT),
                                  as.integer(object$totalNodeCount),
                                  as.integer(object$seed),
                                  as.integer(get.rf.cores()),
                                  as.integer(0),
                                  as.integer(0),
                                  as.integer(NULL),
                                  as.integer(0),
                                  as.integer(NULL),
                                  as.integer(get.type(family, partial.type)),
                                  as.integer(which(xvar.names == partial.xvar)),
                                  as.integer(length(partial.values)),
                                  as.double(partial.values),
                                  as.integer(length(partial.xvar2)),
                                  as.integer(match(partial.xvar2, xvar.names)),
                                  as.double(partial.values2),
                                  as.integer(0),
                                  as.integer(0),
                                  as.double(NULL),
                                  as.double(NULL),
                                  as.double((object$nativeArrayTNDS$tnSURV)),
                                  as.double((object$nativeArrayTNDS$tnMORT)),
                                  as.double((object$nativeArrayTNDS$tnNLSN)),
                                  as.double((object$nativeArrayTNDS$tnCSHZ)),
                                  as.double((object$nativeArrayTNDS$tnCIFN)),
                                  as.double((object$nativeArrayTNDS$tnREGR)),
                                  as.integer((object$nativeArrayTNDS$tnCLAS)))}, error = function(e) {
                                    print(e)
                                    NULL})
  if (is.null(nativeOutput)) {
    stop("An error has occurred in prediction.  Please turn trace on for further analysis.")
  }
  rfsrcOutput <- list(call = match.call(),
                      family = family,
                      partial.time = partial.time)
  if (grepl("surv", family)) {
      partial.time.idx <- sapply(partial.time, function(x) {max(which(event.info$time.interest <= x))})
      if (sum(is.na(partial.time.idx)) > 0) {
          stop("partial.time must be a subset of the time interest vector contained in the model")
      }
  }
  if (family == "surv") {
    if ((partial.type == "rel.freq") || (partial.type == "mort")) {
      mort.names <- list(NULL, NULL)
      survOutput <- (if (!is.null(nativeOutput$partialSurv))
                       array(nativeOutput$partialSurv,
                             c(n, length(partial.values)),
                             dimnames=mort.names) else NULL)
    }
    else if (partial.type == "chf") {
      nlsn.names <- list(NULL, NULL, NULL)
      survOutput <- (if (!is.null(nativeOutput$partialSurv))
                       array(nativeOutput$partialSurv,
                             c(n, length(event.info$time.interest), length(partial.values)),
                             dimnames=nlsn.names) else NULL)
      if (!is.null(survOutput)) {
        survOutput <- survOutput[, partial.time.idx, , drop=FALSE]
      }
    }
      else if (partial.type == "surv") {
        surv.names <- list(NULL, NULL, NULL)
        survOutput <- (if (!is.null(nativeOutput$partialSurv))
                         array(nativeOutput$partialSurv,
                               c(n, length(event.info$time.interest), length(partial.values)),
                               dimnames=surv.names) else NULL)
        if (!is.null(survOutput)) {
          survOutput <- survOutput[, partial.time.idx, , drop=FALSE]
        }
      }
        else {
          stop("Invalid choice for 'partial.type' option:  ", partial.type)
        }
    rfsrcOutput <- c(rfsrcOutput, survOutput = list(survOutput))
  }
  else if (family == "surv-CR") {
    if (partial.type == "years.lost") {
      yrls.names <- list(NULL, NULL, NULL)
      survOutput <- (if (!is.null(nativeOutput$partialSurv))
                       array(nativeOutput$partialSurv,
                             c(n, length(event.info$event.type), length(partial.values)),
                             dimnames=yrls.names) else NULL)
    }
      else if (partial.type == "cif") {
        cifn.names <- list(NULL, NULL, NULL, NULL)
        survOutput <- (if (!is.null(nativeOutput$partialSurv))
                         array(nativeOutput$partialSurv,
                               c(n, length(event.info$time.interest), length(event.info$event.type), length(partial.values)),
                               dimnames=cifn.names) else NULL)
        if (!is.null(survOutput)) {
          survOutput <- survOutput[, partial.time.idx, , , drop=FALSE]
        }
      }
        else if (partial.type == "chf") {
          chfn.names <- list(NULL, NULL, NULL, NULL)
          survOutput <- (if (!is.null(nativeOutput$partialSurv))
                           array(nativeOutput$partialSurv,
                                 c(n, length(event.info$time.interest), length(event.info$event.type), length(partial.values)),
                                 dimnames=chfn.names) else NULL)
          if (!is.null(survOutput)) {
            survOutput <- survOutput[, partial.time.idx, , , drop=FALSE]
          }
        }
          else {
            stop("Invalid choice for 'partial.type' option:  ", partial.type)
          }
    rfsrcOutput <- c(rfsrcOutput, survOutput = list(survOutput))
  }
    else {
      class.index <- which(yvar.types != "R")
      class.count <- length(class.index)
      regr.index <- which(yvar.types == "R")
      regr.count <- length(regr.index)
      if (class.count > 0) {
        classOutput <- vector("list", class.count)
        names(classOutput) <- yvar.names[class.index]
        levels.count <- array(0, class.count)
        levels.names <- vector("list", class.count)
        counter <- 0
        for (i in class.index) {
            counter <- counter + 1
            levels.count[counter] <- yvar.nlevels[i]
            if (yvar.types[i] == "C") {
              levels.names[[counter]] <- yfactor$levels[[which(yfactor$factor == yvar.names[i])]]
            }
              else {
                levels.names[[counter]] <- yfactor$order.levels[[which(yfactor$order == yvar.names[i])]]
              }
        }
        iter.start <- 0
        iter.end   <- 0
        offset <- vector("list", class.count)
        for (p in 1:length(partial.values)) {
          for (k in 1:length(outcome.target.idx)) {
            target.idx <- which (class.index == outcome.target.idx[k])
            if (length(target.idx) > 0) {
              iter.start <- iter.end
              iter.end <- iter.start + ((1 + levels.count[target.idx]) * n)
              offset[[target.idx]] <- c(offset[[target.idx]], (iter.start+1):iter.end)
            }
          }
        }
        for (i in 1:length(outcome.target.idx)) {
          target.idx <- which (class.index == outcome.target.idx[i])
          if (length(target.idx) > 0) {
            ens.names <- list(NULL, c("all", levels.names[[target.idx]]), NULL)
            ensemble <- (if (!is.null(nativeOutput$partialClas))
                             array(nativeOutput$partialClas[offset[[target.idx]]],
                                   c(n, 1 + levels.count[target.idx], length(partial.values)),
                                   dimnames=ens.names) else NULL)
            classOutput[[target.idx]] <- ensemble
            remove(ensemble)
          }
        }
        rfsrcOutput <- c(rfsrcOutput, classOutput = list(classOutput))        
      }
      if (regr.count > 0) {
        regrOutput <- vector("list", regr.count)
        names(regrOutput) <- yvar.names[regr.index]
        iter.start <- 0
        iter.end   <- 0
        offset <- vector("list", regr.count)
        for (p in 1:length(partial.values)) {
          for (k in 1:length(outcome.target.idx)) {
            target.idx <- which (regr.index == outcome.target.idx[k])
            if (length(target.idx) > 0) {
              iter.start <- iter.end
              iter.end <- iter.start + n
              offset[[target.idx]] <- c(offset[[target.idx]], (iter.start+1):iter.end)
            }
          }
        }
        for (i in 1:length(outcome.target.idx)) {
          target.idx <- which (regr.index == outcome.target.idx[i])
          if (length(target.idx) > 0) {
            ens.names <- list(NULL, NULL)
            ensemble <- (if (!is.null(nativeOutput$partialRegr))
                             array(nativeOutput$partialRegr[offset[[target.idx]]],
                                   c(n, length(partial.values)),
                                   dimnames=ens.names) else NULL)
            regrOutput[[target.idx]] <- ensemble
            remove(ensemble)
          }
        }
        rfsrcOutput <- c(rfsrcOutput, regrOutput = list(regrOutput))
      }
    }
  class(rfsrcOutput) <- c("rfsrc", "partial",   family)
  return (rfsrcOutput)
}
partial <- partial.rfsrc
get.partial <- function (partial.length) {
  if (!is.null(partial.length)) {
    if (partial.length > 0) {
      bits <- 2^14
    }
      else if (partial.length == 0) {
        bits <- 0
      }
        else {
          stop("Invalid choice for 'partial.length' option:  ", partial.length)
        }
  }
    else {
      stop("Invalid choice for 'partial.length' option:  ", partial.length)
    }
  return (bits)
}
get.type <- function (family, partial.type) {
  if (family == "surv") {
    if (partial.type == "rel.freq") {
      partial.type <- "mort"
    }
    type <- match(partial.type, c("mort", "chf", "surv"))
  }
    else if (family == "surv-CR") {
      type <- match(partial.type, c("years.lost", "cif", "chf"))
    }
      else {
        type <- 0
      }
  if (is.na(type)) {
    stop("Invalid choice for 'partial.type' option:  ", partial.type)
  }
  return (type)
}
get.oob <- function (oob) {
  if (!is.null(oob)) {
    if (oob == TRUE) {
      oob <- 2^1
    }
      else if (oob == FALSE) {
        oob <- 2^0
      }
        else {
          stop("Invalid choice for 'oob' option:  ", oob)
        }
  }
    else {
      stop("Invalid choice for 'oob' option:  ", oob)
    }
  return (oob)
}
