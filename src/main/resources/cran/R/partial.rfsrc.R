partial.rfsrc <- function(
  object,
  oob = TRUE,
  m.target = NULL,
  partial.type = NULL,
  partial.xvar = NULL,
  partial.values = NULL,
  partial.xvar2 = NULL,
  partial.values2 = NULL,
  partial.time = NULL,
  get.tree = NULL,
  seed = NULL,
  do.trace = FALSE,
  ...)
{
  ## hidden options
  user.option <- list(...)
  terminal.qualts <- is.hidden.terminal.qualts(user.option)
  terminal.quants <- is.hidden.terminal.quants(user.option)
  ## Object consistency
  ## (TBD, TBD, TBD) Version checking is NOT implemented in this mode (TBD, TBD, TBD)
  if (missing(object)) {
    stop("object is missing!")
  }
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(rfsrc, grow)' or '(rfsrc, forest)'")
  }
  ## acquire the forest
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.forest.missing(object)) {
      stop("The forest is empty.  Re-run rfsrc (grow call) with forest=TRUE")
    }
    object <- object$forest
  }
  else {
    ## object is already a forest.
  }
  ## convert oob to ensemble
  if (oob) {
    ensemble <- "oob"
  }
  else {
    ensemble <- "inbag"
  }
  ## multivariate family details
  family <- object$family
  splitrule <- object$splitrule
  ## pull the x-variable and y-outcome names from the grow object
  xvar.names <- object$xvar.names
  yvar.names <- object$yvar.names
  ## verify the x-var
  if (length(which(xvar.names == partial.xvar)) != 1) {
    stop("x-variable specified incorrectly:  ", partial.xvar)
  }
  ## verify the x-var2
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
  ## Caution:  There are no checks on partial.type.  Note that "rel.freq" and "mort"
  ## are equivalent from a native code perspective.
  ## Determine the immutable yvar factor map which is needed for
  ## classification sexp dimensioning.  But, first convert object$yvar
  ## to a data frame which is required for factor processing.
  object$yvar <- as.data.frame(object$yvar)
  colnames(object$yvar) <- yvar.names
  yfactor <- extract.factor(object$yvar)
  m.target.idx <- get.outcome.target(family, yvar.names, m.target)
  ## Get the y-outcome type and number of levels
  yvar.types <- get.yvar.type(family, yfactor$generic.types, yvar.names, object$coerce.factor)
  yvar.nlevels <- get.yvar.nlevels(family, yfactor$nlevels, yvar.names, object$yvar, object$coerce.factor)
  ## Get event information for survival families.
  event.info <- get.event.info(object)
  ## CR.bits assignment.
  cr.bits <- get.cr.bits(family)
  ## Determine the immutable xvar factor map.
  xfactor <- extract.factor(object$xvar)
  ## Get the x-variable type and number of levels.
  xvar.types <- get.xvar.type(xfactor$generic.types, xvar.names, object$coerce.factor)
  xvar.nlevels <- get.xvar.nlevels(xfactor$nlevels, xvar.names, object$xvar, object$coerce.factor)
  ## Initialize the number of trees in the forest.
  ntree <- object$ntree
  ## Use the training data na.action protocol.
  na.action = object$na.action
  ## Data conversion for training data.
  xvar <- as.matrix(data.matrix(object$xvar))
  yvar <- as.matrix(data.matrix(object$yvar))
  ## Set the y dimension.
  r.dim <- ncol(cbind(yvar))
  ## remove row and column names for proper processing by the native code
  ## set the dimensions.
  rownames(xvar) <- colnames(xvar) <- NULL
  n.xvar <- ncol(xvar)
  n <- nrow(xvar)
  sampsize <- round(object$sampsize(n))
  ## There is no test data.
  outcome = "train"
  ## Process the get.tree vector that specifies which trees we want
  get.tree <- get.tree.index(get.tree, ntree)
  ## Initialize the low bits.
  ensemble.bits <- get.ensemble(ensemble)
  bootstrap.bits <- get.bootstrap(object$bootstrap)
  na.action.bits <- get.na.action(na.action)
  ## Initalize the high bits
  samptype.bits <- get.samptype(object$samptype)
  partial.bits <- get.partial(length(partial.values))
  terminal.qualts.bits <- get.terminal.qualts(terminal.qualts, object$terminal.qualts)
  terminal.quants.bits <- get.terminal.quants(terminal.quants, object$terminal.quants)
  seed <- get.seed(seed)
  do.trace <- get.trace(do.trace)
  ## Check that hdim is initialized.  If not, set it zero.
  ## This is necessary for backwards compatibility with 2.3.0
  if (is.null(object$hdim)) {
    hdim <- 0
  }
  else {
    hdim <- object$hdim
  }
  ## Marker for start of native forest topology.  This can change with the outputs requested.
  ## For the arithmetic related to the pivot point, you need to refer to stackOutput.c and in
  ## particular, stackForestOutputObjects().
  ## added from generic.predict.rfsrc.R per UBK instruction 11/04/2019
  if (hdim == 0) {
    offset = 0
    chunk = 0
  } else {
    ## Offset starts at parmID2. We adjust for the presence of interactions.
    offset = 7
    ## A chunk is parmID2, contPT2, contPTR2, mwcpSZ2.
    chunk = 4
    if (!is.null(object$base.learner)) {
      if (object$base.learner$trial.depth > 1) {
        ## Offset with interactions is adjusted.
        ## Adjusted for AUGM_X1, AUGM_X2.
        offset = 9
        ## A chunk is parmID2, contPT2, contPTR2, mwcpSZ2, augmXone2, augmXtwo2.
        chunk = 6
      }
    }
  }
  pivot <- which(names(object$nativeArray) == "treeID")
  nativeOutput <- tryCatch({.Call("rfsrcPredict",
                                  as.integer(do.trace),
                                  as.integer(seed),
                                  as.integer(
                                      ensemble.bits +
                                      bootstrap.bits +
                                      cr.bits), 
                                  as.integer(
                                      samptype.bits +
                                      na.action.bits +
                                        terminal.qualts.bits +
                                          terminal.quants.bits +
                                            partial.bits),
                                  ## >>>> start of maxi forest object >>>>
                                  as.integer(ntree),
                                  as.integer(n),
                                  as.integer(r.dim),
                                  as.character(yvar.types),
                                  as.integer(yvar.nlevels),
                                  as.double(as.vector(yvar)),
                                  as.integer(ncol(xvar)),
                                  as.character(xvar.types),
                                  as.integer(xvar.nlevels),
                                  as.double(xvar),
                                  as.integer(sampsize),
                                  as.integer(object$samp),
                                  as.double(object$case.wt),
                                  list(if(is.null(event.info$time.interest)) as.integer(0) else as.integer(length(event.info$time.interest)),
                                       if(is.null(event.info$time.interest)) NULL else as.double(event.info$time.interest)),
                                  as.integer(object$totalNodeCount),
                                  as.integer(object$seed),
                                  as.integer(hdim),
                                  ## Object containing base learner settings, this is never NULL.
                                  object$base.learner, 
                                  as.integer((object$nativeArray)$treeID),
                                  as.integer((object$nativeArray)$nodeID),
                                  ## This is hc_zero.  It is never NULL.
                                  list(as.integer((object$nativeArray)$parmID),
                                  as.double((object$nativeArray)$contPT),
                                  as.integer((object$nativeArray)$mwcpSZ),
                                  as.integer((object$nativeFactorArray)$mwcpPT)),
                                  ## This slot is hc_zero_aug.  This slot can be NULL.
                                  if (!is.null(object$base.learner)) {
                                      if (object$base.learner$trial.depth > 1) {
                                          list(as.integer((object$nativeArray)$augmXone),
                                               as.integer((object$nativeArray)$augmXtwo))
                                      } else { NULL }
                                  } else { NULL },
                                  ## This slot is hc_one.  This slot can be NULL.                                  
                                  if (hdim > 0) {
                                      list(as.integer((object$nativeArray)$hcDim),
                                      as.double((object$nativeArray)$contPTR))
                                  } else { NULL },
                                  ## See the offset documentation in
                                  ## rfsrc.R after the nativeArray
                                  ## SEXP objects are populated in the
                                  ## post-forest parsing for an
                                  ## explanation of the offsets below.
                                  if (hdim > 1) {
                                      ## parmIDx
                                      lapply(0:(hdim-2), function(x) {as.integer(object$nativeArray[, offset + 1 + (chunk * x)])})
                                  } else { NULL },
                                  if (hdim > 1) {
                                      ## contPTx
                                      lapply(0:(hdim-2), function(x) {as.double(object$nativeArray[ , offset + 2 + (chunk * x)])})
                                  } else { NULL },
                                  if (hdim > 1) {
                                      ## contPTRx
                                      lapply(0:(hdim-2), function(x) {as.double(object$nativeArray[ , offset + 3 + (chunk * x)])})
                                  } else { NULL },
                                  if (hdim > 1) {
                                      ## mwcpSZx
                                      lapply(0:(hdim-2), function(x) {as.integer(object$nativeArray[, offset + 4 + (chunk * x)])})
                                  } else { NULL },
                                  if (hdim > 1) {
                                      ## mwcpPTx
                                      lapply(0:(hdim-2), function(x) {as.integer(object$nativeFactorArray[[x + 1]])})
                                  } else { NULL },
                                  if (hdim > 1) {
                                      if (!is.null(object$base.learner)) {
                                          if (object$base.learner$trial.depth > 1) {
                                              ## augmXonex
                                              lapply(0:(hdim-2), function(x) {as.integer(object$nativeArray[ , offset + 5 + (chunk * x)])})
                                          } else { NULL }
                                      } else { NULL }
                                  } else { NULL },
                                  if (hdim > 1) {
                                      if (!is.null(object$base.learner)) {
                                          if (object$base.learner$trial.depth > 1) {
                                              ## augmXtwox
                                              lapply(0:(hdim-2), function(x) {as.integer(object$nativeArray[ , offset + 6 + (chunk * x)])})
                                          } else { NULL }
                                      } else { NULL }
                                  } else { NULL },
                                  as.integer(object$nativeArrayTNDS$tnRMBR),
                                  as.integer(object$nativeArrayTNDS$tnAMBR),
                                  as.integer(object$nativeArrayTNDS$tnRCNT),
                                  as.integer(object$nativeArrayTNDS$tnACNT),
                                  as.double((object$nativeArrayTNDS$tnSURV)),
                                  as.double((object$nativeArrayTNDS$tnMORT)),
                                  as.double((object$nativeArrayTNDS$tnNLSN)),
                                  as.double((object$nativeArrayTNDS$tnCSHZ)),
                                  as.double((object$nativeArrayTNDS$tnCIFN)),
                                  as.double((object$nativeArrayTNDS$tnREGR)),
                                  as.integer((object$nativeArrayTNDS$tnCLAS)),
                                  ## <<<< end of maxi forest object <<<<
                                  list(if (is.null(m.target.idx)) as.integer(0) else as.integer(length(m.target.idx)),
                                       if (is.null(m.target.idx)) NULL else as.integer(m.target.idx)),
                                  as.integer(0),  ## Pruning disabled
                                    
                                  list(as.integer(0), NULL), ## Importance disabled.
                                  ## Partial variables enabled.  Note the as.integer is needed.
                                  list(as.integer(get.type(family, partial.type)),
                                       as.integer(which(xvar.names == partial.xvar)),
                                       as.integer(length(partial.values)),
                                       as.double(partial.values),
                                       as.integer(length(partial.xvar2)),
                                       if (length(partial.xvar2) == 0) NULL else as.integer(match(partial.xvar2, xvar.names)),
                                       as.double(partial.values2)),
                                  as.integer(0),     ## Subsetting disabled.
                                  as.integer(NULL),  ## Subsetting disabled.
                                  as.integer(0),    ## New data disabled.
                                  as.integer(0),    ## New data disabled.
                                  as.double(NULL),  ## New data disabled.
                                  as.double(NULL),  ## New data disabled.
                                  as.integer(ntree), ## block.size is hard-coded.
                                  list(as.integer(0), NULL, as.double(0)), ## Quantiles disabled.
                                  as.integer(get.tree),
                                  as.integer(get.rf.cores()))}, error = function(e) {
                                    print(e)
                                    NULL})
  ## check for error return condition in the native code
  if (is.null(nativeOutput)) {
    stop("An error has occurred in prediction.  Please turn trace on for further analysis.")
  }
  rfsrcOutput <- list(call = match.call(),
                      family = family,
                      partial.time = partial.time)
  ## Subset the user time vector from the grow time interest vector.
  if (grepl("surv", family)) {
      ## Get the indices of the closest points of the grow time interest vector.
      ## Exact:
      ## partial.time.idx <- match(partial.time, event.info$time.interest)
      ## Closest:
      partial.time.idx <- sapply(partial.time, function(x) {max(which(event.info$time.interest <= x))})
      if (sum(is.na(partial.time.idx)) > 0) {
          stop("partial.time must be a subset of the time interest vector contained in the model")
      }
  }
  if (family == "surv") {
    if ((partial.type == "rel.freq") || (partial.type == "mort")) {
      mort.names <- list(NULL, NULL)
      ## Incoming from the native code:
      ##   type = mort
      ##   -> of dim [length(partial.values)] x [1] x [1] x [n]
      ## Outgoing to the R code:
      ##   -> of dim [n] x [length(partial.values)]  
      survOutput <- (if (!is.null(nativeOutput$partialSurv))
                       array(nativeOutput$partialSurv,
                             c(n, length(partial.values)),
                             dimnames=mort.names) else NULL)
    }
    else if (partial.type == "chf") {
      nlsn.names <- list(NULL, NULL, NULL)
      ## Incoming from the native code:
      ##   type = chf
      ##   -> of dim [length(partial.values)] x [1] x [length(partial.time)] x [n]
      ## Outgoing to the R code:
      ##   -> of dim [n] x [length(partial.time)] x [length(partial.values)]  
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
        ## Incoming from the native code:
        ##   type = surv
        ##   -> of dim [length(partial.values)] x [1] x [length(partial.time)] x [n]
        ## Outgoing to the R code:
        ##   -> of dim [n] x [length(partial.time)] x [length(partial.values)]  
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
      ## Incoming from the native code:
      ##   type = years.lost
      ##   -> of dim [length(partial.values)] x [length(event.info$event.type)] x [1] x [n]
      ## Outgoing to the R code:
      ##   -> of dim [n] x [length(event.info$event.type)] x [length(partial.values)]  
      survOutput <- (if (!is.null(nativeOutput$partialSurv))
                       array(nativeOutput$partialSurv,
                             c(n, length(event.info$event.type), length(partial.values)),
                             dimnames=yrls.names) else NULL)
    }
      else if (partial.type == "cif") {
        cifn.names <- list(NULL, NULL, NULL, NULL)
        ## Incoming from the native code:
        ##   type = cif
        ##   -> of dim [length(partial.values)] x [length(event.info$event.type)] x [length(partial.time)] x [n]
        ## Outgoing to the R code:
        ##   -> of dim [n] x [length(partial.time)] x [length(event.info$event.type)] x [length(partial.values)]
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
          ## Incoming from the native code:
          ##   type = chfn
          ##   -> of dim [length(partial.values)] x [length(event.info$event.type)] x [length(partial.time)] x [n]
          ## Outgoing to the R code:
          ##   -> of dim [n] x [length(partial.time)] x [length(event.info$event.type)] x [length(partial.values)]
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
      ## We consider "R", "I", and "C" outcomes.  The outcomes are grouped
      ## by type and sequential.  That is, the first "C" encountered in the
      ## response type vector is in position [[1]] in the classification output
      ## list, the second "C" encountered is in position [[2]] in the
      ## classification output list, and so on.  The same applies to the
      ## regression outputs.  We also have a mapping from the outcome slot back
      ## to the original response vector type, given by the following:
      ## Given yvar.types = c("R", "C", "R", "C", "R" , "I")
      ## regr.index[1] -> 1
      ## regr.index[2] -> 3
      ## regr.index[3] -> 5
      ## clas.index[1] -> 2
      ## clas.index[2] -> 4
      ## clas.index[3] -> 6
      ## This will pick up all "C" and "I".
      class.index <- which(yvar.types != "R")
      class.count <- length(class.index)
      regr.index <- which(yvar.types == "R")
      regr.count <- length(regr.index)
      if (class.count > 0) {
        ## Create and name the classification outputs.
        classOutput <- vector("list", class.count)
        names(classOutput) <- yvar.names[class.index]
        ## Vector to hold the number of levels in each factor response. 
        levels.count <- array(0, class.count)
        ## List to hold the names of levels in each factor response. 
        levels.names <- vector("list", class.count)
        counter <- 0
        for (i in class.index) {
            counter <- counter + 1
            ## Note that [i] is the actual index of the y-variables and not a sequential iterator.
            ## The sequential iteratior is [counter]
            levels.count[counter] <- yvar.nlevels[i]
            if (yvar.types[i] == "C") {
              ## This an unordered factor.
              ## Here, we don't know the sequence of the unordered factor list, so we identify the factor by name.
              levels.names[[counter]] <- yfactor$levels[[which(yfactor$factor == yvar.names[i])]]
            }
              else {
                ## This in an ordered factor.
                ## Here, we don't know the sequence of the ordered factor list, so we identify the factor by name.
                levels.names[[counter]] <- yfactor$order.levels[[which(yfactor$order == yvar.names[i])]]
              }
        }
        iter.start <- 0
        iter.end   <- 0
        offset <- vector("list", class.count)
        for (p in 1:length(partial.values)) {
          for (k in 1:length(m.target.idx)) {
            target.idx <- which (class.index == m.target.idx[k])
            if (length(target.idx) > 0) {
              iter.start <- iter.end
              iter.end <- iter.start + ((1 + levels.count[target.idx]) * n)
              offset[[target.idx]] <- c(offset[[target.idx]], (iter.start+1):iter.end)
            }
          }
        }
        for (i in 1:length(m.target.idx)) {
          target.idx <- which (class.index == m.target.idx[i])
          if (length(target.idx) > 0) {
            ens.names <- list(NULL, c("all", levels.names[[target.idx]]), NULL)
            ## Incoming from the native code:
            ##   type = NULL
            ##   -> of dim [length(partial.values)] x [length(1 + yvar.nlevels[.]] x [n]
            ## Outgoing to the R code:
            ##   -> of dim [n] x [1 + yvar.nlevels[.]] x [length(partial.values)]
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
        ## Create and name the classification outputs.
        regrOutput <- vector("list", regr.count)
        names(regrOutput) <- yvar.names[regr.index]
        iter.start <- 0
        iter.end   <- 0
        offset <- vector("list", regr.count)
        for (p in 1:length(partial.values)) {
          for (k in 1:length(m.target.idx)) {
            target.idx <- which (regr.index == m.target.idx[k])
            if (length(target.idx) > 0) {
              iter.start <- iter.end
              iter.end <- iter.start + n
              offset[[target.idx]] <- c(offset[[target.idx]], (iter.start+1):iter.end)
            }
          }
        }
        for (i in 1:length(m.target.idx)) {
          target.idx <- which (regr.index == m.target.idx[i])
          if (length(target.idx) > 0) {
            ens.names <- list(NULL, NULL)
            ## Incoming from the native code:
            ##   type = NULL
            ##   -> of dim [length(partial.values)] x [1] x [n]
            ## Outgoing to the R code:
            ##   -> of dim [n] x [length(partial.values)]
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
  ## Convert partial option into native code parameter.
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
    ## The native code interprets "rel.freq" as "mort".  The R-side
    ## handles the difference downstream after native code exit.
    if (partial.type == "rel.freq") {
      partial.type <- "mort"
    }
    ## Warning:  Hard coded in global.h
    type <- match(partial.type, c("mort", "chf", "surv"))
  }
    else if (family == "surv-CR") {
      ## Warning:  Hard coded in global.h
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
