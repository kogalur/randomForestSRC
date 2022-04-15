partial.rfsrc <- function(
  object,
  oob = TRUE,
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
  ## Object consistency
  ## (TBD3) Version checking is NOT implemented in this mode
  if (missing(object)) {
    stop("object is missing!")
  }
  ## coherence checks
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(rfsrc, grow)' or '(rfsrc, forest)'")
  }
  ## anonymous not allowed
  if (inherits(object, "anonymous")) {
    stop("this function does work with anonymous forests")
  }
  ## jitt/data.pass
  jitt <- is.hidden.jitt(user.option, "none", "na.omit", FALSE, partial = TRUE)
  data.pass <- is.hidden.data.pass(user.option)
  ## acquire the forest --> hereafter the object is the forest
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
  ## uniquify the x-var
  partial.values <- sort(unique(partial.values))
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
  ## Determine the immutable yvar factor map
  yfactor <- object$yvar.factor
  ## multivariate family details
  ##m.target.idx <- get.outcome.target(family, yvar.names, NULL)
  m.target.idx <- 1:length(yvar.names)
  ## short cut to get the y-outcome type and number of levels
  yvar.types <- yfactor$types
  yvar.nlevels  <- yfactor$nlevels
  yvar.numeric.levels  <- yfactor$numeric.levels
  ## recover the individual subject identifiers, if they exist.
  subj <- object$subj
  ## Get event information for survival families.
  event.info <- object$event.info
  event.type  <- event.info$event.type
  ## CR.bits assignment.
  cr.bits <- get.cr.bits(family)
  ## Determine the immutable xvar factor map.
  xfactor <- object$xvar.factor
  ## short cut to get the x-variable type and number of levels
  xvar.types <- xfactor$types
  xvar.nlevels <- xfactor$nlevels
  xvar.numeric.levels <- xfactor$numeric.levels
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
  ## Recover the case weights and bootstrap sample. 
  case.wt <- object$case.wt
  samp <- object$samp
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
  partial.bits <- get.partial.bits(length(partial.values))
  terminal.qualts.bits <- get.terminal.qualts.predict(object$terminal.qualts)
  terminal.quants.bits <- get.terminal.quants.predict(object$terminal.quants)
  data.pass.bits <- get.data.pass(data.pass)
  ## Just in Time Tree Topology Restoration.  Note that we do not
  ## support augmented data, SGT, TDC, here, but we do no checks for
  ## these. We defer to the user option for JITT.
  jitt.bits <- get.jitt(jitt)
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
  ## The pivot in this predict related function is different from
  ## the pivot used in the grow related function.  In rfsrc we are
  ## referencing the list nativeOutput[[]].  Here we are referencing
  ## the $nativeArray[[]] object which is a massaged version of
  ## nativeOutput[[]].  Here, pivot points to the record PRIOR to
  ## parmID2.  We adjust for the presence of interactions, and
  ## synthetics.  The default chunk is parmIDx, contPTx, contPTRx,
  ## mwcpSZx, fsrecIDx.
  ## WARNING: Note that the maximum number of slots in the following
  ## foreign function call is 64.  Ensure that this limit is not
  ## exceeded.  Otherwise, the program will error on the call.
  if (hdim == 0) {
      pivot = 0
      chunk = 0
  } else {
      pivot <- which(names(object$nativeArray) == "contPTR")
      chunk = 5
      if (!is.null(object$base.learner)) {
          if (object$base.learner$interact.depth > 1) {
              ## Adjusted for pairCT, augmXone2, augmXtwo2 above, but the chunck size needs to increased.
              ## pivot = pivot + 3
              chunk = chunk + 2
          }
          if (object$base.learner$synthetic.depth > 1) {
              ## Adjusted for sythSZ, augmXS2 above, but the chunck size needs to increased.
              ## pivot = pivot + 2
              chunk = chunk + 1
          }
      }
  }
  nativeOutput <- tryCatch({.Call("rfsrcPredict",
                                  as.integer(do.trace),
                                  as.integer(seed),
                                  as.integer(ensemble.bits +
                                             bootstrap.bits +
                                             cr.bits),
                                  as.integer(samptype.bits +
                                             na.action.bits +
                                             terminal.qualts.bits +
                                             terminal.quants.bits +
                                             partial.bits + 
                                             data.pass.bits +
                                             jitt.bits),
                                  ## >>>> start of maxi forest object >>>>
                                  as.numeric(1),
                                  as.integer(ntree),
                                  as.integer(n),
                                  list(as.integer(length(yvar.types)),
                                       if (is.null(yvar.types)) NULL else as.character(yvar.types),
                                       if (is.null(yvar.types)) NULL else as.integer(yvar.nlevels),
                                       if (is.null(yvar.numeric.levels)) NULL else sapply(1:length(yvar.numeric.levels), function(nn) {as.integer(length(yvar.numeric.levels[[nn]]))}),
                                       if (is.null(subj)) NULL else as.integer(subj),
                                       if (is.null(event.type)) NULL else as.integer(length(event.type)),
                                       if (is.null(event.type)) NULL else as.integer(event.type)),
                                  if (is.null(yvar.numeric.levels)) {
                                      NULL
                                  }
                                  else {
                                      lapply(1:length(yvar.numeric.levels),
                                             function(nn) {as.integer(yvar.numeric.levels[[nn]])})
                                  },
                                  if (is.null(yvar.types)) NULL else as.double(as.vector(yvar)),
                                  list(as.integer(n.xvar),
                                       as.character(xvar.types),
                                       if (is.null(xvar.types)) NULL else as.integer(xvar.nlevels),
                                       if (is.null(xvar.numeric.levels)) NULL else sapply(1:length(xvar.numeric.levels), function(nn) {as.integer(length(xvar.numeric.levels[[nn]]))}),
                                       NULL,
                                       NULL),
                                  if (is.null(xvar.numeric.levels)) {
                                      NULL
                                  }
                                  else {
                                      lapply(1:length(xvar.numeric.levels),
                                             function(nn) {as.integer(xvar.numeric.levels[[nn]])})
                                  },
                                  as.double(as.vector(xvar)),
                                  list(as.integer(length(case.wt)),
                                       if (is.null(case.wt)) NULL else as.double(case.wt),
                                       as.integer(sampsize),
                                       if (is.null(samp)) NULL else as.integer(samp)),
                                  list(if(is.null(event.info$time.interest)) as.integer(0) else as.integer(length(event.info$time.interest)),
                                       if(is.null(event.info$time.interest)) NULL else as.double(event.info$time.interest)),
                                  as.integer(object$totalNodeCount),
                                  as.integer(object$leafCount),
                                  list(as.integer(object$seed),
                                       if (is.null(object$seedVimp)) NULL else as.integer(object$seedVimp),
                                       as.integer(object$optLoGrow)),
                                  as.integer(hdim),
                                  ## Object containing base learner settings, this is never NULL.
                                  object$base.learner, 
                                  as.integer((object$nativeArray)$treeID),
                                  as.integer((object$nativeArray)$nodeID),
                                  as.integer((object$nativeArray)$nodeSZ),
                                  as.integer((object$nativeArray)$brnodeID),
                                  ## This is hc_zero.  It is never NULL.
                                  list(as.integer((object$nativeArray)$parmID),
                                  as.double((object$nativeArray)$contPT),
                                  as.integer((object$nativeArray)$mwcpSZ),
                                  as.integer((object$nativeArray)$fsrecID),
                                  as.integer((object$nativeFactorArray)$mwcpPT)),
                                  ## This slot is hc_one_augm_intr.  This slot can be NULL.
                                  if (!is.null(object$base.learner)) {
                                      if (object$base.learner$interact.depth > 1) {
                                          list(as.integer((object$nativeArray)$pairCT),
                                               as.integer((object$nativeArray)$augmXone),
                                               as.integer((object$nativeArray)$augmXtwo))
                                      } else { NULL }
                                  } else { NULL },
                                  ## This slot is hc_one_augm_syth.  This slot can be NULL.
                                  if (!is.null(object$base.learner)) {
                                      if (object$base.learner$synthetic.depth > 1) {
                                          list(as.integer((object$nativeArray)$sythSZ),
                                               as.integer((object$nativeArray)$augmXS))
                                      } else { NULL }
                                  } else { NULL },
                                  
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL,
                                  
                                   
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
                                  
                                  NULL,
                                  
                                    
                                  list(as.integer(0), NULL), ## Importance disabled.
                                  ## Partial variables enabled.  Note the as.integer is needed.
                                  list(as.integer(get.partial.type(family, partial.type)),
                                       as.integer(which(xvar.names == partial.xvar)),
                                       as.integer(length(partial.values)),
                                       as.double(partial.values),
                                       as.integer(length(partial.xvar2)),
                                       if (length(partial.xvar2) == 0) NULL else as.integer(match(partial.xvar2, xvar.names)),
                                       as.double(partial.values2)),
                                  as.integer(0),    ## New data disabled.
                                  as.integer(0),    ## New data disabled.
                                  as.double(NULL),  ## New data disabled.
                                  as.double(NULL),  ## New data disabled.
                                  as.integer(ntree), ## block.size is hard-coded.
                                  list(as.integer(0),
                                       NULL,
                                       as.double(0)), ## Quantiles disabled.
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
  ## append other useful information here
  rfsrcOutput$partial.values <- partial.values
  rfsrcOutput$yvar.names <- object$yvar.names
  rfsrcOutput$event.info <- object$event.info 
  ## return the goodies
  class(rfsrcOutput) <- c("rfsrc", "partial",   family)
  return (rfsrcOutput)
}
partial <- partial.rfsrc
##################################################################
##
##
## work-horse for converting raw data to useable format for plots
##
##
##################################################################
get.partial.plot.data <- function(o, target, m.target = NULL) {
  ##--------------------------------------------------------------
  ##
  ## coherency checks
  ##
  ##--------------------------------------------------------------
  if (sum(inherits(o, c("rfsrc", "partial"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(rfsrc, partial)'")
  }
  ##--------------------------------------------------------------
  ##
  ## extract the partial data, family specific
  ##
  ##--------------------------------------------------------------
  family <- class(o)[3]
  yvar.names <- o$yvar.names
  ## multivariate families
  if (family == "regr+" | family == "class+" | family == "mix+") {
    ## determine the y-outcome to be used
    if (missing(m.target)) {
      m.target <- yvar.names[1]
    }
    else {
      m.target <- match.arg(m.target, yvar.names)
    }      
    ## pull the data
    pdta  <- unlist(list(o$classOutput, o$regrOutput), recursive = FALSE)[[m.target]]
    ## convert to either regression or classification
    ## classification
    if (length(dim(pdta)) > 2) {
      ## identify target class level
      if (missing(target)) {
        target <- 1
      }
      else {
        target <- match(match.arg(target, yvar.names), yvar.names)
      }
      pdta <- pdta[, 1 + target, ]
      pdta[is.infinite(pdta)] <- NA
      yhat <- colMeans(pdta, na.rm = TRUE)    
    }
    ## regression
    else {
      pdta[is.infinite(pdta)] <- NA
      yhat <- colMeans(pdta, na.rm = TRUE)
    }
  }
  ## regression
  else if (family == "regr") {
    pdta <- o$regrOutput[[1]]
    pdta[is.infinite(pdta)] <- NA
    yhat <- colMeans(pdta, na.rm = TRUE)    
  }
  ## classification
  else if (family == "class") {
    pdta <- o$classOutput[[1]]
    yvar.names <- colnames(pdta[,,1])[-1]
    if (missing(target)) {
      target <- 1
    }
    else {
      target <- match(match.arg(target, yvar.names), yvar.names)
    }
    pdta <- pdta[, 1 + target, ]
    pdta[is.infinite(pdta)] <- NA
    yhat <- colMeans(pdta, na.rm = TRUE)    
  }
  ## survival
  else if (family == "surv") {
    pdta <- o$survOutput
    pdta[is.infinite(pdta)] <- NA
    ## surv, chf
    if (length(dim(pdta)) > 2) {
      yhat <- do.call(cbind, lapply(1:dim(pdta)[2], function(k) {
        colMeans(pdta[, k, ], na.rm = TRUE)
      }))
    }
    ## mortality
    else {
      yhat <- colMeans(pdta, na.rm = TRUE)
    }
  }
  ## competing risk
  else if (family == "surv-CR") {
    pdta <- o$survOutput
    pdta[is.infinite(pdta)] <- NA
    ## identify target cause specific event
    if (missing(target)) {
        target <- 1
      }
      else {
        target <- match(target, o$event.info$event.type)
      }
    ## cif, chf
    if (length(dim(pdta)) > 3) {
      yhat <- do.call(cbind, lapply(1:dim(pdta)[2], function(k) {
        colMeans(pdta[, k, target, ], na.rm = TRUE)
      }))
    }
    ## years-lost
    else {
      yhat <- colMeans(pdta[, target, ], na.rm = TRUE)
    }
  }
  else {
    stop("family specified is not supported")
  }
  ##--------------------------------------------------------------
  ##
  ## return the goodies
  ##
  ##--------------------------------------------------------------
  list(x = o$partial.values, yhat = yhat, partial.time = o$partial.time)
}
##################################################################
##
##
## utility functions
##
##
##################################################################
get.partial.bits <- function (partial.length) {
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
get.partial.type <- function (family, partial.type) {
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
