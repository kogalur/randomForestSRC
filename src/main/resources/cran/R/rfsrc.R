rfsrc <- function(formula, data, ntree = 1000,
                  mtry = NULL, ytry = NULL,
                  nodesize = NULL, nodedepth = NULL,
                  splitrule = NULL, nsplit = 10,
                  importance = c(FALSE, TRUE, "none", "permute", "random", "anti"),
                  block.size = if (importance == "none" || as.character(importance) == "FALSE") NULL else 10,
                  ensemble = c("all", "oob", "inbag"),
                  bootstrap = c("by.root", "by.node", "none", "by.user"),
                  samptype = c("swr", "swor"), sampsize = NULL, samp = NULL, membership = FALSE,
                  na.action = c("na.omit", "na.impute"), nimpute = 1,
                  ntime, cause,
                  proximity = FALSE, distance = FALSE, forest.wt = FALSE,
                  xvar.wt = NULL, yvar.wt = NULL, split.wt = NULL, case.wt = NULL, 
                  forest = TRUE,
                  var.used = c(FALSE, "all.trees", "by.tree"),
                  split.depth = c(FALSE, "all.trees", "by.tree"),
                  seed = NULL,
                  do.trace = FALSE,
                  statistics = FALSE,
                   ...)
{
  univariate.nomenclature = TRUE
  ## get any hidden options
  user.option <- list(...)
  impute.only <- is.hidden.impute.only(user.option)
  terminal.qualts <- is.hidden.terminal.qualts(user.option)
  terminal.quants <- is.hidden.terminal.quants(user.option)
  ## TBD TBD make perf.type visible TBD TBD
  perf.type <- is.hidden.perf.type(user.option)
  rfq <- is.hidden.rfq(user.option)
  htry <- is.hidden.htry(user.option)
  ## verify key options
  ensemble <- match.arg(ensemble, c("all", "oob", "inbag"))  
  bootstrap <- match.arg(bootstrap, c("by.root", "by.node", "none", "by.user"))
  if (bootstrap == "by.node" || bootstrap == "none") {
    ensemble <- "inbag"
  }
  importance <- match.arg(as.character(importance), c(FALSE, TRUE, "none", "permute", "random", "anti"))
  na.action <- match.arg(na.action, c("na.omit", "na.impute"))
  proximity <- match.arg(as.character(proximity), c(FALSE, TRUE, "inbag", "oob", "all"))
  distance <- match.arg(as.character(distance), c(FALSE, TRUE, "inbag", "oob", "all"))
  var.used <- match.arg(as.character(var.used), c("FALSE", "all.trees", "by.tree"))
  split.depth <- match.arg(as.character(split.depth),  c("FALSE", "all.trees", "by.tree"))
  if (var.used == "FALSE") var.used <- FALSE
  if (split.depth == "FALSE") split.depth <- FALSE
  ## data cannot be missing
  if (missing(data)) stop("data is missing")
  ## conduct preliminary formula validation
  if (missing(formula) | (!missing(formula) && is.null(formula))) {
    if (is.null(ytry)) {
      formula <- as.formula("Unsupervised() ~ .")
    }
    else {
      formula <- as.formula(paste("Unsupervised(", ytry, ")~."))
    }
  } 
  formulaPrelim <- parseFormula(formula, data, ytry)
  ## save the call/formula for the return object
  my.call <- match.call()
  my.call$formula <- eval(formula)
  ## conduct preliminary processing of missing data
  ## record whether there's no missing data: efficiency step
  if (any(is.na(data))) {
    data <- parseMissingData(formulaPrelim, data)
    miss.flag <- TRUE
  }
    else {
      miss.flag <- FALSE
    }
  ## finalize the formula based on the pre-processed data
  formulaDetail <- finalizeFormula(formulaPrelim, data)
  ## coherence checks on option parameters
  ntree <- round(ntree)
  if (ntree < 1) stop("Invalid choice of 'ntree'.  Cannot be less than 1.")
  if (!is.null(nodesize) && nodesize < 1) stop("Invalid choice of 'nodesize'. Cannot be less than 1.")
  if (!is.null(nodedepth)) nodedepth = round(nodedepth) else nodedepth = -1
  nimpute <- round(nimpute)
  if (nimpute < 1) stop("Invalid choice of 'nimpute'.  Cannot be less than 1.")
  ## initialize the seed
  seed <- get.seed(seed)
  ## save the family for convenient access
  family <- formulaDetail$family
  ## save the names for convenient access
  xvar.names <- formulaDetail$xvar.names
  yvar.names <- formulaDetail$yvar.names
  ## reality check on x and y
  ## are there any x-variables?  (can happen when for multivariate formula)
  if (length(xvar.names) == 0) {
    stop("something seems wrong: your formula did not define any x-variables")
  }
  ## .. are there any y-variables?  (do not test for the unsupervised case)
  if (family != "unsupv" && length(yvar.names) == 0) {
    stop("something seems wrong: your formula did not define any y-variables")
  }
  ## missing levels are allowed.
  if (family == "class") {
    if (length(setdiff(levels(data[, yvar.names]), unique(data[, yvar.names]))) > 0) {
      warning("empty classes found when implementing classification\n")
    }
  }
  ## mark missing factor levels as NA.
  data <- rm.na.levels(data, xvar.names)
  data <- rm.na.levels(data, yvar.names)
  ## Determine the immutable yvar factor map.
  yfactor <- extract.factor(data, yvar.names)
  ## Determine the immutable xvar factor map.
  xfactor <- extract.factor(data, xvar.names)
  ## get the y-outcome type and nlevels
  yvar.types <- get.yvar.type(family, yfactor$generic.types, yvar.names)
  yvar.nlevels <- get.yvar.nlevels(family, yfactor$nlevels, yvar.names, data)
  ## get the x-variable type and nlevels
  xvar.types <- get.xvar.type(xfactor$generic.types, xvar.names)
  xvar.nlevels <- get.xvar.nlevels(xfactor$nlevels, xvar.names, data)
  ## Convert the data to numeric mode, apply the na.action protocol.
  data <- finalizeData(c(yvar.names, xvar.names), data, na.action, miss.flag)
  ## Save the row and column names for later overlay
  data.row.names <- rownames(data)
  data.col.names <- colnames(data)
  ## Finalize the xvar matrix.
  xvar <- as.matrix(data[, xvar.names, drop = FALSE])
  rownames(xvar) <- colnames(xvar) <- NULL
  ## Initialize sample size
  ## Set mtry
  ## Set the weight matrix for xvar, and the split
  n <- nrow(xvar)
  n.xvar <- ncol(xvar)
  mtry <- get.grow.mtry(mtry, n.xvar, family)
  samptype <- match.arg(samptype, c("swr", "swor"))
  if (bootstrap == "by.root") {
    ## Set the default value for the size of the bootstrap
    if(missing(sampsize) | is.null(sampsize)) {
      if (samptype == "swr") {
        sampsize <- nrow(xvar)
      }
      if (samptype == "swor") {
        sampsize <- round(nrow(xvar) * (1 - exp(-1)))
      }
    }
    else {
      sampsize <- round(sampsize)
      if (sampsize < 1) {
        stop("sampsize must be greater than zero")
      }
      if (samptype == "swr") {
        ## Size is not limited.
      }
      if (samptype == "swor") {
        ## Size is limited by the number of cases.
        sampsize <- min(sampsize, nrow(xvar))
      }
    }
    samp <- NULL
  }
  else if (bootstrap == "by.user") {    
    ## Check for coherent sample.  It will of be dim [n] x [ntree].
    if (is.null(samp)) {
      stop("samp must not be NULL when bootstrapping by user")
    }
    ## ntree should be coherent with the sample provided
    ntree <- ncol(samp)
    sampsize <- colSums(samp)
    if (sum(sampsize == sampsize[1]) != ntree) {
      stop("sampsize must be identical for each tree")
    }
    sampsize <- sampsize[1]
  }
  else {
    ## Override sample size when not bootstrapping by root or user.
    sampsize = nrow(xvar)
    ## Override sample type when not bootstrapping by root or user.
    samptype <- "swr"
  }
  case.wt  <- get.weight(case.wt, n)
  forest.wt <- match.arg(as.character(forest.wt), c(FALSE, TRUE, "inbag", "oob", "all"))
  split.wt <- get.weight(split.wt, n.xvar)
  ## Note that yvar.types is based on the formula:
  ## If the family is unsupervised, yvar.types is NULL
  ## If the family is regression, classification, or
  ## multivariate/mixed, it is of length greater than zero (0).  In the case
  ## of surival and competing risk, it is specifically of length two (2) for
  ## coherency, though it is ignored in the native-code.
  ## Note that ytry is set by the formula.
  ## If the family is unsupervised, ytry assumes the value specified
  ## via the formula.  If the family is regression, classification, or
  ## multivariate/mixed, it defaults to length(yvar.types).  In the case
  ## of survival and competing risk, it assumes the value of two (2) for
  ## coherency, though it is ignored in the native-code.
  if (family == "unspv") {
    ## Override any incoming weights.
    yvar.wt <- NULL
  }
  else {
    yvar.wt <- get.weight(yvar.wt, length(yvar.types))
  }
  xvar.wt  <- get.weight(xvar.wt, n.xvar)
  ## Get the y-outcome
  yvar <- as.matrix(data[, yvar.names, drop = FALSE])
  if(dim(yvar)[2] == 0) {
    ## Override yvar if it is not present.
    yvar <- NULL
  }
  ## Determine the number of missing values
  if (miss.flag) {
    n.miss <- get.nmiss(xvar, yvar)
  }
    else {
      n.miss <- 0
    }
  ## In impute.only, if there is no missing data,
  ## return the data and do nothing else.
  if (impute.only && n.miss == 0) {
    return(data)
  }
  ## Don't need the data anymore.
  remove(data)
  ## Memory management option
  ## (TBD) Not currently implemented.
  big.data <- FALSE
  ## TBD I think we need to separate event related data and dimensioning
  ## of the response which this function combines. TBD 
  ## Get event information and dimensioning for families
  event.info <- get.grow.event.info(yvar, family, ntime = ntime)
  ## Initialize nsplit, noting that it may be been overridden.
  splitinfo <- get.grow.splitinfo(formulaDetail, splitrule, htry, nsplit, event.info$event.type)
  ## Set the cause weights for the split statistic calculation.
  if (family == "surv" || family == "surv-CR") {
    if (length(event.info$event.type) > 1) {
      ## This may or may not be CR, but we do have multiple events.
      ## This requires associated weights.
      if (missing(cause) || is.null(cause)) {
        cause <- NULL
        cause.wt <- rep(1, length(event.info$event.type))
      }
      else {
        if (length(cause) == 1) {
          ## Event specific selection
          if (cause >= 1 && cause <= length(event.info$event.type)) {
            cause.wt <- rep(0, length(event.info$event.type))
            cause.wt[cause] <- 1
          }
          else {
            cause.wt <- rep(1, length(event.info$event.type))
          }
        }
        else {
          ## The user has specified event specific weights
          if (length(cause) == length(event.info$event.type) && all(cause >= 0) && !all(cause == 0)) {
            cause.wt <- cause / sum(cause)
          }
          else {
            cause.wt <- rep(1, length(event.info$event.type))
          }
        }
      }
    }
    else {
      ## We are conducting non-CR, without multiple events.  Send in a dummy value.
      cause <- NULL
      cause.wt = 1
    }
    ## Set the coerced family.
    family <- get.coerced.survival.fmly(family, event.info$event.type, splitinfo$name)
  }
  else {
    ## We pass in NULL for the cause weight for all other families
    cause <- cause.wt <- NULL
  }
  ## Nodesize determination
  nodesize <- get.grow.nodesize(family, nodesize)
  ## Turn ensemble outputs off unless bootstrapping by root or user.
  if ((bootstrap != "by.root") && (bootstrap != "by.user")) {
    importance <- "none"
    perf.type <- "none"
  }
  ## Turn ensemble outputs off for unsupervised mode only
  if (family == "unsupv") {
    importance <- "none"
    perf.type <- "none"
  }
  ## Impute only mode
  ## Some of this may be redundant
  if (impute.only) {
    forest       <- FALSE
    proximity    <- FALSE
    distance     <- FALSE
    var.used     <- FALSE
    split.depth  <- FALSE
    membership   <- FALSE
    perf.type    <- "none"
    importance   <- "none"
    terminal.qualts <- FALSE
    terminal.quants <- FALSE
    ## We have deleted the following line because it was incorrectly
    ## overriding the user specified option na.action in 
    ## impute.rfsrc( generic.impute.rfsrc( rfsrc() ) ).  Now, the user specified option
    ## is respected.  Note that calls from impute.rfsrc()
    ## set the impute.only bit which results in turning of all outputs
    ## other than imputed data.
    ## na.action    <- "na.impute"
  }
  ## Bit dependencies:
  if (terminal.qualts | terminal.quants) {
    forest <- TRUE
  }
  ## Assign low bits for the native code
  ensemble.bits <- get.ensemble(ensemble)
  impute.only.bits <- get.impute.only(impute.only, n.miss)
  var.used.bits <- get.var.used(var.used)
  split.depth.bits <- get.split.depth(split.depth)
  importance.bits <- get.importance(importance)
  bootstrap.bits <- get.bootstrap(bootstrap)
  forest.bits <- get.forest(forest)
  proximity.bits <- get.proximity(TRUE, proximity)
  distance.bits <- get.distance(TRUE, distance)
  membership.bits <-  get.membership(membership)
  statistics.bits <- get.statistics(statistics)
  split.cust.bits <- get.split.cust(splitinfo$cust)
  ## get performance and rfq bits
  perf.type <- get.perf(perf.type, impute.only, family)
  perf.bits <-  get.perf.bits(perf.type)
  rfq <- get.rfq(rfq)
  rfq.bits <- get.rfq.bits(rfq, family)
  ## Assign high bits for the native code
  samptype.bits <- get.samptype(samptype)
  forest.wt.bits <- get.forest.wt(TRUE, bootstrap, forest.wt)
  na.action.bits <- get.na.action(na.action)
  block.size <- get.block.size(block.size, ntree)
  terminal.qualts.bits <- get.terminal.qualts(terminal.qualts, FALSE)
  terminal.quants.bits <- get.terminal.quants(terminal.quants, FALSE)
  ## Set the trace
  do.trace <- get.trace(do.trace)
  nativeOutput <- tryCatch({.Call("rfsrcGrow",
                                  as.integer(do.trace),
                                  as.integer(seed),
                                  as.integer(ensemble.bits +
                                             impute.only.bits +
                                             var.used.bits +
                                             split.depth.bits +
                                             importance.bits +
                                             bootstrap.bits +
                                             forest.bits +
                                             proximity.bits +
                                             perf.bits +
                                             rfq.bits +
                                             statistics.bits),
                                  as.integer(samptype.bits +
                                             forest.wt.bits +
                                             distance.bits + 
                                             na.action.bits +
                                             split.cust.bits +
                                             membership.bits +
                                             terminal.qualts.bits +
                                             terminal.quants.bits),
                                  as.integer(splitinfo$index),
                                  as.integer(splitinfo$nsplit),
                                  as.integer(mtry),
                                  as.integer(htry),
                                  as.integer(formulaDetail$ytry),
                                  as.integer(nodesize),
                                  as.integer(nodedepth),
                                  as.integer(length(cause.wt)),
                                  as.double(cause.wt),
                                  as.integer(ntree),
                                  as.integer(n),
                                  as.integer(length(yvar.types)),
                                  as.character(yvar.types),
                                  as.integer(yvar.nlevels),
                                  as.double(as.vector(yvar)),
                                  as.integer(n.xvar),
                                  as.character(xvar.types),
                                  as.integer(xvar.nlevels),
                                  as.integer(sampsize),
                                  as.integer(samp),
                                  as.double(case.wt),
                                  as.double(split.wt),
                                  as.double(yvar.wt),
                                  as.double(xvar.wt),
                                  as.double(xvar),
                                  as.integer(length(event.info$time.interest)),
                                  as.double(event.info$time.interest),
                                  as.integer(nimpute),
                                  as.integer(block.size),
                                  as.integer(get.rf.cores()))}, error = function(e) {
                                    print(e)
                                    NULL})
  ## check for error return condition in the native code
  if (is.null(nativeOutput)) {
    if (impute.only) {
      ## in impute mode we proceed and return NULL
      return(NULL)
    }
      else {
        ## currently all other modes will error out with the RFSRC error message
        stop("An error has occurred in the grow algorithm.  Please turn trace on for further analysis.")
      }
  }
  ## check if there was missing data, and assign imputed data if so.
  if (n.miss > 0) {
    imputed.data <- matrix(nativeOutput$imputation, nrow = n.miss, byrow = FALSE)
    imputed.indv <- imputed.data[, 1]
    imputed.data <- as.matrix(imputed.data[, -1, drop = FALSE])
    ##if (n.miss == 1) {
    ##  imputed.data <- t(imputed.data)
    ##}
    nativeOutput$imputation <- NULL
    ## fill NA's in original GROW data with multiply imputed values.
    ## This will now serve as the forest data set and will enable
    ## recovery of terminal node membership with the head of the
    ## seed chain
    if (nimpute > 1) {
      ## the forest was grown using overlaid full (all) summary values
      if (grepl("surv", family)) {
        yvar[imputed.indv, 1] <- imputed.data[, 1]
        yvar[imputed.indv, 2] <- imputed.data[, 2]
        xvar[imputed.indv, ] <- imputed.data[, -c(1:2), drop = FALSE]
      }
        else {
          if (!is.null(yvar.types)) {
            yvar[imputed.indv, ] <- imputed.data[, 1:length(yvar.types), drop = FALSE]
            xvar[imputed.indv, ] <- imputed.data[, -c(1:length(yvar.types)), drop = FALSE]
          }
            else {
              xvar[imputed.indv, ] <- imputed.data
            }
        }
      ## remove the imputed data outputs
      imputed.indv <- NULL
      imputed.data <- NULL
      imputedOOBData <- NULL
      ## the missingness action for the training data is formally safed.  This is saved as part of the
      ## forest object, but the setting is now irrelevant as the training data has no missingness.
      na.action = "na.omit"
    }
      else {
        ## add column names to the imputed data outputs in the absence
        ## of multiple imputation.
        colnames(imputed.data) <- c(yvar.names, xvar.names)
        imputed.data <- as.data.frame(imputed.data)
      }
    ## now map the imputed.data columns back to their original order
    ## commented out because of difficulty in maintaining order across different modes
    ## imputed.data <- imputed.data[, data.col.names]
  }
  ## add row and column names to xvar matrix
  xvar <- as.data.frame(xvar)
  rownames(xvar) <- data.row.names
  colnames(xvar) <- xvar.names
  ## map xvar factors back to original values
  xvar <- map.factor(xvar, xfactor)
  ## add column names to response matrix
  ## does not apply for unsupervised mode
  if (family != "unsupv") {
    yvar <- as.data.frame(yvar)
    colnames(yvar) <- yvar.names
  }
  else {
    yvar <- NULL
  }
  ## map response factors back to original values
  ## does not apply for unsupervised mode
  if (family != "unsupv") {
    if (family == "regr+" | family == "class+" | family == "mix+") {
      yvar <- map.factor(yvar, yfactor)
    }
    else {
      yvar <- amatrix.remove.names(map.factor(yvar, yfactor))
    }
  }
  ## class imbalanced processing
  pi.hat <- NULL
  if (family == "class" && rfq) {
    pi.hat <- table(yvar) / length(yvar)
  }
  ## map imputed data factors back to original values
  ## does NOT apply for y under unsupervised mode
  if ((n.miss > 0) & (nimpute < 2)) {
    imputed.data <- map.factor(imputed.data, xfactor)
    if (family != "unsupv") {
      imputed.data <- map.factor(imputed.data, yfactor)
    }
  }
  ## Define the forest.
  if (forest) {
      ## We allocate, in the native code, the native array to its
      ## theoretical maximum.  This is trimmed below to the actual
      ## size.
      nativeArraySize = 0
      if (htry == 0) {
          mwcpCountSummary = rep (0, 1)
          nativeFactorArray <- vector("list", 1)
      }
      else {
          mwcpCountSummary = rep (0, htry)
          nativeFactorArray <- vector("list", htry)
      }
    ## Marker for start of native forest topology.  This can change with the outputs requested.
    ## For the arithmetic related to the pivot point, you need to refer to stackOutput.c and in
    ## particular, stackForestOutputObjects().
    pivot <- which(names(nativeOutput) == "treeID")
    nullO <- lapply(1:ntree, function(b) {
      if (nativeOutput$leafCount[b] > 0) {
        ## The tree was not rejected.  Count the number of internal
        ## and external (terminal) nodes in the forest.
        nativeArraySize <<- nativeArraySize + (2 * nativeOutput$leafCount[b]) - 1
        mwcpCountSummary[1] <<- mwcpCountSummary[1] + nativeOutput$mwcpCT[b]
          if (htry > 1) {
            for (i in 0:(htry-2)) {
              mwcpCountSummary[i+2] <<- mwcpCountSummary[i+2] + nativeOutput[[pivot + 9 + (5 *(htry-1)) + i]][b]
            }
        }
      }
      else {
        ## The tree was rejected.  However, it acts as a
        ## placeholder, being a stump topologically and thus
        ## adds to the total node count.
        nativeArraySize <<- nativeArraySize + 1
      }
      NULL
    })
    rm(nullO)##memory saving device
    nativeArray <- as.data.frame(cbind(nativeOutput$treeID[1:nativeArraySize],
                                       nativeOutput$nodeID[1:nativeArraySize]))
    nativeArrayHeader <- c("treeID", "nodeID")
    nativeArray <- as.data.frame(cbind(nativeArray,
                                       nativeOutput$parmID[1:nativeArraySize],
                                       nativeOutput$contPT[1:nativeArraySize],
                                       nativeOutput$mwcpSZ[1:nativeArraySize]))
    nativeArrayHeader <- c(nativeArrayHeader, "parmID", "contPT", "mwcpSZ") 
    if (mwcpCountSummary[1] > 0) {
      ## This can be NULL if there are no factor splits along this dimension.
      nativeFactorArray[[1]] <- nativeOutput$mwcpPT[1:mwcpCountSummary[1]]
    }
    nativeFactorArrayHeader <- "mwcpPT"
    if (htry > 0) {
      nativeArray <- as.data.frame(cbind(nativeArray,
                                         nativeOutput[[pivot + 7]][1:nativeArraySize]))
      nativeArrayHeader <- c(nativeArrayHeader, "hcDim")
      nativeArray <- as.data.frame(cbind(nativeArray,
                                         nativeOutput[[pivot + 8]][1:nativeArraySize]))
      nativeArrayHeader <- c(nativeArrayHeader, "contPTR")
    }
      if (htry > 1) {
          for (i in 0:(htry-2)) {
              nativeArray <- as.data.frame(cbind(nativeArray,
                                                 nativeOutput[[pivot + 9 + (0 * (htry-1)) + i]][1:nativeArraySize]))
              nativeArrayHeader <- c(nativeArrayHeader, paste("parmID", i+2, sep=""))
              nativeArray <- as.data.frame(cbind(nativeArray,
                                                 nativeOutput[[pivot + 9 + (1 * (htry-1)) + i]][1:nativeArraySize]))
              nativeArrayHeader <- c(nativeArrayHeader, paste("contPT", i+2, sep=""))
              nativeArray <- as.data.frame(cbind(nativeArray,
                                                 nativeOutput[[pivot + 9 + (2 * (htry-1)) + i]][1:nativeArraySize]))
              nativeArrayHeader <- c(nativeArrayHeader, paste("contPTR", i+2, sep=""))
              nativeArray <- as.data.frame(cbind(nativeArray,
                                                 nativeOutput[[pivot + 9 + (3 * (htry-1)) + i]][1:nativeArraySize]))
              nativeArrayHeader <- c(nativeArrayHeader, paste("mwcpSZ", i+2, sep=""))
              if (mwcpCountSummary[i+2] > 0) {
                  ## This can be NULL if there are no factor splits along this dimension.
                  nativeFactorArray[[i+2]] <- nativeOutput[[pivot + 9 + (4 * (htry-1)) + i]][1:mwcpCountSummary[i+2]]
              }
              nativeFactorArrayHeader <- c(nativeFactorArrayHeader, paste("mwcpPT", i+2, sep=""))
          }
      }
    names(nativeArray) <- nativeArrayHeader
    names(nativeFactorArray) <- nativeFactorArrayHeader
    if (terminal.qualts | terminal.quants) {
      ## The terminal node qualitative and quantitative outputs
      ## are allocated to their theoretical maximum.
      ## Each tree-specific segment has an initialized and uninitialized
      ## partition. It will have [1:leafCount] validly populated, and
      ## [leafCount+1 : treeTheoreticalMaximum] uninitialized.  We parse the
      ## valid segments here and concatenate them.
      temp <- 2 * (nodesize - 1)
      if (sampsize  > temp) { 
        treeTheoreticalMaximum <- sampsize - temp;
      }
        else {
          treeTheoreticalMaximum <- 1;
        }
      forestTheoreticalMaximum <- treeTheoreticalMaximum * ntree
      offset <- 0
      valid.mcnt.indices <- unlist(lapply(1:ntree, function(b) {
        if (b > 1) {
          offset <<- offset + treeTheoreticalMaximum
        }
        (offset + 1):(offset + nativeOutput$leafCount[b])
      }))
      if (terminal.quants) {
        ## family specific additions to the grow object
        if (grepl("surv", family)) {
          offset <- 0
          valid.2D.surv.indices <- unlist(lapply(1:ntree, function(b) {
            if (b > 1) {
              offset <<- offset + (treeTheoreticalMaximum * length(event.info$event.type) * length(event.info$time.interest))
            }
            (offset + 1):(offset + nativeOutput$leafCount[b] * length(event.info$event.type) * length(event.info$time.interest))
          }))
          offset <- 0
          valid.1D.surv.indices <- unlist(lapply(1:ntree, function(b) {
            if (b > 1) {
              offset <<- offset + (treeTheoreticalMaximum * length(event.info$time.interest))
            }
            (offset + 1):(offset + nativeOutput$leafCount[b] * length(event.info$time.interest))
          }))
          offset <- 0
          valid.mort.indices <- unlist(lapply(1:ntree, function(b) {
            if (b > 1) {
              offset <<- offset + (treeTheoreticalMaximum * length(event.info$event.type))
            }
            (offset + 1):(offset + nativeOutput$leafCount[b] * length(event.info$event.type))
          }))
        }
        else {
          ## This will pick up all "C" and "I".
          class.index <- which(yvar.types != "R")
          class.count <- length(class.index)
          regr.index <- which(yvar.types == "R")
          regr.count <- length(regr.index)
          if (class.count > 0) {
            ## Vector to hold the number of levels in each factor response. 
            levels.count <- array(0, class.count)
            counter <- 0
            for (i in class.index) {
              counter <- counter + 1
              ## Note that [i] is the actual index of the y-variables and not a sequential iterator.
              ## The sequential iteratior is [counter]
              levels.count[counter] <- yvar.nlevels[i]
            }
            offset <- 0
            valid.clas.indices <- unlist(lapply(1:ntree, function(b) {
              if (b > 1) {
                offset <<- offset + (treeTheoreticalMaximum * sum(levels.count))
              }
              (offset + 1):(offset + nativeOutput$leafCount[b] * sum(levels.count))
            }))
          }
          if (regr.count > 0) {
            offset <- 0
            valid.regr.indices <- unlist(lapply(1:ntree, function(b) {
              if (b > 1) {
                offset <<- offset + (treeTheoreticalMaximum * regr.count)
              }
              (offset + 1):(offset + nativeOutput$leafCount[b] * regr.count)
            }))
          }
        }
      }
      nativeArrayTNDS <- list(if(!is.null(nativeOutput$tnSURV)) nativeOutput$tnSURV[valid.1D.surv.indices] else NULL,
                              if(!is.null(nativeOutput$tnMORT)) nativeOutput$tnMORT[valid.mort.indices] else NULL,
                              if(!is.null(nativeOutput$tnNLSN)) nativeOutput$tnNLSN[valid.1D.surv.indices] else NULL,
                              if(!is.null(nativeOutput$tnCSHZ)) nativeOutput$tnCSHZ[valid.2D.surv.indices] else NULL,
                              if(!is.null(nativeOutput$tnCIFN)) nativeOutput$tnCIFN[valid.2D.surv.indices] else NULL,
                              if(!is.null(nativeOutput$tnREGR)) nativeOutput$tnREGR[valid.regr.indices] else NULL,
                              if(!is.null(nativeOutput$tnCLAS)) nativeOutput$tnCLAS[valid.clas.indices] else NULL,
                              nativeOutput$rmbrMembership,
                              nativeOutput$ambrMembership,
                              nativeOutput$tnRCNT[valid.mcnt.indices],
                              nativeOutput$tnACNT[valid.mcnt.indices]);
      names(nativeArrayTNDS) <- c("tnSURV","tnMORT","tnNLSN","tnCSHZ","tnCIFN","tnREGR","tnCLAS", "tnRMBR", "tnAMBR", "tnRCNT", "tnACNT")
    }
    else {
      nativeArrayTNDS <- NULL
    }
    forest.out <- list(htry = htry,
                       nativeArray = nativeArray,
                       nativeFactorArray = nativeFactorArray,
                       totalNodeCount = dim(nativeArray)[1],
                       nodesize = nodesize,
                       nodedepth = nodedepth,
                       ntree = ntree,
                       family = family,
                       splitrule = splitinfo$name,
                       yvar = yvar,
                       yvar.names = yvar.names,
                       xvar = xvar,
                       xvar.names = xvar.names,
                       seed = nativeOutput$seed,
                       bootstrap = bootstrap,
                       sampsize = sampsize,
                       samptype = samptype,
                       samp     = samp,
                       case.wt  = case.wt,
                       terminal.qualts = terminal.qualts,
                       terminal.quants = terminal.quants,
                       nativeArrayTNDS = nativeArrayTNDS,
                       version = "@PROJECT_VERSION@",
                       na.action = na.action,
                       perf.type = perf.type,
                       rfq = rfq,
                       block.size = block.size)
    ## family specific additions to the forest object
    if (grepl("surv", family)) {
      forest.out$time.interest <- event.info$time.interest
    }
    ## Initialize the default class of the forest.
    class(forest.out) <- c("rfsrc", "forest", family)
    if (big.data) {
      class(forest.out) <- c(class(forest.out), "bigdata")
    }
  }
  ## the forest is NULL (the user has requested not to save the forest)
  else {
    forest.out <- NULL
  }
  ## process the proximity matrix
  if (proximity != FALSE) {
    proximity.out <- matrix(0, n, n)
    count <- 0
    for (k in 1:n) {
      proximity.out[k,1:k] <- nativeOutput$proximity[(count+1):(count+k)]
      proximity.out[1:k,k] <- proximity.out[k,1:k]
      count <- count + k
    }
    nativeOutput$proximity <- NULL
  }
  else {
    proximity.out <- NULL
  }
  ## process the proximity matrix
  if (distance != FALSE) {
    distance.out <- matrix(0, n, n)
    count <- 0
    for (k in 1:n) {
      distance.out[k,1:k] <- nativeOutput$distance[(count+1):(count+k)]
      distance.out[1:k,k] <- distance.out[k,1:k]
      count <- count + k
    }
    nativeOutput$distance <- NULL
  }
  else {
    distance.out <- NULL
  }
  ## process weight matrix
  if (forest.wt != FALSE) {
    forest.wt.out <- matrix(nativeOutput$weight, c(n, n), byrow = TRUE)
    nativeOutput$weight <- NULL
  }
  else {
    forest.wt.out <- NULL
  }
  ## membership
  if (membership) {
    membership.out <- matrix(nativeOutput$nodeMembership, c(n, ntree))
    inbag.out <- matrix(nativeOutput$bootMembership, c(n, ntree))
    nativeOutput$nodeMembership <- NULL
    nativeOutput$bootMembership <- NULL
  }
    else {
      membership.out <- NULL
      inbag.out <- NULL
    }
  ## variables used
  if (var.used != FALSE) {
    if (var.used == "all.trees") {
      var.used.out <- nativeOutput$varUsed
      names(var.used.out) <- xvar.names
    }
      else {
        var.used.out <- matrix(nativeOutput$varUsed, nrow = ntree, byrow = TRUE)
        colnames(var.used.out) <- xvar.names
      }
    nativeOutput$varUsed <- NULL
  }
    else {
      var.used.out <-  NULL
    }
  ## split depth
  if (split.depth != FALSE) {
    if (split.depth == "all.trees") {
      split.depth.out <- array(nativeOutput$splitDepth, c(n, n.xvar))
    }
    else {
      split.depth.out <- array(nativeOutput$splitDepth, c(n, n.xvar, ntree))
    }
    nativeOutput$splitDepth <- NULL
  }
  else {
    split.depth.out <-  NULL
  }
  ## node statistics
  if (statistics) {
    node.stats <- as.data.frame(cbind(nativeOutput$spltST[1:nativeArraySize]))
    colnames(node.stats) <- "spltST"
    node.mtry.stats <- t(array(nativeOutput$mtryST[1:nativeArraySize], c(mtry, forest.out$totalNodeCount)))
    node.mtry.index <- t(array(nativeOutput$mtryID[1:nativeArraySize], c(mtry, forest.out$totalNodeCount)))
    ## In the case of unsupervised splitting, output additional
    ## statistics.
    ## (TBD TBD) We need to output ytry related statistics in
    ## the supervised case as well, following the introduction of ytry in
    ## supervised settings.
    if (family == "unsupv") {
      node.ytry.index <- t(array(nativeOutput$uspvST[1:nativeArraySize], c(formulaDetail$ytry, forest.out$totalNodeCount)))
    }
      else {
        node.ytry.index <- NULL
      }
  }
    else {
      node.stats      <- NULL
      node.mtry.stats <- NULL
      node.mtry.index <- NULL
      node.ytry.index <- NULL
    }
  ## make the output object
  rfsrcOutput <- list(
    call = my.call,
    family = family,
    n = n,
    ntree = ntree,
    nimpute = nimpute,
    mtry = mtry,
    nodesize = nodesize,
    nodedepth = nodedepth,
    nsplit = splitinfo$nsplit,
    yvar = yvar,
    yvar.names = yvar.names,
    xvar = xvar,
    xvar.names = xvar.names,
    xvar.wt = xvar.wt,
    split.wt = split.wt,
    cause.wt = cause.wt,
    leaf.count = nativeOutput$leafCount,
    proximity = proximity.out,
    forest = forest.out,
    forest.wt = forest.wt.out,    
    distance = distance.out,
    membership = membership.out,
    splitrule = splitinfo$name,
    inbag = inbag.out,
    var.used = var.used.out,
    imputed.indv = (if (n.miss > 0) imputed.indv else NULL),
    imputed.data = (if (n.miss > 0) imputed.data else NULL),
    split.depth  = split.depth.out,
    node.stats      = node.stats,
    node.mtry.stats = node.mtry.stats,
    node.mtry.index = node.mtry.index,
    node.ytry.index = node.ytry.index,
    ensemble = ensemble,
    block.size = block.size
  )
  ## memory management
  remove(yvar)
  remove(xvar)
  nativeOutput$leafCount <- NULL
  remove(proximity.out)
  remove(forest.out)
  remove(forest.wt.out)
  remove(distance.out)
  remove(membership.out)
  remove(inbag.out)
  remove(var.used.out)
  if (n.miss > 0) remove(imputed.indv)
  if (n.miss > 0) remove(imputed.data)
  remove(split.depth.out)
  ## save the outputs
  survOutput <- NULL
  classOutput <- NULL
  regrOutput <- NULL
  ## EFFICIENCY EFFICIENCY EFFICIENCY
  ## for efficiency - the following does not need to be executed in impute.only mode
  if (!impute.only) {
    ## family specific additions to the grow object
    if (grepl("surv", family)) {
        if ((length(event.info$event.type) > 1) &&
            (splitinfo$name != "l2.impute") &&
            (splitinfo$name != "logrankscore")) {
        coerced.event.count <- length(event.info$event.type)
      }
        else {
          coerced.event.count <- 1
        }
      if (family == "surv") {
        ## Right Censored names.
        ens.names <- list(NULL, NULL)
        mortality.names <- list(NULL, NULL)
        err.names <- list(NULL, NULL)
        vimp.names <- list(NULL, xvar.names)
      }
        else {
          ## Competing Risk names.
          ens.names <- list(NULL, NULL, c(paste("condCHF.", 1:length(event.info$event.type), sep = "")))
          mortality.names <- list(NULL, paste("event.", 1:length(event.info$event.type), sep = ""))
          cif.names <- list(NULL, NULL, c(paste("CIF.", 1:length(event.info$event.type), sep = "")))
          err.names <- list(c(paste("event.", 1:length(event.info$event.type), sep = "")), NULL)
          vimp.names <- list(paste("event.", 1:length(event.info$event.type), sep = ""), xvar.names)
        }
      ## From the native code:
      ##   "allEnsbCHF"
      ##   "oobEnsbCHF"
      ## -> of dim [length(event.info$event.type)] x [RF_sortedTimeInterestSize] x [n]
      ##    where [length(event.info$event.type)] may be equal to [1].
      ## To the R code:
      ## -> of dim [n] x [RF_sortedTimeInterestSize] x [length(event.info$event.type)]  
      chf <- (if (!is.null(nativeOutput$allEnsbCHF))
                adrop3d.last(array(nativeOutput$allEnsbCHF,
                                   c(n, length(event.info$time.interest), length(event.info$event.type)),
                                   dimnames=ens.names), coerced.event.count) else NULL)
      nativeOutput$allEnsbCHF <- NULL
      survOutput <- list(chf = chf)
      remove(chf)
      chf.oob <- (if (!is.null(nativeOutput$oobEnsbCHF))
                    adrop3d.last(array(nativeOutput$oobEnsbCHF,
                                       c(n, length(event.info$time.interest), length(event.info$event.type)),
                                       dimnames=ens.names), coerced.event.count) else NULL)
      nativeOutput$oobEnsbCHF <- NULL
      survOutput = c(survOutput, chf.oob = list(chf.oob))
      remove(chf.oob)
      ## From the native code:
      ##   "allEnsbMRT"
      ##   "oobEnsbMRT"
      ## -> of dim [length(event.info$event.type)] x [n]
      ## To the R code:
      ## -> of dim [n] x [length(event.info$event.type)] 
      predicted <- (if (!is.null(nativeOutput$allEnsbMRT))
                      adrop2d.last(array(nativeOutput$allEnsbMRT,
                                         c(n, length(event.info$event.type)), dimnames=mortality.names), coerced.event.count) else NULL)
      nativeOutput$allEnsbMRT <- NULL
      survOutput = c(survOutput, predicted = list(predicted))
      remove(predicted)
      predicted.oob <- (if (!is.null(nativeOutput$oobEnsbMRT))
                          adrop2d.last(array(nativeOutput$oobEnsbMRT,
                                             c(n, length(event.info$event.type)), dimnames=mortality.names), coerced.event.count) else NULL)
      nativeOutput$oobEnsbMRT <- NULL
      survOutput <- c(survOutput, predicted.oob = list(predicted.oob))
      remove(predicted.oob)
      ## From the native code:
      ##   "allEnsbSRV"
      ##   "oobEnsbSRV"
      ## -> of dim [RF_sortedTimeInterestSize] x [n]
      ## To the R code:
      ## -> of dim [n] x [RF_sortedTimeInterestSize]
      survival <-  (if (!is.null(nativeOutput$allEnsbSRV))
                      matrix(nativeOutput$allEnsbSRV,
                             c(n, length(event.info$time.interest))) else NULL)
      nativeOutput$allEnsbSRV <- NULL
      survOutput <- c(survOutput, survival = list(survival))
      remove(survival)
      survival.oob <-  (if (!is.null(nativeOutput$oobEnsbSRV))
                          matrix(nativeOutput$oobEnsbSRV,
                                 c(n, length(event.info$time.interest))) else NULL)
      nativeOutput$oobEnsbSRV <- NULL
      survOutput <- c(survOutput, survival.oob = list(survival.oob))
      remove(survival.oob)
      ## From the native code:
      ##   "allEnsbCIF"
      ##   "oobEnsbCIF"
      ##   -> of dim [length(event.info$event.type)] x [RF_sortedTimeInterestSize] x [n]
      ## To the native code:
      ##   -> of dim  [n] x [RF_sortedTimeInterestSize] x [length(event.info$event.type)]
      cif <- (if (!is.null(nativeOutput$allEnsbCIF))
                array(nativeOutput$allEnsbCIF,
                      c(n, length(event.info$time.interest), length(event.info$event.type)),
                      dimnames=cif.names) else NULL)
      nativeOutput$allEnsbCIF <- NULL
      survOutput <- c(survOutput, cif = list(cif))
      remove(cif)
      cif.oob <- (if (!is.null(nativeOutput$oobEnsbCIF))
                    array(nativeOutput$oobEnsbCIF,
                          c(n, length(event.info$time.interest), length(event.info$event.type)),
                          dimnames=cif.names) else NULL)
      nativeOutput$oobEnsbCIF <- NULL
      survOutput = c(survOutput, cif.oob = list(cif.oob))
      remove(cif.oob)
      ## From the native code:
      ##   "perfSurv"
      ##   -> of dim [ntree] x length(event.info$event.type)]
      ## To the R code:
      ##   -> of dim [ntree] x length(event.info$event.type)]
      if (!is.null(nativeOutput$perfSurv)) {
        err.rate <- adrop2d.first(array(nativeOutput$perfSurv,
                                        c(length(event.info$event.type), ntree),
                                        dimnames=err.names),
                                  coerced.event.count)
        nativeOutput$perfSurv <- NULL
        if (family == "surv-CR") {
          survOutput = c(survOutput, err.rate = list(t(err.rate)))
        }
          else {
            survOutput = c(survOutput, err.rate = list(err.rate))
          }
        remove(err.rate)
      }
      ## From the native code:
      ##   "blockSurv"
      ##   -> of dim [block.cnt] x length(event.info$event.type)]
      ## To the R code:
      ##   -> of dim [block.cnt] x length(event.info$event.type)]
      if (!is.null(nativeOutput$blockSurv)) {
        err.block.rate <- adrop2d.first(array(nativeOutput$blockSurv,
                                              c(length(event.info$event.type), floor(ntree/block.size)),
                                              dimnames=err.names),
                                        coerced.event.count)
        nativeOutput$blockSurv <- NULL
        if (family == "surv-CR") {
          survOutput = c(survOutput, err.block.rate = list(t(err.block.rate)))
        }
          else {
            survOutput = c(survOutput, err.block.rate = list(err.block.rate))
          }
        remove(err.block.rate)
      }
      ## From the native code:
      ##   "vimpSurv"
      ##   -> of dim [n.xvar] x length(event.info$event.type)]
      ## To the R code:
      ##   -> of dim length(event.info$event.type)] x [n.xvar]
      if (!is.null(nativeOutput$vimpSurv)) {
        importance <- adrop2d.first(array(nativeOutput$vimpSurv,
                                          c(length(event.info$event.type), n.xvar),
                                          dimnames = vimp.names),
                                    coerced.event.count)
        nativeOutput$vimpSurv <- NULL
        if (family == "surv-CR") {
          survOutput = c(survOutput, importance = list(t(importance)))
        }
          else {
            survOutput = c(survOutput, importance = list(importance))
          }
        remove(importance)
      }
      survOutput = c(
        survOutput, list(
          time.interest = event.info$time.interest,
          ndead = sum(na.omit(event.info$cens) != 0))
      )
      ## When TRUE we revert to univariate nomenclature for all the outputs.
      if(univariate.nomenclature) {
        rfsrcOutput <- c(rfsrcOutput, survOutput)
      }
        else {
          rfsrcOutput <- c(rfsrcOutput, survOutput = list(survOutput))
        }
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
        classOutput <- vector("list", class.count)
        ## Names of the classification outputs.
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
        ## Incoming error rates: T=tree R=response L=level
        ## T1R1L0 T1R1L1 T1R1L2 T1R1L3 T1R2L0 T1R2L1 T1R2L2, T2R1L0 T2R1L1 T2R1L2 T2R1L3 T2R2L0 T2R2L1 T2R2L2, ... 
        ## In GROW mode, all class objects are represented in the tree offset calculation.
        ## In PRED mode, the offsets are dependent on the only those targets that are requested!
        ## Yields tree.offset = c(1, 8, ...) 
        tree.offset <- array(1, ntree)
        if (ntree > 1) {
          tree.offset[2:ntree] <- sum(1 + levels.count)
        }
        tree.offset <-  cumsum(tree.offset)
        ## The block offset calculation mirrors the tree offset calculation, but differs in only the primary dimension.
        block.offset <- array(1, floor(ntree/block.size))
        if (floor(ntree/block.size) > 1) {
          block.offset[2:floor(ntree/block.size)] <- sum(1 + levels.count)
        }
        block.offset <-  cumsum(block.offset)
        ## Incoming vimp rates: V=xvar R=response L=level
        ## V1R1L0 V1R1L1 V1R1L2 V1R1L3 V1R1L0 V1R2L1 V1R2L2, V2R1L0 V2R1L1 V2R1L2 V2R1L3 V2R2L0 V2R2L1 V2R2L2, ... 
        ## Yields vimp.offset = c(1, 8, ...) 
        vimp.offset <- array(1, n.xvar)
        if (n.xvar > 1) {
          vimp.offset[2:n.xvar] <- sum(1 + levels.count)
        }
        vimp.offset <-  cumsum(vimp.offset)
        iter.ensb.start <- 0
        iter.ensb.end   <- 0
        ## From the native code:
        ##   "allEnsbCLS"
        ##   "oobEnsbCLS"
        ## -> of dim [class.count] x [levels.count[]] x [n]
        ##    where this is a ragged array.
        ## From the native code:
        ##   "perfClas"
        ## -> of dim [ntree] x [class.count] x [1 + levels.count[]]
        ## where the slot [.] x [.] x [1] holds the unconditional error rate.
        ## Note that this is a ragged array.
        ## To the R code:
        ## -> of dim [[class.count]] x [ntree] x [1 + levels.count[]] 
        ## From the native code:
        ##   "blockClas"
        ## -> of dim [block.cnt] x [class.count] x [1 + levels.count[]]
        ## where the slot [.] x [.] x [1] holds the unconditional blocked error rate.
        ## Note that this is a ragged array.
        ## To the R code:
        ## -> of dim [[class.count]] x [block.cnt] x [1 + levels.count[]] 
        ## From the native code:
        ##   "vimpClas"
        ## -> of dim [n.xvar] x [class.count] x [1 + levels.count[]]
        ## where the slot [.] x [.] x [1] holds the unconditional vimp.
        ## Note that this is a ragged array.
        ## To the R code:
        ## -> of dim [[class.count]] x [1 + levels.count[]] x [n.xvar] 
        for (i in 1:class.count) {
          iter.ensb.start <- iter.ensb.end
          iter.ensb.end <- iter.ensb.end + (levels.count[i] * n)
          ens.names <- list(NULL, levels.names[[i]])
          err.names <- c("all", levels.names[[i]])
          vimp.names <- list(c("all", levels.names[[i]]), xvar.names)
          predicted <- (if (!is.null(nativeOutput$allEnsbCLS))
                        array(nativeOutput$allEnsbCLS[(iter.ensb.start + 1):iter.ensb.end],
                              c(n, levels.count[i]), dimnames=ens.names) else NULL)
          classOutput[[i]] <- list(predicted = predicted)
          response <- (if (!is.null(predicted)) bayes.rule(predicted, pi.hat) else NULL)
          classOutput[[i]] <- c(classOutput[[i]], class = list(response))
          remove(predicted)
          remove(response)
          predicted.oob <- (if (!is.null(nativeOutput$oobEnsbCLS))
                            array(nativeOutput$oobEnsbCLS[(iter.ensb.start + 1):iter.ensb.end],
                                  c(n, levels.count[i]), dimnames=ens.names) else NULL)
          classOutput[[i]] <- c(classOutput[[i]], predicted.oob = list(predicted.oob))
          response.oob <- (if (!is.null(predicted.oob)) bayes.rule(predicted.oob, pi.hat) else NULL)
          classOutput[[i]] <- c(classOutput[[i]], class.oob = list(response.oob))
          remove(predicted.oob)
          remove(response.oob)
          if (!is.null(nativeOutput$perfClas)) {
            err.rate <- array(0, c(1 + levels.count[i], ntree))
            for (j in 1: (1 + levels.count[i])) {
              err.rate[j, ]  <- nativeOutput$perfClas[tree.offset]
              tree.offset <- tree.offset + 1
            }
            row.names(err.rate) <- err.names
            classOutput[[i]] <- c(classOutput[[i]], err.rate = list(t(err.rate)))
            remove(err.rate)
          }
          if (!is.null(nativeOutput$blockClas)) {
            err.block.rate <- array(0, c(1 + levels.count[i], floor(ntree/block.size)))
            for (j in 1: (1 + levels.count[i])) {
              err.block.rate[j, ]  <- nativeOutput$blockClas[block.offset]
              block.offset <- block.offset + 1
            }
            row.names(err.block.rate) <- err.names
            classOutput[[i]] <- c(classOutput[[i]], err.block.rate = list(t(err.block.rate)))
            remove(err.block.rate)
          }
          if (!is.null(nativeOutput$vimpClas)) {
            importance <- array(0, c(1 + levels.count[i], n.xvar), dimnames=vimp.names)
            for (j in 1: (1 + levels.count[i])) {
              importance[j, ]  <- nativeOutput$vimpClas[vimp.offset]
              vimp.offset <- vimp.offset + 1
            }
            classOutput[[i]] <- c(classOutput[[i]], importance = list(t(importance)))
            remove(importance)
          }
        }
        nativeOutput$allEnsbCLS <- NULL
        nativeOutput$oobEnsbCLS <- NULL
        nativeOutput$perfClas <- NULL
        nativeOutput$vimpClas <- NULL
        nativeOutput$blockClas <- NULL
        ## When TRUE we revert to univariate nomenclature for all the outputs.
        if(univariate.nomenclature) {
          if ((class.count == 1) & (regr.count == 0)) {
            names(classOutput) <- NULL
            rfsrcOutput <- c(rfsrcOutput, unlist(classOutput, recursive=FALSE))
          }
          else {
            rfsrcOutput <- c(rfsrcOutput, classOutput = list(classOutput))
          }
        }
        else {
          rfsrcOutput <- c(rfsrcOutput, classOutput = list(classOutput))
        }
      }
      if (regr.count > 0) {
        regrOutput <- vector("list", regr.count)
        names(regrOutput) <- yvar.names[regr.index]
        ## Incoming: T=tree R=response
        ## T1R1 T1R2, T2R1 T2R2, T3R1 T3R2, ... 
        ## Yields tree.offset = c(1, 3, 5, ...) 
        tree.offset <- array(1, ntree)
        if (ntree > 1) {
          tree.offset[2:ntree] <- length(regr.index)
        }
        tree.offset <-  cumsum(tree.offset)
        ## The block offset calculation mirrors the tree offset calculation, but differs in only the primary dimension.
        block.offset <- array(1, floor(ntree/block.size))
        if (floor(ntree/block.size) > 1) {
          block.offset[2:floor(ntree/block.size)] <- length(regr.index)
        }
        block.offset <-  cumsum(block.offset)
        ## Incoming vimp rates: V=xvar R=response L=level
        ## V1R1 V1R2, V2R1 V2R2, V3R1 V3R2, ... 
        ## Yields vimp.offset = c(1, 3, 5, ...) 
        vimp.offset <- array(1, n.xvar)
        if (n.xvar > 1) {
          vimp.offset[2:n.xvar] <- length(regr.index)
        }
        vimp.offset <-  cumsum(vimp.offset)
        iter.ensb.start <- 0
        iter.ensb.end   <- 0
        ## From the native code:
        ##   "allEnsbRGR"
        ##   "oobEnsbRGR"
        ## -> of dim [regr.count] x [obsSize]
        ## From the native code:
        ##   "perfRegr"
        ## -> of dim [ntree] x [regr.count]
        ## To the R code:
        ## -> of dim [[regr.count]] x [ntree] 
        ## From the native code:
        ##   "blockRegr"
        ## -> of dim [block.cnt] x [regr.count]
        ## To the R code:
        ## -> of dim [[regr.count]] x [block.cnt]
        ## From the native code:
        ##   "vimpRegr"
        ## -> of dim [n.vxar] x [regr.count]
        ## To the R code:
        ## -> of dim  [[regr.count]] x [n.xvar]
        for (i in 1:regr.count) {
          iter.ensb.start <- iter.ensb.end
          iter.ensb.end <- iter.ensb.end + n
          vimp.names <- xvar.names
          predicted <- (if (!is.null(nativeOutput$allEnsbRGR))
                        array(nativeOutput$allEnsbRGR[(iter.ensb.start + 1):iter.ensb.end], n) else NULL)
          regrOutput[[i]] <- list(predicted = predicted)
          remove(predicted)
          predicted.oob <- (if (!is.null(nativeOutput$oobEnsbRGR))
                            array(nativeOutput$oobEnsbRGR[(iter.ensb.start + 1):iter.ensb.end], n) else NULL)
          regrOutput[[i]] <- c(regrOutput[[i]], predicted.oob = list(predicted.oob))
          remove(predicted.oob)
          if (!is.null(nativeOutput$perfRegr)) {
            err.rate <- nativeOutput$perfRegr[tree.offset]
            tree.offset <- tree.offset + 1
            regrOutput[[i]] <- c(regrOutput[[i]], err.rate = list(err.rate))
            remove(err.rate)
          }
          if (!is.null(nativeOutput$blockRegr)) {
            err.block.rate <- nativeOutput$blockRegr[block.offset]
            block.offset <- block.offset + 1
            regrOutput[[i]] <- c(regrOutput[[i]], err.block.rate = list(err.block.rate))
            remove(err.block.rate)
          }
          if (!is.null(nativeOutput$vimpRegr)) {
            importance <- nativeOutput$vimpRegr[vimp.offset]
            names(importance) <- xvar.names
            vimp.offset <- vimp.offset + 1
            regrOutput[[i]] <- c(regrOutput[[i]], importance = list(importance))
            remove(importance)
          }
        }
        nativeOutput$allEnsbRGR <- NULL
        nativeOutput$oobEnsbRGR <- NULL
        nativeOutput$perfRegr <- NULL
        nativeOutput$vimpRegr <- NULL
        nativeOutput$blockRegr <- NULL
        ## When TRUE we revert to univariate nomenclature for all the outputs.
        if(univariate.nomenclature) {
          if ((class.count == 0) & (regr.count == 1)) {
            names(regrOutput) <- NULL
            rfsrcOutput <- c(rfsrcOutput, unlist(regrOutput, recursive=FALSE))
          }
          else {
            rfsrcOutput <- c(rfsrcOutput, regrOutput = list(regrOutput))
          }
        }
        else {
          rfsrcOutput <- c(rfsrcOutput, regrOutput = list(regrOutput))
        }
      }
    }
  }## BLOCK OF CODE IS NOT EXECUTED IN IMPUTE.ONLY
  class(rfsrcOutput) <- c("rfsrc", "grow", family)
  if (big.data) {
    class(rfsrcOutput) <- c(class(rfsrcOutput), "bigdata")
  }
  return(rfsrcOutput)
}
