rfsrc <- function(formula, data, ntree = 500,
                  mtry = NULL, ytry = NULL,
                  nodesize = NULL, nodedepth = NULL,
                  splitrule = NULL, nsplit = NULL,
                  importance = c(FALSE, TRUE, "none", "anti", "permute", "random"),
                  block.size = if (any(is.element(as.character(importance), c("none", "FALSE")))) NULL else 10,
                  bootstrap = c("by.root", "none", "by.user"),
                  samptype = c("swor", "swr"),  samp = NULL, membership = FALSE,
                  sampsize = if (samptype == "swor") function(x){x * .632} else function(x){x},
                  na.action = c("na.omit", "na.impute"), nimpute = 1,
                  ntime = 150, cause,
                  perf.type = NULL,
                  proximity = FALSE, distance = FALSE, forest.wt = FALSE,
                  xvar.wt = NULL, yvar.wt = NULL, split.wt = NULL, case.wt = NULL, case.depth = FALSE, 
                  forest = TRUE,
                  save.memory = FALSE,
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
  terminal.quants <- is.hidden.terminal.quants(user.option, save.memory)
  ## set the big data flag
  big.data <- FALSE
  ## manually over-ride parameters if forest is not requested
  if (!forest) {
    terminal.qualts <- FALSE
    terminal.quants <- FALSE
    big.data <- TRUE
  }
  cse  <- is.hidden.cse(user.option)
  csv  <- is.hidden.csv(user.option)
  presort.xvar  <- is.hidden.presort.xvar(user.option)
  ## experimental option(s)
  experimental <- is.hidden.experimental(user.option)
  ## this is always false by default
  data.pass <- is.hidden.data.pass(user.option)
  ## rfq
  rfq <- is.hidden.rfq(user.option)
  ## quantile regression
  gk.quantile <- is.hidden.gk.quantile(user.option)
  quantile.regr <- is.hidden.quantile.regr(user.option)
  prob <- is.hidden.prob(user.option)
  prob.epsilon <- is.hidden.prob.epsilon(user.option)
  ## lot
  lot <- is.hidden.lot(user.option)    
  hdim <- lot$hdim
  ## mahalanobis sigma matrix
  mahalanobis.sigma <- is.hidden.mahalanobis.sigma(user.option)
  ## misc.
  base.learner <- is.hidden.base.learner(user.option)
  vtry <- is.hidden.vtry(user.option)
  holdout.array <- is.hidden.holdout.array(user.option)
  holdout.specs <- is.hidden.holdout.specs(user.option)
  empirical.risk <- is.hidden.empirical.risk(user.option)
  tdc.rule <- is.hidden.tdc.rule(user.option)
  vimp.threshold  <- is.hidden.vimp.threshold(user.option)
  ## verify key options
  ensemble <- "all"##remove option to select between "all", "oob", "inbag"
  bootstrap <- match.arg(bootstrap, c("by.root", "none", "by.user"))
  if (bootstrap == "none") {
    ensemble <- "inbag"
  }
  importance <- match.arg(as.character(importance), c(FALSE, TRUE, "none", "anti", "permute", "random"))
  na.action <- match.arg(na.action, c("na.omit", "na.impute"))
  proximity <- match.arg(as.character(proximity), c(FALSE, TRUE, "inbag", "oob", "all"))
  distance <- match.arg(as.character(distance), c(FALSE, TRUE, "inbag", "oob", "all"))
  var.used <- match.arg(as.character(var.used), c("FALSE", "all.trees", "by.tree"))
  split.depth <- match.arg(as.character(split.depth),  c("FALSE", "all.trees", "by.tree"))
  if (var.used == "FALSE") var.used <- FALSE
  if (split.depth == "FALSE") split.depth <- FALSE
  ## data cannot be missing
  if (missing(data)) stop("data is missing")
  ## data cannot have Inf or -Inf values
  ## but we are disabling this check for now
  ## if (any(unlist(mclapply(data, function(x) {any(is.infinite(x))})))) stop("data contains Inf or -Inf values")
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
  if (data.pass == TRUE) {
    ## Give the data a complete pass and assume no missing data.
    ## This is a user directive and should be used with extreme caution.
    ## In general, control of this parameter should be left to the algorithm.
    miss.flag = FALSE
  }
  else {
    if (any(is.na(data))) {
      data <- parseMissingData(formulaPrelim, data)
      ## We detect missingness.
      miss.flag <- TRUE
    }
    else {
      ## We detect no missingness.
      miss.flag <- FALSE
    }
  }
  ## Set the C-side protocol for the data pass option. If there is
  ## no missingness, allow the C-side to skip a lot of missingness
  ## checks, and hand over control of many split routines to their
  ## non-missing counterparts.
  if (miss.flag == FALSE) {
    data.pass = TRUE
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
  subj.names <- formulaDetail$subj.names
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
  ## get y-outcome type and nlevels and append to y-immutable map (saved in forest later)
  yfactor <- extract.factor(data, yvar.names)
  yfactor$types <- yvar.types <- get.yvar.type(family, yfactor$generic.types, yvar.names)
  yfactor$nlevels <- yvar.nlevels <- get.yvar.nlevels(family, yfactor$nlevels, yvar.names, data)
  ## Determine the immutable xvar factor map.
  xfactor <- extract.factor(data, xvar.names)
  xfactor$types <- xvar.types <- xfactor$generic.types
  xvar.nlevels <- xfactor$nlevels
  ## Convert the data to numeric mode, apply the na.action protocol.
  data <- finalizeData(c(subj.names, yvar.names, xvar.names), data, na.action, miss.flag)
  ## Save the row and column names for later overlay
  ## data.col.names <- colnames(data)
  data.row.names <- rownames(data)
  ## Finalize the xvar matrix.
  xvar <- as.matrix(data[, xvar.names, drop = FALSE])
  rownames(xvar) <- colnames(xvar) <- NULL
  ## get numeric levels for x-variables and append to x-immutable map (saved in forest)
  xfactor$numeric.levels <- xvar.numeric.levels <- get.numeric.levels(family, xfactor$nlevels, xvar)
  ## initialize sample size, set mtry, samptype
  n <- nrow(xvar)
  n.xvar <- length(xvar.names)
  mtry <- get.grow.mtry(mtry, n.xvar, family, splitrule)
  samptype <- match.arg(samptype, c("swor", "swr"))
  ## Get the indvidual subject identifiers if they exist - do this *after* na.protocol
  ## The default is that we do not support time dependent x-variables.
  subj.unique.count <- n
  subj <- NULL
  xvar.time <- NULL
  subj.time <- NULL
   
  ## bootstrap case
  if (bootstrap == "by.root") {
    ## now applies to Time Dependent Covariates (TDC)
    ## Nominal bootstrap options.
    if (!is.function(sampsize) && !is.numeric(sampsize)) {
      stop("sampsize must be a function or number specifying size of subsampled data")
    }
    if (is.function(sampsize)) {
      sampsize.function <- sampsize
    }
    else {
      ## convert user passed sample size to a function
      sampsize.function <- make.samplesize.function(sampsize / subj.unique.count)
    }
    sampsize <- round(sampsize.function(subj.unique.count))
    if (sampsize < 1) {
      stop("sampsize must be greater than zero")
    }
    ## for swor size is limited by the number of cases
    if (samptype == "swor" && (sampsize > subj.unique.count)) {
      sampsize.function <- function(x){x}
      sampsize <- subj.unique.count
    }
    samp <- NULL
    case.wt  <- get.weight(case.wt, subj.unique.count)
  }
  ## user specified bootstrap
  else if (bootstrap == "by.user") {    
    ## check for coherent sample: it will of be dim [.] x [ntree]
    ## where [.] is the number of rows in the data set or unique subjects in
    ## the case of time dependent covariates.
    if (is.null(samp)) {
      stop("samp must not be NULL when bootstrapping by user")
    }
    ## ntree should be coherent with the sample provided
    ntree <- ncol(samp)
    sampsize <- colSums(samp)
    if (sum(sampsize == sampsize[1]) != ntree) {
      stop("sampsize must be identical for each tree")
    }
    ## set the sample size and function
    sampsize <- sampsize[1]
    sampsize.function <- make.samplesize.function(sampsize[1] / subj.unique.count)
    ## set the case weights
    case.wt  <- get.weight(NULL, subj.unique.count)
  }
  ## This is "none".
  else {
    sampsize <- subj.unique.count
    sampsize.function <- function(x){x}
    case.wt  <- get.weight(case.wt, sampsize)
  }
  ## Set the weight matrix for xvar, split.
  ## Set the flag for forest weight.
  split.wt <- get.weight(split.wt, n.xvar)
  forest.wt <- match.arg(as.character(forest.wt), c(FALSE, TRUE, "inbag", "oob", "all"))
  ## Note that yvar.types is based on the formula:
  ## If the family is unsupervised, yvar.types is NULL
  ## If the family is regression, classification, or
  ## multivariate/mixed, it is of length greater than zero (0).  In the case
  ## of surival and competing risk, it is specifically of length two (2) for
  ## coherency, though it is ignored in the native-code.
  ## Note that ytry is set by the formula:
  ## If the family is unsupervised, ytry assumes the value specified
  ## via the formula.  If the family is regression, classification,
  ## or multivariate/mixed, it defaults to length(yvar.types).  In
  ## the case of survival, competing risk, and time dependent
  ## covariates, it is zero (0) and ignored in the native-code.
  ## With time dependent covariates, the weight y-vector is
  ## implicitly of length 3, one for "start", "stop", and "event".
  ## But "start" and "event" are ignored.  Send in a weight in the
  ## correct slot for "stop" (slot 2 based on the formula), and it
  ## will be appended, in the native code, to xvar.wt.  For example,
  ## sending in something like ywar.wt=c(0,2,0), and xvar.wt =
  ## rep(1, 17) will favour time over indigenous x-vars 2 to 1.
  if (family == "unspv") {
    ## Override any incoming weights.
    yvar.wt <- NULL
  }
  else {
    yvar.wt <- get.weight(yvar.wt, length(yvar.types))
  }
  xvar.wt  <- get.weight(xvar.wt, n.xvar)
  ## Get the y-outcomes
  yvar <- as.matrix(data[, yvar.names, drop = FALSE])
  if(dim(yvar)[2] == 0) {
    ## Override yvar if it is not present.
    yvar <- 
      yvar.nlevels  <-  yvar.numeric.levels <- yfactor <- NULL
  }
  else {
    yfactor$numeric.levels <- yvar.numeric.levels <- get.numeric.levels(family, yfactor$nlevels, yvar)
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
  ## TBD I think we need to separate event related data and dimensioning
  ## of the response which this function combines.
  ## Get event information and dimensioning for families
  event.info <- get.grow.event.info(yvar, family, ntime = ntime)
  ## Initialize nsplit, noting that it may be been overridden.
  splitinfo <- get.grow.splitinfo(formulaDetail, splitrule, hdim, nsplit, event.info)
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
    family <- get.coerced.survival.fmly(family, subj, event.info$event.type, splitinfo$name)
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
    big.data     <- FALSE
    proximity    <- FALSE
    distance     <- FALSE
    var.used     <- FALSE
    split.depth  <- FALSE
    membership   <- FALSE
    perf.type    <- "none"
    importance   <- "none"
    case.depth   <- FALSE
    terminal.qualts <- FALSE
    terminal.quants <- FALSE
    cse  <- FALSE
    csv  <- FALSE
    ## We have deleted the following line because it was incorrectly
    ## overriding the user specified option na.action in 
    ## impute.rfsrc( generic.impute.rfsrc( rfsrc() ) ).  Now, the user specified option
    ## is respected.  Note that calls from impute.rfsrc()
    ## set the impute.only bit which results in turning of all outputs
    ## other than imputed data.
    ## na.action    <- "na.impute"
  }
  ## deal with user specified holdout array
  ## confirm dimension is right
  ## set vtry to 1 to trigger holdout algorithm
  if (!is.null(holdout.array)) {
    if (nrow(holdout.array) != n.xvar | ncol(holdout.array) != ntree) {
      stop("dimension of holdout.array does not conform to p x ntree")
    }
    vtry <- 1##this only needs to be non-zero to trigger holdout
  }
  ## !!! here's where prob and prob.epsilon are set globally !!!
  ## for the user convenience if not set make coherent
  ## global assingments for  prob and prob.epsilon values
  ## takes into account splitrule and whether gk.quantile is requested
  gk.quantile <- get.gk.quantile(gk.quantile)
  prob.assign <- global.prob.assign(prob, prob.epsilon, gk.quantile, quantile.regr, splitinfo$name, n)
  prob <- prob.assign$prob
  prob.epsilon <- prob.assign$prob.epsilon
  ## mahalanobis user specified sigma matrix
  if (splitinfo$name ==  "mahalanobis" && !is.null(mahalanobis.sigma)) {
    ## check that sigma is the correct dimension
    if (length(yvar.names) != nrow(mahalanobis.sigma)
        || length(yvar.names) != ncol(mahalanobis.sigma)) {
          stop("dimension of mahalanobis sigma matrix is incorrect: must match dimension of y response")
    }
  }
  ## Bit dependencies:
  if (terminal.qualts | terminal.quants) {
    forest <- TRUE
  }
  ## get performance and rfq, gk bits
  perf.type <- get.perf(perf.type, impute.only, family)
  perf.bits <-  get.perf.bits(perf.type)
  rfq <- get.rfq(rfq)
  rfq.bits  <- get.rfq.bits(rfq, family)
  gk.quantile.bits <- get.gk.quantile.bits(gk.quantile)
  empirical.risk.bits <- get.empirical.risk.bits(empirical.risk)
  tdc.rule.bits <- get.tdc.rule.bits(tdc.rule)
  ## Assign low bits for the native code
  ensemble.bits <- get.ensemble.bits(ensemble)
  impute.only.bits <- get.impute.only.bits(impute.only, n.miss)
  var.used.bits <- get.var.used.bits(var.used)
  split.depth.bits <- get.split.depth.bits(split.depth)
  importance.value <-  get.importance(importance, perf.type)
  importance.bits <- get.importance.bits(importance, perf.type)
  bootstrap.bits <- get.bootstrap.bits(bootstrap)
  forest.bits <- get.forest.bits(forest)
  proximity.bits <- get.proximity.bits(TRUE, proximity)
  distance.bits <- get.distance.bits(TRUE, distance)
  membership.bits <-  get.membership.bits(membership)
  statistics.bits <- get.statistics.bits(statistics)
  split.cust.bits <- get.split.cust.bits(splitinfo$cust)
  case.depth.bits  <- get.case.depth.bits(case.depth)
  ## Assign high bits for the native code
  samptype.bits <- get.samptype.bits(samptype)
  forest.wt.bits <- get.forest.wt.bits(TRUE, bootstrap, forest.wt)
  na.action.bits <- get.na.action.bits(na.action)
  block.size <- get.block.size.bits(block.size, ntree)
  terminal.qualts.bits <- get.terminal.qualts.bits(terminal.qualts)
  terminal.quants.bits <- get.terminal.quants.bits(terminal.quants)
  cse.bits = get.cse.bits(cse)
  csv.bits = get.csv.bits(csv)
  presort.xvar  <- get.presort.xvar.bits(presort.xvar)
  data.pass.bits <- get.data.pass.bits(data.pass)
  experimental.bits <- get.experimental.bits(experimental)
  ## set the maximum class levels
  max.class.levels <- 0
  ## Set the trace
  do.trace <- get.trace.bits(do.trace)
  ## Start the C external timer.
  ctime.external.start  <- proc.time()
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
                                             gk.quantile.bits +
                                             statistics.bits +
                                             empirical.risk.bits +
                                             case.depth.bits),
                                  as.integer(samptype.bits +
                                             forest.wt.bits +
                                             distance.bits + 
                                             na.action.bits +
                                             split.cust.bits +
                                             membership.bits +
                                             terminal.qualts.bits +
                                             terminal.quants.bits +
                                             tdc.rule.bits +
                                             cse.bits +
                                             csv.bits +
                                             data.pass.bits),
                                  as.integer(experimental.bits),
                                  as.integer(splitinfo$index),
                                  as.integer(splitinfo$nsplit),
                                  as.integer(mtry),
                                  lot, ## object containing lot settings
                                  base.learner, ## Object containing base learner settings.  This is never NULL.
                                  as.integer(vtry),
                                  as.integer(holdout.array),
                                  holdout.specs, ## object containing speculative holdout settings
                                  as.integer(formulaDetail$ytry),
                                  as.integer(nodesize),
                                  as.integer(nodedepth),
                                  as.integer(length(cause.wt)),
                                  as.double(cause.wt),
                                  as.double(vimp.threshold),
                                  as.integer(ntree),
                                  as.integer(n),
                                  list(as.integer(length(yvar.types)),
                                       if (is.null(yvar.types)) NULL else as.character(yvar.types),
                                       if (is.null(yvar.types)) NULL else as.integer (sapply(1:length(yvar.nlevels), function(nn) {
                                         if(yvar.nlevels[nn] > 0) max(max.class.levels, yvar.nlevels[nn]) else 0})),
                                       if (is.null(yvar.numeric.levels)) NULL else sapply(1:length(yvar.numeric.levels), function(nn) {
                                         as.integer(length(yvar.numeric.levels[[nn]]))}),
                                       if (is.null(subj)) NULL else as.integer(subj),
                                       if (is.null(event.info)) NULL else as.integer(length(event.info$event.type)),
                                       if (is.null(event.info)) NULL else as.integer(event.info$event.type)),
                                  if (is.null(yvar.numeric.levels)) {
                                    NULL
                                  }
                                  else {
                                    lapply(1:length(yvar.numeric.levels),
                                           function(nn) {as.integer(yvar.numeric.levels[[nn]])})
                                  },
                                  if (is.null(yvar)) NULL else as.double(as.vector(yvar)),
                                  list(as.integer(n.xvar),
                                       as.character(xvar.types),
                                       if (is.null(xvar.types)) NULL else as.integer(xvar.nlevels),
                                       if (is.null(xvar.numeric.levels)) NULL else sapply(1:length(xvar.numeric.levels), function(nn) {as.integer(length(xvar.numeric.levels[[nn]]))}),
                                       if (is.null(xvar.time)) NULL else as.integer(xvar.time),
                                       if (is.null(subj.time)) NULL else as.integer(subj.time)),
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
                                  as.double(split.wt),
                                  as.double(yvar.wt),
                                  as.double(xvar.wt),
                                  list(if(is.null(event.info$time.interest)) as.integer(0) else as.integer(length(event.info$time.interest)),
                                       if(is.null(event.info$time.interest)) NULL else as.double(event.info$time.interest)),
                                  as.integer(nimpute),
                                  as.integer(block.size),
                                  list(if (is.null(prob)) as.integer(0) else as.integer(length(prob)),
                                       if (is.null(prob)) NULL else as.double(prob),
                                       if (is.null(prob.epsilon)) as.double(0) else as.double(prob.epsilon)),
                                  if (is.null(mahalanobis.sigma)) NULL else as.double(as.vector(mahalanobis.sigma)),
                                  as.double(presort.xvar),
                                  as.integer(get.rf.cores()))}, error = function(e) {
                                    print(e)
                                    NULL})
  ## Stop the C external timer.
  ctime.external.stop <- proc.time()
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
  class.relfrq <- NULL
  if (family == "class" && rfq) {
    class.relfrq <- table(yvar) / length(yvar)
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
    nativeArraySize = 0
    if (hdim == 0) {
      mwcpCountSummary <- rep (0, 1)
      nativeFactorArray <- vector("list", 1)
    }
     
    ## Marker for start of forest topology.  This can
    ## change with the outputs requested.  For the arithmetic
    ## related to the pivot point, refer to
    ## stackOutput.c and in particular,
    ## stackForestOutputObjects().  Do not reference internal.c or
    ## global.h for proxies of these as they may not always map to
    ## stackForestOutputObjects(). The order in the list represented by
    ## nativeOutput[[]] is TOTALLY dependent on the order in
    ## which stackAndProtect() in the C-code is executed.  So,
    ## pay attention to this, when you are not calling the list
    ## elements by their name and instead using the offset nativeOutput[[...]].
    ## The pivot is "one" based.
    pivot <- which(names(nativeOutput) == "treeID")
    ## The offset, when added to the pivot gives the start of the
    ## repeatable chunks representing the hyper-splits, 
    ## including $hcDim, and $contPTR for the case of hdim == 1.
    if (hdim == 0) {
      ## The offset is irrelevant, no calculations are conducted, so we just safe it.
      offset = 0
    }
     
    if (!is.null(base.learner)) {
      if (base.learner$synthetic.depth > 1) {
        ## Generalized MWCP counts, for when we support hdim > 1.
        mwcpCountSummarySyth <- rep (0, 1)
        nullO <- lapply(1:ntree, function(b) {
          ## Add the tree-specific number of mwcp's to the total. 
          mwcpCountSummarySyth[1] <<- mwcpCountSummarySyth[1] + nativeOutput$mwcpCTsyth[b]                    
          NULL
        })
      }
    }
    nullO <- lapply(1:ntree, function(b) {
      if (nativeOutput$leafCount[b] > 0) {
        ## The tree was not rejected.  Count the number of internal
        ## and external (terminal) nodes in the forest.
        nativeArraySize <<- nativeArraySize + (2 * nativeOutput$leafCount[b]) - 1
        mwcpCountSummary[1] <<- mwcpCountSummary[1] + nativeOutput$mwcpCT[b]
                   
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
    ## Note that we don't save the $blnodeID that is output from the C-side.
    nativeArray <- as.data.frame(cbind(nativeOutput$treeID[1:nativeArraySize],
                                       nativeOutput$nodeID[1:nativeArraySize],
                                       nativeOutput$nodeSZ[1:nativeArraySize],
                                       nativeOutput$brnodeID[1:nativeArraySize]))
    nativeArrayHeader <- c("treeID", "nodeID", "nodeSZ", "brnodeID")
    nativeArray <- as.data.frame(cbind(nativeArray,
                                       nativeOutput$parmID[1:nativeArraySize],
                                       nativeOutput$contPT[1:nativeArraySize],
                                       nativeOutput$mwcpSZ[1:nativeArraySize],
                                       nativeOutput$fsrecID[1:nativeArraySize]))
    nativeArrayHeader <- c(nativeArrayHeader, "parmID", "contPT", "mwcpSZ", "fsrecID") 
    if (mwcpCountSummary[1] > 0) {
      ## This can be NULL if there are no factor splits along this dimension.
      nativeFactorArray[[1]] <- nativeOutput$mwcpPT[1:mwcpCountSummary[1]]
    }
    nativeFactorArrayHeader <- "mwcpPT"
     
     
    ## Finally, we parse the synthetic topologies.
    nativeArraySyth <- nativeFactorArraySyth <- NULL
    nodeCountSyth <- NULL
    totalNodeCountSyth = 0
    if (!is.null(base.learner)) {
      if (base.learner$synthetic.depth > 1) {
        if (!is.null(nativeOutput$treeIDsyth)) {
          nativeArraySyth <- as.data.frame(cbind(nativeOutput$treeIDsyth,
                                                 nativeOutput$nodeIDsyth,
                                                 nativeOutput$hcDimsyth,
                                                 nativeOutput$parmIDsyth,
                                                 nativeOutput$contPTsyth,
                                                 nativeOutput$contPTRsyth,
                                                 nativeOutput$mwcpSZsyth))
          nativeArrayHeaderSyth <- c("treeID", "nodeID", "hcDim", "parmID", "contPT", "contPTR", "mwcpSZ")
          names(nativeArraySyth) = nativeArrayHeaderSyth
          totalNodeCountSyth <- length(nativeOutput$treeIDsyth)
          nodeCountSyth <- nativeOutput$nodeCountSyth
          if (mwcpCountSummarySyth[1] > 0) {
            ## This can be NULL if there are no factor splits along this dimension.
            nativeFactorArraySyth <- nativeOutput$mwcpPTsyth[1:mwcpCountSummarySyth[1]]
            nativeFactorArrayHeaderSyth <- "mwcpPT"
            names(nativeFactorArraySyth) = nativeFactorArrayHeaderSyth
          }
        }
      }
    }
     
    names(nativeArray) <- nativeArrayHeader
    names(nativeFactorArray) <- nativeFactorArrayHeader
    if (terminal.qualts | terminal.quants) {
      totalLeafCount <- sum(nativeOutput$leafCount)
      valid.mcnt.indices <- 1:totalLeafCount
      if (terminal.quants) {
        ## family specific additions to the grow object
        if (grepl("surv", family)) {
          valid.2D.surv.indices <- 1:(totalLeafCount * length(event.info$event.type) * length(event.info$time.interest))
          valid.1D.surv.indices <- 1:(totalLeafCount * length(event.info$time.interest))
          valid.mort.indices <- 1:(totalLeafCount * length(event.info$event.type))
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
            valid.clas.indices <- 1:(totalLeafCount * sum(levels.count))
          }
          if (regr.count > 0) {
            valid.regr.indices <- 1: (totalLeafCount * regr.count)
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
                              nativeOutput$ombrMembership,
                              nativeOutput$imbrMembership,
                              nativeOutput$tnRCNT[valid.mcnt.indices],
                              nativeOutput$tnACNT[valid.mcnt.indices],
                              nativeOutput$tnOCNT[valid.mcnt.indices],
                              nativeOutput$tnICNT[valid.mcnt.indices],                             
                              nativeOutput$oobSZ,
                              nativeOutput$ibgSZ);
      names(nativeArrayTNDS) <- c("tnSURV","tnMORT","tnNLSN","tnCSHZ","tnCIFN","tnREGR","tnCLAS", "tnRMBR", "tnAMBR", "tnOMBR", "tnIMBR", "tnRCNT", "tnACNT", "tnOCNT", "tnICNT", "oobSZ", "ibgSZ")
    }
    else {
      nativeArrayTNDS <- NULL
    }
    ## Node statistics are processed here.  They are not part of the forest object proper, but their close
    ## processing relationship to the forest topology makes it convenient to place that code here.  Actual
    ## inclusing in the returned model object occurs later, and not immediately below.
    if (statistics) {
      node.stats <- as.data.frame(cbind(nativeOutput$spltST[1:nativeArraySize],
                                        nativeOutput$dpthST[1:nativeArraySize]))
      names(node.stats) <- c("spltST", "dpthST")
    }
    else {
      node.stats      <- NULL
    }
    forest.out <- list(forest = TRUE,
                       nativeArray = nativeArray,
                       nativeFactorArray = nativeFactorArray,
                       leafCount = nativeOutput$leafCount,
                       totalNodeCount = dim(nativeArray)[1],
                       nativeArraySyth = nativeArraySyth,
                       nativeFactorArraySyth = nativeFactorArraySyth,
                       nodeCountSyth = nodeCountSyth,
                       nodesize = nodesize,
                       nodedepth = nodedepth,
                       ntree = ntree,
                       family = family,
                       n = n,
                       splitrule = splitinfo$name,
                       xvar = xvar,
                       xvar.names = xvar.names,
                       xvar.factor = xfactor,
                       yvar = yvar,
                       yvar.names = yvar.names,
                       yvar.factor = yfactor,
                       base.learner = base.learner,
                       block.size = block.size,
                       bootstrap = bootstrap,
                       case.wt = case.wt,
                       event.info = event.info,
                       gk.quantile = gk.quantile,
                       hdim = hdim,
                       na.action = na.action,
                       nativeArrayTNDS = nativeArrayTNDS,
                       optLoGrow = nativeOutput$optLoGrow,
                       experimental = experimental,
                       data.pass = data.pass,
                       perf.type = perf.type,
                       prob = prob,
                       prob.epsilon = prob.epsilon,
                       quantile.regr = quantile.regr,
                       rfq = rfq,
                       sampsize = sampsize.function,
                       samptype = samptype,
                       samp = samp,
                       subj = subj,
                       subj.names = subj.names,
                       seed = nativeOutput$seed,
                       seedVimp = if (importance == FALSE) NULL else nativeOutput$seedVimp,
                       terminal.qualts = terminal.qualts,
                       terminal.quants = terminal.quants,
                       importance = importance.value,
                       vimp.threshold = vimp.threshold,
                       version = "3.4.0")
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
  ## add basic information needed for downstream niceties like printing
  else {
    node.stats      <- NULL
    forest.out <- list(forest = FALSE,
                       hdim = hdim,
                       base.learner = base.learner,
                       nodesize = nodesize,
                       nodedepth = nodedepth,
                       ntree = ntree,
                       family = family,
                       n = n,
                       splitrule = splitinfo$name,
                       yvar = yvar,
                       yvar.names = yvar.names,
                       #xvar = xvar,
                       xvar.names = xvar.names,
                       subj = subj,
                       subj.names = subj.names,
                       bootstrap = bootstrap,
                       sampsize = sampsize.function,
                       samptype = samptype,
                       samp = samp,
                       case.wt = case.wt,
                       version = "3.4.0",
                       na.action = na.action,
                       perf.type = perf.type,
                       rfq = rfq,
                       gk.quantile = gk.quantile,
                       quantile.regr = quantile.regr,
                       prob = prob,
                       prob.epsilon = prob.epsilon,
                       block.size = block.size)
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
  ## process case.depth
  if (case.depth != FALSE) {
    case.depth.out <- matrix(nativeOutput$caseDepth, c(ntree, n), byrow = TRUE)
    nativeOutput$caseDepth <- NULL
  }
  else {
    case.depth.out <- NULL
  }
  ## membership
  if (membership) {
    membership.out <- matrix(nativeOutput$nodeMembership, c(n, ntree))
    inbag.out <- matrix(nativeOutput$bootMembership, c(n, ntree))
    nativeOutput$nodeMembership <- NULL
    nativeOutput$bootMembership <- NULL
    ## The result of this calculation is a list of dimension
    ## [[ntree]][[n]] in which each element [[b]]][[i]] is a
    ## vector containing the (multiple) terminal node identifiers
    ## for case [[i]] in tree [[b]].  These vectors are of varying
    ## lengths because each case [[i]] in each tree [[b]] can live
    ## in more than one node when splits on time occur.
    if(!is.null(subj)) {
      tdc.membership.cnt <- matrix(nativeOutput$nodeMembershipTDC[[1]], c(n, ntree))
      tdc.membership.out  <- vector("list", ntree)            
      begin.indx  <- 0
      end.indx  <- 0
      for (i in 1:ntree) {
        temp  <- vector("list", n)
        for (j in 1:n) {
          begin.indx  <- end.indx + 1
          end.indx  <- end.indx + tdc.membership.cnt[j,i]
          temp[[j]]  <- nativeOutput$nodeMembershipTDC[[2]][begin.indx:end.indx] 
        }
        tdc.membership.out[[i]]  <- temp
      }
    }
    else {
      tdc.membership.out  <- NULL
    }
  }
  else {
    membership.out <- NULL
    inbag.out <- NULL
    tdc.membership.out  <- NULL        
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
  empr.risk <- NULL
  oob.empr.risk <- NULL
  if (empirical.risk) {
    if (!is.null(nativeOutput$emprRisk)) {
      empr.risk <- array(nativeOutput$emprRisk, c(lot$treesize, ntree))
      nativeOutput$emprRisk <- NULL
    }
    if (!is.null(nativeOutput$oobEmprRisk)) {
      oob.empr.risk <- array(nativeOutput$oobEmprRisk, c(lot$treesize, ntree))
      nativeOutput$oobEmprRisk <- NULL
    }
  }
  if (!is.null(holdout.specs)) {
    holdout.blk <- nativeOutput$holdoutBlk
    nativeOutput$holdoutBlk <- NULL
  }
  else {
    holdout.blk = NULL
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
    xvar = if (big.data) NULL else xvar,
    xvar.names = xvar.names,
    event.info = event.info,
    subj = subj,
    subj.names = subj.names,
    xvar.wt = xvar.wt,
    split.wt = split.wt,
    cause.wt = cause.wt,
    leaf.count = nativeOutput$leafCount,
    proximity = proximity.out,
    forest = forest.out,
    forest.wt = forest.wt.out,
    case.depth = case.depth.out,
    distance = distance.out,
    membership = membership.out,
    tdc.membership = tdc.membership.out,
    splitrule = splitinfo$name,
    inbag = inbag.out,
    var.used = var.used.out,
    imputed.indv = (if (n.miss > 0) imputed.indv else NULL),
    imputed.data = (if (n.miss > 0) imputed.data else NULL),
    split.depth  = split.depth.out,
    node.stats = node.stats,
    ensemble = ensemble,
    holdout.array = holdout.array,
    block.size = block.size,
    holdout.blk = holdout.blk,
    empr.risk = empr.risk,
    oob.empr.risk = oob.empr.risk,
    ctime.internal = nativeOutput$cTimeInternal,
    ctime.external = ctime.external.stop - ctime.external.start
  )
  ## memory management
  remove(yvar)
  remove(xvar)
  remove(proximity.out)
  remove(forest.out)
  remove(forest.wt.out)
  remove(case.depth.out)
  remove(distance.out)
  remove(membership.out)
  remove(inbag.out)
  remove(var.used.out)
  if (n.miss > 0) remove(imputed.indv)
  if (n.miss > 0) remove(imputed.data)
  remove(split.depth.out)
  remove(holdout.array)
  remove(empr.risk)
  remove(oob.empr.risk)
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
      else if (family == "surv-CR") {
        ## Competing Risk names.
        ens.names <- list(NULL, NULL, c(paste("condCHF.", 1:length(event.info$event.type), sep = "")))
        mortality.names <- list(NULL, paste("event.", 1:length(event.info$event.type), sep = ""))
        cif.names <- list(NULL, NULL, c(paste("CIF.", 1:length(event.info$event.type), sep = "")))
        err.names <- list(c(paste("event.", 1:length(event.info$event.type), sep = "")), NULL)
        vimp.names <- list(paste("event.", 1:length(event.info$event.type), sep = ""), xvar.names)
      }
      else {
        ## Time dependent covariates.
        ens.names <- list(NULL, NULL)
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
      ##   "allEnsbKHZ"
      ##   "oobEnsbKHZ"
      ## -> of dim [RF_sortedTimeInterestSize] x [subject count]
      ## To the R code:
      ## -> of dim [n] x [RF_sortedTimeInterestSize]
      hazard <-  (if (!is.null(nativeOutput$allEnsbKHZ))
                    matrix(nativeOutput$allEnsbKHZ,
                           c(subj.unique.count, length(event.info$time.interest))) else NULL)
      nativeOutput$allEnsbKHZ <- NULL
      survOutput <- c(survOutput, hazard = list(hazard))
      remove(hazard)
      hazard.oob <-  (if (!is.null(nativeOutput$oobEnsbKHZ))
                        matrix(nativeOutput$oobEnsbKHZ,
                               c(subj.unique.count, length(event.info$time.interest))) else NULL)
      nativeOutput$oobEnsbKHZ <- NULL
      survOutput <- c(survOutput, hazard.oob = list(hazard.oob))
      remove(hazard.oob)
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
      ## From the native code:
      ##   "holdoutSurv"
      ##   -> of dim [p] x [holdout.blk[.]] x [RF_eventTypeSize]
      ## To the R code:
      ##   -> of dim [[p]] x [RF_eventTypeSize] x [holdout.blk[.]]
      ##              ^^^
      ##    list element can be null
      ##   if no blocks for this x-var  
      ## 2 blocks for X1                                     1 block for X2            0 blocks X3, 1 block X4
      ## --------------------------------------------------- ------------------------- -------------------------
      ## [X1]BX1[1]E1 [X1]BX1[1]E2 [X1]BX1[2]E1 [X1]BX1[2]E2 [X2]BX2[1]E1 [X2]BX2[1]E2 [X4]BX4[1]E1 [X4]BX4[1]E2  
      ## ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------
      ## 1            2            3            4            5            6            7            8
      ## holdout.blk =    c(2, 1, 0, 1, ...)
      ## offset      =    c(4, 2, 0, 2, ...)
      ## offset.sum  = c(0, 4, 6, 6, 8, ...)
      ## X1:  [[1]] is an array of dim 2 x 2
      ## X2:  [[2]] is an array of dim 2 x 1
      ## X3:  [[3]] is NULL
      ## X4:  [[4]] is an array of dim 2 x 1
      if (!is.null(nativeOutput$holdoutSurv)) {
        holdout.vimp <- vector("list", length(rfsrcOutput$holdout.blk))
        names(holdout.vimp) <- xvar.names
        ## These are lengths of of each x-var dimension.    
        holdout.offset <- rfsrcOutput$holdout.blk * length(event.info$event.type) 
        holdout.offset.sum <- c(0, cumsum(holdout.offset))
        for (i in 1:length(holdout.vimp)) {
          if (rfsrcOutput$holdout.blk[i] > 0) {
            if (length(event.info$event.type) > 1) {
              holdout.vimp[[i]] <- array(nativeOutput$holdoutSurv[(holdout.offset.sum[i] + 1) : holdout.offset.sum[i+1]], c(length(event.info$event.type), rfsrcOutput$holdout.blk[i]))
            }
            else {
              holdout.vimp[[i]] <- nativeOutput$holdoutSurv[(holdout.offset.sum[i] + 1) : holdout.offset.sum[i+1]]
            }
          }
          else {
            holdout.vimp[[i]] = NA
          }
        }
        survOutput = c(survOutput, holdout.vimp = list(holdout.vimp))
        remove(holdout.vimp)
      }
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
      resp.clas.count <- length(class.index)
      regr.index <- which(yvar.types == "R")
      resp.regr.count <- length(regr.index)
      if (resp.clas.count > 0) {
        classOutput <- vector("list", resp.clas.count)
        ## Names of the classification outputs.
        names(classOutput) <- yvar.names[class.index]
        ## Vector to hold the number of levels in each factor response. 
        levels.count <- array(0, resp.clas.count)
        ## List to hold the names of levels in each factor response. 
        levels.names <- vector("list", resp.clas.count)
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
        ## -> of dim [resp.clas.count] x [levels.count[]] x [n]
        ##    where this is a ragged array.
        ## From the native code:
        ##   "perfClas"
        ## -> of dim [ntree] x [resp.clas.count] x [1 + levels.count[]]
        ## where the slot [.] x [.] x [1] holds the unconditional error rate.
        ## Note that this is a ragged array.
        ## To the R code:
        ## -> of dim [[resp.clas.count]] x [ntree] x [1 + levels.count[]] 
        ## From the native code:
        ##   "blockClas"
        ## -> of dim [block.cnt] x [resp.clas.count] x [1 + levels.count[]]
        ## where the slot [.] x [.] x [1] holds the unconditional blocked error rate.
        ## Note that this is a ragged array.
        ## To the R code:
        ## -> of dim [[resp.clas.count]] x [block.cnt] x [1 + levels.count[]] 
        ## From the native code:
        ##   "vimpClas"
        ## -> of dim [n.xvar] x [resp.clas.count] x [1 + levels.count[]]
        ## where the slot [.] x [.] x [1] holds the unconditional vimp.
        ## Note that this is a ragged array.
        ## To the R code:
        ## -> of dim [[resp.clas.count]] x [1 + levels.count[]] x [n.xvar] 
        ## From the native code:
        ##   "cseNum"
        ## -> of dim [resp.regr.count] x [n]
        ## To the R code:
        ## -> of dim [[resp.regr.count]] x [n]
        ## From the native code:
        ##   "cseDen"
        ## -> of dim [n]
        ## To the R code:
        ## -> of dim [n]
        ## From the native code:
        ##   "csvNum"
        ## -> of dim [n.xvar] x [resp.regr.count] x [n]
        ## To the R code:
        ## -> of dim [[resp.regr.count]] x [n] x [n.xvar]
        ## From the native code:
        ##   "csvDen"
        ## -> of dim [n.xvar] x [n] 
        ## To the R code:
        ## -> of dim  x [n]  x [n.xvar]
        ## Preparing for holdout classification vimp.  See the details and examples further below for
        ## parsing the output.
        if (!is.null(nativeOutput$holdoutClas)) {
          ## These are lengths of of each x-var dimension.    
          holdout.offset.x <- rfsrcOutput$holdout.blk * (sum(1 + levels.count))
          ## This is the cumulative offset for each x-var dimension.
          holdout.offset.sum.x <- c(0, cumsum(holdout.offset.x))
          ## Relative response offset.
          holdout.offset.r <- 0
        }
        for (i in 1:resp.clas.count) {
          iter.ensb.start <- iter.ensb.end
          iter.ensb.end <- iter.ensb.end + (levels.count[i] * n)
          ens.names <- list(NULL, levels.names[[i]])
          err.names <- c("all", levels.names[[i]])
          vimp.names <- list(c("all", levels.names[[i]]), xvar.names)
          predicted <- (if (!is.null(nativeOutput$allEnsbCLS))
                          array(nativeOutput$allEnsbCLS[(iter.ensb.start + 1):iter.ensb.end],
                                c(n, levels.count[i]), dimnames=ens.names) else NULL)
          classOutput[[i]] <- list(predicted = predicted)
          response <- (if (!is.null(predicted)) get.bayes.rule(predicted, class.relfrq) else NULL)
          classOutput[[i]] <- c(classOutput[[i]], class = list(response))
          remove(predicted)
          remove(response)
          predicted.oob <- (if (!is.null(nativeOutput$oobEnsbCLS))
                              array(nativeOutput$oobEnsbCLS[(iter.ensb.start + 1):iter.ensb.end],
                                    c(n, levels.count[i]), dimnames=ens.names) else NULL)
          classOutput[[i]] <- c(classOutput[[i]], predicted.oob = list(predicted.oob))
          response.oob <- (if (!is.null(predicted.oob)) get.bayes.rule(predicted.oob, class.relfrq) else NULL)
          classOutput[[i]] <- c(classOutput[[i]], class.oob = list(response.oob))
          remove(predicted.oob)
          remove(response.oob)
          ## New case specific arrays, numerator:
          cse.num <- (if (!is.null(nativeOutput$cseClas))
                        array(nativeOutput$cseClas[(iter.ensb.start + 1):iter.ensb.end], n) else NULL)
          classOutput[[i]] <- c(classOutput[[i]], cse.num = list(cse.num))
          remove(cse.num)
          ## New case specific arrays, denominator:
          cse.den <- (if (!is.null(nativeOutput$cseDen))
                        array(nativeOutput$cseDen[1:n], n) else NULL)
          classOutput[[i]] <- c(classOutput[[i]], cse.den = list(cse.den))
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
          ## Incoming csv rates:  V=xvar R=response n=case
          ## In elemental form:
          ## V1R1n1 V1R1n2 V1R1n3, ..., V1R2n1 V1R2n2 V1R2n3, ..., V2R1n1, V2R1n2, V2R1n3, ..., V2R2n1, V2R2n2, V2R2n3, ...
          ## Equivalent in blocks on (n):
          ##       V1R1(n)  V1R2(n)  V1R3(n)  V2R1(n)  V2R2(n)  V2R3(n)  V3R1(n)  V3R2(n)  V3R3(n)
          ## To parse we use:
          ## R1:    iter11                    iter12                     iter13
          ## R2:             iter21                    iter22                     iter23
          ## R3:                      iter31                    iter32                     iter33
          ## All zeros initially and of length for this [[i]] slot in the response list.
          ## Then we make this into an [n] x [n.xvar] array for position [[i]] in the target response.
          csv.idx  <- array(0, n * n.xvar) 
          ## Jesus Christ. For an example, try n = n.xvar = 7, resp.clas.cnt = 2
          ## For response 1, i = 1:   1  2  3  4  5  6  7  8  9 10   21 22 23 24 25 26 27 28 29 30   41 42 43 44 45 46 47 48 49 50
          ## For response 2, i = 2:  11 12 13 14 15 16 17 18 19 20   31 32 33 34 35 36 37 38 39 40   51 52 53 54 55 56 57 58 59 60
          ## Then we make this into an [n] x [n.xvar] array for position [[i]] in the target response.
          for (j in 1:n.xvar) {
            csv.idx[(((j-1) * n) + 1) : (((j-1) * n) + n)]  <- ( ((i-1) * n) + ((j-1) * n * resp.clas.count) + 1 ) : ( ((i-1) * n) + ((j-1) * n * resp.clas.count) + n )
          }
          ## New case specific arrays, numerator:
          csv.num <- (if (!is.null(nativeOutput$csvClas))
                        array(nativeOutput$csvClas[csv.idx], c(n, n.xvar)) else NULL)
          classOutput[[i]] <- c(classOutput[[i]], csv.num = list(csv.num))
          remove(csv.num)
          ## New case specific arrays, denominator:
          csv.den <- (if (!is.null(nativeOutput$csvDen))
                        array(nativeOutput$csvDen, c(n, n.xvar)) else NULL)
          classOutput[[i]] <- c(classOutput[[i]], csv.den = list(csv.den))
          remove(csv.den)
          if (!is.null(nativeOutput$vimpClas)) {
            importance <- array(0, c(1 + levels.count[i], n.xvar), dimnames=vimp.names)
            for (j in 1: (1 + levels.count[i])) {
              importance[j, ]  <- nativeOutput$vimpClas[vimp.offset]
              vimp.offset <- vimp.offset + 1
            }
            classOutput[[i]] <- c(classOutput[[i]], importance = list(t(importance)))
            remove(importance)
          }
          ## From the native code:
          ##   "holdoutClas"
          ##   -> of dim [p] x [holdout.blk[.]] x [resp.clas.count] x [1 + levels.count[]]
          ## To the R code:
          ##   -> of dim [[p]] x [1 + levels.count[]] x [holdout.blk[.]]
          ##              ^^^
          ##    list element can be null
          ##   if no blocks for this x-var  
          ##
          ##   Also note that ONLY the current response [resp.clas.count] is extracted, after its offset is calculated.
          ## Assume two responses:
          ##
          ## X1: 1 of 2 blocks, 1st response (2 class), 2nd response (3 class)
          ##
          ## -------------- -------------- -------------- -------------- -------------- -------------- --------------
          ## [X1]BX1[1]R1L0 [X1]BX1[1]R1L1 [X1]BX1[1]R1L2 [X1]BX1[1]R2L0 [X1]BX1[1]R2L1 [X1]BX1[1]R2L2 [X1]BX1[1]R2L3 
          ## -------------- -------------- -------------- -------------- -------------- -------------- --------------
          ## 1              2              3              4              5              6              7            
          ## X1: 2 of 2 blocks, 1st response (2 class), 2nd response (3 class)
          ##
          ## -------------- -------------- -------------- -------------- -------------- -------------- --------------
          ## [X1]BX1[2]R1L0 [X1]BX1[2]R1L1 [X1]BX1[2]R1L2 [X1]BX1[2]R2L0 [X1]BX1[2]R2L1 [X1]BX1[2]R2L2 [X1]BX1[2]R2L3 
          ## -------------- -------------- -------------- -------------- -------------- -------------- --------------
          ## 8              9              10             11             12             13             14            
          ##
          ## X2: 1 of 1 blocks, 1st response (2 class), 2nd response (3 class)
          ##
          ## -------------- -------------- -------------- -------------- -------------- -------------- --------------
          ## [X2]BX2[1]R1L0 [X2]BX2[1]R1L1 [X2]BX2[1]R1L2 [X2]BX2[1]R2L0 [X2]BX2[1]R2L1 [X2]BX2[1]R2L2 [X2]BX2[1]R2L3 
          ## -------------- -------------- -------------- -------------- -------------- -------------- --------------
          ## 15             16             17             18             19             20             21            
          ##
          ## X3: zero blocks
          ##
          ## X4: 1 of 1 blocks, 1st response (2 class), 2nd response (3 class)
          ##
          ## -------------- -------------- -------------- -------------- -------------- -------------- --------------
          ## [X4]BX4[1]R1L0 [X4]BX4[1]R1L1 [X4]BX4[1]R1L2 [X4]BX4[1]R2L0 [X4]BX4[1]R2L1 [X4]BX4[1]R2L2 [X4]BX4[1]R2L3 
          ## -------------- -------------- -------------- -------------- -------------- -------------- --------------
          ## 22             23             24             25             26             27             28            
          ##
          ## holdout.blk =    c( 2,  1,  0,  1, ...)
          ## offset      =    c(14,  7,  0,  7, ...)
          ## offset.sum  = c(0, 14, 21, 21, 28, ...)
          ## arrays for R1:
          ## X1:  [[1]] is an array of dim [1 + levels.count[1]] x 2 = 3 x 2
          ## X2:  [[2]] is an array of dim [1 + levels.count[1]] x 1 = 3 x 1
          ## X3:  [[3]] is NULL
          ## X4:  [[4]] is an array of dim [1 + levels.count[1]] x 1 = 3 x 1
          ## arrays for R2:        
          ## X1:  [[1]] is an array of dim [1 + levels.count[2]] x 2 = 4 x 2
          ## X2:  [[2]] is an array of dim [1 + levels.count[2]] x 1 = 4 x 1
          ## X3:  [[3]] is NULL
          ## X4:  [[4]] is an array of dim [1 + levels.count[2]] x 1 = 4 x 1
          if (!is.null(nativeOutput$holdoutClas)) {    
            holdout.vimp <- vector("list", length(rfsrcOutput$holdout.blk))
            names(holdout.vimp) <- xvar.names
            ## Update the relative response offset.
            if (i > 1) {
              holdout.offset.r <- holdout.offset.r + (1 + levels.count[i-1])
            }
            ## R1:  0
            ## R2:  3
            ## Loop over each x-var, k.
            for (k in 1:length(holdout.vimp)) {
              ## Is the block count for this x-var non-zero?
              if (rfsrcOutput$holdout.blk[k] > 0) {
                ## Get the location for this x-var.
                offset.x <- holdout.offset.sum.x[k]
                ## x = 1 -> offset.x =  0
                ## x = 2 -> offset.x = 14
                ## x = 3 -> offset.x = NA
                ## x = 4 -> offset.x = 21
                ## Relative block offset.
                offset.b <- 0
                index.m <- NULL
                ## Loop over each block for this x-var
                for (m in 1:rfsrcOutput$holdout.blk[k]) {
                  if (m > 1) {
                    ## Update the relative block offset.                            
                    offset.b <- offset.b +  (sum(1 + levels.count))
                  }
                  index.m <- c(index.m, seq(from = offset.x + offset.b + holdout.offset.r + 1, by = 1, length.out = levels.count[i] + 1))
                  ## blk = 1 -> offset.b =  0
                  ## blk = 2 -> offset.b =  7
                  ## indices for R1, X1, all blocks:  01:03, 08:10
                  ## indices for R1, X2, all blocks:  15:17
                  ## indices for R1, X3, all blocks:  NULL
                  ## indices for R1, X4, all blocks:  22:24
                  ## indices for R2, X1, all blocks:  04:07, 11:14
                  ## indices for R2, X2, all blocks:  18:21
                  ## indices for R2, X3, all blocks:  NULL
                  ## indices for R2, X4, all blocks:  25:28                
                }
                ## holdout.vimp[[k]] = index.m
                ## holdout.vimp[[k]] = nativeOutput$holdoutClas[index.m]
                holdout.vimp[[k]] = array(nativeOutput$holdoutClas[index.m], c(levels.count[i] + 1, rfsrcOutput$holdout.blk[k]))
              }
              else {
                holdout.vimp[[k]] = NA
              }
            }
            classOutput[[i]] = c(classOutput[[i]], holdout.vimp = list(holdout.vimp))
            remove(holdout.vimp)
          }
        }           
        nativeOutput$allEnsbCLS  <- NULL
        nativeOutput$oobEnsbCLS  <- NULL          
        nativeOutput$perfClas    <- NULL
        nativeOutput$blockClas   <- NULL
        nativeOutput$vimpClas    <- NULL
        nativeOutput$holdoutClas <- NULL
        ## When TRUE we revert to univariate nomenclature for all the outputs.
        if(univariate.nomenclature) {
          if ((resp.clas.count == 1) & (resp.regr.count == 0)) {
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
      if (resp.regr.count > 0) {
        regrOutput <- vector("list", resp.regr.count)
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
        ## Incoming vimp rates: V=xvar R=response
        ## V1R1 V1R2, V2R1 V2R2, V3R1 V3R2, ... 
        vimp.offset <- array(1, n.xvar)
        if (n.xvar > 1) {
          vimp.offset[2:n.xvar] <- length(regr.index)
        }
        ## Yields vimp.offset = c(1, 3, 5, ...) 
        vimp.offset <-  cumsum(vimp.offset)
        iter.ensb.start <- 0
        iter.ensb.end   <- 0
        iter.qntl.start <- 0
        iter.qntl.end   <- 0
        ## From the native code:
        ##   "allEnsbRGR"
        ##   "oobEnsbRGR"
        ## -> of dim [resp.regr.count] x [n]
        ## To the R code:
        ## -> of dim [[resp.regr.count]] x [n] 
        ## From the native code:
        ##   "allEnsbQNT"
        ##   "oobEnsbQNT"
        ## -> of dim [resp.regr.count] x [length(prob)] x [n]
        ## To the R code:
        ## -> of dim [[resp.regr.count]] x [n] x [length(prob)]
        ## From the native code:
        ##   "perfRegr"
        ## -> of dim [ntree] x [resp.regr.count]
        ## To the R code:
        ## -> of dim [[resp.regr.count]] x [ntree] 
        ## From the native code:
        ##   "blockRegr"
        ## -> of dim [block.cnt] x [resp.regr.count]
        ## To the R code:
        ## -> of dim [[resp.regr.count]] x [block.cnt]
        ## From the native code:
        ##   "vimpRegr"
        ## -> of dim [n.xvar] x [resp.regr.count]
        ## To the R code:
        ## -> of dim  [[resp.regr.count]] x [n.xvar]
        ## From the native code:
        ##   "cseNum"
        ## -> of dim [resp.regr.count] x [n]
        ## To the R code:
        ## -> of dim [[resp.regr.count]] x [n]
        ## From the native code:
        ##   "cseDen"
        ## -> of dim [n]
        ## To the R code:
        ## -> of dim [n]
        ## From the native code:
        ##   "csvNum"
        ## -> of dim [n.xvar] x [resp.regr.count] x [n]
        ## To the R code:
        ## -> of dim [[resp.regr.count]] x [n] x [n.xvar]
        ## From the native code:
        ##   "csvDen"
        ## -> of dim [n.xvar] x [n] 
        ## To the R code:
        ## -> of dim  x [n]  x [n.xvar]
        ## Preparing for holdout regression vimp.  See the details and examples further below for
        ## parsing the output.
        if (!is.null(nativeOutput$holdoutRegr)) {
          ## These are lengths of of each x-var dimension.    
          holdout.offset.x <- rfsrcOutput$holdout.blk * resp.regr.count
          ## This is the cumulative offset for each x-var dimension.
          holdout.offset.sum.x <- c(0, cumsum(holdout.offset.x))
          ## Relative response offset.
          holdout.offset.r <- 0
        }
        for (i in 1:resp.regr.count) {
          iter.ensb.start <- iter.ensb.end
          iter.ensb.end <- iter.ensb.end + n
          iter.qntl.start <- iter.qntl.end
          iter.qntl.end <- iter.qntl.end + (length(prob) * n)
          vimp.names <- xvar.names
          predicted <- (if (!is.null(nativeOutput$allEnsbRGR))
                          array(nativeOutput$allEnsbRGR[(iter.ensb.start + 1):iter.ensb.end], n) else NULL)
          regrOutput[[i]] <- list(predicted = predicted)
          remove(predicted)
          predicted.oob <- (if (!is.null(nativeOutput$oobEnsbRGR))
                              array(nativeOutput$oobEnsbRGR[(iter.ensb.start + 1):iter.ensb.end], n) else NULL)
          regrOutput[[i]] <- c(regrOutput[[i]], predicted.oob = list(predicted.oob))
          remove(predicted.oob)
          ## New case specific arrays, numerator:
          cse.num <- (if (!is.null(nativeOutput$cseRegr))
                        array(nativeOutput$cseRegr[(iter.ensb.start + 1):iter.ensb.end], n) else NULL)
          regrOutput[[i]] <- c(regrOutput[[i]], cse.num = list(cse.num))
          remove(cse.num)
          ## New case specific arrays, denominator:
          cse.den <- (if (!is.null(nativeOutput$cseDen))
                        array(nativeOutput$cseDen[1:n], n) else NULL)
          regrOutput[[i]] <- c(regrOutput[[i]], cse.den = list(cse.den))
          quantile <- (if (!is.null(nativeOutput$allEnsbQNT))
                         array(nativeOutput$allEnsbQNT[(iter.qntl.start + 1):iter.qntl.end],
                               c(n, length(prob))) else NULL)
          regrOutput[[i]] <- c(regrOutput[[i]], quantile = list(quantile))
          remove(quantile)
          quantile.oob <- (if (!is.null(nativeOutput$oobEnsbQNT))
                             array(nativeOutput$oobEnsbQNT[(iter.qntl.start + 1):iter.qntl.end],
                                   c(n, length(prob))) else NULL)
          regrOutput[[i]] <- c(regrOutput[[i]], quantile.oob = list(quantile.oob))
          remove(quantile.oob)
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
          ## Incoming csv rates:  V=xvar R=response n=case
          ## In elemental form:
          ## V1R1n1 V1R1n2 V1R1n3, ..., V1R2n1 V1R2n2 V1R2n3, ..., V2R1n1, V2R1n2, V2R1n3, ..., V2R2n1, V2R2n2, V2R2n3, ...
          ## Equivalent in blocks on (n):
          ##       V1R1(n)  V1R2(n)  V1R3(n)  V2R1(n)  V2R2(n)  V2R3(n)  V3R1(n)  V3R2(n)  V3R3(n)
          ## To parse we use:
          ## R1:    iter11                    iter12                     iter13
          ## R2:             iter21                    iter22                     iter23
          ## R3:                      iter31                    iter32                     iter33
          ## All zeros initially and of length for this [[i]] slot in the response list.
          ## Then we make this into an [n] x [n.xvar] array for position [[i]] in the target response.
          csv.idx  <- array(0, n * n.xvar) 
          ## Jesus Christ. For an example, try n = n.xvar = 7, resp.regr.cnt = 2
          ## For response 1, i = 1:   1  2  3  4  5  6  7  8  9 10   21 22 23 24 25 26 27 28 29 30   41 42 43 44 45 46 47 48 49 50
          ## For response 2, i = 2:  11 12 13 14 15 16 17 18 19 20   31 32 33 34 35 36 37 38 39 40   51 52 53 54 55 56 57 58 59 60
          ## Then we make this into an [n] x [n.xvar] array for position [[i]] in the target response.
          for (j in 1:n.xvar) {
            csv.idx[(((j-1) * n) + 1) : (((j-1) * n) + n)]  <- ( ((i-1) * n) + ((j-1) * n * resp.regr.count) + 1 ) : ( ((i-1) * n) + ((j-1) * n * resp.regr.count) + n )
          }
          ## New case specific arrays, numerator:
          csv.num <- (if (!is.null(nativeOutput$csvRegr))
                        array(nativeOutput$csvRegr[csv.idx], c(n, n.xvar)) else NULL)
          regrOutput[[i]] <- c(regrOutput[[i]], csv.num = list(csv.num))
          remove(csv.num)
          ## New case specific arrays, denominator:
          csv.den <- (if (!is.null(nativeOutput$csvDen))
                        array(nativeOutput$csvDen, c(n, n.xvar)) else NULL)
          regrOutput[[i]] <- c(regrOutput[[i]], csv.den = list(csv.den))
          remove(csv.den)
          ## From the native code:
          ##   "holdoutRegr"
          ##   -> of dim [p] x [holdout.blk[.]] x [resp.regr.count]
          ## To the R code:
          ##   -> of dim [[p]] x [holdout.blk[.]]
          ##              ^^^
          ##    list element can be null
          ##   if no blocks for this x-var  
          ##
          ##   Also note that ONLY the current response [resp.regr.count] is extracted, after its offset is calculated.
          ## Assume three responses:
          ##
          ## X1: 1 of 2 blocks, 3 responses
          ##
          ## -------------- -------------- -------------- 
          ## [X1]BX1[1]R1   [X1]BX1[1]R2   [X1]BX1[1]R3   
          ## -------------- -------------- -------------- 
          ## 1              2              3              
          ## X1: 2 of 2 blocks, 1st response (2 class), 2nd response (3 class)
          ##
          ## -------------- -------------- -------------- 
          ## [X1]BX1[2]R1   [X1]BX1[2]R2   [X1]BX1[2]R3   
          ## -------------- -------------- -------------- 
          ## 4              5              6             
          ##
          ## X2: 1 of 1 blocks, 3 responses
          ##
          ## -------------- -------------- -------------- 
          ## [X2]BX2[1]R1   [X2]BX2[1]R2   [X2]BX2[1]R3   
          ## -------------- -------------- -------------- 
          ## 7              8              9             
          ##
          ## X3: zero blocks
          ##
          ## X4: 1 of 1 blocks, 3 responses
          ##
          ## -------------- -------------- --------------
          ## [X4]BX4[1]R1   [X4]BX4[1]R2   [X4]BX4[1]R3  
          ## -------------- -------------- --------------
          ## 10             11             12            
          ##
          ## holdout.blk =    c( 2,  1,  0,  1, ...)
          ## offset      =    c(6, 3,  0,  3, ...)
          ## offset.sum  = c(0, 6, 9,  9, 12, ...)
          ## arrays for R1:
          ## X1:  [[1]] is an array of len 2
          ## X2:  [[2]] is an array of len 1
          ## X3:  [[3]] is NULL
          ## X4:  [[4]] is an array of len 1
          ## arrays for R2:        
          ## X1:  [[1]] is an array of len 2
          ## X2:  [[2]] is an array of len 1
          ## X3:  [[3]] is NULL
          ## X4:  [[4]] is an array of len 1
          if (!is.null(nativeOutput$holdoutRegr)) {
            holdout.vimp <- vector("list", length(rfsrcOutput$holdout.blk))
            names(holdout.vimp) <- xvar.names
            ## Update the relative response offset.
            if (i > 1) {
              holdout.offset.r <- holdout.offset.r + 1
            }
            ## R1:  0
            ## R2:  1
            ## R3:  2
            ## Loop over each x-var, k.
            for (k in 1:length(holdout.vimp)) {
              ## Is the block count for this x-var non-zero?
              if (rfsrcOutput$holdout.blk[k] > 0) {
                ## Get the location for this x-var.
                offset.x <- holdout.offset.sum.x[k]
                ## x = 1 -> offset.x =  0
                ## x = 2 -> offset.x =  6
                ## x = 3 -> offset.x = NA
                ## x = 4 -> offset.x =  9
                ## Relative block offset.
                offset.b <- 0
                index.m <- NULL
                ## Loop over each block for this x-var
                for (m in 1:rfsrcOutput$holdout.blk[k]) {
                  if (m > 1) {
                    ## Update the relative block offset.                            
                    offset.b <- offset.b + resp.regr.count
                  }
                  index.m <- c(index.m, offset.x + offset.b + holdout.offset.r + 1)
                  ## blk = 1 -> offset.b =  0
                  ## blk = 2 -> offset.b =  3
                  ## indices for R1, X1, all blocks:  01,03
                  ## indices for R1, X2, all blocks:  07
                  ## indices for R1, X3, all blocks:  NULL
                  ## indices for R1, X4, all blocks:  10
                  ## indices for R2, X1, all blocks:  02,05
                  ## indices for R2, X2, all blocks:  08
                  ## indices for R2, X3, all blocks:  NULL
                  ## indices for R2, X4, all blocks:  11                
                }
                ## holdout.vimp[[k]] = index.m
                ## holdout.vimp[[k]] = nativeOutput$holdoutRegr[index.m]
                holdout.vimp[[k]] <- nativeOutput$holdoutRegr[index.m]
              }
              else {
                holdout.vimp[[k]] <- NA
              }
            }
            regrOutput[[i]] <- c(regrOutput[[i]], holdout.vimp = list(holdout.vimp))        
          }
        }
        nativeOutput$allEnsbRGR  <- NULL
        nativeOutput$oobEnsbRGR  <- NULL
        nativeOutput$allEnsbQNT  <- NULL          
        nativeOutput$oobEnsbQNT  <- NULL          
        nativeOutput$perfRegr    <- NULL
        nativeOutput$blockRegr   <- NULL
        nativeOutput$vimpRegr    <- NULL
        nativeOutput$holdoutRegr <- NULL
        ## When TRUE we revert to univariate nomenclature for all the outputs.
        if(univariate.nomenclature) {
          if ((resp.clas.count == 0) & (resp.regr.count == 1)) {
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
  #if (big.data) {
  #  class(rfsrcOutput) <- c(class(rfsrcOutput), "bigdata")
  #}
  return(rfsrcOutput)
}
