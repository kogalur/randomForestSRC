generic.predict.rfsrc <-
  function(object,
           newdata,
           m.target = NULL,
           importance = FALSE,
           get.tree = NULL,
           block.size = NULL,
           importance.xvar,
           na.action = c("na.omit", "na.impute", "na.random"),
           outcome = c("train", "test"),
           perf.type = NULL,
           proximity = FALSE,
           forest.wt = FALSE,
           ptn.count = 0,
           distance = FALSE,
           var.used = c(FALSE, "all.trees", "by.tree"),
           split.depth = c(FALSE, "all.trees", "by.tree"),
           case.depth = FALSE,
           seed = NULL,
           do.trace = FALSE,
           membership = FALSE,
           marginal.xvar = NULL,
           ...)
{
  univariate.nomenclature <- TRUE
  ## get any hidden options
  user.option <- list(...)
  ## incoming parameter checks: all are fatal
  if (missing(object)) {
    stop("object is missing!")
  }
  ## incoming object must be a grow forest or a forest object
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2)
    stop("this function only works for objects of class `(rfsrc, grow)' or '(rfsrc, forest)'")
  ## grow forests must have true forest information
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.forest.missing(object)) {
      stop("Forest information for prediction is missing.  Re-run rfsrc (grow call) with forest=TRUE")
    }
  }
  ## case specific performance
  cse  <- is.hidden.cse(user.option)
  csv  <- is.hidden.csv(user.option)
  ## verify importance (needed for jitt)
  importance <- match.arg(as.character(importance)[1], c(FALSE, TRUE,
          "none", "anti", "permute", "random", "anti.joint", "permute.joint", "random.joint"))
  if (grepl("joint", importance)) {
    vimp.joint <- TRUE
  }
  else {
    vimp.joint <- FALSE
  }
  ## verify na.action (needed for jitt)
  na.action <- match.arg(na.action, c("na.omit", "na.impute", "na.random"))
  if (inherits(object, "anonymous") && na.action == "na.impute") {
      na.action <- "na.mean"
  }
  ## verify other key options
  outcome <- match.arg(outcome, c("train", "test"))
  proximity <- match.arg(as.character(proximity), c(FALSE, TRUE, "inbag", "oob", "all"))
  forest.wt <- match.arg(as.character(forest.wt), c(FALSE, TRUE, "inbag", "oob", "all"))  
  distance <- match.arg(as.character(distance), c(FALSE, TRUE, "inbag", "oob", "all"))
  ## set restore.mode and the ensemble option
  ## newdata missing --> restore.mode = TRUE
  ## outcome = "test" --> restore.mode = FALSE for R-processing but switched later for C-code function call
  if (missing(newdata)) {
    restore.mode <- TRUE
    outcome <- "train"
    ensemble <- "all"##remove option to select between "all", "oob", "inbag"
  }
  else {##there is test data present
    restore.mode <- FALSE
    ## special treatment for outcome=="test" (which is really restore mode)
    if (outcome == "test") {
      ensemble <- "all"##remove option to select between "all", "oob", "inbag"
    }
    ## standard prediction scenario on new test data - there is no OOB
    else {
      ensemble <- "inbag"
    }
  }
  ## check if this is an anonymous object
  ## coerce values as necessary
  ## graceful return if restore.mode = TRUE which is not allowed for anonymous
  if (inherits(object, "anonymous")) {
    anonymize.bits <- 2^26
    outcome <- "train"
    if (restore.mode) {
      stop("in order to predict with anonymous forests please provide a test data set")
    }
  }
  else {
    anonymize.bits <- 0
  }
  ## jitt
  jitt <- is.hidden.jitt(user.option, importance, na.action, inherits(object, "anonymous"), partial = FALSE)
  ## rfq
  rfq <- is.hidden.rfq(user.option)
  ## quantile regression
  gk.quantile <- is.hidden.gk.quantile(user.option)
  prob <- is.hidden.prob(user.option)
  prob.epsilon <- is.hidden.prob.epsilon(user.option)
  ## vimp
  vimp.threshold  <- is.hidden.vimp.threshold(user.option)
  ## set the family
  family <- object$family
  ## pull the x-variable and y-outcome names from the grow object
  xvar.names <- object$xvar.names
  yvar.names <- object$yvar.names
  subj.names <- object$subj.names
  importance.xvar <- get.importance.xvar(importance.xvar, importance, object)
  importance.xvar.idx <- match(importance.xvar, xvar.names)
  marginal.xvar <- get.marginal.xvar(marginal.xvar, object)
  marginal.xvar.idx <- NULL
  if (length(marginal.xvar) > 0) {
    marginal.xvar.idx <- match(marginal.xvar, xvar.names)
  }
  var.used <- match.arg(as.character(var.used), c("FALSE", "all.trees", "by.tree"))
  if (var.used == "FALSE") var.used <- FALSE
  split.depth <- match.arg(as.character(split.depth),  c("FALSE", "all.trees", "by.tree"))
  if (split.depth == "FALSE") split.depth <- FALSE
  ## Currently we do not support split.depth for test data. TBD2
  split.depth = FALSE
  ## initialize the seed
  seed <- get.seed(seed)
  ## REDUCES THE OBJECT TO THE FOREST -- REDUCTION STARTS HERE
  ## hereafter we only need the forest and reassign "object" to the forest
  ## (TBD3) memory management "big.data" not currently implemented:
  big.data <- FALSE
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
    if (inherits(object, "bigdata")) {
      big.data <- TRUE
    }
    object <- object$forest
  }
  else {
    ## object is already a forest
    if (inherits(object, "bigdata")) {
      big.data <- TRUE
    }
  }
  ## confirm version coherence
  if (is.null(object$version)) {
    cat("\n  This function only works with objects created with the following minimum version of the package:")
    cat("\n    Minimum version:  ")
    cat("2.3.0")
    cat("\n    Your version:     ")
    cat("unknown")
    cat("\n")
    stop()
  }
    else {
      object.version <- as.integer(unlist(strsplit(object$version, "[.]")))
      installed.version <- as.integer(unlist(strsplit("3.4.4", "[.]")))
      minimum.version <- as.integer(unlist(strsplit("2.3.0", "[.]")))
      object.version.adj <- object.version[1] + (object.version[2]/10) + (object.version[3]/100)
      installed.version.adj <- installed.version[1] + (installed.version[2]/10) + (installed.version[3]/100)
      minimum.version.adj <- minimum.version[1] + (minimum.version[2]/10) + (minimum.version[3]/100)
      ## Minimum object version must be satisfied for us to proceed.  This is the only way
      ## terminal node restoration is guaranteed, due to RNG coherency.
      if (object.version.adj >= minimum.version.adj) {
        ## We are okay
      }
        else {
          cat("\n  This function only works with objects created with the following minimum version of the package:")
          cat("\n    Minimum version:  ")
          cat("2.3.0")
          cat("\n    Your version:     ")
          cat(object$version)
          cat("\n")
          stop()
        }
  }
  ## classification specific details related to rfq and perf.type
  class.relfrq <- NULL
  if (family == "class") {
    ## rfq specific details
    if (!is.null(rfq)) {##predict has specified rfq
      if (!rfq) {##predict does not want rfq
        ## nothing 
      }
      else {##predict has requested rfq
        class.relfrq <- prop.table(table(object$yvar))
      }
    }
    if (is.null(rfq)) {##predict  ambivalent about rfq
      if (!object$rfq) {##grow did not use rfq
        ## nothing -> rfq = FALSE
      }
      else {##grow used rfq - use grow spec
        class.relfrq <- prop.table(table(object$yvar))
        rfq <- TRUE
      }
    }
    ## performance details
    if (is.null(perf.type) && !is.null(object$perf.type)) {
      perf.type <- object$perf.type
    }
  }
  ## recover the split rule
  splitrule <- object$splitrule
  ## gk processing
  if (!is.null(gk.quantile) || object$gk.quantile) { 
    if (is.null(gk.quantile)) {##predict ambivalent about gk - use grow spec
      gk.quantile <- object$gk.quantile
    }
  }
  ## !!! here's where prob and prob.epsilon are set globally !!!
  gk.quantile <- get.gk.quantile(gk.quantile)
  prob.assign <- global.prob.assign(if (is.null(prob)) object$prob else prob,
                                    if (is.null(prob.epsilon)) object$prob.epsilon else prob.epsilon,
                                    gk.quantile, object$quantile.regr, splitrule, object$n)
  ## Determine the immutable yvar factor map which is needed for
  ## classification sexp dimensioning.  But, first convert object$yvar
  ## to a data frame which is required for factor processing
  #object$yvar <- as.data.frame(object$yvar)
  #colnames(object$yvar) <- yvar.names
  yfactor <- object$yvar.factor
  ## multivariate family details
  m.target.idx <- get.outcome.target(family, yvar.names, m.target)
  ## short cut to get the y-outcome type and number of levels
  yvar.types <- yfactor$types
  yvar.nlevels  <- yfactor$nlevels
  yvar.numeric.levels  <- yfactor$numeric.levels
  ## recover the individual subject identifiers, if they exist.
  subj <- object$subj
  ## get event information for survival families
  event.info <- object$event.info
  event.type <- event.info$event.type
  ## CR.bits assignment
  cr.bits <- get.cr.bits(family)
  ## determine the immutable xvar factor map
  xfactor <- object$xvar.factor
  any.xvar.factor <-  (length(xfactor$factor) + length(xfactor$order)) > 0
  ## short cut to get the x-variable type and number of levels
  xvar.types <- xfactor$types
  xvar.nlevels <- xfactor$nlevels
  xvar.numeric.levels <- xfactor$numeric.levels
  ## set dimensions
  n.xvar <- length(xvar.names)
  n <- object$n
  r.dim <- event.info$rdim
  ## set flags appropriately for unsupervised forests
  ## there are layers of checks appearing later, so some of these are redundant
  if (family == "unsupv") {
    outcome <- "train"
    perf.type <- "none"
    importance <- "none"
  }
  ## Override ptn.count by family.  Pruning is based on variance,
  ## thus, this is not yet well-defined in [S] settings.
  if (grepl("surv", family)) {
    ptn.count <- 0
  }
  ## ----------------------------------------------------------------
  ## From the native code's perspective, PRED mode can process one
  ## or two data sets.  If one data set is sent in, we assume
  ## we wish to restore the forest with original-training or
  ## pseudo-training data.  If two data sets are sent in, we
  ## assume we are sending in the original-training data set, and a
  ## test data set.  When outcome="test" we make a call to the native
  ## code with only one data set (the test data which becomes pseudo-training
  ## data).  We set the test outcome bit to allow the native
  ## code to distinguish this case from the case of the restore
  ## with the original training data.
  ## ----------------------------------------------------------------
  ##--------------------------------------------------------
  ##
  ## NON-RESTORE MODE PROCESSING (includes outcome=="test")
  ##
  ##--------------------------------------------------------
  if (!restore.mode) {
    ## Filter the test data based on the formula
    newdata <- newdata[, is.element(names(newdata), c(yvar.names, xvar.names)), drop = FALSE]
    ## impute.mean for anonymous forests
    if (inherits(object, "anonymous") && na.action == "na.mean") {
      newdata <- assign.impute.mean(newdata, object$impute.mean)
      na.action <- "na.omit"
    }
    ## Check that test/train factors are the same.  If factor has an
    ## NA in its levels, remove it.  Confirm factor labels overlap.
    newdata <- rm.na.levels(newdata, xvar.names)
    #newdata.xfactor <- extract.factor(newdata, xvar.names)
    #if (!setequal(xfactor$factor, newdata.xfactor$factor)) {
    #  stop("x-variable factors from test data do not match original training data")
    #}
    #if (!setequal(xfactor$order, newdata.xfactor$order)) {
    #  stop("(ordered) x-variable factors from test data do not match original training data")
    #}
    any.outcome.factor <- family == "class"
    if (family == "class+" | family ==  "mix+") {
      if (length(intersect("R", yfactor$generic.types[m.target.idx])) == 0) {
        any.outcome.factor <- TRUE
      }
    }
    ## If the outcomes contain factors we need to check that train/test y-outcomes are compatible
    if (any.outcome.factor) {
      if (sum(is.element(names(newdata), yvar.names)) > 0) {
        newdata <- rm.na.levels(newdata, yvar.names)
        #newdata.yfactor <- extract.factor(newdata, yvar.names)
        #if (!setequal(yfactor$factor, newdata.yfactor$factor)) {
        #  stop("class outcome from test data does not match original training data")
        #}
        #if (!setequal(yfactor$order, newdata.yfactor$order)) {
        #  stop("(ordered) class outcome from test data does not match original training data")
        #}
      }
    }
    ## one final (possibly redundant) check confirming coherence of train/test xvars
    ## previously commented out: but put we this back because error message is helpful
    if (length(xvar.names) != sum(is.element(xvar.names, names(newdata)))) {
      stop("x-variables in test data do not match original training data")
    }
    ## coherence of train/test yvars (assuming test yvars are available)
    yvar.present <- sum(is.element(yvar.names, names(newdata))) > 0
    if (yvar.present && length(yvar.names) != sum(is.element(yvar.names, names(newdata)))) {
      stop("y-variables in test data do not match original training data")
    }
    ## Force test factor levels to equal grow factor levels
    ## this is crucial to ensuring an immutable map
    if (any.xvar.factor) {
      newdata <- check.factor(newdata, xfactor)
    }
    ## If the outcomes contain factors we need to check factor coherence.
    if (any.outcome.factor) {
      if (yvar.present) {
        ### COMMENT THIS OUT IF YOU WANT TO ELIMINATE BAD Y-TEST OUTCOME LABELS
        ### newdata <- check.factor(newdata, yfactor, ignore = FALSE)
        newdata <- check.factor(newdata, yfactor, ignore = TRUE)
      }
    }
    ## Extract test yvar names (if any) and xvar names.
    if (yvar.present) {
      fnames <- c(yvar.names, xvar.names)
    }
      else {
        fnames <- xvar.names
      }
    ## Data conversion to numeric mode
    newdata <- finalizeData(fnames, newdata, na.action)
    ## Extract the test x-matrix and sort the columns as in the original training data.
    ## this accomodates incoming test x-matrix in a different order.
    xvar.newdata  <- as.matrix(newdata[, xvar.names, drop = FALSE])
    n.newdata <- nrow(newdata)
    ## Save the row names for later overlay
    newdata.row.names <- rownames(xvar.newdata)
    ## Process the y-outcomes and set their dimension
    ## note that in unsupervised mode there are no responses
    ## r.dim.newdata is set correctly by get.grow.event.info() in this case
    if (yvar.present) {
      yvar.newdata <- as.matrix(newdata[, yvar.names, drop = FALSE])
      event.info.newdata <- get.grow.event.info(yvar.newdata, family, need.deaths = FALSE)
      r.dim.newdata <- event.info.newdata$r.dim
      ## Survival specific coherency checks
      ## if there are no deaths, turn off performance values and VIMP
      if (grepl("surv", family) && all(na.omit(event.info.newdata$cens) == 0)) {
        perf.type <- "none"
        importance <- "none"
      }
      ## Ensure consistency of event types
      if (grepl("surv", family) &&
          length(setdiff(na.omit(event.info.newdata$cens), na.omit(event.info$cens))) > 1) {
        stop("survival events in test data do not match training data")
      }
    }
      else {
        ## Disallow outcome=TEST without y-outcomes
        if (outcome == "test") {
          stop("outcome=TEST, but the test data has no y values, which is not permitted")
        }
        ## There are no outcomes.
        r.dim.newdata <- 0
        yvar.newdata <-  NULL
        perf.type <- "none"
        importance <- "none"
      }
    ## Remove xvar row and column names for proper processing by the native library
    ## does not apply when outcome = TEST because the xvar TEST data has been made NULL
    if (outcome != "test") {
      rownames(xvar.newdata) <- colnames(xvar.newdata) <- NULL
    }
    ## We don't need the test data anymore
    remove(newdata)
  }
  ##--------------------------------------------------------
  ##
  ## RESTORE MODE PROCESSING (excludes outcome=="test")
  ##
  ##--------------------------------------------------------
  else {
    ## There cannot be test data in restore mode
    ## The native code switches based on n.newdata being zero (0).  Be careful.
    n.newdata <- 0
    r.dim.newdata <- 0
    xvar.newdata <- NULL
    yvar.newdata <-  NULL
    ## Outcome is set to train for the native code
    ## Determine whether performance values are requested
    outcome <- "train"
    if (object$bootstrap == "none" || object$bootstrap == "by.node" || family == "unsupv") {
      importance <- "none"
      perf.type <- "none"
    }
    else {
      ## do nothing.
    }
    ## restore hidden parameters
    if (is.null(user.option$vimp.threshold)) {
      vimp.threshold <- object$vimp.threshold 
    }      
  } ## ends restore.mode check
  ## ------------------------------------------------------------
  ## We have completed the restore/non-restore mode processing
  ## ------------------------------------------------------------
  ## Final processing of xvar and yvar test data
  ## depends on "outcome"
  ## outcome=train
  if (outcome == "train") {
    ## data conversion for y-training data
    yvar <- as.matrix(data.matrix(data.frame(object$yvar)))
    ## respect the training options related to bootstrapping:
    sampsize <- round(object$sampsize(n))
    case.wt <- object$case.wt
    samp <- object$samp
    ## data conversion for x-training data
    if (anonymize.bits == 0) {
      xvar <- as.matrix(data.matrix(object$xvar))
      rownames(xvar) <- colnames(xvar) <- NULL
    }
    else {
      xvar <- NULL
    }
  }
  if (outcome == "test") {
    ## outcome=test
    ## From the native code perspective we are in pseudo-restore mode
    ## From the R side, it is convenient to now pretend we are
    ## in restore mode, so we swap the training data out with the test data
    ## Performance is always requested for this setting
    ## swap the data
    xvar <- xvar.newdata
    yvar <- yvar.newdata
    rownames(xvar) <- colnames(xvar) <- NULL
    restore.mode <- TRUE
    ## set the dimensions
    n <- nrow(xvar)
    r.dim <- r.dim.newdata
    ## pretend there is no test data, but do *not* get rid of it
    ## we need (and use) this data *after* the native code call
    n.newdata <- 0
    r.dim.newdata <- 0
    ## set the sample size
    ## "swor" now handled because we now make sampsize a function of n
    sampsize <- round(object$sampsize(n))
    case.wt <- get.weight(NULL, n)
    samp <- NULL
  }
  ## initialize the number of trees in the forest
  ntree <- object$ntree
  ## process the get.tree vector that specifies which trees we want
  ## to extract from the forest.  This is only relevant to restore mode.
  ## The parameter is ignored in predict mode.
  get.tree <- get.tree.index(get.tree, ntree)
  ## get performance and rfq, gk bits
  perf.type <- get.perf(perf.type, FALSE, family)
  perf.bits <-  get.perf.bits(perf.type)
  rfq <- get.rfq(rfq)
  rfq.bits <- get.rfq.bits(rfq, family)
  gk.quantile.bits <- get.gk.quantile.bits(gk.quantile)
  bootstrap.bits <- get.bootstrap.bits(object$bootstrap)
  ## initialize the low bits
  ensemble.bits <- get.ensemble.bits(ensemble)
  importance.bits <- get.importance.bits(importance, perf.type)
  proximity.bits <- get.proximity.bits(restore.mode, proximity)
  distance.bits <- get.distance.bits(restore.mode, distance)
  split.depth.bits <- get.split.depth.bits(split.depth)
  var.used.bits <- get.var.used.bits(var.used)
  outcome.bits <- get.outcome.bits(outcome)
  case.depth.bits  <- get.case.depth.bits(case.depth)
  ## Initalize the high bits
  samptype.bits <- get.samptype.bits(object$samptype)
  ## forest weights
  forest.wt.bits <- get.forest.wt.bits(restore.mode, object$bootstrap, forest.wt)
  membership.bits <-  get.membership.bits(membership)
  terminal.qualts.bits <- get.terminal.qualts.predict.bits(object$terminal.qualts)
  terminal.quants.bits <- get.terminal.quants.predict.bits(object$terminal.quants)
  cse.bits = get.cse.bits(cse)
  csv.bits = get.csv.bits(csv)
  jitt.bits <- get.jitt.bits(jitt)
  ## set the data.pass flags: we do this here because the restore.mode flag is now finalized
  ## training data.pass acquires the grow data.pass flag 
  data.pass.bits <- get.data.pass.bits(object$data.pass)
  ## testing data.pass is na.action AND restore.mode dependent
  if (restore.mode == FALSE) {
    if (na.action == "na.omit") {
      data.pass.predict.bits  <- get.data.pass.predict.bits(TRUE)
    }
    else {
      data.pass.predict.bits  <- get.data.pass.predict.bits(FALSE)
    }
  }
  else {
    ## we are in restore mode -> we are safe and the test data.pass flag irrelevant
    data.pass.predict.bits  <- get.data.pass.predict.bits(FALSE)
  }
  ## We over-ride block-size in the case that get.tree is user specified
  block.size <- min(get.block.size.bits(block.size, ntree), sum(get.tree))
  ## Turn off partial option.
  partial.bits <- get.partial.bits(0)
  ## na.action bit
  ## in restore mode, we send in the protocol used for the training data.
  if (restore.mode) {
    na.action = object$na.action
  }
  na.action.bits <- get.na.action.bits(na.action)
  do.trace <- get.trace.bits(do.trace)
  ## The pivot in this predict related function is different from
  ## the pivot used in the grow related function.  In rfsrc we are
  ## referencing the list nativeOutput[[]].  Here we are referencing
  ## the $nativeArray[[]] object which is a massaged version of
  ## nativeOutput[[]].
  ## WARNING: Note that the maximum number of slots in the following
  ## foreign function call is 64.  Ensure that this limit is not
  ## exceeded.  Otherwise, the program will error on the call.
      pivot = 0
      chunk = 0
  ## set the maximum class levels
  max.class.levels <- 0
  ## Start the C external timer.
  ctime.external.start  <- proc.time()
  nativeOutput <- tryCatch({.Call("rfsrcPredict",
                                  as.integer(do.trace),
                                  as.integer(seed),
                                  as.integer(ensemble.bits +
                                             importance.bits +
                                             bootstrap.bits +
                                             proximity.bits +
                                             split.depth.bits +
                                             var.used.bits +
                                             outcome.bits +
                                             perf.bits +
                                             rfq.bits +
                                             cr.bits +
                                             gk.quantile.bits +
                                             anonymize.bits +
                                             case.depth.bits),
                                  as.integer(forest.wt.bits +
                                             distance.bits +
                                             samptype.bits +
                                             na.action.bits +
                                             membership.bits +
                                             terminal.qualts.bits +
                                             terminal.quants.bits +
                                             cse.bits +
                                             csv.bits +
                                             data.pass.bits +
                                             data.pass.predict.bits +
                                             jitt.bits),
                                  ## >>>> start of maxi forest object >>>>
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
                                  if (is.null(xvar)) NULL else as.double(as.vector(xvar)),
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
                                  as.integer(0),
                                  ## Deleted base learner slot.
                                  NULL,
                                  as.integer((object$nativeArray)$treeID),
                                  as.integer((object$nativeArray)$nodeID),
                                  as.integer((object$nativeArray)$nodeSZ),
                                  as.integer((object$nativeArray)$brnodeID),
                                  ## This is hc_zero.  It is never NULL.
                                  list(as.integer((object$nativeArray)$parmID),
                                  as.double((object$nativeArray)$contPT),
                                  as.integer((object$nativeArray)$mwcpSZ),
                                  as.integer((object$nativeArray)$fsrecID),
                                  if (is.null((object$nativeFactorArray)$mwcpPT)) NULL else as.integer((object$nativeFactorArray)$mwcpPT)),
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
                                  as.integer(ptn.count),
                                  list(if (is.null(marginal.xvar.idx)) as.integer(0) else as.integer(length(marginal.xvar.idx)),
                                       if (is.null(marginal.xvar.idx)) NULL else as.integer(marginal.xvar.idx)),
                                  list(if (is.null(importance.xvar.idx)) as.integer(0) else as.integer(length(importance.xvar.idx)),
                                       if (is.null(importance.xvar.idx)) NULL else as.integer(importance.xvar.idx)),
                                  ## Partial variables disabled.
                                  list(as.integer(0),
                                       as.integer(0),
                                       as.integer(0),
                                       NULL,
                                       as.integer(0),
                                       NULL,
                                       NULL),
                                  as.integer(n.newdata),
                                  as.integer(r.dim.newdata),
                                  as.double(if (outcome != "test") yvar.newdata else NULL),
                                  as.double(if (outcome != "test") xvar.newdata else NULL),
                                  as.integer(block.size),
                                  list(if (is.null(prob.assign$prob)) as.integer(0) else as.integer(length(prob.assign$prob)),
                                       if (is.null(prob.assign$prob)) NULL else as.double(prob.assign$prob),
                                       if (is.null(prob.assign$prob.epsilon)) as.double(0) else as.double(prob.assign$prob.epsilon)),
                                  as.integer(get.tree),
                                  as.integer(get.rf.cores()))}, error = function(e) {
                                    print(e)
                                    NULL})
  ## Stop the C external timer.
  ctime.external.stop <- proc.time()
  ## check for error return condition in the native code
  if (is.null(nativeOutput)) {
    stop("An error has occurred in prediction.  Please turn trace on for further analysis.")
  }
  ## determine missingness which is used for imputed data
  ## not done for anonymous forests or na.random imputation
  if (!inherits(object, "anonymous") && na.action != "na.random") {
    if (restore.mode) {
      n.miss <- get.nmiss(xvar, yvar)
    }
    else {
      n.miss <- get.nmiss(xvar.newdata, yvar.newdata)
     }
  }
  else {
    n.miss <- 0
  }
  ## extract the imputed data if there was missingness
  if (n.miss > 0) {
    imputed.data <- matrix(nativeOutput$imputation, nrow = n.miss)
    nativeOutput$imputation <- NULL
    imputed.indv <- imputed.data[, 1]
    imputed.data <- as.data.frame(imputed.data[, -1, drop = FALSE])
    if ((r.dim.newdata > 0) | (perf.type != "none")) {
      colnames(imputed.data) <- c(yvar.names, xvar.names)
    }
      else {
        colnames(imputed.data) <- xvar.names
      }
  }
  ## post-process the test data
  ## for restore mode there is no test data *except* when outcome=TEST
  if (!restore.mode | outcome == "test") {
    ## add row and column names to test xvar matrix
    xvar.newdata <- as.data.frame(xvar.newdata)
    rownames(xvar.newdata) <- newdata.row.names
    colnames(xvar.newdata) <- xvar.names
    ## map xvar factors back to original values
    xvar.newdata <- map.factor(xvar.newdata, xfactor)
    if (perf.type != "none") {
      ## add column names to test response matrix
      yvar.newdata <- as.data.frame(yvar.newdata)
      colnames(yvar.newdata) <- yvar.names
      ## map response factors back to original values
      yvar.newdata <- map.factor(yvar.newdata, yfactor)
    }
  }
  ## Map imputed data factors back to original values
  if (n.miss > 0) {
    imputed.data <- map.factor(imputed.data, xfactor)
    if (perf.type != "none") {
      imputed.data <- map.factor(imputed.data, yfactor)
    }
  }
  ## proximity
  if (proximity != FALSE) {
    if (restore.mode) {
      prox.n <- n
    }
      else {
        prox.n <- n.newdata
      }
    proximity.out <- matrix(0, prox.n, prox.n)
    count <- 0
    for (k in 1:prox.n) {
      proximity.out[k,1:k] <- nativeOutput$proximity[(count+1):(count+k)]
      proximity.out[1:k,k] <- proximity.out[k,1:k]
      count <- count + k
    }
    nativeOutput$proximity <- NULL
  }
    else {
      proximity.out <- NULL
    }
  ## distance
  if (distance != FALSE) {
    if (restore.mode) {
      dist.n <- n
    }
      else {
        dist.n <- n.newdata
      }
    distance.out <- matrix(0, dist.n, dist.n)
    count <- 0
    for (k in 1:dist.n) {
      distance.out[k,1:k] <- nativeOutput$distance[(count+1):(count+k)]
      distance.out[1:k,k] <- distance.out[k,1:k]
      count <- count + k
    }
    nativeOutput$distance <- NULL
  }
    else {
      distance.out <- NULL
    }
  ## forest weight matrix
  if (forest.wt != FALSE) {
    if (restore.mode) {
      forest.wt.n <- c(n, n)
    }
      else {
        forest.wt.n <- c(n.newdata, n)
      }
    forest.wt.out <- matrix(nativeOutput$weight, forest.wt.n, byrow = TRUE)
    nativeOutput$weight <- NULL
  }
    else {
      forest.wt.out <- NULL
    }
  ## forest case.depth
  if (case.depth != FALSE) {
    if (restore.mode) {
      case.depth.out <- matrix(nativeOutput$caseDepth, c(ntree, n), byrow = TRUE)
      nativeOutput$caseDepth <- NULL
    }
    else {
      case.depth.out <- matrix(nativeOutput$caseDepth, c(ntree, n.newdata), byrow = TRUE)
    }
    nativeOutput$caseDepth <- NULL
  }
  else {
    case.depth.out <- NULL
  }
  n.observed = if (restore.mode) n else n.newdata
  ## membership
  if (membership) {
    membership.out <- matrix(nativeOutput$nodeMembership, c(n.observed, ntree))
    nativeOutput$nodeMembership <- NULL
    if (restore.mode) {
      inbag.out <- matrix(nativeOutput$bootMembership, c(n.observed, ntree))
      nativeOutput$bootMembership <- NULL
    }
      else {
        inbag.out <- NULL
      }
    if (ptn.count > 0) {
      ptn.membership.out <- matrix(nativeOutput$pstnMembership, c(n.observed, ntree))
      nativeOutput$pstnMembership <- NULL
    }
      else {
        ptn.membership.out <- NULL
      }
  }
    else {
      membership.out <- NULL
      inbag.out <- NULL
      ptn.membership.out <- NULL
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
  ## make the output object
  rfsrcOutput <- list(
    call = match.call(),
    family = family,
    n = n.observed,
    ntree = ntree,
    yvar = (if ((outcome == "train" & restore.mode) | (perf.type != "none")) {
      if (outcome == "train" & restore.mode)
        amatrix.remove.names(object$yvar) else amatrix.remove.names(yvar.newdata)} else NULL),
    yvar.names = yvar.names,
    xvar = (if(outcome != "test" & restore.mode) object$xvar else xvar.newdata),
    xvar.names = xvar.names,
    leaf.count = nativeOutput$leafCount,
    proximity = proximity.out,
    forest = object,
    forest.wt = forest.wt.out,
    case.depth = case.depth.out,
    distance = distance.out,
    ptn.membership = ptn.membership.out,
    membership = membership.out,
    splitrule = splitrule,
    inbag = inbag.out,
    var.used = var.used.out,
    imputed.indv = (if (n.miss>0) imputed.indv else NULL),
    imputed.data = (if (n.miss>0) imputed.data else NULL),
    split.depth  = split.depth.out,
    block.size = block.size,
    perf.type = perf.type,
    ctime.internal = nativeOutput$cTimeInternal,
    ctime.external = ctime.external.stop - ctime.external.start
  )
  ## memory management
  nativeOutput$leafCount <- NULL
  remove(object)
  remove(proximity.out)
  remove(forest.wt.out)
  remove(case.depth.out)
  remove(distance.out)
  remove(ptn.membership.out)
  remove(membership.out)
  remove(inbag.out)
  if (n.miss > 0) remove(imputed.indv)
  if (n.miss > 0) remove(imputed.data)
  remove(var.used.out)
  remove(split.depth.out)
  ## Safe the outputs.
  survOutput <- NULL
  classOutput <- NULL
  regrOutput <- NULL
  if(vimp.joint) {
    vimp.count <- 1
  }
    else {
      vimp.count <- length(importance.xvar)
    }
  ## family specific additions to the predict object
  if (grepl("surv", family)) {
      if ((length(event.info$event.type) > 1) &&
          (splitrule != "l2.impute") &&
          (splitrule != "logrankscore")) {
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
      vimp.names <- list(NULL, if (vimp.joint) "joint" else importance.xvar)
    }
      else {
        ## Competing Risk names.
        ens.names <- list(NULL, NULL, c(paste("condCHF.", 1:length(event.info$event.type), sep = "")))
        mortality.names <- list(NULL, paste("event.", 1:length(event.info$event.type), sep = ""))
        cif.names <- list(NULL, NULL, c(paste("CIF.", 1:length(event.info$event.type), sep = "")))
        err.names <- list(c(paste("event.", 1:length(event.info$event.type), sep = "")), NULL)
        vimp.names <- list(paste("event.", 1:length(event.info$event.type), sep = ""),
                           if(vimp.joint) "joint" else importance.xvar)
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
                                 c(n.observed, length(event.info$time.interest), length(event.info$event.type)),
                                 dimnames=ens.names), length(event.info$event.type)) else NULL)
    nativeOutput$allEnsbCHF <- NULL
    survOutput <- list(chf = chf)
    remove(chf)
    chf.oob <- (if (!is.null(nativeOutput$oobEnsbCHF))
                  adrop3d.last(array(nativeOutput$oobEnsbCHF,
                                     c(n.observed, length(event.info$time.interest), length(event.info$event.type)),
                                     dimnames=ens.names), length(event.info$event.type)) else NULL)
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
                                       c(n.observed, length(event.info$event.type)), dimnames=mortality.names), coerced.event.count) else NULL)
    nativeOutput$allEnsbMRT <- NULL
    survOutput = c(survOutput, predicted = list(predicted))
    remove(predicted)
    predicted.oob <- (if (!is.null(nativeOutput$oobEnsbMRT))
                        adrop2d.last(array(nativeOutput$oobEnsbMRT,
                                           c(n.observed, length(event.info$event.type)), dimnames=mortality.names), coerced.event.count) else NULL)
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
                           c(n.observed, length(event.info$time.interest))) else NULL)
    nativeOutput$allEnsbSRV <- NULL
    survOutput <- c(survOutput, survival = list(survival))
    remove(survival)
    survival.oob <-  (if (!is.null(nativeOutput$oobEnsbSRV))
                        matrix(nativeOutput$oobEnsbSRV,
                               c(n.observed, length(event.info$time.interest))) else NULL)
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
                    c(n.observed, length(event.info$time.interest), length(event.info$event.type)),
                    dimnames=cif.names) else NULL)
    nativeOutput$allEnsbCIF <- NULL
    survOutput <- c(survOutput, cif = list(cif))
    remove(cif)
    cif.oob <- (if (!is.null(nativeOutput$oobEnsbCIF))
                  array(nativeOutput$oobEnsbCIF,
                        c(n.observed, length(event.info$time.interest), length(event.info$event.type)),
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
                                        c(length(event.info$event.type), vimp.count),
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
        ndead = (if (perf.type != "none") sum((if (restore.mode) yvar[, 2] else yvar.newdata[, 2]) !=0 , na.rm=TRUE) else NULL))
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
        ## Yeilds tree.offset = c(1, 8, ...)
        tree.offset <- array(1, ntree)
        ## Sum of all level counts targetted classification.  For example,
        ## if R1 has 3 levels, and R2 has 2 levels, and R3 has 6 levels, and
        ## only R1 and R3 and targetted, then levels.total = 9 + 3
        levels.total <- 0
        if (ntree > 1) {
          ## Iterate over all the target outcomes and map them to the class list.
          for (i in 1:length(m.target.idx)) {
            ## This is the slot in the class list.  The class list spans all the classification
            ## outputs, some of which may be NULL.
            target.idx <- which (class.index == m.target.idx[i])
            ## Is the target a classification y-variable.
            if (length(target.idx) > 0) {
              levels.total <- levels.total + 1 + levels.count[target.idx]
            }
          }
          tree.offset[2:ntree] <- levels.total
        }
        tree.offset <-  cumsum(tree.offset)
        ## The block offset calculation mirrors the tree offset calculation, but differs in only the primary dimension.
        block.offset <- array(1, floor(ntree/block.size))
        if (floor(ntree/block.size) > 1) {
          block.offset[2:floor(ntree/block.size)] <- levels.total
        }
        block.offset <-  cumsum(block.offset)
        ## Incoming vimp rates: V=xvar R=response L=level
        ## V1R1L0 V1R1L1 V1R1L2 V1R1L3 V1R1L0 V1R2L1 V1R2L2, V2R1L0 V2R1L1 V2R1L2 V2R1L3 V2R2L0 V2R2L1 V2R2L2, ...
        ## Yeilds vimp.offset = c(1, 8, ...)
        vimp.offset <- array(1, vimp.count)
        if (vimp.count > 1) {
          vimp.offset[2:vimp.count] <- levels.total
        }
        vimp.offset <-  cumsum(vimp.offset)
        iter.ensb.start <- 0
        iter.ensb.end   <- 0
        ## From the native code:
        ##   "allEnsbCLS"
        ##   "oobEnsbCLS"
        ## -> of dim [[class.count]] x [levels.count[]] x [n]
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
        ## where the slot [.] x [.] x [1] holds the unconditional error rate.
        ## Note that this is a ragged array.
        ## To the R code:
        ## -> of dim [[class.count]] x [block.cnt] x [1 + levels.count[]]
        ## From the native code:
        ##   "vimpClas"
        ## -> of dim [n.xvar] x [class.count] x [1 + levels.count]
        ## where the slot [.] x [.] x [1] holds the unconditional vimp.
        ## Note that this is a ragged array.
        ## To the R code:
        ## -> of dim [[class.count]] x [1 + levels.count] x [n.xvar]
        ## From the native code:
        ##   "cseNum"
        ## -> of dim [resp.class.count] x [n]
        ## To the R code:
        ## -> of dim [[resp.class.count]] x [n]
        ## From the native code:
        ##   "cseDen"
        ## -> of dim [n]
        ## To the R code:
        ## -> of dim [n]
        ## From the native code:
        ##   "csvNum"
        ## -> of dim [n.xvar] x [resp.class.count] x [n]
        ## To the R code:
        ## -> of dim [[resp.class.count]] x [n] x [n.xvar]
        ## From the native code:
        ##   "csvDen"
        ## -> of dim [n.xvar] x [n] 
        ## To the R code:
        ## -> of dim  x [n]  x [n.xvar]
        for (i in 1:length(m.target.idx)) {
          target.idx <- which(class.index == m.target.idx[i])
          if (length(target.idx) > 0) {
            iter.ensb.start <- iter.ensb.end
            iter.ensb.end <- iter.ensb.end + (levels.count[target.idx] * n.observed)
            ens.names <- list(NULL, levels.names[[target.idx]])
            err.names <- c("all", levels.names[[target.idx]])
            vimp.names <- list(c("all", levels.names[[target.idx]]), if(vimp.joint) "joint" else importance.xvar)
            predicted <- (if (!is.null(nativeOutput$allEnsbCLS))
                            array(nativeOutput$allEnsbCLS[(iter.ensb.start + 1):iter.ensb.end],
                                  c(n.observed, levels.count[target.idx]), dimnames=ens.names) else NULL)
            classOutput[[target.idx]] <- list(predicted = predicted)
            response <- (if (!is.null(predicted)) get.bayes.rule(predicted, class.relfrq) else NULL)
            classOutput[[target.idx]] <- c(classOutput[[target.idx]], class = list(response))
            remove(predicted)
            remove(response)
            predicted.oob <- (if (!is.null(nativeOutput$oobEnsbCLS))
                                array(nativeOutput$oobEnsbCLS[(iter.ensb.start + 1):iter.ensb.end],
                                      c(n.observed, levels.count[target.idx]), dimnames=ens.names) else NULL)
            classOutput[[target.idx]] <- c(classOutput[[target.idx]], predicted.oob = list(predicted.oob))
            response.oob <- (if (!is.null(predicted.oob)) get.bayes.rule(predicted.oob, class.relfrq) else NULL)
            classOutput[[target.idx]] <- c(classOutput[[target.idx]], class.oob = list(response.oob))
            remove(predicted.oob)
            remove(response.oob)
            ## New case specific arrays, numerator:
            cse.num <- (if (!is.null(nativeOutput$cseClas))
                            array(nativeOutput$cseClas[(iter.ensb.start + 1):iter.ensb.end], n) else NULL)
            classOutput[[target.idx]] <- c(classOutput[[target.idx]], cse.num = list(cse.num))
            remove(cse.num)
            ## New case specific arrays, denominator:
            cse.den <- (if (!is.null(nativeOutput$cseDen))
                            array(nativeOutput$cseDen[(iter.ensb.start + 1):iter.ensb.end], n) else NULL)
            classOutput[[target.idx]] <- c(classOutput[[target.idx]], cse.den = list(cse.den))
            remove(cse.den)
            if (!is.null(nativeOutput$perfClas)) {
              err.rate <- array(0, c(1 + levels.count[target.idx], ntree))
              for (j in 1: (1 + levels.count[target.idx])) {
                err.rate[j, ]  <- nativeOutput$perfClas[tree.offset]
                tree.offset <- tree.offset + 1
              }
              row.names(err.rate) <- err.names
              classOutput[[target.idx]] <- c(classOutput[[target.idx]], err.rate = list(t(err.rate)))
              remove(err.rate)
            }
            if (!is.null(nativeOutput$blockClas)) {
              err.block.rate <- array(0, c(1 + levels.count[target.idx], floor(ntree/block.size)))
              for (j in 1: (1 + levels.count[target.idx])) {
                err.block.rate[j, ]  <- nativeOutput$blockClas[block.offset]
                block.offset <- block.offset + 1
              }
              row.names(err.block.rate) <- err.names
              classOutput[[target.idx]] <- c(classOutput[[target.idx]], err.block.rate = list(t(err.block.rate)))
              remove(err.block.rate)
            }
            if (!is.null(nativeOutput$vimpClas)) {
              importance <- array(0, c(1 + levels.count[target.idx], vimp.count), dimnames=vimp.names)
              for (j in 1: (1 + levels.count[target.idx])) {
                importance[j, ]  <- nativeOutput$vimpClas[vimp.offset]
                vimp.offset <- vimp.offset + 1
              }
              classOutput[[target.idx]] <- c(classOutput[[target.idx]], importance = list(t(importance)))
              remove(importance)
            }
            ## See GROW call for explanation of arrays:
            csv.idx  <- array(0, n.observed * n.xvar) 
            j2seq <- (target.idx - 1) * n.observed
            for (j in 1:n.xvar) {
              j1seq <- (j - 1) * n.observed
              csv.idx[(j1seq + 1) : (j1seq + n.observed)]  <- ( (j2seq + j1seq * class.count + 1) : (j2seq + j1seq * class.count + n.observed) ) 
            }
            ## New case specific arrays, numerator:
            csv.num <- (if (!is.null(nativeOutput$csvClas))
                          array(nativeOutput$csvClas[csv.idx], c(n.observed, n.xvar)) else NULL)
            classOutput[[target.idx]] <- c(classOutput[[target.idx]], csv.num = list(csv.num))
            remove(csv.num)
            ## New case specific arrays, denominator:
            csv.den <- (if (!is.null(nativeOutput$csvDen))
                          array(nativeOutput$csvDen, c(n.observed, n.xvar)) else NULL)
            classOutput[[target.idx]] <- c(classOutput[[target.idx]], csv.den = list(csv.den))
            remove(csv.den)
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
        ## Yeilds tree.offset = c(1, 3, 5)
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
        ## Yeilds vimp.offset = c(1, 3, 5, ...)
        vimp.offset <- array(1, vimp.count)
        if (vimp.count > 1) {
          vimp.offset[2:vimp.count] <- length(regr.index)
        }
        vimp.offset <-  cumsum(vimp.offset)
        iter.ensb.start <- 0
        iter.ensb.end   <- 0
        iter.qntl.start <- 0
        iter.qntl.end   <- 0
        ## From the native code:
        ##   "allEnsbRGR"
        ## -> of dim [regr.count] x [obsSize]
        ## From the native code:
        ##   "allEnsbQNT"
        ##   "oobEnsbQNT"
        ## -> of dim [regr.count] x [length(prob)] x [n]
        ## From the native code:
        ##   "perfRegr"
        ## -> of dim [regr.count] x [ntree]
        ## To the R code:
        ## -> of dim  [[regr.count]] x [ntree]
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
        ## From the native code:
        ##   "cseNum"
        ## -> of dim [regr.count] x [n]
        ## To the R code:
        ## -> of dim [[regr.count]] x [n]
        ## From the native code:
        ##   "cseDen"
        ## -> of dim [n]
        ## To the R code:
        ## -> of dim [n]
        ## From the native code:
        ##   "csvNum"
        ## -> of dim [n.xvar] x [regr.count] x [n]
        ## To the R code:
        ## -> of dim [[regr.count]] x [n] x [n.xvar]
        ## From the native code:
        ##   "csvDen"
        ## -> of dim [n.xvar] x [n] 
        ## To the R code:
        ## -> of dim  x [n]  x [n.xvar]
        for (i in 1:length(m.target.idx)) {
          target.idx <- which(regr.index == m.target.idx[i])
          if (length(target.idx) > 0) {
            iter.ensb.start <- iter.ensb.end
            iter.ensb.end <- iter.ensb.end + n.observed
            iter.qntl.start <- iter.qntl.end
            iter.qntl.end <- iter.qntl.end + (length(prob) * n.observed)
            vimp.names <- if(vimp.joint) "joint" else importance.xvar
            predicted <- (if (!is.null(nativeOutput$allEnsbRGR))
                            array(nativeOutput$allEnsbRGR[(iter.ensb.start + 1):iter.ensb.end], n.observed) else NULL)
            regrOutput[[target.idx]] <- list(predicted = predicted)
            remove(predicted)
            predicted.oob <- (if (!is.null(nativeOutput$oobEnsbRGR))
                                array(nativeOutput$oobEnsbRGR[(iter.ensb.start + 1):iter.ensb.end], n.observed) else NULL)
            regrOutput[[target.idx]] <- c(regrOutput[[target.idx]], predicted.oob = list(predicted.oob))
            remove(predicted.oob)
            ## New case specific arrays, numerator:
            cse.num <- (if (!is.null(nativeOutput$cseRegr))
                            array(nativeOutput$cseRegr[(iter.ensb.start + 1):iter.ensb.end], n) else NULL)
            regrOutput[[target.idx]] <- c(regrOutput[[target.idx]], cse.num = list(cse.num))
            remove(cse.num)
            ## New case specific arrays, denominator:
            cse.den <- (if (!is.null(nativeOutput$cseDen))
                            array(nativeOutput$cseDen[(iter.ensb.start + 1):iter.ensb.end], n) else NULL)
            regrOutput[[target.idx]] <- c(regrOutput[[target.idx]], cse.den = list(cse.den))
            remove(cse.den)
            quantile <- (if (!is.null(nativeOutput$allEnsbQNT))
                             array(nativeOutput$allEnsbQNT[(iter.qntl.start + 1):iter.qntl.end],
                                   c(n.observed, length(prob))) else NULL)
            regrOutput[[target.idx]] <- c(regrOutput[[target.idx]], quantile = list(quantile))
            remove(quantile)
            quantile.oob <- (if (!is.null(nativeOutput$oobEnsbQNT))
                                 array(nativeOutput$oobEnsbQNT[(iter.qntl.start + 1):iter.qntl.end],
                                       c(n.observed, length(prob))) else NULL)
            regrOutput[[target.idx]] <- c(regrOutput[[target.idx]], quantile.oob = list(quantile.oob))
            remove(quantile.oob)
            if (!is.null(nativeOutput$perfRegr)) {
              err.rate <- nativeOutput$perfRegr[tree.offset]
              tree.offset <- tree.offset + 1
              regrOutput[[target.idx]] <- c(regrOutput[[target.idx]], err.rate = list(err.rate))
              remove(err.rate)
            }
            if (!is.null(nativeOutput$blockRegr)) {
              err.block.rate <- nativeOutput$blockRegr[block.offset]
              block.offset <- block.offset + 1
              regrOutput[[target.idx]] <- c(regrOutput[[target.idx]], err.block.rate = list(err.block.rate))
              remove(err.block.rate)
            }
            if (!is.null(nativeOutput$vimpRegr)) {
              importance <- nativeOutput$vimpRegr[vimp.offset]
              names(importance) <- vimp.names
              vimp.offset <- vimp.offset + 1
              regrOutput[[target.idx]] <- c(regrOutput[[target.idx]], importance = list(importance))
              remove(importance)
            }
            ## See GROW call for explanation of arrays:
            csv.idx  <- array(0, n.observed * n.xvar)
            j2seq <- (target.idx - 1) * n.observed
            for (j in 1:n.xvar) {
              j1seq <- (j - 1) * n.observed
              csv.idx[(j1seq + 1) : (j1seq + n.observed)]  <- ( (j2seq + j1seq * regr.count + 1) : (j2seq + j1seq * regr.count + n.observed) ) 
            }
            ## New case specific arrays, numerator:
            csv.num <- (if (!is.null(nativeOutput$csvRegr))
                          array(nativeOutput$csvRegr[csv.idx], c(n.observed, n.xvar)) else NULL)
            regrOutput[[target.idx]] <- c(regrOutput[[target.idx]], csv.num = list(csv.num))
            remove(csv.num)
            ## New case specific arrays, denominator:
            csv.den <- (if (!is.null(nativeOutput$csvDen))
                          array(nativeOutput$csvDen, c(n.observed, n.xvar)) else NULL)
            regrOutput[[target.idx]] <- c(regrOutput[[target.idx]], csv.den = list(csv.den))
            remove(csv.den)
          }
        }
        nativeOutput$allEnsbRGR <- NULL
        nativeOutput$oobEnsbRGR <- NULL
        nativeOutput$allEnsbQNT <- NULL
        nativeOutput$oobEnsbQNT <- NULL
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
  class(rfsrcOutput) <- c("rfsrc", "predict",   family)
  return(rfsrcOutput)
}
