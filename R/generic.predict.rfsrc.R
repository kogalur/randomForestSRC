##  **********************************************************************
##  **********************************************************************
##  
##    RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
##  
##    This program is free software; you can redistribute it and/or
##    modify it under the terms of the GNU General Public License
##    as published by the Free Software Foundation; either version 3
##    of the License, or (at your option) any later version.
##  
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##  
##    You should have received a copy of the GNU General Public
##    License along with this program; if not, write to the Free
##    Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
##    Boston, MA  02110-1301, USA.
##  
##    ----------------------------------------------------------------
##    Project Partially Funded By: 
##    ----------------------------------------------------------------
##    Dr. Ishwaran's work was funded in part by DMS grant 1148991 from the
##    National Science Foundation and grant R01 CA163739 from the National
##    Cancer Institute.
##  
##    Dr. Kogalur's work was funded in part by grant R01 CA163739 from the 
##    National Cancer Institute.
##    ----------------------------------------------------------------
##    Written by:
##    ----------------------------------------------------------------
##      Hemant Ishwaran, Ph.D.
##      Director of Statistical Methodology
##      Professor, Division of Biostatistics
##      Clinical Research Building, Room 1058
##      1120 NW 14th Street
##      University of Miami, Miami FL 33136
##  
##      email:  hemant.ishwaran@gmail.com
##      URL:    http://web.ccs.miami.edu/~hishwaran
##      --------------------------------------------------------------
##      Udaya B. Kogalur, Ph.D.
##      Adjunct Staff
##      Department of Quantitative Health Sciences
##      Cleveland Clinic Foundation
##      
##      Kogalur & Company, Inc.
##      5425 Nestleway Drive, Suite L1
##      Clemmons, NC 27012
##  
##      email:  ubk@kogalur.com
##      URL:    https://github.com/kogalur/randomForestSRC
##      --------------------------------------------------------------
##  
##  **********************************************************************
##  **********************************************************************


generic.predict.rfsrc <-
  function(object,
           newdata,
           outcome.target = NULL, 
           importance = c(FALSE, TRUE, "none", "permute", "random", "anti", "permute.ensemble", "random.ensemble", "anti.ensemble"),
           importance.xvar,
           subset = NULL,
           na.action = c("na.omit", "na.impute"),
           outcome = c("train", "test"),
           proximity = FALSE,
           var.used = c(FALSE, "all.trees", "by.tree"),
           split.depth = c(FALSE, "all.trees", "by.tree"),
           seed = NULL,
           do.trace = FALSE,
           membership = FALSE,
           tree.err = FALSE,
           statistics = FALSE,
           ...)
{
  univariate.nomenclature = TRUE
  user.option <- list(...)
  terminal.qualts <- is.hidden.terminal.qualts(user.option)
  terminal.quants <- is.hidden.terminal.quants(user.option)
  ptn.count <- is.hidden.ptn.count(user.option)
  if (missing(object)) {
    stop("object is missing!")
  }
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2)  
    stop("this function only works for objects of class `(rfsrc, grow)' or '(rfsrc, forest)'")
  importance <- match.arg(as.character(importance), c(FALSE, TRUE,
                                                      "none", "permute", "random", "anti",
                                                      "permute.ensemble", "random.ensemble", "anti.ensemble",
                                                      "permute.joint", "random.joint", "anti.joint",
                                                      "permute.joint.ensemble", "random.joint.ensemble", "anti.joint.ensemble"))
  if (grepl("joint", importance)) {
    vimp.joint <- TRUE
  }
    else {
      vimp.joint <- FALSE
    }
  xvar.names <- object$xvar.names
  yvar.names <- object$yvar.names
  importance.xvar <- get.importance.xvar(importance.xvar, importance, object)
  importance.xvar.idx <- match(importance.xvar, xvar.names)
  na.action <- match.arg(na.action, c("na.omit", "na.impute"))
  outcome <- match.arg(outcome, c("train", "test"))
  proximity <- match.arg(as.character(proximity), c(FALSE, TRUE, "inbag", "oob", "all"))
  var.used <- match.arg(as.character(var.used), c("FALSE", "all.trees", "by.tree"))
  if (var.used == "FALSE") var.used <- FALSE
  split.depth <- match.arg(as.character(split.depth),  c("FALSE", "all.trees", "by.tree"))
  if (split.depth == "FALSE") split.depth <- FALSE
  seed <- get.seed(seed)
  if (missing(newdata)) {
    outcome <- "train"
    grow.equivalent <- TRUE
  }
    else {
      grow.equivalent <- FALSE
    }
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.null(object$forest)) {
      stop("The forest is empty.  Re-run rfsrc (grow) call with forest=TRUE")
    }
    if (inherits(object, "bigdata")) {
      big.data <- TRUE
    }
      else {
        big.data <- FALSE
      }
    object <- object$forest
  }
    else {
      if (inherits(object, "bigdata")) {
        big.data <- TRUE
      }
        else {
          big.data <- FALSE
        }
    }
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
      installed.version <- as.integer(unlist(strsplit("2.4.2", "[.]")))
      minimum.version <- as.integer(unlist(strsplit("2.3.0", "[.]")))
      object.version.adj <- object.version[1] + (object.version[2]/10) + (object.version[3]/100)
      installed.version.adj <- installed.version[1] + (installed.version[2]/10) + (installed.version[3]/100)
      minimum.version.adj <- minimum.version[1] + (minimum.version[2]/10) + (minimum.version[3]/100) 
      if (object.version.adj >= minimum.version.adj) {
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
  splitrule <- object$splitrule
  object$yvar <- as.data.frame(object$yvar)
  colnames(object$yvar) <- yvar.names
  yfactor <- extract.factor(object$yvar)
  family <- object$family
  outcome.target.idx <- get.outcome.target(family, yvar.names, outcome.target)
  yvar.types <- get.yvar.type(family, yfactor$generic.types, yvar.names, object$coerce.factor)
  yvar.nlevels <- get.yvar.nlevels(family, yfactor$nlevels, yvar.names, object$yvar, object$coerce.factor)
  event.info <- get.event.info(object)
  cr.bits <- get.cr.bits(family)
  xfactor <- extract.factor(object$xvar)
  any.xvar.factor <-  (length(xfactor$factor) + length(xfactor$order)) > 0
  xvar.types <- get.xvar.type(xfactor$generic.types, xvar.names, object$coerce.factor)
  xvar.nlevels <- get.xvar.nlevels(xfactor$nlevels, xvar.names, object$xvar, object$coerce.factor)
  if (family == "unsupv") {
    outcome <- "train"
    perf.flag <- FALSE
    importance <- "none"
  }
  if (grepl("surv", family)) {
    ptn.count <- 0
  }
  if (!grow.equivalent) {
    newdata <- newdata[, is.element(names(newdata),
                                    c(yvar.names, xvar.names)), drop = FALSE]
    newdata <- rm.na.levels(newdata, xvar.names)
    newdata.xfactor <- extract.factor(newdata, xvar.names)
    if (!setequal(xfactor$factor, newdata.xfactor$factor)) {
      stop("x-variable factors from test data do not match original training data")
    }
    if (!setequal(xfactor$order, newdata.xfactor$order)) {
      stop("(ordered) x-variable factors from test data do not match original training data")
    }
    any.outcome.factor <- family == "class"
    if (family == "class+" | family ==  "mix+") {
      if (length(intersect("R", yfactor$generic.types[outcome.target.idx])) == 0) {
        any.outcome.factor <- TRUE
      }
    }
    if (any.outcome.factor) {
      if (sum(is.element(names(newdata), yvar.names)) > 0) {
        newdata <- rm.na.levels(newdata, yvar.names)
        newdata.yfactor <- extract.factor(newdata, yvar.names)
        if (!setequal(yfactor$factor, newdata.yfactor$factor)) {
          stop("class outcome from test data does not match original training data")
        }
        if (!setequal(yfactor$order, newdata.yfactor$order)) {
          stop("(ordered) class outcome from test data does not match original training data")
        }
      }
    }
    if (length(xvar.names) != sum(is.element(xvar.names, names(newdata)))) {
      stop("x-variables in test data do not match original training data")
    }
    yvar.present <- sum(is.element(yvar.names, names(newdata))) > 0
    if (yvar.present && length(yvar.names) != sum(is.element(yvar.names, names(newdata)))) {
      stop("y-variables in test data do not match original training data")
    }
    if (any.xvar.factor) {
      newdata <- check.factor(object$xvar, newdata, xfactor)
    }
    if (any.outcome.factor) {
      if (yvar.present) {
        newdata <- check.factor(object$yvar, newdata, yfactor)
      }
    }
    if (yvar.present) {
      fnames <- c(yvar.names, xvar.names)
    }
      else {
        fnames <- xvar.names
      }
    newdata <- finalizeData(fnames, newdata, na.action)
    xvar.newdata  <- as.matrix(newdata[, xvar.names, drop = FALSE])
    n.newdata <- nrow(newdata)
    newdata.row.names <- rownames(xvar.newdata)
    if (yvar.present) {
      yvar.newdata <- as.matrix(newdata[, yvar.names, drop = FALSE])
      event.info.newdata <- get.grow.event.info(yvar.newdata, family, need.deaths = FALSE)
      r.dim.newdata <- event.info.newdata$r.dim
      perf.flag <- TRUE
      if (grepl("surv", family) && all(na.omit(event.info.newdata$cens) == 0)) {
        perf.flag <- FALSE
        importance <- "none"
      }
      if (grepl("surv", family) &&
          length(setdiff(na.omit(event.info.newdata$cens), na.omit(event.info$cens))) > 1) {
        stop("survival events in test data do not match training data")
      }
    }
      else {
        if (outcome == "test") {
          stop("outcome=TEST, but the test data has no y values, which is not permitted")
        }
        r.dim.newdata <- 0
        yvar.newdata <-  NULL
        perf.flag <- FALSE
        importance <- "none"
      }
    if (outcome != "test") {
      rownames(xvar.newdata) <- colnames(xvar.newdata) <- NULL
    }
    remove(newdata)
  }
    else {
      n.newdata <- 0
      r.dim.newdata <- 0
      xvar.newdata <- NULL
      yvar.newdata <-  NULL
      outcome <- "train"
      if (object$bootstrap != "by.root" | family == "unsupv") {
        importance <- "none"
        perf.flag <- FALSE
      }
        else {
          perf.flag <- TRUE
        }
    } 
  if (outcome == "train") {
    xvar <- as.matrix(data.matrix(object$xvar))
    yvar <- as.matrix(data.matrix(object$yvar))
    sampsize <- object$sampsize
    case.wt <- object$case.wt
    samp <- object$samp
  }
    else {
      xvar <- xvar.newdata
      yvar <- yvar.newdata
      grow.equivalent <- TRUE
      n.newdata <- 0
      r.dim.newdata <- 0
      sampsize <- nrow(xvar)
      case.wt <- get.weight(NULL, nrow(xvar))
      samp <- NULL
    }
  r.dim <- ncol(cbind(yvar))
  rownames(xvar) <- colnames(xvar) <- NULL
  n.xvar <- ncol(xvar)
  n <- nrow(xvar)
  split.null <- object$split.null
  ntree <- object$ntree
  importance.bits <- get.importance(importance)
  proximity.bits <- get.proximity(grow.equivalent, proximity)
  split.null.bits <- get.split.null(split.null)
  split.depth.bits <- get.split.depth(split.depth)
  var.used.bits <- get.var.used(var.used)
  outcome.bits <- get.outcome(outcome)
  perf.bits <-  get.perf.bits(perf.flag)
  statistics.bits <- get.statistics(statistics)
  bootstrap.bits <- get.bootstrap(object$bootstrap)
  samptype.bits <- get.samptype(object$samptype)
  membership.bits <-  get.membership(membership)
  terminal.qualts.bits <- get.terminal.qualts(terminal.qualts, object$terminal.qualts)
  terminal.quants.bits <- get.terminal.quants(terminal.quants, object$terminal.quants)
  tree.err.bits <- get.tree.err(tree.err)
  partial.bits <- get.partial(0)
  if (outcome == "test") {
  }
    else {
      na.action = object$na.action
    }
  na.action.bits <- get.na.action(na.action)
  if (missing(subset) | is.null(subset)) {
    subset <- NULL
  }
    else {
      if (is.logical(subset)) {
        subset <- which(subset)
      }
      subset <- unique(subset[subset >= 1 & subset <= n])
      if (length(subset) == 0) {
        stop("'subset' not set properly")
      }
    }
  do.trace <- get.trace(do.trace)
  nativeOutput <- tryCatch({.Call("rfsrcPredict",
                                  as.integer(do.trace),
                                  as.integer(seed),
                                  as.integer(importance.bits +
                                               bootstrap.bits +
                                                 proximity.bits +
                                                   split.null.bits +
                                                     split.depth.bits +
                                                       var.used.bits +
                                                         outcome.bits +
                                                           perf.bits +
                                                             cr.bits +
                                                               statistics.bits),
                                  as.integer(
                                        samptype.bits +
                                          na.action.bits +
                                            tree.err.bits +
                                              membership.bits +
                                                terminal.qualts.bits +
                                                  terminal.quants.bits),
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
                                  as.integer(sampsize),
                                  as.integer(samp),
                                  as.double(case.wt),
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
                                  as.integer(ptn.count),
                                  as.integer(length(importance.xvar.idx)),
                                  as.integer(importance.xvar.idx),
                                  as.integer(length(subset)),
                                  as.integer(subset),
                                  as.integer(0),
                                  as.integer(0),
                                  as.integer(0),
                                  as.double(NULL),
                                  as.integer(0),
                                  as.integer(NULL),
                                  as.double(NULL),
                                  as.integer(n.newdata),
                                  as.integer(r.dim.newdata),
                                  as.double(if (outcome != "test") yvar.newdata else NULL),
                                  as.double(if (outcome != "test") xvar.newdata else NULL),
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
  if (grow.equivalent) {
    n.miss <- get.nmiss(xvar, yvar)
  }
    else {
      n.miss <- get.nmiss(xvar.newdata, yvar.newdata)
    }
  if (n.miss > 0) {
    imputed.data <- matrix(nativeOutput$imputation, nrow = n.miss)
    nativeOutput$imputation <- NULL
    imputed.indv <- imputed.data[, 1]
    imputed.data <- as.data.frame(imputed.data[, -1, drop = FALSE])
    if (r.dim.newdata > 0 | perf.flag) {
      colnames(imputed.data) <- c(yvar.names, xvar.names)
    }
      else {
        colnames(imputed.data) <- xvar.names
      }
  }
  if (!grow.equivalent | outcome == "test") {
    xvar.newdata <- as.data.frame(xvar.newdata)
    rownames(xvar.newdata) <- newdata.row.names
    colnames(xvar.newdata) <- xvar.names
    xvar.newdata <- map.factor(xvar.newdata, xfactor)
    if (perf.flag) {
      yvar.newdata <- as.data.frame(yvar.newdata)
      colnames(yvar.newdata) <- yvar.names
      yvar.newdata <- map.factor(yvar.newdata, yfactor)
    }
  }
  if (n.miss > 0) {
    imputed.data <- map.factor(imputed.data, xfactor)
    if (perf.flag) {
      imputed.data <- map.factor(imputed.data, yfactor)
    }
  }
  if (proximity != FALSE) {
    if (grow.equivalent) {
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
  n.observed = if (grow.equivalent) n else n.newdata
  if (membership) {
    membership.out <- matrix(nativeOutput$nodeMembership, c(n.observed, ntree))
    nativeOutput$nodeMembership <- NULL
    if (grow.equivalent) {
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
  if (statistics == TRUE) {
    node.stats <- as.data.frame(cbind(nativeOutput$spltST))
    names(node.stats) <- c("spltST")
  }
    else {
      node.stats <- NULL
    }
  rfsrcOutput <- list(
    call = match.call(),
    family = family,
    n = n.observed,
    ntree = ntree,
    yvar = (if ((outcome == "train" & grow.equivalent) | perf.flag) {
      if (outcome == "train" & grow.equivalent)
        amatrix.remove.names(object$yvar) else amatrix.remove.names(yvar.newdata)} else NULL),
    yvar.names = yvar.names,
    xvar = (if(outcome != "test" & grow.equivalent) object$xvar else xvar.newdata),
    xvar.names = xvar.names,
    leaf.count = nativeOutput$leafCount,
    proximity = proximity.out,
    forest = object,
    ptn.membership = ptn.membership.out,
    membership = membership.out,
    splitrule = splitrule,
    inbag = inbag.out,
    var.used = var.used.out,
    imputed.indv = (if (n.miss>0) imputed.indv else NULL),
    imputed.data = (if (n.miss>0) imputed.data else NULL),
    split.depth  = split.depth.out,
    node.stats = node.stats,
    tree.err = tree.err
  )
  nativeOutput$leafCount <- NULL
  remove(object)
  remove(proximity.out)
  remove(ptn.membership.out)
  remove(membership.out)
  remove(inbag.out)
  if (n.miss > 0) remove(imputed.indv)
  if (n.miss > 0) remove(imputed.data)
  remove(var.used.out)
  remove(split.depth.out)
  remove(node.stats)
  survOutput <- NULL
  classOutput <- NULL
  regrOutput <- NULL
  if(vimp.joint) {
    vimp.count <- 1
  }
    else {
      vimp.count <- length(importance.xvar)
    }
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
      ens.names <- list(NULL, NULL)
      mortality.names <- list(NULL, NULL)
      err.names <- list(NULL, NULL)
      vimp.names <- list(NULL, if (vimp.joint) "joint" else importance.xvar)
    }
      else {
        ens.names <- list(NULL, NULL, c(paste("condCHF.", 1:length(event.info$event.type), sep = "")))
        mortality.names <- list(NULL, paste("event.", 1:length(event.info$event.type), sep = ""))
        cif.names <- list(NULL, NULL, c(paste("CIF.", 1:length(event.info$event.type), sep = "")))
        err.names <- list(c(paste("event.", 1:length(event.info$event.type), sep = "")), NULL)
        vimp.names <- list(paste("event.", 1:length(event.info$event.type), sep = ""),
                           if(vimp.joint) "joint" else importance.xvar)
      }
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
        ndead = (if (perf.flag) sum((if (grow.equivalent) yvar[, 2] else yvar.newdata[, 2]) !=0 , na.rm=TRUE) else NULL))
    )
    if(univariate.nomenclature) {
      rfsrcOutput <- c(rfsrcOutput, survOutput)
    }
      else {
        rfsrcOutput <- c(rfsrcOutput, survOutput = list(survOutput))
      }
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
        tree.offset <- array(1, ntree)
        levels.total <- 0 
        if (ntree > 1) {
          for (i in 1:length(outcome.target.idx)) {
            target.idx <- which (class.index == outcome.target.idx[i])
            if (length(target.idx) > 0) {
              levels.total <- levels.total + 1 + levels.count[target.idx]
            }
          }
          tree.offset[2:ntree] <- levels.total
        }
        tree.offset <-  cumsum(tree.offset)
        vimp.offset <- array(1, vimp.count)
        if (vimp.count > 1) {
          vimp.offset[2:vimp.count] <- levels.total
        }
        vimp.offset <-  cumsum(vimp.offset)
        iter.ensb.start <- 0
        iter.ensb.end   <- 0
        for (i in 1:length(outcome.target.idx)) {
          target.idx <- which (class.index == outcome.target.idx[i])
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
            response <- (if (!is.null(predicted)) bayes.rule(predicted) else NULL)
            classOutput[[target.idx]] <- c(classOutput[[target.idx]], class = list(response))
            remove(predicted)
            remove(response)
            predicted.oob <- (if (!is.null(nativeOutput$oobEnsbCLS))
                                array(nativeOutput$oobEnsbCLS[(iter.ensb.start + 1):iter.ensb.end],
                                      c(n.observed, levels.count[target.idx]), dimnames=ens.names) else NULL)
            classOutput[[target.idx]] <- c(classOutput[[target.idx]], predicted.oob = list(predicted.oob))
            response.oob <- (if (!is.null(predicted.oob)) bayes.rule(predicted.oob) else NULL)
            classOutput[[target.idx]] <- c(classOutput[[target.idx]], class.oob = list(response.oob))
            remove(predicted.oob)
            remove(response.oob)
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
            if (!is.null(nativeOutput$vimpClas)) {
              importance <- array(0, c(1 + levels.count[target.idx], vimp.count), dimnames=vimp.names)
              for (j in 1: (1 + levels.count[target.idx])) {
                importance[j, ]  <- nativeOutput$vimpClas[vimp.offset]
                vimp.offset <- vimp.offset + 1
              }
              classOutput[[target.idx]] <- c(classOutput[[target.idx]], importance = list(t(importance)))
              remove(importance)
            }
          }
        }
        nativeOutput$allEnsbCLS <- NULL
        nativeOutput$oobEnsbCLS <- NULL
        nativeOutput$perfClas <- NULL
        nativeOutput$vimpClas <- NULL
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
        tree.offset <- array(1, ntree)
        if (ntree > 1) {
          tree.offset[2:ntree] <- length(regr.index)
        }
        tree.offset <-  cumsum(tree.offset)
        vimp.offset <- array(1, vimp.count)
        if (vimp.count > 1) {
          vimp.offset[2:vimp.count] <- length(regr.index)
        }
        vimp.offset <-  cumsum(vimp.offset)
        iter.ensb.start <- 0
        iter.ensb.end   <- 0
        for (i in 1:length(outcome.target.idx)) {
          target.idx <- which (regr.index == outcome.target.idx[i])
          if (length(target.idx) > 0) {
            iter.ensb.start <- iter.ensb.end
            iter.ensb.end <- iter.ensb.end + n.observed
            vimp.names <- if(vimp.joint) "joint" else importance.xvar
            predicted <- (if (!is.null(nativeOutput$allEnsbRGR))
                            array(nativeOutput$allEnsbRGR[(iter.ensb.start + 1):iter.ensb.end], n.observed) else NULL)
            regrOutput[[target.idx]] <- list(predicted = predicted)
            remove(predicted)
            predicted.oob <- (if (!is.null(nativeOutput$oobEnsbRGR))
                                array(nativeOutput$oobEnsbRGR[(iter.ensb.start + 1):iter.ensb.end], n.observed) else NULL)
            regrOutput[[target.idx]] <- c(regrOutput[[target.idx]], predicted.oob = list(predicted.oob))
            remove(predicted.oob)
            if (!is.null(nativeOutput$perfRegr)) {
              err.rate <- nativeOutput$perfRegr[tree.offset]
              tree.offset <- tree.offset + 1
              regrOutput[[target.idx]] <- c(regrOutput[[target.idx]], err.rate = list(err.rate))
              remove(err.rate)
            }
            if (!is.null(nativeOutput$vimpRegr)) {
              importance <- nativeOutput$vimpRegr[vimp.offset]
              names(importance) <- vimp.names
              vimp.offset <- vimp.offset + 1
              regrOutput[[target.idx]] <- c(regrOutput[[target.idx]], importance = list(importance))
              remove(importance)
            }
          }
        }
        nativeOutput$allEnsbRGR <- NULL
        nativeOutput$oobEnsbRGR <- NULL
        nativeOutput$perfRegr <- NULL
        nativeOutput$vimpRegr <- NULL
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
