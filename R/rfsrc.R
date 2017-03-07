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


rfsrc <- function(formula,
                  data,
                  ntree = 1000,
                  bootstrap = c("by.root", "by.node", "none", "by.user"),
                  mtry = NULL,
                  nodesize = NULL,
                  nodedepth = NULL,
                  splitrule = NULL,
                  nsplit = 0,
                  split.null = FALSE,
                  importance = c(FALSE, TRUE, "none", "permute", "random", "anti", "permute.ensemble", "random.ensemble", "anti.ensemble"),
                  na.action = c("na.omit", "na.impute"),
                  nimpute = 1,
                  ntime,
                  cause,
                  proximity = FALSE,
                  sampsize = NULL,
                  samptype = c("swr", "swor"),
                  samp = NULL,
                  case.wt = NULL,
                  xvar.wt = NULL,
                  forest = TRUE,
                  var.used = c(FALSE, "all.trees", "by.tree"),
                  split.depth = c(FALSE, "all.trees", "by.tree"),
                  seed = NULL,
                  do.trace = FALSE,
                  membership = FALSE,
                  statistics = FALSE,
                  tree.err = FALSE,
                  coerce.factor = NULL,
                  ...)
{
  univariate.nomenclature = TRUE
  user.option <- list(...)
  impute.only <- is.hidden.impute.only(user.option)
  miss.tree <- is.hidden.impute.only(user.option)
  terminal.qualts <- is.hidden.terminal.qualts(user.option)
  terminal.quants <- is.hidden.terminal.quants(user.option)
  bootstrap <- match.arg(bootstrap, c("by.root", "by.node", "none", "by.user"))
  importance <- match.arg(as.character(importance), c(FALSE, TRUE, "none", "permute", "random", "anti", "permute.ensemble", "random.ensemble", "anti.ensemble"))
  na.action <- match.arg(na.action, c("na.omit", "na.impute"))
  proximity <- match.arg(as.character(proximity), c(FALSE, TRUE, "inbag", "oob", "all"))
  var.used <- match.arg(as.character(var.used), c("FALSE", "all.trees", "by.tree"))
  split.depth <- match.arg(as.character(split.depth),  c("FALSE", "all.trees", "by.tree"))
  if (var.used == "FALSE") var.used <- FALSE
  if (split.depth == "FALSE") split.depth <- FALSE
  if (missing(formula) | (!missing(formula) && is.null(formula))) {
    formula <- as.formula("Unsupervised() ~ .")
  }
  if (missing(data)) stop("data is missing")
  formulaPrelim <- parseFormula(formula, data, NULL, coerce.factor)
  coerce.factor <- formulaPrelim$coerce.factor
  if (any(is.na(data))) {
    data <- parseMissingData(formulaPrelim, data)
    miss.flag <- TRUE
  }
    else {
      miss.flag <- FALSE
    }
  formulaDetail <- finalizeFormula(formulaPrelim, data)
  ntree <- round(ntree)
  if (ntree < 1) stop("Invalid choice of 'ntree'.  Cannot be less than 1.")
  if (!is.null(nodesize) && nodesize < 1) stop("Invalid choice of 'nodesize'. Cannot be less than 1.")
  if (!is.null(nodedepth)) nodedepth = round(nodedepth) else nodedepth = -1
  nimpute <- round(nimpute)
  if (nimpute < 1) stop("Invalid choice of 'nimpute'.  Cannot be less than 1.")
  seed <- get.seed(seed)
  family <- formulaDetail$family
  xvar.names <- formulaDetail$xvar.names
  yvar.names <- formulaDetail$yvar.names
  if (length(xvar.names) == 0) {
    stop("something seems wrong: your formula did not define any x-variables")
  }
  if (family != "unsupv" && length(yvar.names) == 0) {
    stop("something seems wrong: your formula did not define any y-variables")
  }
  if (family == "class") {
    if (length(setdiff(levels(data[, yvar.names]), unique(data[, yvar.names]))) > 0) {
      warning("empty classes found when implementing classification\n")
    }
  }
  data <- rm.na.levels(data, xvar.names)
  data <- rm.na.levels(data, yvar.names)
  yfactor <- extract.factor(data, yvar.names)
  xfactor <- extract.factor(data, xvar.names)
  yvar.types <- get.yvar.type(family, yfactor$generic.types, yvar.names, coerce.factor)
  yvar.nlevels <- get.yvar.nlevels(family, yfactor$nlevels, yvar.names, data, coerce.factor)
  xvar.types <- get.xvar.type(xfactor$generic.types, xvar.names, coerce.factor)
  xvar.nlevels <- get.xvar.nlevels(xfactor$nlevels, xvar.names, data, coerce.factor)
  data <- finalizeData(c(yvar.names, xvar.names), data, na.action, miss.flag)
  data.row.names <- rownames(data)
  data.col.names <- colnames(data)
  xvar <- as.matrix(data[, xvar.names, drop = FALSE])
  rownames(xvar) <- colnames(xvar) <- NULL
  n <- nrow(xvar)
  n.xvar <- ncol(xvar)
  mtry <- get.grow.mtry(mtry, n.xvar, family)
  samptype <- match.arg(samptype, c("swr", "swor"))
  if (bootstrap == "by.root") {
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
      }
      if (samptype == "swor") {
        sampsize <- min(sampsize, nrow(xvar))
      }
    }
    samp = NULL
  }
    else if (bootstrap == "by.user") {
      samptype <- "swr"
      if (is.null(samp)) {
        stop("samp must not be NULL when bootstrapping by user")
      }
      sampsize <- apply(samp, 2, sum)
      if (sum(sampsize == sampsize[1]) != ntree) {
        stop("sampsize must identical for each tree")
      }
      sampsize <- sampsize[1]
    }
      else {
        sampsize = nrow(xvar)
        samptype <- "swr"
      }
  case.wt  <- get.weight(case.wt, n)
  xvar.wt  <- get.weight(xvar.wt, n.xvar)
  yvar <- as.matrix(data[, yvar.names, drop = FALSE])
  if(dim(yvar)[2] == 0) {
    yvar <- NULL
  }
  if (miss.flag) {
    n.miss <- get.nmiss(xvar, yvar)
  }
    else {
      n.miss <- 0
    }
  if (impute.only && n.miss == 0) {
    return(data)
  }
  remove(data)
  big.data <- FALSE
  event.info <- get.grow.event.info(yvar, family, ntime = ntime)
  splitinfo <- get.grow.splitinfo(formulaDetail, splitrule, nsplit, event.info$event.type)
  if (family == "surv") {
    if (length(event.info$event.type) > 1) {
      if (missing(cause)) {
        cause.wt <- rep(1, length(event.info$event.type))
      }
        else {
          if (length(cause) == 1) {
            if (cause >= 1 && cause <= length(event.info$event.type)) {
              cause.wt <- rep(0, length(event.info$event.type))
              cause.wt[cause] <- 1
            }
              else {
                cause.wt <- rep(1, length(event.info$event.type))
              }
          }
            else {
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
        cause.wt = 1
      }
    family <- get.coerced.survival.fmly(family, event.info$event.type, splitinfo$name)
  }
    else {
      cause.wt <- NULL
    }
  nodesize <- get.grow.nodesize(family, nodesize)
  perf <- NULL
  if ((bootstrap != "by.root") && (bootstrap != "by.user")) {
    importance <- "none"
    perf <- FALSE
  }
  if (family == "unsupv") {
    importance <- "none"
    perf <- FALSE
  }
  if (impute.only) {
    forest       <- FALSE
    proximity    <- FALSE
    var.used     <- FALSE
    split.depth  <- FALSE
    membership   <- FALSE
    perf         <- FALSE
    importance   <- "none"
    terminal.qualts <- FALSE
    terminal.quants <- FALSE
  }
  if (terminal.qualts | terminal.quants) {
    forest <- TRUE
  }
  impute.only.bits <- get.impute.only(impute.only, n.miss)
  var.used.bits <- get.var.used(var.used)
  split.depth.bits <- get.split.depth(split.depth)
  importance.bits <- get.importance(importance)
  bootstrap.bits <- get.bootstrap(bootstrap)
  forest.bits <- get.forest(forest)
  proximity.bits <- get.proximity(TRUE, proximity)
  split.null.bits <- get.split.null(split.null)
  membership.bits <-  get.membership(membership)
  statistics.bits <- get.statistics(statistics)
  split.cust.bits <- get.split.cust(splitinfo$cust)
  perf.flag <- get.perf(perf, impute.only, family)
  perf.bits <-  get.perf.bits(perf.flag)
  samptype.bits <- get.samptype(samptype)
  na.action.bits <- get.na.action(na.action)
  tree.err.bits <- get.tree.err(tree.err)
  terminal.qualts.bits <- get.terminal.qualts(terminal.qualts, FALSE)
  terminal.quants.bits <- get.terminal.quants(terminal.quants, FALSE)
  do.trace <- get.trace(do.trace)
  nativeOutput <- tryCatch({.Call("rfsrcGrow",
                                  as.integer(do.trace),
                                  as.integer(seed),
                                  as.integer(impute.only.bits +
                                               var.used.bits +
                                                 split.depth.bits +
                                                   importance.bits +
                                                     bootstrap.bits +
                                                       forest.bits +
                                                         proximity.bits +
                                                           split.null.bits +
                                                             perf.bits +
                                                               statistics.bits),
                                  as.integer(samptype.bits +
                                                   na.action.bits +
                                                     split.cust.bits +
                                                       tree.err.bits +
                                                         membership.bits +
                                                           terminal.qualts.bits +
                                                             terminal.quants.bits),
                                  as.integer(splitinfo$index),
                                  as.integer(splitinfo$nsplit),
                                  as.integer(mtry),
                                  as.integer(formulaDetail$ytry),
                                  as.integer(nodesize),
                                  as.integer(nodedepth),
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
                                  as.double(xvar.wt),
                                  as.double(xvar),
                                  as.integer(length(event.info$time.interest)),
                                  as.double(event.info$time.interest),
                                  as.double(miss.tree),
                                  as.integer(nimpute),
                                  as.integer(get.rf.cores()))}, error = function(e) {
                                    print(e)
                                    NULL})
  if (is.null(nativeOutput)) {
    if (impute.only) {
      return(NULL)
    }
      else {
        stop("An error has occurred in the grow algorithm.  Please turn trace on for further analysis.")
      }
  }
  if (n.miss > 0) {
    imputed.data <- matrix(nativeOutput$imputation, nrow = n.miss, byrow = FALSE)
    imputed.indv <- imputed.data[, 1]
    imputed.data <- as.matrix(imputed.data[, -1, drop = FALSE])
    nativeOutput$imputation <- NULL
    if (nimpute > 1) {
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
      imputed.indv <- NULL
      imputed.data <- NULL
      imputedOOBData <- NULL
      na.action = "na.omit"
    }
      else {
        colnames(imputed.data) <- c(yvar.names, xvar.names)
        imputed.data <- as.data.frame(imputed.data)
      }
  }
  xvar <- as.data.frame(xvar)
  rownames(xvar) <- data.row.names
  colnames(xvar) <- xvar.names
  xvar <- map.factor(xvar, xfactor)
  if (family != "unsupv") {
    yvar <- as.data.frame(yvar)
    colnames(yvar) <- yvar.names
  }
  else {
    yvar <- NULL
  }
  if (family != "unsupv") {
    if (family == "regr+" | family == "class+" | family == "mix+") {
      yvar <- map.factor(yvar, yfactor)
    }
    else {
      yvar <- amatrix.remove.names(map.factor(yvar, yfactor))
    }
  }
  if ((n.miss > 0) & (nimpute < 2)) {
    imputed.data <- map.factor(imputed.data, xfactor)
    if (family != "unsupv") {
      imputed.data <- map.factor(imputed.data, yfactor)
    }
  }
  if (forest) {
    nativeArraySize = 0
    mwcpPTSize = 0
    for (b in 1:ntree) {
      if (nativeOutput$leafCount[b] > 0) {
        nativeArraySize = nativeArraySize + (2 * nativeOutput$leafCount[b]) - 1
        mwcpPTSize = mwcpPTSize + nativeOutput$mwcpCount[b]
      }
        else {
          nativeArraySize = nativeArraySize + 1
        }
    }
    nativeArray <- as.data.frame(cbind(nativeOutput$treeID[1:nativeArraySize],
                                       nativeOutput$nodeID[1:nativeArraySize],
                                       nativeOutput$parmID[1:nativeArraySize],
                                       nativeOutput$contPT[1:nativeArraySize],
                                       nativeOutput$mwcpSZ[1:nativeArraySize]))
    names(nativeArray) <- c("treeID", "nodeID", "parmID", "contPT", "mwcpSZ")
    if (mwcpPTSize > 0) {
      nativeFactorArray <- nativeOutput$mwcpPT[1:mwcpPTSize]
    }
      else {
        nativeFactorArray <- NULL
      }
    if (terminal.qualts | terminal.quants) {
      temp <- 2 * (nodesize - 1)
      if (sampsize  > temp) { 
        treeTheoreticalMaximum <- sampsize - temp;
      }
        else {
          treeTheoreticalMaximum <- 1;
        }
      forestTheoreticalMaximum <- treeTheoreticalMaximum * ntree
      offset <- 0
      valid.mcnt.indices <- NULL
      for (b in 1:ntree) {
        valid.mcnt.indices = c(valid.mcnt.indices, (offset + 1):(offset + nativeOutput$leafCount[b]))
        offset <- offset + treeTheoreticalMaximum
      }
      if (terminal.quants) {
        if (grepl("surv", family)) {
          offset <- 0
          valid.2D.surv.indices <- NULL
          for (b in 1:ntree) {
            valid.2D.surv.indices = c(valid.2D.surv.indices, (offset + 1):(offset + nativeOutput$leafCount[b] * length(event.info$event.type) * length(event.info$time.interest)))
            offset <- offset + (treeTheoreticalMaximum * length(event.info$event.type) * length(event.info$time.interest))
          }
          offset <- 0
          valid.1D.surv.indices <- NULL
          for (b in 1:ntree) {
            valid.1D.surv.indices = c(valid.1D.surv.indices, (offset + 1):(offset + nativeOutput$leafCount[b] * length(event.info$time.interest)))
            offset <- offset + (treeTheoreticalMaximum * length(event.info$time.interest))
          }
          offset <- 0
          valid.mort.indices <- NULL
          for (b in 1:ntree) {
            valid.mort.indices = c(valid.mort.indices, (offset + 1):(offset + nativeOutput$leafCount[b] * length(event.info$event.type)))
            offset <- offset + (treeTheoreticalMaximum * length(event.info$event.type))
          }
        }
          else {
            class.index <- which(yvar.types != "R")
            class.count <- length(class.index)
            regr.index <- which(yvar.types == "R")
            regr.count <- length(regr.index)
            if (class.count > 0) {
              levels.count <- array(0, class.count)
              counter <- 0
              for (i in class.index) {
                counter <- counter + 1
                levels.count[counter] <- yvar.nlevels[i]
              }
              offset <- 0
              valid.clas.indices <- NULL
              for (b in 1:ntree) {
                valid.clas.indices = c(valid.clas.indices, (offset + 1):(offset + nativeOutput$leafCount[b] * sum(levels.count)))
                offset <- offset + (treeTheoreticalMaximum * sum(levels.count))
              }
            }
            if (regr.count > 0) {
              offset <- 0
              valid.regr.indices <- NULL
              for (b in 1:ntree) {
                valid.regr.indices = c(valid.regr.indices, (offset + 1):(offset + nativeOutput$leafCount[b] * regr.count))
                offset <- offset + (treeTheoreticalMaximum * regr.count)
              }
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
    forest.out <- list(nativeArray = nativeArray,
                       nativeFactorArray = nativeFactorArray,
                       totalNodeCount = dim(nativeArray)[1],
                       nodesize = nodesize,
                       nodedepth = nodedepth,
                       split.null = split.null,
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
                       version = "2.4.2",
                       na.action = na.action,
                       coerce.factor = coerce.factor)
    if (grepl("surv", family)) {
      forest.out$time.interest <- event.info$time.interest
    }
    class(forest.out) <- c("rfsrc", "forest", family)
    if (big.data) {
      class(forest.out) <- c(class(forest.out), "bigdata")
    }
  }
  else {
    forest.out <- NULL
  }
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
  if (statistics) {
    node.stats <- as.data.frame(cbind(nativeOutput$spltST[1:nativeArraySize]))
    colnames(node.stats) <- "spltST"
    node.mtry.stats <- t(array(nativeOutput$mtryST[1:nativeArraySize], c(mtry, forest.out$totalNodeCount)))
    node.mtry.index <- t(array(nativeOutput$mtryID[1:nativeArraySize], c(mtry, forest.out$totalNodeCount)))
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
  rfsrcOutput <- list(
    call = match.call(),
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
    leaf.count = nativeOutput$leafCount,
    proximity = proximity.out,
    forest = forest.out,
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
    tree.err = tree.err
  )
  if (is.null(coerce.factor$yvar.names)) remove(yvar)
  remove(xvar)
  nativeOutput$leafCount <- NULL
  remove(proximity.out)
  remove(forest.out)
  remove(membership.out)
  remove(inbag.out)
  remove(var.used.out)
  if (n.miss > 0) remove(imputed.indv)
  if (n.miss > 0) remove(imputed.data)
  remove(split.depth.out)
  survOutput <- NULL
  classOutput <- NULL
  regrOutput <- NULL
  if (!impute.only) {
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
        ens.names <- list(NULL, NULL)
        mortality.names <- list(NULL, NULL)
        err.names <- list(NULL, NULL)
        vimp.names <- list(NULL, xvar.names)
      }
        else {
          ens.names <- list(NULL, NULL, c(paste("condCHF.", 1:length(event.info$event.type), sep = "")))
          mortality.names <- list(NULL, paste("event.", 1:length(event.info$event.type), sep = ""))
          cif.names <- list(NULL, NULL, c(paste("CIF.", 1:length(event.info$event.type), sep = "")))
          err.names <- list(c(paste("event.", 1:length(event.info$event.type), sep = "")), NULL)
          vimp.names <- list(paste("event.", 1:length(event.info$event.type), sep = ""), xvar.names)
        }
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
        if (!is.null(coerce.factor$yvar.names)) remove(yvar)
        tree.offset <- array(1, ntree)
        if (ntree > 1) {
          tree.offset[2:ntree] <- sum(1 + levels.count)
        }
        tree.offset <-  cumsum(tree.offset)
        vimp.offset <- array(1, n.xvar)
        if (n.xvar > 1) {
          vimp.offset[2:n.xvar] <- sum(1 + levels.count)
        }
        vimp.offset <-  cumsum(vimp.offset)
        iter.ensb.start <- 0
        iter.ensb.end   <- 0
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
          response <- (if (!is.null(predicted)) bayes.rule(predicted) else NULL)
          classOutput[[i]] <- c(classOutput[[i]], class = list(response))
          remove(predicted)
          remove(response)
          predicted.oob <- (if (!is.null(nativeOutput$oobEnsbCLS))
                            array(nativeOutput$oobEnsbCLS[(iter.ensb.start + 1):iter.ensb.end],
                                  c(n, levels.count[i]), dimnames=ens.names) else NULL)
          classOutput[[i]] <- c(classOutput[[i]], predicted.oob = list(predicted.oob))
          response.oob <- (if (!is.null(predicted.oob)) bayes.rule(predicted.oob) else NULL)
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
        vimp.offset <- array(1, n.xvar)
        if (n.xvar > 1) {
          vimp.offset[2:n.xvar] <- length(regr.index)
        }
        vimp.offset <-  cumsum(vimp.offset)
        iter.ensb.start <- 0
        iter.ensb.end   <- 0
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
  }
  class(rfsrcOutput) <- c("rfsrc", "grow", family)
  if (big.data) {
    class(rfsrcOutput) <- c(class(rfsrcOutput), "bigdata")
  }
  return(rfsrcOutput)
}
