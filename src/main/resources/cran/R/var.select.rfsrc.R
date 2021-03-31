var.select.rfsrc <-
  function(formula,          
           data,
           object,
           cause,
           m.target = NULL,
           method = c("md", "vh", "vh.vimp"),
           conservative = c("medium", "low", "high"),
           ntree = (if (method == "md") 1000 else 500),
           mvars = (if (method != "md") ceiling(ncol(data)/5) else NULL),
           mtry = (if (method == "md") ceiling(ncol(data)/3) else NULL),
           nodesize = 2,
           splitrule = NULL,
           nsplit = 10,
           xvar.wt = NULL,
           refit = (method != "md"),
           fast = FALSE,
           na.action = c("na.omit", "na.impute"),
           always.use = NULL,  
           nrep = 50,        
           K = 5,             
           nstep = 1,         
           prefit =  list(action = (method != "md"), ntree = 100, mtry = 500, nodesize = 3, nsplit = 1),
           verbose = TRUE,
           ...
           )
{
  ## --------------------------------------------------------------
  ##  
  ##  workhorse: variable hunting algorithm
  ##
  ## --------------------------------------------------------------
  rfsrc.var.hunting <- function(train.id, var.pt, nstep) {
    ## ------------------filtering step-----------------------
    if (verbose) cat("\t", paste("selecting variables using", mName), "...\n")
    ## which variables to include
    drop.var.pt <- setdiff(var.columns, var.pt)
    ## family specific checks
    if (grepl("surv", family)) {
      if (sum(data[train.id, 2], na.rm = TRUE) < 2) {
        stop("training data has insufficient deaths: K is probably set too high\n")
      }
    }
    ## filtered forest
    ## over-ride user mtry setting: use an aggressive value
    rfsrc.filter.obj  <- rfsrc(rfsrc.all.f,
                               data=(if (LENGTH(var.pt, drop.var.pt)) data[train.id, -drop.var.pt]
                                       else data[train.id, ]),
                               ntree = ntree,
                               splitrule = splitrule,
                               nsplit = nsplit,
                               mtry = Mtry(var.columns, drop.var.pt),
                               nodesize = nodesize,
                               cause = cause,
                               na.action = na.action,
                               importance=TRUE,
                               block.size = block.size, perf.type = perf.type)
    ## set the target dimension for CR families
    if (rfsrc.filter.obj$family == "surv-CR") {
      target.dim <- max(1, min(cause, max(get.event.info(rfsrc.filter.obj)$event.type)), na.rm = TRUE)
    }
    ## extract vimp
    ## for multivariate families we must manually extract the importance and error rate
    imp <- get.varselect.imp(coerce.multivariate(rfsrc.filter.obj, m.target), target.dim)
    names(imp) <- rfsrc.filter.obj$xvar.names
    ## selection using vimp
    if (method == "vh.vimp") {
      VarStrength <- sort(imp, decreasing = TRUE)
      lower.VarStrength <- min(VarStrength) - 1 #need theoretical lower bound to vimp
      n.lower <- min(2, length(VarStrength))    #n.lower cannot be > no. available variables
      forest.depth <- m.depth <- NA
      sig.vars.old <- names(VarStrength)[1]
    }
    ## selection using minimal depth
      else {
        max.obj <- max.subtree(rfsrc.filter.obj, conservative = (conservative == "high"))
        if (is.null(max.obj$order)) {
          ## maximal information failed; revert to vimp
          VarStrength <- lower.VarStrength <- 0
          forest.depth <- m.depth <- NA
          sig.vars.old <- names(sort(imp, decreasing = TRUE))[1]
        }
          else {
            m.depth <- VarStrength <- max.obj$order[, 1]
            forest.depth <- floor(mean(apply(max.obj$nodes.at.depth, 2, function(x){sum(!is.na(x))}), na.rm=TRUE))
            exact.threshold <- ifelse(conservative == "low", max.obj$threshold.1se, max.obj$threshold)
            n.lower <- max(min(2, length(VarStrength)),#n.lower cannot be > no. available variables
                           sum(VarStrength <= exact.threshold))
            VarStrength <- max(VarStrength) - VarStrength
            names(m.depth) <- names(VarStrength) <- rfsrc.filter.obj$xvar.names
            VarStrength <- sort(VarStrength, decreasing = TRUE)
            lower.VarStrength <- -1 #need theoretical upper bound to first order statistic
            sig.vars.old <- names(VarStrength)[1]
          }
      }
    ## set nstep
    nstep <- max(round(length(rfsrc.filter.obj$xvar.names)/nstep), 1)
    imp.old <- 0
    ## regularized forward selection using joint vimp
    for (b in 1:nstep) {
      if (b == 1) {
        if (sum(VarStrength > lower.VarStrength) == 0) {
          sig.vars <- sig.vars.old
          break
        }
        n.upper <- max(which(VarStrength > lower.VarStrength), n.lower)
        threshold <- unique(round(seq(n.lower, n.upper, length = nstep)))
        if (length(threshold) < nstep) {
          threshold <- c(threshold, rep(max(threshold), nstep - length(threshold)))
        }
      }
      sig.vars <- names(VarStrength)[1:threshold[b]]
      if (!is.null(always.use)) {
        sig.vars <- unique(c(sig.vars, always.use))
      }
      if (length(sig.vars) <= 1) {#break if there is only one variable
        sig.vars <- sig.vars.old
        break
      }
      imp <- coerce.multivariate(vimp(rfsrc.filter.obj, sig.vars, m.target = m.target,
                                      joint = TRUE), m.target)$importance[target.dim]
      ## verbose output
      if (verbose) cat("\t iteration: ", b,
                       "  # vars:",     length(sig.vars),
                       "  joint-vimp:",  round(imp, 3),
                       "\r")
      ## break when joint vimp no longer increases (strict inequality is a safety feature)
      if (imp  <= imp.old) {
        sig.vars <- sig.vars.old
        break
      }
        else {
          var.pt <- var.columns[match(sig.vars, xvar.names)]
          sig.vars.old <- sig.vars
          imp.old <- imp
        }
    }
    ## refit forest and exit
    ## over-ride user specified nodesize: use default mtry as variables 
    ## should now be filtered and default settings should therefore apply
    var.pt <- var.columns[match(sig.vars, xvar.names)]
    drop.var.pt <- setdiff(var.columns, var.pt)
    rfsrc.obj  <- rfsrc(rfsrc.all.f,
                        data=(if (LENGTH(var.pt, drop.var.pt)) data[train.id, -drop.var.pt]
                                else data[train.id, ]),
                        ntree = ntree,
                        splitrule = splitrule,
                        nsplit = nsplit,
                        cause = cause,
                        na.action = na.action,
                        perf.type = perf.type)
    return(list(rfsrc.obj=rfsrc.obj, sig.vars=rfsrc.obj$xvar.names, forest.depth=forest.depth, m.depth=m.depth))
  }
  ## --------------------------------------------------------------
  ##   
  ##   preliminary processing
  ##
  ## --------------------------------------------------------------
  ## Incoming parameter checks.  All are fatal.
  if (!missing(object)) {
    if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2) {
      stop("This function only works for objects of class `(rfsrc, grow)'")
    }
    if (inherits(object, "anonymous")) {
      stop("this function does work with anonymous forests")
    }
    if (is.null(object$forest)) {
      stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'")
    }
    rfsrc.all.f <- object$formula
  }
    else {
      if (missing(formula) || missing(data)) {
        ## allowance for users who overlook the correct way to assign the object
        if (sum(inherits(formula, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
          object <- formula
        }
          else {
            stop("Need to specify 'formula' and 'data' or provide a grow forest object")
          }
      }
      rfsrc.all.f <- formula
    }
  ## If an object is provided, after the above checks, make the m.target coherent.
  if (!missing(object)) {
    ## initialize the target outcome, in case it is NULL
    ## coersion of an object depends on a target outcome
    m.target <- get.univariate.target(object, m.target)
  }
  ## rearrange the data
  ## need to handle unsupervised families carefully: minimal depth is the only permissible method
  ## pull performance parameters 
  if (missing(object)) {
    ## parse the formula
    formulaDetail <- finalizeFormula(parseFormula(rfsrc.all.f, data), data)
    family <- formulaDetail$family
    xvar.names <- formulaDetail$xvar.names
    yvar.names <- formulaDetail$yvar.names
    if (family != "unsupv") {
      data <- cbind(data[, yvar.names, drop = FALSE], data[, match(xvar.names, names(data))])
      yvar <- data[, yvar.names]
      yvar.dim <- ncol(cbind(yvar))
    }
      else {
        data <- data[, match(xvar.names, names(data))]
        yvar.dim <- 0
        method <- "md"
    }
    dots <- list(...)
    block.size <- dots$block.size
    perf.type <- is.hidden.perf.type(dots)
  }
    else {
      ## parse the object
      family <- object$family
      xvar.names <- object$xvar.names
      if (family != "unsupv") {
        yvar.names <- object$yvar.names
        data <- data.frame(object$yvar, object$xvar)
        colnames(data) <- c(yvar.names, xvar.names)
        yvar <- data[, yvar.names]
        yvar.dim <- ncol(cbind(yvar))
      }
        else {
          data <- object$xvar
          yvar.dim <- 0
          method <- "md"
      }
    block.size <- object$block.size
    perf.type <- object$forest$perf.type
    }
  ## specify the default event type for CR
  if (missing(cause)) {
    cause <- 1
  }
  ## verify key options
  method <- match.arg(method, c("md", "vh", "vh.vimp"))
  conservative = match.arg(conservative, c("medium", "low", "high"))
  ## pretty names for method
  mName <- switch(method,
                  "md"      = "Minimal Depth",
                  "vh"      = "Variable Hunting",
                  "vh.vimp" = "Variable Hunting (VIMP)")
  ## simplify the formula: needed when we drop variables
  rfsrc.all.f <- switch(family,
                        "surv"   = as.formula(paste("Surv(",yvar.names[1],",",yvar.names[2],") ~ .")),
                        "surv-CR"= as.formula(paste("Surv(",yvar.names[1],",",yvar.names[2],") ~ .")),
                        "regr"   = as.formula(paste(yvar.names, "~ .")),
                        "class"  = as.formula(paste(yvar.names, "~ .")),
                        "unsupv" = NULL,
                        "regr+"  = as.formula(paste("Multivar(", paste(yvar.names, collapse = ","), paste(") ~ ."), sep = "")),
                        "class+" = as.formula(paste("Multivar(", paste(yvar.names, collapse = ","), paste(") ~ ."), sep = "")),
                        "mix+"   = as.formula(paste("Multivar(", paste(yvar.names, collapse = ","), paste(") ~ ."), sep = ""))
                        )
  ## initialize dimensions
  n <- nrow(data)
  P <- length(xvar.names)
  ## Specify the target event.  Later will be over-written for CR
  target.dim <- 1
  ## make special allowance for always.use x-variables
  var.columns <- (1 + yvar.dim):ncol(data)
  if (!is.null(always.use)) {
    always.use.pt <- var.columns[match(always.use, xvar.names)]
  }
    else {
      always.use.pt <- NULL
    }
  ## if xvar weight is specified, then make sure it is defined correctly
  xvar.wt <- get.weight(xvar.wt, P)
  ## Final checks on option parameters  
  if (!is.null(mtry)) {
    mtry <- round(mtry)
    if (mtry < 1 | mtry > P) mtry <- max(1, min(mtry, P))
  }
  ## prefit forest parameter details
  prefit.masterlist <- list(action = (method != "md"), ntree = 100, mtry = 500, nodesize = 3, nsplit = 1)
  parm.match <- na.omit(match(names(prefit), names(prefit.masterlist)))
  if (length(parm.match) > 0) {
    for (l in 1:length(parm.match)) {
      prefit.masterlist[[parm.match[l]]] <- prefit[[l]]
    }
  }
  prefit <- prefit.masterlist
  prefit.flag  <- prefit$action
  ## --------------------------------------------------------------
  ##  
  ##  minimal depth analysis
  ##
  ## --------------------------------------------------------------
  if (method == "md") {
    ## ------------------------------------------------
    ## run preliminary forest to determine weights for variables
    ## we do this OUTSIDE of the loop
    if (prefit.flag && is.null(xvar.wt) && missing(object)) {
      if (verbose) cat("Using forests to preweight each variable's chance of splitting a node...\n")
      rfsrc.prefit.obj  <- rfsrc(rfsrc.all.f,
                                 data = data,
                                 ntree = prefit$ntree,
                                 nodesize = prefit$nodesize,
                                 mtry = prefit$mtry,
                                 splitrule = splitrule,
                                 nsplit = prefit$nsplit,
                                 cause = cause,
                                 na.action = na.action,
                                 importance = TRUE,
                                 block.size = block.size, perf.type = perf.type)
      ## set the target dimension for CR families
      if (rfsrc.prefit.obj$family == "surv-CR") {
        target.dim <- max(1, min(cause, max(get.event.info(rfsrc.prefit.obj)$event.type)), na.rm = TRUE)
      }
      ## for multivariate families we must manually extract the importance and error rate
      wts <- pmax(get.varselect.imp(coerce.multivariate(rfsrc.prefit.obj, m.target), target.dim), 0)
      if (any(wts > 0)) {
        xvar.wt <- get.weight(wts, P)
      }
      rm(rfsrc.prefit.obj)
    }
    ## ------------------------------------------------
    ## extract minimal depth
    if (!missing(object)) {
    }
    if (!missing(object) && !prefit.flag) {
      if (verbose) cat("minimal depth variable selection ...\n")
      md.obj <- max.subtree(object, conservative = (conservative == "high"))
      ## for multivariate families we must manually extract the importance and error rate
      object <- coerce.multivariate(object, m.target)
      m.target <- object$outcome.target
      pe <- get.varselect.err(object)
      ntree <- object$ntree
      nsplit <- object$nsplit
      mtry <- object$mtry
      nodesize <- object$nodesize
      ## set the target dimension for CR families
      if (family == "surv-CR") {
        target.dim <- max(1, min(cause, max(get.event.info(object)$event.type)), na.rm = TRUE)
      }
      imp <- get.varselect.imp(object, target.dim)
      imp.all <- get.varselect.imp.all(object)
      rm(object)
    }
    ## ------------------------------------------------
    ## otherwise run forests...then extract minimal depth    
      else {
        if (verbose) cat("running forests ...\n")
        rfsrc.obj <- rfsrc(rfsrc.all.f,
                           data,
                           ntree = ntree,
                           mtry = mtry,
                           nodesize = nodesize,
                           splitrule = splitrule,
                           nsplit = nsplit,
                           cause = cause,
                           na.action = na.action,
                           xvar.wt = xvar.wt,
                           importance = TRUE,
                           block.size = block.size, perf.type = perf.type)
        ## set the target dimension for CR families
        if (rfsrc.obj$family == "surv-CR") {
          target.dim <- max(1, min(cause, max(get.event.info(rfsrc.obj)$event.type)), na.rm = TRUE)
        }
        if (verbose) cat("minimal depth variable selection ...\n")
        md.obj <- max.subtree(rfsrc.obj, conservative = (conservative == "high"))
        ## for multivariate families we must manually extract the importance and error rate
        rfsrc.obj <- coerce.multivariate(rfsrc.obj, m.target)
        m.target <- rfsrc.obj$outcome.target
        pe <- get.varselect.err(rfsrc.obj)
        imp <- get.varselect.imp(rfsrc.obj, target.dim)
        imp.all <- get.varselect.imp.all(rfsrc.obj)
        mtry <- rfsrc.obj$mtry
        nodesize <- rfsrc.obj$nodesize
        n <- nrow(rfsrc.obj$xvar)
        family <- rfsrc.obj$family
        rm(rfsrc.obj)#don't need the grow object
      }
    ## parse minimal depth information
    depth <- md.obj$order[, 1]
    threshold <- ifelse(conservative == "low", md.obj$threshold.1se, md.obj$threshold)
    top.var.pt <- (depth <= threshold)
    modelsize <- sum(top.var.pt)
    o.r.m <- order(depth, decreasing = FALSE)
    top.var.pt <- top.var.pt[o.r.m]
    varselect <- as.data.frame(cbind(depth = depth, vimp = imp.all))[o.r.m, ]
    topvars <- unique(c(always.use, rownames(varselect)[top.var.pt]))
    ## fit a forest to final variable list
    ## use default settings for nodesize, mtry due to dimension reduction
    if (refit == TRUE) {
      if (verbose) cat("fitting forests to minimal depth selected variables ...\n")
      var.pt <- var.columns[match(topvars, xvar.names)]
      var.pt <- unique(c(var.pt, always.use.pt))
      drop.var.pt <- setdiff(var.columns, var.pt)
      rfsrc.refit.obj  <- rfsrc(rfsrc.all.f,
                                data=(if (LENGTH(var.pt, drop.var.pt)) data[, -drop.var.pt, drop = FALSE] else data),
                                ntree = ntree,
                                splitrule = splitrule,
                                nsplit = nsplit,
                                na.action = na.action,
                                perf.type = perf.type)
      ## for multivariate families we must manually extract the importance and error rate
      rfsrc.refit.obj <- coerce.multivariate(rfsrc.refit.obj, m.target)
    }
      else {
        rfsrc.refit.obj <- NULL
      }
    ## output: all nicely packaged
    if (verbose) {
      cat("\n\n")
      cat("-----------------------------------------------------------\n")
      cat("family             :", family, "\n")
      if (family == "regr+" | family == "class+" | family == "mix+") {
        cat("no. y-variables    : ", yvar.dim,       "\n", sep="")
        cat("response used      : ", m.target, "\n", sep="")
      }    
      cat("var. selection     :", mName, "\n")
      cat("conservativeness   :", conservative, "\n")
      cat("x-weighting used?  :", !is.null(xvar.wt), "\n")
      cat("dimension          :", P, "\n")
      cat("sample size        :", n, "\n")
      cat("ntree              :", ntree, "\n")
      cat("nsplit             :", nsplit, "\n")
      cat("mtry               :", mtry, "\n")
      cat("nodesize           :", nodesize, "\n")
      cat("refitted forest    :", refit, "\n")
      cat("model size         :", modelsize, "\n")
      cat("depth threshold    :", round(threshold, 4), "\n")
      if (!prefit.flag) {
        cat("PE (true OOB)      :", round(pe, 4), "\n")
      }
        else {
          cat("PE (biased)        :", round(pe, 4), "\n")
        }
      cat("\n\n")
      cat("Top variables:\n")
      print(round(varselect[top.var.pt, ], 3))
      cat("-----------------------------------------------------------\n")
    }
    ## Return the goodies
    return(invisible((list(err.rate=pe,
                           modelsize=modelsize,
                           topvars=topvars,
                           varselect=varselect,
                           rfsrc.refit.obj=rfsrc.refit.obj,
                           md.obj=md.obj
                           ))))
  }  
  ## --------------------------------------------------------------
  ##  
  ##  VH algorithm
  ##
  ## --------------------------------------------------------------
  ## vectors/matrices etc.
  pred.results <- dim.results <- forest.depth <- rep(0, nrep)
  var.signature <- NULL
  var.depth <- matrix(NA, nrep, P)
  ## ------------------------------------------------
  ## run preliminary forest to determine weights for variables
  ## we do this OUTSIDE of the loop
  outside.loop <- FALSE
  if (prefit.flag & is.null(xvar.wt)) {
    if (verbose) cat("Using forests to select a variables likelihood of splitting a node...\n")
    rfsrc.prefit.obj  <- rfsrc(rfsrc.all.f,
                               data = data,
                               ntree = prefit$ntree,
                               mtry = prefit$mtry,
                               nodesize = prefit$nodesize,
                               nsplit = prefit$nsplit,
                               cause = cause,
                               splitrule = splitrule,
                               na.action = na.action,
                               perf.type = perf.type)
    ## set the target dimension for CR families
    if (rfsrc.prefit.obj$family == "surv-CR") {
      target.dim <- max(1, min(cause, max(get.event.info(rfsrc.prefit.obj)$event.type)), na.rm = TRUE)
    }
    ## record that a pre-fit has occurred
    outside.loop <- TRUE
  }
  ## ------------------------------------------------
  ## loop
  for (m in 1:nrep) {
    if (verbose & nrep>1) cat("---------------------  Iteration:", m, "  ---------------------\n")
    ## train/test subsamples
    ## use balanced sampling for CR/multiclass
    all.folds <- switch(family,
                        "surv"     =  balanced.folds(yvar[, 2], K),
                        "surv-CR"  =  balanced.folds(yvar[, 2], K),
                        "class"    =  balanced.folds(yvar, K),
                        "regr"     =  cv.folds(n, K),
                        "class+"   =  balanced.folds(yvar, K),
                        "regr+"    =  cv.folds(n, K),
                        "mix+"     =  cv.folds(n, K)
                        )    
    if (fast == TRUE) {
      train.id <- all.folds[[1]]
      test.id <- all.folds[[2]]
    }
      else {
        test.id <- all.folds[[1]]
        train.id <- setdiff(1:n, test.id)
      }
    ## run preliminary forest to determine weights for variables
    ## we do this INSIDE of the loop
    if (is.null(xvar.wt)) {
      if (!prefit.flag) {
        if (verbose) cat("Using forests to determine variable selection weights...\n")
        rfsrc.prefit.obj  <- rfsrc(rfsrc.all.f,
                                   data = data[train.id,, drop = FALSE],
                                   ntree = prefit$ntree,
                                   mtry = prefit$mtry,
                                   nodesize = prefit$nodesize,                                    
                                   nsplit = prefit$nsplit,
                                   cause = cause,
                                   splitrule = splitrule,
                                   na.action = na.action,
                                   importance = TRUE,
                                   block.size = block.size, perf.type = perf.type)
        ## set the target dimension for CR families
        if (rfsrc.prefit.obj$family == "surv-CR") {
          target.dim <- max(1, min(cause, max(get.event.info(rfsrc.prefit.obj)$event.type)), na.rm = TRUE)
        }
      }
      ## for multivariate families we must manually extract the importance and error rate
      rfsrc.prefit.obj <- coerce.multivariate(rfsrc.prefit.obj, m.target)
      wts <- pmax(get.varselect.imp(rfsrc.prefit.obj, target.dim), 0)
      if (any(wts > 0)) {
        var.pt <- unique(resample(var.columns, mvars, replace = TRUE, prob = wts))
      }
        else {
          var.pt <- var.columns[1:P]
        }
    }
      else {
        var.pt <- var.columns[1:P]
      }
    ## pre-guided gene selection
    if (!is.null(xvar.wt)) {
      var.pt <- unique(resample(var.columns, mvars, replace = TRUE, prob = xvar.wt))
    }
    ## always.use variables 
    if (!is.null(always.use)) {
      var.pt <- unique(c(var.pt, always.use.pt))
    }
    ## RFSRC gene hunting call
    object <- rfsrc.var.hunting(train.id, var.pt, nstep)
    rfsrc.obj <- object$rfsrc.obj
    m.target <- get.univariate.target(rfsrc.obj, m.target)
    sig.vars <- object$sig.vars
    if (method == "vh") {
      forest.depth[m] <- object$forest.depth
      var.depth[m, match(names(object$m.depth), xvar.names)] <- object$m.depth
    }
    ## RFSRC prediction
    ## for multivariate families we must manually extract the importance and error rate
    pred.out <- coerce.multivariate(predict(rfsrc.obj, data[test.id, ]), m.target)
    pred.results[m] <- get.varselect.err(pred.out)[target.dim] 
    dim.results[m] <- length(sig.vars)
    var.signature <- c(var.signature, sig.vars)
    ## nice output
    if (verbose) {
      cat("\t                                                                \r")
      cat("\t PE:", round(pred.results[m], 4), "     dim:", dim.results[m], "\n")
    }
  }
  ## --------------------------------------------------------------
  ## finalize details
  ## remove NA's in PE 
  pred.results <- c(na.omit(pred.results))
  ## frequency
  var.freq.all.temp <- 100 * tapply(var.signature, var.signature, length) / nrep
  freq.pt <- match(names(var.freq.all.temp), xvar.names)
  var.freq.all <- rep(0, P)
  var.freq.all[freq.pt] <- var.freq.all.temp
  ##  package it up for output
  if (method == "vh") {
    var.depth.all <- apply(var.depth, 2, mean, na.rm = T)
    varselect <- cbind(depth = var.depth.all, rel.freq = var.freq.all)
  }
    else {
      varselect <- cbind(rel.freq = var.freq.all)
    }
  o.r.f <- order(var.freq.all, decreasing = TRUE)
  rownames(varselect) <- xvar.names
  varselect <- varselect[o.r.f,, drop = FALSE]
  modelsize <- ceiling(mean(dim.results))  
  topvars <- unique(c(always.use, rownames(varselect)[1:modelsize]))
  ## fit a forest to final variable list
  if (refit == TRUE) {
    if (verbose) cat("fitting forests to final selected variables ...\n")
    var.pt <- var.columns[match(rownames(varselect)[1:modelsize], xvar.names)]
    drop.var.pt <- setdiff(var.columns, var.pt)
    rfsrc.refit.obj  <- rfsrc(rfsrc.all.f,
                              data = (if (LENGTH(var.pt, drop.var.pt)) data[, -drop.var.pt]
                                        else data),
                              na.action = na.action,
                              ntree = ntree,
                              nodesize = nodesize,
                              nsplit = nsplit,
                              cause = cause,
                              splitrule = splitrule,
                              block.size = block.size, perf.type = perf.type)
  }
    else {
      rfsrc.refit.obj <- NULL
    }
  ## output: all nicely packaged
  if (verbose) {
    cat("\n\n")
    cat("-----------------------------------------------------------\n")
    cat("family             :", family, "\n")
    if (family == "regr+" | family == "class+" | family == "mix+") {
      cat("no. y-variables    : ", yvar.dim,              "\n", sep="")
      cat("response used      : ", m.target, "\n", sep="")
    }    
    cat("var. selection     :", mName, "\n")
    cat("conservativeness   :", conservative, "\n")
    cat("dimension          :", P, "\n")
    cat("sample size        :", n, "\n")
    cat("K-fold             :", K, "\n")
    cat("no. reps           :", nrep, "\n")
    cat("nstep              :", nstep, "\n")
    cat("ntree              :", ntree, "\n")
    cat("nsplit             :", nsplit, "\n")
    cat("mvars              :", mvars, "\n")
    cat("nodesize           :", nodesize, "\n")
    cat("refitted forest    :", refit, "\n")
    if (method == "vh") {
      cat("depth ratio        :", round(mean(mvars/(2^forest.depth)), 4), "\n")
    }
    cat("model size         :", round(mean(dim.results), 4), "+/-", round(SD(dim.results), 4), "\n")
    if (outside.loop) {
      cat("PE (K-fold, biased):", round(mean(pred.results), 4), "+/-", round(SD(pred.results), 4), "\n")
    }
      else {
        cat("PE (K-fold)        :", round(mean(pred.results), 4), "+/-", round(SD(pred.results), 4), "\n")
      }
    cat("\n\n")
    cat("Top variables:\n")
    print(round(varselect[1:modelsize,, drop = FALSE], 3))
    cat("-----------------------------------------------------------\n")
  }
  ## return the goodies 
  return(invisible(list(err.rate=pred.results,
                        modelsize=modelsize,
                        topvars=topvars,
                        varselect=varselect,
                        rfsrc.refit.obj=rfsrc.refit.obj,
                        md.obj=NULL
                        )))
}
## --------------------------------------------------------------
##  
## internal functions
##
## --------------------------------------------------------------
get.varselect.imp <- function(f.o, target.dim) {
  if (!is.null(f.o$importance)) {
    c(cbind(f.o$importance)[, target.dim])
  }
    else {
      rep(NA, length(f.o$xvar.names))
    }
}
get.varselect.imp.all <- function(f.o) {
  if (!is.null(f.o$importance)) {
    imp.all <- cbind(f.o$importance)
    if (ncol(imp.all) == 1) {
      colnames(imp.all) <- "vimp"
    }
      else {
        colnames(imp.all) <- paste("vimp.", colnames(imp.all), sep = "")
      }
    imp.all
  }
    else {
      rep(NA, length(f.o$xvar.names))
    }
}
get.varselect.err <- function(f.o) {
  if (!is.null(f.o$err.rate)) {
    if (grepl("surv", f.o$family)) {
      err <- 100 * cbind(f.o$err.rate)[f.o$ntree, ]
    }
      else {
        err <- cbind(f.o$err.rate)[f.o$ntree, ]
      }
  }
    else {
      err = NA
    }
  err
}
SD <- function(x) {
  if (all(is.na(x))) {
    NA
  }
    else {
      sd(x, na.rm = TRUE)
    }
}
LENGTH <- function(x, y) {
  (length(x) > 0 & length(y) > 0)
}
Mtry <- function(x, y) {
  mtry <- round((length(x) - length(y))/3)
  if (mtry == 0) {
    round(length(x)/3)
  }
    else {
      mtry
    }
}
permute.rows <-function(x) {
  n <- nrow(x)
  p <- ncol(x)
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}
balanced.folds <- function(y, nfolds = min(min(table(y)), 10)) {
  y[is.na(y)] <- resample(y[!is.na(y)], size = sum(is.na(y)), replace = TRUE)
  totals <- table(y)
  if (length(totals) < 2) {
    return(cv.folds(length(y), nfolds))
  }
    else {
      fmax <- max(totals)
      nfolds <- min(nfolds, fmax)     
      nfolds <- max(nfolds, 2)
      folds <- as.list(seq(nfolds))
      yids <- split(seq(y), y) 
      bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
      for(i in seq(totals)) {
        if(length(yids[[i]])>1){bigmat[seq(totals[i]), i] <- sample(yids[[i]])}
        if(length(yids[[i]])==1){bigmat[seq(totals[i]), i] <- yids[[i]]}
      }
      smallmat <- matrix(bigmat, nrow = nfolds)
      smallmat <- permute.rows(t(smallmat)) 
      res <- vector("list", nfolds)
      for(j in 1:nfolds) {
        jj <- !is.na(smallmat[, j])
        res[[j]] <- smallmat[jj, j]
      }
      return(res)
    }
}
var.select <- var.select.rfsrc
