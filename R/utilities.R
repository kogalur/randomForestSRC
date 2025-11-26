get.block.size.bits <- function (block.size, ntree) {
  ## Check for user silliness.
  if (!is.null(block.size)) {
    ## for backwards compatibility allow TRUE/FALSE
    if (is.logical(block.size)) {
      if (block.size) {
        block.size <- 1
      }
      else {
        block.size <- ntree
      }
    }
    else if ((block.size < 1) || (block.size > ntree)) {
      block.size <- ntree
    }
    else {
      block.size <- round(block.size)
    }
  }
  else {
    block.size <- ntree
  }
  return (block.size)
}
get.bootstrap.bits <- function (bootstrap) {
  if (bootstrap == "by.root") {
    bootstrap <- 0
  }
  else if (bootstrap == "by.node") {
    bootstrap <- 2^19
  }
  else if (bootstrap == "none") {
    bootstrap <- 2^20
  }
  else if (bootstrap == "by.user") {
    bootstrap <- 2^19 + 2^20
  }
  else {
    stop("Invalid choice for 'bootstrap' option:  ", bootstrap)
  }
  return (bootstrap)
}
get.case.depth.bits <- function (case.depth) {
  if (!is.null(case.depth)) {
    if (case.depth == TRUE) {
      case.depth <- 2^11
    }
    else if (case.depth == FALSE) {
      case.depth <- 0
    }
    else {
      stop("Invalid choice for 'case.depth' option:  ", case.depth)
    }
  }
  else {
    stop("Invalid choice for 'case.depth' option:  ", case.depth)
  }
  return (case.depth)
}
get.cr.bits <- function (fmly) {
  if (fmly == "surv-CR") {
    return(2^21)
  } else {
    return(0)
  }
}
get.cse.bits <- function (cse) {
  if (!is.null(cse)) {
    if (cse == TRUE) {
      cse <- 2^28
    }
    else if (cse == FALSE) {
      cse <- 0
    }
    else {
      stop("Invalid choice for 'cse' option:  ", cse)
    }
  }
  else {
    stop("Invalid choice for 'cse' option:  ", cse)
  }
  return (cse)
}
get.csv.bits <- function (csv) {
  if (!is.null(csv)) {
    if (csv == TRUE) {
      csv <- 2^29
    }
    else if (csv == FALSE) {
      csv <- 0
    }
    else {
      stop("Invalid choice for 'csv' option:  ", csv)
    }
  }
  else {
    stop("Invalid choice for 'csv' option:  ", csv)
  }
  return (csv)
}
get.data.pass.bits <- function (data.pass) {
  if (!is.null(data.pass)) {
    if (data.pass == TRUE) {
      data.pass <- 2^15
    }
    else if (data.pass == FALSE) {
      data.pass <- 0
    }
    else {
      stop("Invalid choice for 'data.pass' option:  ", data.pass)
    }
  }
  else {
    stop("Invalid choice for 'data.pass' option:  ", data.pass)
  }
  return (data.pass)
}
get.data.pass.predict.bits <- function (data.pass) {
  if (!is.null(data.pass)) {
    if (data.pass == TRUE) {
      data.pass <- 2^27
    }
    else if (data.pass == FALSE) {
      data.pass <- 0
    }
    else {
      stop("Invalid choice for 'data.pass' option:  ", data.pass)
    }
  }
  else {
    stop("Invalid choice for 'data.pass' option:  ", data.pass)
  }
  return (data.pass)
}
get.distance.bits <- function (grow.equivalent, distance) {
  ## Convert distance option into native code parameter.
  if (!is.null(distance)) {
    if (distance == FALSE) {
      dist.bits <- 0
    }
    else if (grow.equivalent == TRUE) {
      if (distance == TRUE) {
        dist.bits <- 2^20 + 2^21
      }
      else if (distance == "inbag") {
        dist.bits <- 2^20 + 2^21
      }
      else if (distance == "oob") {
        dist.bits <- 2^20 + 2^22
      }
      else if (distance == "all") {
        dist.bits <- 2^20 + 2^21 + 2^22
      }
      else {
        stop("Invalid choice for 'distance' option:  ", distance)
      }
    }
    else if (grow.equivalent == FALSE) {
      if (distance == TRUE) {
        dist.bits <- 2^20 + 2^21 + 2^22
      }
      else if (distance == "all") {
        dist.bits <- 2^20 + 2^21 + 2^22
      }
      else {
        stop("Invalid choice for 'distance' option:  ", distance)
      }
    }
    else {
      stop("Invalid choice for 'grow.equivalent' in distance:  ", grow.equivalent)
    }
  }
  else {
    stop("Invalid choice for 'distance' option:  ", distance)
  }
  return (dist.bits)
}
## convert ensemble option into native code parameter.
get.ensemble.bits <- function (ensemble) {
  if (ensemble == "oob") {
    ensemble <- 2^1
  }
  else if (ensemble == "inbag") {
    ensemble <- 2^0
  }
  else if (ensemble == "all") {
    ensemble <- 2^0 + 2^1
  }    
  else {
    stop("Invalid choice for 'ensemble' option:  ", ensemble)
  }
  return (ensemble)
}
get.experimental.bits <- function (experimental) {
  ## Convert experimental option into native code parameter.
  bits <- 0
  for (i in 1:length(experimental)) {
    if (experimental[i] == TRUE) {
      bits <- bits + 2^(i-1)
    }
  }
  return (bits)
}
get.forest.bits <- function (forest) {
  ## Convert forest option into native code parameter.
  if (!is.null(forest)) {
    if (forest == TRUE) {
      forest <- 2^5
    }
    else if (forest == FALSE) {
      forest <- 0
    }
    else {
      stop("Invalid choice for 'forest' option:  ", forest)
    }
  }
  else {
    stop("Invalid choice for 'forest' option:  ", forest)
  }
  return (forest)
}
get.forest.wt.bits <- function (grow.equivalent, bootstrap, weight) {
  ## Convert weight option into native code parameter.
  if (!is.null(weight)) {
    if (weight == FALSE) {
      wght.bits <- 0
    }
    else if (grow.equivalent == TRUE) {
      if (weight == TRUE) {
        wght.bits <- 2^0 + 2^1
      }
      else if (weight == "inbag") {
        wght.bits <- 2^0 + 2^1
      }
      else if (weight == "oob") {
        wght.bits <- 2^0 + 2^2
      }
      else if (weight == "all") {
        wght.bits <- 2^0 + 2^1 + 2^2
      }
      else {
        stop("Invalid choice for 'weight' option:  ", weight)
      }
    }
    else if (grow.equivalent == FALSE) {
      if (weight == TRUE) {
        wght.bits <- 2^0 + 2^1 + 2^2
      }
      else if (weight == "all") {
        wght.bits <- 2^0 + 2^1 + 2^2
      }
      else {
        stop("Invalid choice for 'weight' option:  ", weight)
      }
    }
    else {
      stop("Invalid choice for 'grow.equivalent' in weight:  ", grow.equivalent)
    }
  }
  else {
    stop("Invalid choice for 'weight' option:  ", weight)
  }
  return (wght.bits)
}
get.gk.quantile <- function(gk.quantile) {
  if (is.null(gk.quantile)) {
    gk.quantile <- FALSE
  }
  gk.quantile
}
get.gk.quantile.bits <-  function (gk.quantile) {
  if (gk.quantile) {
    return (2^24)
  }
  else {
    return (0)
  }
}
## This has been deprecated!
get.insitu.ensemble.bits <- function (insitu.ensemble) {
  ## Convert insitu.ensemble option into native code parameter.
  if (!is.null(insitu.ensemble)) {
    if (insitu.ensemble == TRUE) {
      insitu.ensemble <- 2^11
    }
    else if (insitu.ensemble == FALSE) {
      insitu.ensemble <- 0
    }
    else {
      stop("Invalid choice for 'insitu.ensemble' option:  ", insitu.ensemble)
    }
  }
  else {
    stop("Invalid choice for 'insitu.ensemble' option:  ", insitu.ensemble)
  }
  return (insitu.ensemble)
}
get.importance <-  function(importance, perf.type = NULL) {
  ## convert importance option into native code parameter
  if (!is.null(importance)) {
    ## use permute with gmean as default
    if (!is.null(perf.type) && perf.type == "gmean" && as.character(importance) == TRUE) {
      importance <- "permute"
    }
    ## default vimp is now anti
    if (importance == TRUE) {
      importance <- "anti"
    }
    else if (importance == FALSE) {
      importance <- "none"
    }
  }
  importance
}
get.importance.bits <-  function(importance, perf.type = NULL) {
  ## convert importance option into native code parameter
  if (!is.null(importance)) {
    importance <- get.importance(importance, perf.type)
    if (importance == "none") {
      importance <- 0
    }
    else if (importance == "anti") {
      importance <- 2^25 + 0
    }
    else if (importance == "permute") {
      importance <- 2^25 + 2^8
    }
    else if (importance == "random") {
      importance <- 2^25 + 2^9
    }
    else if (importance == "anti.joint") {
      importance <- 2^25 + 2^10 + 0
    }
    else if (importance == "permute.joint") {
      importance <- 2^25 + 2^10 + 2^8
    }
    else if (importance == "random.joint") {
      importance <- 2^25 + 2^10 + 2^9
    }
    else {
      stop("Invalid choice for 'importance' option:  ", importance)
    }
  }
  else {
    stop("Invalid choice for 'importance' option:  ", importance)
  }
  return (importance)
}
get.impute.only.bits <-  function (impute.only, nMiss) {
  if (impute.only) {
    if (nMiss > 0) {
      return (2^16)
    }
    else {
      stop("Data has no missing values, using 'impute' makes no sense.")
    }
  }
  else {
    return (0)
  }
}
## Convert jitt option into native code parameter.
get.jitt.bits <- function (jitt) {
  if (!is.null(jitt)) {
    if (jitt == TRUE) {
      jitt <- 2^23
    }
    else if (jitt == FALSE) {
      jitt <- 0
    }
    else {
      stop("Invalid choice for 'jitt' option:  ", jitt)
    }
  }
  else {
    stop("Invalid choice for 'jitt' option:  ", jitt)
  }
  return (jitt)
}
get.membership.bits <- function (membership) {
  ## Convert option into native code parameter.
  bits <- 0
  if (!is.null(membership)) {
    if (membership == TRUE) {
      bits <- 2^6
    }
    else if (membership != FALSE) {
      stop("Invalid choice for 'membership' option:  ", membership)
    }
  }
  else {
    stop("Invalid choice for 'membership' option:  ", membership)
  }
  return (bits)
}
get.na.action.bits <- function (na.action) {
  if (na.action == "na.omit") {
    ## This is the high byte!
    na.action <- 0
  }
  else if (na.action == "na.impute" || na.action == "na.random") {
    ## This is the high byte!
    na.action <- 2^4
    ## To recover the original functionality in which the split
    ## statistic uses missing in-node imputed values, uncomment 
    ## the following statement:
    ## na.action <- 0
  }
  else {
    stop("Invalid choice for 'na.action' option:  ", na.action)
  }
  return (na.action)
}
get.outcome.bits <- function (outcome) {
  ## Convert outcome option into native code parameter.
  if (outcome == "train") {
    outcome <- 0
  }
  else if (outcome == "test") {
    outcome <- 2^17
  }
  else {
    stop("Invalid choice for 'outcome' option:  ", outcome)
  }
  return (outcome)
}
get.presort.xvar.bits <- function (presort.xvar) {
  ## Convert trace into native code parameter.
  if (!is.logical(presort.xvar)) {
    ## Leave it as is.
  }
  else {
    presort.xvar <- 1 * presort.xvar
  }
  return (presort.xvar)
}
get.perf <-  function (perf, impute.only, family) {
  ## first deal with impute.only where there is no performance
  if (impute.only == TRUE) {
    return("none")
  }
  ## now deal with non-classification
  if (family != "class") {
    if (is.null(perf)) {
      return("default")
    }
    perf <- match.arg(perf, c("none", "default", "standard"))##only allowed values
    if (perf == "standard") {
      perf <- "default"
    }
    return(perf)
  }
  ## now deal with classification
  if (is.null(perf)) {
    return("default")
  }
  perf <- match.arg(perf, c("none", "default", "standard", "misclass", "brier", "gmean", "g.mean"))
  if (perf == "standard" || perf == "misclass") {
    perf <- "default"
  }
  if (perf == "g.mean") {
    perf <- "gmean"
  }
  perf
}
get.perf.bits <- function (perf) {
  if (perf == "default") {
    return (2^2)
  }
  else if (perf == "gmean" || perf == "g.mean") {
    return (2^2 + 2^14)
  }
  else if (perf == "brier") {
    return (2^2 + 2^3)
  }
  else {#everything else becomes "none"
    return (0)
  }
}
get.proximity.bits <- function (grow.equivalent, proximity) {
  ## Convert proximity option into native code parameter.
  if (!is.null(proximity)) {
    if (proximity == FALSE) {
      prox.bits <- 0
    }
    else if (grow.equivalent == TRUE) {
      if (proximity == TRUE) {
        prox.bits <- 2^28 + 2^29
      }
      else if (proximity == "inbag") {
        prox.bits <- 2^28 + 2^29
      }
      else if (proximity == "oob") {
        prox.bits <- 2^28 + 2^30
      }
      else if (proximity == "all") {
        prox.bits <- 2^28 + 2^29 + 2^30
      }
      else {
        stop("Invalid choice for 'proximity' option:  ", proximity)
      }
    }
    else if (grow.equivalent == FALSE) {
      if (proximity == TRUE) {
        prox.bits <- 2^28 + 2^29 + 2^30
      }
      else if (proximity == "all") {
        prox.bits <- 2^28 + 2^29 + 2^30
      }
      else {
        stop("Invalid choice for 'proximity' option:  ", proximity)
      }
    }
    else {
      stop("Invalid choice for 'grow.equivalent' in proximity:  ", grow.equivalent)
    }
  }
  else {
    stop("Invalid choice for 'proximity' option:  ", proximity)
  }
  return (prox.bits)
}
get.rfq <- function(rfq) {
  if (is.null(rfq)) {
    rfq <- FALSE
  }
  rfq
}
get.rfq.bits <- function (rfq, family) {
  result <- 0
  if (family == "class") {
    if (rfq) {
      result <- 2^15
    }
  }
  return (result)
}
get.rf.cores <- function () {
    ## PART I:  Two ways for the user to specify cores:
    ## (1) R-option "rf.cores"
    ## (2) Shell-environment-option "RF_CORES"
    if (is.null(getOption("rf.cores", NULL))) {
        if (!is.na(as.numeric(Sys.getenv("RF_CORES")))) {
            options(rf.cores = as.integer(Sys.getenv("RF_CORES")))
        }
    }
    ## If the user has set the cores using either of the two methods, we respect it.
    if (!is.null(getOption("rf.cores", NULL))) {
        return (getOption("rf.cores"))
    }
    ## PART II:  Respect R CMD check limit
    chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
    if (nzchar(chk) && chk != "false") {
        ## under R CMD check --as-cran (CRAN sets this)
        return(2L)
    }
    ## PART III:  Use everything.
    return (-1L)
}
## convert samptype option into native code parameter.
get.samptype.bits <- function (samptype) {
  if (samptype == "swr") {
    bits <- 0
  }
  else if (samptype == "swor") {
    bits <- 2^12
  }
  else {
    stop("Invalid choice for 'samptype' option:  ", samptype)
  }
  return (bits)
}
get.seed <- function (seed) {
  if ((is.null(seed)) || (abs(seed) < 1)) {
    seed <- runif(1,1,1e6)
  }
  seed <- -round(abs(seed))
  return (seed)
}
get.split.cust.bits <- function (split.cust) {
  ## Convert split.cust option into native code parameter.
  if (!is.null(split.cust)) {
    if ((split.cust >= 1) && (split.cust <= 16)) {
      ## Bit shift eight left.
      split.cust <- 256 * (split.cust - 1)
    }
    else {
      stop("Invalid choice for 'split.cust' option:  ", split.cust)
    }
  }
  else {
    split.cust <- 0
  }
  return (split.cust)
}
get.split.depth.bits <- function (split.depth) {
  ## Convert split.depth option into native code parameter.
  if (!is.null(split.depth)) {
    if (split.depth == "all.trees") {
      split.depth <- 2^22
    }
    else if (split.depth == "by.tree") {
      split.depth <- 2^23
    }
    else if (split.depth == FALSE) {
      split.depth <- 0
    }
    else {
      stop("Invalid choice for 'split.depth' option:  ", split.depth)
    }
  }
  else {
    stop("Invalid choice for 'split.depth' option:  ", split.depth)
  }
  return (split.depth)
}
get.split.null.bits <- function (split.null) {
  ## Convert split.null option into native code parameter.
  if (!is.null(split.null)) {
    if (split.null == TRUE) {
      split.null <- 2^18
    }
    else if (split.null == FALSE) {
      split.null <- 0
    }
    else {
      stop("Invalid choice for 'split.null' option:  ", split.null)
    }
  }
  else {
    stop("Invalid choice for 'split.null' option:  ", split.null)
  }
  return (split.null)
}
get.tdc.rule.bits <- function (tdc.rule) {
  if (is.null(tdc.rule)) {
    ## Allow splits on everything.
    result <- 2^24 + 2^25 + 2^26
  }
  else if (tdc.rule == "all") {
    ## Allow splits on everything.
    result <- 2^24 + 2^25 + 2^26
  }
  else if (tdc.rule == "time") {
    ## Allow splits only on time.
    result <- 2^24
  }
  else if (tdc.rule == "tdx") {
    ## Allow splits on time dependent x-variables only, and not on time.
    result <- 2^25
  }
  else if (tdc.rule == "tsx") {
    ## Allow splits on time static x-variables only, and not on time.
    result <- 2^26
  }
  else if (tdc.rule == "time.and.tsx") {
    ## Allow splits on time and time static x-variables only.
    result <- 2^24 + 2^26
  }
  else if (tdc.rule == "time.and.tdx") {
    ## Allow splits on time and time dependent x-variables only.
    result <- 2^24 + 2^25
  }
  return (result)
}
get.terminal.qualts.bits <- function(terminal.qualts) {
  bits <- 0
  if (!is.null(terminal.qualts)) {
    if (terminal.qualts == TRUE) {
      bits <- bits + 2^16
    }
    else if (terminal.qualts != FALSE) {
      stop("Invalid choice for 'terminal.qualts' option:  ", terminal.qualts)
    }
  }
  else {
    stop("Invalid choice for 'terminal.qualts' option:  ", terminal.qualts)
  }
  return(bits)
}
get.terminal.qualts.predict.bits <- function(incoming.flag) {
  bits <- 2^16
  if (!is.null(incoming.flag)) {
    if (incoming.flag == TRUE) {
      bits <- bits + 2^17
    }
  }
  else {
    stop("Invalid choice for 'incoming.flag':  ", incoming.flag)
  }
  return(bits)
}
get.terminal.quants.bits <- function(terminal.quants) {
  bits <- 0
  if (!is.null(terminal.quants)) {
    if (terminal.quants == TRUE) {
      bits <- bits + 2^18
    }
    else if (terminal.quants != FALSE) {
      stop("Invalid choice for 'terminal.quants' option:  ", terminal.quants)
    }
  }
  else {
    stop("Invalid choice for 'terminal.quants' option:  ", terminal.quants)
  }
  return(bits)
}
get.terminal.quants.predict.bits <- function(incoming.flag) {
  bits <- 2^18
  if (!is.null(incoming.flag)) {
    if (incoming.flag == TRUE) {
      bits <- bits + 2^19
    }
  }
  else {
    stop("Invalid choice for 'incoming.flag':  ", incoming.flag)
  }
  return(bits)
}
get.trace.bits <- function (do.trace) {
  ## Convert trace into native code parameter.
  if (!is.logical(do.trace)) {
    if (do.trace >= 1) {
      do.trace <- round(do.trace)
    }
    else {
      do.trace <- 0
    }
  }
  else {
    do.trace <- 1 * do.trace
  }
  return (do.trace)
}
get.tree.index <- function(get.tree, ntree) {
  ## NULL --> default setting
  if (is.null(get.tree)) {
    rep(1, ntree)
  }
  ## the user has specified a subset of trees
  else {
    pt <- get.tree >=1 & get.tree <= ntree
    if (sum(pt) > 0) {
      get.tree <- get.tree[pt]
      get.tree.temp <- rep(0, ntree)
      get.tree.temp[get.tree] <- 1
      get.tree.temp
    }
    else {
      rep(1, ntree)
    }
  }
}
get.var.used.bits <- function (var.used) {
  ## Convert var.used option into native code parameter.
  if (!is.null(var.used)) {
    if (var.used == "all.trees") {
      var.used <- 2^12
    }
    else if (var.used == "by.tree") {
      var.used <- 2^13
    }
    else if (var.used == FALSE) {
      var.used <- 0
    }
    else {
      stop("Invalid choice for 'var.used' option:  ", var.used)
    }
  }
  else {
    stop("Invalid choice for 'var.used' option:  ", var.used)
  }
  return (var.used)
}
get.vimp.only.bits <-  function (vimp.only) {
  ## Convert vimp.only option into native code parameter.
  if (!is.null(vimp.only)) {
    if (vimp.only) {
      return (2^27)
    }
    else if (!vimp.only) {
      return (0)
    }
    else {
      stop("Invalid choice for 'vimp.only' option:  ", vimp.only)
    }
  }
  else {
    stop("Invalid choice for 'vimp.only' option:  ", vimp.only)
  }
}
## Check for presence of forest
is.forest.missing <- function(object) {
  ## for backwards compatability
  if(is.null(object$forest$forest)) {
    is.null(object$forest)
  }
  ## current stealth build moving forward
  else {
    !object$forest$forest
  }
}
is.hidden.cse <-  function(user.option) {
  if (is.null(user.option$cse)) {
    FALSE
  }
  else {
    as.logical(as.character(user.option$cse))
  }
}
is.hidden.csv <-  function(user.option) {
  if (is.null(user.option$csv)) {
    FALSE
  }
  else {
    as.logical(as.character(user.option$csv))
  }
}
is.hidden.data.pass <-  function(user.option) {
  if (is.null(user.option$data.pass)) {
    FALSE
  }
  else {
    as.logical(as.character(user.option$data.pass))
  }
}
is.hidden.do.trace <-  function(user.option) {
  if (is.null(user.option$do.trace)) {
    FALSE
  }
  else {
    user.option$do.trace
  }
}
is.hidden.experimental <-  function(user.option) {
  if (is.null(user.option$experimental)) {
    e <- FALSE
  }
  else {
    e  <- as.logical(as.character(user.option$experimental))
  }
  return (e)
}
is.hidden.gk.quantile <-  function(user.option) {
  if (is.null(user.option$gk.quantile)) {
    NULL
  }
  else {
    as.logical(as.character(user.option$gk.quantile))
  }
}
is.hidden.holdout.array <-  function(user.option) {
  if (is.null(user.option$holdout.array)) {
    NULL
  }
  else {
    user.option$holdout.array
  }
}
is.hidden.holdout.specs <-  function(user.option) {
  ## Default value is NULL
  if (is.null(user.option$holdout.specs)) {
    obj <- NULL
  }
  else {
    obj <- user.option$holdout.specs
  }
  return (obj)
}
is.hidden.impute.only <-  function(user.option) {
  if (is.null(user.option$impute.only)) {
    FALSE
  }
  else {
    as.logical(as.character(user.option$impute.only))
  }
}
## This has been deprecated!
is.hidden.insitu.ensemble <-  function(user.option) {
  if (is.null(user.option$insitu.ensemble)) {
    ## From > 2.14.0 we now assert insitu.ensemble as the default. 
    TRUE
  }
  else {
    as.logical(as.character(user.option$insitu.ensemble))
  }
}
is.hidden.jitt <-  function (user.option, importance, na.action = "na.omit", anonymous = FALSE,  partial = FALSE) {
  ## no missing data
  ## anonymous = TRUE, jitt = TRUE,
  ## anonymous = TRUE, jitt = FALSE,
  ## anonymous = FALSE, jitt = TRUE,
  ## anonymous = FALSE, jitt = FALSE,
  ## missing data
  ## anonymous = TRUE, na.action=na.mean  =>>  no restriction on jitt
  ## anonymous = TRUE, jitt = TRUE   =>>  na.random
  ## anonymous = TRUE, jitt = FALSE  =>>  failure and exit
  ## anonymous = FALSE, jitt = TRUE  =>>  na.random
  ## anonymous = FALSE, jitt = FALSE =>>  na.impute (old-school)
  if (na.action == "na.random") {
    return(TRUE)
  }
  if (na.action == "na.impute" && !anonymous) {
    return(FALSE)
  }
  ## now that the danger is passed, proceed forward ...
  if (is.null(user.option$jitt)) {
    if (partial || importance == "permute" || importance == "permute.joint") {
      ##partial is generally a full restore
      ##permutation importance is synthetic, so default JITT = FALSE
      FALSE
    }
    else {
      TRUE
    }
  }
  else {
    as.logical(as.character(user.option$jitt))
  }
}
is.hidden.mahalanobis.sigma <- function(user.option) {
  if (is.null(user.option$mahalanobis.sigma) & is.null(user.option$sigma)) {
    NULL
  }
  else {
    if (!is.null(user.option$mahalanobis.sigma)) {
      ginverse(user.option$mahalanobis.sigma)
    }
    else {
      ginverse(user.option$sigma)
    }
  }
}
is.hidden.perf.type <-  function(user.option) {
  ## Default value is NULL
  if (is.null(user.option$perf.type)) {
    NULL
  }
  else {
    as.character(user.option$perf.type)
  }
}
is.hidden.presort.xvar <-  function(user.option) {
  if (is.null(user.option$presort.xvar)) {
    FALSE
  }
  else {
    user.option$presort.xvar
  }
}
is.hidden.prob <-  function(user.option) {
  if (is.null(user.option$prob)) {
    NULL
  }
  else {
    prob <- user.option$prob
    sort(prob[prob>0 & prob<1])
  }
}
is.hidden.prob.epsilon <-  function(user.option) {
  if (is.null(user.option$prob.epsilon)) {
    NULL
  }
  else {
    prob.epsilon <- user.option$prob.epsilon
    if ((prob.epsilon <= 0) || (prob.epsilon >= 0.50)) {
      stop("parameter 'prob.epsilon' must be in range (0, 1/2) :  ", prob.epsilon)
    }
    prob.epsilon
  }
}
is.hidden.quantile.regr <-  function(user.option) {
  if (is.null(user.option$quantile.regr)) {
    FALSE
  }
  else {
    as.logical(as.character(user.option$quantile.regr))
  }
}
is.hidden.rfq <-  function(user.option) {
  if (is.null(user.option$rfq)) {
    NULL
  }
  else {
    as.logical(as.character(user.option$rfq))
  }
}
is.hidden.tdc.rule <-  function(user.option) {
  if (is.null(user.option$tdc.rule)) {
    NULL
  }
  else {
    user.option$tdc.rule
  }
}
is.hidden.terminal.qualts <-  function(user.option) {
  ## From > 2.14.0 we now assert terminal qualts and quants as the default. 
  if (is.null(user.option$terminal.qualts)) {
    TRUE
  }
  else {
    as.logical(as.character(user.option$terminal.qualts))
  }
}
is.hidden.terminal.quants <-  function (user.option, save.memory = FALSE) {
  ## From > 2.14.0 we now assert terminal qualts and quants as the default - save.memory override though 
  if (is.null(user.option$terminal.quants)) {
    TRUE & !save.memory
  }
  else {
    as.logical(as.character(user.option$terminal.quants))
  }
}
is.hidden.vtry <-  function(user.option) {
  ## Default value is 0
  if (is.null(user.option$vtry)) {
    0
  }
  else {
    user.option$vtry
  }
}
is.hidden.ytry <-  function(user.option) {
  if (is.null(user.option$ytry)) {
    NULL
  }
  else {
    as.integer(user.option$ytry)
  }
}
is.hidden.vimp.threshold <-  function(user.option) {
  if (is.null(user.option$vimp.threshold)) {
    ## This is value from zero (0) to one (1) that represents the
    ## probability that variables will be noised up for random and
    ## anti vimp.  Thus zero means that vimp is effectively turned off.
    ## One means that the vimp protocol will be respected all the
    ## time.  If the threshold value alpha is such that 0 < alpha
    ## < 1, we draw a random value b. If b <= alpha, we respect
    ## the vimp protocol.
    1.0
  }
  else {
    user.option$vimp.threshold
  }
}
