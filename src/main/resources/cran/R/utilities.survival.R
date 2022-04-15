####################################################################
##
## survival functions
##
####################################################################
get.event.info <- function(obj, subset = NULL) {
  ## survival case
  if (grepl("surv", obj$family)) {
    if (!is.null(obj$yvar)) {
      if (is.null(subset)) {
        subset <- (1:nrow(cbind(obj$yvar)))
      }
        if (is.null(obj$subj)) { 
            r.dim <- 2
            time <- obj$yvar[subset, 1]
            cens <- obj$yvar[subset, 2]
        }
        else {
            r.dim <- 3
            start.time <- obj$yvar[subset, 1]
            time <- obj$yvar[subset, 2]
            cens <- obj$yvar[subset, 3]
        }
      ## censoring must be coded coherently
      if (!all(floor(cens) == abs(cens), na.rm = TRUE)) {
        stop("for survival families censoring variable must be coded as a non-negative integer")
      }
      ## Extract the unique event types.
      event <- na.omit(cens)[na.omit(cens) > 0]
      event.type <- sort(unique(event))
    }
    ##everything else
    else {
      r.dim <- 0
      event <- event.type <- cens <- cens <- time <- NULL
    }
    ## Set grid of time points.
    time.interest <- obj$time.interest
  }
  else {
    ## NULL for other families
    if ((obj$family == "regr+") | (obj$family == "class+")) {
      r.dim <- dim(obj$yvar)[2]
    }
    else {
      r.dim <- 1
    }
    event <- event.type <- cens <- time.interest <- cens <- time <- NULL
  }
  return(list(event = event, event.type = event.type, cens = cens,
              time.interest = time.interest, time = time, r.dim = r.dim))
}
get.grow.event.info <- function(yvar, fmly, need.deaths = TRUE, ntime = NULL) {
  if (grepl("surv", fmly)) {
    ##-----------------------------------------------------------
    ## survival, competing risks, or time dependent covariates
    ##-----------------------------------------------------------
    if (dim(yvar)[2] == 2) {
      ##---------------------------------
      ## survival or competing risks:
      ##---------------------------------
      r.dim <- 2
      time <- yvar[, 1]
      cens <- yvar[, 2]
      start.time <- NULL
      ## censoring must be coded coherently
      if (!all(floor(cens) == abs(cens), na.rm = TRUE)) {
        stop("for survival families censoring variable must be coded as a non-negative integer (perhaps the formula is set incorrectly?)")
      }
      ## check if deaths are available (if user specified)
      if (need.deaths && (all(na.omit(cens) == 0))) {
        stop("no deaths in data!")
      }
      ## Check for event time consistency.
      ## we over-ride this now to allow for negative time (see Stute)
      ##if (!all(na.omit(time) >= 0)) {
      ##  stop("time must be  positive")
      ##}
      ## Extract the unique event types.
      event.type <- unique(na.omit(cens))
      ## Ensure they are all greater than or equal to zero.
      if (sum(event.type >= 0) != length(event.type)) {
        stop("censoring variable must be coded as NA, 0, or greater than 0.")
      }
      ## Discard the censored state, if it exists.
      event <- na.omit(cens)[na.omit(cens) > 0]
      event.type <- unique(event)
      ## Set grid of time points.
      nonMissingOutcome <- which(!is.na(cens) & !is.na(time))
      nonMissingDeathFlag <- (cens[nonMissingOutcome] != 0)
      time.interest <- sort(unique(time[nonMissingOutcome[nonMissingDeathFlag]]))
      ## trim the time points if the user has requested it
      ## we also allow the user to pass requested time points
      if (!is.null(ntime) && !((length(ntime) == 1) && ntime == 0)) {
        if (length(ntime) == 1 && length(time.interest) > ntime) {
          time.interest <- time.interest[
            unique(round(seq.int(1, length(time.interest), length.out = ntime)))]
        }
        if (length(ntime) > 1) {
          time.interest <- unique(sapply(ntime, function(tt) {
            time.interest[max(1, sum(tt >= time.interest, na.rm = TRUE))]
          }))
        }
      }
    }
    ##-------------------------------
    ## time dependent covariates:
    ##-------------------------------
    else {
      r.dim <- 3
      start.time <- yvar[, 1]
      time <- yvar[, 2]
      cens <- yvar[, 3]
      ## censoring must be coded coherently
      if (!all(floor(cens) == abs(cens), na.rm = TRUE)) {
        stop("for survival families censoring variable must be coded as a non-negative integer (perhaps the formula is set incorrectly?)")
      }
      ## check if deaths are available (if user specified)
      if (need.deaths && (all(na.omit(cens) == 0))) {
        stop("no deaths in data!")
      }
      ## Check for event time consistency.
      if (!all(na.omit(time) >= 0)) {
        stop("time must be  positive")
      }
      ## Extract the unique event types.
      event.type <- unique(na.omit(cens))
      ## Ensure they are all greater than or equal to zero.
      if (sum(event.type >= 0) != length(event.type)) {
        stop("censoring variable must be coded as NA, 0, or greater than 0.")
      }
      ## Discard the censored state, if it exists.
      event <- na.omit(cens)[na.omit(cens) > 0]
      event.type <- unique(event)
      ## Set grid of time points.
      nonMissingOutcome <- which(!is.na(cens) & !is.na(time))
      nonMissingDeathFlag <- (cens[nonMissingOutcome] != 0)
      time.interest <- sort(unique(time[nonMissingOutcome[nonMissingDeathFlag]]))
      ## trim the time points if the user has requested it
      ## we also allow the user to pass requested time points
      if (!is.null(ntime) && !((length(ntime) == 1) && ntime == 0)) {
        if (length(ntime) == 1 && length(time.interest) > ntime) {
          ## select evenly spaced values over [0,1] and not event times 
          time.interest <- seq(0,  min(1, max(time[nonMissingOutcome])), length = ntime)
          time.interest <- time.interest[time.interest > 0]
        }
        if (length(ntime) > 1) {
          ## over-ride the default setting and allow the user to specify anything they want between [0,1]
          time.pt <- ntime <= min(1, max(time[nonMissingOutcome])) & ntime > 0
          if (sum(time.pt) == 0) {
            stop("the ntime vector supplied must be between [0,1]:", ntime)
          }
          time.interest <- sort(unique(ntime[time.pt]))
        }
      }
    }
  }
  ##---------------------
  ## other families
  ##---------------------
  else {
    if ((fmly == "regr+") | (fmly == "class+") | (fmly == "mix+")) {
      r.dim <- dim(yvar)[2]
    }
    else {
      if (fmly == "unsupv") {
        r.dim <- 0
      }
      else {
        r.dim <- 1
      }
    }
    event <- event.type <- cens <- time.interest <- cens <- time <- start.time <- NULL
  }
  return(list(event = event, event.type = event.type, cens = cens,
              time.interest = time.interest,
              time = time, start.time = start.time, r.dim = r.dim))
}
## ---------------------------------------------------------------------
##
## brier score
##
## ---------------------------------------------------------------------
## trapezoidal rule
trapz <- function (x, y) {
  idx = 2:length(x)
  return(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2)
}
## returns an index of positions for evaluating a step function at selected times
sIndex <- function(x,y) {sapply(1:length(y), function(j) {sum(x <= y[j])})}
## main brier function
get.brier.survival <- function(o, subset, cens.model = c("km", "rfsrc"), papply = lapply) {
  ## incoming parameter checks
  if (is.null(o)) {
    return(NULL)
  }
  if (o$family != "surv") {
    stop("this function only supports right-censored survival settings")
  }
  if (sum(inherits(o, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(o, c("rfsrc", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("This function only works for objects of class `(rfsrc, grow)' or '(rfsrc, predict)'")
  }
  ## use imputed missing time or censoring indicators
  if (!is.null(o$yvar) && !is.null(o$imputed.indv)) {
    o$yvar[o$imputed.indv, ] <- o$imputed.data[, 1:2]
  }
  ## verify the cens.model option
  cens.model <- match.arg(cens.model, c("km", "rfsrc"))
  ## subsetting: assumes entire data set to be used if not specified
  if (missing(subset) || is.null(subset)) {
    subset <- 1:length(o$predicted)
  }
  else {
    ## convert the user specified subset into a usable form
    if (is.logical(subset)) subset <- which(subset)
    subset <- unique(subset[subset >= 1 & subset <= length(o$predicted)])
    if (length(subset) == 0) {
      stop("'subset' not set properly.")
    }
  }
  ## yvar is used for building the training (grow) censoring distribution
  ## however, there is no guarantee that yvar will exist in predict mode
  ## the forest however always contains the yvar in play and so we use that
  pred.no.y <- is.null(o$yvar)
  yvar <- o$forest$yvar
  o$yvar <- yvar
  event.info <- get.event.info(o)
  ## obtain subset event info, but then put original yvar back
  if (!pred.no.y) {
    o$yvar <- yvar[subset,, drop = FALSE]
    subset.event.info <- get.event.info(o)
    o$yvar <- yvar
  }
  ## use OOB values if available
  if (is.null(o$predicted.oob)) {
    mort <- o$predicted[subset]
    surv.ensb <- t(o$survival[subset,, drop = FALSE])
  }
  else {
    mort <- o$predicted.oob[subset]
    surv.ensb <- t(o$survival.oob[subset,, drop = FALSE])
  }
  ##-------------------------------------------------------------------------------
  ##
  ## KM for training/testing data - for testing, there must be y
  ## match time to grow master list, time.interest
  ##
  ##-------------------------------------------------------------------------------
  if (!pred.no.y) {
    km.obj <- do.call(rbind, papply(1:length(subset.event.info$time.interest), function(j) {
      c(sum(subset.event.info$time >= subset.event.info$time.interest[j], na.rm = TRUE),
        sum(subset.event.info$time[subset.event.info$cens != 0] == subset.event.info$time.interest[j], na.rm = TRUE))
    }))
    Y <- km.obj[, 1]
    d <- km.obj[, 2]
    r <- d / (Y + 1 * (Y == 0))
    surv.aalen <- exp(-cumsum(r))[1 + sIndex(subset.event.info$time.interest, event.info$time.interest)]
  }
  else {
    surv.aalen <- NULL
  }
  ##-------------------------------------------------------------------------------
  ##
  ## censoring distribution estimator for training (grow) data
  ## match time to grow master list, time.interest
  ##
  ##-------------------------------------------------------------------------------
  ## we match the censoring times with the master list time.interest
  ## this unifies all further calculations
  censTime <- sort(unique(event.info$time[event.info$cens == 0]))
  censTime.pt <- c(sIndex(censTime, event.info$time.interest))
  ## check to see if there are censoring cases
  if (length(censTime) > 0) {
    ## KM censoring distribution estimator
    if (cens.model == "km") {
      censModel.obj <- do.call(rbind, papply(1:length(censTime), function(j) {
        c(sum(event.info$time >= censTime[j], na.rm = TRUE),
          sum(event.info$time[event.info$cens == 0] == censTime[j], na.rm = TRUE))
      }))
      Y <- censModel.obj[, 1]
      d <- censModel.obj[, 2]
      r <- d / (Y + 1 * (Y == 0))
      cens.dist <- c(1, exp(-cumsum(r)))[1 + censTime.pt]
    }
    ## rfsrc censoring distribution estimator
    else {
      cens.dta <- data.frame(tm = o$forest$yvar[, 1],
                             cens = 1 * (o$yvar[, 2] == 0),
                             o$forest$xvar[,, drop = FALSE])
      cens.o <- rfsrc(Surv(tm, cens) ~ ., cens.dta,
                      ntree = 100, nsplit = 1)
      cens.dist <- cens.o$survival.oob
      censTime.pt <- c(sIndex(cens.o$time.interest, event.info$time.interest))
      cens.dist <- t(cbind(1, cens.dist)[, 1 + censTime.pt])
    }
  }
  ## no censoring cases; assign a default distribution
  else {
    cens.dist <- rep(1, length(censTime.pt))
  }
  ## brier calculations
  brier.matx <- do.call(rbind, papply(1:ncol(surv.ensb), function(i) {
    tau <-  event.info$time
    event <- event.info$cens
    t.unq <- event.info$time.interest
    cens.pt <- sIndex(t.unq, tau[i])
    if (cens.model == "km") {
      c1 <- 1 * (tau[i] <= t.unq & event[i] != 0)/c(1, cens.dist)[1 + cens.pt]
      c2 <- 1 * (tau[i] > t.unq) / cens.dist
    }
    else {
      c1 <- 1 * (tau[i] <= t.unq & event[i] != 0)/c(1, cens.dist[, i])[1 + cens.pt]
      c2 <- 1 * (tau[i] > t.unq) / cens.dist[, i]
    }
    (1 * (tau[i] > t.unq) - surv.ensb[, i])^2 * (c1 + c2)
  }))
  brier.score <- data.frame(time = event.info$time.interest,
                            brier.score = colMeans(brier.matx, na.rm = TRUE))
  ## crps - continuous rank probability score
  crps <- trapz(brier.score$time, brier.score$brier.score)
  ## return the goodies
  list(brier.matx = brier.matx,
       brier.score = brier.score,
       cens.dist = cens.dist,
       crps = crps,
       crps.std = crps / max(brier.score$time),
       time = event.info$time.interest,
       event.info = event.info,
       subset = subset,
       mort = mort,
       surv.aalen = surv.aalen,
       surv.ensb = surv.ensb)
}
