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
## rmst 
##
## ---------------------------------------------------------------------
get.rmst <- function(o, tau.horizon = NULL, q = .95) {
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
  ## extract time, survival (use OOB values if available)
  time <- o$time.interest
  if (is.null(o$survival.oob)) {
    surv <- o$survival.oob
  }
  else {
    surv <- o$survival
  }
  ## set the time horizon
  if (is.null(tau.horizon)) {
    ## can replace this with maximum
    ## tau.horizon <- max(time, na.rm = TRUE)
    tau.horizon <- quantile(time, probs = q, na.rm = TRUE)
  }
  ## adjustment for when time doesn't include tau.horizon
  etime <- sort(unique(c(time, tau.horizon)))
  surv <- cbind(1, surv)[, 1 + sIndex(time, etime)]
  time <- etime
  ## restrict time to tau horizon
  time.pt <- time <= tau.horizon
  ## calculate rmst for the restricted time
  c(surv[, time.pt, drop = FALSE] %*% diff(c(0, time[time.pt])))
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
## set nodesize
set.nodesize <- function(n, p, nodesize = NULL) {
  if (is.null(nodesize)) {
    if (n <= 300 & p > n) {
      nodesize <- 2
    }
    else if (n <= 300 & p <= n) {
      nodesize <- 5
    }
    else if (n > 300 & n <= 2000) {
      nodesize <- 10
    }
    else {
      nodesize <- n / 200
    }
  }
  nodesize
}
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
      sum(inherits(o, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(o, c("rfsrc", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("This function only works for objects of class `(rfsrc, grow)', '(rfsrc, forest)' or '(rfsrc, predict)'")
  }
  ## special handling if object is a forest
  if (sum(inherits(o, c("rfsrc", "forest"), TRUE) == c(1, 2)) == 2) {
    predO <- predict(o, perf.type = "none")
    o$predicted <- predO$predicted
    o$predicted.oob <- predO$predicted.oob
    o$survival.oob <- predO$survival.oob
    o$forest <- list()
    o$forest$yvar <- o$yvar
    o$forest$xvar <- o$xvar
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
  ## the forest however always contains yvar, so we use that
  ## also see above for special handling of forest
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
    ## rfsrc censoring distribution estimator using random splitting
    else {
      cens.dta <- data.frame(time = o$forest$yvar[, 1],
                             cens = 1 * (o$forest$yvar[, 2] == 0),
                             o$forest$xvar)
      cens.o <- rfsrc(Surv(time, cens) ~ ., cens.dta,                      
                      ntree = 50,
                      nsplit = 1,
                      splitrule = "random",
                      nodesize = set.nodesize(nrow(cens.dta), ncol(o$forest$xvar)),
                      perf.type = "none")
      if (!is.null(o$imputed.indv)) {
        o$xvar[o$imputed.indv, ] <- o$imputed.data[, -(1:2)]
      }
      cens.dist <- predict(cens.o, o$xvar[subset,, drop = FALSE])$survival
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
## ------------------------------------------------------------
## Uno weights
## - training mode: KM or OOB KM
## - test mode: works generically
## - censors happen after deaths at tied times
## ------------------------------------------------------------
## fit censoring KM on training outcomes only
km_censor_fit <- function(time, status) {
  stopifnot(length(time) == length(status))
  ok <- !is.na(time) & !is.na(status)
  time   <- as.numeric(time[ok])
  status <- as.integer(status[ok])
  n <- length(time)
  if (n == 0L) stop("No non-missing training outcomes.")
  ord <- order(time)
  t <- time[ord]
  s <- status[ord]
  times <- numeric(n)
  G     <- numeric(n)
  surv   <- 1.0
  n_risk <- n
  k <- 0L
  i <- 1L
  while (i <= n) {
    ti <- t[i]
    j <- i
    d_death <- 0L
    d_cens  <- 0L
    while (j <= n && t[j] == ti) {
      if (s[j] == 1L) d_death <- d_death + 1L else d_cens <- d_cens + 1L
      j <- j + 1L
    }
    ## Update censoring survival AFTER removing deaths at this time
    n_after_death <- n_risk - d_death
    if (d_cens > 0L) {
      if (n_after_death <= 0L) {
        surv <- 0.0
      } else {
        surv <- surv * (1 - d_cens / n_after_death)
      }
    }
    k <- k + 1L
    times[k] <- ti
    G[k]     <- surv
    n_risk <- n_risk - (d_death + d_cens)
    i <- j
  }
  list(time = times[1L:k], G = G[1L:k])
}
## Generic left-limit step evaluation: returns Ghat(t-)
## for arbitrary t_new using knots/time_knots and post-step values G.
uno_Ghat_minus_predict <- function(time_knots, G, t_new) {
  t_new <- as.numeric(t_new)
  out <- rep(NA_real_, length(t_new))
  ok <- !is.na(t_new)
  if (!any(ok)) return(out)
  if (length(time_knots) == 0L) {
    out[ok] <- 1.0
    return(out)
  }
  idx <- findInterval(t_new[ok], time_knots, left.open = TRUE)
  out[ok] <- ifelse(idx == 0L, 1.0, G[idx])
  out
}
## Effective sample size of positive weights
uno_ess <- function(w) {
  w <- w[is.finite(w) & !is.na(w) & (w > 0)]
  if (length(w) == 0L) return(NA_real_)
  (sum(w)^2) / sum(w^2)
}
## Choose gmin automatically from training event-time Ghat(t-).
##
## Input:  G_event = vector of Ghat(t-) evaluated at event times only.
## Output: list(gmin, ess_target, ess_kept, n_events, n_dropped, wmax_kept)
##
## Rule: drop the largest weights (smallest G) until ESS >= ess_target,
## where ess_target = max(ess_min, ceil(ess_frac * n_events)).
uno_choose_gmin_auto <- function(G_event,
                                 eps = 1e-12,
                                 ess_frac = 0.20,
                                 ess_min  = 20L) {
  g <- as.numeric(G_event)
  g <- g[is.finite(g) & !is.na(g)]
  d <- length(g)
  if (d <= 1L) {
    return(list(gmin = 0.0, ess_target = NA_real_, ess_kept = NA_real_,
                n_events = d, n_dropped = 0L, wmax_kept = NA_real_))
  }
  ## If essentially no censoring (G ~ 1), no trimming needed
  if (min(g) >= 1 - 1e-12) {
    return(list(gmin = 0.0, ess_target = d, ess_kept = d,
                n_events = d, n_dropped = 0L, wmax_kept = 1.0))
  }
  ## weights are monotone in g: smaller g => larger weight
  g_sorted <- sort(g)                       # ascending g
  w_desc   <- 1.0 / pmax(g_sorted, eps)^2   # descending weights
  ## prefix sums (with leading 0)
  p1 <- c(0.0, cumsum(w_desc))
  p2 <- c(0.0, cumsum(w_desc * w_desc))
  ess_target <- max(as.integer(ess_min), as.integer(ceiling(ess_frac * d)))
  ess_target <- min(ess_target, d)
  ## drop k largest weights; must keep at least ess_target events
  best_k <- 0L
  best_ess <- NA_real_
  for (k in 0L:(d - ess_target)) {
    sum_w  <- p1[d + 1L] - p1[k + 1L]
    sum_w2 <- p2[d + 1L] - p2[k + 1L]
    ess_k  <- if (sum_w2 > 0) (sum_w * sum_w) / sum_w2 else NA_real_
    if (is.finite(ess_k) && (ess_k >= ess_target)) {
      best_k <- k
      best_ess <- ess_k
      break
    }
  }
  ## If never hit the target (rare), keep only ess_target events
  if (!is.finite(best_ess)) {
    best_k <- d - ess_target
    sum_w  <- p1[d + 1L] - p1[best_k + 1L]
    sum_w2 <- p2[d + 1L] - p2[best_k + 1L]
    best_ess <- if (sum_w2 > 0) (sum_w * sum_w) / sum_w2 else NA_real_
  }
  gmin <- g_sorted[best_k + 1L]
  wmax <- 1.0 / pmax(gmin, eps)^2
  list(gmin = gmin,
       ess_target = ess_target,
       ess_kept = best_ess,
       n_events = d,
       n_dropped = best_k,
       wmax_kept = wmax)
}
## Train-mode Uno weights
get.uno.weights.train <- function(time, status,
                                  gmin = "auto",
                                  ess_frac = 0.20,
                                  ess_min  = 20L,
                                  eps = 1e-12,
                                  eps_keep = .Machine$double.eps,
                                  drop_if_G0 = FALSE,
                                  return_fit = TRUE) {
  stopifnot(length(time) == length(status))
  ## Fit KM censoring curve on training outcomes
  fit <- km_censor_fit(time, status)
  ## Global Ghat(t-) used for gating and (also) weight magnitude
  G_gate <- uno_Ghat_minus_predict(fit$time, fit$G, time)
  ## Decide gmin (train-once)
  if (is.character(gmin)) {
    gmin <- match.arg(gmin, c("auto", "none"))
    if (gmin == "none") {
      gmin_used <- 0.0
      ginfo <- list(gmin = 0.0)
    } else {
      ## events = non-censored cases (otherwise breaks for CR)
      ev <- !is.na(status) & (as.integer(status) != 0L) & !is.na(G_gate)
      ginfo <- uno_choose_gmin_auto(G_gate[ev], eps = eps,
                                    ess_frac = ess_frac, ess_min = ess_min)
      gmin_used <- ginfo$gmin
    }
  } else {
    gmin_used <- as.numeric(gmin)[1L]
    if (!is.finite(gmin_used) || gmin_used < 0) gmin_used <- 0.0
    ginfo <- list(gmin = gmin_used)
  }
  ## Missing => exclude (weight 0)
  miss <- is.na(time) | is.na(status)
  G_gate[miss] <- NA_real_
  w <- rep(0.0, length(time))
  ok <- !is.na(G_gate)
  if (!drop_if_G0) {
    ## keep-as-comparator always; event contribution trimmed by gmin
    keep <- ok & (G_gate >= gmin_used)
    drop <- ok & !keep
    if (any(keep)) {
      Gsafe <- pmax(G_gate[keep], eps)
      w[keep] <- 1.0 / (Gsafe * Gsafe)
    }
    if (any(drop)) {
      w[drop] <- eps_keep
    }
  } else {
    ## trim additionally when G is essentially 0
    keep <- ok & (G_gate >= gmin_used) & (G_gate > eps)
    drop <- ok & !keep
    if (any(keep)) {
      Gsafe <- pmax(G_gate[keep], eps)
      w[keep] <- 1.0 / (Gsafe * Gsafe)
    }
    if (any(drop)) {
      w[drop] <- eps_keep
    }
  }
  ## Store once; test will reuse automatically
  fit$gmin      <- gmin_used
  fit$eps_keep  <- eps_keep
  fit$gmin_info <- ginfo
  if (return_fit) {
    return(list(weight = w,
                Ghat_minus = G_gate,
                fit = fit))
  }
  w
}
## Test-mode Uno weights
get.uno.weights.test <- function(time_test, fit,
                                 eps = 1e-12,
                                 drop_if_G0 = FALSE) {
  if (is.null(fit$time) || is.null(fit$G))
    stop("fit must be a list with elements $time and $G")
  gmin_used <- if (!is.null(fit$gmin) && is.finite(fit$gmin)) fit$gmin else 0.0
  eps_keep  <- if (!is.null(fit$eps_keep) && is.finite(fit$eps_keep)) fit$eps_keep else .Machine$double.eps
  G_gate <- uno_Ghat_minus_predict(fit$time, fit$G, time_test)
  w <- rep(0.0, length(time_test))
  ok <- !is.na(G_gate)
  if (!drop_if_G0) {
    keep <- ok & (G_gate >= gmin_used)
  } else {
    keep <- ok & (G_gate >= gmin_used) & (G_gate > eps)
  }
  drop <- ok & !keep
  if (any(keep)) {
    Gsafe <- pmax(G_gate[keep], eps)
    w[keep] <- 1.0 / (Gsafe * Gsafe)
  }
  if (any(drop)) {
    w[drop] <- eps_keep
  }
  w
}
## ------------------------------------------------------------
## Convenience helper for test evaluation
## ------------------------------------------------------------
uno.prepare.test <- function(time_test, status_test, fit,
                             eps = 1e-12, drop_if_G0 = FALSE) {
  w <- get.uno.weights.test(time_test, fit, eps = eps, drop_if_G0 = drop_if_G0)
  list(time = time_test, status = status_test, weight = w)
}
## ------------------------------------------------------------
## Convenience one-liner: training weights only
## ------------------------------------------------------------
get.uno.weights <- function(time, status) {
  get.uno.weights.train(time, status, return_fit = FALSE)
}
