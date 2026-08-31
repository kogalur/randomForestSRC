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
## Brier score and time-dependent AUC
##
## ---------------------------------------------------------------------
## trapezoidal rule
trapz <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y)) {
    stop("'x' and 'y' must have the same length.")
  }
  if (length(x) < 2L) {
    return(0)
  }
  idx <- 2:length(x)
  sum((x[idx] - x[idx - 1L]) * (y[idx] + y[idx - 1L]) / 2)
}
## returns an index of positions for evaluating a step function at selected times
sIndex <- function(x, y) {
  findInterval(y, x)
}
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
## exact element extraction for stripped RF-SRC objects
##
## The dollar operator performs partial matching on lists.  This matters for
## rfsrc.fast(..., forest = FALSE), where xvar is absent but xvar.names is
## retained: o$forest$xvar can then partially match xvar.names.  Use exact
## extraction throughout the performance-data path.
.surv.get.element <- function(x, name) {
  if (is.null(x) || !is.list(x)) {
    return(NULL)
  }
  x[[name]]
}
## identify a supported RF-SRC object
.surv.object.type <- function(o, allow.forest = TRUE) {
  if (!inherits(o, "rfsrc")) {
    return(NULL)
  }
  if (inherits(o, "grow")) {
    return("grow")
  }
  if (inherits(o, "predict")) {
    return("predict")
  }
  if (allow.forest && inherits(o, "forest")) {
    return("forest")
  }
  NULL
}
## convert a forest object to a grow-like object with OOB predictions
.surv.prepare.forest <- function(o) {
  pred.o <- predict(o, perf.type = "none")
  o[["predicted"]] <- .surv.get.element(pred.o, "predicted")
  o[["predicted.oob"]] <- .surv.get.element(pred.o, "predicted.oob")
  o[["survival"]] <- .surv.get.element(pred.o, "survival")
  o[["survival.oob"]] <- .surv.get.element(pred.o, "survival.oob")
  if (is.null(.surv.get.element(o, "time.interest"))) {
    o[["time.interest"]] <- .surv.get.element(pred.o, "time.interest")
  }
  o[["forest"]] <- list(yvar = .surv.get.element(o, "yvar"),
                         xvar = .surv.get.element(o, "xvar"))
  o
}
## use imputed values stored on the current object, when available
.surv.apply.imputation <- function(o, yvar, xvar) {
  if (is.null(o$imputed.indv) || is.null(o$imputed.data) ||
      (is.null(yvar) && is.null(xvar))) {
    return(list(yvar = yvar, xvar = xvar))
  }
  n <- if (!is.null(yvar)) {
    nrow(yvar)
  }
  else {
    nrow(xvar)
  }
  if (is.null(n) || length(n) != 1L || n == 0L) {
    return(list(yvar = yvar, xvar = xvar))
  }
  idx <- o$imputed.indv
  if (is.logical(idx)) {
    if (length(idx) != n) {
      return(list(yvar = yvar, xvar = xvar))
    }
    idx <- which(idx)
  }
  else {
    idx <- suppressWarnings(as.integer(idx))
    idx <- idx[is.finite(idx) & idx >= 1L & idx <= n]
  }
  if (length(idx) == 0L) {
    return(list(yvar = yvar, xvar = xvar))
  }
  imp <- data.frame(o$imputed.data, check.names = FALSE)
  if (nrow(imp) == length(idx)) {
    imp.row <- seq_along(idx)
  }
  else if (nrow(imp) == n) {
    imp.row <- idx
  }
  else {
    return(list(yvar = yvar, xvar = xvar))
  }
  ## Survival imputation stores time and status in the first two columns.
  if (!is.null(yvar) && ncol(yvar) >= 2L && ncol(imp) >= 2L) {
    yvar[idx, 1:2] <- imp[imp.row, 1:2, drop = FALSE]
  }
  if (!is.null(xvar)) {
    nx <- ncol(xvar)
    if (nx > 0L && !is.null(yvar) && ncol(imp) >= 2L + nx) {
      x.col <- 2L + seq_len(nx)
    }
    else if (nx > 0L && ncol(imp) >= nx) {
      ## Accommodate prediction objects whose imputation matrix contains
      ## covariates only.
      x.col <- ncol(imp) - nx + seq_len(nx)
    }
    else {
      x.col <- integer(0L)
    }
    if (length(x.col) == nx && nx > 0L) {
      xvar[idx, ] <- imp[imp.row, x.col, drop = FALSE]
    }
  }
  list(yvar = yvar, xvar = xvar)
}
## coerce a survival prediction object to a case-by-time matrix
.surv.prediction.matrix <- function(x, ncase, ntime, name) {
  if (is.null(x)) {
    stop("survival predictions are unavailable in the supplied object.")
  }
  if (is.null(dim(x))) {
    if (ncase == 1L && length(x) == ntime) {
      x <- matrix(x, nrow = 1L)
    }
    else if (ntime == 1L && length(x) == ncase) {
      x <- matrix(x, ncol = 1L)
    }
    else {
      stop(name, " does not have dimensions compatible with the prediction object.")
    }
  }
  else {
    x <- as.matrix(x)
  }
  if (nrow(x) == ncase && ncol(x) == ntime) {
    return(x)
  }
  if (nrow(x) == ntime && ncol(x) == ncase) {
    return(t(x))
  }
  stop(name, " does not have dimensions compatible with the prediction object.")
}
## evaluate case-by-time step functions at requested times
.surv.step.matrix <- function(x, time.grid, times, initial = 1) {
  time.grid <- as.numeric(time.grid)
  times <- as.numeric(times)
  if (length(time.grid) == 0L || is.unsorted(time.grid, strictly = FALSE)) {
    stop("the survival time grid must be nondecreasing.")
  }
  idx <- findInterval(times, time.grid)
  out <- matrix(initial, nrow = nrow(x), ncol = length(times))
  use <- idx > 0L
  if (any(use)) {
    out[, use] <- x[, idx[use], drop = FALSE]
  }
  out
}
## evaluate each row of a case-by-time step function at its own time
.surv.step.row <- function(x, time.grid, times,
                           type = c("right", "left"), initial = 1) {
  type <- match.arg(type)
  times <- as.numeric(times)
  if (length(times) != nrow(x)) {
    stop("row-specific times do not match the number of predicted cases.")
  }
  out <- rep(NA_real_, length(times))
  ok <- is.finite(times)
  if (!any(ok)) {
    return(out)
  }
  idx <- if (type == "left") {
    findInterval(times[ok], time.grid, left.open = TRUE)
  }
  else {
    findInterval(times[ok], time.grid)
  }
  out[ok] <- initial
  use <- idx > 0L
  if (any(use)) {
    row <- which(ok)[use]
    out[row] <- x[cbind(row, idx[use])]
  }
  out
}
## create the event-information list used by the existing plotting code
.surv.event.info <- function(yvar, times) {
  if (is.null(yvar)) {
    return(list(event = NULL, event.type = NULL, cens = NULL,
                time.interest = times, time = NULL, r.dim = 0L))
  }
  if (ncol(yvar) != 2L) {
    stop("survival outcomes must contain exactly two columns: time and censoring.")
  }
  time <- as.numeric(yvar[, 1L])
  cens <- as.numeric(yvar[, 2L])
  valid.cens <- is.na(cens) |
    (is.finite(cens) & cens >= 0 & floor(cens) == cens)
  if (!all(valid.cens)) {
    stop("for survival families censoring variable must be coded as a non-negative integer")
  }
  event <- na.omit(cens)[na.omit(cens) > 0]
  event.type <- sort(unique(event))
  list(event = event, event.type = event.type, cens = cens,
       time.interest = times, time = time, r.dim = 2L)
}
## validate and normalize requested confidence level
.surv.conf.level <- function(conf.int) {
  if (length(conf.int) != 1L || is.na(conf.int)) {
    stop("'conf.int' must be FALSE, TRUE, or a single confidence level.")
  }
  if (is.logical(conf.int)) {
    return(if (conf.int) 0.95 else NULL)
  }
  conf.level <- suppressWarnings(as.numeric(conf.int))
  if (!is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("a numeric 'conf.int' must be strictly between zero and one.")
  }
  conf.level
}
## extract grow and evaluation data without mixing their outcomes
.get.survival.performance.data <- function(o, subset = NULL, times = NULL) {
  if (is.null(o)) {
    return(NULL)
  }
  if (o$family != "surv") {
    stop("this function only supports right-censored survival settings")
  }
  object.type <- .surv.object.type(o, allow.forest = TRUE)
  if (is.null(object.type)) {
    stop("This function only works for objects of class `(rfsrc, grow)', `(rfsrc, forest)' or `(rfsrc, predict)'")
  }
  if (object.type == "forest") {
    o <- .surv.prepare.forest(o)
  }
  time.grid <- as.numeric(.surv.get.element(o, "time.interest"))
  if (length(time.grid) == 0L || any(!is.finite(time.grid)) ||
      is.unsorted(time.grid, strictly = FALSE)) {
    stop("the supplied object does not contain a valid survival time grid.")
  }
  if (is.null(times)) {
    times <- time.grid
  }
  else {
    times <- as.numeric(times)
    times <- sort(unique(times[is.finite(times)]))
    if (length(times) == 0L) {
      stop("'times' does not contain any finite values.")
    }
  }
  ## The full grow data always define the censoring distribution. For a
  ## prediction object, never substitute test outcomes when grow outcomes are
  ## unavailable: scoring must stop rather than estimate G from the test set.
  ## Exact extraction is essential for stripped rfsrc.fast objects, which keep
  ## xvar.names but intentionally omit xvar when forest = FALSE.
  forest <- .surv.get.element(o, "forest")
  forest.yvar <- .surv.get.element(forest, "yvar")
  forest.xvar <- .surv.get.element(forest, "xvar")
  object.yvar <- .surv.get.element(o, "yvar")
  object.xvar <- .surv.get.element(o, "xvar")
  if (object.type == "predict") {
    grow.yvar <- forest.yvar
    grow.xvar <- forest.xvar
  }
  else {
    grow.yvar <- if (!is.null(forest.yvar)) forest.yvar else object.yvar
    grow.xvar <- if (!is.null(forest.xvar)) forest.xvar else object.xvar
  }
  if (!is.null(grow.yvar) &&
      (is.null(dim(grow.yvar)) || ncol(grow.yvar) != 2L)) {
    stop("grow outcomes must contain exactly two columns: time and censoring.")
  }
  if (!is.null(grow.xvar) && is.null(dim(grow.xvar))) {
    stop("grow covariates must be stored in a matrix or data frame.")
  }
  ## Evaluation data and predictions depend on the object type. Exact
  ## extraction also prevents predicted from partially matching predicted.oob,
  ## or survival from partially matching survival.oob, in reduced objects.
  predicted <- .surv.get.element(o, "predicted")
  predicted.oob <- .surv.get.element(o, "predicted.oob")
  survival.full <- .surv.get.element(o, "survival")
  survival.oob <- .surv.get.element(o, "survival.oob")
  if (object.type == "predict") {
    eval.yvar <- object.yvar
    eval.xvar <- object.xvar
    mort <- as.numeric(predicted)
    survival <- survival.full
    prediction <- "test"
  }
  else {
    eval.yvar <- object.yvar
    eval.xvar <- object.xvar
    if (!is.null(predicted.oob) && !is.null(survival.oob)) {
      mort <- as.numeric(predicted.oob)
      survival <- survival.oob
      prediction <- "oob"
    }
    else {
      mort <- as.numeric(predicted)
      survival <- survival.full
      prediction <- "ensemble"
    }
  }
  if (length(mort) == 0L) {
    stop("mortality predictions are unavailable in the supplied object.")
  }
  ncase <- length(mort)
  survival <- .surv.prediction.matrix(survival, ncase, length(time.grid),
                                      "survival predictions")
  if (!is.null(eval.yvar)) {
    if (is.null(dim(eval.yvar)) || ncol(eval.yvar) != 2L) {
      stop("evaluation outcomes must contain exactly two columns: time and censoring.")
    }
    if (nrow(eval.yvar) != ncase) {
      stop("evaluation outcomes do not align with the predicted cases.")
    }
  }
  if (!is.null(eval.xvar)) {
    if (is.null(dim(eval.xvar))) {
      stop("evaluation covariates must be stored in a matrix or data frame.")
    }
    if (nrow(eval.xvar) != ncase) {
      stop("evaluation covariates do not align with the predicted cases.")
    }
  }
  ## Imputation information on a predict object applies to prediction data;
  ## on a grow/forest object it applies to both grow and evaluation data.
  eval.imp <- .surv.apply.imputation(o, eval.yvar, eval.xvar)
  eval.yvar <- eval.imp$yvar
  eval.xvar <- eval.imp$xvar
  if (object.type != "predict") {
    grow.imp <- .surv.apply.imputation(o, grow.yvar, grow.xvar)
    grow.yvar <- grow.imp$yvar
    grow.xvar <- grow.imp$xvar
  }
  if (is.null(subset)) {
    subset <- seq_len(ncase)
  }
  else {
    if (is.logical(subset)) {
      if (length(subset) != ncase) {
        stop("a logical 'subset' must have length equal to the number of predicted cases.")
      }
      subset <- which(subset)
    }
    else if (!is.numeric(subset)) {
      stop("'subset' must be a numeric or logical vector.")
    }
    subset <- as.numeric(subset)
    subset <- unique(subset[is.finite(subset) & subset == floor(subset) &
                              subset >= 1L & subset <= ncase])
    subset <- as.integer(subset)
    if (length(subset) == 0L) {
      stop("'subset' not set properly.")
    }
  }
  survival <- survival[subset, , drop = FALSE]
  survival <- .surv.step.matrix(survival, time.grid, times, initial = 1)
  mort <- mort[subset]
  if (!is.null(eval.yvar)) {
    eval.yvar <- eval.yvar[subset, , drop = FALSE]
  }
  if (!is.null(eval.xvar)) {
    eval.xvar <- eval.xvar[subset, , drop = FALSE]
  }
  list(object = o,
       object.type = object.type,
       prediction = prediction,
       subset = subset,
       time = times,
       prediction.time = time.grid,
       survival = survival,
       mort = mort,
       grow.yvar = grow.yvar,
       grow.xvar = grow.xvar,
       eval.yvar = eval.yvar,
       eval.xvar = eval.xvar,
       grow.event.info = .surv.event.info(grow.yvar, times),
       eval.event.info = .surv.event.info(eval.yvar, times))
}
## Nelson-Aalen survival estimator for the evaluation outcomes
.surv.aalen <- function(time, status, times) {
  ok <- is.finite(time) & is.finite(status)
  time <- as.numeric(time[ok])
  status <- as.numeric(status[ok])
  if (length(time) == 0L) {
    return(NULL)
  }
  event.time <- sort(unique(time[status != 0]))
  if (length(event.time) == 0L) {
    return(rep(1, length(times)))
  }
  increment <- vapply(event.time, function(tt) {
    y <- sum(time >= tt)
    d <- sum(time == tt & status != 0)
    if (y > 0) d / y else 0
  }, numeric(1L))
  surv <- exp(-cumsum(increment))
  idx <- findInterval(times, event.time)
  out <- rep(1, length(times))
  use <- idx > 0L
  out[use] <- surv[idx[use]]
  out
}
## fit the grow-data censoring distribution and evaluate it on eval cases
## G.time is a length-m vector for KM and an n-by-m matrix for an RF model.
.surv.censoring <- function(dat, cens.model) {
  cens.model <- match.arg(cens.model, c("km", "rfsrc"))
  grow.time <- as.numeric(dat$grow.yvar[, 1L])
  grow.status <- as.numeric(dat$grow.yvar[, 2L])
  eval.time <- as.numeric(dat$eval.yvar[, 1L])
  n <- nrow(dat$survival)
  m <- length(dat$time)
  grow.complete <- is.finite(grow.time) & is.finite(grow.status)
  if (!any(grow.complete)) {
    stop("no non-missing grow outcomes are available for estimating censoring.")
  }
  has.censoring <- any(grow.complete & grow.status == 0)
  if (!has.censoring) {
    g.minus <- rep(1, n)
    if (cens.model == "km") {
      g.time <- rep(1, m)
      cens.dist <- g.time
    }
    else {
      g.time <- matrix(1, nrow = n, ncol = m)
      cens.dist <- t(g.time)
    }
    return(list(G.time = g.time,
                G.minus = g.minus,
                cens.dist = cens.dist,
                fit = NULL))
  }
  if (cens.model == "km") {
    fit <- km_censor_fit(grow.time[grow.complete],
                         as.integer(grow.status[grow.complete] != 0))
    g.time <- .km.censor.predict(fit$time, fit$G, dat$time,
                                 type = "right")
    g.minus <- .km.censor.predict(fit$time, fit$G, eval.time,
                                  type = "left")
    return(list(G.time = g.time,
                G.minus = g.minus,
                cens.dist = g.time,
                fit = fit))
  }
  if (is.null(dat$grow.xvar) || is.null(dat$eval.xvar)) {
    stop("'cens.model = \"rfsrc\"' requires stored grow and evaluation ",
         "covariates. They are unavailable in the supplied object. For an ",
         "rfsrc.fast fit, use cens.model = \"km\" or refit with ",
         "forest = TRUE.")
  }
  if (nrow(dat$grow.xvar) != nrow(dat$grow.yvar)) {
    stop("grow covariates do not align with grow outcomes.")
  }
  ## Missing grow outcomes cannot train the censoring forest. Predictor
  ## missingness, if present, is handled by RF-SRC in the usual way.
  grow.use <- which(grow.complete)
  cens.dta <- data.frame(time = grow.time[grow.use],
                         cens = 1 * (grow.status[grow.use] == 0),
                         dat$grow.xvar[grow.use, , drop = FALSE])
  cens.o <- rfsrc(Surv(time, cens) ~ ., cens.dta,
                  ntree = 50,
                  nsplit = 1,
                  splitrule = "random",
                  nodesize = set.nodesize(nrow(cens.dta),
                                           ncol(dat$grow.xvar)),
                  perf.type = "none")
  cens.pred <- predict(cens.o, newdata = dat$eval.xvar)$survival
  cens.pred <- .surv.prediction.matrix(cens.pred, n,
                                       length(cens.o$time.interest),
                                       "censoring survival predictions")
  g.time <- .surv.step.matrix(cens.pred, cens.o$time.interest,
                              dat$time, initial = 1)
  g.minus <- .surv.step.row(cens.pred, cens.o$time.interest,
                            eval.time, type = "left", initial = 1)
  list(G.time = g.time,
       G.minus = g.minus,
       cens.dist = t(g.time),
       fit = cens.o)
}
## Brier score and empirical conditional asymptotic variance
.surv.brier <- function(dat, censoring, conf.level, papply, keep.matrix) {
  time <- as.numeric(dat$eval.yvar[, 1L])
  status <- as.numeric(dat$eval.yvar[, 2L])
  times <- dat$time
  survival <- dat$survival
  n <- nrow(survival)
  m <- ncol(survival)
  ans <- papply(seq_len(n), function(i) {
    loss <- rep(NA_real_, m)
    bad.g <- rep(FALSE, m)
    if (!is.finite(time[i]) || !is.finite(status[i])) {
      return(list(loss = loss, bad.g = bad.g))
    }
    case <- time[i] <= times & status[i] != 0
    control <- time[i] > times
    censored <- time[i] <= times & status[i] == 0
    ## Censored observations contribute zero after their censoring time.
    loss[censored] <- 0
    if (any(case)) {
      jj <- which(case)
      good.s <- is.finite(survival[i, jj])
      good.g <- is.finite(censoring$G.minus[i]) &&
        censoring$G.minus[i] > 0
      if (!good.g) {
        ## A zero censoring probability invalidates only horizons at which
        ## this case otherwise has an evaluable prediction.
        bad.g[jj[good.s]] <- TRUE
      }
      else {
        loss[jj[good.s]] <- survival[i, jj[good.s]]^2 /
          censoring$G.minus[i]
      }
    }
    if (any(control)) {
      jj <- which(control)
      g <- if (is.matrix(censoring$G.time)) {
        censoring$G.time[i, jj]
      }
      else {
        censoring$G.time[jj]
      }
      good.g <- is.finite(g) & g > 0
      good.s <- is.finite(survival[i, jj])
      bad.g[jj[good.s & !good.g]] <- TRUE
      use <- good.g & good.s
      loss[jj[use]] <- (1 - survival[i, jj[use]])^2 / g[use]
    }
    list(loss = loss, bad.g = bad.g)
  })
  brier.matx <- do.call(rbind, lapply(ans, function(z) z$loss))
  bad.g <- do.call(rbind, lapply(ans, function(z) z$bad.g))
  ## A horizon is unsupported when a required censoring probability is zero.
  unsupported <- colSums(bad.g, na.rm = TRUE) > 0L
  if (any(unsupported)) {
    brier.matx[, unsupported] <- NA_real_
  }
  point <- vapply(seq_len(m), function(j) {
    z <- brier.matx[, j]
    z <- z[is.finite(z)]
    if (length(z) > 0L) mean(z) else NA_real_
  }, numeric(1L))
  n.eval <- vapply(seq_len(m), function(j) {
    sum(is.finite(brier.matx[, j]))
  }, integer(1L))
  brier.score <- data.frame(time = times, brier.score = point)
  if (!is.null(conf.level)) {
    std.err <- vapply(seq_len(m), function(j) {
      z <- brier.matx[, j]
      z <- z[is.finite(z)]
      nz <- length(z)
      if (nz > 1L) {
        influence <- z - mean(z)
        sqrt(sum(influence^2) / (nz * (nz - 1L)))
      }
      else {
        NA_real_
      }
    }, numeric(1L))
    crit <- qnorm((1 + conf.level) / 2)
    brier.score$std.err <- std.err
    brier.score$lower <- point - crit * std.err
    brier.score$upper <- point + crit * std.err
    brier.score$n.eval <- n.eval
  }
  crps <- trapz(times, point)
  time.max <- max(times, na.rm = TRUE)
  crps.std <- if (is.finite(time.max) && time.max > 0) {
    crps / time.max
  }
  else {
    NA_real_
  }
  list(brier.matx = if (keep.matrix) brier.matx else NULL,
       brier.score = brier.score,
       crps = crps,
       crps.std = crps.std,
       n.eval = n.eval)
}
## weighted case/control placements for cumulative/dynamic AUC
.surv.auct.placements <- function(case.risk, control.risk,
                                   case.weight, control.weight) {
  if (length(case.risk) == 0L || length(control.risk) == 0L ||
      length(case.risk) != length(case.weight) ||
      length(control.risk) != length(control.weight)) {
    stop("case and control risks must have matching nonempty weight vectors.")
  }
  if (any(!is.finite(case.risk)) || any(!is.finite(control.risk)) ||
      any(!is.finite(case.weight) | case.weight <= 0) ||
      any(!is.finite(control.weight) | control.weight <= 0)) {
    stop("AUCT risks must be finite and AUCT weights must be finite and positive.")
  }
  ## Scale before normalization to avoid overflow when IPCW weights are
  ## highly concentrated. The common scale does not change the AUC.
  p <- case.weight / max(case.weight)
  p <- pmax(p, .Machine$double.xmin)
  p <- p / sum(p)
  q <- control.weight / max(control.weight)
  q <- pmax(q, .Machine$double.xmin)
  q <- q / sum(q)
  ## Control-weight distribution used for case placements.
  ord <- order(control.risk)
  risk <- control.risk[ord]
  weight <- q[ord]
  run <- rle(risk)
  run.end <- cumsum(run$lengths)
  run.start <- run.end - run$lengths + 1L
  control.unique <- run$values
  control.weight.group <- vapply(seq_along(run.start), function(k) {
    sum(weight[run.start[k]:run.end[k]])
  }, numeric(1L))
  control.cum <- cumsum(control.weight.group)
  idx <- findInterval(case.risk, control.unique)
  equal <- idx > 0L
  equal[equal] <- control.unique[idx[equal]] == case.risk[equal]
  lower.idx <- idx - as.integer(equal)
  lower <- rep(0, length(case.risk))
  use <- lower.idx > 0L
  lower[use] <- control.cum[lower.idx[use]]
  same <- rep(0, length(case.risk))
  same[equal] <- control.weight.group[idx[equal]]
  case.place <- lower + 0.5 * same
  ## Case-weight distribution used for control placements.
  ord <- order(case.risk)
  risk <- case.risk[ord]
  weight <- p[ord]
  run <- rle(risk)
  run.end <- cumsum(run$lengths)
  run.start <- run.end - run$lengths + 1L
  case.unique <- run$values
  case.weight.group <- vapply(seq_along(run.start), function(k) {
    sum(weight[run.start[k]:run.end[k]])
  }, numeric(1L))
  case.cum <- cumsum(case.weight.group)
  idx <- findInterval(control.risk, case.unique)
  equal <- idx > 0L
  equal[equal] <- case.unique[idx[equal]] == control.risk[equal]
  leq <- rep(0, length(control.risk))
  use <- idx > 0L
  leq[use] <- case.cum[idx[use]]
  same <- rep(0, length(control.risk))
  same[equal] <- case.weight.group[idx[equal]]
  control.place <- 1 - leq + 0.5 * same
  auct <- sum(p * case.place)
  list(auct = auct,
       case.place = case.place,
       control.place = control.place,
       case.weight = p,
       control.weight = q)
}
## stratified delete-one jackknife for weighted cumulative/dynamic AUC
.surv.auct.jackknife <- function(place) {
  p <- place$case.weight
  q <- place$control.weight
  n.case <- length(p)
  n.control <- length(q)
  n.case.eff <- if (n.case > 0L) 1 / sum(p^2) else NA_real_
  n.control.eff <- if (n.control > 0L) 1 / sum(q^2) else NA_real_
  max.case.weight <- if (n.case > 0L) max(p) else NA_real_
  max.control.weight <- if (n.control > 0L) max(q) else NA_real_
  ## Delete one member of a stratum and renormalize the remaining IPCW
  ## weights. Prefix and suffix sums avoid subtracting a dominant weighted
  ## placement from the full AUC, which can lose precision under heavy
  ## censoring. Each set of replicates is centered at its own replicate mean.
  jackknife.component <- function(weight, placement) {
    n <- length(weight)
    if (n <= 1L) {
      return(NA_real_)
    }
    weighted.placement <- weight * placement
    keep <- seq_len(n - 1L)
    weight.before <- c(0, cumsum(weight)[keep])
    weight.after <- rev(c(0, cumsum(rev(weight))[keep]))
    value.before <- c(0, cumsum(weighted.placement)[keep])
    value.after <- rev(c(0, cumsum(rev(weighted.placement))[keep]))
    remain.weight <- weight.before + weight.after
    remain.value <- value.before + value.after
    if (any(!is.finite(remain.weight)) ||
        any(remain.weight <= 0) ||
        any(!is.finite(remain.value))) {
      return(NA_real_)
    }
    delete.value <- remain.value / remain.weight
    if (any(!is.finite(delete.value))) {
      return(NA_real_)
    }
    delete.mean <- mean(delete.value)
    (n - 1) / n * sum((delete.value - delete.mean)^2)
  }
  var.case <- jackknife.component(p, place$case.place)
  var.control <- jackknife.component(q, place$control.place)
  std.err <- if (is.finite(var.case) && is.finite(var.control)) {
    sqrt(max(0, var.case + var.control))
  }
  else {
    NA_real_
  }
  list(std.err = std.err,
       var.case = var.case,
       var.control = var.control,
       n.case.eff = n.case.eff,
       n.control.eff = n.control.eff,
       max.case.weight = max.case.weight,
       max.control.weight = max.control.weight)
}
## cumulative/dynamic AUC and stratified delete-one jackknife variance
.surv.auct <- function(dat, censoring, conf.level, papply) {
  time <- as.numeric(dat$eval.yvar[, 1L])
  status <- as.numeric(dat$eval.yvar[, 2L])
  times <- dat$time
  risk <- 1 - dat$survival
  empty.result <- function(n.case, n.control) {
    c(auct = NA_real_, std.err = NA_real_,
      n.case = n.case, n.control = n.control,
      n.case.eff = NA_real_, n.control.eff = NA_real_,
      max.case.weight = NA_real_, max.control.weight = NA_real_)
  }
  ans <- papply(seq_along(times), function(j) {
    complete <- is.finite(time) & is.finite(status) & is.finite(risk[, j])
    case <- complete & time <= times[j] & status != 0
    control <- complete & time > times[j]
    n.case <- sum(case)
    n.control <- sum(control)
    if (n.case == 0L || n.control == 0L) {
      return(empty.result(n.case, n.control))
    }
    g.case <- censoring$G.minus[case]
    g.control <- if (is.matrix(censoring$G.time)) {
      censoring$G.time[control, j]
    }
    else {
      rep(censoring$G.time[j], n.control)
    }
    if (any(!is.finite(g.case) | g.case <= 0) ||
        any(!is.finite(g.control) | g.control <= 0)) {
      return(empty.result(n.case, n.control))
    }
    ## Only relative IPCW weights enter AUC. Scaling by the smallest
    ## censoring probability keeps the largest weight at one and avoids
    ## overflow in heavily censored settings.
    case.weight <- pmax(min(g.case) / g.case, .Machine$double.xmin)
    control.weight <- pmax(min(g.control) / g.control,
                           .Machine$double.xmin)
    if (any(!is.finite(case.weight) | case.weight <= 0) ||
        any(!is.finite(control.weight) | control.weight <= 0)) {
      return(empty.result(n.case, n.control))
    }
    place <- .surv.auct.placements(
      case.risk = risk[case, j],
      control.risk = risk[control, j],
      case.weight = case.weight,
      control.weight = control.weight
    )
    ## Weight diagnostics are useful even when confidence intervals are not
    ## requested. The jackknife calculation itself is only needed for an SE.
    p <- place$case.weight
    q <- place$control.weight
    n.case.eff <- 1 / sum(p^2)
    n.control.eff <- 1 / sum(q^2)
    max.case.weight <- max(p)
    max.control.weight <- max(q)
    std.err <- NA_real_
    if (!is.null(conf.level)) {
      jackknife <- .surv.auct.jackknife(place)
      std.err <- jackknife$std.err
    }
    c(auct = place$auct, std.err = std.err,
      n.case = n.case, n.control = n.control,
      n.case.eff = n.case.eff,
      n.control.eff = n.control.eff,
      max.case.weight = max.case.weight,
      max.control.weight = max.control.weight)
  })
  ans <- as.data.frame(do.call(rbind, ans))
  auct.score <- data.frame(time = times, auct = ans$auct)
  if (!is.null(conf.level)) {
    crit <- qnorm((1 + conf.level) / 2)
    auct.score$std.err <- ans$std.err
    auct.score$lower <- ans$auct - crit * ans$std.err
    auct.score$upper <- ans$auct + crit * ans$std.err
  }
  auct.score$n.case <- as.integer(ans$n.case)
  auct.score$n.control <- as.integer(ans$n.control)
  auct.score$n.case.eff <- ans$n.case.eff
  auct.score$n.control.eff <- ans$n.control.eff
  auct.score$max.case.weight <- ans$max.case.weight
  auct.score$max.control.weight <- ans$max.control.weight
  list(auct.score = auct.score)
}
## shared engine for Brier score and time-dependent AUC
.get.survival.performance <- function(o,
                                      subset = NULL,
                                      cens.model = c("km", "rfsrc"),
                                      papply = lapply,
                                      times = NULL,
                                      conf.int = FALSE,
                                      keep.matrix = TRUE,
                                      metrics = c("brier", "auct")) {
  cens.model <- match.arg(cens.model)
  metrics <- unique(match.arg(metrics, c("brier", "auct"),
                              several.ok = TRUE))
  if (!is.function(papply)) {
    stop("'papply' must be a function such as lapply or mclapply.")
  }
  if (length(keep.matrix) != 1L || is.na(keep.matrix) ||
      !is.logical(keep.matrix)) {
    stop("'keep.matrix' must be TRUE or FALSE.")
  }
  conf.level <- .surv.conf.level(conf.int)
  dat <- .get.survival.performance.data(o, subset = subset, times = times)
  if (is.null(dat)) {
    return(NULL)
  }
  if (is.null(dat$eval.yvar)) {
    stop("evaluation outcomes are required to calculate Brier score or time-dependent AUC.")
  }
  if (is.null(dat$grow.yvar)) {
    stop("grow outcomes are required to estimate the censoring distribution.")
  }
  eval.complete <- is.finite(as.numeric(dat$eval.yvar[, 1L])) &
    is.finite(as.numeric(dat$eval.yvar[, 2L]))
  if (!any(eval.complete)) {
    stop("no non-missing evaluation outcomes are available.")
  }
  censoring <- .surv.censoring(dat, cens.model)
  surv.aalen <- .surv.aalen(dat$eval.event.info$time,
                            dat$eval.event.info$cens,
                            dat$time)
  common <- list(cens.dist = censoring$cens.dist,
                 time = dat$time,
                 event.info = dat$grow.event.info,
                 grow.event.info = dat$grow.event.info,
                 eval.event.info = dat$eval.event.info,
                 test.event.info = dat$eval.event.info,
                 subset = dat$subset,
                 mort = dat$mort,
                 surv.aalen = surv.aalen,
                 surv.ensb = t(dat$survival),
                 prediction = dat$prediction,
                 cens.model = cens.model,
                 conf.level = conf.level)
  out <- list()
  if ("brier" %in% metrics) {
    brier <- .surv.brier(dat, censoring, conf.level, papply, keep.matrix)
    out$brier <- c(brier, common)
  }
  if ("auct" %in% metrics) {
    auct <- .surv.auct(dat, censoring, conf.level, papply)
    out$auct <- c(auct, common)
  }
  out
}
## Brier score helper
get.brier.survival <- function(o,
                               subset = NULL,
                               cens.model = c("km", "rfsrc"),
                               papply = lapply,
                               times = NULL,
                               conf.int = FALSE,
                               keep.matrix = TRUE) {
  ans <- .get.survival.performance(o,
                                   subset = subset,
                                   cens.model = cens.model,
                                   papply = papply,
                                   times = times,
                                   conf.int = conf.int,
                                   keep.matrix = keep.matrix,
                                   metrics = "brier")
  if (is.null(ans)) NULL else ans$brier
}
## cumulative/dynamic time-dependent AUC helper
get.auct.survival <- function(o,
                              subset = NULL,
                              cens.model = c("km", "rfsrc"),
                              papply = lapply,
                              times = NULL,
                              conf.int = FALSE) {
  ans <- .get.survival.performance(o,
                                   subset = subset,
                                   cens.model = cens.model,
                                   papply = papply,
                                   times = times,
                                   conf.int = conf.int,
                                   keep.matrix = FALSE,
                                   metrics = "auct")
  if (is.null(ans)) NULL else ans$auct
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
## Generic step evaluation for the censoring distribution.
## The stored values G are post-step values at time_knots.
.km.censor.predict <- function(time_knots, G, t_new,
                               type = c("left", "right")) {
  type <- match.arg(type)
  t_new <- as.numeric(t_new)
  out <- rep(NA_real_, length(t_new))
  ok <- !is.na(t_new)
  if (!any(ok)) return(out)
  if (length(time_knots) == 0L) {
    out[ok] <- 1.0
    return(out)
  }
  idx <- if (type == "left") {
    findInterval(t_new[ok], time_knots, left.open = TRUE)
  }
  else {
    findInterval(t_new[ok], time_knots)
  }
  out[ok] <- c(1.0, G)[1L + idx]
  out
}
## Generic left-limit step evaluation retained for the Uno helpers.
uno_Ghat_minus_predict <- function(time_knots, G, t_new) {
  .km.censor.predict(time_knots, G, t_new, type = "left")
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
