plot.survival.rfsrc <- function (x,
                                 plots.one.page = TRUE,
                                 show.plots = TRUE,
                                 subset, collapse = FALSE,
                                 cens.model = c("km", "rfsrc"),
                                 ...)
{
  ## Incoming parameter checks.  All are fatal.
  if (is.null(x)) {
    stop("object x is empty!")
  }
  if (sum(inherits(x, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rfsrc", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("This function only works for objects of class `(rfsrc, grow)' or '(rfsrc, predict)'.")
  }
  if (x$family != "surv") {
    stop("this function only supports right-censored survival settings")
  }
  ## predict object does not contain OOB values
  if (sum(inherits(x, c("rfsrc", "predict"), TRUE) == c(1, 2)) == 2) {
    pred.flag <- TRUE
  }
  else {
    pred.flag <- FALSE
  }
  ## grow objects under non-standard bootstrapping are OOB devoid
  ## treat the object as if it were predict
  if (is.null(x$predicted.oob)) {
    pred.flag <- TRUE
  }
  ## verify the cens.model option
  cens.model <- match.arg(cens.model, c("km", "rfsrc"))
  ## use imputed missing time or censoring indicators
  if (!is.null(x$yvar) && !is.null(x$imputed.indv)) {
    x$yvar[x$imputed.indv, ]=x$imputed.data[, 1:2]
  }
  ## get the event data
  event.info <- get.event.info(x)
  ## Process the subsetted index
  ## Assumes the entire data set is to be used if not specified
  if (missing(subset)) {
    subset <- 1:x$n
    subset.provided <- FALSE
  }
    else {
      ## convert the user specified subset into a usable form
      if (is.logical(subset)) subset <- which(subset)
      subset <- unique(subset[subset >= 1 & subset <= x$n])
      show.plots <- subset.provided <- TRUE
      if (length(subset) == 0) {
        stop("'subset' not set properly.")
      }
    }
  ## no point in producing plots if sample size is too small
  if (!pred.flag && !subset.provided && (x$n < 2 | x$ndead < 1)) {
    stop("sample size or number of deaths is too small for meaningful analysis")
  }
  ## use OOB values if available
  if (is.null(x$predicted.oob)) {
    mort <- x$predicted[subset]
    surv.ensb <- t(x$survival[subset,, drop = FALSE])
    chf.ensb <- x$chf[subset,, drop = FALSE]
    y.lab <- "Mortality"
    title.1 <- "Survival"
    title.2 <- "Cumulative Hazard"
    title.3 <- "Mortality vs Time"
  }
  else {
    mort <- x$predicted.oob[subset]
    surv.ensb <- t(x$survival.oob[subset,, drop = FALSE])
    chf.ensb <- x$chf.oob[subset,, drop = FALSE]
    y.lab <- "OOB Mortality"
    title.1 <- "OOB Survival"
    title.2 <- "OOB Cumulative Hazard"
    title.3 <- "OOB Mortality vs Time"
  }
  ## mean ensemble survival
  if (!subset.provided) {
    surv.mean.ensb <- rowMeans(surv.ensb, na.rm = TRUE)
  }
  ## collapse across the subset?
  if (subset.provided && collapse) {
    surv.ensb <- rowMeans(surv.ensb, na.rm = TRUE)
    chf.ensb <- rbind(colMeans(chf.ensb, na.rm = TRUE))
  }
  ## ------------------------------------------------------------
  ## 
  ## survival calculations
  ##
  ## -------------------------------------------------------------
  if (!pred.flag && !subset.provided) {
    ## KM estimator
    km.obj <- do.call(rbind, mclapply(1:length(event.info$time.interest), function(j) {
                c(sum(event.info$time >= event.info$time.interest[j], na.rm = TRUE),
                sum(event.info$time[event.info$cens != 0] == event.info$time.interest[j], na.rm = TRUE))
    }))
    Y <- km.obj[, 1]
    d <- km.obj[, 2]
    r <- d / (Y + 1 * (Y == 0))
    surv.aalen <- exp(-cumsum(r))
    ## Estimate the censoring distribution
    sIndex <- function(x,y) {sapply(1:length(y), function(j) {sum(x <= y[j])})}
    censTime <- sort(unique(event.info$time[event.info$cens == 0]))
    censTime.pt <- c(sIndex(censTime, event.info$time.interest))
    ## check to see if there are censoring cases
    if (length(censTime) > 0) {
      ## KM estimator for the censoring distribution
      if (cens.model == "km") {
        censModel.obj <- do.call(rbind, mclapply(1:length(censTime), function(j) {
             c(sum(event.info$time >= censTime[j], na.rm = TRUE),
               sum(event.info$time[event.info$cens == 0] == censTime[j], na.rm = TRUE))
        }))
        Y <- censModel.obj[, 1]
        d <- censModel.obj[, 2]
        r <- d / (Y + 1 * (Y == 0))
        cens.dist <- c(1, exp(-cumsum(r)))[1 + censTime.pt]
      }
      ## rfsrc estimator for the censoring distribution
      else {
        cens.dist <- rfsrc(Surv(tm, cens) ~ .,
                data.frame(tm = x$yvar[, 1], cens = 1 * (x$yvar[, 2] == 0), x$xvar),
                ntime = 250, ntree = 250, nsplit = 1)$survival.oob
        cens.dist <- t(cbind(1, cens.dist)[, 1 + censTime.pt])
      }
    }
    ## no censoring cases; assign a default distribution
    else {
      cens.dist <- rep(1, length(censTime.pt))
    }
    ## ------------------------------------------------------------
    ## 
    ## brier calculations
    ##
    ## -------------------------------------------------------------
    brier.obj <- do.call(rbind, mclapply(1:x$n, function(i) {
      tau <-  event.info$time
      event <- event.info$cens
      t.unq <- event.info$time.interest
      cens.pt <- sIndex(t.unq, tau[i])
      if (cens.model == "km") {
        c1 <- 1 * (tau[i] <= t.unq & event[i] != 0)/c(1, cens.dist)[1 + cens.pt]
        c2 <- 1 * (tau[i] > t.unq)/cens.dist
      }
      else {
        c1 <- 1 * (tau[i] <= t.unq & event[i] != 0)/c(1, cens.dist[, i])[1 + cens.pt]
        c2 <- 1 * (tau[i] > t.unq)/cens.dist[, i]
      }
      (1 * (tau[i] > t.unq) - surv.ensb[, i])^2 * (c1 + c2)
    }))
    brier.score <- matrix(NA, length(event.info$time.interest), 4)
    mort.perc   <- c(min(mort, na.rm = TRUE) - 1e-5, quantile(mort, (1:4)/4, na.rm = TRUE))
    brier.score <- do.call(cbind, lapply(1:4, function(k) {
      mort.pt <- (mort > mort.perc[k]) & (mort <= mort.perc[k+1])
      apply(brier.obj[mort.pt,, drop = FALSE], 2, mean, na.rm = TRUE)
    }))
    brier.score <- data.frame(brier.score, all = colMeans(brier.obj, na.rm = TRUE))
    colnames(brier.score) <- c("bs.q25", "bs.q50", "bs.q75", "bs.q100", "bs.all")
    ## crps - continuous rank probability score
    crps <- do.call(cbind, lapply(1:5, function(k) {
      mclapply(1:length(event.info$time.interest), function(j) {
        trapz(event.info$time.interest[1:j], brier.score[1:j, k]) / diff(range(event.info$time.interest[1:j]))
      })
    }))
    colnames(crps) <- c("crps.q25", "crps.q50", "crps.q75", "crps.q100", "crps.all")
  }
  ## ------------------------------------------------------------
  ## 
  ## plots
  ##
  ## -------------------------------------------------------------
  ## should we display the plots?
  if (show.plots) {
    old.par <- par(no.readonly = TRUE)
    ## plots on one page
    if (plots.one.page) {
      if (pred.flag && !subset.provided) {
        if (!is.null(x$yvar)) {
          ## survival/mortality only
          par(mfrow = c(1,2))
        }
          else {
            ## predict mode but no outcomes: survival only
            par(mfrow = c(1,1))
          }
      }
      ## grow mode
      else {
        if (!subset.provided) {
          par(mfrow = c(2, 2))
        }
        else {
          par(mfrow = c(1, 2))
        }
      }
    }
    ## plots on separate pages
    else {
      par(mfrow=c(1,1))
    }
    par(cex = 1.0)
    ## ----survival plot----
    if (!subset.provided && x$n > 500) {
      r.pt <- sample(1:x$n, 500, replace = FALSE)
      matplot(event.info$time.interest,
              surv.ensb[, r.pt],
              xlab = "Time",
              ylab = title.1,
              type = "l",
              col = 1,
              lty = 3, ...)
    }
    else {
      matplot(event.info$time.interest,
              surv.ensb,
              xlab = "Time",
              ylab = title.1,
              type = "l",
              col = 1,
              lty = 3, ...)
    }
    if (!pred.flag && !subset.provided) {
      lines(event.info$time.interest, surv.aalen, lty = 1, col = 3, lwd = 3)
    }
    if (!subset.provided) {
      lines(event.info$time.interest, surv.mean.ensb, lty = 1, col = 2, lwd = 3)
    }
    rug(event.info$time.interest, ticksize=-0.03)
    if (plots.one.page) {
      title(title.1, cex.main = 1.25)
    }
    ## ----CHF plot----
    if (subset.provided) {
      matplot(event.info$time.interest,
              t(chf.ensb),
              xlab = "Time",
              ylab = title.2,
              type = "l",
              col = 1,
              lty = 3, ...)
      if (plots.one.page) {
        title(title.2, cex.main = 1.25)
      }
    }
    ## ----brier plot----
    if (!pred.flag && !subset.provided) {
      matplot(event.info$time.interest, brier.score,
              xlab = "Time",
              ylab = "OOB Brier Score",
              type = "l",
              lwd  = c(rep(1, 4), 2),
              col  = c(rep(1, 4), 2),
              lty  = c(1:4, 1), ...)
      rng <- range(unlist(brier.score), na.rm = TRUE)
      abline(h = seq(rng[1], rng[2], length = 20), col = gray(.6), lty = 3, lwd = .85)
      point.x=round(length(event.info$time.interest)*c(3,4)/4)
      text(event.info$time.interest[point.x],brier.score[point.x,1],"0-25",col=4)
      text(event.info$time.interest[point.x],brier.score[point.x,2],"25-50",col=4)
      text(event.info$time.interest[point.x],brier.score[point.x,3],"50-75",col=4)
      text(event.info$time.interest[point.x],brier.score[point.x,4],"75-100",col=4)
      rug(event.info$time.interest,ticksize=0.03)
      if (plots.one.page) title("OOB Brier Score",cex.main = 1.25)
      matplot(event.info$time.interest, crps,
              xlab = "Time",
              ylab = "OOB CRPS",
              type = "l",
              lwd  = c(rep(1, 4), 2),
              col  = c(rep(1, 4), 2),
              lty  = c(1:4, 1), ...)
      rng <- range(unlist(crps), na.rm = TRUE)
      abline(h = seq(rng[1], rng[2], length = 20), col = gray(.6), lty = 3, lwd = .85)
      point.x=round(length(event.info$time.interest)*c(3,4)/4)
      text(event.info$time.interest[point.x],crps[point.x,1],"0-25",col=4)
      text(event.info$time.interest[point.x],crps[point.x,2],"25-50",col=4)
      text(event.info$time.interest[point.x],crps[point.x,3],"50-75",col=4)
      text(event.info$time.interest[point.x],crps[point.x,4],"75-100",col=4)
      rug(event.info$time.interest,ticksize=0.03)
      if (plots.one.page) title("OOB CRPS",cex.main = 1.25)
    }
    ## ----mortality plot----
    if (!subset.provided && !is.null(x$yvar)) {
      plot(event.info$time, mort, xlab = "Time", ylab = y.lab, type = "n", ...)
      if (plots.one.page) {
        title(title.3, cex.main = 1.25)
      }
      if (x$n > 500) cex <- 0.5 else cex <- 0.75
      points(event.info$time[event.info$cens != 0], mort[event.info$cens != 0], pch = 16, col = 4, cex = cex)
      points(event.info$time[event.info$cens == 0], mort[event.info$cens == 0], pch = 16, cex = cex)
      if (sum(event.info$cens != 0) > 1)
        lines(supsmu(event.info$time[event.info$cens != 0][order(event.info$time[event.info$cens != 0])],
                     mort[event.info$cens != 0][order(event.info$time[event.info$cens != 0])]),
              lty = 3,
              col = 4,
              cex = cex)
      if (sum(event.info$cens == 0) > 1)
        lines(supsmu(event.info$time[event.info$cens == 0][order(event.info$time[event.info$cens == 0])],
                     mort[event.info$cens == 0][order(event.info$time[event.info$cens == 0])]),
              lty = 3,
              cex = cex)
      rug(event.info$time.interest, ticksize=-0.03)
    }
    ## reset par
    par(old.par)
  }
  ## ------------------------------------------------------------
  ## 
  ## return invisible objects 
  ##
  ## -------------------------------------------------------------
  ## invisibly return the brier score
  if (!pred.flag && !subset.provided) {
    invisible(data.frame(time = event.info$time.interest, brier.score, crps))
  }
}
plot.survival <- plot.survival.rfsrc
## trapezoidal rule
trapz <- function (x, y) {
  idx = 2:length(x)
  return(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2)
}
