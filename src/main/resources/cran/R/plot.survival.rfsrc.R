plot.survival.rfsrc <- function (x,
                                 show.plots = TRUE,
                                 subset, collapse = FALSE,
                                 cens.model = c("km", "rfsrc"),
                                 ...)
{
  ## incoming parameter checks
  if (is.null(x)) {
    stop("object x is empty!")
  }
  if (sum(inherits(x, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rfsrc", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("This function only works for objects of class `(rfsrc, grow)' or '(rfsrc, predict)'")
  }
  if (x$family != "surv") {
    stop("this function only supports right-censored survival settings")
  }
  ## acquire brier score, censoring distribution and other useful quantities
  brier.obj <- get.brier.survival(x, subset, cens.model)
  brier.matx <- brier.obj$brier.matx
  brier.score <- brier.obj$brier.score$brier.score
  mort <- brier.obj$mort
  surv.ensb <- brier.obj$surv.ensb
  surv.aalen <- brier.obj$surv.aalen
  event.info <- brier.obj$event.info
  test.event.info <- brier.obj$test.event.info
  subset <- brier.obj$subset
  ## brier processing
  mort.perc   <- c(min(mort, na.rm = TRUE) - 1e-5, quantile(mort, (1:4)/4, na.rm = TRUE))
  brier.score <- data.frame(cbind(do.call(cbind, lapply(1:4, function(k) {
    mort.pt <- (mort > mort.perc[k]) & (mort <= mort.perc[k+1])
    colMeans(brier.matx[mort.pt,, drop = FALSE], na.rm = TRUE)
  })), brier.score))
  colnames(brier.score) <- c("bs.q25", "bs.q50", "bs.q75", "bs.q100", "bs.all")
  ## crps processing
  crps <- do.call(cbind, lapply(1:5, function(k) {
    lapply(1:length(event.info$time.interest), function(j) {
      trapz(event.info$time.interest[1:j], brier.score[1:j, k]) / diff(range(event.info$time.interest[1:j]))
    })
  }))
  colnames(crps) <- c("crps.q25", "crps.q50", "crps.q75", "crps.q100", "crps.all")
  ## labels and titles
  if (is.null(x$predicted.oob)) {
    ylab.1 <- "Survival"
    ylab.2 <- "Brier Score"
    ylab.3 <- "CRPS"
    ylab.4 <- "Mortality vs Time"
  }
  else {
    ylab.1 <- "OOB Survival"
    ylab.2 <- "OOB Brier"
    ylab.3 <- "OOB CRPS"
    ylab.4 <- "OOB Mortality vs Time"
  }
  ## mean ensemble survival
  surv.mean.ensb <- rowMeans(surv.ensb, na.rm = TRUE)
  ## collapse across the subset?
  if (collapse) {
    surv.ensb <- rowMeans(surv.ensb, na.rm = TRUE)
  }
  ## ------------------------------------------------------------
  ## 
  ## plots
  ##
  ## -------------------------------------------------------------
  ## should we display the plots?
  if (show.plots) {
    old.par <- par(no.readonly = TRUE)
    par(mfrow=c(2,2), cex = 1.0)
    ## ----survival plot----
    if (!collapse && length(subset) > 500) {
      r.pt <- sample(1:length(subset), 500, replace = FALSE)
      matplot(x$time.interest,
              surv.ensb[, r.pt, drop = FALSE],
              xlab = "Time",
              ylab = ylab.1,
              type = "l",
              col = 1,
              lty = 3, ...)
    }
    else {
      matplot(x$time.interest,
              surv.ensb,
              xlab = "Time",
              ylab = ylab.1,
              type = "l",
              col = 1,
              lty = 3, ...)
    }
    if (!is.null(surv.aalen)) {
      lines(event.info$time.interest, surv.aalen, lty = 1, col = 3, lwd = 3)
    }
    lines(event.info$time.interest, surv.mean.ensb, lty = 1, col = 2, lwd = 3)
    rug(event.info$time.interest, ticksize=-0.03)
    ## ----brier plot----
    matplot(event.info$time.interest, brier.score,
            xlab = "Time",
            ylab = ylab.2,
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
    matplot(event.info$time.interest, crps,
            xlab = "Time",
            ylab = ylab.3,
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
    ## ----mortality plot----
    if (!is.null(x$yvar)) {
      yvar <- x$yvar[subset,, drop = FALSE]
      colnames(yvar) <- c("time", "cens")
      plot(yvar$time, mort, xlab = "Time", ylab = ylab.4, type = "n", ...)
      if (x$n > 500) cex <- 0.5 else cex <- 0.75
      points(yvar$time[yvar$cens != 0], mort[yvar$cens != 0], pch = 16, col = 4, cex = cex)
      points(yvar$time[yvar$cens == 0], mort[yvar$cens == 0], pch = 16, cex = cex)
      if (sum(yvar$cens != 0) > 1)
        lines(supsmu(yvar$time[yvar$cens != 0][order(yvar$time[yvar$cens != 0])],
                     mort[yvar$cens != 0][order(yvar$time[yvar$cens != 0])]),
              lty = 3,
              col = 4,
              cex = cex)
      if (sum(yvar$cens == 0) > 1)
        lines(supsmu(yvar$time[yvar$cens == 0][order(yvar$time[yvar$cens == 0])],
                     mort[yvar$cens == 0][order(yvar$time[yvar$cens == 0])]),
              lty = 3,
              cex = cex)
      if (sum(yvar$cens != 0) > 1) {
        rug(yvar$time[yvar$cens != 0], ticksize=-0.03)
      }
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
  invisible(data.frame(time = event.info$time.interest, brier.score, crps))
}
plot.survival <- plot.survival.rfsrc
