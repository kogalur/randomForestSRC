plot.quantreg.rfsrc <- function(x, prbL = .25, prbU = .75,
                    m.target = NULL, crps = TRUE, subset = NULL, ...) {
  ##--------------------------------------------------------------
  ##
  ## prelimary checks
  ##
  ##--------------------------------------------------------------
  ## does not apply to predict objects without y
  if (sum(grepl("predict", class(x))) > 0 && is.null(x$yvar)) {
    stop("no yvar present in quantreg predict object")
  }
  ## check probs
  if (prbL < 0 || prbU > 1) {
    stop("requested probabilities must lie in (0, 1)")
  }
  if (prbL >= prbU) {
    stop("prbL must be less than prbU")
  }
  prbM <- max(prbL, .5)
  if (prbM == prbL) {
    prbM <- (prbL + prbU) / 2
  }
  ##--------------------------------------------------------------
  ##
  ## subset assignment
  ##
  ##--------------------------------------------------------------
  if (is.null(subset)) {
    subset <- 1:x$n
  }
  else {
    if (is.logical(subset)) {
      subset <- which(subset)
    }
    subset <- subset[subset >=1 & subset <= x$n]
  }
  if (length(subset) == 0) {
    stop("requested subset analysis has subset that is empty")
  }
  ##--------------------------------------------------------------
  ##
  ## assemble the quantile data for plotting
  ##
  ##--------------------------------------------------------------
  quant.dat <- get.quantile(x, c(prbL, prbM, prbU), FALSE)
  ## we have a univariate quantile regression object
  if (length(quant.dat) == 1) {
    y.names <- names(quant.dat)[1]
    quant.dat <- quant.dat[[1]][subset,, drop = FALSE]
    y <- x$yvar[subset]
    if (crps) crps.dat <- get.quantile.crps(x, subset = subset)
  }
  ## the quantile regression object is multivariate
  else {
    if (is.null(m.target) || length(intersect(m.target, names(quant.dat))) == 0) {
      y.names <- names(quant.dat)[1]
      quant.dat <- quant.dat[[1]][subset,, drop = FALSE]
      y <- x$yvar[, y.names][subset]
      if (crps) crps.dat <- get.quantile.crps(x, subset = subset)[[1]]
    }
    else {
      y.names <- m.target
      quant.dat <- quant.dat[[y.names]][subset,, drop = FALSE]
      y <- x$yvar[, y.names][subset]
      if (crps) crps.dat <- get.quantile.crps(x, subset = subset)[[y.names]]
    }
  }
  ##--------------------------------------------------------------
  ##
  ## quantile regression plot
  ##
  ##--------------------------------------------------------------
  jitter.y <- jitter(y, 10)
  rng <- range(c(y, quant.dat, jitter.y))
  plot(rng, rng, xlab = y.names, ylab = "Target Quantiles", type = "n")
  points(jitter.y, quant.dat[, 2], pch = 15, col = 4, cex = 0.75)
  segments(jitter.y, quant.dat[, 2], jitter.y, quant.dat[, 1], col = "grey")
  segments(jitter.y, quant.dat[, 2], jitter.y, quant.dat[, 3], col = "grey")
  points(jitter.y, quant.dat[, 1], pch = "-", cex = 1)
  points(jitter.y, quant.dat[, 3], pch = "-", cex = 1)
  abline(0, 1, lty = 2, col = 2)
  ##--------------------------------------------------------------
  ##
  ## inset the CRPS
  ##
  ##--------------------------------------------------------------
  if (crps) {
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    u <- par("usr")
    v <- c(
      grconvertX(u[1:2], "user", "ndc"),
      grconvertY(u[3:4], "user", "ndc")
    )
    v <- c(v[1], (v[1]+v[2]) / 3,  (v[3]+v[4]) / (1.5), v[4])
    figO <- tryCatch({par(fig = v, new = TRUE, mar = c(0,0,0,0), mgp = c(0,0,0))}, error=function(ex){NULL})
    if (!is.null(figO)) {
      plot(crps.dat, yaxt = "n", col = 2, type = "l", lwd = 2, tck = .05, cex.axis = .75, ...)
      axis(4, cex.axis = .75, tck = .05)
      box()
    }
  }
}
plot.quantreg <- plot.quantreg.rfsrc
