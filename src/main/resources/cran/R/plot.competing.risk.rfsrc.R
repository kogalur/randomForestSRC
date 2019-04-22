plot.competing.risk.rfsrc <- function (x, plots.one.page = FALSE, ...) {
  ## Incoming parameter checks.  All are fatal.
  if (is.null(x)) {
    stop("object x is empty!")
  }
  if (sum(inherits(x, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rfsrc", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("This function only works for objects of class `(rfsrc, grow)' or '(rfsrc, predict)'.")
  }
  if (x$family != "surv-CR") {
    stop("this function only supports competing risk settings")
  }
  ## work-horse plotting function
  matPlot <- function(matx, ylab = "", legend = "", pos = 1) {
    m <- dim(cbind(matx))[2]
    if (m > 1) legend <- paste(legend, 1:m, "  ")
    matplot(x$time.interest, matx, xlab = "Time", ylab = ylab, type = "l",
            col = (1:m), lty = 1, lwd = 3)
    legend(c("topright", "bottomright")[pos], legend = legend, col = (1:m), lty = 1, lwd = 3)
  }
  ## save par settings
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  if (plots.one.page) par(mfrow = c(1,1)) else par(mfrow = c(2,2))
  ## acquire the estimators - use OOB whenever possible
  if (!is.null(x$chf.oob)) {
    cschf <- apply(x$chf.oob, c(2, 3), mean, na.rm = TRUE)
    cif <- apply(x$cif.oob, c(2, 3), mean, na.rm = TRUE)
  }
  else {
    cschf <- apply(x$chf, c(2, 3), mean, na.rm = TRUE)
    cif <- apply(x$cif, c(2, 3), mean, na.rm = TRUE)
  }
  cpc <- do.call(cbind, lapply(1:ncol(cif), function(j) {
    cif[, j] / (1 - rowSums(cif[, -j, drop = FALSE]))
  }))
  ## plot the results 
  matPlot(cschf, "Cause-Specific CHF", "CSCHF", pos = 2)
  matPlot(100 * cif, "Probability (%)", "CIF", 2)
  matPlot(100 * cpc, "Probability (%)", "CPC", 2)
}
plot.competing.risk <- plot.competing.risk.rfsrc
