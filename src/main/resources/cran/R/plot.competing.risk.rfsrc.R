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
  if (plots.one.page) par(mfrow = c(1,1)) else par(mfrow = c(1,2))
  matPlot(apply(x$chf, c(2, 3), mean, na.rm = TRUE), "CHF", "CSCHF", pos = 2)
  matPlot(100 * apply(x$cif, c(2, 3), mean, na.rm = TRUE), "Probability (%)", "CIF", 2)
}
plot.competing.risk <- plot.competing.risk.rfsrc
