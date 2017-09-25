\name{quantileReg}
\alias{quantileReg}
\alias{quantileReg.rfsrc}
\title{Quantile Regression Forests}
\description{
  Determines the conditional quantiles for a regression forest.  Applies to both
  univariate and multivariate forests.  Can be used for both training
  and testing purposes.
}
\usage{\method{quantileReg}{rfsrc}(obj, oob = TRUE, prob = (1:10) / 10, newdata = NULL)
}
\arguments{
  \item{obj}{A previously grown forest.}
  \item{oob}{Return OOB (out-of-bag) quantiles?  If false, in-bag values
  are returned.}
  \item{prob}{Target quantile probabilities.}
  \item{newdata}{Test data (optional) over which conditional quantiles
    are evaluated over.}
}
\details{
  Given a regression forest, or a multivariate forest with at least one
  regression outcome, returns the conditional quantiles for the target
  outcomes.  Also returns the conditional density, which can be used to
  calculate conditional moments, such as the mean and standard
  deviation.  
}
\value{
  Conditional quantiles and densitites are returned for each data point
  in the training set, or if a test set is provided, the values are
  returned for each data point in the test data.  If more than one
  target outcome is available, the returned object will be a list of
  length equal to the number of targets.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Meinshausen N. (2006) Quantile Regression Forests, \emph{Journal of
  Machine Learning Research}, 7:983--999.
}
\seealso{
  \command{\link{rfsrc}}
}
\examples{
\dontrun{
## ------------------------------------------------------------
## regression example
## ------------------------------------------------------------

## grow call, followed by quantileReg call
o <- rfsrc(mpg ~ ., mtcars)
qo <- quantileReg(o, prob = c(0.25, 0.5, 0.75))

## calculate the conditional mean, compare to OOB predicted value
c.mean <- qo$density \%*\% qo$yunq
print(data.frame(c.mean = c.mean, pred.oob = o$predicted.oob))

## calculate conditional standard deviation
c.std <- sqrt(qo$density \%*\% qo$yunq^2 - c.mean ^ 2)
quant <- qo$quantile
colnames(quant) <- paste("q", 100 * qo$prob, sep = "")
print(data.frame(quant, c.std))


## ------------------------------------------------------------
## train/test regression example
## ------------------------------------------------------------

## grow call, followed by quantileReg call
o <- rfsrc(mpg ~ ., mtcars[1:20,])
qo <- quantileReg(o, newdata = mtcars[-(1:20),], prob = c(0.25, 0.5, 0.75))

## calculate test set conditional mean and standard deviation
c.mean <- qo$density \%*\% qo$yunq
c.std <- sqrt(qo$density \%*\% qo$yunq^2 - c.mean ^ 2)
quant <- qo$quant
colnames(quant) <- paste("q", 100 * qo$prob, sep = "")
print(data.frame(quant, c.mean, c.std))


## ------------------------------------------------------------
## multivariate mixed outcomes example
## ------------------------------------------------------------

dta <- mtcars
dta$cyl <- factor(dta$cyl)
dta$carb <- factor(dta$carb, ordered = TRUE)
o <- rfsrc(cbind(carb, mpg, cyl, disp) ~., data = dta)
qo <- quantileReg(o)

print(head(qo$mpg$quant))
print(head(qo$disp$quant))


## ------------------------------------------------------------
## quantile regression plot for Boston Housing data
## ------------------------------------------------------------

if (library("mlbench", logical.return = TRUE)) {

  ## apply quantile regression to Boston Housing data
  data(BostonHousing)
  o <- rfsrc(medv ~ ., BostonHousing)
  qo <- quantileReg(o, prob = c(0.25, 0.5, 0.75))

  ## quantile data for plotting 
  quant.dat <- qo$quant
  y <- o$yvar

  ## quantile regression plot
  plot(range(y), range(quant.dat), xlab = "y",
       ylab = ".25-.75 Quantiles", type = "n")
  jitter.y <- jitter(y, 10)
  points(jitter.y, quant.dat[, 2], pch = 15, col = 4, cex = 0.75)
  segments(jitter.y, quant.dat[, 2], jitter.y, quant.dat[, 1], col = "grey")
  segments(jitter.y, quant.dat[, 2], jitter.y, quant.dat[, 3], col = "grey")
  points(jitter.y, quant.dat[, 1], pch = "-", cex = 1)
  points(jitter.y, quant.dat[, 3], pch = "-", cex = 1)
  abline(0, 1, lty = 2, col = 2)

}



}}
\keyword{quantile regression forests}

