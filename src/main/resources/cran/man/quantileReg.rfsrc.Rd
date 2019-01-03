\name{quantileReg}
\alias{quantileReg}
\alias{quantileReg.rfsrc}
\title{Quantile Regression Forests}
\description{
  Grows a univariate or multivariate quantile regression forest and returns
  its conditional quantile and density values.  Can be used for
  both training and testing purposes.
}
\usage{\method{quantileReg}{rfsrc}(formula, data, object, newdata,
  method = "forest", prob = NULL, prob.epsilon = NULL,
  oob = TRUE, fast = FALSE, maxn = 1e3, ...)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit.
    Must be specified unless \code{object} is given.}
  \item{data}{Data frame containing the y-outcome and x-variables in
    the model. Must be specified unless \code{object} is given.}
  \item{object}{(Optional) A previously grown quantile regression
    forest.}
  \item{method}{Method used to calculate quantiles.  Forest weighted
    averaging is used by default.  While this works well for standard
    data, consider using the Greenwald-Khanna algorithm for big data.
    The latter is specified by any one of the following: "gk", "GK",
    "G-K", "g-k".}
  \item{prob}{Target quantile probabilities when training.  If left unspecified,
    uses percentiles (1 through 99) for \code{method = "forest"}, and 
    for Greenwald-Khanna selects equally spaced percentiles optimized
    for accuracy (see below).}
  \item{prob.epsilon}{Greenwald-Khanna allowable error for quantile
    probabilities when training.}
  \item{newdata}{Test data (optional) over which conditional quantiles
    are evaluated over.}
  \item{oob}{Return OOB (out-of-bag) quantiles?  If false, in-bag values
    are returned.}
  \item{fast}{Use fast random forests, \code{rfsrcFast}, in place of
    \code{rfsrc}?  Improves speed but may be less accurate.}
  \item{maxn}{Maximum number of unique y training values used when
    calculating the conditional density.}
  \item{...}{Further arguments to be passed to the \code{rfsrc}
    function used for fitting the quantile regression forest.}  
}
\details{

  Grows a univariate or multivariate quantile regression forest using
  quantile regression splitting using the new splitrule
  \code{quantile.regr} based on the quantile loss function (often called
  the "check function").
  
  The default method for calculating quantiles is \code{method}="forest"
  which uses forest weights as in Meinshausen (2006).  However, because
  quantile regression splitting is used, and not mean-squared error
  splitting (as used by Meinshuasen), results may differ substantially
  from Meinshausen.  We believe quantile regression splitting will
  provide superior performance.

  While calculating quantiles using forest weights works well for
  standard data, a second approach, the Greenwald-Khanna (2001)
  algorithm, will be more appropriate for big data due to its high
  memory efficiency.
  
  The Greenwald-Khanna algorithm is implemented roughly as follows.  To
  form a distribution of values for each case, from which we sample to
  determine quantiles, we create a chain of values for the case as we
  grow the forest.  Every time a case lands in a terminal node, we
  insert all of its co-inhabitants to its chain of values.

  The best case scenario is when tree node size is 1 because each case
  gets only one insert into its chain for that tree.  The worst case
  scenario is when node size is so large that trees stump. This is
  because each case receives insertions for the entire in-bag
  population.  
  
  What the user needs to know is that Greenwald-Khanna can become slow
  in counter-intutive settings such as when node size is large.  The
  easy fix is to change the epsilon quantile approximation that is
  requested.  You will see a significant speed-up just by doubling
  \code{prob.epsilon}.  This is because the chains stay a lot smaller as
  epsilon increases, which is exactly what you want when node sizes are
  large. Both time and space requirements for the algorithm are affected
  by epsilon.

  The best results for Greenwald-Khanna come from setting the number of
  quantiles equal to 2 times the sample size and epsilon to 1 over 2
  times the sample size which is the default values used if left
  unspecified.  This will be slow, especially for big data, and less
  stringent choices should be used if computational speed is of concern.

}
\value{
  
  Returns quantiles for each of the requested probabilities.  Also
  returns the conditional density (and conditional cdf) for unique
  y-values in the training data (or test data if provided).  The
  conditional density can be used to calculate conditional moments, such
  as the mean and standard deviation.  For convenience, the mean is
  returned as the object \code{yhat}.

  For multivariate forests, the returned object will be a list of length
  equal to the number of target outcomes.

}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Greenwald M. and Khanna S. (2001).  Space-efficient online computation of
  quantile summaries. \emph{Proceedings of ACM SIGMOD}, 30(2):58--66.
  
  Meinshausen N. (2006) Quantile regression forests, \emph{Journal of
  Machine Learning Research}, 7:983--999.
}
\seealso{
  \command{\link{rfsrc}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## regression example
## ------------------------------------------------------------

## standard call
o <- quantileReg(mpg ~ ., mtcars)
qo <- o$quantileReg

## calculate the conditional mean, compare to OOB predicted value
## note that the conditional mean is returned as "yhat"
c.mean <- qo$density \%*\% qo$yunq
print(data.frame(c.mean = c.mean, yhat = qo$yhat, pred.oob = o$predicted.oob))

## calculate conditional standard deviation
c.std <- sqrt(qo$density \%*\% qo$yunq^2 - c.mean ^ 2)
quant <- qo$quantile
colnames(quant) <- paste("q", 100 * qo$prob, sep = "")
print(data.frame(quant, c.std))


## ------------------------------------------------------------
## train/test regression example
## ------------------------------------------------------------

## train (grow) call followed by test call
trn <- quantileReg(mpg ~ ., mtcars[1:20,])
test <- quantileReg(object = trn, newdata = mtcars[-(1:20),-1])

## calculate test set conditional mean and standard deviation
qo <- test$quantileReg
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
o <- quantileReg(cbind(carb, mpg, cyl, disp) ~., data = dta)

print(head(o$quantileReg$mpg$quant))
print(head(o$quantileReg$disp$quant))


## ------------------------------------------------------------
## quantile regression plot for Boston Housing data
## ------------------------------------------------------------

if (library("mlbench", logical.return = TRUE)) {

  ## apply quantile regression to Boston Housing data
  data(BostonHousing)
  o <- quantileReg(medv ~ ., BostonHousing, nodesize = 1)
  y <- o$yvar
  qo <- o$quantileReg

  ## pull desired quantiles - nice little wrapper for doing this
  get.quantile <- function(q, target.prob) {
    target.prob <- sort(unique(target.prob))
    q.dta <- do.call(cbind, lapply(target.prob, function(pr) {
      q$quant[, which.min(abs(pr - q$prob))]
    }))
    colnames(q.dta) <-  paste("q.", 100 * target.prob, sep = "")
    q.dta
   }

  ## extract 25,50,75 quantiles
  quant.dat <- get.quantile(qo, c(.25, .50, .75))
  
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

  ## compare 25-75 percentiles to values expected under normality
  c.mean <- qo$density \%*\% qo$yunq
  c.std <- sqrt(qo$density \%*\% qo$yunq^2 - c.mean ^ 2)
  q.25.est <- c.mean + qnorm(.25) * c.std
  q.75.est <- c.mean + qnorm(.75) * c.std
  print(head(data.frame(quant.dat[, -2],  q.25.est, q.75.est)))

  ## compare performance of quantile regression estimator to
  ## standard random forest estimator of averaged tree mean
  cat("quantile regression yhat error:", mean((o$yvar-qo$yhat)^2), "\n")
  cat("RF averaged tree mean error:", mean((o$yvar-o$predicted.oob)^2), "\n")


}

## ------------------------------------------------------------
## example of quantile regression for ordinal data
## ------------------------------------------------------------

 ## use the wine data for illustration
 data(wine, package = "_PROJECT_PACKAGE_NAME_")

 ## run quantile regression
 o <- quantileReg(quality ~ ., wine, ntree = 100)

 ## extract "probabilities" = density values
 qo <- o$quantileReg
 yunq <- qo$yunq
 yvar <- factor(cut(o$yvar, c(-1, yunq), labels = yunq)) 
 qo.dens <- qo$density
 colnames(qo.dens) <- yunq 
 qo.class <- _PROJECT_PACKAGE_NAME_:::bayes.rule(qo.dens)
 qo.confusion <- table(yvar, qo.class)
 qo.err <- 1 - diag(qo.confusion) / rowSums(qo.confusion)
 qo.confusion <- cbind(qo.confusion, qo.err)
 print(qo.confusion)
 cat("Normalized Brier:", 100 * _PROJECT_PACKAGE_NAME_:::brier(yvar, qo.dens), "\n")

}}
\keyword{quantile regression forests}

