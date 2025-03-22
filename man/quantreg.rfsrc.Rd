\name{quantreg.rfsrc}
\alias{quantreg.rfsrc}
\alias{quantreg}
\alias{extract.quantile}
\alias{get.quantile}
\alias{get.quantile.crps}
\alias{get.quantile.stat}
\title{Quantile Regression Forests}
\description{
  Grows a univariate or multivariate quantile regression forest and returns
  its conditional quantile and density values.  Can be used for
  both training and testing purposes.
}
\usage{\method{quantreg}{rfsrc}(formula, data, object, newdata,
  method = "local", splitrule = NULL, prob = NULL, prob.epsilon = NULL,
  oob = TRUE, fast = FALSE, maxn = 1e3, ...)
}
\arguments{

  \item{formula}{A symbolic description of the model to be fit.
    Must be specified unless \code{object} is given.}

  \item{data}{Data frame containing the y-outcome and x-variables in
    the model. Must be specified unless \code{object} is given.}

  \item{object}{(Optional) A previously grown quantile regression
    forest.}

  \item{newdata}{(Optional) Test data frame used for prediction.  Note
    that prediction on test data must always be done with the
    \code{quantreg} function and not the \code{predict} function.  See
    example below.}

  \item{method}{Method used to calculate quantiles.  Three methods are
    provided: (1) A variation of the method used in Meinshausen (2006)
    based on forest weight (\code{method = "forest"}); (2) The
    Greenwald-Khanna algorithm, suited for big data, and specified by
    any one of the following: "gk", "GK", "G-K", "g-k"; (3) The default
    method, \code{method = "local"}, which uses the local adjusted cdf
    approach of Zhang et al. (2019).  This does not rely on forest
    weights and is reasonably fast. See below for further discussion.}

  \item{splitrule}{The default action is local adaptive quantile regression splitting,
    but this can be over-ridden by the user.  Not applicable to
    multivariate forests.  See details below.}
  
  \item{prob}{Target quantile probabilities when training.  If left unspecified,
    uses percentiles (1 through 99) for \code{method = "forest"}, and 
    for Greenwald-Khanna selects equally spaced percentiles optimized
    for accuracy (see below).}
  
  \item{prob.epsilon}{Greenwald-Khanna allowable error for quantile
    probabilities when training.}
  
  \item{oob}{Return OOB (out-of-bag) quantiles?  If false, in-bag values
    are returned.}

  \item{fast}{Use fast random forests, \code{rfsrc.fast}, in place of
    \code{rfsrc}?  Improves speed but may be less accurate.}

  \item{maxn}{Maximum number of unique y training values used when
    calculating the conditional density.}

  \item{...}{Further arguments to be passed to the \code{rfsrc}
    function used for fitting the quantile regression forest.}  
}
\details{
  
  The most common method for calculating RF quantiles uses the method
  described in Meinshausen (2006) using forest weights.  The forest
  weights method employed here (specified using \code{method}="forest"),
  however differs in that quantiles are estimated using a
  weighted local cumulative distribution function estimator.  For this
  reason, results may differ from Meinshausen (2006).  Moreover, results
  may also differ as the default splitting rule uses local adaptive
  quantile regression splitting instead of CART regression mean squared
  splitting which was used by Meinshausen (2006).  Note that local
  adaptive quantile regression splitting is not available for
  multivariate forests which reverts to the default multivariate
  composite splitting rule. In multivariate regression, users however do
  have the option to over-ride this using
  Mahalanobis splitting by setting \code{splitrule="mahalanobis"}

  A second method for estimating quantiles uses the Greenwald-Khanna
  (2001) algorithm (invoked by \code{method}="gk", "GK", "G-K" or
  "g-k").  While this will not be as accurate as forest weights, the
  high memory efficiency of Greenwald-Khanna makes it feasible to
  implement in big data settings unlike forest weights.
  
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

  Finally, the default method, \code{method}="local", implements the
  locally adjusted cdf estimator of Zhang et al. (2019).  This does not
  use forest weights and is reasonably fast and can be used for large
  data.  However, this relies on the assumption of homogeneity of the
  error distribution, i.e. that errors are iid and therefore have equal
  variance.  While this is reasonably robust to departures of homogeneity,
  there are instances where this may perform poorly; see Zhang et
  al. (2019) for details.  If hetereogeneity is suspected we recommend
  \code{method}="forest".
  
}
\value{
  
  Returns the object \code{quantreg} containing quantiles for each of
  the requested probabilities (which can be conveniently extracted using
  \code{get.quantile}).  Also contains the conditional density (and
  conditional cdf) for each case in the training data (or test data if
  provided) evaluated at each of the unique grow y-values.  The
  conditional density can be used to calculate conditional moments, such
  as the mean and standard deviation.  Use \code{get.quantile.stat} as a
  way to conveniently obtain these quantities.

  For multivariate forests, returned values will be a list of length
  equal to the number of target outcomes.

}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Greenwald M. and Khanna S. (2001).  Space-efficient online computation of
  quantile summaries. \emph{Proceedings of ACM SIGMOD}, 30(2):58-66.
  
  Meinshausen N. (2006) Quantile regression forests, \emph{Journal of
  Machine Learning Research}, 7:983-999.

  Zhang H., Zimmerman J., Nettleton D. and Nordman D.J. (2019).  Random
  forest prediction intervals. \emph{The American Statistician}. 4:1-5.
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
o <- quantreg(mpg ~ ., mtcars)

## extract conditional quantiles
print(get.quantile(o))
print(get.quantile(o, c(.25, .50, .75)))

## extract conditional mean and standard deviation
print(get.quantile.stat(o))

## standardized continuous rank probabiliy score (crps) performance
plot(get.quantile.crps(o), type = "l")


## ------------------------------------------------------------
## train/test regression example
## ------------------------------------------------------------

## train (grow) call followed by test call
o <- quantreg(mpg ~ ., mtcars[1:20,])
o.tst <- quantreg(object = o, newdata = mtcars[-(1:20),])

## extract test set quantiles and conditional statistics
print(get.quantile(o.tst))
print(get.quantile.stat(o.tst))


## ------------------------------------------------------------
## quantile regression for Boston Housing using forest method
## ------------------------------------------------------------

if (library("mlbench", logical.return = TRUE)) {

  ## quantile regression with mse splitting
  data(BostonHousing)
  o <- quantreg(medv ~ ., BostonHousing, method = "forest", nodesize = 1)

  ## standardized continuous rank probabiliy score (crps) 
  plot(get.quantile.crps(o), type = "l")

  ## quantile regression plot
  plot.quantreg(o, .05, .95)
  plot.quantreg(o, .25, .75)

  ## (A) extract 25,50,75 quantiles
  quant.dat <- get.quantile(o, c(.25, .50, .75))

  ## (B) values expected under normality
  quant.stat <- get.quantile.stat(o)
  c.mean <- quant.stat$mean
  c.std <- quant.stat$std
  q.25.est <- c.mean + qnorm(.25) * c.std
  q.75.est <- c.mean + qnorm(.75) * c.std

  ## compare (A) and (B)
  print(head(data.frame(quant.dat[, -2],  q.25.est, q.75.est)))


}

## ------------------------------------------------------------
## multivariate mixed outcomes example
## quantiles are only returned for the continous outcomes
## ------------------------------------------------------------

dta <- mtcars
dta$cyl <- factor(dta$cyl)
dta$carb <- factor(dta$carb, ordered = TRUE)
o <- quantreg(cbind(carb, mpg, cyl, disp) ~., data = dta)

plot.quantreg(o, m.target = "mpg")
plot.quantreg(o, m.target = "disp")

## ------------------------------------------------------------
## multivariate regression example using Mahalanobis splitting
## ------------------------------------------------------------

dta <- mtcars
o <- quantreg(cbind(mpg, disp) ~., data = dta, splitrule = "mahal")

plot.quantreg(o, m.target = "mpg")
plot.quantreg(o, m.target = "disp")

## ------------------------------------------------------------
## example of quantile regression for ordinal data
## ------------------------------------------------------------

## use the wine data for illustration
data(wine, package = "randomForestSRC")

## run quantile regression
o <- quantreg(quality ~ ., wine, ntree = 100)

## extract "probabilities" = density values
qo.dens <- o$quantreg$density
yunq <- o$quantreg$yunq
colnames(qo.dens) <- yunq

## convert y to a factor 
yvar <- factor(cut(o$yvar, c(-1, yunq), labels = yunq)) 
 
## confusion matrix
qo.confusion <- get.confusion(yvar, qo.dens)
print(qo.confusion)

## normalized Brier score
cat("Brier:", 100 * get.brier.error(yvar, qo.dens), "\n")


## ------------------------------------------------------------
## example of large data using Greenwald-Khanna algorithm 
## ------------------------------------------------------------

## load the data and do quick and dirty imputation
data(housing, package = "randomForestSRC")
housing <- impute(SalePrice ~ ., housing,
         ntree = 50, nimpute = 1, splitrule = "random")

## Greenwald-Khanna algorithm 
## request a small number of quantiles 
o <- quantreg(SalePrice ~ ., housing, method = "gk",
        prob = (1:20) / 20, prob.epsilon = 1 / 20, ntree = 250)
plot.quantreg(o)
  
## ------------------------------------------------------------
## using mse splitting with local cdf method for large data
## ------------------------------------------------------------

## load the data and do quick and dirty imputation
data(housing, package = "randomForestSRC")
housing <- impute(SalePrice ~ ., housing,
         ntree = 50, nimpute = 1, splitrule = "random")

## use mse splitting and reduce number of trees
o <- quantreg(SalePrice ~ ., housing, splitrule = "mse", ntree = 250)
plot.quantreg(o)
  


}}
\keyword{quantile regression forests}

