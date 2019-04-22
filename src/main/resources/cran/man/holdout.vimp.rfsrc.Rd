\name{holdout.vimp.rfsrc}
\alias{holdout.vimp.rfsrc}
\alias{holdout.vimp}
\title{Hold out variable importance (VIMP)}
\description{
  Hold out VIMP is calculated from the error rate for trees 
  grown with and without a variable.  Applies to all families.
}
\usage{\method{holdout.vimp}{rfsrc}(formula, data,
  ntree = function(p, vtry){1000 * p / vtry},
  ntree.max = 2000,
  ntree.allvars = NULL,
  nsplit = 10,
  ntime = 50,
  mtry = NULL,
  vtry = 1,
  fast = FALSE,
  verbose = TRUE, 
  ...)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit.}
  \item{data}{Data frame containing the y-outcome and x-variables.}
  \item{ntree}{Function specifying requested number of trees used for
    growing the forest.  Inputs are dimension and number of holdout
    variables. The requested number of trees can also be a number.}
  \item{ntree.max}{Maximum number of trees used when calculating prediction
    error for determing hold out VIMP.}
  \item{ntree.allvars}{Grow this many additional trees and use them for
    calculating the baseline error rate.  Ignored if \code{NULL}.}
  \item{nsplit}{Non-negative integer value specifying number of 
    random split points used to split a node (deterministic splitting
    corresponds to the value zero and is much slower).}
  \item{ntime}{Integer value used for survival to
    constrain ensemble calculations to a grid of \code{ntime} time points.}
  \item{mtry}{Number of variables randomly selected as candidates for
    splitting a node.}
  \item{vtry}{Number of variables randomly selected to be held out when
    growning a tree.}
  \item{fast}{Use fast random forests, \code{rfsrc.fast}, in place of
    \code{rfsrc}?  Improves speed but is less accurate.}
  \item{verbose}{Provide verbose output?}
  \item{...}{Further arguments to be passed to \code{\link{rfsrc}}.}
}
\details{

  Prior to growing a tree, a random set of \code{vtry} features are held
  out.  Tree growing proceeds as usual with the remaining features.
  Once the forest is grown, hold out VIMP for a given variable v is
  calculated as follows.  Gather all trees where v was held out and
  calculate OOB prediction error.  Next gather all trees were v was not
  held out and calculate OOB prediction error. Hold out VIMP for v is
  the difference between these two values.  Thus hold out VIMP measures
  the importance of a variable when that variable is truly removed from
  tree growing.

  If \code{ntree.allvars} is set to an integer value, then a total of
  this many trees are grown using all variables.  The above procedure is
  then implemented with the following change.  Determine the error rate
  for these additional trees.  Hold out VIMP for v is the difference
  between this value and the error rate for trees where v was held out.
  Unlike the above procedure, this makes sure that the baseline used for
  calculating holdout VIMP is the same for all v.  This feature is
  probably most useful in low-dimensional settings. 
  
  Note that accuracy of hold out VIMP depends heavily on the size of the
  forest.  If the number of trees is too small, then number of times a
  variable is held out will be small and OOB error may suffer from high
  variance.  Thus, \code{ntree} should be set fairly high - we recommend
  using 1000 times the number of features.  Increasing \code{vtry} is
  another way to increase number of hold out trees.  In particular,
  number of trees needed should decrease linearly with \code{vtry}.  For
  this reason the default \code{ntree} equals 1000 trees for each
  feature divided by \code{vtry}.  Keep in mind that intrepretation of
  holdout VIMP is altered when \code{vtry} is different than 1.

   
  
}
\value{
  Hold out VIMP for each variable.  For multivariate forests, hold out
  VIMP is calculated for each of the target outcomes.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{

  Ishwaran H. (2019).  Holdout variable importance for random forest
  models.
    
  Lu M. and Ishwaran H. (2018). Expert Opinion: A prediction-based
  alternative to p-values in regression models. \emph{J. Thoracic and
    Cardiovascular Surgery}, 155(3), 1130--1136.

  
}

\seealso{
  \command{\link{vimp.rfsrc}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## boston housing example
## ------------------------------------------------------------

if (library("mlbench", logical.return = TRUE)) {

  data(BostonHousing)
  hv <- holdout.vimp(medv ~ ., BostonHousing)
  print(hv)

}

## ------------------------------------------------------------
## iris example illustrating vtry
## ------------------------------------------------------------

print(100 * holdout.vimp(Species ~ ., iris))
print(100 * holdout.vimp(Species ~ ., iris, vtry=2))

## ------------------------------------------------------------
## example involving class imbalanced data
## illustrates the new RFQ classifier
## see the function "imbalanced" for more information about RFQ
## ------------------------------------------------------------

data(breast, package = "_PROJECT_PACKAGE_NAME_")
breast <- na.omit(breast)
f <- as.formula(status ~ .)
hv <- holdout.vimp(f, breast, rfq=TRUE, perf.type="g.mean")
print(10 * hv)


## ------------------------------------------------------------
## multivariate regression analysis example
## ------------------------------------------------------------

print(holdout.vimp(cbind(mpg, cyl) ~., mtcars))

## ------------------------------------------------------------
## white wine classification example
## ------------------------------------------------------------

data(wine, package = "_PROJECT_PACKAGE_NAME_")
wine$quality <- factor(wine$quality)
hv <- holdout.vimp(quality ~ ., wine, vtry = 5)
print(100 * hv)


## ------------------------------------------------------------
## pbc survival example
## ------------------------------------------------------------

data(pbc, package = "_PROJECT_PACKAGE_NAME_")
hv <- holdout.vimp(Surv(days, status) ~ ., pbc, splitrule = "random")
print(100 * hv)

## ------------------------------------------------------------
## WIHS competing risk example
## ------------------------------------------------------------

data(wihs, package = "_PROJECT_PACKAGE_NAME_")
hv <- holdout.vimp(Surv(time, status) ~ ., wihs, ntree = 1000)
print(100 * hv)

}}

\keyword{vimp}
