\name{holdoutvimp}
\alias{holdoutvimp}
\alias{holdoutvimp.rfsrc}
\title{Hold out variable importance (VIMP)}
\description{
  Hold out VIMP is calculated from the error rate for trees 
  grown with and without a variable.  Applies to all families.
}
\usage{\method{holdoutvimp}{rfsrc}(formula, data,
  ntree = 1000 * ncol(data) / vtry,
  ntree.max = 2000,
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
  \item{ntree}{Number of trees used for growing the forest.}
  \item{ntree.max}{Maximum number of trees used when calculating prediction
    error for determing hold out VIMP.}
  \item{nsplit}{Non-negative integer value specifying number of 
    random split points used to split a node (deterministic splitting
    corresponds to the value zero and is much slower).}
  \item{ntime}{Integer value used for survival to
    constrain ensemble calculations to a grid of \code{ntime} time points.}
  \item{mtry}{Number of variables randomly selected as candidates for
    splitting a node.}
  \item{vtry}{Number of variables randomly selected to be held out when
    growning a tree.}
  \item{fast}{Use fast random forests, \code{rfsrcFast}, in place of
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

  Accuracy of hold out VIMP depends heavily on the size of the forest.
  If the number of trees is too small, then number of trees where v is
  held out will be small, and the resulting OOB error will have high
  variance.  Thus, \code{ntree} should be set fairly high - we recommend
  using 1000 times the number of features.  Increasing \code{vtry} is
  another way to increase number of hold out trees.  In particular,
  number of trees needed should decrease linearly with \code{vtry}.
  Keep in mind however that intrepretation of holdout VIMP is altered
  when \code{vtry} is different than 1.  This is likely to be more of a
  concern in low dimensional settings.

  Uses the new \code{get.tree} option in \code{predict} to extract
  specific trees from a forest and the hidden option \code{vtry} in
  \code{rfsrc}.  The latter creates a hidden array \code{holdout.array}
  of zeroes and ones indicating which variable to hold out in a tree
  where number of rows equals number of features and number of columns
  equals number of trees.  The array can also be passed as a hidden
  option but is not checked for coherence so users should be careful
  when doing so.  
}
\value{
  Hold out VIMP for each variable.  For multivariate forests, hold out
  VIMP is calculated for each of the target outcomes.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Lu M. and Ishwaran H. (2018). Expert Opinion: A prediction-based
  alternative to p-values in regression models. \emph{J. Thoracic and
  Cardiovascular Surgery}, 155(3), 1130--1136. 
}

\seealso{
  \command{\link{vimp}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## Boston housing example
## ------------------------------------------------------------

if (library("mlbench", logical.return = TRUE)) {

  data(BostonHousing)
  hv <- holdoutvimp(medv ~ ., BostonHousing)
  print(hv)

}

## ------------------------------------------------------------
## Multivariate regression analysis
## ------------------------------------------------------------

hv <- holdoutvimp(cbind(mpg, cyl) ~., mtcars)
print(hv)

## ------------------------------------------------------------
## White wine classification example
## ------------------------------------------------------------

data(wine, package = "_PROJECT_PACKAGE_NAME_")
wine$quality <- factor(wine$quality)
hv <- holdoutvimp(quality ~ ., wine, vtry = 5)
print(100 * hv)


## ------------------------------------------------------------
## pbc survival example
## ------------------------------------------------------------

data(pbc, package = "_PROJECT_PACKAGE_NAME_")
hv <- holdoutvimp(Surv(days, status) ~ ., pbc, splitrule = "random")
print(100 * hv)

## ------------------------------------------------------------
## WIHS competing risk example
## ------------------------------------------------------------

data(wihs, package = "_PROJECT_PACKAGE_NAME_")
hv <- holdoutvimp(Surv(time, status) ~ ., wihs, ntree = 1000)
print(100 * hv)

}}

\keyword{vimp}
