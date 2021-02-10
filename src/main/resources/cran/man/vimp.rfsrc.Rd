\name{vimp.rfsrc}
\alias{vimp.rfsrc}
\alias{vimp}
\title{VIMP for Single or Grouped Variables}
\description{
  Calculate variable importance (VIMP) for a single variable or group of
  variables for training or test data.
}
\usage{\method{vimp}{rfsrc}(object, xvar.names, m.target = NULL, 
  importance = c("permute", "random", "anti"), block.size = 10,
  joint = FALSE, subset, seed = NULL, do.trace = FALSE, ...)
}
\arguments{
  \item{object}{An object of class \code{(rfsrc, grow)} or
    \code{(rfsrc, forest)}. Requires \option{forest=TRUE} in the
    original \command{rfsrc} call.}
  \item{xvar.names}{Names of the x-variables to be used.  If not
    specified all variables are used.}
  \item{m.target}{Character value for multivariate families
    specifying the target outcome to be used.  If left unspecified, the
    algorithm will choose a default target.}
  \item{importance}{Type of VIMP.}
  \item{block.size}{Specifies number of trees in a block when calculating
    VIMP.}
  \item{joint}{Individual or joint VIMP?}
  \item{subset}{Vector indicating which rows of the grow data to
    restrict VIMP calculations to; i.e. this option yields VIMP which is
    restricted to a specific subset of the data.  Note that the vector
    should correspond to the rows of \code{object$xvar} and not the
    original data passed in the grow call.  All rows used if not
    specified.}
  \item{seed}{Negative integer specifying seed for the random number
    generator.} 
  \item{do.trace}{Number of seconds between updates to the user on
    approximate time to completion.}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Using a previously grown forest, calculate the VIMP for variables
  \code{xvar.names}.  By default, VIMP is calculated for the original
  data, but the user can specify a new test data for the VIMP
  calculation using \code{newdata}.  See \command{rfsrc} for more
  details about how VIMP is calculated.
  
  Joint VIMP is requested using \option{joint} and equals importance for
  a group of variables when the group is perturbed simultaneously.
}
\value{
  An object of class \code{(rfsrc, predict)} containing importance
  values.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran H. (2007).  Variable importance in binary regression
  trees and forests,  \emph{Electronic J. Statist.}, 1:519-537.
}
\seealso{
  \command{\link{holdout.vimp.rfsrc}},
  \command{\link{rfsrc}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## classification example
## showcase different vimp
## ------------------------------------------------------------

iris.obj <- rfsrc(Species ~ ., data = iris)

# Permutation vimp (default)
print(vimp(iris.obj)$importance)

# VIMP using brier prediction error
print(vimp(iris.obj, perf.type = "brier")$importance)

# Random daughter vimp
print(vimp(iris.obj, importance = "random")$importance)

# Joint permutation vimp 
print(vimp(iris.obj, joint = TRUE)$importance)

# Paired vimp
print(vimp(iris.obj, c("Petal.Length", "Petal.Width"), joint = TRUE)$importance)
print(vimp(iris.obj, c("Sepal.Length", "Petal.Width"), joint = TRUE)$importance)



## ------------------------------------------------------------
## regression example
## ------------------------------------------------------------

airq.obj <- rfsrc(Ozone ~ ., airquality)
print(vimp(airq.obj))


## ------------------------------------------------------------
## regression example where vimp is calculated on test data
## ------------------------------------------------------------

set.seed(100080)
train <- sample(1:nrow(airquality), size = 80)
airq.obj <- rfsrc(Ozone~., airquality[train, ])

#training data vimp
print(airq.obj$importance)
print(vimp(airq.obj)$importance)

#test data vimp
print(vimp(airq.obj, newdata = airquality[-train, ])$importance)

## ------------------------------------------------------------
## survival example
## study how vimp depends on tree imputation
## makes use of the subset option
## ------------------------------------------------------------

data(pbc, package = "_PROJECT_PACKAGE_NAME_")

# determine which records have missing values
which.na <- apply(pbc, 1, function(x){any(is.na(x))})

# impute the data using na.action = "na.impute"
pbc.obj <- rfsrc(Surv(days,status) ~ ., pbc, nsplit = 3,
        na.action = "na.impute", nimpute = 1)

# compare vimp based on records with no missing values
# to those that have missing values
# note the option na.action="na.impute" in the vimp() call
vimp.not.na <- vimp(pbc.obj, subset = !which.na, na.action = "na.impute")$importance
vimp.na <- vimp(pbc.obj, subset = which.na, na.action = "na.impute")$importance
print(data.frame(vimp.not.na, vimp.na))
}
}
\keyword{variable selection}
\keyword{predict}
