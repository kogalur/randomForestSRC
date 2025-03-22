\name{vimp.rfsrc}
\alias{vimp.rfsrc}
\alias{vimp}
\title{VIMP for Single or Grouped Variables}
\description{
  Calculate variable importance (VIMP) for a single variable or group of
  variables for training or test data.
}
\usage{\method{vimp}{rfsrc}(object, xvar.names, m.target = NULL, 
  importance = c("anti", "permute", "random"), block.size = 10,
  joint = FALSE, seed = NULL, do.trace = FALSE, ...)
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

  \item{seed}{Negative integer specifying seed for the random number
    generator.} 

  \item{do.trace}{Number of seconds between updates to the user on
    approximate time to completion.}

  \item{...}{Further arguments passed to or from other methods.}
  
}
\details{
  Using a previously trained forest, calculate the VIMP for variables
  \code{xvar.names}.  By default, VIMP is calculated for the original
  data, but the user can specify a new test data for the VIMP
  calculation using \code{newdata}.  See \command{rfsrc} for more
  details about how VIMP is calculated.
  
  \option{joint=TRUE} returns joint VIMP, defined as importance for
  a group of variables when the group is perturbed simultaneously.

  \code{csv=TRUE} return case specific VIMP.  Applies to
  all families except survival families.  See example below.  
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

## anti vimp (default)
print(vimp(iris.obj)$importance)

## anti vimp using brier prediction error
print(vimp(iris.obj, perf.type = "brier")$importance)

## permutation vimp
print(vimp(iris.obj, importance = "permute")$importance)

## random daughter vimp
print(vimp(iris.obj, importance = "random")$importance)

## joint anti vimp 
print(vimp(iris.obj, joint = TRUE)$importance)

## paired anti vimp
print(vimp(iris.obj, c("Petal.Length", "Petal.Width"), joint = TRUE)$importance)
print(vimp(iris.obj, c("Sepal.Length", "Petal.Width"), joint = TRUE)$importance)

## ------------------------------------------------------------
## survival example
## anti versus permute VIMP with different block sizes
## ------------------------------------------------------------

data(pbc, package = "randomForestSRC")
pbc.obj <- rfsrc(Surv(days, status) ~ ., pbc)

print(vimp(pbc.obj)$importance)
print(vimp(pbc.obj, block.size=1)$importance)
print(vimp(pbc.obj, importance="permute")$importance)
print(vimp(pbc.obj, importance="permute", block.size=1)$importance)

## ------------------------------------------------------------
## imbalanced classification example
## see the imbalanced function for more details
## ------------------------------------------------------------

data(breast, package = "randomForestSRC")
breast <- na.omit(breast)
f <- as.formula(status ~ .)
o <- rfsrc(f, breast, ntree = 2000)

## permutation vimp
print(100 * vimp(o, importance = "permute")$importance)

## anti vimp using gmean performance
print(100 * vimp(o, perf.type = "gmean")$importance[, 1])

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

## training data vimp
print(airq.obj$importance)
print(vimp(airq.obj)$importance)

## test data vimp
print(vimp(airq.obj, newdata = airquality[-train, ])$importance)

## ------------------------------------------------------------
## case-specific vimp
## returns VIMP for each case
## ------------------------------------------------------------

o <- rfsrc(mpg~., mtcars)
v <- vimp(o, csv = TRUE)
csvimp <- get.mv.csvimp(v, standardize=TRUE)
print(csvimp)

## ------------------------------------------------------------
## case-specific joint vimp
## returns joint VIMP for each case
## ------------------------------------------------------------

o <- rfsrc(mpg~., mtcars)
v <- vimp(o, joint = TRUE, csv = TRUE)
csvimp <- get.mv.csvimp(v, standardize=TRUE)
print(csvimp)

## ------------------------------------------------------------
## case-specific joint vimp for multivariate regression
## returns joint VIMP for each case, for each outcome
## ------------------------------------------------------------

o <- rfsrc(Multivar(mpg, cyl) ~., data = mtcars)
v <- vimp(o, joint = TRUE, csv = TRUE)
csvimp <- get.mv.csvimp(v, standardize=TRUE)
print(csvimp)

}
}
\keyword{variable selection}
\keyword{predict}
