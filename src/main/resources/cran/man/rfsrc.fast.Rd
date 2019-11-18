\name{rfsrc.fast}
\alias{rfsrc.fast}
\title{Fast Random Forests}
\description{
  Fast approximate random forests using subsampling with forest options
  set to encourage computational speed.  Applies to all families.
}
\usage{rfsrc.fast(formula, data,
  ntree = 500,
  nsplit = 10,
  bootstrap = "by.root",
  ensemble = "oob",
  sampsize = function(x){min(x * .632, max(150, x ^ (3/4)))},
  samptype = "swor",
  samp = NULL,
  ntime = 50,
  forest = FALSE,
  terminal.qualts = FALSE,
  ...)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit. If
    missing, unsupervised splitting is implemented.}
  \item{data}{Data frame containing the y-outcome and x-variables.} 
  \item{ntree}{Number of trees.}
  \item{nsplit}{Non-negative integer value specifying number of 
    random split points used to split a node (deterministic splitting
    corresponds to the value zero and is much slower).}
  \item{bootstrap}{Bootstrap protocol used in growing a tree.}
  \item{ensemble}{Specifies the type of ensemble.  We request only
    out-of-sample which corresponds to "oob".}
  \item{sampsize}{Function specifying size of subsampled data. Can also be a number.}
  \item{samptype}{Type of bootstrap used.}
  \item{samp}{Bootstrap specification when \code{"by.user"} is used.}
  \item{ntime}{Integer value used for survival to
    constrain ensemble calculations to a grid of \code{ntime} time points.}
  \item{forest}{Should the forest object be returned?  Turn this on if
    you want prediction on test data but for big data this can be large.}
  \item{terminal.qualts}{Should terminal node membership information be
    returned?  Ensure this is off in the presence of big data as memory
    warnings can occur otherwise.  In either case, this parameter does not effect
    the ability to restore the model.  When this is \code{FALSE}, the
    training data must be sent down the tree when restoring the model.
    When this is \code{TRUE}, terminal node membership is assigned based
    on the contents of this output during model restoration, resulting
    in a faster CPU times.}
  \item{...}{Further arguments to be passed to \code{\link{rfsrc}}.}
}
\details{
  Calls \code{\link{rfsrc}} under various options (including subsampling) to
  encourage computational speeds.  This will provide a good
  approximation but will not be as good as default settings of
  \code{\link{rfsrc}}.
}
\value{
  An object of class \code{(rfsrc, grow)}.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\seealso{
  \command{\link{rfsrc}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## Iowa housing regression example
## ------------------------------------------------------------

## load the Iowa housing data
data(housing, package = "_PROJECT_PACKAGE_NAME_")

## do quick and *dirty* imputation
housing <- impute(SalePrice ~ ., housing,
         ntree = 50, nimpute = 1, splitrule = "random")

## grow a fast forest
o1 <- rfsrc.fast(SalePrice ~ ., housing)
o2 <- rfsrc.fast(SalePrice ~ ., housing, nodesize = 1)
print(o1)
print(o2)

## grow a fast bivariate forest
o3 <- rfsrc.fast(cbind(SalePrice,Overall.Qual) ~ ., housing)
print(o3)

## ------------------------------------------------------------
## White wine classification example
## ------------------------------------------------------------

data(wine, package = "_PROJECT_PACKAGE_NAME_")
wine$quality <- factor(wine$quality)
o <- rfsrc.fast(quality ~ ., wine)
print(o)


## ------------------------------------------------------------
## pbc survival example
## ------------------------------------------------------------

data(pbc, package = "_PROJECT_PACKAGE_NAME_")
o <- rfsrc.fast(Surv(days, status) ~ ., pbc)
print(o)

## ------------------------------------------------------------
## WIHS competing risk example
## ------------------------------------------------------------

data(wihs, package = "_PROJECT_PACKAGE_NAME_")
o <- rfsrc.fast(Surv(time, status) ~ ., wihs)
print(o)

}}
\keyword{forest}
\keyword{fast}