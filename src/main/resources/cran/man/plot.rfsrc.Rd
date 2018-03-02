\name{plot.rfsrc}
\alias{plot.rfsrc}
\title{Plot Error Rate and Variable Importance from a RF-SRC analysis}
\description{
  Plot out-of-bag (OOB) error rates and variable importance (VIMP)
  from a RF-SRC analysis. This is the default plot method for the package.
}
\usage{\method{plot}{rfsrc}(x, m.target = NULL,
  plots.one.page = TRUE, sorted = TRUE, verbose = TRUE,  ...)}
\arguments{
  \item{x}{An object of class \code{(rfsrc, grow)}, \code{(rfsrc, synthetic)},
    or \code{(rfsrc, predict)}.}
  \item{m.target}{Character value for multivariate families
    specifying the target outcome to be used.  If left unspecified, the
    algorithm will choose a default target.}
  \item{plots.one.page}{Should plots be placed on one page?}
  \item{sorted}{Should variables be sorted by importance values?}
  \item{verbose}{Should VIMP be printed?}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{  
  Plot cumulative OOB error rates as a function of number of trees and
  variable importance (VIMP) if available.  Note that the default
  settings are now such that the error rate is no longer calculated on
  every tree and VIMP is only calculated if requested.  To get OOB error
  rates for ever tree, use the option \code{tree.err = TRUE} when
  growing the forest or restore the model using the option
  \code{tree.err = TRUE}.  Likewise, to view VIMP, use the option
  \code{importance} when growing the forest or restore the forest using
  the option \code{importance}.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.
  
  Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7(2):25-31.
}
\seealso{
  \command{\link{predict.rfsrc}},
  \command{\link{rfsrc}}
}
\examples{
\dontrun{
## ------------------------------------------------------------
## classification example
## ------------------------------------------------------------

iris.obj <- rfsrc(Species ~ ., data = iris,
     tree.err = TRUE, importance = TRUE)
plot(iris.obj)

## ------------------------------------------------------------
## competing risk example
## ------------------------------------------------------------

## use the pbc data from the survival package
## events are transplant (1) and death (2)
if (library("survival", logical.return = TRUE)) {
  data(pbc, package = "survival")
  pbc$id <- NULL
  plot(rfsrc(Surv(time, status) ~ ., pbc, nsplit = 10, tree.err = TRUE))
}

## ------------------------------------------------------------
## multivariate mixed forests
## ------------------------------------------------------------

mtcars.new <- mtcars
mtcars.new$cyl <- factor(mtcars.new$cyl)
mtcars.new$carb <- factor(mtcars.new$carb, ordered = TRUE)
mv.obj <- rfsrc(cbind(carb, mpg, cyl) ~., data = mtcars.new, tree.err = TRUE)
plot(mv.obj, m.target = "carb")
plot(mv.obj, m.target = "mpg")
plot(mv.obj, m.target = "cyl")

}
}
\keyword{plot}
