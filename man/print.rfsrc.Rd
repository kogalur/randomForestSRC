\name{print.rfsrc}
\alias{print.rfsrc}
\title{Print Summary Output of a RF-SRC Analysis}
\description{
  Print summary output from a RF-SRC analysis. This is the default
  print method for the package.
}
\usage{\method{print}{rfsrc}(x, outcome.target = NULL, ...)}
\arguments{
  \item{x}{An object of class \code{(rfsrc, grow)}, \code{(rfsrc, synthetic)},
    or \code{(rfsrc, predict)}.}
  \item{outcome.target}{Character value for multivariate families
    specifying the target outcome to be used. The default is to use the
    first coordinate.}
  \item{...}{Further arguments passed to or from other methods.}
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7/2:25-31.
}
\seealso{
  \command{\link{rfsrc}}, \command{\link{rfsrcSyn}},
  \command{\link{predict.rfsrc}}
}
\examples{
  iris.obj <- rfsrc(Species ~., data = iris, ntree=100)
  print(iris.obj)
}
\keyword{print}
