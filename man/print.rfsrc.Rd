\name{print.rfsrc}
\alias{print.rfsrc}
\title{Print Summary Output of a RF-SRC Analysis}
\description{
  Print summary output from a RF-SRC analysis. This is the default
  print method for the package.
}
\usage{\method{print}{rfsrc}(x, outcome.target = NULL, ...,
                          quantreg.tau = NULL)}
\arguments{
  \item{x}{An object of class \code{(rfsrc, grow)}, 
    or \code{(rfsrc, predict)}.}
  \item{outcome.target}{Character value for multivariate families
    specifying the target outcome to be used. The default is to use the
    first coordinate from the continuous outcomes (otherwise if none,
    the first coordinate from the categorical outcomes).}
  \item{...}{Further arguments passed to or from other methods.}
  \item{quantreg.tau}{Optional numeric vector of finite levels strictly
    between zero and one for the printed pinball losses of a quantile
    regression object. Supply this argument by name. When \code{NULL},
    uses the session option \code{quantreg.tau}, then the former
    \code{rfsrc.pinball.taus} option for compatibility, then
    \code{c(0.1, 0.5, 0.9)}. Explicit levels apply to this print call
    only and do not change the fitted object or session options.
    See \code{\link{quantreg}} for details.}
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7/2:25-31.
}
\examples{
  iris.obj <- rfsrc(Species ~., data = iris, ntree=10)
  print(iris.obj)
}
\keyword{print}
