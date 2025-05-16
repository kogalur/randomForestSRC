\name{plot.quantreg.rfsrc}
\alias{plot.quantreg.rfsrc}
\alias{plot.quantreg}
\title{Plot Quantiles from Quantile Regression Forests}
\description{
  Plots quantiles obtained from a quantile regression forest.
  Additionally insets the continuous rank probability score (crps), a
  useful diagnostic of accuracy.
}
\usage{\method{plot.quantreg}{rfsrc}(x, prbL = .25, prbU = .75,
   m.target = NULL, crps = TRUE, subset = NULL, xlab = NULL, ylab = NULL, ...)
}
\arguments{

  \item{x}{A quantile regression object returned by a call to \command{quantreg}.}

  \item{prbL}{Lower quantile level, typically less than \code{0.5}.}

  \item{prbU}{Upper quantile level, typically greater than \code{0.5}.}

  \item{m.target}{Character string specifying the target outcome for multivariate families. 
  If not provided, a default target is selected automatically.}

  \item{crps}{Logical. If \code{TRUE}, calculates the continuous ranked probability score (CRPS) and adds it to the plot.}

  \item{subset}{Optional vector specifying a subset of the data to be plotted. 
  Defaults to plotting all data points.}

  \item{xlab}{Label for the x-axis.}

  \item{ylab}{Label for the y-axis.}

  \item{...}{Additional arguments passed to or from other methods.}
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\seealso{
  \command{\link{quantreg.rfsrc}}
}
\keyword{plot}
