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
   m.target = NULL, crps = TRUE, subset = NULL, ...)
}
\arguments{
  \item{x}{A quantile regression object obtained from calling \command{quantreg}.}
  \item{prbL}{Lower quantile (preferably < .5).}
  \item{prbU}{Upper quantile (preferably > .5).}
  \item{m.target}{Character value for multivariate families
    specifying the target outcome to be used.  If left unspecified, the
    algorithm will choose a default target.}
  \item{crps}{Calculate crps and inset it?}
  \item{subset}{Restricts plotted values to a subset of the data.
    Default is to use the entire data.}
  \item{...}{Further arguments passed to or from other methods.}
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\seealso{
  \command{\link{quantreg.rfsrc}}
}
\keyword{plot}
