\name{plot.competing.risk}
\alias{plot.competing.risk}
\alias{plot.competing.risk.rfsrc}
\title{Plots for Competing Risks}
\description{
  Plot the ensemble cumulative incidence function (CIF) and
  cause-specific cumulative hazard function (CSCHF) from a competing risk
  analysis.
}
\usage{\method{plot.competing.risk}{rfsrc}(x, plots.one.page = FALSE, ...)}
\arguments{
  \item{x}{An object of class \code{(rfsrc, grow)} or
    \code{(rfsrc, predict)}.}
  \item{plots.one.page}{Should plots be placed on one page?}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Ensemble ensemble CSCHF and CIF functions for each event type.
  Does not apply to right-censored data.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran H., Gerds T.A., Kogalur U.B., Moore R.D., Gange S.J. and Lau
  B.M. (2014). Random survival forests for competing risks.
  \emph{Biostatistics}, 15(4):757-773.
}
\seealso{
  \command{\link{follic}},
  \command{\link{hd}},
  \command{\link{rfsrc}},
  \command{\link{wihs}}
}
\examples{
\dontrun{
## ------------------------------------------------------------
## follicular cell lymphoma
## ------------------------------------------------------------

  data(follic, package = "randomForestSRC")
  follic.obj <- rfsrc(Surv(time, status) ~ ., follic, nsplit = 3, ntree = 100)
  plot.competing.risk(follic.obj)

## ------------------------------------------------------------
## competing risk analysis of pbc data from the survival package
## events are transplant (1) and death (2)
## ------------------------------------------------------------

if (library("survival", logical.return = TRUE)) {
   data(pbc, package = "survival")
   pbc$id <- NULL
   plot.competing.risk(rfsrc(Surv(time, status) ~ ., pbc, nsplit = 10))
}
}
}
\keyword{plot}
