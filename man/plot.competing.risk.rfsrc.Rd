\name{plot.competing.risk.rfsrc}
\alias{plot.competing.risk.rfsrc}
\alias{plot.competing.risk}
\title{Plots for Competing Risks}
\description{
  Plot useful summary curves from a random survival forest competing risk
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

  Given a random survival forest object from a competing risk analysis
  (Ishwaran et al. 2014), plots from top to bottom, left to right: (1)
  cause-specific cumulative hazard function (CSCHF) for each event, (2)
  cumulative incidence function (CIF) for each event, and (3) continuous
  probability curves (CPC) for each event (Pepe and Mori, 1993).
  
  Does not apply to right-censored data.  Whenever possible, out-of-bag
  (OOB) values are displayed.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  
  Ishwaran H., Gerds T.A., Kogalur U.B., Moore R.D., Gange S.J. and Lau
  B.M. (2014). Random survival forests for competing risks.
  \emph{Biostatistics}, 15(4):757-773.

  Pepe, M.S. and Mori, M., (1993). Kaplan-Meier, marginal or conditional
  probability curves in summarizing competing risks failure time
  data? \emph{Statistics in Medicine}, 12(8):737-751.
 
  
}
\seealso{
  \command{\link{follic}},
  \command{\link{hd}},
  \command{\link{rfsrc}},
  \command{\link{wihs}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## follicular cell lymphoma
## ------------------------------------------------------------

  data(follic, package = "randomForestSRC")
  follic.obj <- rfsrc(Surv(time, status) ~ ., follic, nsplit = 3, ntree = 100)
  print(follic.obj)
  plot.competing.risk(follic.obj)

## ------------------------------------------------------------
## Hodgkin's Disease
## ------------------------------------------------------------

  data(hd, package = "randomForestSRC")
  hd.obj <- rfsrc(Surv(time, status) ~ ., hd, nsplit = 3, ntree = 100)
  print(hd.obj)
  plot.competing.risk(hd.obj)

## ------------------------------------------------------------
## competing risk analysis of pbc data from the survival package
## events are transplant (1) and death (2)
## ------------------------------------------------------------

if (library("survival", logical.return = TRUE)) {
   data(pbc, package = "survival")
   pbc$id <- NULL
   plot.competing.risk(rfsrc(Surv(time, status) ~ ., pbc))
}
}
}
\keyword{plot}
