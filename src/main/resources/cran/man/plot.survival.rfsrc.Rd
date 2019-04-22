\name{plot.survival.rfsrc}
\alias{plot.survival.rfsrc}
\alias{plot.survival}
\title{Plot of Survival Estimates}
\description{
  Plot various survival estimates.
}
\usage{\method{plot.survival}{rfsrc}(x, plots.one.page = TRUE, 
  show.plots = TRUE, subset, collapse = FALSE,
  cens.model = c("km", "rfsrc"), ...)
}
\arguments{
  \item{x}{An object of class \code{(rfsrc, grow)} or
    \code{(rfsrc, predict)}.}
  \item{plots.one.page}{Should plots be placed on one page?}
  \item{show.plots}{Should plots be displayed?}
  \item{subset}{Vector indicating which individuals we want estimates for.
    All individuals are used if not specified.}
  \item{collapse}{Collapse the survival and cumulative hazard function
    across the individuals specified by \option{subset}?  Only applies
    when \option{subset} is specified.}
  \item{cens.model}{Method for estimating the censoring distribution
    used in the inverse probability of censoring weights (IPCW) for the
    Brier score:      
    \describe{
      \item{\code{km}:}{Uses the Kaplan-Meier estimator.}
      \item{\code{rfscr}:}{Uses a censoring random survival forest estimator.}
    }
  }
  \item{...}{Further arguments passed to or from other methods.}
}
\details{  
  If \option{subset} is not specified, generates the following 
  plots (going from top to bottom, left to right):
  
  \enumerate{
    
    \item Forest estimated survival function for each individual (thick
    red line is overall ensemble survival, thick green line is
    Nelson-Aalen estimator).
    
    \item Brier score (0=perfect, 1=poor, and 0.25=guessing) stratified
    by ensemble mortality.  Based on the IPCW method described in Gerds
    et al. (2006).  Stratification is into 4 groups corresponding to the
    0-25, 25-50, 50-75 and 75-100 percentile values of mortality.  Red
    line is the overall (non-stratified) Brier score.

    \item Continuous rank probability score (CRPS) equal to the
    integrated Brier score divided by time.
    
    \item Plot of mortality of each individual versus observed time.
    Points in blue correspond to events, black points are censored
    observations.

  }
  
  When \option{subset} is specified, then for each individual in
  \option{subset}, the following three plots are generated:
  
  \enumerate{
    \item Forest estimated survival function.
    
    \item Forest estimated cumulative hazard function (CHF) (displayed
    using black lines).  Blue lines are the CHF from the estimated hazard
    function.
  }

  Note that when the object \code{x} is of class \code{(rfsrc, predict)}
  not all plots will be produced.  In particular, Brier scores are not
  calculated.
  
  Only applies to survival families.  In particular, fails for competing
  risk analyses.  Use \code{plot.competing.risk} in such cases.
  
  Whenever possible, out-of-bag (OOB) values are used.
  
}
\value{
  Invisibly, the conditional and unconditional Brier scores, and the
  integrated Brier score (if they are available).
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Gerds T.A and Schumacher M. (2006).  Consistent estimation of the
  expected Brier score in general survival models with right-censored
  event times, \emph{Biometrical J.}, 6:1029-1040.
  
  Graf E., Schmoor C., Sauerbrei W. and Schumacher M. (1999).
  Assessment and comparison of prognostic classification
  schemes for survival data, \emph{Statist. in Medicine},
  18:2529-2545.
  
  Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7(2):25-31.

}
\seealso{
  \command{\link{plot.competing.risk.rfsrc}},
  \command{\link{predict.rfsrc}},
  \command{\link{rfsrc}}
}
\examples{
\donttest{
## veteran data
data(veteran, package = "_PROJECT_PACKAGE_NAME_") 
plot.survival(rfsrc(Surv(time, status)~ ., veteran), cens.model = "rfsrc")

## pbc data
data(pbc, package = "_PROJECT_PACKAGE_NAME_") 
pbc.obj <- rfsrc(Surv(days, status) ~ ., pbc)

## use subset to focus on specific individuals
plot.survival(pbc.obj, subset = 3)
plot.survival(pbc.obj, subset = c(3, 10))



}
}
\keyword{plot}
