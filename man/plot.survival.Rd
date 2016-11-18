\name{plot.survival}
\alias{plot.survival}
\alias{plot.survival.rfsrc}
\title{Plot of Survival Estimates}
\description{
  Plot various survival estimates.
}
\usage{\method{plot.survival}{rfsrc}(x, plots.one.page = TRUE, 
  show.plots = TRUE, subset, collapse = FALSE,
  haz.model = c("spline", "ggamma", "nonpar", "none"),
  k = 25, span = "cv", cens.model = c("km", "rfsrc"), ...)
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
  \item{haz.model}{Method for estimating the hazard.  See details below.
    Applies only when \option{subset} is specified.}
  \item{k}{The number of natural cubic spline knots used
    for estimating the hazard function.  Applies only when \option{subset}
    is specified.}
  \item{span}{The fraction of the observations in the span of Friedman's
    super-smoother used for estimating the hazard function.  Applies only when
    \option{subset} is specified.}
  \item{cens.model}{Method for estimating the censoring distribution
    used in the inverse probability of censoring weights (IPCW) for the
    Brier score:      
    \describe{
      \item{\code{km}:}{Uses the Kaplan-Meier estimator.}
      
      \item{\code{rfscr}:}{Uses random survival forests.}
    }
  }
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  
  If \option{subset} is not specified, generates the following three
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
    function.  See the next item.
    
    \item A smoothed hazard function derived from the forest estimated
    CHF (or survival function).  The default method,
    \option{haz.model="spline"}, models the log CHF using natural cubic
    splines as described in Royston and Parmar (2002).  The lasso is
    used for model selection, implemented using the \code{glmnet}
    package (this package must be installed for this option to work).
    If \option{haz.model="ggamma"}, a three-parameter generalized gamma
    distribution (using the parameterization described in Cox et al,
    2007) is fit to the smoothed forest survival function, where
    smoothing is imposed using Friedman's supersmoother (implemented by
    \code{supsmu}).  If \option{haz.model="nonpar"}, Friedman's
    supersmoother is applied to the forest estimated hazard function
    (obtained by taking the crude derivative of the smoothed forest
    CHF).  Finally, setting \option{haz.model="none"} suppresses
    hazard estimation and no hazard estimate is provided.
    
    At this time, please note that all hazard estimates are considered
    experimental and users should interpret the results with caution.}
  
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
  Cox C., Chu, H., Schneider, M. F. and Munoz, A. (2007).  Parametric
  survival analysis and taxonomy of hazard functions for the
  generalized gamma distribution.  Statistics in Medicine
  26:4252-4374.
  
  Gerds T.A and Schumacher M. (2006).  Consistent estimation of the
  expected Brier score in general survival models with right-censored
  event times, \emph{Biometrical J.}, 6:1029-1040.
  
  Graf E., Schmoor C., Sauerbrei W. and Schumacher M. (1999).
  Assessment and comparison of prognostic classification
  schemes for survival data, \emph{Statist. in Medicine},
  18:2529-2545.
  
  Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7(2):25-31.
  
  Royston P. and Parmar M.K.B. (2002).  Flexible parametric
  proportional-hazards and proportional-odds models for censored
  survival data, with application to prognostic modelling and
  estimation of treatment effects, \emph{Statist. in Medicine},
  21::2175-2197.
}
\seealso{
  \command{\link{plot.competing.risk}},
  \command{\link{predict.rfsrc}},
  \command{\link{rfsrc}}
}
\examples{
\dontrun{
## veteran data
data(veteran, package = "randomForestSRC") 
plot.survival(rfsrc(Surv(time, status)~ ., veteran), cens.model = "rfsrc")

## pbc data
data(pbc, package = "randomForestSRC") 
pbc.obj <- rfsrc(Surv(days, status) ~ ., pbc, nsplit = 10)

# default spline approach
plot.survival(pbc.obj, subset = 3)
plot.survival(pbc.obj, subset = 3, k = 100)

# three-parameter generalized gamma is approximately the same
# but notice that its CHF estimate (blue line) is not as accurate
plot.survival(pbc.obj, subset = 3, haz.model = "ggamma")

# nonparametric method is too wiggly or undersmooths
plot.survival(pbc.obj, subset = 3, haz.model = "nonpar", span = 0.1)
plot.survival(pbc.obj, subset = 3, haz.model = "nonpar", span = 0.8)

}
}
\keyword{plot}
