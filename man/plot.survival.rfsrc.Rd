\name{plot.survival.rfsrc}
\alias{plot.survival.rfsrc}
\alias{plot.survival}
\alias{get.brier.survival}
\title{Plot of Survival Estimates}
\description{
  Plot various survival estimates.
}
\usage{\method{plot.survival}{rfsrc}(x, show.plots = TRUE, subset,
  collapse = FALSE, cens.model = c("km", "rfsrc"), ...)
}
\arguments{
  \item{x}{An object of class \code{(rfsrc, grow)} or
    \code{(rfsrc, predict)}.}

  \item{show.plots}{Should plots be displayed?}

  \item{subset}{Vector indicating which cases from \code{x} we want
    estimates for.  All cases used if not specified.}

  \item{collapse}{Collapse the survival function?}

  \item{cens.model}{Using the training data, specifies method for
    estimating the censoring distribution used in the inverse
    probability of censoring weights (IPCW) for calculating the Brier
    score:

    \describe{
      \item{\code{km}:}{Uses the Kaplan-Meier estimator.}
      \item{\code{rfscr}:}{Uses a censoring random survival forest estimator.}
    }
  }
  
 \item{...}{Further arguments passed to or from other methods.}

}
\details{  
  Produces the following plots (going from top to bottom, left to right):
  
  \enumerate{
    
    \item Forest estimated survival function for each individual (thick
    red line is overall ensemble survival, thick green line is
    Nelson-Aalen estimator).
    
    \item Brier score (0=perfect, 1=poor, and 0.25=guessing) stratified
    by ensemble mortality.  Based on the IPCW method described in Gerds
    et al. (2006).  Stratification is into 4 groups corresponding to the
    0-25, 25-50, 50-75 and 75-100 percentile values of mortality.  Red
    line is overall (non-stratified) Brier score.

    \item Continuous rank probability score (CRPS) equal to the
    integrated Brier score divided by time.
    
    \item Plot of mortality of each individual versus observed time.
    Points in blue correspond to events, black points are censored
    observations.  Not given for prediction settings lacking
    survival response information.
    
  }
  
  Whenever possible, out-of-bag (OOB) values are used.

  Only applies to survival families.  In particular, fails for competing
  risk analyses.  Use \command{plot.competing.risk} in such cases.

  Mortality (Ishwaran et al., 2008) represents estimated risk for an
  individual calibrated to the scale of number of events (as a specific
  example, if \emph{i} has a mortality value of 100, then if all
  individuals had the same x-values as \emph{i}, we would expect an
  average of 100 events).
  
  The utility function \code{get.brier.survival} can be used to extract
  the Brier score among other useful quantities.
  
}
\value{
  Invisibly, the conditional and unconditional Brier scores, and the
  integrated Brier score.
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

  Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
  (2008).  Random survival forests, \emph{Ann. App. Statist.}, 2:841-860.

}
\seealso{
  \command{\link{plot.competing.risk.rfsrc}},
  \command{\link{predict.rfsrc}},
  \command{\link{rfsrc}}
}
\examples{
\donttest{
## veteran data
data(veteran, package = "randomForestSRC") 
plot.survival(rfsrc(Surv(time, status)~ ., veteran), cens.model = "rfsrc")

## pbc data
data(pbc, package = "randomForestSRC") 
pbc.obj <- rfsrc(Surv(days, status) ~ ., pbc)

## use subset to focus on specific individuals
plot.survival(pbc.obj, subset = 3)
plot.survival(pbc.obj, subset = c(3, 10))
plot.survival(pbc.obj, subset = c(3, 10), collapse = TRUE)

## get.brier.survival function does many nice things!
plot(get.brier.survival(pbc.obj, cens.model="km")$brier.score,type="s", col=2)
lines(get.brier.survival(pbc.obj, cens.model="rfsrc")$brier.score, type="s", col=4)
legend("bottomright", legend=c("cens.model = km", "cens.model = rfsrc"), fill=c(2,4))

}
}
\keyword{plot}
