\name{plot.survival.rfsrc}
\alias{plot.survival.rfsrc}
\alias{plot.survival}
\alias{get.brier.survival}
\alias{get.auct.survival}
\alias{plotBrierAUC}
\title{Survival Prediction Performance and Diagnostic Plots}
\description{
  Plot survival estimates and calculate inverse-probability-of-censoring
  weighted Brier scores and cumulative/dynamic time-dependent AUC.
}
\usage{
\method{plot.survival}{rfsrc}(x, show.plots = TRUE, subset,
  collapse = FALSE, cens.model = c("km", "rfsrc"), ...)

get.brier.survival(o, subset = NULL,
  cens.model = c("km", "rfsrc"), papply = lapply,
  times = NULL, conf.int = FALSE, keep.matrix = TRUE)

get.auct.survival(o, subset = NULL,
  cens.model = c("km", "rfsrc"), papply = lapply,
  times = NULL, conf.int = FALSE)

plotBrierAUC(x, subset = NULL,
  cens.model = c("km", "rfsrc"), papply = lapply,
  times = NULL, conf.int = TRUE,
  plots = c("brier", "auct"), show.plots = TRUE, ...)
}
\arguments{
  \item{x, o}{An object of class \code{(rfsrc, grow)} or
    \code{(rfsrc, predict)}. The two numerical helpers and
    \code{plotBrierAUC} also accept an object of class
    \code{(rfsrc, forest)}.}

  \item{show.plots}{Should plots be displayed?}

  \item{subset}{Vector indicating which predicted cases are to be used.
    All cases are used if not specified.}

  \item{collapse}{Collapse the individual survival functions into their
    ensemble mean?}

  \item{cens.model}{Method used to estimate the censoring distribution
    for inverse probability of censoring weighting. In all cases the
    censoring model is estimated from the full grow data, not from the
    evaluation outcomes:

    \describe{
      \item{\code{km}:}{Kaplan--Meier estimator.}
      \item{\code{rfsrc}:}{Random survival forest estimator of the
        conditional censoring distribution. When censoring is present,
        this option requires the grow and evaluation covariates to be stored
        in the object. Because
        \code{rfsrc.fast} uses \code{forest=FALSE} by default, use
        \code{cens.model="km"} for its default reduced object or refit with
        \code{forest=TRUE}.}
    }
  }

  \item{papply}{Function used in place of \code{lapply} for repeated
    calculations.}

  \item{times}{Optional numeric vector of prediction horizons. The
    forest survival curves are evaluated as right-continuous step
    functions. The forest time grid is used by default.}

  \item{conf.int}{Controls pointwise confidence intervals. Use
    \code{FALSE} for no intervals, \code{TRUE} for 95 percent intervals,
    or a numeric value strictly between zero and one for a different
    confidence level.}

  \item{keep.matrix}{Should the subject-by-time matrix of IPCW Brier
    contributions be retained? This matrix is used by
    \code{plot.survival} for mortality-stratified curves.}

  \item{plots}{One or both of \code{"brier"} and \code{"auct"}.
    \code{"auc"} is also accepted as an alias for \code{"auct"}.}

  \item{...}{Further graphical arguments. For \code{plotBrierAUC},
    \code{col}, \code{lwd}, \code{lty}, and \code{type} control
    the estimated curve; \code{band.col} and \code{band.border} control
    the confidence band.}
}
\details{
  \code{plot.survival} produces the following plots, going from top to
  bottom and left to right:

  \enumerate{
    \item Forest estimated survival function for each individual. The
    thick red line is the ensemble mean survival and the thick green
    line is the marginal Nelson--Aalen survival estimate.

    \item Brier score stratified by ensemble mortality. Stratification
    is into four groups corresponding to the 0--25, 25--50, 50--75 and
    75--100 percentile ranges of mortality. The red line is the overall
    Brier score.

    \item Cumulative integrated Brier score divided by elapsed time,
    labeled CRPS in the plot.

    \item Mortality (Ishwaran et al., 2008) versus observed time. Blue
    points are events and black points are censored observations.
    Mortality is estimated risk calibrated to the scale of the number of
    events. For example, a mortality value of 100 means that if all
    individuals had the same covariate values, an average of 100 events
    would be expected.

  }

  Whenever possible, out-of-bag predictions are used. For a prediction
  object, the predicted test survival curves and test outcomes are used
  for scoring, while the censoring distribution is still estimated from
  the grow data. If a prediction object contains no outcomes,
  \code{plot.survival} displays only its predicted survival curves.
  Brier score, AUC, CRPS, and mortality-versus-time cannot be calculated
  without evaluation outcomes. Direct calls to \code{get.brier.survival},
  \code{get.auct.survival}, or \code{plotBrierAUC} therefore require
  outcomes for the cases being evaluated.

  For subject \eqn{i} at horizon \eqn{t}, the IPCW Brier contribution is
  \deqn{
    L_i(t) =
    \frac{I(T_i \le t, \Delta_i > 0)}{\widehat G(T_i-)}
       \widehat S_i(t)^2
    +
    \frac{I(T_i > t)}{\widehat G(t)}
       \{1-\widehat S_i(t)\}^2,
  }
  with subject-specific values of \eqn{\widehat G} when
  \code{cens.model="rfsrc"}. The Brier score is the sample mean of these
  contributions.

  The time-dependent AUC uses a cumulative/dynamic definition. Cases
  at \eqn{t} are observed events satisfying \eqn{T_i \le t}; controls
  satisfy \eqn{T_i > t}; observations censored at or before \eqn{t} are
  not included in the case-control comparison. The time-specific risk
  score is \eqn{1-\widehat S_i(t)}, and comparisons are weighted by the
  same grow-data censoring distribution.

  When confidence intervals are requested, the Brier standard error is
  obtained by deleting one subject-level IPCW loss at a time while holding
  the survival predictions and censoring weights fixed. 

  The AUC standard error uses a stratified delete-one jackknife. One observed
  case or control is deleted at a time, the remaining IPCW weights within
  that stratum are renormalized, and the resulting case-delete and
  control-delete variance components are added. With equal weights this
  calculation reduces exactly to the ordinary DeLong variance. The
  effective case and control sample sizes and the largest normalized IPCW
  weights are returned to diagnose horizons at which a few observations
  dominate the weighted comparison.

  Both are pointwise conditional standard errors, and the reported
  intervals are pointwise normal approximations rather than simultaneous
  bands. The survival predictions and the estimated censoring distribution
  are treated as fixed. They do not account for uncertainty from training
  the survival forest or estimating the censoring model. For an independent
  prediction sample, the intervals have a direct conditional test-performance
  interpretation. For a grow object, shared out-of-bag fits and the shared
  censoring estimate induce dependence among subject-level quantities, so
  the intervals should be viewed as fixed-fit working intervals rather than
  repeated-training confidence intervals. A score and its interval are
  returned as \code{NA} at a horizon where a required censoring survival
  probability is zero. The AUC standard error is also \code{NA} if fewer
  than two cases or two controls are available, or if deleting one subject
  leaves no positive IPCW weight in its stratum.

  \code{plotBrierAUC} is a base-R plotting helper. It calls a shared
  numerical engine so that Brier score and AUC use the same fitted
  censoring model, then adds the requested pointwise confidence bands.

  Only right-censored survival families are supported. Competing-risk
  analyses should use \command{plot.competing.risk}.


}
\value{
  \code{get.brier.survival} returns a list containing
  \code{brier.score}, the optional \code{brier.matx}, integrated scores
  \code{crps} and \code{crps.std}, the estimated censoring distribution,
  aligned grow and evaluation event information, survival predictions,
  mortality, and the selected subset. When confidence intervals are
  requested, \code{brier.score} also contains \code{std.err},
  \code{lower}, \code{upper}, and \code{n.eval}.

  \code{get.auct.survival} returns a corresponding list whose
  \code{auct.score} contains \code{time}, \code{auct}, \code{n.case},
  \code{n.control}, effective sample sizes \code{n.case.eff} and
  \code{n.control.eff}, and largest normalized weights
  \code{max.case.weight} and \code{max.control.weight}. The columns
  \code{std.err}, \code{lower}, and \code{upper} are included when
  confidence intervals are requested.

  \code{plotBrierAUC} invisibly returns a list containing the requested
  Brier and/or AUC result objects.

  With evaluation outcomes, \code{plot.survival} invisibly returns the
  mortality-stratified and overall Brier and cumulative integrated Brier
  curves. Without evaluation outcomes, it invisibly returns the predicted
  survival curves and their ensemble mean.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Gerds T.A. and Schumacher M. (2006). Consistent estimation of the
  expected Brier score in general survival models with right-censored
  event times, \emph{Biometrical Journal}, 48:1029--1040.

  Graf E., Schmoor C., Sauerbrei W. and Schumacher M. (1999).
  Assessment and comparison of prognostic classification schemes for
  survival data, \emph{Statistics in Medicine}, 18:2529--2545.

  DeLong E.R., DeLong D.M. and Clarke-Pearson D.L. (1988). Comparing the
  areas under two or more correlated receiver operating characteristic
  curves: a nonparametric approach, \emph{Biometrics}, 44:837--845.

  Efron B. and Tibshirani R.J. (1993). \emph{An Introduction to the
  Bootstrap}. Chapman and Hall, New York.

  Heagerty P.J. and Zheng Y. (2005). Survival model predictive accuracy
  and ROC curves, \emph{Biometrics}, 61:92--105.

  Ishwaran H. and Kogalur U.B. (2007). Random survival forests for R,
  \emph{R News}, 7(2):25--31.

  Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008).
  Random survival forests, \emph{Annals of Applied Statistics},
  2:841--860.
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
plot.survival(rfsrc(Surv(time, status) ~ ., veteran),
              cens.model = "rfsrc")

## pbc data
data(pbc, package = "randomForestSRC")
pbc.obj <- rfsrc(Surv(days, status) ~ ., pbc)

## standard survival diagnostics
plot.survival(pbc.obj)
plot.survival(pbc.obj, subset = c(3, 10), collapse = TRUE)

## Brier and AUCT helpers with pointwise intervals
brier.obj <- get.brier.survival(pbc.obj, conf.int = TRUE)
head(brier.obj$brier.score)
auct.obj <- get.auct.survival(pbc.obj, conf.int = TRUE)
head(auct.obj$auct.score)

## compare two grow-data censoring models
brier.km <- get.brier.survival(pbc.obj, cens.model = "km")
brier.rf <- get.brier.survival(pbc.obj, cens.model = "rfsrc")
plot(brier.km$brier.score$time,
     brier.km$brier.score$brier.score, type = "s", col = 2,
     xlab = "Time", ylab = "Brier Score")
lines(brier.rf$brier.score$time,
      brier.rf$brier.score$brier.score, type = "s", col = 4)
legend("bottomright",
       legend = c("cens.model = km", "cens.model = rfsrc"),
       col = c(2, 4), lty = 1)

## Brier and AUCT curves with confidence bands
plotBrierAUC(pbc.obj)
plotBrierAUC(pbc.obj, plots = "auct", conf.int = 0.90)
plotBrierAUC(pbc.obj, plots = "brier", conf.int = 0.90)



}
}
\keyword{plot}
\keyword{survival}
