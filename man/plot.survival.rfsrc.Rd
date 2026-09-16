\name{plot.survival.rfsrc}
\alias{plot.survival.rfsrc}
\alias{plot.survival}
\alias{get.brier.survival}
\alias{get.auct.survival}
\alias{plotBrierAUC}
\alias{get.cindex}
\title{Survival Prediction Performance and Diagnostic Plots}
\description{
  Plot survival estimates and calculate inverse-probability-of-censoring
  weighted Brier scores and cumulative/dynamic time-dependent AUC.
  \code{plotBrierAUC} also compares the forest Brier curve with a
  predictor-free Kaplan--Meier reference estimated from the training data.
  \code{get.cindex} calculates concordance-based prediction error
  from observed outcomes and risk predictions for right-censored
  survival or competing risks.
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
  plots = c("brier", "auct"), show.plots = TRUE,
  brier.null = TRUE, ...)

get.cindex(time, censoring, predicted, weight, fast, do.trace = FALSE)
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

  \item{brier.null}{For \code{plotBrierAUC}, calculate and display the
    Brier score of a predictor-free Kaplan--Meier survival estimate from
    the full grow data? The default is \code{TRUE}. The reference is also
    returned when \code{show.plots=FALSE}. This argument has no effect
    when only AUC is requested.}

  \item{...}{Further graphical arguments. For \code{plotBrierAUC},
    \code{xlim}, \code{ylim}, titles, and axis options apply to each
    requested panel. The arguments \code{col}, \code{lwd}, \code{lty},
    and \code{type} control the forest curve; \code{band.col} and
    \code{band.border} control its confidence band. Use \code{null.col},
    \code{null.lwd}, and \code{null.lty} to style the Brier reference
    (defaults: black, width 2, dashed). Use \code{legend.pos} to position
    its Forest/Null legend (default \code{"topright"}), or
    \code{legend.pos=FALSE} to omit the legend.}
  \item{time}{For \code{get.cindex}, the observed follow-up times,
    aligned with \code{censoring} and the prediction rows.}
  \item{censoring}{For \code{get.cindex}, numeric status codes: zero
    for censoring and one for an event in right-censored survival.
    For competing risks, use event codes \code{1}, \ldots, \code{J}.}
  \item{predicted}{For \code{get.cindex}, risk predictions, with larger
    values representing higher risk. For right-censored survival, a
    vector or one-column matrix or data frame. For competing risks,
    a matrix or data frame with column \code{j} corresponding to event
    code \code{j}, or a list of event-specific prediction vectors.
    The first \code{J} columns or list entries are used.}
  \item{weight}{Optional numeric vector of precomputed concordance
    weights, one per observation. Omit it for the unweighted calculation.
    The helper does not estimate censoring weights from the supplied
    outcomes.}
  \item{fast}{Optional native concordance-algorithm selector. Normally
    leave it unspecified so that the native library chooses the
    implementation. It is not used by the weighted competing-risk branch.}
  \item{do.trace}{For \code{get.cindex}, enable native diagnostic output?}
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
  numerical engine so that Brier score, its null reference, and AUC use
  the same censoring model, then adds the requested pointwise confidence
  bands for the forest curves. The AUC panel retains its horizontal
  reference at 0.5.

  \subsection{Null reference for the Brier score}{
    With \code{brier.null=TRUE}, the Brier panel includes a dashed
    predictor-free reference and a Forest/Null legend. Its survival
    prediction is the same for every evaluated subject:
    \deqn{
      \widehat S_0(t) =
      \prod_{u \le t}\left\{1-\frac{d_{\mathrm{grow}}(u)}
                                      {Y_{\mathrm{grow}}(u)}\right\},
    }
    where \eqn{d_{\mathrm{grow}}(u)} and
    \eqn{Y_{\mathrm{grow}}(u)} are the event count and risk-set size
    at grow-data event time \eqn{u}. All finite grow outcomes are used,
    including their frequencies at tied times. Events are processed
    before censorings at a tied time. This is the event-survival
    estimate; \code{cens.model} separately controls the censoring
    distribution used for scoring.

    The null Brier curve is obtained by replacing \eqn{\widehat S_i(t)}
    with \eqn{\widehat S_0(t)} in the IPCW loss above. It uses the same
    horizons, evaluation subset, censoring weights, and finite forest-loss
    contributions as the forest curve. Censored observations at or before
    a horizon retain their zero contribution in both averages. A reference
    horizon is \code{NA} if the forest horizon is unsupported or a required
    null loss is nonfinite. Confidence bands are drawn only for the forest.

    For test evaluation, the reference is estimated from the saved grow
    outcomes, not from the test outcomes. Changing \code{subset} changes
    the cases being scored, not the grow-data reference. For training
    evaluation, the reference uses the full training sample and is an
    in-sample benchmark, even when the forest predictions are OOB.
    The last Kaplan--Meier estimate is carried forward after its last
    observed time, subject to the existing censoring-support checks.

    A forest curve below the null curve has lower estimated prediction
    error than this predictor-free reference at that horizon. The two
    curves are compared directly on the Brier scale; the score is not
    rescaled and no comparison test is added.
  }

  Only right-censored survival families are supported. Competing-risk
  analyses should use \command{plot.competing.risk}.


  \subsection{Concordance error}{
    \code{get.cindex} accepts aligned outcome and risk-prediction values
    directly. Supply \code{predicted.oob} for OOB evaluation or
    \code{predicted} with the corresponding outcomes for new-data
    evaluation. It uses the supplied values without fitting a forest or
    selecting an OOB component automatically.

    A status code greater than one selects the competing-risk branch.
    For that branch, \code{J} is the largest usable status code and
    prediction columns are indexed by event code, not matched by their
    names. Without weights, the calculation for event \code{j} uses
    observations censored or experiencing event \code{j}; observations
    with other events are excluded. With weights, it uses the native
    event-specific weighted concordance calculation instead.

    Missing times, statuses, predictions, and supplied weights are
    excluded from the applicable calculation. The helper returns the
    native concordance error; it is distinct from the time-specific
    cumulative/dynamic AUC returned by \code{get.auct.survival}.
  }
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
  Brier and/or AUC result objects. When Brier is requested with
  \code{brier.null=TRUE}, the Brier result also contains a \code{null}
  list with \code{brier.score} (columns \code{time} and
  \code{brier.score}), integrated scores \code{crps} and
  \code{crps.std}, and the number of evaluated contributions
  \code{n.eval}. It also contains the common null \code{survival}
  vector on the scoring grid, \code{n.train} (finite grow outcomes),
  \code{method="kaplan-meier"}, and \code{source="grow"}.
  No subject-by-time null loss matrix or null confidence band is returned.
  The integrated null scores use the same trapezoidal rule and
  standardization as the forest scores. Direct numerical-helper calls
  retain their existing return structure.

  With evaluation outcomes, \code{plot.survival} invisibly returns the
  mortality-stratified and overall Brier and cumulative integrated Brier
  curves. Without evaluation outcomes, it invisibly returns the predicted
  survival curves and their ensemble mean.

  \code{get.cindex} returns a numeric concordance-error value for
  right-censored survival, or a numeric vector in event-code order
  \code{1}, \ldots, \code{J} for competing risks. In the unweighted
  competing-risk calculation an event with fewer than two usable
  observations returns \code{NA}.
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

## ------------------------------------------------------------
## Concordance error from the stored OOB risk predictions
## ------------------------------------------------------------
print(get.cindex(pbc.obj$yvar[, 1], pbc.obj$yvar[, 2],
                 pbc.obj$predicted.oob))

## standard survival diagnostics
plot.survival(pbc.obj)
plot.survival(pbc.obj, subset = c(3, 10), collapse = TRUE)

## Brier and AUCT helpers with pointwise intervals
brier.obj <- get.brier.survival(pbc.obj, conf.int = TRUE)
print(head(brier.obj$brier.score))
auct.obj <- get.auct.survival(pbc.obj, conf.int = TRUE)
print(head(auct.obj$auct.score))

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

## Brier and AUCT curves; Brier includes the grow-data KM null reference
perf <- plotBrierAUC(pbc.obj)
print(head(data.frame(
  time = perf$brier$time,
  forest = perf$brier$brier.score$brier.score,
  null = perf$brier$null$brier.score$brier.score
)))
plotBrierAUC(pbc.obj, plots = "auct", conf.int = 0.90)
plotBrierAUC(pbc.obj, plots = "brier", conf.int = 0.90,
             ylim = c(0, .4), null.lty = 3, legend.pos = "topleft")

## Omit the Brier reference, retaining the forest curve and its band
plotBrierAUC(pbc.obj, plots = "brier", brier.null = FALSE)

## Independent test evaluation: the null still uses training outcomes
set.seed(19)
pbc.complete <- na.omit(pbc)
trn <- sample(seq_len(nrow(pbc.complete)),
              size = floor(0.7 * nrow(pbc.complete)))
grow <- rfsrc(Surv(days, status) ~ ., pbc.complete[trn, ], ntree = 100)
test <- predict(grow, newdata = pbc.complete[-trn, ])
test.perf <- plotBrierAUC(test, plots = "brier")
print(c(forest = test.perf$brier$crps.std,
        null = test.perf$brier$null$crps.std))

## Obtain both curves without opening a graphics device
perf <- plotBrierAUC(test, plots = "brier", show.plots = FALSE)
print(head(perf$brier$null$brier.score))



}
}
\keyword{plot}
\keyword{survival}
