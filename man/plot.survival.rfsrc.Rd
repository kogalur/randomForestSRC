\name{plot.survival.rfsrc}
\alias{plot.survival.rfsrc}
\alias{plot.survival}
\alias{get.brier.survival}
\alias{get.auct.survival}
\alias{plotBrierAUC}
\alias{get.cindex}
\title{Survival Prediction Performance and Diagnostic Plots}
\description{
  Display survival predictions and assess their accuracy using Brier
  scores, integrated Brier scores, and time-dependent AUC. Compare the
  forest Brier score with a predictor-free Kaplan--Meier reference,
  or obtain numerical scores for further analysis.

  The separate \code{get.cindex} helper calculates concordance error
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
  \item{x, o}{A fitted or prediction object of class \code{(rfsrc, grow)}
    or \code{(rfsrc, predict)}. The functions \code{get.brier.survival},
    \code{get.auct.survival}, and \code{plotBrierAUC} also accept a saved
    forest of class \code{(rfsrc, forest)}.}

  \item{show.plots}{If \code{TRUE}, draw the requested plots. Use
    \code{FALSE} to return the results without drawing.}

  \item{subset}{Observations to include, specified by positive row
    indices or a logical vector with one entry per predicted observation.
    By default, all observations are included.}

  \item{collapse}{For \code{plot.survival}, show the mean survival curve
    without the individual curves? The marginal Nelson--Aalen curve,
    when available, is also shown.}

  \item{cens.model}{Censoring model for inverse-probability-of-censoring
    weighting (IPCW), estimated from the full training data:
    \describe{
      \item{\code{"km"}:}{A Kaplan--Meier estimate common to all
        observations. This is the default.}
      \item{\code{"rfsrc"}:}{A random survival forest estimate that
        allows the censoring distribution to depend on the predictors.
        When censoring is present, the object must contain both training
        and evaluation predictors. See \bold{Note}.}
    }
  }

  \item{papply}{Function used in place of \code{lapply} for repeated
    calculations.}

  \item{times}{Optional numeric vector of times at which to calculate
    Brier score or AUC. Survival predictions are evaluated as
    right-continuous step functions on this grid. By default, use the
    forest's time grid.}

  \item{conf.int}{Pointwise confidence intervals: \code{FALSE} omits
    intervals, \code{TRUE} requests 95 percent intervals, and a number
    strictly between zero and one specifies another confidence level.
    Defaults to \code{TRUE} for \code{plotBrierAUC} and \code{FALSE}
    for the numerical helpers.}

  \item{keep.matrix}{For \code{get.brier.survival}, include the
    subject-by-time matrix of IPCW Brier contributions in the result?
    Use \code{FALSE} when only the summary scores are needed.}

  \item{plots}{For \code{plotBrierAUC}, select \code{"brier"},
    \code{"auct"}, or both. The name \code{"auc"} is an alias for
    \code{"auct"}.}

  \item{brier.null}{For \code{plotBrierAUC}, compare the forest Brier
    score with a predictor-free Kaplan--Meier estimate from the full
    training data? The default is \code{TRUE}. The reference is also
    returned with \code{show.plots=FALSE}. Ignored when only AUC is
    requested.}

  \item{...}{Further graphical arguments. For \code{plotBrierAUC},
    \code{xlim}, \code{ylim}, titles, and axis options apply to each
    requested panel. The arguments \code{col}, \code{lwd}, \code{lty},
    and \code{type} control the forest curve; \code{band.col} and
    \code{band.border} control its confidence band. Use \code{null.col},
    \code{null.lwd}, and \code{null.lty} for the Brier reference
    (defaults: black, width 2, dashed). The Forest/Null legend is placed
    at \code{legend.pos="topleft"} by default; set \code{legend.pos}
    to another position or \code{FALSE} to omit it.}

  \item{time}{For \code{get.cindex}, observed follow-up times, in the
    same order as \code{censoring} and the prediction rows.}

  \item{censoring}{For \code{get.cindex}, numeric status codes: zero
    for censoring and one for an event in right-censored survival.
    For competing risks, use event codes \code{1}, \ldots, \code{J}.}

  \item{predicted}{For \code{get.cindex}, risk predictions, with larger
    values representing higher risk. For right-censored survival, supply
    a vector or one-column matrix or data frame. For competing risks,
    supply a matrix or data frame with column \code{j} for event code
    \code{j}, or a list of event-specific prediction vectors. The first
    \code{J} columns or list entries are used.}

  \item{weight}{For \code{get.cindex}, optional precomputed concordance
    weights, one per observation. Omit for unweighted concordance.
    This helper does not estimate censoring weights.}

  \item{fast}{For \code{get.cindex}, optional concordance-algorithm
    selector. Usually left unspecified so that the underlying library
    chooses the algorithm. Ignored for weighted competing-risk
    concordance.}

  \item{do.trace}{For \code{get.cindex}, display diagnostic output
    from the underlying concordance calculation?}
}
\details{
  These functions assess survival predictions over follow-up. Brier
  score measures agreement between predicted survival probabilities
  and observed outcomes. Time-dependent AUC measures how well the
  predictions distinguish individuals with an event by a given time
  from those still event-free.

  Use \code{plot.survival} for a four-panel diagnostic display,
  \code{plotBrierAUC} for performance curves with confidence bands,
  or \code{get.brier.survival} and \code{get.auct.survival} for numerical
  results and custom plots.

  These plots and scores support right-censored survival. For
  competing-risk plots, use \code{\link{plot.competing.risk.rfsrc}}.
  The separate \code{get.cindex} helper supports both settings.

  \subsection{Diagnostic panels}{
    \code{plot.survival} draws two rows of panels, from left to right:
    \enumerate{
      \item \bold{Survival curves.} Individual forest predictions, their
      mean across the selected individuals (thick red line), and the
      marginal Nelson--Aalen survival estimate from their outcomes
      (thick green line).

      \item \bold{Brier score.} Prediction error over time for four
      groups defined by the 0--25, 25--50, 50--75, and 75--100 percentile
      ranges of predicted mortality. The red curve is the overall
      Brier score.

      \item \bold{Integrated Brier score.} Cumulative average Brier
      error from the first plotting time to each subsequent time,
      labeled CRPS. The groups and red overall curve match the Brier
      panel.

      \item \bold{Mortality versus observed time.} Predicted mortality
      plotted against observed follow-up time, with blue points for
      events and black points for censored observations. Mortality
      (Ishwaran et al., 2008) expresses risk on the scale of the number
      of events. For example, a value of 100 means that an average of
      100 events would be expected if all individuals had the same
      covariate values.
    }
  }

  \subsection{Training and test evaluation}{
    For a grow object, the functions use out-of-bag (OOB) survival and
    mortality predictions when both are available; otherwise, they use
    full-ensemble predictions. For a prediction object, they use its
    test predictions and corresponding test outcomes.

    The censoring model always uses the full training data, including
    when scoring a test set. Training outcomes must therefore be
    available. The argument \code{subset} selects the observations to
    evaluate, not the data used to estimate censoring.

    If no evaluation outcomes are available, \code{plot.survival}
    displays only predicted survival curves and their mean.
    \code{get.brier.survival}, \code{get.auct.survival}, and
    \code{plotBrierAUC} require outcomes for the observations being
    evaluated.
  }

  \subsection{Brier score and integrated error}{
    The Brier score uses squared prediction errors adjusted for
    censoring by IPCW (Graf et al., 1999; Gerds and Schumacher, 2006).
    Smaller scores indicate better predictions.

    Let \eqn{T_i} be the observed follow-up time, \eqn{\Delta_i} the
    event status (zero for censoring), \eqn{\widehat S_i(t)} the predicted
    survival probability, and \eqn{\widehat G(t)} the estimated censoring
    survival probability. The contribution from individual \eqn{i} is
  \deqn{
    L_i(t) =
    \frac{I(T_i \le t, \Delta_i > 0)}{\widehat G(T_i-)}
       \widehat S_i(t)^2
    +
    \frac{I(T_i > t)}{\widehat G(t)}
       \{1-\widehat S_i(t)\}^2,
  }
    For \code{cens.model="rfsrc"}, \eqn{\widehat G} also depends on the
    individual's predictors. The Brier score averages the finite
    contributions at each time. Individuals censored at or before
    \eqn{t} contribute zero and are included in this average.

    Integrated Brier scores use trapezoidal integration. The CRPS
    panel in \code{plot.survival} divides the integral through each time
    by the elapsed time since the first grid point. Its first value
    is \code{NA} because no time has elapsed.

    The scalar summaries from \code{get.brier.survival} are \code{crps},
    the integral over the evaluation grid, and \code{crps.std}, that
    integral divided by the largest evaluation time. This differs from
    the elapsed-time normalization in the CRPS panel.
  }

  \subsection{Time-dependent AUC}{
    The cumulative/dynamic AUC compares individuals with an observed
    event by time \eqn{t} (cases) with those followed beyond \eqn{t}
    (controls). Individuals censored at or before \eqn{t} are excluded
    from this comparison. Comparisons use the risk score
    \eqn{1-\widehat S_i(t)} and IPCW weights from the training-data
    censoring model.

    A larger AUC indicates better discrimination; \code{plotBrierAUC}
    shows a reference line at 0.5. The effective sample sizes and
    largest normalized weights in the output help identify times
    when a few observations dominate the comparison.
  }

  \subsection{Null reference for the Brier score}{
    By default, \code{plotBrierAUC} compares the forest with a
    Kaplan--Meier estimate from the full training sample. This reference
    assigns the same survival curve to every individual, without using
    predictors. A forest Brier curve below the dashed reference indicates
    lower estimated prediction error.

    The reference survival curve is
    \deqn{
      \widehat S_0(t) =
      \prod_{u \le t}\left\{1-\frac{d_{\mathrm{grow}}(u)}
                                      {Y_{\mathrm{grow}}(u)}\right\},
    }
    where \eqn{d_{\mathrm{grow}}(u)} and
    \eqn{Y_{\mathrm{grow}}(u)} are the event count and risk-set size
    at training event time \eqn{u}. This estimates event survival;
    \code{cens.model} separately determines the censoring distribution.

    The reference substitutes \eqn{\widehat S_0(t)} for
    \eqn{\widehat S_i(t)} in the Brier loss. Both curves use the same
    times, censoring weights, and observations with finite forest
    losses, including zero contributions from individuals already
    censored. Confidence bands apply only to the forest; no test of
    the difference between curves is performed.

    The reference uses training outcomes only, regardless of
    \code{subset}. For training evaluation, it is an in-sample benchmark
    even when the forest predictions are OOB. Set \code{brier.null=FALSE}
    to omit it.
  }

  \subsection{Confidence intervals}{
    Confidence intervals use pointwise normal approximations, treating
    the survival predictions and censoring estimates as fixed. They
    describe evaluation uncertainty, excluding uncertainty from fitting
    either model. The shaded regions in \code{plotBrierAUC} show these
    pointwise intervals, not simultaneous confidence bands.

    \describe{
      \item{Brier score:}{The standard error is the sample standard
        deviation of the finite individual IPCW losses divided by the
        square root of their count. This is equivalent to deleting one
        loss at a time while keeping the predictions and weights fixed.}
      \item{AUC:}{A stratified delete-one jackknife removes one case or
        control at a time and renormalizes the remaining weights within
        that group. The case and control variance components are added.
        With equal weights, this gives the DeLong variance
        (DeLong et al., 1988).}
    }

    For an independent test sample, the intervals describe conditional
    test performance. For a grow object, shared OOB fits and censoring
    estimates create dependence between observations. Interpret these
    as fixed-fit working intervals, not repeated-training confidence
    intervals.
  }

  \subsection{Concordance error}{
    \code{get.cindex} uses supplied outcomes and risk predictions
    directly. Supply \code{predicted.oob} for OOB evaluation, or
    \code{predicted} with the corresponding new-data outcomes for test
    evaluation. Concordance error is distinct from the time-specific
    AUC returned by \code{get.auct.survival}.

    A status code greater than one selects competing-risk concordance.
    Here, \code{J} is the largest usable status code, and prediction
    column \code{j} must correspond to event code \code{j}; column names
    do not determine this match. Without weights, each event-specific
    calculation includes only censored observations and observations
    with that event. With weights, the event-specific weighted
    concordance calculation is used.

    Missing times, statuses, predictions, and supplied weights are
    excluded from the applicable calculation.
  }
}
\value{
  \subsection{\code{plot.survival}}{
    With evaluation outcomes, invisibly returns a data frame with
    \code{time}, mortality-stratified Brier curves \code{bs.q25},
    \code{bs.q50}, \code{bs.q75}, \code{bs.q100}, and the overall
    \code{bs.all}. The corresponding integrated curves are
    \code{crps.q25}, \code{crps.q50}, \code{crps.q75}, \code{crps.q100},
    and \code{crps.all}.

    Without outcomes, invisibly returns a list with \code{time},
    \code{survival} (a time-by-individual matrix), \code{survival.mean},
    and \code{subset}.
  }

  \subsection{\code{get.brier.survival}}{
    Returns a list containing:
    \describe{
      \item{\code{brier.score}}{Data frame with \code{time} and
        \code{brier.score}. With intervals, also includes \code{std.err},
        \code{lower}, \code{upper}, and \code{n.eval}.}
      \item{\code{brier.matx}}{Subject-by-time IPCW losses; \code{NULL}
        when \code{keep.matrix=FALSE}.}
      \item{\code{crps}, \code{crps.std}}{Integrated Brier scores
        defined in \bold{Details}.}
      \item{\code{n.eval}}{Number of finite Brier contributions at each
        time, including zeros from censored observations.}
    }
    Supporting components include the time grid, censoring distribution,
    training and evaluation event information, survival predictions,
    mortality, and selected subset.
  }

  \subsection{\code{get.auct.survival}}{
    Returns a list with \code{auct.score} and the same supporting
    components as the Brier helper. The \code{auct.score} data frame
    contains \code{time}, \code{auct}, counts \code{n.case} and
    \code{n.control}, effective sample sizes \code{n.case.eff} and
    \code{n.control.eff}, and largest normalized weights
    \code{max.case.weight} and \code{max.control.weight}. With intervals,
    also includes \code{std.err}, \code{lower}, and \code{upper}.
  }

  \subsection{\code{plotBrierAUC}}{
    Invisibly returns the requested \code{brier} and/or \code{auct}
    result lists, including with \code{show.plots=FALSE}. The Brier
    result has \code{brier.matx=NULL}.

    With \code{brier.null=TRUE}, \code{brier$null} contains reference
    scores \code{brier.score} (columns \code{time} and \code{brier.score}),
    \code{crps}, \code{crps.std}, and contribution counts \code{n.eval}.
    It also contains the reference \code{survival} vector, \code{n.train}
    (finite training outcomes), \code{method="kaplan-meier"}, and
    \code{source="grow"}. Integration and standardization match the
    forest scores. Individual reference losses and reference confidence
    intervals are not returned.
  }

  \subsection{\code{get.cindex}}{
    Returns a numeric concordance error for right-censored survival or
    a vector in event-code order \code{1}, \ldots, \code{J} for competing
    risks. In unweighted competing-risk calculations, an event with
    fewer than two usable observations has error \code{NA}.
  }
}
\note{
  When censoring is present, \code{cens.model="rfsrc"} requires training
  and evaluation predictors. For a reduced \code{rfsrc.fast} object
  fitted with \code{forest=FALSE}, use \code{cens.model="km"} or refit
  with \code{forest=TRUE}.

  A score and its interval are \code{NA} when a required censoring
  probability is zero or nonfinite. AUC also requires at least one
  case and one control; its standard error requires at least two of
  each and positive remaining weight after each deletion. The Brier
  standard error requires at least two finite losses. The null Brier
  score is \code{NA} if the forest score is unavailable or a required
  reference loss is nonfinite.

  The Kaplan--Meier reference uses all finite training outcomes,
  counts repeated observations at ties, and processes events before
  censoring at a tie. Its final estimate applies beyond the last
  observed time, provided the required censoring probabilities remain
  positive and finite.

  With \code{collapse=FALSE} and more than 500 selected individuals,
  \code{plot.survival} draws 500 randomly selected individual curves.
  The mean curve and scores still use the full selected group.
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

## Show only the forest Brier curve and its confidence band
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
