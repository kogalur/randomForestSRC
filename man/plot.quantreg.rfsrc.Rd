\name{plot.quantreg.rfsrc}
\alias{plot.quantreg.rfsrc}
\alias{plot.quantreg}
\title{Plot Conditional Quantiles and CRPS Diagnostics}
\description{
  Plots observed responses against their predicted conditional quantiles.
  An optional inset compares the forest's continuous ranked probability
  score (CRPS) with a reference that uses the training-response
  distribution without predictors. Numerical pinball losses can also
  be displayed for selected quantile levels.
}
\usage{
\method{plot.quantreg}{rfsrc}(x, prbL = .25, prbU = .75,
  m.target = NULL, crps = TRUE, subset = NULL,
  xlab = NULL, ylab = NULL, ...,
  inset.args = list(), crps.null = TRUE, quantreg.tau = NULL)
}
\arguments{
  \item{x}{A training or prediction object returned by
    \code{\link{quantreg}}, containing observed responses for the
    observations to be plotted. For new-data diagnostics, include
    response columns in the \code{quantreg()} prediction call.}
  \item{prbL, prbU}{Lower and upper quantile levels, each a single
    probability in \eqn{(0,1)}, with \code{prbL < prbU}. The middle
    level is 0.5 when it lies strictly between them; otherwise it is
    their midpoint. Each level is matched to the nearest stored
    probability by \code{get.quantile()}. Include the desired levels
    in \code{prob} when calling \code{quantreg()} for exact matching.}
  \item{m.target}{Name of one continuous response to display. The default
    is the first response with quantile results. For multivariate or
    mixed outcomes, this selects results from \code{x} without
    calculating new predictions.}
  \item{crps}{If \code{TRUE}, display the standardized, finite-grid CRPS
    curve from \code{get.quantile.crps()} in an inset.}
  \item{subset}{Observations to plot and score. Supply positive integer
    row indices or a logical vector with one nonmissing entry per row
    of the quantile output. \code{NULL} uses all rows. Indices refer to
    the returned results after training or prediction preprocessing;
    repeated indices include an observation more than once.}
  \item{xlab, ylab}{Main-panel axis labels. Defaults are the selected
    response name and \code{"Target Quantiles"}.}
  \item{...}{Named graphical arguments for the main panel, passed to
    \code{\link[graphics]{plot.default}}. Use \code{xlim} and
    \code{ylim} for axis limits, \code{main} for the title, and
    \code{cex.axis}, \code{las}, and related options for the axes.
    The arguments \code{pch}, \code{col}, and \code{cex} style the
    middle-quantile points. The function sets the plotted coordinates
    and initial panel \code{type}. These settings do not style the
    CRPS inset.}
  \item{inset.args}{Named list of graphical arguments for the CRPS
    inset, passed to \code{\link[graphics]{plot.default}}. Use
    \code{xlim} and \code{ylim} for inset limits, and \code{col},
    \code{lty}, and \code{lwd} for the forest curve. Axis options such
    as \code{cex.axis}, \code{las}, and \code{tck} also apply.
    \code{axes = FALSE} hides both inset axes; \code{yaxt = "n"}
    hides its right-hand axis. The default \code{list()} uses the
    inset's own settings. Ignored when \code{crps = FALSE}. Supply
    the full argument name.}
  \item{quantreg.tau}{Quantile levels for numerical pinball-loss labels,
    as a vector of finite probabilities strictly between zero and one.
    For each level, the plot displays the mean loss for the selected
    response and subset in the lower-right corner of the main panel.
    The default \code{NULL} omits these labels, even when the session
    option \code{quantreg.tau} is set. The plotted quantiles and CRPS
    inset are unchanged. Supply the full argument name.}
  \item{crps.null}{If \code{TRUE}, include the null reference in the
    CRPS inset as a dashed black curve, with a legend identifying
    \code{Forest} and \code{Null}. Set to \code{FALSE} for the forest
    curve alone. Ignored when \code{crps = FALSE}. Supply the full
    argument name.}
}
\details{
  The plot compares observed responses with the predicted response
  distributions. The main panel shows selected quantiles for individual
  observations; the CRPS inset summarizes distributional prediction
  error across observations. Both displays use the response selected
  by \code{m.target} and the observations selected by \code{subset}.

  \subsection{Quantile panel}{
    Each observation is placed at its observed response on the
    horizontal axis. A point marks its predicted middle quantile, and a
    vertical segment with endpoint marks spans its lower and upper
    quantiles. The dashed diagonal marks equality between observed and
    predicted values. The default levels request the median and an
    interval from the 25th to the 75th percentile.

    Horizontal jitter separates overlapping observations. All scores
    use the original response values, not these display positions.
    Observations with an unavailable response or requested quantile
    are omitted from the main panel with a warning.
  }

  \subsection{Pinball-loss labels}{
    Supply \code{quantreg.tau = c(.1, .5, .9)} to display mean pinball
    losses for the 10th, 50th, and 90th percentile predictions as text
    in the lower-right corner of the main panel. Smaller losses
    indicate better quantile predictions. These levels can differ from
    the quantiles drawn in the panel.

    The values are calculated by
    \code{get.pinball.error(x, tau = quantreg.tau, subset = subset,
    m.target = m.target)} for the selected response. Each loss averages
    over the selected observations with finite responses and quantile
    predictions at that level. Plotting omits these labels unless
    levels are supplied in the call; the session reporting option
    alone does not add them.
  }

  \subsection{CRPS inset}{
    The inset evaluates the predicted CDF across response thresholds.
    At each threshold, it averages the squared difference between the
    predicted CDF and the indicator that the observed response is at
    or below that threshold. It integrates these errors by the
    trapezoidal rule from the first reporting-grid value to each
    successive value, then divides by the integration width.

    The horizontal axis is the response threshold; the right-hand
    vertical axis is standardized CRPS. Lower values indicate less
    error over the corresponding integration interval. This is the
    finite-grid curve returned by \code{get.quantile.crps()}, so the
    grid's range and resolution affect the scores. The first value is
    unavailable because its integration width is zero. An entirely
    unavailable curve is omitted with a warning.
  }

  \subsection{Null reference}{
    The null curve provides a comparison with predictions based on the
    training-response distribution alone. It uses the same empirical
    CDF for every evaluated observation:
    \deqn{F_0(t)=\frac{1}{n_0}\sum_{j=1}^{n_0}
      I(Y_j^{\mathrm{train}}\leq t),}
    where \eqn{n_0} is the number of finite training responses for the
    selected outcome. Every training observation contributes, so
    repeated response values contribute according to their frequencies.
    Test responses are used only to evaluate the reference, not to
    construct it.

    Forest and null curves use the same response grid, subset,
    integration rule, and standardization. At each threshold, both
    exclude observations with an unavailable response or forest CDF
    prediction. The null distribution itself always uses all finite
    training responses, regardless of \code{subset}. For training
    plots, it includes the evaluated responses and is an in-sample
    benchmark rather than an OOB or leave-one-out estimate.

    Training responses come from \code{x$forest$yvar}, or
    \code{x$yvar} for a training object without responses saved in its
    forest. If they are unavailable, the null curve is omitted with
    a warning.
  }

  \subsection{Customizing the display}{
    Use \code{...} for the main panel and \code{inset.args} for the
    inset. For example,
    \code{inset.args = list(ylim = c(0, 0.3), lwd = 1.5)} changes the
    inset's vertical limits and forest-curve line width. Default inset
    limits cover both curves.

    Axis limits change only the displayed region; they do not change
    the observations scored or the CRPS integration limits. After the
    inset is drawn, graphics settings return to the main panel so that
    further drawing applies there. Successive calls can be used in a
    multi-panel layout.
  }
}
\value{
  Used for its graphical side effect. Invisibly returns \code{NULL}.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\seealso{
  \code{\link{quantreg}}, \code{\link{get.quantile}},
  \code{\link{get.quantile.crps}}, \code{\link{get.pinball.error}},
  \code{\link[graphics]{plot.default}}
}
\examples{
\donttest{
## Univariate quantiles, with forest and null CRPS curves.
library(randomForestSRC)
set.seed(19)
dta <- na.omit(airquality)
prob <- c(.25, .50, .75)
q <- quantreg(Temp ~ ., data = dta, prob = prob, ntree = 100)
plot.quantreg(q)

## Add pinball losses for this plot only; the interval and inset are unchanged.
plot.quantreg(q, quantreg.tau = c(.2, .5, .8))
print(get.pinball.error(q, tau = c(.2, .5, .8)))

## Customize the main panel directly, and the inset separately.
plot.quantreg(q, main = "Temperature",
              xlim = c(50, 100), ylim = c(45, 110), pch = 19,
              inset.args = list(xlim = c(60, 95), ylim = c(0, .35),
                                lwd = 1.5))

## Keep the inset, but suppress its null reference.
plot.quantreg(q, crps.null = FALSE)

## Multivariate forest: select each response for plotting only.
mv <- quantreg(cbind(Ozone, Temp) ~ ., data = dta,
               splitrule = "mahalanobis", prob = prob, ntree = 100)
print(names(get.quantile(mv, pretty = FALSE)))
op <- par(mfrow = c(1, 2))
plot.quantreg(mv, m.target = "Ozone")
plot.quantreg(mv, m.target = "Temp", main = "Temperature",
              xlim = c(50, 100), ylim = c(45, 110),
              inset.args = list(ylim = c(0, .35)))
par(op)
}
}
\keyword{plot}
