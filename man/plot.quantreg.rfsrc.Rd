\name{plot.quantreg.rfsrc}
\alias{plot.quantreg.rfsrc}
\alias{plot.quantreg}
\title{Plot Conditional Quantiles and CRPS Diagnostics}
\description{
  Plots estimated conditional response quantiles against observed responses.
  An optional inset compares the forest's continuous ranked probability
  score (CRPS) with the score for the empirical CDF of training responses.
  Mean pinball losses can also be displayed for selected quantiles.
}
\usage{
\method{plot.quantreg}{rfsrc}(x, prbL = .25, prbU = .75,
  m.target = NULL, crps = TRUE, subset = NULL,
  xlab = NULL, ylab = NULL, ...,
  inset.args = list(), crps.null = TRUE, quantreg.tau = NULL)
}
\arguments{
  \item{x}{A training or prediction object returned by
    \code{\link{quantreg}}, with observed responses for the observations
    to be plotted. For new-data diagnostics, include response columns
    in the \code{quantreg()} prediction call.}
  \item{prbL, prbU}{Lower and upper quantile levels: single probabilities
    in \eqn{(0,1)} with \code{prbL < prbU}. The middle level is 0.5 if
    it lies strictly between them, or their midpoint otherwise.
    \code{get.quantile()} matches each level to the nearest stored
    probability. Include the desired levels in \code{prob} when calling
    \code{quantreg()} for exact matching.}
  \item{m.target}{Name of the continuous response to display. Defaults
    to the first response with quantile results. Selects existing
    results without calculating new predictions.}
  \item{crps}{If \code{TRUE}, show the standardized CRPS curve from
    \code{get.quantile.crps()} in an inset. Scores use the finite CDF
    reporting grid.}
  \item{subset}{Rows to plot and score. Supply positive integer indices
    or a logical vector with one nonmissing entry per row of the
    quantile output. \code{NULL} uses all rows. Indices refer to results
    after training or prediction preprocessing; repeated indices count
    a row more than once.}
  \item{xlab, ylab}{Main-panel axis labels. Defaults are the selected
    response name and \code{"Target Quantiles"}.}
  \item{...}{Named graphical arguments passed to
    \code{\link[graphics]{plot.default}} for the main panel, including
    \code{xlim}, \code{ylim}, \code{main}, \code{cex.axis}, and
    \code{las}. Use \code{pch}, \code{col}, and \code{cex} to style
    the middle-quantile points. The function sets the coordinates and
    initial panel \code{type}. These arguments do not style the inset.}
  \item{inset.args}{Named list of graphical arguments passed to
    \code{\link[graphics]{plot.default}} for the CRPS inset. Use
    \code{xlim} and \code{ylim} for limits; \code{col}, \code{lty},
    and \code{lwd} for the forest curve; and \code{cex.axis},
    \code{las}, \code{tck}, or related options for the axes.
    \code{axes = FALSE} hides both axes; \code{yaxt = "n"} hides
    the right-hand axis. The default \code{list()} uses the inset's
    own settings, with limits covering both curves. Ignored when
    \code{crps = FALSE}. Supply the full argument name.}
  \item{quantreg.tau}{Quantile levels for mean pinball-loss labels, as
    a vector of finite probabilities strictly between zero and one.
    Losses use the selected response and rows and appear in the
    lower-right corner of the main panel. \code{NULL} omits labels
    even when the session option \code{quantreg.tau} is set. The plotted
    quantiles and CRPS inset are unchanged. Supply the full argument name.}
  \item{crps.null}{If \code{TRUE}, add a dashed black null curve and a
    \code{Forest}/\code{Null} legend to the inset. \code{FALSE} shows
    the forest curve alone. Training responses for the null distribution
    come from \code{x$forest$yvar}, or from \code{x$yvar} for a training
    object without responses saved in its forest. If unavailable, omit
    the null curve with a warning. Ignored when \code{crps = FALSE}.
    Supply the full argument name.}
}
\details{
  Each observation is plotted at its observed response on the horizontal
  axis. A point marks its predicted middle quantile, and a vertical
  segment with endpoint marks spans its lower and upper quantiles.
  The segment is a prediction interval for an individual response,
  not a confidence interval for the conditional mean.
  The dashed diagonal marks equality between observed and predicted
  values. Defaults request the median and the 25th and 75th percentiles.
  The response selected by \code{m.target} and rows selected by
  \code{subset} apply to the panel and all scores.

  Horizontal jitter separates overlapping observations; scores use the
  original response values. Rows with an unavailable response or plotted
  quantile are omitted from the panel with a warning.

  To add pinball losses, supply levels such as
  \code{quantreg.tau = c(.1, .5, .9)}. They can differ from the levels
  drawn in the panel. Losses are calculated by
  \code{get.pinball.error(x, tau = quantreg.tau, subset = subset,
  m.target = m.target)}. Each mean uses selected rows with finite
  responses and quantile predictions at that level. Smaller losses
  indicate better quantile predictions.

  The inset shows the standardized CRPS curve from
  \code{get.quantile.crps()}. At each response threshold on the
  horizontal axis, the right-hand axis shows prediction error integrated
  from the first grid value to that threshold and divided by the
  integration width. Smaller values indicate less error. Scores depend
  on the grid's range and resolution; see \code{\link{quantreg}} for
  the calculation. The first value is unavailable because its integration
  width is zero. An entirely unavailable curve is omitted with a warning.

  The null distribution is the empirical CDF of all finite training
  responses for the selected outcome and is unchanged by \code{subset}.
  Its CRPS curve uses the same grid, observations, integration, and
  standardization as the forest curve. At each threshold, both scores
  exclude rows with an unavailable response or forest CDF prediction.
  For training plots, the null distribution includes the evaluated
  responses and is an in-sample reference, even for OOB predictions.

  Use \code{...} to customize the main panel and \code{inset.args} for
  the inset. Axis limits change only the displayed region, not the rows
  scored or the CRPS integration limits. After drawing the inset,
  graphics settings return to the main panel so that further drawing
  applies there. Successive calls can be used in a multi-panel layout.
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
