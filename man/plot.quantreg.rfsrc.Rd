\name{plot.quantreg.rfsrc}
\alias{plot.quantreg.rfsrc}
\alias{plot.quantreg}
\title{Plot Conditional Quantiles and CRPS Diagnostics}
\description{
  Plots observed responses against their conditional quantiles. An optional
  inset compares the forest's continuous ranked probability score (CRPS)
  curve with a predictor-free reference based on the training responses.
}
\usage{
\method{plot.quantreg}{rfsrc}(x, prbL = .25, prbU = .75,
  m.target = NULL, crps = TRUE, subset = NULL,
  xlab = NULL, ylab = NULL, ...,
  inset.args = list(), crps.null = TRUE, quantreg.tau = NULL)
}
\arguments{
  \item{x}{A quantile regression object returned by \code{\link{quantreg}}.
    Observed responses must be retained in the object. Predictions made
    without response values cannot be displayed by this diagnostic.}
  \item{prbL, prbU}{Single probabilities in \eqn{(0,1)} specifying the lower
    and upper conditional quantiles. Require \code{prbL < prbU}. The middle
    probability is \code{0.5} when it lies strictly between these bounds;
    otherwise it is their midpoint. Each requested probability is matched
    to the nearest stored probability by \code{get.quantile()}. Include
    the desired probabilities in the original \code{quantreg()} call for
    exact matching.}
  \item{m.target}{Name of one continuous response to display. The default
    is the first response with stored quantiles. For multivariate or mixed
    outcomes, this selects existing results for plotting and does not
    request new predictions.}
  \item{crps}{Logical. If \code{TRUE}, add the standardized, finite-grid
    CRPS curve returned by \code{get.quantile.crps()}.}
  \item{subset}{Positive integer row indices, or a logical vector with one
    nonmissing entry per row of the quantile output. The default uses all
    rows. Repeated integer indices repeat those observations in the plot
    and score calculation. Row indices refer to the stored results, after
    any preprocessing performed during training or prediction.}
  \item{xlab, ylab}{Axis labels for the main quantile panel. Defaults are
    the selected response name and \code{"Target Quantiles"}.}
  \item{...}{Named graphical arguments for the main quantile panel, passed
    to \code{\link[graphics]{plot.default}}. In particular, \code{xlim}
    and \code{ylim} set its axis limits; \code{main}, \code{cex.axis},
    \code{las}, and related options customize its title and axes. Entries
    \code{pch}, \code{col}, and \code{cex} customize the middle quantile
    points. The plotted coordinates and initial empty-panel \code{type}
    are controlled internally. These arguments apply whether or not the
    CRPS inset is displayed.}
  \item{inset.args}{Named list of graphical arguments for the CRPS inset,
    passed to \code{\link[graphics]{plot.default}}. Entries \code{xlim}
    and \code{ylim} set the inset limits; \code{col}, \code{lty}, and
    \code{lwd} control the forest curve. Axis options such as
    \code{cex.axis}, \code{las}, and \code{tck} also apply to the inset.
    Use \code{axes = FALSE} to suppress both inset axes, or
    \code{yaxt = "n"} to suppress its right-hand axis. The default
    \code{list()} uses the inset's own settings, independently of
    \code{...}. Ignored when \code{crps = FALSE}. Must be supplied by its
    full argument name.}
  \item{quantreg.tau}{Optional numeric vector of levels strictly between
    zero and one.  Adds a lower-right annotation of mean pinball losses
    at these levels for the selected response and subset, using
    \code{get.pinball.error()}.  The default \code{NULL} adds no pinball
    annotation, even when a session reporting option is set.  It does
    not change \code{prbL}, \code{prbU}, or the CRPS inset.  Must be
    supplied by its full argument name.}
  \item{crps.null}{Logical. If \code{TRUE}, add the null reference when
    \code{crps = TRUE}. The reference is drawn as a dashed black curve,
    with a legend identifying \code{Forest} and \code{Null}. Set to
    \code{FALSE} to display only the forest curve. Must be supplied by its
    full argument name.}
}
\details{
  The main panel places observed response values on the horizontal axis
  and the requested conditional quantiles on the vertical axis. A point
  marks the middle quantile, and a vertical segment with endpoint marks
  spans the lower and upper quantiles. The dashed diagonal is the identity
  line. Horizontal jitter reduces overplotting; it does not change the
  response values used for CRPS. Rows with unavailable responses or
  requested quantiles are omitted from the main panel with a warning.

  An explicit \code{quantreg.tau} adds pinball losses to the main panel.
  Each loss is averaged over the requested \code{subset}, excluding
  unavailable responses and quantiles separately at that level.  It
  uses the original response values, not their jittered display
  positions.  These are the same values returned by
  \code{get.pinball.error(x, tau = quantreg.tau, subset = subset,
  m.target = m.target)} for the selected response.  The reporting levels
  are independent of the lower and upper quantiles shown in the plot.
  The session option \code{quantreg.tau} sets defaults for printing and
  loss extraction; plotting keeps this annotation off unless levels
  are supplied explicitly.

  The inset shows the average squared CDF error, integrated by the
  trapezoidal rule from the first response-grid value to each successive
  grid value and divided by the width of that interval. Its horizontal
  axis is the response threshold and its right-hand vertical axis is
  standardized CRPS. Lower values indicate less error over the
  corresponding integration interval. This is the finite-grid curve
  computed by \code{get.quantile.crps()}, with the same dependence on the
  reporting grid. The first value is unavailable because its integration
  interval has zero width. An entirely unavailable curve is omitted with
  a warning.

  The null reference uses the empirical CDF of all finite training
  responses for the selected outcome,
  \deqn{F_0(t) = \frac{1}{n_0}\sum_{j=1}^{n_0} I(Y_j^{\mathrm{train}}\leq t),}
  where \eqn{n_0} is the number of finite training responses. Repeated
  response values retain their frequencies. This same predictor-free CDF
  is used for every evaluated observation. Training responses are read
  from \code{x$forest$yvar}, or from \code{x$yvar} for a grow object
  without saved training responses in its forest. Test responses are
  used only for evaluation, not to construct the reference. Distinct
  reporting-grid values are not a replacement for the training sample.

  Both curves use the same response grid, subset, integration rule, and
  standardization. At each threshold, both scores exclude the same rows
  with an unavailable observed response or forest CDF prediction. The
  reference distribution itself is not restricted by \code{subset}.
  In training plots, the full-training empirical reference includes the
  evaluated training response; it is an in-sample benchmark, not an OOB
  or leave-one-out null estimate. If training responses are unavailable,
  the null curve is omitted with a warning.

  Graphical arguments in \code{...} control the main panel. Customize
  the inset separately with \code{inset.args}, for example,
  \code{inset.args = list(ylim = c(0, 0.3), lwd = 1.5)}. Default inset
  limits cover both curves. Limits supplied inside \code{inset.args}
  change the displayed inset window only; they do not recompute or
  truncate either score integral. Main-panel limits also leave both score
  calculations unchanged. Inset graphics settings are restored before
  returning so additional drawing refers to the main panel, and successive
  calls can be used in a multi-panel layout.
}
\value{
  Used for its graphical side effect. Invisibly returns \code{NULL}.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\seealso{
  \code{\link{quantreg}}, \code{get.quantile()},
  \code{get.quantile.crps()}, \code{\link[graphics]{plot.default}}
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
