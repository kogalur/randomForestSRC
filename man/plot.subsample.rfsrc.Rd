\name{plot.subsample.rfsrc}
\alias{plot.subsample.rfsrc}
\alias{plot.subsample}
\title{Plot Subsampling Confidence Intervals for Variable Importance}
\description{
  Display confidence intervals for variable importance (VIMP) from a
  \code{subsample} object. Select predictors, an outcome, and an
  interval method using the saved replicate estimates. Prediction-error
  and joint-VIMP rows can also be displayed when requested during
  subsampling.
}
\usage{\method{plot.subsample}{rfsrc}(x, alpha = .01, xvar.names,
 standardize = TRUE, normal = TRUE, jknife = FALSE, target, m.target = NULL,
 pmax = 75, main = "", sorted = TRUE, show.plots = TRUE, ...)}
\arguments{
  \item{x}{An object returned by \code{\link{subsample.rfsrc}}.}
  \item{alpha}{Significance level for the displayed intervals, whose
    nominal confidence level is \eqn{1-\alpha}. The plotting default
    \code{.01} gives 99 percent intervals. Use \code{alpha = .05} to
    match the default of \code{extract.subsample} and the print method.}
  \item{xvar.names}{Names of predictor rows to display. If omitted,
    all available rows are considered, subject to \code{pmax}.
    This selects stored results and does not regrow the forest.}
  \item{standardize}{For regression outcomes, divide VIMP by the
    variance of the selected response in the full training data.
    The same variance is used for all replicate estimates and any
    error row. Other families are unchanged. Set \code{FALSE} to
    display the original error scale.}
  \item{normal}{Use normal-approximation intervals when \code{TRUE}
    and nonparametric intervals when \code{FALSE}. For double-bootstrap
    objects, these are bootstrap normal and percentile intervals.}
  \item{jknife}{Select the delete-\eqn{d} jackknife standard error for
    normal intervals. Use \code{normal = TRUE, jknife = TRUE} for this
    display. The default normal display uses the subsampling standard
    error. Ignored when \code{normal = FALSE} and for double-bootstrap
    objects.}
  \item{target}{For classification, an integer or class label selecting
    class-specific VIMP; the default is overall VIMP. For competing
    risks, an integer from 1 to \code{J} selecting an event type;
    the default is the first event. This selects a statistic within
    the outcome identified by \code{m.target}.}
  \item{m.target}{Name of one response in a multivariate or mixed-outcome
    forest. If omitted, a default response is selected. This selects
    an outcome in the saved results; it is separate from choosing
    a class or event with \code{target}.}
  \item{pmax}{Maximum number of rows retained for display, selected
    by VIMP.}
  \item{main}{Main plot title.}
  \item{sorted}{Order the displayed rows by importance.}
  \item{show.plots}{Draw the plot. Set \code{FALSE} to obtain its
    plot-data return value without opening or updating a plot.}
  \item{...}{Graphical arguments for the interval display, including
    \code{xlim}, \code{ylim}, \code{xlab}, \code{ylab}, \code{boxfill},
    \code{border}, \code{whisklty}, \code{whisklwd}, and axis controls.
    A named \code{pars} list is also accepted; direct arguments take
    precedence. See Details.}
}
\details{
  \subsection{Choosing an interval display}{
    The default uses normal intervals and the subsampling standard
    error. Use \code{normal = TRUE, jknife = TRUE} for normal intervals
    with the jackknife standard error, or \code{normal = FALSE} for
    nonparametric intervals. All use the estimates already stored
    in \code{x}; changing \code{alpha} or the display method does not
    perform another round of subsampling.

    The subsampling standard error measures replicate dispersion
    around the replicate mean. The jackknife calculation measures
    dispersion around the full-data estimate, incorporating the
    displacement between those two centers. The jackknife standard
    error is not necessarily larger for a finite replicate set.
    The nonparametric interval reverses empirical quantiles of the
    centered subsampling roots. See \code{\link{subsample.rfsrc}}
    for the calculations.

    The display summarizes uncertainty in VIMP, not the distribution
    of the observed predictor values. Its intervals are calculated
    separately for each displayed statistic. A positive lower VIMP
    endpoint identifies a positive-importance result under that
    interval procedure; the display does not apply a multiple-testing
    correction. Prediction-error rows summarize the error itself.
  }
  \subsection{Selecting predictors and outcomes}{
    \code{xvar.names} selects rows and \code{m.target} selects a response.
    For a classification response or competing-risk outcome,
    \code{target} additionally selects its class-specific or
    event-specific statistic. These selections are made after
    resampling, so one saved object can be used for several displays.

    Set \code{alpha} explicitly when comparing results with
    \code{print} or \code{extract.subsample}; their default is .05,
    whereas the plotting default is .01.
  }
  \subsection{Graphical customization}{
    The default display is horizontal. \code{xlim} controls its
    horizontal VIMP/error axis and \code{ylim} controls its vertical
    row positions. These arguments always refer to the displayed axes;
    the wrapper handles the internal limit reversal used by \code{bxp}
    for horizontal boxes. \code{horizontal = FALSE} gives a vertical
    display, with VIMP/error on the vertical axis.

    Use \code{boxfill} (or \code{col}) for the box fill,
    \code{border} for its border, and \code{whisklty} and
    \code{whisklwd} for interval whiskers. By default, boxes are red
    when their lower endpoint is positive and blue otherwise. User
    styles are retained when the intervals are redrawn above the guides.
    \code{outline = FALSE} remains the default. The whiskers are
    confidence limits and the boxes show the inner interval summaries,
    rather than quartiles of the observed predictor values.

    \code{cex.axis}, \code{col.axis}, and \code{las} control axis
    labels. Set \code{xaxt = "n"} or \code{yaxt = "n"} to suppress
    one axis, or \code{axes = FALSE} to suppress both. A scalar
    \code{ylab} is a vertical-axis title. \code{names} supplies one
    row label per selected statistic before sorting and trimming; labels
    follow the same ordering as the statistics. \code{at} supplies
    positions for the final displayed intervals.

    \code{show.plots = FALSE} returns the plot data invisibly without
    opening a graphics device. \code{extract.subsample(x, raw = TRUE)}
    returns the interval matrices and replicate estimates directly.
  }
}
\value{
  Invisibly returns a boxplot-summary list. Its \code{stats} component
  is the selected five-row confidence-interval matrix, in displayed
  order, and \code{names} contains the row labels. The remaining
  components are the boxplot scaffold, not additional subsampling
  confidence intervals. The same list is returned with
  \code{show.plots = FALSE}.
}
\author{Hemant Ishwaran and Udaya B. Kogalur}
\references{
  Ishwaran H. and Lu M. (2019). Standard errors and confidence
  intervals for variable importance in random forest regression,
  classification, and survival. \emph{Statistics in Medicine}, 38,
  558-582.

  Politis, D.N. and Romano, J.P. (1994). Large sample confidence
  regions based on subsamples under minimal assumptions. \emph{The
  Annals of Statistics}, 22(4):2031-2050.

  Shao, J. and Wu, C.J. (1989). A general theory for jackknife variance
  estimation. \emph{The Annals of Statistics}, 17(3):1176-1197.
}
\seealso{
  \code{\link{subsample.rfsrc}}, \code{\link{extract.subsample}},
  \code{\link[graphics]{bxp}}
}
\examples{
\donttest{
## Small settings are for illustration; increase B for final inference.
set.seed(19)
dta <- na.omit(airquality)
o <- rfsrc(Ozone ~ ., data = dta, ntree = 100,
            importance = "permute", block.size = 1)
smp <- subsample(o, B = 25, verbose = FALSE)

## Three interval displays, all at the same confidence level.
plot.subsample(smp, alpha = .05, main = "Subsampling normal intervals")
plot.subsample(smp, alpha = .05, jknife = TRUE,
                 main = "Jackknife normal intervals")
plot.subsample(smp, alpha = .05, normal = FALSE,
                 main = "Nonparametric intervals")

## Restrict the display and customize the main axes and whiskers.
plot.subsample(smp, alpha = .05,
                 xvar.names = c("Solar.R", "Wind", "Temp"),
                 xlim = c(-.1, .8), las = 2, cex.axis = .75,
                 whisklty = 1, whisklwd = 1.5)
plot.data <- plot.subsample(smp, alpha = .05, show.plots = FALSE)
print(plot.data)

## Multivariate regression: one subsample bank, two response displays.
mv <- rfsrc(cbind(Ozone, Temp) ~ ., data = dta, ntree = 100,
             importance = "permute", block.size = 1)
mv.smp <- subsample(mv, B = 25, verbose = FALSE)
plot.subsample(mv.smp, m.target = "Ozone", alpha = .05,
                 main = "Ozone")
plot.subsample(mv.smp, m.target = "Temp", alpha = .05,
                 main = "Temperature")
print(extract.subsample(mv.smp, m.target = "Temp", alpha = .05)$var.sel.Z)
}
}
\keyword{plot}
