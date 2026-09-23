\name{plot.subsample.rfsrc}
\alias{plot.subsample.rfsrc}
\alias{plot.subsample}
\title{Plot Subsampling Confidence Intervals for Variable Importance}
\description{
  Plot confidence intervals for variable importance (VIMP) from a
  \code{subsample} analysis. Compare predictors and assess whether their
  importance is positive under the chosen interval method. Joint VIMP
  and prediction error can also be displayed when included in the
  resampling analysis.
}
\usage{\method{plot.subsample}{rfsrc}(x, alpha = .01, xvar.names,
 standardize = TRUE, normal = TRUE, jknife = FALSE, target, m.target = NULL,
 pmax = 75, main = "", sorted = TRUE, show.plots = TRUE, ...)}
\arguments{
  \item{x}{An object returned by \code{\link{subsample.rfsrc}}, including
    double-bootstrap output.}
  \item{alpha}{Significance level, strictly between zero and one.
    Intervals have nominal confidence level \eqn{1-\alpha}. The default
    \code{.01} gives 99 percent intervals. Use \code{alpha = .05} to
    match the default for extraction and printing.}
  \item{xvar.names}{Names of rows to display: predictor names, or
    \code{"joint"} and \code{"err"} when those statistics were requested.
    If omitted, all rows are considered, subject to \code{pmax}.}
  \item{standardize}{For regression, divide estimates by the selected
    response's variance in the full training data. The same divisor
    applies to all replicates, VIMP, and prediction error. Other
    families are unchanged. Set \code{FALSE} for the original scale.}
  \item{normal}{Use normal-approximation intervals when \code{TRUE}
    and nonparametric intervals when \code{FALSE}. For double-bootstrap
    objects, these are bootstrap normal and percentile intervals.}
  \item{jknife}{Use the delete-\eqn{d} jackknife standard error for
    normal intervals. The default uses the subsampling standard error.
    Ignored when \code{normal = FALSE} or for double-bootstrap objects.}
  \item{target}{For classification, use \code{0} (the default) for the
    overall statistic, or a class label or position \code{1}, \code{2},
    and so on for a class-specific statistic. For competing risks, use
    an event position from \code{1} to \code{J}, where \code{J} is the
    number of event types; the first event is the default. Ignored for
    regression and right-censored survival.}
  \item{m.target}{Name of one response in a multivariate or mixed-outcome
    forest. The default selects the first regression response with
    results, or the first classification response with results if none
    is available for regression. Use \code{target} to select a class
    or event within that response.}
  \item{pmax}{Maximum number of rows to display. Supply a positive
    integer or \code{Inf} for all rows. When this limit is exceeded,
    the rows with the largest reference estimates are displayed in
    decreasing order, even if \code{sorted = FALSE}.}
  \item{main}{Main plot title.}
  \item{sorted}{Order rows by decreasing reference estimate. For VIMP,
    this places the most important predictors first. See \code{pmax}.}
  \item{show.plots}{Draw the plot. Set \code{FALSE} to return its
    numerical summaries without drawing.}
  \item{...}{Named graphical arguments, such as \code{xlim},
    \code{ylim}, \code{xlab}, \code{ylab}, \code{boxfill},
    \code{border}, \code{whisklty}, and \code{whisklwd}. A named
    \code{pars} list is also accepted; direct arguments take precedence.
    See Details.}
}
\details{
  After fitting \code{smp <- subsample(obj)}, use
  \code{plot.subsample(smp)} to compare VIMP estimates and their
  uncertainty. Each row represents a predictor or a requested joint-VIMP
  or prediction-error statistic. Changing the confidence level, interval
  method, or selected results uses the saved replicates without fitting
  more forests.

  \subsection{Interpreting the display}{
    The default plot is horizontal, with VIMP or prediction error on
    the horizontal axis. Whiskers mark the confidence limits; the box
    spans the lower and upper quartiles of the interval summary.
    For subsampling, the center line marks the full-data estimate for
    normal intervals and the interval median for nonparametric
    intervals. For double-bootstrap output, it marks the bootstrap mean.

    Boxes are red when the lower confidence limit exceeds zero and
    blue otherwise. For VIMP, red identifies positive importance under
    the chosen interval method. An \code{err} row summarizes prediction
    error, so its color is not a variable-selection result. Intervals
    are calculated separately, without a multiple-testing adjustment.
    These boxes describe uncertainty in the estimates, not variation
    in the observed predictor values.
  }
  \subsection{Choosing the interval method}{
    Subsampling provides three displays:
    \describe{
      \item{Normal subsampling (default).}{Use
        \code{normal = TRUE, jknife = FALSE}. The standard error measures
        variation about the subsample mean.}
      \item{Normal jackknife.}{Use \code{normal = TRUE, jknife = TRUE}.
        The standard error uses deviations from the full-data estimate,
        including its difference from the subsample mean.}
      \item{Nonparametric subsampling.}{Use \code{normal = FALSE}.
        Limits use empirical quantiles of centered and scaled
        subsample estimates.}
    }
    The jackknife standard error need not be larger than the subsampling
    standard error. See \code{\link{subsample.rfsrc}} for the formulas
    and finite-replicate scaling.

    For double-bootstrap objects, \code{normal = TRUE} uses normal
    intervals centered at the bootstrap mean; \code{normal = FALSE}
    uses percentile intervals. \code{jknife} has no effect.

    Use the same \code{alpha} when comparing plots with printed or
    extracted results. Plotting defaults to \code{.01}; printing and
    extraction default to \code{.05}.
  }
  \subsection{Selecting predictors and outcomes}{
    Use \code{xvar.names} to choose rows, \code{m.target} to choose a
    response, and \code{target} to choose its class or event statistic.
    For example, \code{m.target = "Ozone"} selects results for the
    response named \code{Ozone}, not a predictor to display.

    Competing-risk indexing differs from extraction: the first event
    uses \code{target = 1} here but \code{target = 0} in
    \code{extract.subsample()}.
  }
  \subsection{Customizing the plot}{
    \code{xlim} and \code{ylim} control the horizontal and vertical
    axes as displayed. Use \code{horizontal = FALSE} to place VIMP or
    prediction error on the vertical axis.

    Set \code{boxfill} (or \code{col}) for box colors, \code{border}
    for borders, and \code{whisklty} and \code{whisklwd} for the
    confidence-limit lines. Outliers are hidden by default
    (\code{outline = FALSE}).

    \code{cex.axis}, \code{col.axis}, and \code{las} control axis
    labels. Use \code{xaxt = "n"} or \code{yaxt = "n"} to suppress one
    axis, or \code{axes = FALSE} for both. A single \code{ylab} string
    sets the vertical-axis title. For custom row labels, supply
    \code{names} with one label per selected row before sorting and
    limiting the display to \code{pmax} rows. \code{at} sets positions
    for the final displayed intervals.
  }
  \subsection{Obtaining numerical results}{
    \code{plot.subsample(smp, show.plots = FALSE)} returns the selected
    interval matrix and labels without drawing.
    \code{extract.subsample(smp, raw = TRUE)} returns all interval
    methods and the replicate estimates. Both use the requested
    \code{standardize} setting.
  }
}
\value{
  Invisibly returns a list with the following plotting components:
  \item{stats}{Five-row interval matrix, with one column per displayed
    statistic. Rows contain the lower confidence limit, lower quartile,
    center, upper quartile, and upper confidence limit.}
  \item{names}{Labels corresponding to the columns of \code{stats}.}
  Other components support the boxplot drawing and are not additional
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
