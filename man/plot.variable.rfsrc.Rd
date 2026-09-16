\name{plot.variable.rfsrc}
\alias{plot.variable.rfsrc}
\alias{plot.variable}
\title{Plot Marginal and Partial Dependence of Predictors}
\description{
  Display ensemble predictions against predictors, either at their observed
  values (marginal plots) or after fixing one predictor at a sequence of values
  and averaging over the other predictors (partial plots). Plot data can be
  returned without drawing and reused for custom graphics.
}
\usage{\method{plot.variable}{rfsrc}(x, xvar.names, target,
  m.target = NULL, time, surv.type = c("mort", "rel.freq",
  "surv", "years.lost", "cif", "chf"), class.type =
  c("prob", "bayes"), partial = FALSE, oob = TRUE,
  show.plots = TRUE, plots.per.page = 4, granule = 5, sorted = TRUE,
  nvar, npts = 25, smooth.lines = FALSE, subset, ...)
}
\arguments{
  \item{x}{An object of class \code{(rfsrc, grow)},
    \code{(rfsrc, synthetic)}, or a saved \code{(rfsrc, plot.variable)}
    object. New calculations require the training predictors and a usable
    saved forest. Anonymous and unsupervised forests are not supported.}
  \item{xvar.names}{Character vector naming predictors to include. The
    default is all available predictors. When replotting, select among
    predictors already present in the saved plotting object.}
  \item{target}{For classification, a class label or its position in the
    response levels; the first class is used by default. For competing
    risks, the event index, with 1 selecting the first event.}
  \item{m.target}{Single response name for a multivariate or mixed fit.
    The default selects the first available regression response, or the
    first available classification response when there is no regression
    response. For a factor response, \code{target} selects its class.}
  \item{time}{Single evaluation time for time-dependent survival summaries.
    The default is the median of the forest's stored \code{time.interest}
    grid. Supply a finite scalar no earlier than the first stored time.
    Between stored times, the last time not exceeding the request is used;
    requests beyond the grid use its last point. Ignored for mortality,
    relative mortality, and years lost.}
  \item{surv.type}{Survival summary to display. For right-censored data,
    select \code{"mort"}, \code{"rel.freq"}, \code{"surv"}, or \code{"chf"}. For competing risks, select \code{"years.lost"},
    \code{"cif"}, or \code{"chf"}. Defaults are mortality and years lost,
    respectively. See Details.}
  \item{class.type}{For a factor response, \code{"prob"} displays the
    probability of the selected class; \code{"bayes"} selects the
    class-assignment display. In partial mode, \code{"bayes"} retains the
    numeric training class-code display. Use \code{"prob"} for interpretable
    partial class probabilities; averaging unordered class codes does not
    estimate a class probability.}
  \item{partial}{Logical. \code{FALSE} produces marginal plots;
    \code{TRUE} computes partial dependence summaries.}
  \item{oob}{Logical. Request OOB predictions by default; set to
    \code{FALSE} for the non-OOB ensemble. OOB predictions use trees that
    excluded each observation during training.}
  \item{show.plots}{Logical. Set to \code{FALSE} to compute and return
    plotting data without opening or changing a graphics device.}
  \item{plots.per.page}{Positive integer controlling the panel layout.
    The current layout uses up to this many columns and up to this many
    rows. Thus a value of 4 can place more than four panels on a page.}
  \item{granule}{Numeric threshold for displaying a numeric predictor as
    a boxplot. Factors always use boxplots. For marginal plots the
    comparison uses the number of distinct observed values; for partial
    plots it uses the number of evaluation-grid values.}
  \item{sorted}{Logical. Order predictors by decreasing available variable
    importance before applying \code{nvar}. Unavailable importance values
    follow finite ones. Ties retain predictor order. When no importance is stored,
    retain the predictor order.}
  \item{nvar}{Maximum number of predictors selected for a new calculation.
    The default is all predictors selected by \code{xvar.names}.}
  \item{npts}{Maximum number of evaluation points for a numeric predictor
    in partial mode. Points are selected from sorted distinct values,
    approximately evenly in their index order. Factors are not thinned.}
  \item{smooth.lines}{Logical. Apply \code{lowess} to continuous partial
    curves and, when available, their display bands. The stored partial
    means are unchanged. Continuous marginal plots use \code{lowess}
    regardless of this option.}
  \item{subset}{Numeric row indices or a logical selection relative to the
    fitted object's \code{xvar}. The default uses all stored training rows.
    Use \code{x$xvar}, not row numbers from an original data set that may
    have lost rows during fitting. In partial mode the selection determines
    both the displayed predictor grid and the requested averaging rows.
    Logical selections must have one entry per stored row. Missing and
    out-of-range numeric indices are omitted; fractional indices are rejected.}
  \item{...}{Named graphical arguments for the displayed panels.
    Use \code{main}, \code{xlab}, \code{ylab}, \code{xlim}, and \code{ylim}
    for titles, labels, and limits. \code{col}, \code{lty}, \code{lwd},
    \code{pch}, and \code{cex} customize continuous displays. \code{type}
    controls the fitted series; the partial mean defaults to points joined
    by lines. Partial variability bands are drawn separately.
    Boxplots accept \code{notch}, \code{boxfill}, \code{whisklwd},
    \code{outline}, \code{horizontal}, \code{names}, and related
    \code{bxp} settings. Limits always refer to the displayed axes,
    including horizontal boxes. A named \code{pars} list is also accepted;
    individually supplied arguments take precedence. Graphical settings do
    not change the predictions, means, or standard errors.}

}
\details{
  \subsection{Marginal and partial plots}{
    Marginal plots relate the ensemble's predictions to the observed value
    of a predictor. They preserve the observed combinations of predictors.
    Continuous predictors receive a smoothed curve; factors and numeric
    predictors with few distinct values receive grouped boxplots.

    Partial plots fix the selected predictor at each grid value for every
    observation, retaining the other predictor coordinates, and summarize
    those predictions. For a grid value \eqn{z}, write the observation-level
    partial prediction as \eqn{P_i(z)=\widehat f(z,X_{i,-j})}. With complete
    predictions the curve is \eqn{n^{-1}\sum_i P_i(z)}; missing predictions
    are omitted from the mean. When OOB output is requested, \eqn{P_i(z)}
    uses the OOB trees for row \eqn{i}.

    Partial mode reuses the fitted trees but performs additional prediction
    work at each grid value. Restrict \code{xvar.names} or \code{nvar} and
    use a moderate \code{npts} to control the computation. For explicit
    grids or several fixed predictors, use \code{\link{partial.rfsrc}}.
  }
  \subsection{Response selection}{
    Regression displays the predicted response. Classification displays
    the selected class probability by default. The \code{m.target} argument
    chooses the response for a multivariate or mixed forest; \code{target}
    then chooses a class when that response is a factor.

    Right-censored survival displays mortality, relative-frequency
    mortality, survival probability, or cumulative hazard. Competing risks displays expected
    years lost, cumulative incidence, or cause-specific cumulative hazard
    for the selected event. Survival probability, cumulative incidence,
    and cumulative hazard depend on \code{time}; mortality and years lost
    are scalar summaries and do not use the requested time.

    The existing marginal display multiplies right-censored survival and
    cumulative hazard values by 100; partial displays retain their raw
    scales. Relative mortality also retains its existing mode-specific
    scaling: marginal mortality is divided by the maximum of the selected
    row count and the available full-vector mortality values; partial
    mortality is divided by the available row count at each grid value.
  }
  \subsection{Subsets and missing values}{
    The plotting calculation starts from the stored training predictors
    and overlays saved imputed predictor values when available. It does
    not grow a replacement forest. Define subsets relative to the stored
    row order. In partial mode, grid construction uses nonmissing values
    of the selected predictor within that subset; the same subset is
    passed to the partial-prediction extraction step.

    A subset affects the population being summarized, not which observations
    trained the forest. Likewise, OOB selection changes the contributing
    trees for an observation, not the definition of the plotting subset.
  }
  \subsection{Partial means and standard errors}{
    A single calculation is used for regression predictions, class
    probabilities, and scalar survival or competing-risk predictions.
    At grid value \eqn{z}, let \eqn{I_z} contain the selected rows with
    finite partial predictions, and let \eqn{n_z} be their number. The mean
    and empirical row-averaging standard error are
    \deqn{
      \bar P(z)=\frac{1}{n_z}\sum_{i\in I_z}P_i(z),\qquad
      \mathrm{se}(z)=\left\{
        \frac{\sum_{i\in I_z}[P_i(z)-\bar P(z)]^2}{n_z(n_z-1)}
      \right\}^{1/2}.
    }
    Thus \code{se} is the sample standard deviation divided by
    \eqn{\sqrt{n_z}}. Classification uses the observed variation among
    predicted probabilities, rather than a Bernoulli variance calculated
    from their mean. Both numerator and denominator use the same available
    rows. With no finite predictions the mean is \code{NA}; with fewer
    than two, the standard deviation and standard error are \code{NA}.
    Identical finite predictions have zero standard error when at least
    two are present.

    This standard error describes reference-row averaging with a fixed
    prediction function under independent reference sampling. For the
    training/OOB plots here, it is a descriptive row-averaging scale:
    it does not account for fitting the forest or for dependence between
    overlapping OOB fits. The continuous display uses pointwise bands
    \eqn{\bar P(z)\pm2\mathrm{se}(z)}, not refit-based confidence intervals
    or simultaneous confidence bands. The numeric class-code display
    uses this calculation on its code scale, not a probability scale.
  }
  \subsection{Discrete partial displays}{
    Discrete displays retain boxplots of the transformed values
    \eqn{\bar P(z)+[P_i(z)-\bar P(z)]/\sqrt{n_z}}. Their sample standard
    deviation is exactly \code{se} when \eqn{n_z>1}, so continuous and
    discrete displays use the same row-averaging scale. The boxes and
    whiskers remain ordinary boxplot summaries of these transformed
    values; their endpoints are not the mean plus or minus two standard
    errors. Notches likewise concern these display values, not uncertainty
    from refitting the forest. A single available value remains at its
    mean, with an unavailable standard error.

    The returned \code{plotthis} tables retain these display values.
    \code{partial.summary} gives the means, standard deviations, standard
    errors, and available counts directly for every predictor, including
    factors. For unscaled observation-level predictions, use
    \code{partial()} and \code{get.partial.plot.data(..., granule = TRUE)}.
  }
  \subsection{Returning and reusing plot data}{
    Assign the result of \code{plot.variable(..., show.plots = FALSE)} to
    retain the calculations. For partial mode, \code{plotthis} contains
    one data frame per predictor with columns \code{x} and \code{yhat}.
    Continuous curves normally have one row per grid value; discrete
    displays retain the shrunken values for each selected observation.

    Call \code{plot.variable(saved)} to redraw without recalculating
    predictions. The saved object determines the marginal or partial
    mode, outcome, subset, grid, and display settings. Its
    \code{plots.per.page} and \code{smooth.lines} components can be changed
    before replotting. To change the prediction target, subset, or grid,
    create a new plotting object from the grow object.
  }
}
\value{
  Invisibly returns an object of class \code{(rfsrc, plot.variable, family)}.
  Common components include \code{family}, \code{partial}, \code{event.info},
  \code{target}, \code{ylabel}, \code{n}, \code{xvar.names}, \code{nvar},
  \code{plots.per.page}, \code{granule}, and \code{smooth.lines}.
  \item{pData}{For partial plots, a named list of predictor-specific
    components: \code{xvar.names}, \code{yhat}, \code{yhat.se}, \code{n.x},
    \code{x.uniq}, and the nonmissing predictor values \code{x}.}
  \item{plotthis}{For partial plots, named data frames with columns
    \code{x} and \code{yhat}, suitable for custom displays. These are
    display data, as described in Details.}
  \item{partial.summary}{For partial plots, named data frames with one row
    per evaluated grid value and columns \code{x}, \code{mean}, \code{sd},
    \code{se}, and \code{n.used}. Counts refer to available predictions,
    not an estimated independence-adjusted effective sample size.}
  \item{yhat, xvar}{For marginal plots, the selected predictions and
    predictor data.}
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Friedman J.H. (2001). Greedy function approximation: a gradient
  boosting machine, \emph{Ann. of Statist.}, 5:1189-1232.
  
  Ishwaran H., Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7(2):25-31.
  
  Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
  (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.
  
  Ishwaran H., Gerds T.A., Kogalur U.B., Moore R.D., Gange S.J. and Lau
  B.M. (2014). Random survival forests for competing risks.
  \emph{Biostatistics}, 15(4):757-773.
  
}
\seealso{
  \command{\link{rfsrc}},
  \command{\link{partial.rfsrc}},
  \command{\link{predict.rfsrc}}
}
\examples{
\donttest{
## ------------------------------------------------------------
##
## Regression: marginal versus partial display.
##
## ------------------------------------------------------------
set.seed(19)
dta <- na.omit(airquality)
fit <- rfsrc(Ozone ~ ., data = dta, ntree = 100, importance = TRUE)
plot.variable(fit, xvar.names = c("Wind", "Temp"), sorted = FALSE,
              plots.per.page = 2)
pv <- plot.variable(fit, xvar.names = c("Wind", "Temp"), sorted = FALSE,
                    partial = TRUE, npts = 15, plots.per.page = 2,
                    show.plots = FALSE)
print(head(pv$plotthis$Wind))
print(pv$partial.summary$Wind)
plot.variable(pv, main = "Partial dependence", lwd = 2)

## ------------------------------------------------------------
##
## Inspect the mean and its row-averaging standard error.
##
## ------------------------------------------------------------
ss <- pv$partial.summary$Wind
plot(ss$x, ss$mean, type = "b", xlab = "Wind", ylab = "Partial prediction")
lines(ss$x, ss$mean - 2 * ss$se, lty = 3)
lines(ss$x, ss$mean + 2 * ss$se, lty = 3)

## ------------------------------------------------------------
##
## Customize a stored continuous curve without another forest calculation.
##
## ------------------------------------------------------------
pv$plots.per.page <- 1
pv$smooth.lines <- TRUE
plot.variable(pv, xvar.names = "Temp", main = "Temperature")

## ------------------------------------------------------------
##
## Subset indices refer to the fitted object's stored rows.
##
## ------------------------------------------------------------
use <- fit$xvar$Solar.R < 200
ps <- plot.variable(fit, xvar.names = "Wind", partial = TRUE,
                    subset = use, npts = 15, show.plots = FALSE)
print(head(ps$plotthis$Wind))

## ------------------------------------------------------------
##
## Multivariate regression: choose a response for the plot.
##
## ------------------------------------------------------------
mv <- rfsrc(cbind(Ozone, Temp) ~ ., data = dta, ntree = 100)
mv.temp <- plot.variable(mv, m.target = "Temp", xvar.names = "Wind",
                         partial = TRUE, npts = 15, show.plots = FALSE)
print(mv.temp$partial.summary$Wind)
plot.variable(mv.temp)

## ------------------------------------------------------------
##
## Classification: choose a class label.
##
## ------------------------------------------------------------
iris.fit <- rfsrc(Species ~ ., data = iris, ntree = 100)
plot.variable(iris.fit, target = "versicolor", xvar.names = "Petal.Length",
              partial = TRUE, npts = 15)

## ------------------------------------------------------------
##
## A factor predictor receives a boxplot.
##
## ------------------------------------------------------------
car.data <- mtcars
car.data$cyl <- factor(car.data$cyl)
car.fit <- rfsrc(mpg ~ ., data = car.data, ntree = 100)
pc <- plot.variable(car.fit, xvar.names = "cyl", partial = TRUE,
                     show.plots = FALSE)
boxplot(yhat ~ x, data = pc$plotthis$cyl,
        xlab = "Cylinders", ylab = "Partial display values")

## ------------------------------------------------------------
##
## Mixed outcomes: m.target selects the response, target selects its class.
##
## ------------------------------------------------------------
car.data$am <- factor(car.data$am)
mix <- rfsrc(Multivar(mpg, am) ~ ., data = car.data, ntree = 100)
plot.variable(mix, m.target = "am", target = "1", xvar.names = "wt",
              partial = TRUE, npts = 15)

## ------------------------------------------------------------
##
## Right-censored survival at a stored evaluation time.
##
## ------------------------------------------------------------
data(veteran, package = "randomForestSRC")
vfit <- rfsrc(Surv(time, status) ~ ., data = veteran, ntree = 100)
when <- vfit$time.interest[ceiling(length(vfit$time.interest) / 2)]
plot.variable(vfit, xvar.names = "age", surv.type = "surv", time = when,
              partial = TRUE, npts = 15)

## ------------------------------------------------------------
##
## Competing risks: cumulative incidence for the second event.
##
## ------------------------------------------------------------
data(follic, package = "randomForestSRC")
crfit <- rfsrc(Surv(time, status) ~ ., data = follic, ntree = 100, nsplit = 3)
when <- crfit$time.interest[ceiling(length(crfit$time.interest) / 2)]
plot.variable(crfit, xvar.names = "age", target = 2, surv.type = "cif",
              time = when, partial = TRUE, npts = 15)
}
}
\keyword{plot}
