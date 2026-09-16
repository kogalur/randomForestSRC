\name{partial.rfsrc}
\alias{partial.rfsrc}
\alias{partial}
\alias{get.partial.plot.data}
\title{Compute Partial Dependence Values}
\description{
  Compute observation-level predictions after setting a predictor to specified
  values. The saved forest is reused. Results can be averaged for partial
  dependence curves or retained for custom plots. Regression, classification,
  multivariate and mixed outcomes, right-censored survival, and competing
  risks have family-specific outputs.
}
\usage{
partial.rfsrc(object, oob = TRUE,
  partial.type = NULL, partial.xvar = NULL, partial.values = NULL,
  partial.xvar2 = NULL, partial.values2 = NULL,
  partial.time = NULL, get.tree = NULL, seed = NULL, do.trace = FALSE, ...)

get.partial.plot.data(o, target, m.target = NULL, granule = FALSE)
}
\arguments{
  \item{object}{An object of class \code{(rfsrc, grow)} with a saved forest,
    or an object of class \code{(rfsrc, forest)} retaining its training data.
    Anonymous forests are not supported.}
  \item{oob}{Logical. Request out-of-bag (OOB) partial predictions by default.
    Each observation is predicted using trees for which it was OOB in the
    original training sample. Set to \code{FALSE} to request the non-OOB
    ensemble. This selects predictions within the saved forest; it does not
    refit the forest or select a different set of observations for averaging.}
  \item{partial.type}{For right-censored survival, one of \code{"mort"},
    \code{"rel.freq"}, \code{"surv"}, or \code{"chf"}. For competing risks,
    one of \code{"years.lost"}, \code{"cif"}, or \code{"chf"}. Specify a
    single type for these families. For regression and classification this
    argument is not used and may be \code{NULL}.}
  \item{partial.xvar}{Name of the single predictor to vary.}
  \item{partial.values}{Nonempty numeric vector of values for
    \code{partial.xvar}. Values are sorted and duplicates removed before
    calculation. For a factor, supply its integer codes in the original
    training level order, not its labels.}
  \item{partial.xvar2}{Optional character vector of additional predictors
    to fix while varying \code{partial.xvar}.}
  \item{partial.values2}{One numeric value for each predictor in
    \code{partial.xvar2}, in the same order. Factors require their original
    training integer codes. These are fixed values, not additional grids.}
  \item{partial.time}{Numeric vector of evaluation times for \code{"surv"},
    \code{"cif"}, or \code{"chf"}. Supply at least one time; using points
    from the forest's \code{time.interest} gives exact grid matches.
    Between grid points, the last stored time not exceeding the requested
    time is used. Times must be finite and must not precede the first
    stored time. A request beyond the last stored time uses that last time.
    Missing, empty, or invalid requests are rejected before native prediction.
    Mortality and years-lost outputs have no requested-time dimension and
    ignore this argument.}
  \item{get.tree}{Optional vector of tree indices. By default all trees
    are considered, with OOB eligibility applied when \code{oob = TRUE}.}
  \item{seed}{Negative integer specifying the random-number seed.}
  \item{do.trace}{Number of seconds between progress updates; \code{FALSE}
    suppresses updates.}
  \item{...}{Additional internal prediction options. Plotting arguments
    belong in the subsequent plotting call, not in \code{partial()}.}
  \item{o}{An object returned by \code{partial()}.}
  \item{target}{For classification extraction, a class label or its position
    within the selected response's levels; the first class is used when
    omitted. This applies to univariate, multivariate, and mixed fits. For competing risks, an event
    code from \code{o$event.info$event.type}; omission selects the first
    event. Ignored for regression.}
  \item{m.target}{For multivariate or mixed \code{partial} objects, the
    response name to extract. The default is the first response in
    \code{o$yvar.names}. All response arrays are already in \code{o}.}
  \item{granule}{Logical used by \code{get.partial.plot.data()}. The default
    \code{FALSE} averages the observation-level partial predictions at each
    grid value. With \code{TRUE}, retain them all and repeat each grid value
    once per observation. This is different from the numeric \code{granule}
    argument of \code{plot.variable()}.}
}
\details{
  \subsection{Partial predictions and averaging}{
    For a primary predictor \eqn{X_j}, the calculation replaces its value by
    \eqn{z} for every training observation while leaving the other predictor
    values at that observation's values. Write the resulting prediction as
    \eqn{P_i(z) = \widehat f(z, X_{i,-j})}. The default summary returned by
    \code{get.partial.plot.data()} is the average of the available
    \eqn{P_i(z)} values. With OOB prediction, the prediction function for
    row \eqn{i} uses that row's OOB trees.

    \code{partial()} returns the observation-level arrays and does not
    average them. The extraction helper replaces infinite values by
    \code{NA} and omits missing values when averaging. A grid value with
    no available predictions has an unavailable mean.
  }
  \subsection{Predictor grids and factors}{
    The returned \code{partial.values} gives the actual sorted, unique
    grid. The input order and duplicated values are not retained.
    Reducing the grid reduces the number of partial predictions and the
    size of the result.

    Factor codes refer to the level order saved during training. For
    example, if the original levels are \code{c("low", "middle", "high")},
    the code for \code{"high"} is 3 even when a plotting subset contains
    only \code{"middle"} and \code{"high"}. Use
    \code{match(labels, levels(training.factor))} to obtain these codes.
    The extraction helper returns numeric grid values, so factor labels
    can be added in the custom plotting call.
  }
  \subsection{Fixing additional predictors}{
    \code{partial.xvar2} and \code{partial.values2} fix one or more
    additional coordinates for every observation. For example, vary wind
    over a grid while holding temperature at a specified value. To make a
    two-dimensional display, repeat the call at several fixed temperature
    values. A single call does not expand \code{partial.values2} into a
    Cartesian product of grids.
  }
  \subsection{Responses and survival summaries}{
    Regression outputs contain predicted responses. Classification outputs
    contain an \code{"all"} column followed by columns named for the class
    levels. Multivariate and mixed fits return named response components;
    \code{m.target} in the extraction helper selects among these components.
    The class selector \code{target} is then resolved using that response's
    class labels, not response names. A single observation or grid point
    retains its observation-by-grid meaning during extraction.

    For right-censored survival, \code{"mort"} returns mortality,
    \code{"surv"} returns survival probabilities, and \code{"chf"} returns
    cumulative hazards. At this low-level interface, \code{"rel.freq"}
    requests the same raw output as \code{"mort"}; neither \code{partial()}
    nor \code{get.partial.plot.data()} normalizes this output.

    Competing-risk outputs are expected years lost, cumulative incidence,
    or cause-specific cumulative hazard, selected by \code{partial.type}.
    The time-dependent outputs are calculated on the forest's stored time
    grid and then selected in R for \code{partial.time}. The returned
    \code{partial.time} retains the requested times.
  }
  \subsection{Custom plots}{
    With \code{granule = FALSE}, the extraction helper returns one mean
    for each primary grid value. Time-dependent survival summaries have
    one column per requested time. With \code{granule = TRUE}, values are
    ordered by primary grid value and then observation; the result contains
    raw predictions, without the display shrinkage used for boxplots by
    \code{plot.variable(partial = TRUE)}.

    The automatic plotting wrapper also returns \code{partial.summary},
    containing the partial mean, empirical row-averaging standard error,
    and available row count for each grid value. See
    \code{\link{plot.variable.rfsrc}} for that scale's interpretation.

    Use \code{plot.variable()} for an automatically selected predictor
    grid and built-in displays. Use \code{partial()} when explicit grids,
    several evaluation times, additional fixed predictors, or raw
    observation-level results are needed.
  }
}
\value{
  \code{partial()} returns an object of class \code{(rfsrc, partial, family)}.
  It contains \code{call}, \code{family}, \code{partial.values},
  \code{partial.time}, \code{yvar.names}, and \code{event.info}, together with
  the applicable prediction components below. Let \eqn{n} be the number of
  training rows, \eqn{K} the number of sorted unique primary values,
  \eqn{T} the number of requested times, and \eqn{J} the number of events.
  \item{regrOutput}{Named list of regression-response matrices, each
    \eqn{n \times K}.}
  \item{classOutput}{Named list of classification-response arrays, each
    \eqn{n \times (1 + C) \times K}, where \eqn{C} is that response's number
    of classes. The second dimension is named \code{all} followed by the
    class labels.}
  \item{survOutput}{For right-censored mortality, an \eqn{n \times K}
    matrix. For survival or cumulative hazard, an \eqn{n \times T \times K}
    array. For competing-risk years lost, an \eqn{n \times J \times K}
    array. For competing-risk cumulative incidence or cumulative hazard,
    an \eqn{n \times T \times J \times K} array.}

  \code{get.partial.plot.data()} returns a list with \code{x}, \code{yhat},
  and \code{partial.time}. Normally \code{x} is the primary grid and
  \code{yhat} is a vector of means or a grid-by-time matrix. With
  \code{granule = TRUE}, \code{x} is repeated across observations and
  \code{yhat} contains the corresponding raw values, with one column per
  time for time-dependent summaries.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran H., Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7(2):25-31.
  
  Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
  (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.
}
\seealso{
  \command{\link{plot.variable.rfsrc}}
}
\examples{
\donttest{
## ------------------------------------------------------------
##
## Regression: explicit grid and raw observation-level output.
##
## ------------------------------------------------------------
set.seed(19)
dta <- na.omit(airquality)
fit <- rfsrc(Ozone ~ ., data = dta, ntree = 100)
wind <- unique(as.numeric(quantile(fit$xvar$Wind, seq(0, 1, length.out = 11))))
p <- partial(fit, partial.xvar = "Wind", partial.values = wind)
pd <- get.partial.plot.data(p)
print(head(data.frame(Wind = pd$x, partial = pd$yhat)))
plot(pd$x, pd$yhat, type = "b", xlab = "Wind", ylab = "Partial prediction")
raw <- get.partial.plot.data(p, granule = TRUE)
boxplot(raw$yhat ~ raw$x, xlab = "Wind", ylab = "Individual partial predictions")

## ------------------------------------------------------------
##
## Vary wind while fixing temperature at three values.
##
## ------------------------------------------------------------
temperature <- as.numeric(quantile(fit$xvar$Temp, c(.25, .5, .75)))
slices <- lapply(temperature, function(tt) {
  pp <- partial(fit, partial.xvar = "Wind", partial.values = wind,
                partial.xvar2 = "Temp", partial.values2 = tt)
  get.partial.plot.data(pp)$yhat
})
matplot(wind, do.call(cbind, slices), type = "l", lty = 1,
        xlab = "Wind", ylab = "Partial prediction")
legend("topright", legend = paste("Temp =", temperature),
       col = seq_along(temperature), lty = 1)

## ------------------------------------------------------------
##
## regression: partial effects for three variables simultaneously
## (can be slow, so modify accordingly)
##
## ------------------------------------------------------------

n <- 200
x <- matrix(rnorm(n * 3), ncol = 3)
y <- x[, 1] + x[, 1] * x[, 2] + x[, 1] * x[, 2] * x[, 3]
o <- rfsrc(y ~ ., data = data.frame(y = y, x))

## define target x values
x1 <- seq(-3, 3, length = 40)
x2 <- x3 <- seq(-3, 3, length = 10)

## extract second order partial effects
pdta <- do.call(rbind,
          lapply(x3, function(x3v) {
            cat("outer loop x3 = ", x3v, "\n")
            do.call(rbind,lapply(x2, function(x2v) {
              o <- partial(o,
                      partial.xvar = "X1",
                      partial.values = x1,
                      partial.xvar2 = c("X2", "X3"),
                      partial.values2 = c(x2v, x3v))
              cbind(x1, x2v, x3v, get.partial.plot.data(o)$yhat)
            }))
          }))
pdta <- data.frame(pdta)
colnames(pdta) <- c("x1", "x2", "x3", "effectSize")

## coplot of partial effects
coplot(effectSize ~ x1|x2*x3, pdta, pch = 16, overlap = 0)


## ------------------------------------------------------------
##
## Classification: extract a named class from a univariate response.
##
## ------------------------------------------------------------
iris.fit <- rfsrc(Species ~ ., data = iris, ntree = 100)
sp <- partial(iris.fit, partial.xvar = "Sepal.Length",
              partial.values = seq(4.5, 7.5, length.out = 11))
sp.data <- get.partial.plot.data(sp, target = "versicolor")
plot(sp.data$x, sp.data$yhat, type = "b", ylim = c(0, 1),
     xlab = "Sepal length", ylab = "Partial probability of versicolor")

## ------------------------------------------------------------
##
## Multivariate regression: one partial call, two extracted responses.
##
## ------------------------------------------------------------
mv <- rfsrc(cbind(Ozone, Temp) ~ ., data = dta, ntree = 100)
mp <- partial(mv, partial.xvar = "Wind", partial.values = wind)
print(names(mp$regrOutput))
mp.ozone <- get.partial.plot.data(mp, m.target = "Ozone")
mp.temp <- get.partial.plot.data(mp, m.target = "Temp")
print(head(data.frame(Wind = mp.temp$x, Temp = mp.temp$yhat)))

## ------------------------------------------------------------
##
## Mixed outcomes: extract a named class from the selected response.
##
## ------------------------------------------------------------
car.mix <- mtcars
car.mix$am <- factor(car.mix$am)
mix <- rfsrc(Multivar(mpg, am) ~ ., data = car.mix, ntree = 100)
mix.p <- partial(mix, partial.xvar = "wt", partial.values = c(2, 3, 4))
mix.am <- get.partial.plot.data(mix.p, m.target = "am", target = "1")
print(data.frame(wt = mix.am$x, probability = mix.am$yhat))

## ------------------------------------------------------------
##
## Factors: keep the training level codes when choosing labels.
##
## ------------------------------------------------------------
car.data <- mtcars
car.data$cyl <- factor(car.data$cyl)
car.fit <- rfsrc(mpg ~ ., data = car.data, ntree = 100)
lev <- levels(car.fit$xvar$cyl)
code <- match(lev, levels(car.fit$xvar$cyl))
cp <- partial(car.fit, partial.xvar = "cyl", partial.values = code)
cd <- get.partial.plot.data(cp, granule = TRUE)
boxplot(cd$yhat ~ factor(cd$x, levels = code, labels = lev),
        xlab = "Cylinders", ylab = "Individual partial predictions")

## ------------------------------------------------------------
##
## Right-censored survival: several stored evaluation times.
##
## ------------------------------------------------------------
data(veteran, package = "randomForestSRC")
vfit <- rfsrc(Surv(time, status) ~ ., data = veteran, ntree = 100)
grid <- vfit$time.interest
when <- grid[unique(round(seq(1, length(grid), length.out = 3)))]
karno <- unique(as.numeric(quantile(vfit$xvar$karno)))
vp <- partial(vfit, partial.type = "surv", partial.xvar = "karno",
              partial.values = karno, partial.time = when)
vd <- get.partial.plot.data(vp)
print(vd$partial.time)
matplot(vd$x, vd$yhat, type = "l", lty = 1,
        xlab = "Karnofsky score", ylab = "Partial survival probability")
legend("bottomright", legend = paste("Time =", vd$partial.time),
       col = seq_along(vd$partial.time), lty = 1)

## ------------------------------------------------------------
##
## Survival curves at selected quantiles of Karnofsky score.
## Reuse vfit from the preceding right-censored survival example.
##
## ------------------------------------------------------------
karno.probs <- c(.10, .50, .90)
karno.values <- quantile(vfit$xvar$karno, probs = karno.probs,
                         na.rm = TRUE)
print(karno.values)

## Set Karnofsky score to each requested value for every reference row.
## Keep the other predictors at their observed values and use OOB trees.
surv.partial <- partial(vfit,
                         oob = TRUE,
                         partial.type = "surv",
                         partial.xvar = "karno",
                         partial.values = as.numeric(karno.values),
                         partial.time = vfit$time.interest)
surv.data <- get.partial.plot.data(surv.partial)

## yhat has one row per distinct Karnofsky value and one column per time.
## Transpose it so that each plotted column is a survival curve over time.
## Each curve averages over the same reference rows.
curve.id <- seq_along(surv.data$x)
matplot(surv.data$partial.time, t(surv.data$yhat),
        type = "s", col = curve.id, lty = curve.id, lwd = 2,
        ylim = c(0, 1), xlab = "Time", ylab = "Survival probability",
        main = "Partial survival curves by Karnofsky score")

## Use the returned grid for labels: partial() removes tied quantile values.
legend("topright",
       legend = paste("Karnofsky =", format(surv.data$x, trim = TRUE)),
       col = curve.id, lty = curve.id, lwd = 2, bty = "n")

## ------------------------------------------------------------
##
## Competing risks: select an event code from the saved event map.
##
## ------------------------------------------------------------
data(follic, package = "randomForestSRC")
crfit <- rfsrc(Surv(time, status) ~ ., data = follic, ntree = 100, nsplit = 3)
age <- unique(as.numeric(quantile(crfit$xvar$age)))
crp <- partial(crfit, partial.type = "years.lost",
               partial.xvar = "age", partial.values = age)
event <- crp$event.info$event.type[1]
crd <- get.partial.plot.data(crp, target = event)
plot(crd$x, crd$yhat, type = "b", xlab = "Age", ylab = "Partial years lost")
}
}
\keyword{partial}
