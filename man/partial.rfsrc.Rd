\name{partial.rfsrc}
\alias{partial.rfsrc}
\alias{partial}
\alias{get.partial.plot.data}
\title{Compute Partial Dependence Values}
\description{
  Compute predictions for partial dependence analysis using a fitted forest.
  This advanced interface is the workhorse called by
  \code{plot.variable(..., partial = TRUE)}, which provides automatic grid
  selection and plotting. Use \code{partial()} directly for custom predictor
  grids, additional fixed predictors, several survival evaluation times,
  or observation-level predictions.

  \code{partial()} returns predictions for each training observation.
  \code{get.partial.plot.data()} extracts their means or individual values
  for custom plots. Supported outcomes include regression, classification,
  multivariate and mixed outcomes, right-censored survival, and competing
  risks. See \code{\link{plot.variable.rfsrc}} for the simpler plotting
  interface.
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
    or an object of class \code{(rfsrc, forest)}. The forest must contain
    its training data. Anonymous forests are not supported.}
  \item{oob}{Logical. If \code{TRUE} (default), predict each observation
    using trees that excluded it during training. If \code{FALSE}, use
    the non-OOB ensemble.}
  \item{partial.type}{Prediction summary for survival outcomes. For
    right-censored survival, choose \code{"mort"}, \code{"rel.freq"},
    \code{"surv"}, or \code{"chf"}. For competing risks, choose
    \code{"years.lost"}, \code{"cif"}, or \code{"chf"}. Specify one type
    for these families; see Details for their meanings. Ignored for
    regression and classification.}
  \item{partial.xvar}{Name of the predictor to vary.}
  \item{partial.values}{Nonempty vector of finite numeric values for
    \code{partial.xvar}. Values are sorted and duplicates removed before
    calculation. For a factor, supply integer codes from its training
    level order; see Details.}
  \item{partial.xvar2}{Optional character vector naming additional
    predictors to hold fixed. Names must be distinct and different from
    \code{partial.xvar}.}
  \item{partial.values2}{Numeric vector giving one finite fixed value per
    predictor in \code{partial.xvar2}, in the same order. Factors use
    their training level codes. Each value applies to all observations
    throughout the primary predictor grid.}
  \item{partial.time}{Nonempty vector of finite evaluation times for
    \code{"surv"}, \code{"cif"}, or \code{"chf"}. Times must be at least
    the first value in the forest's \code{time.interest}. See Details for
    matching requests to the time grid. Ignored for \code{"mort"},
    \code{"rel.freq"}, and \code{"years.lost"}.}
  \item{get.tree}{Optional vector of tree indices. By default all trees
    are considered, with OOB eligibility applied when \code{oob = TRUE}.}
  \item{seed}{Negative integer specifying the random-number seed.}
  \item{do.trace}{Number of seconds between progress updates; \code{FALSE}
    suppresses updates.}
  \item{...}{Additional internal prediction options. Supply graphical
    arguments to the subsequent plotting call.}
  \item{o}{An object returned by \code{partial()}.}
  \item{target}{For classification, a class label or its numeric position
    in the selected response's levels. For competing risks, an event code
    from \code{o$event.info$event.type}. The first class or event is used
    when omitted. Ignored for regression and right-censored survival.}
  \item{m.target}{Response name to extract from a multivariate or mixed
    \code{partial} object. The default is the first response in
    \code{o$yvar.names}.}
  \item{granule}{Logical used by \code{get.partial.plot.data()}. If
    \code{FALSE} (default), return the mean prediction at each grid value.
    If \code{TRUE}, return individual predictions with a grid value for
    each observation. This differs from the numeric \code{granule}
    argument of \code{plot.variable()}.}
}
\details{
  \subsection{Partial predictions and averaging}{
    At each value in \code{partial.values}, \code{partial()} sets
    \code{partial.xvar} to that value for every training observation and
    predicts with the fitted forest. The other predictors stay at each
    observation's training values, except for any additional predictors
    fixed by \code{partial.xvar2}. No trees are refitted.

    For a predictor \eqn{X_j}, with no additional predictors fixed, write
    the prediction for observation \eqn{i} at value \eqn{z} as
    \eqn{P_i(z) = \widehat f(z, X_{i,-j})}.
    Averaging these predictions gives the partial dependence at \eqn{z}:
    the model's mean prediction when \eqn{X_j} is set to \eqn{z}, averaged
    over the training values of the other predictors.

    \code{partial()} returns the individual \eqn{P_i(z)} values;
    \code{get.partial.plot.data()} computes their means. With
    \code{oob = TRUE}, each \eqn{P_i(z)} uses the trees for which
    observation \eqn{i} was OOB. This changes the contributing trees,
    while the training observations still provide the values over which
    the mean is calculated.
  }
  \subsection{Choosing predictor values}{
    Choose a grid for the range of interest, such as equally spaced
    values or selected training quantiles. Each grid point requires a
    prediction for every observation, so fewer points reduce computation
    and output size.

    For factors, use the integer codes corresponding to the training
    levels. For example, levels \code{c("low", "middle", "high")} have
    codes 1, 2, and 3. Use \code{match(labels, levels(training.factor))}
    to obtain codes for selected labels. The extraction helper returns
    numeric codes, which can be labeled in the plotting call.
  }
  \subsection{Fixing additional predictors}{
    Use \code{partial.xvar2} and \code{partial.values2} to hold additional
    predictors at specified values for every observation. For example,
    vary wind over a grid while fixing temperature at one value. The
    resulting curve averages predictions over the remaining predictors.

    To compare curves at several temperatures, repeat the call with a
    different fixed temperature each time. This also provides the values
    for a two-dimensional display: one predictor varies along the grid,
    and the other changes between calls.
  }
  \subsection{Selecting a response}{
    \code{get.partial.plot.data()} extracts predicted responses for
    regression and the selected class probability for classification.
    For multivariate and mixed outcomes, \code{partial()} computes all
    responses in one call. Use \code{m.target} in the extraction helper
    to choose a response, and \code{target} to choose a class when that
    response is a factor. Several responses or classes can be extracted
    from the same partial object without recalculating predictions.

    For competing risks, \code{target} selects an event code from
    \code{o$event.info$event.type}. This differs from \code{plot.variable()},
    where \code{target} selects the event's position in that vector.
  }
  \subsection{Survival summaries and evaluation times}{
    For right-censored survival, \code{"mort"} returns mortality,
    \code{"surv"} returns survival probabilities, and \code{"chf"} returns
    cumulative hazards. In \code{partial()} and its extraction helper,
    \code{"rel.freq"} returns the same unnormalized mortality as
    \code{"mort"}.

    For competing risks, \code{"years.lost"} returns expected years lost,
    \code{"cif"} returns cumulative incidence, and \code{"chf"} returns
    cause-specific cumulative hazard. The extraction helper selects the
    event of interest.

    Supply \code{partial.time} for \code{"surv"}, \code{"cif"}, or
    \code{"chf"}. Predictions use the forest's \code{time.interest} grid.
    Each request uses the last grid time less than or equal to it;
    a request after the final grid time uses that final time. Using
    values from \code{time.interest} gives exact matches. The returned
    \code{partial.time} contains the requested times, which can differ
    from the grid times used for prediction.

    Mortality and years lost are scalar summaries and do not use
    \code{partial.time}.
  }
  \subsection{Extracting results for custom plots}{
    Pass the partial object to \code{get.partial.plot.data()} to obtain
    predictor values \code{x} and corresponding predictions \code{yhat}.
    The \code{granule} argument controls the summary:
    \describe{
      \item{\code{granule = FALSE}}{Return the mean at each grid value.
        Only finite predictions contribute to the mean; if none are
        finite, the mean is \code{NA}.}
      \item{\code{granule = TRUE}}{Return the individual predictions,
        grouped by grid value and then in training-observation order.
        Each grid value is repeated once per observation. Nonfinite
        predictions become \code{NA}.}
    }
    For time-dependent survival summaries, \code{yhat} has one column
    per requested time. The individual values are on the original
    prediction scale. Unlike the boxplot values from
    \code{plot.variable(..., partial = TRUE)}, they are not shrunk toward
    their mean.

    For automatic plots and numerical summaries, use
    \code{plot.variable(..., partial = TRUE)}. Its \code{partial.summary}
    component contains means, standard deviations, row-averaging standard
    errors, and finite prediction counts. See
    \code{\link{plot.variable.rfsrc}} for the interpretation of its
    standard errors and display bands.
  }
}
\value{
  \code{partial()} returns an object of class \code{(rfsrc, partial, family)}.
  It includes \code{call}, \code{family}, \code{partial.values} (the grid
  used), \code{partial.time} (the requested times), \code{yvar.names}, and
  \code{event.info}, together with the prediction components below.
  Let \eqn{n} be the number of training observations, \eqn{K} the number
  of predictor grid values, \eqn{T} the number of requested times, and
  \eqn{J} the number of events.
  \item{regrOutput}{Named list with one \eqn{n \times K} matrix per
    regression response. Rows correspond to observations and columns
    to predictor grid values.}
  \item{classOutput}{Named list with one
    \eqn{n \times (1 + C) \times K} array per classification response,
    where \eqn{C} is its number of classes. The second dimension is named
    \code{"all"} followed by the class labels; the class-labeled columns
    contain class probabilities.}
  \item{survOutput}{For right-censored mortality (\code{"mort"} or
    \code{"rel.freq"}), an \eqn{n \times K} matrix. For survival or
    cumulative hazard, an \eqn{n \times T \times K} array. For
    competing-risk years lost, an \eqn{n \times J \times K} array.
    For competing-risk cumulative incidence or cumulative hazard, an
    \eqn{n \times T \times J \times K} array.}

  \code{get.partial.plot.data()} returns a list with:
  \item{x}{Predictor grid values. With \code{granule = TRUE}, each value
    is repeated \eqn{n} times to match the individual predictions.}
  \item{yhat}{Mean predictions when \code{granule = FALSE}, or individual
    predictions when \code{granule = TRUE}. Scalar summaries give a
    vector of length \eqn{K} or \eqn{nK}, respectively. Time-dependent
    summaries give a matrix with \eqn{K} or \eqn{nK} rows and \eqn{T}
    columns.}
  \item{partial.time}{The evaluation times supplied to \code{partial()}.}
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
