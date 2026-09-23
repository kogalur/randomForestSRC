\name{quantreg.rfsrc}
\alias{quantreg.rfsrc}
\alias{quantreg}
\alias{extract.quantile}
\alias{get.quantile}
\alias{get.quantile.stat}
\alias{get.quantile.crps}
\alias{get.pinball.error}
\title{Quantile Regression Forests}
\description{
  Estimates conditional quantiles of continuous responses using random
  forests. Returns predicted percentiles, cumulative distribution
  functions (CDFs), and distribution summaries for training or new
  observations. Supports univariate, multivariate, and mixed-outcome
  forests, with quantile estimates for each continuous response.
}
\usage{\method{quantreg}{rfsrc}(formula, data, object, newdata,
  method = "exchangeable", splitrule = NULL, prob = NULL, prob.epsilon = NULL,
  oob = TRUE, fast = FALSE, maxn = 1e3, ...)

extract.quantile(o)
get.quantile(o, target.prob = NULL, pretty = TRUE)
get.quantile.stat(o, pretty = TRUE)
get.quantile.crps(o, pretty = TRUE, subset = NULL, standardize = TRUE)
get.pinball.error(o, tau = NULL, subset = NULL, m.target = NULL,
                 pretty = TRUE)
}
\arguments{
  \item{formula}{Model formula. Required when growing a new forest and
    ignored when \code{object} is supplied. At least one response must
    be continuous.}
  \item{data}{Training data containing the responses and predictors.
    Required when growing a new forest. Converted to a plain data frame.}
  \item{object}{The original fitted object returned by \code{quantreg},
    including its saved forest, training responses, and residuals.
    Use this object for subsequent prediction calls.}
  \item{newdata}{Data for new observations, containing the same predictor
    variables as the training data. Converted to a plain data frame.
    Response columns are optional and are used for performance evaluation.
    Use \code{quantreg(object = fit, newdata = testdata)} to obtain
    quantile output. Without \code{newdata}, or with \code{newdata = NULL},
    results are calculated for the training observations.}
  \item{method}{How to estimate the response distribution:
    \describe{
      \item{\code{"exchangeable"}}{Add quantiles of the equally weighted
        training-residual distribution to each forest mean prediction
        (default).}
      \item{\code{"forest"}}{Use forest weights to estimate a
        predictor-specific residual distribution, then add its quantiles
        to the forest mean prediction.}
      \item{\code{"gk"}}{Compute quantiles with the Greenwald--Khanna
        algorithm. Also accepts \code{"GK"}, \code{"G-K"}, and
        \code{"g-k"}.}
    }
    The former name \code{"local"} is an alias for
    \code{"exchangeable"}. The method applies to the current call;
    specify it again in later calls to use a nondefault method.}
  \item{splitrule}{Splitting rule used to grow the forest. The default
    for univariate regression is \code{"la.quantile.regr"}, local
    adaptive quantile regression splitting. Multivariate forests use
    their default rule when either \code{"la.quantile.regr"} or
    \code{"quantile.regr"} is specified. Other applicable rules include
    \code{"mse"} for regression and \code{"mahalanobis"} for
    multivariate regression.}
  \item{prob}{Quantile levels to calculate, as a nonempty numeric vector
    of finite probabilities strictly between zero and one. Can be
    supplied during fitting or prediction; levels are sorted, including
    any duplicates. With \code{NULL}, the residual methods calculate
    percentiles 1 through 99 during fitting, and an \code{object} call
    uses the probabilities saved during fitting. GK training defaults
    also depend on \code{prob.epsilon}; supply \code{prob} to control
    the number of quantiles.}
  \item{prob.epsilon}{Greenwald--Khanna probability approximation
    tolerance, a finite number in \eqn{(0,1)}. Larger values allow
    coarser summaries. A prediction-time \code{NULL} uses the saved
    training tolerance. When both \code{prob} and \code{prob.epsilon}
    are \code{NULL} during GK fitting, the default requests \eqn{2n}
    quantiles, where \eqn{n} is the processed training sample size.}
  \item{oob}{Use out-of-bag (OOB) predictions for training observations?
    During fitting, \code{TRUE} also selects OOB training residuals for
    the residual methods; \code{FALSE} selects residuals from
    full-ensemble predictions. Later calls reuse that residual selection.
    New observations always use full-ensemble predictions. See
    \emph{Training residuals and OOB predictions} in Details for the
    distinction between the saved residuals and current predictions.}
  \item{fast}{If \code{TRUE}, grow a new forest with \code{rfsrc.fast}
    instead of \code{rfsrc}. Ignored when \code{object} is supplied.}
  \item{maxn}{Maximum number of response values in the CDF reporting
    grid, as a positive integer or \code{Inf}. Controls the resolution
    of \code{cdf}, \code{density}, and grid-integrated CRPS. For a fixed
    forest and prediction call, it does not change the quantiles or the
    residual-method means and standard deviations.}
  \item{o}{A training or prediction object returned by \code{quantreg}.
    CRPS and pinball scoring require observed responses in \code{o$yvar}.}
  \item{target.prob}{Quantile levels to extract with \code{get.quantile}.
    \code{NULL} returns all stored levels. Otherwise, supply finite
    probabilities in \eqn{[0,1]}; each is matched to the nearest stored
    level.}
  \item{pretty}{If \code{TRUE}, simplify a helper result for one
    continuous response to a matrix, data frame, or numeric vector, as
    described in Value. If \code{FALSE}, return a response-named list.}
  \item{subset}{Rows to evaluate for CRPS or pinball loss. Supply positive
    integer indices or a logical vector with one nonmissing entry per
    row of the returned object. \code{NULL} uses all rows. The selection
    must be nonempty; repeated indices include a row more than once.}
  \item{standardize}{If \code{TRUE}, divide each cumulative CRPS integral
    by its integration width.}
  \item{tau}{Quantile levels at which to calculate pinball loss, as
    finite probabilities strictly between zero and one. Scores follow
    the supplied order, including repeated levels. \code{NULL} uses
    the session option \code{quantreg.tau}, with the defaults described
    under \emph{Printed performance}. Supplied levels apply to this
    call only.}
  \item{m.target}{For \code{get.pinball.error}, the name of one continuous
    response to evaluate. \code{NULL} evaluates all available continuous
    responses. Selects results from \code{o} without a prediction call.}
  \item{...}{Additional options for \code{rfsrc} when growing a forest,
    or for its prediction method when \code{object} is supplied.
    \code{quantreg} sets the quantile method and any required
    forest-weight options.}
}
\details{
  Quantile regression describes the distribution of a response given its
  predictors. Whereas a regression forest gives a mean prediction,
  \code{quantreg()} estimates lower, middle, and upper quantiles. These
  estimates describe variation in response values and can be used to
  form prediction intervals. For example,
  \code{prob = c(.05, .50, .95)} requests the estimated 5th percentile,
  median, and 95th percentile for each observation.

  The default method combines each forest mean prediction with a common
  distribution of training residuals. The \code{"forest"} method allows
  this residual distribution to vary with the predictors, and
  \code{"gk"} computes quantiles with the Greenwald--Khanna algorithm.
  Each method returns quantiles and a CDF for every evaluated observation.

  \subsection{Fitting and prediction}{
    Supply \code{formula} and \code{data} to fit a forest. Save the
    resulting training object, for example as \code{fit}, and use
    \code{quantreg(object = fit, newdata = testdata)} to predict new
    observations. Use \code{quantreg(object = fit, prob = c(.1, .5, .9))}
    to calculate different quantiles for the training observations
    without growing new trees. Subsequent calls require the original
    training object, rather than a prediction object.

    New data can contain predictors alone. Include response columns to
    calculate performance scores or use \code{plot.quantreg()}.
    These responses are used for evaluation; the saved training data
    continue to define the residual distributions and CDF reporting grid.
    Training and new-data prediction can also be requested in a single
    call by supplying \code{formula}, \code{data}, and \code{newdata}.
  }

  \subsection{Training residuals and OOB predictions}{
    With the default \code{oob = TRUE}, the residual methods use
    \deqn{R_j = Y_j - \widehat m_j^{\mathrm{OOB}},}
    where \eqn{\widehat m_j^{\mathrm{OOB}}} averages predictions from
    trees whose training samples exclude observation \eqn{j}. With
    \code{oob = FALSE} during fitting, residuals instead use
    full-ensemble training predictions. This choice is saved during
    fitting and used in all subsequent calls.

    For training observations, \code{oob} also selects the current mean
    predictions, forest weights, or GK quantiles. Changing \code{oob}
    in a later call changes these predictions, but not the saved
    residual selection. New observations use full-ensemble predictions
    together with the previously selected residuals. Thus a fit made
    with the default settings uses OOB training residuals even when
    predicting new observations. When fitting and new-data prediction
    are combined in one call, \code{oob} first selects the training
    residuals before the new observations are predicted.

    If the requested OOB output component is absent, full-ensemble
    output is used instead. Missing entries within an existing OOB
    component stay missing; nonfinite residuals are excluded from the
    residual distribution. The common exchangeable distribution
    includes each training observation's finite OOB residual, including
    when calculating that observation's quantiles. GK quantiles are
    computed directly by the forest and do not use residuals.
  }

  \subsection{Exchangeable residual method}{
    The default \code{method = "exchangeable"} assigns equal weight to
    all finite training residuals. Their empirical CDF is
    \deqn{\widehat F_R(u)
      = \frac{1}{N_R}\sum_{j\in\mathcal I_R} I\{R_j\leq u\},}
    where \eqn{\mathcal I_R} indexes the \eqn{N_R} finite residuals.
    Let \eqn{\widehat Q_R(\tau)} be its quantile at probability
    \eqn{\tau}, and \eqn{\widehat m(x)} the current forest mean
    prediction. The predicted response quantile and CDF are
    \deqn{\widehat Q(\tau\mid x)
      = \widehat m(x) + \widehat Q_R(\tau),\qquad
      \widehat F(s\mid x) = \widehat F_R\{s-\widehat m(x)\}.}

    Every observation therefore has the same residual-distribution
    shape and spread, shifted by its own forest prediction.
    Exchangeability is the working assumption behind pooling residuals
    in this way; OOB describes how the default residuals are calculated.
    This is the locally adjusted residual-CDF construction associated
    with Zhang et al. (2019).
  }

  \subsection{Forest-weighted residual method}{
    With \code{method = "forest"}, each observation has a weighted
    distribution of the saved training residuals:
    \deqn{\widehat F(s\mid x)
      = \sum_j w_j(x) I\{R_j\leq s-\widehat m(x)\}.}
    The forest weights \eqn{w_j(x)} sum to one over the usable residuals.
    Because the weights depend on \eqn{x}, both the shape and spread of
    the residual distribution can vary with the predictors. The
    weights are applied to residuals, rather than directly to training
    responses as in Meinshausen (2006).

    Training observations use OOB weights when \code{oob = TRUE} and
    full-ensemble weights otherwise. New observations use full-ensemble
    weights. The weight matrix has one row per evaluated observation
    and one column per training observation, which can require
    substantial memory for large data sets.
  }

  \subsection{Residual quantiles and distribution summaries}{
    Both residual methods calculate a quantile by ordering the residuals
    and selecting the first whose cumulative weight reaches the
    requested probability. They then add the current forest prediction.
    This is the lower generalized inverse; equal weights give the
    type-1 sample quantile convention. Predicted quantiles can extend
    beyond the observed training-response range.

    The mean and standard deviation describe the same predicted
    distribution. With \eqn{\bar R(x)=\sum_j w_j(x)R_j}, they are
    \deqn{\widehat\mu(x)=\widehat m(x)+\bar R(x),\qquad
      \widehat\sigma(x)=\left\{\sum_j w_j(x)
      [R_j-\bar R(x)]^2\right\}^{1/2}.}
    The exchangeable method uses equal weights. Residuals are used
    without recentering, so the distribution mean can differ from the
    forest mean prediction. These calculations use the full usable
    residual distribution, independently of the CDF reporting grid.

    Nonfinite residuals are excluded together with their weights, and
    the remaining forest weights are normalized for each observation.
    Quantiles and summaries are \code{NA} when there are no usable
    residuals, the current prediction is nonfinite, or the weights for
    finite residuals contain nonfinite values or have no positive sum.
    Negative weights and incompatible weight-matrix dimensions cause
    an error.
  }

  \subsection{Greenwald--Khanna method}{
    The \code{"gk"} method computes quantiles with the native
    Greenwald--Khanna algorithm, without requesting an explicit
    forest-weight matrix. The returned quantiles come directly from
    this algorithm. A step CDF is constructed from their values and
    probability levels: at each response-grid value, it uses the largest
    requested probability whose quantile is no greater than that value.
    It is zero below all requested quantiles and reaches one at the
    largest finite training response.

    For this method, \code{prob} controls both the requested quantiles
    and the resolution of the reconstructed CDF. The argument
    \code{maxn} controls only the grid at which that CDF is evaluated.
  }

  \subsection{CDF reporting grid and storage}{
    The reporting grid, \code{yunq}, contains sorted unique finite
    training responses, separately for each continuous response. When
    there are more than \code{maxn} values, the function selects
    approximately equally spaced positions in this ordered set.
    Use \code{maxn = Inf} to include all values.

    For the residual methods, \code{cdf} evaluates the predicted CDF at
    these grid values. Its last value may be less than one when the
    shifted residual distribution extends above the grid. Quantiles,
    means, and standard deviations still use the entire usable residual
    distribution, including values outside the grid.

    The component \code{density} contains interval probabilities, not
    density values divided by interval widths. Its first column is the
    CDF at the first grid value; subsequent columns are successive CDF
    differences. The row sum equals the final reported CDF value,
    without adjustment for probability above the grid.

    Quantile storage grows with the number of evaluated observations
    times \code{length(prob)}. CDF and interval-probability storage grow
    with that observation count times the grid size. Reducing
    \code{maxn} reduces these latter matrices, but not the forest-weight
    matrix or the number of GK quantiles.
  }

  \subsection{Extracting results}{
    \describe{
      \item{\code{extract.quantile(o)}}{Returns all quantile components
        as a named list, with one entry per continuous response.}
      \item{\code{get.quantile(o)}}{Returns the stored quantiles. Use
        \code{target.prob} to select particular levels. The helper
        matches each requested level to the nearest stored probability;
        column labels show the requested levels. To calculate additional
        levels, call \code{quantreg()} on the original training object
        with those levels in \code{prob}.}
      \item{\code{get.quantile.stat(o)}}{Returns \code{mean},
        \code{median}, and \code{std} for each observation. The residual
        methods use the distribution summaries above. GK means and
        standard deviations are calculated from the reporting-grid
        CDF increments and depend on that grid and its represented
        probability mass. The median uses the stored level nearest
        0.5; include 0.5 in \code{prob} to obtain that exact level.}
    }
    Quantile matrices have observations in rows and probability levels
    in columns, including when there is only one row or level.
  }

  \subsection{CRPS: evaluating the predicted distribution}{
    \code{get.quantile.crps(o)} assesses the predicted CDF using the
    observed responses. It returns a finite-grid continuous ranked
    probability score (CRPS) curve, with columns \code{y} and
    \code{crps}. Smaller scores indicate less prediction error.

    At each reporting-grid value \eqn{s_k}, it averages
    \eqn{[I\{Y_i\leq s_k\}-\widehat F(s_k\mid X_i)]^2} over the
    selected observations with finite response and CDF values. It then
    integrates these errors by the trapezoidal rule from the first grid
    value to \eqn{s_k}. With \code{standardize = TRUE}, each integral is
    divided by its integration width.

    The last entry summarizes the entire reporting grid. Integration
    covers only the range of this grid. The first
    standardized entry is \code{NA} because its integration width is
    zero; with \code{standardize = FALSE}, the first entry is zero.
    Grid resolution and range therefore affect this score.
  }

  \subsection{Pinball loss: evaluating individual quantiles}{
    \code{get.pinball.error(o)} evaluates predictions at specified
    quantile levels. For response \eqn{Y_i} and predicted quantile
    \eqn{\widehat Q_i(\tau)}, the loss is
    \deqn{\rho_\tau\{Y_i-\widehat Q_i(\tau)\},\qquad
      \rho_\tau(u)=u\{\tau-I(u<0)\}.}
    The helper returns the mean loss at each requested level; smaller
    values indicate better quantile predictions. Each mean uses the
    selected observations with finite responses and quantiles at that
    level. If none are available, the result is \code{NA}.

    For \code{"exchangeable"} and \code{"forest"}, the helper uses
    exactly matching stored quantiles when available and otherwise
    calculates them from the saved residual distribution. Thus
    \code{tau} need not be included in \code{prob}, and scoring is
    independent of the reporting grid. GK scoring uses interpolation of
    the reporting CDF. Scoring uses the object's current predictions
    and saved residual selection, without refitting or a prediction call.
  }

  \subsection{Printed performance}{
    For a univariate quantile object with observed responses and
    performance output, \code{print.rfsrc} reports raw and standardized
    grid-integrated CRPS and mean pinball losses. For a multivariate
    object, use \code{print(object, outcome.target = "response.name")}
    to display quantile scores for one continuous response. Without
    \code{outcome.target}, the multivariate display gives mean and
    response-specific regression or classification errors. Selecting a
    response for printing or extraction does not recalculate predictions.
    The OOB label for quantile scores follows that response's quantile
    predictions; the regression error is labeled according to the
    predictions used for that error.

    Set levels for one print call with
    \code{print(object, quantreg.tau = c(.2, .5, .8))}, or set defaults
    for the session with \code{options(quantreg.tau = c(.2, .5, .8))}.
    Supplied levels override the session option. When neither is set,
    the former option \code{rfsrc.pinball.taus} is used if available,
    followed by the default \code{c(.1, .5, .9)}. These settings also
    determine the default \code{tau} for \code{get.pinball.error()}.
    They control reporting only and do not modify the fitted object.
    Invalid levels cause an error.

    To display numerical pinball losses on a plot, supply levels in
    the plotting call, for example,
    \code{plot.quantreg(object, quantreg.tau = c(.2, .5, .8))}. The
    losses appear as text in the lower-right corner of the main panel.
    The plotted quantile interval and CRPS inset are unchanged.
    Without this argument, plotting omits the pinball labels even when
    the session option is set.
  }
}
\value{
  \code{quantreg()} returns an \code{rfsrc} training (grow) or prediction
  object with the additional class \code{quantreg}. For a univariate
  regression response, \code{o$quantreg} contains:
  \item{quantiles}{Matrix with evaluated observations in rows and
    requested probability levels in columns.}
  \item{prob}{Sorted probability levels corresponding to \code{quantiles}.}
  \item{cdf}{CDF matrix with evaluated observations in rows and
    \code{yunq} grid values in columns.}
  \item{density}{CDF increments, with the same dimensions as \code{cdf}.}
  \item{yunq}{Ordered training-response grid used to report the CDF.}
  \item{mean, std}{Predicted distribution mean and standard deviation
    vectors for the residual methods. \code{NULL} for GK; use
    \code{get.quantile.stat()} for its grid-based summaries.}
  \item{method}{Method used in this call: \code{"exchangeable"},
    \code{"forest"}, or \code{"gk"}. The former name \code{"local"}
    is reported as \code{"exchangeable"}.}
  \item{oob}{Whether the current mean predictions or native GK
    quantiles use OOB output.}
  \item{prediction}{Current forest mean predictions used to shift the
    residual distribution. \code{NULL} for GK.}
  \item{residual.source}{Source of the saved training residuals:
    \code{"oob"} or \code{"ensemble"} for the residual methods.
    \code{NULL} for GK. A new-data result can have
    \code{residual.source = "oob"} and \code{oob = FALSE}: the
    training residuals are OOB, but the new predictions use the full
    ensemble.}

  For multivariate and mixed-outcome forests, \code{o$quantreg} is a
  named list of these components, one per returned continuous response.
  Classification responses have no quantile component.

  \subsection{Training data and residuals}{
    The usual \code{o$xvar} and \code{o$yvar} describe the observations
    being predicted, with \code{yvar} available when responses are
    supplied for evaluation. The saved forest's \code{xvar} and
    \code{yvar} describe the training data. The CDF grid and null CRPS
    reference use training responses; performance scores use the
    responses for the observations being evaluated.

    For a univariate training object,
    \code{o$residual = o$yvar - o$predicted} and
    \code{o$residual.oob = o$yvar - o$predicted.oob}, when the respective
    predictions are available. For multivariate and mixed outcomes,
    these vectors are in \code{o$regrOutput[[response.name]]}.
    An unavailable prediction component gives a \code{NULL} residual
    component; missing entries within OOB predictions stay missing in
    the corresponding residuals.

    The saved forest also contains both training-residual components
    for future prediction calls. They are vectors for univariate
    regression or matrices with training observations in rows and named
    continuous responses in columns. An unavailable response has an
    \code{NA} column in a partly available matrix; a wholly unavailable
    residual component is \code{NULL}.
    \code{o$forest$quantreg$residual.source} identifies the original
    selection for each response. Prediction objects include these
    training residuals under \code{forest}, separately from their
    current predictions.
  }

  \subsection{Helper functions}{
    \code{extract.quantile()} always returns a response-named list of
    quantile components. The other helpers return the following for one
    continuous response when \code{pretty = TRUE}, and a response-named
    list otherwise:
    \describe{
      \item{\code{get.quantile}}{An observation-by-probability matrix,
        with column names such as \code{q.50}.}
      \item{\code{get.quantile.stat}}{A data frame with columns
        \code{mean}, \code{median}, and \code{std}, one row per observation.}
      \item{\code{get.quantile.crps}}{A data frame with columns \code{y}
        and \code{crps}, one row per reporting-grid value.}
      \item{\code{get.pinball.error}}{A numeric vector of mean losses,
        with names such as \code{"tau=0.2"}, one entry per requested level.}
    }
  }
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Greenwald M. and Khanna S. (2001).  Space-efficient online computation
  of quantile summaries. \emph{Proceedings of ACM SIGMOD}, 30(2):58--66.

  Meinshausen N. (2006).  Quantile regression forests.
  \emph{Journal of Machine Learning Research}, 7:983--999.

  Zhang H., Zimmerman J., Nettleton D. and Nordman D.J. (2019).
  Random forest prediction intervals. \emph{The American Statistician}.
}
\seealso{
  \code{\link{rfsrc}}, \code{\link{predict.rfsrc}},
  \code{\link{plot.quantreg}}
}
\examples{
## ------------------------------------------------------------
## A basic analysis using default settings
## ------------------------------------------------------------
o <- quantreg(Temp ~ ., data = na.omit(airquality))
plot.quantreg(o)

\donttest{
## ------------------------------------------------------------
## Exchangeable OOB residual quantiles for wine alcohol
## ------------------------------------------------------------
data(wine, package = "randomForestSRC")
set.seed(17)
prob <- c(.05, .25, .50, .75, .95)
o <- quantreg(alcohol ~ ., data = wine, method = "exchangeable",
              oob = TRUE, prob = prob, ntree = 100)
print(head(get.quantile(o)))
print(head(get.quantile(o, c(.25, .50, .75))))
print(head(get.quantile.stat(o)))

## The common residual quantiles are added to each OOB mean prediction.
r.oob <- as.numeric(o$residual.oob)
r.oob <- r.oob[is.finite(r.oob)]
q.res <- as.numeric(quantile(r.oob, probs = prob, type = 1))
q.oob <- outer(as.numeric(o$predicted.oob), q.res, "+")
print(all.equal(unname(get.quantile(o)), unname(q.oob)))
print(o$quantreg$residual.source)

## Inspect the stored levels and reporting-grid probability mass.
print(o$quantreg$prob)
print(summary(rowSums(o$quantreg$density)))
crps <- get.quantile.crps(o)
print(crps)
plot(crps, type = "l")

## Pinball losses at nondefault levels: no refitting is needed.
print(o, quantreg.tau = c(.2, .5, .8))
print(get.pinball.error(o, tau = c(.2, .5, .8)))
plot.quantreg(o, quantreg.tau = c(.2, .5, .8))

## Optionally use the same reporting levels throughout the session.
op <- options(quantreg.tau = c(.2, .5, .8))
print(getOption("quantreg.tau"))
print(o)
options(op)

## ------------------------------------------------------------
## Forest-weighted residuals and predictor-only deployment
## ------------------------------------------------------------
set.seed(23)
train <- sample.int(nrow(wine), floor(.7 * nrow(wine)))
o <- quantreg(alcohol ~ ., data = wine[train, ],
              method = "forest", prob = prob, ntree = 100)
o.test <- quantreg(object = o, newdata = wine[-train, ],
                   method = "forest")
print(head(get.quantile(o.test)))
print(tail(get.quantile.crps(o.test, standardize = FALSE), 1))

x.test <- wine[-train, setdiff(names(wine), "alcohol"), drop = FALSE]
o.predict <- quantreg(object = o, newdata = x.test, method = "forest")
print(head(get.quantile(o.predict)))

## Exchangeable prediction uses the same saved OOB residual bank.
o.exchangeable <- quantreg(object = o, newdata = x.test,
                            method = "exchangeable")
print(head(get.quantile(o.exchangeable)))
print(identical(o.exchangeable$forest$residual.oob, o$residual.oob))
print(o.exchangeable$quantreg$residual.source)  # "oob": training residuals
print(o.exchangeable$quantreg$oob)              # FALSE: new-row predictions
print(c(evaluated.rows = o.exchangeable$n,
        training.rows = length(o.exchangeable$forest$residual.oob)))

## Calculate a single level using the original grow object.
o.median <- quantreg(object = o, newdata = x.test,
                     method = "forest", prob = .5)
print(head(get.quantile(o.median)))

## ------------------------------------------------------------
## Grid resolution does not determine residual quantiles or moments
## ------------------------------------------------------------
coarse <- quantreg(object = o, newdata = x.test, method = "exchangeable",
                  prob = prob, maxn = 2, seed = 37)
fine <- quantreg(object = o, newdata = x.test, method = "exchangeable",
                prob = prob, maxn = Inf, seed = 37)
print(all.equal(get.quantile(coarse), get.quantile(fine)))
print(all.equal(get.quantile.stat(coarse), get.quantile.stat(fine)))

## ------------------------------------------------------------
## Multivariate and mixed outcomes
## ------------------------------------------------------------
dta <- na.omit(airquality)
mv <- quantreg(cbind(Ozone, Temp) ~ ., data = dta,
               splitrule = "mahalanobis", prob = prob, ntree = 100)
q.mv <- get.quantile(mv, pretty = FALSE)
print(names(q.mv))
print(head(q.mv$Ozone))
print(head(q.mv$Temp))
print(head(mv$regrOutput$Temp$residual.oob))
print(identical(mv$regrOutput$Temp$residual.oob,
                as.numeric(mv$forest$residual.oob[, "Temp"])))

## Restore all responses; select Temp only for extraction and printing.
mv.restore <- quantreg(object = mv, prob = .5)
q.restore <- get.quantile(mv.restore, pretty = FALSE)
print(names(q.restore))
print(head(q.restore$Temp))
print(get.mv.error(mv.restore))
print(mv.restore, outcome.target = "Temp")

dta$Month <- factor(dta$Month)
mixed <- quantreg(cbind(Ozone, Temp, Month) ~ ., data = dta,
                  prob = prob, ntree = 100)
print(names(extract.quantile(mixed)))

## ------------------------------------------------------------
## Native GK quantiles and alternative regression splitting
## ------------------------------------------------------------
gk <- quantreg(alcohol ~ ., data = wine, method = "gk",
               prob = c(.05, .25, .50, .75, .95),
               prob.epsilon = .01, maxn = 100, ntree = 100)
print(head(get.quantile(gk)))

mse <- quantreg(alcohol ~ ., data = wine, method = "exchangeable",
                splitrule = "mse", prob = c(.05, .50, .95), ntree = 100)
print(head(get.quantile(mse)))

## ------------------------------------------------------------
## Larger data set; iowa housing
## ------------------------------------------------------------

data(housing,  package = "randomForestSRC")

## the original data contains lots of missing data; use fast imputation
iowa <- housing
iowa$PID <- NULL
iowa$SalePrice <- log(iowa$SalePrice)
iowa <- impute(SalePrice ~. , iowa, splitrule = "random", nimpute = 1)

## use fewer trees and shallow trees for speed
o <- quantreg(SalePrice ~., iowa, ntree = 50, nodesize = 20)
plot.quantreg(o, prbL=.05, prbU=.95)

}}
\keyword{quantile regression forests}
