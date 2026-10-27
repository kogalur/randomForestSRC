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
  Estimates conditional response distributions and their quantiles
  using random forests. The main method estimates a distribution of
  training residuals and shifts it by the forest mean prediction.
  The default
  \code{method = "exchangeable"} gives residuals equal weight;
  \code{method = "forest"} uses forest weights. An optional
  Greenwald--Khanna approximation is available for large data sets.
  Supports continuous responses in univariate, multivariate, and
  mixed-outcome forests.
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
  \item{object}{The original training object returned by \code{quantreg},
    including its saved forest, training responses, and residuals.}
  \item{newdata}{Data for new observations, containing the same predictor
    variables as the training data. Converted to a plain data frame.
    Response columns are optional and are used for performance evaluation.
    If omitted or \code{NULL}, calculate results for the training
    observations.}
  \item{method}{The residual method adds residual quantiles to the
    forest mean prediction. Choose how the residuals are weighted:
    \describe{
      \item{\code{"exchangeable"}}{Equal weights (default).}
      \item{\code{"forest"}}{Forest weights that depend on the
        predictor values.}
    }
    For large data sets, \code{"gk"} uses the Greenwald--Khanna quantile
    approximation instead of the residual method. It also accepts
    \code{"GK"}, \code{"G-K"}, and \code{"g-k"}.
    The former name \code{"local"} is an alias for
    \code{"exchangeable"}. To reuse a nondefault method, specify it in
    each call to \code{quantreg(object = fit, ...)}.}
  \item{splitrule}{Splitting rule used to grow the forest. The default
    for univariate regression is \code{"la.quantile.regr"}, local
    adaptive quantile regression splitting. Multivariate forests use
    their default rule when either \code{"la.quantile.regr"} or
    \code{"quantile.regr"} is specified. Other applicable rules include
    \code{"mse"} for regression and \code{"mahalanobis"} for
    multivariate regression.}
  \item{prob}{Quantile levels to calculate, as a nonempty numeric vector
    of finite probabilities strictly between zero and one. Can be
    supplied during fitting or prediction; levels are sorted and
    duplicates retained. With \code{NULL}, the residual methods calculate
    percentiles 1 through 99 when fitting. With \code{object} supplied,
    \code{prob = NULL} reuses the saved levels. GK training defaults
    also depend on \code{prob.epsilon}; supply \code{prob} to control
    the number of quantiles.}
  \item{prob.epsilon}{Greenwald--Khanna approximation tolerance on the
    probability scale, a finite number in \eqn{(0,1)}. Larger values allow
    coarser summaries. With \code{object} supplied, \code{NULL} reuses
    the saved tolerance. When both \code{prob} and \code{prob.epsilon}
    are \code{NULL} during GK fitting, the default requests \eqn{2n}
    quantiles, where \eqn{n} is the processed training sample size.}
  \item{oob}{Use out-of-bag (OOB) predictions for training observations
    if \code{TRUE} (default), or predictions from all trees if
    \code{FALSE}. When fitting, this also selects the training residuals
    for the residual methods. Calls with \code{object} reuse the saved
    residual choice. Predictions for \code{newdata} always use all trees.
    See \emph{OOB predictions} in Details.}
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
    under \emph{Evaluating predictions}. Supplied levels apply to this
    call only.}
  \item{m.target}{For \code{get.pinball.error}, the name of one continuous
    response to evaluate. \code{NULL} evaluates all available continuous
    responses. Extracts results from \code{o} without recalculating
    predictions.}
  \item{...}{Additional options for \code{rfsrc} when growing a forest,
    or for its prediction method when \code{object} is supplied.
    \code{quantreg} sets the quantile method and any required
    forest-weight options.}
}
\details{
  
  \code{quantreg()} estimates the conditional cumulative distribution
  function (CDF) of a response \eqn{Y} given predictor values \eqn{X=x},
  \deqn{F(s\mid x)=P(Y\leq s\mid X=x),} and obtains conditional
  quantiles from this distribution.  This is different than
  nonparametric regression which targets the conditional mean
  \eqn{E(Y\mid X=x)}. Intervals formed from the conditional quantiles
  describe individual response values; they are not confidence intervals
  for the conditional mean.

  The residual method estimates this conditional distribution using
  forest mean predictions and the empirical distribution of training
  residuals. For training observation \eqn{j}, the residual is
  \deqn{R_j = Y_j - \widehat m_j,}
  where \eqn{\widehat m_j} is the training prediction selected when
  fitting the forest. For an observation with predictor values
  \eqn{x}, adding each residual to the forest mean prediction
  \eqn{\widehat m(x)} gives adjusted response values
  \eqn{\widehat m(x)+R_j}. Their weighted empirical distribution
  estimates \eqn{F(s\mid x)}.

  Let \eqn{\mathcal I_R} index the finite training residuals, and let
  \eqn{w_j(x)} be their nonnegative weights, summing to one. The
  residual CDF is
  \deqn{\widehat F_R(u\mid x)
    = \sum_{j\in\mathcal I_R} w_j(x) I\{R_j\leq u\}.}
  Its quantile \eqn{\widehat Q_R(\tau\mid x)} is the smallest residual
  whose cumulative weight is at least \eqn{\tau}. Adding the forest
  mean prediction gives the response quantile and shifts the CDF:
  \deqn{\widehat Q(\tau\mid x)
    = \widehat m(x) + \widehat Q_R(\tau\mid x),\qquad
    \widehat F(s\mid x) = \widehat F_R\{s-\widehat m(x)\mid x\}.}
  Predicted quantiles can fall outside the range of the training responses.

  The default \code{method = "exchangeable"} uses
  \eqn{w_j(x)=1/N_R}, where \eqn{N_R} is the number of finite
  training residuals. This uses a common residual distribution for
  all predictor values. The estimated conditional response distributions
  therefore have the same shape and spread, with location determined
  by \eqn{\widehat m(x)}. With the default OOB residuals, prediction
  intervals follow the signed OOB-error construction in Zhang et al.
  (2020, Section 2.3). Equal weights give the same residual quantiles
  as \code{quantile(..., type = 1)}. When calculating quantiles for a
  training observation, its own finite residual remains in the pool.

  With \code{method = "forest"}, the same residuals receive forest
  weights that depend on \eqn{x}. This allows both the shape and
  spread of the predicted distribution to vary with the predictors.

  Meinshausen (2006) estimates the conditional CDF by weighting the
  training responses directly:
  \deqn{\widehat F_{\mathrm{QRF}}(s\mid x)
    = \sum_{j=1}^n v_j(x) I\{Y_j\leq s\},}
  where \eqn{v_j(x)} are forest weights. The \code{"forest"} option
  here instead weights residuals and then shifts their distribution
  by \eqn{\widehat m(x)}. It assigns probability to the adjusted
  responses \eqn{Y_j-\widehat m_j+\widehat m(x)}.
  The two estimators generally differ even if the same weights are used.

  \subsection{OOB predictions}{
    When fitting a forest, \code{oob} determines which trees supply
    the training predictions used to calculate residuals.
    With \code{oob = TRUE} (default),
    \eqn{\widehat m_j} averages predictions from trees whose training
    samples exclude observation \eqn{j}. With \code{oob = FALSE},
    it uses all trees. The fitted forest retains this choice of
    training residuals.

    When the original fitted object is supplied as \code{object},
    \code{oob} selects OOB or full-ensemble mean predictions, forest
    weights, or GK quantiles for the training observations. It does
    not change the saved residuals. For example,
    \code{quantreg(object = fit, oob = FALSE)} uses all trees for the
    training predictions while keeping the residuals selected when
    \code{fit} was created.

    For \code{newdata}, predictions and weights always use all trees.
    The residual methods reuse the saved training residuals, so a
    default fit uses OOB residuals for new observations as well.
    When fitting and prediction are combined in one call, \code{oob}
    selects the training residuals before predicting the new data.

    If an entire OOB prediction vector, weight matrix, or GK quantile
    matrix is unavailable (\code{NULL}), the corresponding result from
    all trees is used. Individual missing entries are not replaced.
  }

  \subsection{Greenwald--Khanna method for large data}{
    For large data sets, \code{method = "gk"} uses the
    Greenwald--Khanna (2001) algorithm for compact approximate quantile
    summaries. This option obtains quantiles from the forest without
    using training residuals or requesting a forest-weight matrix.
    Specify \code{prob} to control the number of returned quantiles.
    The argument \code{prob.epsilon} controls the approximation tolerance
    on the probability scale; larger values allow coarser summaries.

    The CDF is reconstructed from the returned quantiles.
    At each grid value, it takes the largest requested probability
    whose quantile is no greater than that value. The CDF is zero below
    all requested quantiles and reaches one at the largest finite
    training response. Thus \code{prob} controls the quantiles and the
    resolution of the reconstructed CDF, while \code{maxn} controls
    where that CDF is evaluated. GK means and standard deviations use
    the CDF increments and depend on the grid and the probability mass
    it represents.
  }

  \subsection{Fitting and results}{
    Supply \code{formula} and \code{data} to fit a forest, and save the
    result as \code{fit}. The argument \code{prob} specifies the quantile
    levels; for example, \code{prob = c(.05, .50, .95)} requests the
    5th percentile, median, and 95th percentile. Use
    \code{quantreg(object = fit, newdata = testdata)} to predict new
    observations, or \code{quantreg(object = fit, prob = c(.1, .5, .9))}
    to calculate different quantiles for the training observations.
    These calls reuse the fitted trees and require the original
    training object.

    To fit and predict in one call, supply \code{formula}, \code{data},
    and \code{newdata} together. New data can contain predictors alone.
    Include response columns to calculate performance scores or use
    \code{plot.quantreg()}. These responses are used only for evaluation;
    residual distributions and the CDF reporting grid use training data.

    The residual methods calculate the mean and standard deviation of
    the predicted distribution from the same residuals and weights.
    With \eqn{\bar R(x)=\sum_{j\in\mathcal I_R} w_j(x)R_j}, they are
    \deqn{\widehat\mu(x)=\widehat m(x)+\bar R(x),\qquad
      \widehat\sigma(x)=\left\{\sum_{j\in\mathcal I_R} w_j(x)
      [R_j-\bar R(x)]^2\right\}^{1/2}.}
    The exchangeable method uses equal weights. Residuals are used
    without recentering, so the distribution mean can differ from the
    forest mean prediction.

    To report the CDF, the function constructs a grid \code{yunq} from
    the sorted unique finite training responses, separately for each
    continuous response. When there are more than \code{maxn} values,
    it selects approximately equally spaced positions in this ordered
    set. Use \code{maxn = Inf} to include all values.

    For the residual methods, \code{cdf} evaluates the predicted CDF at
    these grid values. Its last value may be less than one when the
    shifted residual distribution extends above the grid. Quantiles,
    means, and standard deviations use the full residual distribution,
    including values outside the grid, and do not depend on \code{maxn}.
    The component \code{density} contains interval probabilities:
    its first column is the CDF at the first grid value, and subsequent
    columns are successive CDF differences. These probabilities are
    not divided by interval widths. Their row sum equals the last
    reported CDF value, without adjustment for probability above the
    grid.

    Use \code{get.quantile(o)} to extract quantiles and
    \code{get.quantile.stat(o)} to obtain the mean, median, and standard
    deviation for each observation. In \code{get.quantile()}, each
    requested \code{target.prob} is matched to the nearest stored
    probability, while column labels show the requested levels.
    To calculate new levels, use
    \code{quantreg(object = fit, prob = ...)} with the original training
    object. The median returned by \code{get.quantile.stat()} uses
    the stored level nearest 0.5; include 0.5 in \code{prob} to obtain
    that exact level. Use \code{extract.quantile(o)} to obtain all
    quantile components in a named list, one entry per continuous
    response.

    With \code{method = "forest"}, the weight matrix has one row per
    evaluated observation and one column per training observation,
    which can require substantial memory. Reducing \code{maxn}
    reduces the \code{cdf} and \code{density} matrices; it does not
    reduce the forest-weight matrix or the number of GK quantiles.
  }

  \subsection{Evaluating predictions}{
    \code{get.quantile.crps(o)} compares the predicted CDF with observed
    responses using the continuous ranked probability score (CRPS).
    At each reporting-grid value \eqn{s_k}, it averages
    \eqn{[I\{Y_i\leq s_k\}-\widehat F(s_k\mid X_i)]^2} over the
    selected observations with finite response and CDF values. It
    integrates these errors by the trapezoidal rule from the first
    grid value to \eqn{s_k}, returning a curve with columns \code{y}
    and \code{crps}. Smaller scores indicate less prediction error.

    The last entry integrates over the full grid range, so the score
    depends on both the grid range and resolution. With
    \code{standardize = TRUE}, each integral is divided by its
    integration width. The first entry is then \code{NA} because its
    width is zero; with \code{standardize = FALSE}, the first entry
    is zero.

    \code{get.pinball.error(o)} evaluates individual quantiles. For
    response \eqn{Y_i} and predicted quantile \eqn{\widehat Q_i(\tau)},
    the pinball loss is
    \deqn{\rho_\tau\{Y_i-\widehat Q_i(\tau)\},\qquad
      \rho_\tau(u)=u\{\tau-I(u<0)\}.}
    The helper returns the mean loss at each requested level, using
    the selected observations with finite responses and quantiles.
    Smaller values indicate better quantile predictions. If no
    observations are available at a level, its mean loss is \code{NA}.

    For the residual methods, pinball scoring uses stored quantiles
    at exactly matching levels, or calculates them from the saved
    residual distribution. Levels in \code{tau} need not appear in
    \code{prob}, and the score does not depend on the CDF reporting
    grid. GK scoring interpolates the reported CDF. Scoring uses the
    predictions and residual choice saved in \code{o}, without
    refitting the forest or recalculating its mean predictions.

    For a univariate quantile object with observed responses and
    performance output, \code{print.rfsrc} reports raw and standardized
    grid-integrated CRPS and mean pinball losses. For a multivariate
    object, use \code{print(object, outcome.target = "response.name")}
    to display quantile scores for one continuous response. Without
    \code{outcome.target}, the display gives mean and response-specific
    regression or classification errors. Selecting a response for
    printing or extraction does not recalculate predictions. OOB labels
    refer to the predictions used in each calculation.

    Set pinball levels for one print call with
    \code{print(object, quantreg.tau = c(.2, .5, .8))}, or for the
    session with \code{options(quantreg.tau = c(.2, .5, .8))}.
    Explicit levels override the session option. Otherwise the defaults
    are taken, in order, from \code{quantreg.tau}, the former option
    \code{rfsrc.pinball.taus}, or \code{c(.1, .5, .9)}.
    These defaults also apply to \code{get.pinball.error()}.
    Invalid levels cause an error. These settings control reporting
    without modifying the fitted object.

    To display pinball losses on a plot, use
    \code{plot.quantreg(object, quantreg.tau = c(.2, .5, .8))}.
    The losses appear in the lower-right corner of the main panel;
    the quantile interval and CRPS inset are unchanged. Pinball labels
    are omitted unless \code{quantreg.tau} is supplied in the plotting
    call, even when the session option is set.
  }
}
\value{
  \code{quantreg()} returns an \code{rfsrc} training (grow) or prediction
  object with the additional class \code{quantreg}. For a univariate
  regression response, \code{o$quantreg} contains:
  \item{quantiles}{Matrix with evaluated observations in rows and
    requested probability levels in columns, including when there is
    only one observation or level.}
  \item{prob}{Sorted probability levels corresponding to \code{quantiles}.}
  \item{cdf}{CDF matrix with evaluated observations in rows and
    \code{yunq} grid values in columns.}
  \item{density}{CDF increments, with the same dimensions as \code{cdf}.}
  \item{yunq}{Ordered training-response grid used to report the CDF.}
  \item{mean, std}{Mean and standard deviation of the estimated
    conditional response distribution, as vectors for the residual
    methods. The standard deviation describes response variability,
    not the standard error of a mean estimate.
    \code{NULL} for GK; use \code{get.quantile.stat()} for its grid-based
    summaries.}
  \item{method}{Method used in this call: \code{"exchangeable"},
    \code{"forest"}, or \code{"gk"}. The former name \code{"local"}
    is reported as \code{"exchangeable"}.}
  \item{oob}{Whether the reported mean predictions or GK quantiles
    use OOB output.}
  \item{prediction}{Forest mean predictions used to shift the
    residual distribution. \code{NULL} for GK.}
  \item{residual.source}{Source of the saved training residuals:
    \code{"oob"} or \code{"ensemble"} for the residual methods.
    \code{NULL} for GK. A new-data result can have
    \code{residual.source = "oob"} and \code{oob = FALSE}: the
    training residuals are OOB and the new predictions use all trees.}

  For multivariate and mixed-outcome forests, \code{o$quantreg} is a
  named list of these components, one per returned continuous response.
  Classification responses have no quantile component.

  \subsection{Training data and residuals}{
    \code{o$xvar} contains predictors for the observations being
    predicted; \code{o$yvar} contains their responses, if supplied.
    The saved forest's \code{xvar} and \code{yvar} contain the training
    data. The CDF grid and null CRPS reference use training responses;
    performance scores use responses for the observations being
    evaluated.

    For a univariate training object,
    \code{o$residual = o$yvar - o$predicted} and
    \code{o$residual.oob = o$yvar - o$predicted.oob}, when the respective
    predictions are available. For multivariate and mixed outcomes,
    these vectors are in \code{o$regrOutput[[response.name]]}.
    If \code{predicted} or \code{predicted.oob} is \code{NULL}, its
    residual vector is also \code{NULL}. Missing OOB predictions give
    missing residuals.

    Both sets of training residuals are saved in
    \code{o$forest$residual} and \code{o$forest$residual.oob}.
    They are vectors for univariate regression and matrices for
    multivariate or mixed outcomes, with training observations in rows
    and named continuous responses in columns. If residuals are
    available for only some responses, the remaining columns contain
    \code{NA}. If none are available, the residual matrix is \code{NULL}.
    \code{o$forest$quantreg$residual.source} records which set was
    selected for each response. Prediction objects retain these
    training residuals under \code{forest}.

    Nonfinite residuals are excluded together with their weights, and
    the remaining forest weights are normalized for each observation.
    Quantiles and summaries are \code{NA} if no finite residuals remain,
    the forest mean prediction is nonfinite, or the remaining weights
    include nonfinite values or their sum is not positive.
    Negative weights and incompatible weight-matrix dimensions cause
    an error.
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
  of quantile summaries. \emph{Proceedings of the 2001 ACM SIGMOD
  International Conference on Management of Data}, 58--66.
  \url{https://sanjeevkhanna.org/papers/sigmod01_quantiles.pdf}.

  Meinshausen N. (2006).  Quantile regression forests.
  \emph{Journal of Machine Learning Research}, 7:983--999.
  \url{https://jmlr.org/papers/v7/meinshausen06a.html}.

  Zhang H., Zimmerman J., Nettleton D. and Nordman D.J. (2020).
  Random forest prediction intervals. \emph{The American Statistician},
  74(4):392--406. \doi{10.1080/00031305.2019.1585288}.
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

## Check the stored levels and probability covered by the CDF grid.
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
## Forest-weighted residuals and prediction for new observations
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

## Exchangeable prediction reuses the saved OOB residuals.
o.exchangeable <- quantreg(object = o, newdata = x.test,
                            method = "exchangeable")
print(head(get.quantile(o.exchangeable)))
print(identical(o.exchangeable$forest$residual.oob, o$residual.oob))
print(o.exchangeable$quantreg$residual.source)  # "oob": training residuals
print(o.exchangeable$quantreg$oob)              # FALSE: new observations
print(c(evaluated.rows = o.exchangeable$n,
        training.rows = length(o.exchangeable$forest$residual.oob)))

## Calculate a single level using the original training object.
o.median <- quantreg(object = o, newdata = x.test,
                     method = "forest", prob = .5)
print(head(get.quantile(o.median)))

## ------------------------------------------------------------
## Changing maxn leaves residual quantiles and summaries unchanged
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

## Recalculate all responses; select Temp for extraction and printing.
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
## GK quantiles and alternative regression splitting
## ------------------------------------------------------------
gk <- quantreg(alcohol ~ ., data = wine, method = "gk",
               prob = c(.05, .25, .50, .75, .95),
               prob.epsilon = .01, maxn = 100, ntree = 100)
print(head(get.quantile(gk)))

mse <- quantreg(alcohol ~ ., data = wine, method = "exchangeable",
                splitrule = "mse", prob = c(.05, .50, .95), ntree = 100)
print(head(get.quantile(mse)))

## ------------------------------------------------------------
## Larger data set: Iowa housing
## ------------------------------------------------------------

data(housing,  package = "randomForestSRC")

## Impute missing values using random splitting.
iowa <- housing
iowa$PID <- NULL
iowa$SalePrice <- log(iowa$SalePrice)
iowa <- impute(SalePrice ~. , iowa, splitrule = "random", nimpute = 1)

## Use fewer trees and larger terminal nodes for speed.
o <- quantreg(SalePrice ~., iowa, ntree = 50, nodesize = 20)
plot.quantreg(o, prbL=.05, prbU=.95)

}}
\keyword{quantile regression forests}
