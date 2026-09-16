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
  Estimates conditional quantiles for continuous responses using a
  forest prediction and the training-residual distribution, or the
  native Greenwald--Khanna quantile algorithm.  Returns quantiles,
  CDF values, and distribution summaries for training or new observations.
  Univariate, multivariate, and mixed-outcome forests are supported.
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
  \item{formula}{Model formula.  Required when growing a new forest and
    ignored when \code{object} is supplied.  At least one response must
    be continuous.}
  \item{data}{Training data containing the response and predictor
    variables.  Required when growing a new forest.  Input is converted
    to a plain data frame before formula processing.}
  \item{object}{The original \code{quantreg} grow object, including its
    saved forest, training responses, and residuals.  Retain this object
    for subsequent prediction calls.}
  \item{newdata}{Optional data containing the training predictors.
    Input is converted to a plain data frame.  Response columns can also
    be supplied for performance evaluation.  Use \code{quantreg}, rather
    than a direct \code{predict} call, to obtain quantile output for new
    observations.  Omitted or \code{NULL} \code{newdata} restores
    predictions for the training observations.}
  \item{method}{Quantile calculation method.  The default
    \code{"exchangeable"} treats the training residuals as exchangeable,
    pools their finite values with equal weights, and adds quantiles
    of this common residual distribution to the forest prediction.
    With the default \code{oob = TRUE} during growth, this is an
    OOB exchangeable residual distribution.  \code{"forest"} uses a
    forest-weighted residual distribution for each evaluated row.
    \code{"gk"} requests native Greenwald--Khanna quantiles;
    \code{"GK"}, \code{"G-K"}, and \code{"g-k"} are also accepted.
    The former name \code{"local"} is retained as an alias for
    \code{"exchangeable"}.  The method applies to the current call.
    Specify it again in later calls to retain a nondefault method.}
  \item{splitrule}{Splitting rule used to grow the forest.  The default
    for univariate regression is \code{"la.quantile.regr"}, local
    adaptive quantile regression splitting.  Multivariate forests use
    their default splitting rule when either \code{"la.quantile.regr"}
    or \code{"quantile.regr"} is specified.  Other applicable rules can
    be supplied, including \code{"mse"} for regression and
    \code{"mahalanobis"} for multivariate regression.}
  \item{prob}{Nonempty numeric vector of finite probabilities strictly
    between zero and one.  Levels are sorted; duplicates are retained.
    Can be supplied at training or prediction time.  With \code{NULL},
    the residual methods use percentiles 1 through 99 during training,
    and an object call uses the saved training probabilities.  For GK,
    the training defaults also depend on \code{prob.epsilon}; specify
    \code{prob} explicitly to control the number of native quantiles.}
  \item{prob.epsilon}{Finite probability approximation tolerance in
    \eqn{(0,1)} for Greenwald--Khanna.  Larger values permit coarser
    summaries.  A prediction-time \code{NULL} uses the saved training
    tolerance.  When both probability arguments are \code{NULL} during
    GK growth, the default requests \eqn{2n} quantiles, where \eqn{n} is
    the processed training sample size.}
  \item{oob}{Selects the residual distribution and the current training-row
    predictions.  During growth, \code{TRUE} selects residuals computed
    as the observed response minus its OOB mean prediction, when that
    component is available.  These residuals define the default
    exchangeable distribution and are also used by the forest method.
    \code{FALSE} during growth selects full-ensemble residuals instead.
    Both numerical quantities are retained as \code{residual} and
    \code{residual.oob} when their predictions are available.  For
    training-row output, \code{oob} also selects the current mean
    predictions, forest weights, or native GK quantiles.  Changing it
    in a later \code{object} call leaves the original residual selection
    unchanged.  Supplied \code{newdata} uses ordinary full-ensemble
    predictions together with the selected training residuals.  See the
    OOB residual discussion below.}
  \item{fast}{Use \code{rfsrc.fast} instead of \code{rfsrc} when growing
    a new forest?  Does not regrow an existing object.}
  \item{maxn}{Positive integer, or \code{Inf}, limiting the number of
    ordered, unique training-response values in the CDF reporting grid.
    Controls the resolution of \code{cdf}, \code{density}, and the
    grid-integrated CRPS curve.  For a fixed forest and prediction call,
    it does not change the returned quantiles or the residual-method
    means and standard deviations.}
  \item{o}{A grow or prediction object returned by \code{quantreg}.
    CRPS and pinball scoring also require observed responses in
    \code{o$yvar}.}
  \item{target.prob}{Probabilities to extract with \code{get.quantile}.
    \code{NULL} returns all stored levels.  Otherwise, supply finite
    values in \eqn{[0,1]}; each is matched to the nearest stored level.}
  \item{pretty}{Simplify the helper result when there is one continuous
    response?  \code{TRUE} returns a matrix, data frame, or numeric
    vector, as described in Value.  \code{FALSE} retains a response-named
    list.}
  \item{subset}{Rows to evaluate for CRPS or pinball loss.  Supply positive
    integer indices or a nonmissing logical vector with one entry per
    row of the returned object.  \code{NULL} uses all rows.  The selection
    must be nonempty; repeated indices include a row repeatedly.}
  \item{standardize}{Divide each cumulative CRPS integral by its
    integration width?  The default is \code{TRUE}.}
  \item{tau}{Finite pinball-loss levels strictly between zero and one.
    Order and duplicates are retained.  \code{NULL} uses the session
    option \code{quantreg.tau}, with the fallbacks described under
    Printed performance.  Explicit levels apply to this call only.}
  \item{m.target}{For \code{get.pinball.error}, the name of one continuous
    response to evaluate.  \code{NULL} evaluates all available continuous
    responses.  Selects stored results without making a prediction call.}
  \item{...}{Additional forest options.  Training options are selected
    from those recognized by \code{rfsrc}; object calls forward options
    to the prediction method.  Options controlled by \code{quantreg},
    including the quantile method and associated forest weights, take
    precedence.}
}
\details{
  \subsection{Training and prediction}{
    Supply \code{formula} and \code{data} to grow a forest, or
    \code{object} to reuse its saved forest.  \code{newdata} can be
    supplied with either form of call.  Training responses and residuals
    are saved once and reused by subsequent calls.  Response values in
    \code{newdata} are used for evaluation; they do not replace this
    training information.
  }

  \subsection{OOB residuals and their exchangeable distribution}{
    With \code{oob = TRUE} during growth and OOB predictions available,
    the saved residual for training observation \eqn{j} is
    \deqn{R_j^{\mathrm{OOB}} = Y_j - \widehat m_j^{\mathrm{OOB}},}
    where \eqn{\widehat m_j^{\mathrm{OOB}}} averages predictions from
    trees whose training samples exclude observation \eqn{j}.  Thus,
    each residual compares the observed response with a prediction from
    trees that did not use that observation to grow.

    The exchangeable method pools these residuals into the empirical CDF
    \deqn{\widehat F_R^{\mathrm{OOB}}(u)
      = \frac{1}{N_R}\sum_{j\in\mathcal I_R}
        I\{R_j^{\mathrm{OOB}}\leq u\},}
    where \eqn{\mathcal I_R} indexes the finite residuals and
    \eqn{N_R=|\mathcal I_R|}.  This is the OOB exchangeable residual
    distribution: each finite residual receives equal weight, and the
    same residual distribution is used at every predictor row.
    Exchangeability is the working assumption motivating this common
    distribution; OOB specifies how its residuals are constructed.
    Its lower-inverse quantile is denoted by
    \eqn{\widehat Q_R^{\mathrm{OOB}}(\tau)}.

    For training-row output using OOB mean predictions, the quantile is
    \deqn{\widehat Q_j^{\mathrm{OOB}}(\tau)
      = \widehat m_j^{\mathrm{OOB}}
        + \widehat Q_R^{\mathrm{OOB}}(\tau).}
    The common residual distribution includes each observation's finite
    OOB residual.  For a new predictor row \eqn{x}, the quantile is
    \deqn{\widehat Q(\tau\mid x)
      = \widehat m(x) + \widehat Q_R^{\mathrm{OOB}}(\tau),}
    where \eqn{\widehat m(x)} is the full-ensemble prediction from the
    saved forest.  The OOB residual distribution is retained even though
    the new-row mean prediction uses the full ensemble.

    With \code{oob = FALSE} during growth, the selected residuals
    instead use full-ensemble training predictions.  Both residual
    quantities are retained when available, but the original selection
    is reused in subsequent calls.  Changing \code{oob} in an
    \code{object} call changes the current training-row output, not
    that residual selection.
    When growth and \code{newdata} are requested in the same call,
    \code{oob} first selects the training residual bank and then ordinary
    predictions are used for the new rows.

    If an OOB output component is absent, the existing full-ensemble
    fallback is used.  Missing values within an available OOB component
    remain missing; nonfinite residuals are excluded from the empirical
    distribution.  Each residual-method quantile component records
    \code{residual.source = "oob"} or \code{"ensemble"} separately
    from its \code{oob} setting for the current predictions.  A
    new-data result can therefore have \code{residual.source = "oob"}
    and \code{oob = FALSE}.  Its numerical training OOB residuals
    remain in \code{forest$residual.oob}; the suffix \code{.oob}
    identifies the OOB values, not a logical flag.

    The forest method uses the same saved residual bank but assigns
    row-specific forest weights.  For training rows it requests OOB
    weights when \code{oob = TRUE} and full-ensemble weights otherwise.
    For GK, \code{oob} selects the native OOB or full-ensemble quantile
    output; GK quantiles are not calculated from residuals.
  }

  \subsection{Response data and residuals in the returned object}{
    The usual \code{xvar} and \code{yvar} components describe the
    observations for which predictions were returned.  In a grow object
    they contain training data.  In a new-data prediction object they
    contain the retained prediction data, with \code{yvar} available
    when responses were retained for evaluation.  The saved forest's
    \code{xvar} and \code{yvar} continue to describe the training data.
    Quantile reporting grids and the null CRPS reference use these
    training responses, whereas performance scores use the current
    object's responses.

    For a univariate grow object \code{o}, \code{o$residual} is
    \code{o$yvar - o$predicted}, and \code{o$residual.oob} is
    \code{o$yvar - o$predicted.oob}.  For multivariate and mixed
    outcomes, these vectors are stored beside the corresponding
    predictions in \code{o$regrOutput[[response.name]]}.  No residual
    is substituted for an unavailable OOB prediction.

    The saved forest retains \code{residual} and \code{residual.oob}
    for use with future prediction rows.  These are vectors for
    univariate regression and matrices for multivariate or mixed
    outcomes, with training observations in rows and named continuous
    responses in columns.  An unavailable component is \code{NULL}
    in the grow output; its column in a partially available saved bank
    is \code{NA}.  A wholly unavailable bank is \code{NULL}.

    The named setting \code{forest$quantreg$residual.source} records
    the original residual selection for each continuous response.
    Residual-method quantile components also report their
    \code{residual.source}.  Neither setting replaces a numerical
    \code{residual.oob} component.  A prediction object retains the
    training residuals under \code{forest}; it does not place these
    training-length vectors beside its current predictions.  Retain
    the original grow object for further \code{quantreg} calls.
  }

  \subsection{Residual quantiles}{
    For one continuous response, write \eqn{R_j} for its saved
    residuals, constructed as described above, and
    \eqn{\widehat m(x)} for the current forest mean prediction.

    The default \code{"exchangeable"} method adds quantiles of the
    common empirical residual distribution to the current prediction:
    \deqn{\widehat Q(\tau\mid x)
      = \widehat m(x) + \widehat Q_R(\tau).}
    Its CDF is
    \deqn{\widehat F(s\mid x)
      = \frac{1}{N_R}\sum_{j\in\mathcal I_R}
        I\{R_j \leq s-\widehat m(x)\},}
    where \eqn{\mathcal I_R} indexes its \eqn{N_R} finite residuals.
    Each row uses the same residual-distribution shape, translated by
    its own forest
    prediction.  This is the locally adjusted residual-CDF construction
    associated with Zhang et al. (2019).

    The \code{"forest"} method weights the residuals using forest
    weights for the evaluated row:
    \deqn{\widehat F(s\mid x)
      = \sum_j w_j(x) I\{R_j \leq s-\widehat m(x)\}.}
    These weights allow the residual distribution to vary with the
    predictor row.  They are applied to residuals, rather than directly
    to training responses as in Meinshausen (2006).

    Both methods use the lower generalized inverse.  After ordering
    residuals and their weights, the residual quantile is the first
    residual whose cumulative weight reaches \eqn{\tau}.  The current
    forest prediction is then added.  Equal weights give the inverse
    empirical-CDF convention, corresponding to a type-1 sample quantile.
    Quantiles can extend beyond the observed training-response range.

    The mean and standard deviation are computed from the same residual
    distribution.  With \eqn{\bar R(x)=\sum_j w_j(x)R_j}, they are
    \deqn{\widehat\mu(x)=\widehat m(x)+\bar R(x),\qquad
      \widehat\sigma(x)=\left\{\sum_j w_j(x)
      [R_j-\bar R(x)]^2\right\}^{1/2}.}
    Equal weights are used for \code{"exchangeable"}.  The residuals are not
    recentered to force their weighted mean to zero, so the distribution
    mean need not equal the forest prediction.

    Nonfinite residuals and their associated weight columns are removed
    together.  Remaining forest weights are normalized separately for
    each row.  A row with a nonfinite current prediction, nonfinite
    retained weights, or no positive usable weight has unavailable
    quantiles and summaries.  No usable residuals also yields unavailable
    output.  Negative weights and incompatible weight dimensions are
    rejected.
  }

  \subsection{Greenwald--Khanna}{
    The \code{"gk"} method requests quantiles from the native
    Greenwald--Khanna algorithm.  It avoids requesting the explicit
    forest-weight matrix.  The returned quantiles are retained directly;
    a reporting CDF is then reconstructed from their values and
    probability levels.  At each grid value this CDF uses the largest
    requested probability whose quantile is at or below that value,
    and reaches one at the upper training-response endpoint.
    Consequently, \code{prob} affects both native quantiles and the
    resolution of this reconstructed CDF.  \code{maxn} only controls
    its evaluation grid.
  }

  \subsection{CDF reporting grid}{
    The grid, \code{yunq}, consists of sorted unique finite training
    responses, separately for each continuous response.  When there are
    more than \code{maxn} values, an approximately equally spaced
    selection of their order positions is retained.  \code{maxn = Inf}
    retains all values.

    For the residual methods, \code{cdf} evaluates the residual-based
    CDF at these grid values.  These evaluations do not determine the
    primary quantiles or moments.  In particular, the last CDF value can
    be less than one when a shifted residual lies above the reporting
    grid; the primary quantiles and moments still use that residual.

    \code{density} contains CDF increments.  Its first column is the
    CDF at the first grid value, and later columns are successive
    differences.  Each increment is probability in a reporting interval,
    rather than a density divided by grid spacing.  Its row sum equals
    the last reported CDF value.  These increments are not renormalized
    to represent omitted upper-tail mass.

    Quantile storage is proportional to the number of evaluated rows
    times \code{length(prob)}.  CDF and increment storage are
    proportional to that row count times the grid size.  Forest weights
    additionally require one column per training observation.  Reducing
    \code{maxn} does not reduce the forest-weight matrix or the number
    of native GK quantiles requested.
  }

  \subsection{Extracting results}{
    \code{extract.quantile(o)} returns a named list with one
    component per returned continuous response, including univariate
    objects.

    \code{get.quantile(o)} extracts stored quantiles.  With
    \code{target.prob = NULL}, all stored levels are returned.  Requested
    probabilities are sorted and deduplicated, then matched to the
    nearest stored levels.  Column labels use the requested probabilities.
    Quantiles are not recomputed by this helper.  To calculate additional
    levels, call \code{quantreg} on the original grow object with those
    levels in \code{prob}.

    \code{get.quantile.stat(o)} returns columns \code{mean},
    \code{median}, and \code{std}.  For the residual methods, the mean
    and standard deviation use the direct residual summaries described
    above.  The median retains the nearest-stored lookup at probability
    0.5; include 0.5 in \code{prob} to obtain that exact level.  For GK,
    the mean and standard deviation retain the reporting-grid moment
    calculation from the CDF increments, so these summaries depend on
    the grid and its represented probability mass.

    Matrices retain their row and column dimensions when only one row,
    probability, or grid value is present.
  }

  \subsection{CRPS}{
    \code{get.quantile.crps(o)} returns a score curve with columns
    \code{y} and \code{crps}.  Observed responses must be present.

    At each grid value \eqn{s_k}, the helper averages
    \eqn{[I\{Y_i\leq s_k\}-\widehat F(s_k\mid X_i)]^2} over the
    selected rows with finite response and CDF values.  It accumulates
    the trapezoidal integral from the first grid value to \eqn{s_k}.
    With \code{standardize = TRUE}, each integral is divided by its
    integration width.  Zero-width standardized entries are \code{NA};
    the corresponding raw integrals are zero.  The last entry summarizes
    the integral over the retained grid.  Integration outside that grid
    is not included.
  }

  \subsection{Pinball loss}{
    \code{get.pinball.error(o)} returns the mean pinball loss at each
    requested level; smaller values are better.  For observed response
    \eqn{Y_i} and predicted quantile \eqn{\widehat Q_i(\tau)}, the loss is
    \deqn{\rho_\tau\{Y_i-\widehat Q_i(\tau)\},\qquad
      \rho_\tau(u)=u\{\tau-I(u<0)\}.}
    Averaging uses the selected rows with finite responses and quantiles,
    separately at each level; no usable rows gives \code{NA}.

    For \code{"exchangeable"} and \code{"forest"}, scoring uses exact
    stored levels or quantiles calculated from the saved residual
    distribution, independently of the reporting grid.  Thus \code{tau}
    need not have been included in \code{prob}.  GK retains interpolation
    of its reporting CDF.  Scoring preserves the object's prediction mode
    and residual selection without refitting or calling prediction.
  }

  \subsection{Printed performance}{
    For univariate quantile output with observed responses and performance
    output, \code{print.rfsrc} reports raw and standardized grid-integrated
    CRPS and pinball losses.  A selected continuous response from a
    multivariate object can be printed using
    \code{print(object, outcome.target = "response.name")}.
    To extract that response, use its name in the list returned by
    \code{get.quantile(object, pretty = FALSE)}.  These selections act on
    the returned results; restoration and prediction retain all responses.
    Without \code{outcome.target}, the multivariate print displays the
    existing performance summary, with the mean error followed by the
    response-specific errors, rather than response-specific quantile scores.
    Quantile-score OOB labels follow the quantile calculation for that
    response.  The requested regression error retains its own provenance.

    Pinball levels can be supplied for one print call using
    \code{print(object, quantreg.tau = c(.2, .5, .8))}, or set for the
    session using \code{options(quantreg.tau = c(.2, .5, .8))}.
    Explicit levels take precedence over the session option.  When
    neither is supplied, the former \code{rfsrc.pinball.taus} option is
    accepted for compatibility, followed by the default
    \code{c(0.1, 0.5, 0.9)}.  These are reporting settings, not
    arguments to \code{quantreg}, and do not modify the fitted object.
    Levels must be finite and lie strictly between zero and one; they
    need not have been included in \code{prob}.  Invalid levels give
    an explanatory error rather than silently suppressing the scores.

    \code{plot.quantreg(object, quantreg.tau = c(.2, .5, .8))} adds the
    same losses to the main plot as an annotation.  This optional
    annotation leaves the plotted interval and the CRPS inset unchanged.
  }
}
\value{
  An \code{rfsrc} grow or prediction object with the additional class
  \code{quantreg}.  Existing classes are retained.  For a univariate
  regression response, the component \code{quantreg} contains:
  \item{quantiles}{Matrix with evaluated observations in rows and
    requested probability levels in columns.}
  \item{prob}{Sorted probability levels corresponding to \code{quantiles}.}
  \item{cdf}{CDF matrix evaluated at \code{yunq}.}
  \item{density}{Successive CDF increments, with the same dimensions
    as \code{cdf}.}
  \item{yunq}{Ordered training-response reporting grid.}
  \item{mean, std}{Direct distribution mean and standard deviation
    vectors for the residual methods.  \code{NULL} for GK, whose moment
    helper uses the reporting grid.}
  \item{method}{Canonical method name for this call: \code{"exchangeable"},
    \code{"forest"}, or \code{"gk"}.  Calls using the former
    \code{"local"} name return \code{"exchangeable"}.}
  \item{oob}{Whether the current mean predictions or native GK
    quantiles came from the OOB component.  This setting does not
    change the saved residual selection.}
  \item{prediction}{Current mean predictions used to shift the residuals.
    \code{NULL} for GK.}
  \item{residual.source}{The selected training residual distribution:
    \code{"oob"} or \code{"ensemble"} for the residual methods.
    \code{NULL} for GK, whose quantiles do not use residuals.}

  For multivariate and mixed-outcome forests, \code{quantreg} is a named
  list of these components for the returned continuous responses.
  Classification responses have no quantile component.

  Grow objects also provide numerical \code{residual} and
  \code{residual.oob} components, at the top level for univariate
  regression and within each continuous response's \code{regrOutput}
  component otherwise.  The saved forest retains the training residual
  vectors or named matrices and the residual-source settings described
  in Details.  The usual \code{xvar}, \code{yvar}, and other forest
  output components retain their meanings.

  \subsection{Helper functions}{
    \code{extract.quantile} always returns a response-named list of the
    quantile components described above.  The other helpers return the
    following for one continuous response when \code{pretty = TRUE},
    and a response-named list otherwise:
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
  \code{\link{rfsrc}}, \code{\link{predict.rfsrc}}
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
