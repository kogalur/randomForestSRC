\name{imbalanced.rfsrc}
\alias{imbalanced.rfsrc}
\alias{imbalanced}
\alias{get.imbalanced.performance}
\alias{get.imbalanced.optimize}
\alias{get.rfq.threshold}
\alias{print.imbalanced.performance}
\title{Random Forests for Imbalanced Two-Class Classification}

\description{
  Fits a random forest for a two-class outcome using the random forest
  quantile classifier (RFQ), balanced random forest (BRF), or standard
  classification.  RFQ changes the probability threshold used to assign
  class labels; BRF changes the sampling used to grow the trees.  The
  splitting rule, performance measure, and variable importance can be
  specified separately.

  Companion helpers summarize performance, select a probability
  threshold, or extract the class-frequency threshold from response
  labels.
}

\usage{\method{imbalanced}{rfsrc}(formula, data, ntree = 3000, 
  method = c("rfq", "brf", "standard"), splitrule = "auc",
  perf.type = NULL, block.size = NULL, fast = FALSE,
  ratio = NULL, ...)

get.imbalanced.performance(obj, prob = NULL, threshold = NULL,
                           confusion = FALSE, robust = FALSE)
get.imbalanced.optimize(obj, prob = NULL, newdata = NULL,
                        measure = c("gmean", "F1", "F1mod", "F1modgmean"),
                        ngrid = 1000, plot.it = TRUE)
get.rfq.threshold(y)

\method{print}{imbalanced.performance}(x, digits = 4,
                                      show.confusion = TRUE, ...)
}

\arguments{
  \item{formula}{A model formula with one factor response having two
    levels and one or more predictors.}

  \item{data}{A data frame containing the response and predictors.
    Data-frame subclasses are converted to a plain data frame.}

  \item{ntree}{Number of trees to grow.  The default is 3000.}

  \item{method}{Classification method.  \code{"rfq"}, the default, uses
    the RFQ decision rule of O'Brien and Ishwaran (2019).
    \code{"brf"} uses class-weighted sampling with replacement and a
    reduced tree sample size.  \code{"standard"} uses ordinary random
    forest classification without enabling the RFQ rule.  See Details.}

  \item{splitrule}{Splitting rule passed to the forest grow function.
    The default is \code{"auc"}.  Other choices include \code{"gini"}
    and \code{"entropy"}; see \code{\link{rfsrc}} for supported
    classification rules.  This setting is separate from
    \code{perf.type}.}

  \item{perf.type}{Performance measure used by the forest and its
    variable importance calculations.  The default is \code{"gmean"}
    for RFQ and BRF, and \code{"default"} (misclassification error) for
    standard classification.  Choices accepted by this wrapper are
    \code{"none"}, \code{"default"}, \code{"standard"},
    \code{"misclass"}, \code{"brier"}, \code{"gmean"}, and the legacy
    spelling \code{"g.mean"}.  The latter is converted to
    \code{"gmean"}.  Use \code{"none"} to omit forest performance
    calculations.}

  \item{block.size}{Block size passed to the selected forest grow
    function for performance and VIMP calculations; see
    \code{\link{rfsrc}}.  This wrapper passes \code{NULL} explicitly
    when the argument is not set.  Specify the same block size when
    comparing VIMP across methods.}

  \item{fast}{Logical.  Use \code{\link{rfsrc.fast}} for the RFQ or
    standard method when \code{TRUE}.  The BRF branch always calls
    \code{rfsrc} and does not use this argument.}

  \item{ratio}{Optional, experimental majority-class undersampling for
    the RFQ and standard methods.  Supply a number between zero and one.
    The majority class is sampled without replacement, while the
    minority class is bootstrapped with replacement.  Small requested
    fractions may be increased by the sampling helper.  This is a
    majority-class sampling fraction, rather than a requested final
    minority prevalence.  Ignored for BRF.  See Details for OOB support.}

  \item{obj}{For the performance and optimization helpers, a grow or
    prediction object with a single two-class outcome and observed
    responses, or a two-level response factor when \code{prob} is
    supplied.}

  \item{prob}{Optional probabilities aligned with the responses in
    \code{obj}.  A numeric vector or one-column matrix or data frame
    contains minority-class probabilities.  A two-column matrix or data
    frame contains both class probabilities; named columns are matched
    to factor levels, and unnamed columns follow factor-level order.
    \code{NULL} uses the object's predictions as described in Details.}

  \item{threshold}{Cutoff for minority-class probabilities in
    \code{get.imbalanced.performance}.  \code{NULL}, the default, uses
    the minority proportion among the observed evaluation responses.
    \code{"auto"} uses that prevalence for an RFQ object and \eqn{0.5}
    otherwise.  \code{"bayes"} always uses \eqn{0.5}.  A single finite
    numeric value specifies the cutoff directly; values outside
    \eqn{[0,1]} are clipped to that interval.  With direct response and
    probability inputs, \code{"auto"} uses \eqn{0.5} because no fitted
    classifier is supplied.}

  \item{confusion}{Return a list including the confusion matrix instead
    of the default named vector of performance measures?}

  \item{robust}{Add one to the diagonal counts when calculating the
    four classification rates and their derived summaries?  See Details.}

  \item{newdata}{Optional threshold-selection data for
    \code{get.imbalanced.optimize}, supplied with a grow object and
    \code{prob = NULL}.  Must contain the predictors and observed
    response.  \code{NULL} uses predictions already in \code{obj}.
    Ignored when \code{prob} is supplied.}

  \item{measure}{Performance measure to maximize during threshold
    selection: \code{"gmean"} (the default), \code{"F1"},
    \code{"F1mod"}, or \code{"F1modgmean"}.  See Details for definitions.}

  \item{ngrid}{Number of equally spaced thresholds from zero to one,
    including both endpoints.  Must be an integer of at least two.}

  \item{plot.it}{Draw the four threshold-search diagnostic panels?}

  \item{y}{Class-label vector for \code{get.rfq.threshold}, typically a
    two-level factor.  Missing responses are omitted when counting classes.}

  \item{...}{Additional arguments passed to \code{\link{rfsrc}} or,
    when selected, \code{\link{rfsrc.fast}}.  Examples include
    \code{mtry}, \code{nodesize}, \code{nsplit}, \code{importance}, and
    \code{na.action}.  BRF sets its own \code{case.wt}, \code{sampsize},
    and \code{samptype}.  When \code{ratio} is supplied, the wrapper
    constructs \code{samp} and sets \code{bootstrap = "by.user"}.
    Additional arguments to the print method are currently unused.}
  \item{x}{For the print method, the result of
    \code{get.imbalanced.performance}.}
  \item{digits}{Number of decimal places used in the performance display.}
  \item{show.confusion}{Print model and expected random-reference
    confusion tables when a confusion matrix is present in \code{x}?
    Obtain that matrix with \code{confusion = TRUE} in
    \code{get.imbalanced.performance}. This option does not recompute it.}
}

\details{
  \subsection{Choice of classifier}{
    Let \eqn{\pi} be the training proportion of the minority class and
    \eqn{\widehat p(x)}{p-hat(x)} its predicted probability.  RFQ assigns
    an observation to the minority class when
    \eqn{\widehat p(x) \geq \pi}{p-hat(x) >= pi}.  The class-frequency
    threshold can identify minority observations whose predicted
    probability is below one half.  The standard method instead uses
    ordinary class assignment based on the largest predicted probability.
    RFQ here is a classification rule, distinct from the continuous-response
    quantile estimation performed by \code{\link{quantreg}}.

    BRF calls the grow function with replacement sampling, a tree sample
    size of twice the smaller observed class count, and case weights that
    give the two classes equal total sampling weight.  The wrapper does
    not construct separate fixed-size samples from each class.  Both
    \code{fast} and \code{ratio} are ignored in this branch.
  }
  \subsection{Performance and variable importance}{
    G-mean is the geometric mean of minority sensitivity and majority
    specificity, \eqn{\sqrt{\mathrm{sensitivity}\,
    \mathrm{specificity}}}{sqrt(sensitivity * specificity)}.  Larger
    values indicate better balance between the two class-specific rates.
    Misclassification and Brier scores are error measures for which smaller
    values are better.  Changing \code{perf.type} selects the forest
    performance measure; it does not select the RFQ decision rule.

    Request permutation VIMP explicitly with
    \code{importance = "permute"}.  For comparisons, use a common
    \code{perf.type}, splitting rule, and \code{block.size}, unless the
    comparison is intended to study those settings.  Increase \code{ntree}
    as needed to assess the stability of OOB performance and VIMP.
  }
  \subsection{Sampling and missing values}{
    With \code{ratio = NULL}, the RFQ and standard branches pass
    missing-data handling to the selected grow function.  BRF and calls
    using \code{ratio} first apply complete-case omission to the supplied
    data frame.  This omission currently includes columns that are not
    selected by the formula and occurs before the grow function handles
    \code{na.action}.

    With \code{ratio} supplied, each tree uses a bootstrap sample of the
    minority observations and a without-replacement sample of majority
    observations.  At \code{ratio = 1}, every majority observation is
    in-bag for every tree, so majority-class OOB performance is unavailable.
    A single minority observation is also always in-bag in this sampling
    scheme.  Use a separate evaluation sample when OOB support is absent.
  }
  \subsection{Detailed performance summaries}{
    \code{\link{get.imbalanced.performance}} accepts a fitted or predicted
    forest object, or a two-level factor supplied as \code{obj} together
    with \code{prob}.  For a forest object, it prefers
    \code{predicted.oob} when that component exists and is not entirely
    missing; otherwise it uses \code{predicted}.  A partially missing OOB
    component is not filled with full-ensemble predictions.

    A numeric vector or one-column matrix or data frame represents
    minority-class probabilities, with one value per response.  No column
    name is needed for these inputs.  A two-column matrix or data frame
    contains both class probabilities.  Named columns are matched to the
    response factor levels; unnamed columns follow factor-level order.
    The response and probability rows must already be aligned.

    Finite probabilities are restricted to \eqn{[0,1]} using the existing
    clipping rule.  Finite two-column rows must sum to one within numerical
    tolerance.  Rows with a missing response or a nonfinite probability are
    excluded from scoring; both columns must be finite for a two-column
    input.  Unavailable scores are returned as \code{NA}.  Class counts,
    minority identity, and the default threshold are determined from all
    observed evaluation responses before excluding unavailable predictions.

    The helper determines the minority class and its proportion from the
    response supplied for evaluation.  If the class counts tie, the first
    factor level is selected.  The default \code{threshold = NULL} uses
    this evaluation proportion for all forest methods.

    With \code{threshold = "auto"}, the saved \code{forest$rfq} flag
    selects the cutoff: \code{TRUE} uses the same evaluation prevalence
    as \code{NULL}, and otherwise the cutoff is \eqn{0.5}.  This applies
    to both grow and prediction objects.  For RFQ test predictions the
    prevalence is therefore calculated from the evaluation responses,
    not substituted from the training responses.  Splitting rules and
    \code{perf.type} do not determine whether the classifier is RFQ.
    Direct response/probability inputs have no fitted-classifier flag,
    so \code{"auto"} uses \eqn{0.5}.  Use \code{threshold = "bayes"}
    to request \eqn{0.5} explicitly for any classifier.

    All settings threshold the probabilities rather than reading the
    stored \code{class} component.  Supply a single finite numeric
    \code{threshold} to examine a different cutoff; finite values outside
    \eqn{[0,1]} are clipped to that interval.  Threshold ties are assigned
    to the minority class.  A supplied threshold changes the cutoff, not
    how the minority class is identified or how the forest was fitted.

    Besides the four rates \code{sens}, \code{spec}, \code{prec}, and
    \code{npv}, the summary includes \code{misclass}, \code{brier},
    \code{brier.norm}, \code{auc}, \code{logloss}, and \code{pr.auc}.
    For complementary binary probabilities, \code{brier} is the mean
    squared error of the minority probability and \code{brier.norm} is
    four times this value.  The current log-loss helper excludes infinite
    losses rather than clipping zero probabilities.

    \code{F1} is the harmonic mean of precision and sensitivity;
    \code{F1mod} is the harmonic mean of sensitivity, specificity,
    precision, and negative predictive value.  \code{F1gmean} averages
    F1 and G-mean; \code{F1modgmean} averages F1mod and G-mean.  A
    harmonic mean is zero when all its component rates are defined and at
    least one is zero.  It is \code{NA} when any required rate is undefined.

    \code{robust = TRUE} adds one to the two diagonal counts when
    calculating the four rates and their derived summaries.  It does not
    change the returned confusion counts, misclassification error, or
    probability-based scores.  Use \code{confusion = TRUE} to return the
    confusion matrix in addition to the metrics.  Its rows are observed
    classes and columns are predicted classes, with \code{0} denoting the
    majority and \code{1} the minority.
  }
  \subsection{Random-score reference}{
    The performance summary has a \code{"rand"} attribute containing
    reference metrics and expected confusion counts.  This reference uses
    a score \eqn{U} uniformly distributed on \eqn{[0,1]}, independent of
    class, and predicts the minority class when \eqn{U} reaches the chosen
    threshold.  For threshold \eqn{t}, its sensitivity is \eqn{1-t},
    specificity is \eqn{t}, and ROC AUC is \eqn{1/2}.  The separately
    returned PR-AUC reference is the proportion of positive observations
    among the rows scored by the PR helper.  Its unnormalized Brier, normalized Brier,
    and log-loss references are \eqn{1/3}, \eqn{4/3}, and \eqn{1},
    respectively.  These differ from a predictor that always returns
    probability \eqn{1/2}.

    Nonlinear rate summaries are computed from the reference rates (or
    smoothed expected counts under \code{robust = TRUE}); they are not
    finite-sample expectations of those nonlinear statistics.
    The performance print method compares the model and reference with
    \code{Delta} and percentage \code{Gain}, oriented so that positive
    values indicate improvement.  The fitted-forest printer has its own
    probability reference; it need not use the same baseline.
  }
  \subsection{Threshold selection}{
    \code{get.imbalanced.optimize} searches \code{ngrid} equally spaced
    thresholds from zero to one, using G-mean, F1, F1mod, or F1modgmean as
    selected by \code{measure}.  It uses the first grid point attaining
    the largest available value.  It returns one row of performance
    metrics at the selected threshold and does not change the forest or
    its stored class predictions.  \code{plot.it = FALSE} suppresses its
    four diagnostic panels.  \code{ngrid} must be an integer of at least
    two.  When the chosen measure is unavailable at every grid point, the
    function stops with an explanatory error.  A valid optimum of zero is
    retained.  The caller's graphics settings are restored after plotting.

    A grow object ordinarily supplies OOB probabilities for this search.
    Supplying \code{newdata} searches using that data's observed responses
    and predictions, making it a threshold-selection sample.  For an
    independent assessment of the chosen threshold, apply it unchanged
    to a separate test sample.  The examples select a threshold from
    training OOB predictions and carry it forward to the test data.
  }
  \subsection{Printing performance}{
    Printing a performance result displays two metric tables. The first
    contains classification rates and their F1 and G-mean combinations;
    the second contains overall performance, including G-mean, AUC,
    PR-AUC, misclassification, Brier scores, and log loss. The same
    G-mean row appears in both tables; F1-related combinations appear
    only in the first.

    \code{Delta} and percentage \code{Gain} compare the model with its
    random-score reference, with positive values indicating improvement.
    Larger rates, G-mean, and AUC values are better; smaller errors and
    losses are better. \code{digits} controls the display precision,
    and \code{show.confusion = FALSE} omits the optional confusion tables.
    These settings do not change the stored performance values.
  }
}

\value{
  \code{imbalanced} returns an object from the selected forest grow
  function.  Usual components include the class-probability matrices
  \code{predicted}
  and \code{predicted.oob}, class assignments \code{class} and
  \code{class.oob}, and requested performance and importance values.
  Availability depends on the selected grow interface and sampling
  options; see \code{\link{rfsrc}} and \code{\link{rfsrc.fast}}.

  \subsection{Performance and threshold helpers}{
    \describe{
      \item{\code{get.imbalanced.performance}}{A named numeric vector of
        class counts, imbalance ratio, resolved numeric threshold, and
        performance measures.  With \code{confusion = TRUE}, a named list also
        includes \code{confusion}, the confusion matrix with a
        \code{class.error} column.  Both forms have class
        \code{"imbalanced.performance"} and a \code{"rand"} attribute
        containing the random-reference metrics and expected confusion
        counts.}
      \item{\code{get.imbalanced.optimize}}{A one-row data frame of
        performance measures at the selected threshold.  Extract the
        cutoff from the \code{threshold} column; the row name is its
        grid index.}
      \item{\code{get.rfq.threshold}}{The smaller class count divided by
        the total observed class count, as a numeric scalar.  Returns
        \code{NULL} unless the class-frequency table has exactly two
        entries.  No forest or probability predictions are required.}
    }
    The performance and optimization helpers return \code{NULL} when
    the evaluation response is not a two-level factor.
  }

  \code{print.imbalanced.performance} displays the summaries and returns
  \code{x} invisibly.
}

\author{Hemant Ishwaran and Udaya B. Kogalur}

\references{
  Chen, C., Liaw, A. and Breiman, L. (2004). Using random forest to learn
  imbalanced data. University of California, Berkeley, Technical Report
  110.

  Kubat, M., Holte, R. and Matwin, S. (1997). Learning when negative
  examples abound. \emph{Machine Learning}, ECML-97: 146-153.

  O'Brien, R. and Ishwaran, H. (2019). A random forests quantile
  classifier for class imbalanced data. \emph{Pattern Recognition},
  90, 232-249.
}

\seealso{
  \code{\link{rfsrc}}, \code{\link{rfsrc.fast}},
  \code{\link{get.imbalanced.performance}}, \code{\link{subsample}}
}

\examples{
\donttest{
## ------------------------------------------------------------
## RFQ, BRF, and standard classification
## ------------------------------------------------------------
data(breast, package = "randomForestSRC")
dta <- na.omit(breast)
f <- status ~ .

set.seed(19)
rfq <- imbalanced(f, dta, ntree = 300)
print(rfq)
print(get.imbalanced.performance(rfq, confusion = TRUE))

## Extract the class-frequency threshold from the observed responses.
print(get.rfq.threshold(rfq$yvar))

brf <- imbalanced(f, dta, ntree = 300, method = "brf")
standard <- imbalanced(f, dta, ntree = 300, method = "standard")

## Compare all three using the helper's default evaluation rule.
metrics <- c("sens", "spec", "gmean", "misclass", "brier", "auc")
comparison <- rbind(
  RFQ = as.numeric(get.imbalanced.performance(rfq)[metrics]),
  BRF = as.numeric(get.imbalanced.performance(brf)[metrics]),
  Standard = as.numeric(get.imbalanced.performance(standard)[metrics])
)
colnames(comparison) <- metrics
print(comparison)

## Supply a cutoff explicitly to inspect another probability rule.
print(get.imbalanced.performance(standard, threshold = .5))

## Let the saved classifier setting select prevalence or 0.5.
print(get.imbalanced.performance(rfq, threshold = "auto"))
print(get.imbalanced.performance(standard, threshold = "auto"))
print(get.imbalanced.performance(brf, threshold = "auto"))

## Request the 0.5 Bayes cutoff explicitly, including for RFQ.
print(get.imbalanced.performance(rfq, threshold = "bayes"))

## Vectors, one-column objects, and named two-column matrices.
class.frequency <- table(rfq$yvar)
minority <- names(class.frequency)[which.min(class.frequency)]
p <- rfq$predicted.oob[, minority]
print(get.imbalanced.performance(rfq$yvar, p))
print(get.imbalanced.performance(rfq$yvar, data.frame(p = p)))
print(get.imbalanced.performance(
  rfq$yvar, rfq$predicted.oob[, rev(levels(rfq$yvar)), drop = FALSE]
))

## ------------------------------------------------------------
## Train/test assessment with a training-selected threshold
## ------------------------------------------------------------
set.seed(23)
trn.id <- sort(unlist(lapply(
  split(seq_len(nrow(dta)), dta$status),
  function(idx) idx[sample.int(length(idx),
    size = floor(2 * length(idx) / 3))]
)))
trn <- dta[trn.id, , drop = FALSE]
tst <- dta[-trn.id, , drop = FALSE]
print(rbind(training = table(trn$status), test = table(tst$status)))

fit <- imbalanced(f, trn, ntree = 300)
selected <- get.imbalanced.optimize(fit, measure = "gmean",
                                     ngrid = 501, plot.it = FALSE)
print(selected)
threshold.selected <- selected$threshold

pred <- predict(fit, newdata = tst)
print(get.imbalanced.performance(pred))
print(get.imbalanced.performance(pred, threshold = threshold.selected))

## Predicting labels for new rows does not require their outcomes.
class.frequency <- table(trn$status)
minority <- names(class.frequency)[which.min(class.frequency)]
majority <- setdiff(levels(trn$status), minority)
pred.x <- predict(fit, newdata = tst[, fit$xvar.names, drop = FALSE])
selected.class <- factor(
  ifelse(pred.x$predicted[, minority] >= threshold.selected,
         minority, majority),
  levels = levels(trn$status)
)
print(head(selected.class))

## ------------------------------------------------------------
## Explicit permutation VIMP and blocking
## ------------------------------------------------------------
rfq.vimp <- imbalanced(f, dta, ntree = 300,
                       importance = "permute", block.size = 1)
rfq.block <- imbalanced(f, dta, ntree = 300,
                        importance = "permute", block.size = 10)
print(cbind(unblocked = rfq.vimp$importance[, 1],
            blocked = rfq.block$importance[, 1]))

## Small B is for illustration; increase it for final inference.
smp <- subsample(rfq.block, B = 25, verbose = FALSE)
print(extract.subsample(smp)$var.sel.Z)
plot.subsample(smp, cex.axis = .7, alpha = .5)

## ------------------------------------------------------------
## Optional sampling and grow-interface choices
## ------------------------------------------------------------
rfq.sample <- imbalanced(f, dta, ntree = 300, ratio = .25)
print(get.imbalanced.performance(rfq.sample))
rfq.fast <- imbalanced(f, dta, ntree = 300, fast = TRUE)
print(get.imbalanced.performance(rfq.fast))
}
}
\keyword{classif}
\keyword{tree}
