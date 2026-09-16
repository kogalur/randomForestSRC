\name{classification.performance}
\alias{classification.performance}
\alias{get.confusion}
\alias{get.misclass.error}
\alias{get.auc}
\alias{get.brier.error}
\alias{get.logloss}
\alias{get.pr.auc}
\alias{get.pr.curve}
\alias{get.bayes.rule}
\title{Classification Performance Metrics}
\description{
  Evaluate classification predictions using confusion matrices,
  class-specific misclassification rates, ROC area under the curve (AUC),
  Brier score, log loss, and binary precision-recall summaries.  The
  helpers take observed responses and predicted classes or probabilities
  directly.  A class-assignment helper is also provided.
}
\usage{
get.confusion(y, class.or.prob)

get.misclass.error(y, yhat)

get.auc(y, prob)

get.brier.error(y, prob, normalized = TRUE, vector = FALSE)

get.logloss(y, prob, robust = TRUE)

get.pr.auc(truth, yhat)

get.pr.curve(truth, yhat)

get.bayes.rule(prob, class.relfrq = NULL)
}
\arguments{
  \item{y}{Observed class labels, normally a factor.  Keep the same
    factor levels and ordering as the probability columns.  For
    \code{get.auc()}, a nonfactor response is also accepted, with
    probability columns in the order of \code{sort(unique(y))}.
    \code{get.misclass.error()} also accepts vectors of class labels.}
  \item{prob}{Numeric matrix with one row per observation and one column
    per class.  Arrange columns in response-level order and name them
    with the corresponding class labels.  \code{get.auc()} uses column
    positions even when names are present.  \code{get.logloss()} and
    \code{get.bayes.rule()} require class column names;
    \code{get.brier.error()} can supply missing names from
    \code{levels(y)}.}
  \item{class.or.prob}{For \code{get.confusion()}, either a factor of
    predicted labels with the same levels as \code{y}, or a matrix of
    class probabilities as described for \code{prob}.  Missing matrix
    column names are supplied from \code{levels(y)}.}
  \item{yhat}{For \code{get.misclass.error()}, predicted class labels
    aligned with \code{y}.  For the precision-recall helpers, a numeric
    score vector, a one-column matrix or data frame, or a two-column
    matrix or data frame.  A vector or single column contains scores for
    the positive class; larger scores favor that class.  Scores need not
    lie in \eqn{[0,1]}.  See Details for the class and column conventions.}
  \item{truth}{Binary observed responses for precision-recall summaries:
    a vector coded \code{0}/\code{1}, or a factor with two levels.
    Code \code{1} is positive.  For factors with labels other than
    \code{0} and \code{1}, the less frequent class is positive.}
  \item{normalized}{If \code{TRUE}, use the normalized Brier scale
    described below, on which equal probabilities for all classes give
    score one.  If \code{FALSE}, average the squared probability errors
    over classes.}
  \item{vector}{If \code{TRUE}, return the scaled Brier contribution
    for each class, rather than their sum.  These are class-level
    contributions, not observation-level losses.}
  \item{robust}{For \code{get.logloss()}, remove infinite log losses
    before averaging.  The default is \code{TRUE}.  This does not clip
    probabilities or modify confusion counts.}
  \item{class.relfrq}{Optional vector of class relative frequencies in
    probability-column order, for a two-class problem.  With \code{NULL},
    \code{get.bayes.rule()} selects the largest class probability.
    Otherwise it applies the minority-prevalence cutoff described below.}
}
\details{
  \subsection{Supplying predictions}{
    These helpers score the predictions supplied to them; they do not fit
    a forest or automatically choose OOB predictions.  For a grow object,
    use \code{o$yvar} with \code{o$predicted.oob} for OOB probability
    scores, or \code{o$class.oob} for its stored class assignments.  For
    new-data evaluation, use the prediction object's responses and
    \code{predicted} or \code{class} components.

    Keep responses and predictions aligned by row.  For a common scoring
    sample, restrict to observed responses and finite probabilities before
    calling the helpers.  In particular, \code{get.logloss()} and
    \code{get.misclass.error()} require nonmissing observed labels.
    The examples illustrate this preparation for OOB predictions.

    For a combined binary summary with sensitivity, specificity, F1,
    G-mean, and random-reference comparisons, use
    \code{\link{get.imbalanced.performance}}.
  }
  \subsection{Class assignments and misclassification}{
    With \code{class.relfrq = NULL}, \code{get.bayes.rule()} assigns each
    row to a class with the largest probability, sampling among tied
    maximizers.  A row whose probabilities are all missing receives
    \code{NA}.  With two class relative frequencies supplied, it assigns
    the minority class when that class's probability is at least its
    supplied relative frequency.  Otherwise it assigns the majority
    class.  If the frequencies tie, the first probability column is
    treated as the minority.  Supply finite probabilities for this
    frequency-based rule.

    \code{get.confusion()} has observed classes in rows and predicted
    classes in columns.  Its final \code{class.error} column is one minus
    the fraction correctly classified within each observed class, rounded
    to four decimal places.  Missing response/prediction pairs do not
    enter the counts.  When supplied a probability matrix, it uses the
    largest-probability rule.  To evaluate RFQ or another thresholded
    rule, supply that rule's predicted factor instead, for example
    \code{o$class.oob}.

    \code{get.misclass.error()} returns one misclassification rate per
    observed class, in \code{sort(unique(y))} order.  It returns
    class-specific rates rather than an overall error rate.  Missing
    predicted labels make the corresponding class rate unavailable.
  }
  \subsection{ROC area under the curve}{
    \code{get.auc()} uses the pairwise multiclass construction of Hand
    and Till.  For each pair of classes, it calculates a rank-based AUC
    using each class's probability column as its score and averages the
    two directions.  It then averages over available class pairs.
    For complementary two-class probabilities, the two directions agree.
    Larger values indicate better discrimination.

    Missing observed labels and nonfinite scores are excluded from the
    relevant rank calculations.  Each direction requires at least two
    finite scores in each of the two classes.  Unavailable pairs are
    omitted; the result is \code{NA} when no pair is available.
  }
  \subsection{Brier score}{
    For \eqn{J} classes, let
    \deqn{b_j=\mathrm{mean}_i\left[(I(y_i=j)-p_{ij})^2\right].}
    The scalar returned by \code{get.brier.error()} is
    \deqn{\frac{J}{J-1}\sum_{j=1}^{J} b_j}
    with \code{normalized = TRUE}, and
    \deqn{\frac{1}{J}\sum_{j=1}^{J} b_j}
    with \code{normalized = FALSE}.  Smaller values are better.  Equal
    probabilities \eqn{p_{ij}=1/J} give normalized score one.  For
    complementary binary probabilities, the unnormalized score is the
    usual mean squared error of either class probability, and the
    normalized score is four times that value.

    With \code{vector = TRUE}, each returned element is \eqn{b_j}
    multiplied by the selected scaling constant.  Missing losses are
    omitted separately within each class.  The scalar sums available
    class contributions and is \code{NA} when none are available.
  }
  \subsection{Log loss}{
    With every response level represented, \code{get.logloss()} averages
    \eqn{-\log(p_{i,y_i})}, using the natural logarithm.  Smaller values
    are better.  With \code{robust = TRUE}, infinite losses from zero
    probabilities are excluded, rather than replaced using a probability
    floor.  Use \code{robust = FALSE} to retain these infinite losses.
    Missing losses are excluded in either case.  An unused factor level
    contributes one zero term to the helper's average.
  }
  \subsection{Precision-recall summaries}{
    These helpers apply to two classes.  For truth coded \code{0}/\code{1},
    including a factor with those two labels, \code{1} is positive.
    For other two-level factors, the less frequent observed class is
    positive; ties select the first factor level.  This choice is made
    before excluding missing scores.

    A score vector or single column always refers to the positive class.
    With two columns, named columns are matched to the original response
    levels; unnamed columns follow factor-level order, or \code{0},
    \code{1} order for a nonfactor response.  The positive-class column
    is then selected.  Rows with a missing response or nonfinite selected
    score are excluded.

    \code{get.pr.auc()} returns the precision-recall area and a random
    reference area equal to the positive-class proportion among the
    scored observations.  The area uses the analytic precision-recall
    interpolation implemented by the package.  \code{get.pr.curve()}
    supplies recall, precision, and threshold values for plotting; its
    rows run from high to low recall.  Both classes must have a usable
    score.  Larger precision-recall areas are better.
  }
}
\value{
  \describe{
    \item{\code{get.confusion}}{A numeric matrix of class counts with
      an additional \code{class.error} column.}
    \item{\code{get.misclass.error}}{An unnamed numeric vector of
      class-specific error rates, in \code{sort(unique(y))} order.}
    \item{\code{get.auc}}{A scalar AUC, or \code{NA} when unavailable.}
    \item{\code{get.brier.error}}{A scalar Brier score by default, or an
      unnamed length-\eqn{J} vector of scaled class contributions in
      probability-column order when \code{vector = TRUE}.}
    \item{\code{get.logloss}}{A scalar mean log loss.  It can be infinite
      with \code{robust = FALSE}, or \code{NaN} when no losses remain
      for averaging.}
    \item{\code{get.pr.auc}}{An unnamed numeric vector of length two:
      the model area followed by the random reference area.  Both
      entries are \code{NA} when the calculation is unavailable.}
    \item{\code{get.pr.curve}}{A numeric matrix with columns
      \code{recall}, \code{precision}, and \code{threshold}, or
      \code{NULL} when the calculation is unavailable.}
    \item{\code{get.bayes.rule}}{A factor of predicted classes with
      levels given by the probability-column names.}
  }
}
\seealso{
  \code{\link{rfsrc}}, \code{\link{predict.rfsrc}},
  \code{\link{imbalanced.rfsrc}},
  \code{\link{get.imbalanced.performance}},
  \code{\link{get.brier.survival}}, \code{\link{get.auct.survival}}
}
\examples{
## ------------------------------------------------------------
## A basic calculation from observed labels and probabilities
## ------------------------------------------------------------
y <- factor(c("no", "no", "no", "no", "yes", "yes"),
            levels = c("no", "yes"))
p <- c(.10, .20, .65, .35, .40, .85)
prob <- cbind(no = 1 - p, yes = p)
print(get.confusion(y, prob))
print(get.auc(y, prob))
print(get.brier.error(y, prob))

\donttest{
## ------------------------------------------------------------
## Class-specific errors and probability losses
## ------------------------------------------------------------
yhat <- get.bayes.rule(prob)
print(setNames(get.misclass.error(y, yhat), levels(y)))
print(get.brier.error(y, prob, normalized = FALSE))
print(setNames(get.brier.error(y, prob, vector = TRUE), colnames(prob)))
print(get.logloss(y, prob))

## Use supplied class frequencies for the binary RFQ decision rule.
class.frq <- as.numeric(prop.table(table(y)))
yhat.rfq <- get.bayes.rule(prob, class.relfrq = class.frq)
print(get.confusion(y, yhat.rfq))

## ------------------------------------------------------------
## Precision-recall from a score vector
## ------------------------------------------------------------
truth <- as.integer(y == "yes")
pr.auc <- get.pr.auc(truth, p)
print(setNames(pr.auc, c("model", "random")))
pr <- get.pr.curve(truth, p)
plot(pr[, "recall"], pr[, "precision"], type = "l",
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Recall", ylab = "Precision")
abline(h = pr.auc[2], lty = 2)

## A single score column is also accepted.
print(get.pr.auc(truth, matrix(p, ncol = 1)))

## ------------------------------------------------------------
## Multiclass forest: OOB and test-data performance
## ------------------------------------------------------------
set.seed(17)
train <- c(1:35, 51:85, 101:135)
o <- rfsrc(Species ~ ., data = iris[train, ], ntree = 100)

## Select one common set of observed responses and finite OOB scores.
p.oob <- o$predicted.oob
keep <- !is.na(o$yvar) & rowSums(!is.finite(p.oob)) == 0
print(get.confusion(o$yvar[keep], o$class.oob[keep]))
print(get.auc(o$yvar[keep], p.oob[keep, , drop = FALSE]))

## Use the current responses and probabilities for test-data scoring.
p.test <- predict(o, newdata = iris[-train, ])
print(c(
  auc = get.auc(p.test$yvar, p.test$predicted),
  brier = get.brier.error(p.test$yvar, p.test$predicted),
  logloss = get.logloss(p.test$yvar, p.test$predicted)
))
}
}
\keyword{classif}
\keyword{models}
