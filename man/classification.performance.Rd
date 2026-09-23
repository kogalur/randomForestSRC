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
  Evaluate predicted classes and probabilities using confusion matrices,
  class-specific misclassification rates, ROC area under the curve (AUC),
  Brier score, log loss, and binary precision-recall summaries.
  \code{get.bayes.rule()} assigns classes from predicted probabilities.
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
  \item{y}{Observed class labels, usually a factor with levels in
    probability-column order.  \code{get.auc()} also accepts nonfactor
    labels, with columns in \code{sort(unique(y))} order.
    \code{get.misclass.error()} also accepts vectors of class labels.}
  \item{prob}{Numeric matrix with one row per observation and one column
    per class, in response-level order.  Name columns with the class
    labels.  \code{get.auc()} uses column positions;
    \code{get.logloss()} and \code{get.bayes.rule()} require column names.
    \code{get.confusion()} and \code{get.brier.error()} supply missing
    names from \code{levels(y)}.}
  \item{class.or.prob}{A factor of predicted labels with the same levels
    as \code{y}, or a probability matrix as described for \code{prob}.}
  \item{yhat}{For \code{get.misclass.error()}, predicted class labels
    aligned with \code{y}.  For precision-recall summaries, a numeric
    score vector or a matrix or data frame with one or two columns.
    A vector or single column gives positive-class scores; larger values
    favor that class.  Scores need not lie in \eqn{[0,1]}.
    See Details for two-column inputs.}
  \item{truth}{Binary responses coded \code{0}/\code{1}, or a two-level
    factor.  For \code{0}/\code{1} labels, \code{1} is positive.
    For other factors, the less frequent class is positive, with ties
    resolved by the first factor level.  The positive class is selected
    before excluding missing scores.}
  \item{normalized}{If \code{TRUE}, scale the Brier score so that equal
    class probabilities give score one.  If \code{FALSE}, average the
    squared probability errors over classes.  See Details.}
  \item{vector}{If \code{TRUE}, return scaled Brier contributions for
    each class instead of their sum.}
  \item{robust}{For \code{get.logloss()}, exclude infinite losses before
    averaging.  Probabilities are unchanged.}
  \item{class.relfrq}{Optional two-class relative frequencies, in
    probability-column order.  \code{NULL} selects the largest-probability
    rule; supplying frequencies selects the minority-prevalence rule
    described below.}
}
\details{
  \subsection{Supplying predictions}{
    Supply responses and predictions in matching row order.  For OOB
    evaluation of a grow object \code{o}, use \code{o$yvar} with
    \code{o$predicted.oob} for probabilities or \code{o$class.oob} for
    predicted classes.  For new data, use the prediction object's
    \code{yvar} with \code{predicted} or \code{class}.

    To compare metrics on the same observations, retain observed labels
    and finite probabilities.  \code{get.logloss()} and
    \code{get.misclass.error()} require nonmissing observed labels.

    For a combined binary summary with sensitivity, specificity, F1,
    G-mean, and random-reference comparisons, use
    \code{\link{get.imbalanced.performance}}.
  }
  \subsection{Class assignments and misclassification}{
    With \code{class.relfrq = NULL}, \code{get.bayes.rule()} selects the
    class with the largest probability, breaking ties at random.  Rows
    with all probabilities missing receive \code{NA}.  With two class
    frequencies supplied, it assigns the minority class when its
    probability is at least its relative frequency, and the majority
    class otherwise.  If frequencies tie, the class in the first
    probability column is treated as the minority.  Supply finite
    probabilities for this rule.

    \code{get.confusion()} tabulates observed classes in rows and
    predicted classes in columns.  The final \code{class.error} column
    gives each observed class's misclassification rate, rounded to four
    decimal places.  Missing response/prediction pairs are excluded.
    Probability matrices use the largest-probability rule.  For RFQ or
    another thresholded rule, supply the predicted factor instead, such
    as \code{o$class.oob}.

    \code{get.misclass.error()} returns one error rate per observed
    class, in \code{sort(unique(y))} order.  A missing predicted label
    makes the corresponding class rate unavailable.
  }
  \subsection{ROC area under the curve}{
    \code{get.auc()} computes the multiclass AUC of Hand and Till (2001).
    For each class pair, it averages two rank-based AUCs, one using each
    class's probability as the score, then averages across available
    pairs.  For complementary binary probabilities, the two directions
    agree.  Larger values indicate better discrimination.

    Missing labels and nonfinite scores are excluded from each rank
    calculation.  Each direction requires at least two finite scores
    per class.  Unavailable pairs are omitted; the result is \code{NA}
    when none remain.
  }
  \subsection{Brier score}{
    The Brier score measures squared probability error (Brier, 1950).
    For \eqn{J} classes, let
    \deqn{b_j=\mathrm{mean}_i\left[(I(y_i=j)-p_{ij})^2\right].}
    \code{get.brier.error()} returns
    \deqn{\frac{J}{J-1}\sum_{j=1}^{J} b_j}
    for \code{normalized = TRUE}, or
    \deqn{\frac{1}{J}\sum_{j=1}^{J} b_j}
    for \code{normalized = FALSE}.  Smaller values are better.
    Equal probabilities \eqn{p_{ij}=1/J} give normalized score one.
    For complementary binary probabilities, the unnormalized score is
    the mean squared error of either class probability; the normalized
    score is four times that value.

    With \code{vector = TRUE}, each element is \eqn{b_j} multiplied by
    the selected scaling constant.  Missing losses are omitted within
    each class.  The scalar sums available contributions, returning
    \code{NA} when none remain.
  }
  \subsection{Log loss}{
    With every response level represented, \code{get.logloss()} averages
    the negative log probability of the observed class,
    \eqn{-\log(p_{i,y_i})}, using the natural logarithm (see Gneiting and
    Raftery, 2007).  Smaller values are better.  Missing losses are
    omitted.  With \code{robust = TRUE}, infinite losses from zero
    probabilities are also omitted; \code{robust = FALSE} retains them.
    Probabilities are not clipped.  An unused factor level contributes
    one zero term to the average.
  }
  \subsection{Precision-recall summaries}{
    For two-column scores, named columns match the original response
    levels; unnamed columns follow factor-level order,
    or \code{0}, \code{1} order for nonfactor responses.  The
    positive-class column is selected.  Rows with a missing response
    or nonfinite selected score are excluded; both classes must remain.

    \code{get.pr.auc()} returns the precision-recall area, calculated
    using analytic precision-recall interpolation, and a random reference
    area equal to the positive-class proportion among scored observations.
    \code{get.pr.curve()} returns recall, precision, and thresholds in
    decreasing recall order.  Larger precision-recall areas are better.
    See Davis and Goadrich (2006) for background on precision-recall
    curves.
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
\references{
  Brier, G.W. (1950). Verification of forecasts expressed in terms of
  probability. \emph{Monthly Weather Review}, 78, 1-3.
  \doi{10.1175/1520-0493(1950)078<0001:VOFEIT>2.0.CO;2}

  Davis, J. and Goadrich, M. (2006). The relationship between
  precision-recall and ROC curves. \emph{Proceedings of the 23rd
  International Conference on Machine Learning}, 233-240.
  \doi{10.1145/1143844.1143874}

  Gneiting, T. and Raftery, A.E. (2007). Strictly proper scoring rules,
  prediction, and estimation. \emph{Journal of the American Statistical
  Association}, 102, 359-378.
  \doi{10.1198/016214506000001437}

  Hand, D.J. and Till, R.J. (2001). A simple generalisation of the area
  under the ROC curve for multiple class classification problems.
  \emph{Machine Learning}, 45, 171-186.
  \doi{10.1023/A:1010920819831}
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
