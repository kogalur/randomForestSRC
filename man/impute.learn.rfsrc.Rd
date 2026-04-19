\name{impute.learn.rfsrc}
\alias{impute.learn.rfsrc}
\alias{impute.learn}
\alias{save.impute.learn.rfsrc}
\alias{save.impute.learn}
\alias{load.impute.learn.rfsrc}
\alias{load.impute.learn}
\alias{predict.impute.learn.rfsrc}
\alias{predict.impute.learn}
\alias{impute.ood.rfsrc}
\alias{impute.ood}
\alias{print.impute.learn.rfsrc}
\alias{print.impute.learn}
\title{Learn a predictive imputer for test-time imputation and OOD scoring}
\description{
  Learns a predictive imputer from training data for later use on new data.

  If the training data contain missing values, the function first
  imputes them using \code{impute}. It then fits one saved full-sweep
  learner per selected target on the completed training data and reuses
  those learners later to update missing values in new data without
  refitting on the test set.

  The same saved learner bank can also be used to score new data for
  out-of-distribution (OOD) behavior. Note that OOD scores are available
  even when new data have missing values.  Each selected target is
  reconstructed from its saved conditional learner and compared with the
  observed value.  Target-wise discrepancies are calibrated against a
  training reference calculated from out-of-bag predictions computed
  during the training.

  If the training data are complete and \code{target.mode = "all"},
  the initial training-data imputation step is skipped and the
  full-sweep learners are fit directly from the complete training data.
}
\usage{
impute.learn.rfsrc(formula, data,
  ntree = 100, nodesize = 1, nsplit = 10,
  nimpute = 2, fast = FALSE, blocks,
  mf.q, max.iter = 10, eps = 0.01,
  ytry = NULL, always.use = NULL, verbose = TRUE,
  ...,
  full.sweep.options = list(ntree = 100, nsplit = 10),
  target.mode = c("missing.only", "all"),
  deployment.xvars = NULL,
  anonymous = TRUE,
  learner.prefix = "impute.learner.",
  learner.root = "learners",
  out.dir = NULL,
  wipe = TRUE,
  keep.models = is.null(out.dir),
  keep.ximp = FALSE,
  save.on.fit = !is.null(out.dir),
  save.ood = TRUE,
  weight = NULL)

save.impute.learn.rfsrc(object, path, wipe = TRUE, verbose = TRUE)

load.impute.learn.rfsrc(path, targets = NULL, lazy = TRUE, verbose = TRUE)

\method{predict}{impute.learn.rfsrc}(object, newdata,
  max.predict.iter = 3L,
  eps = 1e-3,
  targets = NULL,
  restore.integer = TRUE,
  cache.learners = c("session", "none", "all"),
  verbose = TRUE,
  ...)

impute.ood.rfsrc(object, newdata,
  targets = NULL,
  max.predict.iter = 3L,
  eps = 1e-3,
  cache.learners = c("all", "session", "none"),
  weight = NULL,
  aggregate = c("bounded.product",  "weighted.mean",
    "weighted.lp", "weighted.lp.log", "top.k"),
  aggregate.args = list(),
  return.details = FALSE,
  verbose = TRUE,
  ...)
}
\arguments{
  \item{formula}{A symbolic model description. Can be omitted. The same
    interpretation as in \code{impute} is used for the initial
    training-data imputation stage. The saved full-sweep learner bank is
    controlled by \code{deployment.xvars}, not by \code{formula}.}

  \item{data}{Training data. Variables that are not real-valued are
    coerced to factors before fitting when possible; otherwise fitting
    stops with an error. Rows and columns that are entirely missing are
    dropped before the training begins.}

  \item{ntree, nodesize, nsplit, nimpute, fast, blocks,
    max.iter, ytry, always.use, verbose}{Arguments passed to
    \code{impute} for the initial training-data imputation. The
    argument \code{full.sweep} is controlled internally and should not
    be supplied. The same \code{verbose}, \code{max.predict.iter},
    \code{eps}, and \code{cache.learners} controls are also used by
    \code{predict.impute.learn} and \code{impute.ood}.}

  \item{mf.q}{Controls the imputation engine used by \code{impute}. If
    \code{mf.q = 1}, training uses standard \code{missForest}. If
    \code{mf.q > 1}, training uses the multivariate \code{missForest}
    generalization. If \code{mf.q} is omitted, the training imputation
    follows the default behavior of \code{impute}.}

  \item{eps}{Convergence threshold. In \code{impute.learn} this
    controls the initial training-data imputation. In
    \code{predict.impute.learn} and \code{impute.ood} it controls early
    stopping for the prediction-time sweep.}

  \item{\dots}{For \code{impute.learn}, additional arguments passed to
    \code{impute}. For \code{predict.impute.learn} and
    \code{impute.ood}, additional arguments are currently ignored.}

  \item{full.sweep.options}{A \code{list} of options used when fitting
    the full sweep after the training data have been imputed.
    Recognized entries include \code{ntree}, \code{nodesize},
    \code{nsplit}, \code{mtry}, \code{splitrule}, \code{bootstrap},
    \code{sampsize}, \code{samptype}, \code{perf.type}, \code{rfq},
    \code{save.memory}, \code{importance}, and \code{proximity}.
    Unknown entries are ignored with a warning.}

  \item{target.mode}{Determines which variables receive a saved
    full-sweep learner. The default \code{"missing.only"} saves
    learners only for variables that were missing in the training data.
    The option \code{"all"} saves a learner for every variable. If the
    training data are complete, \code{target.mode = "all"} must be
    used. For the broadest OOD coverage, \code{target.mode = "all"} is
    recommended so every deployment-time variable can be reconstructed.}

  \item{deployment.xvars}{Controls which predictors are assumed to be
    available later when the saved imputer is used on new data. If
    \code{NULL}, all columns except the target are used. If a character
    vector, the same predictor set is used for all targets. If a named
    list, each target can have its own predictor set. By default, all
    non-target columns are eligible predictors, so users should exclude
    outcomes, future information, identifiers, or any variables that
    will not be available at deployment time.}

  \item{anonymous}{If \code{TRUE}, uses \code{rfsrc.anonymous} when
    fitting the full sweep. This usually reduces the size of the saved
    object.}

  \item{learner.prefix, learner.root}{Names used when writing saved
    full-sweep learners to disk.}

  \item{out.dir}{Optional output directory. If supplied and
    \code{save.on.fit = TRUE}, the manifest and the saved full-sweep
    learners are written to this directory during fitting. This requires
    the \pkg{fst} package because learners are serialized with
    \code{fast.save}.}

  \item{wipe}{If \code{TRUE}, removes an existing output directory
    before writing a new one.}

  \item{keep.models}{If \code{TRUE}, keeps the fitted full-sweep
    learners in memory in the returned object. At least one storage mode
    must be enabled: either \code{keep.models = TRUE} or
    \code{out.dir} with \code{save.on.fit = TRUE}.}

  \item{keep.ximp}{If \code{TRUE}, keeps the completed training data in
    the returned object. This is not required for later prediction.}

  \item{save.on.fit}{If \code{TRUE} and \code{out.dir} is supplied,
    writes the imputer to disk during fitting.}

  \item{save.ood}{If \code{TRUE}, computes and stores an OOD reference.
    The reference is built during training from target-wise out-of-bag
    reconstruction discrepancies, their target-wise calibrated training
    scores, and a default row-level weighted mean using the saved OOD
    target weights. If no \code{weight} is supplied at fit time, equal
    target weights are used. The saved target-wise training-score matrix
    allows \code{impute.ood} to rebuild a calibrated row-level
    percentile for arbitrary target subsets, test-time weight overrides,
    and alternate row aggregates besides the weighted mean.}

  \item{object}{An object returned by \code{impute.learn} or
    \code{load.impute.learn}.}

  \item{path}{Directory containing a saved imputer. Save and load
    operations require the \pkg{fst} package because learners are read
    and written with \code{fast.save} and \code{fast.load}.}

  \item{targets}{Optional subset of target variables to load, update,
    or score. Unknown names are ignored with a warning. For
    \code{impute.ood}, row-level percentile calibration is rebuilt for
    the requested target subset from the saved target-wise training OOD
    scores whenever those scores are available in the manifest.}

  \item{lazy}{If \code{TRUE}, saved learners are loaded only when they
    are needed. If \code{FALSE}, all saved learners are loaded at once.}

  \item{newdata}{New data to be imputed or scored. Missing columns are
    added and extra columns are dropped to match the training schema.
    Unseen factor levels are converted to \code{NA} for harmonization,
    but they are also tracked row-wise. In \code{impute.ood}, any row
    containing an unseen factor level is flagged and its row-level OOD
    score is set to the maximum value.}

  \item{max.predict.iter}{Maximum number of full-sweep passes applied to
    \code{newdata} before the saved learner bank is used for OOD
    reconstruction or returned prediction-time imputations.}

  \item{restore.integer}{If \code{TRUE}, integer columns in the
    returned data are rounded and restored as integers. Factor columns
    are always conformed back to the training schema. The package
    operates on real-valued and factor variables; inputs that are not
    real-valued are coerced to factors during preprocessing when
    possible, otherwise an error is raised.}

  \item{cache.learners}{How saved learners are reused during
    prediction or OOD scoring. For \code{predict.impute.learn}, the
    default \code{"session"} loads each needed learner once per call.
    The option \code{"none"} reloads a learner every time it is needed.
    The option \code{"all"} loads all requested learners before work
    starts. For \code{impute.ood}, \code{"all"} is the default because
    the saved learner bank is typically reused once for predictor-side
    completion and again for target reconstruction.}

  \item{weight}{Optional nonnegative target weights used for row-level OOD
    aggregation. In \code{impute.learn}, these weights define the
    default row-level OOD weighting scheme stored in the manifest. In
    \code{impute.ood}, they define the active row-level weighting scheme
    for the current scoring call. If supplied as a named vector, entries
    are matched to targets by name; omitted targets are set to zero, and
    extra names are ignored. If omitted in \code{impute.ood}, the saved
    training-time OOD weights are used automatically. Because the fit
    stores target-wise training OOD scores, \code{score.percentile} can
    be recalibrated automatically for test-time weight overrides rather
    than being limited to the original training-time weights.}

  \item{aggregate}{Row-level aggregation metric used by
    \code{impute.ood} to combine calibrated target-wise OOD scores.  The
    default \code{"bounded.product"} applies a weighted product of
    the form \eqn{1-\prod_j (1-u_j+\varepsilon)^{\tilde w_j}} where
    \eqn{u_j} are row level calibrated scores for target feature
    \eqn{j}.  \code{"weighted.mean"} is the weighted
    average. \code{"weighted.lp"} applies a weighted Minkowski \eqn{L_p}
    aggregation to the calibrated target scores.
    \code{"weighted.lp.log"} first applies the tail-stretching transform
    \eqn{-\log(1-u_j+\varepsilon)} to each calibrated target score
    \eqn{u_j} and then applies the weighted \eqn{L_p}
    aggregation. \code{"top.k"} averages only the \eqn{k} largest
    scoreable target scores among the positively weighted targets.  When
    the fit stores target-wise training OOD scores,
    \code{score.percentile} is rebuilt for the requested aggregate as
    well as the requested targets and weights.}

  \item{aggregate.args}{Optional list of tuning arguments for
    \code{aggregate}. Recognized entries are \code{p} for
    \code{"weighted.lp"} and \code{"weighted.lp.log"}, \code{k}
    (or \code{top.k}) for \code{"top.k"}, and \code{eps} for
    \code{"weighted.lp.log"} and \code{"bounded.product"}. The
    default values are \code{p = 2}, \code{k = 1}, and
    \code{eps = 1e-12}. Unknown entries are ignored with a warning.}

  \item{return.details}{If \code{TRUE}, \code{impute.ood} returns the
    per-target discrepancy and calibrated target-score matrices and a
    row-by-variable unseen-level mask in the diagnostics.}
}
\details{
  A predictive imputer is calculated in two stages.

  The training data are first normalized to a data frame. Variables
  that are not real-valued are coerced to factors when possible;
  otherwise fitting stops with an error. Rows and columns that are
  entirely missing are removed before the training schema is stored.

  If the resulting training data contain missing values, the first
  stage uses \code{impute} to complete the training data. The
  imputation engine is chosen in exactly the same way as for
  \code{impute} itself. In particular, \code{mf.q = 1} gives standard
  \code{missForest}, \code{mf.q > 1} gives the multivariate
  \code{missForest} generalization, and if \code{mf.q} is omitted the
  default \code{impute} behavior is used. If the training data are
  already complete and \code{target.mode = "all"}, this initial
  imputation step is skipped.

  In the second stage, a full sweep is fit on the completed training
  data. For each target selected by \code{target.mode}, rows where that
  target was observed are used to fit a forest with that target on the
  left-hand side and the selected deployment predictors on the
  right-hand side. The saved learner bank therefore depends on
  \code{deployment.xvars}. The \code{formula} argument affects the
  initial training-data imputation step, but it does not define the
  saved predictor bank for the later test-time sweep.

  By default, \code{deployment.xvars = NULL} allows every non-target
  column to be used as a predictor. This is convenient, but it can also
  introduce leakage if the training data include outcomes, future-only
  variables, identifiers, or any fields that will not be available when
  the learned imputer is applied to new data. Restrict
  \code{deployment.xvars} when that is a concern.

  When the imputer is saved to disk, each full-sweep learner is written
  separately using \code{fast.save}. Loading uses \code{fast.load}. In
  practice this gives a small manifest plus a directory of saved
  learners. The \pkg{fst} package is therefore required for save and
  load operations. The explicit save method can write learners either
  from memory or by reloading them from an attached saved path.

  Prediction starts by matching \code{newdata} to the training schema,
  filling missing values with training means or modes, and then
  applying one or more full-sweep passes. Only the targets selected by
  \code{target.mode} are updated by saved learners.

  If \code{target.mode = "missing.only"}, a variable that was complete
  in training but missing in new data is initialized from the training
  fit but does not receive a model-based update. Use
  \code{target.mode = "all"} if missing values may appear later in any
  variable. Complete training data also require
  \code{target.mode = "all"}, because otherwise there are no missing
  variables from which to determine the saved targets.

  If \code{save.ood = TRUE}, the fit also stores an OOD reference in the
  manifest. For each saved target, the out-of-bag prediction from the
  fitted learner is compared with the observed training value to form a
  target-wise reconstruction discrepancy. Continuous and integer targets
  use absolute reconstruction error. Factor targets prefer the negative
  log predictive probability assigned to the observed class. When class
  probabilities are unavailable, unordered factors fall back to a 0/1
  mismatch score and ordered factors fall back to a scaled rank
  distance.

  The row-level OOD calibration stored at fit time is built by
  aggregating the target-wise training scores with a weighted mean using
  \code{weight}. If \code{weight} is omitted at fit time, all saved OOD
  targets receive weight 1. If a named vector is supplied, entries are
  matched by target name, omitted saved OOD targets receive weight 0,
  and the resulting weighting scheme is carried forward in the manifest
  for later deployment-time scoring.

  \code{impute.ood} first completes the predictor side of \code{newdata}
  using the same harmonization, initialization, and iterative sweep
  logic used by \code{predict.impute.learn}. It then reconstructs each
  requested target directly from its saved learner and compares the
  reconstruction with the observed value. Raw target discrepancies are
  converted to target-wise OOD scores using the saved target-specific
  training references, which places continuous and factor targets on a
  common scale.

  The row-level OOD score combines those calibrated target-wise scores
  over the targets that are both observed and scoreable for that row. By
  default, \code{impute.ood} uses a bounded product rule, but the row
  aggregate can be changed to weighted mean, a weighted \eqn{L_p} rule,
  a log-tail weighted \eqn{L_p} rule, or a top-\eqn{k} rule. This makes
  it possible to explore row scores that are more sensitive to sparse
  but severe coordinate shifts.  By default, \code{impute.ood} reuses
  the same OOD weights saved during \code{impute.learn}, so a pipeline
  can fix its weighting scheme once upstream and carry it forward
  automatically.

  A second component, \code{score.percentile}, is obtained by rebuilding
  the row-level training reference from the saved target-wise training
  OOD scores using the requested target subset, the active weight
  vector, and the active row aggregate. This means percentile
  calibration remains available when the user leaves the saved weights
  in place, overrides them at test time, scores only a subset of the
  saved OOD targets, or experiments with alternate row aggregates.

  Unseen factor levels are tracked row-wise during harmonization.
  Because such values are immediate anomalies relative to the training
  schema, \code{impute.ood} flags those rows and assigns them the
  maximum row-level score. If the unseen level occurs in a scored target
  itself, the corresponding target-level discrepancy is also treated as
  maximal.
}
\value{
  \code{impute.learn} returns an object of class
  \code{c("impute.learn.rfsrc", "impute.learn")}. The object
  contains a manifest, optionally the fitted full-sweep learners,
  optionally the completed training data, and optionally a path to the
  saved imputer on disk. If \code{save.ood = TRUE}, the manifest also
  contains an \code{ood} component storing compact target-wise OOD
  references, the saved row-by-target training OOD score matrix used for
  later percentile recalibration, and the default OOD aggregation
  weights.

  \code{load.impute.learn} returns an object of the same class.

  \code{predict.impute.learn} returns a data frame with imputed values
  overlaid. An attribute named \code{"impute.learn.info"} contains
  prediction-time diagnostics such as the number of sweep passes,
  pass-difference history, caching mode, disk-load counts, schema
  harmonization details, row-wise unseen-factor flags, and any targets
  skipped because a learner was unavailable or a prediction failed.

  \code{impute.ood} returns an object of class
  \code{c("impute.ood.rfsrc", "impute.ood")}. It is a list with the
  following components:
  \itemize{
    \item \code{score}: the row-level aggregate of calibrated
      target-wise OOD scores under the requested \code{aggregate} and
      \code{weight}. Larger values indicate greater
      out-of-distribution behavior. For \code{aggregate = "weighted.lp.log"},
      this raw score is on a positive unbounded scale.
    \item \code{score.percentile}: the percentile of \code{score}
      relative to a row-level training reference rebuilt from the saved
      target-wise training OOD scores for the requested targets,
      weights, and row aggregate. For legacy fitted objects that do not
      contain those saved
      training scores, the original saved row-level reference is used
      when possible; otherwise \code{NA}.
    \item \code{targets.used}: the number of weighted targets that
      contributed to each row-level score.
    \item \code{target.score}: optional matrix of target-wise calibrated
      OOD scores, returned when \code{return.details = TRUE}.
    \item \code{target.delta}: optional matrix of raw target-wise
      reconstruction discrepancies, returned when
      \code{return.details = TRUE}.
    \item \code{info}: a list of diagnostics including harmonization
      details, row-wise unseen-factor flags, learner-loading
      information, the active row aggregate and its arguments, whether
      the saved row-level calibration was used, and any target-specific
      issues.
  }
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Stekhoven D.J. and Buhlmann P. (2012). MissForest--non-parametric
  missing value imputation for mixed-type data.
  \emph{Bioinformatics}, 28(1):112--118.

  Tang F. and Ishwaran H. (2017). Random forest missing data
  algorithms. \emph{Statistical Analysis and Data Mining}, 10:363--377.
}
\seealso{
  \code{\link{impute.rfsrc}}, \code{\link{rfsrc}}, and
  \code{\link{predict.rfsrc}}.
}
\examples{
\donttest{
## ------------------------------------------------------------
## small data example: uses missForest for impute engine
## ------------------------------------------------------------

set.seed(101)
aq <- airquality[, c("Ozone", "Solar.R", "Wind", "Temp", "Month")]
aq$Month <- factor(aq$Month)

id <- sample(1:nrow(aq), 100)
train <- aq[id, ]
test <- aq[-id, ]

fit <- impute.learn(
  data = train,
  ntree = 25,
  mf.q = 1,
  max.iter = 5,
  full.sweep.options = list(ntree = 25, nsplit = 5)
)

test.imp <- predict(fit, test, max.predict.iter = 2, verbose = FALSE)
head(test.imp)

## OOD scoring is most informative when every deployment-time
## variable can be reconstructed, so target.mode = "all" is recommended.
## Optional named OOD weights can also be supplied here. Any omitted
## targets receive weight 0, and the saved weights are reused
## automatically later by impute.ood().
ood.fit <- impute.learn(
  data = train,
  ntree = 25,
  mf.q = 1,
  max.iter = 5,
  target.mode = "all",
  save.ood = TRUE,
  full.sweep.options = list(ntree = 25, nsplit = 5),
  verbose = FALSE
)

ood <- impute.ood(ood.fit, test, return.details = TRUE, verbose = FALSE)
head(ood$score)
head(ood$score.percentile)

## try a more spike-sensitive row aggregate
ood.lp <- impute.ood(ood.fit, test,
                     aggregate = "weighted.lp",
                     aggregate.args = list(p = 4),
                     verbose = FALSE)
head(ood.lp$score.percentile)
}

\dontrun{
## ------------------------------------------------------------
## Save the learned imputer to disk and load it later.
## This explicit save example writes learners kept in memory.
## Uses missForest for the impute engine.
## ------------------------------------------------------------

bundle.dir <- file.path(tempdir(), "aq.imputer")

fit <- impute.learn(
  data = train,
  ntree = 25,
  mf.q = 1,
  max.iter = 5,
  full.sweep.options = list(ntree = 25, nsplit = 5),
  keep.models = TRUE,
  verbose = FALSE
)

save.impute.learn(fit, bundle.dir, verbose = FALSE)
imp <- load.impute.learn(bundle.dir, lazy = TRUE, verbose = FALSE)
test.imp <- predict(imp, test, max.predict.iter = 2, verbose = FALSE)

unlink(bundle.dir, recursive = TRUE)



## ------------------------------------------------------------
## Challenging example with factors, uses save/reload
## ------------------------------------------------------------

## load pbc, convert everything to factors
data(pbc, package = "randomForestSRC")
dta <- data.frame(lapply(pbc, factor))
dta$days <- pbc$days
dta$status <- dta$status

## split the data into unbalanced train/test data (25/75)
## the train/test data have the same levels, but different labels
idx <- sample(1:nrow(dta), round(nrow(dta) * .25))
train <- dta[idx,]
test <- dta[-idx,]

## even harder ... factor level not previously encountered in training
levels(test$stage) <- c(levels(test$stage), "fake")
test$stage[sample(seq_len(nrow(test)), 10)] <- "fake"

## train forest
fit <- suppressWarnings(
  impute.learn(Surv(days, status) ~ ., train,
               target.mode = "all",
               save.ood = TRUE,
               keep.models = TRUE)
)

## save/reload
bundle.dir <- file.path(tempdir(), "pbc.imputer")
save.impute.learn(fit, bundle.dir, verbose = FALSE)
imp <- load.impute.learn(bundle.dir, lazy = TRUE, verbose = FALSE)
test.imp <- predict(imp, test, max.predict.iter = 2, verbose = FALSE)
ood <- impute.ood(imp, test, return.details = TRUE, verbose = FALSE)
which(ood$info$unseen.rows)
print(summary(test.imp))
unlink(bundle.dir, recursive = TRUE)
}
}
\keyword{missing data}
