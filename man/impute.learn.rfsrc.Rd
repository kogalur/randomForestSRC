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
  even when new data have missing values. Each selected target is
  reconstructed from its saved conditional learner and compared with the
  observed value. Target-wise discrepancies are calibrated against a
  training reference calculated from out-of-bag predictions computed
  during training.

  If the training data are complete and \code{target.mode = "all"},
  the initial training-data imputation step is skipped and the
  full-sweep learners are fit directly from the complete training data.
  
  If \code{supervised.formula} is supplied, the function also fits an
  internal supervised forest from the training data.  The supervised
  forest plays no role in the training imputation but is used downstream
  to create auxiliary variables which are appended as predictors which
  can improve test-time imputation and OOD scoring. Users who want only
  the basic unsupervised learned imputer and do not want test time
  contamination from outcomes should leave \code{supervised.formula}
  unspecified.
  
}
\usage{
impute.learn.rfsrc(formula, data,
  ntree = 100, nodesize = 1, nsplit = 10,
  nimpute = 2, fast = FALSE, blocks,
  mf.q, max.iter = 10, eps = 0.01,
  ytry = NULL, always.use = NULL, verbose = TRUE,
  ...,
  supervised.formula = NULL,
  supervised.args = list(),
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
  return.reconstruction = FALSE,
  verbose = TRUE,
  ...)
}
\arguments{

  \item{formula}{An optional symbolic model description passed to
    \code{impute} for the initial training-data imputation when
    \code{supervised.formula} is not supplied. It follows the same
    guidelines as in \code{impute}. This formula does not specify the
    full-sweep learner bank, the OOD targets, or the supervised
    auxiliary forest. When \code{supervised.formula} is supplied, the
    raw predictor block is taken from the right-hand side of
    \code{supervised.formula}, and \code{formula} is not used.}

  \item{data}{Training data. Variables that are not real-valued are
    coerced to factors before fitting when possible; otherwise fitting
    stops with an error. Rows and columns that are entirely missing are
    dropped before training begins. If \code{supervised.formula} is
    supplied, \code{data} should contain both the supervised response
    columns and the raw predictor columns. The learned imputer is then
    built on the raw predictor block defined by the right-hand side of
    \code{supervised.formula}.}

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

  \item{supervised.formula}{Optional supervised learning formula used to
    augment the learned imputer for improved OOD detection in supervised
    settings.  The left-hand side defines the supervised response and
    the right-hand side defines the raw predictor block to be learned by
    \code{impute.learn}. The supervised forest is fit internally after
    the raw predictor block has been completed. Its training out-of-bag
    \code{predicted} values and test-time \code{predicted} values are
    appended internally as auxiliary predictors.  These auxiliary
    predictors are internal and are not expected in \code{newdata}.}
  
  \item{supervised.args}{Optional named list of arguments passed to the
    internal supervised \code{rfsrc} fit. Entries for \code{formula},
    \code{data}, and \code{forest} are controlled internally and are
    ignored if supplied. All other entries override the defaults used
    for this supervised fit, including defaults inherited from
    top-level arguments such as \code{fast} and internally supplied
    defaults such as \code{perf.type = "none"}.}

  \item{full.sweep.options}{A \code{list} of options used when fitting
    the full sweep after the training data have been imputed.
    Recognized entries include \code{ntree}, \code{nodesize},
    \code{nsplit}, \code{mtry}, \code{splitrule}, \code{bootstrap},
    \code{sampsize}, \code{samptype}, \code{perf.type}, \code{rfq},
    \code{save.memory}, \code{importance}, and \code{proximity}.
    Unknown entries are ignored with a warning.}

  \item{target.mode}{Determines which raw variables receive a saved
    full-sweep learner. The default \code{"missing.only"} saves
    learners only for variables that were missing in the training data.
    The option \code{"all"} saves a learner for every retained raw
    variable. If the training data are complete, \code{target.mode =
    "all"} must be used. When \code{supervised.formula} is supplied,
    the default is promoted internally to \code{"all"} so that every
    retained raw predictor can later be updated and scored. For the
    broadest OOD coverage, \code{target.mode = "all"} is recommended so
    every deployment-time variable can be reconstructed.}

  \item{deployment.xvars}{Controls which raw predictors are assumed to
    be available later when the saved imputer is used on new data. If
    \code{NULL}, all raw non-target columns are used. If a character
    vector, the same raw predictor set is used for all targets. If a
    named list, names must be target variables and each entry gives the
    raw predictor set for that target; targets omitted from the list
    fall back to all raw non-target columns. Unnamed or duplicated list
    entries are not allowed. Entries named for non-target variables, and
    predictor names not found in the training data, are ignored with a
    warning. When \code{supervised.formula} is supplied, the internal
    auxiliary \code{predicted.*} columns are added automatically on top
    of these raw predictor sets. In other words,
    \code{deployment.xvars} restricts the raw predictors only.}

  \item{anonymous}{If \code{TRUE}, uses \code{rfsrc.anonymous} when
    fitting the saved target-wise full-sweep learner bank. The internal
    supervised forest, when requested through \code{supervised.formula},
    is fit with \code{rfsrc} because it must be retained for later
    prediction-time auxiliary variables. Thus \code{anonymous} usually
    reduces the size of the target-wise learner bank, but not
    necessarily the size of the supervised forest.}

  \item{learner.prefix, learner.root}{Names used when writing saved
    full-sweep learners to disk. If \code{supervised.formula} is
    supplied, the internal supervised forest is saved under the same
    learner root.}

  \item{out.dir}{Optional output directory. If supplied and
    \code{save.on.fit = TRUE}, the manifest and the saved full-sweep
    learners are written to this directory during fitting. This requires
    the \pkg{fst} package because learners are serialized with
    \code{fast.save}. If \code{supervised.formula} is supplied, the
    internal supervised forest is also saved there.}

  \item{wipe}{If \code{TRUE}, removes an existing output directory
    before writing a new one.}

  \item{keep.models}{If \code{TRUE}, keeps the fitted full-sweep
    learners in memory in the returned object. At least one storage mode
    must be enabled: either \code{keep.models = TRUE} or
    \code{out.dir} with \code{save.on.fit = TRUE}. If
    \code{supervised.formula} is supplied, the internal supervised
    forest is kept in memory under the same rule.}

  \item{keep.ximp}{If \code{TRUE}, keeps the completed raw training
    predictor table in the returned object. This is not required for
    later prediction.}

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
    are needed. If \code{FALSE}, all saved learners are loaded at once.
    If \code{supervised.formula} was used at fit time, the internal
    supervised forest follows the same lazy versus eager loading rule.}

  \item{newdata}{New data to be imputed or scored. Missing columns are
    added and extra columns are dropped to match the training schema.
    Unseen factor levels are converted to \code{NA} for harmonization,
    but they are also tracked row-wise. In supervised mode,
    \code{newdata} should contain only the raw predictor columns learned
    by the imputer. The internal auxiliary columns are created
    automatically and should not be supplied by the user. If supervised
    response columns are supplied in \code{newdata} and are not part of
    the learned raw predictor block, they are treated as extra columns
    and dropped; they are not used to compute auxiliary predictors,
    imputations, or OOD scores. In \code{impute.ood},
    observed raw target values are used to compute target-wise
    discrepancies. Missing raw target values are completed for
    predictor-side reconstruction but do not themselves contribute a
    finite target-wise OOD score. Any row containing an unseen factor
    level is flagged and its row-level OOD score is set to the maximum
    value.}

  \item{max.predict.iter}{Maximum number of full-sweep passes applied to
    \code{newdata} before the saved learner bank is used for OOD
    reconstruction or returned prediction-time imputations.}

  \item{restore.integer}{If \code{TRUE}, integer-supported variables
    in the returned data are restored to the integer grid. This includes
    columns stored as R integers and numeric columns whose observed
    finite training values were all integer-valued. Numeric columns
    remain numeric; returned values are rounded to integer support.
    Columns stored as R integers are returned as integer vectors.
    Factor columns are always conformed back to the training schema.
    The package operates on real-valued and factor variables;
    inputs that are not real-valued are coerced to factors during
    preprocessing when possible, otherwise an error is raised.}

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
    \code{impute.ood} to combine calibrated target-wise OOD scores. The
    default \code{"bounded.product"} applies a weighted product of
    the form \eqn{1-\prod_j (1-u_j+\varepsilon)^{\tilde w_j}} where
    \eqn{u_j} are row level calibrated scores for target feature
    \eqn{j}. \code{"weighted.mean"} is the weighted
    average. \code{"weighted.lp"} applies a weighted Minkowski \eqn{L_p}
    aggregation to the calibrated target scores.
    \code{"weighted.lp.log"} first applies the tail-stretching transform
    \eqn{-\log(1-u_j+\varepsilon)} to each calibrated target score
    \eqn{u_j} and then applies the weighted \eqn{L_p}
    aggregation. \code{"top.k"} averages only the \eqn{k} largest
    scoreable target scores among the positively weighted targets. When
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
    richer \code{info} list.}
  
  \item{return.reconstruction}{If \code{TRUE}, \code{impute.ood} returns
    the target-wise reconstructed values used during OOD scoring.  The
    reconstructed values are returned on the raw target scale, with one
    column per scored target.  Typically these values are computed
    internally and discarded.}


}
\details{
  A predictive imputer is calculated in two stages.

  The training data are first processed. Variables that are not
  real-valued are coerced to factors when possible; otherwise fitting
  stops with an error. Rows and columns that are entirely missing are
  removed before the training schema is stored.

  If the resulting training data contain missing values, the first stage
  uses \code{impute} as the imputation engine to complete the training
  data.  Options are specified exactly as in \code{impute}.  In
  particular, \code{mf.q = 1} gives standard \code{missForest},
  \code{mf.q > 1} gives the multivariate \code{missForest}
  generalization, and if \code{mf.q} is omitted, on the fly imputation
  is used if \code{formula} is specified, otherwise default unsupervised
  imputation is used. If the training data are already complete and
  \code{target.mode = "all"}, this initial imputation step is skipped.

  In the second stage, a full forest sweep is fit on the completed
  training data. For each target selected by \code{target.mode}, rows
  where that target was observed are used to fit a forest with that
  target on the left-hand side and the selected deployment predictors on
  the right-hand side are specified by \code{deployment.xvars}.

  By default, \code{deployment.xvars = NULL} allows every non-target
  column to be used as a predictor. This is convenient, but it can also
  introduce leakage if the training data include outcomes, future-only
  variables, identifiers, or any fields that will not be available when
  the learned imputer is applied to new data. Restrict
  \code{deployment.xvars} when that is a concern.

  If \code{supervised.formula} is supplied, the function also fits an
  internal supervised forest from the training data. Out-of-bag
  \code{predicted} values for the training sample, and \code{predicted}
  values for new data, are appended internally as auxiliary
  predictors. These auxiliary variables are not used as imputation
  targets, however, they are available as predictors for the saved
  forests and can therefore influence both test-time imputation and OOD
  scoring. Users who want only the basic unsupervised learned imputer
  should leave \code{supervised.formula} unspecified.
  
  When \code{supervised.formula} is used, \code{deployment.xvars}
  continues to restrict only the raw predictors. The internally created
  auxiliary variables are added automatically to the predictor sets for
  the saved raw targets. This lets the learned imputer benefit from
  supervised signal without asking the user to prepare and pass those
  auxiliary columns manually.

  When the imputer is saved to disk, each full-sweep learner is written
  separately using \code{fast.save}. Loading uses \code{fast.load}. In
  practice this gives a small manifest plus a directory of saved
  learners. The \pkg{fst} package is therefore required for save and
  load operations. The explicit save method can write learners either
  from memory or by reloading them from an attached saved path. If
  \code{supervised.formula} is used, the internal supervised forest is
  saved alongside the target-wise learner bank.

  Prediction starts by matching \code{newdata} to the training schema,
  filling missing values with training means or modes, and then applying
  one or more full-sweep passes. Only the targets selected by
  \code{target.mode} are updated by saved learners. In supervised mode,
  this same test-time preparation first initializes the raw predictor
  table, then computes the internal auxiliary variables from the saved
  supervised forest, and then runs the forest test-time sweep.
  Supervised response columns, if present in \code{newdata} but not in
  the learned raw predictor block, are dropped before this step and play
  no role at prediction or scoring time.

  If \code{target.mode = "missing.only"}, a variable that was complete
  in training but missing in new data is initialized from the training
  fit but does not receive a model-based update. Use
  \code{target.mode = "all"} if missing values may appear later in any
  raw variable. Complete training data also require
  \code{target.mode = "all"}, because otherwise there are no missing
  variables from which to determine the saved targets.

  If \code{save.ood = TRUE}, the fit also stores an OOD reference in the
  manifest. For each saved target, the out-of-bag prediction from the
  fitted learner is compared with the observed training value to form a
  target-wise reconstruction discrepancy. Continuous and integer targets
  use absolute reconstruction error. Factor targets use negative log
  predictive probabilities.
  
  The row-level OOD calibration stored at fit time is built by
  aggregating the target-wise training scores with a weighted mean using
  \code{weight}. If \code{weight} is omitted at fit time, all saved OOD
  targets receive weight 1. If a named vector is supplied, entries are
  matched by target name, omitted saved OOD targets receive weight 0,
  and the resulting weighting scheme is carried forward for later
  deployment-time scoring.

  \code{impute.ood} first completes the predictor side of \code{newdata}
  using the same harmonization, initialization, and iterative sweep
  logic used by \code{predict.impute.learn}. It then reconstructs each
  requested raw target directly from its saved learner and compares the
  reconstruction with the observed raw target value supplied in
  \code{newdata}. If a scored raw target is missing in a row, that
  target does not contribute to that row's OOD score. Raw target
  discrepancies are converted to target-wise OOD scores using the saved
  target-specific training references. In supervised mode, auxiliary
  variables participate as predictors in the reconstruction, but the
  supervised response itself is not used at deployment time.

  The row-level OOD score combines those calibrated target-wise scores
  over the targets that are both observed and scoreable for that row. By
  default, \code{impute.ood} uses a bounded product rule, but the row
  aggregate can be changed to weighted mean, a weighted \eqn{L_p} rule,
  a log-tail weighted \eqn{L_p} rule, or a top-\eqn{k} rule. This makes
  it possible to explore row scores that are more sensitive to sparse
  but severe coordinate shifts. By default, \code{impute.ood} reuses
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
  optionally the internal supervised forest when
  \code{supervised.formula} is used, optionally the completed raw
  training predictor table, and optionally a path to the saved imputer
  on disk. If \code{save.ood = TRUE}, the manifest also contains an
  \code{ood} component storing compact target-wise OOD references, the
  saved row-by-target training OOD score matrix used for later
  percentile recalibration, and the default OOD aggregation weights.
  When supervised mode is active, the manifest also records the
  supervised family, response names, and the internally created
  auxiliary \code{predicted.*} variable names.

  \code{load.impute.learn} returns an object of the same class.

  \code{predict.impute.learn} returns a data frame with imputed values
  overlaid on the raw predictor table. An attribute named
  \code{"impute.learn.info"} contains prediction-time diagnostics such
  as the number of sweep passes, pass-difference history, caching mode,
  disk-load counts, schema harmonization details, dropped supervised
  response columns when present, row-wise unseen-factor flags,
  supervised-auxiliary diagnostics when present, and any targets
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
      contain those saved training scores, the original saved row-level
      reference is used when possible; otherwise \code{NA}.
    \item \code{targets.used}: the number of weighted targets that
      contributed to each row-level score.
    \item \code{target.score}: optional matrix of target-wise calibrated
      OOD scores, returned when \code{return.details = TRUE}.
    \item \code{target.delta}: optional matrix of raw target-wise
      reconstruction discrepancies, returned when
      \code{return.details = TRUE}.
    \item \code{info}: a list of diagnostics including harmonization
      details, dropped supervised response columns when present,
      row-wise unseen-factor flags, learner-loading information,
      supervised-auxiliary diagnostics when present, the active row
      aggregate and its arguments, whether the saved row-level
      calibration was used, and any target-specific issues.
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
  \code{\link{impute.rfsrc}},
  \code{\link{rfsrc}},
  \code{\link{predict.rfsrc}}.
}
\examples{
## ------------------------------------------------------------
## small data example: uses missForest for impute engine
## ------------------------------------------------------------

set.seed(101)
aq <- airquality[, c("Ozone", "Solar.R", "Wind", "Temp", "Month")]
aq$Month <- factor(aq$Month)

id <- sample(1:nrow(aq), 100)
train <- aq[id, ]
test <- aq[-id, ]

## training the imputer
fit <- impute.learn(
  data = train,
  ntree = 25,
  mf.q = 1,
  max.iter = 5,
  full.sweep.options = list(ntree = 25, nsplit = 5)
)

## test time imputation
test.imp <- predict(fit, test, max.predict.iter = 2, verbose = FALSE)
print(head(test.imp))

\donttest{

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
print(head(ood$score))
print(head(ood$score.percentile))

## try a more spike-sensitive row aggregate
ood.lp <- impute.ood(ood.fit, test,
                     aggregate = "weighted.lp",
                     aggregate.args = list(p = 4),
                     verbose = FALSE)
print(head(ood.lp$score.percentile))

## ------------------------------------------------------------
## supervised OOD example: regression benchmark
## the user supplies only the raw x-columns at test time
## ------------------------------------------------------------

friedman1_sim <- function(n = 150, p = 10, sigma = 1) {
  X <- matrix(runif(n * p), nrow = n)
  y <- 10 * sin(pi * X[, 1] * X[, 2]) +
       20 * (X[, 3] - 0.5)^2 +
       10 * X[, 4] + 5 * X[, 5] +
       rnorm(n, sd = sigma)
  list(X = X, y = y)
}


trn <- data.frame(friedman1_sim())
tst <- data.frame(friedman1_sim())
xvars <- setdiff(names(trn), "y")


## impute data using missForest, construct a supervised forest
## - supervised forests are used to create auxiliary variables
## - improves test time OOD in supervised problems

sup.fit <- impute.learn(
  data = trn,
  mf.q = 1,
  supervised.formula = y ~ .,
  supervised.args = list(ntree = 50, nsplit = 5),
  full.sweep.options = list(ntree = 25, nsplit = 5),
  save.ood = TRUE,
  verbose = FALSE
)

## add some missing values to the test data
xnew <- tst[, xvars, drop = FALSE]
xnew[sample(seq_len(nrow(xnew)), 5), xvars[1]] <- NA
xnew[sample(seq_len(nrow(xnew)), 5), xvars[2]] <- NA

## imputation
xnew.imp <- predict(sup.fit, xnew, max.predict.iter = 2, verbose = FALSE)
print(head(xnew.imp))

## OOD score
ood.sup <- impute.ood(sup.fit, xnew, verbose = FALSE)
print(head(ood.sup$score.percentile))

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
print(which(ood$info$unseen.rows))
print(summary(test.imp))
unlink(bundle.dir, recursive = TRUE)
}
}
\keyword{missing data}
