\name{rfsrc}
\alias{rfsrc}
\alias{rfsrc.cart}
\alias{randomForestSRC}
\title{Fast Unified Random Forests for Survival, Regression, and Classification (RF-SRC)}
\description{

Fast OpenMP-parallel implementation of random forests (Breiman, 2001)
for regression, classification, survival analysis (Ishwaran et al.,
2008), competing risks (Ishwaran et al., 2012), multivariate outcomes
(Segal and Xiao, 2011), unsupervised learning (Mantero and Ishwaran,
2020), quantile regression (Meinshausen, 2006; Zhang et al., 2019;
Greenwald and Khanna, 2001), and imbalanced q-classification (O'Brien
and Ishwaran, 2019).
  

The package supports both deterministic and randomized splitting rules
(Geurts et al., 2006; Ishwaran, 2015) across all families. Multiple
types of variable importance (VIMP) are available, including holdout
VIMP and confidence regions (Ishwaran and Lu, 2019), for both individual
and grouped variables. Variable selection can be performed using minimal
depth (Ishwaran et al., 2010, 2011). Fast interfaces for missing data
imputation are provided using several forest-based algorithms (Tang and
Ishwaran, 2017).


\strong{Highlighted updates:}

\enumerate{
  \item For variable selection, we recommend using \pkg{VarPro}, an R
  package for model-independent variable selection using rule-based
  variable priority. It supports regression, classification, survival
  analysis, and includes a new mode for unsupervised learning. See
  \url{https://www.varprotools.org} for more information.
  

  \item For computational speed, the default VIMP method has changed
  from "permute" (Breiman-Cutler permutation) to "anti"
  (\code{importance = "anti"} or \code{importance = TRUE}). While
  faster, this may be less accurate in settings such as highly
  imbalanced classification. To revert to permutation VIMP, use
  \code{importance = "permute"}.
}


This is the main entry point to the \pkg{randomForestSRC}
package. For more information on OpenMP support and the package as a
whole, see \code{package?randomForestSRC}.

}
\usage{
rfsrc(formula, data, ntree = 500,
  mtry = NULL, ytry = NULL,
  nodesize = NULL, nodedepth = NULL,
  splitrule = NULL, nsplit = NULL,
  importance = c(FALSE, TRUE, "none", "anti", "permute", "random"),
  block.size = if (any(is.element(as.character(importance),
                     c("none", "FALSE")))) NULL else 10,
  bootstrap = c("by.root", "none", "by.user"),
  samptype = c("swor", "swr"), samp = NULL, membership = FALSE,
  sampsize = if (samptype == "swor") function(x){x * .632} else function(x){x},
  na.action = c("na.omit", "na.impute"), nimpute = 1,
  ntime = 150, cause,
  perf.type = NULL,
  proximity = FALSE, distance = FALSE, forest.wt = FALSE,
  xvar.wt = NULL, yvar.wt = NULL, split.wt = NULL, case.wt  = NULL,
  case.depth = FALSE,
  forest = TRUE,
  save.memory = FALSE,
  var.used = c(FALSE, "all.trees", "by.tree"),
  split.depth = c(FALSE, "all.trees", "by.tree"),
  seed = NULL,
  do.trace = FALSE,
  ...)

## convenient interface for growing a CART tree
rfsrc.cart(formula, data, ntree = 1, mtry = ncol(data),
   bootstrap = "none", nsplit = 0, ...)
}

\arguments{

  \item{formula}{A formula describing the model to fit. Interaction terms are not supported. If missing, unsupervised splitting is used.}

  \item{data}{Data frame containing the response and predictor variables.}

  \item{ntree}{Number of trees to grow.}

  \item{mtry}{Number of candidate variables randomly selected at each split. Defaults: regression uses \code{p/3}, others use \code{sqrt(p)}; rounded up.}

  \item{ytry}{Number of pseudo-response variables randomly selected for unsupervised splitting. Default is 1.}

  \item{nodesize}{Minimum terminal node size. Defaults: survival/competing risks (15), regression (5), classification (1), mixed/unsupervised (3).}

  \item{nodedepth}{Maximum tree depth. Ignored by default.}

  \item{splitrule}{Splitting rule. See Details.}

  \item{nsplit}{Number of random split points per variable. \code{0} uses all values (deterministic). Default is 10.}

  \item{importance}{Variable importance (VIMP) method. Choices: \code{FALSE}, \code{TRUE}, \code{"none"}, \code{"anti"}, \code{"permute"}, \code{"random"}. Default is \code{"none"}. VIMP can be computed later using \code{vimp} or \code{predict}.}

  \item{block.size}{Controls frequency of cumulative error/VIMP updates. Default is \code{10} if importance is requested; otherwise \code{NULL}. See Details.}

  \item{bootstrap}{Bootstrap method. Options: \code{"by.root"} (default), \code{"by.user"}, or \code{"none"} (no OOB error possible).}

  \item{samptype}{Sampling type for \code{by.root} bootstrap. Options: \code{"swor"} (without replacement, default), \code{"swr"} (with replacement).}

  \item{samp}{Bootstrap weights (only for \code{bootstrap="by.user"}). A matrix of size \code{n x ntree} giving in-bag counts per tree.}

  \item{membership}{Return inbag and terminal node membership?}

  \item{sampsize}{Bootstrap sample size (used when \code{bootstrap="by.root"}). Defaults: 0.632 x n for \code{swor}; n for \code{swr}. Can also be numeric.}

  \item{na.action}{Missing data handling. \code{"na.omit"} (default) removes rows with any \code{NA}; \code{"na.impute"} performs fast internal imputation. See also \code{impute}.}

  \item{nimpute}{Number of iterations for internal imputation. If >1, OOB error rates may be optimistic.}

  \item{ntime}{For survival models: number or grid of time points used in ensemble estimation. If \code{NULL} or \code{0}, uses all event times.}

  \item{cause}{For competing risks: event of interest (1 to \code{J}), or a vector of weights over all \code{J} events. Defaults to an average over all events.}

  \item{perf.type}{Optional performance metric for prediction, VIMP, and error. Defaults to the family-specific metric. \code{"none"} disables performance. See Details.}

  \item{proximity}{Compute proximity matrix? Options: \code{"inbag"}, \code{"oob"}, \code{"all"}, \code{TRUE} (inbag), or \code{FALSE}.}

  \item{distance}{Compute pairwise distances between cases? Similar options as \code{proximity}. See Details.}

  \item{forest.wt}{Return forest weight matrix? Same options as \code{proximity}. Default is \code{TRUE} (inbag).}

  \item{xvar.wt}{Optional weights on x-variables for sampling at splits. Does not need to sum to 1. Defaults to uniform.}

  \item{yvar.wt}{Weights on response variables (for multivariate regression). Used when \code{y} is high-dimensional.}

  \item{split.wt}{Weights applied to each variable's split statistic. Higher weight increases likelihood of splitting.}

  \item{case.wt}{Sampling weights for cases in the bootstrap. Higher values increase selection probability. See class imbalance example.}

  \item{case.depth}{Return matrix recording depth of first split for each case? Default is \code{FALSE}.}

  \item{forest}{Save forest object for future prediction? Set \code{FALSE} if prediction is not needed.}

  \item{save.memory}{Reduce memory usage by avoiding storage of prediction quantities. Recommended for large survival or competing risk forests.}

  \item{var.used}{Return variable usage statistics? Options: \code{FALSE}, \code{"all.trees"}, \code{"by.tree"}.}

  \item{split.depth}{Return minimal depth of splits for each variable? Options: \code{FALSE}, \code{"all.trees"}, \code{"by.tree"}. See Details.}

  \item{seed}{Integer seed for reproducibility (negative values only).}

  \item{do.trace}{Print progress updates every \code{do.trace} seconds.}

  \item{...}{Additional arguments passed to or from other methods.}
}

\details{
  \enumerate{

  \item \emph{Types of forests}

  The type of forest is automatically inferred from the outcome and formula. Supported forest types include:

  \itemize{
    \item Regression forests for continuous outcomes.
    \item Classification forests for factor outcomes.
    \item Multivariate forests for continuous, categorical, or mixed outcomes.
    \item Unsupervised forests when no outcome is specified.
    \item Survival forests for right-censored time-to-event data.
    \item Competing risk forests for multi-event survival settings.
  }

  \item \emph{Splitting}

  \enumerate{
    \item Splitting rules are set using the \code{splitrule} option.
    \item Random splitting is invoked via \code{splitrule = "random"}.
    \item Use \code{nsplit} to enable randomized splitting and improve speed; see \emph{Improving computational speed}.
  }

  \item \emph{Available splitting rules}

  \itemize{

    \item \strong{Regression}
    \enumerate{
      \item \code{"mse"} (default): weighted mean squared error (Breiman et al., 1984).
      \item \code{"quantile.regr"}: quantile regression via check-loss; see \code{quantreg.rfsrc}.
      \item \code{"la.quantile.regr"}: local adaptive quantile regression.
    }

    \item \strong{Classification}
    \enumerate{
      \item \code{"gini"} (default): Gini index.
      \item \code{"auc"}: AUC-based splitting; appropriate for imbalanced data.
      \item \code{"entropy"}: entropy-based splitting.
    }

    \item \strong{Survival}
    \enumerate{
      \item \code{"logrank"} (default): log-rank splitting.
      \item \code{"bs.gradient"}: Brier score gradient splitting. Uses 90th percentile of observed times by default, or set \code{prob}.
      \item \code{"logrankscore"}: log-rank score splitting.
    }

    \item \strong{Competing risks} (see Ishwaran et al., 2014)
    \enumerate{
      \item \code{"logrankCR"} (default): Gray's test-based weighted log-rank splitting.
      \item \code{"logrank"}: cause-specific weighted log-rank; use \code{cause} to target specific events.
    }

    \item \strong{Multivariate}
    \enumerate{
      \item Default: normalized composite splitting (Tang and Ishwaran, 2017).
      \item \code{"mahalanobis"}: Mahalanobis splitting with optional covariance matrix; for multivariate regression.
    }

    \item \strong{Unsupervised}
    Splitting uses pseudo-outcomes and the composite rule. See \code{sidClustering} for advanced unsupervised analysis.

    \item \strong{Custom splitting}
    Custom rules can be defined using \code{splitCustom.c}. Up to 16 rules per family are allowed. Use \code{"custom"}, \code{"custom1"}, etc. Compilation required.
  }

  \item \emph{Improving computational speed}

  See \code{\link{rfsrc.fast}}. Strategies include:

  \itemize{
    \item Increase \code{nodesize}.
    \item Set \code{save.memory = TRUE} for large survival or competing risk models.
    \item Set \code{block.size = NULL} to avoid repeated cumulative error computation.
    \item Use \code{perf.type = "none"} to disable VIMP and C-index calculations.
    \item Set \code{nsplit} to a small integer (e.g., 1-10).
    \item Reduce bootstrap size with \code{sampsize}, \code{samptype}.
    \item Set \code{ntime} to a coarse grid (e.g., 50) for survival models.
    \item Pre-filter variables; use \code{max.subtree} for fast variable selection.
  }

  \item \emph{Prediction Error}

  Error is computed using OOB data:
  \itemize{
    \item Regression: mean squared error.
    \item Classification: misclassification rate, Brier score, AUC.
    \item Survival: C-error = 1 - Harrell's concordance index.
  }

  If \code{bootstrap = "none"}, OOB error is unavailable. Use \code{predict.rfsrc} for cross-validation error instead.

  \item \emph{Variable Importance (VIMP)}

  VIMP methods:
  \itemize{
    \item \code{"permute"}: permutation VIMP (Breiman-Cutler).
    \item \code{"random"}: randomized left/right assignment.
    \item \code{"anti"} (default): anti-split assignment.
  }

  The \code{block.size} option controls granularity. For confidence intervals, see \code{subsampling}. Also see \code{holdout.vimp} for a more conservative variant.

  \item \emph{Multivariate Forests}

  Use:
  \preformatted{rfsrc(Multivar(y1, ..., yd) ~ ., data)}
  or
  \preformatted{rfsrc(cbind(y1, ..., yd) ~ ., data)}

  Use \code{get.mv.formula}, \code{get.mv.predicted}, \code{get.mv.error} for multivariate extraction.

  \item \emph{Unsupervised Forests}

  Use:
  \preformatted{rfsrc(data = X)}
  or
  \preformatted{rfsrc(Unsupervised(ytry) ~ ., data = X)}

  Random subsets of \code{ytry} pseudo-responses are used for each \code{mtry} variable. No performance metrics are computed.

  \item \emph{Survival, Competing Risks}

  \itemize{
    \item Survival: use \code{Surv(time, status) ~ .}. Status must be 0 (censored) or 1 (event).
    \item Competing risks: status = 0 (censored), 1-J (event types). Use \code{cause} to target specific events.
    \item Larger \code{nodesize} is typically needed for competing risks.
  }

  \item \emph{Missing data imputation}

  Use \code{na.action = "na.impute"}. Iteration with \code{nimpute > 1} replaces missing values using OOB predictions. Observations or variables with all missing values are removed.

  \item \emph{Allowable data types and factors}

  Variables must be numeric, integer, factor, or logical. Non-factors are coerced to numeric. For unordered factors, all complementary subsets are considered for splits.

  Factor levels are mapped to ensure consistency across training/test data. Consider converting factors to numeric for high-dimensional settings.
}}

\value{
An object of class \code{(rfsrc, grow)} with the following components:

\describe{
  \item{call}{The original call to \code{rfsrc}.}
  \item{family}{The family used in the analysis.}
  \item{n}{Sample size after applying \code{na.action}.}
  \item{ntree}{Number of trees grown.}
  \item{mtry}{Number of variables randomly selected at each node.}
  \item{nodesize}{Minimum terminal node size.}
  \item{nodedepth}{Maximum depth allowed for each tree.}
  \item{splitrule}{Splitting rule used.}
  \item{nsplit}{Number of random split points.}
  \item{yvar}{Response values.}
  \item{yvar.names}{Character vector of response variable names.}
  \item{xvar}{Data frame of predictor variables.}
  \item{xvar.names}{Character vector of predictor variable names.}
  \item{xvar.wt}{Non-negative weights specifying the selection probability of each variable.}
  \item{split.wt}{Non-negative weights adjusting each variable's split statistic.}
  \item{cause.wt}{Weights for composite competing risk splitting.}
  \item{leaf.count}{Number of terminal nodes per tree. A value of 0 indicates a rejected tree (may occur with missing data); a value of 1 indicates a stump.}
  \item{proximity}{Proximity matrix indicating how often case pairs fall in the same terminal node.}
  \item{forest}{Forest object, returned if \code{forest=TRUE}. Required for prediction and most wrappers.}
  \item{forest.wt}{Forest weight matrix.}
  \item{membership}{Terminal node membership matrix (rows: trees; columns: cases).}
  \item{inbag}{Inbag count matrix (rows: trees; columns: cases).}
  \item{var.used}{Number of times each variable is used to split a node.}
  \item{imputed.indv}{Indices of individuals with missing values.}
  \item{imputed.data}{Imputed dataset with responses followed by predictors.}
  \item{split.depth}{Matrix or array recording minimal split depth of variables by case and tree.}
  \item{err.rate}{Cumulative OOB error rate.}
  \item{err.block.rate}{Cumulative error per ensemble block (size defined by \code{block.size}). If \code{block.size = 1}, error per tree.}
  \item{importance}{Variable importance (VIMP) for each predictor.}
  \item{predicted}{In-bag predicted values.}
  \item{predicted.oob}{Out-of-bag (OOB) predicted values.}

  \item{class}{(Classification) In-bag predicted class labels.}
  \item{class.oob}{(Classification) OOB predicted class labels.}

  \item{regrOutput}{(Multivariate) List of performance results for continuous outcomes.}
  \item{clasOutput}{(Multivariate) List of performance results for categorical outcomes.}

  \item{survival}{(Survival) In-bag survival functions.}
  \item{survival.oob}{(Survival) OOB survival functions.}
  \item{chf}{(Survival or competing risks) In-bag cumulative hazard function.}
  \item{chf.oob}{(Survival or competing risks) OOB cumulative hazard function.}
  \item{time.interest}{(Survival or competing risks) Unique sorted event times.}
  \item{ndead}{(Survival or competing risks) Total number of observed events.}

  \item{cif}{(Competing risks) In-bag cumulative incidence function by cause.}
  \item{cif.oob}{(Competing risks) OOB cumulative incidence function by cause.}
}
}

\note{
Values returned by the forest depend on the family:

\itemize{
  \item \strong{Regression}: \code{predicted} and \code{predicted.oob} are vectors of predicted values.
  
  \item \strong{Classification}: \code{predicted} and \code{predicted.oob} are matrices of class probabilities. VIMP and performance metrics are returned as a matrix with \code{J+1} columns (J = number of classes). The first column ("all") gives unconditional results; remaining columns give class-conditional results.
  
  \item \strong{Survival}: \code{predicted} contains mortality estimates (Ishwaran et al., 2008). These are calibrated to the number of expected events under identical covariate profiles. Also returned are matrices of the survival function and CHF for each individual over \code{time.interest}.
  
  \item \strong{Competing risks}: \code{predicted} contains expected life years lost by cause (Ishwaran et al., 2013). Also returned are three-dimensional arrays for CIF and CSCHF indexed by case, time, and event type.
  
  \item \strong{Multivariate}: Predictions, VIMP, and error rates are returned in \code{regrOutput} and \code{clasOutput}. Use \code{get.mv.predicted}, \code{get.mv.vimp}, and \code{get.mv.error} to extract results.
}
}

  
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}

\references{
  Breiman L., Friedman J.H., Olshen R.A. and Stone C.J. (1984).
  \emph{Classification and Regression Trees}, Belmont, California.
  
  Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.
  
  Cutler A. and Zhao G. (2001). PERT-Perfect random tree ensembles.
  \emph{Comp. Sci. Statist.}, 33: 490-497.

  Dietterich, T. G. (2000). An experimental comparison of three methods for
  constructing ensembles of decision trees: bagging, boosting, and randomization.
  \emph{Machine Learning}, 40, 139-157.
  
  Gray R.J. (1988).  A class of k-sample tests for comparing the
  cumulative incidence of a competing risk, \emph{Ann. Statist.},
  16: 1141-1154.
  
  Geurts, P., Ernst, D. and Wehenkel, L., (2006). Extremely randomized
  trees. \emph{Machine learning}, 63(1):3-42.

  Greenwald M. and Khanna S. (2001).  Space-efficient online computation of
  quantile summaries. \emph{Proceedings of ACM SIGMOD}, 30(2):58-66.

  Harrell et al. F.E. (1982).  Evaluating the yield of medical tests,
  \emph{J. Amer. Med. Assoc.}, 247:2543-2546.
  
  Hothorn T. and Lausen B. (2003). On the exact distribution of maximally selected
  rank statistics, \emph{Comp. Statist. Data Anal.}, 43:121-137.
  
  Ishwaran H. (2007).  Variable importance in binary regression
  trees and forests,  \emph{Electronic J. Statist.}, 1:519-537.
  
  Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7(2):25-31.
  
  Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
  (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.
  
  Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and
  Lauer M.S. (2010).  High-dimensional variable selection for survival
  data.  \emph{J. Amer. Statist. Assoc.}, 105:205-217.
  
  Ishwaran H., Kogalur U.B., Chen X. and Minn A.J. (2011). Random survival
  forests for high-dimensional data. \emph{Stat. Anal. Data Mining}, 4:115-132
  
  Ishwaran H., Gerds T.A., Kogalur U.B., Moore R.D., Gange S.J. and Lau
  B.M. (2014). Random survival forests for competing risks.
  \emph{Biostatistics}, 15(4):757-773.
  
  Ishwaran H. and Malley J.D. (2014). Synthetic learning
  machines. \emph{BioData Mining}, 7:28.

  Ishwaran H. (2015).  The effect of splitting on random forests.
  \emph{Machine Learning}, 99:75-118.

  Lin, Y. and Jeon, Y. (2006). Random forests and adaptive nearest
  neighbors. \emph{J. Amer. Statist. Assoc.}, 101(474), 578-590.
  
  Lu M., Sadiq S., Feaster D.J. and Ishwaran H. (2018). Estimating
  individual treatment effect in observational data using random forest
  methods. \emph{J. Comp. Graph. Statist}, 27(1), 209-219

  Ishwaran H. and Lu M.  (2019).  Standard errors and confidence
  intervals for variable importance in random forest regression,
  classification, and survival. \emph{Statistics in Medicine}, 38,
  558-582.
  
  LeBlanc M. and Crowley J. (1993).  Survival trees by goodness of split,
  \emph{J. Amer. Statist. Assoc.}, 88:457-467.
  
  Loh W.-Y and Shih Y.-S (1997).  Split selection methods for
  classification trees, \emph{Statist. Sinica}, 7:815-840.

  Mantero A. and Ishwaran H. (2021).  Unsupervised random forests.
  \emph{Statistical Analysis and Data Mining}, 14(2):144-167.

  Meinshausen N. (2006) Quantile regression forests, \emph{Journal of
  Machine Learning Research}, 7:983-999.

  Mogensen, U.B, Ishwaran H. and Gerds T.A. (2012). Evaluating random
  forests for survival analysis using prediction error curves,
  \emph{J.  Statist. Software}, 50(11): 1-23.
  
  O'Brien R. and Ishwaran H. (2019).  A random forests quantile
  classifier for class imbalanced data. \emph{Pattern Recognition},
  90, 232-249

  Segal M.R. (1988).  Regression trees for censored data,
  \emph{Biometrics},  44:35-47.

  Segal M.R. and Xiao Y. Multivariate random
  forests. (2011). \emph{Wiley Interdisciplinary Reviews: Data Mining
  and Knowledge Discovery}. 1(1):80-87.
  
  Tang F. and Ishwaran H. (2017).  Random forest missing data
  algorithms.  \emph{Statistical Analysis and Data Mining}, 10:363-377.

  Zhang H., Zimmerman J., Nettleton D. and Nordman D.J. (2019).  Random
  forest prediction intervals. \emph{The American Statistician}. 4:1-5.
  

}
\seealso{

  \command{\link{get.tree.rfsrc}},

  \command{\link{holdout.vimp.rfsrc}},

  \command{\link{imbalanced.rfsrc}},
  \command{\link{impute.rfsrc}},

  \command{\link{max.subtree.rfsrc}},

  \command{\link{partial.rfsrc}},
  \command{\link{plot.competing.risk.rfsrc}},
  \command{\link{plot.rfsrc}},
  \command{\link{plot.survival.rfsrc}},
  \command{\link{plot.variable.rfsrc}},
  \command{\link{predict.rfsrc}},
  \command{\link{print.rfsrc}},

  \command{\link{quantreg.rfsrc}},

  \command{\link{rfsrc}},
  \command{\link{rfsrc.anonymous}},
  \command{\link{rfsrc.cart}},
  \command{\link{rfsrc.fast}},

  
  \command{\link{sidClustering.rfsrc}},

  \command{\link{subsample.rfsrc}},
  
  \command{\link{tune.rfsrc}},

  \command{\link{vimp.rfsrc}}
}
\examples{
\donttest{
##------------------------------------------------------------
## survival analysis
##------------------------------------------------------------

## veteran data
## randomized trial of two treatment regimens for lung cancer
data(veteran, package = "randomForestSRC")
v.obj <- rfsrc(Surv(time, status) ~ ., data = veteran, block.size = 1)

## plot tree number 3
plot(get.tree(v.obj, 3))

## print results of trained forest
print(v.obj)

## plot results of trained forest
plot(v.obj)

## plot survival curves for first 10 individuals -- direct way
matplot(v.obj$time.interest, 100 * t(v.obj$survival.oob[1:10, ]),
    xlab = "Time", ylab = "Survival", type = "l", lty = 1)

## plot survival curves for first 10 individuals
## using function "plot.survival" 
plot.survival(v.obj, subset = 1:10)

## obtain Brier score using KM and RSF censoring distribution estimators
bs.km <- get.brier.survival(v.obj, cens.model = "km")$brier.score
bs.rsf <- get.brier.survival(v.obj, cens.model = "rfsrc")$brier.score

## plot the brier score
plot(bs.km, type = "s", col = 2)
lines(bs.rsf, type ="s", col = 4)
legend("topright", legend = c("cens.model = km", "cens.model = rfsrc"), fill = c(2,4))

## plot CRPS (continuous rank probability score) as function of time
## here's how to calculate the CRPS for every time point
trapz <- randomForestSRC:::trapz
time <- v.obj$time.interest
crps.km <- sapply(1:length(time), function(j) {
  trapz(time[1:j], bs.km[1:j, 2] / diff(range(time[1:j])))
})
crps.rsf <- sapply(1:length(time), function(j) {
  trapz(time[1:j], bs.rsf[1:j, 2] / diff(range(time[1:j])))
})
plot(time, crps.km, ylab = "CRPS", type = "s", col = 2)
lines(time, crps.rsf, type ="s", col = 4)
legend("bottomright", legend=c("cens.model = km", "cens.model = rfsrc"), fill=c(2,4))


## fast nodesize optimization for veteran data
## optimal nodesize in survival is larger than other families
## see the function "tune" for more examples
tune.nodesize(Surv(time,status) ~ ., veteran)


## Primary biliary cirrhosis (PBC) of the liver
data(pbc, package = "randomForestSRC")
pbc.obj <- rfsrc(Surv(days, status) ~ ., pbc)
print(pbc.obj)


## save.memory example for survival
## growing many deep trees creates memory issue without this option!
data(pbc, package = "randomForestSRC")
print(rfsrc(Surv(days, status) ~ ., pbc, splitrule = "random",
            ntree = 25000, nodesize = 1, save.memory = TRUE))



##------------------------------------------------------------
## trees can be plotted for any family 
## see get.tree for details and more examples
##------------------------------------------------------------

## survival where factors have many levels
data(veteran, package = "randomForestSRC")
vd <- veteran
vd$celltype=factor(vd$celltype)
vd$diagtime=factor(vd$diagtime)
vd.obj <- rfsrc(Surv(time,status)~., vd, ntree = 100, nodesize = 5)
plot(get.tree(vd.obj, 3))

## classification
iris.obj <- rfsrc(Species ~., data = iris)
plot(get.tree(iris.obj, 25, class.type = "bayes"))
plot(get.tree(iris.obj, 25, target = "setosa"))
plot(get.tree(iris.obj, 25, target = "versicolor"))
plot(get.tree(iris.obj, 25, target = "virginica"))

## ------------------------------------------------------------
## simple example of VIMP using iris classification
## ------------------------------------------------------------

## directly from trained forest
print(rfsrc(Species~.,iris,importance=TRUE)$importance)

## VIMP (and performance) use misclassification error by default
## but brier prediction error can be requested
print(rfsrc(Species~.,iris,importance=TRUE,perf.type="brier")$importance)

## example using vimp function (see vimp help file for details)
iris.obj <- rfsrc(Species ~., data = iris)
print(vimp(iris.obj)$importance)
print(vimp(iris.obj,perf.type="brier")$importance)

## example using hold out vimp (see holdout.vimp help file for details)
print(holdout.vimp(Species~.,iris)$importance)
print(holdout.vimp(Species~.,iris,perf.type="brier")$importance)

## ------------------------------------------------------------
## confidence interval for vimp using subsampling
## compare with holdout vimp
## ------------------------------------------------------------

## new York air quality measurements
o <- rfsrc(Ozone ~ ., data = airquality)
so <- subsample(o)
plot(so)

## compare with holdout vimp
print(holdout.vimp(Ozone ~ ., data = airquality)$importance)


##------------------------------------------------------------
## example of imputation in survival analysis
##------------------------------------------------------------

data(pbc, package = "randomForestSRC")
pbc.obj2 <- rfsrc(Surv(days, status) ~ ., pbc, na.action = "na.impute")

## same as above but iterate the missing data algorithm
pbc.obj3 <- rfsrc(Surv(days, status) ~ ., pbc,
         na.action = "na.impute", nimpute = 3)

## fast way to impute data (no inference is done)
## see impute for more details
pbc.imp <- impute(Surv(days, status) ~ ., pbc, splitrule = "random")

##------------------------------------------------------------
## compare RF-SRC to Cox regression
## Illustrates C-error and Brier score measures of performance
## assumes "pec" and "survival" libraries are loaded
##------------------------------------------------------------

if (library("survival", logical.return = TRUE)
    & library("pec", logical.return = TRUE)
    & library("prodlim", logical.return = TRUE))

{
  ##prediction function required for pec
  predictSurvProb.rfsrc <- function(object, newdata, times, ...){
    ptemp <- predict(object,newdata=newdata,...)$survival
    pos <- sindex(jump.times = object$time.interest, eval.times = times)
    p <- cbind(1,ptemp)[, pos + 1]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop("Prediction failed")
    p
  }

  ## data, formula specifications
  data(pbc, package = "randomForestSRC")
  pbc.na <- na.omit(pbc)  ##remove NA's
  surv.f <- as.formula(Surv(days, status) ~ .)
  pec.f <- as.formula(Hist(days,status) ~ 1)

  ## run cox/rfsrc models
  ## for illustration we use a small number of trees
  cox.obj <- coxph(surv.f, data = pbc.na, x = TRUE)
  rfsrc.obj <- rfsrc(surv.f, pbc.na, ntree = 150)

  ## compute bootstrap cross-validation estimate of expected Brier score
  ## see Mogensen, Ishwaran and Gerds (2012) Journal of Statistical Software
  set.seed(17743)
  prederror.pbc <- pec(list(cox.obj,rfsrc.obj), data = pbc.na, formula = pec.f,
                        splitMethod = "bootcv", B = 50)
  print(prederror.pbc)
  plot(prederror.pbc)

  ## compute out-of-bag C-error for cox regression and compare to rfsrc
  rfsrc.obj <- rfsrc(surv.f, pbc.na)
  cat("out-of-bag Cox Analysis ...", "\n")
  cox.err <- sapply(1:100, function(b) {
    if (b\%\%10 == 0) cat("cox bootstrap:", b, "\n")
    train <- sample(1:nrow(pbc.na), nrow(pbc.na), replace = TRUE)
    cox.obj <- tryCatch({coxph(surv.f, pbc.na[train, ])}, error=function(ex){NULL})
    if (!is.null(cox.obj)) {
      get.cindex(pbc.na$days[-train], pbc.na$status[-train], predict(cox.obj, pbc.na[-train, ]))
    } else NA
  })
  cat("\n\tOOB error rates\n\n")
  cat("\tRSF            : ", rfsrc.obj$err.rate[rfsrc.obj$ntree], "\n")
  cat("\tCox regression : ", mean(cox.err, na.rm = TRUE), "\n")
}

##------------------------------------------------------------
## competing risks
##------------------------------------------------------------

## WIHS analysis
## cumulative incidence function (CIF) for HAART and AIDS stratified by IDU

data(wihs, package = "randomForestSRC")
wihs.obj <- rfsrc(Surv(time, status) ~ ., wihs, nsplit = 3, ntree = 100)
plot.competing.risk(wihs.obj)
cif <- wihs.obj$cif.oob
Time <- wihs.obj$time.interest
idu <- wihs$idu
cif.haart <- cbind(apply(cif[,,1][idu == 0,], 2, mean),
                   apply(cif[,,1][idu == 1,], 2, mean))
cif.aids  <- cbind(apply(cif[,,2][idu == 0,], 2, mean),
                   apply(cif[,,2][idu == 1,], 2, mean))
matplot(Time, cbind(cif.haart, cif.aids), type = "l",
        lty = c(1,2,1,2), col = c(4, 4, 2, 2), lwd = 3,
        ylab = "Cumulative Incidence")
legend("topleft",
       legend = c("HAART (Non-IDU)", "HAART (IDU)", "AIDS (Non-IDU)", "AIDS (IDU)"),
       lty = c(1,2,1,2), col = c(4, 4, 2, 2), lwd = 3, cex = 1.5)


## illustrates the various splitting rules
## illustrates event specific and non-event specific variable selection
if (library("survival", logical.return = TRUE)) {

  ## use the pbc data from the survival package
  ## events are transplant (1) and death (2)
  data(pbc, package = "survival")
  pbc$id <- NULL

  ## modified Gray's weighted log-rank splitting
  ## (equivalent to cause=c(1,1) and splitrule="logrankCR")
  pbc.cr <- rfsrc(Surv(time, status) ~ ., pbc)
 
  ## log-rank cause-1 specific splitting and targeted VIMP for cause 1
  pbc.log1 <- rfsrc(Surv(time, status) ~ ., pbc, 
              splitrule = "logrankCR", cause = c(1,0), importance = TRUE)

  ## log-rank cause-2 specific splitting and targeted VIMP for cause 2
  pbc.log2 <- rfsrc(Surv(time, status) ~ ., pbc, 
              splitrule = "logrankCR", cause = c(0,1), importance = TRUE)

  ## extract VIMP from the log-rank forests: event-specific
  ## extract minimal depth from the Gray log-rank forest: non-event specific
  var.perf <- data.frame(md = max.subtree(pbc.cr)$order[, 1],
                         vimp1 = 100 * pbc.log1$importance[ ,1],
                         vimp2 = 100 * pbc.log2$importance[ ,2])
  print(var.perf[order(var.perf$md), ], digits = 2)

}

## ------------------------------------------------------------
## regression analysis
## ------------------------------------------------------------

## new York air quality measurements
airq.obj <- rfsrc(Ozone ~ ., data = airquality, na.action = "na.impute")

# partial plot of variables (see plot.variable for more details)
plot.variable(airq.obj, partial = TRUE, smooth.lines = TRUE)

## motor trend cars
mtcars.obj <- rfsrc(mpg ~ ., data = mtcars)

## ------------------------------------------------------------
## regression with custom bootstrap
## ------------------------------------------------------------

ntree <- 25
n <- nrow(mtcars)
s.size <- n / 2
swr <- TRUE
samp <- randomForestSRC:::make.sample(ntree, n, s.size, swr)
o <- rfsrc(mpg ~ ., mtcars, bootstrap = "by.user", samp = samp)

## ------------------------------------------------------------
## classification analysis
## ------------------------------------------------------------

## iris data
iris.obj <- rfsrc(Species ~., data = iris)

## wisconsin prognostic breast cancer data
data(breast, package = "randomForestSRC")
breast.obj <- rfsrc(status ~ ., data = breast, block.size=1)
plot(breast.obj)

## ------------------------------------------------------------
## big data set, reduce number of variables using simple method
## ------------------------------------------------------------

## use Iowa housing data set
data(housing, package = "randomForestSRC")

## original data contains lots of missing data, use fast imputation
## however see impute for other methods
housing2 <- impute(data = housing, fast = TRUE)

## run shallow trees to find variables that split any tree
xvar.used <- rfsrc(SalePrice ~., housing2, ntree = 250, nodedepth = 4,
                   var.used="all.trees", mtry = Inf, nsplit = 100)$var.used

## now fit forest using filtered variables
xvar.keep  <- names(xvar.used)[xvar.used >= 1]
o <- rfsrc(SalePrice~., housing2[, c("SalePrice", xvar.keep)])
print(o)

## ------------------------------------------------------------
## imbalanced classification data
## see the "imbalanced" function for further details
##
## (a) use balanced random forests with undersampling of the majority class
## Specifically let n0, n1 be sample sizes for majority, minority
## cases.  We sample 2 x n1 cases with majority, minority cases chosen
## with probabilities n1/n, n0/n where n=n0+n1
##
## (b) balanced random forests using "imbalanced"
##
## (c) q-classifier (RFQ) using "imbalanced" 
##
## ------------------------------------------------------------

## Wisconsin breast cancer example
data(breast, package = "randomForestSRC")
breast <- na.omit(breast)

## balanced random forests - brute force
y <- breast$status
obdirect <- rfsrc(status ~ ., data = breast, nsplit = 10,
            case.wt = randomForestSRC:::make.wt(y),
            sampsize = randomForestSRC:::make.size(y))
print(obdirect)
print(get.imbalanced.performance(obdirect))

## balanced random forests - using "imbalanced" 
ob <- imbalanced(status ~ ., data = breast, method = "brf")
print(ob)
print(get.imbalanced.performance(ob))

## q-classifier (RFQ) - using "imbalanced" 
oq <- imbalanced(status ~ ., data = breast)
print(oq)
print(get.imbalanced.performance(oq))

## q-classifier (RFQ) - with auc splitting
oqauc <- imbalanced(status ~ ., data = breast, splitrule = "auc")
print(oqauc)
print(get.imbalanced.performance(oqauc))

## ------------------------------------------------------------
## unsupervised analysis
## ------------------------------------------------------------

## two equivalent ways to implement unsupervised forests
mtcars.unspv <- rfsrc(Unsupervised() ~., data = mtcars)
mtcars2.unspv <- rfsrc(data = mtcars)


## illustration of sidClustering for the mtcars data
## see sidClustering for more details
mtcars.sid <- sidClustering(mtcars, k = 1:10)
print(split(mtcars, mtcars.sid$cl[, 3]))
print(split(mtcars, mtcars.sid$cl[, 10]))


## ------------------------------------------------------------
## bivariate regression using Mahalanobis splitting
## also illustrates user specified covariance matrix
## ------------------------------------------------------------

if (library("mlbench", logical.return = TRUE)) {

  ## load boston housing data, specify the bivariate regression
  data(BostonHousing)
  f <- formula("Multivar(lstat, nox) ~.")
  
  ## Mahalanobis splitting  
  bh.mreg <- rfsrc(f, BostonHousing, importance = TRUE, splitrule = "mahal")
  
  ## performance error and vimp
  vmp <- get.mv.vimp(bh.mreg)
  pred <- get.mv.predicted(bh.mreg)

  ## standardized error and vimp
  err.std <- get.mv.error(bh.mreg, standardize = TRUE)
  vmp.std <- get.mv.vimp(bh.mreg, standardize = TRUE)

  ## same analysis, but with user specified covariance matrix
  sigma <- cov(BostonHousing[, c("lstat","nox")])
  bh.mreg2 <- rfsrc(f, BostonHousing, splitrule = "mahal", sigma = sigma)
  
}

## ------------------------------------------------------------
## multivariate mixed forests (nutrigenomic study)
## study effects of diet, lipids and gene expression for mice
## diet, genotype and lipids used as the multivariate y
## genes used for the x features
## ------------------------------------------------------------

## load the data (data is a list)
data(nutrigenomic, package = "randomForestSRC")

## assemble the multivariate y data
ydta <- data.frame(diet = nutrigenomic$diet,
                   genotype = nutrigenomic$genotype,
                   nutrigenomic$lipids)

## multivariate mixed forest call
## uses "get.mv.formula" for conveniently setting formula
mv.obj <- rfsrc(get.mv.formula(colnames(ydta)),
                data.frame(ydta, nutrigenomic$genes),
		importance=TRUE, nsplit = 10)

## print results for diet and genotype y values	    
print(mv.obj, outcome.target = "diet")
print(mv.obj, outcome.target = "genotype")

## extract standardized VIMP
svimp <- get.mv.vimp(mv.obj, standardize = TRUE)

## plot standardized VIMP for diet, genotype and lipid for each gene
boxplot(t(svimp), col = "bisque", cex.axis = .7, las = 2,
        outline = FALSE,
        ylab = "standardized VIMP",
        main = "diet/genotype/lipid VIMP for each gene")
	

## ------------------------------------------------------------
## illustrates yvar.wt which sets the probability of selecting 
## the response variables in multivariate regression
## ------------------------------------------------------------

## use mtcars: add fake responses
mult.mtcars  <- cbind(mtcars, mtcars$mpg, mtcars$mpg)
names(mult.mtcars) = c(names(mtcars), "mpg2", "mpg3")

## noise up the fake responses
mult.mtcars$mpg2 <- sample(mtcars$mpg)
mult.mtcars$mpg3 <- sample(mtcars$mpg)

formula = as.formula(Multivar(mpg, mpg2, mpg3)  ~ .)

## select 2 of the 3 responses randomly at each split with an associated weight vector.
## choose the noisy y responses which should degrade performance
yvar.wt = c(0.000001, 0.5, 0.5)
ytry = 2

mult.grow <- rfsrc(formula = formula, data = mult.mtcars, ytry = ytry, yvar.wt = yvar.wt)

print(mult.grow)
print(get.mv.error(mult.grow))

## Also, compare the following two results, as they should be similar:
yvar.wt = c(1.0, 00000.1, 00000.1)
ytry = 1

result1 = rfsrc(formula = formula, data = mult.mtcars, ytry = ytry, yvar.wt = yvar.wt)
result2 = rfsrc(mpg ~ ., mtcars)

print(get.mv.error(result1))
print(get.mv.error(result2))

## ------------------------------------------------------------
## custom splitting using the pre-coded examples
## ------------------------------------------------------------

## motor trend cars
mtcars.obj <- rfsrc(mpg ~ ., data = mtcars, splitrule = "custom")

## iris analysis
iris.obj <- rfsrc(Species ~., data = iris, splitrule = "custom1")

## WIHS analysis
wihs.obj <- rfsrc(Surv(time, status) ~ ., wihs, nsplit = 3,
                  ntree = 100, splitrule = "custom1")

}
}
\keyword{forest}
