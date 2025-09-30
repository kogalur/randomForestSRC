\name{predict.rfsrc}
\alias{predict.rfsrc}
\title{Prediction for Random Forests for Survival, Regression, and Classification}
\description{
  Obtain predicted values using a forest.  Also returns performance
  values if the test data contains y-outcomes.
}
\usage{\method{predict}{rfsrc}(object,
  newdata,
  importance = c(FALSE, TRUE, "none", "anti", "permute", "random"),
  get.tree = NULL,
  block.size = if (any(is.element(as.character(importance),
                     c("none", "FALSE")))) NULL else 10,
  na.action = c("na.omit", "na.impute", "na.random"),
  outcome = c("train", "test"),
  perf.type = NULL,
  proximity = FALSE,
  forest.wt = FALSE,
  ptn.count = 0,
  distance = FALSE,
  var.used = c(FALSE, "all.trees", "by.tree"),
  split.depth = c(FALSE, "all.trees", "by.tree"),
  case.depth = FALSE,
  seed = NULL,
  do.trace = FALSE, membership = FALSE,
  marginal.xvar = NULL,  ...)
}

\arguments{
  \item{object}{An object of class \code{(rfsrc, grow)} or \code{(rfsrc, forest)}.}
  
  \item{newdata}{Test data. If omitted, the original training data is used.}
  
  \item{importance}{Method for computing variable importance (VIMP). See \code{\link{vimp}} for additional options including joint importance. See \code{\link{holdout.vimp}} for an alternative importance measure.}
  
  \item{get.tree}{Vector of integers specifying which trees to use for ensemble calculations. Defaults to all trees. Useful for extracting ensembles, VIMP, or proximity from specific trees. If specified, \code{block.size} is overridden to match the number of trees. See examples for per-tree VIMP extraction.}
  
  \item{block.size}{Controls the granularity of error rate and VIMP calculation. If \code{NULL}, error is reported only for the final tree. Set to an integer \code{k} to compute error every \code{k} trees. For VIMP, calculations are done in blocks of size \code{block.size}, balancing between tree-level and forest-level assessments.}
  
  \item{na.action}{Action to take when missing values are present. Options are \code{"na.omit"} (default), \code{"na.random"} for fast random imputation, or \code{"na.impute"} to use the imputation method in \code{rfsrc}.}
  
  \item{outcome}{Specifies whether predicted values should be based on the outcomes from the training data (\code{"train"}, default) or test data. Ignored if \code{newdata} is missing or if test outcomes are unavailable.}
  
  \item{perf.type}{Optional metric for prediction, VIMP, and error. Currently used for classification and multivariate classification. Choices: \code{"misclass"} (default), \code{"brier"}, and \code{"gmean"}.}
  
  \item{proximity}{Whether to compute the proximity matrix for test observations. Options include \code{"inbag"}, \code{"oob"}, \code{"all"}, \code{TRUE}, or \code{FALSE}. Not all options are valid in all contexts; \code{TRUE} is the safest choice.}
  
  \item{distance}{Whether to compute the distance matrix. Options are the same as for \code{proximity}.}
  
  \item{forest.wt}{Whether to compute the forest weight matrix. Options are the same as for \code{proximity}.}
  
  \item{ptn.count}{If nonzero, each tree is pruned to have this many terminal nodes. Only the terminal node membership is returned; no prediction is made. Default is \code{ptn.count = 0}.}
  
  \item{var.used}{If \code{TRUE}, records how many times each variable was used for splitting.}
  
  \item{split.depth}{If \code{TRUE}, returns minimal depth of each variable per case.}
  
  \item{case.depth}{If \code{TRUE}, returns a matrix of the depth at which each case first splits in each tree.}
  
  \item{seed}{Negative integer used to set the random seed.}
  
  \item{do.trace}{Number of seconds between progress updates during execution.}
  
  \item{membership}{If \code{TRUE}, returns terminal node membership and in-bag information.}
  
  \item{marginal.xvar}{Vector of variable names to marginalize over when calculating weights or proximity. If a variable is marginalized, its split does not partition the data; all cases are passed to both daughters. When all splits involve marginalized variables, terminal nodes contain the full dataset. When no marginalized variables are used, membership is unchanged. Default is \code{NULL} (no marginalization).}
  
  \item{...}{Additional arguments passed to or from other methods.}
}


\details{

Predicted values are obtained by "dropping" the test data down the trained forest-i.e., the forest grown using the training data. If the test data includes y-outcome values, performance metrics are also returned. Variable importance (VIMP), including joint VIMP, is returned if requested.

If no test data is supplied, the function uses the original training data and enters "restore" mode. This allows users to extract outputs from the trained forest that were not requested during the original grow call.

If \code{outcome = "test"}, predictions are computed using y-outcomes from the test data (which must include outcome values). Terminal node statistics are recalculated using these outcomes, while the tree topology remains fixed from training. Error rates and VIMP are then computed by bootstrapping the test set and applying out-of-bagging to maintain unbiased estimates.

Set \code{csv = TRUE} to return case-specific VIMP, and \code{cse = TRUE} to return case-specific error rates. These apply to all families except survival. Both options can also be used at training time.
}

\value{
An object of class \code{(rfsrc, predict)}, which is a list with the following components:

\describe{
  \item{call}{The original grow call to \code{rfsrc}.}
  \item{family}{The family used in the analysis.}
  \item{n}{Sample size of the test data (after handling missing values).}
  \item{ntree}{Number of trees in the trained forest.}
  \item{yvar}{Y-outcome values from the test data or original grow data (if \code{newdata} is missing).}
  \item{yvar.names}{Character vector of response variable names.}
  \item{xvar}{Data frame of test set predictor variables.}
  \item{xvar.names}{Character vector of predictor variable names.}
  \item{leaf.count}{Vector of length \code{ntree} giving the number of terminal nodes per tree.}
  \item{proximity}{Proximity matrix computed on the test data.}
  \item{forest}{The trained forest object.}
  \item{forest.wt}{Forest weight matrix for test cases.}
  \item{ptn.membership}{Matrix of pruned terminal node membership. Only returned if \code{ptn.count > 0}.}
  \item{membership}{Matrix of terminal node membership for test cases. Each column corresponds to one tree.}
  \item{inbag}{Matrix indicating how many times each case appears in the bootstrap sample for each tree.}
  \item{var.used}{Number of times each variable was used in splitting.}
  \item{imputed.indv}{Indices of test observations with missing values.}
  \item{imputed.data}{Imputed version of the test data. Columns are ordered with responses first, followed by predictors.}
  \item{split.depth}{Matrix or array recording minimal depth of each variable for each case, optionally by tree.}
  \item{err.rate}{Cumulative out-of-bag (OOB) error rate, if y-outcomes are present.}
  \item{importance}{Variable importance (VIMP) for the test data. May be \code{NULL}.}
  \item{predicted}{Predicted values for the test data.}
  \item{predicted.oob}{OOB predicted values. May be \code{NULL} depending on context.}
  \item{quantile}{Estimated quantile values at the requested probabilities (quantile regression only).}
  \item{quantile.oob}{OOB quantile values. May be \code{NULL}.}

  \item{class}{(Classification only) Predicted class labels.}
  \item{class.oob}{(Classification only) OOB predicted class labels.}

  \item{regrOutput}{(Multivariate only) List of performance measures for multivariate regression outcomes.}
  \item{clasOutput}{(Multivariate only) List of performance measures for multivariate categorical outcomes.}

  \item{chf}{(Survival or competing risks) Cumulative hazard function (CHF).}
  \item{chf.oob}{(Survival or competing risks) OOB CHF. May be \code{NULL}.}
  \item{survival}{(Survival only) Survival function estimates.}
  \item{survival.oob}{(Survival only) OOB survival function. May be \code{NULL}.}
  \item{time.interest}{(Survival or competing risks) Sorted unique event times.}
  \item{ndead}{(Survival or competing risks) Number of deaths observed.}

  \item{cif}{(Competing risks only) Cumulative incidence function (CIF) for each event type.}
  \item{cif.oob}{(Competing risks only) OOB CIF. May be \code{NULL}.}
  \item{chf}{(Competing risks only) Cause-specific cumulative hazard function (CSCHF).}
  \item{chf.oob}{(Competing risks only) OOB CSCHF. May be \code{NULL}.}
}
}


\note{

The dimensions and contents of returned objects depend on the forest family and whether y-outcomes are available in the test data. In particular, performance-related components (e.g., error rate, VIMP) will be \code{NULL} if y-outcomes are missing.

For multivariate families, predicted values, VIMP, error rates, and performance metrics are stored in the lists \code{regrOutput} and \code{clasOutput}. These can be accessed using the helper functions \code{\link{get.mv.predicted}}, \code{\link{get.mv.vimp}}, and \code{\link{get.mv.error}}.

}


\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.
  
  Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
  (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.
  
  Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7(2):25-31.
}	
\seealso{
  \command{\link{holdout.vimp.rfsrc}},
  \command{\link{plot.competing.risk.rfsrc}},
  \command{\link{plot.rfsrc}},
  \command{\link{plot.survival.rfsrc}},
  \command{\link{plot.variable.rfsrc}},
  \command{\link{rfsrc}},
  \command{\link{rfsrc.fast}},
  \command{\link{vimp.rfsrc}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## typical train/testing scenario
## ------------------------------------------------------------

data(veteran, package = "randomForestSRC")
train <- sample(1:nrow(veteran), round(nrow(veteran) * 0.80))
veteran.grow <- rfsrc(Surv(time, status) ~ ., veteran[train, ]) 
veteran.pred <- predict(veteran.grow, veteran[-train, ])
print(veteran.grow)
print(veteran.pred)


## ------------------------------------------------------------
## restore mode
## - if predict is called without specifying the test data
##   the original training data is used and the forest is restored
## ------------------------------------------------------------

## first train the forest
airq.obj <- rfsrc(Ozone ~ ., data = airquality)

## now we restore it and compare it to the original call
## they are identical
predict(airq.obj)
print(airq.obj)

## we can retrieve various outputs that were not asked for in
## in the original call

## here we extract the proximity matrix
prox <- predict(airq.obj, proximity = TRUE)$proximity
print(prox[1:10,1:10])

## here we extract the number of times a variable was used to grow
## the grow forest
var.used <- predict(airq.obj, var.used = "by.tree")$var.used
print(head(var.used))

## ------------------------------------------------------------
## prediction when test data has missing values
## ------------------------------------------------------------

data(pbc, package = "randomForestSRC")
trn <- pbc[1:312,]
tst <- pbc[-(1:312),]
o <- rfsrc(Surv(days, status) ~ ., trn)

## default imputation method used by rfsrc
print(predict(o, tst, na.action = "na.impute"))

## random imputation
print(predict(o, tst, na.action = "na.random"))

## ------------------------------------------------------------
## requesting different performance for classification
## ------------------------------------------------------------

## default performance is misclassification
o <- rfsrc(Species~., iris)
print(o)

## get (normalized) brier performance
print(predict(o, perf.type = "brier"))

## ------------------------------------------------------------
## vimp for each tree: illustrates get.tree 
## ------------------------------------------------------------

## regression analysis but no VIMP
o <- rfsrc(mpg~., mtcars)

## now extract VIMP for each tree using get.tree
vimp.tree <- do.call(rbind, lapply(1:o$ntree, function(b) {
     predict(o, get.tree = b, importance = TRUE)$importance
}))

## boxplot of tree VIMP
boxplot(vimp.tree, outline = FALSE, col = "cyan")
abline(h = 0, lty = 2, col = "red")

## summary information of tree VIMP
print(summary(vimp.tree))

## extract tree-averaged VIMP using importance=TRUE
## remember to set block.size to 1
print(predict(o, importance = TRUE, block.size = 1)$importance)

## use direct call to vimp() for tree-averaged VIMP
print(vimp(o, block.size = 1)$importance)

## ------------------------------------------------------------
## vimp for just a few trees
## illustrates how to get vimp if you have a large data set
## ------------------------------------------------------------

## survival analysis but no VIMP
data(pbc, package = "randomForestSRC")
o <- rfsrc(Surv(days, status) ~ ., pbc, ntree = 2000)

## get vimp for a small number of trees
print(predict(o, get.tree=1:250, importance = TRUE)$importance)


## ------------------------------------------------------------
## case-specific vimp
## returns VIMP for each case
## ------------------------------------------------------------

o <- rfsrc(mpg~., mtcars)
op <- predict(o, importance = TRUE, csv = TRUE)
csvimp <- get.mv.csvimp(op, standardize=TRUE)
print(csvimp)

## ------------------------------------------------------------
## case-specific error rate
## returns tree-averaged error rate for each case
## ------------------------------------------------------------

o <- rfsrc(mpg~., mtcars)
op <- predict(o, importance = TRUE, cse = TRUE)
cserror <- get.mv.cserror(op, standardize=TRUE)
print(cserror)


## ------------------------------------------------------------
## predicted probability and predicted class labels are returned
## in the predict object for classification analyses
## ------------------------------------------------------------

data(breast, package = "randomForestSRC")
breast.obj <- rfsrc(status ~ ., data = breast[(1:100), ])
breast.pred <- predict(breast.obj, breast[-(1:100), ])
print(head(breast.pred$predicted))
print(breast.pred$class)


## ------------------------------------------------------------
## unique feature of randomForestSRC
## cross-validation can be used when factor labels differ over
## training and test data
## ------------------------------------------------------------

## first we convert all x-variables to factors
data(veteran, package = "randomForestSRC")
veteran2 <- data.frame(lapply(veteran, factor))
veteran2$time <- veteran$time
veteran2$status <- veteran$status

## split the data into unbalanced train/test data (25/75)
## the train/test data have the same levels, but different labels
train <- sample(1:nrow(veteran2), round(nrow(veteran2) * .25))
summary(veteran2[train,])
summary(veteran2[-train,])

## train the forest and use this to predict on test data
o.grow <- rfsrc(Surv(time, status) ~ ., veteran2[train, ]) 
o.pred <- predict(o.grow, veteran2[-train , ])
print(o.grow)
print(o.pred)

## even harder ... factor level not previously encountered in training
veteran3 <- veteran2[1:3, ]
veteran3$celltype <- factor(c("newlevel", "1", "3"))
o2.pred <- predict(o.grow, veteran3)
print(o2.pred)
## the unusual level is treated like a missing value but is not removed
print(o2.pred$xvar)

## ------------------------------------------------------------
## example illustrating the flexibility of outcome = "test"
## illustrates restoration of forest via outcome = "test"
## ------------------------------------------------------------

## first we train the forest
data(pbc, package = "randomForestSRC")
pbc.grow <- rfsrc(Surv(days, status) ~ ., pbc)

## use predict with outcome = TEST
pbc.pred <- predict(pbc.grow, pbc, outcome = "test")

## notice that error rates are the same!!
print(pbc.grow)
print(pbc.pred)

## note this is equivalent to restoring the forest
pbc.pred2 <- predict(pbc.grow)
print(pbc.grow)
print(pbc.pred)
print(pbc.pred2)

## similar example, but with na.action = "na.impute"
airq.obj <- rfsrc(Ozone ~ ., data = airquality, na.action = "na.impute")
print(airq.obj)
print(predict(airq.obj))
## ... also equivalent to outcome="test" but na.action = "na.impute" required
print(predict(airq.obj, airquality, outcome = "test", na.action = "na.impute"))

## classification example
iris.obj <- rfsrc(Species ~., data = iris)
print(iris.obj)
print(predict.rfsrc(iris.obj, iris, outcome = "test"))

## ------------------------------------------------------------
## another example illustrating outcome = "test"
## unique way to check reproducibility of the forest
## ------------------------------------------------------------

## training step
set.seed(542899)
data(pbc, package = "randomForestSRC")
train <- sample(1:nrow(pbc), round(nrow(pbc) * 0.50))
pbc.out <- rfsrc(Surv(days, status) ~ .,  data=pbc[train, ])

## standard prediction call
pbc.train <- predict(pbc.out, pbc[-train, ], outcome = "train")
##non-standard predict call: overlays the test data on the grow forest
pbc.test <- predict(pbc.out, pbc[-train, ], outcome = "test")

## check forest reproducibilility by comparing "test" predicted survival
## curves to "train" predicted survival curves for the first 3 individuals
Time <- pbc.out$time.interest
matplot(Time, t(pbc.train$survival[1:3,]), ylab = "Survival", col = 1, type = "l")
matlines(Time, t(pbc.test$survival[1:3,]), col = 2)

## ------------------------------------------------------------
## multivariate forest example
## ------------------------------------------------------------

## train the forest
trn <- 1:20
o <- rfsrc(cbind(mpg, disp)~.,mtcars[trn,])

## print training results for each outcome
print(o, outcome.target="mpg")
print(o, outcome.target="disp")

## print test results for each outcome
p <- predict(o, mtcars[-trn,])
print(p, outcome.target="mpg")
print(p, outcome.target="disp")


}
}
\keyword{predict}
\keyword{forest}

