\name{predict.rfsrc}
\alias{predict.rfsrc}
\title{Prediction for Random Forests for Survival, Regression, and Classification}
\description{
  Obtain predicted values using a forest.  Also returns performance
  values if the test data contains y-outcomes.
}
\usage{\method{predict}{rfsrc}(object,
  newdata,
  m.target = NULL,
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
  do.trace = FALSE, membership = FALSE, statistics = FALSE,
   
  ...)
}
\arguments{
  \item{object}{An object of class \code{(rfsrc, grow)} or \code{(rfsrc,
	    forest)}.}

  \item{newdata}{Test data. If missing, the original grow (training)
    data is used.}
  

  \item{m.target}{Character vector for multivariate families
    specifying the target outcomes to be used. The default uses all
    coordinates.}
  
  \item{importance}{Method used for variable importance (VIMP).  Also
    see \command{vimp} for more flexibility, including joint vimp
    calculations.  See \command{holdoutvimp} for an alternate importance
    measure.}
  
  \item{get.tree}{Vector of integer(s) identifying trees over which the
    ensembles are calculated over.  By default, uses all trees in the
    forest.  As an example, the user can extract the ensemble, the VIMP
    , or proximity from a single tree (or several trees).  Note that
    \code{block.size} will be over-ridden so that it is no larger than
    the requested number of trees.  See example below illustrating how
    to extract VIMP for each tree.}
  
  \item{block.size}{Should the error rate be calculated on every tree?
    When \code{NULL}, it will only be calculated on the last tree.  To
    view the error rate on every nth tree, set the value to an integer
    between \code{1} and \code{ntree}.  If importance is requested, VIMP
    is calculated in "blocks" of size equal to \code{block.size}, thus
    resulting in a compromise between ensemble and permutation VIMP.}
  
  \item{na.action}{Missing value action. The default \code{na.omit}
    removes the entire record if any entry is \code{NA}.  
    Selecting \option{na.random} uses fast random imputation, while
    \option{na.impute} uses the imputation method described in
    \command{rfsrc}.}
  
  \item{outcome}{Determines whether the y-outcomes from the training
    data or the test data are used to calculate the predicted value.
    The default and natural choice is \code{train} which uses the
    original training data.  Option is ignored when \code{newdata} is
    missing as the training data is used for the test data in such
    settings.  The option is also ignored whenever the test data is
    devoid of y-outcomes.  See the details and examples below for more
    information.}

  \item{perf.type}{Optional character value for requesting metric used
    for predicted value, variable importance (VIMP) and error rate.  If not
    specified, values returned are calculated by the default action used
    for the family. Currently applicable only to classification and
    multivariate classification; allowed values are
    \code{perf.type="misclass"} (default), \code{perf.type="brier"} and
    \code{perf.type="gmean"}.}

  \item{proximity}{Should proximity between test observations
    be calculated?  Possible choices are \code{"inbag"}, \code{"oob"},
    \code{"all"}, \code{TRUE}, or \code{FALSE} --- but some options may
    not be valid and will depend on the context of the predict call.
    The safest choice is \code{TRUE} if proximity is desired.}
  
  \item{distance}{Should distance between test observations
    be calculated?  Possible choices are \code{"inbag"}, \code{"oob"},
    \code{"all"}, \code{TRUE}, or \code{FALSE} --- but some options may
    not be valid and will depend on the context of the predict call.
    The safest choice is \code{TRUE} if distance is desired.}
  
  \item{forest.wt}{Should the forest weight matrix for test observations
    be calculated?  Choices are the same as proximity.}

  \item{ptn.count}{The number of terminal nodes that each tree in the
    grow forest should be pruned back to.  The terminal node membership
    for the pruned forest is returned but no other action is taken.  The
    default is \code{ptn.count=0} which does no pruning.}

  \item{var.used}{Record the number of times a variable is split?}
  
  \item{split.depth}{Return minimal depth for each variable for each case?}

  \item{case.depth}{Return a matrix recording the depth at which a case
    first splits in a tree. Default is \code{FALSE}.}
  
  \item{seed}{Negative integer specifying seed for the random number
    generator.} 

  \item{do.trace}{Number of seconds between updates to the user on
    approximate time to completion.}

  \item{membership}{Should terminal node membership and inbag
    information be returned?}

  \item{statistics}{Should split statistics be returned?  Values can be
    parsed using \command{stat.split}.}

   
  

  \item{...}{Further arguments passed to or from other methods.}

}
\details{

  Predicted values are obtained by "dropping" test data down the trained
  forest (forest calculated using training data).  Performance values are
  returned if test data contains y-outcome values.  Single as well as
  joint VIMP are also returned if requested.
  
  If no test data is provided, the original training data is used, and
  the code reverts to restore mode allowing the user to restore the
  original trained forest.  This feature allows extracting outputs from
  the forest not asked for in the original grow call.
  
  If \option{outcome="test"}, the predictor is calculated by using
  y-outcomes from the test data (outcome information must be present).
  Terminal nodes from the trained forest are recalculated using y-outcomes
  from the test set.  This yields a modified predictor in which the
  topology of the forest is based solely on the training data, but where
  predicted values are obtained from test data.  Error rates and VIMP
  are calculated by bootstrapping the test data and using out-of-bagging
  to ensure unbiased estimates.  

  \code{csv=TRUE} returns case specific VIMP;  \code{cse=TRUE} returns
  case specific error rates.  Applies to all families except survival.
  These options can also be applied while training.
  
}    
\value{
  An object of class \code{(rfsrc, predict)}, which is a list with the
  following components:
  \item{call}{The original grow call to \command{rfsrc}.}
  \item{family}{The family used in the analysis.}
  \item{n}{Sample size of test data (depends upon \code{NA} values).}
  \item{ntree}{Number of trees in the grow forest.}
  \item{yvar}{Test set y-outcomes or original grow y-outcomes if none.}
  \item{yvar.names}{A character vector of the y-outcome names.}
  \item{xvar}{Data frame of test set x-variables.}
  \item{xvar.names}{A character vector of the x-variable names.}
  \item{leaf.count}{Number of terminal nodes for each tree in the
    grow forest.  Vector of length \code{ntree}.}
  \item{proximity}{Symmetric proximity matrix of the test data.}
  \item{forest}{The grow forest.}
  
                    
  
  \item{membership}{Matrix recording terminal node membership for the
    test data where each column contains the node number that a
    case falls in for that tree.}
  \item{inbag}{Matrix recording inbag membership for the test data
    where each column contains the number of times that a case
    appears in the bootstrap sample for that tree.}
  \item{var.used}{Count of the number of times a variable was used in
    growing the forest.}
  \item{imputed.indv}{Vector of indices of records in test data with
    missing values.}
  \item{imputed.data}{Data frame comprising imputed test data.  The first
    columns are the y-outcomes followed by the x-variables.}
  \item{split.depth}{Matrix (i,j) or array (i,j,k) recording the
    minimal depth for variable j for case i, either averaged over
    the forest, or by tree k.}
  \item{node.stats}{Split statistics returned when
    \code{statistics=TRUE} which can be parsed using \command{stat.split}.}
  \item{err.rate}{Cumulative OOB error rate for the test data if
    y-outcomes are present.}
  \item{importance}{Test set variable importance (VIMP).  Can be
    \code{NULL}.}
  \item{predicted}{Test set predicted value.}
  \item{predicted.oob}{OOB predicted value (\code{NULL} unless
    \option{outcome="test"}).}\cr
  \item{quantile}{Quantile value at probabilities requested.}
  \item{quantile.oob}{OOB quantile value at probabilities requested (\code{NULL} unless
    \option{outcome="test"}).}\cr

  
  \item{++++++++}{for classification settings, additionally ++++++++} \cr     
  \item{class}{In-bag predicted class labels.}
  \item{class.oob}{OOB predicted class labels (\code{NULL} unless \option{outcome="test"}).}\cr
  
  \item{++++++++}{for multivariate settings, additionally ++++++++} \cr
  
  \item{regrOutput}{List containing performance values for
	  test multivariate regression responses (applies only in
	  multivariate settings).}
	\item{clasOutput}{List containing performance values for
	  test multivariate categorical (factor) responses (applies only in
	  multivariate settings).}
  
  \item{++++++++}{for survival settings, additionally ++++++++} \cr
  
  \item{chf}{Cumulative hazard function (CHF).}
  \item{chf.oob}{OOB CHF (\code{NULL} unless \option{outcome="test"}).}
  \item{survival}{Survival function.}
  \item{survival.oob}{OOB survival function (\code{NULL} unless \option{outcome="test"}).}
  \item{time.interest}{Ordered unique death times.}
  \item{ndead}{Number of deaths.}\cr
  
  \item{++++++++}{for competing risks, additionally ++++++++} \cr
  
  \item{chf}{Cause-specific cumulative hazard function (CSCHF)
    for each event.}
  \item{chf.oob}{OOB CSCHF for each event (\code{NULL} unless \option{outcome="test"}).}
  \item{cif}{Cumulative incidence function (CIF) for each event.}
  \item{cif.oob}{OOB CIF (\code{NULL} unless \option{outcome="test"}).}
}
\note{

  The dimensions and values of returned objects depend heavily on the
  underlying family and whether y-outcomes are present in the test data.
  In particular, items related to performance will be \code{NULL} when
  y-outcomes are not present.  For multivariate families, predicted
  values, VIMP, error rate, and performance values are stored in the
  lists \code{regrOutput} and \code{clasOutput} which can be extracted
  using functions \code{get.mv.error}, \code{get.mv.predicted} and
  \code{get.mv.vimp}.
  
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
  \command{\link{stat.split.rfsrc}},
  \command{\link{synthetic.rfsrc}},
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
## ... just for _fun_ ...
## survival analysis using mixed multivariate outcome analysis 
## compare the predicted value to RSF
## ------------------------------------------------------------

## train survival forest using pbc data
data(pbc, package = "randomForestSRC")
rsf.obj <- rfsrc(Surv(days, status) ~ ., pbc)
yvar <- rsf.obj$yvar

## fit a mixed outcome forest using days and status as y-variables
pbc.mod <- pbc
pbc.mod$status <- factor(pbc.mod$status)
mix.obj <- rfsrc(Multivar(days, status) ~., pbc.mod)

## compare oob predicted values
rsf.pred <- rsf.obj$predicted.oob
mix.pred <- mix.obj$regrOutput$days$predicted.oob
plot(rsf.pred, mix.pred)

## compare C-error rate
rsf.err <- get.cindex(yvar$days, yvar$status, rsf.pred)
mix.err <- 1 - get.cindex(yvar$days, yvar$status, mix.pred)
cat("RSF                :", rsf.err, "\n")
cat("multivariate forest:", mix.err, "\n")

}
}
\keyword{predict}
\keyword{forest}

