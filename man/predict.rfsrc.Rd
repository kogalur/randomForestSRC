\name{predict.rfsrc}
\alias{predict.rfsrc}
\title{Prediction for Random Forests for Survival, Regression, and Classification}
\description{
  Obtain predicted values using a forest.  Also returns performance
  values if the test data contains y-outcomes.
}
\usage{\method{predict}{rfsrc}(object, newdata, outcome.target = NULL,
  importance = c(FALSE, TRUE, "none", "permute", "random", "anti",
                 "permute.ensemble", "random.ensemble", "anti.ensemble")[1],
  na.action = c("na.omit", "na.impute"), outcome = c("train", "test"),
  proximity = FALSE,

  

  var.used = c(FALSE, "all.trees", "by.tree"),
  split.depth = c(FALSE, "all.trees", "by.tree"), seed = NULL,
  do.trace = FALSE, membership = FALSE, statistics = FALSE, ...)
}
\arguments{
  \item{object}{An object of class \code{(rfsrc, grow)} or \code{(rfsrc,
	    forest)}.}
  \item{newdata}{Test data. If missing, the original grow
    (training) data is used.}
  \item{outcome.target}{Character vector for multivariate families
    specifying the target outcomes to be used. The default is to use all
    coordinates.}
  \item{importance}{Method for computing variable importance (VIMP).
    See \command{rfsrc} for details. Only applies when the test data
    contains y-outcome values.}
  \item{na.action}{Missing value action. The default
    \code{na.omit} removes the entire record if
    even one of its entries is \code{NA}.  Selecting \option{na.impute} 
    imputes the test data.}
  \item{outcome}{Determines whether the y-outcomes from the training
    data or the test data are used to calculate the predicted value.
    The default and natural choice is \code{train} which uses the
    original training data.  Option is ignored when \code{newdata} is
    missing as the training data is used for the test data in such
    settings.  The option is also ignored whenever the test data is
    devoid of y-outcomes.  See the details and examples below for more
    information.}
  \item{proximity}{Should proximity between test observations
    be calculated?  Possible choices are \code{"inbag"}, \code{"oob"},
    \code{"all"}, \code{TRUE}, or \code{FALSE} --- but some options may
    not be valid and will depend on the context of the predict call.
    The safest choice is \code{TRUE} if proximity is desired.}
  

  
  
  \item{var.used}{Record the number of times a variable is split?}
  \item{split.depth}{Return minimal depth for each variable for each case?}
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
  Predicted values are obtained by dropping test data down the grow
  forest (the forest grown using the training data).  The overall error
  rate and VIMP are also returned if the test data contains y-outcome
  values.  Single as well as joint VIMP measures can be requested.  Note
  that calculating VIMP can be computationally expensive (especially
  when the dimension is high), thus if VIMP is not needed, computational
  times can be significantly improved by setting
  \option{importance="none"} which turns VIMP off.
  
  Setting \option{na.action="na.impute"} imputes missing test data
  (x-variables and/or y-outcomes).  Imputation uses the grow-forest and
  only training data is used to impute test data to avoid biasing error
  rates and VIMP (Ishwaran et al. 2008).  See the \command{rfsrc}
  help file for details.
  
  If no test data is provided, then the original training data is used
  and the code reverts to restore mode allowing the user to restore the
  original grow forest.  This is useful, because it gives the user the
  ability to extract outputs from the forest that were not asked for in
  the original grow call.  
  
  If \option{outcome="test"}, the predictor is calculated by using
  y-outcomes from the test data (outcome information must be present).
  In this case, the terminal nodes from the grow-forest are recalculated
  using the y-outcomes from the test set.  This yields a modified
  predictor in which the topology of the forest is based solely on the
  training data, but where the predicted value is based on the test
  data.  Error rates and VIMP are calculated by bootstrapping the test
  data and using out-of-bagging to ensure unbiased estimates.  See the
  examples for illustration.
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
  \item{split.depth}{Matrix [i][j] or array [i][j][k] recording the
    minimal depth for variable [j] for case [i], either averaged over
    the forest, or by tree [k].}
  \item{node.stats}{Split statistics returned when
    \code{statistics=TRUE} which can be parsed using \command{stat.split}.}
  \item{err.rate}{Cumulative OOB error rate for the test data if
    y-outcomes are present.}
  \item{importance}{Test set variable importance (VIMP).  Can be
    \code{NULL}.}
  \item{predicted}{Test set predicted value.}
  \item{predicted.oob}{OOB predicted value (\code{NULL} unless
    \option{outcome="test"}).}\cr
  
  \item{++++++++}{for classification settings, additionally ++++++++} \cr     
  \item{class}{In-bag predicted class labels.}
  \item{class.oob}{OOB predicted class labels (\code{NULL} unless \option{outcome="test"}).}\cr
  
  \item{++++++++}{for multivariate settings, additionally ++++++++} \cr
  
  \item{regrOutput}{List containing performance values for
	  test regression responses (if such values are present).}
	\item{clasOutput}{List containing performance values for
	  test categorical (factor) responses (if such values are present).}
  
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
  \item{cif.oob}{OOB CIF for each event (\code{NULL} unless \option{outcome="test"}).}
  \item{time.interest}{Ordered unique event times.}
  \item{ndead}{Number of events.}
}
\note{
  The dimensions and values of returned objects depend heavily on the
  underlying family and whether y-outcomes are present in the test data.
  In particular, items related to performance will be \code{NULL} when
  y-outcomes are not present.  For multivariate families, predicted
  values, VIMP, error rate, and performance values are stored in the
  lists \code{regrOutput} and \code{clasOutput}.
  
  For detailed definitions of returned values (such as \code{predicted})
  see the \command{rfsrc} help file.
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
  \command{\link{plot.competing.risk}},
  \command{\link{plot.rfsrc}},
  \command{\link{plot.survival}},
  \command{\link{plot.variable}},
  \command{\link{rfsrc}},
  \command{\link{stat.split}},
  \command{\link{vimp}}
}
\examples{
## ------------------------------------------------------------
## typical train/testing scenario
## ------------------------------------------------------------

data(veteran, package = "randomForestSRC")
train <- sample(1:nrow(veteran), round(nrow(veteran) * 0.80))
veteran.grow <- rfsrc(Surv(time, status) ~ ., veteran[train, ], ntree = 100) 
veteran.pred <- predict(veteran.grow, veteran[-train , ])
print(veteran.grow)
print(veteran.pred)

\dontrun{
## ------------------------------------------------------------
## predicted probability and predicted class labels are returned
## in the predict object for classification analyses
## ------------------------------------------------------------

data(breast, package = "randomForestSRC")
breast.obj <- rfsrc(status ~ ., data = breast[(1:100), ], nsplit = 10)
breast.pred <- predict(breast.obj, breast[-(1:100), ])
print(head(breast.pred$predicted))
print(breast.pred$class)

## ------------------------------------------------------------
## example illustrating restore mode
## if predict is called without specifying the test data
## the original training data is used and the forest is restored
## ------------------------------------------------------------

# first we make the grow call
airq.obj <- rfsrc(Ozone ~ ., data = airquality)

# now we restore it and compare it to the original call
# they are identical
predict(airq.obj)
print(airq.obj)

# we can retrieve various outputs that were not asked for in
# in the original call

#here we extract the proximity matrix
prox <- predict(airq.obj, proximity = TRUE)$proximity
print(prox[1:10,1:10])

#here we extract the number of times a variable was used to grow
#the grow forest
var.used <- predict(airq.obj, var.used = "by.tree")$var.used
print(head(var.used))

## ------------------------------------------------------------
## unique feature of randomForestSRC
## cross-validation can be used when factor labels differ over
## training and test data
## ------------------------------------------------------------

# first we convert all x-variables to factors
data(veteran, package = "randomForestSRC")
veteran.factor <- data.frame(lapply(veteran, factor))
veteran.factor$time <- veteran$time
veteran.factor$status <- veteran$status

# split the data into unbalanced train/test data (5/95)
# the train/test data have the same levels, but different labels
train <- sample(1:nrow(veteran), round(nrow(veteran) * .05))
summary(veteran.factor[train,])
summary(veteran.factor[-train,])

# grow the forest on the training data and predict on the test data
veteran.f.grow <- rfsrc(Surv(time, status) ~ ., veteran.factor[train, ]) 
veteran.f.pred <- predict(veteran.f.grow, veteran.factor[-train , ])
print(veteran.f.grow)
print(veteran.f.pred)

## ------------------------------------------------------------
## example illustrating the flexibility of outcome = "test"
## illustrates restoration of forest via outcome = "test"
## ------------------------------------------------------------

# first we make the grow call
data(pbc, package = "randomForestSRC")
pbc.grow <- rfsrc(Surv(days, status) ~ ., pbc, nsplit = 10)

# now use predict with outcome = TEST
pbc.pred <- predict(pbc.grow, pbc, outcome = "test")

# notice that error rates are the same!!
print(pbc.grow)
print(pbc.pred)

# note this is equivalent to restoring the forest
pbc.pred2 <- predict(pbc.grow)
print(pbc.grow)
print(pbc.pred)
print(pbc.pred2)

# similar example, but with na.action = "na.impute"
airq.obj <- rfsrc(Ozone ~ ., data = airquality, na.action = "na.impute")
print(airq.obj)
print(predict(airq.obj))
# ... also equivalent to outcome="test" but na.action = "na.impute" required
print(predict(airq.obj, airquality, outcome = "test", na.action = "na.impute"))

# classification example
iris.obj <- rfsrc(Species ~., data = iris)
print(iris.obj)
print(predict.rfsrc(iris.obj, iris, outcome = "test"))

## ------------------------------------------------------------
## another example illustrating outcome = "test"
## unique way to check reproducibility of the forest
## ------------------------------------------------------------

# primary call
set.seed(542899)
data(pbc, package = "randomForestSRC")
train <- sample(1:nrow(pbc), round(nrow(pbc) * 0.50))
pbc.out <- rfsrc(Surv(days, status) ~ .,  data=pbc[train, ],
        nsplit = 10)

# standard predict call
pbc.train <- predict(pbc.out, pbc[-train, ], outcome = "train")
#non-standard predict call: overlays the test data on the grow forest
pbc.test <- predict(pbc.out, pbc[-train, ], outcome = "test")

# check forest reproducibilility by comparing "test" predicted survival
# curves to "train" predicted survival curves for the first 3 individuals
Time <- pbc.out$time.interest
matplot(Time, t(exp(-pbc.train$chf)[1:3,]), ylab = "Survival", col = 1, type = "l")
matlines(Time, t(exp(-pbc.test$chf)[1:3,]), col = 2)

## ------------------------------------------------------------
## survival analysis using mixed multivariate outcome analysis 
## compare the predicted value to RSF
## ------------------------------------------------------------

# fit the pbc data using RSF
data(pbc, package = "randomForestSRC")
rsf.obj <- rfsrc(Surv(days, status) ~ ., pbc, nsplit = 10)
yvar <- rsf.obj$yvar

# fit a mixed outcome forest using days and status as y-variables
pbc.mod <- pbc
pbc.mod$status <- factor(pbc.mod$status)
mix.obj <- rfsrc(Multivar(days, status) ~., pbc.mod, nsplit = 10)

# compare oob predicted values
rsf.pred <- rsf.obj$predicted.oob
mix.pred <- mix.obj$regrOutput$days$predicted.oob
plot(rsf.pred, mix.pred)

# compare C-index error rate
rsf.err <- randomForestSRC:::cindex(yvar$days, yvar$status, rsf.pred)
mix.err <- 1 - randomForestSRC:::cindex(yvar$days, yvar$status, mix.pred)
cat("RSF                :", rsf.err, "\n")
cat("multivariate forest:", mix.err, "\n")

}
}
\keyword{predict}
\keyword{forest}

