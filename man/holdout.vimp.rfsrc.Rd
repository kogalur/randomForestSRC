\name{holdout.vimp.rfsrc}
\alias{holdout.vimp.rfsrc}
\alias{holdout.vimp}
\title{Hold out variable importance (VIMP)}

\description{
  Hold out VIMP is calculated from the error rate of mini ensembles of
  trees (blocks of trees) grown with and without a variable.  Applies to
  all families.
}

\usage{\method{holdout.vimp}{rfsrc}(formula, data,
  ntree = function(p, vtry){1000 * p / vtry},
  nsplit = 10,
  ntime = 50,
  sampsize = function(x){x * .632},
  samptype = "swor",
  block.size = 10,
  vtry = 1,
  ...)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit.}

  \item{data}{Data frame containing the y-outcome and x-variables.}

  \item{ntree}{Specifies the number of trees used to grow the forest. Can be a function of 
  data dimension and number of holdout variables, or a fixed numeric value.}

  \item{nsplit}{Non-negative integer specifying the number of random split points used to split a node. 
  A value of zero corresponds to deterministic splitting, which is significantly slower.}

  \item{ntime}{Integer value used for survival settings to constrain ensemble 
  calculations to a grid of \code{ntime} time points.}

  \item{sampsize}{Specifies the size of the subsampled data. Can be either a function or a numeric value.}

  \item{samptype}{Type of bootstrap used when subsampling.}

  \item{vtry}{Number of variables randomly selected to be held out when
  growing a tree.  Can also be a list for targeted holdout variable
  importance analysis. See \code{details} for more information.}

  \item{block.size}{Specifies the number of trees in a block when calculating holdout variable importance.}

  \item{...}{Further arguments passed to \code{\link{rfsrc}}.}
}

\details{

  Holdout variable importance (holdout VIMP) measures the importance of a variable 
  by comparing prediction error between two forests (blocks of trees): 
  one in which selected variables are held out during tree growing (the \emph{holdout forest}) 
  and one in which no variables are held out (the \emph{baseline forest}).

  For each variable-block combination, the bootstrap samples used to grow the trees 
  are the same in both forests. The difference in out-of-bag (OOB) prediction error 
  between the holdout and baseline forests gives the holdout VIMP for that variable-block pair. 
  The final holdout VIMP for a variable is the average of these differences over all blocks 
  in which the variable was held out.

  The option \code{vtry} controls how many variables are held out per tree. 
  The default is one, meaning a single variable is held out per tree. 
  Larger values of \code{vtry} increase the number of times each variable is held out, 
  reducing the required total number of trees. However, interpretation of holdout VIMP 
  changes when \code{vtry} exceeds one, and this option should be used cautiously.

  High accuracy requires a sufficiently large number of trees. 
  As a general guideline, we recommend using \code{ntree = 1000 * p / vtry}, 
  where \code{p} is the number of features. Accuracy also depends on \code{block.size}, 
  which determines how many trees comprise a block. Smaller values yield better accuracy 
  but are computationally more demanding. The most accurate setting is \code{block.size = 1}. 
  Ensure that \code{block.size} does not exceed \code{ntree / p}, 
  otherwise insufficient trees may be available for certain variables.

  Targeted holdout VIMP analysis can be requested by specifying \code{vtry} as a list 
  with two components: a vector of variable indices (\code{xvar}) and a logical flag 
  \code{joint} indicating whether to compute joint VIMP. For example, to compute holdout VIMP 
  only for variables 1, 4, and 5 individually:

  \code{vtry = list(xvar = c(1, 4, 5), joint = FALSE)}

  To compute the joint effect of removing these three variables together:

  \code{vtry = list(xvar = c(1, 4, 5), joint = TRUE)}

  Targeted analysis is useful when the user has prior knowledge of variables of interest 
  and can significantly reduce computation. Joint VIMP quantifies the combined importance 
  of specific groups of variables. See the Iris example below for illustration.

}

\value{
  Invisibly a list with the following components (which themselves can be lists):
  \item{importance}{Holdout VIMP.}
  \item{baseline}{Prediction error for the baseline forest.}
  \item{holdout}{Prediction error for the holdout forest.}
}

\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{

  Lu M. and Ishwaran H. (2018). Expert Opinion: A prediction-based
  alternative to p-values in regression models. \emph{J. Thoracic and
    Cardiovascular Surgery}, 155(3), 1130--1136.

  
}

\seealso{
  \command{\link{vimp.rfsrc}}
}
\examples{
\donttest{


## ------------------------------------------------------------
## regression analysis
## ------------------------------------------------------------

## new York air quality measurements
airq.obj <- holdout.vimp(Ozone ~ ., data = airquality, na.action = "na.impute")
print(airq.obj$importance)

## ------------------------------------------------------------
## classification analysis
## ------------------------------------------------------------

## iris data
iris.obj <- holdout.vimp(Species ~., data = iris)
print(iris.obj$importance)

## iris data using brier prediction error
iris.obj <- holdout.vimp(Species ~., data = iris, perf.type = "brier")
print(iris.obj$importance)

## ------------------------------------------------------------
## illustration of targeted holdout vimp analysis
## ------------------------------------------------------------

## iris data - only interested in variables 3 and 4
vtry <- list(xvar = c(3, 4), joint = FALSE)
print(holdout.vimp(Species ~., data = iris, vtry = vtry)$impor)

## iris data - joint importance of variables 3 and 4
vtry <- list(xvar = c(3, 4), joint = TRUE)
print(holdout.vimp(Species ~., data = iris, vtry = vtry)$impor)

## iris data - joint importance of variables 1 and 2
vtry <- list(xvar = c(1, 2), joint = TRUE)
print(holdout.vimp(Species ~., data = iris, vtry = vtry)$impor)


## ------------------------------------------------------------
## imbalanced classification (using RFQ)
## ------------------------------------------------------------

if (library("caret", logical.return = TRUE)) {

  ## experimental settings
  n <- 400
  q <- 20
  ir <- 6
  f <- as.formula(Class ~ .)
 
  ## simulate the data, create minority class data
  d <- twoClassSim(n, linearVars = 15, noiseVars = q)
  d$Class <- factor(as.numeric(d$Class) - 1)
  idx.0 <- which(d$Class == 0)
  idx.1 <- sample(which(d$Class == 1), sum(d$Class == 1) / ir , replace = FALSE)
  d <- d[c(idx.0,idx.1),, drop = FALSE]

  ## VIMP for RFQ with and without blocking
  vmp1 <- imbalanced(f, d, importance = TRUE, block.size = 1)$importance[, 1]
  vmp10 <- imbalanced(f, d, importance = TRUE, block.size = 10)$importance[, 1]

  ## holdout VIMP for RFQ with and without blocking
  hvmp1 <- holdout.vimp(f, d, rfq =  TRUE,
               perf.type = "g.mean", block.size = 1)$importance[, 1]
  hvmp10 <- holdout.vimp(f, d, rfq =  TRUE,
               perf.type = "g.mean", block.size = 10)$importance[, 1]
  
  ## compare VIMP values
  imp <- 100 * cbind(vmp1, vmp10, hvmp1, hvmp10)
  legn <- c("vimp-1", "vimp-10","hvimp-1", "hvimp-10")
  colr <- rep(4,20+q)
  colr[1:20] <- 2
  ylim <- range(c(imp))
  nms <- 1:(20+q)
  par(mfrow=c(2,2))
  barplot(imp[,1],col=colr,las=2,main=legn[1],ylim=ylim,names.arg=nms)
  barplot(imp[,2],col=colr,las=2,main=legn[2],ylim=ylim,names.arg=nms)
  barplot(imp[,3],col=colr,las=2,main=legn[3],ylim=ylim,names.arg=nms)
  barplot(imp[,4],col=colr,las=2,main=legn[4],ylim=ylim,names.arg=nms)

}

## ------------------------------------------------------------
## multivariate regression analysis
## ------------------------------------------------------------
mtcars.mreg <- holdout.vimp(Multivar(mpg, cyl) ~., data = mtcars,
                                    vtry = 3,
                                    block.size = 1,
                                    samptype = "swr",
                                    sampsize = dim(mtcars)[1])
print(mtcars.mreg$importance)

## ------------------------------------------------------------
## mixed outcomes analysis
## ------------------------------------------------------------

mtcars.new <- mtcars
mtcars.new$cyl <- factor(mtcars.new$cyl)
mtcars.new$carb <- factor(mtcars.new$carb, ordered = TRUE)
mtcars.mix <- holdout.vimp(cbind(carb, mpg, cyl) ~., data = mtcars.new,
                                   ntree = 100,
                                   block.size = 2,
                                   vtry = 1)
print(mtcars.mix$importance)

##------------------------------------------------------------
## survival analysis
##------------------------------------------------------------

## Primary biliary cirrhosis (PBC) of the liver
data(pbc, package = "randomForestSRC")
pbc.obj <- holdout.vimp(Surv(days, status) ~ ., pbc,
                                nsplit = 10,
                                ntree = 1000,
                                na.action = "na.impute")
print(pbc.obj$importance)

##------------------------------------------------------------
## competing risks
##------------------------------------------------------------

## WIHS analysis
## cumulative incidence function (CIF) for HAART and AIDS stratified by IDU

data(wihs, package = "randomForestSRC")
wihs.obj <- holdout.vimp(Surv(time, status) ~ ., wihs,
                                 nsplit = 3,
                                 ntree = 100)
print(wihs.obj$importance)

}}

\keyword{vimp}
