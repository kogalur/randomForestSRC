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
  \item{ntree}{Function specifying requested number of trees used for
    growing the forest.  Inputs are dimension and number of holdout
    variables. The requested number of trees can also be a number.}
  \item{nsplit}{Non-negative integer value specifying number of 
    random split points used to split a node (deterministic splitting
    corresponds to the value zero and is much slower).}
  \item{ntime}{Integer value used for survival to
    constrain ensemble calculations to a grid of \code{ntime} time points.}
  \item{sampsize}{Function specifying size of subsampled data. Can also
    be a number.}
  \item{samptype}{Type of bootstrap used.}
  \item{vtry}{Number of variables randomly selected to be held out when
    growing a tree.  This can also be set to a list for a targeted hold
    out VIMP analysis.  See details below for more information.}
  \item{block.size}{Specifies number of trees in a block when calculating
    holdout variable importance.}
  \item{...}{Further arguments to be passed to \code{\link{rfsrc}}.}
}
\details{

  Holdout variable importance (holdout VIMP) is based on comparing error
  performance of two mini forests of trees (blocks of trees): the first in
  which a random set of \code{vtry} features are held out (the holdout
  forest), and the second in which no features are held out (the
  baseline forest).  

  To summarize, holdout VIMP measures the importance of a variable when
  that variable is truly removed from the tree growing process.

  Specifically, if a feature is held out in a block of trees, we refer
  to this as the (feature, block) pair.  The bootstrap for the trees in
  a (feature, block) pair are identical in both forests. That is, the
  holdout block is grown by holding out the feature, and the baseline
  block is grown over the same trees, with the same bootstrap, but
  without holding out any features.  \code{vtry} controls how many
  features are held out in every tree. If set to one (default), only one
  variable is held out in every tree.  Once a (feature, block) of trees
  has been grown, holdout VIMP for a given variable v is calculated as
  follows.  Gather the block of trees where the feature was held out
  (from the holdout forest) and calculate OOB prediction error.  Next
  gather the corresponding block of trees where v was not held out (from
  the baseline forest) and calculate OOB prediction error.  Holdout VIMP
  for the (feature, block) pair is the difference between these two
  values.  The final holdout VIMP estimate for a feature v is obtained
  by averaging holdout VIMP for (feature=v, block) over all blocks.

  Accuracy of hold out VIMP depends critically on total number of trees.
  If total number of trees is too small, then number of times a variable
  is held out will be small and OOB error can suffer from high variance.
  Therefore, \code{ntree} should be set fairly high---we recommend using
  1000 times the number of features.  Increasing \code{vtry} is another
  way to increase number of times a variable is held out and therefore
  reduces the burden of growing a large number of trees.  In particular,
  total number of trees needed decreases linearly with \code{vtry}.  The
  default \code{ntree} equals 1000 trees for each feature divided by
  \code{vtry}.  Keep in mind intrepretation of holdout VIMP is altered
  when \code{vtry} is different than one.  Thus this option should be
  used with caution.

  Accuracy also depends on the value of \code{block.size}.  Smaller
  values generally produce better results but are more computationally
  demanding.  The most computationally demanding, but most accurate, is
  \code{block.size=1}.  This is similar to how \code{block.size} is used
  for usual variable importance: see the help file for \command{rfsrc}
  for details.  Note the value of \code{block.size} should not exceed
  \code{ntree} divided by number of features, otherwise there may not be
  enough trees to satisify the target block size for a feature and
  missing values will result.

  A targeted hold out VIMP analysis can be requested by setting
  \code{vtry} to a list with two entries.  The first entry is a vector
  of integer values specifying the variables of interest.  The second
  entry is a boolean logical flag indicating whether individual or joint
  VIMP should be calculated.  For example, suppose variables 1, 4 and 5
  are our variables of interest.  To calculate holdout VIMP for these
  variables, and these variables only, \code{vtry} would be specified by

   vtry = list(xvar = c(1, 4, 5), joint = FALSE)

   On the other hand, if we are interested in the joint effect when we
   remove the three variables simultaneously, then

   vtry = list(xvar = c(1, 4, 5), joint = TRUE)

   The benefits of a targeted analysis is that the user may have a
   pre-conceived idea of which variables are interesting.  Only VIMP for
   these variables will be calculated which greatly reduces
   computational time.  Another benefit is that when joint VIMP is
   requested, this provides the user with a way to assess importance of
   specific groups of variables.  See the iris example below for
   illustration.
     

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
