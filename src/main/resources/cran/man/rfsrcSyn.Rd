\name{rfsrcSyn}
\alias{rfsrcSyn}
\alias{rfsrcSyn.rfsrc}
\title{Synthetic Random Forests}
\description{
  Grows a synthetic random forest (RF) using RF machines as synthetic
  features.  Applies only to regression and classification settings.
}
\usage{\method{rfsrcSyn}{rfsrc}(formula, data, object, newdata,
  ntree = 1000,
  mtry = NULL,
  mtrySeq = NULL,
  nodesize = 5,
  nodesizeSeq = c(1:10,20,30,50,100),
  nsplit = 0,
  min.node = 3,
  use.org.features = TRUE,
  
  na.action = c("na.omit", "na.impute"),
  
   
  oob = TRUE,
  verbose = TRUE,
  ...)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit.
    Must be specified unless \code{object} is given.}
  \item{data}{Data frame containing the y-outcome and x-variables in
    the model. Must be specified unless \code{object} is given.}
  \item{object}{An object of class \code{(rfsrc, synthetic)}.
    Not required when \code{formula} and \code{data} are supplied.}
  \item{newdata}{Test data used for prediction (optional).}
  \item{ntree}{Number of trees.}
  \item{mtry}{mtry value for synthetic forest.}
  \item{mtrySeq}{Sequence of mtry values used for fitting the
    collection of RF machines.  If \code{NULL}, set to the default
    value \code{p}/3.}
  \item{nodesize}{Nodesize value for the synthetic forest.}
  \item{nodesizeSeq}{Sequence of nodesize values used for the fitting the
    collection of RF machines.}
  \item{nsplit}{If non-zero, nsplit-randomized splitting is used which can
    significantly increase speed.}
  \item{min.node}{Minimum forest averaged number of nodes a RF machine
    must exceed in order to be used as a synthetic feature.}
  \item{use.org.features}{In addition to synthetic features, should
    the original features be used when fitting synthetic forests?}
  \item{na.action}{Missing value action. The default \code{na.omit}
    removes the entire record if even one of its entries is \code{NA}.
    The action \code{na.impute} pre-imputes the data using fast
    imputation via \command{impute.rfsrc}.}
  \item{oob}{Preserve "out-of-bagness" so that error rates and
    VIMP are honest?  Default is yes (\option{oob=TRUE}).}
  \item{verbose}{Set to \code{TRUE} for verbose output.}
  \item{...}{Further arguments to be passed to the \code{rfsrc}
    function used for fitting the synthetic forest.}
}
\details{
  A collection of random forests are fit using different nodesize
  values.  The predicted values from these machines are then used as
  synthetic features (called RF machines) to fit a synthetic random
  forest (the original features are also used in constructing the
  synthetic forest).  Currently only implemented for regression and
  classification settings (univariate and multivariate).

  Synthetic features are calculated using out-of-bag (OOB) data to avoid
  over-using training data.  However, to guarantee that performance
  values such as error rates and VIMP are honest, bootstrap draws are
  fixed across all trees used in the construction of the synthetic
  forest and its synthetic features.  The option \option{oob=TRUE}
  ensures that this happens.  Change this option at your own peril.
  
  If values for \code{mtrySeq} are given, RF machines are constructed
  for each combination of nodesize and mtry values specified by
  \code{nodesizeSeq} \code{mtrySeq}.
}
\value{
  A list with the following components:
  \item{rfMachines}{RF machines used to construct the synthetic
    features.}
  \item{rfSyn}{The (grow) synthetic RF built over training data.}
  \item{rfSynPred}{The predict synthetic RF built over test data (if available).}
  \item{synthetic}{List containing the synthetic features.}
  \item{opt.machine}{Optimal machine: RF machine with smallest
    OOB error rate.}
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran H. and Malley J.D. (2014).  Synthetic learning machines.
  \emph{BioData Mining}, 7:28.
}
\seealso{
  \command{\link{rfsrc}},
  \command{\link{impute.rfsrc}}
}
\examples{
\dontrun{
## ------------------------------------------------------------
## compare synthetic forests to regular forest (classification)
## ------------------------------------------------------------

## rfsrc and rfsrcSyn calls
if (library("mlbench", logical.return = TRUE)) {

  ## simulate the data 
  ring <- data.frame(mlbench.ringnorm(250, 20))

  ## classification forests
  ringRF <- rfsrc(classes ~., data = ring)

  ## synthetic forests:
  ## 1 = nodesize varied
  ## 2 = nodesize/mtry varied
  ringSyn1 <- rfsrcSyn(classes ~., data = ring)
  ringSyn2 <- rfsrcSyn(classes ~., data = ring, mtrySeq = c(1, 10, 20))

  ## test-set performance
  ring.test <- data.frame(mlbench.ringnorm(500, 20))
  pred.ringRF <- predict(ringRF, newdata = ring.test)
  pred.ringSyn1 <- rfsrcSyn(object = ringSyn1, newdata = ring.test)$rfSynPred
  pred.ringSyn2 <- rfsrcSyn(object = ringSyn2, newdata = ring.test)$rfSynPred


  print(pred.ringRF)
  print(pred.ringSyn1)
  print(pred.ringSyn2)

}

## ------------------------------------------------------------
## compare synthetic forest to regular forest (regression)
## ------------------------------------------------------------

## simulate the data
n <- 250
ntest <- 1000
N <- n + ntest
d <- 50
std <- 0.1
x <- matrix(runif(N * d, -1, 1), ncol = d)
y <- 1 * (x[,1] + x[,4]^3 + x[,9] + sin(x[,12]*x[,18]) + rnorm(n, sd = std)>.38)
dat <- data.frame(x = x, y = y)
test <- (n+1):N

## regression forests
regF <- rfsrc(y ~ ., data = dat[-test, ], )
pred.regF <- predict(regF, dat[test, ], importance = "none")

## synthetic forests
## we pass both the training and testing data
## but this can be split into separate commands as in the
## previous classification example
synF1 <- rfsrcSyn(y ~ ., data = dat[-test, ],
  newdata = dat[test, ])
synF2 <- rfsrcSyn(y ~ ., data = dat[-test, ],
  newdata = dat[test, ], mtrySeq = c(1, 10, 20, 30, 40, 50))

## standardized MSE performance
mse <- c(tail(pred.regF$err.rate, 1),
         tail(synF1$rfSynPred$err.rate, 1),
         tail(synF2$rfSynPred$err.rate, 1)) / var(y[-test])
names(mse) <- c("forest", "synthetic1", "synthetic2")
print(mse)

## ------------------------------------------------------------
## multivariate synthetic forests
## ------------------------------------------------------------

mtcars.new <- mtcars
mtcars.new$cyl <- factor(mtcars.new$cyl)
mtcars.new$carb <- factor(mtcars.new$carb, ordered = TRUE)
trn <- sample(1:nrow(mtcars.new), nrow(mtcars.new)/2)
mvSyn <- rfsrcSyn(cbind(carb, mpg, cyl) ~., data = mtcars.new[trn,])
mvSyn.pred <- rfsrcSyn(object = mvSyn, newdata = mtcars.new[-trn,])
}
}
\keyword{forest}
\keyword{predict}