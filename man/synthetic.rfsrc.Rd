\name{synthetic}
\alias{synthetic.rfsrc}
\alias{synthetic}
\title{Synthetic Random Forests}

\description{
  Grows a synthetic random forest (RF) using RF machines as synthetic
  features.  Applies only to regression and classification settings.
}

\usage{\method{synthetic}{rfsrc}(formula, data, object, newdata,
  ntree = 1000, mtry = NULL, nodesize = 5, nsplit = 10,
  mtrySeq = NULL, nodesizeSeq = c(1:10,20,30,50,100),
  min.node = 3,
  fast = TRUE,
  use.org.features = TRUE,
  na.action = c("na.omit", "na.impute"),
  oob = TRUE,
  verbose = TRUE,
  ...)
}

\arguments{

  \item{formula}{Model to be fit. Must be specified unless \code{object} is provided.}

  \item{data}{Data frame containing the y-outcome and x-variables. Must be specified unless \code{object} is provided.}

  \item{object}{An object of class \code{(rfsrc, synthetic)}. Used to bypass the fitting step. Not required if \code{formula} and \code{data} are supplied.}

  \item{newdata}{Optional test data for prediction. If omitted, the training data is used.}

  \item{ntree}{Number of trees used in each RF machine.}

  \item{mtry}{\code{mtry} value used in the final synthetic forest.}

  \item{nodesize}{\code{nodesize} value used in the final synthetic forest.}

  \item{nsplit}{Number of random splits used in randomized splitting. Increases speed when set to a small positive integer.}

  \item{mtrySeq}{Sequence of \code{mtry} values used to train the ensemble of RF machines. If \code{NULL}, defaults to \code{ceiling(p / 3)}, where \code{p} is the number of variables.}

  \item{nodesizeSeq}{Sequence of \code{nodesize} values used to train the ensemble of RF machines.}

  \item{min.node}{Minimum forest-averaged number of terminal nodes required for an RF machine to be retained as a synthetic feature.}

  \item{fast}{Use \code{rfsrc.fast} instead of \code{rfsrc} to fit base learners? Improves speed at the cost of accuracy.}

  \item{use.org.features}{Should original features be included alongside synthetic features in the final synthetic forest?}

  \item{na.action}{Action to be taken on missing data. The default, \code{"na.omit"}, removes records with any missing values. Set to \code{"na.impute"} to pre-impute data using \command{impute.rfsrc}.}

  \item{oob}{Preserve out-of-bag (OOB) estimation for error rates and VIMP? Defaults to \code{TRUE}.}

  \item{verbose}{Display detailed output of the fitting process? Defaults to \code{FALSE}.}

  \item{...}{Additional arguments passed to \code{rfsrc} for training the synthetic forest.}
}

\details{

  A collection of random forests are trained using different values of
  \code{nodesize} (and optionally \code{mtry}). The out-of-bag (OOB)
  predicted values from these forests are then used as synthetic
  features (referred to as RF machines) to train a final synthetic random
  forest. The original features can optionally be included in the final
  model.
  

  This approach is currently implemented for regression and
  classification settings (both univariate and multivariate).
  

  Synthetic features are generated using OOB predictions to prevent
  overfitting. To ensure that performance metrics (such as error rates
  and VIMP) remain valid, the same bootstrap samples are reused across
  all trees for both the synthetic forest and its constituent RF
  machines. This behavior is controlled by the \code{oob=TRUE}
  option. Disabling this may yield misleading performance estimates and
  should be done with caution.
  

  If values for \code{mtrySeq} are provided, RF machines are constructed
  for every combination of \code{nodesizeSeq} and \code{mtrySeq}.
  

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
  \command{\link{rfsrc.fast}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## compare synthetic forests to regular forest (classification)
## ------------------------------------------------------------

## rfsrc and synthetic calls
if (library("mlbench", logical.return = TRUE)) {

  ## simulate the data 
  ring <- data.frame(mlbench.ringnorm(250, 20))

  ## classification forests
  ringRF <- rfsrc(classes ~., ring)

  ## synthetic forests
  ## 1 = nodesize varied
  ## 2 = nodesize/mtry varied
  ringSyn1 <- synthetic(classes ~., ring)
  ringSyn2 <- synthetic(classes ~., ring, mtrySeq = c(1, 10, 20))

  ## test-set performance
  ring.test <- data.frame(mlbench.ringnorm(500, 20))
  pred.ringRF <- predict(ringRF, newdata = ring.test)
  pred.ringSyn1 <- synthetic(object = ringSyn1, newdata = ring.test)$rfSynPred
  pred.ringSyn2 <- synthetic(object = ringSyn2, newdata = ring.test)$rfSynPred


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
regF <- rfsrc(y ~ ., dat[-test, ], )
pred.regF <- predict(regF, dat[test, ])

## synthetic forests using fast rfsrc
synF1 <- synthetic(y ~ ., dat[-test, ], newdata = dat[test, ])
synF2 <- synthetic(y ~ ., dat[-test, ],
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
mvSyn <- synthetic(cbind(carb, mpg, cyl) ~., mtcars.new[trn,])
mvSyn.pred <- synthetic(object = mvSyn, newdata = mtcars.new[-trn,])
}
}
\keyword{forest}
\keyword{predict}
