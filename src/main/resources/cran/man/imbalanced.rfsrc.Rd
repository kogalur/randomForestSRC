\name{imbalanced.rfsrc}
\alias{imbalanced.rfsrc}
\alias{imbalanced}
\alias{get.imbalanced.performance}
\title{Imbalanced Two Class Problems}
\description{
  Implements various solutions to the two-class imbalanced problem,
  including the newly proposed quantile-classifier approach of
  O'Brien and Ishwaran (2017).  Also includes Breiman's balanced random
  forests undersampling of the majority class.  Performance is assesssed
  using the G-mean, but misclassification error can be requested.
}
\usage{\method{imbalanced}{rfsrc}(formula, data, ntree = 3000, 
  method = c("rfq", "brf", "standard"),
  block.size = NULL, perf.type = NULL, fast = FALSE,
  ratio = NULL, optimize = FALSE, ngrid = 1e4,
  ...)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit.}
  
  \item{data}{Data frame containing the two-class y-outcome and
    x-variables.}

  \item{ntree}{Number of trees.}

  \item{method}{Method used for fitting the classifier.  The default is
    \code{rfq} which is the random forests quantile-classifer (RFQ)
    approach of O'Brien and Ishwaran (2017).  The method \code{brf}
    implements the balanced random forest (BRF) method of Chen et
    al. (2004) which undersamples the majority class so that its
    cardinality matches that of the minority class.  The method
    \code{standard} implements a standard random forest analysis.}

  \item{perf.type}{Measure used for assessing performance (and all
    downstream calculations based on it such as variable importance).
    The default for \code{rfq} and \code{brf} is to use the G-mean
    (Kubat et al., 1997).  For standard random forests, the default is
    misclassification error.  Users can over-ride the default performance
    measure by manually selecting either \code{g.mean} for the G-mean,
    \code{misclass} for misclassification error, or \code{brier} for the
    normalized Brier score. See the examples below.}

    \item{block.size}{Should the cumulative error rate be calculated on
    every tree?  When \code{NULL}, it will only be calculated on the
    last tree. If importance is requested, VIMP is calculated in
    "blocks" of size equal to \code{block.size}.}

  \item{fast}{Use fast random forests, \code{rfsrc.fast}, in place of
    \code{rfsrc}?  Improves speed but is less accurate.  Only applies to
    RFQ.}

  \item{ratio}{This is an optional parameter for expert users 
    and included only for experimental purposes.  Used to specify the
    ratio (between 0 and 1) for undersampling the majority class.
    Sampling is without replacement.  Option is ignored for BRF.} 

  \item{optimize}{Calculate the G-mean under various threshold values?
    Returns out-of-bag G-mean values for each tested threshold value.
    See examples below for illustration.}

  \item{ngrid}{Number of threshold values attempted when \code{optimize}
    is requested}

  \item{...}{Further arguments to be passed to the \command{rfsrc}
    function to specify random forest parameters.}

}
\details{
  Imbalanced data, or the so-called imbalanced minority class problem,
  refers to classification settings involving two-classes where the
  ratio of the majority class to the minority class is much larger than
  one. Two solutions to the two-class imbalanced problem are provided
  here, including the newly proposed random forests quantile-classifier
  (RFQ) of O'Brien and Ishwaran (2017), and the balanced random forests
  (BRF) undersampling approach of Chen et al. (2004).  The default
  performance metric is the G-mean (Kubat et al., 1997).

  Currently, missing values cannot be handled for BRF or when the
  \code{ratio} option is used; in these cases, missing data is removed
  prior to the analysis.

  We recommend setting \code{ntree} to a relatively large value when
  dealing with imbalanced data to ensure convergence of the performance
  value -- this is especially true for the G-mean.  Consider using 5 times
  the usual number of trees.
}   
\value{
  A two-class random forest fit under the requested method and
  performance value.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{

  Chen, C., Liaw, A. and Breiman, L. (2004). Using random forest to learn
  imbalanced data. University of California, Berkeley, Technical Report
  110.
  
  Kubat, M., Holte, R. and Matwin, S. (1997). Learning when negative
  examples abound. \emph{Machine Learning}, ECML-97: 146-153.

  O'Brien R. and Ishwaran H. (2019).  A random forests quantile
  classifier for class imbalanced data. \emph{Pattern Recognition},
  90, 232-249
  
}
\seealso{
  \command{\link{rfsrc}},
  \command{\link{rfsrc.fast}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## use the breast data for illustration
## ------------------------------------------------------------

data(breast, package = "_PROJECT_PACKAGE_NAME_")
breast <- na.omit(breast)
f <- as.formula(status ~ .)

##----------------------------------------------------------------
## default RFQ call
##----------------------------------------------------------------

o.rfq <- imbalanced(f, breast)
print(o.rfq)

## equivalent to:
## rfsrc(f, breast, rfq =  TRUE, ntree = 3000, perf.type = "g.mean") 

##----------------------------------------------------------------
## detailed output using customized performance function
##----------------------------------------------------------------

print(get.imbalanced.performance(o.rfq))

##----------------------------------------------------------------
## RFQ with AUC splitting
##----------------------------------------------------------------

print(o.rfq <- imbalanced(f, breast, splitrule = "auc"))
print(get.imbalanced.performance(o.rfq))

##----------------------------------------------------------------
## RFQ call with fast rfsrc
##----------------------------------------------------------------

o.rfq <- imbalanced(f, breast, fast = TRUE)
print(o.rfq)

## equivalent to:
## rfsrc.fast(f, breast, ntree = 3000, rfq =  TRUE, perf.type = "g.mean") 

##-----------------------------------------------------------------
## standard RF (uses misclassification)
## ------------------------------------------------------------

o.std <- imbalanced(f, breast, method = "stand")

##-----------------------------------------------------------------
## standard RF using G-mean performance
## ------------------------------------------------------------

o.std <- imbalanced(f, breast, method = "stand", perf.type = "g.mean")

## equivalent to:
## rfsrc(f, breast, ntree = 3000, perf.type = "g.mean")

##----------------------------------------------------------------
## default BRF call 
##----------------------------------------------------------------

o.brf <- imbalanced(f, breast, method = "brf")

## equivalent to:
## imbalanced(f, breast, method = "brf", perf.type = "g.mean")

##----------------------------------------------------------------
## BRF call with misclassification performance 
##----------------------------------------------------------------

o.brf <- imbalanced(f, breast, method = "brf", perf.type = "misclass")

##----------------------------------------------------------------
## RFQ with optimized threshold
##----------------------------------------------------------------

o.rfq.opt <- imbalanced(f, breast, optimize = TRUE)
plot(o.rfq.opt$gmean, type = "l")

##----------------------------------------------------------------
## train/test example
##----------------------------------------------------------------

trn <- sample(1:nrow(breast), size = nrow(breast) / 2)
o.trn <- imbalanced(f, breast[trn,], importance = TRUE)
o.tst <- predict(o.trn, breast[-trn,], importance = TRUE)
print(o.trn)
print(o.tst)
print(100 * cbind(o.trn$impo[, 1], o.tst$impo[, 1]))

##----------------------------------------------------------------
## simulation example using the caret R-package
## creates imbalanced data by randomly sampling the class 1 data
##
## uses SMOTE from "imbalance" package to oversample the minority
## illustrates RFQ with and without SMOTE
##
##----------------------------------------------------------------

if (library("caret", logical.return = TRUE) &
    library("imbalance", logical.return = TRUE)) {

  ## experimental settings
  n <- 5000
  q <- 20
  ir <- 6
  f <- as.formula(Class ~ .)
 
  ## simulate the data, create minority class data
  d <- twoClassSim(n, linearVars = 15, noiseVars = q)
  d$Class <- factor(as.numeric(d$Class) - 1)
  idx.0 <- which(d$Class == 0)
  idx.1 <- sample(which(d$Class == 1), sum(d$Class == 1) / ir , replace = FALSE)
  d <- d[c(idx.0,idx.1),, drop = FALSE]
  d <- d[sample(1:nrow(d)), ]

  ## define train/test split
  trn <- sample(1:nrow(d), size = nrow(d) / 2, replace = FALSE)

  ## now make SMOTE training data
  newd.50 <- mwmote(d[trn, ], numInstances = 50, classAttr = "Class")
  newd.500 <- mwmote(d[trn, ], numInstances = 500, classAttr = "Class")

  ## fit RFQ with and without SMOTE
  o.with.50 <- imbalanced(f, rbind(d[trn, ], newd.50)) 
  o.with.500 <- imbalanced(f, rbind(d[trn, ], newd.500))
  o.without <- imbalanced(f, d[trn, ])
  
  ## compare performance on test data
  print(predict(o.with.50, d[-trn, ]))
  print(predict(o.with.500, d[-trn, ]))
  print(predict(o.without, d[-trn, ]))
  
}

##----------------------------------------------------------------
## simulation example using the caret R-package similar to above
##
## illustrates effectiveness of blocked VIMP
##
##----------------------------------------------------------------

if (library("caret", logical.return = TRUE)) {

  ## experimental settings
  n <- 1000
  q <- 20
  ir <- 6
  f <- as.formula(Class ~ .)
 
  ## simulate the data, create minority class data
  d <- twoClassSim(n, linearVars = 15, noiseVars = q)
  d$Class <- factor(as.numeric(d$Class) - 1)
  idx.0 <- which(d$Class == 0)
  idx.1 <- sample(which(d$Class == 1), sum(d$Class == 1) / ir , replace = FALSE)
  d <- d[c(idx.0,idx.1),, drop = FALSE]

  ## VIMP for BRF with and without blocking
  ## blocked VIMP is a hybrid of Breiman-Cutler/Ishwaran-Kogalur VIMP
  brf <- imbalanced(f, d, method = "brf", importance = TRUE, block.size = 1)
  brfB <- imbalanced(f, d, method = "brf", importance = TRUE, block.size = 10)

  ## VIMP for RFQ with and without blocking
  rfq <- imbalanced(f, d, importance = TRUE, block.size = 1)
  rfqB <- imbalanced(f, d, importance = TRUE, block.size = 10)

  ## compare VIMP values
  imp <- 100 * cbind(brf$importance[, 1], brfB$importance[, 1],
                     rfq$importance[, 1], rfqB$importance[, 1])
  legn <- c("BRF", "BRF-block", "RFQ", "RFQ-block")
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

##----------------------------------------------------------------
##
## confidence intervals for G-mean VIMP using subsampling
##
##----------------------------------------------------------------

if (library("caret", logical.return = TRUE)) {

  ## experimental settings
  n <- 1000
  q <- 20
  ir <- 6
  f <- as.formula(Class ~ .)
 
  ## simulate the data, create minority class data
  d <- twoClassSim(n, linearVars = 15, noiseVars = q)
  d$Class <- factor(as.numeric(d$Class) - 1)
  idx.0 <- which(d$Class == 0)
  idx.1 <- sample(which(d$Class == 1), sum(d$Class == 1) / ir , replace = FALSE)
  d <- d[c(idx.0,idx.1),, drop = FALSE]

  ## q-classifier
  oq <- imbalanced(Class ~ ., d, splitrule = "auc",
              importance = TRUE, block.size = 10)

  ## subsample the q-classifier
  smp.oq <- subsample(oq, B = 100)
  plot(smp.oq, cex.axis = .7)

}


}}
\keyword{imbalanced two-class data}

