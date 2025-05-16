\name{imbalanced.rfsrc}
\alias{imbalanced.rfsrc}
\alias{imbalanced}
\title{Imbalanced Two Class Problems}

\description{
  Implements various solutions to the two-class imbalanced problem,
  including the newly proposed quantile-classifier approach of
  O'Brien and Ishwaran (2017).  Also includes Breiman's balanced random
  forests undersampling of the majority class.  Performance is assesssed
  using the G-mean, but misclassification error can be requested.
}

\usage{\method{imbalanced}{rfsrc}(formula, data, ntree = 3000, 
  method = c("rfq", "brf", "standard"), splitrule = "auc",
  perf.type = NULL, block.size = NULL, fast = FALSE,
  ratio = NULL, ...)
}

\arguments{

  \item{formula}{A symbolic description of the model to be fit.}

  \item{data}{A data frame containing the two-class y-outcome and x-variables.}

  \item{ntree}{Number of trees to grow.}

  \item{method}{Method used to fit the classifier. The default is \code{"rfq"}, which implements the random forest quantile classifier (RFQ) of O'Brien and Ishwaran (2017). The option \code{"brf"} applies the balanced random forest (BRF) approach of Chen et al. (2004), which undersamples the majority class to match the minority class size. The option \code{"standard"} performs a standard random forest analysis.}

  \item{splitrule}{Splitting rule used to grow trees. The default is \code{"auc"}, which optimizes G-mean performance. Other supported options are \code{"gini"} and \code{"entropy"}.}

  \item{perf.type}{Performance metric used for evaluating the classifier and computing downstream quantities such as VIMP. Defaults depend on the method: \code{"gmean"} for RFQ and BRF; \code{"misclass"} (misclassification error) for standard random forests. Users may override this by specifying \code{"gmean"}, \code{"misclass"}, or \code{"brier"} (normalized Brier score). See examples for usage.}

  \item{block.size}{Controls how the cumulative error rate is computed. If \code{NULL}, it is calculated only once for the final tree. If set to an integer, cumulative error and VIMP are computed in blocks of that size. If unspecified, uses the default in \code{\link{rfsrc}}.}

  \item{fast}{Logical. If \code{TRUE}, uses the fast random forest implementation via \code{\link{rfsrc.fast}} instead of \code{rfsrc}. Improves speed at the cost of accuracy. Applies only to RFQ.}

  \item{ratio}{Optional and experimental. Specifies the proportion (between 0 and 1) of majority class cases to sample during RFQ training. Sampling is without replacement. Ignored for BRF.}

  \item{...}{Additional arguments passed to \code{\link{rfsrc}} to control random forest behavior.}

}

\details{

Imbalanced data, also known as the minority class problem, refers to two-class classification settings where the majority class significantly outnumbers the minority class. This function supports two approaches to address class imbalance:

\itemize{
  \item The random forests quantile classifier (RFQ) proposed by O'Brien and Ishwaran (2017).
  \item The balanced random forest (BRF) undersampling method of Chen et al. (2004).
}

By default, the performance metric used is the G-mean (Kubat et al., 1997), which balances sensitivity and specificity.

\strong{Handling of missing values:} Missing data are not supported for BRF or when the \code{ratio} option is specified. In these cases, records with missing values are removed prior to analysis.

\strong{Variable importance:} Permutation-based VIMP is used by default in this setting, in contrast to anti-VIMP which is the default for other families. Empirical results suggest that permutation VIMP is more reliable in highly imbalanced settings.

\strong{Tree count recommendation:} We recommend using a relatively
large value for \code{ntree} in imbalanced problems to ensure stable
performance estimation, especially for G-mean. As a general guideline,
use at least five times the usual number of trees.

\strong{Performance metrics:} A helper function, \code{\link{get.imbalanced.performance}}, is provided for extracting classification performance summaries. The metric names are self-explanatory in most cases. Some key metrics include:
\itemize{
  \item \code{F1}: The harmonic mean of precision and recall.
  \item \code{F1mod}: The harmonic mean of sensitivity, specificity, precision, and negative predictive value.
  \item \code{F1gmean}: The average of F1 and G-mean.
  \item \code{F1modgmean}: The average of F1mod and G-mean.
}
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

data(breast, package = "randomForestSRC")
breast <- na.omit(breast)
f <- as.formula(status ~ .)

##----------------------------------------------------------------
## default RFQ call
##----------------------------------------------------------------

o.rfq <- imbalanced(f, breast)
print(o.rfq)

## equivalent to:
## rfsrc(f, breast, rfq =  TRUE, ntree = 3000,
##       perf.type = "gmean", splitrule = "auc") 

##----------------------------------------------------------------
## detailed output using customized performance function
##----------------------------------------------------------------

print(get.imbalanced.performance(o.rfq))

##-----------------------------------------------------------------
## RF using misclassification error with gini splitting
## ------------------------------------------------------------

o.std <- imbalanced(f, breast, method = "stand", splitrule = "gini")

##-----------------------------------------------------------------
## RF using G-mean performance with AUC splitting
## ------------------------------------------------------------

o.std <- imbalanced(f, breast, method = "stand", perf.type = "gmean")

## equivalent to:
## rfsrc(f, breast, ntree = 3000, perf.type = "gmean", splitrule = "auc")

##----------------------------------------------------------------
## default BRF call 
##----------------------------------------------------------------

o.brf <- imbalanced(f, breast, method = "brf")

## equivalent to:
## imbalanced(f, breast, method = "brf", perf.type = "gmean")

##----------------------------------------------------------------
## BRF call with misclassification performance 
##----------------------------------------------------------------

o.brf <- imbalanced(f, breast, method = "brf", perf.type = "misclass")

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
##
##  illustrates how to optimize threshold on training data
##  improves Gmean for RFQ in many situations
##
##----------------------------------------------------------------

if (library("caret", logical.return = TRUE)) {

  ## experimental settings
  n <- 2 * 5000 
  q <- 20
  ir <- 6
  f <- as.formula(Class ~ .)
 
  ## simulate the data, create minority class data
  d <- twoClassSim(n, linearVars = 15, noiseVars = q)
  d$Class <- factor(as.numeric(d$Class) - 1)
  idx.0 <- which(d$Class == 0)
  idx.1 <- sample(which(d$Class == 1), sum(d$Class == 1) / ir , replace = FALSE)
  d <- d[c(idx.0,idx.1),, drop = FALSE]

  ## split data into train and test
  trn.pt <- sample(1:nrow(d), size = nrow(d) / 2)
  trn <- d[trn.pt, ]
  tst <- d[setdiff(1:nrow(d), trn.pt), ]

  ## run rfq on training data
  o <- imbalanced(f, trn)

  ## (1) default threshold (2) directly optimized gmean threshold
  th.1 <- get.imbalanced.performance(o)["threshold"]
  th.2 <- get.imbalanced.optimize(o)["threshold"]

  ## training performance
  cat("-------- train performance ---------\n")
  print(get.imbalanced.performance(o, thresh=th.1))
  print(get.imbalanced.performance(o, thresh=th.2))

  ## test performance
  cat("-------- test performance ---------\n")
  pred.o <- predict(o, tst)
  print(get.imbalanced.performance(pred.o, thresh=th.1))
  print(get.imbalanced.performance(pred.o, thresh=th.2))
 
} 

##----------------------------------------------------------------
##  illustrates RFQ with and without SMOTE
## 
## - simulation example using the caret R-package
## - creates imbalanced data by randomly sampling the class 1 data
## - use SMOTE from "imbalance" package to oversample the minority
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

  ## permutation VIMP for BRF with and without blocking
  ## blocked VIMP is a hybrid of Breiman-Cutler/Ishwaran-Kogalur VIMP
  brf <- imbalanced(f, d, method = "brf", importance = "permute", block.size = 1)
  brfB <- imbalanced(f, d, method = "brf", importance = "permute", block.size = 10)

  ## permutation VIMP for RFQ with and without blocking
  rfq <- imbalanced(f, d, importance = "permute", block.size = 1)
  rfqB <- imbalanced(f, d, importance = "permute", block.size = 10)

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
## confidence intervals for G-mean permutation VIMP using subsampling
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

  ## RFQ
  o <- imbalanced(Class ~ ., d, importance = "permute", block.size = 10)

  ## subsample RFQ
  smp.o <- subsample(o, B = 100)
  plot(smp.o, cex.axis = .7)

}


}}
\keyword{imbalanced two-class data}

