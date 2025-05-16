\name{subsample.rfsrc}
\alias{subsample}
\alias{subsample.rfsrc}
\title{Subsample Forests for VIMP Confidence Intervals}

\description{
  Use subsampling to calculate confidence intervals and standard
  errors for VIMP (variable importance). Applies to all families.
}

\usage{\method{subsample}{rfsrc}(obj,
  B = 100,
  block.size = 1,
  importance,
  subratio = NULL,
  stratify = TRUE,
  performance = FALSE,
  performance.only = FALSE,
  joint = FALSE,
  xvar.names = NULL,
  bootstrap = FALSE,
  verbose = TRUE)
}

\arguments{

  \item{obj}{A forest grow object of class \code{(rfsrc, grow)}.}

  \item{B}{Number of subsamples (or bootstrap iterations, if \code{bootstrap = TRUE}).}

  \item{block.size}{Number of trees in each block used when calculating VIMP. If VIMP is already included in the original grow object, that setting is used instead.}

  \item{importance}{Type of variable importance (VIMP) to compute. Choices are \code{"anti"}, \code{"permute"}, or \code{"random"}. If not specified, the default importance setting from the original grow call is used (if available).}

  \item{subratio}{Subsample size as a proportion of the original sample size. The default is approximately the inverse square root of the sample size.}

  \item{stratify}{Logical. If \code{TRUE}, uses stratified subsampling to preserve class balance. See \code{Details} for more information.}

  \item{performance}{Logical. If \code{TRUE}, calculates generalization error along with standard error and confidence intervals.}

  \item{performance.only}{Logical. If \code{TRUE}, only generalization error and its uncertainty are returned; VIMP is not computed.}

  \item{joint}{Logical. If \code{TRUE}, joint VIMP is computed for all variables. To calculate joint VIMP for a subset of variables, use \code{xvar.names}.}

  \item{xvar.names}{Character vector specifying variables to be used for joint VIMP. If omitted, all variables are included.}

  \item{bootstrap}{Logical. If \code{TRUE}, uses the double bootstrap instead of subsampling. This is typically slower but may provide more accurate uncertainty estimates.}

  \item{verbose}{Logical. If \code{TRUE}, prints progress updates during computation.}

}

\details{

This function applies subsampling (or optional double bootstrapping) to a previously trained forest to estimate standard errors and construct confidence intervals for variable importance (VIMP), as described in Ishwaran and Lu (2019). It also supports inference for the out-of-bag (OOB) prediction error via the \code{performance = TRUE} option. Joint VIMP for selected or all variables can be obtained using \code{joint} and \code{xvar.names}.

If the original forest does not include VIMP, it will be computed prior to subsampling. For repeated calls to \code{subsample}, it is recommended that VIMP be requested in the original \code{rfsrc} call. This not only avoids redundant computation, but also ensures consistency of the importance type (e.g., anti, permute, or random) and related parameters, which may otherwise be unclear. Note that permutation importance is \emph{not} the default for most families.

Subsampled forests are constructed using the same tuning parameters as the original forest. While most settings are automatically recovered, certain advanced configurations (e.g., custom sampling schemes) may not be fully supported.

Both subsampled variance estimates (Politis and Romano, 1994) and delete-\(d\) jackknife variance estimates (Shao and Wu, 1989) are returned. The jackknife estimator tends to produce larger standard errors, offering a conservative bias correction, particularly for signal variables.

By default, stratified subsampling is used for classification, survival, and competing risk families:
\itemize{
  \item For classification, strata correspond to class labels.
  \item For survival and competing risks, strata include event type and censoring.
}
Stratification helps ensure representation of key outcome types and is especially important for small sample sizes. Overriding this behavior is discouraged. Note that stratification is \emph{not} available for multivariate families, and caution should be exercised when subsampling in that context.

The function \code{\link{extract.subsample}} can be used to retrieve detailed information from the subsample object. By default, returned VIMP values are standardized: for regression families, VIMP is divided by the variance of the response; for other families, no transformation is applied. To obtain raw (unstandardized) values, set \code{standardize = FALSE}. For expert users, the option \code{raw = TRUE} returns detailed internal output, including VIMP from each individual subsampled forest (constructed on a smaller sample size), which is used internally by \code{\link{plot.subsample.rfsrc}} to generate confidence intervals.

Printed and plotted outputs also standardize VIMP by default. This behavior can be disabled via \code{standardize}. The \code{alpha} option controls the confidence level and is preset in wrapper functions but can be adjusted by the user.
}

\value{

  A list with the following key components:
  \item{rf}{Original forest grow object.}
  \item{vmp}{Variable importance values for grow forest.}
  \item{vmpS}{Variable importance subsampled values.}
  \item{subratio}{Subratio used.}  

}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{

  Ishwaran H. and Lu M.  (2019).  Standard errors and confidence
  intervals for variable importance in random forest regression,
  classification, and survival. \emph{Statistics in Medicine}, 38,
  558-582.

  Politis, D.N. and Romano, J.P. (1994). Large sample confidence
  regions based on subsamples under minimal assumptions. \emph{The
  Annals of Statistics}, 22(4):2031-2050.

  Shao, J. and Wu, C.J. (1989). A general theory for jackknife variance
  estimation. \emph{The Annals of Statistics}, 17(3):1176-1197.

}
\seealso{
  \command{\link{holdout.vimp.rfsrc}}
  \command{\link{plot.subsample.rfsrc}},
  \command{\link{rfsrc}},
  \command{\link{vimp.rfsrc}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## regression
## ------------------------------------------------------------

## training the forest
reg.o <- rfsrc(Ozone ~ ., airquality)

## default subsample call
reg.smp.o <- subsample(reg.o)

## plot confidence regions
plot.subsample(reg.smp.o)

## summary of results
print(reg.smp.o)

## joint vimp and confidence region for generalization error
reg.smp.o2 <- subsample(reg.o, performance = TRUE,
           joint = TRUE, xvar.names = c("Day", "Month"))
plot.subsample(reg.smp.o2)

## now try the double bootstrap (slower)
reg.dbs.o <- subsample(reg.o, B = 25, bootstrap = TRUE)
print(reg.dbs.o)
plot.subsample(reg.dbs.o)

## standard error and confidence region for generalization error only
gerror <- subsample(reg.o, performance.only = TRUE)
plot.subsample(gerror)

## ------------------------------------------------------------
## classification
## ------------------------------------------------------------

## 3 non-linear, 15 linear, and 5 noise variables 
if (library("caret", logical.return = TRUE)) {
  d <- twoClassSim(1000, linearVars = 15, noiseVars = 5)

  ## VIMP based on (default) misclassification error
  cls.o <- rfsrc(Class ~ ., d)
  cls.smp.o <- subsample(cls.o, B = 100)
  plot.subsample(cls.smp.o, cex.axis = .7)

  ## same as above, but with VIMP defined using normalized Brier score
  cls.o2 <- rfsrc(Class ~ ., d, perf.type = "brier")
  cls.smp.o2 <- subsample(cls.o2, B = 100)
  plot.subsample(cls.smp.o2, cex.axis = .7)
}

## ------------------------------------------------------------
## class-imbalanced data using RFQ classifier with G-mean VIMP
## ------------------------------------------------------------

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

  ## RFQ classifier
  oq <- imbalanced(Class ~ ., d, importance = TRUE, block.size = 10)

  ## subsample the RFQ-classifier
  smp.oq <- subsample(oq, B = 100)
  plot.subsample(smp.oq, cex.axis = .7)

}

## ------------------------------------------------------------
## survival 
## ------------------------------------------------------------

data(pbc, package = "randomForestSRC")
srv.o <- rfsrc(Surv(days, status) ~ ., pbc)
srv.smp.o <- subsample(srv.o, B = 100)
plot(srv.smp.o)

## ------------------------------------------------------------
## competing risks
## target event is death (event = 2)
## ------------------------------------------------------------

if (library("survival", logical.return = TRUE)) {
  data(pbc, package = "survival")
  pbc$id <- NULL
  cr.o <- rfsrc(Surv(time, status) ~ ., pbc, splitrule = "logrankCR", cause = 2)
  cr.smp.o <- subsample(cr.o, B = 100)
  plot.subsample(cr.smp.o, target = 2)
}

## ------------------------------------------------------------
## multivariate 
## ------------------------------------------------------------

if (library("mlbench", logical.return = TRUE)) {
  ## simulate the data 
  data(BostonHousing)
  bh <- BostonHousing
  bh$rm <- factor(round(bh$rm))
  o <- rfsrc(cbind(medv, rm) ~ ., bh)
  so <- subsample(o)
  plot.subsample(so)
  plot.subsample(so, m.target = "rm")
  ##generalization error
  gerror <- subsample(o, performance.only = TRUE)
  plot.subsample(gerror, m.target = "medv")
  plot.subsample(gerror, m.target = "rm")
}

## ------------------------------------------------------------
## largish data example - use rfsrc.fast for fast forests
## ------------------------------------------------------------

if (library("caret", logical.return = TRUE)) {
  ## largish data set
  d <- twoClassSim(1000, linearVars = 15, noiseVars = 5)

  ## use a subsampled forest with Brier score performance
  ## remember to set forest=TRUE for rfsrc.fast
  o <- rfsrc.fast(Class ~ ., d, ntree = 100,
           forest = TRUE, perf.type = "brier")
  so <- subsample(o, B = 100)
  plot.subsample(so, cex.axis = .7)
}


}
}
\keyword{subsampling}
\keyword{vimp}
\keyword{confidence interval}

