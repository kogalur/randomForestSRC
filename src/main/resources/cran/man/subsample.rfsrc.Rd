\name{subsample.rfsrc}
\alias{subsample.rfsrc}
\alias{subsample}
\alias{extract.subsample}
\alias{extract.bootsample}
\title{Subsample Forests for VIMP Confidence Intervals}
\description{
  Use subsampling to calculate confidence intervals and standard
  errors for VIMP (variable importance). Applies to all families.
}
\usage{\method{subsample}{rfsrc}(obj,
  B = 100,
  block.size = 1,
  subratio = NULL,
  stratify = TRUE,
  joint = FALSE,
  bootstrap = FALSE,
  verbose = TRUE)
}
\arguments{
  \item{obj}{A forest grow object.}
  \item{B}{Number of subsamples (or number of bootstraps).}
   \item{block.size}{Specifies number of trees in a block when calculating
    VIMP.  This is over-ridden if VIMP is present in the original grow
    call in which case the grow value is used.} 
  \item{subratio}{Ratio of subsample size to original sample size.  The
    default is the inverse square root of the sample size.}
  \item{stratify}{Use stratified subsampling?  See details below.}
  \item{joint}{Include the VIMP for all variables jointly perturbed?
    This is useful reference problems where one might be suspicious that
    many (or all) variables are noise.}  
  \item{bootstrap}{Use double bootstrap approach in place of
    subsampling?  Much slower, but potentially more accurate.}
  \item{verbose}{Provide verbose output?}
}
\details{
  Given a forest object, subsamples the forest to obtain standard errors
  and confidence intervals for VIMP (Ishwaran and Lu, 2019).  If
  bootstrapping is requested, then the double bootstrap is applied in
  place of subsampling.

  If VIMP is not present in the original forest object, the algorithm
  will first need to calculate VIMP.  Therefore, if the user plans to
  make repeated calls to \code{subsample}, it is advisable to include
  VIMP in the original grow call.  Subsampled forests are calculated
  using the same tuning parameters as the original forest.  While a
  sophisticated algorithm is utilized to acquire as many of these
  parameters as possible, keep in mind there are some conditions where
  this will fail: for example there are certain settings where the user
  has specified non-standard sampling in the grow forest.
  
  Delete-d jackknife estimators (Shao and Wu, 1989) are returned along
  with subsampling estimators (Politis and Romano, 1994).  While these
  two methods are closely related, standard errors for delete-d
  estimators are generally larger than the subsampled estimates, which
  is a form of bias correction, which occurs primarily for variables
  with true signal.  Confidence interval coverage is generally better
  under delete-d estimators.  Note that undercoverage for strong
  variables and overcoverage for noise variables exhibited by both
  estimators may be beneficial if the goal is variable selection
  (Ishwaran and Lu, 2019).

  By default, stratified subsampling is used for classification,
  survival, and competing risk families.  For classification,
  stratification is on the class label, while for survival and competing
  risk, stratification is on the event type and censoring.  Users are
  discouraged from over-riding this option, especially in small sample
  settings, as this could lead to error due to subsampled data not
  having full representation of class labels in classification settings,
  and in survival settings, subsampled data may be devoid of deaths
  and/or have reduced number of competing risks.  Note also that
  stratified sampling is not available for multivariate families --
  users should especially exercise caution when selecting subsampling
  rates here.

  The function \command{extract.subsample} conveniently extracts summary
  information from the subsampled object. It parses objects \code{rf},
  \code{vmp} and \code{vmpS} returned by the function (see below for
  what these are).

  When printing and or plotting results, the default setting is to
  standardize VIMP, where for regression families, VIMP is standardized
  by dividing by the variance and multiplying by 100.  For all other
  families, VIMP is scaled by 100.  This can be turned off using the
  option \code{standardize} in those wrappers.
  
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
## regression example
## ------------------------------------------------------------

## grow the forest - request VIMP
reg.o <- rfsrc(mpg ~ ., mtcars)

## very small sample size so need largish subratio
reg.smp.o <- subsample(reg.o, B = 100, subratio = .5)

## plot confidence regions
plot.subsample(reg.smp.o)

## summary of results
print(reg.smp.o)

## now try the double bootstrap (slow!!)
reg.dbs.o <- subsample(reg.o, B = 100, bootstrap = TRUE)
print(reg.dbs.o)
plot.subsample(reg.dbs.o)

## ------------------------------------------------------------
## classification example
## ------------------------------------------------------------

## 3 non-linear, 15 linear, and 5 noise variables 
if (library("caret", logical.return = TRUE)) {
  d <- twoClassSim(1000, linearVars = 15, noiseVars = 5)

  ## VIMP based on (default) misclassification error
  cls.o <- rfsrc(Class ~ ., d)
  cls.smp.o <- subsample(cls.o, B = 100)
  plot.subsample(cls.smp.o, cex = .7)

  ## same as above, but with VIMP defined using normalized Brier score
  cls.o2 <- rfsrc(Class ~ ., d, perf.type = "brier")
  cls.smp.o2 <- subsample(cls.o2, B = 100)
  plot.subsample(cls.smp.o2, cex = .7)
}

## ------------------------------------------------------------
## survival example
## ------------------------------------------------------------

data(pbc, package = "_PROJECT_PACKAGE_NAME_")
srv.o <- rfsrc(Surv(days, status) ~ ., pbc)
srv.smp.o <- subsample(srv.o, B = 100)
plot.subsample(srv.smp.o)

## ------------------------------------------------------------
## competing risk example
## target event is death (event = 2)
## ------------------------------------------------------------

if (library("survival", logical.return = TRUE)) {
  data(pbc, package = "survival")
  pbc$id <- NULL
  cr.o <- rfsrc(Surv(time, status) ~ ., pbc, splitrule = "logrank", cause = 2)
  cr.smp.o <- subsample(cr.o, B = 100)
  plot.subsample(cr.smp.o, target = 2)
}

## ------------------------------------------------------------
## multivariate family
## ------------------------------------------------------------

if (library("mlbench", logical.return = TRUE)) {
  ## simulate the data 
  data(BostonHousing)
  bh <- BostonHousing
  bh$rm <- factor(round(bh$rm))
  o <- rfsrc(cbind(medv, rm) ~ ., bh)
  so <- subsample(o)
  plot(so)
  plot(so, m.target = "rm")
}

## ------------------------------------------------------------
## largish data example - use rfsrc.fast for fast forests
## ------------------------------------------------------------

if (library("caret", logical.return = TRUE)) {
  ## largish data set
  d <- twoClassSim(1000, linearVars = 15, noiseVars = 5)

  ## use a subsampled forest with Brier score performance
  o <- rfsrc.fast(Class ~ ., d, ntree = 100, perf.type = "brier")
  so <- subsample(o, B = 100)
  plot.subsample(so, cex = .7)
}


}
}
\keyword{subsampling}
\keyword{vimp}
\keyword{confidence interval}

