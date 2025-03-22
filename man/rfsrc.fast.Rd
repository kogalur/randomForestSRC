\name{rfsrc.fast}
\alias{rfsrc.fast}
\title{Fast Random Forests}
\description{
  Fast approximate random forests using subsampling with forest options
  set to encourage computational speed.  Applies to all families.
}
\usage{rfsrc.fast(formula, data,
  ntree = 500,
  nsplit = 10,
  bootstrap = "by.root",
  sampsize = function(x){min(x * .632, max(150, x ^ (3/4)))},
  samptype = "swor",
  samp = NULL,
  ntime = 50,
  forest = FALSE,
  save.memory = TRUE,
  ...)
}
\arguments{

  \item{formula}{Model to be fit. If missing, unsupervised splitting is
    implemented.}

  \item{data}{Data frame containing the y-outcome and x-variables.} 

  \item{ntree}{Number of trees.}

  \item{nsplit}{Non-negative integer value specifying number of 
    random split points used to split a node (deterministic splitting
    corresponds to the value zero and can be slower).}

  \item{bootstrap}{Bootstrap protocol used in growing a tree.}

  \item{sampsize}{Function specifying size of subsampled data. Can also be a number.}

  \item{samptype}{Type of bootstrap used.}

  \item{samp}{Bootstrap specification when \code{"by.user"} is used.}

  \item{ntime}{Integer value used for survival to
    constrain ensemble calculations to a grid of \code{ntime} time points.}

  \item{forest}{Save key forest values? Turn this on if you want prediction on test data.}

  \item{save.memory}{Save memory?  Setting this to \code{FALSE} stores
    terminal node quantities used for prediction on test data.  This
    yields rapid prediction but can be memory intensive for big data,
    especially competing risks and survival models.}
  
  \item{...}{Further arguments to be passed to \code{\link{rfsrc}}.}
  
}
\details{
  Calls \code{\link{rfsrc}} by choosing options (like subsampling) to
  encourage computational speeds.  This will provide a good
  approximation but will not be as good as default settings of
  \code{\link{rfsrc}}.
}
\value{
  An object of class \code{(rfsrc, grow)}.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\seealso{
  \command{\link{rfsrc}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## regression 
## ------------------------------------------------------------

## load the Iowa housing data
data(housing, package = "randomForestSRC")

## do quick and *dirty* imputation
housing <- impute(SalePrice ~ ., housing,
         ntree = 50, nimpute = 1, splitrule = "random")

## grow a fast forest
o1 <- rfsrc.fast(SalePrice ~ ., housing)
o2 <- rfsrc.fast(SalePrice ~ ., housing, nodesize = 1)
print(o1)
print(o2)

## grow a fast bivariate forest
o3 <- rfsrc.fast(cbind(SalePrice,Overall.Qual) ~ ., housing)
print(o3)

## ------------------------------------------------------------
## classification 
## ------------------------------------------------------------

data(wine, package = "randomForestSRC")
wine$quality <- factor(wine$quality)
o <- rfsrc.fast(quality ~ ., wine)
print(o)

## ------------------------------------------------------------
## grow fast random survival forests without C-calculation
## use brier score to assess model performance
## compare pure random splitting to logrank splitting
## ------------------------------------------------------------

data(peakVO2, package = "randomForestSRC")
f <- as.formula(Surv(ttodead, died)~.)
o1 <- rfsrc.fast(f, peakVO2, perf.type = "none")
o2 <- rfsrc.fast(f, peakVO2, perf.type = "none", splitrule = "random")
bs1 <- get.brier.survival(o1, cens.model = "km")
bs2 <- get.brier.survival(o2, cens.model = "km")
plot(bs2$brier.score, type = "s", col = 2)
lines(bs1$brier.score, type = "s", col = 4)
legend("bottomright", legend = c("random", "logrank"), fill = c(2,4))

## ------------------------------------------------------------
## competing risks
## ------------------------------------------------------------

data(wihs, package = "randomForestSRC")
o <- rfsrc.fast(Surv(time, status) ~ ., wihs)
print(o)

## ------------------------------------------------------------
## class imbalanced data using gmean performance
## ------------------------------------------------------------

data(breast, package = "randomForestSRC")
breast <- na.omit(breast)
f <- as.formula(status ~ .)
o <- rfsrc.fast(f, breast, perf.type = "gmean")
print(o)

## ------------------------------------------------------------
## class imbalanced data using random forests quantile-classifer (RFQ)
## fast=TRUE => rfsrc.fast
## see imbalanced function for further details
## ------------------------------------------------------------

data(breast, package = "randomForestSRC")
breast <- na.omit(breast)
f <- as.formula(status ~ .)
o <- imbalanced(f, breast, fast = TRUE)
print(o)

}}
\keyword{forest}
\keyword{fast}