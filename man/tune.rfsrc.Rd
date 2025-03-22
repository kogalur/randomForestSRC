\name{tune.rfsrc}
\alias{tune.rfsrc}
\alias{tune}
\alias{tune.nodesize.rfsrc}
\alias{tune.nodesize}
\title{Tune Random Forest for the optimal mtry and nodesize parameters}
\description{
  Finds the optimal mtry and nodesize tuning parameter for a random
  forest using out-of-sample error.  Applies to all families. 
}
\usage{

\method{tune}{rfsrc}(formula, data,
  mtryStart = ncol(data) / 2,
  nodesizeTry = c(1:9, seq(10, 100, by = 5)), ntreeTry = 100,
  sampsize = function(x){min(x * .632, max(150, x ^ (3/4)))},
  nsplit = 1, stepFactor = 1.25, improve = 1e-3, strikeout = 3, maxIter = 25,
  trace = FALSE, doBest = FALSE, ...)

\method{tune.nodesize}{rfsrc}(formula, data,
  nodesizeTry = c(1:9, seq(10, 150, by = 5)), ntreeTry = 100,
  sampsize = function(x){min(x * .632, max(150, x ^ (4/5)))},
  nsplit = 1, trace = TRUE, ...)

}
\arguments{
  \item{formula}{A symbolic description of the model to be fit.}
  \item{data}{Data frame containing the y-outcome and x-variables.} 
  \item{mtryStart}{Starting value of mtry.}
  \item{nodesizeTry}{Values of nodesize optimized over.}
  \item{ntreeTry}{Number of trees used for the tuning step.}
  \item{sampsize}{Function specifying requested size of subsampled data.
    Can also be passed in as a number.}
  \item{nsplit}{Number of random splits used for splitting.}
  \item{stepFactor}{At each iteration, mtry is inflated (or deflated) by
    this value.}
  \item{improve}{The (relative) improvement in out-of-sample error must be by this
    much for the search to continue.}
  \item{strikeout}{The search is discontinued when the relative
    improvement in OOB error is negative.  However \code{strikeout}
    allows for some tolerance in this.  If a negative improvement is
    noted a total of \code{strikeout} times, the search is stopped.
    Increase this value only if you want an exhaustive search.}
  \item{maxIter}{The maximum number of iterations allowed for each mtry
    bisection search.}
  \item{trace}{Print the progress of the search?}
  \item{doBest}{Return a forest fit with the optimal mtry and nodesize parameters?}
  \item{...}{Further options to be passed to \code{\link{rfsrc.fast}}.}
}
\details{
  \code{tune} returns a matrix whose first and second
  columns contain the nodesize and mtry values searched and whose third
  column is the corresponding out-of-sample error.  Uses standardized error
  and in the case of multivariate forests it is the averaged
  standardized rror over the outcomes and for competing risks it is
  the averaged standardized error over the event types.

  If \code{doBest=TRUE}, also returns a forest object fit using the
  optimal \code{mtry} and \code{nodesize} values.

  All calculations (including the final optimized forest) are based on
  the fast forest interface \code{\link{rfsrc.fast}} which utilizes
  subsampling.  However, while this yields a fast optimization strategy,
  such a solution can only be considered approximate.  Users may wish to
  tweak various options to improve accuracy.  Increasing the default
  \code{sampsize} will definitely help. Increasing \code{ntreeTry}
  (which is set to 100 for speed) may also help.  It is also useful to
  look at contour plots of the out-of-sample error as a function of
  \code{mtry} and \code{nodesize} (see example below) to identify
  regions of the parameter space where error rate is small.

  \code{tune.nodesize} returns the optimal nodesize where optimization is
  over \code{nodesize} only.  

}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\seealso{
  \command{\link{rfsrc.fast}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## White wine classification example
## ------------------------------------------------------------

## load the data
data(wine, package = "randomForestSRC")
wine$quality <- factor(wine$quality)

## set the sample size manually
o <- tune(quality ~ ., wine, sampsize = 100)

## here is the optimized forest 
print(o$rf)

## visualize the nodesize/mtry OOB surface
if (library("interp", logical.return = TRUE)) {

  ## nice little wrapper for plotting results
  plot.tune <- function(o, linear = TRUE) {
    x <- o$results[,1]
    y <- o$results[,2]
    z <- o$results[,3]
    so <- interp(x=x, y=y, z=z, linear = linear)
    idx <- which.min(z)
    x0 <- x[idx]
    y0 <- y[idx]
    filled.contour(x = so$x,
                   y = so$y,
                   z = so$z,
                   xlim = range(so$x, finite = TRUE) + c(-2, 2),
                   ylim = range(so$y, finite = TRUE) + c(-2, 2),
                   color.palette =
                     colorRampPalette(c("yellow", "red")),
                   xlab = "nodesize",
                   ylab = "mtry",
                   main = "error rate for nodesize and mtry",
                   key.title = title(main = "OOB error", cex.main = 1),
                   plot.axes = {axis(1);axis(2);points(x0,y0,pch="x",cex=1,font=2);
                                points(x,y,pch=16,cex=.25)})
  }

  ## plot the surface
  plot.tune(o)

}

## ------------------------------------------------------------
## tuning for class imbalanced data problem
## - see imbalanced function for details
## - use rfq and perf.type = "gmean" 
## ------------------------------------------------------------

data(breast, package = "randomForestSRC")
breast <- na.omit(breast)
o <- tune(status ~ ., data = breast, rfq = TRUE, perf.type = "gmean")
print(o)


## ------------------------------------------------------------
## tune nodesize for competing risk - wihs data 
## ------------------------------------------------------------

data(wihs, package = "randomForestSRC")
plot(tune.nodesize(Surv(time, status) ~ ., wihs, trace = TRUE)$err)

}}
\keyword{forest}
\keyword{tune}


