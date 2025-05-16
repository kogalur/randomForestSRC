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
  trace = FALSE, doBest = TRUE, ...)

\method{tune.nodesize}{rfsrc}(formula, data,
  nodesizeTry = c(1:9, seq(10, 150, by = 5)), ntreeTry = 100,
  sampsize = function(x){min(x * .632, max(150, x ^ (4/5)))},
  nsplit = 1, trace = TRUE, ...)

}

\arguments{
  \item{formula}{A symbolic formula describing the model to be fit.}
  
  \item{data}{A data frame containing the response variable and predictor variables.}
  
  \item{mtryStart}{Initial value of \code{mtry} used to start the tuning search.}
  
  \item{nodesizeTry}{Vector of \code{nodesize} values over which tuning is performed.}
  
  \item{ntreeTry}{Number of trees used during the tuning step.}
  
  \item{sampsize}{Function specifying the size of the subsample. Can also be a numeric value.}
  
  \item{nsplit}{Number of random split points considered when splitting a node.}
  
  \item{stepFactor}{Multiplicative factor used to adjust \code{mtry} at each iteration.}
  
  \item{improve}{Minimum relative improvement in out-of-sample error required to continue the search.}
  
  \item{strikeout}{Number of consecutive non-improving steps (negative improvement) allowed before stopping the search. Increase to allow a more exhaustive search.}
  
  \item{maxIter}{Maximum number of iterations allowed for the \code{mtry} bisection search.}
  
  \item{trace}{If \code{TRUE}, prints progress during the search.}
  
  \item{doBest}{If \code{TRUE}, fits and returns a forest using the optimal \code{mtry} and \code{nodesize}.}
  
  \item{...}{Additional arguments passed to \code{\link{rfsrc.fast}}.}
}

\details{

\code{tune} returns a matrix with three columns: the first and second columns contain the \code{nodesize} and \code{mtry} values evaluated during the tuning process, and the third column contains the corresponding out-of-sample error.

The error is standardized. For multivariate forests, it is averaged over the outcomes; for competing risks, it is averaged over the event types.

If \code{doBest = TRUE}, the function also returns a forest object fit using the optimal \code{mtry} and \code{nodesize} values.

All tuning calculations, including the final optimized forest, are performed using the fast forest interface \code{\link{rfsrc.fast}}, which relies on subsampling. This makes the procedure computationally efficient but approximate. Users seeking more accurate tuning results may wish to adjust parameters such as:

\itemize{
  \item Increasing \code{sampsize}, which controls the size of the subsample used for tuning.
  \item Increasing \code{ntreeTry}, which defaults to 100 for speed.
}

It is also helpful to visualize the out-of-sample error surface as a function of \code{mtry} and \code{nodesize} using a contour plot (see example below) to identify regions of low error.

The function \code{tune.nodesize} performs a simplified search by optimizing only over \code{nodesize}.
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


