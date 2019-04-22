\name{tune.rfsrc}
\alias{tune.rfsrc}
\alias{tune}
\alias{tune.nodesize.rfsrc}
\alias{tune.nodesize}
\title{Tune Random Forest for the optimal mtry and nodesize parameters}
\description{
  Finds the optimal mtry and nodesize tuning parameter for a random
  forest using out-of-bag (OOB) error.  Applies to all families. 
}
\usage{

\method{tune}{rfsrc}(formula, data,
  mtryStart = ncol(data) / 2,
  nodesizeTry = c(1:9, seq(10, 100, by = 5)), ntreeTry = 50,
  sampsize = function(x){min(x * .632, max(150, x ^ (3/4)))},
  nsplit = 10, stepFactor = 1.25, improve = 1e-3, strikeout = 3, maxIter = 25,
  trace = FALSE, doBest = TRUE, ...)

\method{tune.nodesize}{rfsrc}(formula, data,
  nodesizeTry = c(1:9, seq(10, 150, by = 5)),
  sampsize = function(x){min(x * .632, max(150, x ^ (4/5)))},
  nsplit = 1, trace = FALSE, ...)

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
  \item{improve}{The (relative) improvement in OOB error must be by this
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
  column is the corresponding OOB error.  Uses standardized OOB error
  and in the case of multivariate forests it is the averaged
  standardized OOB error over the outcomes and for competing risks it is
  the averaged standardized OOB error over the event types.

  If \code{doBest=TRUE}, also returns a forest object fit using the
  optimal \code{mtry} and \code{nodesize} values.

  All calculations (including the final optimized forest) are based on
  the fast forest interface \code{\link{rfsrc.fast}} which utilizes
  subsampling.  However, while this yields a fast optimization strategy,
  such a solution can only be considered approximate.  Users may wish to
  tweak various options to improve accuracy.  Increasing the default
  \code{sampsize} will definitely help. Increasing \code{ntreeTry}
  (which is set to 50 for speed) may also help.  It is also useful to
  look at contour plots of the OOB error as a function of \code{mtry}
  and \code{nodesize} (see example below) to identify regions of the
  parameter space where error rate is small.

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
data(wine, package = "_PROJECT_PACKAGE_NAME_")
wine$quality <- factor(wine$quality)

## default tuning call
o <- tune(quality ~ ., wine)

## here is the optimized forest 
print(o$rf)

## visualize the nodesize/mtry OOB surface
if (library("akima", logical.return = TRUE)) {

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
                   main = "OOB error for nodesize and mtry",
                   key.title = title(main = "OOB error", cex.main = 1),
                   plot.axes = {axis(1);axis(2);points(x0,y0,pch="x",cex=1,font=2);
                                points(x,y,pch=16,cex=.25)})
  }

  ## plot the surface
  plot.tune(o)

}

## ------------------------------------------------------------
## tune nodesize for competing risk - wihs data 
## ------------------------------------------------------------

data(wihs, package = "_PROJECT_PACKAGE_NAME_")
plot(tune.nodesize(Surv(time, status) ~ ., wihs, trace = TRUE)$err)

}}
\keyword{forest}
\keyword{tune}


