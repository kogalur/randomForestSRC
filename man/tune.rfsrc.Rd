\name{tune.rfsrc}
\alias{tune.rfsrc}
\alias{tune}
\alias{tune.nodesize.rfsrc}
\alias{tune.nodesize}
\title{Tune Random Forest for optimal \code{mtry} and \code{nodesize}}

\description{

  Finds the optimal \code{mtry} and \code{nodesize} for a random forest
using out-of-bag (OOB) error.  Two search strategies are supported: a
grid-based search and a golden-section search with noise control.  Works
for all response families supported by \code{rfsrc.fast}.

}

\usage{

\method{tune}{rfsrc}(formula, data,
  mtry.start = ncol(data) / 2,
  nodesize.try = c(1:9, seq(10, 100, by = 5)), ntree.try = 100,
  sampsize = function(x) { min(x * .632, max(150, x^(3/4))) },
  nsplit = 1, step.factor = 1.25, improve = 1e-3, strikeout = 3, max.iter = 25,
  method = c("grid", "golden"),
  final.window = 5, reps.initial = 2, reps.final = 3,
  trace = FALSE, do.best = TRUE, seed = NULL, ...)

\method{tune.nodesize}{rfsrc}(formula, data,
  nodesize.try = c(1:9, seq(10, 150, by = 5)), ntree.try = 100,
  sampsize = function(x) { min(x * .632, max(150, x^(4/5))) },
  nsplit = 1, method = c("grid", "golden"),
  final.window = 5, reps.initial = 2, reps.final = 3, max.iter = 50,
  trace = TRUE, seed = NULL, ...)
}

\arguments{
  \item{formula}{A model formula.}
  \item{data}{A data frame with response and predictors.}
  \item{mtry.start}{Initial \code{mtry} for \code{tune}.}
  \item{nodesize.try}{Candidate \code{nodesize} values. Only values \eqn{\le} \code{floor(sampsize(n)/2)} are used.}
  \item{ntree.try}{Number of trees grown at each tuning evaluation.}
  \item{sampsize}{Function or numeric giving the per-tree subsample size. During tuning a single numeric size \code{ssize} is computed and passed to \code{rfsrc.fast}. If a vector is supplied (e.g., class specific), its total is used for \code{ssize}.}
  \item{nsplit}{Number of random split points to consider at each node.}
  \item{step.factor}{Multiplicative step-out factor over \code{mtry} for grid search in \code{tune}.}
  \item{improve}{Minimum relative improvement required to continue a search step in \code{tune}.}
  \item{strikeout}{Maximum number of consecutive non-improving steps allowed in \code{tune}.}
  \item{max.iter}{Maximum number of iterations for the step-out search in \code{tune} or the coordinate loop when \code{method = "golden"}.}
  \item{method}{Search strategy: \code{"grid"} (default) or \code{"golden"}.}
  \item{final.window}{For golden search, the terminal bracket width for the one-dimensional line search.}
  \item{reps.initial}{Replicates averaged at interior evaluations during golden iterations.}
  \item{reps.final}{Replicates averaged for each candidate during the final local sweep in golden search.}
  \item{trace}{If \code{TRUE}, prints progress.}
  \item{do.best}{If \code{TRUE}, \code{tune} fits and returns a forest at the optimal pair.}
  \item{seed}{Optional integer for reproducible tuning. The holdout split (when used) and all tuning fits become deterministic for a given seed.}
  \item{...}{Additional arguments passed to \code{\link{rfsrc.fast}}. Arguments that control tuning itself (\code{perf.type}, \code{forest}, \code{save.memory}, \code{ntree}, \code{mtry}, \code{nodesize}, \code{sampsize}, \code{nsplit}) are managed internally.}
}

\details{

\emph{Error estimate.} If \code{2 * ssize < n}, a disjoint holdout of
size \code{ssize} is used for evaluation; otherwise OOB error is
used.


\emph{Subsample used during tuning.} Both functions derive a single
integer \code{ssize} from \code{sampsize} and pass it to
\code{\link{rfsrc.fast}} for all tuning fits. This improves stability
and comparability across candidates. When \code{do.best = TRUE} in
\code{tune}, the final forest is fit with the user-supplied
\code{sampsize} exactly as provided.


\emph{Grid search.} \code{tune} performs a step-out search over
\code{mtry} for each \code{nodesize} in \code{nodesize.try}, using
\code{step.factor}, \code{improve}, \code{strikeout}, and
\code{max.iter}. \code{tune.nodesize} evaluates the supplied
\code{nodesize.try} grid directly.


\emph{Golden search.} Uses a guarded golden-section line search with
noise control. For each one-dimensional search (over \code{nodesize} or
\code{mtry}), the routine probes a small left-anchor grid \code{1:9},
iterates golden shrinkage until the bracket width is at most
\code{final.window}, then runs a short local sweep with
\code{reps.final} replicates. In \code{tune} the searches over
\code{nodesize} and \code{mtry} alternate in a simple coordinate loop,
with \code{improve} and \code{strikeout} as stopping controls.


}

\value{
For \code{tune}:
\itemize{
  \item \code{results}: matrix with columns \code{nodesize}, \code{mtry}, \code{err}.
  \item \code{optimal}: named numeric vector \code{c(nodesize = ..., mtry = ...)}.
  \item \code{rf}: fitted forest at the optimum if \code{do.best = TRUE}.
}

For \code{tune.nodesize}:
\itemize{
  \item \code{nsize.opt}: optimal \code{nodesize}.
  \item \code{err}: data frame with columns \code{nodesize} and \code{err}.
}
}

\author{
Hemant Ishwaran and Udaya B. Kogalur
}

\seealso{
\code{\link{rfsrc.fast}}
}

\examples{
\donttest{
## ------------------------------------------------------------
## White wine classification example
## ------------------------------------------------------------
data(wine, package = "randomForestSRC")
wine$quality <- factor(wine$quality)

## Fixed seed makes tuning reproducible
set.seed(1)

## Full tuner over nodesize and mtry (grid)
o1 <- tune(quality ~ ., wine, sampsize = 100, method = "grid")
print(o1$optimal)

## Golden search alternative
o2 <- tune(quality ~ ., wine, sampsize = 100, method = "golden",
           reps.initial = 2, reps.final = 3, seed = 1)
print(o2$optimal)

## visualize the nodesize/mtry surface
if (library("interp", logical.return = TRUE)) {

  plot.tune <- function(o, linear = TRUE) {
    x <- o$results[, 1]
    y <- o$results[, 2]
    z <- o$results[, 3]
    so <- interp(x = x, y = y, z = z, linear = linear)
    idx <- which.min(z)
    x0 <- x[idx]; y0 <- y[idx]
    filled.contour(x = so$x, y = so$y, z = so$z,
                   xlim = range(so$x, finite = TRUE) + c(-2, 2),
                   ylim = range(so$y, finite = TRUE) + c(-2, 2),
                   color.palette = colorRampPalette(c("yellow", "red")),
                   xlab = "nodesize", ylab = "mtry",
                   main = "error rate for nodesize and mtry",
                   key.title = title(main = "OOB error", cex.main = 1),
                   plot.axes = {
                     axis(1); axis(2)
                     points(x0, y0, pch = "x", cex = 1, font = 2)
                     points(x, y, pch = 16, cex = .25)
                   })
  }

  plot.tune(o1)
  plot.tune(o2)
}

## ------------------------------------------------------------
## nodesize only: grid vs golden
## ------------------------------------------------------------
o3 <- tune.nodesize(quality ~ ., wine, sampsize = 100, method = "grid",
                    trace = TRUE, seed = 1)
o4 <- tune.nodesize(quality ~ ., wine, sampsize = 100, method = "golden",
                    reps.initial = 2, reps.final = 3, trace = TRUE, seed = 1)
plot(o3$err, type = "s", xlab = "nodesize", ylab = "error")

## ------------------------------------------------------------
## Tuning for class imbalance (rfq with geometric mean performance)
## ------------------------------------------------------------
data(breast, package = "randomForestSRC")
breast <- na.omit(breast)
o5 <- tune(status ~ ., data = breast, rfq = TRUE, perf.type = "gmean",
           method = "golden", seed = 1)
print(o5$optimal)

## ------------------------------------------------------------
## Competing risks example (nodesize only)
## ------------------------------------------------------------
data(wihs, package = "randomForestSRC")
plot(tune.nodesize(Surv(time, status) ~ ., wihs, trace = TRUE)$err, type = "s")
}}

\keyword{forest}
\keyword{tune}
