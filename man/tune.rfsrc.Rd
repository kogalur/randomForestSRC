\name{tune.rfsrc}
\alias{tune.rfsrc}
\alias{tune}
\alias{tune.nodesize.rfsrc}
\alias{tune.nodesize}
\title{Tune Random Forest Parameters}

\description{

  Select \code{mtry} and \code{nodesize} for a supervised random forest
  using holdout or out-of-bag (OOB) prediction error. \code{tune()} searches
  over both parameters; \code{tune.nodesize()} searches over
  \code{nodesize} only. Both functions offer grid and golden-section
  search and fit the tuning forests with \code{\link{rfsrc.fast}}.

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
  \item{formula}{A model formula specifying the response and predictors.}
  \item{data}{A data frame containing the response and predictors.}
  \item{mtry.start}{Starting value of \code{mtry} for \code{tune()}.}
  \item{nodesize.try}{Candidate terminal-node sizes for grid search,
    subject to the limits in \bold{Note}. Not used for golden search.}
  \item{ntree.try}{Number of trees in each tuning forest, excluding the
    final forest returned by \code{tune()}.}
  \item{sampsize}{Per-tree sample size, supplied as a function of the
    sample size or as a number. For a numeric vector, such as class-specific
    sizes, the total is used during tuning.}
  \item{nsplit}{Number of random split points considered at each node.}
  \item{step.factor}{Factor for increasing or decreasing \code{mtry}
    during grid search in \code{tune()}. Must be greater than one.}
  \item{improve}{Relative-error reduction threshold for \code{tune()}.
    Must be nonnegative; \code{1e-3} corresponds to a 0.1 percent decrease.}
  \item{strikeout}{Number of additional steps allowed after error
    increases in grid search, counted across each search direction.
    For golden search in \code{tune()}, the number of consecutive
    alternating-search iterations without sufficient improvement before
    stopping.}
  \item{max.iter}{Maximum iterations per search direction (grid) or
    per one-dimensional search (golden). For \code{tune()} with golden
    search, also limits the number of alternating searches.}
  \item{method}{Search method: \code{"grid"} (default) or
    \code{"golden"}. See \bold{Details}.}
  \item{final.window}{Target interval width for golden search. Also
    sets the distance searched on either side of the best value from
    the initial small-value comparisons.}
  \item{reps.initial}{Number of replicate forests per setting during
    golden-section search; their errors are averaged.}
  \item{reps.final}{Number of replicate forests per setting in the final
    integer comparisons for golden search; their errors are averaged.}
  \item{trace}{If \code{TRUE}, print tuning progress and errors.}
  \item{do.best}{If \code{TRUE}, \code{tune()} also fits a forest on all
    available data using the selected settings.}
  \item{seed}{Optional integer seed for reproducible tuning, including
    the holdout sample when used.}
  \item{...}{Additional options for \code{\link{rfsrc.fast}}, such as
    \code{splitrule} and \code{perf.type}. In \code{tune()},
    \code{ntree}, \code{mtry}, \code{nodesize}, \code{forest}, and
    \code{save.memory} are controlled by the tuner; values supplied
    here are ignored.}
}

\details{

  Tuning compares forest settings to find those with the lowest
  prediction error. The parameter \code{mtry} controls how many
  predictors are considered at each split, and \code{nodesize} controls
  terminal-node size. A forest is grown for each candidate setting and
  its prediction error is calculated. The selected setting has the
  smallest error in the returned results table.

  \subsection{Choosing a search method}{

    \describe{
      \item{\code{method = "grid"}:}{
        For each eligible \code{nodesize.try} value, \code{tune()}
        starts at \code{mtry.start} and searches toward smaller and
        larger \code{mtry} values. \code{step.factor} sets the spacing;
        \code{improve}, \code{strikeout}, and \code{max.iter} control
        when to stop. The grid specifies \code{nodesize} values, while
        \code{mtry} values are chosen during the search.
        \code{tune.nodesize()} evaluates the \code{nodesize.try} grid
        directly, subject to the sample-size limit in \bold{Note}.
      }
      \item{\code{method = "golden"}:}{
        Comparisons within an interval progressively narrow the search
        toward values with smaller error. \code{tune()} alternates
        between searching over \code{nodesize} with \code{mtry} fixed
        and searching over \code{mtry} with \code{nodesize} fixed.
        \code{tune.nodesize()} searches over \code{nodesize} only.

        Each one-dimensional search also checks \code{1:9} within its
        limits. It finishes with individual integer comparisons near
        the best small value or within the remaining interval. Averaging
        replicate forests, controlled by \code{reps.initial} and
        \code{reps.final}, reduces the effect of random variation in
        forest fits. \code{nodesize.try} and \code{step.factor} are not
        used.
      }
    }

  }

  \subsection{Comparing prediction errors}{

    Both functions use one per-tree sample size, \code{ssize}, throughout
    tuning. Let \code{n} be the sample size after preprocessing. When
    \code{2 * ssize < n}, a random holdout of \code{ssize} observations
    is set aside. Forests are grown on the remaining data, and all
    candidates are evaluated on this same holdout. Otherwise, forests
    use all available data and are evaluated by OOB error.

    Supply \code{perf.type} to choose the tuning error measure.
    \code{tune()} averages the available errors from
    \code{\link{get.mv.error}} with \code{standardize = TRUE}:
    regression errors are divided by the response variance, and errors
    for multivariate responses are averaged.

  }

  \subsection{Fitting the final forest}{

    With \code{do.best = TRUE}, \code{tune()} fits a new forest on all
    available data, using the selected \code{mtry} and \code{nodesize},
    the original \code{sampsize} argument, and the specified
    \code{nsplit}. The forest is returned as \code{rf}.

    This final fit uses the default tree count and performance measure
    of \code{rfsrc.fast}. To choose these yourself, set
    \code{do.best = FALSE} and fit a forest separately with the returned
    \code{optimal} values. \code{tune.nodesize()} returns its selected
    node size and error table without a final forest fit.

  }

}

\value{
  \code{tune()} returns a list with components:
  \describe{
    \item{results}{A numeric matrix with columns \code{nodesize},
      \code{mtry}, and \code{err} (prediction error). For golden search,
      this table contains the final integer comparisons; repeated
      parameter pairs are combined by averaging their errors.}
    \item{optimal}{A named numeric vector giving the \code{nodesize}
      and \code{mtry} with the smallest error in \code{results}.}
    \item{rf}{The final fitted forest when \code{do.best = TRUE};
      otherwise \code{NULL}.}
  }

  \code{tune.nodesize()} returns a list with components:
  \describe{
    \item{nsize.opt}{The selected \code{nodesize}.}
    \item{err}{A data frame with columns \code{nodesize} and \code{err},
      giving the candidate values and their prediction errors.}
  }
}

\note{
  The tuning sample size \code{ssize} is obtained by evaluating
  \code{sampsize(n)}, using a numeric scalar directly, or summing a
  numeric vector. In \code{tune()}, it is rounded and bounded between
  2 and \code{n}.

  For \code{tune()}, \code{mtry} ranges from 1 to the number of
  predictors, and the upper \code{nodesize} limit is
  \code{max(10, floor(ssize / 2))}. Grid values above this limit are
  excluded. For \code{tune.nodesize()}, the upper grid limit is half
  the per-tree sample size, rounded down.
}

\author{
Hemant Ishwaran and Udaya B. Kogalur
}

\seealso{
  \code{\link{rfsrc.fast}},
  \code{\link{rfsrc}},
  \code{\link{get.mv.error}}
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
                   key.title = title(main = "Tuning error", cex.main = 1),
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
