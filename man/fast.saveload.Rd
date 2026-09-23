\name{fast.saveload}
\alias{fast.saveload}
\alias{fast.save}
\alias{fast.load}
\alias{fast.save.list}
\alias{fast.load.list}
\title{Fast Saving and Loading of Random Forests}
\description{
  Save and load a minimal forest object for prediction.  These functions
  are designed to use substantially less disk space and load large forests
  faster than saving the full fitted object with \code{\link[base]{saveRDS}}
  and loading it with \code{\link[base]{readRDS}}.  Functions are available
  for a single forest or a list of forests.
}
\usage{
fast.save(o, path = NULL, testing = TRUE, units = "Mb", compact = FALSE)

fast.load(directory, path = NULL, testing = FALSE, units = "Mb")

fast.save.list(o, path = NULL, testing = FALSE, units = "Mb", compact = FALSE)

fast.load.list(directory, path = NULL, testing = FALSE, units = "Mb")
}
\arguments{
  \item{o}{For \code{fast.save()}, a fitted object with class beginning
    \code{c("rfsrc", "grow")} and a \code{forest} component.
    Supply the fitted object, not \code{o$forest}.  For
    \code{fast.save.list()}, a nonempty list of these objects.}
  \item{path}{For saving, the destination directory, including its name.
    The default is \code{file.path(getwd(), "forest")}.
    An existing destination and all its contents are deleted.
    For loading, the parent directory containing \code{directory};
    \code{NULL} uses the current working directory.}
  \item{directory}{Name of the saved directory relative to \code{path}.
    For \code{fast.load.list()}, this directory contains the individual
    forest subdirectories.}
  \item{testing}{Print component memory sizes using \code{lsos()}.
    Defaults to \code{TRUE} for \code{fast.save()} and \code{FALSE}
    for the other functions.}
  \item{units}{Units for the memory-size display, passed to
    \code{format(object.size(...), units = units)}.
    The default is \code{"Mb"}.}
  \item{compact}{If \code{TRUE}, omit terminal-node information from the
    saved forest and recompute it during prediction.  Requires a
    non-anonymous forest containing the training predictors and responses.
    If \code{FALSE} (default), save the forest with its available
    terminal-node information.  See \emph{Compact saving} below.}
}
\details{
  \subsection{Saving and loading a forest}{
    \code{fast.save()} addresses the large files and slow loading that
    can result from saving a full fitted forest with \code{saveRDS()}.
    It saves only \code{o$forest}, the minimal object used for prediction,
    including its tree information and any training data it contains.
    Predictions and performance summaries outside this component are
    excluded.  \code{fast.load()} returns the forest for use with
    \code{\link{predict.rfsrc}}.

    Large tabular components are saved efficiently with
    \code{fst::write_fst()}, and the remaining forest information with
    \code{saveRDS(..., compress = FALSE)}.  Saving requires the \pkg{fst}
    and \pkg{data.table} packages; loading requires \pkg{fst}.
    All files have a \file{.rda} extension despite the different formats.
    Keep the directory's files together and load the forest with
    \code{fast.load()}.
  }
  \subsection{Compact saving}{
    With \code{compact = TRUE}, the saved forest includes the tree
    structure, factor-split information, training data, seeds, and other
    forest information, but omits stored terminal-node information.
    Prediction recomputes this information from the saved forest and
    training data without growing new trees.  The fitted object in memory
    is unchanged.

    The saved copy has \code{terminal.qualts = FALSE} and
    \code{terminal.quants = FALSE}, with \code{nativeArrayTNDS} removed.
    Use the same \code{fast.load()} or \code{fast.load.list()} call
    for either saving mode.  Compact saving can further reduce disk space,
    but prediction may take longer and peak memory use may not decrease.

    Use \code{compact = FALSE} for anonymous forests.  Compact saving
    includes the training data.  Eligibility checks run before the
    destination is deleted; for lists, every forest is checked first.
  }
  \subsection{Directories}{
    \strong{Both save functions delete an existing destination directory
    and all its contents.}  Use a dedicated directory for each forest or
    list.  The destination and any missing parent directories are then
    created.

    When saving, \code{path} is the complete destination.  When loading,
    \code{directory} is its name and \code{path} is its parent.  For a
    destination \code{save.path}, use
    \code{fast.load(basename(save.path), path = dirname(save.path))}.
    For the default destination, use \code{fast.load("forest")}.
  }
  \subsection{Lists of forests}{
    \code{fast.save.list()} saves forests in subdirectories
    \file{forest1}, \file{forest2}, and so on, inside \code{path}.
    List names are not saved.  \code{fast.load.list()} loads every entry
    listed by \code{list.files()} in that directory, so it should contain
    only saved forest subdirectories.

    The loaded list is unnamed and follows directory-listing order,
    which can differ from the original order: \file{forest10} can precede
    \file{forest2}.  To load one forest, call \code{fast.load()} with its
    subdirectory name.
  }
}
\value{
  \describe{
    \item{\code{fast.save}}{Saves the forest to disk and returns the
      memory-usage matrix produced by \code{gc(FALSE)}.}
    \item{\code{fast.load}}{The saved forest component with its original
      class, usually \code{c("rfsrc", "forest", family)}.}
    \item{\code{fast.save.list}}{An unnamed list of the memory-usage
      matrices returned by the individual \code{fast.save()} calls.}
    \item{\code{fast.load.list}}{An unnamed list of forest objects in
      directory-listing order.}
  }
}
\note{
  The sampling-size function is saved as
  \code{sampfrac = o$forest$sampsize(1)} and loaded as
  \code{function(x) x * sampfrac}.  Use \code{saveRDS()} to save a custom
  function that is not proportional to sample size.

  Compact saving checks that \code{round(sampsize(n))} equals
  \code{round(n * sampsize(1))} at the training sample size \code{n}.
  A mismatch stops saving before the destination is deleted.
  This check applies only at the training sample size.
}
\seealso{
  \code{\link{rfsrc}}, \code{\link{predict.rfsrc}},
  \code{\link[base]{saveRDS}}, \code{\link[base]{readRDS}}
}
\examples{
## ------------------------------------------------------------
## Regression: save, load, and compare predictions
## ------------------------------------------------------------
if (requireNamespace("fst", quietly = TRUE) &&
    requireNamespace("data.table", quietly = TRUE)) {
  o <- rfsrc(mpg ~ ., data = mtcars)
  print(o)
  save.path <- tempfile("rfsrc-forest-")
  fast.save(o, path = save.path, testing = FALSE)
  oo <- fast.load(basename(save.path), path = dirname(save.path))
  p <- predict(o)
  pp <- predict(oo)
  print(summary(p$predicted - pp$predicted))
  print(summary(p$predicted.oob - pp$predicted.oob))
  unlink(save.path, recursive = TRUE)
}

\donttest{
## ------------------------------------------------------------
## Compact saving: the same forest, fewer saved components
## ------------------------------------------------------------
if (requireNamespace("fst", quietly = TRUE) &&
    requireNamespace("data.table", quietly = TRUE)) {
  o <- rfsrc(mpg ~ ., data = mtcars, ntree = 100)
  set.seed(19)
  reference <- predict(o, seed = -19)
  save.path <- tempfile("rfsrc-compact-")
  fast.save(o, path = save.path, compact = TRUE, testing = FALSE)
  oo <- fast.load(basename(save.path), path = dirname(save.path))
  print(oo[c("terminal.qualts", "terminal.quants")])
  print(list.files(save.path, pattern = "^nativeArrayTDNS_"))
  set.seed(19)
  restored <- predict(oo, seed = -19)
  print(all.equal(reference$predicted, restored$predicted))
  print(all.equal(reference$predicted.oob, restored$predicted.oob))
  unlink(save.path, recursive = TRUE)
}

## ------------------------------------------------------------
## Regression: a list of forests with different node sizes
## ------------------------------------------------------------
if (requireNamespace("fst", quietly = TRUE) &&
    requireNamespace("data.table", quietly = TRUE)) {
  o1 <- rfsrc(mpg ~ ., data = mtcars, nodesize = 1)
  o2 <- rfsrc(mpg ~ ., data = mtcars, nodesize = 10)
  print(o1)
  print(o2)
  models <- list(o1, o2)
  save.path <- tempfile("rfsrc-forest-list-")
  invisible(fast.save.list(models, path = save.path))
  oo <- fast.load.list(basename(save.path), path = dirname(save.path))
  print(predict(oo[[1]]))
  print(predict(oo[[2]]))
  unlink(save.path, recursive = TRUE)
}

## ------------------------------------------------------------
## RFQ for imbalanced classification
## ------------------------------------------------------------
## Use matching prediction seeds when comparing class labels.
if (requireNamespace("fst", quietly = TRUE) &&
    requireNamespace("data.table", quietly = TRUE)) {
  data(breast, package = "randomForestSRC")
  dta <- na.omit(breast)
  o <- imbalanced(status ~ ., data = dta, ntree = 100)
  print(o)
  save.path <- tempfile("rfsrc-forest-")
  fast.save(o, path = save.path, testing = FALSE)
  oo <- fast.load(basename(save.path), path = dirname(save.path))
  set.seed(19)
  p <- predict(o, seed = -19)
  set.seed(19)
  pp <- predict(oo, seed = -19)
  print(summary(p$predicted - pp$predicted))
  print(summary(p$predicted.oob - pp$predicted.oob))
  print(all.equal(as.character(p$class), as.character(pp$class)))
  print(all.equal(as.character(p$class.oob),
                  as.character(pp$class.oob)))
  unlink(save.path, recursive = TRUE)
}

## ------------------------------------------------------------
## Binary classification with rfq = TRUE
## ------------------------------------------------------------
if (requireNamespace("fst", quietly = TRUE) &&
    requireNamespace("data.table", quietly = TRUE)) {
  data(breast, package = "randomForestSRC")
  dta <- na.omit(breast)
  o <- rfsrc(status ~ ., data = dta, rfq = TRUE, ntree = 100,
             perf.type = "gmean", splitrule = "auc")
  print(o)
  save.path <- tempfile("rfsrc-forest-")
  fast.save(o, path = save.path, testing = FALSE)
  oo <- fast.load(basename(save.path), path = dirname(save.path))
  set.seed(19)
  p <- predict(o, seed = -19)
  set.seed(19)
  pp <- predict(oo, seed = -19)
  print(summary(p$predicted - pp$predicted))
  print(summary(p$predicted.oob - pp$predicted.oob))
  print(all.equal(as.character(p$class), as.character(pp$class)))
  print(all.equal(as.character(p$class.oob),
                  as.character(pp$class.oob)))
  unlink(save.path, recursive = TRUE)
}

## ------------------------------------------------------------
## Anonymous RFQ: supply the same prediction data to both forests
## ------------------------------------------------------------
if (requireNamespace("fst", quietly = TRUE) &&
    requireNamespace("data.table", quietly = TRUE)) {
  data(breast, package = "randomForestSRC")
  dta <- na.omit(breast)
  o <- rfsrc.anonymous(status ~ ., data = dta, rfq = TRUE,
                       ntree = 100, perf.type = "gmean", splitrule = "auc")
  print(o)
  save.path <- tempfile("rfsrc-forest-")
  fast.save(o, path = save.path, testing = FALSE)
  oo <- fast.load(basename(save.path), path = dirname(save.path))
  set.seed(19)
  p <- predict(o, newdata = dta, seed = -19)
  set.seed(19)
  pp <- predict(oo, newdata = dta, seed = -19)
  print(summary(p$predicted - pp$predicted))
  print(all.equal(as.character(p$class), as.character(pp$class)))
  unlink(save.path, recursive = TRUE)
}

## ------------------------------------------------------------
## Survival
## ------------------------------------------------------------
if (requireNamespace("fst", quietly = TRUE) &&
    requireNamespace("data.table", quietly = TRUE)) {
  data(pbc, package = "randomForestSRC")
  o <- rfsrc(Surv(days, status) ~ ., data = pbc, ntree = 100)
  print(o)
  save.path <- tempfile("rfsrc-forest-")
  fast.save(o, path = save.path, testing = FALSE)
  oo <- fast.load(basename(save.path), path = dirname(save.path))
  set.seed(19)
  p <- predict(o, seed = -19)
  set.seed(19)
  pp <- predict(oo, seed = -19)
  print(summary(p$predicted - pp$predicted))
  print(summary(p$predicted.oob - pp$predicted.oob))
  unlink(save.path, recursive = TRUE)
}

## ------------------------------------------------------------
## Survival with save.memory = TRUE
## ------------------------------------------------------------
if (requireNamespace("fst", quietly = TRUE) &&
    requireNamespace("data.table", quietly = TRUE)) {
  data(pbc, package = "randomForestSRC")
  o <- rfsrc(Surv(days, status) ~ ., data = pbc,
             ntree = 100, save.memory = TRUE)
  print(o)
  save.path <- tempfile("rfsrc-forest-")
  fast.save(o, path = save.path, testing = FALSE)
  oo <- fast.load(basename(save.path), path = dirname(save.path))
  set.seed(19)
  p <- predict(o, seed = -19)
  set.seed(19)
  pp <- predict(oo, seed = -19)
  print(summary(p$predicted - pp$predicted))
  print(summary(p$predicted.oob - pp$predicted.oob))
  unlink(save.path, recursive = TRUE)
}

## ------------------------------------------------------------
## Competing risks
## ------------------------------------------------------------
if (requireNamespace("fst", quietly = TRUE) &&
    requireNamespace("data.table", quietly = TRUE)) {
  data(wihs, package = "randomForestSRC")
  o <- rfsrc(Surv(time, status) ~ ., data = wihs, nsplit = 3, ntree = 100)
  print(o)
  save.path <- tempfile("rfsrc-forest-")
  fast.save(o, path = save.path, testing = FALSE)
  oo <- fast.load(basename(save.path), path = dirname(save.path))
  set.seed(19)
  p <- predict(o, seed = -19)
  set.seed(19)
  pp <- predict(oo, seed = -19)
  print(summary(p$predicted - pp$predicted))
  print(summary(p$predicted.oob - pp$predicted.oob))
  unlink(save.path, recursive = TRUE)
}

## ------------------------------------------------------------
## Multivariate regression and classification
## ------------------------------------------------------------
if (requireNamespace("fst", quietly = TRUE) &&
    requireNamespace("data.table", quietly = TRUE)) {
  data(nutrigenomic, package = "randomForestSRC")
  ydta <- data.frame(diet = nutrigenomic$diet,
                     genotype = nutrigenomic$genotype,
                     nutrigenomic$lipids)
  o <- rfsrc(get.mv.formula(colnames(ydta)),
             data = data.frame(ydta, nutrigenomic$genes),
             ntree = 100, importance = TRUE, nsplit = 10)
  print(o)
  save.path <- tempfile("rfsrc-forest-")
  fast.save(o, path = save.path, testing = FALSE)
  oo <- fast.load(basename(save.path), path = dirname(save.path))
  set.seed(19)
  p <- predict(o, seed = -19)
  set.seed(19)
  pp <- predict(oo, seed = -19)
  print(summary(get.mv.predicted(p, oob = FALSE) -
                get.mv.predicted(pp, oob = FALSE)))
  print(summary(get.mv.predicted(p) - get.mv.predicted(pp)))
  for (yn in names(p$classOutput)) {
    print(yn)
    print(all.equal(as.character(p$classOutput[[yn]]$class),
                    as.character(pp$classOutput[[yn]]$class)))
    print(all.equal(as.character(p$classOutput[[yn]]$class.oob),
                    as.character(pp$classOutput[[yn]]$class.oob)))
  }
  unlink(save.path, recursive = TRUE)
}
}

\dontrun{

## ------------------------------------------------------------
## Classification: optional alzheimers data from varPro
## ------------------------------------------------------------
if (requireNamespace("fst", quietly = TRUE) &&
    requireNamespace("data.table", quietly = TRUE)) {
  data(alzheimers, package = "varPro")
  o <- rfsrc(Diagnosis ~ ., data = alzheimers)
  print(o)
  save.path <- tempfile("rfsrc-forest-")
  fast.save(o, path = save.path, testing = FALSE)
  oo <- fast.load(basename(save.path), path = dirname(save.path))
  set.seed(19)
  p <- predict(o, seed = -19)
  set.seed(19)
  pp <- predict(oo, seed = -19)
  print(summary(p$predicted - pp$predicted))
  print(summary(p$predicted.oob - pp$predicted.oob))
  print(all.equal(as.character(p$class), as.character(pp$class)))
  print(all.equal(as.character(p$class.oob),
                  as.character(pp$class.oob)))
  unlink(save.path, recursive = TRUE)
}


## ------------------------------------------------------------
## Optional memory-intensive anonymous survival test
## ------------------------------------------------------------
## This test repeats each PBC row 250 times and can require substantial memory.
if (requireNamespace("fst", quietly = TRUE) &&
    requireNamespace("data.table", quietly = TRUE)) {
  data(pbc, package = "randomForestSRC")
  dta <- pbc[rep(seq_len(nrow(pbc)), each = 250), ]
  o <- rfsrc.anonymous(Surv(days, status) ~ ., data = dta)
  print(o)
  save.path <- tempfile("rfsrc-forest-")
  fast.save(o, path = save.path, testing = FALSE)
  oo <- fast.load(basename(save.path), path = dirname(save.path))
  set.seed(19)
  p <- predict(o, newdata = dta, seed = -19)
  set.seed(19)
  pp <- predict(oo, newdata = dta, seed = -19)
  print(summary(p$predicted - pp$predicted))
  unlink(save.path, recursive = TRUE)
}
}
}
\keyword{file}
