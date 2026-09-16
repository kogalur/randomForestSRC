\name{fast.saveload}
\alias{fast.saveload}
\alias{fast.save}
\alias{fast.load}
\alias{fast.save.list}
\alias{fast.load.list}
\title{Fast Saving and Loading Objects}
\description{
  Save and load the forest component of a fitted random forest using a
  directory of files.  The loaded forest can be passed to
  \code{predict()}.  Corresponding helpers save and load lists of forests.
}
\usage{
fast.save(o, path = NULL, testing = TRUE, units = "Mb")

fast.load(directory, path = NULL, testing = FALSE, units = "Mb")

fast.save.list(o, path = NULL, testing = FALSE, units = "Mb")

fast.load.list(directory, path = NULL, testing = FALSE, units = "Mb")
}
\arguments{
  \item{o}{For \code{fast.save()}, a grow object with class beginning
    \code{c("rfsrc", "grow")} and a retained \code{forest} component.
    Supply the grow object, not \code{o$forest}.  For
    \code{fast.save.list()}, a nonempty list of such grow objects.}
  \item{path}{For saving, the destination directory, including its name.
    The default is \code{file.path(getwd(), "forest")}.  An existing
    destination and all its contents are removed before saving.
    For loading, the parent directory containing \code{directory};
    \code{NULL} uses the current working directory.  See Details for
    matching save and load paths.}
  \item{directory}{For loading, the name of the saved forest directory
    relative to \code{path}, or to the current working directory when
    \code{path = NULL}.  For \code{fast.load.list()}, the name of the
    directory containing the saved forest subdirectories.}
  \item{testing}{Print component-size diagnostics using \code{lsos()}.
    Defaults to \code{TRUE} for \code{fast.save()} and \code{FALSE}
    for the other helpers.  This is a size display, not a save/load
    comparison or a timing benchmark.}
  \item{units}{Character string specifying the units for the diagnostic
    size display, passed to \code{format(object.size(...), units = units)}.
    The default is \code{"Mb"}.  This does not affect the saved format.}
}
\details{
  \subsection{Saving and loading a forest}{
    \code{fast.save()} saves \code{o$forest}, including its retained
    training data and tree information.  \code{fast.load()} returns the
    reconstructed forest object, which can be used directly with
    \code{\link{predict.rfsrc}}.

    Components outside \code{o$forest}, such as the grow object's
    prediction and performance summaries, are not saved.  To preserve
    the complete R object instead, use \code{\link[base]{saveRDS}} and
    \code{\link[base]{readRDS}}.

    Saving requires the \pkg{fst} and \pkg{data.table} packages;
    loading requires \pkg{fst}.  Large tabular components are written
    with \code{fst::write_fst()}, and the remaining forest metadata are
    stored with \code{saveRDS(..., compress = FALSE)}.  The files use a
    \file{.rda} suffix despite these different formats.  Keep the entire
    directory together and read it using \code{fast.load()}.
  }
  \subsection{Directories}{
    \strong{Both save helpers remove the destination directory and all
    its contents before writing.}  Use a dedicated directory for the
    forest or forest list, rather than a directory containing other work.
    The destination is then created, including parent directories as
    needed.

    Saving takes the full destination in \code{path}.  Loading takes
    its name in \code{directory} and its parent in \code{path}.  For a
    destination stored in \code{save.path}, the matching load call is
    \code{fast.load(basename(save.path), path = dirname(save.path))}.
    With the default save destination, use \code{fast.load("forest")}.
  }
  \subsection{Lists of forests}{
    \code{fast.save.list()} saves the elements in subdirectories
    \file{forest1}, \file{forest2}, and so on, inside \code{path}.
    Original list names are not stored.

    \code{fast.load.list()} loads every entry returned by
    \code{list.files()} in the list directory and returns an unnamed
    list in that order.  Keep only saved forest subdirectories there.
    Directory-name ordering can differ from the original numeric index
    order: for example, \file{forest10} can precede \file{forest2}.
    Use \code{fast.load()} on a named subdirectory when selecting a
    particular saved element.
  }
  \subsection{Sampling-size function}{
    The saved sampling-size information is
    \code{sampfrac = o$forest$sampsize(1)}.  Loading reconstructs the
    function as \code{function(x) x * sampfrac}.  This preserves a
    proportional sampling-size rule; it does not preserve an arbitrary
    nonlinear custom function.  Use \code{saveRDS()} when that original
    function must be retained.
  }
}
\value{
  \describe{
    \item{\code{fast.save}}{Used for its file-writing side effect.
      Its return value is the memory-usage matrix from the final
      \code{gc(FALSE)} call, not the saved forest.}
    \item{\code{fast.load}}{The reconstructed forest component, retaining
      its saved class, ordinarily \code{c("rfsrc", "forest", family)}.
      It is a forest object rather than the original grow object.}
    \item{\code{fast.save.list}}{An unnamed list of the return values
      from the individual \code{fast.save()} calls.  The forests are
      saved in the destination subdirectories.}
    \item{\code{fast.load.list}}{An unnamed list of reconstructed forest
      objects in the directory-listing order described above.}
  }
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
    cat("check coherence for outcome:", yn, "\n")
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
