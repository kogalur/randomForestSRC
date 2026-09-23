## Helpers for the learned imputer and its OOD reconstruction scores.
##
## Raw schema: retained column names/types learned before imputation.
## Target: a raw column with a saved conditional learner.
## Auxiliary: supervised prediction columns, used as predictors only.
## Harmonized data: deployment data before initialization/model updates.
## Working data: initialized raw data, optionally augmented by auxiliaries.
##------------------------------------------------------------------
## Small utilities and persistence prerequisites
##------------------------------------------------------------------
## Use a default only for NULL; zero-length and missing-valued inputs are retained.
`%||%` <- function(x, y) if (is.null(x)) y else x
## Require the serialization dependency before any fast.save/fast.load operation.
.check.fst <- function() {
  if (!requireNamespace("fst", quietly = TRUE)) {
    stop("Package 'fst' is required for save/load support.", call. = FALSE)
  }
  invisible(TRUE)
}
## Format local clock time for progress messages.
.timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}
## Emit one timestamped progress line only when verbose is exactly TRUE.
.msg <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) {
    cat(sprintf("[%s] ", .timestamp()), ..., "\n", sep = "")
    flush.console()
  }
}
## Create missing parent directories; return the supplied path invisibly.
.safe.dir.create <- function(path) {
  if (!dir.exists(path) &&
      !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create directory: ", path, call. = FALSE)
  }
  invisible(path)
}
## Remove an existing directory tree; the caller decides whether wiping is allowed.
.safe.unlink.dir <- function(path) {
  if (file.exists(path) || dir.exists(path)) {
    status <- unlink(path, recursive = TRUE, force = TRUE, expand = FALSE)
    if (status != 0L || file.exists(path) || dir.exists(path)) {
      stop("Could not remove path: ", path, call. = FALSE)
    }
  }
  invisible(path)
}
##------------------------------------------------------------------
## Staged persistence
##------------------------------------------------------------------
## Resolve existing ancestors as well as relative/nonexistent destinations.
## This makes path containment checks meaningful before directory creation.
.impute.canonical.path <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(trimws(path))) {
    stop("A nonempty directory path is required.", call. = FALSE)
  }
  path <- path.expand(path)
  if (file.exists(path) || dir.exists(path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  parent <- dirname(path)
  if (identical(parent, path)) {
    stop("Cannot resolve directory path: ", path, call. = FALSE)
  }
  parent <- .impute.canonical.path(parent)
  leaf <- basename(path)
  if (leaf == ".") return(parent)
  if (leaf == "..") return(dirname(parent))
  file.path(parent, leaf)
}
.impute.path.within <- function(path, parent) {
  ## Conservative comparisons on platforms with commonly case-insensitive
  ## filesystems; false overlap is safer than deleting an aliased source.
  if (.Platform$OS.type == "windows" || identical(Sys.info()[["sysname"]], "Darwin")) {
    path <- tolower(path)
    parent <- tolower(parent)
  }
  identical(path, parent) || startsWith(path, paste0(sub("/+$", "", parent), "/"))
}
.impute.path.part <- function(x, arg, nested = FALSE) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !nzchar(trimws(x)) || grepl("\\", x, fixed = TRUE) ||
      grepl(":", x, fixed = TRUE) || startsWith(x, "~") ||
      (!nested && grepl("/", x, fixed = TRUE))) {
    stop("'", arg, "' must be a relative ",
         if (nested) "directory path" else "file-name component", ".", call. = FALSE)
  }
  parts <- strsplit(paste0(x, "/."), "/", fixed = TRUE)[[1L]]
  parts <- parts[-length(parts)]
  if (any(!nzchar(parts) | parts %in% c(".", ".."))) {
    stop("'", arg, "' cannot contain empty, '.' or '..' path components.",
         call. = FALSE)
  }
  invisible(x)
}
.impute.check.layout <- function(manifest) {
  .impute.path.part(manifest$learner.root %||% "learners", "learner.root", nested = TRUE)
  .impute.path.part(manifest$learner.prefix %||% "impute.learner.", "learner.prefix")
  saved.names <- character(0)
  for (target in manifest$targets) {
    info <- manifest$learners[[target]]
    if (!identical(info$status, "ok")) next
    .impute.path.part(info$learner.name, "learner.name")
    saved.names <- c(saved.names, info$learner.name)
  }
  if (.is.supervised.manifest(manifest)) {
    .impute.path.part(manifest$supervised$model.name, "supervised model.name")
    saved.names <- c(saved.names, manifest$supervised$model.name)
  }
  if (anyDuplicated(saved.names)) {
    stop("Saved learner names must be unique.", call. = FALSE)
  }
  invisible(manifest)
}
.impute.save.destination <- function(path, source.path = NULL) {
  dest <- .impute.canonical.path(path)
  protected <- unique(c(.impute.canonical.path(getwd()),
                        .impute.canonical.path(path.expand("~"))))
  if (identical(dirname(dest), dest) ||
      any(vapply(protected, .impute.path.within, logical(1), parent = dest))) {
    stop("Use a dedicated imputer directory, not a root, home, working directory ",
         "or an ancestor of those directories: ", dest, call. = FALSE)
  }
  if (file.exists(dest) && !dir.exists(dest)) {
    stop("Save destination is not a directory: ", dest, call. = FALSE)
  }
  if (!is.null(source.path)) {
    source.path <- .impute.canonical.path(source.path)
    same <- .impute.path.within(dest, source.path) &&
      .impute.path.within(source.path, dest)
    if (!same && (.impute.path.within(dest, source.path) ||
                  .impute.path.within(source.path, dest))) {
      stop("Source and destination imputer directories must not contain one another.",
           call. = FALSE)
    }
  }
  dest
}
.impute.begin.save <- function(path, wipe = TRUE, source.path = NULL,
                                targets = NULL) {
  .impute.flag(wipe, "wipe")
  dest <- .impute.save.destination(path, source.path)
  if (!is.null(source.path)) {
    source.path <- .impute.canonical.path(source.path)
    same <- .impute.path.within(dest, source.path) &&
      .impute.path.within(source.path, dest)
    previous.manifest <- file.path(dest, "manifest.rds")
    if (same && file.exists(previous.manifest)) {
      previous <- readRDS(previous.manifest)
      if (!setequal(previous$targets, targets)) {
        stop("A target-subset object cannot replace its source imputer in place. ",
             "Save the subset to a different directory.", call. = FALSE)
      }
    }
  }
  ## Stage beside the destination so directory renames stay on one filesystem.
  .safe.dir.create(dirname(dest))
  stage <- tempfile(pattern = paste0(".", basename(dest), ".stage-"),
                    tmpdir = dirname(dest))
  .safe.dir.create(stage)
  ready <- FALSE
  on.exit(if (!ready) try(.safe.unlink.dir(stage), silent = TRUE), add = TRUE)
  ## wipe = FALSE preserves unrelated existing files without modifying them
  ## while new learners are written. Copies are private to the staging bundle.
  if (!isTRUE(wipe) && dir.exists(dest)) {
    entries <- list.files(dest, all.files = TRUE, no.. = TRUE, full.names = TRUE)
    if (length(entries) && !all(file.copy(entries, stage, recursive = TRUE,
                                          copy.mode = TRUE, copy.date = TRUE))) {
      stop("Could not stage the existing destination directory.", call. = FALSE)
    }
  }
  .safe.unlink.dir(file.path(stage, "manifest.rds"))
  ready <- TRUE
  list(path = dest, stage = stage)
}
.impute.write.model <- function(model, name, root, label) {
  .impute.path.part(name, "learner.name")
  path <- file.path(root, name)
  .safe.unlink.dir(path)
  status <- fast.save(model, path, testing = FALSE)
  if (identical(status, FALSE)) stop("Failed to save ", label, ".", call. = FALSE)
  ## Verify each serialized learner while staging, without retaining the whole
  ## bank in memory. A failed read leaves the existing destination untouched.
  restored <- .fast.load.named(name, root, label = label, strict = TRUE)
  if (is.null(restored) || !is.list(restored)) {
    stop("Saved ", label, " could not be read back as a model.", call. = FALSE)
  }
  invisible(NULL)
}
.impute.commit.save <- function(transaction, manifest) {
  .impute.check.layout(manifest)
  saveRDS(manifest, file.path(transaction$stage, "manifest.rds"))
  if (!is.list(readRDS(file.path(transaction$stage, "manifest.rds")))) {
    stop("Staged manifest could not be read back.", call. = FALSE)
  }
  dest <- transaction$path
  backup <- NULL
  installed <- FALSE
  on.exit({
    if (!installed && !is.null(backup) && dir.exists(backup)) {
      restored <- !file.exists(dest) && tryCatch(file.rename(backup, dest),
                                                  error = function(e) FALSE)
      if (!isTRUE(restored)) {
        warning("The previous imputer is preserved at: ", backup,
                call. = FALSE)
      }
    }
  }, add = TRUE)
  if (dir.exists(dest)) {
    backup <- tempfile(pattern = paste0(".", basename(dest), ".backup-"),
                       tmpdir = dirname(dest))
    if (!file.rename(dest, backup)) {
      stop("Could not move the previous imputer to a backup; it is unchanged.",
           call. = FALSE)
    }
  }
  if (!file.rename(transaction$stage, dest)) {
    stop("Could not install the staged imputer; restoring the previous destination.",
         call. = FALSE)
  }
  installed <- TRUE
  if (!is.null(backup)) {
    tryCatch(.safe.unlink.dir(backup), error = function(e) {
      warning("Imputer saved, but its backup could not be removed: ", backup,
              call. = FALSE)
    })
  }
  invisible(dest)
}
## Make a filename component from a target label; indices distinguish collisions.
.safe.name <- function(x) {
  x <- gsub("[^[:alnum:].]+", ".", x)
  x <- gsub("\\.+", ".", x)
  x <- gsub("^\\.+|\\.+$", "", x)
  if (!nzchar(x)) x <- "var"
  x
}
##------------------------------------------------------------------
## Raw column normalization, schema and initialization
##------------------------------------------------------------------
## Identify integer/double columns, excluding the date/time classes handled as factors.
.is.real.valued <- function(x) {
  is.null(dim(x)) && (is.double(x) || is.integer(x)) &&
    !inherits(x, c("Date", "POSIXt", "difftime"))
}
## Retain factors and real-valued inputs; otherwise attempt factor conversion.
.coerce.supported.column <- function(x, nm) {
  if (is.factor(x) || .is.real.valued(x)) {
    return(x)
  }
  out <- tryCatch(
    factor(x),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    stop("Column `", nm, "` is neither real-valued nor factor and could not be coerced to factor: ",
         conditionMessage(out),
         call. = FALSE)
  }
  out
}
## Normalize the container and column types before constructing the missingness mask.
.normalize.training.data <- function(data) {
  data <- .impute.data.frame(data)
  for (nm in names(data)) {
    data[[nm]] <- .coerce.supported.column(data[[nm]], nm)
  }
  data
}
## Return the retained raw table plus original-row and original-column keep masks.
.drop.all.na.train <- function(data) {
  which.na <- is.na(data)
  keep.rows <- rowSums(which.na) < ncol(data)
  keep.cols <- colSums(which.na) < nrow(data)
  list(
    data = data[keep.rows, keep.cols, drop = FALSE],
    keep.rows = keep.rows,
    keep.cols = keep.cols
  )
}
## Record factor levels/order and pre-imputation numeric storage/support by column.
.build.schema <- function(data) {
  out <- lapply(names(data), function(nm) {
    x <- data[[nm]]
    integer.info <- .integer.support.info(x)
    list(
      class = class(x),
      is.factor = is.factor(x),
      ordered = is.ordered(x),
      levels = if (is.factor(x)) levels(x) else NULL,
      is.integer = is.integer(x),
      is.numeric = .is.real.valued(x),
      integer.support = isTRUE(integer.info$restore),
      integer.storage = isTRUE(integer.info$storage.integer),
      integer.frac = integer.info$integer.frac
    )
  })
  names(out) <- names(data)
  out
}
## Supply integer metadata defaults for older stored schemas without refitting.
.normalize.schema.integer.support <- function(schema) {
  if (is.null(schema) || length(schema) == 0L) {
    return(schema)
  }
  for (nm in names(schema)) {
    sc <- schema[[nm]]
    if (!is.list(sc)) next
    if (is.null(sc$integer.storage)) {
      sc$integer.storage <- isTRUE(sc$is.integer)
    }
    if (is.null(sc$integer.support)) {
      sc$integer.support <- isTRUE(sc$is.integer)
    }
    if (is.null(sc$integer.frac)) {
      sc$integer.frac <- if (isTRUE(sc$is.integer)) 1 else NA_real_
    }
    schema[[nm]] <- sc
  }
  schema
}
## Test whether the saved schema requests integer-grid restoration.
.schema.restores.integer <- function(sc) {
  is.list(sc) && (isTRUE(sc$integer.support) || isTRUE(sc$is.integer))
}
## Translate a schema entry to the shared integer-restoration helper format.
.schema.integer.info <- function(sc) {
  list(
    restore = .schema.restores.integer(sc),
    storage.integer = isTRUE(sc$is.integer) || isTRUE(sc$integer.storage),
    integer.frac = sc$integer.frac %||% (if (isTRUE(sc$is.integer)) 1 else NA_real_)
  )
}
## Convert factor labels, rather than their internal codes, to numeric values.
## Unparseable labels become NA with conversion warnings suppressed.
.as.numeric.safe <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  suppressWarnings(as.numeric(x))
}
## Numeric conversion entry point; schema is currently accepted but not consulted.
.as.numeric.from.schema <- function(x, schema = NULL) {
  .as.numeric.safe(x)
}
## Return a most frequent observed label, or NA when no observed values remain.
.mode.value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA)
  tb <- sort(table(x), decreasing = TRUE)
  names(tb)[1]
}
## Use the mean of finite observed values; return NA_real_ if none are available.
.numeric.init <- function(x) {
  x <- .as.numeric.safe(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}
## Compute one factor mode or numeric mean per supplied table column.
.compute.init <- function(data, schema) {
  out <- lapply(names(data), function(nm) {
    x <- data[[nm]]
    sc <- schema[[nm]]
    if (isTRUE(sc$is.factor)) {
      .mode.value(as.character(x))
    }
    else {
      .numeric.init(.as.numeric.from.schema(x, sc))
    }
  })
  names(out) <- names(data)
  out
}
## Compute numeric SDs for pass-change scaling; use one for factors or zero/invalid SDs.
.compute.scale <- function(data, schema) {
  out <- lapply(names(data), function(nm) {
    sc <- schema[[nm]]
    if (isTRUE(sc$is.factor)) {
      1
    }
    else {
      x <- .as.numeric.from.schema(data[[nm]], sc)
      s <- stats::sd(x, na.rm = TRUE)
      if (!is.finite(s) || s <= 0) s <- 1
      s
    }
  })
  names(out) <- names(data)
  out
}
##------------------------------------------------------------------
## Target selection and raw predictor maps
##------------------------------------------------------------------
## Select targets from the original missingness mask or select every retained column.
.resolve.targets <- function(which.na, target.mode = c("missing.only", "all")) {
  target.mode <- match.arg(target.mode)
  miss.count <- colSums(which.na)
  if (target.mode == "missing.only") {
    names(miss.count)[miss.count > 0]
  }
  else {
    colnames(which.na)
  }
}
## Build named raw-predictor sets, always excluding each target itself.
## NULL/list omissions use all non-target columns; character(0) selects none.
.resolve.predictor.map <- function(targets, all.names, deployment.xvars = NULL) {
  if (is.null(deployment.xvars)) {
    out <- lapply(targets, function(y) setdiff(all.names, y))
    names(out) <- targets
    return(out)
  }
  if (is.character(deployment.xvars)) {
    unknown <- setdiff(deployment.xvars, all.names)
    if (length(unknown) > 0L) {
      warning("Ignoring deployment.xvars predictors not found in training data: ",
              paste(unknown, collapse = ", "),
              call. = FALSE)
    }
    xvars <- intersect(all.names, deployment.xvars)
    out <- lapply(targets, function(y) setdiff(xvars, y))
    names(out) <- targets
    return(out)
  }
  if (is.list(deployment.xvars)) {
    nms <- names(deployment.xvars)
    if (is.null(nms) || length(nms) != length(deployment.xvars) ||
        any(is.na(nms) | !nzchar(nms))) {
      stop("'deployment.xvars' supplied as a list must be named by target variable.",
           call. = FALSE)
    }
    dup <- unique(nms[duplicated(nms)])
    if (length(dup) > 0L) {
      stop("'deployment.xvars' contains duplicated target names: ",
           paste(dup, collapse = ", "),
           call. = FALSE)
    }
    extra.targets <- setdiff(nms, targets)
    if (length(extra.targets) > 0L) {
      warning("Ignoring deployment.xvars entries for non-target variables: ",
              paste(extra.targets, collapse = ", "),
              call. = FALSE)
    }
    bad.entries <- names(deployment.xvars)[vapply(deployment.xvars, function(x) {
      !(is.null(x) || is.character(x))
    }, logical(1))]
    if (length(bad.entries) > 0L) {
      stop("Each 'deployment.xvars' list entry must be NULL or a character vector. ",
           "Problem entries: ", paste(bad.entries, collapse = ", "),
           call. = FALSE)
    }
    listed.predictors <- unique(unlist(deployment.xvars[intersect(nms, targets)],
                                       use.names = FALSE))
    unknown <- setdiff(listed.predictors, all.names)
    if (length(unknown) > 0L) {
      warning("Ignoring deployment.xvars predictors not found in training data: ",
              paste(unknown, collapse = ", "),
              call. = FALSE)
    }
    out <- lapply(targets, function(y) {
      xvars <- deployment.xvars[[y]]
      if (is.null(xvars)) {
        xvars <- setdiff(all.names, y)
      }
      setdiff(intersect(all.names, xvars), y)
    })
    names(out) <- targets
    return(out)
  }
  stop("'deployment.xvars' must be NULL, a character vector, or a named list.",
       call. = FALSE)
}
##------------------------------------------------------------------
## Manifest compatibility and supervised forest support
##------------------------------------------------------------------
## Supply legacy defaults for supervised, raw-predictor and integer-schema metadata.
## This normalizes known fields; it is not a full manifest integrity/version validator.
.normalize.impute.learn.manifest <- function(manifest) {
  if (is.null(manifest$spec.version)) {
    manifest$spec.version <- 1L
  }
  if (is.null(manifest$supervised) || !is.list(manifest$supervised)) {
    manifest$supervised <- list(enabled = FALSE)
  }
  if (is.null(manifest$supervised$enabled)) {
    manifest$supervised$enabled <- FALSE
  }
  manifest$schema <- .normalize.schema.integer.support(manifest$schema)
  if (is.null(manifest$supervised$response.schema)) {
    manifest$supervised$response.schema <- NULL
  }
  else {
    manifest$supervised$response.schema <- .normalize.schema.integer.support(
      manifest$supervised$response.schema
    )
  }
  if (is.null(manifest$predictor.map.raw)) {
    manifest$predictor.map.raw <- manifest$predictor.map
  }
  manifest
}
## Identify manifests containing an enabled supervised auxiliary forest.
.is.supervised.manifest <- function(manifest) {
  is.list(manifest$supervised) && isTRUE(manifest$supervised$enabled)
}
## Resolve a package helper from the current environment or randomForestSRC namespace.
.get.rfsrc.internal <- function(name) {
  fn <- get0(name, mode = "function")
  if (!is.null(fn)) {
    return(fn)
  }
  ns <- tryCatch(
    asNamespace("randomForestSRC"),
    error = function(e) NULL
  )
  if (is.null(ns) || !exists(name, envir = ns, inherits = FALSE)) {
    return(NULL)
  }
  get(name, envir = ns, inherits = FALSE)
}
## Create a terms object while identifying the offending formula argument on errors.
.parse.formula.terms <- function(formula, data, arg.name = "formula") {
  tt <- tryCatch(
    stats::terms(formula, data = data, keep.order = TRUE),
    error = function(e) e
  )
  if (inherits(tt, "error")) {
    stop("Failed to parse '", arg.name, "': ", conditionMessage(tt),
         call. = FALSE)
  }
  tt
}
## Retain the user formula environment and make supported response constructors available.
.make.supervised.formula.env <- function(parent = environment()) {
  if (is.null(parent)) {
    parent <- baseenv()
  }
  env <- new.env(parent = parent)
  surv.fun <- get0("Surv", envir = parent, mode = "function", inherits = TRUE)
  if (is.null(surv.fun)) {
    surv.fun <- tryCatch(
      get("Surv", envir = asNamespace("survival"), inherits = FALSE),
      error = function(e) NULL
    )
  }
  if (!is.null(surv.fun)) {
    env$Surv <- surv.fun
  }
  multivar.fun <- .get.rfsrc.internal("Multivar")
  if (!is.null(multivar.fun)) {
    env$Multivar <- multivar.fun
  }
  env$cbind <- base::cbind
  env
}
## Restrict the supervised RHS to raw columns retained after all-NA preprocessing.
.raw.supervised.predictors <- function(tt, data, arg.name) {
  labels <- attr(tt, "term.labels") %||% character(0)
  if (length(labels) == 0L) {
    stop("'", arg.name, "' must include at least one predictor.", call. = FALSE)
  }
  vars <- lapply(labels, function(label) {
    expr <- tryCatch(parse(text = label, keep.source = FALSE),
                     error = function(e) NULL)
    if (length(expr) == 1L && is.symbol(expr[[1L]])) {
      as.character(expr[[1L]])
    } else NA_character_
  })
  vars <- unlist(vars, use.names = FALSE)
  bad <- is.na(vars) | !(vars %in% names(data))
  if (any(bad) || length(attr(tt, "offset")) > 0L) {
    stop("'", arg.name, "' must reference raw predictor columns in 'data'. ",
         "Derived terms, interactions, and transformations are not supported. ",
         "Problem terms: ", paste(c(labels[bad],
           if (length(attr(tt, "offset"))) "offset(...)"), collapse = ", "),
         call. = FALSE)
  }
  vars
}
.project.supervised.formula <- function(supervised.formula, data, train.names,
                                        arg.name = "supervised.formula") {
  tt <- .parse.formula.terms(supervised.formula, data = data, arg.name = arg.name)
  xvar.names <- .raw.supervised.predictors(tt, data, arg.name)
  xvar.names <- intersect(xvar.names, train.names)
  if (length(xvar.names) == 0L) {
    stop("No retained predictor columns remain for '", arg.name,
         "' after training preprocessing.", call. = FALSE)
  }
  env <- .make.supervised.formula.env(environment(supervised.formula))
  out <- .impute.formula(supervised.formula[[2L]], xvar.names, env = env)
  if (identical(attr(tt, "intercept"), 0L)) {
    out[[3L]] <- call("-", out[[3L]], 1)
  }
  out
}
## Resolve supervised response and raw predictor columns before raw-table imputation.
## The raw RHS is checked separately from the package-specific response parser.
.parse.supervised.spec <- function(supervised.formula, data) {
  if (is.null(supervised.formula)) {
    return(list(enabled = FALSE))
  }
  if (!inherits(supervised.formula, "formula")) {
    stop("'supervised.formula' must be a formula.", call. = FALSE)
  }
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  tt <- .parse.formula.terms(supervised.formula, data = data,
                             arg.name = "supervised.formula")
  response.pos <- attr(tt, "response") %||% 0L
  if (response.pos != 1L) {
    stop("'supervised.formula' must include a response on the left-hand side.",
         call. = FALSE)
  }
  xvar.names <- .raw.supervised.predictors(tt, data, "supervised.formula")
  parser <- .get.rfsrc.internal("parseFormula")
  parsed <- NULL
  if (!is.null(parser)) {
    parsed <- tryCatch(
      parser(supervised.formula, data = data),
      error = function(e) e
    )
    if (inherits(parsed, "error")) {
      parsed <- NULL
    }
  }
  if (is.null(parsed)) {
    yvar.names <- all.vars(supervised.formula[[2L]])
    subj.names <- character(0)
    family.guess <- NA_character_
  }
  else {
    yvar.names <- parsed$yvar.names %||% character(0)
    subj.names <- parsed$subj.names %||% character(0)
    family.guess <- parsed$family %||% NA_character_
  }
  response.names <- unique(c(subj.names, yvar.names))
  if (length(response.names) == 0L) {
    stop("Unable to resolve supervised response variables from 'supervised.formula'.",
         call. = FALSE)
  }
  if (!all(response.names %in% names(data))) {
    bad <- setdiff(response.names, names(data))
    stop("The supervised response variables were not found in 'data': ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  response.predictor.overlap <- intersect(response.names, xvar.names)
  if (length(response.predictor.overlap) > 0L) {
    stop("'supervised.formula' must not include supervised response variables ",
         "on the right-hand side. Problem variables: ",
         paste(response.predictor.overlap, collapse = ", "),
         call. = FALSE)
  }
  y.data <- data[, response.names, drop = FALSE]
  for (nm in names(y.data)) {
    y.data[[nm]] <- .coerce.supported.column(y.data[[nm]], nm)
  }
  list(
    enabled = TRUE,
    formula = supervised.formula,
    family.guess = family.guess,
    xvar.names = xvar.names,
    yvar.names = yvar.names,
    subj.names = subj.names,
    response.names = response.names,
    x.data = data[, xvar.names, drop = FALSE],
    y.data = y.data
  )
}
## Identify retained rows with complete supervised response information.
.supervised.usable.rows <- function(ydata) {
  out <- tryCatch(
    stats::complete.cases(ydata),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    ydf <- tryCatch(
      as.data.frame(ydata, stringsAsFactors = FALSE),
      error = function(e) e
    )
    if (inherits(ydf, "error")) {
      stop("Failed to determine which supervised responses are observed: ",
           conditionMessage(out), call. = FALSE)
    }
    out <- stats::complete.cases(ydf)
  }
  out
}
## Validate named overrides and remove arguments controlled by the supervised wrapper.
.parse.supervised.args <- function(supervised.args = list()) {
  supervised.args <- .impute.named.options(supervised.args, "supervised.args")
  blocked <- intersect(names(supervised.args), c("formula", "data", "forest"))
  if (length(blocked) > 0L) {
    warning("Ignoring supervised.args entries controlled internally: ",
            paste(blocked, collapse = ", "),
            call. = FALSE)
    supervised.args[blocked] <- NULL
  }
  supervised.args
}
## Recover family and variable names from the fit, forest or parsed formula metadata.
.resolve.supervised.fit.meta <- function(supervised.fit, fallback = NULL) {
  if (is.null(supervised.fit) || !is.list(supervised.fit)) {
    stop("The supervised forest fit did not return a valid randomForestSRC object.",
         call. = FALSE)
  }
  forest <- .resolve.predict.forest(supervised.fit)
  family <- supervised.fit$family %||%
    (if (!is.null(forest)) forest$family else NULL) %||%
    (if (!is.null(fallback)) fallback$family.guess else NULL) %||%
    NA_character_
  xvar.names <- supervised.fit$xvar.names %||%
    (if (!is.null(forest)) forest$xvar.names else NULL) %||%
    (if (!is.null(fallback)) fallback$xvar.names else NULL)
  yvar.names <- supervised.fit$yvar.names %||%
    (if (!is.null(forest)) forest$yvar.names else NULL) %||%
    (if (!is.null(fallback)) fallback$yvar.names else NULL)
  if (is.null(xvar.names) || length(xvar.names) == 0L) {
    stop("The supervised forest fit did not return 'xvar.names'.",
         call. = FALSE)
  }
  if (is.null(yvar.names) || length(yvar.names) == 0L) {
    stop("The supervised forest fit did not return 'yvar.names'.",
         call. = FALSE)
  }
  list(
    family = family,
    xvar.names = xvar.names,
    yvar.names = yvar.names
  )
}
## Infer prediction row count for constructing response metadata when yvar is absent.
.infer.supervised.pred.n <- function(x) {
  n <- x$n %||% NULL
  if (!is.null(n) && length(n) == 1L && is.finite(n)) {
    return(as.integer(n))
  }
  if (!is.null(x$xvar)) {
    xdf <- tryCatch(as.data.frame(x$xvar, check.names = FALSE,
                                  stringsAsFactors = FALSE),
                    error = function(e) NULL)
    if (!is.null(xdf)) {
      return(nrow(xdf))
    }
  }
  pred <- x$predicted %||% x$predicted.oob %||% NULL
  if (!is.null(pred)) {
    return(NROW(pred))
  }
  for (slot in c("classOutput", "regrOutput")) {
    obj <- x[[slot]]
    if (!is.null(obj) && length(obj) > 0L) {
      for (nm in names(obj)) {
        comp <- obj[[nm]]
        if (!is.null(comp$predicted)) {
          return(NROW(comp$predicted))
        }
        if (!is.null(comp$predicted.oob)) {
          return(NROW(comp$predicted.oob))
        }
      }
    }
  }
  0L
}
## Supply multivariate response metadata needed by get.mv.predicted().
## Missing response columns are placeholders, not deployment-time observed outcomes.
.patch.supervised.prediction.object <- function(x, family,
                                                response.schema = NULL,
                                                yvar.names = NULL) {
  if (!(family %in% c("regr+", "class+", "mix+"))) {
    return(x)
  }
  n <- .infer.supervised.pred.n(x)
  if (n <= 0L) {
    return(x)
  }
  y <- x$yvar %||% NULL
  if (is.null(y)) {
    if (is.null(response.schema) || length(response.schema) == 0L) {
      return(x)
    }
    y <- .na.data.from.schema(response.schema, n)
  }
  else {
    y <- tryCatch(
      as.data.frame(y, check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(y)) {
      if (is.null(response.schema) || length(response.schema) == 0L) {
        return(x)
      }
      y <- .na.data.from.schema(response.schema, n)
    }
    else if (!is.null(response.schema) && length(response.schema) > 0L) {
      wanted <- names(response.schema)
      extra.cols <- setdiff(names(y), wanted)
      miss.cols <- setdiff(wanted, names(y))
      if (length(extra.cols) > 0L) {
        y <- y[, setdiff(names(y), extra.cols), drop = FALSE]
      }
      if (length(miss.cols) > 0L) {
        y.miss <- .na.data.from.schema(response.schema[miss.cols], nrow(y))
        y <- data.frame(y, y.miss, check.names = FALSE, stringsAsFactors = FALSE)
      }
      y <- y[, wanted, drop = FALSE]
      y <- .restore.schema(y, response.schema, restore.integer = TRUE)
    }
  }
  if (!is.null(yvar.names) && length(yvar.names) > 0L) {
    have <- intersect(yvar.names, names(y))
    if (length(have) == length(yvar.names)) {
      y <- y[, yvar.names, drop = FALSE]
    }
  }
  x$yvar <- y
  if (!is.null(yvar.names) && length(yvar.names) > 0L) {
    x$yvar.names <- yvar.names
  }
  x
}
## Extract family-specific predictions into a numeric auxiliary block.
## Training requests OOB slots; deployment requests ordinary prediction slots.
.supervised.materialize.aux <- function(x, family, oob = FALSE,
                                        response.schema = NULL,
                                        yvar.names = NULL) {
  x <- .patch.supervised.prediction.object(
    x,
    family = family,
    response.schema = response.schema,
    yvar.names = yvar.names
  )
  if (isTRUE(oob)) {
    if (family %in% c("regr+", "class+", "mix+")) {
      getter <- .get.rfsrc.internal("get.mv.predicted")
      if (is.null(getter)) {
        stop("get.mv.predicted() is required to materialize supervised multivariate predictions.",
             call. = FALSE)
      }
      pred <- getter(x, oob = TRUE)
    }
    else {
      pred <- x$predicted.oob
    }
  }
  else {
    if (family %in% c("regr+", "class+", "mix+")) {
      getter <- .get.rfsrc.internal("get.mv.predicted")
      if (is.null(getter)) {
        stop("get.mv.predicted() is required to materialize supervised multivariate predictions.",
             call. = FALSE)
      }
      pred <- getter(x)
    }
    else {
      pred <- x$predicted
    }
  }
  out <- tryCatch(
    data.frame(predicted = pred, check.names = FALSE),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    stop("Failed to materialize supervised predicted values: ",
         conditionMessage(out), call. = FALSE)
  }
  out
}
## Convert auxiliary columns to numeric and impose stored names by column position.
.align.supervised.aux <- function(aux, expected.names = NULL) {
  aux <- as.data.frame(aux, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(aux) == 0L) {
    stop("The supervised prediction block has no columns.", call. = FALSE)
  }
  for (nm in names(aux)) {
    aux[[nm]] <- .as.numeric.safe(aux[[nm]])
  }
  if (is.null(expected.names)) {
    return(aux)
  }
  if (length(expected.names) != ncol(aux)) {
    stop("The supervised prediction block produced ", ncol(aux),
         " column(s), but ", length(expected.names),
         " were learned during training.", call. = FALSE)
  }
  names(aux) <- expected.names
  aux[, expected.names, drop = FALSE]
}
## Replace nonfinite auxiliary values with the corresponding stored initialization value.
.fill.supervised.aux <- function(aux, init) {
  aux <- as.data.frame(aux, check.names = FALSE, stringsAsFactors = FALSE)
  if (length(aux) == 0L) {
    return(aux)
  }
  for (nm in names(aux)) {
    x <- .as.numeric.safe(aux[[nm]])
    bad <- !is.finite(x)
    if (any(bad)) {
      x[bad] <- as.numeric(init[[nm]] %||% NA_real_)
    }
    aux[[nm]] <- x
  }
  aux
}
## Append all auxiliary predictors to each target's selected raw predictors.
.augment.predictor.map <- function(targets, predictor.map.raw, aux.names = NULL) {
  predictor.map.raw <- predictor.map.raw[targets]
  if (is.null(aux.names) || length(aux.names) == 0L) {
    return(predictor.map.raw)
  }
  out <- lapply(targets, function(target) {
    unique(c(predictor.map.raw[[target]], aux.names))
  })
  names(out) <- targets
  out
}
##------------------------------------------------------------------
## Serialized objects and learner-specific predictor conformity
##------------------------------------------------------------------
## Load a named serialized object; optionally return an error message instead of stopping.
.fast.load.named <- function(name, root, label = "saved object", strict = TRUE) {
  .check.fst()
  out <- tryCatch(
    fast.load(name, root),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    msg <- paste0("Failed to load ", label, " from ",
                  file.path(root, name), ": ",
                  conditionMessage(out))
    if (isTRUE(strict)) {
      stop(msg, call. = FALSE)
    }
    return(list(model = NULL, error = msg))
  }
  if (isTRUE(strict)) {
    return(out)
  }
  list(model = out, error = NULL)
}
## Accept either a grow object containing a forest or a forest object directly.
.resolve.predict.forest <- function(model) {
  if (is.null(model)) return(NULL)
  if (is.list(model) && !is.null(model$forest) &&
      is.list(model$forest) && !is.null(model$forest$xvar.names)) {
    return(model$forest)
  }
  if (is.list(model) && !is.null(model$xvar.names)) {
    return(model)
  }
  NULL
}
## Match prediction columns and factor encodings to the particular target forest.
## The bank-wide schema alone is insufficient when a learner saw only some rows.
.conform.x.to.forest <- function(x, model, ignore.levels = TRUE) {
  x <- data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  forest <- .resolve.predict.forest(model)
  if (is.null(forest) || is.null(forest$xvar.names)) {
    return(x)
  }
  wanted <- forest$xvar.names
  ## add/drop/reorder to match the learner exactly
  miss.cols <- setdiff(wanted, names(x))
  extra.cols <- setdiff(names(x), wanted)
  if (length(extra.cols) > 0) {
    x <- x[, setdiff(names(x), extra.cols), drop = FALSE]
  }
  if (length(miss.cols) > 0) {
    for (nm in miss.cols) x[[nm]] <- rep(NA, nrow(x))
  }
  x <- x[, wanted, drop = FALSE]
  ## restore unordered and ordered factors using forest metadata
  if (!is.null(forest$xvar.factor)) {
    x <- check.factor(x, forest$xvar.factor, ignore = ignore.levels)
    gtypes <- forest$xvar.factor$generic.types
    if (!is.null(gtypes) && length(gtypes) == ncol(x)) {
      real.idx <- which(gtypes == "R")
      if (length(real.idx) > 0) {
        x[real.idx] <- lapply(x[real.idx], .as.numeric.safe)
      }
    }
  }
  x
}
##------------------------------------------------------------------
## Bank fitting options and internal names
##------------------------------------------------------------------
## Separate bank-level forest options from their structural defaults.
## These settings are independent of the initial imputation forest settings.
.parse.full.sweep.options <- function(full.sweep.options = NULL) {
  fs <- .impute.named.options(full.sweep.options, "full.sweep.options")
  .impute.forest.options(fs)
  allowed <- c(
    "mtry", "splitrule", "bootstrap", "sampsize", "samptype",
    "perf.type", "rfq", "save.memory", "importance", "proximity",
    "terminal.qualts", "terminal.quants"
  )
  structural <- c("ntree", "nodesize", "nsplit")
  unknown <- setdiff(names(fs), c(structural, allowed))
  if (length(unknown) > 0L) {
    warning("Ignoring unknown full.sweep.options entries: ",
            paste(unknown, collapse = ", "),
            call. = FALSE)
  }
  list(
    ntree = fs$ntree %||% 100L,
    nodesize = fs$nodesize,
    nsplit = fs$nsplit %||% 10L,
    dots = fs[names(fs) %in% allowed]
  )
}
## Combine an index and sanitized target name into a stable saved-model name.
.make.learner.name <- function(i, target, prefix = "impute.learner.") {
  sprintf("%s%03d.%s", prefix, i, .safe.name(target))
}
## Choose a temporary response name that does not collide with predictor columns.
.make.response.name <- function(existing.names = character()) {
  .impute.response.name(existing.names)
}
## Rebuild factors in stored training order; unknown labels become NA.
.coerce.factor.levels <- function(x, levels, ordered = FALSE) {
  x <- as.character(x)
  x[!(x %in% levels)] <- NA_character_
  factor(x, levels = levels, ordered = ordered)
}
##------------------------------------------------------------------
## Deployment schema, initialization and restoration
##------------------------------------------------------------------
## Match raw deployment columns and types to the manifest before any imputation.
## Keep unseen-label masks separately so OOD scoring can distinguish novelty
## from ordinary NA values. Supplied supervised response columns are discarded.
.harmonize.newdata <- function(newdata, manifest, verbose = TRUE) {
  newdata <- .impute.data.frame(newdata, "newdata", check.finite = FALSE)
  required.cols <- manifest$columns
  added.cols <- setdiff(required.cols, names(newdata))
  extra.cols <- setdiff(names(newdata), required.cols)
  response.cols <- character(0)
  if (.is.supervised.manifest(manifest) &&
      !is.null(manifest$supervised$response.schema)) {
    response.cols <- intersect(names(manifest$supervised$response.schema), extra.cols)
  }
  other.extra.cols <- setdiff(extra.cols, response.cols)
  if (length(added.cols) > 0L) {
    for (nm in added.cols) {
      newdata[[nm]] <- rep(NA, nrow(newdata))
    }
    .msg("Added missing columns to newdata: ", paste(added.cols, collapse = ", "),
         verbose = verbose)
  }
  if (length(response.cols) > 0L) {
    .msg("Dropping supervised response columns from newdata: ",
         paste(response.cols, collapse = ", "),
         verbose = verbose)
  }
  if (length(other.extra.cols) > 0L) {
    .msg("Dropping extra columns from newdata: ", paste(other.extra.cols, collapse = ", "),
         verbose = verbose)
  }
  newdata <- newdata[, required.cols, drop = FALSE]
  conversion.issues <- list()
  unseen.levels <- vector("list", length(required.cols))
  names(unseen.levels) <- required.cols
  unseen.mask <- matrix(FALSE, nrow = nrow(newdata), ncol = length(required.cols),
                        dimnames = list(NULL, required.cols))
  for (nm in required.cols) {
    sc <- manifest$schema[[nm]]
    x <- newdata[[nm]]
    if (isTRUE(sc$is.factor)) {
      xchr <- as.character(x)
      bad <- !is.na(xchr) & !(xchr %in% sc$levels)
      unseen <- unique(xchr[bad])
      unseen <- stats::na.omit(unseen)
      if (length(unseen) > 0L) {
        unseen.levels[[nm]] <- as.character(unseen)
      }
      if (nrow(newdata) > 0L) {
        unseen.mask[, nm] <- bad
      }
      newdata[[nm]] <- .coerce.factor.levels(xchr, sc$levels, ordered = sc$ordered)
    }
    else {
      converted <- .as.numeric.from.schema(x, sc)
      if (any(!is.na(converted) & !is.finite(converted))) {
        stop("Numeric values in newdata column `", nm,
             "` must be finite or missing.", call. = FALSE)
      }
      bad <- !is.na(x) & is.na(converted)
      if (any(bad)) conversion.issues[[nm]] <- which(bad)
      newdata[[nm]] <- converted
    }
  }
  if (length(conversion.issues) > 0L) {
    warning("Nonmissing values could not be converted to numeric and were treated ",
            "as missing in columns: ", paste(names(conversion.issues), collapse = ", "),
            ". Row indices are recorded in conversion.issues.", call. = FALSE)
  }
  list(
    data = newdata,
    conversion.issues = conversion.issues,
    added.cols = added.cols,
    extra.cols = extra.cols,
    response.cols = response.cols,
    unseen.levels = unseen.levels,
    unseen.mask = as.data.frame(unseen.mask, optional = TRUE, stringsAsFactors = FALSE),
    unseen.rows = if (nrow(newdata) == 0L) logical(0) else rowSums(unseen.mask) > 0L
  )
}
## Fill currently missing cells from stored training means/modes, without model updates.
.apply.init <- function(data, init, schema) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  for (nm in names(data)) {
    miss <- is.na(data[[nm]])
    if (!any(miss)) next
    sc <- schema[[nm]]
    val <- init[[nm]]
    if (isTRUE(sc$is.factor)) {
      x <- as.character(data[[nm]])
      x[miss] <- as.character(val)
      data[[nm]] <- .coerce.factor.levels(x, sc$levels, ordered = sc$ordered)
    }
    else {
      x <- .as.numeric.from.schema(data[[nm]], sc)
      x[miss] <- as.numeric(val)
      data[[nm]] <- x
    }
  }
  data
}
## Restore factor encodings and the integer support of generated numeric cells.
## Observed numeric values are retained both as predictors and in returned data.
.restore.schema <- function(data, schema, restore.integer = TRUE,
                              generated = NULL) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  generated <- .impute.generated.mask(data, generated)
  for (nm in names(data)) {
    sc <- schema[[nm]]
    if (isTRUE(sc$is.factor)) {
      data[[nm]] <- .coerce.factor.levels(data[[nm]], sc$levels, ordered = sc$ordered)
    }
    else if (.schema.restores.integer(sc)) {
      x <- .as.numeric.from.schema(data[[nm]], sc)
      data[[nm]] <- .restore.integer.vector(
        x,
        .schema.integer.info(sc),
        restore.integer = restore.integer,
        generated = generated[, nm]
      )
    }
    else {
      data[[nm]] <- .as.numeric.from.schema(data[[nm]], sc)
    }
  }
  data
}
## Allocate a typed missing-value table for response metadata or reconstructions.
.na.data.from.schema <- function(schema, nrow) {
  if (is.null(schema) || length(schema) == 0L) {
    return(data.frame())
  }
  out <- setNames(vector("list", length(schema)), names(schema))
  for (nm in names(schema)) {
    sc <- schema[[nm]]
    if (isTRUE(sc$is.factor)) {
      out[[nm]] <- factor(rep(NA_character_, nrow),
                          levels = sc$levels,
                          ordered = sc$ordered)
    }
    else if (isTRUE(sc$is.integer) || isTRUE(sc$integer.storage)) {
      out[[nm]] <- rep(NA_integer_, nrow)
    }
    else {
      out[[nm]] <- rep(NA_real_, nrow)
    }
  }
  data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
}
##------------------------------------------------------------------
## Prediction extraction and target reconstruction discrepancies
##------------------------------------------------------------------
## Extract regression predictions or class labels for imputation/reconstruction.
## Numeric integer-grid restoration is controlled by the caller.
.target.prediction.vector <- function(x) {
  if (is.data.frame(x)) {
    if (ncol(x) != 1L) stop("Expected one prediction column.", call. = FALSE)
    x <- x[[1L]]
  }
  if (!is.null(dim(x))) {
    if (length(dim(x)) > 2L || (length(dim(x)) == 2L && ncol(x) != 1L)) {
      stop("Expected one prediction per row.", call. = FALSE)
    }
    dim(x) <- NULL
  }
  if (!is.null(x) && !is.atomic(x)) {
    stop("Predictions must be an atomic vector.", call. = FALSE)
  }
  x
}
.validate.target.values <- function(values, n, schema) {
  values <- .target.prediction.vector(values)
  if (is.null(values) || length(values) != n) {
    stop("Prediction length (", length(values), ") does not match requested rows (",
         n, ").", call. = FALSE)
  }
  if (isTRUE(schema$is.factor)) {
    valid <- !is.na(values) & as.character(values) %in% schema$levels
  } else {
    if (!is.numeric(values)) stop("Numeric predictions are required.", call. = FALSE)
    valid <- is.finite(values)
  }
  list(values = values, valid = valid)
}
.extract.prediction <- function(pred, family, target.schema,
                                restore.integer = TRUE) {
  if (is.null(pred)) return(NULL)
  if (identical(family, "regr")) {
    out <- .target.prediction.vector(pred$predicted)
    if (!is.null(out) && !is.numeric(out)) {
      stop("Numeric predictions are required.", call. = FALSE)
    }
    if (.schema.restores.integer(target.schema)) {
      out <- .restore.integer.vector(
        out,
        .schema.integer.info(target.schema),
        restore.integer = restore.integer
      )
    }
    return(out)
  }
  out <- pred$class
  if (is.null(out) && !is.null(pred$predicted)) {
    out <- pred$predicted
  }
  out <- .target.prediction.vector(out)
  if (isTRUE(target.schema$is.factor)) {
    out <- .coerce.factor.levels(out, target.schema$levels,
                                 ordered = target.schema$ordered)
  }
  out
}
## Extract the classification probability matrix; recover labels from schema if possible.
.extract.class.prob <- function(pred, target.schema = NULL) {
  if (is.null(pred)) return(NULL)
  prob <- pred$predicted
  if (is.data.frame(prob)) {
    prob <- as.matrix(prob)
  }
  if (!is.matrix(prob)) {
    return(NULL)
  }
  if (is.null(colnames(prob)) && !is.null(target.schema$levels) &&
      ncol(prob) == length(target.schema$levels)) {
    colnames(prob) <- target.schema$levels
  }
  prob
}
## Compute absolute reconstruction errors, rounding integer-supported predictions.
## Observed missing target values propagate to missing discrepancies.
.ood.delta.regression <- function(observed, pred, target.schema = NULL) {
  if (is.null(pred) || is.null(pred$predicted)) {
    return(rep(NA_real_, length(observed)))
  }
  values <- .extract.prediction(pred, "regr", target.schema)
  checked <- .validate.target.values(values, length(observed), target.schema)
  obs <- .as.numeric.from.schema(observed, target.schema)
  out <- rep(NA_real_, length(obs))
  ok <- !is.na(obs) & checked$valid
  out[ok] <- abs(obs[ok] - checked$values[ok])
  out
}
## Use negative log probability of the observed class when probabilities are present.
## Otherwise use misclassification or scaled ordered-level distance.
.ood.delta.classification <- function(observed, pred, target.schema,
                                      prob.floor = 1e-12) {
  if (!is.numeric(prob.floor) || length(prob.floor) != 1L ||
      !is.finite(prob.floor) || prob.floor <= 0 || prob.floor >= 1) {
    stop("'prob.floor' must satisfy 0 < prob.floor < 1.", call. = FALSE)
  }
  y <- as.character(observed)
  out <- rep(NA_real_, length(y))
  missing.y <- is.na(y)
  prob <- .extract.class.prob(pred, target.schema)
  if (!is.null(prob)) {
    if (!is.numeric(prob) || nrow(prob) != length(y) ||
        is.null(colnames(prob)) || anyNA(colnames(prob)) ||
        anyDuplicated(colnames(prob))) {
      stop("Class probabilities must have one row per observation and named classes.",
           call. = FALSE)
    }
    lev <- colnames(prob)
    pos <- match(y, lev)
    ok <- !missing.y & !is.na(pos)
    if (any(ok)) {
      p.true <- prob[cbind(which(ok), pos[ok])]
      ## Unavailable/invalid probabilities remain missing reference values.
      ## A genuine finite zero probability is still scored using the floor.
      valid <- is.finite(p.true) & p.true >= 0 & p.true <= 1
      rows <- which(ok)[valid]
      out[rows] <- -log(pmax(p.true[valid], prob.floor))
    }
    if (any(!missing.y & is.na(pos))) {
      out[!missing.y & is.na(pos)] <- Inf
    }
    return(out)
  }
  cls <- pred$class
  if (is.null(cls) && !is.null(pred$predicted)) {
    cls <- pred$predicted
  }
  if (is.null(cls)) return(out)
  cls <- as.character(.target.prediction.vector(cls))
  if (length(cls) != length(y)) {
    stop("Class prediction length does not match observed rows.", call. = FALSE)
  }
  ok <- !missing.y & !is.na(cls) & cls %in% target.schema$levels
  if (!any(ok)) {
    return(out)
  }
  if (isTRUE(target.schema$ordered)) {
    obs.code <- match(y[ok], target.schema$levels)
    cls.code <- match(cls[ok], target.schema$levels)
    denom <- max(1L, length(target.schema$levels) - 1L)
    out[ok] <- abs(obs.code - cls.code) / denom
  }
  else {
    out[ok] <- as.numeric(y[ok] != cls[ok])
  }
  out
}
## Dispatch target discrepancy calculation using the raw target schema.
.compute.ood.delta <- function(observed, pred, target.schema) {
  if (isTRUE(target.schema$is.factor)) {
    .ood.delta.classification(observed, pred, target.schema = target.schema)
  }
  else {
    .ood.delta.regression(observed, pred, target.schema = target.schema)
  }
}
##------------------------------------------------------------------
## Compact OOD references and row aggregation
##------------------------------------------------------------------
## Compress finite training discrepancies into an ordered type-8 quantile grid.
## n records the number of finite values before compression.
.make.ood.reference <- function(x, probs = seq(0, 1, length.out = 257)) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  probs <- unique(sort(pmin(1, pmax(0, probs))))
  if (length(probs) < 2L) {
    probs <- c(0, 1)
  }
  if (length(x) == 0L) {
    return(list(
      probs = c(0, 1),
      quantiles = c(0, 0),
      n = 0L
    ))
  }
  q <- as.numeric(stats::quantile(x, probs = probs, names = FALSE,
                                  na.rm = TRUE, type = 8))
  q <- cummax(q)
  list(
    probs = probs,
    quantiles = q,
    n = length(x)
  )
}
## Map discrepancies to stored quantile-grid probability positions.
## Constant references use a separate two-level mapping.
.eval.ood.reference <- function(x, ref) {
  out <- rep(NA_real_, length(x))
  if (length(out) == 0L || is.null(ref) || is.null(ref$quantiles) ||
      (!is.null(ref$n) && (!is.finite(ref$n) || ref$n <= 0))) {
    return(out)
  }
  q <- as.numeric(ref$quantiles)
  p <- as.numeric(ref$probs)
  if (length(q) == 0L || length(q) != length(p) ||
      any(!is.finite(q)) || any(!is.finite(p)) ||
      is.unsorted(q) || is.unsorted(p) || any(p < 0 | p > 1)) {
    return(out)
  }
  ok <- !is.na(x)
  if (!any(ok)) return(out)
  ## Preserve the established constant-reference convention, including equality.
  if (length(q) == 1L || all(q == q[1L])) {
    out[ok] <- ifelse(x[ok] <= q[1L], p[1L], p[length(p)])
    return(out)
  }
  ## Preserve interior ties and equality at the largest quantile. Strict upper
  ## tails (including +Inf) use the last probability; NA/NaN remain unscored.
  idx <- findInterval(x[ok], q, rightmost.closed = TRUE, all.inside = FALSE)
  idx <- pmax(1L, pmin(length(p), idx))
  out[ok] <- p[idx]
  out
}
## Resolve a row-aggregation name and its applicable tuning arguments.
.canonicalize.ood.aggregate <- function(aggregate = "weighted.mean",
                                        aggregate.args = NULL) {
  choices <- c(
    "bounded.product",
    "weighted.mean",
    "weighted.lp",
    "weighted.lp.log",
    "top.k"
  )
  if (is.null(aggregate) || length(aggregate) == 0L) {
    aggregate <- "weighted.mean"
  }
  aggregate <- match.arg(aggregate, choices = choices)
  aggregate.args <- .impute.named.options(aggregate.args, "aggregate.args")
  used <- character(0)
  args <- list()
  if (aggregate %in% c("weighted.lp", "weighted.lp.log")) {
    p <- aggregate.args$p %||% 2
    if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p < 1) {
      stop("The row aggregate power 'p' must be a finite scalar >= 1.",
           call. = FALSE)
    }
    args$p <- as.numeric(p)
    used <- c(used, "p")
  }
  if (aggregate == "weighted.lp.log") {
    agg.eps <- aggregate.args$eps %||% 1e-12
    if (!is.numeric(agg.eps) || length(agg.eps) != 1L ||
        !is.finite(agg.eps) || agg.eps <= 0 || agg.eps >= 1) {
      stop("The row aggregate 'eps' must satisfy 0 < eps < 1.",
           call. = FALSE)
    }
    args$eps <- as.numeric(agg.eps)
    used <- c(used, "eps")
  }
  if (aggregate == "top.k") {
    k <- aggregate.args$k %||% aggregate.args$top.k %||% 1L
    if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 1 || k > .Machine$integer.max) {
      stop("The top-k row aggregate requires 1 <= k <= .Machine$integer.max.",
           call. = FALSE)
    }
    args$k <- as.integer(max(1L, round(k)))
    used <- c(used, "k", "top.k")
  }
  if (aggregate == "bounded.product") {
    agg.eps <- aggregate.args$eps %||% 1e-12
    if (!is.numeric(agg.eps) || length(agg.eps) != 1L ||
        !is.finite(agg.eps) || agg.eps <= 0 || agg.eps >= 1) {
      stop("The bounded-product row aggregate requires 0 < eps < 1.",
           call. = FALSE)
    }
    args$eps <- as.numeric(agg.eps)
    used <- c(used, "eps")
  }
  extra <- setdiff(names(aggregate.args), unique(used))
  if (length(extra) > 0L) {
    warning("Ignoring unknown aggregate.args entries: ",
            paste(extra, collapse = ", "),
            call. = FALSE)
  }
  list(name = aggregate, args = args)
}
## Compare canonical aggregation specifications for legacy row-reference reuse.
.same.ood.aggregate <- function(aggregate, aggregate.args = NULL,
                                ref.aggregate, ref.aggregate.args = NULL,
                                tolerance = sqrt(.Machine$double.eps)) {
  ax <- tryCatch(
    .canonicalize.ood.aggregate(aggregate, aggregate.args),
    error = function(e) NULL
  )
  ay <- tryCatch(
    .canonicalize.ood.aggregate(ref.aggregate, ref.aggregate.args),
    error = function(e) NULL
  )
  if (is.null(ax) || is.null(ay)) {
    return(FALSE)
  }
  if (!identical(ax$name, ay$name)) {
    return(FALSE)
  }
  if (!setequal(names(ax$args), names(ay$args))) {
    return(FALSE)
  }
  if (length(ax$args) == 0L) {
    return(TRUE)
  }
  isTRUE(all(vapply(names(ax$args), function(nm) {
    isTRUE(all.equal(ax$args[[nm]], ay$args[[nm]],
                     tolerance = tolerance,
                     check.attributes = FALSE))
  }, logical(1))))
}
## Combine one row's finite target scores using positive weights only.
## The denominator renormalizes weights over targets scoreable in this row.
.ood.aggregate.vector <- function(x, weight, aggregate.spec) {
  ok <- is.finite(x) & is.finite(weight) & weight > 0
  if (!any(ok)) return(NA_real_)
  x <- as.numeric(x[ok])
  weight <- as.numeric(weight[ok])
  ## Top-k selection precedes normalization: weights outside the selected set
  ## must not cause very small selected weights to underflow to zero.
  if (aggregate.spec$name == "top.k") {
    k <- min(length(x), aggregate.spec$args$k)
    sel <- order(x, decreasing = TRUE, na.last = NA)[seq_len(k)]
    x <- x[sel]
    weight <- weight[sel]
  }
  log.weight <- log(weight)
  log.weight <- log.weight - max(log.weight)
  log.weight <- log.weight - log(sum(exp(log.weight)))
  weight <- weight / max(weight)
  weight <- weight / sum(weight)
  ## Use log weights and max scaling for the norm. This also retains the
  ## influence of very small weights when a large power magnifies a tail score.
  lp <- function(z, p) {
    scale <- max(abs(z))
    if (scale == 0) return(0)
    positive <- abs(z) > 0
    terms <- log.weight[positive] + p * (log(abs(z[positive])) - log(scale))
    shift <- max(terms)
    log.norm <- shift + log(sum(exp(terms - shift)))
    min(scale, exp(log(scale) + min(0, log.norm) / p))
  }
  switch(
    aggregate.spec$name,
    weighted.mean = sum(weight * x),
    weighted.lp = lp(x, aggregate.spec$args$p),
    weighted.lp.log = {
      u <- pmin(1, pmax(0, x))
      z <- pmin(-log1p(-u), -log(aggregate.spec$args$eps))
      lp(z, aggregate.spec$args$p)
    },
    top.k = sum(weight * x),
    bounded.product = {
      u <- pmin(1, pmax(0, x))
      log.surv <- pmax(log1p(-u), log(aggregate.spec$args$eps))
      pmin(1, pmax(0, -expm1(sum(weight * log.surv))))
    },
    stop("Unknown OOD row aggregate: ", aggregate.spec$name, call. = FALSE)
  )
}
## Return the sentinel maximum used by the row-wide unseen-factor override.
.max.ood.aggregate.value <- function(aggregate = "weighted.mean",
                                    aggregate.args = NULL) {
  aggregate.spec <- .canonicalize.ood.aggregate(aggregate, aggregate.args)
  switch(
    aggregate.spec$name,
    weighted.lp.log = -log(aggregate.spec$args$eps),
    1
  )
}
## Apply the chosen aggregate row by row; return NA for rows with no usable targets.
.aggregate.ood.row <- function(mat, weight = NULL,
                               aggregate = "weighted.mean",
                               aggregate.args = NULL) {
  mat <- as.matrix(mat)
  if (ncol(mat) == 0L) {
    return(rep(NA_real_, nrow(mat)))
  }
  if (is.null(weight)) {
    weight <- rep(1, ncol(mat))
  }
  if (length(weight) != ncol(mat)) {
    stop("'weight' must have one entry per target.", call. = FALSE)
  }
  weight <- as.numeric(weight)
  aggregate.spec <- .canonicalize.ood.aggregate(aggregate, aggregate.args)
  out <- rep(NA_real_, nrow(mat))
  for (i in seq_len(nrow(mat))) {
    out[i] <- .ood.aggregate.vector(mat[i, ], weight, aggregate.spec)
  }
  out
}
## Recover a numeric training-score matrix and its target-column labels.
.as.ood.train.score.matrix <- function(x, targets = NULL) {
  if (is.null(x)) {
    return(NULL)
  }
  mat <- tryCatch(as.matrix(x), error = function(e) NULL)
  if (is.null(mat)) {
    return(NULL)
  }
  storage.mode(mat) <- "double"
  if (is.null(colnames(mat)) && !is.null(targets) && ncol(mat) == length(targets)) {
    colnames(mat) <- targets
  }
  mat
}
## Reaggregate saved training scores using the active targets, weights and rule.
## Return a reason when rebuilding is unavailable so callers can assess legacy fallback.
.rebuild.ood.row.reference <- function(train.target.score, targets, weight,
                                       saved.targets = NULL,
                                       aggregate = "weighted.mean",
                                       aggregate.args = NULL) {
  mat <- .as.ood.train.score.matrix(train.target.score, targets = saved.targets)
  if (is.null(mat)) {
    return(list(
      reference = NULL,
      row.score = NULL,
      n.finite = 0L,
      reason = paste(
        "No saved training OOD target scores are available.",
        "Refit with the updated imputer to enable row-level percentile",
        "recalibration for arbitrary target subsets, weights, and row aggregates."
      )
    ))
  }
  if (length(targets) == 0L || ncol(mat) == 0L) {
    return(list(
      reference = NULL,
      row.score = rep(NA_real_, nrow(mat)),
      n.finite = 0L,
      reason = "No saved training OOD target scores are available for the requested targets."
    ))
  }
  if (is.null(colnames(mat))) {
    return(list(
      reference = NULL,
      row.score = NULL,
      n.finite = 0L,
      reason = paste(
        "Saved training OOD target scores lack target names.",
        "Refit with the updated imputer to enable row-level percentile",
        "recalibration for arbitrary target subsets, weights, and row aggregates."
      )
    ))
  }
  missing.targets <- setdiff(targets, colnames(mat))
  if (length(missing.targets) > 0L) {
    return(list(
      reference = NULL,
      row.score = NULL,
      n.finite = 0L,
      reason = paste0(
        "Saved training OOD target scores are missing requested targets: ",
        paste(missing.targets, collapse = ", "), "."
      )
    ))
  }
  row.score <- .aggregate.ood.row(
    mat[, targets, drop = FALSE],
    weight = weight,
    aggregate = aggregate,
    aggregate.args = aggregate.args
  )
  n.finite <- sum(is.finite(row.score))
  if (n.finite == 0L) {
    return(list(
      reference = NULL,
      row.score = row.score,
      n.finite = 0L,
      reason = paste(
        "No finite row-level training OOD scores were available for the",
        "requested targets, weights, and row aggregates."
      )
    ))
  }
  list(
    reference = .make.ood.reference(row.score),
    row.score = row.score,
    n.finite = n.finite,
    reason = NULL
  )
}
## Align optional named weights to targets; unmatched/omitted named targets receive zero.
## Require finite nonnegative weights and at least one positive active target.
.resolve.ood.weight <- function(targets, weight = NULL, default = NULL,
                                warn.extra = TRUE) {
  targets <- as.character(targets %||% character(0))
  if (length(targets) == 0L) {
    return(setNames(numeric(0), character(0)))
  }
  validate.named.weight <- function(x, what = "'weight'") {
    if (is.null(names(x))) {
      return(invisible(NULL))
    }
    dup <- unique(names(x)[duplicated(names(x))])
    if (length(dup) > 0L) {
      stop(what, " contains duplicated target names: ",
           paste(dup, collapse = ", "),
           call. = FALSE)
    }
    invisible(NULL)
  }
  align.named.weight <- function(x, what = "'weight'",
                                 warn.extra = TRUE,
                                 require.match = TRUE,
                                 fill = 0) {
    validate.named.weight(x, what = what)
    out <- setNames(rep(fill, length(targets)), targets)
    extra.targets <- setdiff(names(x), targets)
    if (isTRUE(warn.extra) && length(extra.targets) > 0L) {
      warning("Ignoring ", what, " entries for unknown targets: ",
              paste(extra.targets, collapse = ", "),
              call. = FALSE)
    }
    matched.targets <- intersect(targets, names(x))
    if (length(matched.targets) == 0L) {
      if (isTRUE(require.match)) {
        stop(what, " did not match any requested targets.", call. = FALSE)
      }
      return(out)
    }
    out[matched.targets] <- as.numeric(x[matched.targets])
    out
  }
  if (is.null(weight) || length(weight) == 0L) {
    if (is.null(default) || length(default) == 0L) {
      weight <- setNames(rep(1, length(targets)), targets)
    }
    else if (is.null(names(default))) {
      if (length(default) != length(targets)) {
        stop("'default' must have one entry per target.", call. = FALSE)
      }
      weight <- setNames(as.numeric(default), targets)
    }
    else {
      weight <- align.named.weight(
        x = default,
        what = "'default'",
        warn.extra = FALSE,
        require.match = FALSE,
        fill = 0
      )
    }
  }
  else if (is.null(names(weight))) {
    if (length(weight) != length(targets)) {
      stop("'weight' must have one entry per target.", call. = FALSE)
    }
    weight <- setNames(as.numeric(weight), targets)
  }
  else {
    weight <- align.named.weight(
      x = weight,
      what = "'weight'",
      warn.extra = warn.extra,
      require.match = TRUE,
      fill = 0
    )
  }
  if (any(!is.finite(weight) | weight < 0)) {
    stop("'weight' must contain finite nonnegative values.", call. = FALSE)
  }
  if (!any(weight > 0)) {
    stop("At least one target weight must be positive.", call. = FALSE)
  }
  weight
}
## Compare target-aligned weights when deciding whether an old row reference is reusable.
.same.ood.weight <- function(targets, x, y,
                             tolerance = sqrt(.Machine$double.eps)) {
  if (is.null(x) || is.null(y)) {
    return(FALSE)
  }
  wx <- tryCatch(
    .resolve.ood.weight(targets, x, warn.extra = FALSE),
    error = function(e) NULL
  )
  wy <- tryCatch(
    .resolve.ood.weight(targets, y, warn.extra = FALSE),
    error = function(e) NULL
  )
  if (is.null(wx) || is.null(wy)) {
    return(FALSE)
  }
  isTRUE(all.equal(unname(wx), unname(wy), tolerance = tolerance,
                   check.attributes = FALSE))
}
##------------------------------------------------------------------
## Prediction-pass diagnostics and shared deployment preparation
##------------------------------------------------------------------
## Average target-wise changes over originally missing cells.
## Factors use disagreement rates; numerics use RMS changes scaled by training SD.
.compute.pass.diff <- function(old.data, new.data, missing.mask, schema, scale, targets) {
  diffs <- sapply(targets, function(y) {
    idx <- missing.mask[[y]]
    if (!any(idx)) return(NA_real_)
    xo <- old.data[[y]][idx]
    xn <- new.data[[y]][idx]
    sc <- schema[[y]]
    if (isTRUE(sc$is.factor)) {
      return(sum(as.character(xn) != as.character(xo), na.rm = TRUE) /
               (0.001 + length(xn)))
    }
    xo <- .as.numeric.from.schema(xo, sc)
    xn <- .as.numeric.from.schema(xn, sc)
    sy <- scale[[y]] %||% 1
    if (!is.finite(sy) || sy <= 0) sy <- 1
    sqrt(mean((xn - xo)^2, na.rm = TRUE) / (0.001 + sy^2))
  })
  mean(diffs, na.rm = TRUE)
}
## Shared prediction/OOD preparation. Preserve a pre-fill harmonized table and
## maintain separate raw and augmented working tables through ordered updates.
## Only originally missing cells in the requested targets receive forest updates.
.prepare.impute.learn.newdata <- function(object, newdata,
                                          targets = NULL,
                                          max.predict.iter = 3L,
                                          eps = 1e-3,
                                          restore.integer = TRUE,
                                          cache.learners = c("session", "none", "all"),
                                          verbose = TRUE) {
  if (!inherits(object, "impute.learn.rfsrc")) {
    stop("'object' must inherit from class 'impute.learn.rfsrc'.", call. = FALSE)
  }
  if (missing(newdata)) stop("'newdata' is missing.", call. = FALSE)
  .impute.scalar(max.predict.iter, "max.predict.iter", lower = 0, whole = TRUE,
                 upper = .Machine$integer.max)
  .impute.scalar(eps, "eps", lower = 0)
  .impute.flag(restore.integer, "restore.integer")
  .impute.flag(verbose, "verbose")
  if (!is.null(targets) && (!is.character(targets) || anyNA(targets))) {
    stop("'targets' must be NULL or a character vector without missing names.",
         call. = FALSE)
  }
  object$manifest <- .normalize.impute.learn.manifest(object$manifest)
  if (!is.null(targets)) {
    bad.targets <- setdiff(targets, object$manifest$targets)
    if (length(bad.targets) > 0L) {
      warning("Ignoring unknown targets: ",
              paste(bad.targets, collapse = ", "),
              call. = FALSE)
    }
  }
  use.targets <- if (is.null(targets)) object$manifest$targets else {
    intersect(object$manifest$targets, targets)
  }
  if (length(use.targets) == 0L) {
    stop("No valid targets requested.", call. = FALSE)
  }
  cache.learners <- match.arg(cache.learners)
  harmonized <- .harmonize.newdata(newdata, object$manifest, verbose = verbose)
  data <- harmonized$data
  ## Capture missingness in every raw column before initialization, including
  ## columns outside the requested update targets. Their initialization values
  ## are generated too. harmonized$data remains the pre-fill OOD comparison table.
  original.missing <- as.data.frame(is.na(data))
  names(original.missing) <- names(data)
  data <- .apply.init(data, object$manifest$init, object$manifest$schema)
  ## Restore generated initialization values before they become predictors.
  ## Observed fractions are not rounded to the saved training integer grid.
  data <- .restore.schema(data, object$manifest$schema,
                          restore.integer = restore.integer,
                          generated = original.missing)
  working.data <- data
  supervised.info <- list(enabled = FALSE)
  ## Compute supervised auxiliaries once from initialized raw data. They are
  ## held fixed during all subsequent prediction-time target sweeps.
  if (.is.supervised.manifest(object$manifest)) {
    if (nrow(data) == 0L) {
      aux <- as.data.frame(matrix(numeric(0), nrow = 0L,
                                  ncol = length(object$manifest$supervised$aux.names)))
      names(aux) <- object$manifest$supervised$aux.names
      sup.info <- list(loaded.from.disk = FALSE)
    } else {
      sup.info <- .predict.get.supervised.model(object)
      if (is.null(sup.info$model)) {
        stop(sup.info$error %||% "The supervised forest could not be loaded.",
             call. = FALSE)
      }
      aux <- tryCatch(
        {
          pred.newdata <- .conform.x.to.forest(
            data,
            sup.info$model,
            ignore.levels = TRUE
          )
          pred <- predict(sup.info$model, pred.newdata)
          aux <- .supervised.materialize.aux(
            pred,
            family = object$manifest$supervised$family,
            oob = FALSE,
            response.schema = object$manifest$supervised$response.schema,
            yvar.names = object$manifest$supervised$yvar.names
          )
          aux <- .align.supervised.aux(
            aux,
            expected.names = object$manifest$supervised$aux.names
          )
          aux <- .fill.supervised.aux(aux, object$manifest$supervised$aux.init)
          if (nrow(aux) != nrow(data) || any(!is.finite(as.matrix(aux)))) {
            stop("Supervised auxiliaries must be finite and have one row per observation.",
                 call. = FALSE)
          }
          aux
        },
        error = function(e) e
      )
      if (inherits(aux, "error")) {
        stop("Failed to compute supervised auxiliary predictors for newdata: ",
             conditionMessage(aux), call. = FALSE)
      }
    }
    working.data <- data.frame(data, aux, check.names = FALSE)
    supervised.info <- list(
      enabled = TRUE,
      family = object$manifest$supervised$family,
      aux.names = names(aux),
      n.aux = ncol(aux),
      response.names = names(object$manifest$supervised$response.schema %||% list()),
      dropped.response.columns = harmonized$response.cols,
      loaded.from.disk = isTRUE(sup.info$loaded.from.disk)
    )
  }
  pass.history <- numeric(0)
  pass.updated.cells <- integer(0)
  pass.failed.cells <- integer(0)
  converged <- FALSE
  stopping.reason <- if (nrow(data) == 0L) "empty.newdata" else "no.missing"
  sweep.order <- object$manifest$sweep.order
  sweep.order <- sweep.order[sweep.order %in% use.targets]
  ## Cache lifetime is this preparation/scoring call, not the R session.
  ## OOD scoring receives this environment and can reuse models already loaded.
  cache.env <- if (identical(cache.learners, "none")) NULL else new.env(parent = emptyenv())
  disk.load.targets <- character(0)
  n.disk.loads <- 0L
  target.issues <- setNames(vector("list", length(use.targets)), use.targets)
  record.issue <- function(target, message) {
    current <- target.issues[[target]]
    if (is.null(current)) current <- character(0)
    if (!(message %in% current)) {
      target.issues[[target]] <<- c(current, message)
    }
    invisible(NULL)
  }
  if (nrow(data) > 0L && identical(cache.learners, "all")) {
    .msg("Preloading learner bank...", verbose = verbose)
    for (target in use.targets) {
      info <- object$manifest$learners[[target]]
      if (!identical(info$status, "ok")) next
      mdl.info <- .predict.get.model(object, target, cache.env = cache.env)
      if (isTRUE(mdl.info$loaded.from.disk)) {
        disk.load.targets <- c(disk.load.targets, target)
        n.disk.loads <- n.disk.loads + 1L
      }
      if (is.null(mdl.info$model) && !is.null(mdl.info$error)) {
        record.issue(target, mdl.info$error)
      }
    }
  }
  any.target.missing <- length(sweep.order) > 0L &&
    any(as.matrix(original.missing[, sweep.order, drop = FALSE]))
  if (isTRUE(any.target.missing)) {
    stopping.reason <- if (max.predict.iter == 0L) "initialization.only" else "max.iter"
    .msg("Starting prediction-time sweep...", verbose = verbose)
    for (iter in seq_len(max.predict.iter)) {
      old.data <- data
      n.updated <- 0L
      n.failed <- 0L
      .msg("  prediction pass ", iter, "/", max.predict.iter, verbose = verbose)
      ## Ordered updates: later targets use earlier targets' current imputations
      ## from working.data. Auxiliary predictor columns are not recomputed.
      for (target in sweep.order) {
        miss.idx <- which(original.missing[[target]])
        if (length(miss.idx) == 0L) next
        info <- object$manifest$learners[[target]]
        if (!identical(info$status, "ok")) {
          n.failed <- n.failed + length(miss.idx)
          msg <- paste0("No trained learner is available (status = ",
                        info$status %||% "unknown", ").")
          record.issue(target, msg)
          .msg("    skipping `", target, "` (", msg, ")", verbose = verbose)
          next
        }
        mdl.info <- .predict.get.model(object, target, cache.env = cache.env)
        mdl <- mdl.info$model
        if (isTRUE(mdl.info$loaded.from.disk)) {
          disk.load.targets <- c(disk.load.targets, target)
          n.disk.loads <- n.disk.loads + 1L
        }
        if (is.null(mdl)) {
          n.failed <- n.failed + length(miss.idx)
          msg <- mdl.info$error %||% "learner could not be loaded"
          record.issue(target, msg)
          .msg("    skipping `", target, "` (", msg, ")", verbose = verbose)
          next
        }
        ## Conformity, prediction, extraction and shape checks are one
        ## target operation. Invalid entries retain their previous fill.
        result <- tryCatch({
          xvars <- object$manifest$predictor.map[[target]]
          pred.df <- working.data[miss.idx, xvars, drop = FALSE]
          pred.df <- .conform.x.to.forest(pred.df, mdl)
          pred <- predict(mdl, pred.df)
          values <- .extract.prediction(
            pred, info$family, object$manifest$schema[[target]],
            restore.integer = restore.integer
          )
          .validate.target.values(values, length(miss.idx),
                                   object$manifest$schema[[target]])
        }, error = function(e) e)
        if (inherits(result, "error")) {
          msg <- paste0("Prediction failed: ", conditionMessage(result))
          record.issue(target, msg)
          .msg("    ", target, ": ", msg, verbose = verbose)
          n.failed <- n.failed + length(miss.idx)
        } else {
          invalid <- sum(!result$valid)
          if (invalid > 0L) {
            record.issue(target, paste0(invalid,
              " unavailable prediction(s); previous fills retained."))
          }
          rows <- miss.idx[result$valid]
          if (length(rows) > 0L) {
            data[rows, target] <- result$values[result$valid]
            working.data[rows, target] <- result$values[result$valid]
          }
          n.updated <- n.updated + length(rows)
          n.failed <- n.failed + invalid
        }
        if (identical(cache.learners, "none") && is.null(object$models[[target]])) {
          rm(mdl)
          gc()
        }
      }
      ## A pass with no valid model updates is not evidence of convergence.
      diff.err <- if (n.updated == 0L) NA_real_ else {
        .compute.pass.diff(old.data, data, original.missing,
                            object$manifest$schema, object$manifest$scale, sweep.order)
      }
      pass.history <- c(pass.history, diff.err)
      pass.updated.cells <- c(pass.updated.cells, n.updated)
      pass.failed.cells <- c(pass.failed.cells, n.failed)
      .msg("    pass diff = ", format(diff.err, digits = 4), verbose = verbose)
      if (n.updated == 0L) {
        stopping.reason <- "no.valid.updates"
        .msg("    no valid model updates; retaining previous fills.", verbose = verbose)
        break
      }
      if (is.finite(diff.err) && diff.err < eps) {
        converged <- n.failed == 0L
        stopping.reason <- if (converged) "converged" else "stable.with.failures"
        .msg(if (converged) "    convergence criterion met; stopping early." else
               "    available updates are stable; some cells retain fallback fills.",
             verbose = verbose)
        break
      }
    }
  }
  else {
    .msg("No missing values were found among requested targets; ",
         "skipping iterative sweep.", verbose = verbose)
  }
  data <- .restore.schema(data, object$manifest$schema,
                          restore.integer = restore.integer,
                          generated = original.missing)
  working.data[names(data)] <- data
  target.issues <- target.issues[lengths(target.issues) > 0L]
  list(
    data = data,
    working.data = working.data,
    use.targets = use.targets,
    harmonized = harmonized,
    cache.learners = cache.learners,
    cache.env = cache.env,
    info = list(
      n.passes = length(pass.history),
      pass.diff = pass.history,
      pass.updated.cells = pass.updated.cells,
      pass.failed.cells = pass.failed.cells,
      converged = converged,
      stopping.reason = stopping.reason,
      conversion.issues = harmonized$conversion.issues,
      targets = use.targets,
      added.columns = harmonized$added.cols,
      dropped.extra.columns = harmonized$extra.cols,
      dropped.response.columns = harmonized$response.cols,
      unseen.levels = harmonized$unseen.levels,
      unseen.rows = harmonized$unseen.rows,
      cache.learners = cache.learners,
      n.disk.loads = n.disk.loads,
      disk.load.targets = unique(disk.load.targets),
      supervised = supervised.info,
      target.issues = target.issues
    )
  )
}
##------------------------------------------------------------------
## Lazy learner retrieval
##------------------------------------------------------------------
## Load one target learner using its manifest entry, with optional nonfatal errors.
.fast.load.learner <- function(target, info, learner.root, strict = TRUE) {
  .check.fst()
  out <- tryCatch(
    fast.load(info$learner.name, learner.root),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    msg <- paste0("Failed to load learner for `", target, "` from ",
                  file.path(learner.root, info$learner.name), ": ",
                  conditionMessage(out))
    if (isTRUE(strict)) {
      stop(msg, call. = FALSE)
    }
    return(list(model = NULL, error = msg))
  }
  if (isTRUE(strict)) {
    return(out)
  }
  list(model = out, error = NULL)
}
## Resolve a target learner from this call's cache, the object, or its saved path.
## A successful disk read is cached only when a cache environment is supplied.
.predict.get.model <- function(object, target, cache.env = NULL) {
  if (!inherits(object, "impute.learn.rfsrc")) {
    stop("'object' must inherit from class 'impute.learn.rfsrc'.", call. = FALSE)
  }
  object$manifest <- .normalize.impute.learn.manifest(object$manifest)
  if (!is.null(cache.env) && exists(target, envir = cache.env, inherits = FALSE)) {
    return(list(
      model = get(target, envir = cache.env, inherits = FALSE),
      loaded.from.disk = FALSE,
      cached = TRUE,
      error = NULL
    ))
  }
  if (!is.null(object$models[[target]])) {
    mdl <- object$models[[target]]
    if (!is.null(cache.env)) assign(target, mdl, envir = cache.env)
    return(list(model = mdl, loaded.from.disk = FALSE,
                cached = !is.null(cache.env), error = NULL))
  }
  info <- object$manifest$learners[[target]]
  if (is.null(info) || !identical(info$status, "ok")) {
    return(list(model = NULL, loaded.from.disk = FALSE,
                cached = FALSE, error = NULL))
  }
  if (is.null(object$path)) {
    return(list(
      model = NULL,
      loaded.from.disk = FALSE,
      cached = FALSE,
      error = paste0("Learner for `", target, "` is not available in memory ",
                     "and no saved imputer path is attached to 'object'.")
    ))
  }
  learner.root <- file.path(object$path, object$manifest$learner.root)
  loaded <- .fast.load.learner(target, info, learner.root, strict = FALSE)
  mdl <- loaded$model
  if (!is.null(cache.env) && !is.null(mdl)) {
    assign(target, mdl, envir = cache.env)
  }
  list(model = mdl,
       loaded.from.disk = !is.null(mdl),
       cached = !is.null(cache.env) && !is.null(mdl),
       error = loaded$error)
}
## Resolve the optional supervised forest from memory or its saved learner directory.
.predict.get.supervised.model <- function(object) {
  if (!inherits(object, "impute.learn.rfsrc")) {
    stop("'object' must inherit from class 'impute.learn.rfsrc'.", call. = FALSE)
  }
  object$manifest <- .normalize.impute.learn.manifest(object$manifest)
  if (!.is.supervised.manifest(object$manifest)) {
    return(list(
      model = NULL,
      loaded.from.disk = FALSE,
      enabled = FALSE,
      error = NULL
    ))
  }
  if (!is.null(object$supervised.model)) {
    return(list(
      model = object$supervised.model,
      loaded.from.disk = FALSE,
      enabled = TRUE,
      error = NULL
    ))
  }
  if (is.null(object$path)) {
    return(list(
      model = NULL,
      loaded.from.disk = FALSE,
      enabled = TRUE,
      error = "The supervised forest is not available in memory and no saved imputer path is attached to 'object'."
    ))
  }
  learner.root <- file.path(object$path, object$manifest$learner.root)
  loaded <- .fast.load.named(
    name = object$manifest$supervised$model.name,
    root = learner.root,
    label = "supervised forest",
    strict = FALSE
  )
  list(
    model = loaded$model,
    loaded.from.disk = !is.null(loaded$model),
    enabled = TRUE,
    error = loaded$error
  )
}
