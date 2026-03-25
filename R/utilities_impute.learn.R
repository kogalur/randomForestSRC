`%||%` <- function(x, y) if (is.null(x)) y else x
.check.fst <- function() {
  if (!requireNamespace("fst", quietly = TRUE)) {
    stop("Package 'fst' is required for save/load support.", call. = FALSE)
  }
  invisible(TRUE)
}
.timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}
.msg <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) {
    cat(sprintf("[%s] ", .timestamp()), ..., "\n", sep = "")
    flush.console()
  }
}
.safe.dir.create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}
.safe.unlink.dir <- function(path) {
  if (dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
  }
  invisible(path)
}
.safe.name <- function(x) {
  x <- gsub("[^[:alnum:].]+", ".", x)
  x <- gsub("\\.+", ".", x)
  x <- gsub("^\\.+|\\.+$", "", x)
  if (!nzchar(x)) x <- "var"
  x
}
.is.real.valued <- function(x) {
  (is.double(x) || is.integer(x)) &&
    !inherits(x, c("Date", "POSIXt", "difftime"))
}
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
.normalize.training.data <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  for (nm in names(data)) {
    data[[nm]] <- .coerce.supported.column(data[[nm]], nm)
  }
  data
}
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
.build.schema <- function(data) {
  out <- lapply(names(data), function(nm) {
    x <- data[[nm]]
    list(
      class = class(x),
      is.factor = is.factor(x),
      ordered = is.ordered(x),
      levels = if (is.factor(x)) levels(x) else NULL,
      is.integer = is.integer(x),
      is.numeric = .is.real.valued(x)
    )
  })
  names(out) <- names(data)
  out
}
.as.numeric.safe <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  suppressWarnings(as.numeric(x))
}
.as.numeric.from.schema <- function(x, schema = NULL) {
  .as.numeric.safe(x)
}
.mode.value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA)
  tb <- sort(table(x), decreasing = TRUE)
  names(tb)[1]
}
.numeric.init <- function(x) {
  x <- .as.numeric.safe(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}
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
.resolve.predictor.map <- function(targets, all.names, deployment.xvars = NULL) {
  if (is.null(deployment.xvars)) {
    out <- lapply(targets, function(y) setdiff(all.names, y))
    names(out) <- targets
    return(out)
  }
  if (is.character(deployment.xvars)) {
    xvars <- intersect(all.names, deployment.xvars)
    out <- lapply(targets, function(y) setdiff(xvars, y))
    names(out) <- targets
    return(out)
  }
  if (is.list(deployment.xvars)) {
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
    for (nm in miss.cols) x[[nm]] <- NA
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
.parse.full.sweep.options <- function(full.sweep.options = NULL) {
  fs <- full.sweep.options %||% list()
  allowed <- c(
    "mtry", "splitrule", "bootstrap", "sampsize", "samptype",
    "perf.type", "rfq", "save.memory", "importance", "proximity"
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
.make.learner.name <- function(i, target, prefix = "impute.learner.") {
  sprintf("%s%03d.%s", prefix, i, .safe.name(target))
}
.make.response.name <- function(existing.names = character()) {
  nm <- ".impute.learn.response."
  while (nm %in% existing.names) {
    nm <- paste0(nm, ".")
  }
  nm
}
.coerce.factor.levels <- function(x, levels, ordered = FALSE) {
  x <- as.character(x)
  x[!(x %in% levels)] <- NA_character_
  factor(x, levels = levels, ordered = ordered)
}
.harmonize.newdata <- function(newdata, manifest, verbose = TRUE) {
  newdata <- as.data.frame(newdata, stringsAsFactors = FALSE)
  required.cols <- manifest$columns
  added.cols <- setdiff(required.cols, names(newdata))
  extra.cols <- setdiff(names(newdata), required.cols)
  if (length(added.cols) > 0L) {
    for (nm in added.cols) {
      newdata[[nm]] <- NA
    }
    .msg("Added missing columns to newdata: ", paste(added.cols, collapse = ", "),
         verbose = verbose)
  }
  if (length(extra.cols) > 0L) {
    .msg("Dropping extra columns from newdata: ", paste(extra.cols, collapse = ", "),
         verbose = verbose)
  }
  newdata <- newdata[, required.cols, drop = FALSE]
  unseen.levels <- vector("list", length(required.cols))
  names(unseen.levels) <- required.cols
  for (nm in required.cols) {
    sc <- manifest$schema[[nm]]
    x <- newdata[[nm]]
    if (isTRUE(sc$is.factor)) {
      xchr <- as.character(x)
      unseen <- unique(stats::na.omit(xchr[!(xchr %in% sc$levels)]))
      if (length(unseen) > 0L) {
        unseen.levels[[nm]] <- unseen
      }
      newdata[[nm]] <- .coerce.factor.levels(xchr, sc$levels, ordered = sc$ordered)
    }
    else {
      newdata[[nm]] <- .as.numeric.from.schema(x, sc)
    }
  }
  list(
    data = newdata,
    added.cols = added.cols,
    extra.cols = extra.cols,
    unseen.levels = unseen.levels
  )
}
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
.restore.schema <- function(data, schema, restore.integer = TRUE) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  for (nm in names(data)) {
    sc <- schema[[nm]]
    if (isTRUE(sc$is.factor)) {
      data[[nm]] <- .coerce.factor.levels(data[[nm]], sc$levels, ordered = sc$ordered)
    }
    else if (isTRUE(sc$is.integer)) {
      x <- .as.numeric.from.schema(data[[nm]], sc)
      if (isTRUE(restore.integer)) {
        xi <- as.integer(round(x))
        xi[is.na(x)] <- NA_integer_
        data[[nm]] <- xi
      }
      else {
        data[[nm]] <- x
      }
    }
    else {
      data[[nm]] <- .as.numeric.from.schema(data[[nm]], sc)
    }
  }
  data
}
.extract.prediction <- function(pred, family, target.schema) {
  if (is.null(pred)) return(NULL)
  if (identical(family, "regr")) {
    out <- pred$predicted
    if (isTRUE(target.schema$is.integer)) {
      out <- as.integer(round(out))
    }
    return(out)
  }
  out <- pred$class
  if (is.null(out) && !is.null(pred$predicted)) {
    out <- pred$predicted
  }
  if (isTRUE(target.schema$is.factor)) {
    out <- .coerce.factor.levels(out, target.schema$levels,
                                 ordered = target.schema$ordered)
  }
  out
}
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
.predict.get.model <- function(object, target, cache.env = NULL) {
  if (!inherits(object, "impute.learn.rfsrc")) {
    stop("'object' must inherit from class 'impute.learn.rfsrc'.", call. = FALSE)
  }
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
