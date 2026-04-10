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
      newdata[[nm]] <- .as.numeric.from.schema(x, sc)
    }
  }
  list(
    data = newdata,
    added.cols = added.cols,
    extra.cols = extra.cols,
    unseen.levels = unseen.levels,
    unseen.mask = as.data.frame(unseen.mask, stringsAsFactors = FALSE),
    unseen.rows = if (nrow(newdata) == 0L) logical(0) else rowSums(unseen.mask) > 0L
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
.ood.delta.regression <- function(observed, pred, target.schema = NULL) {
  if (is.null(pred) || is.null(pred$predicted)) {
    return(rep(NA_real_, length(observed)))
  }
  out <- pred$predicted
  if (isTRUE(target.schema$is.integer)) {
    out <- as.integer(round(out))
  }
  obs <- .as.numeric.from.schema(observed, target.schema)
  est <- .as.numeric.from.schema(out, target.schema)
  abs(obs - est)
}
.ood.delta.classification <- function(observed, pred, target.schema,
                                      prob.floor = 1e-12) {
  y <- as.character(observed)
  out <- rep(NA_real_, length(y))
  missing.y <- is.na(y)
  prob <- .extract.class.prob(pred, target.schema)
  if (!is.null(prob) && nrow(prob) == length(y)) {
    lev <- colnames(prob)
    pos <- match(y, lev)
    ok <- !missing.y & !is.na(pos)
    if (any(ok)) {
      p.true <- prob[cbind(which(ok), pos[ok])]
      p.true[!is.finite(p.true)] <- 0
      out[ok] <- -log(pmax(p.true, prob.floor))
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
  cls <- as.character(cls)
  ok <- !missing.y & !is.na(cls)
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
.compute.ood.delta <- function(observed, pred, target.schema) {
  if (isTRUE(target.schema$is.factor)) {
    .ood.delta.classification(observed, pred, target.schema = target.schema)
  }
  else {
    .ood.delta.regression(observed, pred, target.schema = target.schema)
  }
}
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
.eval.ood.reference <- function(x, ref) {
  out <- rep(NA_real_, length(x))
  if (length(out) == 0L || is.null(ref) || is.null(ref$quantiles)) {
    return(out)
  }
  ok <- is.finite(x)
  if (!any(ok)) {
    return(out)
  }
  q <- as.numeric(ref$quantiles)
  p <- as.numeric(ref$probs)
  if (length(q) == 0L || length(p) == 0L) {
    return(out)
  }
  if (length(q) == 1L || all(q == q[1L])) {
    out[ok] <- ifelse(x[ok] <= q[1L], p[1L], p[length(p)])
    return(out)
  }
  idx <- findInterval(x[ok], q, rightmost.closed = TRUE, all.inside = TRUE)
  idx <- pmax(1L, pmin(length(p), idx))
  out[ok] <- p[idx]
  out
}
.aggregate.ood.row <- function(mat, weight = NULL) {
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
  out <- rep(NA_real_, nrow(mat))
  for (i in seq_len(nrow(mat))) {
    ok <- is.finite(mat[i, ]) & is.finite(weight) & weight > 0
    if (!any(ok)) next
    out[i] <- sum(weight[ok] * mat[i, ok]) / sum(weight[ok])
  }
  out
}
.resolve.ood.weight <- function(targets, weight = NULL, default = NULL) {
  if (is.null(weight)) {
    weight <- default
  }
  if (is.null(weight)) {
    weight <- setNames(rep(1, length(targets)), targets)
  }
  if (is.null(names(weight))) {
    if (length(weight) != length(targets)) {
      stop("'weight' must have one entry per target.", call. = FALSE)
    }
    weight <- setNames(as.numeric(weight), targets)
  }
  else {
    missing.targets <- setdiff(targets, names(weight))
    if (length(missing.targets) > 0L) {
      stop("Weights are missing for targets: ",
           paste(missing.targets, collapse = ", "),
           call. = FALSE)
    }
    weight <- as.numeric(weight[targets])
    names(weight) <- targets
  }
  if (any(!is.finite(weight) | weight < 0)) {
    stop("'weight' must contain finite nonnegative values.", call. = FALSE)
  }
  if (!any(weight > 0)) {
    stop("At least one target weight must be positive.", call. = FALSE)
  }
  weight
}
.same.ood.weight <- function(targets, x, y,
                             tolerance = sqrt(.Machine$double.eps)) {
  if (is.null(x) || is.null(y)) {
    return(FALSE)
  }
  wx <- tryCatch(.resolve.ood.weight(targets, x), error = function(e) NULL)
  wy <- tryCatch(.resolve.ood.weight(targets, y), error = function(e) NULL)
  if (is.null(wx) || is.null(wy)) {
    return(FALSE)
  }
  isTRUE(all.equal(unname(wx), unname(wy), tolerance = tolerance,
                   check.attributes = FALSE))
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
  original.missing <- as.data.frame(is.na(data[, use.targets, drop = FALSE]))
  names(original.missing) <- use.targets
  data <- .apply.init(data, object$manifest$init, object$manifest$schema)
  ## makes working copy of `data` look more like the training x-schema before iterative sweep
  data <- .restore.schema(data, object$manifest$schema, restore.integer = TRUE)
  pass.history <- numeric(0)
  sweep.order <- object$manifest$sweep.order
  sweep.order <- sweep.order[sweep.order %in% use.targets]
  cache.env <- if (identical(cache.learners, "none")) NULL else new.env(parent = emptyenv())
  disk.load.targets <- character(0)
  target.issues <- setNames(vector("list", length(use.targets)), use.targets)
  record.issue <- function(target, message) {
    current <- target.issues[[target]]
    if (is.null(current)) current <- character(0)
    if (!(message %in% current)) {
      target.issues[[target]] <<- c(current, message)
    }
    invisible(NULL)
  }
  if (identical(cache.learners, "all")) {
    .msg("Preloading learner bank...", verbose = verbose)
    for (target in use.targets) {
      info <- object$manifest$learners[[target]]
      if (!identical(info$status, "ok")) next
      mdl.info <- .predict.get.model(object, target, cache.env = cache.env)
      if (isTRUE(mdl.info$loaded.from.disk)) {
        disk.load.targets <- c(disk.load.targets, target)
      }
      if (is.null(mdl.info$model) && !is.null(mdl.info$error)) {
        record_issue <- mdl.info$error
        record.issue(target, record_issue)
      }
    }
  }
  any.target.missing <- length(sweep.order) > 0L &&
    any(as.matrix(original.missing[, sweep.order, drop = FALSE]))
  if (isTRUE(any.target.missing)) {
    .msg("Starting prediction-time sweep...", verbose = verbose)
    for (iter in seq_len(max.predict.iter)) {
      old.data <- data
      .msg("  prediction pass ", iter, "/", max.predict.iter, verbose = verbose)
      for (target in sweep.order) {
        miss.idx <- which(original.missing[[target]])
        if (length(miss.idx) == 0L) next
        info <- object$manifest$learners[[target]]
        if (!identical(info$status, "ok")) {
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
        }
        if (is.null(mdl)) {
          msg <- mdl.info$error %||% "learner could not be loaded"
          record.issue(target, msg)
          .msg("    skipping `", target, "` (", msg, ")", verbose = verbose)
          next
        }
        xvars <- object$manifest$predictor.map[[target]]
        pred.df <- data[miss.idx, xvars, drop = FALSE]
        pred.df <- .conform.x.to.forest(pred.df, mdl)
        pred <- tryCatch(
          predict(mdl, pred.df),
          error = function(e) e
        )
        if (inherits(pred, "error")) {
          msg <- paste0("Prediction failed: ", conditionMessage(pred))
          record.issue(target, msg)
          .msg("    prediction failed for `", target, "`: ", conditionMessage(pred),
               verbose = verbose)
          next
        }
        values <- .extract.prediction(pred, info$family, object$manifest$schema[[target]])
        data[miss.idx, target] <- values
        if (identical(cache.learners, "none") && is.null(object$models[[target]])) {
          rm(mdl)
          gc()
        }
      }
      diff.err <- .compute.pass.diff(old.data, data, original.missing,
                                     object$manifest$schema,
                                     object$manifest$scale,
                                     sweep.order)
      pass.history <- c(pass.history, diff.err)
      .msg("    pass diff = ", format(diff.err, digits = 4), verbose = verbose)
      if (is.finite(diff.err) && diff.err < eps) {
        .msg("    convergence criterion met; stopping early.", verbose = verbose)
        break
      }
    }
  }
  else {
    .msg("No missing values were found among requested targets; ",
         "skipping iterative sweep.", verbose = verbose)
  }
  data <- .restore.schema(data, object$manifest$schema,
                          restore.integer = restore.integer)
  target.issues <- target.issues[lengths(target.issues) > 0L]
  list(
    data = data,
    use.targets = use.targets,
    harmonized = harmonized,
    cache.learners = cache.learners,
    cache.env = cache.env,
    info = list(
      n.passes = length(pass.history),
      pass.diff = pass.history,
      targets = use.targets,
      added.columns = harmonized$added.cols,
      dropped.extra.columns = harmonized$extra.cols,
      unseen.levels = harmonized$unseen.levels,
      unseen.rows = harmonized$unseen.rows,
      cache.learners = cache.learners,
      n.disk.loads = length(disk.load.targets),
      disk.load.targets = unique(disk.load.targets),
      target.issues = target.issues
    )
  )
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
