impute.learn.rfsrc <- function(formula, data,
                               ntree = 100, nodesize = 1, nsplit = 10,
                               nimpute = 2, fast = FALSE, blocks,
                               mf.q, max.iter = 10, eps = 0.01,
                               ytry = NULL, always.use = NULL, verbose = TRUE,
                               ...,
                               full.sweep.options = list(ntree = 100, nsplit = 10),
                               target.mode = c("missing.only", "all"),
                               deployment.xvars = NULL,
                               anonymous = TRUE,
                               learner.prefix = "impute.learner.",
                               learner.root = "learners",
                               out.dir = NULL,
                               wipe = TRUE,
                               keep.models = is.null(out.dir),
                               keep.ximp = FALSE,
                               save.on.fit = !is.null(out.dir)) {
  target.mode <- match.arg(target.mode)
  persist.on.fit <- !is.null(out.dir) && isTRUE(save.on.fit)
  if (isTRUE(save.on.fit) && is.null(out.dir)) {
    stop("'save.on.fit = TRUE' requires 'out.dir' to be supplied.",
         call. = FALSE)
  }
  if (!isTRUE(keep.models) && !persist.on.fit) {
    stop("At least one trained-learner storage mode must be enabled. ",
         "Use keep.models = TRUE or supply out.dir with save.on.fit = TRUE.",
         call. = FALSE)
  }
  if (persist.on.fit) {
    .check.fst()
  }
  ## make sure data is a data frame
  if (missing(data)) {
    stop("'data' is missing.", call. = FALSE)
  }
  if (is.character(data) && length(data) == 1L && is.null(dim(data))) {
    stop("'data' must be a data frame-like object, not a character string. ",
         "Did you mean data = ", data, " rather than data = \"", data, "\"?",
         call. = FALSE)
  }
  if (is.atomic(data) && is.null(dim(data)) && !is.list(data)) {
    stop("'data' must be a data frame, matrix, or other tabular object. ",
         "A bare vector is not allowed.",
         call. = FALSE)
  }
  train.input <- .normalize.training.data(data)
  dropped <- .drop.all.na.train(train.input)
  train.data <- dropped$data
  if (nrow(train.data) == 0L || ncol(train.data) == 0L) {
    stop("No usable rows or columns remain after removing all-NA rows or columns.",
         call. = FALSE)
  }
  which.na <- is.na(train.data)
  has.missing <- any(which.na)
  schema <- .build.schema(train.data)
  fs <- .parse.full.sweep.options(full.sweep.options)
  if (!has.missing) {
    if (identical(target.mode, "missing.only")) {
      stop("Training data have no missing values after preprocessing. ",
           "Use target.mode = \"all\" to fit predictive imputation learners ",
           "from complete training data.",
           call. = FALSE)
    }
    train.imputation <- "none; training data are complete"
    fit.seconds <- 0
    ximp <- as.data.frame(train.data, stringsAsFactors = FALSE)
    .msg("Training data are complete. Skipping training-time imputation and fitting ",
         "full-sweep learners directly.", verbose = verbose)
  }
  else {
    .msg("Starting training-time imputation...", verbose = verbose)
    fit.start <- proc.time()[[3]]
    impute.args <- c(
      list(
        data = train.data,
        ntree = ntree,
        nodesize = nodesize,
        nsplit = nsplit,
        nimpute = nimpute,
        fast = fast,
        max.iter = max.iter,
        eps = eps,
        ytry = ytry,
        always.use = always.use,
        verbose = verbose,
        full.sweep = FALSE
      ),
      list(...)
    )
    if (!missing(formula)) impute.args$formula <- formula
    if (!missing(blocks)) impute.args$blocks <- blocks
    mf.q.missing <- missing(mf.q)
    if (!mf.q.missing) impute.args$mf.q <- mf.q
    train.imputation <- if (mf.q.missing) {
      "selected by impute.rfsrc defaults"
    } else if (length(mf.q) == 1L && !is.na(mf.q) && identical(as.integer(mf.q), 1L)) {
      "missForest"
    } else {
      "multivariate missForest"
    }
    ximp <- tryCatch(
      do.call(impute.rfsrc, impute.args),
      error = function(e) e
    )
    if (inherits(ximp, "error")) {
      stop("Training-time imputation failed: ", conditionMessage(ximp),
           call. = FALSE)
    }
    fit.seconds <- proc.time()[[3]] - fit.start
    .msg("Training-time imputation finished in ", format(fit.seconds, digits = 4),
         " seconds.", verbose = verbose)
    ximp <- as.data.frame(ximp, stringsAsFactors = FALSE)
  }
  init <- .compute.init(ximp, schema)
  scale <- .compute.scale(ximp, schema)
  targets <- .resolve.targets(which.na, target.mode = target.mode)
  predictor.map <- .resolve.predictor.map(targets, names(train.data), deployment.xvars)
  bad.targets <- targets[lengths(predictor.map[targets]) == 0L]
  if (length(bad.targets) > 0L) {
    stop("Some targets have no deployment predictors: ",
         paste(bad.targets, collapse = ", "), call. = FALSE)
  }
  miss.frac <- colMeans(which.na)
  sweep.order <- targets[order(miss.frac[targets], decreasing = FALSE)]
  if (persist.on.fit) {
    if (isTRUE(wipe)) .safe.unlink.dir(out.dir)
    .safe.dir.create(out.dir)
    .safe.dir.create(file.path(out.dir, learner.root))
  }
  engine <- if (isTRUE(anonymous)) rfsrc.anonymous else rfsrc
  models <- setNames(vector("list", length(targets)), targets)
  learners <- setNames(vector("list", length(targets)), targets)
  .msg("Training final-sweep learner bank...", verbose = verbose)
  sweep.start <- proc.time()[[3]]
  for (i in seq_along(sweep.order)) {
    yname <- sweep.order[[i]]
    trn <- which(!which.na[, yname])
    tst <- which(which.na[, yname])
    xvars <- predictor.map[[yname]]
    learner.name <- .make.learner.name(i, yname, prefix = learner.prefix)
    learners[[yname]] <- list(
      learner.name = learner.name,
      predictors = xvars,
      n.obs = length(trn),
      n.missing.train = length(tst),
      status = "pending",
      error = NULL,
      family = NA_character_
    )
    .msg("  [", i, "/", length(sweep.order), "] target = `", yname,
         "`  predictors = ", length(xvars), "  observed rows = ", length(trn),
         verbose = verbose)
    if (length(trn) == 0L) {
      learners[[yname]]$status <- "skipped.all.missing"
      next
    }
    yy <- ximp[trn, yname]
    xtrain <- ximp[trn, xvars, drop = FALSE]
    response.name <- .make.response.name(names(xtrain))
    fit.df <- data.frame(xtrain, check.names = FALSE)
    fit.df[[response.name]] <- yy
    fit.df <- fit.df[, c(response.name, xvars), drop = FALSE]
    fit.args <- c(
      list(
        formula = stats::as.formula(paste0("`", response.name, "` ~ .")),
        data = fit.df,
        ntree = fs$ntree,
        nodesize = fs$nodesize,
        nsplit = fs$nsplit,
        perf.type = fs$dots$perf.type %||% "none",
        fast = fast
      ),
      fs$dots[names(fs$dots) != "perf.type"]
    )
    grow <- tryCatch(
      do.call(engine, fit.args),
      error = function(e) e
    )
    if (inherits(grow, "error")) {
      learners[[yname]]$status <- "error"
      learners[[yname]]$error <- conditionMessage(grow)
      .msg("      fit failed for `", yname, "`: ", conditionMessage(grow),
           verbose = verbose)
      next
    }
    learners[[yname]]$status <- "ok"
    learners[[yname]]$family <- grow$family
    if (isTRUE(keep.models)) {
      models[[yname]] <- grow
    }
    if (persist.on.fit) {
      learner.path <- file.path(out.dir, learner.root, learner.name)
      .msg("      saving learner to ", learner.path, verbose = verbose)
      fast.save(grow, learner.path, testing = FALSE)
    }
    if (!isTRUE(keep.models)) {
      rm(grow)
      gc()
    }
  }
  sweep.seconds <- proc.time()[[3]] - sweep.start
  .msg("Final-sweep learner bank finished in ", format(sweep.seconds, digits = 4),
       " seconds.", verbose = verbose)
  manifest <- list(
    version = "1.1",
    created.at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    formula = if (missing(formula)) NULL else paste(deparse(formula), collapse = ""),
    formula.scope = "initial imputation stage only",
    train.imputation = train.imputation,
    columns = names(train.data),
    schema = schema,
    init = init,
    scale = scale,
    targets = targets,
    sweep.order = sweep.order,
    predictor.map = predictor.map,
    deployment.xvars = deployment.xvars,
    learners = learners,
    learner.root = learner.root,
    learner.prefix = learner.prefix,
    target.mode = target.mode,
    anonymous = anonymous,
    fast = fast,
    full.sweep.options = full.sweep.options,
    train.missing.count = colSums(which.na),
    train.missing.frac = colMeans(which.na),
    n.train = nrow(train.data),
    p.train = ncol(train.data),
    dropped.all.na.rows = sum(!dropped$keep.rows),
    dropped.all.na.cols = names(train.input)[!dropped$keep.cols],
    fit.seconds = fit.seconds,
    sweep.seconds = sweep.seconds,
    call = match.call()
  )
  object <- list(
    manifest = manifest,
    models = if (isTRUE(keep.models)) models else setNames(vector("list", length(targets)), targets),
    ximp = if (isTRUE(keep.ximp)) ximp else NULL,
    path = if (persist.on.fit) normalizePath(out.dir, mustWork = FALSE) else NULL
  )
  class(object) <- c("impute.learn.rfsrc", "impute.learn")
  if (persist.on.fit) {
    saveRDS(manifest, file.path(out.dir, "manifest.rds"))
    .msg("Wrote manifest: ", file.path(out.dir, "manifest.rds"), verbose = verbose)
  }
  object
}
impute.learn <- impute.learn.rfsrc
save.impute.learn.rfsrc <- function(object, path, wipe = TRUE, verbose = TRUE) {
  if (!inherits(object, "impute.learn.rfsrc")) {
    stop("'object' must inherit from class 'impute.learn.rfsrc'.", call. = FALSE)
  }
  .check.fst()
  learner.root <- object$manifest$learner.root %||% "learners"
  source.path <- if (is.null(object$path)) NULL else normalizePath(object$path, mustWork = FALSE)
  dest.path <- normalizePath(path, mustWork = FALSE)
  if (!is.null(source.path) && identical(source.path, dest.path) && isTRUE(wipe)) {
    .msg("Save path matches the existing imputer path; leaving directory in place.",
         verbose = verbose)
    wipe <- FALSE
  }
  if (isTRUE(wipe)) .safe.unlink.dir(path)
  .safe.dir.create(path)
  .safe.dir.create(file.path(path, learner.root))
  saveRDS(object$manifest, file.path(path, "manifest.rds"))
  .msg("Saved manifest to ", file.path(path, "manifest.rds"), verbose = verbose)
  source.root <- if (is.null(source.path)) NULL else file.path(source.path, learner.root)
  for (target in object$manifest$targets) {
    info <- object$manifest$learners[[target]]
    if (is.null(info) || !identical(info$status, "ok")) next
    mdl <- object$models[[target]]
    if (is.null(mdl)) {
      if (is.null(source.root)) {
        stop("Learner for `", target, "` is not available in memory and no saved ",
             "learner path is attached to 'object'.",
             call. = FALSE)
      }
      .msg("Loading learner for `", target, "` from attached path before saving.",
           verbose = verbose)
      mdl <- .fast.load.learner(target, info, source.root, strict = TRUE)
    }
    learner.path <- file.path(path, learner.root, info$learner.name)
    .msg("Saving learner for `", target, "` to ", learner.path, verbose = verbose)
    fast.save(mdl, learner.path, testing = FALSE)
    if (is.null(object$models[[target]])) {
      rm(mdl)
      gc()
    }
  }
  invisible(load.impute.learn.rfsrc(path, lazy = TRUE, verbose = FALSE))
}
save.impute.learn <- save.impute.learn.rfsrc
load.impute.learn.rfsrc <- function(path, targets = NULL, lazy = TRUE, verbose = TRUE) {
  .check.fst()
  manifest.path <- file.path(path, "manifest.rds")
  if (!file.exists(manifest.path)) {
    stop("Manifest not found: ", manifest.path, call. = FALSE)
  }
  manifest <- readRDS(manifest.path)
  all.targets <- manifest$targets
  if (!is.null(targets)) {
    bad.targets <- setdiff(targets, all.targets)
    if (length(bad.targets) > 0L) {
      warning("Ignoring unknown targets: ", paste(bad.targets, collapse = ", "),
              call. = FALSE)
    }
  }
  use.targets <- if (is.null(targets)) all.targets else intersect(all.targets, targets)
  if (length(use.targets) == 0L) {
    stop("No requested targets were found in manifest.", call. = FALSE)
  }
  manifest$targets <- use.targets
  manifest$sweep.order <- manifest$sweep.order[manifest$sweep.order %in% use.targets]
  manifest$predictor.map <- manifest$predictor.map[use.targets]
  manifest$learners <- manifest$learners[use.targets]
  models <- setNames(vector("list", length(use.targets)), use.targets)
  object <- list(
    manifest = manifest,
    models = models,
    ximp = NULL,
    path = normalizePath(path, mustWork = TRUE)
  )
  class(object) <- c("impute.learn.rfsrc", "impute.learn")
  if (!isTRUE(lazy)) {
    learner.root <- file.path(path, manifest$learner.root)
    .msg("Loading learner bank into memory...", verbose = verbose)
    for (target in use.targets) {
      info <- manifest$learners[[target]]
      if (!identical(info$status, "ok")) next
      .msg("  loading `", target, "`", verbose = verbose)
      object$models[[target]] <- .fast.load.learner(target, info, learner.root,
                                                    strict = TRUE)
    }
  }
  object
}
load.impute.learn <- load.impute.learn.rfsrc
predict.impute.learn.rfsrc <- function(object, newdata,
                                       max.predict.iter = 3L,
                                       eps = 1e-3,
                                       targets = NULL,
                                       restore.integer = TRUE,
                                       cache.learners = c("session", "none", "all"),
                                       verbose = TRUE,
                                       ...) {
  if (!inherits(object, "impute.learn.rfsrc")) {
    stop("'object' must inherit from class 'impute.learn.rfsrc'.", call. = FALSE)
  }
  harmonized <- .harmonize.newdata(newdata, object$manifest, verbose = verbose)
  data <- harmonized$data
  if (!is.null(targets)) {
    bad.targets <- setdiff(targets, object$manifest$targets)
    if (length(bad.targets) > 0L) {
      warning("Ignoring unknown prediction targets: ",
              paste(bad.targets, collapse = ", "),
              call. = FALSE)
    }
  }
  use.targets <- if (is.null(targets)) object$manifest$targets else {
    intersect(object$manifest$targets, targets)
  }
  if (length(use.targets) == 0L) {
    stop("No valid targets requested for prediction-time imputation.", call. = FALSE)
  }
  cache.learners <- match.arg(cache.learners)
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
  any.target.missing <- length(sweep.order) > 0L &&
    any(as.matrix(original.missing[, sweep.order, drop = FALSE]))
  if (isTRUE(any.target.missing) && identical(cache.learners, "all")) {
    .msg("Preloading learner bank for prediction...", verbose = verbose)
    for (target in sweep.order) {
      info <- object$manifest$learners[[target]]
      if (!identical(info$status, "ok")) next
      mdl.info <- .predict.get.model(object, target, cache.env = cache.env)
      if (isTRUE(mdl.info$loaded.from.disk)) {
        disk.load.targets <- c(disk.load.targets, target)
      }
      if (is.null(mdl.info$model) && !is.null(mdl.info$error)) {
        record.issue(target, mdl.info$error)
      }
    }
  }
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
        ## conform step: makes sure factors are treated properly
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
    .msg("No missing values were found among prediction-time targets; ",
         "skipping iterative sweep.", verbose = verbose)
  }
  data <- .restore.schema(data, object$manifest$schema,
                          restore.integer = restore.integer)
  target.issues <- target.issues[lengths(target.issues) > 0L]
  attr(data, "impute.learn.info") <- list(
    n.passes = length(pass.history),
    pass.diff = pass.history,
    targets = use.targets,
    added.columns = harmonized$added.cols,
    dropped.extra.columns = harmonized$extra.cols,
    unseen.levels = harmonized$unseen.levels,
    cache.learners = cache.learners,
    n.disk.loads = length(disk.load.targets),
    disk.load.targets = unique(disk.load.targets),
    target.issues = target.issues
  )
  data
}
predict.impute.learn <- predict.impute.learn.rfsrc
print.impute.learn.rfsrc <- function(x, ...) {
  cat("Predictive imputer (randomForestSRC)\n")
  cat("  version:       ", x$manifest$version, "\n", sep = "")
  cat("  imputation:    ", x$manifest$train.imputation %||% "<unknown>", "\n", sep = "")
  cat("  training rows: ", x$manifest$n.train, "\n", sep = "")
  cat("  training cols: ", x$manifest$p.train, "\n", sep = "")
  cat("  targets:       ", length(x$manifest$targets), "\n", sep = "")
  cat("  learner root:  ", x$manifest$learner.root, "\n", sep = "")
  cat("  path:          ", x$path %||% "<memory>", "\n", sep = "")
  invisible(x)
}
print.impute.learn <- print.impute.learn.rfsrc
