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
                               save.on.fit = !is.null(out.dir),
                               save.ood = TRUE) {
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
  ood.delta <- if (isTRUE(save.ood)) setNames(vector("list", length(targets)), targets) else NULL
  ood.issues <- if (isTRUE(save.ood)) setNames(vector("list", length(targets)), targets) else NULL
  record.ood.issue <- function(target, message) {
    if (!isTRUE(save.ood)) return(invisible(NULL))
    current <- ood.issues[[target]]
    if (is.null(current)) current <- character(0)
    if (!(message %in% current)) {
      ood.issues[[target]] <<- c(current, message)
    }
    invisible(NULL)
  }
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
      record.ood.issue(yname, "No observed training rows were available for OOD reference.")
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
      record.ood.issue(yname, paste0("Learner fit failed: ", conditionMessage(grow)))
      .msg("      fit failed for `", yname, "`: ", conditionMessage(grow),
           verbose = verbose)
      next
    }
    learners[[yname]]$status <- "ok"
    learners[[yname]]$family <- grow$family
    if (isTRUE(save.ood)) {
      pred.oob <- list(
        predicted = grow$predicted.oob,
        class = grow$class.oob
      )
      delta.oob <- tryCatch(
        .compute.ood.delta(yy, pred.oob, schema[[yname]]),
        error = function(e) e
      )
      if (inherits(delta.oob, "error")) {
        record.ood.issue(yname, paste0("Failed to compute OOD reference: ",
                                       conditionMessage(delta.oob)))
      }
      else if (length(delta.oob) != length(trn)) {
        record.ood.issue(yname, "OOB reference length did not match observed rows.")
      }
      else {
        tmp <- rep(NA_real_, nrow(train.data))
        tmp[trn] <- as.numeric(delta.oob)
        ood.delta[[yname]] <- tmp
      }
    }
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
  ood <- NULL
  if (isTRUE(save.ood)) {
    target.reference <- list()
    valid.ood.targets <- character(0)
    for (yname in targets) {
      delta.y <- ood.delta[[yname]]
      if (is.null(delta.y) || !any(is.finite(delta.y))) {
        if (is.null(ood.issues[[yname]]) || length(ood.issues[[yname]]) == 0L) {
          record.ood.issue(yname, "No finite OOD reference values were available.")
        }
        next
      }
      target.reference[[yname]] <- .make.ood.reference(delta.y)
      valid.ood.targets <- c(valid.ood.targets, yname)
    }
    if (length(valid.ood.targets) > 0L) {
      target.score.train <- do.call(cbind, lapply(valid.ood.targets, function(yname) {
        .eval.ood.reference(ood.delta[[yname]], target.reference[[yname]])
      }))
      colnames(target.score.train) <- valid.ood.targets
      default.weight <- setNames(rep(1, length(valid.ood.targets)), valid.ood.targets)
      row.score.train <- .aggregate.ood.row(target.score.train, default.weight)
      row.reference <- .make.ood.reference(row.score.train)
    }
    else {
      default.weight <- setNames(numeric(0), character(0))
      row.reference <- NULL
    }
    ood <- list(
      version = "1.0",
      reference = "oob",
      aggregate = "weighted.mean",
      target.metric = c(
        numeric = "absolute.error",
        factor = "negative.log.probability",
        fallback = "misclass.or.rank.distance"
      ),
      targets = valid.ood.targets,
      weight = default.weight,
      target.reference = target.reference,
      row.reference = row.reference,
      issues = ood.issues
    )
  }
  manifest <- list(
    version = "1.2",
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
    save.ood = save.ood,
    ood = ood,
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
  cache.learners <- match.arg(cache.learners)
  prep <- .prepare.impute.learn.newdata(
    object = object,
    newdata = newdata,
    targets = targets,
    max.predict.iter = max.predict.iter,
    eps = eps,
    restore.integer = restore.integer,
    cache.learners = cache.learners,
    verbose = verbose
  )
  attr(prep$data, "impute.learn.info") <- prep$info
  prep$data
}
predict.impute.learn <- predict.impute.learn.rfsrc
impute.ood.rfsrc <- function(object, newdata,
                             targets = NULL,
                             max.predict.iter = 3L,
                             eps = 1e-3,
                             cache.learners = c("all", "session", "none"),
                             weight = NULL,
                             return.details = FALSE,
                             verbose = TRUE,
                             ...) {
  cache.learners <- match.arg(cache.learners)
  if (!inherits(object, "impute.learn.rfsrc")) {
    stop("'object' must inherit from class 'impute.learn.rfsrc'.", call. = FALSE)
  }
  ood.ref <- object$manifest$ood
  if (is.null(ood.ref) || length(ood.ref$targets %||% character(0)) == 0L) {
    stop("No saved OOD reference was found in 'object$manifest$ood'. ",
         "Refit with save.ood = TRUE to enable OOD scoring.",
         call. = FALSE)
  }
  if (!is.null(targets)) {
    bad.targets <- setdiff(targets, object$manifest$targets)
    if (length(bad.targets) > 0L) {
      warning("Ignoring unknown OOD targets: ",
              paste(bad.targets, collapse = ", "),
              call. = FALSE)
    }
  }
  score.targets <- if (is.null(targets)) object$manifest$targets else {
    intersect(object$manifest$targets, targets)
  }
  if (length(score.targets) == 0L) {
    stop("No valid targets requested for OOD scoring.", call. = FALSE)
  }
  prep <- .prepare.impute.learn.newdata(
    object = object,
    newdata = newdata,
    targets = NULL,
    max.predict.iter = max.predict.iter,
    eps = eps,
    restore.integer = TRUE,
    cache.learners = cache.learners,
    verbose = verbose
  )
  use.targets <- intersect(score.targets, ood.ref$targets %||% character(0))
  missing.ref.targets <- setdiff(score.targets, use.targets)
  if (length(missing.ref.targets) > 0L) {
    warning("Skipping targets without a saved OOD reference: ",
            paste(missing.ref.targets, collapse = ", "),
            call. = FALSE)
  }
  if (length(use.targets) == 0L) {
    stop("No requested targets have a saved OOD reference.", call. = FALSE)
  }
  weight <- .resolve.ood.weight(use.targets, weight, default = ood.ref$weight)
  completed.data <- prep$data
  n <- nrow(completed.data)
  target.delta <- matrix(NA_real_, nrow = n, ncol = length(use.targets),
                         dimnames = list(NULL, use.targets))
  target.score <- matrix(NA_real_, nrow = n, ncol = length(use.targets),
                         dimnames = list(NULL, use.targets))
  target.issues <- prep$info$target.issues
  if (is.null(target.issues)) {
    target.issues <- setNames(vector("list", length(use.targets)), use.targets)
  }
  for (nm in setdiff(use.targets, names(target.issues))) {
    target.issues[[nm]] <- character(0)
  }
  record.issue <- function(target, message) {
    current <- target.issues[[target]]
    if (is.null(current)) current <- character(0)
    if (!(message %in% current)) {
      target.issues[[target]] <<- c(current, message)
    }
    invisible(NULL)
  }
  disk.load.targets <- prep$info$disk.load.targets %||% character(0)
  unseen.mask <- prep$harmonized$unseen.mask
  unseen.rows <- prep$harmonized$unseen.rows
  .msg("Scoring OOD targets...", verbose = verbose)
  for (target in use.targets) {
    info <- object$manifest$learners[[target]]
    if (!identical(info$status, "ok")) {
      msg <- paste0("No trained learner is available (status = ",
                    info$status %||% "unknown", ").")
      record.issue(target, msg)
      next
    }
    mdl.info <- .predict.get.model(object, target, cache.env = prep$cache.env)
    mdl <- mdl.info$model
    if (isTRUE(mdl.info$loaded.from.disk)) {
      disk.load.targets <- c(disk.load.targets, target)
    }
    if (is.null(mdl)) {
      record.issue(target, mdl.info$error %||% "learner could not be loaded")
      next
    }
    xvars <- object$manifest$predictor.map[[target]]
    pred.df <- completed.data[, xvars, drop = FALSE]
    pred.df <- .conform.x.to.forest(pred.df, mdl)
    pred <- tryCatch(
      predict(mdl, pred.df),
      error = function(e) e
    )
    if (inherits(pred, "error")) {
      record.issue(target, paste0("Prediction failed: ", conditionMessage(pred)))
      if (identical(prep$cache.learners, "none") && is.null(object$models[[target]])) {
        rm(mdl)
        gc()
      }
      next
    }
    observed <- prep$harmonized$data[[target]]
    delta <- .compute.ood.delta(observed, pred,
                                object$manifest$schema[[target]])
    target.unseen <- if (target %in% names(unseen.mask)) {
      unseen.mask[[target]]
    } else {
      rep(FALSE, n)
    }
    if (length(target.unseen) > 0L && any(target.unseen)) {
      delta[target.unseen] <- Inf
    }
    score.j <- .eval.ood.reference(delta, ood.ref$target.reference[[target]])
    if (length(target.unseen) > 0L && any(target.unseen)) {
      score.j[target.unseen] <- 1
    }
    target.delta[, target] <- delta
    target.score[, target] <- score.j
    if (identical(prep$cache.learners, "none") && is.null(object$models[[target]])) {
      rm(mdl)
      gc()
    }
  }
  score <- .aggregate.ood.row(target.score[, use.targets, drop = FALSE], weight)
  weight.mask <- matrix(rep(weight > 0, each = n), nrow = n)
  targets.used <- rowSums(is.finite(target.score[, use.targets, drop = FALSE]) & weight.mask)
  score.percentile <- rep(NA_real_, n)
  row.reference.used <- FALSE
  row.reference.reason <- NULL
  same.targets <- setequal(use.targets, ood.ref$targets %||% character(0)) &&
    length(use.targets) == length(ood.ref$targets %||% character(0))
  same.weight <- .same.ood.weight(use.targets, weight, ood.ref$weight)
  if (!is.null(ood.ref$row.reference) && isTRUE(same.targets) && isTRUE(same.weight)) {
    score.percentile <- .eval.ood.reference(score, ood.ref$row.reference)
    row.reference.used <- TRUE
  }
  else {
    if (is.null(ood.ref$row.reference)) {
      row.reference.reason <- "No saved row-level OOD reference is available."
    }
    else if (!isTRUE(same.targets)) {
      row.reference.reason <- "Saved row-level OOD calibration requires scoring the original target set."
    }
    else if (!isTRUE(same.weight)) {
      row.reference.reason <- "Saved row-level OOD calibration requires the default target weights."
    }
  }
  if (length(unseen.rows) > 0L && any(unseen.rows)) {
    score[unseen.rows] <- 1
    if (isTRUE(row.reference.used)) {
      score.percentile[unseen.rows] <- 1
    }
  }
  target.issues <- target.issues[lengths(target.issues) > 0L]
  out <- list(
    score = score,
    score.percentile = score.percentile,
    targets.used = targets.used,
    target.score = if (isTRUE(return.details)) target.score[, use.targets, drop = FALSE] else NULL,
    target.delta = if (isTRUE(return.details)) target.delta[, use.targets, drop = FALSE] else NULL,
    info = list(
      targets = use.targets,
      weight = weight,
      added.columns = prep$info$added.columns,
      dropped.extra.columns = prep$info$dropped.extra.columns,
      unseen.levels = prep$info$unseen.levels,
      unseen.rows = unseen.rows,
      maxed.rows = unseen.rows,
      cache.learners = prep$cache.learners,
      n.disk.loads = length(unique(disk.load.targets)),
      disk.load.targets = unique(disk.load.targets),
      row.reference.used = row.reference.used,
      row.reference.reason = row.reference.reason,
      target.issues = target.issues
    )
  )
  if (isTRUE(return.details)) {
    out$info$unseen.mask <- unseen.mask
  }
  class(out) <- c("impute.ood.rfsrc", "impute.ood")
  out
}
impute.ood <- impute.ood.rfsrc
print.impute.learn.rfsrc <- function(x, ...) {
  cat("Predictive imputer (randomForestSRC)\n")
  cat("  version:       ", x$manifest$version, "\n", sep = "")
  cat("  imputation:    ", x$manifest$train.imputation %||% "<unknown>", "\n", sep = "")
  cat("  training rows: ", x$manifest$n.train, "\n", sep = "")
  cat("  training cols: ", x$manifest$p.train, "\n", sep = "")
  cat("  targets:       ", length(x$manifest$targets), "\n", sep = "")
  cat("  ood targets:   ", length(x$manifest$ood$targets %||% character(0)), "\n", sep = "")
  cat("  learner root:  ", x$manifest$learner.root, "\n", sep = "")
  cat("  path:          ", x$path %||% "<memory>", "\n", sep = "")
  invisible(x)
}
print.impute.learn <- print.impute.learn.rfsrc
