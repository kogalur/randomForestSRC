impute.learn.rfsrc <- function(formula, data,
                               ntree = 100, nodesize = 1, nsplit = 10,
                               nimpute = 2, fast = FALSE, blocks,
                               mf.q, max.iter = 10, eps = 0.01,
                               ytry = NULL, always.use = NULL, verbose = TRUE,
                               ...,
                               supervised.formula = NULL,
                               supervised.args = list(),
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
                               save.ood = TRUE,
                               weight = NULL) {
  target.mode.missing <- missing(target.mode)
  supervised.enabled <- !missing(supervised.formula) && !is.null(supervised.formula)
  formula.supplied <- !missing(formula)
  formula.used <- isTRUE(formula.supplied) && !isTRUE(supervised.enabled)
  if (isTRUE(formula.supplied) && isTRUE(supervised.enabled)) {
    warning("'formula' is ignored when 'supervised.formula' is supplied. ",
            "The raw predictor block is defined by the right-hand side of ",
            "'supervised.formula'.", call. = FALSE)
  }
  if (isTRUE(supervised.enabled) && isTRUE(target.mode.missing)) {
    target.mode <- "all"
  }
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
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (isTRUE(supervised.enabled)) {
    supervised.spec <- .parse.supervised.spec(supervised.formula, data)
    train.input <- .normalize.training.data(supervised.spec$x.data)
    response.data <- supervised.spec$y.data
    supervised.args <- .parse.supervised.args(supervised.args)
  }
  else {
    supervised.spec <- list(enabled = FALSE)
    train.input <- .normalize.training.data(data)
    response.data <- NULL
    supervised.args <- list()
  }
  dropped <- .drop.all.na.train(train.input)
  train.data <- dropped$data
  if (nrow(train.data) == 0L || ncol(train.data) == 0L) {
    stop("No usable rows or columns remain after removing all-NA rows or columns.",
         call. = FALSE)
  }
  if (isTRUE(supervised.enabled)) {
    response.data <- response.data[dropped$keep.rows, , drop = FALSE]
    response.schema <- .build.schema(response.data)
  }
  else {
    response.schema <- NULL
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
    if (isTRUE(formula.used)) {
      impute.args$formula <- formula
    }
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
  if (anyNA(ximp)) {
    stop("Training-time imputation did not return a fully completed raw predictor table.",
         call. = FALSE)
  }
  init <- .compute.init(ximp, schema)
  scale <- .compute.scale(ximp, schema)
  engine <- if (isTRUE(anonymous)) rfsrc.anonymous else rfsrc
  supervised.model <- NULL
  supervised.info <- list(enabled = FALSE)
  augmented.ximp <- ximp
  if (isTRUE(supervised.enabled)) {
    response.rows <- .supervised.usable.rows(response.data)
    if (length(response.rows) != nrow(train.data)) {
      stop("Internal supervised-response bookkeeping failed after training preprocessing.",
           call. = FALSE)
    }
    if (!any(response.rows)) {
      stop("No retained training rows have an observed supervised response.",
           call. = FALSE)
    }
    fit.supervised.formula <- .project.supervised.formula(
      supervised.formula,
      data = data,
      train.names = names(train.data),
      arg.name = "supervised.formula"
    )
    fit.sup.data <- data.frame(
      response.data[response.rows, , drop = FALSE],
      ximp[response.rows, , drop = FALSE],
      check.names = FALSE
    )
    supervised.fit.args <- list(
      formula = fit.supervised.formula,
      data = fit.sup.data,
      forest = TRUE,
      perf.type = "none",
      fast = fast
    )
    if (length(supervised.args) > 0L) {
      supervised.fit.args[names(supervised.args)] <- supervised.args
    }
    .msg("Fitting supervised forest...", verbose = verbose)
    supervised.fit <- tryCatch(
      do.call(rfsrc, supervised.fit.args),
      error = function(e) e
    )
    if (inherits(supervised.fit, "error")) {
      stop("Supervised forest fit failed: ", conditionMessage(supervised.fit),
           call. = FALSE)
    }
    supervised.meta <- .resolve.supervised.fit.meta(supervised.fit, fallback = supervised.spec)
    if (!setequal(supervised.meta$xvar.names, names(train.data))) {
      stop("The supervised forest returned predictor names that do not match the retained training predictors.",
           call. = FALSE)
    }
    aux.used <- .align.supervised.aux(
      .supervised.materialize.aux(
        supervised.fit,
        family = supervised.meta$family,
        oob = TRUE,
        response.schema = response.schema,
        yvar.names = supervised.meta$yvar.names
      )
    )
    if (nrow(aux.used) != sum(response.rows)) {
      stop("The supervised OOB prediction block produced ", nrow(aux.used),
           " row(s), but ", sum(response.rows),
           " training row(s) with observed supervised response were expected.",
           call. = FALSE)
    }
    aux.train <- as.data.frame(
      matrix(NA_real_, nrow = nrow(train.data), ncol = ncol(aux.used)),
      check.names = FALSE
    )
    names(aux.train) <- names(aux.used)
    aux.train[response.rows, ] <- aux.used
    response.omit <- which(!response.rows)
    if (length(response.omit) > 0L) {
      aux.omit.newdata <- .conform.x.to.forest(
        ximp[response.omit, , drop = FALSE],
        supervised.fit,
        ignore.levels = TRUE
      )
      aux.omit.pred <- tryCatch(
        predict(supervised.fit, aux.omit.newdata),
        error = function(e) e
      )
      if (inherits(aux.omit.pred, "error")) {
        stop("Failed to predict supervised auxiliary values for rows with missing response: ",
             conditionMessage(aux.omit.pred), call. = FALSE)
      }
      aux.omit <- .align.supervised.aux(
        .supervised.materialize.aux(
          aux.omit.pred,
          family = supervised.meta$family,
          oob = FALSE,
          response.schema = response.schema,
          yvar.names = supervised.meta$yvar.names
        ),
        expected.names = names(aux.train)
      )
      if (nrow(aux.omit) != length(response.omit)) {
        stop("The supervised prediction block for rows with missing response produced ",
             nrow(aux.omit), " row(s), but ", length(response.omit),
             " were expected.", call. = FALSE)
      }
      aux.train[response.omit, ] <- aux.omit
    }
    aux.invalid <- if (ncol(aux.train) > 0L) {
      rowSums(!is.finite(as.matrix(aux.train))) > 0L
    } else {
      rep(FALSE, nrow(aux.train))
    }
    if (any(aux.invalid)) {
      aux.fill.newdata <- .conform.x.to.forest(
        ximp[aux.invalid, , drop = FALSE],
        supervised.fit,
        ignore.levels = TRUE
      )
      aux.fill.pred <- tryCatch(
        predict(supervised.fit, aux.fill.newdata),
        error = function(e) e
      )
      if (!inherits(aux.fill.pred, "error")) {
        aux.fill <- .align.supervised.aux(
          .supervised.materialize.aux(
            aux.fill.pred,
            family = supervised.meta$family,
            oob = FALSE,
            response.schema = response.schema,
            yvar.names = supervised.meta$yvar.names
          ),
          expected.names = names(aux.train)
        )
        aux.train[aux.invalid, ] <- aux.fill
      }
    }
    if (length(intersect(names(aux.train), names(train.data))) > 0L) {
      stop("Supervised auxiliary names collide with raw predictor names: ",
           paste(intersect(names(aux.train), names(train.data)), collapse = ", "),
           call. = FALSE)
    }
    aux.schema <- .build.schema(aux.train)
    aux.init <- .compute.init(aux.train, aux.schema)
    aux.scale <- .compute.scale(aux.train, aux.schema)
    aux.train <- .fill.supervised.aux(aux.train, aux.init)
    augmented.ximp <- data.frame(ximp, aux.train, check.names = FALSE)
    supervised.model <- supervised.fit
    supervised.info <- list(
      enabled = TRUE,
      formula = paste(deparse(supervised.formula), collapse = ""),
      fit.formula = paste(deparse(fit.supervised.formula), collapse = ""),
      formula.scope = "supervised auxiliary forest",
      family = supervised.meta$family,
      xvar.names = supervised.meta$xvar.names,
      yvar.names = supervised.meta$yvar.names,
      aux.names = names(aux.train),
      aux.role = "predictor.only",
      aux.init = aux.init,
      aux.scale = aux.scale,
      fit.args = supervised.args,
      response.policy = "fit on rows with observed supervised response",
      response.n.used = sum(response.rows),
      response.n.omit = sum(!response.rows),
      subj.names = supervised.spec$subj.names %||% character(0),
      response.schema = response.schema,
      train.source = "predicted.oob then prediction fill",
      predict.source = "predicted",
      recompute.during.sweep = FALSE,
      model.name = .make.learner.name(0L, "supervised", prefix = learner.prefix)
    )
  }
  targets <- .resolve.targets(which.na, target.mode = target.mode)
  predictor.map.raw <- .resolve.predictor.map(targets, names(train.data), deployment.xvars)
  predictor.map <- if (isTRUE(supervised.info$enabled)) {
    .augment.predictor.map(targets, predictor.map.raw, supervised.info$aux.names)
  } else {
    predictor.map.raw
  }
  bad.targets <- targets[lengths(predictor.map[targets]) == 0L]
  if (length(bad.targets) > 0L) {
    stop("Some targets have no available predictors after applying deployment restrictions: ",
         paste(bad.targets, collapse = ", "), call. = FALSE)
  }
  miss.frac <- colMeans(which.na)
  sweep.order <- targets[order(miss.frac[targets], decreasing = FALSE)]
  if (persist.on.fit) {
    if (isTRUE(wipe)) .safe.unlink.dir(out.dir)
    .safe.dir.create(out.dir)
    .safe.dir.create(file.path(out.dir, learner.root))
  }
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
  if (isTRUE(supervised.info$enabled) && isTRUE(persist.on.fit)) {
    learner.path <- file.path(out.dir, learner.root, supervised.info$model.name)
    .msg("Saving supervised forest to ", learner.path, verbose = verbose)
    fast.save(supervised.model, learner.path, testing = FALSE)
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
      predictors.raw = predictor.map.raw[[yname]],
      predictors.supervised = if (isTRUE(supervised.info$enabled)) supervised.info$aux.names else character(0),
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
    xtrain <- augmented.ximp[trn, xvars, drop = FALSE]
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
      saved.weight <- .resolve.ood.weight(
        valid.ood.targets,
        weight = weight,
        warn.extra = FALSE
      )
      target.score.train <- do.call(cbind, lapply(valid.ood.targets, function(yname) {
        .eval.ood.reference(ood.delta[[yname]], target.reference[[yname]])
      }))
      colnames(target.score.train) <- valid.ood.targets
      row.score.train <- .aggregate.ood.row(
        target.score.train,
        weight = saved.weight,
        aggregate = "weighted.mean",
        aggregate.args = list()
      )
      row.reference <- .make.ood.reference(row.score.train)
      row.reference.targets <- valid.ood.targets
      row.reference.weight <- saved.weight
      row.reference.aggregate <- "weighted.mean"
      row.reference.aggregate.args <- list()
    }
    else {
      saved.weight <- setNames(numeric(0), character(0))
      target.score.train <- matrix(numeric(0), nrow = nrow(train.data), ncol = 0L,
                                   dimnames = list(NULL, character(0)))
      row.reference <- NULL
      row.reference.targets <- character(0)
      row.reference.weight <- setNames(numeric(0), character(0))
      row.reference.aggregate <- "weighted.mean"
      row.reference.aggregate.args <- list()
    }
    ood <- list(
      reference = "oob",
      aggregate = "weighted.mean",
      aggregate.args = list(),
      target.metric = c(
        numeric = "absolute.error",
        factor = "negative.log.probability",
        fallback = "misclass.or.rank.distance"
      ),
      targets = valid.ood.targets,
      weight = saved.weight,
      target.reference = target.reference,
      train.target.score = target.score.train,
      row.reference = row.reference,
      row.reference.targets = row.reference.targets,
      row.reference.weight = row.reference.weight,
      row.reference.aggregate = row.reference.aggregate,
      row.reference.aggregate.args = row.reference.aggregate.args,
      issues = ood.issues
    )
  }
  manifest <- list(
    spec.version = 2L,
    created.at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    formula = if (isTRUE(formula.used)) paste(deparse(formula), collapse = "") else NULL,
    formula.scope = if (isTRUE(formula.used)) {
      "initial imputation stage only"
    } else if (isTRUE(formula.supplied) && isTRUE(supervised.enabled)) {
      "ignored because supervised.formula defines the raw predictor block"
    } else {
      NULL
    },
    supervised = supervised.info,
    train.imputation = train.imputation,
    columns = names(train.data),
    schema = schema,
    init = init,
    scale = scale,
    targets = targets,
    sweep.order = sweep.order,
    predictor.map.raw = predictor.map.raw,
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
  if (!isTRUE(keep.models) && isTRUE(supervised.info$enabled)) {
    rm(supervised.model)
    gc()
    supervised.model <- NULL
  }
  object <- list(
    manifest = manifest,
    models = if (isTRUE(keep.models)) models else setNames(vector("list", length(targets)), targets),
    supervised.model = if (isTRUE(keep.models) && isTRUE(supervised.info$enabled)) supervised.model else NULL,
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
  object$manifest <- .normalize.impute.learn.manifest(object$manifest)
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
  if (.is.supervised.manifest(object$manifest)) {
    sup.mdl <- object$supervised.model
    if (is.null(sup.mdl)) {
      if (is.null(source.root)) {
        stop("The supervised forest is not available in memory and no saved ",
             "imputer path is attached to 'object'.",
             call. = FALSE)
      }
      .msg("Loading supervised forest from attached path before saving.",
           verbose = verbose)
      sup.mdl <- .fast.load.named(
        object$manifest$supervised$model.name,
        source.root,
        label = "supervised forest",
        strict = TRUE
      )
    }
    learner.path <- file.path(path, learner.root, object$manifest$supervised$model.name)
    .msg("Saving supervised forest to ", learner.path, verbose = verbose)
    fast.save(sup.mdl, learner.path, testing = FALSE)
    if (is.null(object$supervised.model)) {
      rm(sup.mdl)
      gc()
    }
  }
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
  manifest <- .normalize.impute.learn.manifest(readRDS(manifest.path))
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
  if (!is.null(manifest$predictor.map.raw)) {
    manifest$predictor.map.raw <- manifest$predictor.map.raw[use.targets]
  }
  manifest$learners <- manifest$learners[use.targets]
  models <- setNames(vector("list", length(use.targets)), use.targets)
  supervised.model <- NULL
  object <- list(
    manifest = manifest,
    models = models,
    supervised.model = supervised.model,
    ximp = NULL,
    path = normalizePath(path, mustWork = TRUE)
  )
  class(object) <- c("impute.learn.rfsrc", "impute.learn")
  if (!isTRUE(lazy)) {
    learner.root <- file.path(path, manifest$learner.root)
    if (.is.supervised.manifest(manifest)) {
      .msg("Loading supervised forest into memory...", verbose = verbose)
      object$supervised.model <- .fast.load.named(
        manifest$supervised$model.name,
        learner.root,
        label = "supervised forest",
        strict = TRUE
      )
    }
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
  object$manifest <- .normalize.impute.learn.manifest(object$manifest)
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
                             aggregate = c("bounded.product",
                                           "weighted.mean",
                                           "weighted.lp",
                                           "weighted.lp.log",
                                           "top.k"),
                             aggregate.args = list(),
                             return.details = FALSE,
                             return.reconstruction = FALSE,
                             verbose = TRUE,
                             ...) {
  cache.learners <- match.arg(cache.learners)
  if (!inherits(object, "impute.learn.rfsrc")) {
    stop("'object' must inherit from class 'impute.learn.rfsrc'.", call. = FALSE)
  }
  object$manifest <- .normalize.impute.learn.manifest(object$manifest)
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
  aggregate.spec <- .canonicalize.ood.aggregate(aggregate, aggregate.args)
  completed.data <- prep$working.data
  n <- nrow(completed.data)
  target.delta <- matrix(NA_real_, nrow = n, ncol = length(use.targets),
                         dimnames = list(NULL, use.targets))
  target.score <- matrix(NA_real_, nrow = n, ncol = length(use.targets),
                         dimnames = list(NULL, use.targets))
  ## Optional: store the model-based reconstruction used by OOD scoring.
  ## For each target, this is the learner prediction "target ~ other variables"
  ## evaluated on the completed deployment data.  This is distinct from
  ## prep$data / prep$working.data, which are completed/imputed copies of
  ## newdata rather than the per-target OOD reconstructions.
  target.reconstruction <- if (isTRUE(return.reconstruction)) {
    .na.data.from.schema(object$manifest$schema[use.targets], n)
  } else {
    NULL
  }
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
    if (isTRUE(return.reconstruction)) {
      reconstructed <- .extract.prediction(
        pred,
        info$family,
        object$manifest$schema[[target]]
      )
      if (!is.null(reconstructed) && length(reconstructed) == n) {
        target.reconstruction[[target]] <- reconstructed
      } else {
        record.issue(target, "Failed to extract target reconstruction from prediction object.")
      }
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
  score <- .aggregate.ood.row(
    target.score[, use.targets, drop = FALSE],
    weight = weight,
    aggregate = aggregate.spec$name,
    aggregate.args = aggregate.spec$args
  )
  weight.mask <- matrix(rep(weight > 0, each = n), nrow = n)
  targets.used <- rowSums(is.finite(target.score[, use.targets, drop = FALSE]) & weight.mask)
  score.percentile <- rep(NA_real_, n)
  row.reference.used <- FALSE
  row.reference.mode <- NULL
  row.reference.reason <- NULL
  row.reference.n.train <- 0L
  rebuilt.row.reference <- .rebuild.ood.row.reference(
    train.target.score = ood.ref$train.target.score,
    targets = use.targets,
    weight = weight,
    saved.targets = ood.ref$targets %||% character(0),
    aggregate = aggregate.spec$name,
    aggregate.args = aggregate.spec$args
  )
  if (!is.null(rebuilt.row.reference$reference)) {
    score.percentile <- .eval.ood.reference(score, rebuilt.row.reference$reference)
    row.reference.used <- TRUE
    row.reference.mode <- "recomputed.from.train.target.score"
    row.reference.n.train <- rebuilt.row.reference$n.finite %||% 0L
  }
  else {
    row.reference.reason <- rebuilt.row.reference$reason
    row.reference.targets <- ood.ref$row.reference.targets
    row.reference.weight <- ood.ref$row.reference.weight
    row.reference.aggregate <- ood.ref$row.reference.aggregate %||% ood.ref$aggregate %||% "weighted.mean"
    row.reference.aggregate.args <- ood.ref$row.reference.aggregate.args %||% ood.ref$aggregate.args %||% list()
    has.row.reference.meta <- !is.null(row.reference.targets) && !is.null(row.reference.weight)
    same.targets <- isTRUE(has.row.reference.meta) &&
      setequal(use.targets, row.reference.targets %||% character(0)) &&
      length(use.targets) == length(row.reference.targets %||% character(0))
    same.weight <- isTRUE(has.row.reference.meta) &&
      .same.ood.weight(use.targets, weight, row.reference.weight)
    same.aggregate <- .same.ood.aggregate(
      aggregate = aggregate.spec$name,
      aggregate.args = aggregate.spec$args,
      ref.aggregate = row.reference.aggregate,
      ref.aggregate.args = row.reference.aggregate.args
    )
    if (!is.null(ood.ref$row.reference) && isTRUE(same.targets) &&
        isTRUE(same.weight) && isTRUE(same.aggregate)) {
      score.percentile <- .eval.ood.reference(score, ood.ref$row.reference)
      row.reference.used <- TRUE
      row.reference.mode <- "saved.row.reference"
      row.reference.n.train <- ood.ref$row.reference$n %||% 0L
      row.reference.reason <- NULL
    }
    else if (is.null(row.reference.reason)) {
      if (is.null(ood.ref$row.reference)) {
        row.reference.reason <- "No saved row-level OOD reference is available."
      }
      else if (!isTRUE(has.row.reference.meta)) {
        row.reference.reason <- paste(
          "Saved row-level OOD calibration is legacy and lacks the",
          "training target-score matrix needed for arbitrary test-time",
          "weight and row-aggregate recalibration. Refit with the updated",
          "imputer to enable score.percentile for all target subsets,",
          "weights, and row aggregates."
        )
      }
      else if (!isTRUE(same.targets)) {
        row.reference.reason <- "Saved row-level OOD calibration requires scoring the original target set."
      }
      else if (!isTRUE(same.weight)) {
        row.reference.reason <- "Saved row-level OOD calibration requires the saved row-reference target weights."
      }
      else if (!isTRUE(same.aggregate)) {
        row.reference.reason <- "Saved row-level OOD calibration requires the saved row-reference aggregate."
      }
    }
  }
  if (length(unseen.rows) > 0L && any(unseen.rows)) {
    score[unseen.rows] <- .max.ood.aggregate.value(
      aggregate = aggregate.spec$name,
      aggregate.args = aggregate.spec$args
    )
    if (isTRUE(row.reference.used)) {
      score.percentile[unseen.rows] <- 1
    }
  }
  target.issues <- target.issues[lengths(target.issues) > 0L]
  reconstructed.data <- NULL
  if (isTRUE(return.reconstruction)) {
    ## User-facing reconstructed row: start from the harmonized raw deployment
    ## data and replace scored targets by their OOD reconstruction.  This keeps
    ## the original column layout and excludes supervised auxiliary columns.
    reconstructed.data <- prep$harmonized$data[, object$manifest$columns, drop = FALSE]
    for (target in use.targets) {
      reconstructed.data[[target]] <- target.reconstruction[[target]]
    }
    reconstructed.data <- .restore.schema(
      reconstructed.data,
      object$manifest$schema,
      restore.integer = TRUE
    )
  }
  out <- list(
    score = score,
    score.percentile = score.percentile,
    targets.used = targets.used,
    target.score = if (isTRUE(return.details)) target.score[, use.targets, drop = FALSE] else NULL,
    target.delta = if (isTRUE(return.details)) target.delta[, use.targets, drop = FALSE] else NULL,
    target.reconstruction = if (isTRUE(return.reconstruction)) target.reconstruction[, use.targets, drop = FALSE] else NULL,
    reconstructed.data = reconstructed.data,
    completed.data = if (isTRUE(return.details) && isTRUE(return.reconstruction)) {
      prep$data[, object$manifest$columns, drop = FALSE]
    } else {
      NULL
    },
    info = list(
      targets = use.targets,
      weight = weight,
      aggregate = aggregate.spec$name,
      aggregate.args = aggregate.spec$args,
      added.columns = prep$info$added.columns,
      dropped.extra.columns = prep$info$dropped.extra.columns,
      dropped.response.columns = prep$info$dropped.response.columns,
      unseen.levels = prep$info$unseen.levels,
      unseen.rows = unseen.rows,
      maxed.rows = unseen.rows,
      cache.learners = prep$cache.learners,
      n.disk.loads = length(unique(disk.load.targets)),
      disk.load.targets = unique(disk.load.targets),
      row.reference.used = row.reference.used,
      row.reference.mode = row.reference.mode,
      row.reference.n.train = row.reference.n.train,
      row.reference.reason = row.reference.reason,
      supervised = prep$info$supervised,
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
  x$manifest <- .normalize.impute.learn.manifest(x$manifest)
  cat("Predictive imputer (randomForestSRC)\n")
  cat("  imputation:    ", x$manifest$train.imputation %||% "<unknown>", "\n", sep = "")
  cat("  training rows: ", x$manifest$n.train, "\n", sep = "")
  cat("  training cols: ", x$manifest$p.train, "\n", sep = "")
  cat("  targets:       ", length(x$manifest$targets), "\n", sep = "")
  cat("  ood targets:   ", length(x$manifest$ood$targets %||% character(0)), "\n", sep = "")
  cat("  supervised:    ", if (.is.supervised.manifest(x$manifest)) "yes" else "no", "\n", sep = "")
  if (.is.supervised.manifest(x$manifest)) {
    cat("  aux columns:   ", length(x$manifest$supervised$aux.names %||% character(0)), "\n", sep = "")
    cat("  sup family:    ", x$manifest$supervised$family %||% "<unknown>", "\n", sep = "")
  }
  cat("  learner root:  ", x$manifest$learner.root, "\n", sep = "")
  cat("  path:          ", x$path %||% "<memory>", "\n", sep = "")
  invisible(x)
}
print.impute.learn <- print.impute.learn.rfsrc
