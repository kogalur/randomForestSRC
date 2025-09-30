tune.rfsrc <- function(formula, data,
      mtry.start = ncol(data) / 2,
      nodesize.try = c(1:9, seq(10, 100, by = 5)), ntree.try = 100,
      sampsize = function(x){min(x * .632, max(150, x ^ (3/4)))},
      nsplit = 1, step.factor = 1.25, improve = 1e-3, strikeout = 3, max.iter = 25,
      method = c("grid","golden"),
      ## golden controls (ignored for grid)
      final.window = 5L, reps.initial = 2L, reps.final = 3L,
      trace = FALSE, do.best = TRUE, seed = NULL, ...) {
  method <- match.arg(method)
  ## input checks 
  if (improve < 0) stop("improve must be non-negative.")
  if (step.factor <= 1) stop("step.factor must be great than 1.")
  if (missing(formula)) stop("a formula must be supplied (only supervised forests allowed).")
  ## capture/sanitize dots (avoid duplicate formals)
  dots <- list(...)
  perf.type <- if ("perf.type" %in% names(dots)) dots$perf.type else "default"
  drop.from.dots <- c("perf.type","forest","save.memory","ntree","mtry","nodesize",
                      "sampsize","nsplit","data","formula")
  if (length(dots)) dots <- dots[setdiff(names(dots), drop.from.dots)]
  ## reproducibility and a base seed for CRN
  have.old.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (have.old.seed) old.seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (exists("old.seed", inherits = FALSE)) assign(".Random.seed", old.seed, envir = .GlobalEnv)
  }, add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  base.seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else seed
  ## clean model frame via stump (handles NA & factor levels consistently)
  stump <- rfsrc(formula, data, nodedepth = 0, perf.type = "none",
                 save.memory = TRUE, ntree = 1, splitrule = "random")
  n <- stump$n
  yvar.names <- stump$yvar.names
  data <- data.frame(stump$yvar, stump$xvar)
  colnames(data)[seq_along(yvar.names)] <- yvar.names
  p <- ncol(data) - length(yvar.names)  ## number of predictors
  rm(stump)
  ## effective numeric ssize for tuning and holdout
  if (is.function(sampsize)) {
    ssize <- sampsize(n)
  } else if (is.numeric(sampsize) && length(sampsize) == 1L) {
    ssize <- sampsize
  } else if (is.numeric(sampsize) && length(sampsize) > 1L) {
    ssize <- sum(sampsize)
  } else {
    stop("'sampsize' must be a function or numeric (scalar or vector).")
  }
  ssize <- max(2L, min(n, as.integer(round(ssize))))
  if (trace) cat("effective sampsize:", ssize, "of", n, "\n")
  ## holdout equal to ssize if possible; else rely on OOB
  use.holdout <- (2L * ssize) < n
  if (use.holdout) {
    tst <- sample.int(n, ssize)
    trn <- setdiff(seq_len(n), tst)
    newdata <- data[tst, , drop = FALSE]
  } else {
    trn <- seq_len(n)
    newdata <- NULL
  }
  ## admissible nodesize bound
  U.calc <- max(10L, ssize %/% 2L)
  ## evaluation kernel for a pair (nodesize, mtry) with CRN
  eval.pair <- function(nsz, mtry, reps) {
    vals <- numeric(reps)
    for (j in seq_len(reps)) {
      s <- as.integer((base.seed + 10007L * as.integer(nsz) + 101L * as.integer(mtry) + j) %% .Machine$integer.max)
      if (s <= 0L) s <- j
      set.seed(s)
      if (is.null(newdata)) {
        tr.args <- c(list(formula = formula,
                          data = data[trn, , drop = FALSE],
                          ntree = ntree.try,
                          mtry = mtry,
                          nodesize = nsz,
                          sampsize = ssize,      # numeric subsample during tuning
                          nsplit = nsplit,
                          perf.type = perf.type,
                          forest = FALSE,
                          save.memory = TRUE),
                     dots)
        fit <- tryCatch(do.call(rfsrc.fast, tr.args), error = function(e) e)
        if (inherits(fit, "error")) return(NA_real_)
        vals[j] <- tryCatch(mean(get.mv.error(fit, TRUE), na.rm = TRUE), error = function(e) NA_real_)
      } else {
        tr.args <- c(list(formula = formula,
                          data = data[trn, , drop = FALSE],
                          ntree = ntree.try,
                          mtry = mtry,
                          nodesize = nsz,
                          sampsize = ssize,
                          nsplit = nsplit,
                          perf.type = "none",
                          forest = TRUE,
                          save.memory = FALSE),
                     dots)
        fit <- tryCatch(do.call(rfsrc.fast, tr.args), error = function(e) e)
        if (inherits(fit, "error")) return(NA_real_)
        pr <- tryCatch(predict(fit, newdata = newdata, perf.type = perf.type), error = function(e) e)
        if (inherits(pr, "error")) return(NA_real_)
        vals[j] <- tryCatch(mean(get.mv.error(pr, TRUE), na.rm = TRUE), error = function(e) NA_real_)
      }
    }
    if (trace && is.finite(mean(vals))) {
      cat("nodesize =", nsz, " mtry =", mtry, " reps =", reps,
          " error =", paste0(100 * round(mean(vals), 4), "%"), "\n")
    }
    mean(vals)
  }
  ## ---------------------------------------
  ## method == "grid": original method
  ## ---------------------------------------
  if (method == "grid") {
    grid.cand <- unique(sort(as.integer(nodesize.try[nodesize.try <= U.calc])))
    if (!length(grid.cand)) stop("No admissible 'nodesize.try' after filtering by ssize/2.")
    res <- list()
    counter1 <- 0
    for (nsz in grid.cand) {
      counter1 <- counter1 + 1
      ## acquire starting error at mtry.start
      if (is.null(newdata)) {
        o <- tryCatch(do.call(rfsrc.fast,
                              c(list(formula = formula, data = data,
                                     ntree = ntree.try, mtry = mtry.start, nodesize = nsz,
                                     sampsize = ssize, nsplit = nsplit, perf.type = perf.type),
                                dots)),
                      error = function(e) e)
        if (inherits(o, "error")) stop(conditionMessage(o))
        error.old <- mean(get.mv.error(o, TRUE), na.rm = TRUE)
      } else {
        o <- tryCatch(do.call(rfsrc.fast,
                              c(list(formula = formula, data = data[trn, , drop = FALSE],
                                     ntree = ntree.try, mtry = mtry.start, nodesize = nsz,
                                     sampsize = ssize, nsplit = nsplit, perf.type = "none",
                                     forest = TRUE, save.memory = FALSE),
                                dots)),
                      error = function(e) e)
        if (inherits(o, "error")) stop(conditionMessage(o))
        error.old <- mean(get.mv.error(predict(o, newdata, perf.type = perf.type), TRUE), na.rm = TRUE)
      }
      mtry.start <- o$mtry
      mtry.max <- p
      if (trace) {
        cat("nodesize = ", nsz, " mtry =", mtry.start,
            " error =", paste(100 * round(error.old, 4), "%", sep = ""), "\n")
      }
      if (!is.null(o$leaf.count) && mean(o$leaf.count) <= 2) break
      oob.error <- list()
      oob.error[[1]] <- error.old
      names(oob.error)[1] <- mtry.start
      for (direction in c("left", "right")) {
        if (trace) cat("Searching", direction, "...\n")
        Improve <- 1.1 * improve
        mtry.best <- mtry.start
        mtry.cur <- mtry.start
        counter2 <- 1
        strikes <- 0
        while (counter2 <= max.iter && (Improve >= improve || (Improve < 0 && strikes < strikeout))) {
          counter2 <- counter2 + 1
          if (Improve < 0) strikes <- strikes + 1
          mtry.old <- mtry.cur
          if (direction == "left") {
            mtry.cur <- max(1, min(ceiling(mtry.cur / step.factor), mtry.cur - 1))
          } else {
            mtry.cur <- min(mtry.max, max(floor(mtry.cur * step.factor), mtry.cur + 1))
          }
          if (mtry.cur == mtry.old) break
          error.cur <- eval.pair(nsz, mtry.cur, reps = 1L)
          if (trace) {
            cat("nodesize = ", nsz, " mtry =", mtry.cur,
                " error =", paste(100 * round(error.cur, 4), "%", sep = ""), "\n")
          }
          oob.error[[as.character(mtry.cur)]] <- error.cur
          Improve <- 1 - error.cur / error.old
          if (trace) cat(Improve, improve, "\n")
          if (Improve > improve) {
            error.old <- error.cur
            mtry.best <- mtry.cur
          }
        }
      }
      mtry.vec <- sort(as.numeric(names(oob.error)))
      err.vec <- unlist(oob.error[as.character(mtry.vec)])
      res[[counter1]] <- cbind(nodesize = nsz, mtry = mtry.vec, err = err.vec)
    }
    if (is.null(res)) stop("NULL results - check tuning settings.")
    res <- do.call(rbind, res)
    res <- res[order(res[, 1], res[, 2]), , drop = FALSE]
    rownames(res) <- seq_len(nrow(res))
    opt.idx <- which.min(res[, 3])
    rf <- NULL
    if (do.best) {
      rf <- do.call(rfsrc.fast,
                    c(list(formula = formula, data = data,
                           mtry = res[opt.idx, 2], nodesize = res[opt.idx, 1],
                           sampsize = sampsize, nsplit = nsplit, forest = TRUE),
                      dots))
    }
    return(list(results = res, optimal = res[opt.idx, -3], rf = rf))
  }
  ## ---------------------------------------
  ## method == "golden": guarded coordinate line search
  ## ---------------------------------------
  ## 1D golden on nodesize with left guard; final local sweep
  golden.search.nodesize <- function(mfix) {
    L <- 1L; U <- U.calc
    if (U <= L) stop("Admissible nodesize range is empty; increase 'ssize' or check data.")
    left.max <- min(9L, U)
    left.cand <- seq.int(1L, left.max)
    left.err <- vapply(left.cand, function(nsz) eval.pair(nsz, mfix, reps.initial), numeric(1))
    best.left.idx <- which.min(left.err); best.left.nsz <- left.cand[best.left.idx]
    best.left.err <- left.err[best.left.idx]
    phi <- (1 + sqrt(5)) / 2; invphi <- 1 / phi
    xL2 <- max(L, min(U, floor(U - invphi * (U - L))))
    xR2 <- max(L, min(U, ceiling(L + invphi * (U - L))))
    if (xL2 == xR2) xR2 <- min(U, xL2 + 1L)
    fL2 <- eval.pair(xL2, mfix, reps.initial)
    fR2 <- eval.pair(xR2, mfix, reps.initial)
    if (is.finite(best.left.err) && (best.left.err <= min(fL2, fR2, na.rm = TRUE))) {
      Lf <- max(L, best.left.nsz - final.window)
      Uf <- min(U, best.left.nsz + final.window)
      cand <- sort(unique(c(seq.int(Lf, Uf), left.cand)))
      err <- vapply(cand, function(nsz) eval.pair(nsz, mfix, reps.final), numeric(1))
      return(list(best = cand[which.min(err)], grid = data.frame(nodesize = cand, mtry = mfix, err = err)))
    }
    iter <- 0L
    while ((U - L) > final.window && iter < max.iter) {
      if (!is.finite(fL2) || !is.finite(fR2)) break
      if (fR2 < fL2) {
        L <- xL2; xL2 <- xR2; fL2 <- fR2
        xR2 <- max(L, min(U, ceiling(L + invphi * (U - L))))
        if (xR2 == xL2) xR2 <- min(U, xR2 + 1L)
        fR2 <- eval.pair(xR2, mfix, reps.initial)
      } else {
        U <- xR2; xR2 <- xL2; fR2 <- fL2
        xL2 <- max(L, min(U, floor(U - invphi * (U - L))))
        if (xL2 == xR2) xL2 <- max(L, xL2 - 1L)
        fL2 <- eval.pair(xL2, mfix, reps.initial)
      }
      iter <- iter + 1L
    }
    cand <- sort(unique(c(seq.int(L, U), left.cand)))
    err <- vapply(cand, function(nsz) eval.pair(nsz, mfix, reps.final), numeric(1))
    list(best = cand[which.min(err)], grid = data.frame(nodesize = cand, mtry = mfix, err = err))
  }
  ## 1D golden on mtry with left guard; final local sweep
  golden.search.mtry <- function(nfix) {
    L <- 1L; U <- p
    if (U <= L) stop("Admissible mtry range is empty; check data.")
    left.max <- min(9L, U)
    left.cand <- seq.int(1L, left.max)
    left.err <- vapply(left.cand, function(m) eval.pair(nfix, m, reps.initial), numeric(1))
    best.left.idx <- which.min(left.err); best.left.m <- left.cand[best.left.idx]
    best.left.err <- left.err[best.left.idx]
    phi <- (1 + sqrt(5)) / 2; invphi <- 1 / phi
    xL2 <- max(L, min(U, floor(U - invphi * (U - L))))
    xR2 <- max(L, min(U, ceiling(L + invphi * (U - L))))
    if (xL2 == xR2) xR2 <- min(U, xL2 + 1L)
    fL2 <- eval.pair(nfix, xL2, reps.initial)
    fR2 <- eval.pair(nfix, xR2, reps.initial)
    if (is.finite(best.left.err) && (best.left.err <= min(fL2, fR2, na.rm = TRUE))) {
      Lf <- max(L, best.left.m - final.window)
      Uf <- min(U, best.left.m + final.window)
      cand <- sort(unique(c(seq.int(Lf, Uf), left.cand)))
      err <- vapply(cand, function(m) eval.pair(nfix, m, reps.final), numeric(1))
      return(list(best = cand[which.min(err)], grid = data.frame(nodesize = nfix, mtry = cand, err = err)))
    }
    iter <- 0L
    while ((U - L) > final.window && iter < max.iter) {
      if (!is.finite(fL2) || !is.finite(fR2)) break
      if (fR2 < fL2) {
        L <- xL2; xL2 <- xR2; fL2 <- fR2
        xR2 <- max(L, min(U, ceiling(L + invphi * (U - L))))
        if (xR2 == xL2) xR2 <- min(U, xR2 + 1L)
        fR2 <- eval.pair(nfix, xR2, reps.initial)
      } else {
        U <- xR2; xR2 <- xL2; fR2 <- fL2
        xL2 <- max(L, min(U, floor(U - invphi * (U - L))))
        if (xL2 == xR2) xL2 <- max(L, xL2 - 1L)
        fL2 <- eval.pair(nfix, xL2, reps.initial)
      }
      iter <- iter + 1L
    }
    cand <- sort(unique(c(seq.int(L, U), left.cand)))
    err <- vapply(cand, function(m) eval.pair(nfix, m, reps.final), numeric(1))
    list(best = cand[which.min(err)], grid = data.frame(nodesize = nfix, mtry = cand, err = err))
  }
  ## coordinate descent: alternate golden on nodesize and mtry
  res.grid <- list()
  res.row <- 0L
  add.grid <- function(df) {
    if (is.null(df) || !nrow(df)) return()
    res.row <<- res.row + 1L
    res.grid[[res.row]] <<- df
  }
  mtry.cur <- as.integer(max(1L, min(p, round(mtry.start))))
  best.err <- Inf
  strikes <- 0L
  for (iter.cd in seq_len(max.iter)) {
    nsres <- golden.search.nodesize(mtry.cur)
    add.grid(nsres$grid)
    nsz.cur <- nsres$best
    mtres <- golden.search.mtry(nsz.cur)
    add.grid(mtres$grid)
    mtry.new <- mtres$best
    cur.err <- eval.pair(nsz.cur, mtry.new, reps.final)
    if (trace) {
      cat("CD iter", iter.cd, ": nodesize =", nsz.cur, " mtry =", mtry.new,
          " error =", paste0(100 * round(cur.err, 4), "%"), "\n")
    }
    if (is.finite(best.err)) {
      rel.imp <- 1 - cur.err / best.err
      if (rel.imp <= improve) {
        strikes <- strikes + 1L
      } else {
        strikes <- 0L
      }
      if (rel.imp <= improve && strikes >= strikeout) break
    }
    if (cur.err < best.err) best.err <- cur.err
    mtry.cur <- mtry.new
  }
  if (!length(res.grid)) stop("No evaluations; check settings.")
  res <- do.call(rbind, res.grid)
  res <- stats::aggregate(err ~ nodesize + mtry, data = res, FUN = mean)  # average if duplicates
  res <- res[order(res$nodesize, res$mtry), , drop = FALSE]
  colnames(res) <- c("nodesize", "mtry", "err")
  opt.idx <- which.min(res$err)
  opt.nsz <- res$nodesize[opt.idx]
  opt.mtry <- res$mtry[opt.idx]
  rf <- NULL
  if (do.best) {
    rf <- do.call(rfsrc.fast,
                  c(list(formula = formula, data = data,
                         mtry = opt.mtry, nodesize = opt.nsz,
                         sampsize = sampsize, nsplit = nsplit, forest = TRUE),
                    dots))
  }
  list(results = as.matrix(res), optimal = c(nodesize = opt.nsz, mtry = opt.mtry), rf = rf)
}
tune <- tune.rfsrc
