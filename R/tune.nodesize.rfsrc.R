## Finalized: uses numeric 'ssize' during tuning (OOB/holdout), supports method = "grid" or "golden"
tune.nodesize.rfsrc <- function(formula, data,
                                nodesize.try = c(1:9, seq(10, 150, by = 5)),
                                ntree.try = 100,
                                sampsize = function(x) {min(x * 0.632, max(150, x^(4/5)))},
                                nsplit = 1,
                                method = c("grid", "golden"),
                                ## golden-search controls (ignored for grid)
                                final.window = 5L,
                                reps.initial = 2L,
                                reps.final = 3L,
                                max.iter = 50L,
                                trace = TRUE,
                                seed = NULL,
                                ...) {
  method <- match.arg(method)
  ## capture and sanitize dots (avoid duplicate formal arguments)
  dots <- list(...)
  perf.type <- if ("perf.type" %in% names(dots)) dots$perf.type else "default"
  drop.from.dots <- c("perf.type","forest","save.memory","ntree","nodesize",
                      "sampsize","nsplit","data","formula")
  if (length(dots)) {
    dots <- dots[setdiff(names(dots), drop.from.dots)]
  }
  ## reproducibility: save/restore RNG; choose base seed even if user did not pass one
  have.old.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (have.old.seed) old.seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (exists("old.seed", inherits = FALSE)) assign(".Random.seed", old.seed, envir = .GlobalEnv)
  }, add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  base.seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else seed
  ## build stump to reconstruct a clean model frame (handles NA, factors, Surv, etc.)
  stump <- rfsrc(formula, data, nodedepth = 0, perf.type = "none",
                 save.memory = TRUE, ntree = 1, splitrule = "random")
  n <- stump$n
  yvar.names <- stump$yvar.names
  data <- data.frame(stump$yvar, stump$xvar)
  colnames(data)[seq_along(yvar.names)] <- yvar.names
  rm(stump)
  ## compute effective subsample size and coerce to integer bounds
  if (is.function(sampsize)) {
    ssize <- sampsize(n)
  } else if (is.numeric(sampsize) && length(sampsize) == 1L) {
    ssize <- sampsize
  } else if (is.numeric(sampsize) && length(sampsize) > 1L) {
    ## class-specific vector: use its total count per tree during tuning
    ssize <- sum(sampsize)
  } else {
    stop("'sampsize' must be a function or numeric (scalar or vector).")
  }
  ssize <- max(2L, min(n, as.integer(round(ssize))))
  if (trace) cat("effective sampsize:", ssize, "of", n, "\n")
  ## holdout equal to ssize if feasible; otherwise rely on OOB
  use.holdout <- (2L * ssize) < n
  if (use.holdout) {
    tst <- sample.int(n, ssize)
    trn <- setdiff(seq_len(n), tst)
    newdata <- data[tst, , drop = FALSE]
  } else {
    trn <- seq_len(n)
    newdata <- NULL
  }
  ## admissible upper bound on nodesize from the tree subsample size
  U.calc <- max(10L, ssize %/% 2L)
  ## filter the provided grid to admissible values
  grid.cand <- unique(sort(as.integer(nodesize.try[nodesize.try <= U.calc])))
  if (!length(grid.cand) && method == "grid") {
    stop("No admissible 'nodesize.try' after filtering by ssize/2. Increase 'ssize' or adjust the grid.")
  }
  ## evaluation kernel; uses numeric 'ssize' during tuning
  ## seeds are deterministic for reproducibility when 'seed' is supplied
  eval.nodesize <- function(nsz, reps) {
    vals <- numeric(reps)
    for (j in seq_len(reps)) {
      set.seed(base.seed + j)  ## shared across candidates within replicate j
      if (is.null(newdata)) {
        tr.args <- c(list(formula = formula,
                          data = data[trn, , drop = FALSE],
                          ntree = ntree.try,
                          nodesize = nsz,
                          sampsize = ssize,   ## numeric ssize during tuning
                          nsplit = nsplit,
                          perf.type = perf.type,
                          forest = FALSE,
                          save.memory = TRUE),
                     dots)
        fit <- tryCatch(do.call(rfsrc.fast, tr.args), error = function(e) e)
        if (inherits(fit, "error")) return(NA_real_)
        vals[j] <- tryCatch(mean(get.mv.error(fit, TRUE), na.rm = TRUE),
                            error = function(e) NA_real_)
      } else {
        tr.args <- c(list(formula = formula,
                          data = data[trn, , drop = FALSE],
                          ntree = ntree.try,
                          nodesize = nsz,
                          sampsize = ssize,   ## numeric ssize during tuning
                          nsplit = nsplit,
                          perf.type = "none",
                          forest = TRUE,
                          save.memory = FALSE),
                     dots)
        fit <- tryCatch(do.call(rfsrc.fast, tr.args), error = function(e) e)
        if (inherits(fit, "error")) return(NA_real_)
        pr <- tryCatch(predict(fit, newdata = newdata, perf.type = perf.type),
                       error = function(e) e)
        if (inherits(pr, "error")) return(NA_real_)
        vals[j] <- tryCatch(mean(get.mv.error(pr, TRUE), na.rm = TRUE),
                            error = function(e) NA_real_)
      }
    }
    if (trace && is.finite(mean(vals))) {
      cat("nodesize =", nsz, " reps =", reps, " error =", paste0(100 * round(mean(vals), 4), "%"), "\n")
    }
    mean(vals)
  }
  ## ----------------------------
  ## method = "grid"
  ## ----------------------------
  if (method == "grid") {
    err <- vapply(grid.cand, function(nsz) eval.nodesize(nsz, reps = 1L), numeric(1))
    if (all(is.na(err))) {
      msg <- if (is.null(newdata)) "All OOB errors are NA; check settings (especially 'sampsize')."
             else "All holdout errors are NA; check settings (especially 'sampsize')."
      warning(msg)
      return(list(nsize.opt = NA_integer_,
                  err = data.frame(nodesize = grid.cand, err = err)))
    }
    nsize.opt <- grid.cand[which.min(err)]
    if (trace) cat("optimal nodesize:", nsize.opt, "\n")
    return(list(nsize.opt = nsize.opt,
                err = data.frame(nodesize = grid.cand, err = err)))
  }
  ## ----------------------------
  ## method = "golden" (guarded)
  ## ----------------------------
  ## always probe the left boundary region (1..min(9, U.calc))
  L <- 1L; U <- U.calc
  if (U <= L) stop("Admissible nodesize range is empty; increase 'ssize' or check data.")
  left.max <- min(9L, U)
  left.cand <- seq.int(1L, left.max)
  left.err <- vapply(left.cand, function(nsz) eval.nodesize(nsz, reps.initial), numeric(1))
  best.left.idx <- which.min(left.err)
  best.left.nsz <- left.cand[best.left.idx]
  best.left.err <- left.err[best.left.idx]
  ## interior points for golden shrinkage
  phi <- (1 + sqrt(5)) / 2
  invphi <- 1 / phi
  xL2 <- max(L, min(U, floor(U - invphi * (U - L))))
  xR2 <- max(L, min(U, ceiling(L + invphi * (U - L))))
  if (xL2 == xR2) xR2 <- min(U, xL2 + 1L)
  fL2 <- eval.nodesize(xL2, reps.initial)
  fR2 <- eval.nodesize(xR2, reps.initial)
  ## if left boundary looks best, jump to a local sweep
  if (is.finite(best.left.err) && (best.left.err <= min(fL2, fR2, na.rm = TRUE))) {
    Lf <- max(L, best.left.nsz - final.window)
    Uf <- min(U, best.left.nsz + final.window)
    cand <- sort(unique(c(seq.int(Lf, Uf), left.cand)))
    err <- vapply(cand, function(nsz) eval.nodesize(nsz, reps.final), numeric(1))
    nsize.opt <- cand[which.min(err)]
    if (trace) cat("optimal nodesize:", nsize.opt, "\n")
    return(list(nsize.opt = nsize.opt, err = data.frame(nodesize = cand, err = err)))
  }
  ## otherwise proceed with golden shrinkage
  iter <- 0L
  while ((U - L) > final.window && iter < max.iter) {
    if (!is.finite(fL2) || !is.finite(fR2)) break
    if (fR2 < fL2) {
      L <- xL2; xL2 <- xR2; fL2 <- fR2
      xR2 <- max(L, min(U, ceiling(L + invphi * (U - L))))
      if (xR2 == xL2) xR2 <- min(U, xR2 + 1L)
      fR2 <- eval.nodesize(xR2, reps.initial)
    } else {
      U <- xR2; xR2 <- xL2; fR2 <- fL2
      xL2 <- max(L, min(U, floor(U - invphi * (U - L))))
      if (xL2 == xR2) xL2 <- max(L, xL2 - 1L)
      fL2 <- eval.nodesize(xL2, reps.initial)
    }
    iter <- iter + 1L
  }
  ## final local sweep
  cand <- sort(unique(c(seq.int(L, U), left.cand)))
  err <- vapply(cand, function(nsz) eval.nodesize(nsz, reps.final), numeric(1))
  if (all(is.na(err))) {
    msg <- if (is.null(newdata)) "All OOB errors are NA; check 'sampsize' and other settings."
           else "All holdout errors are NA; check 'sampsize' and other settings."
    warning(msg)
    return(list(nsize.opt = NA_integer_, err = data.frame(nodesize = cand, err = err)))
  }
  nsize.opt <- cand[which.min(err)]
  if (trace) cat("optimal nodesize:", nsize.opt, "\n")
  list(nsize.opt = nsize.opt, err = data.frame(nodesize = cand, err = err))
}
tune.nodesize <- tune.nodesize.rfsrc
