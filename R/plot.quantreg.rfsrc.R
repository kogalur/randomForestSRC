plot.quantreg.rfsrc <- function(x, prbL = .25, prbU = .75,
                                m.target = NULL, crps = TRUE, subset = NULL,
                                xlab = NULL, ylab = NULL, ...,
                                inset.args = list(), crps.null = TRUE,
                                quantreg.tau = NULL) {
  ##--------------------------------------------------------------
  ## Validate plotting arguments and select a stored response.
  ## m.target selects a panel only; no forest prediction is requested.
  ##--------------------------------------------------------------
  if (is.null(x$yvar)) {
    stop("Observed responses are required for a quantile diagnostic plot.",
         call. = FALSE)
  }
  valid.prob <- function(p) {
    is.numeric(p) && is.null(dim(p)) && length(p) == 1L &&
      is.finite(p) && p > 0 && p < 1
  }
  if (!valid.prob(prbL) || !valid.prob(prbU)) {
    stop("'prbL' and 'prbU' must be single probabilities in (0, 1).",
         call. = FALSE)
  }
  if (prbL >= prbU) {
    stop("'prbL' must be less than 'prbU'.", call. = FALSE)
  }
  if (!is.logical(crps) || length(crps) != 1L || is.na(crps)) {
    stop("'crps' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(crps.null) || length(crps.null) != 1L || is.na(crps.null)) {
    stop("'crps.null' must be TRUE or FALSE.", call. = FALSE)
  }
  ## Pinball reporting is opt-in for the plot. Its levels are distinct from
  ## prbL/prbU and are never passed through as graphical arguments.
  if (!is.null(quantreg.tau)) {
    quantreg.tau <- .quantreg.resolve.tau(quantreg.tau)
  }
  ## The middle request must lie inside the requested interval. Use the
  ## median when it is interior, otherwise use the interval's midpoint.
  prbM <- if (prbL < .5 && prbU > .5) .5 else prbL + (prbU - prbL) / 2
  requested.prob <- c(prbL, prbM, prbU)
  quant.list <- get.quantile(x, target.prob = requested.prob, pretty = FALSE)
  targets <- names(quant.list)
  if (!length(targets)) {
    stop("No continuous response quantiles are available.", call. = FALSE)
  }
  if (is.null(m.target)) {
    target <- targets[1L]
  } else {
    if (!is.character(m.target) || length(m.target) != 1L ||
        is.na(m.target) || !(m.target %in% targets)) {
      stop("'m.target' must name one available continuous response: ",
           paste(targets, collapse = ", "), ".", call. = FALSE)
    }
    target <- m.target
  }
  ## get.quantile() sorts and deduplicates its requests. Match the columns
  ## explicitly, including the rare case of a rounded midpoint at a bound.
  quant.dat <- quant.list[[target]][,
    match(requested.prob, sort(unique(requested.prob))), drop = FALSE]
  n <- nrow(quant.dat)
  if (is.data.frame(x$yvar) || is.matrix(x$yvar)) {
    if (!(target %in% colnames(x$yvar))) {
      stop("The selected response is missing from 'x$yvar'.", call. = FALSE)
    }
    y <- x$yvar[, target, drop = TRUE]
  } else {
    y <- x$yvar
  }
  if (!is.numeric(y) || length(y) != n) {
    stop("Observed responses and quantile rows are not aligned.",
         call. = FALSE)
  }
  y <- as.numeric(y)
  ##--------------------------------------------------------------
  ## Select rows once for both the main panel and the CRPS calculation.
  ## Repeated numeric indices represent repeated inclusion of those rows.
  ##--------------------------------------------------------------
  if (is.null(subset)) subset <- seq_len(n)
  if (is.logical(subset)) {
    if (!is.null(dim(subset)) || length(subset) != n || anyNA(subset)) {
      stop("A logical 'subset' must have one nonmissing entry per row.",
           call. = FALSE)
    }
    subset <- which(subset)
  } else if (!is.numeric(subset) || !is.null(dim(subset)) ||
             any(!is.finite(subset)) || any(subset != floor(subset)) ||
             any(subset < 1 | subset > n)) {
    stop("'subset' must contain positive integer indices of evaluated rows.",
         call. = FALSE)
  }
  if (!length(subset)) stop("requested subset is empty", call. = FALSE)
  pinball.error <- if (!is.null(quantreg.tau)) {
    get.pinball.error(x, tau = quantreg.tau, subset = subset, m.target = target)
  } else {
    NULL
  }
  y <- y[subset]
  quant.dat <- quant.dat[subset, , drop = FALSE]
  usable <- is.finite(y) & rowSums(!is.finite(quant.dat)) == 0L
  if (!any(usable)) {
    stop("No rows have both a finite response and finite requested quantiles.",
         call. = FALSE)
  }
  if (any(!usable)) {
    warning("Omitting ", sum(!usable),
            " row(s) with unavailable responses or quantiles from the main panel.",
            call. = FALSE)
  }
  y <- y[usable]
  quant.dat <- quant.dat[usable, , drop = FALSE]
  ##--------------------------------------------------------------
  ## Keep the two graphics interfaces separate.
  ## ... customizes the main panel; inset.args customizes the CRPS inset.
  ## Merge defaults by name so col, lwd, type, etc. can be overridden.
  ##--------------------------------------------------------------
  check.args <- function(args, label, reserved) {
    if (!is.list(args) || (length(args) &&
        (is.null(names(args)) || anyNA(names(args)) ||
         any(!nzchar(names(args))) || anyDuplicated(names(args))))) {
      stop(label, " must contain uniquely named graphical arguments.",
           call. = FALSE)
    }
    bad <- intersect(names(args), reserved)
    if (length(bad)) {
      stop(label, " cannot replace ", paste(bad, collapse = ", "), ".",
           call. = FALSE)
    }
  }
  dots <- list(...)
  check.args(dots, "'...'", c("x", "y", "type"))
  inset.plot.args <- NULL
  crps.dat <- null.dat <- NULL
  if (crps) {
    check.args(inset.args, "'inset.args'", c("x", "y"))
    crps.dat <- get.quantile.crps(x, pretty = FALSE, subset = subset)[[target]]
    if (crps.null) {
      ## The marginal training-response CDF provides a predictor-free reference.
      ## Its score uses the same subset and finite-CDF rows as the forest score.
      null.dat <- .plot.quantreg.null.crps(x, target, subset)
      if (is.null(null.dat)) {
        warning("Training responses are unavailable; omitting the null CRPS curve.",
                call. = FALSE)
      }
    }
    finite.crps <- is.finite(crps.dat$y) & is.finite(crps.dat$crps)
    if (!any(finite.crps)) {
      warning("No finite CRPS values are available; omitting the inset.",
              call. = FALSE)
    } else {
      null.values <- numeric(0)
      if (!is.null(null.dat)) {
        finite.null <- is.finite(null.dat$y) & is.finite(null.dat$crps)
        null.values <- null.dat$crps[finite.null]
        null.dat$y[!finite.null] <- NA_real_
        null.dat$crps[!finite.null] <- NA_real_
        if (!length(null.values)) null.dat <- NULL
      }
      inset.plot.args <- list(type = if (sum(finite.crps) > 1L) "l" else "p",
                         col = 2, lty = 1, lwd = 2, tck = .05, cex.axis = .75,
                         xlab = "", ylab = "",
                         xlim = range(crps.dat$y[is.finite(crps.dat$y)]),
                         ylim = range(c(crps.dat$crps[finite.crps], null.values)))
      inset.plot.args[names(inset.args)] <- inset.args
      ## Preserve gaps in a curve rather than joining across unavailable values.
      crps.dat$y[!finite.crps] <- NA_real_
      crps.dat$crps[!finite.crps] <- NA_real_
    }
  }
  ##--------------------------------------------------------------
  ## Observed responses on the horizontal axis, conditional quantiles
  ## vertically. Jitter changes display positions only, not CRPS inputs.
  ##--------------------------------------------------------------
  jitter.y <- jitter(y, factor = 10)
  rng <- range(c(y, quant.dat, jitter.y))
  main.args <- list(xlim = rng, ylim = rng, xlab = target,
                    ylab = "Target Quantiles", pch = 15, col = 4, cex = .75)
  main.args[names(dots)] <- dots
  if (!is.null(xlab)) main.args$xlab <- xlab
  if (!is.null(ylab)) main.args$ylab <- ylab
  do.call(graphics::plot.default, c(list(x = rng, y = rng, type = "n"), main.args))
  segments(jitter.y, quant.dat[, 1L], jitter.y, quant.dat[, 3L], col = "grey")
  points(jitter.y, quant.dat[, 1L], pch = "-", cex = 1)
  points(jitter.y, quant.dat[, 3L], pch = "-", cex = 1)
  point.args <- main.args[intersect(names(main.args),
                                   c("pch", "col", "cex", "bg", "lwd", "xpd"))]
  do.call(points, c(list(x = jitter.y, y = quant.dat[, 2L]), point.args))
  abline(0, 1, lty = 2, col = 2, untf = TRUE)
  if (!is.null(pinball.error)) {
    ## The requested score levels annotate the main panel only. The displayed
    ## quantile interval and the CRPS inset retain their separate controls.
    graphics::legend("bottomright", title = "Pinball loss",
      legend = paste0(names(pinball.error), ": ",
                      formatC(unname(pinball.error), digits = 4, format = "f")),
      bty = "n", cex = .75)
  }
  ##--------------------------------------------------------------
  ## Draw the inset within the current panel's plotting region.
  ## Changing plt (not fig) retains mfrow/mfcol/layout bookkeeping.
  ## Restore main-panel coordinates even when inset plotting fails.
  ##--------------------------------------------------------------
  if (!is.null(inset.plot.args)) {
    main.par <- par(c("pty", "plt", "mgp", "xlog", "ylog", "usr",
                       "xaxp", "yaxp", "new"))
    on.exit(par(main.par), add = TRUE)
    p <- main.par$plt
    inset.plt <- c(p[1L] + c(.02, .35) * diff(p[1:2]),
                   p[3L] + c(.68, .98) * diff(p[3:4]))
    par(pty = "m", plt = inset.plt, new = TRUE, mgp = c(0, 0, 0))
    ## Let the plot call draw the bottom axis and frame. Draw the y-axis
    ## at the right, honoring the same annotation options on both axes.
    show.right <- !identical(inset.plot.args$axes, FALSE) &&
      !identical(inset.plot.args$yaxt, "n")
    inset.plot.args$yaxt <- "n"
    do.call(graphics::plot.default,
            c(list(x = crps.dat$y, y = crps.dat$crps), inset.plot.args))
    if (!is.null(null.dat)) {
      ## A dashed black reference remains distinct when the forest color changes.
      ## Use points for a one-value curve and retain NA gaps in longer curves.
      graphics::lines(null.dat$y, null.dat$crps,
            type = if (sum(is.finite(null.dat$crps)) > 1L) "l" else "p",
            col = 1, lty = 2, lwd = 2, pch = 1)
      model.col <- if (is.null(inset.plot.args$col)) par("col") else inset.plot.args$col[1L]
      model.lty <- if (is.null(inset.plot.args$lty)) par("lty") else inset.plot.args$lty[1L]
      model.lwd <- if (is.null(inset.plot.args$lwd)) par("lwd") else inset.plot.args$lwd[1L]
      model.pch <- if (is.null(inset.plot.args$pch)) par("pch") else inset.plot.args$pch[1L]
      model.type <- inset.plot.args$type
      if (is.null(model.type)) model.type <- "p"
      model.line <- model.type %in% c("l", "b", "c", "o", "s", "S", "h")
      null.line <- sum(is.finite(null.dat$crps)) > 1L
      ## Mixing character and numeric line types can create an invalid "2"
      ## pattern. Keep the legend vector consistently named or numeric.
      legend.lty <- if (is.character(model.lty)) {
        c(if (model.line) model.lty else "blank",
          if (null.line) "dashed" else "blank")
      } else {
        c(if (model.line) model.lty else 0, if (null.line) 2 else 0)
      }
      graphics::legend("topleft", legend = c("Forest", "Null"),
             col = c(model.col, 1), lty = legend.lty,
             lwd = c(model.lwd, 2),
             pch = c(if (model.type %in% c("p", "b", "o")) model.pch else NA,
                     if (sum(is.finite(null.dat$crps)) > 1L) NA else 1),
             cex = .6, bty = "n", seg.len = 1.5, y.intersp = .9)
    }
    if (show.right) {
      axis.args <- inset.plot.args[intersect(names(inset.plot.args),
        c("cex.axis", "col.axis", "font.axis", "las", "tck", "tcl",
          "lwd", "lwd.ticks", "col.ticks", "gap.axis", "mgp", "fg"))]
      do.call(axis, c(list(side = 4L), axis.args))
    }
  }
  invisible(NULL)
}
plot.quantreg <- plot.quantreg.rfsrc
## Null CRPS on the stored reporting grid. Retain training multiplicities;
## using unique grid values as the reference sample would change the null law.
## This is a full-training empirical reference, including the evaluated row in
## training plots. It is not an OOB or leave-one-out estimate of the null CDF.
.plot.quantreg.null.crps <- function(x, target, subset) {
  response <- function(y) {
    if (is.data.frame(y) || is.matrix(y)) {
      if (!(target %in% colnames(y))) return(NULL)
      y <- y[, target, drop = TRUE]
    }
    if (!is.numeric(y)) return(NULL)
    as.numeric(y)
  }
  ## Read training responses through the same path as quantile construction.
  ## Never substitute the current prediction outcomes for the training sample.
  train <- response(.quantreg.training.yvar(x))
  train <- sort(train[is.finite(train)])
  if (!length(train)) return(NULL)
  q <- extract.quantile(x)[[target]]
  y <- response(x$yvar)
  if (is.null(y) || length(y) != nrow(q$cdf)) {
    stop("Observed responses and CDF rows are not aligned.", call. = FALSE)
  }
  y <- y[subset]
  grid <- q$yunq
  null.cdf <- findInterval(grid, train) / length(train)
  brier <- vapply(seq_along(grid), function(k) {
    ## Match get.quantile.crps()'s evaluation mask at every threshold,
    ## including repeated subset indices and missing forest predictions.
    ok <- is.finite(y) & is.finite(q$cdf[subset, k])
    if (!any(ok)) return(NA_real_)
    mean((as.numeric(y[ok] <= grid[k]) - null.cdf[k])^2)
  }, numeric(1))
  integral <- c(0, cumsum(diff(grid) *
                          (utils::head(brier, -1L) + utils::tail(brier, -1L)) / 2))
  width <- grid - grid[1L]
  value <- rep(NA_real_, length(grid))
  positive <- width > 0
  value[positive] <- integral[positive] / width[positive]
  data.frame(y = grid, crps = value)
}
