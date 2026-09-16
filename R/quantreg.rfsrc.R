quantreg.rfsrc <- function(formula, data, object, newdata,
                           method = "exchangeable", splitrule = NULL,
                           prob = NULL, prob.epsilon = NULL,
                           oob = TRUE, fast = FALSE, maxn = 1e3, ...)
{
  ## ------------------------------------------------------------------
  ## Resolve the current operation and options before calling the forest.
  ## NULL newdata has the same meaning as an omitted prediction dataset.
  ## ------------------------------------------------------------------
  grow <- missing(object)
  test <- !missing(newdata) && !is.null(newdata)
  ## Use the descriptive method name in all newly returned components.
  ## Keep the former name as an alias for existing user calls.
  method <- match.arg(method, c("exchangeable", "forest", "local",
                                "gk", "GK", "G-K", "g-k"))
  if (method == "local") method <- "exchangeable"
  if (!(method %in% c("exchangeable", "forest"))) method <- "gk"
  .quantreg.check.flag(oob, "oob")
  .quantreg.check.flag(fast, "fast")
  if (!is.numeric(maxn) || length(maxn) != 1L || is.na(maxn) ||
      maxn < 1 || (is.finite(maxn) && maxn != floor(maxn))) {
    stop("'maxn' must be a positive integer or Inf.", call. = FALSE)
  }
  if (!is.null(prob)) {
    .quantreg.check.prob(prob)
    prob <- sort(prob)
  }
  if (!is.null(prob.epsilon) &&
      (!is.numeric(prob.epsilon) || length(prob.epsilon) != 1L ||
       !is.finite(prob.epsilon) || prob.epsilon <= 0 || prob.epsilon >= 1)) {
    stop("'prob.epsilon' must be a finite number in (0, 1).", call. = FALSE)
  }
  dots.all <- list(...)
  if (!is.null(dots.all[["outcome"]])) {
    outcome <- match.arg(dots.all[["outcome"]], c("train", "test"))
    if (outcome != "train") {
      stop("'outcome = \"test\"' is not supported by quantreg(); ",
           "supply newdata for ordinary prediction.", call. = FALSE)
    }
  }
  dots.all$outcome <- NULL
  if (test && !identical(class(newdata), "data.frame")) {
    newdata <- as.data.frame(newdata, optional = TRUE)
  }
  ## ------------------------------------------------------------------
  ## Grow the forest, retain full-ensemble and OOB training residuals,
  ## and select the residual bank once. With oob = TRUE the selected bank
  ## uses OOB mean predictions when that output component is available.
  ## The exchangeable method pools its finite residuals with equal weights.
  ## Later restore/newdata calls preserve this original bank selection.
  ## ------------------------------------------------------------------
  if (grow) {
    if (missing(formula)) stop("'formula' is missing.", call. = FALSE)
    if (missing(data)) stop("'data' is missing.", call. = FALSE)
    if (!identical(class(data), "data.frame")) {
      data <- as.data.frame(data, optional = TRUE)
    }
    fmly <- parseFormula(formula, data)$family
    if (!(fmly %in% c("regr", "regr+", "mix+"))) {
      stop("This function requires at least one continuous outcome.", call. = FALSE)
    }
    if (is.null(splitrule)) splitrule <- "la.quantile.regr"
    if (fmly %in% c("regr+", "mix+") &&
        splitrule %in% c("la.quantile.regr", "quantile.regr")) {
      splitrule <- NULL
    }
    rfnames <- setdiff(get.rfnames(hidden = TRUE), c("formula", "data", "splitrule"))
    dots <- dots.all[names(dots.all) %in% rfnames]
    dots$prob <- prob
    dots$prob.epsilon <- prob.epsilon
    dots$gk.quantile <- method == "gk"
    dots$quantile.regr <- TRUE
    if (method == "forest") dots$forest.wt <- if (oob) "oob" else TRUE
    engine <- if (fast) "rfsrc.fast" else "rfsrc"
    object <- do.call(engine, c(list(formula = formula, data = data,
                                    splitrule = splitrule), dots))
    if (!is.list(object) || is.null(object$yvar) || is.null(object$n)) {
      stop("The grow interface did not retain the responses and predictions ",
           "required by quantreg().", call. = FALSE)
    }
  } else {
    if (!inherits(object, "quantreg") || !inherits(object, "rfsrc") ||
        !inherits(object, "grow")) {
      stop("'object' must be the original quantreg grow object.", call. = FALSE)
    }
  }
  if (!(object$family %in% c("regr", "regr+", "mix+"))) {
    stop("This function requires at least one continuous outcome.", call. = FALSE)
  }
  ## The grow object exposes actual residual/residual.oob values. The saved
  ## forest retains both banks and the original source selection for reuse.
  object <- .quantreg.save.residuals(object, oob, initialize = grow)
  if (is.null(object$forest$prob)) {
    stop("The forest did not retain its quantile probability settings.", call. = FALSE)
  }
  prob.grow <- object$forest$prob
  .quantreg.check.prob(prob.grow)
  ## ------------------------------------------------------------------
  ## Restore the training rows or predict new rows using the saved trees.
  ## Resolve probabilities before dispatch so native columns and labels agree.
  ## ------------------------------------------------------------------
  if (!grow || test) {
    if (!isTRUE(object$forest$forest)) {
      stop("Forest information is missing; grow with forest = TRUE.", call. = FALSE)
    }
    ## New rows use ordinary mean predictions (or native GK quantiles).
    ## This does not change the residuals already saved during growth.
    if (test) oob <- FALSE
    dots <- dots.all
    dots$outcome <- "train"
    dots$gk.quantile <- method == "gk"
    settings <- global.prob.assign(
      if (is.null(prob)) prob.grow else prob,
      if (is.null(prob.epsilon)) object$forest$prob.epsilon else prob.epsilon,
      method == "gk", TRUE, object$forest$splitrule, object$forest$n
    )
    prob.grow <- dots$prob <- settings$prob
    dots$prob.epsilon <- settings$prob.epsilon
    if (method == "forest") dots$forest.wt <- if (oob) "oob" else TRUE
    ## The general predictor allows classification-only target requests; the
    ## quantile interface requires at least one requested continuous response.
    if (!is.null(dots$m.target)) {
      continuous <- if (object$family == "regr") object$yvar.names else names(object$regrOutput)
      if (!length(intersect(dots$m.target, continuous))) {
        stop("'m.target' must select at least one continuous outcome.", call. = FALSE)
      }
    }
    ## The predictor carries the saved forest, including training data,
    ## residual banks and their source settings, into the result.
    args <- c(list(object = object), dots)
    if (test) args$newdata <- newdata
    object <- do.call(predict, args)
  }
  prob <- prob.grow
  ## The prediction workhorse leaves NULL slots for unrequested multivariate
  ## outcomes. Process only the continuous output components actually returned.
  if (object$family == "regr") {
    ynames <- object$yvar.names
  } else {
    ynames <- names(object$regrOutput)[!vapply(object$regrOutput, is.null, logical(1))]
  }
  if (!length(ynames)) {
    stop("No continuous-outcome predictions were returned.", call. = FALSE)
  }
  ## ------------------------------------------------------------------
  ## Summarize each response distribution. Quantiles and residual moments
  ## use the full residual bank; yunq is only a CDF/reporting grid.
  ## ------------------------------------------------------------------
  rO <- setNames(vector("list", length(ynames)), ynames)
  for (yn in ynames) {
    component <- if (object$family == "regr") object else object$regrOutput[[yn]]
    y <- .quantreg.response(.quantreg.training.yvar(object), yn)
    grid <- .quantreg.grid(y, maxn)
    selected <- .quantreg.prediction(component, oob, quantile = method == "gk")
    n <- object$n
    mean <- std <- prediction <- residual.source <- NULL
    if (method == "gk") {
      quant <- selected$value
      if (!is.numeric(quant) || length(quant) != n * length(prob) ||
          (!is.null(dim(quant)) &&
           !identical(as.integer(dim(quant)), c(as.integer(n), length(prob))))) {
        stop("Native quantile dimensions do not match the evaluated rows and probabilities.",
             call. = FALSE)
      }
      quant <- matrix(quant, nrow = n, ncol = length(prob))
      cdf <- .quantreg.gk.cdf(quant, prob, grid, max(y[is.finite(y)]))
      ## Keep the public quantile matrix in one location, as before.
      if (object$family == "regr") {
        object$quantile <- object$quantile.oob <- NULL
      } else {
        object$regrOutput[[yn]]$quantile <- NULL
        object$regrOutput[[yn]]$quantile.oob <- NULL
      }
    } else {
      residual.source <- object$forest[["quantreg"]][["residual.source"]][[yn]]
      residual <- .quantreg.training.residual(object, yn, residual.source)
      if (length(residual) != length(y)) {
        stop("Saved responses and residuals are not aligned.", call. = FALSE)
      }
      prediction <- .quantreg.vector(selected$value, n, "Current predictions")
      weight <- if (method == "forest") object$forest.wt else NULL
      if (method == "forest" && is.null(weight)) {
        stop("The forest method requires the requested forest-weight matrix.", call. = FALSE)
      }
      summary <- .quantreg.residual(residual, prediction, prob, grid, weight)
      quant <- summary$quantiles
      cdf <- summary$cdf
      mean <- summary$mean
      std <- summary$std
    }
    ## Keep source metadata separate from the numerical OOB residuals.
    ## oob describes the current prediction; residual.source identifies the
    ## selected training bank. GK uses no residual bank for its quantiles.
    rO[[yn]] <- list(quantiles = quant, prob = prob, cdf = cdf,
                     density = .quantreg.increments(cdf), yunq = grid,
                     mean = mean, std = std, method = method,
                     oob = selected$oob, prediction = prediction,
                     residual.source = residual.source)
  }
  object$quantreg <- if (object$family == "regr") rO[[1L]] else rO
  if (!inherits(object, "quantreg")) class(object) <- c(class(object), "quantreg")
  object
}
quantreg <- quantreg.rfsrc
