tune.rfsrc <- function(formula, data,
      mtryStart = ncol(data) / 2,
      nodesizeTry = c(1:9, seq(10, 100, by = 5)), ntreeTry = 100,
      sampsize = function(x){min(x * .632, max(150, x ^ (3/4)))},
      nsplit = 1, stepFactor = 1.25, improve = 1e-3, strikeout = 3, maxIter = 25,
      trace = FALSE, doBest = FALSE, ...) 
{
  ## inital checks - all are fatal
  if (improve < 0) {
    stop("improve must be non-negative.")
  }
  if (stepFactor <= 1) {
    stop("stepFactor must be great than 1.")
  }
  if (missing(formula)) {
    stop("a formula must be supplied (only supervised forests allowed).")
  }
  ## intialize outer loop values
  res <- list()
  counter1 <- 0
  ## loop over nodesize
  for (nsz in nodesizeTry) {
    counter1 <- counter1 + 1
    ## acquire the starting error rate
    o <- rfsrc.fast(formula, data, ntree = ntreeTry, mtry = mtryStart, nodesize = nsz,
                    sampsize = sampsize, nsplit = nsplit, ...)
    mtryStart <- o$mtry
    mtryMax <- length(o$xvar.names)
    errorOld <- mean(get.mv.error(o, TRUE), na.rm = TRUE)
    if (trace) {
      cat("nodesize = ", nsz,
          " mtry =", mtryStart,
          " OOB error =", paste(100 * round(errorOld, 4), "%", sep = ""), "\n")
    }
    ## terminate if we are stumping
    if (mean(o$leaf.count) <= 2) {
      break
    }
    ## intialize inner loop values
    oobError <- list()
    oobError[[1]] <- errorOld
    names(oobError)[1] <- mtryStart
    ## search over mtry
    for (direction in c("left", "right")) {
      if (trace) 
        cat("Searching", direction, "...\n")
      Improve <- 1.1 * improve
      mtryBest <- mtryStart
      mtryCur <- mtryStart
      counter2 <- 1
      strikes <- 0  
      while (counter2 <= maxIter && (Improve >= improve || (Improve < 0 && strikes < strikeout))) {
        counter2 <- counter2 + 1
        if (Improve < 0) {
          strikes <- strikes + 1
        }
        mtryOld <- mtryCur
        if (direction == "left") {
          mtryCur <- max(1, min(ceiling(mtryCur / stepFactor), mtryCur - 1)) 
        }
        else {
          mtryCur <- min(mtryMax, max(floor(mtryCur * stepFactor), mtryCur + 1))
        }
        if (mtryCur == mtryOld) {
          break
        }
        errorCur <- mean(get.mv.error(rfsrc.fast(formula, data, ntree = ntreeTry, mtry = mtryCur,
                     nodesize = nsz, sampsize = sampsize, nsplit = nsplit, ...), TRUE), na.rm = TRUE)
        if (trace) {
          cat("nodesize = ", nsz,
              " mtry =", mtryCur,
              " OOB error =", paste(100 * round(errorCur, 4), "%", sep = ""), "\n")
        }
        oobError[[as.character(mtryCur)]] <- errorCur
        Improve <- 1 - errorCur / errorOld
        if (trace) {
          cat(Improve, improve, "\n")
        }
        if (Improve > improve) {
          errorOld <- errorCur
          mtryBest <- mtryCur
        }
      }
    }## finished inner loop
    ## save the results by nodesize
    mtry <- sort(as.numeric(names(oobError)))
    err <- unlist(oobError[as.character(mtry)])
    res[[counter1]] <- cbind(nodesize = nsz, mtry = mtry, err = err)
  }## finished outer loop
  ## check that results are not NULL
  if (is.null(res)) {
    stop("NULL results - something is wrong, check parameter (tuning) settings\n")
  }
  ## deparse the results and sort along nodesize and mtry
  res <- do.call(rbind, res)
  res <- res[order(res[, 1], res[, 2]), ]
  rownames(res) <- 1:nrow(res)
  ## determine the optimal value
  opt.idx <- which.min(res[, 3])
  ## fit the optimized forest?
  rf <- NULL
  if (doBest) {
    rf <- rfsrc.fast(formula, data, mtry = res[opt.idx, 2], nodesize = res[opt.idx, 1],
                     sampsize = sampsize, nsplit = nsplit)
  }
  ## return the goodies
  list(results = res, optimal = res[opt.idx, -3], rf = rf)
}
tune <- tune.rfsrc
