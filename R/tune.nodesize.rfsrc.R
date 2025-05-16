tune.nodesize.rfsrc <- function(formula, data, 
                            nodesizeTry = c(1:9, seq(10, 150, by = 5)), ntreeTry = 100,
                            sampsize = function(x){min(x * .632, max(150, x ^ (4/5)))},
                            nsplit = 1, trace = TRUE, ...) 
{
  ## re-define the original data in case there are missing values
  stump <- rfsrc(formula, data, nodedepth = 0, perf.type = "none", save.memory = TRUE, ntree = 1, splitrule = "random")
  n <- stump$n
  yvar.names <- stump$yvar.names
  data <- data.frame(stump$yvar, stump$xvar)
  colnames(data)[1:length(yvar.names)] <- yvar.names
  rm(stump)
  ## sample size 
  if (is.function(sampsize)) {
    ssize <- sampsize(n)
  }
  else {
    ssize <- sampsize
  }
  ## set performance type
  if (is.null(list(...)$perf.type)) {
    perf.type <- "default"
  }
  else {
    perf.type <- list(...)$perf.type
  }
  ## now hold out a test data set equal to the tree sample size (if possible)
  if ((2 * ssize)  < n)  {
    tst <- sample(1:n, size = ssize, replace = FALSE)
    trn <- setdiff(1:n, tst)
    newdata <- data[tst,, drop = FALSE]
  }
  else {
    trn <- 1:n
    newdata <- NULL
  }
  ## restrict nodesize to values less than or equal to sampsize / 2
  nodesizeTry <- nodesizeTry[nodesizeTry <= max(10, ssize / 2)]
  ## loop over nodesize acquiring the error rate
  err <- sapply(nodesizeTry, function(nsz) {
    ## pull the error rate for each candidate nodesize value
    if (is.null(newdata)) {
      err.nsz <- tryCatch({mean(get.mv.error(rfsrc.fast(formula, data,
            ntree = ntreeTry, nodesize = nsz,
            sampsize = sampsize, nsplit = nsplit, ...), TRUE), na.rm = TRUE)}, 
            error=function(ex){NA})
    }
    else {
     err.nsz <- tryCatch({mean(get.mv.error(predict(rfsrc.fast(formula, data[trn,, drop = FALSE],
            ntree = ntreeTry, nodesize = nsz,
            sampsize = sampsize, nsplit = nsplit, forest = TRUE,
            perf.type="none", save.memory = FALSE, ...), newdata, perf.type = perf.type), TRUE), na.rm = TRUE)}, 
            error=function(ex){NA})
    }
    ## verbose output
    if (trace) {
      cat("nodesize = ", nsz,
          "   error =", paste(100 * round(err.nsz, 4), "%", sep = ""), "\n")
    }
    err.nsz 
  })
  ## is there OOB error?
  if (all(is.na(err))) {
    warning("OOB error is NA: check forest settings, especially sampsize")
    return(data.frame(nodesize = nodesizeTry, err = err))
  }
  ## identify the optimal nodesize
  bestidx <- which.min(err)
  if (length(bestidx) > 0) {
    nsize.opt <- nodesizeTry[bestidx]
  }
  if (trace) {
    cat("optimal nodesize:", nsize.opt, "\n")
  }
  return(list(nsize.opt = nsize.opt,
              err = data.frame(nodesize = nodesizeTry,
              err = err)))
}
tune.nodesize <- tune.nodesize.rfsrc
