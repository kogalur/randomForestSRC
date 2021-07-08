tune.nodesize.rfsrc <- function(formula, data, 
                            nodesizeTry = c(1:9, seq(10, 150, by = 5)), ntreeTry = 100,
                            sampsize = function(x){min(x * .632, max(150, x ^ (4/5)))},
                            nsplit = 1, trace = TRUE, ...) 
{
  ## restrict nodesize to values less than or equal to sampsize / 2
  if (is.function(sampsize)) {
    n <- sampsize(nrow(data)) / 2
  }
  else {
    n <- sampsize / 2
  }
  n <- max(n, 10)
  nodesizeTry <- nodesizeTry[nodesizeTry <= n]
  ## loop over nodesize acquiring the error rate
  err <- sapply(nodesizeTry, function(nsz) {
    ## pull the error rate for each candidate nodesize value
    err.nsz <- tryCatch({mean(get.mv.error(rfsrc.fast(formula, data,
            ntree = ntreeTry, nodesize = nsz,
            sampsize = sampsize, nsplit = nsplit, ...), TRUE), na.rm = TRUE)}, 
            error=function(ex){NA})
    if (trace) {
      cat("nodesize = ", nsz,
          " OOB error =", paste(100 * round(err.nsz, 4), "%", sep = ""), "\n")
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
