stat.split.rfsrc <- function(object, ...)
{
  if (sum(inherits(object, c("rfsrc", "synthetic"), TRUE) == c(1, 2)) == 2) {
    object <- object$rfSyn
  }
  if (is.null(object)) stop("Object is empty!")
  if ((sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2) &&
      (sum(inherits(object, c("rfsrc", "predict"), TRUE) == c(1, 2)) != 2))
    stop("This function only works for objects of class `(rfsrc, grow)' or `(rfsrc, pred)',")
  inbag <- object$inbag
  if (is.null(object$node.stats)) {    
    stop("RF-SRC statistics content is NULL.  Please re-run grow call with 'statistics=TRUE'")
  }
    else {
      extendedNativeArray <- object$node.stats
    }
  if (is.null(object$pstn.membership)) {
    prune <- FALSE
  }
    else {
      prune <- TRUE
    }
  if (is.null(object$forest)) {
    stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'.")
  }
  object <- object$forest  
  extendedNativeArray <- cbind(object$nativeArray, extendedNativeArray)
  xvar.names <- object$xvar.names
  if (is.null(xvar.names)) {
    stop("RFSRC xvar.names content is NULL.  Please ensure the object is valid.")
  }
  if (is.null(object$xvar)) {
    stop("RFSRC xvar content is NULL.  Please ensure the object is valid.")
  }
  numTree <- length(as.vector(unique(extendedNativeArray$treeID)))
  numParm <- length(xvar.names)
  xvar <- object$xvar
  splitTreeObj <- mclapply(1:numTree, function(b) {
    inbag.local <- inbag[, b]
    replicates <-  NULL
    while (sum(inbag.local) > 0) {
      replicates <- c(replicates, which(inbag.local > 0))
      inbag.local <-  sapply(1:length(inbag.local), function(k) { max(inbag.local[k] - 1, 0)})
    }
    splitInfo <- vector("list", numParm)
    names(splitInfo) <- xvar.names
    recursiveObject <- list(offset     = min(which(extendedNativeArray$treeID == b)),
                            splitInfo  = splitInfo,
                            terminal   = FALSE)
    return(spParseTree(recursiveObject, extendedNativeArray, xvar,
                       b,
                       depth=0,
                       membership = replicates,
                       ptnFlag = FALSE,
                       prune = prune))
  })
  result <- vector("list", numTree)
  for (b in 1:numTree) {
    result[[b]] <- splitTreeObj[[b]]$splitInfo
  }
  invisible(result)
}
spParseTree <- function(recursiveObject, extendedNativeArray, xvar,
                        b,
                        depth,
                        membership,
                        ptnFlag,
                        prune = prune) {
  if(b != extendedNativeArray$treeID[recursiveObject$offset]) {
    stop("Invalid nativeArray input record (treeID) at ", recursiveObject$offset, ".  Please contact Technical Support.")
  }
  treeID <- extendedNativeArray$treeID[recursiveObject$offset]
  nodeID <- extendedNativeArray$nodeID[recursiveObject$offset]  
  parmID <- extendedNativeArray$parmID[recursiveObject$offset]
  mwcpSZ <- extendedNativeArray$mwcpSZ[recursiveObject$offset]
  contPT <- extendedNativeArray$contPT[recursiveObject$offset]
  spltST <- extendedNativeArray$spltST[recursiveObject$offset]
  if (parmID == 0) {
    recursiveObject$terminal <- TRUE    
  }
    else {
      recursiveObject$terminal <- FALSE
    }
  if (prune) {
    if (ptnFlag == FALSE) {
      if (spltST == 1) {
        ptnFlag = TRUE
      }
    }
  }
  updateFlag <- TRUE
  if ((parmID == 0) || (ptnFlag)) {
    updateFlag <- FALSE    
  }
  if (updateFlag) {
    if(mwcpSZ != 0) {
    }
    split.vec <- sort(unique(xvar[unique(membership), parmID]))
    split.idx <- which(split.vec == contPT)
    spltEC <-  min(length(split.vec) - 1 - split.idx, split.idx - 1) / (length(split.vec) - 1) 
    left.membership <- membership[which(xvar[membership, parmID] <= contPT)]
    rght.membership <- membership[which(xvar[membership, parmID] >  contPT)]    
    localInfo <- c(treeID, nodeID, parmID, contPT, mwcpSZ, depth, 0, spltEC, spltST)
    if (is.null(recursiveObject$splitInfo[[parmID]])) {
      recursiveObject$splitInfo[[parmID]] <- rbind(recursiveObject$splitInfo[[parmID]], localInfo, deparse.level = 0)
      colnames(recursiveObject$splitInfo[[parmID]]) <- c("treeID", "nodeID", "parmID", "contPT", "mwcpSZ", "dpthID",  "spltTY", "spltEC", "spltST")  
    }
      else {
        recursiveObject$splitInfo[[parmID]] <- rbind(recursiveObject$splitInfo[[parmID]], localInfo, deparse.level = 0)
      }
    col.idx <-  which(colnames(recursiveObject$splitInfo[[parmID]]) == "spltTY")    
    row.idx <-  dim(recursiveObject$splitInfo[[parmID]])[1]
  }
  recursiveObject$offset <- recursiveObject$offset + 1
  if (!recursiveObject$terminal) {
    depth <- depth + 1
    recursiveObject <- spParseTree(recursiveObject, extendedNativeArray, xvar,
                                   b, depth, left.membership, ptnFlag, prune)
    if(updateFlag) {    
      split.type <- recursiveObject$splitInfo[[parmID]][row.idx, ]
      split.type[col.idx] <- split.type[col.idx] + 2 
      recursiveObject$splitInfo[[parmID]][row.idx, ] = split.type
      recursiveObject$terminal <-  FALSE          
    }
    recursiveObject <- spParseTree(recursiveObject, extendedNativeArray, xvar,
                                   b, depth, rght.membership, ptnFlag, prune)
    if(updateFlag) {    
      split.type <- recursiveObject$splitInfo[[parmID]][row.idx, ]
      split.type[col.idx] <- split.type[col.idx] + 1 
      recursiveObject$splitInfo[[parmID]][row.idx, ] = split.type
      recursiveObject$terminal <-  FALSE          
    }
  }
  return(recursiveObject)
}
stat.split <- stat.split.rfsrc
