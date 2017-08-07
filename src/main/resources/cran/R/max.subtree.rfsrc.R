max.subtree.rfsrc <- function(object,
                              max.order = 2,
                              sub.order = FALSE,
                              conservative = FALSE,
                              ...)
{
  if (is.null(object)) stop("Object is empty!")
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rfsrc, grow)' or '(rfsrc, forest)'")
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.null(object$forest)) 
      stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'.")
    object <- object$forest
  }
  nativeArray <- object$nativeArray
  if (is.null(nativeArray)) {
    stop("RFSRC nativeArray content is NULL.  Please ensure the object is valid.")
  }
  xvar.names <- object$xvar.names
  if (is.null(xvar.names)) {
    stop("RFSRC xvar.names content is NULL.  Please ensure the object is valid.")
  }
  if (is.null(object$xvar)) {
    stop("RFSRC xvar content is NULL.  Please ensure the object is valid.")
  }
  max.order <- floor(max.order)
  if (max.order < 0) {
    stop("RFSRC 'max.order' requested for distance order statistic must be an integer greater than zero (0).")
  }
  if (max.order == 0) {
    conservative <- FALSE
  }
  MAX.DEPTH <- 10000
  numTree <- length(as.vector(unique(nativeArray$treeID)))
  numParm <- length(xvar.names)
  numSamp <- nrow(object$xvar)
  subtree.obj <- mclapply(1:numTree, function(b) {
    subtree <- vector("list", 8)
    names(subtree) <- c("count",
                        "order",
                        "meanSum",
                        "depth",
                        "terminalDepthSum",
                        "subOrder",
                        "subOrderDiag",
                        "nodesAtDepth")
    recursiveObject <- list(offset = min(which(nativeArray$treeID == b)),
                            subtree = subtree,
                            diagnostic = 0,
                            diagnostic2 = 0)
    recursiveObject$subtree$nodesAtDepth <- rep(NA, MAX.DEPTH)
    recursiveObject$subtree$meanSum <- rep(NA, numParm)
    recursiveObject$subtree$order <- matrix(NA, nrow=numParm, ncol=max(max.order, 1))
    if (sub.order) {
      recursiveObject$subtree$subOrder <- matrix(1.0, nrow=numParm, ncol=numParm)
      recursiveObject$subtree$subOrderDiag <- rep(NA, numParm)
    }
    recursiveObject$subtree$depth <- 0
    recursiveObject$subtree$terminalDepthSum <- 0
    recursiveObject$subtree$count <- rep(0, numParm)
    rootParmID <- nativeArray$parmID[recursiveObject$offset] 
    offsetMark <- recursiveObject$offset
    stumpCnt <- 0
    recursiveObject <- rfsrcParseTree(
      recursiveObject,
      max(max.order, 1),
      sub.order,
      nativeArray,
      b,
      distance=0,
      subtreeFlag=rep(FALSE, numParm))
    if (rootParmID != 0) {
      index <- which(recursiveObject$subtree$count == 0)
      recursiveObject$subtree$meanSum[index] <- recursiveObject$subtree$depth
      forestMeanSum <- recursiveObject$subtree$meanSum
      index <- which(is.na(recursiveObject$subtree$order))
      recursiveObject$subtree$order[index] <- recursiveObject$subtree$depth
      orderTree <- recursiveObject$subtree$order
      subtreeCountSum <- (recursiveObject$subtree$count / ((recursiveObject$offset - offsetMark + 1) / 4))
      terminalDepth <- recursiveObject$subtree$terminalDepthSum / ((recursiveObject$offset - offsetMark + 1) / 2)
      if (sub.order) {
        index <- which(recursiveObject$subtree$count > 0)
        diag(recursiveObject$subtree$subOrder)[index] <- recursiveObject$subtree$subOrderDiag[index]
        index <- which(recursiveObject$subtree$count == 0)
        diag(recursiveObject$subtree$subOrder)[index] <- recursiveObject$subtree$depth
        diag(recursiveObject$subtree$subOrder) <-  diag(recursiveObject$subtree$subOrder) / recursiveObject$subtree$depth
        subOrderSum <- recursiveObject$subtree$subOrder
      }
        else {
          subOrderSum <- NULL
        }
      nodesAtDepthMatrix <- recursiveObject$subtree$nodesAtDepth
    }
      else {
        stumpCnt <- 1
        forestMeanSum <- orderTree <- subtreeCountSum <-
          terminalDepth <- subOrderSum <- nodesAtDepthMatrix <- NULL
      }
    return(list(forestMeanSum = forestMeanSum,
                orderTree = orderTree,
                subtreeCountSum = subtreeCountSum,
                terminalDepth = terminalDepth,
                subOrderSum = subOrderSum,
                stumpCnt = stumpCnt,
                nodesAtDepthMatrix = nodesAtDepthMatrix))
  })
  forestMeanSum <- rep(0, numParm)
  order.tree <- matrix(NA, nrow=numParm, ncol=numTree)
  if (max.order > 0) {
    orderSum <- matrix(0, nrow=numParm, ncol=max.order)
  }
  subtreeCountSum <- rep(0, numParm)
  terminalDepth <- rep(0, numTree)
  subOrderSum <- matrix(0, nrow=numParm, ncol=numParm)
  nodesAtDepthMatrix <- matrix(NA, nrow = MAX.DEPTH, ncol = numTree)
  stumpCnt <- 0
  for (b in 1:numTree) {
    local.obj <- subtree.obj[[b]]
    if (local.obj$stumpCnt == 0) {
      forestMeanSum <- forestMeanSum <- forestMeanSum + local.obj$forestMeanSum
      if (max.order > 0) {
        orderSum   <- orderSum + local.obj$orderTree
        order.tree[ , b] <- local.obj$orderTree[, 1]
      }
        else {
          order.tree[ , b] <- local.obj$orderTree
        }
      subtreeCountSum <- subtreeCountSum + local.obj$subtreeCountSum
      terminalDepth[b] <- local.obj$terminalDepth
      if (sub.order == TRUE) {
        subOrderSum <- subOrderSum + local.obj$subOrderSum
      }
      nodesAtDepthMatrix[, b] <- local.obj$nodesAtDepthMatrix
    }
      else {
        stumpCnt <- stumpCnt + 1
      }
  }
  nameVector <- c("mean",
                  "order",
                  "count",
                  "terminal",
                  "nodes.at.depth",
                  "sub.order",
                  "threshold",
                  "threshold.1se",
                  "topvars",
                  "topvars.1se",
                  "percentile",
                  "density")
  result <- vector("list", length(nameVector))
  names(result) <- nameVector
  if(numTree != stumpCnt) {
    result$terminal <- terminalDepth
    result$mean <- forestMeanSum / (numTree - stumpCnt)
    names(result$mean) <- xvar.names
    minDepthVarTree <- order.tree
    if (max.order > 0) {
      result$order <- orderSum / (numTree - stumpCnt)
      rownames(result$order) <- xvar.names
    }
      else {
        result$order <- minDepthVarTree
        rownames(result$order) <- xvar.names
      }
    result$count <- subtreeCountSum / (numTree - stumpCnt)
    names(result$count) <- xvar.names
    result$nodes.at.depth <- nodesAtDepthMatrix
    if (sub.order == TRUE) {      
      result$sub.order   <- subOrderSum / (numTree - stumpCnt)
      rownames(result$sub.order) <- xvar.names
      colnames(result$sub.order) <- xvar.names
    }
  }
  if (conservative == TRUE) {
    parseDepthOrderObj <- parseDepthOrder(result, max.order)
    if (!is.null(parseDepthOrderObj)) {
      result$threshold <- parseDepthOrderObj$first.moment
      result$threshold.1se <- result$threshold + parseDepthOrderObj$sd / sqrt(numSamp)
      result$density <- parseDepthOrderObj$prob
      result$second.order.threshold <- parseDepthOrderObj$second.order.moment
    }
  }
    else {
      parseDepthOrderObj <- parseDepthOrder(result, 0)
      if (!is.null(parseDepthOrderObj)) {
        threshold.con <- parseDepthOrder(result, 1)$first.moment
        threshold.lib <- mean(parseDepthOrderObj$first.moment, na.rm = TRUE)
        result$threshold <- ifelse(threshold.lib - 0.5 > threshold.con, threshold.lib - 0.5, threshold.lib)
        result$threshold.1se <- result$threshold + mean(parseDepthOrderObj$sd, na.rm = TRUE) / sqrt(numSamp)    
        result$density <- apply(parseDepthOrderObj$prob, 1, mean, na.rm = TRUE)
        result$second.order.threshold <- mean(parseDepthOrderObj$second.order.moment, na.rm = TRUE)
      }
    }
  if (!is.null(result$density)) {
    names(result$density) <- paste(0:(length(result$density) - 1))
  }
  if (max.order > 0) {
    minDepthVar <- result$order[, 1]
  }
    else {
      if (!is.null(result$order)) {
        minDepthVar <- apply(minDepthVarTree, 1, mean, na.rm = TRUE)
      }
        else {
          minDepthVar <- NULL
        }
    }
  top.pt <- minDepthVar <= result$threshold
  top.pt[is.na(top.pt)] <- FALSE
  if (sum(top.pt) > 0) {
    result$topvars <- xvar.names[top.pt]
  }
  top.1se.pt <- minDepthVar <= result$threshold.1se
  top.1se.pt[is.na(top.1se.pt)] <- FALSE
  if (sum(top.1se.pt) > 0) {
    result$topvars.1se <- xvar.names[top.1se.pt]
  }
  if (!is.null(result$order)) {
    cdf <- cumsum(result$density)
    if (conservative == TRUE) {
      result$percentile <- cdf[pmax(round(minDepthVar), 1)]      
    }
      else {
        result$percentile <- apply(rbind(sapply(1:ncol(minDepthVarTree), function(b) {
          cdf[pmax(minDepthVarTree[, b], 1)]
        })), 1, mean, na.rm = TRUE)
      }
    names(result$percentile) <- xvar.names
  }
  result$terminal <- result$mean <- NULL
  invisible(result)
}
rfsrcParseTree <- function(recursiveObject,
                           max.order,
                           sub.order,
                           nativeArray,
                           b,
                           distance,
                           subtreeFlag) {
  recursiveObject$diagnostic <- recursiveObject$diagnostic + 1
  if(b != nativeArray$treeID[recursiveObject$offset]) {
    stop("Invalid nativeArray input record (treeID) at ", recursiveObject$offset, ".  Please contact Technical Support.")
  }
  if (distance > 0) {
    if (distance <= length(recursiveObject$subtree$nodesAtDepth)) {    
      if (is.na(recursiveObject$subtree$nodesAtDepth[distance])) {
        recursiveObject$subtree$nodesAtDepth[distance] <- 1
      }
        else {
          recursiveObject$subtree$nodesAtDepth[distance] <- recursiveObject$subtree$nodesAtDepth[distance] + 1
        }
    }
  }
  splitParameter <- nativeArray$parmID[recursiveObject$offset]
  if (splitParameter == 0) {
    terminalFlag <- TRUE
  }
    else if (splitParameter != 0) {
      terminalFlag <- FALSE
    }
  if (!terminalFlag) {
    if (subtreeFlag[splitParameter] == FALSE) {
      recursiveObject$subtree$count[splitParameter] <- recursiveObject$subtree$count[splitParameter] + 1
      if (is.na(recursiveObject$subtree$meanSum[splitParameter])) {
        recursiveObject$subtree$meanSum[splitParameter] <- distance
      }
        else {
          recursiveObject$subtree$meanSum[splitParameter] <- recursiveObject$subtree$meanSum[splitParameter] + distance
        }
      if (max.order > 0) {
        orderVector <- c(recursiveObject$subtree$order[splitParameter, ], distance, NA)
        index <- which.max(is.na(orderVector))
        orderVector[index] <- distance
        sortedVector <- sort(orderVector[1:index])
        if (index <= max.order) {
          orderVector <- c(sortedVector, rep(NA, max.order-index))
        }
          else {
            orderVector <- sortedVector[1:max.order]
          }
        recursiveObject$subtree$order[splitParameter, ] <- orderVector
      }
        else {
          if (is.na(recursiveObject$subtree$order[splitParameter])) {
            recursiveObject$subtree$order[splitParameter] <- distance
          }
            else {
              recursiveObject$subtree$order[splitParameter] <- min(recursiveObject$order[splitParameter], distance, na.rm = TRUE)
            }
        }
      subtreeFlag[splitParameter] <- TRUE
      if (sub.order == TRUE) {
        if (is.na(recursiveObject$subtree$subOrderDiag[splitParameter])) {
          recursiveObject$subtree$subOrderDiag[splitParameter] <- distance
        }
          else {
            recursiveObject$subtree$subOrderDiag[splitParameter] <- min(recursiveObject$subtree$subOrderDiag[splitParameter], distance, na.rm = TRUE)
          }
        recursive2Object <- list(offset = recursiveObject$offset,
                                 depth = 0,
                                 minimumVector = rep(NA, dim(recursiveObject$subtree$subOrder)[2]),
                                 diagnostic = recursiveObject$diagnostic2)
        subtree2Flag <- rep(FALSE, dim(recursiveObject$subtree$subOrder)[2])        
        subtree2Flag[splitParameter] <- TRUE
        recursive2Object <- rfsrcParse2Tree(recursive2Object,
                                            nativeArray,
                                            b,
                                            distance=0,
                                            subtreeFlag=subtree2Flag)
        recursiveObject$diagnostic2 <- recursiveObject$diagnostic2 + recursive2Object$diagnostic
        recursive2Object$minimumVector[splitParameter] <- recursive2Object$depth
        recursive2Object$minimumVector[which(is.na(recursive2Object$minimumVector))] <- recursive2Object$depth
        recursive2Object$minimumVector <- recursive2Object$minimumVector / recursive2Object$depth
        recursiveObject$subtree$subOrder[splitParameter, ] <- pmin(recursiveObject$subtree$subOrder[splitParameter, ], recursive2Object$minimumVector)
      }  
    }  
  }  
    else if (distance > 0) {
      recursiveObject$subtree$nodesAtDepth[distance] <- recursiveObject$subtree$nodesAtDepth[distance] - 1
    }
  recursiveObject$subtree$depth <- max(recursiveObject$subtree$depth, distance)
  recursiveObject$offset <- recursiveObject$offset + 1
  if (terminalFlag == FALSE) {
    distance <- distance + 1
    recursiveObject <- rfsrcParseTree(recursiveObject, max.order, sub.order, nativeArray, b, distance, subtreeFlag)
    recursiveObject <- rfsrcParseTree(recursiveObject, max.order, sub.order, nativeArray, b, distance, subtreeFlag)
  }
    else {
      recursiveObject$subtree$terminalDepthSum <- recursiveObject$subtree$terminalDepthSum + distance
    }
  return(recursiveObject)
}
rfsrcParse2Tree <- function(recursiveObject,
                            nativeArray,
                            b,
                            distance,
                            subtreeFlag) {
  recursiveObject$diagnostic = recursiveObject$diagnostic + 1
  if(b != nativeArray$treeID[recursiveObject$offset]) {
    stop("Invalid nativeArray input record (treeID) at ", recursiveObject$offset, ".  Please contact Technical Support.")
  }
  splitParameter = nativeArray$parmID[recursiveObject$offset]
  if (splitParameter == 0) {
    terminalFlag = TRUE
  }
    else if (splitParameter != 0) {
      terminalFlag = FALSE
    }
  if (splitParameter != 0) {
    if (subtreeFlag[splitParameter] == FALSE) {
      if (is.na(recursiveObject$minimumVector[splitParameter])) {
        recursiveObject$minimumVector[splitParameter] = distance
      }
        else {
          recursiveObject$minimumVector[splitParameter] = min(recursiveObject$minimumVector[splitParameter], distance, na.rm = TRUE)
        }
      subtreeFlag[splitParameter] = TRUE
    }
  }
  recursiveObject$depth = max(recursiveObject$depth, distance)
  distance = distance + 1
  recursiveObject$offset = recursiveObject$offset + 1
  if (terminalFlag == FALSE) {
    recursiveObject = rfsrcParse2Tree(recursiveObject, nativeArray, b, distance, subtreeFlag)
    recursiveObject = rfsrcParse2Tree(recursiveObject, nativeArray, b, distance, subtreeFlag)
  }
  return(recursiveObject)
}
minDepthProb <- function(p, D, l) {
  if (!is.null(l)) Ld <- 0
  prob <- rep(0, D+1)
  nullObj <- sapply(0:(D-1), function(d) {
    if (is.null(l)) {
      Ld <- 2^d-1
      ld <- 2^d
    }
      else{
        ld <- l[d+1]
        if (d > 0) Ld <- Ld + l[d] 
      }
    prob.1 <- Ld * ifelse(p > 1, log(1 - 1 / p), 0)
    prob.2 <- ld * ifelse(p > 1, log(1 - 1 / p), 0)
    prob[d+1] <<- exp(prob.1) * (1 - exp(prob.2))
    NULL
  })
  prob[D+1] <- 1 - sum(prob[1:D])
  if (prob[D+1] < 0) {
    prob[D+1] <- 0
    prob <- prob / sum(prob)
  }
  prob
}
secondOrderDepthProb <- function(p, D, l) {
  md.prob <- rep(0, D+1)
  Ld <- 0
  if (is.null(l)) {
    l <- sapply(0:(D-1), function(d) 2^d)
  }
  nullObj <- sapply(0:(D-1), function(d) {
    ld <- l[d+1]
    if (d > 0) Ld <- Ld + l[d] 
    md.prob.1 <- Ld * ifelse(p > 1, log(1 - 1 / p), 0)
    md.prob.2 <- ld * ifelse(p > 1, log(1 - 1 / p), 0)
    md.prob[d+1] <<- exp(md.prob.1) * (1 - exp(md.prob.2))
    NULL
  })
  md.prob[D+1] <- 1 - sum(md.prob[1:D])
  if (md.prob[D+1] < 0) {
    md.prob[D+1] <- 0
    md.prob <- md.prob / sum(md.prob)
  }
  prob <- rep(0, D+1)
  if (D >= 2) {
    nullObj <- sapply(1:(D-1), function(d) {
      prob.d <- l[1:d] * ifelse(p > 1, log(1 - 1 / p), 0)
      prob[d+1] <<- md.prob[d+1] * sum(exp(-prob.d) - 1)
      NULL
    })
  }
  prob[D+1] <- 1 - sum(prob[1:D])
  if (prob[D+1] < 0) {
    prob[D+1] <- 0
    prob <- prob/sum(prob)
  }
  prob
}
minDepthStat <- function(p, D=NULL, l=NULL) {
  if (is.null(D) & is.null(l)) stop("set D or l")
  if (!is.null(l)) {
    D <- length(l)
  }
  D.support <- c(0:D)
  prob <- minDepthProb(p, D=D, l=l)
  prob.robust <- prob
  prob.robust[length(prob)] <- 0
  prob.robust <- prob.robust / sum(prob.robust, na.rm = TRUE) 
  first.moment <- sum(D.support * prob, na.rm = TRUE)
  second.moment <- sum((D.support^2) * prob, na.rm = TRUE)
  sd.robust <- sqrt(sum((D.support^2) * prob.robust, na.rm = TRUE)
                    - (sum(D.support * prob.robust, na.rm = TRUE))^2)
  return(list(prob = prob, first.moment = first.moment, second.moment = second.moment, sd.robust = sd.robust))
}
secondOrderDepthStat <- function(p, D=NULL, l=NULL) {
  if (is.null(D) & is.null(l)) stop("set D or l")
  if (!is.null(l)) {
    D <- length(l)
  }
  D.support <- c(0:D)
  prob <- secondOrderDepthProb(p, D=D, l=l)
  sum(D.support * prob, na.rm = TRUE)
}
parseDepthOrder <- function(v, max.order) {
  if (is.null(v$mean) & max.order > 0) {
    return(NULL)
  }
  if (is.null(v$order) & max.order == 0) {
    return(NULL)
  }
  nodes.at.depth <- rbind(1, v$nodes.at.depth)
  treeHeight <- apply(nodes.at.depth, 2, function(l) {sum(!is.na(l))})
  avgTreeHeight <- mean(treeHeight, na.rm=TRUE)
  maxTreeHeight <- max(treeHeight, na.rm=TRUE)
  ntree <- length(v$terminal)
  if (max.order > 0) {
    p <- length(v$mean)
    nodes.at.depth.avg <- apply(nodes.at.depth, 1, mean, na.rm = TRUE)
    l <- nodes.at.depth.avg[1:max(avgTreeHeight - 1, 1)]
    minDepthStatObj <- minDepthStat(p, l = l)
    first.moment <- minDepthStatObj$first.moment
    second.moment <- minDepthStatObj$second.moment
    sd.robust <- minDepthStatObj$sd.robust
    prob <- minDepthStatObj$prob
    second.order.moment <- secondOrderDepthStat(p, l = l)
  }
    else {
      p <- nrow(v$order)
      prob <- sapply(1:ntree, function(tree) {
        l <- nodes.at.depth[, tree][1:sum(nodes.at.depth[, tree] > 0, na.rm = TRUE)]
        c(minDepthStat(p, l = l)$prob, rep(0, maxTreeHeight - treeHeight[tree] * (1 + 1 * (treeHeight[tree] == 1))))
      })
      first.moment <- unlist(sapply(1:ntree, function(tree) {
        l <- nodes.at.depth[, tree][1:sum(nodes.at.depth[, tree] > 0, na.rm = TRUE)]
        minDepthStat(p, l = l)$first.moment
      }))
      second.moment <- unlist(sapply(1:ntree, function(tree) {
        l <- nodes.at.depth[, tree][1:sum(nodes.at.depth[, tree] > 0, na.rm = TRUE)]
        minDepthStat(p, l = l)$second.moment
      }))
      sd.robust  <- unlist(sapply(1:ntree, function(tree) {
        l <- nodes.at.depth[, tree][1:sum(nodes.at.depth[, tree] > 0, na.rm = TRUE)]
        minDepthStat(p, l = l)$sd.robust
      }))
      second.order.moment <- unlist(sapply(1:ntree, function(tree) {
        l <- nodes.at.depth[, tree][1:sum(nodes.at.depth[, tree] > 0, na.rm = TRUE)]
        secondOrderDepthStat(p, l = l)
      }))
    }
  return(list(prob = prob,
              first.moment = first.moment,
              sd = sqrt(second.moment - first.moment^2),
              sd.robust = sd.robust,
              second.order.moment = second.order.moment))
}   
max.subtree <- max.subtree.rfsrc
