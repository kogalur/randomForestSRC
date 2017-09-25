max.subtree.rfsrc <- function(object,
                              max.order = 2,
                              sub.order = FALSE,
                              conservative = FALSE,
                              ...)
{
  ## Incoming parameter checks.  All are fatal.
  if (is.null(object)) stop("Object is empty!")
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rfsrc, grow)' or '(rfsrc, forest)'")
  ## Acquire the forest 
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
  ## verify the max.order option
  max.order <- floor(max.order)
  if (max.order < 0) {
    stop("RFSRC 'max.order' requested for distance order statistic must be an integer greater than zero (0).")
  }
  ## max.order=0 defaults to non-conservative thresholding
  if (max.order == 0) {
    conservative <- FALSE
  }
  ## Maximum depth monitored for nodes at depth counts
  ## TBD Increased from 1000 to 10000 but need a better way to set this value. TBD 
  MAX.DEPTH <- 10000
  ## Count the number of trees in the forest.
  numTree <- length(as.vector(unique(nativeArray$treeID)))
  ## Count the number of parameters in the data set.
  numParm <- length(xvar.names)
  ## Determine the sample size of the data set.
  numSamp <- nrow(object$xvar)
  ## ------------------------------------------------------------------
  ## The maximal distance from the root node of a tree or sub-tree is
  ## zero-based.  This means that the root node represents a distance
  ## of zero (0).
  ##
  ## $count   - of length [numParm].
  ##          - integer count of the number of maximal subtrees for each
  ##            xvar in each tree.
  ##          - this will be zero (0) if the xvar is not used in
  ##            the forest.
  ## $order   - if max.order > 0, matrix of dim [numParm] x [max.order]
  ##            representing the order statistic for the maximal
  ##            subtree distance for each xvar
  ##          - if max.order == 0, matrix of [numParm] x [numTree]
  ##            representing the first order statistic for the maximal
  ##            subtree distance for each xvar by tree
  ##          - rows represent the distribution for the order
  ##            statistic for a xvar
  ##          - thus [2][1] gives the minimum distance for xvar 2.
  ##          - elements will have the value (maximum + 1), the
  ##            penalty, if the xvar is not used in the tree, or
  ##            the statistic does not exist
  ## $meanSum - vector of length [numParm]
  ##          - the pre-mean mimimum maximal subtree distance for each
  ##            xvar in each tree.
  ##          - if the xvar is not used in the tree, this will be
  ##            the maximum depth for the tree, the penalty.
  ## $depth   - integer depth (defined as the maximum terminal node
  ##            distance from the root node) possible in the tree
  ##            topology.  This will be always greater than zero (0)
  ##            for non-trivial trees
  ## $terminalDepthSum
  ##          - integer sum of the depths of the terminal
  ##            nodes of the tree
  ##          - used in calculating the average depth of the
  ##            nodes for tree
  ## $subOrder
  ##          - matrix of [numParm] x [numParm] representing the
  ##            normalized minimum maximal w-subtree distance for
  ##            parameter v's maximal subtree.
  ##          - thus if v=2, w=1, consider a maximal 2-subtree with
  ##            the "root" node of this subtree associated with
  ##            relative depth = 0.  Then consider a 1-subtree of the
  ##            2-subtree.  Then, [2][1] gives the normalized minimum
  ##            maximal subtree distance of all 1-subtrees within the
  ##            2-subtree.  If the maximal w-subtree does not exist,
  ##            the distance is set to be one (1).  The normalized
  ##            distance is calculated as follows.  First, the
  ##            maximum depth of the terminal nodes of the v-subtree
  ##            is determined.  Then the relative depth of the
  ##            w-subtree (from the root node) of the v-subtree is
  ##            determined.  This latter quantity is divided by the
  ##            prior maximum terminal node depth quantity to give a
  ##            number between (0,1].
  ##          - the diagonal [i][i] is the normalized minimum maximal
  ##            v-subtree distance for the parameter i.  If it is not
  ##            split on in the tree, it is assigned the value one
  ##            (1).
  ##          - the resulting matrix is NA-free.
  ## $subOrderDiag
  ##          - interim vector of length [numParm] representing the 
  ##            minimum maximal v-subtree distance for the tree that
  ##            will constitute the diagonal of $subOrder.
  ##          - variables that are not split on are assigned the
  ##            penalty value of the maximum depth.
  ## $nodesAtDepth
  ##          - number of nodes at each depth, where [1] = 2 for
  ##            a non-trivial tree, [2] = 4 for a balanced tree.
  ##            In general, for a balanced tree, at depth d, the
  ##            count is 2^d.
  ##
  ## --------------------------------------------------------------
  ## Loop through all trees
  ## This is a LOCAL operation
  ## Global dependencies:  (xvar.names, forest)
  ##
  ## It is inefficient to use a do-loop, therefore we use lapply
  ## or mclapply.  Later we parse the list to extract the necessary
  ## objects
  ## --------------------------------------------------------------
  subtree.obj <- mclapply(1:numTree, function(b) {
    ## Create the (local) subtree object for recursion.  This is NOT the
    ## output object, but is closely related to it.
    subtree <- vector("list", 8)
    names(subtree) <- c("count",
                        "order",
                        "meanSum",
                        "depth",
                        "terminalDepthSum",
                        "subOrder",
                        "subOrderDiag",
                        "nodesAtDepth")
    ## Create the recursive output object.  This would be unnecessary 
    ## if it was possible to declare global variables in a package.
    recursiveObject <- list(offset = min(which(nativeArray$treeID == b)),
                            subtree = subtree,
                            diagnostic = 0,
                            diagnostic2 = 0)
    ## Reset the nodes at depth count.
    recursiveObject$subtree$nodesAtDepth <- rep(NA, MAX.DEPTH)
    ## Reset the mean maximal subtree distance.
    ## This will be NA if the xvar is not used in this tree.
    recursiveObject$subtree$meanSum <- rep(NA, numParm)
    ## Reset the order statistic matrix representing maximal distance.
    recursiveObject$subtree$order <- matrix(NA, nrow=numParm, ncol=max(max.order, 1))
    if (sub.order) {
      ## Reset the sub-order statistic matrix representing the minimum
      ## maximal w-subtree distance for v != w.  
      recursiveObject$subtree$subOrder <- matrix(1.0, nrow=numParm, ncol=numParm)
      ## Reset the subOrder diagonals to NA.
      recursiveObject$subtree$subOrderDiag <- rep(NA, numParm)
    }
    ## Reset the maximum maximal subtree distance.
    recursiveObject$subtree$depth <- 0
    ## Reset the terminal node depth sums.
    recursiveObject$subtree$terminalDepthSum <- 0
    ## Reset the count of the maximal subtrees in the tree.
    recursiveObject$subtree$count <- rep(0, numParm)
    ## Identify the root node split parameter.
    rootParmID <- nativeArray$parmID[recursiveObject$offset] 
    ## Save the previous value of the offset.  This is used in
    ## determining the number of terminal nodes in the tree.
    offsetMark <- recursiveObject$offset
    ## Set the stump count.  Gets updated when the tree is a stump.
    stumpCnt <- 0
    ## ----------------------------
    ##  PRIMARY RECURSION PRIMARY
    ## ----------------------------
    ## Recursively parse the tree in the primary protocol.
    recursiveObject <- rfsrcParseTree(
      recursiveObject,
      max(max.order, 1),
      sub.order,
      nativeArray,
      b,
      distance=0,
      subtreeFlag=rep(FALSE, numParm))
    ## Check that the current tree is not a stump.
    if (rootParmID != 0) {
      ## Update the forest sum ...
      ## Determine which xvar _were_not_ split on, in the current tree.
      index <- which(recursiveObject$subtree$count == 0)
      ## Penalize these unused xvar by making their mean distance the worst.
      recursiveObject$subtree$meanSum[index] <- recursiveObject$subtree$depth
      forestMeanSum <- recursiveObject$subtree$meanSum
      ## This is tree specific.  Here we view the matrix as a vector.
      index <- which(is.na(recursiveObject$subtree$order))
      ## Penalize these order statistics making their values the worst.      
      recursiveObject$subtree$order[index] <- recursiveObject$subtree$depth
      orderTree <- recursiveObject$subtree$order
      ## Explanation of demominator: The total number of nodes N
      ## (internal and external) in a tree is given by the number of
      ## records in a tree.  In this case, offset - offsetMark = N,
      ## due to the nature of the pointer incrementation.  The number
      ## of external (terminal) nodes N(TERM) is given by (N + 1) / 2.
      ## The maximum number of subtrees for the tree is given by
      ## N(TERM) / 2 which is the denominator.
      ## We normalize the actual number of maximal subtrees found in
      ## the tree by the maximum number maximal subtrees possible.
      subtreeCountSum <- (recursiveObject$subtree$count / ((recursiveObject$offset - offsetMark + 1) / 4))
      ## Get the average depth of the terminal nodes of this tree.
      ## The denominator is just the number of terminal nodes.
      terminalDepth <- recursiveObject$subtree$terminalDepthSum / ((recursiveObject$offset - offsetMark + 1) / 2)
      if (sub.order) {
        ## Determine which xvar _were_ split on, in the current tree.
        index <- which(recursiveObject$subtree$count > 0)
        ## Over-ride those diagonal elements with the minimum depth encountered.
        diag(recursiveObject$subtree$subOrder)[index] <- recursiveObject$subtree$subOrderDiag[index]
        ## Determine which xvar _were_not_ split on, in the
        ## current tree.  This can also be determined by examining the
        ## diagonal vector, and noting those elements with the value
        ## NA
        index <- which(recursiveObject$subtree$count == 0)
        ## Over-ride those missing diagonal elements with penalty value.
        diag(recursiveObject$subtree$subOrder)[index] <- recursiveObject$subtree$depth
        ## Normalize the depths.
        diag(recursiveObject$subtree$subOrder) <-  diag(recursiveObject$subtree$subOrder) / recursiveObject$subtree$depth
        ## Add the result to the forest sum.
        subOrderSum <- recursiveObject$subtree$subOrder
      }
        else {
          subOrderSum <- NULL
        }
      ## Update the nodes at depth count.
      nodesAtDepthMatrix <- recursiveObject$subtree$nodesAtDepth
    }
      else {
        ## The tree is a stump.  It will be excluded.
        stumpCnt <- 1
        ## All other parameters are set to NULL
        forestMeanSum <- orderTree <- subtreeCountSum <-
          terminalDepth <- subOrderSum <- nodesAtDepthMatrix <- NULL
      }
    ## return the local subtree objects for parsing
    return(list(forestMeanSum = forestMeanSum,
                orderTree = orderTree,
                subtreeCountSum = subtreeCountSum,
                terminalDepth = terminalDepth,
                subOrderSum = subOrderSum,
                stumpCnt = stumpCnt,
                nodesAtDepthMatrix = nodesAtDepthMatrix))
  })
  ## ---------------------------------------------------------------------
  ## Parse the local object
  ## Pre-mean maximal subtree distance by parameter over the forest.
  forestMeanSum <- rep(0, numParm)
  ## This is by tree over the forest.
  ## Minimum maximal subree distance by parameter by tree over the forest.
  order.tree <- matrix(NA, nrow=numParm, ncol=numTree)
  ## Order statistics by parameter of maximal subtree distance.
  if (max.order > 0) {
    orderSum <- matrix(0, nrow=numParm, ncol=max.order)
  }
  ## Normalized maximal subtree count sum over the forest.  Linked to $count.
  subtreeCountSum <- rep(0, numParm)
  ## Average depth of each tree.  Linked to $terminalDepthSum 
  terminalDepth <- rep(0, numTree)
  ## Sub-order interim sum.
  subOrderSum <- matrix(0, nrow=numParm, ncol=numParm)
  ## Final list of count of nodes at depth by tree.
  nodesAtDepthMatrix <- matrix(NA, nrow = MAX.DEPTH, ncol = numTree)
  ## Initialize the stump count
  stumpCnt <- 0
  ## iterate over trees
  for (b in 1:numTree) {
    local.obj <- subtree.obj[[b]]
    ## assignments only apply when the tree is not a stump
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
    ## the tree is a stump
      else {
        stumpCnt <- stumpCnt + 1
      }
  }
  ## ---------------------------------------------------------------------
  ## Output preparation  
  ##
  ##
  ## Prepare the ensemble values for the output object.
  ##
  ## $mean - vector of length [numParm]
  ##       - mean minimum maximal subtree distance by parameter over the
  ##         forest.
  ## $order - if max.order > 0, matrix of dim [numParm] x [max.order]
  ##            representing the order statistic for the maximal
  ##            subtree distance for each xvar over the forest
  ##          - if max.order == 0, matrix of [numParm] x [numTree]
  ##            representing the first order statistic for the maximal
  ##            subtree distance for each xvar by tree
  ## $count - vector of length [numParm]
  ##        - normalized maximal subtree count by parameter over the
  ##          forest
  ##        - normalization is by the number of topologically possible
  ##          subtrees in each tree
  ## $sub.order - matrix of [numParm] x [numParm] representing the
  ##             normalized minimum maximal w-subtree distance for
  ##             parameter v's maximal subtree averaged over all trees.
  ##           - see the tree specific description for details.
  ##           - output if sub.order=TRUE
  ## $terminal - vector of length [numTree]
  ##           - average terminal node depth of each tree.
  ## $nodes.at.depth - matrix of [MAX.DEPTH] x [numTree]
  ##               - number of nodes at each depth, where [1] = 2 for
  ##                 a non-trivial tree, [2] = 4 for a balanced tree.
  ##                 In general, for a balanced tree, at depth d, the
  ##                 count is [d] = 2^d.
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
  ## Precautionary check for all stumps.
  if(numTree != stumpCnt) {
    ## Attach the average terminal node depth to the output object.
    result$terminal <- terminalDepth
    ## Determine the forest average maximal subtree distance.
    result$mean <- forestMeanSum / (numTree - stumpCnt)
    names(result$mean) <- xvar.names
    ## Minimal depth by variable by tree
    minDepthVarTree <- order.tree
    if (max.order > 0) {
      ## Summarize the order statistics containing the maximal subtree distances. 
      result$order <- orderSum / (numTree - stumpCnt)
      rownames(result$order) <- xvar.names
    }
      else {
        result$order <- minDepthVarTree
        rownames(result$order) <- xvar.names
      }
    ## Determine the forest average subtree count.
    result$count <- subtreeCountSum / (numTree - stumpCnt)
    names(result$count) <- xvar.names
    ## Copy the matrix for the number of nodes at each depth.
    result$nodes.at.depth <- nodesAtDepthMatrix
    if (sub.order == TRUE) {      
      ## Calculate the forest average of the minimum maximal w-subtree distances.
      result$sub.order   <- subOrderSum / (numTree - stumpCnt)
      rownames(result$sub.order) <- xvar.names
      colnames(result$sub.order) <- xvar.names
    }
  }
  ## print (recursiveObject$diagnostic)
  ## print (recursiveObject$diagnostic2)
  ## Minimal depth for thresholding weak variables
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
  ## variable selection
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
  ## minimal depth percentile
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
  ## discard values no longer needed
  result$terminal <- result$mean <- NULL
  ## exit; invisible return
  invisible(result)
}
## -----------------------------------------------------------------------
##
##  PRIMARY RECURSION
##
## -----------------------------------------------------------------------
## Recursive function to determine first order maximal distances.
rfsrcParseTree <- function(recursiveObject,
                           max.order,
                           sub.order,
                           nativeArray,
                           b,
                           distance,
                           subtreeFlag) {
  ## Diagnostic count of calls.
  recursiveObject$diagnostic <- recursiveObject$diagnostic + 1
  ## Weak consistency check to ensure that the iteration matches the treeID in the nativeArray record.
  if(b != nativeArray$treeID[recursiveObject$offset]) {
    stop("Invalid nativeArray input record (treeID) at ", recursiveObject$offset, ".  Please contact Technical Support.")
  }
  ## Check if we are at the root node.  Otherwise increment the node at depth count.
  if (distance > 0) {
    ## We are not at the root node.
    if (distance <= length(recursiveObject$subtree$nodesAtDepth)) {    
      if (is.na(recursiveObject$subtree$nodesAtDepth[distance])) {
        ## This is the first split at this depth, so initialize the nodes at depth count.
        recursiveObject$subtree$nodesAtDepth[distance] <- 1
      }
        else {
          ## Increment the count of nodes at this depth.
          recursiveObject$subtree$nodesAtDepth[distance] <- recursiveObject$subtree$nodesAtDepth[distance] + 1
        }
    }
  }
  ## Read the current nativeArray split parameter.
  splitParameter <- nativeArray$parmID[recursiveObject$offset]
  ## Determine whether this is a terminal node.  
  if (splitParameter == 0) {
    terminalFlag <- TRUE
  }
    else if (splitParameter != 0) {
      terminalFlag <- FALSE
    }
  if (!terminalFlag) {
    ## Update the maximal subtree information for the parameter if it
    ## has not already been encountered in the tree.
    if (subtreeFlag[splitParameter] == FALSE) {
      ## Increment the subtree count for the parameter.
      recursiveObject$subtree$count[splitParameter] <- recursiveObject$subtree$count[splitParameter] + 1
      ## Update the (pre) mean distance.
      if (is.na(recursiveObject$subtree$meanSum[splitParameter])) {
        recursiveObject$subtree$meanSum[splitParameter] <- distance
      }
        else {
          recursiveObject$subtree$meanSum[splitParameter] <- recursiveObject$subtree$meanSum[splitParameter] + distance
        }
      if (max.order > 0) {
        ## Update distance encountered for the order statistic.
        ## The last element, NA, in the temporary vector is a filler to handle which.max() silliness. 
        orderVector <- c(recursiveObject$subtree$order[splitParameter, ], distance, NA)
        ## Identify the first NA occurence.
        index <- which.max(is.na(orderVector))
        ## The first (relevant) NA is filled with the current distance.
        orderVector[index] <- distance
        ## The vector is then sorted.
        sortedVector <- sort(orderVector[1:index])
        ## Restore the length of the order vector, if necessary, by apppending NA's.
        if (index <= max.order) {
          orderVector <- c(sortedVector, rep(NA, max.order-index))
        }
          else {
            orderVector <- sortedVector[1:max.order]
          }
        ## Update the order statistic.
        recursiveObject$subtree$order[splitParameter, ] <- orderVector
      }
        else {
          ## Update the minimum distance encountered for the parameter.
          if (is.na(recursiveObject$subtree$order[splitParameter])) {
            recursiveObject$subtree$order[splitParameter] <- distance
          }
            else {
              recursiveObject$subtree$order[splitParameter] <- min(recursiveObject$order[splitParameter], distance, na.rm = TRUE)
            }
        }
      ## Indicate that the parameter has been split on.
      subtreeFlag[splitParameter] <- TRUE
      if (sub.order == TRUE) {
        ## Update the diagonal element with the minimum distance encountered for the parameter.
        if (is.na(recursiveObject$subtree$subOrderDiag[splitParameter])) {
          recursiveObject$subtree$subOrderDiag[splitParameter] <- distance
        }
          else {
            recursiveObject$subtree$subOrderDiag[splitParameter] <- min(recursiveObject$subtree$subOrderDiag[splitParameter], distance, na.rm = TRUE)
          }
        ## Create the recursive2 list object.  Note that the relative depth starts at zero (0) for the subtree.  Also note
        ## that the object returns a vector containing the unnormalized minimum maximal w-subtree distance.  
        recursive2Object <- list(offset = recursiveObject$offset,
                                 depth = 0,
                                 minimumVector = rep(NA, dim(recursiveObject$subtree$subOrder)[2]),
                                 diagnostic = recursiveObject$diagnostic2)
        ## Initialize the split flags for the sub-recursion.  Ignore
        ## the diagonal element by setting the split flag.
        subtree2Flag <- rep(FALSE, dim(recursiveObject$subtree$subOrder)[2])        
        subtree2Flag[splitParameter] <- TRUE
        ## Do the sub-recursion.
        recursive2Object <- rfsrcParse2Tree(recursive2Object,
                                            nativeArray,
                                            b,
                                            distance=0,
                                            subtreeFlag=subtree2Flag)
        recursiveObject$diagnostic2 <- recursiveObject$diagnostic2 + recursive2Object$diagnostic
        ## Set the element corresponding to the diagonal element temporarily to the penalized value, to avoid being manipulated.
        recursive2Object$minimumVector[splitParameter] <- recursive2Object$depth
        ## Determine which w-subtrees in this v-subtree do not exist and set them to the penalized value.
        recursive2Object$minimumVector[which(is.na(recursive2Object$minimumVector))] <- recursive2Object$depth
        ## Normalize the relative distances by the maximum terminal node depth encountered for this v-subtree.
        recursive2Object$minimumVector <- recursive2Object$minimumVector / recursive2Object$depth
        ## Update the minimum w-subtree distance encountered in the tree.
        recursiveObject$subtree$subOrder[splitParameter, ] <- pmin(recursiveObject$subtree$subOrder[splitParameter, ], recursive2Object$minimumVector)
      }  ## if (sub.order == TRUE) ...
    }  ## if (subtreeFlag[splitParameter] == FALSE) ...
  }  ## if (terminalFlag and distance > 0) then adjust the node count by - 1
    else if (distance > 0) {
      recursiveObject$subtree$nodesAtDepth[distance] <- recursiveObject$subtree$nodesAtDepth[distance] - 1
    }
  ## Update the maximum depth encountered for this tree.
  recursiveObject$subtree$depth <- max(recursiveObject$subtree$depth, distance)
  ## Increment the offset.
  recursiveObject$offset <- recursiveObject$offset + 1
  ## Parse left and then right, if this is not a terminal node.
  if (terminalFlag == FALSE) {
    ## Increment the (parsed) tree distance. 
    distance <- distance + 1
    ## Parse left:
    recursiveObject <- rfsrcParseTree(recursiveObject, max.order, sub.order, nativeArray, b, distance, subtreeFlag)
    ## Parse right:
    recursiveObject <- rfsrcParseTree(recursiveObject, max.order, sub.order, nativeArray, b, distance, subtreeFlag)
  }
    else {
      ## Update the terminal node depth sum. 
      recursiveObject$subtree$terminalDepthSum <- recursiveObject$subtree$terminalDepthSum + distance
    }
  return(recursiveObject)
}
## -----------------------------------------------------------------------
##
##  SECONDARY RECURSION
##
## -----------------------------------------------------------------------
## Determine the second order maximal subtree distance.
## Recursive Object Definition:
##
## recursiveObject =
##   list(offset,
##        depth = 0,
##        minimumVector = rep(NA, numParm),
##        diagnostic)
##
## Inputs:
##   distance
##     - distance (depth) from the root of the (primary)
##       v-subtree.  Thus, this always starts at zero (0).
##     
##
## Outputs:  
##   recursiveObject$minimumDistance
##     - vector of length [numParm]
##     - second order unnormalized minimum maximum w-subtree distances
##       relative to the subtree defined by the initial root call.
##   recursiveObject$depth
##     - maximum terminal node depth, relative to the v-subtree.
##
## Algorithm Notes:  
##
## We use the terminology (v,w) for a w-maximal subtree within
## (relative to) a v-maximal subtree.
## 
## STEP-1:
## 
## For a given v, and a given v-max subtree, we search for all w-max
## subtrees within this v-max subtree, where v != w.  We single out
## the w-subtree with minimum depth -- that is, the w-split closest to
## the v-split.  Thus, when we talk about depths of w, these are
## relative depths, relative to the v-max subtree.  So the depth at
## root of the v-max subtree is defined to be zero (0).
## 
## For a given v, and a given v-max subtree, we also determine the
## maximum terminal node depth.  Let's call this the penalty depth.
## Again, this is relative to the v-max subtree.
## 
## If a particular w-max subtree does not exist for a given v and
## v-max subtree, it is assigned the penalty depth.
## 
## Now, there are p-1 minimum w-max subtree depths, given by analyzing
## all w != v.  These distances are normalized (divided by) by the
## penalty depth for the v-max subtree, so that all p-1 values are in
## (0,1].
## 
## 
## STEP-2:
## Consider v as fixed.  We do STEP-1 for all v-max subtrees.  Say
## there are count(v) of these.  We then take the MINIMUM normalized depth across the count(v)
## vectors of length p-1 that give the minimum w-max subtree depths by
## dividing by count(v).
## 
## STEP-3:
## Currently, if the v-max subtree does not exist, the row
## corresponding to v in the $subOrder matrix, namely [v, ] will be
## one (1) except for the element corresponding to the diagonal [v,v].
## It will be zero (0).
## 
## STEP-4:
## 
## All of the above are repeated for each tree in the forest such that
## the [p]x[p] matrix for each tree is normalized (divided by) the the
## number of valid (non-stumped) trees in the forest.
## 
## There is the rare possibility that a xvar v will never be
## split upon within the forest.  The entire row for this xvar,
## given by [v, ], will be NA.
rfsrcParse2Tree <- function(recursiveObject,
                            nativeArray,
                            b,
                            distance,
                            subtreeFlag) {
  ## Diagnostic count of calls.
  recursiveObject$diagnostic = recursiveObject$diagnostic + 1
  ## Weak consistency check to ensure that the iteration matches the treeID in the nativeArray record.
  if(b != nativeArray$treeID[recursiveObject$offset]) {
    stop("Invalid nativeArray input record (treeID) at ", recursiveObject$offset, ".  Please contact Technical Support.")
  }
  ## Read the current nativeArray split parameter.
  splitParameter = nativeArray$parmID[recursiveObject$offset]
  ## Determine whether this is a terminal node.  
  if (splitParameter == 0) {
    terminalFlag = TRUE
  }
    else if (splitParameter != 0) {
      terminalFlag = FALSE
    }
  if (splitParameter != 0) {
    ## Check that the parameter has not already been encountered along this path.
    if (subtreeFlag[splitParameter] == FALSE) {
      ## Update the minimum maximal subtree distance encountered for the parameter.
      if (is.na(recursiveObject$minimumVector[splitParameter])) {
        recursiveObject$minimumVector[splitParameter] = distance
      }
        else {
          recursiveObject$minimumVector[splitParameter] = min(recursiveObject$minimumVector[splitParameter], distance, na.rm = TRUE)
        }
      ## Indicate that the parameter has been split on
      subtreeFlag[splitParameter] = TRUE
    }
  }
  ## Update the (relative) maximum depth encountered for this tree.
  recursiveObject$depth = max(recursiveObject$depth, distance)
  ## Increment the (parsed) tree distance. 
  distance = distance + 1
  ## Increment the offset.
  recursiveObject$offset = recursiveObject$offset + 1
  ## Parse left and then right, if this is not a terminal node.
  if (terminalFlag == FALSE) {
    ## Parse left:
    recursiveObject = rfsrcParse2Tree(recursiveObject, nativeArray, b, distance, subtreeFlag)
    ## Parse right:
    recursiveObject = rfsrcParse2Tree(recursiveObject, nativeArray, b, distance, subtreeFlag)
  }
  return(recursiveObject)
}
## -----------------------------------------------------------------------
##
##  Supporting Functions
##
## -----------------------------------------------------------------------
## minimal depth density
minDepthProb <- function(p, D, l) {
  if (!is.null(l)) Ld <- 0
  prob <- rep(0, D+1)
  nullObj <- sapply(0:(D-1), function(d) {##trick to avoid a for-loop
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
## second order depth density
secondOrderDepthProb <- function(p, D, l) {
  md.prob <- rep(0, D+1)
  Ld <- 0
  if (is.null(l)) {
    l <- sapply(0:(D-1), function(d) 2^d)
  }
  nullObj <- sapply(0:(D-1), function(d) {##trick to avoid a for-loop
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
  ## md.prob is the minimal depth density
  ## we now use it to define the second order depth density
  prob <- rep(0, D+1)
  if (D >= 2) {
    nullObj <- sapply(1:(D-1), function(d) {##trick to avoid a for-loop
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
## calculates minimal depth density under the null
## calculates first/second moments and robust estimate of s.d
## treat stumps as equivalent to trees with tree depth = 1.
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
## calculates second order depth density under the null
## returns its first moment
## treat stumps as equivalent to trees with tree depth = 1.
secondOrderDepthStat <- function(p, D=NULL, l=NULL) {
  if (is.null(D) & is.null(l)) stop("set D or l")
  if (!is.null(l)) {
    D <- length(l)
  }
  D.support <- c(0:D)
  prob <- secondOrderDepthProb(p, D=D, l=l)
  sum(D.support * prob, na.rm = TRUE)
}
## parses minimal depth and second order depth
## returns the md density
## returns the md first (*threshold*) and second moment
## returns a robust estimate of standard deviation (gets rid of atom in the right tail)
## returns (new) threshold for second order statistic
## v is the max.subtree object
parseDepthOrder <- function(v, max.order) {
  ## preliminary checks: can this calculation even be done?
  if (is.null(v$mean) & max.order > 0) {
    return(NULL)
  }
  if (is.null(v$order) & max.order == 0) {
    return(NULL)
  }
  ## root node is assigned a value of 1 nodes at depth
  ## although not technically correct, it is relatively benign
  nodes.at.depth <- rbind(1, v$nodes.at.depth)
  ## determine the tree height (height = 1 + depth)
  ## determine the average tree height
  ## determine the maximum tree height
  treeHeight <- apply(nodes.at.depth, 2, function(l) {sum(!is.na(l))})
  avgTreeHeight <- mean(treeHeight, na.rm=TRUE)
  maxTreeHeight <- max(treeHeight, na.rm=TRUE)
  ## determine the number of trees in the forest
  ntree <- length(v$terminal)
  ## mean minimal depth
  ## break this up into two cases
  if (max.order > 0) {
    ## set the dimension
    p <- length(v$mean)
    ## case I: substitute forest averaged values for parameters
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
      ## case II: tree-averaged (by-tree) analysis
      ## set the dimension
      p <- nrow(v$order)
      ## extract the tree density in a matrix format
      ## extract each trees depth
      prob <- sapply(1:ntree, function(tree) {
        l <- nodes.at.depth[, tree][1:sum(nodes.at.depth[, tree] > 0, na.rm = TRUE)]
        c(minDepthStat(p, l = l)$prob, rep(0, maxTreeHeight - treeHeight[tree] * (1 + 1 * (treeHeight[tree] == 1))))
      })
      ## extract the tree minimal threshold value
      first.moment <- unlist(sapply(1:ntree, function(tree) {
        l <- nodes.at.depth[, tree][1:sum(nodes.at.depth[, tree] > 0, na.rm = TRUE)]
        minDepthStat(p, l = l)$first.moment
      }))
      ## extract the tree second moment
      second.moment <- unlist(sapply(1:ntree, function(tree) {
        l <- nodes.at.depth[, tree][1:sum(nodes.at.depth[, tree] > 0, na.rm = TRUE)]
        minDepthStat(p, l = l)$second.moment
      }))
      ## extract the tree robust standard deviation
      sd.robust  <- unlist(sapply(1:ntree, function(tree) {
        l <- nodes.at.depth[, tree][1:sum(nodes.at.depth[, tree] > 0, na.rm = TRUE)]
        minDepthStat(p, l = l)$sd.robust
      }))
      ## extract the tree second order threshold value
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
