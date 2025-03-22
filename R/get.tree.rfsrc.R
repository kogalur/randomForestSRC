get.tree.rfsrc <- function(object,
                           tree.id,
                           target,
                           m.target = NULL,
                           time,
                           surv.type = c("mort", "rel.freq", "surv", "years.lost", "cif", "chf"),
                           class.type = c("bayes", "rfq", "prob"),
                           ensemble = FALSE,
                           oob = TRUE,
                           show.plots = TRUE,
                           do.trace = FALSE)
{
  ##----------------------------------------------------------------
  ##
  ## tolerance value for dealing with float issues related to
  ## splitting left at an exact split value
  ##
  ##----------------------------------------------------------------
  tolerance <- sqrt(.Machine$double.eps)
  ##----------------------------------------------------------------
  ##
  ## The following two utilities were copied from the BMS CRAN
  ## package.  Thanks to Martin Feldkircher and Stefan Zeugne for
  ## these little quickies.
  ##
  ##----------------------------------------------------------------
  hex2bin <- function (hexcode) 
  {
      if (!is.character(hexcode)) 
          stop("please input a character like '0af34c'")
      hexcode <- paste("0", tolower(hexcode), sep = "")
      hexobj <- .hexcode.binvec.convert(length(hexcode) * 16L)
      return(hexobj$as.binvec(hexcode))
  }
  .hexcode.binvec.convert <- function (length.of.binvec) 
  {
      if (length(length.of.binvec) > 1) 
          length.of.binvec = length(length.of.binvec)
      addpositions = 4 - length.of.binvec%%4
      positionsby4 = (length.of.binvec + addpositions)/4
      hexvec = c(0:9, "a", "b", "c", "d", "e", "f")
      hexcodelist = list(`0` = numeric(4),
                         `1` = c(0, 0, 0, 1), 
                         `2` = c(0, 0, 1, 0),
                         `3` = c(0, 0, 1, 1),
                         `4` = c(0, 1, 0, 0),
                         `5` = c(0, 1, 0, 1),
                         `6` = c(0, 1, 1, 0), 
                         `7` = c(0, 1, 1, 1),
                         `8` = c(1, 0, 0, 0),
                         `9` = c(1, 0, 0, 1),
                         a   = c(1, 0, 1, 0),
                         b   = c(1, 0, 1, 1),
                         c   = c(1, 1, 0, 0),
                         d   = c(1, 1, 0, 1),
                         e   = c(1, 1, 1, 0),
                         f   = c(1, 1, 1, 1))
      return(list(
          as.hexcode = function(binvec) {
              incl = c(numeric(addpositions), binvec)
              dim(incl) = c(4, positionsby4)
              return(paste(hexvec[crossprod(incl, 2L^(3:0)) + 1], collapse = ""))
          },
          as.binvec = function(hexcode) {
              return(unlist(
                  hexcodelist[unlist(strsplit(hexcode, "", fixed = TRUE),
                                     recursive = FALSE, use.names = FALSE)], 
                  recursive = FALSE, use.names = FALSE)[-(1:addpositions)])
          }))
  }
  ##----------------------------------------------------------------
  ##
  ## coherence checks
  ##
  ##----------------------------------------------------------------
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2) {
    stop("This function only works for objects of class `(rfsrc, grow)'")
  }
  if (is.forest.missing(object)) {
    stop("forest is missing.  Re-run rfsrc (grow call) with forest=TRUE")
  }
  if (inherits(object, "anonymous")) {
    ## anonymous <- TRUE
    stop("get.tree does not currently work with anonymous forests\n")
  }
  else {
    anonymous <- FALSE
  }
  ## only one tree allowed
  tree.id <- tree.id[1]
  ## must be integer from 1...ntree
  if (tree.id < 1 || tree.id > object$ntree) {
    stop("tree id must be integer from 1 to ntree")
  }
  ## ensure coherency of the multivariate target
  m.target <- get.univariate.target(object, m.target)
  ##----------------------------------------------------------------
  ##
  ##
  ## take care of some preliminary processing for prediction
  ## this will be used for labels on terminal nodes
  ##
  ##
  ##----------------------------------------------------------------
  ## process the object depending on the underlying family
  family <- object$forest$family
  ## coherence for prediction
  predict.flag <- TRUE
  if (family == "unsupv" || family == "surv-TDC") {
    predict.flag <- FALSE
  }
  if (predict.flag) {
    ## survival and competing risk families.
    if (grepl("surv", family)) {
      ## extract event information
      event.info <- object$event.info
      cens <- event.info$cens
      event.type <- event.info$event.type
      ## assign time if missing
      if (missing(time)) {
        time <- median(event.info$time.interest, na.rm = TRUE)
      }
      ## check for single time point
      if (length(time) > 1) {
        stop("time must be a single value:  ", time)
      }
      ## competing risk
      if (family == "surv-CR") {
        if (missing(target)) {
          target <- 1
        }
        else {
          if (target < 1 || target > max(event.type, na.rm = TRUE)) {
            stop("'target' is specified incorrectly")
          }
        }
        ## set the surv.type
        surv.type <- setdiff(surv.type, c("mort", "rel.freq", "surv"))[1]
        pred.type <- match.arg(surv.type, c("years.lost", "cif", "chf"))
      }
      ## survival
      else {
        target <- 1
        ## set the surv.type
        surv.type <- setdiff(surv.type, c("years.lost", "cif", "chf"))[1]
        pred.type <- match.arg(surv.type, c("rel.freq", "mort", "chf", "surv"))
      }
    }
    ## univariate and multivariate families.
    else {
      ## assign a null time value
      event.info <- time <- NULL
      ## factor outcome
      if(is.factor(coerce.multivariate(object, m.target)$yvar)) {
        object.yvar.levels <- levels(coerce.multivariate(object, m.target)$yvar)
        pred.type <- match.arg(class.type, c("bayes", "rfq", "prob"))
        if (missing(target)) {
          target <- object.yvar.levels[1]
        }
        if (is.character(target)) {
          target <- match(match.arg(target, object.yvar.levels), object.yvar.levels)
        }
        else {
          if ((target > length(object.yvar.levels)) | (target < 1)) {
            stop("target is specified incorrectly:", target)
          }
        }
      }
      ## regression or unsupervised
      else {
        pred.type <- "y"
        target <- NULL
      }
    }
  }
  ##----------------------------------------------------------------
  ##
  ## define the target subset of cases
  ## if ensemble=TRUE  --> this is the entire data
  ## if ensemble=FALSE --> this is inbag or oob cases
  ##
  ##----------------------------------------------------------------
  if (ensemble) {
    subset <- 1:object$n
  }
  ## restrict the object to the target
  ## oob not allowed
  else {
    object <- predict(object, get.tree = tree.id, membership = TRUE, do.trace = do.trace)
    subset <- which(object$inbag[, tree.id] != 0)
    oob <- FALSE
  }
  ##----------------------------------------------------------------
  ##
  ## extract x data for the target subset
  ## convert the data to numeric mode, apply the na.action protocol
  ## missing data not allowed
  ##
  ##----------------------------------------------------------------
  xvar.names <- object$forest$xvar.names
  xvar.factor <- object$forest$xvar.factor
  if (!anonymous) {
    x.data <- object$forest$xvar
    if (any(is.na(x.data))) {
      stop("missing data not allowed")
    }
    x.data <- finalizeData(xvar.names, x.data, miss.flag = FALSE)
    x.data <- x.data[subset,, drop = FALSE]
  }
  ##----------------------------------------------------------------
  ##
  ## now acquire the predicted values for the tree labels
  ##
  ##----------------------------------------------------------------
  if (predict.flag) {
    if (ensemble) {
      object <- predict.rfsrc(object, m.target = m.target, do.trace = do.trace)
    }
    yhat <- extract.pred(object, pred.type, subset, time, m.target, target, oob = oob)
  }
  ##----------------------------------------------------------------
  ##
  ## get tree data
  ##
  ##----------------------------------------------------------------
  ## pull the arrays
  native.array <- object$forest$nativeArray
  native.f.array <- object$forest$nativeFactorArray[[1]]
  ## added processing needed for factors
  f.ctr <- 0
  factor.flag <- FALSE
  if (!is.null(native.f.array)) {
    pt.f <- which(native.array$mwcpSZ != 0)
    factPT <- lapply(pt.f, function(j) {
      f.ctr <<- f.ctr + 1
      step <- native.array$mwcpSZ[j] - 1
      mwcpPT <- native.f.array[f.ctr:(f.ctr+step)]
      mwcpPT <- paste0(sapply(mwcpPT, function(mwc) {
        format(as.hexmode(mwcpPT), 8)
      }))      
      mwcpSZ <- hex2bin(mwcpPT)
      paste(mwcpSZ, collapse = "")
    })
    native.array$contPT[pt.f] <- factPT
    factor.flag <- TRUE
  }
  ## define the display tree
  display.tree <- native.array[native.array$treeID == tree.id,, drop = FALSE]
  ## check to see if any factors are left
  ## store relevant informatio for later split-inequality encodings
  if (factor.flag) {
    pt.f <- display.tree$mwcpSZ !=0
    if (sum(pt.f) > 0) {
      f.names <- unique(xvar.names[display.tree$parmID[pt.f]])
    }
    else {
      factor.flag <- FALSE
    }
  }
  ##----------------------------------------------------------------
  ##
  ## prepare the tree to be converted into a network
  ##
  ##----------------------------------------------------------------
  ## conversion
  converted.tree <- display.tree
  vars.id <- data.frame(var = c("<leaf>", xvar.names), parmID = 0:length(xvar.names), stringsAsFactors = FALSE)
  converted.tree$var <- vars.id$var[match(display.tree$parmID, vars.id$parmID)]
  ## special symbol to be used for encoding the counter for variables (see note below)
  special <- "999_999"
  # note: we append a counter to the variables, because the data.tree package has trouble when
  # nodes are not unique.
  var.count <- 1:nrow(converted.tree)
  lapply(unique(converted.tree$var), function(vv) {
    pt <- converted.tree$var == vv
    var.count[which(pt)] <<- 1:sum(pt)
  })
  converted.tree$var_count <- var.count
  converted.tree$var_conc <- paste0(converted.tree$var, special, converted.tree$var_count)
  ##----------------------------------------------------------------
  ##
  ## convert the tree to a network data frame, using the fact that the
  ## nativeArray output is a pre-order traversal 
  ##
  ##----------------------------------------------------------------
  ## preliminary
  from_node <- ""
  network <- data.frame()
  num.children <- data.frame(converted.tree, children = 0)
  num.children <- num.children[num.children$var != "<leaf>",, drop = FALSE]
  num.children <- num.children[!duplicated(num.children$var_conc),, drop = FALSE]
  num_children <- as.list(rep(0, nrow(num.children)))
  names(num_children) <- num.children$var_conc
  ## loop (using lapply)
  lapply(1:nrow(converted.tree), function(i) {
    rowi <- converted.tree[i, ]
    xs <- converted.tree$contPT[converted.tree$var_conc == from_node]
    if(i == 1){
      from_node <<- rowi$var_conc
    }
    else{
      ## develop the split encoding
      if(num_children[[from_node]] == 0) {#left split
        side <- "<="
        contPT.pretty <- round(as.numeric(xs), 3)
        split_ineq_pretty <- paste0(side, contPT.pretty)
      }
      else {#right split
        side <- ">"
        split_ineq_pretty <- ""
      }
      ## both numeric and factors are encoded as <= > but factors are secretely in hex notation
      ## !!ADD MACHINE TOLERANCE WHEN SPLIT IS ON ACTUAL X-VALUE !!
      if (is.numeric(xs)) {
        xs <- xs + tolerance
      }
      split_ineq <- paste0(side, xs)
      ## update the network
      to_node <- rowi$var_conc
      new_node <- list(from = from_node, to = to_node, split = split_ineq, split.pretty = split_ineq_pretty)
      network <<- data.frame(rbind(network, new_node, stringsAsFactors = FALSE))
      num_children[[from_node]] <<- num_children[[from_node]] + 1
      if(rowi$var != "<leaf>")
        from_node <<- to_node
      else{
        if(i != nrow(converted.tree)){
          while(num_children[[from_node]] == 2){
            from_node <<- network$from[network$to == from_node]
          }
        }
      }
    }
  })
  ##----------------------------------------------------------------
  ##
  ## process network for factors - clean up the split encoding
  ## encode as complementary pair string
  ##
  ##----------------------------------------------------------------
  if (factor.flag) {
    ## identify which splits need to be cleaned up
    from.names <- network$from
    pt.f <- sort(unique(unlist(lapply(f.names, function(ptrn) {
      grep(ptrn, from.names)
    }))))
    ## identify levels of the factor
    ## otherwise we would have huge sets full of levels that aren't used
    fs <- gsub(paste0(special,".*"),"",from.names[pt.f])
    fs.levels <- sapply(fs, function(fsn) {
      #length(levels(x.org.data[, fsn]))
      xvar.factor$nlevels[which(fsn == xvar.names)]
    })
    ## clean the splits up and encode as complementary pair sets
    split.str <- lapply(1:length(pt.f), function(j) {
      str <- network$split[pt.f[j]]
      ## left split
      if (grepl("<=", str)) {
        str <- sub("<=", "", str)
        str <- strsplit(str, "")[[1]]
        cpr <- 1 + length(str) - which(str != "0")
        cpr <- cpr[cpr <= fs.levels[j]]
        paste0("{", paste(cpr, collapse = ","), "}")
      }
      ## right split
      else {
        str <- sub(">", "", str)
        str <- strsplit(str, "")[[1]]
        cpr <- 1 + length(str) - which(str == "0")
        cpr <- cpr[cpr <= fs.levels[j]]
        paste0("{", paste(cpr, collapse = ","), "}")
      }
    })
    ## pretty complementary pairs
    split.str.pretty <- lapply(1:length(pt.f), function(j) {
      str <- network$split[pt.f[j]]
      ## left split
      if (grepl("<=", str)) {
        str <- sub("<=", "", str)
        str <- strsplit(str, "")[[1]]
        cpr <- 1 + length(str) - which(str != "0")
        cpr <- cpr[cpr <= fs.levels[j]]
        paste0("{", paste(cpr, collapse = ","), "}")
      }
      else {
        ""
      }
    })
    ## replace the previous fake encoding with the now correct "set encoding"
    network$split[pt.f] <- split.str
    network$split.pretty[pt.f] <- split.str.pretty
  }
  ##----------------------------------------------------------------
  ##
  ## create a tree object, see the data.tree package
  ##
  ##----------------------------------------------------------------
  data.tree.network <- data.tree::FromDataFrameNetwork(network, "split")
  ## INTERNAL NODES
  ## label the edges with the splits, and relabel the nodes so that the appended counters are not visible
  data.tree.network$Get(function(node) {
    data.tree::SetEdgeStyle(node, color = "black", label = node$split.pretty, fontcolor = "black")
    data.tree::SetNodeStyle(node, color = "black", fontcolor = "black", penwidth = 3,
                            label = strsplit(node$name, special)[[1]][1])
  })
  ## TERMINAL NODES
  ## loop through the leaves of the tree, each time getting the path down to the leaf, and using this path
  ## to establish a filtering condition to obtain all cases that match the splitting on the way to this leaf.
  ## we filter out cases to correspond to the leaf, and then modify the leaf display to include the number
  ## of cases, followed by the user requested predicted value
  lapply(data.tree.network$leaves, function(node) {
    path_list <- node$path
    var_list <- sapply(path_list, function(x){strsplit(x, special)[[1]][1]})
    var_list[length(var_list)] <- ""
    node_iter <- data.tree.network
    ## make boolean string operatore
    call <- lapply(2:(length(path_list)), function(i) {
      node_iter <<- node_iter[[path_list[[i]]]]
      str <- node_iter$split
      ## numeric boolean operator
      if (!any(grepl("\\{", str))) {
        str <- paste0("x.data$", paste0(var_list[i-1], str))
      }
      ## complementary pair boolean operator
      else {
        str <- gsub("\\{", "", str)
        str <- gsub("\\}", "", str)
        str <- strsplit(str, ",")[[1]]
        str <- paste("==", str, sep = "")
        str <- paste0("(",paste(paste0("x.data$", var_list[i-1], str), collapse = "|"),")")
      }
      str
    })
    call <- paste(call, collapse = " & ")
    ## evaluate the boolean operator
    ## this yields the id's for the cases in the node
    if (!anonymous) {
      pt <- which(eval(parse(text=call)))
      n.cases <- length(pt)
    }
    ## set the edgelabel
    edge.label <- node$split.pretty
    ## set the node label
    ## extract the predicted value in the node
    if (predict.flag && !anonymous) {
      if (!is.factor(yhat)) {
        yhat <- round(mean(yhat[pt], na.rm = TRUE), 2)
        node.label <- paste0("n=", n.cases, "\n", yhat)
      }
      else {
        frqtable <- tapply(yhat[pt], yhat[pt], length)
        pred <- names(frqtable)[which.max(frqtable)]
        node.label <- paste0("n=", n.cases, "\n", pred)
      }
    }
    ## unsupervised family --> no predicted value
    else if (!predict.flag && !anonymous) {
      node.label <- paste0("n=", n.cases)
    }
    ## anonymous
    else {
      node.label <- NULL
    }
    ## set styles
    data.tree::SetGraphStyle(node, rankdir = "TB")
    data.tree::SetEdgeStyle(node, arrowhead = "vee", color = "grey35",
                            penwidth = 3, label = edge.label)
    data.tree::SetNodeStyle(node, style = "filled,rounded", shape = "box",
                 fillcolor = "GreenYellow", penwidth = 3, 
                 fontname = "helvetica", tooltip = data.tree::GetDefaultTooltip,
                 label = node.label)
  })
  invisible(data.tree.network)
}
get.tree <- get.tree.rfsrc
