find.interaction.rfsrc <- function(
  object, 
  xvar.names,
  cause,
  m.target = NULL,
  importance = c("permute", "random", "anti", "permute.ensemble", "random.ensemble", "anti.ensemble"),
  method = c("maxsubtree", "vimp"),
  sorted = TRUE,
  nvar = NULL, 
  nrep  = 1,
  na.action = c("na.omit", "na.impute", "na.random"),
  seed = NULL,
  do.trace = FALSE,
  verbose = TRUE,
  ...)
{
  ## incoming object must be a grow forest or a forest object
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2)
    stop("this function only works for objects of class `(rfsrc, grow)' or '(rfsrc, forest)'")
  ## grow forests must have true forest information
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.forest.missing(object)) {
      stop("Forest information for prediction is missing.  Re-run rfsrc (grow call) with forest=TRUE")
    }
  }
  ## specify the default event type for CR
  if (missing(cause)) {
    cause <- 1
  }
  ## unsupervised family: minimal depth is the only permissible method
  if (object$family == "unsupv") {
    method <- "maxsubtree"
  }
  ## Verify key options
  method <- match.arg(method,  c("maxsubtree", "vimp"))
  importance <- match.arg(importance, c("permute", "random", "anti",
                                        "permute.ensemble", "random.ensemble", "anti.ensemble"))
  ## get the event data
  event.info <- get.event.info(object)
  n.event <- max(1, length(event.info$event.type))
  ## acquire the target outcome (if there is one)
  m.target <- get.univariate.target(object, m.target)
  ## extract the importance
  if (n.event > 1) {
    object.imp <- NULL
    interact.imp.list.names <- paste("event.", 1:length(event.info$event.type), sep = "")
  }
    else {
      object.imp <- cbind(coerce.multivariate(object, m.target)$importance)[, cause, drop = FALSE]
    }
  ## pull the original xvariable names
  xvar.org.names <- object$xvar.names
  ## has the user provided a subset of variables to focus on?
  if (!missing(xvar.names)) {
    if (sum(is.element(xvar.org.names, xvar.names)) == 0) {
      stop("Variables do not match original analysis:", xvar.names)
    }
    xvar.names <- unique(xvar.names[is.element(xvar.names, xvar.org.names)])
    if (length(xvar.names) == 1) 
      stop("Pairwise comparisons require more than one candidate variable.")
  }
    else {
      xvar.names <- xvar.org.names
    }
  ## determine the number of remaining variables
  n.interact <- length(xvar.names)
  ## sort the variables by VIMP?
  if (sorted) {
    if (!is.null(object.imp)) {
      o.r <- order(object.imp[xvar.names, ], decreasing = TRUE)
      xvar.names <- xvar.names[o.r]
    }
  }
  ## restrict attention to top nvar variables?
  if (!missing(nvar)) {
    n.interact <- min(n.interact, max(round(nvar), 1))
    xvar.names <- xvar.names[1:n.interact]
  }
  if (n.interact == 1) {
    stop("Pairwise comparisons require more than one candidate variable.")
  }
  ## VIMP approach
  if (method == "vimp") {
    yvar.dim <- ncol(object$yvar)
    ## Save the original family.
    family.org <- object$family
    if (n.event > 1) {
      interact.imp.list <- vector("list", n.event)
      names(interact.imp.list) <- interact.imp.list.names
    }
    for (j in 1:n.event) {
      rownames.interact.imp <- interact.imp <- NULL
      target.dim <- ifelse(n.event > 1, j, 1)
      if (verbose && n.event > 1) {
        cat("--> event", j, "\n")
      }
      for (k in 1:(n.interact-1)) {
        n.joint.var <- n.interact - k
        imp <- rep(0 , 1 + n.joint.var)
        imp.joint <- rep(0, n.joint.var)
        for (l in (k+1):n.interact) {
          if (verbose) {
            cat("Pairing",xvar.names[k],"with",xvar.names[l],"\n")
          }
          for (m in 1:nrep) {
            imp.indv.m <- c(cbind(coerce.multivariate(vimp(object, xvar.names[c(k,l)],
                                                           m.target = m.target,
                                                           importance = importance, joint = FALSE, 
                                                           na.action = na.action, seed = seed, 
                                                           do.trace = do.trace), m.target)$importance)[, target.dim])
            imp.joint.m <- coerce.multivariate(vimp(object, xvar.names[c(k,l)], m.target = m.target,
                                                    importance = importance, joint = TRUE,
                                                    na.action = na.action, seed = seed,
                                                    do.trace = do.trace), m.target)$importance[target.dim]
            imp[1] <- imp[1] + imp.indv.m[1]
            imp[l-k+1] <- imp[l-k+1] + imp.indv.m[2]
            imp.joint[l-k] <- imp.joint[l-k] + imp.joint.m
          }
        }
        imp[1] <- imp[1] / n.joint.var
        imp <- imp/nrep
        imp.joint <- imp.joint/nrep
        interact.imp <- rbind(interact.imp,
                              cbind(imp[1], imp[-1], imp.joint, (imp[1] + imp)[-1], imp.joint - (imp[1] + imp)[-1]))
        rownames.interact.imp <- c(rownames.interact.imp,
                                   paste(xvar.names[k],":",xvar.names[(k+1):n.interact],
                                         sep=""))
      }
      colnames(interact.imp) <- c("Var 1", "Var 2","Paired","Additive","Difference")
      rownames(interact.imp) <- rownames.interact.imp
      if (n.event > 1) {
        interact.imp.list[[j]] <- interact.imp
      }
    }
    ## CR details
    if (n.event > 1) {
      interact.imp <- interact.imp.list
    }
    ## output table
    if (verbose) {
      cat("\n")
      cat("                              Method: ", method,                       "\n", sep="")
      cat("                    No. of variables: ", n.interact,                   "\n", sep="")
      if (family.org == "regr+" | family.org == "class+" | family.org == "mix+") {
        cat("              Total no. of responses: ", yvar.dim,                 "\n", sep="")
        cat("         User has requested response: ", m.target,           "\n", sep="")
      }
      cat("           Variables sorted by VIMP?: ", sorted,                       "\n", sep="")
      cat("   No. of variables used for pairing: ", n.interact,                   "\n", sep="")
      cat("    Total no. of paired interactions: ", length(rownames.interact.imp),"\n", sep="")
      cat("            Monte Carlo replications: ", nrep,                         "\n", sep="")
      cat("    Type of noising up used for VIMP: ", importance,                   "\n", sep="")
      cat("\n")
      if (n.event == 1) print(round(interact.imp, 4)) else print(interact.imp)
    }
    ## return the goodies
    invisible(interact.imp)
  }
    else {
      ## maximal subtree approach
      max.obj <- max.subtree(object, sub.order = TRUE, max.order = 1)
      sub.order <- max.obj$sub.order
      if (sorted) {
        o.r <- order(diag(sub.order), decreasing = FALSE)
        sub.order <- sub.order[o.r, o.r]
      }
      xvar.pt <- is.element(colnames(sub.order), xvar.names[1:n.interact])
      sub.order <- sub.order[xvar.pt, xvar.pt]
      ## output table
      if (verbose) {
        cat("\n")
        cat("                              Method: ", method,              "\n", sep="")
        cat("                    No. of variables: ", n.interact,          "\n", sep="")
        cat("  Variables sorted by minimal depth?: ", sorted,              "\n", sep="")
        cat("\n")
        print(round(sub.order, 2))
      }
      ## return the goodies
      invisible(sub.order)
    }
}
find.interaction <- find.interaction.rfsrc
