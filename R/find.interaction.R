##  **********************************************************************
##  **********************************************************************
##  
##    RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
##  
##    This program is free software; you can redistribute it and/or
##    modify it under the terms of the GNU General Public License
##    as published by the Free Software Foundation; either version 3
##    of the License, or (at your option) any later version.
##  
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##  
##    You should have received a copy of the GNU General Public
##    License along with this program; if not, write to the Free
##    Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
##    Boston, MA  02110-1301, USA.
##  
##    ----------------------------------------------------------------
##    Project Partially Funded By: 
##    ----------------------------------------------------------------
##    Dr. Ishwaran's work was funded in part by DMS grant 1148991 from the
##    National Science Foundation and grant R01 CA163739 from the National
##    Cancer Institute.
##  
##    Dr. Kogalur's work was funded in part by grant R01 CA163739 from the 
##    National Cancer Institute.
##    ----------------------------------------------------------------
##    Written by:
##    ----------------------------------------------------------------
##      Hemant Ishwaran, Ph.D.
##      Director of Statistical Methodology
##      Professor, Division of Biostatistics
##      Clinical Research Building, Room 1058
##      1120 NW 14th Street
##      University of Miami, Miami FL 33136
##  
##      email:  hemant.ishwaran@gmail.com
##      URL:    http://web.ccs.miami.edu/~hishwaran
##      --------------------------------------------------------------
##      Udaya B. Kogalur, Ph.D.
##      Adjunct Staff
##      Department of Quantitative Health Sciences
##      Cleveland Clinic Foundation
##      
##      Kogalur & Company, Inc.
##      5425 Nestleway Drive, Suite L1
##      Clemmons, NC 27012
##  
##      email:  ubk@kogalur.com
##      URL:    https://github.com/kogalur/randomForestSRC
##      --------------------------------------------------------------
##  
##  **********************************************************************
##  **********************************************************************


find.interaction.rfsrc <- function(
  object, 
  xvar.names,
  cause,
  outcome.target = NULL,
  importance = c("permute", "random", "anti", "permute.ensemble", "random.ensemble", "anti.ensemble"),
  method = c("maxsubtree", "vimp"),
  sorted = TRUE,
  nvar = NULL, 
  nrep  = 1,
  subset,
  na.action = c("na.omit", "na.impute"),
  seed = NULL,
  do.trace = FALSE,
  verbose = TRUE,
  ...)
{
  if (is.null(object)) stop("Object is empty!")
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rfsrc, grow)' or '(rfsrc, forest)'.")
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.null(object$forest)) 
      stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'.")
  }
  if (missing(cause)) {
    cause <- 1
  }
  if (object$family == "unsupv") {
    method <- "maxsubtree"
  }
  method <- match.arg(method,  c("maxsubtree", "vimp"))
  importance <- match.arg(importance, c("permute", "random", "anti",
                                        "permute.ensemble", "random.ensemble", "anti.ensemble"))
  event.info <- get.event.info(object)
  n.event <- max(1, length(event.info$event.type))
  outcome.target <- get.univariate.target(object, outcome.target)
  if (n.event > 1) {
    object.imp <- NULL
    interact.imp.list.names <- paste("event.", 1:length(event.info$event.type), sep = "")
  }
    else {
      object.imp <- cbind(coerce.multivariate(object, outcome.target)$importance)[, cause, drop = FALSE]
    }
  xvar.org.names <- object$xvar.names
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
  n.interact <- length(xvar.names)
  if (sorted) {
    if (!is.null(object.imp)) {
      o.r <- order(object.imp[xvar.names, ], decreasing = TRUE)
      xvar.names <- xvar.names[o.r]
    }
  }
  if (!missing(nvar)) {
    n.interact <- min(n.interact, max(round(nvar), 1))
    xvar.names <- xvar.names[1:n.interact]
  }
  if (n.interact == 1) {
    stop("Pairwise comparisons require more than one candidate variable.")
  }
  if (method == "vimp") {
    yvar.dim <- ncol(object$yvar)
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
                                                           outcome.target = outcome.target,
                                                           importance = importance, joint = FALSE, 
                                                           na.action = na.action, subset = subset, seed = seed, 
                                                           do.trace = do.trace), outcome.target)$importance)[, target.dim])
            imp.joint.m <- coerce.multivariate(vimp(object, xvar.names[c(k,l)], outcome.target = outcome.target,
                                                    importance = importance, joint = TRUE,
                                                    na.action = na.action, subset = subset, seed = seed,
                                                    do.trace = do.trace), outcome.target)$importance[target.dim]
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
    if (n.event > 1) {
      interact.imp <- interact.imp.list
    }
    if (verbose) {
      cat("\n")
      cat("                              Method: ", method,                       "\n", sep="")
      cat("                    No. of variables: ", n.interact,                   "\n", sep="")
      if (family.org == "regr+" | family.org == "class+" | family.org == "mix+") {
        cat("              Total no. of responses: ", yvar.dim,                 "\n", sep="")
        cat("         User has requested response: ", outcome.target,           "\n", sep="")
      }
      cat("           Variables sorted by VIMP?: ", sorted,                       "\n", sep="")
      cat("   No. of variables used for pairing: ", n.interact,                   "\n", sep="")
      cat("    Total no. of paired interactions: ", length(rownames.interact.imp),"\n", sep="")
      cat("            Monte Carlo replications: ", nrep,                         "\n", sep="")
      cat("    Type of noising up used for VIMP: ", importance,                   "\n", sep="")
      cat("\n")
      if (n.event == 1) print(round(interact.imp, 4)) else print(interact.imp)
    }
    invisible(interact.imp)
  }
    else {
      max.obj <- max.subtree(object, sub.order = TRUE, max.order = 1)
      sub.order <- max.obj$sub.order
      if (sorted) {
        o.r <- order(diag(sub.order), decreasing = FALSE)
        sub.order <- sub.order[o.r, o.r]
      }
      xvar.pt <- is.element(colnames(sub.order), xvar.names[1:n.interact])
      sub.order <- sub.order[xvar.pt, xvar.pt]
      if (verbose) {
        cat("\n")
        cat("                              Method: ", method,              "\n", sep="")
        cat("                    No. of variables: ", n.interact,          "\n", sep="")
        cat("  Variables sorted by minimal depth?: ", sorted,              "\n", sep="")
        cat("\n")
        print(round(sub.order, 2))
      }
      invisible(sub.order)
    }
}
find.interaction <- find.interaction.rfsrc
