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


generic.impute.rfsrc <- function(formula,
                                 data,
                                 ntree = 250,
                                 nimpute = 1,
                                 bootstrap = c("by.root", "by.node", "none"),
                                 mtry = NULL,
                                 nodesize = NULL,
                                 splitrule = NULL,
                                 nsplit = 1,
                                 na.action = c("na.impute"),
                                 xvar.wt = NULL,
                                 seed = NULL,
                                 do.trace = FALSE)
{
  bootstrap <- match.arg(bootstrap, c("by.root", "by.node", "none"))
  importance <- "none"
  na.action <- match.arg(na.action, c("na.impute"))
  forest <- FALSE
  proximity <- FALSE
  var.used <- FALSE
  split.depth <- FALSE
  impute.only <- TRUE
  membership <- FALSE
  miss.tree <- FALSE
  c.names <- colnames(data)
  r.names <- rownames(data)
  object <- rfsrc(formula = formula,
                  data = data,
                  ntree = ntree,
                  bootstrap = bootstrap,
                  mtry = mtry,
                  nodesize = nodesize,
                  splitrule = splitrule,
                  nsplit = nsplit,
                  nimpute = nimpute,
                  xvar.wt = xvar.wt,
                  seed = seed,
                  do.trace = do.trace,
                  importance = importance,
                  na.action = na.action,
                  forest = forest,
                  proximity = proximity,
                  var.used = var.used,
                  split.depth = split.depth,
                  membership = membership,
                  impute.only = impute.only,
                  miss.tree = miss.tree)
  if (is.null(object)) {
    return(NULL)
  }
  rm(data)
  if (is.data.frame(object)) {
    return(invisible(list(data = object, missing = row.col.deleted(object, r.names, c.names)))) 
  }
  if(is.null(object$yvar.names)) {
    imputed.result <- object$xvar
  }
    else {
      imputed.result <- cbind(object$yvar, object$xvar)
    }
  colnames(imputed.result) <- c(object$yvar.names, object$xvar.names)
  if (nimpute == 1) {
    imputed.result[object$imputed.indv, ] <- object$imputed.data
  }
  rm(object)
  invisible(list(data = imputed.result, missing = row.col.deleted(imputed.result, r.names, c.names)))
}
