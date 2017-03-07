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


check.factor <- function(train, test, gfactor) {
  if (!is.null(gfactor)) {
    if (length(gfactor$factor) > 0) {
      test[, match(gfactor$factor, colnames(test))] <- data.frame(
        mclapply(1:length(gfactor$factor),
                 function(k) {
                   fk.train <- train[, colnames(train) == gfactor$factor[k]]
                   fk.test <- test[ , colnames(test) == gfactor$factor[k]]
                   if (length(setdiff(levels(na.omit(fk.test)), levels(fk.train))) > 0) {
                     stop("Factors in test and training data do not match...\n")
                   }
                   factor(as.character(fk.test), levels = gfactor$levels[[k]], exclude = NULL)
                 }))
    }
    if (length(gfactor$order) > 0) {
      test[, match(gfactor$order, colnames(test))] <- data.frame(
        mclapply(1:length(gfactor$order),
                 function(k) {
                   fk.train <- train[, colnames(train) == gfactor$order[k]]
                   fk.test <- test[ , colnames(test) == gfactor$order[k]]
                   if (length(setdiff(levels(na.omit(fk.test)), levels(fk.train))) > 0) {
                     stop("Factors in test and training data do not match...\n")
                   }
                   factor(as.character(fk.test), levels = gfactor$order.levels[[k]], ordered = TRUE)
                 }))
    }
  }
  test
}
is.factor.not.ordered <- function(x) {is.factor(x) && !is.ordered(x)}
extract.factor <- function (dat, generic.names = NULL) {
  generic.types  <- gfactor <- gfactor.order <- gfactor.levels <- gfactor.order.levels <- NULL
  if (is.null(generic.names)) {
    target.names <- names(dat)
  }
  else {
    target.names <- generic.names
  }
  nlevels <- rep(0, length(target.names))
  gfactor <- names(dat)[unlist(lapply(dat, is.factor.not.ordered))]
  gfactor.order <- names(dat)[unlist(lapply(dat, is.ordered))]
  if (!is.null(generic.names)) gfactor <- intersect(gfactor , generic.names)
  if (!is.null(generic.names)) gfactor.order <- intersect(gfactor.order , generic.names)
  if (length(gfactor) > 0) {
    gfactor.levels <- mclapply(cbind(1:(1+length(gfactor))), function(k) {
      if (k <= length(gfactor)) {
        levels(dat[ , names(dat) == gfactor[k]])
      }
        else {
          NULL
        }
    })
    gfactor.levels <- gfactor.levels[-(1+length(gfactor))]
    nlevels[match(gfactor, target.names)] <- unlist(lapply(gfactor.levels, function(g){length(g)}))
  }
  if (length(gfactor.order) > 0 ) {
    gfactor.order.levels <- mclapply(1:(1+length(gfactor.order)), function(k) {
      if (k <= length(gfactor.order)) {
        levels(dat[ , names(dat) == gfactor.order[k]])
      }
        else {
          NULL
        }
    })
    gfactor.order.levels <- gfactor.order.levels[-(1+length(gfactor.order))]
    nlevels[match(gfactor.order, target.names)] <- unlist(lapply(gfactor.order.levels, function(g){length(g)}))
  }
  if (!is.null(generic.names)) {
    generic.types <- rep("R", length(generic.names))
    if (length(gfactor) > 0) {
      generic.types[match(gfactor, generic.names)] <- "C"
    }
    if (length(gfactor.order) > 0) {
      generic.types[match(gfactor.order, generic.names)] <- "I"
    }
  }
    else {
      generic.types <- rep("R", ncol(dat))
      if (length(gfactor) > 0) {
        generic.types[match(gfactor, names(dat))] <- "C"
      }
      if (length(gfactor.order) > 0) {
        generic.types[match(gfactor.order, names(dat))] <- "I"
      }
    }
  if (is.null(gfactor) & is.null(gfactor.order)) {
    return(NULL)
  }
    else {
      return(list(factor=gfactor,
                  order=gfactor.order,
                  levels=gfactor.levels,
                  order.levels=gfactor.order.levels,
                  nlevels = nlevels,
                  generic.types=generic.types))
    }
}
map.factor <- function (gvar, gfactor) {
  if (!is.null(gfactor)) {
    if (length(gfactor$factor) > 0) {
      gvar[, match(gfactor$factor, colnames(gvar))] <- data.frame(
        mclapply(1:length(gfactor$factor),
                 function(k) {            
                   ptk <- (colnames(gvar) == gfactor$factor[k])
                   factor.k <- gfactor$levels[[k]][gvar[ , ptk ]]
                   labels.k <- gfactor$levels[[k]][sort(unique(gvar[ , ptk ]))]
                   gk <- factor(factor.k, labels = labels.k, levels = labels.k, exclude = NULL)  
                   if (length(setdiff(gfactor$levels[[k]], labels.k)) > 0) {
                     gk <- factor(as.character(gk), levels = gfactor$levels[[k]])
                   }
                   gk
                 }))
    }
    if (length(gfactor$order) > 0) {
      gvar[, match(gfactor$order, colnames(gvar))] <- data.frame(
        mclapply(1:length(gfactor$order),
                 function(k) {
                   ptk <- (colnames(gvar) == gfactor$order[k])
                   factor.k <- gfactor$order.levels[[k]][gvar[ , ptk ]]
                   labels.k <- gfactor$order.levels[[k]][sort(unique(gvar[ , ptk ]))]
                   gk <- factor(factor.k, labels = labels.k, levels = labels.k,
                                exclude = NULL, ordered = TRUE)
                   if (length(setdiff(gfactor$order.levels[[k]], labels.k)) > 0) {
                     gk <- factor(as.character(gk), levels = gfactor$order.levels[[k]], ordered = TRUE)
                   }
                   gk
                 }))
    }
  }
  gvar
}
rm.na.levels <- function(dat, xvar.names=NULL) {
  factor.names <- names(dat)[unlist(lapply(dat, is.factor))]
  if (!is.null(xvar.names)) factor.names <- intersect(factor.names , xvar.names)
  if (length(factor.names) > 0) {
    levels.na.pt <- unlist(lapply(1:length(factor.names), function(k) {
      any(levels(dat[ , names(dat) == factor.names[k]]) == "NA")
    }))
    if (any(levels.na.pt)) {
      factor.names <- factor.names[levels.na.pt]
      dat[, match(factor.names, names(dat))] <- data.frame(mclapply(1:length(factor.names),
                                                                    function(k) {
                                                                      x <- dat[ , names(dat) == factor.names[k]]
                                                                      levels(x)[levels(x) == "NA"]  <- NA
                                                                      x
                                                                    }))
    }
  }
  dat
}
