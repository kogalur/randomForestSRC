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


print.rfsrc <- function(x, outcome.target = NULL, ...) {
  if (sum(inherits(x, c("rfsrc", "forest"), TRUE) == c(1, 2)) == 2) {
    print.default(x)
    return()
  }
  if(sum(inherits(x, c("rfsrc", "plot.variable"), TRUE) == c(1, 2)) == 2) {
    print.default(x)
    return()
  }
  if (sum(inherits(x, c("rfsrc", "partial"), TRUE) == c(1, 2)) == 2) {
    print.default(x)
    return()
  }
  sf.flag <- FALSE
  if (sum(inherits(x, c("rfsrc", "synthetic"), TRUE) == c(1, 2)) == 2) {
    if (sum(inherits(x, c("rfsrc", "synthetic", "oob"), TRUE) == c(1, 2, 3)) != 3) {
      sf.flag <- TRUE
      sf.message <- "OOB was not used for synthetic forests, error rates/VIMP will be unreliable"
    }
    x <- x$rfSyn
  }
  if (sum(inherits(x, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rfsrc", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("This function only works for objects of class `(rfsrc, grow)' or '(rfsrc, predict)'.")
  }
  if (sum(inherits(x, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
    grow.mode <- TRUE
  }
    else {
      grow.mode <- FALSE
    }
  family.org <- x$family
  yvar.dim <- ncol(x$yvar)
  outcome.target <- get.univariate.target(x, outcome.target)
  x <- coerce.multivariate(x, outcome.target)
  if (grepl("surv", x$family)) {
    event <- get.event.info(x)$event
    n.event <- 1
    if (!is.null(event)) {
      n.event <- length(unique(event))
      if (length(event) > 0) {
        event.freq <- paste(tapply(event, event, length), collapse = ", ")
      }
        else {
          event.freq <- 0
        }
    }
  }
  if (x$family == "class") {
    if (!is.null(x$yvar)) {
      event.freq <- paste(tapply(x$yvar, x$yvar, length), collapse = ", ")
    }
      else {
        event.freq <- NA
      }
    if (!is.null(x$err.rate)) {
      conf.matx <- table(x$yvar, if(!is.null(x$class.oob) && !all(is.na(x$class.oob))) x$class.oob else x$class)
      conf.matx <- cbind(conf.matx,  class.error = round(1 - diag(conf.matx)/rowSums(conf.matx, na.rm = TRUE), 4))
      names(dimnames(conf.matx)) <- c("  observed", "predicted")
      brierS <- brier(x$yvar, if(!is.null(x$predicted.oob) && !all(is.na(x$predicted.oob))) x$predicted.oob else x$predicted)
    }
      else {
        conf.matx <- brierS <- NULL
      }
  }
  if (!is.null(x$err.rate)) {
    err.rate <- cbind(x$err.rate)
    if (grepl("surv", x$family)) {
      err.rate <- paste(round(100 * err.rate[nrow(err.rate), ], 2), "%", collapse=", ", sep = "")
    }
      else if (x$family == "class") {
        overall.err.rate <- paste(round(100 * err.rate[nrow(err.rate), 1], 2), "%", sep = "")
        err.rate <- paste(round(err.rate[nrow(err.rate), ], 2), collapse=", ", sep = "")
        brierS <- round(100 * brierS, 2)
      }
        else if (x$family == "regr") {
          per.var <- round(100 * (1 - err.rate[nrow(err.rate), ] / var(x$yvar, na.rm = TRUE)), 2)
          err.rate <- round(err.rate[nrow(err.rate), ], 2)
        }
          else {
            err.rate <- NULL
          }
  }
    else {
      err.rate <- NULL
    }
  if (is.null(x$nsplit)) {
    x$nsplit <- 0
  }
  if (grow.mode) {
    cat("                         Sample size: ", x$n,                 "\n", sep="")
    if (grepl("surv", x$family)) {
      if (n.event > 1) {
        cat("                    Number of events: ", event.freq,  "\n", sep="")
      }
        else {
          cat("                    Number of deaths: ", x$ndead,   "\n", sep="")
        }
    }
    if (x$family == "class") {
      cat("           Frequency of class labels: ", event.freq,          "\n", sep="")
    }
    if (!is.null(x$imputed.indv)) {
      cat("                    Was data imputed: ", "yes",               "\n", sep="")
                                        #cat("                         Missingness: ",
                                        #    round(100*length(x$imputed.indv)/x$n,2), "%\n", sep="")      
    }
    cat("                     Number of trees: ", x$ntree,                "\n",sep="")
    cat("          Minimum terminal node size: ", x$nodesize,             "\n", sep="")
    cat("       Average no. of terminal nodes: ", mean(x$leaf.count),     "\n", sep="")
    cat("No. of variables tried at each split: ", x$mtry,                 "\n", sep="")
    cat("              Total no. of variables: ", length(x$xvar.names),   "\n", sep="")
    if (!x$univariate) { 
      cat("              Total no. of responses: ", yvar.dim,   "\n", sep="")
      cat("         User has requested response: ", outcome.target,        "\n", sep="")
    }
    cat("                            Analysis: ", family.pretty(family.org),"\n", sep="")
    cat("                              Family: ", family.org,               "\n", sep="")
    if (x$nsplit > 0 & x$splitrule != "random") {
      cat("                      Splitting rule: ", paste(x$splitrule,"*random*"),"\n", sep="")
      cat("       Number of random split points: ", x$nsplit                   ,  "\n", sep="")
    }
      else {
        cat("                      Splitting rule: ", x$splitrule,         "\n", sep="")
      } 
    if (!is.null(err.rate)) {
      if (x$family == "regr") {
        cat("                % variance explained: ", per.var, "\n", sep="")
      }
      if (x$family == "class" && !is.null(brierS)) {
        cat("              Normalized Brier score:", brierS, "\n")
      } 
      cat("                          Error rate: ", err.rate,            "\n\n", sep="")
    }
    if (x$family == "class" && !is.null(conf.matx)) {
      if (!is.null(x$predicted.oob) && any(is.na(x$predicted.oob))) {
        cat("Confusion matrix (cases with missing OOB predicted values have been removed):\n\n")
      }
        else {
          cat("Confusion matrix:\n\n")
        }
      print(conf.matx)
      cat("\n\tOverall error rate:", overall.err.rate, "\n")
    }
  }
    else {
      cat("  Sample size of test (predict) data: ", x$n,  "\n", sep="")
      if (grepl(x$family, "surv") && !is.null(event)) {
        if (n.event > 1) {
          cat("       Number of events in test data: ", event.freq,  "\n", sep="")
        }
          else {
            cat("       Number of deaths in test data: ", unlist(event.freq),   "\n", sep="")
          }
      }
      if (!is.null(x$imputed.data)) {
        cat("               Was test data imputed: ", "yes",               "\n", sep="")
                                        #cat("                         Missingness: ",
                                        #    round(100*length(x$imputed.indv)/x$n,2), "%\n", sep="")      
      }
      cat("                Number of grow trees: ", x$ntree,             "\n",sep="")
      cat("  Average no. of grow terminal nodes: ", mean(x$leaf.count),  "\n", sep="")
      cat("         Total no. of grow variables: ", length(x$xvar.names), "\n", sep="")  
      if (!x$univariate) { 
        cat("         Total no. of grow responses: ", yvar.dim,   "\n", sep="")
        cat("         User has requested response: ", outcome.target,        "\n", sep="")
      }
      cat("                            Analysis: ", family.pretty(family.org),"\n", sep="")
      cat("                              Family: ", family.org,               "\n", sep="")
      if (!is.null(err.rate)) {
        if (x$family == "regr") {
          cat("                % variance explained: ", per.var, "\n", sep="")
        }
        if (x$family == "class" && !is.null(brierS)) {
          cat("       Test set Normalized Brier score:", brierS, "\n")
        }
        cat("                 Test set error rate: ", err.rate, "\n\n", sep="")
      }
      if (x$family == "class" && !is.null(conf.matx)) {
        if (!is.null(x$predicted.oob) && any(is.na(x$predicted.oob))) {
          cat("Confusion matrix (cases with missing OOB predicted values have been removed):\n\n")
        }
          else {
            cat("Confusion matrix:\n\n")
          }
        print(conf.matx)
        cat("\n\tOverall error rate:", overall.err.rate, "\n")
      }
    }
  if (sf.flag) {
    message(sf.message)
  }
}
