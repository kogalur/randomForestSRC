sidClustering.rfsrc <- function(data,
                                method = "sid",
                                k = NULL,
                                reduce = TRUE,
                                ntree.reduce = function(p, vtry){100 * p / vtry},
                                fast = FALSE,
                                x.no.sid = NULL,
                                use.sid.for.x = TRUE,
                                x.only = NULL, y.only = NULL,
                                dist.sharpen = TRUE,
                                ...)
{
  ##------------------------------------------------------------------
  ##
  ## preliminary processing
  ##
  ##------------------------------------------------------------------
  ## make sure we have a data frame
  if (!is.data.frame(data)) {
    data <- data.frame(data)
  }
  ## decide which method to use
  method <- match.arg(method, c("sid",
                                "unsupv", 
                                "shi-horvath",  "sh",  "SH",  "sh1", "SH1", "sh-1", "SH-1",
                                "shi-horvath2", "sh2", "SH2", "sh-2", "SH-2"))
  if (grepl("sid", method)) {
    methodClass <- 1
  }
  else if (grepl("sh", method)) {
    methodClass <- 2
  }
  else if (grepl("un", method)) {
    methodClass <- 3
  }
  else {
    stop("method not set correctly:", method)
  }
  ## set k if user has not provided a value
  if (is.null(k)) {
    k <- 1:max(1, min(10, nrow(data) - 1))
  }
  else {
    k <- k[k<=nrow(data)]
    k <- k[k != 0]
    if (length(k) == 0) {
      k <- 1:max(1, min(10, nrow(data) - 1))
    }
  }
  ## determine the grow interface - rfsrc or rfsrc.fast?
  if (!fast) {
    rfsrc.grow <- "rfsrc"
  }
  else {
    rfsrc.grow <- "rfsrc.fast"
  }
  ##------------------------------------------------------------------
  ##
  ## run SH-1 to acquire cheap but effective VIMP for dimension reduction
  ##
  ## here's the rationale for Zcut and eps ...
  ## Let I=importance and s=standard error
  ## The left edge of the CI for a strong variable satisfies:
  ##
  ## I - s Zcut > eps' > 0  (eps' quantifies how far left CI is away from 0)
  ##
  ## With some rearrangment:
  ##
  ## I / (eps'/Zcut + s) > Zcut
  ##
  ## Let eps=eps'/Zcut, we have
  ##
  ## I / (eps + s) > Zcut
  ##
  ## As illustration, eps = .25, Zcut = 2 ---> eps'=.5
  ##
  ##------------------------------------------------------------------
  if (methodClass == 1 && reduce && nrow(data) > 2) {
    ## pull unnamed parameters to match against rfsrc formal arguments
    sv.f <- as.formula("classes ~ .")
    dots <- list(...)
    dots$importance <- TRUE
    dots$splitrule <- NULL ##needed to allow Mahalanobis splitting for sidClustering
    rfnames <- names(formals(rfsrc))
    rfnames <- rfnames[rfnames != "distance" & rfnames != "proximity"]
    ## remove any column with less than two unique values
    void.var <- sapply(data, function(x){length(unique(x, na.rm = TRUE)) < 2})
    if (sum(void.var) > 0) {
      data[, which(void.var)] <- NULL
    }
    ## ... the remaining part of the code block works if there is only 1 variable left
    ## ... there needs to be enough unique values always 
    ## run the supervised classifier and acquire VIMP
    dots.sv <- dots
    dots.sv$mtry <- dots.sv$ntree <- NULL
    o <- do.call(holdout.vimp, c(list(formula = sv.f, data = make.sh(data, 1), ntree = ntree.reduce),
                     dots.sv[names(dots.sv) %in% rfnames]))
    ## identify significant variables
    Z <- o$importance[, 1]
    ## we need at least TWO (2) variables to pass filtering (otherwise there will be no SID)
    if (sum(Z > 0) > 1) {
      xvar.filter <- names(Z)[Z > 0]
    }
    else {
      ## there were less than two variables
      ## careful - it's possible there's only one variable
      xvar.filter <- names(Z)[order(Z, decreasing = T)[1:min(length(Z), 2)]]
    }
    ## subset the data
    data <- data[, xvar.filter, drop = FALSE]
  }
  ##-----------------------------------------------------------
  ##
  ## sid clustering method
  ##
  ##-----------------------------------------------------------
  if (methodClass == 1) {
    ## terminate if there is one case or less than two variables
    terminate <- FALSE
    ## did the user provide protected x-values - convert to data frame if so
    if (!is.null(x.no.sid) && !is.character(x.no.sid)) {
      x.no.sid <- data.frame(x.no.sid)
    }
    ## did the user provide protected x-values as a character vector?
    if (!is.null(x.no.sid) && is.character(x.no.sid)) {
      drops <- (names(data) %in% x.no.sid)
      if (sum(drops) > 0) {##check that x.no.sid is still viable - reduce (filtering) could remove it
        x.no.sid <- data[, drops, drop = FALSE]
        if (sum(!drops) > 0) {
          data <- data[ , !drops, drop = FALSE]
        }
        else {
          terminate <- TRUE
        }
      }
      else {
        x.no.sid <- NULL
      }
    }
    ## last chance to fail 
    if (terminate || ncol(data) < 2 || nrow(data) < 2) {
      rO <- list(clustering = rep(1, nrow(data)), rf = NULL, dist = NULL, sid = NULL, outcome.names = NULL)
      class(rO) <- c("rfsrc-failed", "sidClustering", method)
      return(invisible(rO))
    }
    ## sidIFY the data
    sid.o <- make.sid(data)
    rf.dat <- data.frame(sid.o$y, sid.o$x)
    if (use.sid.for.x) {
      outcome.names <- colnames(sid.o$y)
    }
    else {
      outcome.names <- colnames(sid.o$x)
    }
    ## if there are protected x-values, add them now to the data
    if (!is.null(x.no.sid) && is.data.frame(x.no.sid) && sum(names(data) %in% names(x.no.sid)) == 0) {
      rf.dat <- data.frame(rf.dat, x.no.sid)
      outcome.names <- c(outcome.names, colnames(x.no.sid))
    }
    ##remove the original data
    rm(data)
    ## make the mv formula
    mv.f <- as.formula(paste("Multivar(", paste(outcome.names, collapse = ","), paste(") ~ ."), sep = ""))
    ## pull unnamed parameters to match against rfsrc formal arguments
    ### extract the distance (oob by default)
    dots <- list(...)
    if (is.null(dots$distance)) {
      dots$distance <- "oob"
    }
    rfnames <- names(formals(rfsrc))
    ## mvRF call - in extreme cases (n=1,2 sample size) this can fail
    rf <- tryCatch({do.call(rfsrc.grow,
       c(list(formula = mv.f, data = rf.dat), dots[names(dots) %in% rfnames]))}, error=function(ex){NULL})
    if (is.null(rf)) {## call failed
      rO <- list(clustering = rep(1, nrow(rf.dat)), rf = NULL, dist = NULL, sid = NULL, outcome.names = NULL)
      class(rO) <- c("rfsrc-failed", "sidClustering", method)
      return(invisible(rO))
    }
    d <- rf$distance
    ## euclidean distance sharpening (depends on distance.sharpening option - this is slow!!!!)
    if (dist.sharpen) {
      d <- distance(d)
    }
    ## OOB distance may contain NA's - set these to 1
    d[is.na(d)] <- 1
  }##end RF-SID
  ##-----------------------------------------------------------
  ##
  ## Shi-Horvath (mode 1 and 2) clustering
  ## converts problem to classification using artifical class data
  ##
  ##-----------------------------------------------------------
  if (methodClass == 2) {
    ## set NULL values
    sid.o <- outcome.names <- NULL
    ## make SH synthetic class data - mode 1 or 2
    if (grepl("2", method)) {
      mode <- 2
      method <- "sh2"
    }
    else {
      mode <- 1
      method <- "sh1"
    }
    hv.dat <- make.sh(data, mode)
    ## make the supervised formula
    sv.f <- as.formula("classes ~ .")
    ## pull unnamed parameters to match against rfsrc formal arguments
    ## extract the distance (oob by default)
    dots <- list(...)
    if (is.null(dots$proximity)) {
      dots$proximity <- "oob"
    }
    rfnames <- names(formals(rfsrc))
    ## run the supervised classifier and extract the proximity distance
    rf <- do.call(rfsrc.grow, c(list(formula = sv.f, data = hv.dat), dots[names(dots) %in% rfnames]))
    d <- 1 - rf$proximity[rf$yvar == 1, rf$yvar == 1]
    ## OOB proximity may contain NA's - set these to 1
    d[is.na(d)] <- 1
  }
  ##-----------------------------------------------------------
  ##
  ## plain vanilla unsupervised splitting
  ## multivariate splitting x against x
  ##
  ##-----------------------------------------------------------
  if (methodClass == 3) {
    ## set NULL values
    sid.o <- outcome.names <- NULL
    ## scale the y-values
    y.dat <- scale(data.matrix(data))
    ## make mv formula 
    ## duplicate the data
    ## watch out for x.only and y.only requests
    outcome.values <- colnames(y.dat[, setdiff(colnames(y.dat), x.only), drop = FALSE])
    mv.f <- get.mv.formula(outcome.values)
    rf.data <- data.frame(y.dat, data[, setdiff(colnames(data), y.only), drop = FALSE])
    ## remove data
    rm(data)
    ## pull unnamed parameters to match against rfsrc formal arguments
    ### extract the distance (oob by default)
    dots <- list(...)
    if (is.null(dots$distance)) {
      dots$distance <- "oob"
    }
    rfnames <- names(formals(rfsrc))
    ## mvRF call - in extreme cases (n=1,2 sample size) this can fail
    rf <- tryCatch({do.call(rfsrc.grow,
                            c(list(formula = mv.f, data = rf.data), 
                              dots[names(dots) %in% rfnames]))}, error=function(ex){NULL})
    if (is.null(rf)) {## call failed
      rO <- list(clustering = rep(1, nrow(rf.data)), rf = NULL, dist = NULL)
      #class(rO) <- c("rfsrc-failed", "sgreedy.unsupv")
      return(invisible(rO))
    }
    d <- rf$distance
    ## euclidean distance sharpening (depends on distance.sharpening option - this is slow!!!!)
    if (dist.sharpen) {
      d <- distance(d)
    }
    ## OOB distance may contain NA's - set these to 1
    d[is.na(d)] <- 1
  }    
  ##-----------------------------------------------------------
  ##
  ## hierarchical clustering
  ##
  ##-----------------------------------------------------------
  clus <- cutree(hclust(as.dist(d)), k = k)
  ##-----------------------------------------------------------
  ##
  ## return the goodies
  ##
  ##-----------------------------------------------------------
  rO <- list(clustering = clus, rf = rf, dist = d, sid = sid.o, outcome.names = outcome.names)
  class(rO) <- c("rfsrc", "sidClustering", method)
  invisible(rO)
}
sidClustering <- sidClustering.rfsrc
