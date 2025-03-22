##  make SH data (modes 1 and 2)
make.sh <- function(dat, mode = 1) {
  ## extract sample size dimension
  nr <- dim(dat)[[1]]
  nc <- dim(dat)[[2]]
  if (nc == 0) {
    stop("can't make SH data ... not enough unique values\n")
  }
  ## coerce to data frame format
  if (!is.data.frame(dat)) {
    dat <- data.frame(dat)
  }
  ## mode 1
  if (mode == 1) {
    data.frame(classes = factor(c(rep(1, nr), rep(2, nr))),
      rbind(dat, data.frame(mclapply(dat, sample, replace = TRUE))))
  }
  ## mode 2
  else {
    data.frame(classes = factor(c(rep(1, nr), rep(2, nr))),
      rbind(dat, data.frame(mclapply(dat, function(x) {
        if (is.factor(x)) {
          resample(x, replace = TRUE)
        }
        else {
          runif(nr, min(x, na.rm = TRUE), max(x, na.rm = TRUE))
        }
      }))))
  }
}
##  make sid (staggered interaction data)
make.sid <- function(dat, order.by.range = TRUE, delta = NULL) {
  ## coerce to data frame format
  if (!is.data.frame(dat)) {
    dat <- data.frame(dat)
  }
  ## remove any column with less than two unique values
  void.var <- sapply(dat, function(x) {length(unique(x, na.rm = TRUE)) < 2})
  if (sum(void.var) > 0) {
    dat[, which(void.var)] <- NULL
  }
  ## there might be nothing left (small sample size issue)
  if (ncol(dat) == 0) {
    stop("can't make sid data ... not enough unique values\n")
  }
  ## order columns by range of values: noted improvement Alex 04/30/2018
  ## we do this before positivity and translation - but we could just as
  ## well have done this afterwards
  if (order.by.range) {
    order.range <- order(sapply(dat, function(x) {
    if (is.factor(x)) {
      1##changed from -Inf to 1 which is the correct value of a factor noted by Alex 06/05/2018
    }
    else {
      diff(range(x, na.rm = TRUE))
    }
  }), decreasing = TRUE)
  dat <- dat[, order.range, drop = FALSE]
  }
  ## make continuous variables positive
  posdat <- dat
  lapply(1:ncol(posdat), function(i) {
    if(!is.factor(posdat[, i]) && min(posdat[, i], na.rm = TRUE) < 0) {
      posdat[, i] <<- abs(min(posdat[, i], na.rm = TRUE)) + posdat[, i]
    }
    NULL
  })
  ## define the delta offset value (delta = 1 by default)
  ## if ordering is in effect, translate features to have the same maximum value
  ## this was observed in counter-example of theorem of revised paper 05/04/2018 
  if (is.null(delta)) {
    delta <- 1
  }
  if (order.by.range) {
    cont.feature <- sapply(posdat, function(x) {!is.factor(x)})
    ## acquire the maximum value
    if (sum(cont.feature) > 1) {
      maxV <- max(c(0, sapply(which(cont.feature), function(i) {
        max(posdat[, i], na.rm = TRUE)
      })), na.rm = TRUE)
      ## translate the features to have the same max
      lapply(which(cont.feature), function(i) {
        posdat[, i] <<- (maxV - max(posdat[, i], na.rm = TRUE)) + posdat[, i]
        NULL
      })	
    }
  }		
  ## two level factors are converted to numeric
  binary.fac <- sapply(posdat, function(x) {is.factor(x) & length(levels(x)) == 2})
  if (sum(binary.fac) > 0) {
    lapply(which(binary.fac), function(i) {
      names(posdat)[i] <<- names(posdat)[i]
      posdat[, i] <<- as.numeric(posdat[, i])
      NULL 
    })
  }
  ## the following list makes it possible to map SID back to original variable names 
  org.names <- list()
  ## make anova data for remaining factors - i.e. make everything numeric
  ## updated to make column names more appealing for binary factors (03/14/2018)
  counter <- 0
  numdat <- data.frame(lapply(1:ncol(posdat), function(j) {
    m.j <- model.matrix(~.-1, posdat[,j, drop = FALSE])
    ##SID main effect names mapped back to original names
    org.names[(counter + 1):(counter + ncol(m.j))] <<- colnames(posdat)[j]
    counter <<- counter + ncol(m.j)
    m.j
  }))
  ## stagger the positive data using delta
  ## use integer values when staggering to minimize creating double precision values
  ## the latter is accomplished using ceiling
  numvar <- dim(numdat)[2]
  staggerdat <- numdat
  staggerdat[, 1] <- staggerdat[, 1] + delta
  if (numvar > 1) {##handles pathological case of one column design matrix - maybe add this as a failure check?
    lapply(2:numvar, function(i) {
      staggerdat[, i] <<- staggerdat[, i] + ceiling(max(staggerdat[, (i-1)], na.rm = TRUE)) + delta
    })
  }
  ## make interactions
  ## -- do not create interactions WITHIN levels of the same factor -- noted by Alex M. 04/25/17
  ##    such instances are caught using unstaggered numerical data and checking x[,i]*x[,j]=constant
  ## -- interactions BETWEEN distinct factors (two or more levels) must always produce a factor
  ##    we catch this by checking the number of distinct values of each variable
  ##    all factors at this point have <= 2 unique values (whether they are dummy anova or numerical)
  intdat <- staggerdat
  if (numvar > 1) {##interactions do not apply if only one variable is supplied - maybe add this as a failure check?
    ##SID interaction names mapped back to original names
    ##we do this first inside of an lapply and not in next mclapply
    counter <- numvar
    lapply(1:(numvar-1), function(i) {
      for (j in (i + 1):numvar) {
        if (length(unique(numdat[,i] * numdat[,j])) > 1) {
          counter <<- counter + 1
          org.names[[counter]] <<- c(org.names[[i]], org.names[[j]])
        }
      }
      NULL
    })
    ## suggested by Alex M. 02/19/2019
    ## make interaction matrix after creating interactions - faster on big p
    ints <- mclapply(1:(numvar - 1), function(i) {
      d <- data.frame(rep(NA, dim(dat)[1]))
      d <- d[, -1]
      counter.ints <- 1
      for (j in (i + 1):numvar) {
        if (length(unique(numdat[,i] * numdat[,j])) > 1) {
          if (length(unique(numdat[, i])) <= 2 & length(unique(numdat[, j])) <= 2) {
            holder <- factor(staggerdat[, i] * staggerdat[, j])
            ## suggested by Alex M. 11/12/17
            levels(holder) <- c('FF','FT','TF','TT') #new factor levels, should always come out in same order
            d <- cbind(d, holder)
          }
          else {
            d <- cbind(d, staggerdat[, i] * staggerdat[, j])
          }
          names(d)[counter.ints] <- paste(names(intdat)[i], "_", names(intdat)[j], sep = "")
          counter.ints <- counter.ints + 1
        }
      }
      d
    })
    intdat <- data.frame(intdat, ints)
  }
  ## suggested by Alex M. 11/12/17
  ## final processing of numeric 0/1 variables to convert them to a binary factor with new levels
  binary.fac2 <- sapply(intdat, function(x){is.numeric(x) & length(unique(x)) == 2})
  if (sum(binary.fac2) > 0) {
    lapply(which(binary.fac2), function(i) {
      intdat[, i] <<- as.factor(intdat[, i])
      levels(intdat[, i]) <<- c('F','T') #new levels, similarly tested
      NULL
    })
  }
  ## pull the "x" and "y" features
  ## y=staggered main effects
  ## x=staggered interactions
  y <- intdat[, 1:numvar, drop = FALSE]
  names(org.names) <- colnames(intdat)
  y.names <- org.names[1:numvar]
  if (numvar > 1) {
    x <- intdat[, -(1:numvar), drop = FALSE]
    x.names <- org.names[-(1:numvar)]
  }
  else {
    x <- x.names <- NULL
  }
  ## return the goodies
  list(y = y, x = x, y.names = y.names, x.names = x.names, delta = delta)
}
#weighted gini/entropy performance metric
#mode is equal to either 'gini' or 'entropy'
sid.perf.metric <- function(truth,cluster,mode=c("entropy", "gini")){
  ## verify mode option
  mode <- match.arg(mode, c("entropy", "gini"))
  ## confusion matrix
  k=length(unique(truth))
  tab=table(truth,cluster)
  clustersizes=colSums(tab)
  clustersizesnorm=clustersizes/sum(clustersizes)
  tabprop=tab
  lapply(1:(dim(tab)[2]),function(i){
    tabprop[,i]<<-tabprop[,i]/clustersizes[i]
    NULL
  })
  if(mode=="entropy"){
    measure=0
    maxmeasure=0
    lapply(1:(dim(tabprop)[2]),function(i){
      clustermeasure=0
      maxclustermeasure=0
      lapply(1:(dim(tabprop)[1]),function(j){
        if(tabprop[j,i]!=0){
          clustermeasure<<-clustermeasure+-tabprop[j,i]*log2(tabprop[j,i])
        }
        maxclustermeasure<<-maxclustermeasure+-1/k*log2(1/k)
        NULL
      })
      measure<<-measure+clustersizesnorm[i]*clustermeasure
      maxmeasure<<-maxmeasure+clustersizesnorm[i]*maxclustermeasure
    })
  }
  if(mode=="gini"){
    measure=0
    maxmeasure=0
    lapply(1:(dim(tabprop)[2]),function(i){
      clustermeasure=1
      maxclustermeasure=1
      lapply(1:(dim(tabprop)[1]),function(j){
        if(tabprop[j,i]!=0){
          clustermeasure<<-clustermeasure+(-tabprop[j,i]^2)
        }
        maxclustermeasure<<-maxclustermeasure+(-1/k^2)
        NULL
      })
      measure<<-measure+clustersizesnorm[i]*clustermeasure
      maxmeasure<<-maxmeasure+clustersizesnorm[i]*maxclustermeasure
      NULL
    })
  }
  list(result=measure,measure=mode,normalized_measure=measure/maxmeasure)
}
