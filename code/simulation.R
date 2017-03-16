###########################################################################
# Survival Simulation (Cox)
###########################################################################

design.cov <- function(p, cov.x, equalCorr=F) {
  if (!equalCorr){
    x <- matrix(0,p,p)
    for(i in 1:p){
      x[i,(i:p)] <- cov.x^(0:(p - i))
      x[i,(1:(i - 1))] <- x[(1:(i - 1)),i]
    }
  }
  else{
    x <- matrix(cov.x,p,p)
    diag(x) <- 1
  }
  return(x)
}

mvngeneration <- function(n, p, varcov) {
  ## function used to generate mvn dist with given mean vector and variance
  ## covariance matrix (here mean vector is all zeros)
  ## n is the sample size
  ## p the dimension
  choleski <- chol(varcov)
  tcholeski <- t(choleski)
  z <- matrix(rnorm(p * n, 0, 1), p, n)
  tranx <- tcholeski %*% z
  x <- t(tranx)
  x
}

simmOneCox <- function (n, p, cov.x, s, signal) {
  x <- mvngeneration(n, p, design.cov(p, cov.x))
  beta.true <- rep(0, p)
  x.center <- round((p-s)/2)
  beta.true[x.center+(1:s)] <- signal
  linPred <- x %*% beta.true
  Time <- round(sapply(linPred, function(lp){rexp(n=1,exp(-lp))}), 2)
  C <- round(rexp(n, exp(-mean(linPred))), 2)
  dat <- data.frame(time=pmax(1e-3, pmin(Time,C)), event=as.numeric(Time<=C), x)
  return(dat)
}


###########################################################################
# Regression & Classification Simulation
###########################################################################

simmOne <- function (n, p) {
  x <- matrix(data = rnorm(n*p), nrow=n, ncol=p)
  coeff <- c(1:p)
  resp <- x %*% coeff
  dat <- data.frame(outcome = resp, x)
  return(dat)
}


get.data <- function() {

    ## Define the family.
    family = c("surv", "regr", "class")[1]

    ## IF CLASSIFICATION:
    ## Number of levels in classification response
    c          <- 4

    ## IF SURVIVAL:
    ## Survival Simulation Constants
    s <- 5
    cov.x <- 0
    signal <- c(0, 1)[2]

    ## DIMENSIONS:
    p          <- 300
    n          <- 6000

    if (family == "surv") {

        sim.data.org  <- simmOneCox(n, p, cov.x, s, signal)

        ## f <- as.formula(Surv(time, event) ~ .)
    } else {

        sim.data.org  <- simmOne(n, p)

        if (family == "class") {
            outcome.cut <- cut(sim.data.org$outcome, breaks = c)
            sim.data.org$outcome <- outcome.cut
        }

        ## f <- as.formula(outcome ~ .)
    }

    return (sim.data.org)

}
