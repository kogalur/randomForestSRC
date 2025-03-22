\name{impute.rfsrc}
\alias{impute.rfsrc}
\alias{impute}
\title{Impute Only Mode}
\description{
  Fast imputation mode.  A random forest is grown and used to impute
  missing data.  No ensemble estimates or error rates are calculated.
}
\usage{\method{impute}{rfsrc}(formula, data,
  ntree = 100, nodesize = 1, nsplit = 10,
  nimpute = 2, fast = FALSE, blocks, 
  mf.q, max.iter = 10, eps = 0.01, 
  ytry = NULL, always.use = NULL, verbose = TRUE,
  ...)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit.  Can be
    left unspecified if there are no outcomes or we don't care to
    distinguish between y-outcomes and x-variables in the imputation.
    Ignored when using multivariate missForest imputation.}
  \item{data}{Data frame containing the data to be imputed.}
  \item{ntree}{Number of trees to grow.}
  \item{nodesize}{Forest average terminal node size.}
  \item{nsplit}{Non-negative integer value used to specify random splitting.}
  \item{nimpute}{Number of iterations of the missing data algorithm.
    Ignored for multivariate missForest; in which case the algorithm
    iterates until a convergence criteria is achieved (users can
    however enforce a maximum number of iterations with the option
    \code{max.iter}).}
  \item{fast}{Use fast random forests, \code{rfsrcFast}, in place of
    \code{rfsrc}?  Improves speed but is less accurate.}
  \item{blocks}{Integer value specifying the number of blocks the data
    should be broken up into (by rows).  This can improve computational
    efficiency when the sample size is large but imputation efficiency
    decreases.  By default, no action is taken if left unspecified.}
  \item{mf.q}{Use this to turn on missForest (which is off by default).
    Specifies fraction of variables (between 0 and 1) used as responses
    in multivariate missForest imputation.  When set to 1 this
    corresponds to missForest, otherwise multivariate missForest is
    used. Can also be an integer, in
    which case this equals the number of multivariate responses.}
  \item{max.iter}{Maximum number of iterations used when implementing
    multivariate missForest imputation.}
  \item{eps}{Tolerance value used to determine convergence of
    multivariate missForest imputation.}
  \item{ytry}{Number of variables used as pseudo-responses in
    unsupervised forests.  See details below.}
  \item{always.use}{Character vector of variable names to always
    be included as a response in multivariate missForest imputation.
    Does not apply for other imputation methods.}
  \item{verbose}{Send verbose output to terminal (only applies to
    multivariate missForest imputation).}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  \enumerate{
    
    \item Grow a forest and use this to impute data.  All external
    calculations such as ensemble calculations, error rates, etc. are
    turned off.  Use this function if your only interest is imputing the
    data.

    \item Split statistics are calculated using non-misssing data only.
    If a node splits on a variable with missing data, the variable's
    missing data is imputed by randomly drawing values from non-missing
    in-bag data.  The purpose of this is to make it possible to assign
    cases to daughter nodes based on the split.
    
    \item If no formula is specified, unsupervised splitting is
    implemented using a \code{ytry} value of sqrt(\code{p}) where
    \code{p} equals the number of variables.  More precisely,
    \code{mtry} variables are selected at random, and for each of these
    a random subset of \code{ytry} variables are selected and defined as
    the multivariate pseudo-responses.  A multivariate composite
    splitting rule of dimension \code{ytry} is then applied to each of
    the \code{mtry} multivariate regression problems and the node split
    on the variable leading to the best split (Tang and Ishwaran, 2017).
    
    \item If \code{mf.q} is specified, a multivariate version of
    missForest imputation (Stekhoven and Buhlmann, 2012) is applied.
    Specifically, a fraction \code{mf.q} of variables are used as
    multivariate responses and split by the remaining variables using
    multivariate composite splitting (Tang and Ishwaran, 2017).  Missing
    data for responses are imputed by prediction.  The process is
    repeated using a new set of variables for responses (mutually
    exclusive to the previous fit), until all variables have been
    imputed.  This is one iteration.  The entire process is repeated,
    and the algorithm iterated until a convergence criteria is met
    (specified using options \code{max.iter} and \code{eps}).  Integer
    values for \code{mf.q} are allowed and interpreted as a request that
    \code{mf.q} variables be selected for the multivariate response.  If
    \code{mf.q=1}, the algorithm reverts to the original missForest
    procedure. This is generally the most accurate of all the
    imputation procedures, but also the most computationally demanding.
    See examples below for strategies to increase speed.

    \item Prior to imputation, the data is processed and records with
    all values missing are removed, as are variables having all missing
    values.
    
    \item If there is no missing data, either before or after processing
    of the data, the algorithm returns the processed data and no
    imputation is performed.

    \item All options are the same as \command{rfsrc} and the user should
    consult the \command{rfsrc} help file for details.
  }
}
\value{
  Invisibly, the data frame containing the orginal data with imputed
  data overlaid.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
  (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.
  
  Stekhoven D.J. and Buhlmann P. (2012). MissForest--non-parametric
  missing value imputation for mixed-type data.
  \emph{Bioinformatics}, 28(1):112-118.
  
  Tang F. and Ishwaran H. (2017).  Random forest missing data
  algorithms.  \emph{Statistical Analysis and Data Mining}, 10:363-377.
  
}
\seealso{
  \command{\link{rfsrc}},
  \command{\link{rfsrc.fast}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## example of survival imputation
## ------------------------------------------------------------

## default everything - unsupervised splitting
data(pbc, package = "randomForestSRC")
pbc1.d <- impute(data = pbc)

## imputation using outcome splitting
f <- as.formula(Surv(days, status) ~ .)
pbc2.d <- impute(f, data = pbc, nsplit = 3)

## random splitting can be reasonably good
pbc3.d <- impute(f, data = pbc, splitrule = "random", nimpute = 5)

## ------------------------------------------------------------
## example of regression imputation
## ------------------------------------------------------------

air1.d <- impute(data = airquality, nimpute = 5)
air2.d <- impute(Ozone ~ ., data = airquality, nimpute = 5)
air3.d <- impute(Ozone ~ ., data = airquality, fast = TRUE)

## ------------------------------------------------------------
## multivariate missForest imputation
## ------------------------------------------------------------

data(pbc, package = "randomForestSRC")

## missForest algorithm - uses 1 variable at a time for the response
pbc.d <- impute(data = pbc, mf.q = 1)

## multivariate missForest - use 10 percent of variables as responses
## i.e. multivariate missForest
pbc.d <- impute(data = pbc, mf.q = .01)

## missForest but faster by using random splitting
pbc.d <- impute(data = pbc, mf.q = 1, splitrule = "random")

## missForest but faster by increasing nodesize
pbc.d <- impute(data = pbc, mf.q = 1, nodesize = 20, splitrule = "random")

## missForest but faster by using rfsrcFast
pbc.d <- impute(data = pbc, mf.q = 1, fast = TRUE)

## ------------------------------------------------------------
## another example of multivariate missForest imputation
## (suggested by John Sheffield)
## ------------------------------------------------------------

test_rows <- 1000

set.seed(1234)

a <- rpois(test_rows, 500)
b <- a + rnorm(test_rows, 50, 50)
c <- b + rnorm(test_rows, 50, 50)
d <- c + rnorm(test_rows, 50, 50)
e <- d + rnorm(test_rows, 50, 50)
f <- e + rnorm(test_rows, 50, 50)
g <- f + rnorm(test_rows, 50, 50)
h <- g + rnorm(test_rows, 50, 50)
i <- h + rnorm(test_rows, 50, 50)

fake_data <- data.frame(a, b, c, d, e, f, g, h, i)

fake_data_missing <- data.frame(lapply(fake_data, function(x) {
  x[runif(test_rows) <= 0.4] <- NA
  x
}))
  
imputed_data <- impute(
  data = fake_data_missing,
  mf.q = 0.2,
  ntree = 100,
  fast = TRUE,
  verbose = TRUE
)

par(mfrow=c(3,3))
o=lapply(1:ncol(imputed_data), function(j) {
  pt <- is.na(fake_data_missing[, j])
  x <- fake_data[pt, j]
  y <- imputed_data[pt, j]
  plot(x, y, pch = 16, cex = 0.8, xlab = "raw data",
    ylab = "imputed data", col = 2)
  points(x, y, pch = 1, cex = 0.8, col = gray(.9))
  lines(supsmu(x, y, span = .25), lty = 1, col = 4, lwd = 4)
  mtext(colnames(imputed_data)[j])
  NULL
})

}
}
\keyword{missing data}

