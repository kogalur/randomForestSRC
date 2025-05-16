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
  \item{formula}{A symbolic model description. Can be omitted if
    outcomes are unspecified or if distinction between outcomes and
    predictors is unnecessary. Ignored for multivariate missForest.}
  
  \item{data}{A data frame containing variables to be imputed.}
  
  \item{ntree}{Number of trees grown for each imputation.}
  
  \item{nodesize}{Minimum terminal node size in each tree.}
  
  \item{nsplit}{Non-negative integer for specifying random splitting.}
  
  \item{nimpute}{Number of iterations for the missing data
  algorithm. Ignored for multivariate missForest, which iterates to
  convergence unless capped by \code{max.iter}.}

  \item{fast}{If \code{TRUE}, uses \code{rfsrcFast} instead of
    \code{rfsrc}. Increases speed but may reduce accuracy.}
   
  \item{blocks}{Number of row-wise blocks to divide the data into. May
  improve speed for large data, but can reduce imputation accuracy. No
  action if unspecified.}
  
  \item{mf.q}{Enables missForest. Either a fraction (between 0 and 1) of
  variables treated as responses, or an integer indicating number of
  response variables. \code{mf.q = 1} corresponds to standard
  missForest.}
  
  \item{max.iter}{Maximum number of iterations for multivariate missForest.}
  
  \item{eps}{Convergence threshold for multivariate missForest (change in imputed values).}
  
  \item{ytry}{Number of variables used as pseudo-responses in
    unsupervised forests. See Details.}
  
  \item{always.use}{Character vector of variables always included as
    responses in multivariate missForest. Ignored by other methods.}
  
  \item{verbose}{If \code{TRUE}, prints progress during multivariate
    missForest imputation.}
  
  \item{...}{Additional arguments passed to or from methods.}
  
}

\details{
  \enumerate{

  \item Before imputation, observations and variables with all values
  missing are removed.
  

  \item A forest is grown and used solely for imputation. No ensemble
  statistics (e.g., error rates) are computed. Use this function when
  imputation is the only goal.
  

  \item For standard imputation (not missForest), splits are based only
  on non-missing data. If a split variable has missing values, they are
  temporarily imputed by randomly drawing from in-bag, non-missing
  values to allow node assignment.
  

  \item If \code{mf.q} is specified, multivariate missForest imputation
  is applied (Stekhoven and B\"uhlmann, 2012). A fraction (or integer
  count) of variables are selected as multivariate responses, predicted
  using the remaining variables with multivariate composite
  splitting. Each round imputes a disjoint set of variables, and the
  full cycle is repeated until convergence, controlled by
  \code{max.iter} and \code{eps}. Setting \code{mf.q = 1} reverts to
  standard missForest. This method is typically the most accurate, but
  also the most computationally intensive.
  

  \item If no formula is provided, unsupervised splitting is used. The
  default \code{ytry} is \code{sqrt(p)}, where \code{p} is the number of
  variables. For each of \code{mtry} candidate variables, a random
  subset of \code{ytry} variables is selected as pseudo-responses. A
  multivariate composite splitting rule is applied, and the split is
  made on the variable yielding the best result (Tang and Ishwaran,
  2017).
  

  \item If no missing values remain after preprocessing, the function
  returns the processed data without further action.
  

  \item All standard \code{rfsrc} options apply; see examples below for
  illustration.
  

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

