\name{impute}
\alias{impute.rfsrc}
\alias{impute}
\title{Impute Only Mode}
\description{
  Fast imputation mode.  A random forest is grown and used to impute
  missing data.  No ensemble estimates or error rates are calculated.
}
\usage{\method{impute}{rfsrc}(formula, data, ntree = 500, mtry = NULL,
  xvar.wt = NULL, nodesize = 1, splitrule = NULL, nsplit = 10,
  
  na.action = c("na.impute"),
  
   
  nimpute = 2, mf.q, blocks,
  always.use = NULL, max.iter = 10, eps = 0.01, verbose = TRUE,
  do.trace = FALSE, ...)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit.  Can be
    left unspecified if there are no outcomes or we don't care to
    distinguish between y-outcomes and x-variables in the imputation.}
  \item{data}{Data frame containing the data to be imputed.}
  \item{ntree}{Number of trees to grow.}
  \item{mtry}{Number of variables randomly sampled at each split.}
  \item{nodesize}{Forest average terminal node size.}
  \item{splitrule}{Splitting rule used to grow trees.  Do not use for
    unsupervised or multivariate splitting unless
    \code{splitrule}="random" (used as a strategy for speeding up
    calculations).}
  \item{nsplit}{Non-negative integer value used to specify random splitting.}
  \item{na.action}{Missing value action. See details below.}
  \item{nimpute}{Number of iterations of the missing data algorithm.
	  Ignored for multivariate missForest; in which case the algorithm
	  iterates until a convergence criteria is achieved (users can
    however enforce a maximum number of iterations with the option
    \code{max.iter}).}
  \item{mf.q}{Fraction of variables (between 0 and 1) used as responses
    in multivariate missForest imputation.  By default, multivariate
    missForest imputation is not performed if left unspecifed.  Can be
    an integer, in which case this equals the number of multivariate
    responses.}
  \item{blocks}{Integer value specifying the number of blocks the data
    should be broken up into (by rows).  This can improve computational
    efficiency when the sample size is large but imputation efficiency
    decreases.  By default, no action is taken if left unspecified.}
  \item{always.use}{Character vector of variable names to always
    be included as a response in multivariate missForest imputation.
    Does not apply for other imputation methods.}
  \item{xvar.wt}{Weights for selecting variables for splitting on.}
  \item{max.iter}{Maximum number of iterations used when implementing
    multivariate missForest imputation.}
  \item{eps}{Tolerance value used to determine convergence of
    multivariate missForest imputation.}
  \item{verbose}{Send verbose output to terminal (only applies to
    multivariate missForest imputation).}
  \item{do.trace}{Number of seconds between updates to the user on
    approximate time to completion.}
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
    missForest imputation (Stekhoven and Buhlmann, 2012) is applied.  A
    fraction \code{mf.q} of variables are used as multivariate responses
    and split by the remaining variables using multivariate composite
    splitting (Tang and Ishwaran, 2017).  Missing data for responses are
    imputed by prediction.  The process is repeated using a new set of
    variables for responses (mutually exclusive to the previous fit),
    until all variables have been imputed.  This is one iteration.  The
    entire process is repeated, and the algorithm iterated until a
    convergence criteria is met (specified using options \code{max.iter}
    and \code{eps}).  Integer values for \code{mf.q} are allowed and
    interpreted as a request that \code{mf.q} variables be selected for
    the multivariate response.  This is generally the most accurate of
    all the imputation procedures, but also the most computationally
    demanding.
    
    \item Prior to imputation, the data is processed and records with
    all values missing are removed, as are variables having all missing
    values.
    
    \item If there is no missing data, either before or after processing
    of the data, the algorithm returns the processed data and no
    imputation is performed.
    
    \item The default choice \command{nimpute=2} is chosen for coherence
    with the default missing data algorithm implemented in grow mode.
    Thus, if the user imputes data with \command{nimpute=2} and runs a
    grow forest using this imputed data, then performance values such as
    VIMP and error rates will coincide with those obtained by running a
    grow forest on the original non-imputed data using \code{na.action =
      "na.impute"}.  Ignored for multivariate missForest.
    
    \item All options are the same as \command{rfsrc} and the user should
    consult the \command{rfsrc} help file for details.
  }
}   
\value{
  Invisibly, the data frame containing the orginal data with imputed
  data overlayed.
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
  algorithms.  \emph{Statistical Analysis and Data Mining}, 10, 363-377.
}
\seealso{
  \command{\link{rfsrc}}
}
\examples{
\dontrun{
## ------------------------------------------------------------
## example of survival imputation
## ------------------------------------------------------------

## unsupervised splitting
data(pbc, package = "_PROJECT_PACKAGE_NAME_")
pbc1.d <- impute(data = pbc, nimpute = 5)

#imputation using outcome splitting
f <- as.formula(Surv(days, status) ~ .)
pbc2.d <- impute(f, data = pbc, nsplit = 3)

#random splitting can be reasonably good
pbc4.d <- impute(f, data = pbc, splitrule = "random", nimpute = 5)

## ------------------------------------------------------------
## example of regression imputation
## ------------------------------------------------------------

air1.d <- impute(data = airquality, nimpute = 5)
air2.d <- impute(Ozone ~ ., data = airquality, nimpute = 5)
air3.d <- impute(Ozone ~ ., data = airquality, nimpute = 5,
           splitrule = "random", nodesize = 1)

## ------------------------------------------------------------
## multivariate missForest imputation
## ------------------------------------------------------------

data(pbc, package = "_PROJECT_PACKAGE_NAME_")

## use 10 percent of variables as responses
## i.e. multivariate missForest
pbc.d <- impute(data = pbc, mf.q = .01)

## use 1 variable as the response
## i.e. original missForest algorithm
pbc.d <- impute(data = pbc, mf.q = 1)

## faster call by using random splitting
pbc.d <- impute(data = pbc, mf.q = 1, splitrule = "random")

## faster call by increasing nodesize
pbc.d <- impute(data = pbc, mf.q = 1, nodesize = 20, splitrule = "random")


}
}
\keyword{missing data}

