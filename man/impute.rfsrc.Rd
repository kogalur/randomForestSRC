\name{impute.rfsrc}
\alias{impute.rfsrc}
\alias{impute}
\title{Impute Only Mode}
\description{
  Fast imputation mode. A random forest is grown and used to impute
  missing data. No ensemble estimates or error rates are calculated.
  Optionally, a final sweep can be performed to re-fit each variable
  that had original missingness on the final covariates and overwrite
  only its originally-missing entries.
}
\usage{\method{impute}{rfsrc}(formula, data,
  ntree = 100, nodesize = 1, nsplit = 10,
  nimpute = 2, fast = FALSE, blocks,
  mf.q, max.iter = 10, eps = 0.01,
  ytry = NULL, always.use = NULL, verbose = TRUE,
  full.sweep = FALSE,
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

  \item{full.sweep}{If \code{TRUE}, performs a final sweep after the
    main imputation (both standard and missForest). For each variable
    that had any \emph{original} missingness, a forest is fit on rows
    where the variable was observed using the \emph{final} imputed
    covariates; predictions are then written back only to the originally
    missing cells. This can improve self-consistency across variables at
    the cost of extra computation.}

  \item{...}{Additional arguments passed to or from methods.
    Recognized advanced options include \code{full.sweep.options} (a
    \code{list}) controlling the final sweep hyperparameters:
    \code{ntree} (default \code{500}), \code{nodesize} (default
    \code{NULL}), \code{nsplit} (default \code{10}); and the standard
    \code{rfsrc} controls \code{mtry}, \code{splitrule}, \code{bootstrap},
    \code{sampsize}, \code{samptype} that apply to the sweep.}
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

  \item \emph{Optional final sweep}: if \code{full.sweep = TRUE}, a
  post-imputation sweep is performed for every variable with original
  missingness. Each such variable is re-fit on its observed rows using
  the final imputed covariates, and predictions overwrite only the
  originally missing entries. Defaults for the sweep are
  \code{ntree = 500}, \code{nodesize = NULL}, \code{nsplit = 10}, and
  can be customized via \code{full.sweep.options} passed through
  \code{...}. This applies to both standard and missForest modes.

  }
}

\value{
  Invisibly, the data frame containing the original data with imputed
  data overlaid.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
  (2008). Random survival forests, \emph{Ann. Appl. Stat.}, 2:841--860.

  Stekhoven D.J. and Buhlmann P. (2012). MissForest--non-parametric
  missing value imputation for mixed-type data.
  \emph{Bioinformatics}, 28(1):112--118.

  Tang F. and Ishwaran H. (2017). Random forest missing data
  algorithms. \emph{Statistical Analysis and Data Mining}, 10:363--377.
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

## optional final sweep (standard imputation)
pbc3.fs <- impute(f, data = pbc, splitrule = "random", nimpute = 5,
                  full.sweep = TRUE)

## ------------------------------------------------------------
## example of regression imputation
## ------------------------------------------------------------

air1.d <- impute(data = airquality, nimpute = 5)
air2.d <- impute(Ozone ~ ., data = airquality, nimpute = 5)
air3.d <- impute(Ozone ~ ., data = airquality, fast = TRUE)

## final sweep with custom options (e.g., larger forest)
air3.fs <- impute(Ozone ~ ., data = airquality, nimpute = 5,
                  full.sweep = TRUE,
                  full.sweep.options = list(ntree = 1000, nodesize = 5, nsplit = 0,
                                            mtry = 3, splitrule = "random"))

## ------------------------------------------------------------
## multivariate missForest imputation
## ------------------------------------------------------------

data(pbc, package = "randomForestSRC")

## missForest algorithm - uses 1 variable at a time for the response
pbc.d <- impute(data = pbc, mf.q = 1)

## multivariate missForest - use 10 percent of variables as responses
pbc.mv <- impute(data = pbc, mf.q = .10)

## missForest but faster by using random splitting
pbc.fast <- impute(data = pbc, mf.q = 1, splitrule = "random")

## missForest + final sweep
pbc.fast.fs <- impute(data = pbc, mf.q = 1, splitrule = "random",
                      full.sweep = TRUE)

}
}
\keyword{missing data}
