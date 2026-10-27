\name{impute.rfsrc}
\alias{impute.rfsrc}
\alias{impute}
\title{Random Forest Imputation}
\description{
  Imputes missing values using random forests. The default is
  on-the-fly imputation (OTFI), with supervised or unsupervised
  splitting, including an optional autoencoder. Alternatively,
  missForest predicts missing values one variable or group of variables
  at a time. Returns an imputed data frame.
}
\usage{\method{impute}{rfsrc}(formula, data,
  ntree = 100, nodesize = 1, nsplit = 10,
  nimpute = 2, fast = FALSE, blocks,
  mf.q, max.iter = 10, eps = 0.01,
  ytry = NULL, always.use = NULL, verbose = TRUE,
  full.sweep = FALSE, restore.integer = TRUE,
  ...)
}
\arguments{
  
  \item{formula}{Model formula specifying the responses and predictors
    for OTFI. Omit it for unsupervised splitting, when no outcome is
    designated. Ignored for missForest.}

  \item{data}{Data to be imputed. Matrices, tibbles, and
    \code{data.table} objects are converted to a plain data frame.
    Column names must be unique and nonempty, each column must be a
    vector, and numeric values must be finite or missing. Bare vectors
    and nested matrix or list columns are not supported.}

  \item{ntree}{Number of trees grown for each imputation forest.}

  \item{nodesize}{Minimum terminal node size in each tree.}

  \item{nsplit}{Non-negative integer for specifying random splitting.}

  \item{nimpute}{Number of OTFI iterations. Ignored for missForest,
    whose iterations are controlled by \code{max.iter} and \code{eps}.}

  \item{fast}{If \code{TRUE}, uses fast forest settings to reduce
    computation. May reduce imputation accuracy.}

  \item{blocks}{Number of row blocks for the main imputation fits.
    May improve speed for large data, but can reduce imputation
    accuracy. If omitted, fits use all rows. See Details.}

  \item{mf.q}{Selects missForest when supplied. With \code{mf.q = 1}
    and \code{always.use = NULL}, variables are imputed one at a time.
    Other positive values select multivariate missForest and specify an
    approximate fraction (below 1) or number (above 1) of variables with
    missing values to include in each response group. Takes precedence
    over \code{splitrule = "auto"}. See Details.}

  \item{max.iter}{Maximum number of missForest iterations. Must be a
    positive integer.}

  \item{eps}{Minimum decrease in the change between successive
    missForest imputations required to continue iterating. Must be
    finite and nonnegative. See Details.}

  \item{ytry}{Number of response variables used for unsupervised OTFI.
    With the default splitting method, pseudo-responses are selected
    for each candidate splitting variable; the default is
    \code{min(p - 1, max(25, ceiling(sqrt(p))))}. With
    \code{splitrule = "auto"}, response copies are selected at each
    node; the default is \code{p}, and a smaller positive integer can
    be supplied. Here \code{p} is the number of variables remaining
    after all-missing columns are removed. Used when \code{formula}
    and \code{mf.q} are omitted.}

  \item{always.use}{Names of variables to include in every response
    group for multivariate missForest. With \code{mf.q} supplied, any
    non-\code{NULL} value selects the multivariate version, even if
    \code{mf.q = 1}. Empty vectors and names not present in the data
    add no responses but still select this version. Ignored for OTFI.}

  \item{verbose}{If \code{TRUE}, prints progress during missForest
    imputation and the optional final sweep.}

  \item{full.sweep}{If \code{TRUE}, performs one additional pass after
    either imputation method. Each variable with missing values in the
    original data is predicted from the current imputed values of the
    other variables. Only originally missing entries are replaced.
    See Details.}

  \item{restore.integer}{If \code{TRUE}, rounds imputed values to
    integers for columns stored as R integers or whose observed numeric
    values are all integer-valued. Observed values are preserved.
    Numeric columns retain numeric storage. Original integer columns
    are returned as integer vectors only when all nonmissing completed
    values are integral and within R's integer range; otherwise they
    retain numeric storage. Set to \code{FALSE} to keep unrounded
    imputations. See Details for the integer-detection tolerance.}

  \item{...}{Additional forest controls for the main imputation:
    \code{mtry}, \code{splitrule}, \code{bootstrap}, \code{sampsize},
    and \code{samptype}. Use \code{splitrule = "auto"} for autoencoder
    splitting when \code{formula} and \code{mf.q} are omitted.
    Supplying \code{mf.q} selects missForest and ignores the
    \code{"auto"} setting. Controls for the optional final sweep are
    supplied separately in the named list \code{full.sweep.options}.
    This list accepts the same controls, plus \code{ntree} (default
    \code{500}), \code{nodesize} (default \code{NULL}), and
    \code{nsplit} (default \code{10}).}
}

\details{
  Two approaches are available. The default, on-the-fly imputation
  (OTFI), handles missing values during tree construction and obtains
  final imputations from terminal nodes. Supplying \code{mf.q} selects
  missForest, which repeatedly uses the other variables to impute each
  variable or group of variables with missing values. Both approaches
  support numeric and categorical data. The function returns the
  imputed data frame; fitted forests and performance summaries are not
  returned.

  \subsection{On-the-fly imputation (OTFI)}{

     OTFI is the approach underlying \code{na.action = "na.impute"} in
    \code{\link{rfsrc}} and \code{\link{predict.rfsrc}}.  Tang and
    Ishwaran (2017, Sections 2.3--2.4) describe the method for both the
    supervised and unsupervised case.

    Split statistics are calculated from the observed values of the
    variables involved. If an observation is missing the value of the
    chosen splitting variable, a value is drawn from that variable's
    observed inbag values in the node to determine its daughter node.
    This draw is temporary: it is discarded after assignment, and the
    entry is treated as missing again at the next node. Final
    replacements are obtained separately from observed values in the
    terminal nodes and combined across trees.

    With a response formula, OTFI uses supervised splitting guided by
    the specified outcome. For example, \code{Ozone ~ .} uses Ozone
    to guide the splits while imputing missing responses and
    predictors.

    With no formula, OTFI uses unsupervised splitting to learn
    relationships among the variables. At each node, each of \code{mtry}
    candidate splitting variables is evaluated against \code{ytry}
    randomly selected pseudo-responses using a multivariate composite
    splitting rule. The best split is selected.  

    Additionally, for unsupervised data, setting \code{splitrule =
    "auto"} implements a autoencoder type method which regresses
    variables against themselves. Temporary copies serve as multivariate
    responses, and the original variables serve as predictors. Missing
    entries remain missing in both copies for OTFI to process. All
    response variables are used by default (\code{ytry = p}). The
    temporary responses are discarded after imputation; only the
    completed original variables are returned or used in the optional
    final sweep. Other \code{splitrule} settings are ignored when the
    formula is omitted.

    \code{nimpute} controls the number of OTFI iterations. With a
    formula, \code{splitrule = "random"} can be used for completely
    random splitting.
  }

  \subsection{missForest}{
    missForest is an unsupervised imputation procedure in the sense
    that no outcome needs to be designated: every variable can be
    imputed using the others. Its individual forests use supervised
    regression or classification. The \code{formula} argument is
    ignored; an outcome column included in \code{data} is treated like
    any other variable.

    The implementation first fills missing values using a small forest
    with completely random splitting. These initial imputations supply
    values for the predictors in subsequent fits.

    With \code{mf.q = 1} and \code{always.use = NULL}, the method
    follows the variable-by-variable approach of Stekhoven and
    Buehlmann (2012). For each variable with missing values, a forest
    is fit using the rows where that variable was observed. Its
    missing entries are replaced by predictions from this forest.
    Later fits use the imputations already updated for earlier
    variables, and the process repeats.

    The multivariate extension groups variables with missing values
    and imputes each group using the remaining variables (Tang and
    Ishwaran, 2017, Section 2.5). Group sizes are approximate and
    controlled by \code{mf.q}. Within each group, response entries
    that were originally missing are set back to missing, while the
    current imputations remain available in the predictors. A forest
    with multivariate composite splitting imputes the responses.
    Groups are formed anew at each iteration. Variables named in
    \code{always.use} are added to every response group.

    For both versions, iteration is controlled by \code{max.iter} and
    \code{eps}. After each pass, changes at originally missing
    entries are summarized by factor disagreement and scaled numeric
    changes, averaged over the variables. Iteration stops when this
    measure fails to decrease by at least \code{eps}, or when
    \code{max.iter} is reached. The last completed pass is retained.
  }

  \subsection{Data handling and optional refinements}{
    Rows and columns with no observed values are removed before
    imputation. If no missing values remain, the data are returned
    without fitting a forest. If the entire data frame is missing,
    it is returned without imputation.

    For large data sets, \code{blocks} restricts the main imputation
    fits to separate groups of rows. This can reduce computation but
    limits the information available to each fit. Autoencoder fits use
    only response variables with observed values in that block and limit
    \code{ytry} accordingly. The missForest initialization and optional
    final sweep use all retained rows.

    Setting \code{full.sweep = TRUE} adds a final pass after either
    method. Each variable that originally had missing values is fit
    on its observed rows using the current imputed predictors. Only
    its originally missing entries are replaced, and later fits use
    earlier updates from the sweep. Sweep controls are supplied in
    \code{full.sweep.options}; see \code{...} above.

    With \code{restore.integer = TRUE}, integer support is determined
    from the original data before imputation. Numeric values are
    compared with the nearest integer using a fixed absolute
    floating-point tolerance. Rounding applies only to imputed
    entries: after the missForest initialization, before and after
    the optional final sweep, and before the data are returned.
    Values are not rounded between the individual missForest updates.
  }
}

\value{
  A data frame with missing entries replaced by imputed values and
  observed values preserved. All-missing rows and columns are removed
  as described in Details.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Stekhoven D.J. and Buehlmann P. (2012). MissForest--non-parametric
  missing value imputation for mixed-type data.
  \emph{Bioinformatics}, 28(1):112--118.
  \doi{10.1093/bioinformatics/btr597}.

  Tang F. and Ishwaran H. (2017). Random forest missing data
  algorithms. \emph{Statistical Analysis and Data Mining}, 10(6):363--377.
  \doi{10.1002/sam.11348}.
}
\seealso{
  \command{\link{impute.learn.rfsrc}},
  \command{\link{rfsrc}},
  \command{\link{rfsrc.fast}},
  \command{\link{predict.rfsrc}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## example of survival imputation
## ------------------------------------------------------------

## default OTFI: unsupervised splitting
data(pbc, package = "randomForestSRC")
pbc1.d <- impute(data = pbc)

## autoencoder OTFI: regress the variables against themselves
pbc.auto <- impute(data = pbc, splitrule = "auto")

## optionally use fewer response variables at each node
pbc.auto.sub <- impute(data = pbc, splitrule = "auto", ytry = 3)

## OTFI with splits guided by the survival outcome
f <- as.formula(Surv(days, status) ~ .)
pbc2.d <- impute(f, data = pbc, nsplit = 3)

## random splitting can be reasonably good
pbc3.d <- impute(f, data = pbc, splitrule = "random", nimpute = 5)

## optional final sweep after OTFI
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

## multivariate missForest - use about 10 percent of variables
## with missing values as responses in each group
pbc.mv <- impute(data = pbc, mf.q = .10)

## missForest but faster by using random splitting
pbc.fast <- impute(data = pbc, mf.q = 1, splitrule = "random")

## missForest + final sweep
pbc.fast.fs <- impute(data = pbc, mf.q = 1, splitrule = "random",
                      full.sweep = TRUE)

}
}
\keyword{missing data}
