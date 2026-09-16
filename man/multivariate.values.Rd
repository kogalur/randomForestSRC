\name{multivariate.values}
\alias{multivariate.values}
\alias{get.mv.predicted}
\alias{get.mv.error}
\alias{get.mv.error.block}
\alias{get.mv.vimp}
\alias{get.mv.cserror}
\alias{get.mv.csvimp}
\alias{get.mv.formula}
\title{Extracting Multivariate Values}
\description{
  Extract predictions, performance errors, variable importance (VIMP),
  and case-specific values from fitted random forests.  The helpers
  collect response-specific results from multivariate regression,
  multivariate classification, and mixed-outcome forests.  They also
  provide a common interface for univariate results.  A helper for
  constructing a multivariate formula is included.
}
\usage{
get.mv.predicted(obj, oob = TRUE)

get.mv.error(obj, standardize = FALSE, pretty = TRUE, block = FALSE)

get.mv.error.block(obj, standardize = FALSE)

get.mv.vimp(obj, standardize = FALSE, pretty = TRUE)

get.mv.cserror(obj, standardize = FALSE)

get.mv.csvimp(obj, standardize = FALSE)

get.mv.formula(ynames)
}
\arguments{
  \item{obj}{An object returned by \code{rfsrc()} or
    \code{predict.rfsrc()} containing the values to extract.  For
    multivariate and mixed outcomes, retain the response columns in
    \code{obj$yvar}, including factor levels: these identify each
    response's type.  VIMP and case-specific quantities must already be
    present in the object.  See Details for survival and unavailable
    components.}
  \item{oob}{If \code{TRUE}, extract \code{predicted.oob} for each
    response when that component is present; otherwise use
    \code{predicted}.  If \code{FALSE}, use \code{predicted}.  Selection
    is made for the entire response component, not separately for each
    observation.}
  \item{standardize}{If \code{TRUE}, divide regression errors and VIMP
    by the variance of the corresponding response in \code{obj$yvar},
    calculated with \code{na.rm = TRUE}.  Classification and survival
    values are unchanged.  The treatment of a zero or unavailable
    variance is described under Standardization.}
  \item{pretty}{If \code{TRUE}, return a compact error vector or VIMP
    matrix, retaining only the overall \code{all} entry for each
    classification response.  If \code{FALSE}, return a response-named
    list and retain the class-specific entries.  Competing-risk
    event-specific entries are retained in either form.  This argument
    does not round the values.}
  \item{block}{If \code{TRUE}, extract the stored sequence of block
    errors from \code{err.block.rate}, rather than the final entry of
    \code{err.rate}.  This forces \code{pretty = FALSE}.
    \code{get.mv.error.block()} is a shorthand for this choice.}
  \item{ynames}{Character vector of response names for
    \code{get.mv.formula()}.  These should name columns in the data
    supplied to the subsequent forest fit.}
}
\details{
  \subsection{Predictions}{
    \code{get.mv.predicted()} combines the stored predictions into a
    matrix with observations in rows.  Each regression response supplies
    one column.  Each classification response supplies one probability
    column per class, named \code{response.class}.  The responses follow
    \code{obj$yvar.names}; class columns retain their stored order.

    The default uses OOB predictions when available.  An existing OOB
    component is retained even if some or all of its entries are missing.
    The fallback to \code{predicted} occurs only when that entire
    component is \code{NULL}.  Use \code{oob = FALSE} to explicitly
    select full-ensemble or new-data predictions.

    For right-censored survival, the helper extracts the stored mortality
    prediction.  For competing risks, it extracts the event-specific
    \code{predicted} values, with names of the form
    \code{response.event}.  It does not extract the time-indexed survival,
    cumulative hazard, or cumulative incidence arrays.
  }
  \subsection{Performance errors and variable importance}{
    \code{get.mv.error()} extracts the final stored performance error
    for each response.  The error measure is the one used when the object
    was fitted or predicted; the helper does not recalculate it or change
    \code{perf.type}.  With \code{pretty = FALSE}, classification output
    includes \code{all} and the class-specific errors.  Survival output
    retains its stored error columns, including event-specific columns
    for competing risks.

    \code{get.mv.error.block()} returns the full stored block-error
    sequence for each response.  It is equivalent to
    \code{get.mv.error(obj, standardize = standardize, block = TRUE)}.
    It does not choose a block size or form new blocks.

    \code{get.mv.vimp()} collects the stored \code{importance}
    components.  With \code{pretty = TRUE}, predictors index rows and
    results are combined across responses in columns.  With
    \code{pretty = FALSE}, each response has its own matrix, preserving
    classification \code{all} and class-specific columns and
    competing-risk event columns.  Request importance during fitting
    or use \code{\link{vimp}} before extracting it.

    To select a response or predictor, subset the returned vector,
    matrix, or response-named list.  These extraction calls do not fit
    trees, run prediction, or average errors or VIMP across responses.
  }
  \subsection{Case-specific values}{
    \code{get.mv.cserror()} obtains the case-specific error as
    \code{cse.num / cse.den} from the stored response component.
    \code{get.mv.csvimp()} similarly obtains case-specific VIMP as
    \code{csv.num / csv.den}.  These quantities must have been retained
    during fitting or prediction; the helpers do not generate them from
    the ensemble predictions.  In particular, case-specific error is
    not obtained by applying a loss to \code{predicted.oob} and the
    observed response.

    With one response, the case-specific values are returned directly;
    with multiple responses, they are returned in a response-named list.
    Case-specific VIMP has observations in rows and the stored VIMP
    variables in columns.  Variable names are attached when available
    from the \code{importance} component.  Both helpers return
    \code{NULL} for right-censored survival and competing risks.
  }
  \subsection{Standardization}{
    For a regression response \eqn{Y}, \code{standardize = TRUE}
    divides each extracted error or importance value by
    \code{var(Y, na.rm = TRUE)}.  This uses the response values on the
    supplied object: training responses for a grow result and evaluation
    responses for a test-prediction result.  It leaves predictions and
    all classification and survival values unchanged.

    \code{get.mv.error()} and \code{get.mv.vimp()} divide by this
    variance directly; a zero or unavailable variance can therefore
    produce nonfinite values.  The two case-specific helpers instead
    use divisor one when that variance is zero or \code{NA}.  No
    standardization is applied by default.
  }
  \subsection{Availability of stored results}{
    Supply an object containing the requested output for its responses.
    Missing predictions and numerical \code{NA} values remain as stored.
    \code{get.mv.error()} returns \code{NULL} when all response errors
    are absent.  Otherwise its compact output uses \code{NA} for an
    absent response error, while list output retains a \code{NULL} entry.

    \code{get.mv.vimp()}, \code{get.mv.cserror()}, and
    \code{get.mv.csvimp()} return \code{NULL} when the first response
    has no corresponding component, even if a later response has one.
    When a response-named list is returned, absent later components
    remain \code{NULL}.
  }
  \subsection{Constructing a multivariate formula}{
    \code{get.mv.formula(ynames)} returns a formula of the form
    \code{Multivar(y1, y2, ...) ~ .}.  The named responses may be
    continuous, factors, or a mixture; the forest fit determines their
    types from the supplied data.  The dot denotes the remaining data
    columns.  This helper constructs the formula only.
  }
}
\value{
  \describe{
    \item{\code{get.mv.predicted}}{A numeric matrix with one row per
      observation and columns for the response predictions, class
      probabilities, or event-specific predictions.  A single prediction
      column remains a matrix.}
    \item{\code{get.mv.error}}{A named numeric vector by default, or a
      response-named list when \code{pretty = FALSE}.  Without blocking,
      each entry is the final stored error or final error row.  With
      \code{block = TRUE}, list entries contain the stored block-error
      sequences.  Returns \code{NULL} if all error components are absent.}
    \item{\code{get.mv.error.block}}{A response-named list of block-error
      vectors or matrices, or \code{NULL} when all are absent.}
    \item{\code{get.mv.vimp}}{A numeric matrix by default, or a
      response-named list of matrices when \code{pretty = FALSE}.
      Returns \code{NULL} under the availability rule described above.}
    \item{\code{get.mv.cserror}}{Case-specific error values for a single
      response, retaining their stored vector or array dimensions, or a
      response-named list for multiple responses.  Returns \code{NULL}
      for survival families or under the availability rule above.}
    \item{\code{get.mv.csvimp}}{A case-by-variable matrix for a single
      response, or a response-named list of these matrices for multiple
      responses.  Returns \code{NULL} for survival families or under
      the availability rule above.}
    \item{\code{get.mv.formula}}{An R formula with the supplied response
      names on the left of \code{~} and a dot on the right.}
  }
}
\seealso{
  \code{\link{rfsrc}}, \code{\link{predict.rfsrc}},
  \code{\link{vimp}}, \code{\link{subsample}},
  \code{\link{classification.performance}}
}
\examples{
## ------------------------------------------------------------
## A basic multivariate analysis
## ------------------------------------------------------------
o <- rfsrc(cbind(Ozone, Temp) ~ ., data = na.omit(airquality))
print(head(get.mv.predicted(o)))
print(get.mv.error(o))

\donttest{
## ------------------------------------------------------------
## Select a response from the stored results
## ------------------------------------------------------------
pred.oob <- get.mv.predicted(o)
print(head(pred.oob[, "Temp", drop = FALSE]))
print(get.mv.error(o)["Temp"])
print(get.mv.error(o, standardize = TRUE))
print(head(get.mv.predicted(o, oob = FALSE)))

## ------------------------------------------------------------
## Formula construction, VIMP, and block errors
## ------------------------------------------------------------
f <- get.mv.formula(c("Ozone", "Temp"))
print(f)

o.vimp <- rfsrc(f, data = na.omit(airquality), ntree = 100,
                importance = "permute", block.size = 10)
print(get.mv.vimp(o.vimp))
print(get.mv.vimp(o.vimp, standardize = TRUE))
print(get.mv.vimp(o.vimp, pretty = FALSE)[["Temp"]])
print(head(get.mv.error.block(o.vimp)[["Temp"]]))

## Optional case-specific values are NULL when not saved in the object.
print(get.mv.cserror(o.vimp))
print(get.mv.csvimp(o.vimp))

## ------------------------------------------------------------
## Mixed outcomes: retain class-specific entries
## ------------------------------------------------------------
f.mix <- get.mv.formula(c("Sepal.Length", "Species"))
mix <- rfsrc(f.mix, data = iris, ntree = 100,
             importance = "permute", block.size = 10)
print(colnames(get.mv.predicted(mix)))
print(get.mv.error(mix))
print(get.mv.error(mix, pretty = FALSE)[["Species"]])
print(get.mv.vimp(mix, pretty = FALSE)[["Species"]])

## ------------------------------------------------------------
## Extract predictions and errors for held-out observations
## ------------------------------------------------------------
dta <- na.omit(airquality)
set.seed(17)
train <- sample(seq_len(nrow(dta)), floor(.7 * nrow(dta)))
fit <- rfsrc(f, data = dta[train, ], ntree = 100)
p.test <- predict(fit, newdata = dta[-train, ])
print(head(get.mv.predicted(p.test, oob = FALSE)))
print(get.mv.error(p.test))
}
}
\keyword{models}
