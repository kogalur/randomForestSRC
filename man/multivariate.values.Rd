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
  and case-specific errors and VIMP from random forest fits and
  predictions.  The functions combine response-specific results from
  multivariate regression, multivariate classification, and mixed-outcome
  forests.  They also work with univariate forests.
  \code{get.mv.formula()} constructs a multivariate formula from
  response names.
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
    \code{predict.rfsrc()}.  For multivariate and mixed-outcome objects,
    \code{obj$yvar} must contain the response columns with their original
    types and factor levels.}
  \item{oob}{If \code{TRUE}, use \code{predicted.oob} for each response,
    falling back to \code{predicted} when that response's OOB component
    is \code{NULL}.  If \code{FALSE}, use \code{predicted}.}
  \item{standardize}{If \code{TRUE}, divide regression errors and VIMP
    by the corresponding response variance in \code{obj$yvar}.
    Classification and survival values are unchanged.  See
    Standardization for the treatment of zero or missing variances.}
  \item{pretty}{If \code{TRUE}, return errors as a named vector and
    VIMP as a matrix, using only the overall \code{all} value for each
    classification response.  If \code{FALSE}, return a list named by
    response, including class-specific values.  Competing-risk
    event-specific values are included in both formats.  Values are
    not rounded.}
  \item{block}{If \code{TRUE}, return the block-error sequence from
    \code{err.block.rate} instead of the final error from
    \code{err.rate}.  This sets \code{pretty = FALSE}.
    \code{get.mv.error.block()} is shorthand for this choice.}
  \item{ynames}{Character vector of response column names for
    \code{get.mv.formula()}.  These columns must be in the data used
    to fit the forest.}
}
\details{
  The extraction functions use results already computed during fitting
  or prediction.  Results are reported separately for each response.
  Subset the returned vector, matrix, or list to select responses or
  predictors.

  \subsection{Predictions}{
    \code{get.mv.predicted()} returns a matrix with observations in
    rows and predictions in columns:
    \itemize{
      \item Regression: one column per response.
      \item Classification: one probability column per class, named
      \code{response.class}.
      \item Right-censored survival: the mortality prediction.
      \item Competing risks: one prediction column per event, named
      \code{response.event}.
    }
    Responses follow \code{obj$yvar.names}; class and event columns
    follow their order in the object.  Time-indexed survival, cumulative
    hazard, and cumulative incidence arrays are not included.

    By default, OOB predictions are used for each response when its
    \code{predicted.oob} component is present.  Missing values within
    that component remain missing, even if every value is missing.
    The function uses \code{predicted} only when the entire OOB component
    is \code{NULL}.  Set \code{oob = FALSE} to use full-ensemble or
    new-data predictions.
  }
  \subsection{Performance errors and variable importance}{
    \code{get.mv.error()} returns the final performance error for each
    response, using the error measure selected during fitting or
    prediction.  With \code{pretty = FALSE}, classification results
    include \code{all} and the class-specific errors.  Survival results
    include all error columns, including event-specific errors for
    competing risks.

    \code{get.mv.error.block()} returns the full block-error sequence
    computed during fitting or prediction.  It is equivalent to
    \code{get.mv.error(obj, standardize = standardize, block = TRUE)};
    the block size is determined when those errors are computed.

    \code{get.mv.vimp()} extracts variable importance.  Request
    importance when fitting the forest or use \code{\link{vimp}}
    before calling this function.  With \code{pretty = TRUE}, the
    result is a matrix with predictors in rows and response results
    in columns.  With \code{pretty = FALSE}, each response has its own
    matrix, including classification \code{all} and class-specific
    columns and competing-risk event columns.
  }
  \subsection{Case-specific values}{
    \code{get.mv.cserror()} calculates case-specific error as
    \code{cse.num / cse.den}; \code{get.mv.csvimp()} calculates
    case-specific VIMP as \code{csv.num / csv.den}.  These numerator
    and denominator components must be saved in the fitted or prediction
    object.  Case-specific error is calculated from these components,
    not by applying a loss to \code{predicted.oob} and the observed
    response.

    With one response, values are returned directly; with multiple
    responses, they are returned in a list named by response.
    Case-specific VIMP has observations in rows and VIMP variables in
    columns, with variable names taken from \code{importance} when
    available.  Both functions return \code{NULL} for right-censored
    survival and competing risks.
  }
  \subsection{Standardization}{
    For a regression response \eqn{Y}, \code{standardize = TRUE}
    divides each error or VIMP value by \code{var(Y, na.rm = TRUE)}.
    The variance uses the response values in \code{obj$yvar}: training
    responses for a fitted object and evaluation responses for a
    test-prediction object.  Classification and survival values are
    unchanged.  Standardization is off by default and does not apply
    to predictions.

    \code{get.mv.error()}, \code{get.mv.error.block()}, and
    \code{get.mv.vimp()} divide by the variance directly, so a zero
    or missing variance can produce nonfinite values.  The case-specific
    functions use a divisor of one when the variance is zero or
    \code{NA}.
  }
  \subsection{Missing results}{
    A \code{NULL} component means that the corresponding result is
    absent from \code{obj}.  An \code{NA} entry is a missing value
    within an existing result.
    \itemize{
      \item \code{get.mv.error()} returns \code{NULL} when no response
      has an error result.  Otherwise, a response without an error is
      represented by \code{NA} in vector output or \code{NULL} in list
      output.  The same list behavior applies to
      \code{get.mv.error.block()}.
      \item \code{get.mv.vimp()}, \code{get.mv.cserror()}, and
      \code{get.mv.csvimp()} return \code{NULL} if the first response
      lacks the corresponding result, even when later responses have
      results.  If the first response has a result, missing later
      results are \code{NULL} in list output.
    }
  }
  \subsection{Constructing a multivariate formula}{
    \code{get.mv.formula(ynames)} returns
    \code{Multivar(y1, y2, ...) ~ .}.  Responses may be continuous,
    factors, or a mixture; their types are determined from the data
    when the forest is fitted.  The dot specifies the remaining data
    columns as predictors.
  }
}
\value{
  \describe{
    \item{\code{get.mv.predicted}}{A numeric matrix with one row per
      observation and columns for response predictions, class
      probabilities, or event-specific predictions.  A single-column
      result is also a matrix.}
    \item{\code{get.mv.error}}{A named numeric vector by default, or
      a list named by response when \code{pretty = FALSE}.  Entries
      contain the final error or error row.  With \code{block = TRUE},
      the result is a list of block-error sequences.  Returns
      \code{NULL} when no response has an error result.}
    \item{\code{get.mv.error.block}}{A list of block-error vectors or
      matrices named by response, or \code{NULL} when no response has
      block errors.}
    \item{\code{get.mv.vimp}}{A numeric matrix by default, or a list of
      matrices named by response when \code{pretty = FALSE}.  Returns
      \code{NULL} if the first response has no importance result.}
    \item{\code{get.mv.cserror}}{Case-specific error values for one
      response, with their original vector or array dimensions, or a
      list named by response for multiple responses.  Returns
      \code{NULL} for survival families or if the first response has
      no case-specific error result.}
    \item{\code{get.mv.csvimp}}{A case-by-variable matrix for one
      response, or a list of these matrices named by response for
      multiple responses.  Returns \code{NULL} for survival families
      or if the first response has no case-specific VIMP result.}
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
## Mixed outcomes: include class-specific entries
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
