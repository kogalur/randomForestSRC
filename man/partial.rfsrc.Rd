\name{partial.rfsrc}
\alias{partial.rfsrc}
\title{Acquire Partial Effect of a Variable}
\description{
  Acquire the partial effect of a variable on the ensembles.
}
\usage{
partial.rfsrc(object, outcome.target = NULL,
  partial.type = NULL, partial.xvar = NULL, partial.values = NULL,
  partial.xvar2 = NULL, partial.values2 = NULL,
  partial.time = NULL, oob = TRUE, seed = NULL, do.trace = FALSE, ...)
}
\arguments{
  \item{object}{An object of class \code{(rfsrc, grow)}.}
  \item{outcome.target}{Character value for multivariate families
    specifying the target outcome to be used. The default is to use the
    first coordinate.}
  \item{partial.type}{Character value of the type of predicted value.
    See details below.}
  \item{partial.xvar}{Character value specifying the single primary
    partial x-variable to be used.}
  \item{partial.values}{Vector of values that the primary partialy
    x-variable will assume.}
  \item{partial.xvar2}{Vector of character values specifying the second order x-variables to be
    used.}
  \item{partial.values2}{Vector of values that the second order
    x-variables will assume.  Each second order x-variable can only
    assume a single value.  This the length of \code{partial.xvar2} and
    \code{partial.values2} will be the same.  In addition, the user must
    do the appropriate conversion for factors, and represent a value as
    a numeric element.}
  \item{partial.time}{For survival families, the time at which the predicted
    survival value is evaluated at (depends on \code{partial.type}).}
  \item{oob}{OOB (TRUE)  or in-bag (FALSE) predicted values.}
  \item{seed}{Negative integer specifying seed for the random number
    generator.}
  \item{do.trace}{Number of seconds between updates to the user on
    approximate time to completion.}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{

  \enumerate{
    \item For regression, the predicted response is used.
    
    \item For survival, the choices are:
    \itemize{
      \item Relative frequency of mortality (\code{rel.freq}) or
      mortality (\code{mort}) is of dim \code{[n] x
        [length(partial.values)]}.
      \item The cumulative hazard function (\code{chf})
      is of dim \code{[n] x [length(partial.time)] x
        [length(partial.values)]}.
      \item The survival function (\code{surv}) is of dim \code{[n] x
        [length(partial.time)] x [length(partial.values)]}.
    }
    
    \item For competing risks, the choices are:
    \itemize{
	    \item The expected number of life years lost (\code{years.lost})
      is of dim \code{[n] x [length(event.info$event.type)] x
        [length(partial.values)]}.
	    \item The cumulative incidence function (\code{cif}) is of dim
	    \code{[n] x [length(partial.time)] x
	      [length(event.info$event.type)] x
	      [length(partial.values)]}.
	    \item The cumulative hazard function (\code{chf}) is of dim
      \code{[n] x [length(partial.time)] x [length(event.info$event.type)]
        x [length(partial.values)]}.
    }

    \item For regression, it is of dim \code{[n] x [length(partial.values)]}.
    \item For classification, it is of dim \code{[n] x [1 + yvar.nlevels[.]] x [length(partial.values)]}.
  }
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran H., Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7(2):25-31.
  
  Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
  (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.
  
 Ishwaran H., Gerds T.A., Kogalur U.B., Moore R.D., Gange S.J. and Lau
  B.M. (2014). Random survival forests for competing risks.
  \emph{Biostatistics}, 15(4):757-773.
}
\seealso{
  \command{\link{plot.variable.rfsrc}}
}
\examples{
\dontrun{
## ------------------------------------------------------------
## survival/competing risk
## ------------------------------------------------------------

## survival
data(veteran, package = "randomForestSRC")
v.obj <- rfsrc(Surv(time,status)~., veteran, nsplit = 10, ntree = 100)
partial.obj <- partial.rfsrc(v.obj,
  partial.type = "rel.freq",
  partial.xvar = "age",
  partial.values = v.obj$xvar[, "age"],
  partial.time = v.obj$time.interest)

## competing risks
data(follic, package = "randomForestSRC")
follic.obj <- rfsrc(Surv(time, status) ~ ., follic, nsplit = 3, ntree = 100)
partial.obj <- partial.rfsrc(follic.obj,
  partial.type = "cif",
  partial.xvar = "age",
  partial.values = follic.obj$xvar[, "age"],
  partial.time = follic.obj$time.interest,
  oob = TRUE)

## regression
airq.obj <- rfsrc(Ozone ~ ., data = airquality)
partial.obj <- partial.rfsrc(airq.obj,
  partial.xvar = "Wind",
  partial.values = airq.obj$xvar[, "Wind"],
  oob = FALSE)

## classification
iris.obj <- rfsrc(Species ~., data = iris)
partial.obj <- partial.rfsrc(iris.obj,
  partial.xvar = "Sepal.Length",
  partial.values = iris.obj$xvar[, "Sepal.Length"])

## multivariate mixed outcomes
mtcars2 <- mtcars
mtcars2$carb <- factor(mtcars2$carb)
mtcars2$cyl <- factor(mtcars2$cyl)
mtcars.mix <- rfsrc(Multivar(carb, mpg, cyl) ~ ., data = mtcars2)
partial.obj <- partial.rfsrc(mtcars.mix,
  partial.xvar = "disp",
  partial.values = mtcars.mix$xvar[, "disp"])

## second order variable specification
mtcars.obj <- rfsrc(mpg ~., data = mtcars)
partial.obj <- partial.rfsrc(mtcars.obj,
  partial.xvar = "cyl",
  partial.values = c(4, 8),
  partial.xvar2 = c("gear", "disp", "carb"),
  partial.values2 = c(4, 200, 3))

}
}
\keyword{partial}
