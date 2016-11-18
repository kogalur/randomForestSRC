\name{plot.variable}
\alias{plot.variable}
\alias{plot.variable.rfsrc}
\title{Plot Marginal Effect of Variables}
\description{
  Plot the marginal effect of an x-variable on the class probability
  (classification), response (regression), mortality (survival), or
  the expected years lost (competing risk) from a RF-SRC analysis.  Users
  can select between marginal (unadjusted, but fast) and partial plots
  (adjusted, but slow).
}
\usage{\method{plot.variable}{rfsrc}(x, xvar.names, which.class,
  outcome.target = NULL, time, surv.type = c("mort", "rel.freq",
  "surv", "years.lost", "cif", "chf"), class.type =
  c("prob", "bayes"), partial = FALSE, oob = TRUE, show.plots = TRUE,
  plots.per.page = 4, granule = 5, sorted = TRUE, nvar, npts = 25,
  smooth.lines = FALSE, subset, ...)
}
\arguments{
  \item{x}{An object of class \code{(rfsrc, grow)}, \code{(rfsrc, synthetic)},
    \code{(rfsrc, predict)}, or \code{(rfsrc, plot.variable)}.  See the
    examples below for illustration of the latter.}
  \item{xvar.names}{Names of the x-variables to be used.}
  \item{which.class}{For classification families, an integer or
    character value specifying the class to focus on (defaults to the
    first class).  For competing risk families, an integer value between
    1 and \code{J} indicating the event of interest, where \code{J} is
    the number of event types.  The default is to use the first event
    type.}
  \item{outcome.target}{Character value for multivariate families
    specifying the target outcome to be used. The default is to use the
    first coordinate.}
  \item{time}{For survival families, the time at which the predicted
    survival value is evaluated at (depends on \code{surv.type}).}
  \item{surv.type}{For survival families, specifies the predicted value.
    See details below.}
  \item{class.type}{For classification families, specifies the predicted value.
    See details below.}
  \item{partial}{Should partial plots be used?}
  \item{oob}{OOB (TRUE)  or in-bag (FALSE) predicted values.}
  \item{show.plots}{Should plots be displayed?}
  \item{plots.per.page}{Integer value controlling page layout.}
  \item{granule}{Integer value controlling whether a plot for a
    specific variable should be treated as a factor and therefore given
    as a boxplot. Larger values coerce boxplots.}
  \item{sorted}{Should variables be sorted by importance values.}
  \item{nvar}{Number of variables to be plotted. Default is all.}
  \item{npts}{Maximum number of points used when generating partial
    plots for continuous variables.}
  \item{smooth.lines}{Use lowess to smooth partial plots.}
  \item{subset}{Vector indicating which rows of the x-variable matrix
    \code{x$xvar} to use. All rows are used if not specified.}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  The vertical axis displays the ensemble predicted value, while
  x-variables are plotted on the horizontal axis.
  
  \enumerate{
    \item For regression, the predicted response is used.
    
    \item For classification, it is the predicted class probability
    specified by \option{which.class}, or the class of maximum
    probability depending on \option{class.type}.
    
    \item For multivariate families, it is the predicted value of the
    outcome specified by \option{outcome.target} and if that is a
    classification outcome, by \option{which.class}.
    
    \item For survival, the choices are:
    \itemize{
	    \item  Mortality (\code{mort}).
	    \item  Relative frequency of mortality (\code{rel.freq}).
	    \item Predicted survival (\code{surv}), where the predicted
      survival is for the time point specified using
      \code{time} (the default is the median follow up time).
    }
    \item For competing risks, the choices are:
    \itemize{
	    \item  The expected number of life years lost (\code{years.lost}).
	    \item  The cumulative incidence function (\code{cif}).
	    \item  The cumulative hazard function (\code{chf}).
    }
    In all three cases, the predicted value is for the event type
    specified by \option{which.class}.  For \code{cif} and
    \code{chf} the quantity is evaluated at the time point specified
    by \code{time}.
  }
  
  For partial plots use \option{partial=TRUE}.  Their interpretation are
  different than marginal plots.  The y-value for a variable \eqn{X},
  evaluated at \eqn{X=x}, is
  \deqn{
    \tilde{f}(x) = \frac{1}{n} \sum_{i=1}^n \hat{f}(x, x_{i,o}),
  }
  where \eqn{x_{i,o}} represents the value for all other variables
  other than \eqn{X} for individual \eqn{i} and \eqn{\hat{f}} is the
  predicted value. Generating partial plots can be very slow.
  Choosing a small value for \code{npts} can speed up computational
  times as this restricts the number of distinct \eqn{x} values used
  in computing \eqn{\tilde{f}}.
  
  For continuous variables, red points are used to indicate partial
  values and dashed red lines indicate a smoothed error bar of +/- two
  standard errors.  Black dashed line are the partial values.  Set
  \option{smooth.lines=TRUE} for lowess smoothed lines.  For discrete
  variables, partial values are indicated using boxplots with whiskers
  extending out approximately two standard errors from the mean.
  Standard errors are meant only to be a guide and should be
  interpreted with caution.
  
  Partial plots can be slow.  Setting \option{npts} to a smaller number
  can help.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Friedman J.H. (2001). Greedy function approximation: a gradient
  boosting machine, \emph{Ann. of Statist.}, 5:1189-1232.
  
  Ishwaran H., Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7(2):25-31.
  
  Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
  (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.
  
  Ishwaran H., Gerds T.A., Kogalur U.B., Moore R.D., Gange S.J. and Lau
  B.M. (2014). Random survival forests for competing risks.  To appear
  in \emph{Biostatistics}.
}
\seealso{
  \command{\link{rfsrc}}, \command{\link{rfsrcSyn}},
  \command{\link{predict.rfsrc}}
}
\examples{
\dontrun{
## ------------------------------------------------------------
## survival/competing risk
## ------------------------------------------------------------

## survival
data(veteran, package = "randomForestSRC")
v.obj <- rfsrc(Surv(time,status)~., veteran, nsplit = 10, ntree = 100)
plot.variable(v.obj, plots.per.page = 3)
plot.variable(v.obj, plots.per.page = 2, xvar.names = c("trt", "karno", "age"))
plot.variable(v.obj, surv.type = "surv", nvar = 1, time = 200)
plot.variable(v.obj, surv.type = "surv", partial = TRUE, smooth.lines = TRUE)
plot.variable(v.obj, surv.type = "rel.freq", partial = TRUE, nvar = 2)

## example of plot.variable calling a pre-processed plot.variable object
p.v <- plot.variable(v.obj, surv.type = "surv", partial = TRUE, smooth.lines = TRUE)
plot.variable(p.v)
p.v$plots.per.page <- 1
p.v$smooth.lines <- FALSE
plot.variable(p.v)

## competing risks
data(follic, package = "randomForestSRC")
follic.obj <- rfsrc(Surv(time, status) ~ ., follic, nsplit = 3, ntree = 100)
plot.variable(follic.obj, which.class = 2)

## ------------------------------------------------------------
## regression
## ------------------------------------------------------------

## airquality
airq.obj <- rfsrc(Ozone ~ ., data = airquality)
plot.variable(airq.obj, partial = TRUE, smooth.lines = TRUE)

## motor trend cars
mtcars.obj <- rfsrc(mpg ~ ., data = mtcars)
plot.variable(mtcars.obj, partial = TRUE, smooth.lines = TRUE)

## ------------------------------------------------------------
## classification
## ------------------------------------------------------------

## iris
iris.obj <- rfsrc(Species ~., data = iris)
plot.variable(iris.obj, partial = TRUE)

## motor trend cars: predict number of carburetors
mtcars2 <- mtcars
mtcars2$carb <- factor(mtcars2$carb,
   labels = paste("carb", sort(unique(mtcars$carb))))
mtcars2.obj <- rfsrc(carb ~ ., data = mtcars2)
plot.variable(mtcars2.obj, partial = TRUE)

## ------------------------------------------------------------
## multivariate regression
## ------------------------------------------------------------
mtcars.mreg <- rfsrc(Multivar(mpg, cyl) ~., data = mtcars)
plot.variable(mtcars.mreg, outcome.target = "mpg", partial = TRUE, nvar = 1)
plot.variable(mtcars.mreg, outcome.target = "cyl", partial = TRUE, nvar = 1)


## ------------------------------------------------------------
## multivariate mixed outcomes
## ------------------------------------------------------------
mtcars2 <- mtcars
mtcars2$carb <- factor(mtcars2$carb)
mtcars2$cyl <- factor(mtcars2$cyl)
mtcars.mix <- rfsrc(Multivar(carb, mpg, cyl) ~ ., data = mtcars2)
plot.variable(mtcars.mix, outcome.target = "cyl", which.class = "4", partial = TRUE, nvar = 1)
plot.variable(mtcars.mix, outcome.target = "cyl", which.class = 2, partial = TRUE, nvar = 1)


}
}
\keyword{plot}
