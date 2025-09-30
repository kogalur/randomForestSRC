\name{plot.variable.rfsrc}
\alias{plot.variable.rfsrc}
\alias{plot.variable}
\title{Plot Marginal Effect of Variables}

\description{
  Plot the marginal effect of an x-variable on the class probability
  (classification), response (regression), mortality (survival), or the
  expected years lost (competing risk).  Users can select between
  marginal (unadjusted, but fast) and partial plots (adjusted, but
  slower).
}

\usage{\method{plot.variable}{rfsrc}(x, xvar.names, target,
  m.target = NULL, time, surv.type = c("mort", "rel.freq",
  "surv", "years.lost", "cif", "chf"), class.type =
  c("prob", "bayes"), partial = FALSE, oob = TRUE,
  show.plots = TRUE, plots.per.page = 4, granule = 5, sorted = TRUE,
  nvar, npts = 25, smooth.lines = FALSE, subset, ...)
}

\arguments{

  \item{x}{An object of class \code{(rfsrc, grow)}, \code{(rfsrc, synthetic)}, or \code{(rfsrc, plot.variable)}.}

  \item{xvar.names}{Character vector of x-variable names to include. If not specified, all variables are used.}

  \item{target}{For classification, an integer or character specifying the class of interest (default is the first class). For competing risks, an integer between 1 and \code{J} indicating the event of interest, where \code{J} is the number of event types. Default is the first event type.}

  \item{m.target}{Character value for multivariate families specifying the target outcome. If unspecified, a default is automatically chosen.}

  \item{time}{(Survival only) Time point at which the predicted survival value is evaluated, depending on \code{surv.type}.}

  \item{surv.type}{(Survival only) Type of predicted survival value to compute. See \code{\link{plot.variable}} details.}

  \item{class.type}{(Classification only) Type of predicted classification value to use. See \code{\link{plot.variable}} details.}

  \item{partial}{Logical. If \code{TRUE}, partial dependence plots are generated.}

  \item{oob}{Logical. If \code{TRUE}, out-of-bag predictions are used; otherwise, in-bag predictions are used.}

  \item{show.plots}{Logical. If \code{TRUE}, plots are displayed on the screen.}

  \item{plots.per.page}{Integer controlling the number of plots displayed per page.}

  \item{granule}{Integer controlling the coercion of continuous variables to factors (used to generate boxplots). Larger values increase coercion.}

  \item{sorted}{Logical. If \code{TRUE}, variables are sorted by variable importance.}

  \item{nvar}{Number of variables to plot. Defaults to all available variables.}

  \item{npts}{Maximum number of points used when generating partial plots for continuous variables.}

  \item{smooth.lines}{Logical. If \code{TRUE}, applies \code{lowess} smoothing to partial plots.}

  \item{subset}{Vector indicating which rows of \code{x$xvar} to use. Defaults to all rows. Important: do not define \code{subset} based on the original dataset (which may have been altered due to missing data or other processing); define it relative to \code{x$xvar}.}

  \item{...}{Additional arguments passed to or from other methods.}

}

\details{

The vertical axis displays the ensemble-predicted value, while x-variables are plotted along the horizontal axis.

\enumerate{

  \item For regression, the predicted response is plotted.

  \item For classification, the plotted value is the predicted class
  probability for the class specified by \code{target}, or the most
  probable class (Bayes rule) depending on whether \code{class.type} is
  set to \code{"prob"} or \code{"bayes"}.

  \item For multivariate families, the prediction corresponds to the outcome specified by \code{m.target}. If this is a classification outcome, \code{target} may also be used to indicate the class of interest.

  \item For survival, the vertical axis shows the predicted value determined by \code{surv.type}, with the following options:
  
  \itemize{
    \item \code{mort}: Mortality (Ishwaran et al., 2008), interpreted as the expected number of events for an individual with the same covariates.
    
    \item \code{rel.freq}: Relative frequency of mortality.

    \item \code{surv}: Predicted survival probability at a specified time point (default is the median follow-up time), controlled via \code{time}.
  }

  \item For competing risks, the vertical axis shows one of the following quantities, depending on \code{surv.type}:
  
  \itemize{
    \item \code{years.lost}: Expected number of life-years lost.

    \item \code{cif}: Cumulative incidence function for the specified event.

    \item \code{chf}: Cause-specific cumulative hazard function.
  }

  In all competing risks settings, the event of interest is specified using \code{target}, and \code{cif} and \code{chf} are evaluated at the time point given by \code{time}.
}

To generate partial dependence plots, set \code{partial = TRUE}. These differ from marginal plots in that they isolate the effect of a single variable \eqn{X} on the predicted value by averaging over all other covariates:

\deqn{
  \tilde{f}(x) = \frac{1}{n} \sum_{i=1}^n \hat{f}(x, x_{i,o}),
}

where \eqn{x_{i,o}} denotes the observed values of all covariates other than \eqn{X} for individual \eqn{i}, and \eqn{\hat{f}} is the prediction function. Generating partial plots can be computationally expensive; use a smaller value for \code{npts} to reduce the number of grid points evaluated for \eqn{x}.

\emph{Plot display conventions}:

\itemize{
  \item For continuous variables: red points indicate partial values; dashed red lines represent an error band of two standard errors. Black dashed lines show the raw partial values. Use \code{smooth.lines = TRUE} to overlay a \code{lowess} smoothed line.

  \item For discrete (factor) variables: boxplots are used, with whiskers extending approximately two standard errors from the mean.

  \item Standard errors are provided only as rough indicators and should be interpreted cautiously.
}

Partial plots can be slow to compute. Setting \code{npts} to a small value can improve performance.

For additional flexibility and speed, consider using \code{\link{partial.rfsrc}}, which directly computes partial plot data and allows for greater customization.
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
  B.M. (2014). Random survival forests for competing risks.
  \emph{Biostatistics}, 15(4):757-773.
  
}
\seealso{
  \command{\link{rfsrc}},
  \command{\link{partial.rfsrc}},
  \command{\link{predict.rfsrc}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## survival/competing risk
## ------------------------------------------------------------

## survival
data(veteran, package = "randomForestSRC")
v.obj <- rfsrc(Surv(time,status)~., veteran, ntree = 100)
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

## example using a pre-processed plot.variable to define custom plots
p.v <- plot.variable(v.obj, surv.type = "surv", partial = TRUE, show.plots = FALSE)
plotthis <- p.v$plotthis
plot(plotthis[["age"]], xlab = "age", ylab = "partial effect", type = "b")
boxplot(yhat ~ x, plotthis[["trt"]], xlab = "treatment", ylab = "partial effect")


## competing risks
data(follic, package = "randomForestSRC")
follic.obj <- rfsrc(Surv(time, status) ~ ., follic, nsplit = 3, ntree = 100)
plot.variable(follic.obj, target = 2)

## ------------------------------------------------------------
## regression
## ------------------------------------------------------------

## airquality 
airq.obj <- rfsrc(Ozone ~ ., data = airquality)
plot.variable(airq.obj, partial = TRUE, smooth.lines = TRUE)
plot.variable(airq.obj, partial = TRUE, subset = airq.obj$xvar$Solar.R < 200)

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
plot.variable(mtcars.mreg, m.target = "mpg", partial = TRUE, nvar = 1)
plot.variable(mtcars.mreg, m.target = "cyl", partial = TRUE, nvar = 1)


## ------------------------------------------------------------
## multivariate mixed outcomes
## ------------------------------------------------------------
mtcars2 <- mtcars
mtcars2$carb <- factor(mtcars2$carb)
mtcars2$cyl <- factor(mtcars2$cyl)
mtcars.mix <- rfsrc(Multivar(carb, mpg, cyl) ~ ., data = mtcars2)
plot.variable(mtcars.mix, m.target = "cyl", target = "4", partial = TRUE, nvar = 1)
plot.variable(mtcars.mix, m.target = "cyl", target = 2, partial = TRUE, nvar = 1)


}
}
\keyword{plot}
