\name{partial.rfsrc}
\alias{partial.rfsrc}
\alias{partial}
\alias{get.partial.plot.data}
\title{Acquire Partial Effect of a Variable}
\description{
  Direct, fast inferface for partial effect of a variable.  Works for all families.
}
\usage{
partial.rfsrc(object, oob = TRUE, 
  partial.type = NULL, partial.xvar = NULL, partial.values = NULL,
  partial.xvar2 = NULL, partial.values2 = NULL,
  partial.time = NULL, get.tree = NULL, seed = NULL, do.trace = FALSE, ...)
}
\arguments{
  \item{object}{An object of class \code{(rfsrc, grow)}.}
  \item{oob}{By default out-of-bag values are returned, but inbag
    values can be requested by setting this option to \code{FALSE}.}
  \item{partial.type}{Character vector specifying type of predicted
    value requested.  See details below.}
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
  \item{get.tree}{Vector of integer(s) identifying trees over which the
    partial values are calculated over.  By default, uses all trees in the
    forest.}
  \item{seed}{Negative integer specifying seed for the random number
    generator.}
  \item{do.trace}{Number of seconds between updates to the user on
    approximate time to completion.}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{

  Used for direct, efficient call to obtain partial plot effects.
  This function is intended primarily for experts. 
  
  Out-of-bag (OOB) values are returned by default.

  For factors, the partial value should be encoded as a positive integer
  reflecting the level number of the factor.  The actual label of the
  factor should not be used.

  The utility function \code{get.partial.plot.data} is supplied for
  processing returned raw partial effects in a format more convenient
  for plotting.  Options are specified as in \code{plot.variable}.  See
  examples for illustration.
  
  Raw partial plot effects data is returned either as an array or a list
  of length equal to the number of outcomes (length is one for
  univariate families) with entries depending on the underlying family:
  
  \enumerate{
    \item For regression, partial plot data is returned as a list
    in \code{regrOutput} with dim \code{[n] x  [length(partial.values)]}.
    
    \item For classification, partial plot data is returned as a list
    in \code{classOutput} of dim \code{[n] x [1 + yvar.nlevels[.]] x
      [length(partial.values)]}.
    
    \item For mixed multivariate regression, values are returned in 
    list format both in \code{regrOutput} and \code{classOutput}
    
    \item For survival, values are returned as either a matrix or array
    in \code{survOutput}. Depending on partial type specified this can be:
    \itemize{
      \item For partial type \code{surv} returns the survival function
         of dim \code{[n] x [length(partial.time)] x [length(partial.values)]}.
      \item For partial type \code{mort} returns mortality
        of dim \code{[n] x [length(partial.values)]}.
      \item For partial type \code{chf} returns the cumulative hazard function
        of dim \code{[n] x [length(partial.time)] x [length(partial.values)]}.
    }
    
    \item For competing risks, values are returned as either a matrix or
    array in \code{survOutput}.  Depending on the options specified this
    can be:
    \itemize{
      \item For partial type \code{years.lost} returns the expected number of life years lost 
      of dim \code{[n] x [length(event.info$event.type)] x
        [length(partial.values)]}.
      \item For partial type \code{cif} returns the cumulative incidence function of dim
	    \code{[n] x [length(partial.time)] x
	      [length(event.info$event.type)] x
	      [length(partial.values)]}.
      \item For partial type \code{chf} returns the cumulative hazard function of dim
            \code{[n] x [length(partial.time)] x [length(event.info$event.type)]
            x [length(partial.values)]}.
    }

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
}
\seealso{
  \command{\link{plot.variable.rfsrc}}
}
\examples{
\donttest{

## ------------------------------------------------------------
##
## regression
##
## ------------------------------------------------------------

airq.obj <- rfsrc(Ozone ~ ., data = airquality)

## partial effect for wind
partial.obj <- partial(airq.obj,
                  partial.xvar = "Wind",
                  partial.values = airq.obj$xvar$Wind)
pdta <- get.partial.plot.data(partial.obj)

## plot partial values
plot(pdta$x, pdta$yhat, type = "b", pch = 16,
      xlab = "wind", ylab = "partial effect of wind")


## example where we display all the partial effects
## instead of averaging - use the granule=TRUE option
pdta <- get.partial.plot.data(partial.obj, granule = TRUE)
boxplot(pdta$yhat ~ pdta$x, xlab = "Wind", ylab = "partial effect")

## ------------------------------------------------------------
##
## regression: partial effects for two variables simultaneously
##
## ------------------------------------------------------------

airq.obj <- rfsrc(Ozone ~ ., data = airquality)

## specify wind and temperature values of interest
wind <- sort(unique(airq.obj$xvar$Wind))
temp <- sort(unique(airq.obj$xvar$Temp))

## partial effect for wind, for a given temp
pdta <- do.call(rbind, lapply(temp, function(x2) {
  o <- partial(airq.obj,
         partial.xvar = "Wind", partial.xvar2 = "Temp",
         partial.values = wind, partial.values2 = x2)
  cbind(wind, x2, get.partial.plot.data(o)$yhat)
}))
pdta <- data.frame(pdta)
colnames(pdta) <- c("wind", "temp", "effectSize")

## coplot of partial effect of wind and temp 
coplot(effectSize ~ wind|temp, pdta, pch = 16, overlap = 0)


## ------------------------------------------------------------
##
## regression: partial effects for three variables simultaneously
## (can be slow, so modify accordingly)
##
## ------------------------------------------------------------

n <- 1000
x <- matrix(rnorm(n * 3), ncol = 3)
y <- x[, 1] + x[, 1] * x[, 2] + x[, 1] * x[, 2] * x[, 3]
o <- rfsrc(y ~ ., data = data.frame(y = y, x))

## define target x values
x1 <- seq(-3, 3, length = 40)
x2 <- x3 <- seq(-3, 3, length = 10)

## extract second order partial effects
pdta <- do.call(rbind,
          lapply(x3, function(x3v) {
            cat("outer loop x3 = ", x3v, "\n")
            do.call(rbind,lapply(x2, function(x2v) {
              o <- partial(o,
                      partial.xvar = "X1",
                      partial.values = x1,
                      partial.xvar2 = c("X2", "X3"),
                      partial.values2 = c(x2v, x3v))
              cbind(x1, x2v, x3v, get.partial.plot.data(o)$yhat)
            }))
          }))
pdta <- data.frame(pdta)
colnames(pdta) <- c("x1", "x2", "x3", "effectSize")

## coplot of partial effects
coplot(effectSize ~ x1|x2*x3, pdta, pch = 16, overlap = 0)


## ------------------------------------------------------------
##
## classification
##
## ------------------------------------------------------------

iris.obj <- rfsrc(Species ~., data = iris)

## partial effect for sepal length
partial.obj <- partial(iris.obj,
                  partial.xvar = "Sepal.Length",
                  partial.values = iris.obj$xvar$Sepal.Length)

## extract partial effects for each species outcome
pdta1 <- get.partial.plot.data(partial.obj, target = "setosa")
pdta2 <- get.partial.plot.data(partial.obj, target = "versicolor")
pdta3 <- get.partial.plot.data(partial.obj, target = "virginica")

## plot the results
par(mfrow=c(1,1))
plot(pdta1$x, pdta1$yhat, type="b", pch = 16,
     xlab = "sepal length", ylab = "adjusted probability", 
     ylim = range(pdta1$yhat,pdta2$yhat,pdta3$yhat))
points(pdta2$x, pdta2$yhat, col = 2, type = "b", pch = 16)
points(pdta3$x, pdta3$yhat, col = 4, type = "b", pch = 16)
legend("topleft", legend=levels(iris.obj$yvar), fill = c(1, 2, 4))

## ------------------------------------------------------------
##
## survival
##
## ------------------------------------------------------------

data(veteran, package = "randomForestSRC")
v.obj <- rfsrc(Surv(time,status)~., veteran, nsplit = 10, ntree = 100)

## partial effect of age on mortality
partial.obj <- partial(v.obj,
  partial.type = "mort",
  partial.xvar = "age",
  partial.values = v.obj$xvar$age,
  partial.time = v.obj$time.interest)
pdta <- get.partial.plot.data(partial.obj)

plot(lowess(pdta$x, pdta$yhat, f = 1/3),
   type = "l", xlab = "age", ylab = "adjusted mortality")

## example where x is discrete - partial effect of age on mortality
## we use the granule=TRUE option
partial.obj <- partial(v.obj,
       partial.type = "mort",
       partial.xvar = "trt",
       partial.values = v.obj$xvar$trt,
       partial.time = v.obj$time.interest)
pdta <- get.partial.plot.data(partial.obj, granule = TRUE)
boxplot(pdta$yhat ~ pdta$x, xlab = "treatment", ylab = "partial effect")


## partial effects of karnofsky score on survival
karno <- quantile(v.obj$xvar$karno)
partial.obj <- partial(v.obj,
  partial.type = "surv",
  partial.xvar = "karno",
  partial.values = karno,
  partial.time = v.obj$time.interest)
pdta <- get.partial.plot.data(partial.obj)

matplot(pdta$partial.time, t(pdta$yhat), type = "l", lty = 1,
     xlab = "time", ylab = "karnofsky adjusted survival")
legend("topright", legend = paste0("karnofsky = ", karno), fill = 1:5)


## ------------------------------------------------------------
##
## competing risk
##
## ------------------------------------------------------------

data(follic, package = "randomForestSRC")
follic.obj <- rfsrc(Surv(time, status) ~ ., follic, nsplit = 3, ntree = 100)

## partial effect of age on years lost
partial.obj <- partial(follic.obj,
  partial.type = "years.lost",
  partial.xvar = "age",
  partial.values = follic.obj$xvar$age,
  partial.time = follic.obj$time.interest)
pdta1 <- get.partial.plot.data(partial.obj, target = 1)
pdta2 <- get.partial.plot.data(partial.obj, target = 2)

par(mfrow=c(2,2))
plot(lowess(pdta1$x, pdta1$yhat),
   type = "l", xlab = "age", ylab = "adjusted years lost relapse")
plot(lowess(pdta2$x, pdta2$yhat),
   type = "l", xlab = "age", ylab = "adjusted years lost death")

## partial effect of age on cif
partial.obj <- partial(follic.obj,
  partial.type = "cif",
  partial.xvar = "age",
  partial.values = quantile(follic.obj$xvar$age),
  partial.time = follic.obj$time.interest)
pdta1 <- get.partial.plot.data(partial.obj, target = 1)
pdta2 <- get.partial.plot.data(partial.obj, target = 2)

matplot(pdta1$partial.time, t(pdta1$yhat), type = "l", lty = 1,
     xlab = "time", ylab = "age adjusted cif for relapse")
matplot(pdta2$partial.time, t(pdta2$yhat), type = "l", lty = 1,
     xlab = "time", ylab = "age adjusted cif for death")


## ------------------------------------------------------------
##
## multivariate mixed outcomes
##
## ------------------------------------------------------------

mtcars2 <- mtcars
mtcars2$carb <- factor(mtcars2$carb)
mtcars2$cyl <- factor(mtcars2$cyl)
mtcars.mix <- rfsrc(Multivar(carb, mpg, cyl) ~ ., data = mtcars2)

## partial effect of displacement for each the three-outcomes
partial.obj <- partial(mtcars.mix,
                  partial.xvar = "disp",
                  partial.values = mtcars.mix$xvar$disp)
pdta1 <- get.partial.plot.data(partial.obj, m.target = "carb")
pdta2 <- get.partial.plot.data(partial.obj, m.target = "mpg")
pdta3 <- get.partial.plot.data(partial.obj, m.target = "cyl")

par(mfrow=c(2,2))
plot(lowess(pdta1$x, pdta1$yhat), type = "l", xlab="displacement", ylab="carb")
plot(lowess(pdta2$x, pdta2$yhat), type = "l", xlab="displacement", ylab="mpg")
plot(lowess(pdta3$x, pdta3$yhat), type = "l", xlab="displacement", ylab="cyl")


}
}
\keyword{partial}
