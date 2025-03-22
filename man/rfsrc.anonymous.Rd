\name{rfsrc.anonymous}
\alias{rfsrc.anonymous}
\title{Anonymous Random Forests}
\description{

  Anonymous random forests applies random forests but is carefully
  modified so as not to save the original training data.  This allows
  users to share their forest with other researchers but without having
  to share their original data.

}
\usage{rfsrc.anonymous(formula, data, forest = TRUE, ...)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit. If
    missing, unsupervised splitting is implemented.}
  \item{data}{Data frame containing the y-outcome and x-variables.} 
  \item{forest}{Should the forest object be returned?  Used for
    prediction on new data and required by many of the package
    functions.}
  \item{...}{Further arguments as in \code{\link{rfsrc}}.  See the
   \code{rfsrc} help file for details.}
}
\details{

  Calls \code{\link{rfsrc}} and returns an object with the training data
  removed so that users can share their forest while maintaining privacy
  of their data.

  In order to predict on test data, it is however necessary for certain
  minimal information to be saved from the training data.  This includes
  the names of the original variables, and if factor variables are
  present, the levels of the factors.  The mean value and maximal class
  value for real and factor variables in the training data are also
  stored for the purposes of imputation on test data (see below).  The
  topology of grow trees is also saved, which includes among other
  things, the split values used for splitting tree nodes.

  For the most privacy, we recommend that variable names be made
  non-identifiable and that data be coerced to real values.  If factors
  are required, the user should consider using non-identifiable factor
  levels.  However, in all cases, it is the users responsibility to
  de-identify their data and to check that data privacy holds.  We
  provide NO GUARANTEES of this.

  Missing data is especially delicate with anonymous forests.  Training
  data cannot be imputed and the option \code{na.action="na.impute"}
  simply reverts to \code{na.action="na.omit"}.  Therefore if you have
  training data with missing values consider using pre-imputing the data
  using \command{impute}.  It is however possible to impute on test
  data.  The option \code{na.action="na.impute"} in the prediction call
  triggers a rough and fast imputation method where the value of missing
  test data are replaced by the mean (or maximal class) value from the
  training data.  A second option \code{na.action="na.random"} uses a
  fast random imputation method.
  
  In general, it is important to keep in mind that while anonymous
  forests tries to play nice with other functions in the package, it
  only works with calls that do not specifically require training data.
  
    
}
\value{
  An object of class \code{(rfsrc, grow, anonymous)}.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\seealso{
  \command{\link{rfsrc}}
}
\examples{
\donttest{

## ------------------------------------------------------------
## regression
## ------------------------------------------------------------
print(rfsrc.anonymous(mpg ~ ., mtcars))

## ------------------------------------------------------------
## plot anonymous regression tree (using get.tree)
## TBD CURRENTLY NOT IMPLEMENTED 
## ------------------------------------------------------------
## plot(get.tree(rfsrc.anonymous(mpg ~ ., mtcars), 10))

## ------------------------------------------------------------
## classification
## ------------------------------------------------------------
print(rfsrc.anonymous(Species ~ ., iris))

## ------------------------------------------------------------
## survival
## ------------------------------------------------------------
data(veteran, package = "randomForestSRC")
print(rfsrc.anonymous(Surv(time, status) ~ ., data = veteran))

## ------------------------------------------------------------
## competing risks
## ------------------------------------------------------------
data(wihs, package = "randomForestSRC")
print(rfsrc.anonymous(Surv(time, status) ~ ., wihs, ntree = 100))

## ------------------------------------------------------------
## unsupervised forests
## ------------------------------------------------------------
print(rfsrc.anonymous(data = iris))

## ------------------------------------------------------------
## multivariate regression
## ------------------------------------------------------------
print(rfsrc.anonymous(Multivar(mpg, cyl) ~., data = mtcars))

## ------------------------------------------------------------
## prediction on test data with missing values using pbc data
## cases 1 to 312 have no missing values
## cases 313 to 418 having missing values
## ------------------------------------------------------------
data(pbc, package = "randomForestSRC")
pbc.obj <- rfsrc.anonymous(Surv(days, status) ~ ., pbc)
print(pbc.obj)

## mean value imputation
print(predict(pbc.obj, pbc[-(1:312),], na.action = "na.impute"))

## random imputation
print(predict(pbc.obj, pbc[-(1:312),], na.action = "na.random"))

## ------------------------------------------------------------
## train/test setting but tricky because factor labels differ over
## training and test data
## ------------------------------------------------------------

# first we convert all x-variables to factors
data(veteran, package = "randomForestSRC")
veteran.factor <- data.frame(lapply(veteran, factor))
veteran.factor$time <- veteran$time
veteran.factor$status <- veteran$status

# split the data into train/test data (25/75)
# the train/test data have the same levels, but different labels
train <- sample(1:nrow(veteran), round(nrow(veteran) * .5))
summary(veteran.factor[train, ])
summary(veteran.factor[-train, ])

# grow the forest on the training data and predict on the test data
v.grow <- rfsrc.anonymous(Surv(time, status) ~ ., veteran.factor[train, ]) 
v.pred <- predict(v.grow, veteran.factor[-train, ])
print(v.grow)
print(v.pred)



}}
\keyword{forest}
\keyword{anonymous}