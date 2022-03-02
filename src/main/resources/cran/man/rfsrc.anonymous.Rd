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
  present, the levels of the factors.  The topology of grow trees is
  also saved, which includes among other things, the split values used
  for splitting tree nodes.

  For the most privacy, we recommend that variable names be made
  non-identifiable and that data be coerced to real values.  If factors
  are required, the user should consider using non-identifiable factor
  levels.  However, in all cases, it is the users responsibility to
  de-identify their data and to check that data privacy holds.  We
  provide no guarantees of this.

  While anonymous random forests works similar to random forests, there
  are caveats to keep in mind.  First, no missing data is allowed since
  missing data imputation requires training data.  Second, while
  anonymous forest tries to play nice with the functions in the package,
  it only works with functions that specifically do not require training
  data.  Thus users are advised to keep this in mind if they decide to
  go this route.
  
    
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
## ------------------------------------------------------------
## illustrates minimal information saved by the forest
## TBD TBD currently not yet implemented TBD TBD
## plot(get.tree(rfsrc.anonymous(mpg ~ ., mtcars), 10))

## ------------------------------------------------------------
## classification
## ------------------------------------------------------------
print(rfsrc.anonymous(Species ~ ., iris))

## ------------------------------------------------------------
## survival
## ------------------------------------------------------------
data(veteran, package = "_PROJECT_PACKAGE_NAME_")
print(rfsrc.anonymous(Surv(time, status) ~ ., data = veteran))

## ------------------------------------------------------------
## competing risk
## ------------------------------------------------------------
data(wihs, package = "_PROJECT_PACKAGE_NAME_")
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
## train/test setting but tricky because factor labels differ over
## training and test data
## ------------------------------------------------------------

# first we convert all x-variables to factors
data(veteran, package = "_PROJECT_PACKAGE_NAME_")
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