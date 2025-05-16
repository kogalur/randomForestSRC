\name{rfsrc.anonymous}
\alias{rfsrc.anonymous}
\title{Anonymous Random Forests}

\description{

  Anonymous random forests is carefully modified to ensure that the
  original training data is not retained. This enables users to share
  the trained forest with others without disclosing the underlying data.

}

\usage{rfsrc.anonymous(formula, data, forest = TRUE, ...)}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. If missing, unsupervised splitting is performed.}

  \item{data}{A data frame containing the y-outcome and x-variables.}

  \item{forest}{Logical. Should the forest object be returned? Required for prediction on new data and by many other package functions.}

  \item{...}{Additional arguments passed to \code{\link{rfsrc}}. See the \code{rfsrc} help file for full details.}
}


\details{

This function calls \code{\link{rfsrc}} and returns a forest object with the original training data removed. This enables users to share their forest while preserving the privacy of their data.

To enable prediction on new (test) data, certain minimal information from the training data must still be retained. This includes:
\itemize{
  \item Names of the original variables.
  \item For factor variables, the levels of each factor.
  \item Summary statistics used for imputation: the mean for continuous variables and the most frequent class for factors.
  \item Tree topology, including split points used to grow the trees.
}

For maximal privacy, users are strongly encouraged to replace variable names with non-identifiable labels and convert all variables to continuous format when possible. If factor variables are used, their levels should also be anonymized. However, the user is solely responsible for de-identifying the data and verifying that privacy is maintained. \strong{We provide NO GUARANTEES regarding data confidentiality.}

\emph{Missing data handling:} Anonymous forests do not support imputation of training data. The option \code{na.action = "na.impute"} is automatically downgraded to \code{"na.omit"}. If training data contain missing values, we recommend pre-imputing them using \code{\link{impute}}.

Test data, however, \emph{can} be imputed at prediction time:
\itemize{
  \item \code{na.action = "na.impute"} performs a fast imputation by replacing missing values with the training mean (for numeric variables) or most frequent class (for factors).
  \item \code{na.action = "na.random"} uses a fast random draw from training distributions for imputation.
}

Although anonymous forests are compatible with many package functions, they are only guaranteed to work with functions that do not explicitly require access to the original training data.

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
