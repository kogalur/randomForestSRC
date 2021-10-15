\name{housing}
\docType{data}
\alias{housing}
\title{Ames Iowa Housing Data}

\description{
  Data from the Ames Assessor's Office used in assessing values of
  individual residential properties sold in Ames, Iowa from 2006 to
  2010.  This is a regression problem and the goal is to predict 
  "SalePrice" which records the price of a home in thousands of dollars.
}
\references{
  De Cock, D., (2011). Ames, Iowa: Alternative to the Boston housing
  data as an end of semester regression project. \emph{Journal of Statistics
  Education}, 19(3), 1--14. 
}
\examples{
\donttest{
## load the data
data(housing, package = "_PROJECT_PACKAGE_NAME_")

## the original data contains lots of missing data
## here's a fast but reasonably accurate way to impute the data
housing2 <- impute(data = housing, mf.q = .5, fast = TRUE)

## even faster, but potentially less acurate
housing3 <- impute(SalePrice~., housing, splitrule = "random", nimpute = 1)


}}
\keyword{datasets}
