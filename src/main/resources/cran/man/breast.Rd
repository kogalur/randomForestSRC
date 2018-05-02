\name{breast}
\docType{data}
\alias{breast}
\title{Wisconsin Prognostic Breast Cancer Data}
\description{
  Recurrence of breast cancer from 198 breast cancer patients, all of
  which exhibited no evidence of distant metastases at the time of
  diagnosis.  The first 30 features of the data describe characteristics
  of the cell nuclei present in the digitized image of a fine needle
  aspirate (FNA) of the breast mass.
}
\source{
  The data were obtained from the UCI machine learning repository, see
  \url{_HTTP_PREFIX_archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+(Prognostic)}.
}
 
\examples{
\donttest{
## ------------------------------------------------------------
## Standard analysis
## ------------------------------------------------------------

data(breast, package = "_PROJECT_PACKAGE_NAME_")
breast <- na.omit(breast)
o <- rfsrc(status ~ ., data = breast, nsplit = 10)
print(o)

## ------------------------------------------------------------
## The data is imbalanced so we use balanced random forests
## with undersampling of the majority class
##
## Specifically let n0, n1 be sample sizes for majority, minority
## cases.  We sample 2 x n1 cases with majority, minority cases chosen
## with probabilities n1/n, n0/n where n=n0+n1
## ------------------------------------------------------------

y <- breast$status
o <- rfsrc(status ~ ., data = breast, nsplit = 10,
            case.wt = _PROJECT_PACKAGE_NAME_:::make.wt(y),
            sampsize = _PROJECT_PACKAGE_NAME_:::make.size(y))
print(o)

}}
\keyword{datasets}
