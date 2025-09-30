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
  \url{http://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+(Prognostic)}.
}
\examples{
\donttest{
## ------------------------------------------------------------
## Standard analysis
## ------------------------------------------------------------

data(breast, package = "randomForestSRC")
breast <- na.omit(breast)
o <- rfsrc(status ~ ., data = breast, nsplit = 10)
print(o)

}}
\keyword{datasets}
