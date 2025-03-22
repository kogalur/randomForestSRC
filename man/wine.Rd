\name{wine}
\docType{data}
\alias{wine}
\title{White Wine Quality Data}
\description{
  The inputs include objective tests (e.g. PH values) and the output is
  based on sensory data (median of at least 3 evaluations made by wine
  experts) of white wine. Each expert graded the wine quality between 0
  (very bad) and 10 (very excellent).
}
\references{
  Cortez, P., Cerdeira, A., Almeida, F., Matos T. and Reis, J. (2009).
  Modeling wine preferences by data mining from physicochemical properties.
  In \emph{Decision Support Systems}, Elsevier, 47(4):547-553.
}
\examples{
## load wine and convert to a multiclass problem
data(wine, package = "randomForestSRC")
wine$quality <- factor(wine$quality)
}
\keyword{datasets}
