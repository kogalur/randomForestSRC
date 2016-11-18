\name{follic}
\docType{data}
\alias{follic}
\title{Follicular Cell Lymphoma}
\description{
  Competing risk data set involving follicular cell lymphoma. 
}
\format{A data frame containing:
  \tabular{ll}{
    age      \tab age\cr
	  hgb      \tab hemoglobin (g/l)\cr
    clinstg	 \tab clinical stage: 1=stage I, 2=stage II\cr
	  ch       \tab chemotherapy\cr
	  rt       \tab radiotherapy\cr
    time     \tab first failure time\cr
    status   \tab censoring status: 0=censored, 1=relapse, 2=death
  }
}
\source{Table 1.4b, \emph{Competing Risks: A Practical Perspective}.}
\references{
  Pintilie M., (2006) \emph{Competing Risks: A
    Practical Perspective.} West Sussex: John Wiley and Sons.
}
\examples{
\dontrun{
data(follic, package = "randomForestSRC")
follic.obj <- rfsrc(Surv(time, status) ~ ., follic, nsplit = 3, ntree = 100)
}
}
\keyword{datasets}
