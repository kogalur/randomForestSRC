\name{hd}
\docType{data}
\alias{hd}
\title{Hodgkin's Disease}
\description{
  Competing risk data set involving Hodgkin's disease.
}
\format{A data frame containing:
  \tabular{ll}{
    age      \tab age\cr
    sex      \tab gender\cr
	  trtgiven \tab treatment: RT=radition, CMT=Chemotherapy and radiation\cr
	  medwidsi \tab mediastinum involvement: N=no, S=small, L=Large\cr
	  extranod \tab extranodal disease: Y=extranodal disease, N=nodal disease\cr
    clinstg	 \tab clinical stage: 1=stage I, 2=stage II\cr
    time     \tab first failure time\cr
    status   \tab censoring status: 0=censored, 1=relapse, 2=death
  }
}
\source{
  Table 1.6b, \emph{Competing Risks: A Practical Perspective}.
}
\references{
  Pintilie M., (2006) \emph{Competing Risks: A
    Practical Perspective.} West Sussex: John Wiley and Sons.
}
\examples{data(hd, package = "randomForestSRC")}
\keyword{datasets}
