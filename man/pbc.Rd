\name{pbc}
\docType{data}
\alias{pbc}
\title{Primary Biliary Cirrhosis (PBC) Data}
\description{
  Data from the Mayo Clinic trial in primary biliary cirrhosis (PBC) of
  the liver conducted between 1974 and 1984.  A total of 424 PBC
  patients, referred to Mayo Clinic during that ten-year interval, met
  eligibility criteria for the randomized placebo controlled trial of
  the drug D-penicillamine.  The first 312 cases in the data set
  participated in the randomized trial and contain largely complete
  data.  
}
\source{Flemming and Harrington, 1991, Appendix D.1.}
\references{
  Flemming T.R and Harrington D.P., (1991) \emph{Counting Processes
    and Survival Analysis.} New York: Wiley.
}
\examples{
\dontrun{
data(pbc, package = "randomForestSRC")
pbc.obj <- rfsrc(Surv(days, status) ~ ., pbc, nsplit = 3)
}
}
\keyword{datasets}
