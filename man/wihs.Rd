\name{wihs}
\docType{data}
\alias{wihs}
\title{Women's Interagency HIV Study (WIHS)}
\description{
  Competing risk data set involving AIDS in women.
}
\format{A data frame containing:
  \tabular{ll}{
    time \tab time to event\cr
    status \tab censoring status: 0=censoring, 1=HAART initiation, 2=AIDS/Death before HAART\cr
    ageatfda \tab age in years at time of FDA approval of first protease inhibitor\cr
    idu \tab history of IDU: 0=no history, 1=history\cr
    black \tab race:  0=not African-American; 1=African-American\cr
    cd4nadir \tab CD4 count (per 100 cells/ul)
  }
}
\source{
  Study included 1164 women enrolled in WIHS, who were alive, infected
  with HIV, and free of clinical AIDS on December, 1995, when the first
  protease inhibitor (saquinavir mesylate) was approved by the Federal
  Drug Administration. Women were followed until the first of the
  following occurred: treatment initiation, AIDS diagnosis, death, or
  administrative censoring (September, 2006). Variables included history
  of injection drug use at WIHS enrollment, whether an individual was
  African American, age, and CD4 nadir prior to baseline.}
\references{
  Bacon M.C, von Wyl V., Alden C., et al. (2005). The Women's Interagency
  HIV Study: an observational cohort brings clinical sciences to the
  bench, \emph{Clin Diagn Lab Immunol}, 12(9):1013-1019.
}
\examples{
\dontrun{
data(wihs, package = "randomForestSRC")
wihs.obj <- rfsrc(Surv(time, status) ~ ., wihs, nsplit = 3, ntree = 100)
}
}
\keyword{datasets}
