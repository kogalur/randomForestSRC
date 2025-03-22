\name{peakVO2}
\docType{data}
\alias{peakVO2}
\title{Systolic Heart Failure Data}

\description{

  The data involve 2231 patients with systolic heart failure who
  underwent cardiopulmonary stress testing at the Cleveland Clinic.  The
  primary end point was all-cause death.  In total, 39 variables were
  measured for each patient, including baseline clinical values and
  exercise stress test results.  A key variable of interest is
  peak VO2 (mL/kg per min), the peak respiratory exchange ratio.
  More details regarding the data can be found in Hsich et al. (2011).
  
}
\references{
  
  Hsich E., Gorodeski E.Z.,Blackstone E.H., Ishwaran H. and Lauer
  M.S. (2011). Identifying important risk factors for survival in
  systolic heart failure patients using random survival
  forests. Circulation: Cardio. Qual. Outcomes, 4(1), 39-45.
  
}
\examples{
\donttest{
## load the data
data(peakVO2, package = "randomForestSRC")

## random survival forest analysis
o <- rfsrc(Surv(ttodead, died)~., peakVO2)
print(o)

## partial effect of peak V02 on mortality
partial.o <- partial(o,
       partial.type = "mort",
       partial.xvar = "peak.vo2",
       partial.values = o$xvar$peak.vo2,
       partial.time = o$time.interest)
pdta.m <- get.partial.plot.data(partial.o)


## partial effect of peak V02 on survival
pvo2 <- quantile(o$xvar$peak.vo2)
partial.o <- partial(o,
       partial.type = "surv",
       partial.xvar = "peak.vo2",
       partial.values = pvo2,
       partial.time = o$time.interest)
pdta.s <- get.partial.plot.data(partial.o)
     

## compare the two plots
par(mfrow=c(1,2))    

plot(lowess(pdta.m$x, pdta.m$yhat, f = 2/3),
     type = "l", xlab = "peak VO2", ylab = "adjusted mortality")
rug(o$xvar$peak.vo2)

matplot(pdta.s$partial.time, t(pdta.s$yhat), type = "l", lty = 1,
          xlab = "years", ylab = "peak VO2 adjusted survival")
legend("bottomleft", legend = paste0("peak VO2 = ", pvo2),
       bty = "n", cex = .75, fill = 1:5)


}}
\keyword{datasets}
