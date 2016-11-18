\name{nutrigenomic}
\docType{data}
\alias{nutrigenomic}
\title{Nutrigenomic Study}
\description{
  Study the effects of five diet treatments on 21 liver lipids and 120
  hepatic gene expression in wild-type and PPAR-alpha deficient mice.  We
  use a multivariate mixed random forest analysis by regressing gene
  expression, diet and genotype (the x-variables) on lipid expressions
  (the multivariate y-responses).
}
\references{
  Martin P.G. et al. (2007). Novel aspects of PPAR-alpha-mediated
  regulation of lipid and xenobiotic metabolism revealed through a
  nutrigenomic study. \emph{Hepatology}, 45(3), 767--777.
}
\examples{
\dontrun{
## ------------------------------------------------------------
## multivariate mixed forests
## lipids used as the multivariate y-responses
## ------------------------------------------------------------

## load the data
data(nutrigenomic, package = "randomForestSRC")

## nice wrapper for making multivariate formula
mvrfsrc.f <- function(ynames, dat) {
  as.formula(paste("Multivar(", paste(ynames, collapse = ","),paste(") ~ ."), sep = ""))
}

## multivariate mixed forest call
mv.obj <- rfsrc(mvrfsrc.f(colnames(nutrigenomic$lipids)),
            data.frame(do.call(cbind, nutrigenomic)),
            importance=TRUE, nsplit = 10)


## ------------------------------------------------------------
## nice wrappers to extract standardized performance values
## works for all forests - including multivariate forests
## ------------------------------------------------------------

## pull the standardized error from a forest object
get.error <- function(obj) {
  100 * c(sapply(obj$yvar.names, function(nn) {
    o.coerce <- randomForestSRC:::coerce.multivariate(obj, nn)
    if (o.coerce$family == "class") {
      tail(o.coerce$err.rate[, 1], 1)
    }
    else {
      tail(o.coerce$err.rate, 1) / var(o.coerce$yvar, na.rm = TRUE)
    }
  }))
}

## pull the standardized VIMP from a forest object
get.vimp <- function(obj) {
  vimp <- 100 * do.call(cbind, lapply(obj$yvar.names, function(nn) {
    o.coerce <- randomForestSRC:::coerce.multivariate(obj, nn)
    if (o.coerce$family == "class") {
      o.coerce$importance[, 1]
    }
    else {
      o.coerce$importance / var(o.coerce$yvar, na.rm = TRUE)
    }
  }))
  colnames(vimp) <- obj$yvar.names
  vimp
}

## ------------------------------------------------------------
## plot the standarized performance and VIMP values
## ------------------------------------------------------------

## acquire the error rate for each of the 21-coordinates 
## standardize to allow for comparison across coordinates
serr <- get.error(mv.obj)

## acquire standardized VIMP 
svimp <- get.vimp(mv.obj)

par(mfrow = c(1,2))
plot(serr, xlab = "Lipids", ylab = "Standardized Performance")
matplot(svimp, xlab = "Genes/Diet/Genotype", ylab = "Standardized VIMP")


}}
\keyword{datasets}
