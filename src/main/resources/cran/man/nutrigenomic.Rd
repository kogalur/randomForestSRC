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
data(nutrigenomic, package = "_PROJECT_PACKAGE_NAME_")

## nice wrapper for making multivariate formula
get.mv.f <- function(ynames) {
  as.formula(paste("Multivar(", paste(ynames, collapse = ","),paste(") ~ ."), sep = ""))
}

## multivariate mixed forest call
mv.obj <- rfsrc(get.mv.f(colnames(nutrigenomic$lipids)),
            data.frame(do.call(cbind, nutrigenomic)),
            importance=TRUE, nsplit = 10)

## ------------------------------------------------------------
## plot the standarized performance and VIMP values
## ------------------------------------------------------------

## acquire the error rate for each of the 21-coordinates 
## standardize to allow for comparison across coordinates
serr <- get.mv.error(mv.obj, std = TRUE)

## acquire standardized VIMP 
svimp <- get.mv.vimp(mv.obj, std = TRUE)

par(mfrow = c(1,2))
plot(serr, xlab = "Lipids", ylab = "Standardized Performance")
matplot(svimp, xlab = "Genes/Diet/Genotype", ylab = "Standardized VIMP")


}}
\keyword{datasets}
