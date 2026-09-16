\name{nutrigenomic}
\docType{data}
\alias{nutrigenomic}
\title{Nutrigenomic Study}
\description{
  Investigates the effects of five dietary treatments on 21 liver lipids
  and 120 hepatic gene expressions in wild-type and PPAR-alpha deficient
  mice. The examples fit a multivariate regression forest with lipid
  expression profiles as the responses and gene expression, diet, and
  genotype as predictors.
}
\details{
  The examples use \code{\link{get.mv.formula}} to construct the response
  formula, and \code{\link{get.mv.error}} and \code{\link{get.mv.vimp}}
  to extract response-specific performance errors and variable importance
  (VIMP). For these regression responses, \code{standardize = TRUE}
  divides each error or importance value by the variance of the
  corresponding response. See
  Extracting Multivariate Values (\code{\link{multivariate.values}})
  for the helper interfaces and returned formats.
}
\references{
  Martin P.G. et al. (2007). Novel aspects of PPAR-alpha-mediated
  regulation of lipid and xenobiotic metabolism revealed through a
  nutrigenomic study. \emph{Hepatology}, 45(3), 767--777.
}
\seealso{
  \code{\link{rfsrc}}, \code{\link{predict.rfsrc}},
  Extracting Multivariate Values (\code{\link{multivariate.values}})
}
\examples{
\donttest{
## ------------------------------------------------------------
## multivariate regression forests using Mahalanobis splitting
## lipids (all real values) used as the multivariate y
## ------------------------------------------------------------

## load the data
data(nutrigenomic, package = "randomForestSRC")

## parse into y and x data
ydta <- nutrigenomic$lipids
xdta <- data.frame(nutrigenomic$genes,
                   diet = nutrigenomic$diet,
                   genotype = nutrigenomic$genotype)

## multivariate mixed forest call
obj <- rfsrc(get.mv.formula(colnames(ydta)),
             data.frame(ydta, xdta),
             importance=TRUE, nsplit = 10,
             splitrule = "mahalanobis")
print(obj)

## ------------------------------------------------------------
## plot the standarized performance and VIMP values
## ------------------------------------------------------------

## acquire the error rate for each of the 21-coordinates 
## standardize to allow for comparison across coordinates
serr <- get.mv.error(obj, standardize = TRUE)

## acquire standardized VIMP 
svimp <- get.mv.vimp(obj, standardize = TRUE)

par(mfrow = c(1,2))
plot(serr, xlab = "Lipids", ylab = "Standardized Performance")
matplot(svimp, xlab = "Genes/Diet/Genotype", ylab = "Standardized VIMP")


## ------------------------------------------------------------
## plot some trees
## ------------------------------------------------------------

plot(get.tree(obj, 1))
plot(get.tree(obj, 2))
plot(get.tree(obj, 3))


## ------------------------------------------------------------
##
## Compare above to (1) user specified covariance matrix
##                  (2) default composite (independent) splitting
##
## ------------------------------------------------------------

## user specified sigma matrix
obj2 <- rfsrc(get.mv.formula(colnames(ydta)),
              data.frame(ydta, xdta),
              importance = TRUE, nsplit = 10,
              splitrule = "mahalanobis",
              sigma = cov(ydta))
print(obj2)

## default independence split rule
obj3 <- rfsrc(get.mv.formula(colnames(ydta)),
              data.frame(ydta, xdta),
              importance=TRUE, nsplit = 10)
print(obj3)

## compare vimp
imp <- data.frame(mahalanobis  = rowMeans(get.mv.vimp(obj,  standardize = TRUE)),
                  mahalanobis2 = rowMeans(get.mv.vimp(obj2, standardize = TRUE)),
                  default      = rowMeans(get.mv.vimp(obj3, standardize = TRUE)))

print(head(100 * imp[order(imp$mahalanobis, decreasing = TRUE), ], 15))

}}
\keyword{datasets}
