\name{sidClustering.rfsrc}
\alias{sidClustering.rfsrc}
\alias{sidClustering}
\title{sidClustering using SID (Staggered Interaction Data) for Unsupervised Clustering}

\description{
  Clustering of unsupervised data using SID (Mantero and Ishwaran,
  2021).  Also implements the artificial two-class approach of
  Breiman (2003).
}

\usage{\method{sidClustering}{rfsrc}(data,
  method = "sid",
  k = NULL,
  reduce = TRUE,
  ntree = 500,
  ntree.reduce = function(p, vtry){100 * p / vtry},
  fast = FALSE,
  x.no.sid = NULL,
  use.sid.for.x = TRUE,
  x.only = NULL, y.only = NULL,
  dist.sharpen = TRUE, ...)
}

\arguments{

  \item{data}{A data frame containing the unsupervised data.}

  \item{method}{Clustering method. Default is \code{"sid"}, which implements SID clustering using Staggered Interaction Data (Mantero and Ishwaran, 2021). An alternative approach reformulates the problem as a two-class supervised learning task using artificial data, per Breiman (2003) and Shi-Horvath (2006). Mode 1 is specified via \code{"sh"}, \code{"SH"}, \code{"sh1"}, or \code{"SH1"}; Mode 2 via \code{"sh2"} or \code{"SH2"}. A third method, \code{"unsupv"}, uses a plain unsupervised forest where the data act as both features and responses, split using the multivariate rule. This is faster than SID but may be less accurate.}

  \item{k}{Requested number of clusters. Can be a single integer or a vector. If a scalar, returns a vector assigning each observation to a cluster. If a vector, returns a matrix with one column per requested value of \code{k}, each containing a clustering assignment.}

  \item{reduce}{Logical. If \code{TRUE}, applies a variable reduction step via holdout VIMP. This is conservative and computationally intensive but has strong false discovery control. Applies only when \code{method = "sid"}.}

  \item{ntree}{Number of trees used in the main SID clustering analysis.}

  \item{ntree.reduce}{Number of trees used in the holdout VIMP step during variable reduction. See \code{\link{holdout.vimp}} for details.}

  \item{fast}{Logical. If \code{TRUE}, uses the fast implementation \code{\link{rfsrc.fast}} instead of \code{\link{rfsrc}}. Improves speed at the cost of accuracy.}

  \item{x.no.sid}{Variables to exclude from SID transformation. Can be either a separate data frame (not overlapping with \code{data}) or a character vector of variable names from \code{data}. These variables will be included in the final design matrix without SID processing. Applies only when \code{method = "sid"}.}

  \item{use.sid.for.x}{Logical. If \code{FALSE}, reverses the roles of features and responses in the SID process. Staggered interactions are applied to the outcome rather than to features. This option is slower and generally less effective. Included for legacy compatibility. Applies only when \code{method = "sid"}.}

  \item{x.only}{Character vector specifying which variables to use as features. Applies only when \code{method = "unsupv"}.}

  \item{y.only}{Character vector specifying which variables to use as multivariate responses. Applies only when \code{method = "unsupv"}.}

  \item{dist.sharpen}{Logical. If \code{TRUE} (default), applies Euclidean distance to the forest distance matrix to improve clustering ("distance sharpening"). The resulting distance matrix will not be bounded between 0 and 1. Turning this off speeds up computation but may reduce clustering quality. Applies only when \code{method = "sid"} or \code{"unsupv"}.}

  \item{...}{Additional arguments passed to \code{\link{rfsrc}} to control forest construction.}

}

\details{

Given an unsupervised dataset, random forests is used to compute a distance matrix measuring dissimilarity between all pairs of observations. By default, hierarchical clustering is applied to this distance matrix, although users may apply any other clustering algorithm. See the examples below for alternative workflows.

The default method, \code{method = "sid"}, implements SID clustering (sidClustering). The algorithm begins by enhancing the original feature space using Staggered Interaction Data (SID). This transformation creates:

\itemize{
  \item SID main features: shifted and staggered versions of the original features that are made strictly positive and mutually non-overlapping in range;
  \item SID interaction features: pairwise multiplicative interactions formed between all SID main features.
}

A multivariate random forest is trained to predict SID main features
using the SID interaction features as predictors. The rationale is that
if a feature is informative for distinguishing clusters, it will exhibit
systematic variation across the data space. Because each interaction
feature is uniquely defined by the features it is formed from, node
splits on interaction terms are able to capture and separate such
variation, thus effectively identifying the clusters. See Mantero and
Ishwaran (2021) for further details.

Since SID includes all pairwise interactions, the dimensionality of the feature space grows quadratically with the number of original variables (or worse when factor variables are present). As such, the reduction step using holdout variable importance (VIMP) is strongly recommended (enabled by default). This step can be disabled using \code{reduce = FALSE}, but only when the original feature space is of manageable size.

A second approach, proposed by Breiman (2003) and refined by Shi and Horvath (2006), transforms the unsupervised task into a two-class supervised classification problem. The first class consists of the original data, while the second class is generated artificially. The goal is to separate real data from synthetic data. A proximity matrix is constructed from this supervised model, and the proximity values for the original class are extracted and converted into a distance matrix (distance = 1 - proximity) for clustering. 

Artificial data can be generated using two modes:
\itemize{
  \item \code{mode 1} (default): draws random values from the empirical distribution of each feature;
  \item \code{mode 2}: draws uniformly between the observed minimum and maximum of each feature.
}

This method is invoked by setting \code{method = "sh"}, \code{"sh1"}, or \code{"sh2"}. Mantero and Ishwaran (2021) found that while this approach works in certain settings, it can fail when clusters exist in lower-dimensional subspaces (e.g., when defined by interactions or involving both factors and continuous variables). Among the two modes, mode 1 is generally more robust.

The third method, \code{method = "unsupv"}, trains a multivariate forest using the data both as predictors and as responses. The multivariate splitting rule is applied at each node. This method is fast and simple but may be less accurate compared to SID clustering.

The package includes a helper function \code{sid.perf.metric} for evaluating clustering performance using a normalized score; smaller values indicate better performance. See Mantero and Ishwaran (2021) for theoretical background and empirical benchmarking.
}

\value{

A list with the following components:

\item{clustering}{A vector or matrix assigning each observation to a cluster. If multiple values of \code{k} were specified, this is a matrix with one column per clustering solution.}

\item{rf}{The trained random forest object used in the clustering procedure. This is typically a multivariate forest (for \code{method = "sid"} or \code{"unsupv"}) or a classification forest (for Breiman-style methods).}

\item{dist}{The distance matrix computed from the forest. Used for clustering. For \code{method = "sid"}, this is based on the forest dissimilarity; for Breiman/SH methods, this is one minus the proximity matrix.}

\item{sid}{The SID-transformed data used in the clustering (applies only to \code{method = "sid"}). Provided as a list with separate components for the staggered features and their interactions, corresponding to outcomes and predictors in the multivariate forest.}
}

  
\author{Hemant Ishwaran and Udaya B. Kogalur}

\references{

  Breiman, L. (2003). \emph{Manual on setting up, using and
  understanding random forest, V4.0}.  University of California
  Berkeley, Statistics Department, Berkeley.
  
  Mantero A. and Ishwaran H. (2021).  Unsupervised random forests.
  \emph{Statistical Analysis and Data Mining}, 14(2):144-167.
  
  Shi, T. and Horvath, S. (2006). Unsupervised learning with random forest
  predictors. \emph{Journal of Computational and Graphical Statistics},
  15(1):118-138.  
}

\seealso{
  \command{\link{rfsrc}},
  \command{\link{rfsrc.fast}}
}

\examples{
\donttest{
## ------------------------------------------------------------
## mtcars example
## ------------------------------------------------------------

## default SID method 
o1 <- sidClustering(mtcars)
print(split(mtcars, o1$cl[, 10]))

## using artifical class approach
o1.sh <- sidClustering(mtcars, method = "sh")
print(split(mtcars, o1.sh$cl[, 10]))


## ------------------------------------------------------------
## glass data set
## ------------------------------------------------------------

if (library("mlbench", logical.return = TRUE)) {

  ## this is a supervised problem, so we first strip the class label
  data(Glass)
  glass <- Glass
  y <- Glass$Type
  glass$Type <- NULL

  ## default SID call 
  o2 <- sidClustering(glass, k = 6)
  print(table(y, o2$cl))
  print(sid.perf.metric(y, o2$cl))

  ## compare with Shi-Horvath mode 1 
  o2.sh <- sidClustering(glass, method = "sh1", k = 6)
  print(table(y, o2.sh$cl))
  print(sid.perf.metric(y, o2.sh$cl))

  ## plain-vanilla unsupervised analysis
  o2.un <- sidClustering(glass, method = "unsupv", k = 6)
  print(table(y, o2.un$cl))
  print(sid.perf.metric(y, o2.un$cl))

}

## ------------------------------------------------------------
## vowel data set
## ------------------------------------------------------------

if (library("mlbench", logical.return = TRUE) &&
    library("cluster", logical.return = TRUE)) {

  ## strip the class label
  data(Vowel)
  vowel <- Vowel
  y <- Vowel$Class
  vowel$Class <- NULL

  ## SID 
  o3 <- sidClustering(vowel, k = 11)
  print(table(y, o3$cl))
  print(sid.perf.metric(y, o3$cl))

  ## compare to Shi-Horvath which performs poorly in
  ## mixed variable settings
  o3.sh <- sidClustering(vowel, method = "sh1", k = 11)
  print(table(y, o3.sh$cl))
  print(sid.perf.metric(y, o3.sh$cl))

  ## Shi-Horvath improves with PAM clustering
  ## but still not as good as SID
  o3.sh.pam <- pam(o3.sh$dist, k = 11)$clustering
  print(table(y, o3.sh.pam))
  print(sid.perf.metric(y, o3.sh.pam))

  ## plain-vanilla unsupervised analysis
  o3.un <- sidClustering(vowel, method = "unsupv", k = 11)
  print(table(y, o3.un$cl))
  print(sid.perf.metric(y, o3.un$cl))

}

## ------------------------------------------------------------
##  two-d V-shaped cluster (y=x, y=-x) sitting in 12-dimensions 
##  illustrates superiority of SID to Breiman/Shi-Horvath
## ------------------------------------------------------------

p <- 10
m <- 250
n <- 2 * m
std <- .2

x <- runif(n, 0, 1)
noise <- matrix(runif(n * p, 0, 1), n)
y <- rep(NA, n)
y[1:m] <- x[1:m] + rnorm(m, sd = std)
y[(m+1):n] <- -x[(m+1):n] + rnorm(m, sd = std)
vclus <- data.frame(clus = c(rep(1, m), rep(2,m)), x = x, y = y, noise)

## SID
o4 <- sidClustering(vclus[, -1], k = 2)
print(table(vclus[, 1], o4$cl))
print(sid.perf.metric(vclus[, 1], o4$cl))

## Shi-Horvath
o4.sh <- sidClustering(vclus[, -1], method = "sh1", k = 2)
print(table(vclus[, 1], o4.sh$cl))
print(sid.perf.metric(vclus[, 1], o4.sh$cl))

## plain-vanilla unsupervised analysis
o4.un <- sidClustering(vclus[, -1], method = "unsupv", k = 2)
print(table(vclus[, 1], o4.un$cl))
print(sid.perf.metric(vclus[, 1], o4.un$cl))


## ------------------------------------------------------------
##  two-d V-shaped cluster using fast random forests
## ------------------------------------------------------------

o5 <- sidClustering(vclus[, -1], k = 2, fast = TRUE)
print(table(vclus[, 1], o5$cl))
print(sid.perf.metric(vclus[, 1], o5$cl))


}}
\keyword{unsupervised}
\keyword{clustering}

