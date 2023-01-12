\name{sidClustering.rfsrc}
\alias{sidClustering.rfsrc}
\alias{sidClustering}
\alias{sid.perf.metric}
\title{sidClustering using SID (Staggered Interaction Data) for Unsupervised Clustering}
\description{
  Clustering of unsupervised data using SID (Mantero and Ishwaran,
  2020).  Also implements the artificial two-class approach of
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
  \item{data}{Data frame containing the unsupervised data.}

  \item{method}{The method used for unsupervised clustering.  Default is
    "sid" which implements sidClustering using SID (Staggered
    Interaction Data; see Mantero and Ishwaran, 2020).  A second
    approach transforms the unsupervised learning problem into a
    two-class supervised problem (Breiman, 2003) using artificial data
    created using mode 1 or mode 2 of Shi-Horvath (2006).  This approach
    is specified by any one of the following: "sh", "SH", "sh1", "SH1"
    for mode 1, or "sh2", "SH2" for mode 2.  Finally, a third approach
    is a plain vanilla method where the data are used both as features
    and response with splitting implemented using the multivariate
    splitting rule.  This is faster than sidClustering but potentially
    less accurate.  This method is specified using "unsupv".}
  \item{k}{Requested number of clusters.  Can be a number or a vector.
    If a fixed number, returns a vector recording clustering of data.
    If a vector, returns a matrix of clusters with each column recording
    the clustering of the data for the specified number of clusters.}
  \item{reduce}{Apply dimension reduction? Uses holdout vimp which is
    computationally intensive and conservative but has good false
    discovery properties.  Only applies to \code{method="sid"}.}
  \item{ntree}{Number of trees used by sidClustering in the main analysis.}
  \item{ntree.reduce}{Number of trees used by holdout vimp in the
    reduction step.  See \command{holdout.vimp} for details.}
  \item{fast}{Use fast random forests, \code{rfsrcFast}, in place of
    \code{rfsrc}?  Improves speed but is less accurate.}
  \item{x.no.sid}{Features not to be "sid-ified": meaning that these
      features are to be included in the final design matrix without SID
      processing.  Can be either a data frame (should not overlap with
      \code{data}), or a character vector containing the names of
      features from the original data that the user wishes to protect
      from sidification.  Applies only to \code{method="sid"}.}
  \item{use.sid.for.x}{If FALSE, reverses features and outcomes in the
    SID analysis.  Thus, staggered interactions are used for the
    outcomes rather than staggered features.  This is much slower and is
    generally much less effective.  This option is only retained for
    legacy reasons. Applies only to \code{method="sid"}.}
  \item{x.only}{Use only these variables for the features.  Applies only
    to \code{method="unsupv"}.}
  \item{y.only}{Use only these variables for the multivariate outcomes.  Applies only
    to \code{method="unsupv"}.}
  \item{dist.sharpen}{By default, distance sharpening is requested,
    which applies Euclidean distance to the random forest distance
    matrix to sharpen it.  Because of this, the returned distance matrix
    will not have values between 0 and 1 (as for random forests distance) when
    this option is in effect.  Distance sharpening is a useful, but slow
    step. Set this option to \code{FALSE} to improve computational
    times, however clustering performance will not be as good.  Applies
    only when \code{method="sid"} or \code{method="unsupv"}.}
  \item{...}{Further arguments to be passed to the \command{rfsrc}
    function to specify random forest parameters.}
}
\details{

  Given an unsupervised data set, random forests is used to calculate
  the distance between all pairs of data points.  The distance matrix is
  used for clustering the unsupervised data where the default is to use
  hierarchcial clustering.  Users can apply other clustering procedures
  to the distance matrix.  See the examples below.

  The default method, \code{method="sid"}, implements sidClustering.
  The sidClustering algorithm begins by first creating an enhanced SID
  (Staggered Interaction Data) feature space by sidification of the
  original variables.  Sidification results in: (a) SID main features
  which are the original features that have been shifted in order to
  make them strictly positive and staggered so all of their ranges are
  mutually exclusive; and (b) SID interaction features which are the
  multiplicative interactions formed between every pair of SID main
  features.  Multivariate random forests are then trained to predict the
  main SID features using the interaction SID features as predictors.
  The basic premise is if features are informative for clusters, then
  they will vary over the space in a systematic manner, and because each
  SID interaction feature is uniquely determined by the original feature
  values used to form the interaction, cuts along the SID interaction
  feature will be able to find the regions where the informative
  features vary by cluster, thereby not only reducing impurity, but also
  separating the clusters which are dependent on those features.  See
  Mantero and Ishwaran (2020) for details.

  Because SID uses all pairwise interactions, the dimension of the
  feature space is proportional to the square of the number of original
  features (or even larger if factors are present).  Thus it is helpful
  to reduce the feature space. The reduction step (applied by default)
  utilizes holdout VIMP to accomplish this.  It is recommended this step
  be skipped only when the dimension is reasonably small.  For very
  large data sets this step may be slow.
  
  A second approach (Breiman, 2003; Shi-Horvath, 2006) transforms the
  unsupervised learning problem into a two class supervised problem.
  The first class consists of the original observations, while the
  second class is artificially created.  The idea is that in detecting
  the first class out of the second, the model will generate the random
  forest proximity between observations of which those for the original
  class can be extracted and used for clustering.  Note in this approach
  the distance matrix is defined to equal one minus the proximity.  This
  is unlike the distance matrix from SID which is not proximity based.
  Artificial data is created using "mode 1" or "mode 2" of Shi-Horvath
  (2006).  Mode 1 randomly draws from each set of observed features.
  Mode 2 draws a uniform value from the minimum and maximum values of a
  feature.
  
  Mantero and Ishwaran (2020) studied both methods and found SID worked
  well in all settings, whereas Breiman/Shi-Horvath was sensitive to
  cluster structure.  Performance was poor when clusters were hidden in
  lower dimensional subspaces; for example when interactions were
  present or in mixed variable settings (factors/continuous variables).
  See the V-shaped cluster example below.  Generally Shi-Horvath mode 1
  outperforms mode 2.

  Finally, a third method where the data is used for both the features
  and outcome is implemented using \code{method="unsupv"}.  Tree nodes
  are split using the multivariate splitting rule.  This is much faster
  than sidClustering but potentially less accurate.

  There is an internal function \command{sid.perf.metric} for evaluating
  performance of the procedures using a normalized measure score.
  Smaller values indicate better performance.  See Mantero and Ishwaran
  (2020) for details.
  
}

\value{ A list with the following components:
  \item{clustering}{Vector or matrix containing indices mapping data
    points to their clusters.}
  \item{rf}{Random forest object (either a multivariate forest or RF-C object).}
  \item{dist}{Distance matrix.}
  \item{sid}{The "sid-ified" data. Conveniently broken up into
    separate values for outcomes and features used by the
    multivariate forest.}
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

