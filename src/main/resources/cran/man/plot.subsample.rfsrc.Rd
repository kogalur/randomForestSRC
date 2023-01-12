\name{plot.subsample.rfsrc}
\alias{plot.subsample.rfsrc}
\alias{plot.subsample}
\title{Plot Subsampled VIMP Confidence Intervals}
\description{
  Plots VIMP (variable importance) confidence regions obtained from
  subsampling a forest.
}
\usage{\method{plot.subsample}{rfsrc}(x, alpha = .01, xvar.names,
 standardize = TRUE, normal = TRUE, jknife = TRUE,
 target, m.target = NULL, pmax = 75, main = "", sorted = TRUE, ...)}
\arguments{

  \item{x}{An object obtained from calling \command{subample}.}

  \item{alpha}{Desired level of significance.}

  \item{xvar.names}{Names of the x-variables to be used.  If not
  specified all variables used.}
  
  \item{standardize}{Standardize VIMP?  For regression families, VIMP is
    standardized by dividing by the variance and then multipled by 100.
    For all other families, VIMP is scaled by 100.}

  \item{normal}{Use parametric normal confidence regions or
    nonparametric regions?  Generally, parametric regions perform better.}

  \item{jknife}{Use the delete-d jackknife variance estimator?}

  \item{target}{For classification families, an integer or
    character value specifying the class VIMP will be conditioned on
    (default is to use unconditional VIMP).
    For competing risk families, an integer value between
    1 and \code{J} indicating the event VIMP is requested, where \code{J} is
    the number of event types.  The default is to use the first event.}

  \item{m.target}{Character value for multivariate families
    specifying the target outcome to be used.  If left unspecified, the
    algorithm will choose a default target.}

  \item{pmax}{Trims the data to this number of variables (sorted by VIMP).}

  \item{main}{Title used for plot.}
  
  \item{sorted}{Should variables be sorted by importance values?}

  \item{...}{Further arguments that can be passed to \command{bxp}.}
  
}
\details{
  Most of the options used by the R function bxp will work here and can
  be used for customization of plots.  Currently the following
  parameters will work:

  "xaxt", "yaxt", "las", "cex.axis", 
  "col.axis", "cex.main",
  "col.main", "sub", "cex.sub", "col.sub", 
  "ylab", "cex.lab", "col.lab"
		 
}
\value{
  Invisibly, returns the boxplot data that is plotted.
}

\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran H. and Lu M.  (2017).  Standard errors and confidence
  intervals for variable importance in random forest regression,
  classification, and survival.

  Politis, D.N. and Romano, J.P. (1994). Large sample confidence
  regions based on subsamples under minimal assumptions. \emph{The
  Annals of Statistics}, 22(4):2031-2050.

  Shao, J. and Wu, C.J. (1989). A general theory for jackknife variance
  estimation. \emph{The Annals of Statistics}, 17(3):1176-1197.
}
\seealso{
  \command{\link{subsample.rfsrc}}
}
\examples{
\donttest{
o <- rfsrc(Ozone ~ ., airquality)
oo <- subsample(o)
plot.subsample(oo)
plot.subsample(oo, xvar.names = o$xvar.names[1:3])
plot.subsample(oo, jknife = FALSE)
plot.subsample(oo, alpha = .01)
plot(oo,cex.axis=.5)
}
}
\keyword{plot}
