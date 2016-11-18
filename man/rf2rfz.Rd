\name{rf2rfz}
\alias{rf2rfz}
\title{Save RF-SRC in .rfz Compressed Format}
\description{
  \code{rf2rfz} saves a RF-SRC object as a \code{.rfz} compressed file
  that is readable by the \pkg{randomForestSRC} Java plugin that is
  capable of visualizing the trees in the forest.
  
}
\usage{
rf2rfz(object, forestName = NULL, ...)
}
\arguments{
  \item{object}{An object of class \code{(rfsrc, grow)} or \code{(rfsrc,
	    forest)}. Requires \option{forest=TRUE} in the
    original \command{rfsrc} call.}
  \item{forestName}{The desired prefix name for forest as a string.}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  An \code{.rfz} compressed file is actually a \code{.zip} file
  consisting of three files. The first is an ASCII file of type
  \code{.txt} containing the \code{$nativeArray} component of the
  forest.  The second is an ASCII file of type \code{.factor.txt}
  containing the \code{$nativefactorArray} component of the forest.  The
  third is an ASCII file of type \code{.xml} containing the \code{PMML
    DataDictionary} component.
  
  PMML or the Predictive Model Markup Language is an XML based language
  which provides a way for applications to define statistical and data
  mining models and to share models between PMML compliant applications.
  More information about PMML and the Data Mining Group can be found at
  http:
  
  The function \command{rf2rfz} is used to import the geometry of the
  forest to the RF-SRC Java plugin that is capable of visualizing the
  trees in the forest.
  
  The geometry of the forest is saved as a file called
  \code{forestName.rfz} in the users working directory.  This file can
  then be read by the \pkg{randomForestSRC} Java plugin.
  
  Contact the authors on downloading the Java plugin.
}    
\value{
  None.
}
\note{
  Contact the authors on downloading the Java plugin.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  http:
}	

\seealso{
  \command{\link{rfsrc}}
}
\examples{
\dontrun{
# Example 1:  Growing a forest, saving it as a \emph{.rfz} file ready
# for import into the Java plugin.

library("XML")

data(veteran, package = "randomForestSRC")
v.obj <- rfsrc(Surv(time, status)~., data = veteran)
rf2rfz(v.obj$forest, forestName = "veteran")
}
}
\keyword{forest}
