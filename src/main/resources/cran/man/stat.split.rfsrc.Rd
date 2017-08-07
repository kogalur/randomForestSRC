\name{stat.split}
\alias{stat.split}
\alias{stat.split.rfsrc}
\title{Acquire Split Statistic Information}
\description{
  Extract split statistic information from the forest.  The function
  returns a list of length \code{ntree}, in which each element
  corresponds to a tree.  The element [[b]] is itself a vector of length
  xvar.names identified by its x-variable name.  Each element [[b]]$xvar
  contains the complete list of splits on xvar with associated
  identifying information.  The information is as follows:
  
  \enumerate{
    
    \item \emph{treeID} Tree identifier.
    \item \emph{nodeID} Node identifier.
    \item \emph{parmID} Variable indentifier.
    \item \emph{contPT} Value node was split in the case of a
    continuous variable.
    \item \emph{mwcpSZ} Size of the multi-word complementary pair
    in the case of a factor split.
    \item \emph{dpthID} Zero (0) based depth of split.
    \item \emph{spltTY} Split type for parent node:
    
    \tabular{lll}{
      bit 1 \tab  bit 0  \tab meaning\cr
      ----- \tab  -----  \tab ------- \cr
      0   \tab    0    \tab 0 = both daughters have valid splits\cr 
      0   \tab    1    \tab 1 = only the right daughter is terminal\cr 
      1   \tab    0    \tab 2 = only the left daughter is terminal\cr 
      1   \tab    1    \tab 3 = both daughters are terminal\cr
    }
    
	  
    \item \emph{spltEC} End cut statistic for real valued variables
    between [0,0.5] that is small when the split is towards the edge and
    large when the split is towards the middle.  Subtracting this value
    from 0.5 yields the end cut statistic studied in Ishwaran (2014) and
    is a way to identify ECP behavior (end cut preference behavior).
    
    \item \emph{spltST} Split statistic:
    
    \enumerate{
      \item For objects of class (rfsrc, grow), this is the split
      statistic that resulted in the variable being choosen for
      the split.
      \item For an object of class (rfsrc, pred) this is the
      variance of the response within the node for the test data.
      This value is relevant only for real valued responses.  In
      classification and survival, it is not relevant.
    }
  }           
}
\usage{\method{stat.split}{rfsrc}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{(rfsrc, grow)},
    \code{(rfsrc, synthetic)}  or \code{(rfsrc,
	    predict)}}
  \item{...}{Further arguments passed to or from other methods.}
}
\value{
  Invisibly, a list with the following components:
  \item{...}{...}
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran H. (2015).  The effect of splitting on random forests.
  \emph{Machine Learning}, 99:75-118.
}
\examples{
\dontrun{
## run a forest, then make a call to stat.split
grow.obj <- rfsrc(mpg ~., data = mtcars, membership=TRUE, statistics=TRUE)
stat.obj <- stat.split(grow.obj)

## nice wrapper to extract split-statistic for desired variable
## for continuous variables plots ECP data
get.split <- function(splitObj, xvar, inches = 0.1, ...) {
  which.var <- which(names(splitObj[[1]]) == xvar)
  ntree <- length(splitObj)
  stat <- data.frame(do.call(rbind, sapply(1:ntree, function(b) {
    splitObj[[b]][which.var]})))
  dpth <- stat$dpthID
  ecp <- 1/2 - stat$spltEC
  sp <- stat$contPT
  if (!all(is.na(sp))) {
    fgC <- function(x) {
      as.numeric(as.character(cut(x, breaks = c(-1, 0.2, 0.35, 0.5),
      labels = c(1, 4, 2))))
    }
    symbols(jitter(sp), jitter(dpth), ecp, inches = inches, bg = fgC(ecp),
      xlab = xvar, ylab = "node depth", ...)
    legend("topleft", legend = c("low ecp", "med ecp", "high ecp"),
      fill = c(1, 4, 2))
   }
  invisible(stat)
}

## use get.split to investigate ECP behavior of variables
get.split(stat.obj, "disp")
}
}
\keyword{splitting behavior}