\name{find.interaction.rfsrc}
\alias{find.interaction}
\alias{find.interaction.rfsrc}
\title{Find Interactions Between Pairs of Variables}
\description{
  Find pairwise interactions between variables.
}
\usage{\method{find.interaction}{rfsrc}(object, xvar.names, cause, m.target,
  importance = c("permute", "random", "anti",
                 "permute.ensemble", "random.ensemble", "anti.ensemble"),
  method = c("maxsubtree", "vimp"), sorted = TRUE, nvar, nrep = 1, 
  na.action = c("na.omit", "na.impute", "na.random"),
  seed = NULL, do.trace = FALSE, verbose = TRUE, ...)
}
\arguments{
  \item{object}{An object of class \code{(rfsrc, grow)} or
    \code{(rfsrc, forest)}.}  
  \item{xvar.names}{Character vector of names of target x-variables.
    Default is to use all variables.}
  \item{cause}{For competing risk families, integer value between 1
    and \code{J} indicating the event of interest, where \code{J} is
    the number of event types. The default is to use the first event
    type.}
  \item{m.target}{Character value for multivariate families
    specifying the target outcome to be used.  If left unspecified, the
    algorithm will choose a default target.}
  \item{importance}{Type of variable importance (VIMP). See
    \command{rfsrc} for details.}
  \item{method}{Method of analysis: maximal subtree or VIMP.  See details
    below.}
  \item{sorted}{Should variables be sorted by VIMP?  Does not apply for
    competing risks.}
  \item{nvar}{Number of variables to be used.}
  \item{nrep}{Number of Monte Carlo replicates when \option{method="vimp"}.}
  \item{na.action}{Action to be taken if the data contains \code{NA}
    values.  Applies only when \option{method="vimp"}.}
  \item{seed}{Seed for random number generator.  Must be a negative
    integer.}
  \item{do.trace}{Number of seconds between updates to the user on
    approximate time to completion.}
  \item{verbose}{Set to \code{TRUE} for verbose output.}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Using a previously grown forest, identify pairwise interactions for
  all pairs of variables from a specified list.  There are two
  distinct approaches specified by the option \option{method}.
  
  \enumerate{
    \item \option{method="maxsubtree"}
    
    This invokes a maximal subtree analysis.  In this case, a matrix
    is returned where entries [i][i] are the normalized minimal depth
    of variable [i] relative to the root node (normalized wrt the
    size of the tree) and entries [i][j] indicate the normalized
    minimal depth of a variable [j] wrt the maximal subtree for
    variable [i] (normalized wrt the size of [i]'s maximal subtree).
    Smaller [i][i] entries indicate predictive variables.  Small
    [i][j] entries having small [i][i] entries are a sign of an
    interaction between variable i and j (note: the user should scan
    rows, not columns, for small entries).  See Ishwaran et
    al. (2010, 2011) for more details.
    
    \item \option{method="vimp"}
    
    This invokes a joint-VIMP approach.  Two variables are paired and
    their paired VIMP calculated (refered to as 'Paired' importance).
    The VIMP for each separate variable is also calculated.  The sum
    of these two values is refered to as 'Additive' importance.  A
    large positive or negative difference between 'Paired' and
    'Additive' indicates an association worth pursuing if the
    univariate VIMP for each of the paired-variables is reasonably
    large.  See Ishwaran (2007) for more details.
  }
  
  Computations might be slow depending upon the size of the data
  and the forest.  In such cases, consider setting \option{nvar} to a
  smaller number.  If \option{method="maxsubtree"}, consider using a
  smaller number of trees in the original grow call.
  
  If \option{nrep} is greater than 1, the analysis is repeated
  \code{nrep} times and results averaged over the replications
  (applies only when \option{method="vimp"}).
}
\value{
  Invisibly, the interaction table (a list for competing risk data)
  or the maximal subtree matrix.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran H. (2007).  Variable importance in binary regression
  trees and forests,  \emph{Electronic J. Statist.}, 1:519-537.
  
  Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and 
  Lauer M.S. (2010).  High-dimensional variable selection for survival
  data.  \emph{J. Amer. Statist. Assoc.}, 105:205-217.
  
  Ishwaran H., Kogalur U.B., Chen X. and Minn A.J. (2011).  Random
  survival forests for high-dimensional data. \emph{Statist. Anal. Data
    Mining}, 4:115-132.
}
\seealso{
  \command{\link{holdout.vimp.rfsrc}},
  \command{\link{max.subtree.rfsrc}},
  \command{\link{vimp.rfsrc}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## find interactions, survival setting
## ------------------------------------------------------------

data(pbc, package = "randomForestSRC") 
pbc.obj <- rfsrc(Surv(days,status) ~ ., pbc, importance = TRUE)
find.interaction(pbc.obj, method = "vimp", nvar = 8)

## ------------------------------------------------------------
## find interactions, competing risks
## ------------------------------------------------------------

data(wihs, package = "randomForestSRC")
wihs.obj <- rfsrc(Surv(time, status) ~ ., wihs, nsplit = 3, ntree = 100,
                       importance = TRUE)
find.interaction(wihs.obj)
find.interaction(wihs.obj, method = "vimp")

## ------------------------------------------------------------
## find interactions, regression setting
## ------------------------------------------------------------

airq.obj <- rfsrc(Ozone ~ ., data = airquality, importance = TRUE)
find.interaction(airq.obj, method = "vimp", nrep = 3)
find.interaction(airq.obj)

## ------------------------------------------------------------
## find interactions, classification setting
## ------------------------------------------------------------

iris.obj <- rfsrc(Species ~., data = iris, importance = TRUE)
find.interaction(iris.obj, method = "vimp", nrep = 3)
find.interaction(iris.obj)

## ------------------------------------------------------------
## interactions for multivariate mixed forests
## ------------------------------------------------------------

mtcars2 <- mtcars
mtcars2$cyl <- factor(mtcars2$cyl)
mtcars2$carb <- factor(mtcars2$carb, ordered = TRUE)
mv.obj <- rfsrc(cbind(carb, mpg, cyl) ~., data = mtcars2, importance = TRUE)
find.interaction(mv.obj, method = "vimp", outcome.target = "carb")
find.interaction(mv.obj, method = "vimp", outcome.target = "mpg")
find.interaction(mv.obj, method = "vimp", outcome.target = "cyl")
}
}
\keyword{variable selection}

