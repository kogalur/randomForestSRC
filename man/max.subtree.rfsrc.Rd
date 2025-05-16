\name{max.subtree.rfsrc}
\alias{max.subtree.rfsrc}
\alias{max.subtree}
\title{Acquire Maximal Subtree Information}

\description{
    Extract maximal subtree information from a RF-SRC object.  Used for
    variable selection and identifying interactions between variables.
}

\usage{\method{max.subtree}{rfsrc}(object,
  max.order = 2, sub.order = FALSE, conservative = FALSE, ...)
}

\arguments{

  \item{object}{An object of class \code{(rfsrc, grow)} or \code{(rfsrc, forest)}.}

  \item{max.order}{Non-negative integer specifying the maximum interaction order for which minimal depth is calculated. Defaults to 2. Set \code{max.order=0} to return first-order depths only. When \code{max.order=0}, \code{conservative} is automatically set to \code{FALSE}.}

  \item{sub.order}{Logical. If \code{TRUE}, returns the minimal depth of each variable conditional on every other variable. Useful for investigating variable interdependence. See \code{Details}.}

  \item{conservative}{Logical. If \code{TRUE}, uses a conservative threshold for selecting variables based on the marginal minimal depth distribution (Ishwaran et al., 2010). If \code{FALSE}, uses the tree-averaged distribution, which is less conservative and typically identifies more variables in high-dimensional settings.}

  \item{...}{Additional arguments passed to or from other methods.}
}

\details{
  The maximal subtree for a variable \emph{x} is the largest subtree in which the root node splits on \emph{x}. The largest possible maximal subtree is the full tree (root node), though multiple maximal subtrees may exist for a variable. A variable may also have no maximal subtree if it is never used for splitting. See Ishwaran et al. (2010, 2011) for further discussion.

  The minimal depth of a maximal subtree-called the \emph{first-order depth}-quantifies the predictive strength of a variable. It is defined as the distance from the root node to the parent of the closest maximal subtree for \emph{x}. Smaller values indicate stronger predictive impact. A variable is flagged as strong if its minimal depth is below the mean of the minimal depth distribution.

  The \emph{second-order depth} is the distance from the root to the second-closest maximal subtree of \emph{x}. To request depths beyond first order, use the \code{max.order} option (e.g., \code{max.order = 2} returns both first and second-order depths). Set \code{max.order = 0} to retrieve first-order depths for each variable in each tree.

  Set \code{sub.order = TRUE} to obtain the relative minimal depth of
  each variable \emph{j} within the maximal subtree of another variable
  \emph{i}. This returns a \code{p x p} matrix (with \code{p} the number
  of variables) whose entry (i,j) is the normalized relative depth of
  \emph{j} in \emph{i}'s subtree. Entry (i,i) gives the depth of
  \emph{i} relative to the root. Read the matrix across rows to assess
  inter-variable relationships: small (i,j) entries suggest interactions
  between variables \emph{i} and \emph{j}. See \code{find.interaction}
  for further details.

  For competing risks, all analyses are unconditional (non-event specific).
}


\value{

  Invisibly returns a list with the following components:

  \item{order}{Matrix of order depths for each variable up to
  \code{max.order}, averaged over trees. The matrix has \code{p} rows
  and \code{max.order} columns, where \code{p} is the number of
  variables. If \code{max.order = 0}, returns a matrix of dimension
  \code{p x ntree} containing first-order depths for each variable by
  tree.}

  \item{count}{Average number of maximal subtrees per variable, normalized by tree size.}

  \item{nodes.at.depth}{List of vectors recording the number of non-terminal nodes at each depth level for each tree.}

  \item{sub.order}{Matrix of average minimal depths of each variable relative to others (i.e., conditional minimal depth matrix). \code{NULL} if \code{sub.order = FALSE}.}

  \item{threshold}{Threshold value for selecting strong variables based on the mean of the minimal depth distribution.}

  \item{threshold.1se}{Conservative threshold equal to the mean minimal depth plus one standard error.}

  \item{topvars}{Character vector of selected variable names using the \code{threshold} criterion.}

  \item{topvars.1se}{Character vector of selected variable names using the \code{threshold.1se} criterion.}

  \item{percentile}{Percentile value of minimal depth for each variable.}

  \item{density}{Estimated density of the minimal depth distribution.}

  \item{second.order.threshold}{Threshold used for selecting strong second-order depth variables.}

}

\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and 
  Lauer M.S. (2010).  High-dimensional variable selection for survival
  data.  \emph{J. Amer. Statist. Assoc.}, 105:205-217.
  
  Ishwaran H., Kogalur U.B., Chen X. and Minn A.J. (2011).  Random
  survival forests for high-dimensional data. \emph{Statist. Anal. Data
    Mining}, 4:115-132.
}
\seealso{
  \command{\link{holdout.vimp.rfsrc}},
  \command{\link{vimp.rfsrc}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## survival analysis
## first and second order depths for all variables
## ------------------------------------------------------------

data(veteran, package = "randomForestSRC")
v.obj <- rfsrc(Surv(time, status) ~ . , data = veteran)
v.max <- max.subtree(v.obj)

# first and second order depths
print(round(v.max$order, 3))

# the minimal depth is the first order depth
print(round(v.max$order[, 1], 3))

# strong variables have minimal depth less than or equal
# to the following threshold
print(v.max$threshold)

# this corresponds to the set of variables
print(v.max$topvars)

## ------------------------------------------------------------
## regression analysis
## try different levels of conservativeness
## ------------------------------------------------------------

mtcars.obj <- rfsrc(mpg ~ ., data = mtcars)
max.subtree(mtcars.obj)$topvars
max.subtree(mtcars.obj, conservative = TRUE)$topvars
}
}
\keyword{variable selection}
