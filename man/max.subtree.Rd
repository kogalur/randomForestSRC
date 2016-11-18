\name{max.subtree}
\alias{max.subtree}
\alias{max.subtree.rfsrc}
\title{Acquire Maximal Subtree Information}
\description{
    Extract maximal subtree information from a RF-SRC object.  Used for
    variable selection and identifying interactions between variables.
}
\usage{\method{max.subtree}{rfsrc}(object,
  max.order = 2, sub.order = FALSE, conservative = FALSE, ...)
}
\arguments{
  \item{object}{An object of class \code{(rfsrc, grow)} or \code{(rfsrc,
	    forest)}}. 
  \item{max.order}{Non-negative integer specifying the target number
    of order depths.  Default is to return the first and second order
    depths.  Used to identify predictive variables.  Setting
    \option{max.order=0} returns the first order depth for each
    variable by tree.  A side effect is that \option{conservative} is
    automatically set to \code{FALSE}.}
  \item{sub.order}{Set this value to \code{TRUE} to return the
    minimal depth of each variable relative to another variable.
    Used to identify interrelationship between variables.  See
    details below.}
  \item{conservative}{If \code{TRUE}, the threshold value for selecting
    variables is calculated using a conservative marginal
    approximation to the minimal depth distribution (the method used
    in Ishwaran et al. 2010).  Otherwise, the minimal depth
    distribution is the tree-averaged distribution.  The latter method
    tends to give larger threshold values and discovers more
    variables, especially in high-dimensions.}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  The maximal subtree for a variable \emph{x} is the largest subtree
  whose root node splits on \emph{x}.  Thus, all parent nodes of
  \emph{x}'s maximal subtree have nodes that split on variables other
  than \emph{x}.  The largest maximal subtree possible is the root
  node.  In general, however, there can be more than one maximal
  subtree for a variable. A maximal subtree may also not exist if
  there are no splits on the variable.  See Ishwaran et al. (2010,
  2011) for details.
  
  The minimal depth of a maximal subtree (the first order depth)
  measures predictiveness of a variable \emph{x}.  It equals the
  shortest distance (the depth) from the root node to the parent node
  of the maximal subtree (zero is the smallest value possible).  The
  smaller the minimal depth, the more impact \emph{x} has on
  prediction. The mean of the minimal depth distribution is used as
  the threshold value for deciding whether a variable's minimal depth
  value is small enough for the variable to be classified as strong.
  
  The second order depth is the distance from the root node to the
  second closest maximal subtree of \emph{x}.  To specify the target
  order depth, use the \code{max.order} option (e.g., setting
  \option{max.order=2} returns the first and second order depths).
  Setting \option{max.order=0} returns the first order depth for each
  variable for each tree.
  
  Set \option{sub.order=TRUE} to obtain the minimal depth of a
  variable relative to another variable.  This returns a
  \code{p}x\code{p} matrix, where \code{p} is the number of variables,
  and entries [i][j] are the normalized relative minimal depth of a
  variable [j] within the maximal subtree for variable [i], where
  normalization adjusts for the size of [i]'s maximal subtree.  Entry
  [i][i] is the normalized minimal depth of i relative to the root
  node.  The matrix should be read by looking across rows (not down
  columns) and identifies interrelationship between variables.  Small
  [i][j] entries indicate interactions.  See
  \command{find.interaction} for related details.
  
  For competing risk data, maximal subtree analyses are unconditional
  (i.e., they are non-event specific).
}
\value{
  Invisibly, a list with the following components:
  \item{order}{Order depths for a given variable up to \code{max.order}
    averaged over a tree and the forest.  Matrix of dimension
    \code{p}x\code{max.order}.  If \option{max.order=0}, a matrix of
    \code{p}x\code{ntree} is returned containing the first order depth
    for each variable by tree.}
  \item{count}{Averaged number of  maximal subtrees, normalized by
    the size of a tree, for each variable.}
  \item{nodes.at.depth}{Number of non-terminal nodes by depth for each tree.}  
  \item{sub.order}{Average minimal depth of a variable relative to another
    variable.  Can be \code{NULL}.}
  \item{threshold}{Threshold value (the mean minimal depth) used to
    select variables.}
  \item{threshold.1se}{Mean minimal depth plus one standard error.}
  \item{topvars}{Character vector of names of the final selected
    variables.}  
  \item{topvars.1se}{Character vector of names of the final selected
    variables using the 1se threshold rule.}  
  \item{percentile}{Minimal depth percentile for each variable.}
  \item{density}{Estimated minimal depth density.}
  \item{second.order.threshold}{Threshold for second order depth.}
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
  \command{\link{find.interaction}},
  \command{\link{var.select}},
  \command{\link{vimp}}
}
\examples{
\dontrun{
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
