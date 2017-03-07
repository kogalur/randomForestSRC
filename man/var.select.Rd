\name{var.select}
\alias{var.select}
\alias{var.select.rfsrc}
\title{Variable Selection}
\description{
  Variable selection using minimal depth.
}
\usage{\method{var.select}{rfsrc}(formula,
  data,
  object,
  cause,
  outcome.target,
  method = c("md", "vh", "vh.vimp"),
  conservative = c("medium", "low", "high"),
  ntree = (if (method == "md") 1000 else 500),
  mvars = (if (method != "md") ceiling(ncol(data)/5) else NULL),
  mtry = (if (method == "md") ceiling(ncol(data)/3) else NULL),
  nodesize = 2, splitrule = NULL, nsplit = 10, xvar.wt = NULL,
  refit = (method != "md"), fast = FALSE,
  
  na.action = c("na.omit", "na.impute"),
  
  
  always.use = NULL, nrep = 50, K = 5, nstep = 1,
  prefit =  list(action = (method != "md"), ntree = 100,
  mtry = 500, nodesize = 3, nsplit = 1),
  do.trace = 0, verbose = TRUE, ...)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit. 
    Must be specified unless \code{object} is given.}
  \item{data}{Data frame containing the y-outcome and x-variables in
    the model. Must be specified unless \code{object} is given.}
  \item{object}{An object of class \code{(rfsrc, grow)}. 
    Not required when \code{formula} and \code{data} are supplied.}
  \item{cause}{Integer value between 1 and \code{J} indicating
	  the event of interest for competing risks, where \code{J} is
	  the number of event types (this option applies only to
	  competing risk families).  The default is to use the first
	  event type.}
  \item{outcome.target}{Character vector for multivariate families
    specifying the target outcomes to be used when VIMP is utilized.
    The default is to use the first coordinate.} 
  \item{method}{Variable selection method:
    \describe{  
      \item{\code{md}:}{minimal depth (default).}
      \item{\code{vh}:}{variable hunting.}
      \item{\code{vh.vimp}:}{variable hunting with VIMP (variable
        importance).}
    }
  }
  \item{conservative}{Level of conservativeness of the thresholding
    rule used in minimal depth selection:
    \describe{
      \item{\code{high}:}{Use the most conservative threshold.}
      \item{\code{medium}:}{Use the default less conservative tree-averaged
        threshold.}
      \item{\code{low}:}{Use the more liberal one standard error rule.}
    }
  }
  \item{ntree}{Number of trees to grow.}
  \item{mvars}{Number of randomly selected variables used in the
    variable hunting algorithm (ignored when \option{method="md"}).}
  \item{mtry}{The mtry value used.}
  \item{nodesize}{Minimum number of unique cases in a terminal node.}
  \item{splitrule}{Splitting rule used.}
  \item{nsplit}{If non-zero, the specified tree splitting rule is
    randomized which significantly increases speed.}
  \item{xvar.wt}{Vector of non-negative weights specifying the
    probability of selecting a variable for splitting a node.  Must be of
    dimension equal to the number of variables.  Default (\code{NULL})
    invokes uniform weighting or a data-adaptive method depending on
    \code{prefit$action}.}
  \item{refit}{Should a forest be refit using the selected variables?}
  \item{fast}{Speeds up the cross-validation used for variable hunting
    for a faster analysis.  See miscellanea below.}
  \item{na.action}{Action to be taken if the data contains \code{NA} values.}
  \item{always.use}{Character vector of variable names to always
    be included in the model selection procedure and in the final
    selected model.}
  \item{nrep}{Number of Monte Carlo iterations of the variable hunting algorithm.}
  \item{K}{Integer value specifying the \code{K}-fold size used in the variable hunting
    algorithm.}
  \item{nstep}{Integer value controlling the step size used in the
    forward selection process of the variable hunting algorithm.
    Increasing this will encourage more variables to be selected.}
  \item{prefit}{List containing parameters used in preliminary forest
    analysis for determining weight selection of variables.  Users can
    set all or some of the following parameters:
    \describe{
      \item{\code{action}:}{Determines how (or if) the preliminary forest is
        fit.  See details below.}
      \item{\code{ntree}:}{Number of trees used in the preliminary analysis.}
      \item{\code{mtry}:}{mtry used in the preliminary analysis.}
      \item{\code{nodesize}:}{nodesize used in the preliminary analysis.}
      \item{\code{nsplit}:}{nsplit value used in the preliminary analysis.}
    }
  }  
  \item{do.trace}{Number of seconds between updates to the user on
    approximate time to completion.}
  \item{verbose}{Set to \code{TRUE} for verbose output.}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  This function implements random forest variable selection using
  tree minimal depth methodology (Ishwaran et al., 2010).  The option
  \option{method} allows for two different approaches:
  
  \enumerate{
    \item \option{method="md"}
    
    Invokes minimal depth variable selection.  Variables are selected
    using minimal depth variable selection.  Uses all data and all
    variables simultaneously.  This is basically a front-end to the
    \command{max.subtree} wrapper.  Users should consult the
    \command{max.subtree} help file for details.
    
    Set \option{mtry} to larger values in high-dimensional problems.
    
    \item \option{method="vh"} or \option{method="vh.vimp"}
    
    Invokes variable hunting.  Variable hunting is used for problems
    where the number of variables is substantially larger than the
    sample size (e.g., p/n is greater than 10).  It is always prefered
    to use \option{method="md"}, but to find more variables, or when
    computations are high, variable hunting may be preferred.
    
    When \option{method="vh"}: Using training data from a stratified
    K-fold subsampling (stratification based on the y-outcomes), a
    forest is fit using \code{mvars} randomly selected variables
    (variables are chosen with probability proportional to weights
    determined using an initial forest fit; see below for more
    details).  The \code{mvars} variables are ordered by increasing
    minimal depth and added sequentially (starting from an initial
    model determined using minimal depth selection) until joint VIMP
    no longer increases (signifying the final model).  A forest is
    refit to the final model and applied to test data to estimate
    prediction error.  The process is repeated \code{nrep} times.
    Final selected variables are the top P ranked variables, where P
    is the average model size (rounded up to the nearest integer) and
    variables are ranked by frequency of occurrence.
    
    The same algorithm is used when \option{method="vh.vimp"}, but
    variables are ordered using VIMP.  This is faster, but not as
    accurate.
  }
  \emph{Miscellanea}
  \enumerate{
	  \item When variable hunting is used, a preliminary forest is run
    and its VIMP is used to define the probability of selecting a
    variable for splitting a node.  Thus, instead of randomly
    selecting \code{mvars} at random, variables are selected with
    probability proportional to their VIMP (the probability is zero
    if VIMP is negative).  A preliminary forest is run once prior
    to the analysis if \code{prefit$action=TRUE}, otherwise it is
    run prior to each iteration (this latter scenario can be slow).
    When \option{method="md"}, a preliminary forest is fit only if
    \code{prefit$action=TRUE}.  Then instead of randomly selecting
    \code{mtry} variables at random, \code{mtry} variables are
    selected with probability proportional to their VIMP.  In all
    cases, the entire option is overridden if \code{xvar.wt} is
    non-null.
	  
	  \item If \code{object} is supplied and \option{method="md"},
	  the grow forest from \code{object} is parsed for minimal depth
	  information.  While this avoids fitting another forest, thus
	  saving computational time, certain options no longer apply.  In
	  particular, the value of \code{cause} plays no role in the
	  final selected variables as minimal depth is extracted from the
	  grow forest, which has already been grown under a preselected
	  \code{cause} specification.  Users wishing to specify
	  \code{cause} should instead use the formula and data interface.
	  Also, if the user requests a prefitted forest via
	  \code{prefit$action=TRUE}, then \code{object} is not used and a
	  refitted forest is used in its place for variable selection.
	  Thus, the effort spent to construct the original grow forest is
	  not used in this case.
    
	  \item If \option{fast=TRUE}, and variable hunting is used, the
	  training data is chosen to be of size n/K, where n=sample size
	  (i.e., the size of the training data is swapped with the test
	  data).  This speeds up the algorithm.  Increasing K also helps.
	  
	  \item Can be used for competing risk data.  When
	  \option{method="vh.vimp"}, variable selection based on VIMP is
	  confined to an event specific cause specified by \code{cause}.
	  However, this can be unreliable as not all y-outcomes can be
	  guaranteed when subsampling (this is true even when stratifed
	  subsampling is used as done here).
  }
}
\value{
  Invisibly, a list with the following components:
  \item{err.rate}{Prediction error for the forest (a vector of
    length \code{nrep} if variable hunting is used).}
  \item{modelsize}{Number of variables selected.}
  \item{topvars}{Character vector of names of the final selected variables.}  
  \item{varselect}{Useful output summarizing the final selected variables.}
  \item{rfsrc.refit.obj}{Refitted forest using the final set of selected variables  
    (requires \option{refit=TRUE}).}
  \item{md.obj}{Minimal depth object.  \code{NULL} unless \option{method="md"}.}
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
  \command{\link{max.subtree}},
  \command{\link{vimp}}
}
\examples{
\dontrun{
## ------------------------------------------------------------
## Minimal depth variable selection
## survival analysis
## ------------------------------------------------------------

data(pbc, package = "randomForestSRC")
pbc.obj <- rfsrc(Surv(days, status) ~ ., pbc, nsplit = 10, importance = TRUE)

# default call corresponds to minimal depth selection
vs.pbc <- var.select(object = pbc.obj)
topvars <- vs.pbc$topvars

# the above is equivalent to
max.subtree(pbc.obj)$topvars

# different levels of conservativeness
var.select(object = pbc.obj, conservative = "low")
var.select(object = pbc.obj, conservative = "medium")
var.select(object = pbc.obj, conservative = "high")

## ------------------------------------------------------------
## Minimal depth variable selection
## competing risk analysis
## ------------------------------------------------------------

## competing risk data set involving AIDS in women
data(wihs, package = "randomForestSRC")
vs.wihs <- var.select(Surv(time, status) ~ ., wihs, nsplit = 3, 
                      ntree = 100, importance = TRUE)

## competing risk analysis of pbc data from survival package
## implement cause-specific variable selection 
if (library("survival", logical.return = TRUE)) {
  data(pbc, package = "survival")
  pbc$id <- NULL
  var.select(Surv(time, status) ~ ., pbc, nsplit = 10, cause = 1)
  var.select(Surv(time, status) ~ ., pbc, nsplit = 10, cause = 2)
}

## ------------------------------------------------------------
## Minimal depth variable selection
## classification analysis
## ------------------------------------------------------------

vs.iris <- var.select(Species ~ ., iris)

## ------------------------------------------------------------
## Minimal depth variable selection
## Regression analysis
## ------------------------------------------------------------

#Variable hunting (overkill for low dimensions)
vh.air <- var.select(Ozone ~., airquality, method = "vh", nrep = 10, mvars = 5)

#better analysis
vs.air <- var.select(Ozone ~., airquality)

## ------------------------------------------------------------
## Minimal depth high-dimensional example
## van de Vijver microarray breast cancer survival data
## predefined weights for *selecting* a gene for node splitting
## determined from a preliminary forest analysis
## ------------------------------------------------------------

data(vdv, package = "randomForestSRC")
md.breast <- var.select(Surv(Time, Censoring) ~ ., vdv,
  prefit = list(action = TRUE))

## same analysis, but with customization for the preliminary forest fit
## note the large mtry and small nodesize values used
md.breast.custom <- var.select(Surv(Time, Censoring) ~ ., vdv, 
  prefit = list(action = TRUE, mtry = 500, nodesize = 1))

## ------------------------------------------------------------
## Minimal depth high-dimensional example
## van de Vijver microarray breast cancer survival data
## predefined weights for genes for *splitting* tree nodes
## weights defined in terms of cox p-values
## ------------------------------------------------------------

if (library("survival", logical.return = TRUE) 
    & library("parallel", logical.return = TRUE))
{
  cox.weights <- function(rfsrc.f, rfsrc.data) {
    event.names <- all.vars(rfsrc.f)[1:2]
    p <- ncol(rfsrc.data) - 2
    event.pt <- match(event.names, names(rfsrc.data))
    xvar.pt <- setdiff(1:ncol(rfsrc.data), event.pt)
    unlist(mclapply(1:p, function(j) {
      cox.out <- coxph(rfsrc.f, rfsrc.data[, c(event.pt, xvar.pt[j])])
      pvalue <- summary(cox.out)$coef[5]
      if (is.na(pvalue)) 1.0 else 1/(pvalue + 1e-100)
    }))
  }       
  data(vdv, package = "randomForestSRC")
  rfsrc.f <- as.formula(Surv(Time, Censoring) ~ .)
  cox.wts <- cox.weights(rfsrc.f, vdv)
  breast.obj <- rfsrc(rfsrc.f, vdv, nsplit = 10, xvar.wt = cox.wts,
                      importance = TRUE)
  md.breast.splitwt <- var.select(object = breast.obj) 
}


## ------------------------------------------------------------
## Variable hunting high-dimensional example
## van de Vijver microarray breast cancer survival data
## nrep is small for illustration; typical values are nrep = 100
## ------------------------------------------------------------

data(vdv, package = "randomForestSRC")
vh.breast <- var.select(Surv(Time, Censoring) ~ ., vdv,
      method = "vh", nrep = 10, nstep = 5)

# plot top 10 variables
plot.variable(vh.breast$rfsrc.refit.obj,
  xvar.names = vh.breast$topvars[1:10])
plot.variable(vh.breast$rfsrc.refit.obj,
  xvar.names = vh.breast$topvars[1:10], partial = TRUE)

## similar analysis, but using weights from univarate cox p-values
if (library("survival", logical.return = TRUE))
{
  cox.weights <- function(rfsrc.f, rfsrc.data) {
    event.names <- all.vars(rfsrc.f)[1:2]
    p <- ncol(rfsrc.data) - 2
    event.pt <- match(event.names, names(rfsrc.data))
    xvar.pt <- setdiff(1:ncol(rfsrc.data), event.pt)
    sapply(1:p, function(j) {
      cox.out <- coxph(rfsrc.f, rfsrc.data[, c(event.pt, xvar.pt[j])])
      pvalue <- summary(cox.out)$coef[5]
      if (is.na(pvalue)) 1.0 else 1/(pvalue + 1e-100)
    })
  }       
  data(vdv, package = "randomForestSRC")
  rfsrc.f <- as.formula(Surv(Time, Censoring) ~ .)
  cox.wts <- cox.weights(rfsrc.f, vdv)
  vh.breast.cox <- var.select(rfsrc.f, vdv, method = "vh", nstep = 5,
    nrep = 10, xvar.wt = cox.wts)
}

## ------------------------------------------------------------
## variable selection for multivariate mixed forests
## ------------------------------------------------------------

mtcars.new <- mtcars
mtcars.new$cyl <- factor(mtcars.new$cyl)
mtcars.new$carb <- factor(mtcars.new$carb, ordered = TRUE)
mv.obj <- rfsrc(cbind(carb, mpg, cyl) ~., data = mtcars.new,
              importance = TRUE)
var.select(mv.obj, method = "vh.vimp", nrep = 10)

}
}
\keyword{variable selection}

