\name{get.tree.rfsrc}
\alias{get.tree.rfsrc}
\alias{get.tree}
\title{Extract a Single Tree from a Forest and plot it on your browser}
\description{
  Extracts a single tree from a forest which can then be plotted on the users
  browser.  Works for all families.  Missing data not permitted.}
\usage{\method{get.tree}{rfsrc}(object, tree.id, target,
  m.target = NULL, time, surv.type = c("mort", "rel.freq",
  "surv", "years.lost", "cif", "chf"), class.type =
  c("prob", "bayes"), oob = TRUE, show.plots = TRUE)
}
\arguments{
  \item{object}{An object of class \code{(rfsrc, grow)}.}
  \item{tree.id}{Integer value specifying the tree to be extracted.}
  \item{target}{For classification, an integer or
    character value specifying the class to focus on (defaults to the
    first class).  For competing risks, an integer value between
    1 and \code{J} indicating the event of interest, where \code{J} is
    the number of event types.  The default is to use the first event
    type.}
  \item{m.target}{Character value for multivariate families
    specifying the target outcome to be used.  If left unspecified, the
    algorithm will choose a default target.}
   \item{time}{For survival, the time at which the predicted
    survival value is evaluated at (depends on \code{surv.type}).}
  \item{surv.type}{For survival, specifies the predicted value.
    See details below.}
  \item{class.type}{For classification, specifies the predicted value.
    See details below.}
  \item{oob}{OOB (TRUE) or in-bag (FALSE) predicted values.}
   \item{show.plots}{Should plots be displayed?}
}
\details{

  Extracts a specified tree from a forest and converts the tree to a
  hierarchical structure suitable for use with the "data.tree" package.
  Plotting the object will conveniently render the tree on the users
  browser. Left tree splits are displayed.  For continuous values, left
  split is displayed as an inequality with right split equal to the
  reversed inequality.  For factors, split values are described in terms
  of the levels of the factor. In this case, the left daughter split is
  a set consisting of all levels that are assigned to the left daughter
  node.  The right daughter split is the complement of this set.

  Terminal nodes are highlighted by color and display the sample size
  and predicted value.  Here predicted value equals the forest
  ensemble value and not the tree predicted value.  This is done as it
  allows one to visualize the ensemble predictor over a given tree and
  therefore a given partition of the feature space.  The predicted value
  displayed is chosen according to the following rules and options:

  
  \enumerate{
    \item For regression, the predicted response is used.
    
    \item For classification, it is the predicted class probability
    specified by \option{target}, or the class of maximum probability
    depending on whether \option{class.type} is set to "prob" or
    "bayes".
    
    
    \item For multivariate families, it is the predicted value of the
    outcome specified by \option{m.target} and if that is a
    classification outcome, by \option{target}.
    
    \item For survival, the choices are:
    \itemize{
	    \item  Mortality (\code{mort}).
	    \item  Relative frequency of mortality (\code{rel.freq}).
	    \item Predicted survival (\code{surv}), where the predicted
      survival is for the time point specified using
      \code{time} (the default is the median follow up time).
    }
    \item For competing risks, the choices are:
    \itemize{
	    \item  The expected number of life years lost (\code{years.lost}).
	    \item  The cumulative incidence function (\code{cif}).
	    \item  The cumulative hazard function (\code{chf}).
    }
    In all three cases, the predicted value is for the event type
    specified by \option{target}.  For \code{cif} and
    \code{chf} the quantity is evaluated at the time point specified
    by \code{time}.
  }
  
}
\value{
  Invisibly, returns an object with hierarchical structure formatted for use with
  the data.tree package.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur

  Many thanks to @dbarg1 on GitHub for the initial prototype of this
  function
  
  
}
\examples{
\donttest{
## ------------------------------------------------------------
## survival/competing risk
## ------------------------------------------------------------

## survival - veteran data set but with factors
## note that diagtime has many levels
data(veteran, package = "_PROJECT_PACKAGE_NAME_")
vd <- veteran
vd$celltype=factor(vd$celltype)
vd$diagtime=factor(vd$diagtime)
vd.obj <- rfsrc(Surv(time,status)~., vd, ntree = 100, nodesize = 5)
plot(get.tree(vd.obj, 3))

## competing risks
data(follic, package = "_PROJECT_PACKAGE_NAME_")
follic.obj <- rfsrc(Surv(time, status) ~ ., follic, nsplit = 3, ntree = 100)
plot(get.tree(follic.obj, 2))

## ------------------------------------------------------------
## regression
## ------------------------------------------------------------

airq.obj <- rfsrc(Ozone ~ ., data = airquality)
plot(get.tree(airq.obj, 10))

## ------------------------------------------------------------
## classification
## ------------------------------------------------------------

iris.obj <- rfsrc(Species ~., data = iris)
plot(get.tree(iris.obj, 25))
plot(get.tree(iris.obj, 25, class.type = "bayes"))
plot(get.tree(iris.obj, 25, target = "setosa"))
plot(get.tree(iris.obj, 25, target = "versicolor"))
plot(get.tree(iris.obj, 25, target = "virginica"))


## ------------------------------------------------------------
## multivariate regression
## ------------------------------------------------------------

mtcars.mreg <- rfsrc(Multivar(mpg, cyl) ~., data = mtcars)
plot(get.tree(mtcars.mreg, 10, m.target = "mpg"))
plot(get.tree(mtcars.mreg, 10, m.target = "cyl"))


## ------------------------------------------------------------
## multivariate mixed outcomes
## ------------------------------------------------------------

mtcars2 <- mtcars
mtcars2$carb <- factor(mtcars2$carb)
mtcars2$cyl <- factor(mtcars2$cyl)
mtcars.mix <- rfsrc(Multivar(carb, mpg, cyl) ~ ., data = mtcars2)
plot(get.tree(mtcars.mix, 5, m.target = "cyl"))
plot(get.tree(mtcars.mix, 5, m.target = "carb"))

## ------------------------------------------------------------
## unsupervised analysis
## ------------------------------------------------------------

mtcars.unspv <- rfsrc(data = mtcars)
plot(get.tree(mtcars.unspv, 5))



}
}
\keyword{plot}
