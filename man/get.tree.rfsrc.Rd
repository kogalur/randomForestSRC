\name{get.tree.rfsrc}
\alias{get.tree.rfsrc}
\alias{get.tree}
\title{Extract a Single Tree from a Forest and plot it on your browser}
\description{
  Extracts a single tree from a forest which can then be plotted on the users
  browser.  Works for all families.  Missing data not permitted.}
\usage{\method{get.tree}{rfsrc}(object, tree.id, target, m.target = NULL,
   time, surv.type = c("mort", "rel.freq", "surv", "years.lost", "cif", "chf"),
   class.type = c("bayes", "rfq", "prob"),
   ensemble = FALSE, oob = TRUE, show.plots = TRUE, do.trace = FALSE)
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
  \item{ensemble}{Use the ensemble (of all trees) for prediction, or use
    the requested tree for prediction (this is the default).}
  \item{oob}{OOB (TRUE) or in-bag (FALSE) predicted values.  Only
    applies when \code{ensemble=TRUE}.}
  \item{show.plots}{Should plots be displayed?}
  \item{do.trace}{Number of seconds between updates to the user on
    approximate time to completion.}

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
  and predicted value.  By default, predicted value equals the tree
  predicted value and sample size are terminal node inbag sample sizes.
  If \code{ensemble=TRUE}, then the predicted value equals the forest
  ensemble value which could be useful as it allows one to visualize the
  ensemble predictor over a given tree and therefore for a given
  partition of the feature space.  In this case, sample sizes are for
  all cases and not the tree specific inbag cases.

  The predicted value displayed is as follows:

  
  \enumerate{
    \item For regression, the mean of the response.
    
    \item For classification, for the target class specified by
    \option{target}, either the class with most votes if
    \code{class.type="bayes"}; or in a two-class problem the classifier
    using the RFQ quantile threshold if \code{class.type="bayes"}
    (see \link{imbalanced} for more details); or the relative class
    frequency when \code{class.type="prob"}.
    
    \item For multivariate families, the predicted value of the outcome
    specified by \option{m.target}.  This being the value for regression
    or classification described above, depending on whether the outcome
    is real valued or a factor.
    
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
data(veteran, package = "randomForestSRC")
vd <- veteran
vd$celltype=factor(vd$celltype)
vd$diagtime=factor(vd$diagtime)
vd.obj <- rfsrc(Surv(time,status)~., vd, ntree = 100, nodesize = 5)
plot(get.tree(vd.obj, 3))

## competing risks
data(follic, package = "randomForestSRC")
follic.obj <- rfsrc(Surv(time, status) ~ ., follic, nsplit = 3, ntree = 100)
plot(get.tree(follic.obj, 2))

## ------------------------------------------------------------
## regression
## ------------------------------------------------------------

airq.obj <- rfsrc(Ozone ~ ., data = airquality)
plot(get.tree(airq.obj, 10))

## ------------------------------------------------------------
## two-class imbalanced data (see imbalanced function)
## ------------------------------------------------------------

data(breast, package = "randomForestSRC")
breast <- na.omit(breast)
f <- as.formula(status ~ .)
breast.obj <- imbalanced(f, breast)

## compare RFQ to Bayes Rule
plot(get.tree(breast.obj, 1, class.type = "rfq", ensemble = TRUE))
plot(get.tree(breast.obj, 1, class.type = "bayes", ensemble = TRUE))

## ------------------------------------------------------------
## classification
## ------------------------------------------------------------

iris.obj <- rfsrc(Species ~., data = iris, nodesize = 10)

## equivalent
plot(get.tree(iris.obj, 25))
plot(get.tree(iris.obj, 25, class.type = "bayes"))

## predicted probability displayed for terminal nodes
plot(get.tree(iris.obj, 25, class.type = "prob", target = "setosa"))
plot(get.tree(iris.obj, 25, class.type = "prob", target = "versicolor"))
plot(get.tree(iris.obj, 25, class.type = "prob", target = "virginica"))


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
