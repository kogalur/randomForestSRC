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

  \item{tree.id}{Integer specifying the tree to extract.}

  \item{target}{For classification: integer or character indicating the class of interest (defaults to the first class). For competing risks: integer between 1 and \code{J} (number of event types) specifying the event of interest (default is the first event type).}

  \item{m.target}{Character string specifying the target outcome for multivariate families. If unspecified, a default is selected.}

  \item{time}{For survival: time point at which the predicted value is evaluated (depends on \code{surv.type}).}

  \item{surv.type}{For survival: specifies the type of predicted value returned. See \code{Details}.}

  \item{class.type}{For classification: specifies the type of predicted value. See \code{Details}.}

  \item{ensemble}{Logical. If \code{TRUE}, prediction is based on the ensemble of all trees. If \code{FALSE} (default), prediction is based on the specified tree.}

  \item{oob}{Logical. Use OOB predicted values (\code{TRUE}) or in-bag values (\code{FALSE}). Only applies when \code{ensemble=TRUE}.}

  \item{show.plots}{Logical. Should plots be displayed?}

  \item{do.trace}{Number of seconds between progress updates.}
}

\details{

  Extracts a specified tree from a forest and converts it into a
  hierarchical structure compatible with the \pkg{data.tree}
  package. Plotting the resulting object renders an interactive tree
  visualization in the user's web browser.
  

  Left-hand splits are shown. For continuous variables, the left split
  is displayed as an inequality (e.g., \code{x < value}); the right
  split is the reverse. For factor variables, the left daughter node is
  defined by a set of levels assigned to it; the right daughter is its
  complement.
  

  Terminal nodes are highlighted with color and display both sample size
  and predicted value. By default, the predicted value corresponds to
  the prediction from the selected tree, and the sample size refers to
  the in-bag cases reaching the terminal node. If \code{ensemble =
  TRUE}, the predicted value equals the forest ensemble prediction,
  allowing visualization of the full forest predictor over the selected
  tree's partition. In this case, sample sizes refer to all observations
  (not just in-bag cases).

  Predicted values displayed in terminal nodes are defined as follows:

  \enumerate{
    \item For regression: the mean of the response.

    \item For classification: depends on the \code{class.type} argument and target class:
      \itemize{
        \item If \code{class.type = "bayes"}, the predicted class with the most votes, or the RFQ classifier threshold in two-class problems.
        \item If \code{class.type = "prob"}, the class probability for the target class.
      }

    \item For multivariate families: the predicted value for the outcome specified by \code{m.target}, using the logic above depending on whether the outcome is continuous or categorical.

    \item For survival:
      \itemize{
        \item \code{mort}: estimated mortality (Ishwaran et al., 2008).
        \item \code{rel.freq}: relative frequency of mortality.
        \item \code{surv}: predicted survival probability at the specified time (\code{time}).
      }

    \item For competing risks:
      \itemize{
        \item \code{years.lost}: expected number of life years lost.
        \item \code{cif}: cumulative incidence function.
        \item \code{chf}: cause-specific cumulative hazard function.
      }
      For \code{cif} and \code{chf}, predictions are evaluated at the time point given by \code{time}, and all metrics are specific to the event type indicated by \code{target}.
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
