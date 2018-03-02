\name{_PROJECT_PACKAGE_NAME_-package}
\alias{_PROJECT_PACKAGE_NAME_-package}
\alias{cindex}
\alias{extract.bootsample}
\alias{extract.subsample}
\alias{get.mv.error}
\alias{get.mv.predicted}
\alias{get.mv.vimp}
\alias{perf.metric}
\docType{package}
\title{
  Random Forests for Survival, Regression, and Classification (RF-SRC)
}
\description{
  This package provides a unified treatment of Breiman's random forests
  (Breiman 2001) for a variety of data settings.  Regression and
  classification forests are grown when the response is numeric or
  categorical (factor), while survival and competing risk forests
  (Ishwaran et al. 2008, 2012) are grown for right-censored survival
  data.  Multivariate regression and classification responses as well as
  mixed outcomes (regression/classification responses) are also handled.
  Also includes unsupervised forests and quantile regression forests.
  Different splitting rules invoked under deterministic or random
  splitting are available for all families.  Variable predictiveness can
  be assessed using variable importance (VIMP) measures for single, as
  well as grouped variables.  Variable selection is implemented using
  minimal depth variable selection (Ishwaran et al. 2010). Missing data
  (for x-variables and y-outcomes) can be imputed on both training and
  test data.  The underlying code is based on Ishwaran and Kogalur's now
  retired \pkg{randomSurvivalForest} package (Ishwaran and Kogalur
  2007), and has been significantly refactored for improved
  computational speed.
}

\section{Package Overview}{
  
  This package contains many useful functions and users should read the
  help file in its entirety for details.  However, we briefly mention
  several key functions that may make it easier to navigate and
  understand the layout of the package.
  
  \enumerate{
    \item \code{\link{rfsrc}}
    
    This is the main entry point to the package.  It grows a random forest
    using user supplied training data.  We refer to the resulting object
    as a RF-SRC grow object.  Formally, the resulting object has class
    \code{(rfsrc, grow)}.
    
    \item \code{\link{predict.rfsrc}}, \code{predict}
    
    Used for prediction.  Predicted values are obtained by dropping the
    user supplied test data down the grow forest.  The resulting object
    has class \code{(rfsrc, predict)}.
    
    
    \item \code{\link{max.subtree}}, \code{\link{var.select}}
    
    Used for variable selection.  The function \code{max.subtree}
    extracts maximal subtree information from a RF-SRC object which is
    used for selecting variables by making use of minimal depth variable
    selection.  The function \code{var.select} provides
    an extensive set of variable selection options and is a wrapper to
    \code{max.subtree}.
    
    \item \code{\link{impute.rfsrc}}
    
    Fast imputation mode for RF-SRC.  Both \code{rfsrc} and
    \code{predict.rfsrc} are capable of imputing missing data.
    However, for users whose only interest is imputing data, this function
    provides an efficient and fast interface for doing so.

    \item \code{\link{partial.rfsrc}}
    
    Used to extract the partial effects of a variable or variables on the ensembles.

  }
}

\section{Source Code, Beta Builds and Bug Reporting}{Regular stable releases of this package
  are available on CRAN at
  \url{_HTTPS_PREFIX_cran.r-project.org/package=randomForestSRC}.
  Interim unstable development builds with bug fixes and sometimes additional
  functionality are available at
  \url{_HTTPS_PREFIX_github.com/kogalur/randomForestSRC}.
  Bugs may be reported via \url{_HTTPS_PREFIX_github.com/kogalur/randomForestSRC/issues/new}.
  Please provide the accompanying
  information with any reports:

  \enumerate{
    \item \code{sessionInfo()}
    \item A minimal reproducible example consisting of the following items:
    \itemize{
      \item a minimal dataset, necessary to reproduce the error
      \item the minimal runnable code necessary to reproduce the error, 
      which can be run on the given dataset
      \item the necessary information on the used packages, R version and system it is run on
      \item in the case of random processes, a seed (set by \code{set.seed()}) for reproducibility
    }
  }
}

\section{OpenMP Parallel Processing -- Installation}{This package
  implements OpenMP shared-memory parallel programming if the
  target architecture and operating system support it.  This is the
  default mode of execution.

  Additional instructions for configuring OpenMP parallel processing
  are available at
  \url{_HTTPS_PREFIX_kogalur.github.io/randomForestSRC/building.html}.
  
  An understanding of resource utilization (CPU and RAM) is necessary
  when running the package using OpenMP and Open MPI parallel execution.
  Memory usage is greater when running with OpenMP
  enabled.  Diligence should be used not to overtax the hardware available.

}

\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.
  
  Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7(2):25-31.
  
  Ishwaran H. (2007).  Variable importance in binary regression
  trees and forests,  \emph{Electronic J. Statist.}, 1:519-537.

  Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
  (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.
  
  Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and
  Lauer M.S. (2010).  High-dimensional variable selection for survival
  data.  \emph{J. Amer. Statist. Assoc.}, 105:205-217.

  Ishwaran H., Kogalur U.B., Chen X. and Minn A.J. (2011). Random survival
  forests for high-dimensional data. \emph{Stat. Anal. Data Mining}, 4:115-132

  Ishwaran H., Gerds T.A., Kogalur U.B., Moore R.D., Gange S.J. and Lau
  B.M. (2014). Random survival forests for competing risks.
  \emph{Biostatistics}, 15(4):757-773.

  Ishwaran H. and Malley J.D. (2014). Synthetic learning
  machines. \emph{BioData Mining}, 7:28.
  
  Ishwaran H. (2015).  The effect of splitting on random forests.
  \emph{Machine Learning}, 99:75-118.
  
  Ishwaran H. and Lu M.  (2017).  Standard errors and confidence
  intervals for variable importance in random forest regression,
  classification, and survival.

  Mantero A. and Ishwaran H. (2017). Unsupervised random forests.

  O'Brien R. and Ishwaran H. (2017).  A random forests quantile
  classifier for class imbalanced data.

  Tang F. and Ishwaran H. (2017).  Random forest missing data
  algorithms.  \emph{Statistical Analysis and Data Mining}, 10, 363-377.
}
\seealso{
  \command{\link{find.interaction}},
   
  \command{\link{impute}},
  \command{\link{max.subtree}},
   
  \command{\link{plot.competing.risk}},
  \command{\link{plot.rfsrc}},
  \command{\link{plot.survival}},
  \command{\link{plot.variable}},
  \command{\link{predict.rfsrc}},
  \command{\link{print.rfsrc}},
  \command{\link{quantileReg}},
  \command{\link{rfsrc}},
  \command{\link{rfsrcSyn}},
   
  \command{\link{stat.split}},
  \command{\link{var.select}},
  \command{\link{vimp}}
}
\keyword{package}
