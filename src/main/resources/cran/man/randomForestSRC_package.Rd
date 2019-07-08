\name{_PROJECT_PACKAGE_NAME_-package}
\alias{_PROJECT_PACKAGE_NAME_-package}
\docType{package}
\title{
  Fast Unified Random Forests for Survival, Regression, and Classification (RF-SRC)
}
\description{

  Fast OpenMP parallel computing for unified Breiman random forests
  (Breiman 2001) for regression, classification, survival analysis,
  competing risks, multivariate, unsupervised, quantile regression, and
  class imbalanced q-classification.  Missing data imputation includes
  missForest and multivariate missForest. New fast subsampling random
  forests. Confidence intervals for variable importance.  New holdout
  vimp.

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

    \item \code{\link{rfsrc.fast}}
    
    A fast implementation of \code{rfsrc} using subsampling.

    \item \code{\link{quantreg.rfsrc}}, \code{\link{quantreg}}

    Univariate and multivariate quantile regression forest for training
    and testing.  Different methods available including the
    Greenwald-Khanna (2001) algorithm, which is especially suitable for
    big data due to its high memory efficiency.

    \item \code{\link{predict.rfsrc}}, \code{predict}
    
    Used for prediction.  Predicted values are obtained by dropping the
    user supplied test data down the grow forest.  The resulting object
    has class \code{(rfsrc, predict)}.
    
    
    \item \code{\link{vimp.rfsrc}}, \link{subsample.rfsrc},
    \code{\link{holdout.vimp.rfsrc}}
    
    
    Used for variable selection:
    \enumerate{

      \item \code{vimp} calculates variable imporance (VIMP) from a
      RF-SRC grow/predict object by noising up the variable (for example
      by permutation).  Note that grow/predict calls can always directly
      request VIMP.

      \item \code{subsample} calculates VIMP confidence itervals via
      subsampling.

      \item \code{holdout.vimp} measures the importance of a variable
      when it is removed from the model.

    }
    
    \item \code{\link{impute.rfsrc}}
    
    Fast imputation mode for RF-SRC.  Both \code{rfsrc} and
    \code{predict.rfsrc} are capable of imputing missing data.
    However, for users whose only interest is imputing data, this function
    provides an efficient and fast interface for doing so.

    \item \code{\link{partial.rfsrc}}
    
    Used to extract the partial effects of a variable or variables on the ensembles.

  }
  
}

\section{Source Code, Beta Builds and Bug Reporting}{

  \enumerate{

    \item Regular stable releases of this package are available on CRAN
    at \url{https://cran.r-project.org/package=randomForestSRC}

    \item Interim unstable development builds with bug fixes and sometimes additional
    functionality are available at
    \url{https://github.com/kogalur/randomForestSRC}

    \item Bugs may be reported via \url{https://github.com/kogalur/randomForestSRC/issues/new}.
    Please provide the accompanying information with any reports:
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
  
}
\section{OpenMP Parallel Processing -- Installation}{This package
  implements OpenMP shared-memory parallel programming if the
  target architecture and operating system support it.  This is the
  default mode of execution.

  Additional instructions for configuring OpenMP parallel processing
  are available at
  \url{https://kogalur.github.io/randomForestSRC/building.html}.
  
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
  
  Lu M., Sadiq S., Feaster D.J. and Ishwaran H. (2018). Estimating
  individual treatment effect in observational data using random forest
  methods. \emph{J. Comp. Graph. Statist}, 27(1), 209-219

  Ishwaran H. and Lu M.  (2019).  Standard errors and confidence
  intervals for variable importance in random forest regression,
  classification, and survival. \emph{Statistics in Medicine}, 38,
  558-582.

  O'Brien R. and Ishwaran H. (2019).  A random forests quantile
  classifier for class imbalanced data. \emph{Pattern Recognition},
  90, 232-249

  Tang F. and Ishwaran H. (2017).  Random forest missing data
  algorithms.  \emph{Statistical Analysis and Data Mining},
  10, 363-377.
  
}
\seealso{
  \command{\link{find.interaction.rfsrc}},

  \command{\link{holdout.vimp.rfsrc}},

  \command{\link{imbalanced.rfsrc}},
  \command{\link{impute.rfsrc}},

  \command{\link{max.subtree.rfsrc}},

  \command{\link{partial.rfsrc}},
  \command{\link{plot.competing.risk.rfsrc}},
  \command{\link{plot.rfsrc}},
  \command{\link{plot.survival.rfsrc}},
  \command{\link{plot.variable.rfsrc}},
  \command{\link{predict.rfsrc}},
  \command{\link{print.rfsrc}},

  \command{\link{quantreg.rfsrc}},

  \command{\link{rfsrc}},
  \command{\link{rfsrc.cart}},
  \command{\link{rfsrc.fast}},

   

  \command{\link{stat.split.rfsrc}},
  \command{\link{subsample.rfsrc}},
  \command{\link{synthetic.rfsrc}},
  
  \command{\link{tune.rfsrc}},

  \command{\link{var.select.rfsrc}},
  \command{\link{vimp.rfsrc}}
}
\keyword{package}
