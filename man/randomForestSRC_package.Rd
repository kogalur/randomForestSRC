\name{randomForestSRC-package}
\alias{randomForestSRC-package}
\docType{package}

\title{
  Fast Unified Random Forests for Survival, Regression, and Classification (RF-SRC)
}

\description{

Fast OpenMP-parallel implementation of Breiman's random forests
(Breiman, 2001) for regression, classification, survival analysis
(Ishwaran, 2008), competing risks (Ishwaran, 2012), multivariate
outcomes (Segal and Xiao, 2011), unsupervised learning (Mantero and
Ishwaran, 2020), quantile regression (Meinshausen, 2006; Zhang et al.,
2019; Greenwald and Khanna, 2001), and imbalanced q-classification
(O'Brien and Ishwaran, 2019).


Supports deterministic and randomized splitting rules (Geurts et al.,
2006; Ishwaran, 2015) across all families. Variable importance (VIMP),
holdout VIMP, and confidence regions (Ishwaran and Lu, 2019) can be
computed for single and grouped variables. Includes minimal depth
variable selection (Ishwaran et al., 2010, 2011) and a fast interface
for missing data imputation using multiple forest-based methods (Tang
and Ishwaran, 2017).


Tree structures can be visualized in Safari or Chrome for any family;
see \code{\link{get.tree}}.

}

\section{Package Overview}{

This package contains many useful functions. Users are encouraged to read the help files in full for detailed guidance. Below is a brief overview of key functions to help navigate the package.

\enumerate{

  \item \code{\link{rfsrc}}

  The main entry point to the package. Builds a random forest using user-supplied training data. The returned object is of class \code{(rfsrc, grow)}.

  \item \code{\link{rfsrc.fast}}

  A computationally efficient version of \code{rfsrc} using subsampling.

  \item \code{\link{quantreg.rfsrc}}, \code{\link{quantreg}}

  Univariate and multivariate quantile regression forests for training and testing. Includes methods such as the Greenwald-Khanna (2001) algorithm, ideal for large data due to its memory efficiency.

  \item \code{\link{predict.rfsrc}}, \code{predict}

  Predicts outcomes by dropping test data down the trained forest. Returns an object of class \code{(rfsrc, predict)}.

  \item \code{\link{sidClustering.rfsrc}}, \code{sidClustering}

  Unsupervised clustering using SID (Staggered Interaction Data). Also includes Breiman's artificial two-class method (Breiman, 2003).

  \item \code{\link{vimp}}, \code{\link{subsample}}, \code{\link{holdout.vimp}}

  Functions for variable selection and importance assessment:
  \enumerate{
    \item \code{vimp}: Computes variable importance (VIMP) by perturbing each variable (e.g., via permutation). Can also be computed directly in \code{rfsrc} and \code{predict.rfsrc}.
    
    \item \code{subsample}: Computes confidence intervals for VIMP using subsampling.
    
    \item \code{holdout.vimp}: Measures the effect of removing a variable from the model.
    
    \item \strong{VarPro} (\pkg{VarPro} package): For advanced model-independent variable selection using rule-based variable priority. Supports regression, classification, survival, and unsupervised data. See \url{https://www.varprotools.org}.
  }

  \item \code{\link{imbalanced.rfsrc}}, \code{\link{imbalanced}}

  Implements q-classification and G-mean-based VIMP for class-imbalanced data.

  \item \code{\link{impute.rfsrc}}, \code{\link{impute}}

  A fast interface for missing data imputation. While \code{rfsrc} and \code{predict.rfsrc} can handle missing data internally, this provides a dedicated, efficient solution for imputation tasks.

  \item \code{\link{partial.rfsrc}}, \code{\link{partial}}

  Computes partial dependence functions to assess the marginal effect of one or more variables on the forest ensemble.

}
}

\section{Home page, Vignettes, Discussions, Bug Reporting, Source Code, Beta Builds}{

\enumerate{

  \item The package home page, with vignettes, manuals, GitHub links, and additional documentation, is available at:
  \url{https://www.randomforestsrc.org/index.html}

  \item Questions, comments, and general usage discussions (non-bug-related) can be posted at:
  \url{https://github.com/kogalur/randomForestSRC/discussions/}

  \item Bug reports should be submitted at:
  \url{https://github.com/kogalur/randomForestSRC/issues/}

  Please use this only for bugs, and include the following with your report:
  \itemize{
    \item Output from \code{sessionInfo()}.
    \item A minimal reproducible example including:
      \itemize{
        \item A minimal dataset required to reproduce the error.
        \item The smallest runnable code needed to reproduce the issue.
        \item Version details of R and all relevant packages.
        \item A random seed (via \code{set.seed()}) if randomness is involved.
      }
  }

  \item The latest stable release of the package is available on CRAN:
  \url{https://cran.r-project.org/package=randomForestSRC/}

  \item Development builds (unstable) with bug fixes and new features are hosted on GitHub:
  \url{https://github.com/kogalur/randomForestSRC/}
}
}

\section{OpenMP Parallel Processing -- Installation}{

This package supports OpenMP shared-memory parallel programming on systems where the architecture and operating system permit it. OpenMP is enabled by default.

Detailed instructions for configuring OpenMP parallel processing can be found at:
\url{https://www.randomforestsrc.org/articles/installation.html}

Note that running the package with OpenMP (or Open MPI) may increase memory (RAM) usage. Users are advised to understand their system's hardware limits and to monitor resource consumption to avoid overtaxing CPU and memory capacity.
}

\section{Reproducibility}{

Model reproducibility is determined by three components: the random seed, the forest topology (i.e., the structure of trees), and terminal node membership for the training data. These elements together allow the model and its terminal node statistics to be faithfully restored.

Other outputs, such as variable importance (VIMP) and performance metrics, rely on additional internal randomization and are not considered part of the model definition. As a result, such statistics are subject to Monte Carlo variability and may differ across runs, even with the same seed.
}


\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.

  Geurts, P., Ernst, D. and Wehenkel, L., (2006). Extremely randomized
  trees. \emph{Machine learning}, 63(1):3-42.

  Greenwald M. and Khanna S. (2001).  Space-efficient online computation of
  quantile summaries. \emph{Proceedings of ACM SIGMOD}, 30(2):58-66.

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
  
  Ishwaran H. and Lu M.  (2019).  Standard errors and confidence
  intervals for variable importance in random forest regression,
  classification, and survival. \emph{Statistics in Medicine}, 38,
  558-582.

  Lu M., Sadiq S., Feaster D.J. and Ishwaran H. (2018). Estimating
  individual treatment effect in observational data using random forest
  methods. \emph{J. Comp. Graph. Statist}, 27(1), 209-219

  Mantero A. and Ishwaran H. (2021).  Unsupervised random forests.
  \emph{Statistical Analysis and Data Mining}, 14(2):144-167.

  Meinshausen N. (2006) Quantile regression forests, \emph{Journal of
  Machine Learning Research}, 7:983-999.

  O'Brien R. and Ishwaran H. (2019).  A random forests quantile
  classifier for class imbalanced data. \emph{Pattern Recognition},
  90, 232-249

  Segal M.R. and Xiao Y. Multivariate random
  forests. (2011). \emph{Wiley Interdisciplinary Reviews: Data Mining
  and Knowledge Discovery}. 1(1):80-87.

  Tang F. and Ishwaran H. (2017).  Random forest missing data
  algorithms.  \emph{Statistical Analysis and Data Mining}, 10:363-377.

  Zhang H., Zimmerman J., Nettleton D. and Nordman D.J. (2019).  Random
  forest prediction intervals. \emph{The American Statistician}. 4:1-5.

  
}
\seealso{
  \command{\link{find.interaction.rfsrc}},

  \command{\link{get.tree.rfsrc}},

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

  \command{\link{sidClustering.rfsrc}},

  \command{\link{stat.split.rfsrc}},
  \command{\link{subsample.rfsrc}},
  \command{\link{synthetic.rfsrc}},
  
  \command{\link{tune.rfsrc}},

  \command{\link{vimp.rfsrc}}
}
\keyword{package}
