\name{randomForestSRC-package}
\alias{randomForestSRC-package}
\docType{package}
\title{
  Random Forests for Survival, Regression and Classification (RF-SRC)
}
\description{
  This package provides a unified treatment of Breiman's random forests
  (Breiman 2001) for a variety of data settings.  Regression and
  classification forests are grown when the response is numeric or
  categorical (factor), while survival and competing risk forests
  (Ishwaran et al. 2008, 2012) are grown for right-censored survival
  data.  Multivariate regression and classification responses as well as
  mixed outcomes (regression/classification responses) are also handled
  as are unsupervised forests.  Different splitting rules invoked under
  deterministic or random splitting are available for all families.
  Variable predictiveness can be assessed using variable importance
  (VIMP) measures for single, as well as grouped variables.  Variable
  selection is implemented using minimal depth variable selection
  (Ishwaran et al. 2010). Missing data (for x-variables and y-outcomes)
  can be imputed on both training and test data.  The underlying code is
  based on Ishwaran and Kogalur's now retired \pkg{randomSurvivalForest}
  package (Ishwaran and Kogalur 2007), and has been significantly
  refactored for improved computational speed.
}

\section{OpenMP Parallel Processing -- Installation}{This package
  implements OpenMP shared-memory parallel programming.  However, the
  default installation will only execute serially.  To utilize OpenMP,
  the target architecture and operating system must support it.

  There are THREE strategies for installing the OpenMP capable version
  of the CRAN package.  \bold{Method 1} relies on having full package
  development prerequisites, the build environment necessary for
  creating R packages from source.  This is the preferred and
  comprehensive way to guarantee natively compiled, compatible and
  optimized binaries for your system.  \bold{Method 2} relies on having
  partial package development prerequisites.  This makes some
  assumptions about your system, that are not entirely platform
  independent, but will usually work.  \bold{Method 3} is the easiest
  and does not require package development prerequisites.  Method 3
  relies on pre-built binaries, and is intended for users not interested
  in investing the time necessary to have the capability to natively
  build packages. We do not recommend this method as it does not
  guarantee OpenMP execution, nor does it guarantee that our binaries
  will even be compatible with your system.  However, we are providing
  them as a convenience.
  
  \tabular{cccc}{
    \command{\bold{METHOD |}} \tab \command{\bold{| R Development Toolset |}} 
    \tab \command{\bold{| Difficulty |}} 
    \tab \command{\bold{| OpenMP Execution}} \cr
    \command{1} \tab \command{Full} \tab \command{High} \tab \command{Guaranteed} \cr
    \command{2} \tab \command{Partial} \tab \command{Medium} \tab
    \command{Good Success} \cr
    \command{3} \tab \command{None} \tab \command{Low} \tab \command{Moderate Success} \cr
  }

}

\section{OpenMP Parallel Processing -- Method 1}{

  The core software development utilities required for R package
  development vary by operating system.  The difficulty of installing
  this build environment also varies by operating system.  Unix-based
  systems are the friendliest, followed by Mac OS X, followed lastly by
  Windows.  Detailed descriptions on how this is achieved are available
  on a number of sites online and will not be reproduced here.  Once the
  R package development environment is in place, it is possible to build
  our package natively on your platform using the following steps:
  
  \enumerate{
    \item Download the package source code
    \command{randomForestSRC_X.x.x.tar.gz} from CRAN at
    \url{https://cran.r-project.org/package=randomForestSRC}.  The
    X's indicate the version posted. Do not download the binary!
    
    \item Open a console, navigate to the directory containing the
    tarball, and untar it using the command

    \code{tar -xvf randomForestSRC_X.x.x.tar.gz}
    
    \item This will create a directory structure with the root directory
    of the package named \command{randomForestSRC}.  Change into the root
    directory of the package using the command

    \code{cd randomForestSRC}
    
    \item Run autoconf using the command

    \code{autoconf}
    
    \item Change back to your working directory using the command

    \code{cd ..}
    
    \item From your working directory, execute the command

    \code{R CMD INSTALL --preclean --clean randomForestSRC}

    on the modified
    package.  Ensure that you do not target the unmodified tarball, but
    instead act on the directory structure you just modified.
  }
}

\section{OpenMP Parallel Processing -- Method 2}{
  This method hard codes some OpenMP compiler directives.  On Windows systems
  this method is generally easier than Method 1.  The
  instruction for this method follow:
  
  \enumerate{
    \item Download the package source code \command{randomForestSRC_X.x.x.tar.gz}
    from CRAN at \url{https://cran.r-project.org/package=randomForestSRC}.  The
    X's indicate the version posted. Do not
    download the binary!
    
    \item Open a console, navigate to the directory containing the
    tarball, and untar it using the command

    \code{tar -xvf randomForestSRC_X.x.x.tar.gz}
    
    \item This will create a directory structure with the root directory
    of the package named \command{randomForestSRC}.

    \item Download the Makevars file containing the custom compiler
    directives from
    \url{http://www.ccs.miami.edu/~hishwaran/rfsrc/Makevars}.  Copy it
    into the directory \command{randomForestSRC/src}. On Windows
    systems, take the additional step of renaming it to
    \command{Makevars.win}.
    
    \item From your working directory, execute the command

    \code{R CMD INSTALL --preclean --clean randomForestSRC}

    on the modified package.  Ensure that you do not target the
    unmodified tarball, but instead act on the directory structure you
    just modified.
  }
}

\section{OpenMP Parallel Processing -- Method 3}{
  The easiest but least recommended way to install the OpenMP version of
  the package is to download our pre-compiled binaries and attempt to
  use them on your system. Successful execution is not always
  guaranteed.  The instruction for this method follow:
  
  \enumerate{
    \item Download the platform specific binary file
    \command{randomForestSRC_X.x.x.<extension>} from
    \url{http://www.ccs.miami.edu/~hishwaran/rfsrc.html}.
    
    \item From the R GUI, navigate to "Install Packages" (or similar)
    and select "From Local Archive" (or similar).  Then navigate to the
    downloaded binary and click "Install".  If you are not using the R
    GUI, navigate to the directory where you downloaded the compressed
    binary file and execute the command

    \code{R CMD INSTALL randomForestSRC_X.x.x.zip}

  }
}

\section{OpenMP Parallel Processing -- Setting the Number of CPUs}{
  There are several ways to control the number of CPU cores that the
  package accesses during OpenMP parallel execution.  First, you will
  need to determine the number of cores on your local machine.  Do this
  by starting an R session and issuing the command
  \command{detectCores()}.  You will require the \pkg{parallel} package
  for this.
  
  Then you can do the following:
  
  At the start of every R session, you can set the number of cores
  accessed during OpenMP parallel execution by issuing the command
  \command{options(rf.cores = x)}, where \command{x} is the number of
  cores.  If \command{x} is a negative number, the package will access
  the maximum number of cores on your machine.  The options command can
  also be placed in the users .Rprofile file for convenience.  You can,
  alternatively, initialize the environment variable \command{RF_CORES}
  in your shell environment.
  
  The default value for rf.cores is -1 (-1L), if left unspecified, which
  uses all available cores, with a minimum of two.
}
\section{R-side Parallel Processing -- Setting the Number of CPUs}{
  The package also implements R-side parallel processing by replacing
  the R function \command{lapply} with \command{mclapply} found in the
  \pkg{parallel} package.  You can set the number of cores accessed by
  \command{mclapply} by issuing the command

  \code{options(mc.cores = x)}

  where \command{x} is the number of cores.  The options command
  can also be placed in the users .Rprofile file for convenience.  You
  can, alternatively, initialize the environment variable
  \command{MC_CORES} in your shell environment.  See the help files in
  \pkg{parallel} for more information.
  
  The default value for \command{mclapply} on non-Windows systems is
  two (2L) cores.  On Windows systems, the default value is one (1L)
  core.
}
\section{Example: Setting the Number of CPUs}{
  As an example, issuing the following options command uses all
  available cores for both OpenMP and R-side processing:
  
  \command{options(rf.cores=detectCores(), mc.cores=detectCores())}
  
  As stated above, this option command can be placed in the users
  .Rprofile file.
}
\section{CAUTIONARY NOTE}{
  
  Regarding C-side threading (accessed via OpenMP compilation) versus
  R-side forking (accessed via \command{mclapply} in package
  \pkg{parallel}).
  
  \enumerate{
    \item Once the package has been compiled with OpenMP enabled, trees
    will be grown in parallel using the \command{rf.cores} option.
    Independently of this, we also utilize \command{mclapply} to
    parallelize loops in R-side pre-processing and post-processing
    of the forest.  This is always available and independent of
    whether the user chooses to compile the package with the OpenMP
    option enabled.
    
    \item It is important NOT to write programs that fork R processes
    containing OpenMP threads.  That is, one should not use
    \command{mclapply} around the functions \command{rfsrc},
    \command{predict.rfsrc}, \command{vimp.rfsc},
    \command{var.select.rfsrc}, and
    \command{find.interaction.rfsrc}.  In such a scenario, program
    execution is not guaranteed.
    
    \item Note that \command{options(rf.cores=0)} disables C-side
    threading, and \command{options(mc.cores=1)} disables R-side
    forking.  Therefore, setting \command{options(rf.cores=0)}, is
    one means to wrap  \command{mclapply} around the functions
	  listed above in 2.
  }
}

\section{Beta Builds and Bug Reporting}{Regular releases of this package
  are available on CRAN at
  \url{https://cran.r-project.org/package=randomForestSRC}.
  Interim development builds with bug fixes and sometimes additional
  functionality are available at
  \url{https://github.com/kogalur/randomForestSRC}.
  Bugs may be reported via this page.  Please provide the accompanying
  information with any reports:

  \enumerate{
    \item \command{sessionInfo()}
    \item A minimal reproducible example consisting of the following items:
    \enumerate{
      \item a minimal dataset, necessary to reproduce the error
      \item the minimal runnable code necessary to reproduce the error, 
      which can be run on the given dataset
      \item the necessary information on the used packages, R version and system it is run on
      \item in the case of random processes, a seed (set by set.seed()) for reproducibility
    }
  }
}

\section{Package Overview}{
  
  This package contains many useful functions and users should read the
  help file in its entirety for details.  However, we briefly mention
  several key functions that may make it easier to navigate and
  understand the layout of the package.
  
  \enumerate{
    \item \command{\link{rfsrc}}
    
    This is the main entry point to the package.  It grows a random forest
    using user supplied training data.  We refer to the resulting object
    as a RF-SRC grow object.  Formally, the resulting object has class
    \code{(rfsrc, grow)}.
    
    \item \command{\link{predict.rfsrc}} (\command{predict})
    
    Used for prediction.  Predicted values are obtained by dropping the
    user supplied test data down the grow forest.  The resulting object
    has class \code{(rfsrc, predict)}.
    
    
    \item \command{\link{max.subtree}}, \command{\link{var.select}}
    
    Used for variable selection.  The function \command{max.subtree}
    extracts maximal subtree information from a RF-SRC object which is
    used for selecting variables by making use of minimal depth variable
    selection.  The function \command{var.select} provides
    an extensive set of variable selection options and is a wrapper to
    \command{max.subtree}.
    
    \item \command{\link{impute.rfsrc}}
    
    Fast imputation mode for RF-SRC.  Both \command{rfsrc} and
    \command{predict.rfsrc} are capable of imputing missing data.
    However, for users whose only interest is imputing data, this function
    provides an efficient and fast interface for doing so.
  }
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.
  
  Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7(2):25-31.
  
  Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
  (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.
  
  Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and
  Lauer M.S. (2010).  High-dimensional variable selection for survival
  data.  \emph{J. Amer. Statist. Assoc.}, 105:205-217.
  
  Ishwaran H., Gerds T.A., Kogalur U.B., Moore R.D., Gange S.J. and Lau
  B.M. (2014). Random survival forests for competing risks.
  \emph{Biostatistics}, 15(4):757-773.
  
  Ishwaran H. (2015).  The effect of splitting on random forests.
  \emph{Machine Learning}, 99:75-118.
  
}
\seealso{
  \command{\link{find.interaction}},
  \command{\link{impute.rfsrc}},
  \command{\link{max.subtree}},
  \command{\link{plot.competing.risk}},
  \command{\link{plot.rfsrc}},
  \command{\link{plot.survival}},
  \command{\link{plot.variable}},
  \command{\link{predict.rfsrc}},
  \command{\link{print.rfsrc}},
  \command{\link{rf2rfz}},
  \command{\link{rfsrcSyn}},
  \command{\link{stat.split}}
  \command{\link{var.select}},
  \command{\link{vimp}}
}
\keyword{package}


