\name{rfsrc}
\alias{rfsrc}
\alias{rfsrc.cart}
\alias{get.auc}
\alias{get.bayes.rule}
\alias{get.brier.error}
\alias{get.cindex}
\alias{get.confusion}
\alias{get.logloss}
\alias{get.misclass.error}
\alias{get.mv.cserror}
\alias{get.mv.csvimp}
\alias{get.mv.error}
\alias{get.mv.error.block}
\alias{get.mv.formula}
\alias{get.mv.predicted}
\alias{get.mv.vimp}
\alias{randomForestSRC}
\title{Fast Unified Random Forests for Survival, Regression, and Classification (RF-SRC)}
\description{

  Fast OpenMP parallel computing of random forests (Breiman 2001) for
  regression, classification, survival analysis (Ishwaran et al. 2008),
  competing risks (Ishwaran et al. 2012), multivariate (Segal and Xiao
  2011), unsupervised (Mantero and Ishwaran 2020), quantile regression
  (Meinhausen 2006, Zhang et al. 2019, Greenwald-Khanna 2001), and class
  imbalanced q-classification (O'Brien and Ishwaran 2019).  Different
  splitting rules invoked under deterministic or random splitting
  (Geurts et al. 2006, Ishwaran 2015) are available for all families.
  Different types of variable importance (VIMP), holdout VIMP, as well
  as confidence regions (Ishwaran and Lu 2019) can be calculated for
  single and grouped variables.  Minimal depth variable selection
  (Ishwaran et al. 2010, 2011).  Fast interface for missing data
  imputation using a variety of different random forest methods (Tang
  and Ishwaran 2017).
  

  Highlighted new items:
  \enumerate{


    \item For computational speed, the default VIMP is no longer
    "permute" (Breiman-Cutler permutation importance) and has been
    switched to "anti" (\code{importance="anti"},
    \code{importance=TRUE}; see below for details). Be aware in some
    situations, such as highly imbalanced classification, that
    permutation VIMP may perform better. Permutation VIMP is
    obtained using \code{importance="permute"}.


    \item \code{save.memory} can be used for big data to save memory;
    especially useful for survival and competing risks.
    
    \item Mahalanobis splitting for multivariate regression with
    correlated y-outcomes (\code{splitrule="mahalanobis"}).  Now allows
    for a user specified covariance matrix.

    \item Visualize trees on your Safari or Google Chrome
    browser (works for all families).  See \command{\link{get.tree}}.

  }
  
  This is the main entry point to the \pkg{randomForestSRC}
  package.  For more information about this package and OpenMP parallel
  processing, use the command \code{package?randomForestSRC}.
  
}
\usage{
rfsrc(formula, data, ntree = 500,
  mtry = NULL, ytry = NULL,
  nodesize = NULL, nodedepth = NULL,
  splitrule = NULL, nsplit = NULL,
  importance = c(FALSE, TRUE, "none", "anti", "permute", "random"),
  block.size = if (any(is.element(as.character(importance),
                     c("none", "FALSE")))) NULL else 10,
  bootstrap = c("by.root", "none", "by.user"),
  samptype = c("swor", "swr"), samp = NULL, membership = FALSE,
  sampsize = if (samptype == "swor") function(x){x * .632} else function(x){x},
  na.action = c("na.omit", "na.impute"), nimpute = 1,
  ntime = 150, cause,
  perf.type = NULL,
  proximity = FALSE, distance = FALSE, forest.wt = FALSE,
  xvar.wt = NULL, yvar.wt = NULL, split.wt = NULL, case.wt  = NULL,
  case.depth = FALSE,
  forest = TRUE,
  save.memory = FALSE,
  var.used = c(FALSE, "all.trees", "by.tree"),
  split.depth = c(FALSE, "all.trees", "by.tree"),
  seed = NULL,
  do.trace = FALSE,
  statistics = FALSE,
  ...)

## convenient interface for growing a CART tree
rfsrc.cart(formula, data, ntree = 1, mtry = ncol(data), bootstrap = "none", ...)
}

\arguments{

  \item{formula}{Object of class 'formula' describing the model to fit.
    Interaction terms are not supported. If missing, unsupervised
    splitting is implemented.}

  \item{data}{Data frame containing the y-outcome and x-variables.}

  \item{ntree}{Number of trees.}

  \item{mtry}{Number of variables to possibly split at each node.
    Default is number of variables divided by 3 for regression.  For all
    other families (including unsupervised settings), the square root
    of number of variables.  Values are rounded up.}
  
  \item{ytry}{The number of randomly selected pseudo-outcomes for
    unsupervised families (see details below).  Default is
    \code{ytry}=1.}

  \item{nodesize}{Minumum size of terminal node.  The defaults are:
    survival (15), competing risk (15), regression (5), classification
    (1), mixed outcomes (3), unsupervised (3).  It is recommended to
    experiment with different \code{nodesize} values.}
  
  \item{nodedepth}{Maximum depth to which a tree should be grown.
    Parameter is ignored by default.}

  \item{splitrule}{Splitting rule (see below).}

  \item{nsplit}{Non-negative integer specifying number of random splits
    for splitting a variable.  When zero, all split values are
    used (deterministic splitting), which can be slower.  By default 
    10 is used.}

 \item{importance}{Method for computing variable importance (VIMP); see
    below.  Default action is \code{importance="none"} but VIMP can
    be recovered later using \command{vimp} or \command{predict}.}

  \item{block.size}{Determines how cumulative error rate is calculated.
    When \code{NULL}, it is calculated once for the entire forest,
    resulting in a flat cumulative error rate plot. To view the
    cumulative error rate for every nth tree, set \code{block.size} to
    an integer between \code{1} and \code{ntree}. If importance is
    requested, VIMP is calculated in blocks of size \code{block.size},
    balancing between ensemble and tree VIMP. The default is 10 trees.}
  
  \item{bootstrap}{Bootstrap protocol. Default is \code{by.root}, which
    bootstraps the data by sampling with or without replacement (default
    is without replacement; see \code{samptype} below).  If \code{none},
    the data is not bootstrapped (OOB ensembles or prediction error
    cannot be returned).  If \code{by.user}, the bootstrap specified by
    \code{samp} is used.}

  \item{samptype}{Type of bootstrap used when \code{by.root} is in
    effect.  Choices are \code{swor} (sampling without replacement; 
    the default) and \code{swr} (sampling with replacement).}

  \item{samp}{Bootstrap specification when \code{by.user} is in
    effect.  Array of dim \code{n x ntree} specifying how many
    times each record appears inbag in the bootstrap for each tree.}

  \item{membership}{Should terminal node membership and inbag
    information be returned?}


  \item{sampsize}{Specifies bootstrap size when \code{by.root} is in
    effect. For sampling without replacement, it is the requested size
    of the sample, defaulting to .632 times the sample size. For
    sampling with replacement, it is the sample size. Can also be
    specified as a number.}

  \item{na.action}{Action taken if the data contains
    \code{NA}s. Possible values are \code{na.omit} or
    \code{na.impute}. The default, \code{na.omit}, removes the entire
    record if any entry is \code{NA}. Selecting \code{na.impute} imputes
    the data (see below for details). Also see the function
    \command{impute} for fast imputation.}

  \item{nimpute}{Number of iterations of the missing data algorithm.
    Performance measures such as out-of-bag (OOB) error rates are
    optimistic if \code{nimpute} is greater than 1.}



  \item{ntime}{Integer value for survival to constrain ensemble
    calculations to an \code{ntime} grid of time points over the
    observed event times. If a vector of values with length greater than
    one is supplied, it is assumed these are the time points to be used
    (adjusted to match closest observed event times). Setting
    \code{ntime} to zero (or NULL) uses all observed event times.}

  \item{cause}{Integer value between 1 and \code{J} indicating the event
    of interest for splitting a node for competing risks, where \code{J}
    is the number of event types. If not specified, the default is a
    composite splitting rule averaging over all event types.  Can also
    be a vector of non-negative weights of length \code{J} specifying
    weights for each event (a vector of ones reverts to the default
    composite split statistic). Estimates for all event types are
    returned regardless of how \code{cause} is specified.}

  \item{perf.type}{Optional character value specifying the metric used for
    predicted value, variable importance (VIMP), and error rate.
    Defaults to the family metric if not specified.
    \code{perf.type="none"} turns off performance entirely, useful for
    turning off C-error rate calculations for large survival data (which can
    be expensive). Values allowed for univariate/multivariate
    classification are: \code{perf.type="misclass"} (default),
    \code{perf.type="brier"}, and \code{perf.type="gmean"}.}

  
  \item{proximity}{Proximity of cases as measured by the frequency of
    sharing the same terminal node.  This is an \code{n}x\code{n}
    matrix, which can be large.  Choices are \code{inbag}, \code{oob},
    \code{all}, \code{TRUE}, or \code{FALSE}.  Setting \code{proximity =
      TRUE} is equivalent to \code{proximity = "inbag"}.}
  
  \item{distance}{Distance between cases as measured by the ratio of the
    sum of the count of edges from each case to their immediate common
    ancestor node to the sum of the count of edges from each case to the
    root node.  If the cases are co-terminal for a tree, this measure is
    zero and reduces to 1 - the proximity measure.  This is an
    \code{n}x\code{n} matrix, which can be large.  Choices are
    \code{inbag}, \code{oob}, \code{all}, \code{TRUE}, or \code{FALSE}.
    Setting \code{distance = TRUE} is equivalent to \code{distance =
    "inbag"}.}
  
  \item{forest.wt}{Calculate the forest weight matrix?  Creates an
    \code{n}x\code{n} matrix which can be used for prediction and
    constructing customized estimators.  Choices are similar to
    proximity: \code{inbag}, \code{oob}, \code{all}, \code{TRUE}, or
    \code{FALSE}.  The default is \code{TRUE} which is equivalent to
    \code{inbag}.}

  \item{xvar.wt}{Vector of non-negative weights (does not have to sum
    to 1) representing the probability of selecting a variable for
    splitting.  Default is uniform weights.}

  
  \item{yvar.wt}{Vector of non-negative weights (does not have to sum
   to 1) representing the probability of selecting a response variable
   in multivariate regression.  Used when the number of responses are
   high-dimensional.  See the example below.}

  \item{split.wt}{Vector of non-negative weights used for multiplying
    the split statistic for a variable. A large value encourages the
    node to split on a specific variable. Default is uniform
    weights.}

  \item{case.wt}{Vector of non-negative weights (does not have to sum to
    1) for sampling cases.  Observations with larger weights will be
    selected with higher probability in the bootstrap (or subsampled)
    samples.  It is generally better to use real weights rather than
    integers. See the breast data example below illustrating its use
    for class imbalanced data.}
  
  \item{case.depth}{Return a matrix recording the depth at which a case
    first splits in a tree. Default is \code{FALSE}.}

  \item{forest}{Save key forest values?  Used for prediction on new data
    and required by many of the package functions. Turn this off if you
    are only interested in training a forest.}

  \item{save.memory}{Save memory?  Default is to store terminal node
    quantities used for prediction on test data.  This yields rapid
    prediction but can be memory intensive for big data, especially
    competing risks and survival models.  Turn this flag off in those
    cases.}
  
  \item{var.used}{Return statistics on number of times a variable split?
    Default is \code{FALSE}.  Possible values are \code{all.trees} which
    returns total number of splits of each variable, and \code{by.tree}
    which returns a matrix of number a splits for each variable for each
    tree.}
  
  \item{split.depth}{Records the minimal depth for each variable.
    Default is \code{FALSE}.  Possible values are \code{all.trees} which
    returns a matrix of the average minimal depth for a variable
    (columns) for a specific case (rows), and \code{by.tree} which
    returns a three-dimensional array recording minimal depth for a
    specific case (first dimension) for a variable (second dimension)
    for a specific tree (third dimension).}
  
  \item{seed}{Negative integer specifying seed for the random number
    generator.}

  \item{do.trace}{Number of seconds between updates to the user on
    approximate time to completion.}

  \item{statistics}{Should split statistics be returned?  Values can be
    parsed using \command{stat.split}.}

  \item{...}{Further arguments passed to or from other methods.}
  
}
\details{
  \enumerate{
    
    \item \emph{Types of forests}

    There is no need to set the type of forest as the package
    automagically determines the underlying random forest requested from
    the type of outcome and the formula supplied.  There are
    several possible scenarios:
    
    \enumerate{
      \item Regression forests for continuous outcomes.
      \item Classification forests for factor outcomes.
      \item Multivariate forests for continuous and/or factor outcomes
            and for mixed (both type) of outcomes.	    
      \item Unsupervised forests when there is no outcome.
      \item Survival forests for right-censored survival.
      \item Competing risk survival forests for competing risk.
    }
    
    \item \emph{Splitting}

    \enumerate{

      \item Splitting rules are specified by the option \code{splitrule}.
      
      \item For all families, pure random splitting can be invoked by setting
            \code{splitrule="random"}.

      \item For all families, computational speed can be increased using
            randomized splitting invoked by the option \code{nsplit}.
            See Improving Computational Speed.
    }

    \item \emph{Available splitting rules}

    \itemize{
      
      \item Regression analysis:
      \enumerate{

        \item \code{splitrule="mse"} (default split rule):  weighted
	mean-squared error splitting (Breiman et al. 1984, Chapter 8.4).
	
        \item \code{splitrule="quantile.regr"}: quantile regression splitting via
	the "check-loss" function.  Requires specifying the target
	quantiles.  See \code{quantreg.rfsrc} for further details.

        \item \code{la.quantile.regr}: local adaptive quantile
        regression splitting. See \code{quantreg.rfsrc}.
		
      }
      
       \item Classification analysis:
      \enumerate{

	\item \code{splitrule="gini"} (default splitrule): Gini
        index splitting (Breiman et al. 1984, Chapter 4.3).
	
        \item \code{splitrule="auc"}: AUC (area under the ROC curve) splitting
	for both two-class and multiclass setttings.  AUC splitting is
	appropriate for imbalanced data.  See \code{imbalanced} for
	more information.

	\item \code{splitrule="entropy"}: entropy splitting (Breiman et
	al. 1984, Chapter 2.5, 4.3).
	
      }
      
      \item Survival analysis:
      \enumerate{
	
        \item \code{splitrule="logrank"} (default splitrule):
	      log-rank splitting (Segal, 1988; Leblanc and Crowley, 1993).

	\item \code{splitrule="bs.gradient"}: gradient-based (global non-quantile)
	brier score splitting.  The time horizon used for the Brier
	score is set to the 90th percentile of the observed event times.
	This can be over-ridden by the option \code{prob}, which must be
	a value between 0 and 1 (set to .90 by default).
	      
        \item \code{splitrule="logrankscore"}:  log-rank score splitting (Hothorn and
        Lausen, 2003).
	
	
      }
      
      \item Competing risk analysis (for details see Ishwaran et al., 2014):
      \enumerate{
        
        \item \code{splitrule="logrankCR"} (default splitrule): modified
	weighted log-rank splitting rule modeled after Gray's test
	(Gray, 1988).  Use this to find *all* variables 
	that are informative and when the goal is long term
	prediction. 
	      
        
        \item \code{splitrule="logrank"}: weighted log-rank splitting where each
	event type is treated as the event of interest and all
	other events are treated as censored.  The split rule is
	the weighted value of each of log-rank statistics,
	standardized by the variance.  Use this to find variables
	that affect a *specific* cause of interest and when
	the goal is a targeted analysis of a specific cause.
	However in order for this to be effective, remember to set
	the \code{cause} option to the targeted cause of interest.
	See examples below.
	      
      }
	    
      \item Multivariate analysis:
      \enumerate{

         \item Default is the multivariate normalized composite split
	 rule using mean-squared error and Gini index (Tang and
	 Ishwaran, 2017).

	 \item \code{splitrule="mahalanobis"}: Mahalanobis splitting
	 that adjusts for correlation (also allows for a user specified
	 covariance matrix, see example below).  Only works for
	 multivariate regression (all outcomes must be real).

      }
       
      \item Unsupervised analysis: In settings where there is no
            outcome, unsupervised splitting that uses pseudo-outcomes is
            applied using the default multivariate splitting rule (see
            below for details) Also see \command{sidClustering} for a
            more sophisticated method for unsupervised analysis (Mantero
            and Ishwaran, 2020).
            

      \item Custom splitting: All families except unsupervised are
            available for user defined custom splitting.  Some basic
            C-programming skills are required.  The harness for defining
            these rules is in \code{splitCustom.c}.  In this file we
            give examples of how to code rules for regression,
            classification, survival, and competing risk.  Each family
            can support up to sixteen custom split rules.  Specifying
            \code{splitrule="custom"} or \code{splitrule="custom1"} will
            trigger the first split rule for the family defined by the
            training data set. Multivariate families will need a custom
            split rule for both regression and classification.  In the
            examples, we demonstrate how the user is presented with the
            node specific membership.  The task is then to define a
            split statistic based on that membership.  Take note of the
            instructions in \code{splitCustom.c} on how to
            \emph{register} the custom split rules.  It is suggested
            that the existing custom split rules be kept in place for
            reference and that the user proceed to develop
            \code{splitrule="custom2"} and so on. The package must be
            recompiled and installed for the custom split rules to
            become available.

    }
    
    
    \item \emph{Improving computational speed}
    
    See the function \command{\link{rfsrc.fast}} for a fast
    implementation of \code{rfsrc}.  Key methods for
    increasing speed are as follows:
    
      \itemize{

	\item \emph{Nodesize}

	Increasing \code{nodesize} has the greatest effect in speeding
	calculations.  In some big data settings this can also lead to
	better prediction performance.


	\item \emph{Save memory}

	Use option \code{save.memory="TRUE"} for big data competing risk
	and survival models.  By default the package stores terminal
	node quantities to be used in prediction for test data but this
	can be memory intensive for big data.

	
        \item \emph{Block size}

	Make sure \code{block.size="NULL"} (or set to number of trees)
	so that the cumulative error is calculated only once.

	\item \emph{Turn off performace}

	The C-error rate calculation can be very expensive for big
	survival data.  Set \code{perf.type="none"} to turn this off and
	all other performance calculations (then consider using the
	function \code{get.brier.survival} as a fast way to get survival
	performance).

	\code{perf.type="none"} applies to all other families as well.
	
	\item \emph{Randomized splitting rules}

	Set \code{nsplit} to a small non-zero integer value.  Then a
        maximum of \code{nsplit} split points are chosen randomly for
        each of the candidate splitting variables when splitting a tree
        node, thus significantly reducing computational costs.

        For more details about randomized splitting see Loh and Shih
        (1997), Dietterich (2000), and Lin and Jeon (2006).  Geurts et
        al. (2006) introduced extremely randomized trees using the
        extra-trees algorithm.  This algorithm corresponds to
        \code{nsplit}=1.  In our experience however this may be too low
        for general use (Ishwaran, 2015).

        For completely randomized (pure random) splitting use
        \code{splitrule="random"}.  In pure splitting, nodes are split
        by randomly selecting a variable and randomly selecting its
        split point (Cutler and Zhao, 2001).

        \item \emph{Subsampling}

        Reduce the size of the bootstrap using \code{sampsize} and
        \code{samptype}.  See \command{rfsrc.fast} for a fast forest
        implementation using subsampling.

        \item \emph{Unique time points}

	Setting \code{ntime} to a reasonably small value such as 50
	constrains survival ensemble calculations to a restricted grid
	of time points and significantly improves computational times.
      
	\item \emph{Large number of variables}

	Try filtering variables ahead of time.  Make sure not to request
	VIMP (variable importance can always be recovered later using
	\command{vimp} or \command{predict}).  Also if variable
	selection is desired, but is too slow, consider using
	\command{max.subtree} which calculates minimal depth, a measure
	of the depth that a variable splits, and yields fast variable
	selection (Ishwaran, 2010).
	
      }
    
    \item \emph{Prediction Error}
    
    Prediction error is calculated using OOB data.  The metric used is
    mean-squared-error for regression, and misclassification error for
    classification.  A normalized Brier score (relative to a coin-toss)
    and the AUC (area under the ROC curve) is also provided upon
    printing a classification forest.  Performance for Brier score can
    be specified using \code{perf.type="brier"}.  G-mean performance is
    also available, see the function \command{imbalanced} for more
    details.
    
    For survival, prediction error is measured by the C-error rate
    defined as 1-C, where C is Harrell's (Harrell et al., 1982)
    concordance index.  The C-error is between 0 and 1, and measures
    error in ranking (classifying) two random individuals in terms of
    survival.  A value of 0.5 is no better than random guessing. A value
    of 0 is perfect.
    
    When bootstrapping is by \code{none}, a coherent OOB subset is not
    available to assess prediction error.  Thus, all outputs dependent
    on this are suppressed.  In such cases, prediction error is only
    available via classical cross-validation (the user will need to use
    the \command{predict.rfsrc} function).
    
    \item \emph{Variable Importance (VIMP)}

     VIMP is calculated using OOB data in several ways.
     \code{importance="permute"} yields permutation VIMP (Breiman-Cutler
     importance) by permuting OOB cases. \code{importance="random"} uses
     random left/right assignments whenever a split is encountered for
     the target variable.  The default \code{importance="anti"}
     (equivalent to \code{importance=TRUE}) assigns cases to the anti
     (opposite) split.

    VIMP depends upon \code{block.size}, an integer value between 1 and
    \code{ntree}, specifying number of trees in a block used for VIMP.
    When \code{block.size}=1, VIMP is calculated for each tree.  When
    \code{block.size="ntree"}, VIMP is calculated for the entire forest
    by comparing the perturbed OOB forest ensemble (using all trees) to
    the unperturbed OOB forest ensemble (using all trees).  This yields
    ensemble VIMP, which does not measure the tree average effect of a
    variable, but rather its overall forest effect.  

    A useful compromise between tree VIMP and ensemble VIMP can be
    obtained by setting \code{block.size} to a value between 1 and
    \code{ntree}.  Smaller values generally gives better accuracy,
    however computational times will be higher because VIMP is
    calculated over more blocks.  However, see \command{imbalanced} for
    imbalanced classification data where larger \code{block.size}
    often works better (O'Brien and Ishwaran, 2019).

    See \command{vimp} for a user interface for extracting VIMP and
    \command{subsampling} for calculating confidence intervals for VIMP.

    Also see \command{holdout.vimp} for holdout VIMP, which calculates
    importance by holding out variables.  This is more conservative, but
    with good false discovery properties.

    For classification, VIMP is returned as a matrix with J+1
    columns where J is the number of classes.  The first column "all" is
    the unconditional VIMP, while the remaining columns are conditional VIMP
    calculated using only OOB cases with the class label.
      
    
    \item \emph{Multivariate Forests}
    
    Multivariate forests can be specified in two ways:
    
    rfsrc(Multivar(y1, y2, ..., yd) ~ . , my.data, ...)
    
    rfsrc(cbind(y1, y2, ..., yd) ~ . , my.data, ...)

    By default, a multivariate normalized composite splitting rule is
    used to split nodes (for multivariate regression, users have
    the option to use Mahalanobis splitting).

    The nature of the outcomes informs the code as to what type of
    multivariate forest is grown; i.e. whether it is real-valued,
    categorical, or a combination of both (mixed). Performance measures
    (when requested) are returned for all outcomes.

    Helper functions \command{get.mv.formula},
    \command{get.mv.predicted}, \command{get.mv.error} can be used for
    defining the multivariate forest formula and extracting predicted
    values (all outcomes) and VIMP (all variables, all outcomes;
    assuming importance was requested in the call).  The latter two
    functions also work for univariate (regular) forests.  Both
    functions return standardized values (dividing by the variance for
    regression, or multiplying by 100, otherwise) using option
    \code{standardize="TRUE"}.
	 
    \item \emph{Unsupervised Forests and sidClustering}

    See sidClustering \command{sidClustering} for a more sophisticated
    method for unsupervised analysis.

    Otherwise a more direct (but naive) way to proceed is to use the
    unsupervised splitting rule.  The following are equivalent ways to
    grow an unsupervised forest via unsupervised splitting:
    
    rfsrc(data = my.data)
    
    rfsrc(Unsupervised() ~ ., data = my.data)
    
    In unsupervised mode, features take turns acting as target
    y-outcomes and x-variables for splitting.  Specifically, \code{mtry}
    x-variables are randomly selected for splitting the node.  Then for
    each \code{mtry} feature, \code{ytry} variables are selected
    from the remaining features to act as the target pseduo-outcomes.
    Splitting uses the multivariate normalized composite splitting rule.

    The default value of \code{ytry} is 1 but can be increased.  As
    illustration, the following equivalent unsupervised calls set
    \code{mtry=10} and \code{ytry=5}:

    rfsrc(data = my.data, ytry = 5, mtry = 10)
    
    rfsrc(Unsupervised(5) ~ ., my.data, mtry = 10)
    
    Note that all performance values (error rates, VIMP, prediction) are
    turned off in unsupervised mode.
    
    \item \emph{Survival, Competing Risks}
    
    \enumerate{
      \item Survival settings require a time and censoring variable which
      should be identifed in the formula as the outcome using the standard
      \code{Surv} formula specification.  A typical formula call looks like:
      
      Surv(my.time, my.status) ~ .
      
      where \code{my.time} and \code{my.status} are the variables names for
      the event time and status variable in the users data set.
      
      \item For survival forests (Ishwaran et al. 2008), the censoring
      variable must be coded as a non-negative integer with 0 reserved for
      censoring and (usually) 1=death (event). 
      
      \item For competing risk forests (Ishwaran et al., 2013), the
      implementation is similar to survival, but with the following
      caveats:
      
      \itemize{
        \item Censoring must be coded as a non-negative integer, where 0
        indicates right-censoring, and non-zero values indicate different
        event types.  While \code{0,1,2,..,J} is standard, and recommended,
        events can be coded non-sequentially, although 0 must always be used
        for censoring.
        
        \item Setting the splitting rule to \code{logrankscore} will result
        in a survival analysis in which all events are treated as if they
        are the same type (indeed, they will coerced as such).
        
        \item Generally, competing risks requires a larger \code{nodesize} than
        survival settings.
      }
    }
    
    \item \emph{Missing data imputation}
    
    \code{na.action="na.impute"} imputes missing data (both x and
    y-variables) using the missing data algorithm of Ishwaran et
    al. (2008).  But also see the \command{impute} for an
    alternate way to do fast and accurate imputation.
    
    The missing data algorithm can be iterated by setting \code{nimpute}
    to a positive integer greater than 1.  When iterated, at the
    completion of each iteration, missing data is imputed using OOB
    non-missing terminal node data which is then used as input to grow a
    new forest.  A side effect of iteration is that missing values in
    the returned objects \code{xvar}, \code{yvar} are replaced by
    imputed values.  In other words the incoming data is overlaid with
    the missing data.  Also, performance measures such as error rates
    and VIMP become optimistically biased.
    
    Records in which all outcome and x-variable information are missing
    are removed from the forest analysis.  Variables having all missing
    values are also removed.
    
   \item \emph{Allowable data types and factors}

    Data types must be real valued, integer, factor or logical --
    however all except factors are coerced and treated as if real
    valued.  For ordered x-variable factors, splits are similar to real
    valued variables.  For unordered factors, a split will move a subset
    of the levels in the parent node to the left daughter, and the
    complementary subset to the right daughter.  All possible
    complementary pairs are considered and apply to factors with an
    unlimited number of levels.  However, there is an optimization check
    to ensure number of splits attempted is not greater than number of
    cases in a node or the value of \code{nsplit}.

    For coherence, an immutable map is applied to each factor that
    ensures factor levels in the training data are consistent with the
    factor levels in any subsequent test data.  This map is applied to
    each factor before and after the native C library is executed.
    Because of this, if all x-variables all factors, then computational
    time will be long in high dimensional problems.  Consider converting
    factors to real if this is the case.


    
}}
\value{
  An object of class \code{(rfsrc, grow)} with the following
  components:
  
  \item{call}{The original call to \command{rfsrc}.}
  \item{family}{The family used in the analysis.}
  \item{n}{Sample size of the data (depends upon \code{NA}'s, see \code{na.action}).}
  \item{ntree}{Number of trees grown.}
  \item{mtry}{Number of variables randomly selected for splitting at each node.}
  \item{nodesize}{Minimum size of terminal nodes.}
  \item{nodedepth}{Maximum depth allowed for a tree.}
  \item{splitrule}{Splitting rule used.}
  \item{nsplit}{Number of randomly selected split points.}
  \item{yvar}{y-outcome values.}
  \item{yvar.names}{A character vector of the y-outcome names.}
  \item{xvar}{Data frame of x-variables.}
  \item{xvar.names}{A character vector of the x-variable names.}
  \item{xvar.wt}{Vector of non-negative weights specifying the
    probability used to select a variable for splitting a node.}
  \item{split.wt}{Vector of non-negative weights specifying
    multiplier by which the split statistic for a covariate is adjusted.}
  \item{cause.wt}{Vector of weights used for the composite competing
    risk splitting rule.}
  \item{leaf.count}{Number of terminal nodes for each tree in the
    forest. Vector of length \code{ntree}.  A value of zero indicates
    a rejected tree (can occur when imputing missing data).
    Values of one indicate tree stumps.}
  \item{proximity}{Proximity matrix recording the frequency of pairs of data points
    occur within the same terminal node.}
  \item{forest}{If \code{forest=TRUE}, the forest object is returned.
    This object is used for prediction with new test data
	  sets and is required for other R-wrappers.}
  \item{forest.wt}{Forest weight matrix.}
  \item{membership}{Matrix recording terminal node membership where
    each column records node mebership for a case for a tree (rows).}
  \item{splitrule}{Splitting rule used.}
  \item{inbag}{Matrix recording inbag membership where each column
    contains the number of times that a case appears in the bootstrap
    sample for a tree (rows).}
  \item{var.used}{Count of the number of times a variable is used in
    growing the forest.}
  \item{imputed.indv}{Vector of indices for cases with missing
    values.}
  \item{imputed.data}{Data frame of the imputed data. The first
	  column(s) are reserved for the y-outcomes, after which the
	  x-variables are listed.}
  \item{split.depth}{Matrix (i,j) or array (i,j,k) recording the
    minimal depth for variable j for case i, either averaged over
    the forest, or by tree k.}
  \item{node.stats}{Split statistics returned when
    \code{statistics=TRUE} which can be parsed using \command{stat.split}.}

  \item{err.rate}{Tree cumulative OOB error rate.}
  \item{err.block.rate}{When \code{importance=TRUE}, vector of the
  cumulative error rate for each ensemble block comprised of
  \code{block.size} trees.  So with \code{block.size = 10}, entries are
  the cumulative error rate for the first 10 trees, the first 20 trees,
  30 trees, and so on.  As another exmple, if \code{block.size = 1},
  entries are the error rate for each tree.}

  
    
  \item{importance}{Variable importance (VIMP) for each x-variable.}
  \item{predicted}{In-bag predicted value.}
  \item{predicted.oob}{OOB predicted value.}\cr
  
  \item{++++++++}{for classification settings, additionally ++++++++} \cr
  \item{class}{In-bag predicted class labels.}
  \item{class.oob}{OOB predicted class labels.}\cr
  
  
  \item{++++++++}{for multivariate settings, additionally ++++++++} \cr
  
  \item{regrOutput}{List containing performance values for
	  multivariate regression outcomes (applies only in
	  multivariate settings).}
	\item{clasOutput}{List containing performance values for
	  multivariate categorical (factor) outcomes (applies only in
	  multivariate settings).}\cr
  
  \item{++++++++}{for survival settings, additionally ++++++++} \cr
  
  \item{survival}{In-bag survival function.}
  \item{survival.oob}{OOB survival function.}
  \item{chf}{In-bag cumulative hazard function (CHF).}
  \item{chf.oob}{OOB CHF.}
  \item{time.interest}{Ordered unique death times.}
  \item{ndead}{Number of deaths.}\cr
  
  \item{++++++++}{for competing risks, additionally ++++++++} \cr
  
  \item{chf}{In-bag cause-specific cumulative hazard function (CSCHF)
    for each event.}
  \item{chf.oob}{OOB CSCHF.}
  \item{cif}{In-bag cumulative incidence function (CIF) for each event.}
  \item{cif.oob}{OOB CIF.}
}
\note{

  Values returned depend heavily on the family.  In particular,
  predicted values from the forest (\code{predicted} and
  \code{predicted.oob}) are as follows:
    
    \enumerate{
      \item For regression, a vector of predicted y-outcomes.
      
      \item For classification, a matrix with columns containing the
      estimated class probability for each class.  Performance values
      and VIMP for classification are reported as a matrix with J+1
      columns where J is the number of classes.  The first column "all"
      is the unconditional value for performance (VIMP), while the
      remaining columns are performance (VIMP) conditioned on cases
      corresponding to that class label.
      
      \item For survival, a vector of mortality values (Ishwaran et al.,
      2008) representing estimated risk for each individual calibrated
      to the scale of the number of events (as a specific example, if
      \emph{i} has a mortality value of 100, then if all individuals had
      the same x-values as \emph{i}, we would expect an average of 100
      events).  Also returned are matrices containing
      the CHF and survival function.  Each row corresponds to an
      individual's ensemble CHF or survival function evaluated at each
      time point in \code{time.interest}.  
      
      \item For competing risks, a matrix with one column for each event
      recording the expected number of life years lost due to the event
      specific cause up to the maximum follow up (Ishwaran et al.,
      2013).  Also returned are the cause-specific cumulative hazard
      function (CSCHF) and the cumulative incidence function (CIF) for
      each event type.  These are encoded as a three-dimensional array,
      with the third dimension used for the event type, each time point
      in \code{time.interest} making up the second dimension (columns),
      and the case (individual) being the first dimension (rows).
      
      \item For multivariate families, predicted values (and other
      performance values such as VIMP and error rates) are stored in the
      lists \code{regrOutput} and \code{clasOutput} which can be
      extracted using functions \code{get.mv.error},
      \code{get.mv.predicted} and \code{get.mv.vimp}.
      
    }
  }
  
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}

\references{
  Breiman L., Friedman J.H., Olshen R.A. and Stone C.J. (1984).
  \emph{Classification and Regression Trees}, Belmont, California.
  
  Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.
  
  Cutler A. and Zhao G. (2001). PERT-Perfect random tree ensembles.
  \emph{Comp. Sci. Statist.}, 33: 490-497.

  Dietterich, T. G. (2000). An experimental comparison of three methods for
  constructing ensembles of decision trees: bagging, boosting, and randomization.
  \emph{Machine Learning}, 40, 139-157.
  
  Gray R.J. (1988).  A class of k-sample tests for comparing the
  cumulative incidence of a competing risk, \emph{Ann. Statist.},
  16: 1141-1154.
  
  Geurts, P., Ernst, D. and Wehenkel, L., (2006). Extremely randomized
  trees. \emph{Machine learning}, 63(1):3-42.

  Greenwald M. and Khanna S. (2001).  Space-efficient online computation of
  quantile summaries. \emph{Proceedings of ACM SIGMOD}, 30(2):58-66.

  Harrell et al. F.E. (1982).  Evaluating the yield of medical tests,
  \emph{J. Amer. Med. Assoc.}, 247:2543-2546.
  
  Hothorn T. and Lausen B. (2003). On the exact distribution of maximally selected
  rank statistics, \emph{Comp. Statist. Data Anal.}, 43:121-137.
  
  Ishwaran H. (2007).  Variable importance in binary regression
  trees and forests,  \emph{Electronic J. Statist.}, 1:519-537.
  
  Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7(2):25-31.
  
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

  Lin, Y. and Jeon, Y. (2006). Random forests and adaptive nearest
  neighbors. \emph{J. Amer. Statist. Assoc.}, 101(474), 578-590.
  
  Lu M., Sadiq S., Feaster D.J. and Ishwaran H. (2018). Estimating
  individual treatment effect in observational data using random forest
  methods. \emph{J. Comp. Graph. Statist}, 27(1), 209-219

  Ishwaran H. and Lu M.  (2019).  Standard errors and confidence
  intervals for variable importance in random forest regression,
  classification, and survival. \emph{Statistics in Medicine}, 38,
  558-582.
  
  LeBlanc M. and Crowley J. (1993).  Survival trees by goodness of split,
  \emph{J. Amer. Statist. Assoc.}, 88:457-467.
  
  Loh W.-Y and Shih Y.-S (1997).  Split selection methods for
  classification trees, \emph{Statist. Sinica}, 7:815-840.

  Mantero A. and Ishwaran H. (2021).  Unsupervised random forests.
  \emph{Statistical Analysis and Data Mining}, 14(2):144-167.

  Meinshausen N. (2006) Quantile regression forests, \emph{Journal of
  Machine Learning Research}, 7:983-999.

  Mogensen, U.B, Ishwaran H. and Gerds T.A. (2012). Evaluating random
  forests for survival analysis using prediction error curves,
  \emph{J.  Statist. Software}, 50(11): 1-23.
  
  O'Brien R. and Ishwaran H. (2019).  A random forests quantile
  classifier for class imbalanced data. \emph{Pattern Recognition},
  90, 232-249

  Segal M.R. (1988).  Regression trees for censored data,
  \emph{Biometrics},  44:35-47.

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
  \command{\link{rfsrc.anonymous}},
  \command{\link{rfsrc.cart}},
  \command{\link{rfsrc.fast}},

  
  \command{\link{sidClustering.rfsrc}},

  \command{\link{stat.split.rfsrc}},
  \command{\link{subsample.rfsrc}},
  \command{\link{synthetic.rfsrc}},
  
  \command{\link{tune.rfsrc}},

  \command{\link{var.select.rfsrc}},
  \command{\link{vimp.rfsrc}}
}
\examples{
\donttest{
##------------------------------------------------------------
## survival analysis
##------------------------------------------------------------

## veteran data
## randomized trial of two treatment regimens for lung cancer
data(veteran, package = "randomForestSRC")
v.obj <- rfsrc(Surv(time, status) ~ ., data = veteran, block.size = 1)

## plot tree number 3
plot(get.tree(v.obj, 3))

## print results of trained forest
print(v.obj)

## plot results of trained forest
plot(v.obj)

## plot survival curves for first 10 individuals -- direct way
matplot(v.obj$time.interest, 100 * t(v.obj$survival.oob[1:10, ]),
    xlab = "Time", ylab = "Survival", type = "l", lty = 1)

## plot survival curves for first 10 individuals
## using function "plot.survival" 
plot.survival(v.obj, subset = 1:10)

## obtain Brier score using KM and RSF censoring distribution estimators
bs.km <- get.brier.survival(v.obj, cens.model = "km")$brier.score
bs.rsf <- get.brier.survival(v.obj, cens.model = "rfsrc")$brier.score

## plot the brier score
plot(bs.km, type = "s", col = 2)
lines(bs.rsf, type ="s", col = 4)
legend("topright", legend = c("cens.model = km", "cens.model = rfsrc"), fill = c(2,4))

## plot CRPS (continuous rank probability score) as function of time
## here's how to calculate the CRPS for every time point
trapz <- randomForestSRC:::trapz
time <- v.obj$time.interest
crps.km <- sapply(1:length(time), function(j) {
  trapz(time[1:j], bs.km[1:j, 2] / diff(range(time[1:j])))
})
crps.rsf <- sapply(1:length(time), function(j) {
  trapz(time[1:j], bs.rsf[1:j, 2] / diff(range(time[1:j])))
})
plot(time, crps.km, ylab = "CRPS", type = "s", col = 2)
lines(time, crps.rsf, type ="s", col = 4)
legend("bottomright", legend=c("cens.model = km", "cens.model = rfsrc"), fill=c(2,4))


## fast nodesize optimization for veteran data
## optimal nodesize in survival is larger than other families
## see the function "tune" for more examples
tune.nodesize(Surv(time,status) ~ ., veteran)


## Primary biliary cirrhosis (PBC) of the liver
data(pbc, package = "randomForestSRC")
pbc.obj <- rfsrc(Surv(days, status) ~ ., pbc)
print(pbc.obj)


## save.memory example for survival
## growing many deep trees creates memory issue without this option!
data(pbc, package = "randomForestSRC")
print(rfsrc(Surv(days, status) ~ ., pbc, splitrule = "random",
            ntree = 25000, nodesize = 1, save.memory = TRUE))



##------------------------------------------------------------
## trees can be plotted for any family 
## see get.tree for details and more examples
##------------------------------------------------------------

## survival where factors have many levels
data(veteran, package = "randomForestSRC")
vd <- veteran
vd$celltype=factor(vd$celltype)
vd$diagtime=factor(vd$diagtime)
vd.obj <- rfsrc(Surv(time,status)~., vd, ntree = 100, nodesize = 5)
plot(get.tree(vd.obj, 3))

## classification
iris.obj <- rfsrc(Species ~., data = iris)
plot(get.tree(iris.obj, 25, class.type = "bayes"))
plot(get.tree(iris.obj, 25, target = "setosa"))
plot(get.tree(iris.obj, 25, target = "versicolor"))
plot(get.tree(iris.obj, 25, target = "virginica"))

## ------------------------------------------------------------
## simple example of VIMP using iris classification
## ------------------------------------------------------------

## directly from trained forest
print(rfsrc(Species~.,iris,importance=TRUE)$importance)

## VIMP (and performance) use misclassification error by default
## but brier prediction error can be requested
print(rfsrc(Species~.,iris,importance=TRUE,perf.type="brier")$importance)

## example using vimp function (see vimp help file for details)
iris.obj <- rfsrc(Species ~., data = iris)
print(vimp(iris.obj)$importance)
print(vimp(iris.obj,perf.type="brier")$importance)

## example using hold out vimp (see holdout.vimp help file for details)
print(holdout.vimp(Species~.,iris)$importance)
print(holdout.vimp(Species~.,iris,perf.type="brier")$importance)

## ------------------------------------------------------------
## confidence interval for vimp using subsampling
## compare with holdout vimp
## ------------------------------------------------------------

## new York air quality measurements
o <- rfsrc(Ozone ~ ., data = airquality)
so <- subsample(o)
plot(so)

## compare with holdout vimp
print(holdout.vimp(Ozone ~ ., data = airquality)$importance)


##------------------------------------------------------------
## example of imputation in survival analysis
##------------------------------------------------------------

data(pbc, package = "randomForestSRC")
pbc.obj2 <- rfsrc(Surv(days, status) ~ ., pbc, na.action = "na.impute")

## same as above but iterate the missing data algorithm
pbc.obj3 <- rfsrc(Surv(days, status) ~ ., pbc,
         na.action = "na.impute", nimpute = 3)

## fast way to impute data (no inference is done)
## see impute for more details
pbc.imp <- impute(Surv(days, status) ~ ., pbc, splitrule = "random")

##------------------------------------------------------------
## compare RF-SRC to Cox regression
## Illustrates C-error and Brier score measures of performance
## assumes "pec" and "survival" libraries are loaded
##------------------------------------------------------------

if (library("survival", logical.return = TRUE)
    & library("pec", logical.return = TRUE)
    & library("prodlim", logical.return = TRUE))

{
  ##prediction function required for pec
  predictSurvProb.rfsrc <- function(object, newdata, times, ...){
    ptemp <- predict(object,newdata=newdata,...)$survival
    pos <- sindex(jump.times = object$time.interest, eval.times = times)
    p <- cbind(1,ptemp)[, pos + 1]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop("Prediction failed")
    p
  }

  ## data, formula specifications
  data(pbc, package = "randomForestSRC")
  pbc.na <- na.omit(pbc)  ##remove NA's
  surv.f <- as.formula(Surv(days, status) ~ .)
  pec.f <- as.formula(Hist(days,status) ~ 1)

  ## run cox/rfsrc models
  ## for illustration we use a small number of trees
  cox.obj <- coxph(surv.f, data = pbc.na, x = TRUE)
  rfsrc.obj <- rfsrc(surv.f, pbc.na, ntree = 150)

  ## compute bootstrap cross-validation estimate of expected Brier score
  ## see Mogensen, Ishwaran and Gerds (2012) Journal of Statistical Software
  set.seed(17743)
  prederror.pbc <- pec(list(cox.obj,rfsrc.obj), data = pbc.na, formula = pec.f,
                        splitMethod = "bootcv", B = 50)
  print(prederror.pbc)
  plot(prederror.pbc)

  ## compute out-of-bag C-error for cox regression and compare to rfsrc
  rfsrc.obj <- rfsrc(surv.f, pbc.na)
  cat("out-of-bag Cox Analysis ...", "\n")
  cox.err <- sapply(1:100, function(b) {
    if (b\%\%10 == 0) cat("cox bootstrap:", b, "\n")
    train <- sample(1:nrow(pbc.na), nrow(pbc.na), replace = TRUE)
    cox.obj <- tryCatch({coxph(surv.f, pbc.na[train, ])}, error=function(ex){NULL})
    if (!is.null(cox.obj)) {
      get.cindex(pbc.na$days[-train], pbc.na$status[-train], predict(cox.obj, pbc.na[-train, ]))
    } else NA
  })
  cat("\n\tOOB error rates\n\n")
  cat("\tRSF            : ", rfsrc.obj$err.rate[rfsrc.obj$ntree], "\n")
  cat("\tCox regression : ", mean(cox.err, na.rm = TRUE), "\n")
}

##------------------------------------------------------------
## competing risks
##------------------------------------------------------------

## WIHS analysis
## cumulative incidence function (CIF) for HAART and AIDS stratified by IDU

data(wihs, package = "randomForestSRC")
wihs.obj <- rfsrc(Surv(time, status) ~ ., wihs, nsplit = 3, ntree = 100)
plot.competing.risk(wihs.obj)
cif <- wihs.obj$cif.oob
Time <- wihs.obj$time.interest
idu <- wihs$idu
cif.haart <- cbind(apply(cif[,,1][idu == 0,], 2, mean),
                   apply(cif[,,1][idu == 1,], 2, mean))
cif.aids  <- cbind(apply(cif[,,2][idu == 0,], 2, mean),
                   apply(cif[,,2][idu == 1,], 2, mean))
matplot(Time, cbind(cif.haart, cif.aids), type = "l",
        lty = c(1,2,1,2), col = c(4, 4, 2, 2), lwd = 3,
        ylab = "Cumulative Incidence")
legend("topleft",
       legend = c("HAART (Non-IDU)", "HAART (IDU)", "AIDS (Non-IDU)", "AIDS (IDU)"),
       lty = c(1,2,1,2), col = c(4, 4, 2, 2), lwd = 3, cex = 1.5)


## illustrates the various splitting rules
## illustrates event specific and non-event specific variable selection
if (library("survival", logical.return = TRUE)) {

  ## use the pbc data from the survival package
  ## events are transplant (1) and death (2)
  data(pbc, package = "survival")
  pbc$id <- NULL

  ## modified Gray's weighted log-rank splitting
  ## (equivalent to cause=c(1,1) and splitrule="logrankCR")
  pbc.cr <- rfsrc(Surv(time, status) ~ ., pbc)
 
  ## log-rank cause-1 specific splitting and targeted VIMP for cause 1
  pbc.log1 <- rfsrc(Surv(time, status) ~ ., pbc, 
              splitrule = "logrankCR", cause = c(1,0), importance = TRUE)

  ## log-rank cause-2 specific splitting and targeted VIMP for cause 2
  pbc.log2 <- rfsrc(Surv(time, status) ~ ., pbc, 
              splitrule = "logrankCR", cause = c(0,1), importance = TRUE)

  ## extract VIMP from the log-rank forests: event-specific
  ## extract minimal depth from the Gray log-rank forest: non-event specific
  var.perf <- data.frame(md = max.subtree(pbc.cr)$order[, 1],
                         vimp1 = 100 * pbc.log1$importance[ ,1],
                         vimp2 = 100 * pbc.log2$importance[ ,2])
  print(var.perf[order(var.perf$md), ], digits = 2)

}

## ------------------------------------------------------------
## regression analysis
## ------------------------------------------------------------

## new York air quality measurements
airq.obj <- rfsrc(Ozone ~ ., data = airquality, na.action = "na.impute")

# partial plot of variables (see plot.variable for more details)
plot.variable(airq.obj, partial = TRUE, smooth.lines = TRUE)

## motor trend cars
mtcars.obj <- rfsrc(mpg ~ ., data = mtcars)

## ------------------------------------------------------------
## regression with custom bootstrap
## ------------------------------------------------------------

ntree <- 25
n <- nrow(mtcars)
s.size <- n / 2
swr <- TRUE
samp <- randomForestSRC:::make.sample(ntree, n, s.size, swr)
o <- rfsrc(mpg ~ ., mtcars, bootstrap = "by.user", samp = samp)

## ------------------------------------------------------------
## classification analysis
## ------------------------------------------------------------

## iris data
iris.obj <- rfsrc(Species ~., data = iris)

## wisconsin prognostic breast cancer data
data(breast, package = "randomForestSRC")
breast.obj <- rfsrc(status ~ ., data = breast, block.size=1)
plot(breast.obj)

## ------------------------------------------------------------
## big data set, reduce number of variables using simple method
## ------------------------------------------------------------

## use Iowa housing data set
data(housing, package = "randomForestSRC")

## original data contains lots of missing data, use fast imputation
## however see impute for other methods
housing2 <- impute(data = housing, fast = TRUE)

## run shallow trees to find variables that split any tree
xvar.used <- rfsrc(SalePrice ~., housing2, ntree = 250, nodedepth = 4,
                   var.used="all.trees", mtry = Inf, nsplit = 100)$var.used

## now fit forest using filtered variables
xvar.keep  <- names(xvar.used)[xvar.used >= 1]
o <- rfsrc(SalePrice~., housing2[, c("SalePrice", xvar.keep)])
print(o)

## ------------------------------------------------------------
## imbalanced classification data
## see the "imbalanced" function for further details
##
## (a) use balanced random forests with undersampling of the majority class
## Specifically let n0, n1 be sample sizes for majority, minority
## cases.  We sample 2 x n1 cases with majority, minority cases chosen
## with probabilities n1/n, n0/n where n=n0+n1
##
## (b) balanced random forests using "imbalanced"
##
## (c) q-classifier (RFQ) using "imbalanced" 
##
## ------------------------------------------------------------

## Wisconsin breast cancer example
data(breast, package = "randomForestSRC")
breast <- na.omit(breast)

## balanced random forests - brute force
y <- breast$status
obdirect <- rfsrc(status ~ ., data = breast, nsplit = 10,
            case.wt = randomForestSRC:::make.wt(y),
            sampsize = randomForestSRC:::make.size(y))
print(obdirect)
print(get.imbalanced.performance(obdirect))

## balanced random forests - using "imbalanced" 
ob <- imbalanced(status ~ ., data = breast, method = "brf")
print(ob)
print(get.imbalanced.performance(ob))

## q-classifier (RFQ) - using "imbalanced" 
oq <- imbalanced(status ~ ., data = breast)
print(oq)
print(get.imbalanced.performance(oq))

## q-classifier (RFQ) - with auc splitting
oqauc <- imbalanced(status ~ ., data = breast, splitrule = "auc")
print(oqauc)
print(get.imbalanced.performance(oqauc))

## ------------------------------------------------------------
## unsupervised analysis
## ------------------------------------------------------------

## two equivalent ways to implement unsupervised forests
mtcars.unspv <- rfsrc(Unsupervised() ~., data = mtcars)
mtcars2.unspv <- rfsrc(data = mtcars)


## illustration of sidClustering for the mtcars data
## see sidClustering for more details
mtcars.sid <- sidClustering(mtcars, k = 1:10)
print(split(mtcars, mtcars.sid$cl[, 3]))
print(split(mtcars, mtcars.sid$cl[, 10]))


## ------------------------------------------------------------
## bivariate regression using Mahalanobis splitting
## also illustrates user specified covariance matrix
## ------------------------------------------------------------

if (library("mlbench", logical.return = TRUE)) {

  ## load boston housing data, specify the bivariate regression
  data(BostonHousing)
  f <- formula("Multivar(lstat, nox) ~.")
  
  ## Mahalanobis splitting  
  bh.mreg <- rfsrc(f, BostonHousing, importance = TRUE, splitrule = "mahal")
  
  ## performance error and vimp
  vmp <- get.mv.vimp(bh.mreg)
  pred <- get.mv.predicted(bh.mreg)

  ## standardized error and vimp
  err.std <- get.mv.error(bh.mreg, standardize = TRUE)
  vmp.std <- get.mv.vimp(bh.mreg, standardize = TRUE)

  ## same analysis, but with user specified covariance matrix
  sigma <- cov(BostonHousing[, c("lstat","nox")])
  bh.mreg2 <- rfsrc(f, BostonHousing, splitrule = "mahal", sigma = sigma)
  
}

## ------------------------------------------------------------
## multivariate mixed forests (nutrigenomic study)
## study effects of diet, lipids and gene expression for mice
## diet, genotype and lipids used as the multivariate y
## genes used for the x features
## ------------------------------------------------------------

## load the data (data is a list)
data(nutrigenomic, package = "randomForestSRC")

## assemble the multivariate y data
ydta <- data.frame(diet = nutrigenomic$diet,
                   genotype = nutrigenomic$genotype,
                   nutrigenomic$lipids)

## multivariate mixed forest call
## uses "get.mv.formula" for conveniently setting formula
mv.obj <- rfsrc(get.mv.formula(colnames(ydta)),
                data.frame(ydta, nutrigenomic$genes),
		importance=TRUE, nsplit = 10)

## print results for diet and genotype y values	    
print(mv.obj, outcome.target = "diet")
print(mv.obj, outcome.target = "genotype")

## extract standardized VIMP
svimp <- get.mv.vimp(mv.obj, standardize = TRUE)

## plot standardized VIMP for diet, genotype and lipid for each gene
boxplot(t(svimp), col = "bisque", cex.axis = .7, las = 2,
        outline = FALSE,
        ylab = "standardized VIMP",
        main = "diet/genotype/lipid VIMP for each gene")
	

## ------------------------------------------------------------
## illustrates yvar.wt which sets the probability of selecting 
## the response variables in multivariate regression
## ------------------------------------------------------------

## use mtcars: add fake responses
mult.mtcars  <- cbind(mtcars, mtcars$mpg, mtcars$mpg)
names(mult.mtcars) = c(names(mtcars), "mpg2", "mpg3")

## noise up the fake responses
mult.mtcars$mpg2 <- sample(mtcars$mpg)
mult.mtcars$mpg3 <- sample(mtcars$mpg)

formula = as.formula(Multivar(mpg, mpg2, mpg3)  ~ .)

## select 2 of the 3 responses randomly at each split with an associated weight vector.
## choose the noisy y responses which should degrade performance
yvar.wt = c(0.000001, 0.5, 0.5)
ytry = 2

mult.grow <- rfsrc(formula = formula, data = mult.mtcars, ytry = ytry, yvar.wt = yvar.wt)

print(mult.grow)
print(get.mv.error(mult.grow))

## Also, compare the following two results, as they should be similar:
yvar.wt = c(1.0, 00000.1, 00000.1)
ytry = 1

result1 = rfsrc(formula = formula, data = mult.mtcars, ytry = ytry, yvar.wt = yvar.wt)
result2 = rfsrc(mpg ~ ., mtcars)

print(get.mv.error(result1))
print(get.mv.error(result2))

## ------------------------------------------------------------
## custom splitting using the pre-coded examples
## ------------------------------------------------------------

## motor trend cars
mtcars.obj <- rfsrc(mpg ~ ., data = mtcars, splitrule = "custom")

## iris analysis
iris.obj <- rfsrc(Species ~., data = iris, splitrule = "custom1")

## WIHS analysis
wihs.obj <- rfsrc(Surv(time, status) ~ ., wihs, nsplit = 3,
                  ntree = 100, splitrule = "custom1")

}
}
\keyword{forest}
