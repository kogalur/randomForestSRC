\name{rfsrc}
\alias{rfsrc}
\alias{randomForestSRC}
\title{Random Forests for Survival, Regression and Classification (RF-SRC)}
\description{
  A random forest (Breiman, 2001) is grown using user supplied training
  data.  Applies when the response (outcome) is numeric, categorical
  (factor), or right-censored (including competing risk), and yields
  regression, classification, and survival forests, respectively.  The
  resulting forest, informally referred to as a RF-SRC object, contains
  many useful values which can be directly extracted by the user and/or
  parsed using additional functions (see the examples below).  This is
  the main entry point to the \pkg{randomForestSRC} package.
  
  Note that the package now handles multivariate regression and
  classification responses as well as mixed outcomes
  (regression/classification responses).  In such cases, a multivariate
  forest is grown, informally referred to as an mRF-SRC object.
  Unsupervised forests are also available.
  
  The package implements OpenMP shared-memory parallel programming.
  However, the default installation will only execute serially.  Users
  should consult the randomForestSRC-package help file for details on
  installing the OpenMP version of the package.  The help file is
  readily accessible via the command
  \code{package?randomForestSRC}.
}
\usage{
rfsrc(formula, data, ntree = 1000,
  bootstrap = c("by.root", "by.node", "none", "by.user"),
  mtry = NULL,
  
  nodesize = NULL,
  nodedepth = NULL,
  splitrule = NULL,
  nsplit = 0,
  split.null = FALSE,
  importance = c(FALSE, TRUE, "none", "permute", "random", "anti", 
                 "permute.ensemble", "random.ensemble",  "anti.ensemble"),
  
  na.action = c("na.omit", "na.impute"),
  
  
  nimpute = 1,
  ntime,
  cause,
  proximity = FALSE,

  sampsize = NULL,
  samptype = c("swr", "swor"),
  samp = NULL,
  case.wt  = NULL,
   

  xvar.wt = NULL,  
  forest = TRUE,
  var.used = c(FALSE, "all.trees", "by.tree"),
  split.depth = c(FALSE, "all.trees", "by.tree"),
  seed = NULL,
  do.trace = FALSE,
  membership = FALSE,
  statistics = FALSE,
  tree.err = FALSE,
  coerce.factor = NULL,
  ...)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit.  If
    missing, unsupervised splitting is implemented.}
  \item{data}{Data frame containing the y-outcome and x-variables.}
  \item{ntree}{Number of trees in the forest.}
  \item{bootstrap}{Bootstrap protocol.  The default is \code{by.root}
    which bootstraps the data by sampling with replacement at the root
    node before growing the tree.  If \code{by.node} is choosen, the
    data is bootstrapped at each node during the grow process.  If
    \code{none} is chosen, the data is not bootstrapped at all.  if
    \code{by.user} is choosen, the bootstrap specified by \code{samp} is
    used. Note the details below on prediction error when the default choice
    is not in effect.}
  \item{mtry}{Number of variables randomly selected as candidates for
    each node split. The default is sqrt(\code{p}), except for
    regression families where \code{p}/3 is used, where \code{p}
    equals the number of variables.  Values are rounded up.}
  
  \item{nodesize}{Minimum number of unique cases (data points) in a
    terminal node.  The defaults are: survival (3), competing
    risk (6), regression (5), classification (1), mixed outcomes (3).
    It is recommended to experiment with different \code{nodesize}
    values.}
  \item{nodedepth}{Maximum depth to which a tree should be grown.  The
    default behaviour is that this parameter is ignored.}
  \item{splitrule}{Splitting rule used to grow trees.  See below for
    details.}
  \item{nsplit}{Non-negative integer value.  When zero,
    deterministic splitting for an x-variable is in effect.  When
    non-zero, a maximum of nsplit split points are randomly chosen
    among the possible split points for the x-variable. This can
    significantly increase speed over deterministic splitting.  In
    the case of pure random splitting, a value of one is used as the
    default, since deterministic splitting is not a compatible
    concept in that scenario.  However, values greater than one are
    accepted, as with the other split rules.}
  \item{split.null}{Set this value to TRUE when testing the null
    hypothesis.  In particular, this assumes there is no relationship
    between y and x.}
  \item{importance}{Method for computing variable importance (VIMP).
    Calculating VIMP can be computationally expensive when the
    number of variables is high.  The defalut action is
    \code{importance="none"}.  See the details below for more
    about VIMP.}
  \item{na.action}{Action taken if the data contains \code{NA}'s.
    Possible values are \code{na.omit} or \code{na.impute}.
    The default \code{na.omit} removes the entire record if
    even one of its entries is \code{NA} (for x-variables this applies only
    to those specifically listed in 'formula').  Selecting
	  \code{na.impute} imputes the data.
	  See below for more details regarding missing data imputation.}
  \item{nimpute}{Number of iterations of the missing data algorithm.
	  Performance measures such as out-of-bag (OOB) error rates tend
	  to become optimistic if \code{nimpute} is greater than 1.}
  \item{ntime}{Integer value used for survival families to
    constrain ensemble calculations to a grid of time values of no more
    than \code{ntime} time points.  Alternatively if a vector of values
    of length greater than one is supplied, it is assumed these are the
    time points to be used to constrain the calculations (note that the
    constrained time points used will be the observed event times
    closest to the user supplied time points).  If no value is
    specified, the default action is to use all observed event times.}
  \item{cause}{Integer value between 1 and
    \code{J} indicating the event of interest for competing risks, where
    \code{J} is
	  the number of event types (this option applies only to
	  competing risks and is ignored otherwise).  While growing a
	  tree, the splitting of a node is restricted to the event type
	  specified by \code{cause}.  If not specified, the default is
	  to use a composite splitting rule which is an average over the
	  entire set of event types (a democratic approach).  Users can
	  also pass a vector of non-negative weights of length \code{J}
	  if they wish to use a customized composite split statistic
	  (for example, passing a vector of ones reverts to the default
	  composite split statistic).  In all instances when
	  \code{cause} is set incorrectly, splitting reverts to the
	  default.  Finally, note that regardless of how \code{cause} is
	  specified, the returned forest object always provides
	  estimates for all event types.}
  \item{proximity}{Proximity of cases as measured by the frequency of
    sharing the same terminal node.  This is an \code{n}x\code{n}
    matrix, which can be large.  Choices are \code{"inbag"},
    \code{"oob"}, \code{"all"}, \code{TRUE}, or \code{FALSE}.  
    Setting \code{proximity = TRUE} is equivalent to 
    \code{proximity = "inbag"}.}
  
  \item{sampsize}{Requested size of bootstrap when \code{"by.root"} is in
    effect (if missing the default action is the usual bootstrap).}
  \item{samptype}{Type of bootstrap when \code{"by.root"} is in
    effect.  Choices are \code{swr} (sampling with replacement,
    the default action) and \code{swor} (sampling without replacement).}
  \item{samp}{Bootstrap specification when \code{"by.user"} is in
    effect.  This is a array of dim \code{n x ntree} specifying how many
    times each record appears inbag in the bootstrap for each tree.}
  \item{case.wt}{Vector of non-negative weights where entry
    \code{k}, after normalizing, is the probability of selecting
    case \code{k} as a candidate for the bootstrap.  Default is to
    use uniform weights.  Vector must be of dimension \code{n},
    where \code{n} equals the number of cases in the processed data
    set (missing values may be removed, thus altering the original
    sample size).  It is generally better to use real weights rather than
    integers.  With larger sizes of \code{n}, the slightly different
    sampling algorithms used in the two scenarios can result in
    dramatically different execution times.}
   
  
  \item{xvar.wt}{Vector of non-negative weights where entry
    \code{k}, after normalizing, is the probability of selecting
    variable \code{k} as a candidate for splitting a node.  Default is to
    use uniform weights.  Vector must be of dimension \code{p},
    where \code{p} equals the number of variables, otherwise the
	  default is invoked. It is generally better to use real weights rather than
    integers.  With larger sizes of \code{p}, the slightly different
    sampling algorithms used in the two scenarios can result in
    dramatically different execution times.}
  \item{forest}{Should the forest object be returned?  Used for
    prediction on new data and required by many of the functions
	  used to parse the RF-SRC object.  It is recommended not to
    change the default setting.}
  \item{var.used}{Return variables used for splitting?  Default is
    \code{FALSE}.  Possible values are \code{all.trees} and
    \code{by.tree}.}
  \item{split.depth}{Records the minimal depth for each variable.
    Default is \code{FALSE}.  Possible values are \code{all.trees} and
    \code{by.tree}. Used for variable selection.}
  \item{seed}{Negative integer specifying seed for the random number
    generator.}
  \item{do.trace}{Number of seconds between updates to the user on
    approximate time to completion.}
  \item{membership}{Should terminal node membership and inbag
    information be returned?}
  \item{statistics}{Should split statistics be returned?  Values can be
    parsed using \command{stat.split}.}
  \item{tree.err}{Should the error rate be calculated on every tree?
	  When \code{FALSE}, it will only be calculated on the last tree.  In
    such situations, plot of the error rate will result in a flat line.
    To view the error rate over all trees, restore the model with
    option set to \code{TRUE}.}
  \item{coerce.factor}{Names of variables (can either be x-variables or
    y-outcomes) to be analyzed as factors.
    The variables must be coded as sequential positive integers
    1,2,... This option is useful in
    high-dimensional problems with a large number of factors.}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  \enumerate{
    
    \item \emph{Families}
    
    The package automagically determines the underlying random forest
    family to be fit from the following eight choices:
    
    \code{regr}, \code{regr+}, \code{class}, \code{class+}, \code{mix+},
    \code{unsupv}, \code{surv}, and \code{surv-CR}.
    
    \enumerate{
      \item Regression forests (\code{regr}) for continuous responses.
      \item Multivariate regression forests (\code{regr+})
      for multivariate continuous responses.
      \item Classification forests (\code{class}) for factor responses.
      \item Multivariate classification forests (\code{class+})
      for multivariate factor responses.
      \item Multivariate mixed forests (\code{mix+}) for mixed continuous
      and factor responses.
      \item Unsupervised forests (\code{unsupv}) when there is no response.
      \item Survival forest (\code{surv}) for right-censored survival settings.
      \item Competing risk survival forests (\code{surv-CR})
      for competing risk scenarios.
    }
    
    See below for how to code the response in the two different survival
    scenarios and for responses for multivariate forests.
    
    \item \emph{Splitrules}
    
    Splitrules are set according to the option \code{splitrule} as follows:
    \itemize{
      
      \item Regression analysis:
      \enumerate{
        \item The default rule is weighted mean-squared error splitting
        \code{mse} (Breiman et al. 1984, Chapter 8.4).
        \item Unweighted and heavy weighted mean-squared error splitting
        rules can be invoked using splitrules \code{mse.unwt} and
	      \code{mse.hvwt}.  Generally \code{mse} works best, but see
	      Ishwaran (2015) for details.
      }
      
      
      \item Multivariate regression analysis: For multivariate regression
      responses, a composite normalized mean-squared error splitting
      rule is used.
      
      \item Classification analysis:
      \enumerate{
        \item The default rule is Gini index splitting \code{gini} (Breiman
        et al. 1984, Chapter 4.3).
        \item Unweighted and heavy weighted Gini index splitting rules can
   	    be invoked using splitrules \code{gini.unwt} and
	      \code{gini.hvwt}.  Generally \code{gini} works best, but see
	      Ishwaran (2015) for details.
      }	  
      
      \item Multivariate classification analysis: For multivariate
      classification responses, a composite normalized Gini index
      splitting rule is used.
      
	    
      \item Mixed outcomes analysis: When both regression and
	    classification responses are detected, a multivariate
	    normalized composite split rule of mean-squared error and Gini
	    index splitting is invoked.
	    
      
      \item Unsupervised analysis: In settings where there is no outcome,
	    unsupervised splitting is invoked.  In this case, the mixed
	    outcome splitting rule (above) is applied.
      
      \item Survival analysis:
      \enumerate{
        \item The default rule is \code{logrank} which implements
	      log-rank splitting (Segal, 1988; Leblanc and Crowley, 1993).
        \item \code{logrankscore} implements log-rank score splitting
        (Hothorn and Lausen, 2003).
      }
      
      \item Competing risk analysis:
      \enumerate{
        
        \item The default rule is \code{logrankCR} which implements a
	      modified weighted log-rank splitting rule modeled after Gray's
	      test (Gray, 1988).
        
        \item \code{logrank} implements weighted log-rank splitting where
	      each event type is treated as the event of interest and all
	      other events are treated as censored.  The split rule is the
	      weighted value of each of log-rank statistics, standardized by
	      the variance.  For more details see Ishwaran et al. (2014).
      }

      \item Custom splitting: All families except unsupervised are
      available for user defined custom splitting.  Some basic
      C-programming skills are required.  The harness for defining these
      rules is in \code{splitCustom.c}.  In this file we give examples
      of how to code rules for regression, classification, survival, and
      competing risk.  Each family can support up to sixteen custom
      split rules.  Specifying \code{splitrule="custom"} or
      \code{splitrule="custom1"} will trigger the first split rule for
      the family defined by the training data set. Multivariate families
      will need a custom split rule for both regression and
      classification.  In the examples, we demonstrate how the user is
      presented with the node specific membership.  The task is then to
      define a split statistic based on that membership.  Take note of
      the instructions in \code{splitCustom.c} on how to \emph{register} the
      custom split rules.  It is suggested that the existing custom split
      rules be kept in place, for reference, and that the user proceed
      to develop \code{splitrule="custom2"} and so on. The package must
      be recompiled and installed for the custom split rules to become
      available.
    }
    
    
    \item \emph{Allowable data types}
    
    Data types must be real valued, integer, factor or logical -- however
    all except factors are coerced and treated as if real valued.  For
    ordered factors, splits are similar to real valued variables.  If the
    factor is unordered, a split will move a subset of the levels in the
    parent node to the left daughter, and the complementary subset to the
    right daughter.  All possible complementary pairs are considered and
    apply to factors with an unlimited number of levels.  However, there
    is an optimization check to ensure that the number of splits attempted
    is not greater than the number of cases in a node (this internal check
    will override the \code{nsplit} value in random splitting mode if
    \code{nsplit} is large enough; see below for information about
    \code{nsplit}).
    
    
    
    \item \emph{Improving computational speed}
    
    \itemize{
      
      \item \emph{Randomized Splitting Rules}
      
      Trees tend to favor splits on continuous variables and factors with
      large numbers of levels (Loh and Shih, 1997).  To mitigate this
      bias, and considerably improve computational speed, a randomize
      version of a splitting rule can be invoked using \option{nsplit}.
      If \option{nsplit} is set to a non-zero positive integer, then a
      maximum of \option{nsplit} split points are chosen randomly for each
      of the \option{mtry} variables within a node.  The splitting rule is
      applied to the random split points and the node is split on that
      variable and random split point yielding the best value (as measured
      by the splitting rule).  Pure random splitting can be invoked by
      setting \option{splitrule="random"}.  For each node, a variable is
      randomly selected and the node is split using a random split point
      (Cutler and Zhao, 2001; Lin and Jeon, 2006).
      
      The value of \option{nsplit} has a significant impact on the time
      taken to grow a forest.  This is because non-random splitting
      (i.e. the default option \option{nsplit=0}), considers all possible
      split points for each of the \option{mtry} variables, and iterating
      over each split point can be CPU intensive, especially in large
      sample size settings.
      
      In general, regardless of computational speed, it is good practice
      to use the \option{nsplit} when the data contains a mix of
      continuous and discrete variables.  Using a reasonably small value
      mitigates bias and may not compromise accuracy.
      
      
      \item \emph{Large sample size}
      
      For increased efficiency for survival families, users should consider
      setting \option{ntime} to a relatively small value when the sample
      size (number of rows of the data) is large.  This constrains ensemble
      calculations such as survival functions to a restricted grid of time
      points of length no more than \option{ntime} which can considerably
      reduce computational times.
      
      \item \emph{Large number of variables}

      For increased efficiency when the number of variables are large,
      use the defalut setting of \code{importance="none"} which turns
      off variable importance (VIMP) calculations and can considerably
      reduce computational times (see below for more details about
      variable importance).  Note that variable importance calculations
      can always be recovered later from the grow forest using the
      function \command{vimp} and/or \command{predict}.  Alternatively
      if VIMP is desired in grow mode, consider using \option{ensemble}
      VIMP which can be considerably faster, especially for survival
      families.
      
      \item \emph{Factors}
      
      For coherence, an immutable map is applied to each factor that
      ensures that factor levels in the training data set are consistent
      with the factor levels in any subsequent test data set.  This map is
      applied to each factor before and after the native C library is
      executed.  Because of this, if x-variables are all factors, then
      computational times may become very long in high dimensional
      problems.

    }
    
    \item \emph{Prediction Error}
    
    Prediction error is calculated using OOB data.  Performance is
    measured in terms of mean-squared-error for regression, and
    misclassification error for classification.  A normalized Brier score
    (relative to a coin-toss) is also provided upon printing a
    classification forest.
    
    For survival, prediction error is measured by 1-C, where C is
    Harrell's (Harrell et al., 1982) concordance index.  Prediction error
    is between 0 and 1, and measures how well the predictor correctly
    ranks (classifies) two random individuals in terms of survival.  A
    value of 0.5 is no better than random guessing. A value of 0 is
    perfect.
    
    When bootstrapping is \code{by.node} or \code{none}, a coherent OOB
    subset is not available to assess prediction error.  Thus, all outputs
    dependent on this are suppressed.  In such cases, prediction error is
    only available via classical cross-validation (the user will need to
    use the \command{predict.rfsrc} function).
    
    \item \emph{Variable Importance (VIMP)}
    
    The option \option{importance} allows several ways to calculate VIMP.
    The default \code{permute} returns Breiman-Cutler permutation VIMP as
    described in Breiman (2001).  For each tree, the prediction error on
    the out-of-bag (OOB) data is recorded.  Then for a given variable
    \emph{x}, OOB cases are randomly permuted in \emph{x} and the
    prediction error is recorded.  The VIMP for \emph{x} is defined as the
    difference between the perturbed and unperturbed error rate, averaged
    over all trees.  If \code{random} is used, then \emph{x} is not
    permuted, but rather an OOB case is assigned a daughter node randomly
    whenever a split on \emph{x} is encountered in the in-bag tree.  If
    \code{anti} is used, then \emph{x} is assigned to the opposite node
    whenever a split on \emph{x} is encountered in the in-bag tree.
    
    If \option{importance} is selected from one of the \option{ensemble}
    choices, VIMP is calculated by comparing the error rate for the
    perturbed OOB forest ensemble to the unperturbed OOB forest
    ensemble, where the perturbed ensemble is obtained by either
    permuting \emph{x}, or by random daughter node assignment, or by
    anti-splitting on \emph{x}.  Thus, unlike Breiman-Cutler VIMP,
    ensemble VIMP does not measure the tree average effect of \emph{x},
    but rather its overall forest effect (Ishwaran et al. 2008).
    Ensemble VIMP is faster to compute than Breiman-Cutler VIMP and
    therefore may be preferable in large scale problems (especially
    true for survival).
    
    Finally, the option \code{none} turns VIMP off entirely.
    
    Note that the function \command{vimp} provides a friendly user
    interface for extracting VIMP.
    
    
    \item \emph{Multivariate Forests}
    
    Multivariate forests are specified by using the multivariate formula
    interface.  Such a call can take one of two forms:
    
    rfsrc(Multivar(y1, y2, ..., yd) ~ . , my.data, ...)
    
    rfsrc(cbind(y1, y2, ..., yd) ~ . , my.data, ...)
    
    The nature of the outcomes will inform the code as to what type of
    multivariate forest to grow; i.e. whether it is real-valued,
    categorical, or a combination of both (mixed).  Note that performance
    measures such as VIMP and error rates are returned for all outcomes.
    For faster speed, VIMP can be turned off and the \command{predict}
    function used later to extract these values.
    
    
    \item \emph{Unsupervised Splitting}
    
    In the case where no y-outcomes are present, unsupervised splitting can
    be invoked by one of two means:
    
    rfsrc(Unsupervised() ~ ., data = my.data)
    
    rfsrc(data = my.data)
    
    In this case, a random subset of \code{mtry} pairs of variables are
    selected at each tree node.  For each pair, one is chosen at random to
    be the response (called the pseudo-response) and the remaining
    variable is treated as the candidate variable to be split on.  The
    best split (measured in terms of the pseudo-response) over the
    \code{mtry} pairs is used to split the node.
    
    Multivariate unsupervised splitting can also be implemented.  A
    typical call looks like:
    
    rfsrc(Unsupervised(5) ~ ., my.data, mtry = 10)
    
    In the above, \code{mtry=10} variables are selected at random, and for
    each of these a random subset of 5 variables are selected and defined
    as the multivariate pseudo-responses.  In essence this creates a
    collection of 10 multivariate regression problems (the number of
    regressions, 10, is defined by \code{mtry} and the number of
    pseudo-responses in each regression, 5, is determined by the
    unsupervised formula; informally we call this value "ytry").  A
    multivariate normalized composite splitting rule is then applied to
    each of the 10 multivariate regression problems and the node split on
    the variable leading to the best split.
    
    Note that all performance values (error rates, VIMP, prediction) are
    turned off in unsupervised mode.
    
    \item \emph{Survival, Competing Risks}
    
    \enumerate{
      \item Survival settings require a time and censoring variable which
      should be identifed in the formula as the response using the standard
      \code{Surv} formula specification.  A typical formula call looks like:
      
      Surv(my.time, my.status) ~ .
      
      where \code{my.time} and \code{my.status} are the variables names for
      the event time and status variable in the users data set.
      
      \item For survival forests (Ishwaran et al. 2008), the censoring
      variable must be coded as a non-negative integer with 0 reserved for
      censoring and (usually) 1=death (event).  The event time must be
      non-negative.
      
      \item For competing risk forests (Ishwaran et al., 2013), the
      implementation is similar to survival, but with the following
      caveats:
      
      \itemize{
        \item Censoring must be coded as a non-negative integer, where 0
        indicates right-censoring, and non-zero values indicate different
        event types.  While {0,1,2,..,J} is standard, and recommended,
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
    
    Setting \option{na.action="na.impute"} imputes missing data (both x
    and y-variables) using a modification of the missing data algorithm of
    Ishwaran et al. (2008).  Prior to splitting a node, missing data for a
    variable is imputed by randomly drawing values from non-missing in-bag
    data.  The purpose of this imputed data is to make it possible to
    assign cases to daughter nodes in the event the node is split on a
    variable with missing data.  Imputed data is however not used to
    calculate the split-statistic which uses non-missing data only.
    Following a node split, imputed data are reset to missing and the
    process is repeated until terminal nodes are reached.  Missing data in
    terminal nodes are imputed using in-bag non-missing terminal node data.
    For integer valued variables and censoring indicators, imputation uses
    a maximal class rule, whereas continuous variables and survival time
    use a mean rule. 
    
    The missing data algorithm can be iterated by setting \code{nimpute}
    to a positive integer greater than 1.  Using only a few iterations are
    needed to improve accuracy.  When the algorithm is iterated, at the
    completion of each iteration, missing data is imputed using OOB
    non-missing terminal node data which is then used as input to grow a
    new forest.  Note that when the algorithm is iterated, a side effect
    is that missing values in the returned objects \code{xvar},
    \code{yvar} are replaced by imputed values.  Further, imputed objects
    such as \code{imputed.data} are set to \code{NULL}.  Also, keep in
    mind that if the algorithm is iterated, performance measures such as
    error rates and VIMP become optimistically biased.
    
    Finally, records in which all outcome and x-variable information are
    missing are removed from the forest analysis.  Variables having all
    missing values are also removed.
    
    See the function \command{impute.rfsrc} for a fast impute interface.
    
    
    
    \item \emph{Miscellanea}
    \enumerate{
      \item Setting \option{var.used="all.trees"} returns a vector of size
      \code{p} where each element is a count of the number of times a split
      occurred on a variable.  If \option{var.used="by.tree"}, a matrix of
      size \code{ntree}x\code{p} is returned.  Each element [i][j] is the
      count of the number of times a split occurred on variable [j] in tree
      [i].
      
      \item Setting \option{split.depth="all.trees"} returns a matrix of
      size \code{n}x\code{p} where entry [i][j] is the minimal depth for
      variable [j] for case [i] averaged over the forest.  That is, for
      case [i], the entry [i][j] records the first time case [i] splits on
      variable [j] averaged over the forest.  If
      \option{split.depth="by.tree"}, a three-dimensional array is
      returned where the third dimension [k] records the tree and the
      first two coordinates [i][j] record the case and the variable.  Thus
      entry [i][j][k] is the minimal depth for case [i] for variable [j]
      for tree [k].
    }
}}
\value{
  An object of class \code{(rfsrc, grow)} with the following
  components:
  
  \item{call}{The original call to \command{rfsrc}.}
  \item{formula}{The formula used in the call.}
  \item{family}{The family used in the analysis.}
  \item{n}{Sample size of the data (depends upon \code{NA}'s, see \option{na.action}).}
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
  \item{split.wt}{Vector of non-negative weights where entry
    \code{k}, after normalizing, is the multiplier by which the
	  split statistic for a covariate is adjusted.}
  \item{leaf.count}{Number of terminal nodes for each tree in the
    forest. Vector of length \option{ntree}.  A value of zero indicates
    a rejected tree (can occur when imputing missing data).
    Values of one indicate tree stumps.}
  \item{proximity}{Proximity matrix recording the frequency of pairs of data points
    occur within the same terminal node.}
  \item{forest}{If \option{forest=TRUE}, the forest object is returned.
    This object is used for prediction with new test data
	  sets and is required for other R-wrappers.}
  
  \item{membership}{Matrix recording terminal node membership where
    each column contains the node number that a case falls in for that tree.}
  \item{inbag}{Matrix recording inbag membership where each column
    contains the number of times that a case appears in the bootstrap
    sample for that tree.}
  \item{var.used}{Count of the number of times a variable is used in
    growing the forest.}
  \item{imputed.indv}{Vector of indices for cases with missing
    values.}
  \item{imputed.data}{Data frame of the imputed data. The first
	  column(s) are reserved for the y-responses, after which the
	  x-variables are listed.}
  \item{split.depth}{Matrix [i][j] or array [i][j][k] recording the
    minimal depth for variable [j] for case [i], either averaged over
    the forest, or by tree [k].}
  \item{node.stats}{Split statistics returned when
    \code{statistics=TRUE} which can be parsed using \command{stat.split}.}
  \item{err.rate}{Tree cumulative OOB error rate.}
  \item{importance}{Variable importance (VIMP) for each x-variable.}
  \item{predicted}{In-bag predicted value.}
  \item{predicted.oob}{OOB predicted value.}\cr
  
  \item{++++++++}{for classification settings, additionally ++++++++} \cr
  \item{class}{In-bag predicted class labels.}
  \item{class.oob}{OOB predicted class labels.}\cr
  
  
  \item{++++++++}{for multivariate settings, additionally ++++++++} \cr
  
  \item{regrOutput}{List containing performance values for
	  regression responses (if they are present).}
	\item{clasOutput}{List containing performance values for
	  categorical (factor) responses (if they are present).}
  
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
  \item{time.interest}{Ordered unique event times.}
  \item{ndead}{Number of events.}
}
\note{
  \enumerate{
    \item Values returned depend heavily on the family.  In
    particular, \code{predicted} and \code{predicted.oob} are the
    following values calculated using in-bag and OOB data:
    
    \enumerate{
      \item For regression, a vector of predicted y-responses.
      
      \item For classification, a matrix with columns containing the
      estimated class probability for each class.
      
      \item For survival, a vector of mortality values (Ishwaran et al.,
      2008) representing estimated risk for each individual calibrated
      to the scale of the number of events (as a specific example, if
      \emph{i} has a mortality value of 100, then if all individuals had
      the same x-values as \emph{i}, we would expect an average of 100
      events).  Also included in the grow object are matrices containing
      the CHF and survival function.  Each row corresponds to an
      individual's ensemble CHF or survival function evaluated at each
      time point in \code{time.interest}.
      
      \item For competing risks, a matrix with one column for each event
      recording the expected number of life years lost due to the event
      specific cause up to the maximum follow up (Ishwaran et al.,
      2013).  The grow object also contains the cause-specific
      cumulative hazard function (CSCHF) and the cumulative incidence
      function (CIF) for each event type.  These are encoded as a
      three-dimensional array, with the third dimension used for the
      event type, each time point in \code{time.interest} making up the
      second dimension (columns), and the case (individual) being the
      first dimension (rows).
      
      \item For multivariate families, predicted values (and other
      performance values such as VIMP and error rates) are stored in 
      the lists \code{regrOutput} and \code{clasOutput}.
      
    }
    
    \item Different R-wrappers are provided to aid in parsing the grow
    object.
  }
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Breiman L., Friedman J.H., Olshen R.A. and Stone C.J.
  \emph{Classification and Regression Trees}, Belmont, California, 1984.
  
  Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.
  
  Cutler A. and Zhao G. (2001). Pert-Perfect random tree ensembles.
  \emph{Comp. Sci. Statist.}, 33: 490-497.
  
  Gray R.J. (1988).  A class of k-sample tests for comparing the
  cumulative incidence of a competing risk, \emph{Ann. Statist.},
  16: 1141-1154.
  
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
  
  Ishwaran H., Gerds T.A., Kogalur U.B., Moore R.D., Gange S.J. and Lau
  B.M. (2014). Random survival forests for competing risks.
  \emph{Biostatistics}, 15(4):757-773.
  
  Ishwaran H. (2015).  The effect of splitting on random forests.
  \emph{Machine Learning}, 99:75-118.
  
  Lin Y. and Jeon Y. (2006).  Random forests and adaptive nearest
  neighbors, \emph{J. Amer. Statist. Assoc.}, 101:578-590.
  
  LeBlanc M. and Crowley J. (1993).  Survival trees by goodness of split,
  \emph{J. Amer. Statist. Assoc.}, 88:457-467.
  
  Loh W.-Y and Shih Y.-S (1997).  Split selection methods for
  classification trees, \emph{Statist. Sinica}, 7:815-840.
  
  Mogensen, U.B, Ishwaran H. and Gerds T.A. (2012). Evaluating random
  forests for survival analysis using prediction error curves,
  \emph{J.  Statist. Software}, 50(11): 1-23.
  
  Segal M.R. (1988).  Regression trees for censored data,
  \emph{Biometrics},  44:35-47.
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
  \command{\link{stat.split}},
  \command{\link{var.select}},
  \command{\link{vimp}}
}
\examples{
##------------------------------------------------------------
## Survival analysis
##------------------------------------------------------------

## veteran data
## randomized trial of two treatment regimens for lung cancer
data(veteran, package = "randomForestSRC")
v.obj <- rfsrc(Surv(time, status) ~ ., data = veteran, ntree = 100, tree.err=TRUE)

# print and plot the grow object
print(v.obj)
plot(v.obj)

# plot survival curves for first 10 individuals: direct way
matplot(v.obj$time.interest, 100 * t(v.obj$survival[1:10, ]),
    xlab = "Time", ylab = "Survival", type = "l", lty = 1)

# plot survival curves for first 10 individuals
# indirect way: using plot.survival (also generates hazard plots)
plot.survival(v.obj, subset = 1:10, haz.model = "ggamma")

\dontrun{
## Primary biliary cirrhosis (PBC) of the liver

data(pbc, package = "randomForestSRC")
pbc.obj <- rfsrc(Surv(days, status) ~ ., pbc, nsplit = 10)
print(pbc.obj)


##------------------------------------------------------------
## Example of imputation in survival analysis
##------------------------------------------------------------

data(pbc, package = "randomForestSRC")
pbc.obj2 <- rfsrc(Surv(days, status) ~ ., pbc,
           nsplit = 10, na.action = "na.impute")


# here's a nice wrapper to combine original data + imputed data
combine.impute <- function(object) {
 impData <- cbind(object$yvar, object$xvar)
 if (!is.null(object$imputed.indv)) {
   impData[object$imputed.indv, ] <- object$imputed.data
 }
 impData
}

# combine original data + imputed data
pbc.imp.data <- combine.impute(pbc.obj2)

# same as above but we iterate the missing data algorithm
pbc.obj3 <- rfsrc(Surv(days, status) ~ ., pbc, nsplit=10,
         na.action = "na.impute", nimpute = 3)
pbc.iterate.imp.data <- combine.impute(pbc.obj3)

# fast way to impute the data (no inference is done)
# see impute.rfsc for more details
pbc.fast.imp.data <- impute.rfsrc(data = pbc, nsplit = 10, nimpute = 5)

##------------------------------------------------------------
## Compare RF-SRC to Cox regression
## Illustrates C-index and Brier score measures of performance
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
  cox.obj <- coxph(surv.f, data = pbc.na)
  rfsrc.obj <- rfsrc(surv.f, pbc.na, nsplit = 10, ntree = 150)

  ## compute bootstrap cross-validation estimate of expected Brier score
  ## see Mogensen, Ishwaran and Gerds (2012) Journal of Statistical Software
  set.seed(17743)
  prederror.pbc <- pec(list(cox.obj,rfsrc.obj), data = pbc.na, formula = pec.f,
                        splitMethod = "bootcv", B = 50)
  print(prederror.pbc)
  plot(prederror.pbc)

  ## compute out-of-bag C-index for cox regression and compare to rfsrc
  rfsrc.obj <- rfsrc(surv.f, pbc.na, nsplit = 10)
  cat("out-of-bag Cox Analysis ...", "\n")
  cox.err <- sapply(1:100, function(b) {
    if (b\%\%10 == 0) cat("cox bootstrap:", b, "\n")
    train <- sample(1:nrow(pbc.na), nrow(pbc.na), replace = TRUE)
    cox.obj <- tryCatch({coxph(surv.f, pbc.na[train, ])}, error=function(ex){NULL})
    if (is.list(cox.obj)) {
      randomForestSRC:::cindex(pbc.na$days[-train],
                                      pbc.na$status[-train],
                                      predict(cox.obj, pbc.na[-train, ]))
    } else NA
  })
  cat("\n\tOOB error rates\n\n")
  cat("\tRSF            : ", rfsrc.obj$err.rate[rfsrc.obj$ntree], "\n")
  cat("\tCox regression : ", mean(cox.err, na.rm = TRUE), "\n")
}

##------------------------------------------------------------
## Competing risks
##------------------------------------------------------------

## WIHS analysis
## cumulative incidence function (CIF) for HAART and AIDS stratified by IDU

data(wihs, package = "randomForestSRC")
wihs.obj <- rfsrc(Surv(time, status) ~ ., wihs, nsplit = 3, ntree = 100)
plot.competing.risk(wihs.obj)
cif <- wihs.obj$cif
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
  pbc.cr <- rfsrc(Surv(time, status) ~ ., pbc, nsplit = 10)

  ## log-rank event-one specific splitting
  pbc.log1 <- rfsrc(Surv(time, status) ~ ., pbc, nsplit = 10,
              splitrule = "logrank", cause = c(1,0), importance="permute")

  ## log-rank event-two specific splitting
  pbc.log2 <- rfsrc(Surv(time, status) ~ ., pbc, nsplit = 10,
              splitrule = "logrank", cause = c(0,1), importance="permute")

  ## extract VIMP from the log-rank forests: event-specific
  ## extract minimal depth from the Gray log-rank forest: non-event specific
  var.perf <- data.frame(md = max.subtree(pbc.cr)$order[, 1],
                         vimp1 = 100 * pbc.log1$importance[ ,1],
                         vimp2 = 100 * pbc.log2$importance[ ,2])
  print(var.perf[order(var.perf$md), ])

}



## ------------------------------------------------------------
## Regression analysis
## ------------------------------------------------------------

## New York air quality measurements
airq.obj <- rfsrc(Ozone ~ ., data = airquality, na.action = "na.impute")

# partial plot of variables (see plot.variable for more details)
plot.variable(airq.obj, partial = TRUE, smooth.lines = TRUE)

## motor trend cars
mtcars.obj <- rfsrc(mpg ~ ., data = mtcars)

# minimal depth variable selection via max.subtree
md.obj <- max.subtree(mtcars.obj)
cat("top variables:\n")
print(md.obj$topvars)

# equivalent way to select variables
# see var.select for more details
vs.obj <- var.select(object = mtcars.obj)


## ------------------------------------------------------------
## Classification analysis
## ------------------------------------------------------------

## Edgar Anderson's iris data
iris.obj <- rfsrc(Species ~., data = iris)

## Wisconsin prognostic breast cancer data
data(breast, package = "randomForestSRC")
breast.obj <- rfsrc(status ~ ., data = breast, nsplit = 10, tree.err=TRUE)
plot(breast.obj)

## ------------------------------------------------------------
## Unsupervised analysis
## ------------------------------------------------------------

# two equivalent ways to implement unsupervised forests
mtcars.unspv <- rfsrc(Unsupervised() ~., data = mtcars)
mtcars2.unspv <- rfsrc(data = mtcars)

#minimal depth variable selection applies!
var.select(mtcars2.unspv)

## ------------------------------------------------------------
## Multivariate regression analysis
## ------------------------------------------------------------
mtcars.mreg <- rfsrc(Multivar(mpg, cyl) ~., data = mtcars, tree.err=TRUE)
print(mtcars.mreg, outcome.target = "mpg")
print(mtcars.mreg, outcome.target = "cyl")
plot(mtcars.mreg, outcome.target = "mpg")
plot(mtcars.mreg, outcome.target = "cyl")


## ------------------------------------------------------------
## Mixed outcomes analysis
## ------------------------------------------------------------
mtcars.new <- mtcars
mtcars.new$cyl <- factor(mtcars.new$cyl)
mtcars.new$carb <- factor(mtcars.new$carb, ordered = TRUE)
mtcars.mix <- rfsrc(cbind(carb, mpg, cyl) ~., data = mtcars.new, tree.err=TRUE)
print(mtcars.mix, outcome.target = "mpg")
print(mtcars.mix, outcome.target = "cyl")
plot(mtcars.mix, outcome.target = "mpg")
plot(mtcars.mix, outcome.target = "cyl")


## ------------------------------------------------------------
## Custom splitting using the pre-coded examples.
## ------------------------------------------------------------
## motor trend cars
mtcars.obj <- rfsrc(mpg ~ ., data = mtcars, splitrule="custom")
## Edgar Anderson's iris data
iris.obj <- rfsrc(Species ~., data = iris, splitrule="custom1")
## WIHS analysis
wihs.obj <- rfsrc(Surv(time, status) ~ ., wihs, nsplit = 3,
                  ntree = 100, splitrule="custom1")


}
}
\keyword{forest}
