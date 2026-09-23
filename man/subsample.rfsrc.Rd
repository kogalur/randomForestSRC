\name{subsample.rfsrc}
\alias{subsample}
\alias{subsample.rfsrc}
\alias{extract.subsample}
\alias{extract.bootsample}
\alias{print.subsample}
\alias{print.subsample.rfsrc}
\alias{print.bootsample}
\alias{print.bootsample.rfsrc}
\title{Subsampling Inference for Variable Importance and Prediction Error}
\description{
  Estimate standard errors and confidence intervals for variable
  importance (VIMP) and out-of-bag prediction error by fitting forests
  to repeated subsamples of the training data. Extract, print, or plot
  the results at different confidence levels from the same subsample
  estimates. An optional double-bootstrap procedure is also available.
}
\usage{\method{subsample}{rfsrc}(obj,
  B = 100,
  block.size = 1,
  importance,
  subratio = NULL,
  stratify = TRUE,
  performance = FALSE,
  performance.only = FALSE,
  joint = FALSE,
  xvar.names = NULL,
  bootstrap = FALSE,
  verbose = TRUE)

extract.subsample(obj, alpha = .05, target = 0, m.target = NULL,
                  standardize = TRUE, raw = FALSE)

extract.bootsample(obj, alpha = .05, target = 0, m.target = NULL,
                   standardize = TRUE, raw = FALSE)

\method{print}{subsample}(x, alpha = .05, standardize = TRUE, ...)

\method{print}{bootsample}(x, alpha = .05, standardize = TRUE, ...)
}

\arguments{
  \item{obj}{For \code{subsample}, a forest grow object of class
    \code{(rfsrc, grow)} with saved forest information and training data.
    Grow the forest with \code{forest = TRUE} and
    \code{bootstrap = "by.root"}. Anonymous forests and nonuniform
    case weights are not supported. For \code{extract.subsample}, an
    object returned by \code{subsample}, including double-bootstrap
    output. For \code{extract.bootsample}, an object returned by
    \code{subsample(..., bootstrap = TRUE)}.}
  \item{B}{Number of subsamples, or double-bootstrap replicates when
    \code{bootstrap = TRUE}. Each replicate requires a new forest fit.
    Supply a positive integer; at least two replicates are needed to
    calculate a sample standard deviation.}
  \item{block.size}{Number of trees per block for calculating VIMP
    when it must first be computed for \code{obj}. If \code{obj} already
    contains VIMP, its block size is used instead.}
  \item{importance}{VIMP perturbation method, such as \code{"anti"},
    \code{"permute"}, or \code{"random"}. For subsampling, existing VIMP
    in \code{obj} determines the method. If VIMP is absent, it is
    computed using this argument, or \code{importance = TRUE} when
    omitted. Ignored when \code{performance.only = TRUE}.}
  \item{subratio}{Fraction of training observations to draw for each
    subsample. Supply a finite number greater than zero and at most one;
    for inference, choose a fraction below one that gives enough
    observations to fit each forest. The default is
    \eqn{\exp(-1)} for \eqn{n\leq1000} and
    \eqn{\exp(-1)\sqrt{1000/n}} for larger \eqn{n}.
    This controls the data set for each forest, not its per-tree
    \code{sampsize}. Ignored when \code{bootstrap = TRUE}.}
  \item{stratify}{Sample separately within response classes for
    univariate classification, or within event-status groups, including
    censoring, for survival and competing risks. Regression,
    multivariate, and mixed-outcome forests use unstratified subsampling.
    Ignored when \code{bootstrap = TRUE}.}
  \item{performance}{Include uncertainty estimates for OOB prediction
    error as well as VIMP. The error measure is determined by the
    original forest's \code{perf.type}.}
  \item{performance.only}{Compute uncertainty for OOB prediction error
    only. Sets \code{performance = TRUE} and turns off individual and
    joint VIMP.}
  \item{joint}{Include joint VIMP for the predictor group specified by
    \code{xvar.names}, in addition to individual VIMP. Joint VIMP uses
    the default anti perturbation of \code{vimp}, independently of the
    method used for individual VIMP.}
  \item{xvar.names}{Predictors to perturb together for joint VIMP.
    The default \code{NULL} uses all predictors. Individual VIMP is
    still computed for all predictors. To select results for display,
    use \code{xvar.names} in \code{plot.subsample}.}
  \item{bootstrap}{Use the double-bootstrap procedure instead of
    subsampling without replacement. See Details.}
  \item{verbose}{Print progress information during resampling.}
  \item{alpha}{Significance level for extraction or printing, strictly
    between zero and one. Intervals have nominal confidence level
    \eqn{1-\alpha}; the default \code{.05} gives 95 percent intervals.}
  \item{target}{Statistic to extract within the selected response.
    Supply a nonnegative integer: \code{0} for regression or
    right-censored survival; \code{0} for overall classification or
    \code{1}, \code{2}, and so on for its classes; and \code{0},
    \code{1}, and so on for competing-risk events. Classes and events
    follow the order in the forest. For competing risks, the first event
    uses \code{target = 0} here but \code{target = 1} in
    \code{plot.subsample}.}
  \item{m.target}{Name of one response to extract from a multivariate
    or mixed-outcome forest. The default selects the first regression
    response with results, or the first classification response with
    results if none is available for regression. Not needed for
    univariate forests.}
  \item{standardize}{For regression, divide estimates by the response
    variance in the full training data. This applies to VIMP, prediction
    error, and their uncertainty summaries. Other families are unchanged.
    Set \code{FALSE} for the original scale.}
  \item{raw}{Include replicate estimates and interval matrices in the
    extraction result. This controls the amount of output, not its
    scale; use \code{standardize = FALSE} for unscaled values.}
  \item{x}{An object returned by \code{subsample}. Use
    \code{print.subsample} for subsampling and \code{print.bootsample}
    for double-bootstrap output.}
  \item{...}{Additional arguments to the print methods; currently unused.}
}
\details{
  Use \code{subsample()} to add standard errors and confidence intervals
  to a forest's VIMP estimates. The original forest provides the
  full-data estimate for each predictor. Forests grown on repeated
  subsamples show how these estimates vary when the training
  observations change. The same approach can estimate uncertainty in
  OOB prediction error.

  The analysis has two stages: \code{subsample()} fits the replicate
  forests, then \code{extract.subsample()}, \code{print()}, or
  \code{plot.subsample()} uses their saved estimates to construct
  intervals. Once resampling is complete, you can change the confidence
  level, interval method, or response without fitting more forests.

  \subsection{Choosing the statistics}{
    By default, the analysis computes uncertainty for individual VIMP.
    Specify the desired VIMP method in the original grow call, for
    example \code{importance = "permute"}. When \code{obj} already
    contains VIMP, subsampling uses that method and block size, even
    if different values are supplied to \code{subsample()}. Otherwise,
    VIMP is first calculated on the saved trees using the requested
    \code{importance} and \code{block.size}.

    Two additional statistics can be included:
    \describe{
      \item{Joint VIMP.}{\code{joint = TRUE} adds a row named
        \code{joint}, measuring the effect of perturbing the predictors
        in \code{xvar.names} together. Individual VIMP is also included.}
      \item{Prediction error.}{\code{performance = TRUE} adds a row
        named \code{err} for OOB prediction error. Use
        \code{performance.only = TRUE} to compute only this row.}
    }
  }
  \subsection{How subsampling works}{
    Each replicate draws training observations without replacement and
    grows a new forest on that subset. \code{B} controls the number of
    replicates; \code{subratio} controls their nominal sample size.
    Increase \code{B} for more stable interval endpoints, especially at
    high confidence levels.

    Subsampling the data and sampling observations for individual trees
    are separate operations. \code{subratio} determines the data set
    available to each new forest. The original forest's \code{sampsize}
    function is then applied to that data set's size to determine its
    per-tree sample size. The forest settings used for these fits are
    listed in the Note.

    With \code{stratify = TRUE}, univariate classification samples
    within classes, and survival and competing risks sample within
    event-status groups. Integer allocation can leave small groups
    unrepresented; missing group labels are excluded. Other families
    sample without stratification.

    For survival and competing risks, a subsample must contain at least
    one observed event after missing-value omission. Otherwise the
    procedure stops and asks for a larger \code{subratio}; it does not
    draw a replacement subsample.
  }
  \subsection{Extracting and displaying results}{
    After \code{smp <- subsample(obj)},
    \code{extract.subsample(smp, alpha = .05)} returns all three
    subsampling interval summaries for one response: nonparametric,
    normal subsampling, and normal jackknife. For multivariate or mixed
    outcomes, use \code{m.target} to choose the response. For
    classification or competing risks, \code{target} selects its
    overall, class-specific, or event-specific statistic as described
    under Arguments.

    \code{print(smp, alpha = .05)} displays intervals for every response
    and statistic column, rounded to three decimal places.
    \code{plot.subsample(smp, alpha = .05)} draws the normal subsampling
    intervals; its \code{normal} and \code{jknife} arguments select the
    other methods. Use extraction for numerical results at full
    precision, and \code{raw = TRUE} for the replicate estimates and
    interval matrices listed under Value.

    Extraction and printing default to 95 percent intervals
    (\code{alpha = .05}); plotting defaults to 99 percent
    (\code{alpha = .01}). Supply the same \code{alpha} when comparing
    them. Intervals are calculated separately for each statistic,
    without a multiple-testing adjustment.
  }
  \subsection{Standardization}{
    Extraction, printing, and plotting use \code{standardize = TRUE}
    by default. For regression, the full-data and replicate estimates
    are divided by the response variance in the full training data.
    This same divisor is used for every replicate and for both VIMP
    and prediction error. Confidence limits and standard errors are
    then on the corresponding standardized scale.

    Use \code{standardize = FALSE} for the original scale. Classification
    and survival statistics are unchanged by this option. Regression
    standardization requires a finite, positive response variance.
  }
  \subsection{Subsampling intervals}{
    For one statistic, let \eqn{\widehat\theta} be its full-data estimate
    and \eqn{\theta_b^*} its estimate in subsample \eqn{b}. Let \eqn{n}
    be the training sample size, \eqn{r} the value of \code{subratio},
    and \eqn{m=nr} the nominal subsample size. Write \eqn{\bar\theta_*}
    and \eqn{s_*} for the mean and sample standard deviation of the
    available subsample estimates.

    \describe{
      \item{Normal subsampling intervals.}{The subsampling standard
        error measures dispersion about the subsample mean:
        \deqn{\mathrm{SE}_{\mathrm{sub}}=\sqrt{r}\,s_*.}
        It is returned as \code{se.Z}. The interval is
        \eqn{\widehat\theta\pm z_{1-\alpha/2}\mathrm{SE}_{\mathrm{sub}}}
        and is reported in \code{var.sel.Z}.}
      \item{Normal jackknife intervals.}{A subsample of size \eqn{m}
        also deletes \eqn{d=n-m} observations. The delete-\eqn{d}
        jackknife calculation measures dispersion about the full-data
        estimate:
        \deqn{\mathrm{SE}_{\mathrm{jk}}=
          \left\{r\,\mathrm{mean}_b
          (\theta_b^*-\widehat\theta)^2\right\}^{1/2}.}
        It is returned as \code{se.jk.Z}. The interval is
        \eqn{\widehat\theta\pm z_{1-\alpha/2}\mathrm{SE}_{\mathrm{jk}}}
        and is reported in \code{var.jk.sel.Z}.}
      \item{Nonparametric subsampling intervals.}{These use the
        empirical distribution of the centered roots
        \eqn{Z_b^*=\sqrt{m}(\theta_b^*-\widehat\theta)}. If \eqn{q_p}
        is their empirical quantile, computed with the default
        \code{quantile} interpolation, the endpoints are
        \deqn{\left[\widehat\theta-q_{1-\alpha/2}/\sqrt{n},
          \widehat\theta-q_{\alpha/2}/\sqrt{n}\right].}
        The interval is reported in \code{var.sel}. Centering, scaling,
        and reversing the root quantiles converts subsample variation
        to an interval for the full-data statistic; this is not a
        percentile interval of the unscaled VIMP replicates.}
    }

    The jackknife calculation includes both subsample variation and
    the difference between the subsample mean and the full-data estimate:
    \deqn{\mathrm{mean}_b(\theta_b^*-\widehat\theta)^2
      =\mathrm{mean}_b(\theta_b^*-\bar\theta_*)^2
      +(\bar\theta_*-\widehat\theta)^2.}
    The final term is the nonnegative bias adjustment discussed by
    Ishwaran and Lu (2019), equations (10) and (11). The Note explains
    the scaling and finite-replicate denominators used here.
  }
  \subsection{Interpreting the summary tables}{
    Each table gives interval endpoints and a \code{signif} column
    indicating whether the lower endpoint exceeds zero. For VIMP,
    this identifies positive importance under the selected interval
    method. An \code{err} row instead describes uncertainty in prediction
    error; its positive lower endpoint is not a variable-selection
    result. The normal tables report the full-data estimate in
    \code{mean}; the nonparametric table reports an interval median.

    Normal-table p-values are one-sided upper-tail tests of zero VIMP.
    A positive lower endpoint of a two-sided \eqn{1-\alpha} normal
    interval corresponds to a one-sided level of \eqn{\alpha/2}.
    The nonparametric p-value is the proportion of available roots
    satisfying \eqn{Z_b^*>\sqrt{n}\,\widehat\theta}. Ties are excluded,
    and no Monte Carlo adjustment is added. Because interval endpoints
    use interpolated quantiles, the empirical p-value and interval can
    give different decisions at finite-sample boundaries.
  }
  \subsection{Double-bootstrap inference}{
    Set \code{bootstrap = TRUE} to draw \eqn{n} observations with
    replacement for each outer replicate, then generate the per-tree
    samples from that bootstrap data set. With \code{samptype = "swor"},
    trees sample without replacement from its rows; an original
    observation can still appear more than once because the bootstrap
    data set can contain duplicates.

    For these objects, \code{extract.subsample()} calls
    \code{extract.bootsample()}. It returns percentile intervals and
    normal intervals centered at the bootstrap mean, using the
    bootstrap standard deviation. Printing displays both methods.
    The percentile table has no p-value column; the centered-root
    p-value above applies only to subsampling without replacement.
  }
}
\note{
  \subsection{Forest settings and missing data}{
    Subsample forests use the original formula, number of trees,
    \code{mtry}, node-size and depth settings, splitting rule,
    \code{nsplit}, predictor and split weights, sampling-size function,
    sampling type, and performance settings. Survival forests also use
    the saved time grid and choice of Uno weighting. Specialized grow
    options outside these settings are not reproduced.

    Each new fit omits rows with missing values rather than repeating
    the original imputation. Missing-value omission and integer
    allocation within strata can reduce the actual subsample size.
    Interval scaling uses the nominal size \eqn{m=nr}, where \eqn{r}
    is \code{subratio}; actual sample counts are not saved in the result.
  }
  \subsection{Variance scaling}{
    \code{se.jk.Z} uses the factor \eqn{m/n}; equation (11) of Ishwaran
    and Lu (2019) uses \eqn{m/(n-m)}. These factors are asymptotically
    equivalent when \eqn{m/n\to0}. With \eqn{K} available estimates,
    \code{se.Z} uses the sample-variance denominator \eqn{K-1}, whereas
    \code{se.jk.Z} averages squared deviations with denominator \eqn{K}.
    Consequently, the returned jackknife standard error can be smaller
    than the subsampling standard error for a finite set of replicates.
  }
}
\value{
  \code{subsample()} returns the reference forest and replicate
  estimates. Use extraction, printing, or plotting to obtain intervals.
  Its class includes \code{subsample} for subsampling or
  \code{bootsample} for double bootstrapping, together with the
  original forest classes.
  \item{rf}{Reference forest with full-data statistics. If VIMP was
    initially absent, it has been calculated on the saved trees.}
  \item{vmp}{List of full-data statistic matrices, one per outcome.
    Rows correspond to predictors and any requested \code{joint} or
    \code{err} statistic. Classification and competing risks can have
    several columns. Unavailable appended statistics are \code{NA}
    where their outcome and column can be matched. These values are
    not standardized.}
  \item{vmpS}{For subsampling, a length-\code{B} list of replicate
    results, each organized by outcome. The tables use the outcomes,
    row names, and column names in \code{vmp}; unavailable entries are
    \code{NA}.}
  \item{vmpB}{Double-bootstrap replicate estimates, in place of
    \code{vmpS}.}
  \item{subratio}{Nominal subsample fraction, or \code{NULL} for
    double-bootstrap objects.}
  \item{performance.only}{Whether only prediction error was requested.}

  \subsection{Extraction helpers}{
    Each helper returns a list for one response and statistic column.
    Summary tables have one row per predictor or requested joint/error
    statistic. All extracted values use the requested
    \code{standardize} setting.
    \describe{
      \item{\code{extract.subsample()}}{For subsampling, returns
        \code{se.Z} and \code{se.jk.Z}, the subsampling and jackknife
        standard errors; \code{var.sel}, the nonparametric table;
        \code{var.sel.Z}, the normal subsampling table; and
        \code{var.jk.sel.Z}, the normal jackknife table. For a
        double-bootstrap object, it returns the
        \code{extract.bootsample()} result.}
      \item{\code{extract.bootsample()}}{Returns \code{se}, the
        bootstrap standard errors; \code{var.sel}, the percentile
        table; and \code{var.sel.Z}, the normal table. Both tables
        contain \code{lower}, \code{mean}, \code{upper}, and
        \code{signif}; only \code{var.sel.Z} includes \code{pvalue}.}
    }
    With \code{raw = TRUE}, the helpers also return:
    \describe{
      \item{\code{vmp}}{Full-data estimates for subsampling, or
        bootstrap means for double bootstrapping.}
      \item{\code{vmpS}}{Matrix with statistics in rows and replicates
        in columns, for either procedure.}
      \item{\code{ci}, \code{ci.Z}, \code{ci.jk.Z}}{Five-row interval
        matrices for the nonparametric, normal subsampling, and normal
        jackknife methods. Rows give the lower endpoint, lower
        quartile, center, upper quartile, and upper endpoint.
        Double-bootstrap output includes only \code{ci} and
        \code{ci.Z}, both with the bootstrap mean in the middle row.}
      \item{\code{boxplot.dta}}{Centered and rescaled replicate values
        used for subsampling plots. Not returned for double bootstrapping.}
    }
  }

  The print methods display interval summaries. Use the extraction
  helpers for numerical results. \code{print.subsample.rfsrc} and
  \code{print.bootsample.rfsrc} have the same arguments as
  \code{print.subsample} and \code{print.bootsample}, respectively.
}
\author{Hemant Ishwaran and Udaya B. Kogalur}
\references{
  Ishwaran H. and Lu M. (2019). Standard errors and confidence
  intervals for variable importance in random forest regression,
  classification, and survival. \emph{Statistics in Medicine}, 38,
  558-582.

  Geyer, C.J. (2013). 5601 Notes: The Subsampling Bootstrap.
  July 11, 2013. Section 5, equation (12).

  Politis, D.N. and Romano, J.P. (1994). Large sample confidence
  regions based on subsamples under minimal assumptions. \emph{The
  Annals of Statistics}, 22(4):2031-2050.

  Shao, J. and Wu, C.J. (1989). A general theory for jackknife variance
  estimation. \emph{The Annals of Statistics}, 17(3):1176-1197.
}
\seealso{
  \code{\link{plot.subsample.rfsrc}},
  \code{\link{rfsrc}},
  \code{\link{vimp.rfsrc}},
  \code{\link{holdout.vimp.rfsrc}}
}
\examples{
\donttest{
## Small settings keep these examples manageable. Increase B for more
## stable interval endpoints, especially at high confidence levels.
set.seed(19)
dta <- na.omit(airquality)

## ------------------------------------------------------------
## Regression: grow once, then reuse the saved replicate estimates
## ------------------------------------------------------------
reg.o <- rfsrc(Ozone ~ ., data = dta, ntree = 100,
               importance = "permute", block.size = 1)
reg.smp <- subsample(reg.o, B = 25, verbose = FALSE)

ans <- extract.subsample(reg.smp, alpha = .05)
print(ans$var.sel.Z)
print(ans$var.jk.sel.Z)
print(reg.smp, alpha = .05)
plot.subsample(reg.smp, alpha = .05)
plot.subsample(reg.smp, alpha = .05, normal = FALSE)

## Detailed output and unstandardized VIMP are separate requests.
raw <- extract.subsample(reg.smp, alpha = .05,
                         standardize = FALSE, raw = TRUE)
print(raw$ci.Z)
print(dim(raw$vmpS))

## ------------------------------------------------------------
## Joint VIMP and prediction-error uncertainty
## ------------------------------------------------------------
with.error <- subsample(reg.o, B = 25, performance = TRUE,
                         joint = TRUE, xvar.names = c("Day", "Month"),
                         verbose = FALSE)
ans <- extract.subsample(with.error, alpha = .05)
print(ans$var.sel.Z[c("joint", "err"), , drop = FALSE])

error.only <- subsample(reg.o, B = 25, performance.only = TRUE,
                         verbose = FALSE)
print(extract.subsample(error.only, alpha = .05)$var.sel.Z)

## ------------------------------------------------------------
## Multivariate regression: select an outcome after subsampling
## ------------------------------------------------------------
mv <- rfsrc(cbind(Ozone, Temp) ~ ., data = dta, ntree = 100,
             importance = "permute", block.size = 1)
mv.smp <- subsample(mv, B = 25, verbose = FALSE)
print(extract.subsample(mv.smp, m.target = "Temp", alpha = .05)$var.sel.Z)
plot.subsample(mv.smp, m.target = "Ozone", alpha = .05)
plot.subsample(mv.smp, m.target = "Temp", alpha = .05,
                jknife = TRUE)

## ------------------------------------------------------------
## Classification with a different performance measure
## ------------------------------------------------------------
if (requireNamespace("caret", quietly = TRUE)) {
  cls.data <- caret::twoClassSim(500, linearVars = 5, noiseVars = 5)
  cls <- rfsrc(Class ~ ., data = cls.data, ntree = 100,
                importance = "permute", block.size = 1,
                perf.type = "brier")
  cls.smp <- subsample(cls, B = 25, verbose = FALSE)
  print(extract.subsample(cls.smp, alpha = .05)$var.sel.Z)
  plot.subsample(cls.smp, alpha = .05, cex.axis = .7)
}

## ------------------------------------------------------------
## Survival
## ------------------------------------------------------------
data(pbc, package = "randomForestSRC")
pbc.complete <- na.omit(pbc)
srv <- rfsrc(Surv(days, status) ~ ., data = pbc.complete,
               ntree = 100, importance = "permute", block.size = 1)
srv.smp <- subsample(srv, B = 25, verbose = FALSE)
plot.subsample(srv.smp, alpha = .05)

## ------------------------------------------------------------
## Optional double bootstrap (separate resampling procedure)
## ------------------------------------------------------------
reg.boot <- subsample(reg.o, B = 25, bootstrap = TRUE, verbose = FALSE)
print(reg.boot, alpha = .05)
print(extract.bootsample(reg.boot, alpha = .05)$var.sel.Z)
}
}
\keyword{subsampling}
\keyword{vimp}
\keyword{confidence interval}
