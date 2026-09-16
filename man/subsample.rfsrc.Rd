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
  Grow forests on repeated subsamples of the training data to estimate
  uncertainty in variable importance (VIMP) and out-of-bag prediction
  error. The saved subsample estimates can be extracted, printed, and
  plotted at different confidence levels without repeating the forest
  fits. An optional double-bootstrap procedure is also available.
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
    Grow the forest with \code{forest = TRUE}. Anonymous forests and
    nonuniform case weights are not supported.
    For \code{extract.subsample}, an object returned by \code{subsample},
    including a double-bootstrap object. For \code{extract.bootsample},
    an object returned by \code{subsample(..., bootstrap = TRUE)}.}
  \item{B}{Number of outer subsamples, or double-bootstrap replicates
    when \code{bootstrap = TRUE}. Supply a positive integer. Values
    greater than one are needed for sample-standard-deviation estimates.
    Each replicate requires another forest fit.}
  \item{block.size}{Number of trees per block for calculating VIMP
    when it must first be computed for \code{obj}. If \code{obj} already
    contains VIMP, its block size is reused instead.}
  \item{importance}{VIMP perturbation method, such as \code{"anti"},
    \code{"permute"}, or \code{"random"}. With ordinary subsampling,
    existing VIMP in \code{obj} determines the method and is reused.
    Otherwise this argument is passed to \code{vimp}; when omitted,
    \code{subsample} forwards \code{importance = TRUE}.
    Ignored when \code{performance.only = TRUE}.}
  \item{subratio}{Nominal fraction of training observations used in
    each outer subsample. Values must be greater than zero and at most
    one. Choose a fraction below one giving a sufficiently large subsample
    for inference. The default is
    \eqn{\exp(-1)} for \eqn{n\leq1000} and
    \eqn{\exp(-1)\sqrt{1000/n}} for larger \eqn{n}.
    This is separate from the per-tree \code{sampsize} setting.
    Not used by the double-bootstrap procedure.}
  \item{stratify}{Use stratified outer subsampling for univariate
    classification, survival, and competing risks. Classification
    strata are response classes; survival strata are the event-status
    groups, including censoring. Multivariate and mixed-outcome
    forests use unstratified subsampling. Not used by the
    double-bootstrap procedure.}
  \item{performance}{Append the OOB prediction-error statistic to
    VIMP, allowing its uncertainty to be summarized by the same
    extraction and plotting functions. Its error measure is determined
    by the original forest's \code{perf.type}.}
  \item{performance.only}{Compute uncertainty for the prediction-error
    statistic only. Sets \code{performance = TRUE} and suppresses
    individual and joint VIMP.}
  \item{joint}{Append joint VIMP for the predictor group specified by
    \code{xvar.names}. Individual VIMP values are retained. Joint VIMP
    uses the default anti perturbation of \code{vimp}, independently
    of the method used for individual VIMP.}
  \item{xvar.names}{Predictor names defining the group for joint VIMP.
    The default \code{NULL} uses all predictors. This argument does
    not restrict the individual VIMP rows in \code{subsample}. To
    select rows for display, use \code{xvar.names} in
    \code{plot.subsample}.}
  \item{bootstrap}{Use double-bootstrap resampling instead of
    subsampling without replacement. This uses a separate replicate
    generator and interval extractor.}
  \item{verbose}{Print progress information during resampling.}
  \item{alpha}{Significance level for extraction or printing, strictly between
    zero
    and one. Intervals have nominal confidence level \eqn{1-\alpha};
    the default \code{.05} gives 95 percent intervals.}
  \item{target}{For extraction, a nonnegative integer selecting statistic
    column \code{1 + target} within the selected outcome. Use \code{0}
    for regression or right-censored survival. For classification,
    \code{0} selects the overall statistic and \code{1}, \code{2},
    and so on select class-specific columns in stored order.
    For competing risks, \code{0}, \code{1}, and so on select event
    columns in stored order. Thus the first event uses \code{target = 0}
    here, but \code{target = 1} in \code{plot.subsample}.}
  \item{m.target}{Name of one response to extract from a multivariate or
    mixed-outcome object. The default selects the first available
    regression response, or the first available classification response
    if no regression response is available. Not needed for univariate
    objects.}
  \item{standardize}{For a regression response, divide the extracted or
    printed statistics by its variance in the full training data. The same
    divisor applies to every replicate and any error row. Other families
    retain their original scale. Set \code{FALSE} for unscaled values.}
  \item{raw}{Include the replicate estimates, interval matrices, and
    other detailed extraction components in addition to the default
    summaries. This does not disable \code{standardize}.}
  \item{x}{For the print methods, an object returned by \code{subsample}.
    Use \code{print.subsample} for subsampling and \code{print.bootsample}
    for double-bootstrap output.}
  \item{...}{Additional arguments to the print methods; currently unused.}
}
\details{
  \subsection{What is resampled}{
    For ordinary subsampling, \code{subsample} draws a smaller set of
    observations without replacement and grows a new forest on that
    data set. It repeats this operation \code{B} times. The original
    forest supplies the full-data reference VIMP and error statistics;
    the new forests supply their subsample counterparts.

    The new grow calls reuse the stored formula, number of trees,
    \code{mtry}, node-size and depth settings, splitting rule,
    \code{nsplit}, predictor and split weights, sampling-size function,
    sampling type, and performance settings collected from \code{obj}.
    For survival families, the stored time grid and the choice of Uno
    weighting are also used. The procedure requires original
    root-level sampling, \code{bootstrap = "by.root"}.

    Outer subsampling and within-forest tree sampling have different
    roles. \code{subratio} determines the outer data set. The stored
    \code{sampsize} function is then applied to the size of that data
    set when growing its trees.

    Subsample refits use the grow function's default missing-value
    omission. They do not repeat the original forest's imputation,
    avoiding the cost of additional imputation within each replicate.
    When stored training tables contain missing values, the effective
    refit sample can therefore be smaller than the drawn subsample.
    The settings listed above define the refit; specialized grow
    options outside that list are not recovered.
  }
  \subsection{Choosing and reusing VIMP}{
    Requesting VIMP in the original grow call avoids a separate
    calculation on the saved trees. For example, specify
    \code{importance = "permute"} there to obtain uncertainty for
    permutation VIMP. Existing VIMP takes precedence over a different
    \code{importance} or \code{block.size} supplied to ordinary
    \code{subsample}.

    \code{joint = TRUE} adds a row named \code{joint}; it measures the
    effect of perturbing the selected predictors together.
    \code{performance = TRUE} adds a row named \code{err};
    \code{performance.only = TRUE} retains only that error row.
    These rows use the same interval machinery as the individual VIMP
    rows. A positive lower bound for an error rate is not a
    variable-selection result.
  }
  \subsection{Subsample size and stratification}{
    The default fraction is constant through \eqn{n=1000}; it decreases
    in proportion to \eqn{n^{-1/2}} only above that threshold.
    With stratification, integer sample sizes are allocated separately
    within each stratum. Small strata can receive no observations;
    stratification does not guarantee that every class or event type
    appears in each new forest. Missing stratum labels are excluded
    from the allocation.

    For survival and competing risks, each drawn subsample must retain
    at least one observed event after the usual missing-value omission.
    If it does not, the procedure stops with a message recommending
    a larger \code{subratio}. No class-presence check is added for
    classification. The allocation rule is unchanged; the procedure
    does not force one observation into every stratum or redraw a
    rejected subsample.

    The extractor uses the nominal size
    \eqn{m=nr}, where \eqn{r} is \code{subratio}, for scaling. The actual integer sample
    count, which can be smaller after stratification or missing-data
    processing, is not recorded in the returned object.
  }
  \subsection{Extraction, confidence levels, and standardization}{
    For an ordinary subsampling object, \code{extract.subsample} selects
    one outcome and returns all three interval summaries: nonparametric
    subsampling, normal intervals using the subsampling standard error,
    and normal intervals using
    the jackknife standard error. For multivariate or mixed outcomes,
    \code{m.target} selects an outcome in the stored results. It is not
    a new forest-prediction request.

    The extractor and print method default to \code{alpha = .05};
    \code{plot.subsample} defaults to \code{alpha = .01}. Set the same
    \code{alpha} explicitly to compare printed and plotted intervals.
    Each interval has nominal confidence level \eqn{1-\alpha}; the
    calculations are made separately for each statistic, without a
    multiple-testing adjustment.

    Extraction, printing, and plotting use \code{standardize = TRUE}
    by default. For a regression outcome, both full-data and subsample
    statistics are divided by the variance of that outcome in the
    full training data. The same variance is used for every replicate,
    and error rows are scaled in the same way as VIMP rows.
    Classification and survival statistics retain their original scale.
    Use \code{standardize = FALSE} for unscaled values. Regression
    standardization is defined for a finite, positive response variance.

    \code{raw = TRUE} in \code{extract.subsample} requests additional
    output, including the replicate estimates and interval matrices.
    It does not turn off standardization.
  }
  \subsection{The three subsampling intervals}{
    The same subsample fits provide two variance estimators: the
    \emph{subsampling estimator}, which centers the subsample estimates
    at their own mean, and the \emph{delete-\eqn{d} jackknife estimator},
    which centers them at the full-data estimate. A subsample retaining
    \eqn{m} observations from \eqn{n} can also be viewed as deleting
    \eqn{d=n-m} observations. The three interval summaries are normal
    intervals using either variance estimator and nonparametric
    intervals using the centered subsampling distribution.

    For one statistic, let \eqn{\widehat\theta} be the full-data
    estimate and \eqn{\theta_b^*} its estimate in subsample \eqn{b}.
    Write \eqn{r=m/n}, let \eqn{\bar\theta_*} be the mean of the
    available subsample estimates, and let \eqn{s_*} be their sample
    standard deviation.

    The subsampling standard error, returned as \code{se.Z}, is
    \deqn{\mathrm{SE}_{\mathrm{sub}}=\sqrt{r}\,s_*.}
    It measures dispersion about \eqn{\bar\theta_*}. The jackknife
    standard error, returned as \code{se.jk.Z}, uses squared deviations
    from \eqn{\widehat\theta}:
    \deqn{\mathrm{SE}_{\mathrm{jk}}=
      \left\{r\,\mathrm{mean}_b
      (\theta_b^*-\widehat\theta)^2\right\}^{1/2}.}
    This full-data centering retains both the subsample dispersion and
    the squared difference between the subsample mean and the full-data
    estimate:
    \deqn{\mathrm{mean}_b(\theta_b^*-\widehat\theta)^2
      =\mathrm{mean}_b(\theta_b^*-\bar\theta_*)^2
      +(\bar\theta_*-\widehat\theta)^2.}
    The squared mean displacement supplies the nonnegative bias
    adjustment discussed in Ishwaran and Lu (2019), equations (10)
    and (11).

    Both normal intervals use
    \eqn{\widehat\theta\pm z_{1-\alpha/2}\mathrm{SE}}.
    The subsampling version is returned in \code{var.sel.Z}; the
    jackknife version is returned in \code{var.jk.sel.Z}.

    The formulas above describe the implementation's scaling:
    \code{se.jk.Z} uses \eqn{m/n}, whereas the delete-\eqn{d} jackknife
    variance in equation (11) uses \eqn{m/(n-m)}. These scaling factors
    are asymptotically equivalent when \eqn{m/n\to0}. Also,
    for \eqn{K} available estimates, \code{se.Z} uses the
    sample-variance denominator \eqn{K-1}, while \code{se.jk.Z}
    averages squared deviations with denominator \eqn{K}. With these
    finite-replicate conventions, the returned jackknife standard
    error need not exceed the returned subsampling standard error.

    For the nonparametric interval, define the roots
    \eqn{Z_b^*=\sqrt{m}(\theta_b^*-\widehat\theta)} and let
    \eqn{q_p} be their empirical quantile, calculated using the default
    \code{quantile} interpolation. Its endpoints are
    \deqn{\left[\widehat\theta-q_{1-\alpha/2}/\sqrt{n},
      \widehat\theta-q_{\alpha/2}/\sqrt{n}\right].}
    Thus the nonparametric interval uses reversed, centered-root
    quantiles rather than the unscaled percentiles of the replicate
    VIMP values. These formulas describe the subsampling path;
    double-bootstrap objects use their own extractor.
  }
  \subsection{Double-bootstrap output}{
    With \code{bootstrap = TRUE}, each outer sample draws \eqn{n}
    observations with replacement. The per-tree samples are then generated
    from that bootstrap multiset. Under \code{samptype = "swor"},
    sampling is without replacement from multiset positions, so repeated
    original observation identities remain possible.

    \code{extract.subsample} dispatches double-bootstrap objects to
    \code{extract.bootsample}. This returns percentile intervals and
    normal intervals centered at the bootstrap mean, with the ordinary
    bootstrap standard deviation. The nonparametric bootstrap table has
    no p-value column. The centered-root p-value above is specific to
    subsampling without replacement.
  }
  \subsection{Reading extracted summaries}{
    The default extraction returns \code{se.Z}, \code{se.jk.Z},
    \code{var.sel}, \code{var.sel.Z}, and \code{var.jk.sel.Z}.
    The three tables contain interval endpoints and a \code{signif}
    column indicating that the lower endpoint exceeds zero.
    \code{var.sel} reports an interval median; the normal tables report
    the full-data estimate in their \code{mean} column.

    The normal tables report one-sided upper-tail p-values for zero
    VIMP. Their positive-lower-bound rule therefore corresponds to
    a one-sided level of \eqn{\alpha/2}, not \eqn{\alpha}.
    The nonparametric \code{pvalue} uses the same centered roots as
    its interval. It is the proportion of available roots satisfying
    \eqn{Z_b^*>\sqrt{n}\,\widehat\theta}, a one-sided upper-tail
    comparison for zero importance. Ties are excluded by the strict
    inequality. No Monte Carlo adjustment is added. The interval still
    uses the default \code{quantile} interpolation, so empirical
    p-values and interpolated interval endpoints need not yield identical
    decisions at every finite-sample boundary.

    With \code{raw = TRUE}, the extractor also returns \code{vmp},
    \code{vmpS}, \code{ci}, \code{ci.Z}, \code{ci.jk.Z}, and
    \code{boxplot.dta}. The replicate matrix \code{vmpS} has one row
    per statistic and one column per replicate. The five rows of each
    interval matrix contain its lower endpoint, lower quartile,
    center, upper quartile, and upper endpoint. \code{boxplot.dta}
    contains separately centered and rescaled replicate values.
  }
  \subsection{Printing interval summaries}{
    \code{print(x, alpha = .05, standardize = TRUE)} displays interval
    matrices for every stored outcome and statistic column, without
    fitting more forests. Subsampling output includes the nonparametric,
    normal subsampling, and normal jackknife intervals. Double-bootstrap
    output includes its nonparametric and normal intervals. The displayed
    values are rounded to three decimal places.

    \code{print.subsample.rfsrc} and \code{print.bootsample.rfsrc} are
    equivalent to \code{print.subsample} and \code{print.bootsample},
    respectively, with the same arguments. Use the extraction helpers
    to select one outcome or to obtain numerical results at full precision.
  }
}
\value{
  \code{subsample} returns a list containing the full-data reference and
  the replicate estimates.
  Interval summaries are obtained by extraction, printing, or plotting.
  \item{rf}{Reference forest object and full-data statistics. If VIMP
    was initially absent, this includes the result of calculating VIMP
    on the saved forest.}
  \item{vmp}{A list of full-data statistic matrices, indexed by outcome.
    Rows are individual predictors and any requested \code{joint} or
    \code{err} statistics. A requested but unavailable appended statistic
    is represented by \code{NA} when its reference layout is available. Classification and competing-risk outcomes
    can have multiple statistic columns. These stored values are not
    standardized.}
  \item{vmpS}{For subsampling, a length-\code{B} list of outcome-indexed
    subsample statistic tables aligned to the outcome, row, and column
    names of \code{vmp}. Extra replicate rows are discarded; unavailable
    entries retain \code{NA}.}
  \item{vmpB}{For double bootstrapping, the replicate statistics returned
    by the double-bootstrap generator, in place of \code{vmpS}.}
  \item{subratio}{Nominal outer subsample fraction. \code{NULL} for
    double-bootstrap objects.}
  \item{performance.only}{Whether only the performance statistic was
    requested.}
  The returned class includes \code{subsample} for subsampling or
  \code{bootsample} for double bootstrapping, together with the
  retained forest classes.

  \subsection{Extraction helpers}{
    Both helpers return a list for one outcome and statistic column.
    Summary tables have one row per predictor or requested joint/error
    statistic. They reuse the saved replicates without fitting forests.
    \describe{
      \item{\code{extract.subsample}}{For ordinary subsampling, returns
        \code{se.Z}, \code{se.jk.Z}, \code{var.sel}, \code{var.sel.Z},
        and \code{var.jk.sel.Z}. With \code{raw = TRUE}, also returns
        \code{vmp}, \code{vmpS}, \code{ci}, \code{ci.Z},
        \code{ci.jk.Z}, and \code{boxplot.dta}, as described in Details.
        For a double-bootstrap object, returns the
        \code{extract.bootsample} result below.}
      \item{\code{extract.bootsample}}{Returns \code{se}, the bootstrap
        standard errors; \code{var.sel}, the percentile-interval table;
        and \code{var.sel.Z}, the normal-interval table. Both tables
        contain \code{lower}, \code{mean}, \code{upper}, and \code{signif};
        only \code{var.sel.Z} includes a \code{pvalue} column.
        With \code{raw = TRUE}, also returns the bootstrap means
        \code{vmp}, replicate matrix \code{vmpS}, and five-row interval
        matrices \code{ci} and \code{ci.Z}. Both interval matrices use
        the bootstrap mean in their middle row.}
    }
    The \code{vmpS} matrix has statistics in rows and replicates in
    columns for either procedure. All extracted estimates use the
    requested \code{standardize} setting.
  }

  The print methods are used for their display side effect. Use
  \code{extract.subsample} or \code{extract.bootsample} for the numerical
  interval summaries.
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
