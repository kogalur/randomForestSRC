\name{rfsrc.RNGversion}
\alias{rfsrc.RNGversion}
\title{Select the RNG Version for New Forests}
\description{
  Select or query the random number generator (RNG) seed initialization
  protocol used to grow new forests with \code{\link{rfsrc}}.
  Most users do not need to call this function. New forests use the
  recommended wide generator by default. Set the RNG version only to
  reproduce or test legacy behavior.
}
\usage{
rfsrc.RNGversion(version)
}
\arguments{
  \item{version}{A character string: \code{"3.2.3"} selects the legacy
    protocol and \code{"3.10.0"} selects the wide protocol. Omit this
    argument to query the current selection. Supply \code{NULL} to remove
    the override and restore the default, \code{"3.10.0"}. Other values
    are rejected.}
}
\details{
  Use of the legacy generator is not recommended for new analyses.
  Select \code{"3.2.3"} only when reproducing or testing legacy behavior;
  otherwise use the default, \code{"3.10.0"}.

  The selection is stored in the R session option \code{rf.rng.version}.
  When that option is unset, new forests use \code{"3.10.0"}. Functions
  that grow forests through \code{rfsrc()} in the same R session inherit
  this setting. Independent R sessions, including separate worker or
  test processes, must set their own selection.

  Each newly grown forest records its selection in
  \code{fit$forest$rng.version}. Restoration, prediction, partial
  dependence, and variable importance use the saved forest's protocol,
  regardless of the current session setting. Forests without RNG
  metadata use the legacy protocol, \code{"3.2.3"}.

  These identifiers select the package's seed-head protocol. Selecting
  a protocol does not set a seed, change R's \code{\link[base]{RNGkind}},
  or alter \code{.Random.seed}. Use the \code{seed} argument to
  \code{rfsrc()} or \code{\link[base]{set.seed}} to control the seed.

  The legacy identifier marks the earliest supported RNG compatibility
  baseline, including variable importance. Selecting it does not restore
  other historical package defaults or guarantee identical results after
  unrelated algorithm changes.
}
\value{
  With no argument, the current effective RNG version as a character
  string. When setting or resetting the selection, the previous stored
  option is returned invisibly. This is \code{NULL} when no override
  was set, allowing the previous setting to be restored with another
  call to \code{rfsrc.RNGversion()}.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\seealso{
  \code{\link{rfsrc}},
  \code{\link[base]{RNGversion}},
  \code{\link[base]{set.seed}}
}
\examples{
## Only for reproducing or testing legacy behavior.
local({
  previous <- rfsrc.RNGversion("3.2.3")
  on.exit(rfsrc.RNGversion(previous))

  fit <- rfsrc(Sepal.Length ~ ., iris, ntree = 25, seed = -1234)
  print(fit$forest$rng.version)

  ## Return to the recommended wide generator for subsequent fits.
  rfsrc.RNGversion("3.10.0")
  print(rfsrc.RNGversion())

  ## The existing forest keeps its original protocol.
  print(fit$forest$rng.version)

  ## Remove the session override.
  rfsrc.RNGversion(NULL)
  print(rfsrc.RNGversion())
})
}
\keyword{utilities}
