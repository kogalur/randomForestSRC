rfsrc.RNGversion <- function(version) {
  previous <- getOption("rf.rng.version")
  if (missing(version)) {
    current <- if (is.null(previous)) "3.10.0" else previous
    get.rng.version.code(current)
    return(current)
  }
  if (!is.null(version)) {
    get.rng.version.code(version)
  }
  options(rf.rng.version = version)
  invisible(previous)
}
