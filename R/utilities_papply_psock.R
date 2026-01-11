#' @keywords internal
.resolve_papply <- function(papply,
                            cores = getOption("mc.cores", 1L),
                            cl = NULL,
                            packages = NULL,
                            export = NULL,
                            envir = parent.frame()) {
  ## Accept either a function (lapply / mclapply) or a string ("lapply" / "mclapply")
  if (is.character(papply)) {
    papply <- match.arg(papply, c("lapply", "mclapply"))
    papply <- if (papply == "lapply") base::lapply else parallel::mclapply
  } else {
    papply <- match.fun(papply)
  }
  cores <- as.integer(cores)
  if (is.na(cores) || cores < 1L) cores <- 1L
  ## Always return something with lapply-like signature (X, FUN, ...)
  if (identical(papply, base::lapply) || cores < 2L) {
    return(function(X, FUN, ...) base::lapply(X, FUN, ...))
  }
  ## Legacy request for mclapply => use PSOCK instead
  if (identical(papply, parallel::mclapply)) {
    ## In package code, this helps because workers start as "empty" R sessions:
    ## load the calling package namespace (which also loads Imports).
    pkg <- utils::packageName()
    force(cl); force(packages); force(export); force(envir)
    return(function(X, FUN, ...,
                    mc.cores = cores,
                    mc.preschedule = TRUE,
                    mc.set.seed = TRUE,
                    mc.seed = NULL) {
      mc.cores <- as.integer(mc.cores)
      if (is.na(mc.cores) || mc.cores < 2L) {
        return(base::lapply(X, FUN, ...))
      }
      ## Create a cluster per call unless a cluster is supplied for reuse
      local_cluster <- is.null(cl)
      cl_use <- cl
      if (local_cluster) {
        cl_use <- parallel::makePSOCKcluster(mc.cores)
        on.exit(parallel::stopCluster(cl_use), add = TRUE)
      }
      ## Ensure your package (and its Imports) are available on workers
      if (!is.null(pkg) && nzchar(pkg)) {
        parallel::clusterCall(
          cl_use,
          function(p) { loadNamespace(p); NULL },
          pkg
        )
      }
      ## Optionally attach extra packages on workers (if you rely on search path)
      if (!is.null(packages) && length(packages)) {
        parallel::clusterCall(
          cl_use,
          function(pkgs) {
            for (p in pkgs) suppressPackageStartupMessages(
              require(p, character.only = TRUE)
            )
            NULL
          },
          packages
        )
      }
      ## Optional explicit exports (useful mainly for objects in globalenv())
      if (!is.null(export) && length(export)) {
        parallel::clusterExport(cl_use, export, envir = envir)
      }
      ## Optional reproducible parallel RNG
      if (isTRUE(mc.set.seed)) {
        if (is.null(mc.seed)) mc.seed <- sample.int(.Machine$integer.max, 1L)
        parallel::clusterSetRNGStream(cl_use, iseed = mc.seed)
      }
      ## Scheduling: map mclapply's mc.preschedule to parLapply vs parLapplyLB
      if (isTRUE(mc.preschedule)) {
        parallel::parLapply(cl_use, X, FUN, ...)
      } else {
        parallel::parLapplyLB(cl_use, X, FUN, ...)
      }
    })
  }
  ## Otherwise, user provided some other apply-like function; respect it
  if (is.function(papply)) return(papply)
  stop("papply must be 'lapply', 'mclapply', or a function like lapply().")
}
