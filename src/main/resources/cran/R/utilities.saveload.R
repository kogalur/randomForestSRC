## o must be a list of forests!
fast.save.list <- function(o, path=NULL, testing=FALSE, units="Mb") {
  if (is.null(path)) {
    path <- paste0(getwd(), "/forest")
  }
  unlink(path, recursive=TRUE)
  dir.create(path, showWarnings=FALSE, recursive=TRUE)
  lapply(1:length(o), function(j) {
    fast.save(o[[j]], paste0(path, "/forest", j), testing=testing, units=units)
  })
}
## o must be a list of forests!
fast.load.list <- function(directory, path = NULL, testing = FALSE, units="Mb") {
  if (is.null(path)) {
    path <- paste0(getwd(), "/", directory)
  }
  else {
    path <- paste0(path, "/", directory)
  }
  files <- list.files(path)
  lapply(1:length(files), function(j) {
    fast.load(files[j], path, testing=testing, units=units)
  })
}
fast.save <- function(o, path=NULL, testing=TRUE, units="Mb") {
  if (is.null(path)) {
    path <- paste0(getwd(), "/forest")
  }
  unlink(path, recursive=TRUE)
  dir.create(path, showWarnings=FALSE, recursive=TRUE)
  ## coherence check
  if (sum(inherits(o, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2)
    stop("this function only works for objects of class `(rfsrc, grow)'")
  ## extract the forest - hereafter this is what we work with
  o <- o$forest
  ## nativeArrayTNDS
  check <- !sapply(o$nativeArrayTNDS, is.null)
  if (sum(check) > 0) {
    lapply(names(o$nativeArrayTNDS)[check], function(nm) {
      fst::write_fst(data.table::data.table(o$nativeArrayTNDS[[nm]]), paste0(path, "/nativeArrayTDNS_", nm, ".rda"))
    })
    o$nativeArrayTNDS[names(o$nativeArrayTNDS)[check]] <- NULL
  }
  ## nativeArray
  fst::write_fst(data.table::data.table(o$nativeArray), paste0(path, "/nativeArray.rda"))
  o$nativeArray <- NULL
  ## x and y
  if (!is.null(o$xvar)) {
    fst::write_fst(data.table::data.table(o$xvar), paste0(path, "/xvar.rda"))
    o$xvar <- NULL
  }
  if (!is.null(o$yvar)) {
    fst::write_fst(data.frame(o$yvar), paste0(path, "/yvar.rda"))
    o$yvar <- NULL
  }
  ## misc.
  fst::write_fst(data.frame(o$case.wt), paste0(path, "/case.wt.rda"))
  fst::write_fst(data.frame(o$leafCount), paste0(path, "/leafCount.rda"))
  fst::write_fst(data.frame(o$seed), paste0(path, "/seed.rda"))
  o$case.wt <- o$yvar <- o$leafCount <- o$seed <- NULL
  ## convert sample size to a number
  o$sampfrac <- o$sampsize(1)
  o$sampsize <- NULL
  ## survival
  if (o$family == "surv" || o$family == "surv-CR") {
    check <- !sapply(o$event.info, is.null)
    lapply(names(o$event.info)[check], function(nm) {
      fst::write_fst(data.table::data.table(o$event.info[[nm]]), paste0(path, "/event.info_", nm, ".rda"))
    })
    o$event.info[names(o$event.info[check])] <- NULL
  }
  ## output size (used for testing)
  if (testing) {
    print(lsos(o, units=units, n=length(o)))
  }
  ## forest
  saveRDS(o, file=paste0(path, "/forest.rda"), compress=FALSE)
  gc(FALSE)
  #cat("finished\n")
}
fast.load <- function(directory, path = NULL, testing = FALSE, units="Mb") {
  if (is.null(path)) {
    path <- paste0(getwd(), "/", directory)
  }
  else {
    path <- paste0(path, "/", directory)
  }
  files <- list.files(path)
  ## forest
  o <- readRDS(paste0(path, "/forest.rda"))
  gc()
  ## nativeArrayTNDS
  if (any(grepl("nativeArrayTDNS", files))) {
    target <- files[grepl("nativeArrayTDNS", files)]
    target <- gsub(".rda", "", target)
    target <- gsub("nativeArrayTDNS_", "", target)
    lapply(target, function(nm) {
      o$nativeArrayTNDS[[nm]] <<- fst::read_fst(paste0(path, "/nativeArrayTDNS_", nm, ".rda"))[[1]]
    })
  }
  ## nativeArrayTNDS
  o$nativeArray <- fst::read_fst(paste0(path, "/nativeArray.rda"))
  ## x and y
  if (any(grepl("xvar.rda", files))) {
    o$xvar <- fst::read_fst(paste0(path, "/xvar.rda"))
  }
  if (any(grepl("yvar.rda", files))) {
    o$yvar <- fst::read_fst(paste0(path, "/yvar.rda"))
    if (!is.data.frame(o$yvar)) {
      o$yvar <- o$yvar[[1]]
    }
  }
  ## misc.
  o$case.wt <- fst::read_fst(paste0(path, "/case.wt.rda"))[[1]]
  o$leafCount <- fst::read_fst(paste0(path, "/leafCount.rda"))[[1]]
  o$seed <- fst::read_fst(paste0(path, "/seed.rda"))[[1]]
  ## sample size --- hard coded
  o$sampsize <- function(x){x * o$sampfrac} 
  ## survival
  if (any(grepl("event.info", files))) {
    target <- files[grepl("event.info", files)]
    target <- gsub(".rda", "", target)
    target <- gsub("event.info_", "", target)
    lapply(target, function(nm) {
      o$event.info[[nm]] <<- fst::read_fst(paste0(path, "/event.info_", nm, ".rda"))[[1]]
    })
  }
  ## output size (used for testing)
  if (testing) {
    print(lsos(o, units=units, n=length(o)))
  }
  ##return the object
  o
}
## improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5,
                         units="Gb") {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
                           format(object.size(x), units = units)})
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}
## shorthand for improved list of objects
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}
