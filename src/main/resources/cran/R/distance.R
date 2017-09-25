distance <- function (x, 
                      method = "euclidean",
                      rowI = NULL,
                      rowJ = NULL,
                      do.trace = FALSE)
{
    n = dim(x)[1]
    p = dim(x)[2]
    if (n < 2) {
        stop("matrix must have more than one (1) row")
    }
    if (length(rowI) > 0 ) {
        if (length(rowI) != length(rowJ)) {
            stop("rowI and rowJ identifiers must have the same length")
        }
    }
    ## Legal methods.  Not all are implemented. They are placeholders.
    method.names <- c("euclidean", ## 1
                      "canberra",  ## 2
                      "maximum")   ## 3
    ## Default method is euclidean.
    if(is.null(method)) {
        method.idx <- which(method.names == "euclidean")
    }
    else {
        method.idx <- which(method.names == method)
    }
    ## Check for coherent distance method.
    if (length(method.idx) != 1) {
        stop("distance metric invalid")
    }
    ## Jump to native code.  Note that if rowI and rowJ are NULL, we return the distance matrix
    ## in it's entirity.  Only the lower-diagonal entries are returned.  If rowI and rowJ are non-NULL
    ## we return the distance for only the cell pairs (rowI[], rowJ[]).  This allows for master/slave hybrid
    ## OpenMP/MPI cluster processing.
    nativeOutput <- .Call("rfsrcDistance",
                          as.integer(method.idx),
                          as.integer(n),
                          as.integer(p),
                          as.double(x),
                          as.integer(length(rowI)),
                          as.integer(rowI),
                          as.integer(rowJ),
                          as.integer(get.rf.cores()),
                          as.integer(do.trace))
    ## check for error return condition in the native code
    if (is.null(nativeOutput)) {
        stop("An error has occurred in rfsrcDistance.  Please turn trace on for further analysis.")
    }
    if (length(rowI) > 0) {
        ## Return only the cell pairs (rowI[], rowJ[]) for processing in the master/slave scripts. 
        result <- list(rowI = rowI, rowJ = rowJ, distance = nativeOutput$distance)
    }
    else {
        ## Return the matrix in its entirety.
        result <- matrix(0, n, n)
        count <- 0
        for (k in 2:n) {
            result[k,1:(k-1)] <- nativeOutput$distance[(count + 1):(count + k - 1)]
            result[1:(k-1),k] <- result[k,1:(k-1)]
            count <- count + (k-1)
        }
    }
    return (result)
}
