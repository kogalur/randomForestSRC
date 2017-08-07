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
    method.names <- c("euclidean", 
                      "canberra",  
                      "maximum")   
    if(is.null(method)) {
        method.idx <- which(method.names == "euclidean")
    }
    else {
        method.idx <- which(method.names == method)
    }
    if (length(method.idx) != 1) {
        stop("distance metric invalid")
    }
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
    if (is.null(nativeOutput)) {
        stop("An error has occurred in rfsrcDistance.  Please turn trace on for further analysis.")
    }
    if (length(rowI) > 0) {
        result <- list(rowI = rowI, rowJ = rowJ, distance = nativeOutput$distance)
    }
    else {
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
