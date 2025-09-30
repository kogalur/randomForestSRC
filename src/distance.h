#ifndef RF_DISTANCE_H
#define RF_DISTANCE_H
SEXP rfsrcDistance(SEXP sexp_metric,
                   SEXP sexp_n,
                   SEXP sexp_p,
                   SEXP sexp_x,
                   SEXP sexp_sizeIJ,
                   SEXP sexp_rowI,
                   SEXP sexp_rowJ,
                   SEXP sexp_numThreads,
                   SEXP sexp_traceFlag);
double euclidean(uint n, uint p, uint i, uint j, double **x);
#endif
