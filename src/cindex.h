#ifndef RF_CINDEX_H
#define RF_CINDEX_H
SEXP rfsrcCIndex(SEXP sexp_traceFlag,
                 SEXP sexp_size,
                 SEXP sexp_time,
                 SEXP sexp_censoring,
                 SEXP sexp_predicted,
                 SEXP sexp_denom);
SEXP rfsrcCIndexNew(SEXP sexp_traceFlag,
                    SEXP sexp_size,
                    SEXP sexp_time,
                    SEXP sexp_censoring,
                    SEXP sexp_predicted,
                    SEXP sexp_denom);
SEXP rfsrcTestSEXP(SEXP sexp_size);
#endif
