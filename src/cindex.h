#ifndef RF_CINDEX_H
#define RF_CINDEX_H
SEXP rfsrcCIndex(SEXP sexp_traceFlag,
                 SEXP sexp_fastFlag,
                 SEXP sexp_size,
                 SEXP sexp_time,
                 SEXP sexp_censoring,
                 SEXP sexp_predicted,
                 SEXP sexp_denom,
                 SEXP sexp_weight);
SEXP rfsrcCIndexFenwick(SEXP sexp_traceFlag,
                        SEXP sexp_eventType,
                        SEXP sexp_size,
                        SEXP sexp_time,
                        SEXP sexp_status,
                        SEXP sexp_predicted,
                        SEXP sexp_denom,
                        SEXP sexp_weight);
SEXP rfsrcTestSEXP(SEXP sexp_size);
#endif
