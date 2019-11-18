#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> 
#include <R_ext/Rdynload.h>


// ---------------------------//
//     TO UPDATE THIS CODE    //
// ---------------------------//

// 1. Build the source package.
// 2. Launch R, and setwd() to the top-level package directory.
// 3. library(tools)
// 4. tools::package_native_routine_registration_skeleton("randomForestSRC")
// 5. Copy and Paste the output here:


// >>>>>>>>>> Changes Below >>>>>>>>>> //

/* .Call calls */
extern SEXP   rfsrcCIndex(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rfsrcDistance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP     rfsrcGrow(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP);
extern SEXP  rfsrcPredict(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"rfsrcCIndex",   (DL_FUNC) &rfsrcCIndex,    6},
    {"rfsrcDistance", (DL_FUNC) &rfsrcDistance,  9},
    {"rfsrcGrow",     (DL_FUNC) &rfsrcGrow,     42},
    {"rfsrcPredict",  (DL_FUNC) &rfsrcPredict,  59},
    {NULL, NULL, 0}
};

void R_init_randomForestSRC(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
