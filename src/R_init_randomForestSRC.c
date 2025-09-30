
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> 
#include <R_ext/Rdynload.h>
extern SEXP   rfsrcCIndex(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP   rfsrcCIndexNew(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rfsrcDistance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP     rfsrcGrow(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP  rfsrcPredict(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP);
static const R_CallMethodDef CallEntries[] = {
    {"rfsrcCIndex",   (DL_FUNC) &rfsrcCIndex,    6},
    {"rfsrcCIndexNew",(DL_FUNC) &rfsrcCIndex,    6},
    {"rfsrcDistance", (DL_FUNC) &rfsrcDistance,  9},
    {"rfsrcGrow",     (DL_FUNC) &rfsrcGrow,     38},
    {"rfsrcPredict",  (DL_FUNC) &rfsrcPredict,  62},
    {NULL, NULL, 0}
};
void R_init_randomForestSRC(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
