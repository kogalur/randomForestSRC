
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "distance.h"
#include "stackOutput.h"
#include "nativeUtil.h"
#include "nrutil.h"
#include "error.h"
SEXP rfsrcDistance(SEXP sexp_metricType,
                   SEXP sexp_n,
                   SEXP sexp_p,
                   SEXP sexp_x,
                   SEXP sexp_sizeIJ,
                   SEXP sexp_rowI,
                   SEXP sexp_rowJ,
                   SEXP sexp_numThreads,
                   SEXP sexp_traceFlag) {
  uint    traceFlag   = INTEGER(sexp_traceFlag)[0];
  setUserTraceFlag(traceFlag);
  setNativeGlobalEnv(&RF_nativeIndex, &RF_stackCount);
  uint    metricType  = INTEGER(sexp_metricType)[0];
  uint    n           = INTEGER(sexp_n)[0];
  uint    p           = INTEGER(sexp_p)[0];
  double *x           = REAL(sexp_x);
  uint    sizeIJ      = INTEGER(sexp_sizeIJ)[0];
  RF_numThreads       = INTEGER(sexp_numThreads)[0];
  uint    *rowI;
  uint    *rowJ;
  double **xMatrix;
  double *dist;
  uint size;
  uint i, j, k;
  char  *sexpString[3] = {
    "",              
    "",              
    "distance"       
  };
  if (metricType != RF_DISTANCE_EUCLIDEAN) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Parameter verification failed.");
    RF_nativeError("\nRF-SRC:  Distance metric is invalid:  %10d \n", metricType);
    RF_nativeExit();
  }
  if (n < 2) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Parameter verification failed.");
    RF_nativeError("\nRF-SRC:  Matrix must have more than one (1) row:  %10d \n", n);
    RF_nativeExit();
  }
#ifdef _OPENMP
  if (RF_numThreads < 0) {
    RF_numThreads = omp_get_max_threads();
  }
  else {
    RF_numThreads = (RF_numThreads < omp_get_max_threads()) ? (RF_numThreads) : (omp_get_max_threads());
  }
#endif
  ${trace.token}  setTraceFlag(traceFlag, 0);
  ${memor.token}  setTraceFlag(traceFlag, 0);
  ${stack.token}  setMaxMemoryAllocation(0);
  ${stack.token}  setMinMemoryAllocation(0);
  ${trace.token}  if (getTraceFlag(0) & SUMM_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\n\nNative code rfsrc() nominal entry:  @PROJECT_BUILD@ \n");
  ${trace.token}  }
  if (sizeIJ > 0) {
    rowI        = (uint*) INTEGER(sexp_rowI); rowI--;
    rowJ        = (uint*) INTEGER(sexp_rowJ); rowJ--;
    size = sizeIJ;
  }
  else {
    size = (n * (n-1)) >> 1;
    rowI = uivector(1, size);
    rowJ = uivector(1, size);
    k = 0;
    for (i = 1; i <= n; i++) {
      for (j = 1; j < i; j++) {
        k++;
        rowI[k] = i;
        rowJ[k] = j;
      }
    }
  }
  RF_stackCount = 1;
  initProtect(RF_stackCount);
  stackAuxiliaryInfoList(&RF_snpAuxiliaryInfoList, RF_stackCount);
  dist = (double*) stackAndProtect(RF_GROW, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, 2, size, 0, sexpString[2], NULL, 1, size);
  dist --;
  xMatrix = (double **) new_vvector(1, p, NRUTIL_DPTR);
  for (i = 1; i <= p; i++) {
    xMatrix[i] = (x + ((i-1) * n) - 1);
  }
  ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nIncoming matrix:  ");
  ${trace.token}    RF_nativePrint("\n           p");
  ${trace.token}    RF_nativePrint("          n->");
  ${trace.token}    for (i = 1; i <= p - 1; i++) {
  ${trace.token}      RF_nativePrint(" %12s", " ");
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n            ");
  ${trace.token}    for (i = 1; i <= p; i++) {
  ${trace.token}      RF_nativePrint(" %12d", i);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}    for (j = 1; j <= n; j++) {
  ${trace.token}      RF_nativePrint("%12d", j);
  ${trace.token}      for (i = 1; i <= p; i++) {
  ${trace.token}        RF_nativePrint(" %12.4f", xMatrix[i][j]);
  ${trace.token}      }
  ${trace.token}      RF_nativePrint("\n");
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\nAssumed or Incoming [rowI] [rowJ]:  ");
  ${trace.token}    RF_nativePrint("\n       index         rowI         rowJ");
  ${trace.token}    for (k = 1; k <= size; k++) {
  ${trace.token}      RF_nativePrint("\n  %12d %12d %12d", k, rowI[k], rowJ[k]);
  ${trace.token}    }
  ${trace.token}  }
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
  for (k = 1; k <= size; k++) {
    dist[k] = euclidean(n, p, rowI[k], rowJ[k], xMatrix);
  }
  free_new_vvector(xMatrix, 1, p, NRUTIL_DPTR);
  if (sizeIJ > 0) {
  }
  else {
    free_uivector(rowI, 1, size);
    free_uivector(rowJ, 1, size);
  }
  unstackAuxiliaryInfoAndList(FALSE, RF_snpAuxiliaryInfoList, RF_stackCount);
  ${stack.token}  if (getTraceFlag(0) & SUMM_DEF_TRACE) {
  ${stack.token}    memoryCheck();
  ${stack.token}  }
  R_ReleaseObject(RF_sexpVector[RF_OUTP_ID]);
  R_ReleaseObject(RF_sexpVector[RF_STRG_ID]);
  ${trace.token}  if (getTraceFlag(0) & SUMM_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\n\nNative code rfsrc() nominal exit:  @PROJECT_BUILD@ \n");
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}  }
  return RF_sexpVector[RF_OUTP_ID];
}
double euclidean(uint n, uint p, uint i, uint j, double **x) {
  double result;
  double difference;
  uint   k;
  result = 0.0;
  for (k = 1; k <= p; k++) {
    difference = x[k][i] - x[k][j];
    result += (difference * difference);
  }
  result = sqrt(result);
  return result;
}
