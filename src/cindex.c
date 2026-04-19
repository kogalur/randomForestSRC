
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "cindex.h"
#include "stackOutput.h"
#include "nativeUtil.h"
#include "survivalE.h"
#include "error.h"
${trace.token} #include "error.h"
SEXP rfsrcCIndex(SEXP sexp_traceFlag,
                 SEXP sexp_fastFlag,
                 SEXP sexp_size,
                 SEXP sexp_time,
                 SEXP sexp_censoring,
                 SEXP sexp_predicted,
                 SEXP sexp_denom,
                 SEXP sexp_weight) {
  uint    traceFlag   = INTEGER(sexp_traceFlag)[0];
  setUserTraceFlag(traceFlag);
  setNativeGlobalEnv(&RF_nativeIndex, &RF_stackCount);
  char    fastFlag    = (char) INTEGER(sexp_fastFlag)[0];
  uint    size        = (uint) INTEGER(sexp_size)[0];
  double *time        = REAL(sexp_time); time--;
  double *censoring   = REAL(sexp_censoring); censoring--;
  double *predicted   = REAL(sexp_predicted); predicted--;
  double *denom       = REAL(sexp_denom); denom--;
  double *weight;
  double *v;
  char  *sexpString[3] = {
    "",              
    "",              
    "err"            
  };
  ${trace.token}  setTraceFlag(traceFlag, 0);
  ${memor.token}  setTraceFlag(traceFlag, 0);
  ${stack.token}  setMaxMemoryAllocation(0);
  ${stack.token}  setMinMemoryAllocation(0);
  ${trace.token}  if (getTraceFlag(0) & SUMM_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\n\nNative code rfsrc() nominal entry:  @PROJECT_BUILD@ \n");
  ${trace.token}  }
  if (sexp_weight != R_NilValue) {
    weight = REAL(sexp_weight); weight--;
  }
  else {
    weight = NULL;
  }
  RF_stackCount = 1;
  initProtect(RF_stackCount);
  stackAuxiliaryInfoList(&RF_snpAuxiliaryInfoList, RF_stackCount);
  v = (double*) stackAndProtect(RF_GROW,
                                &RF_nativeIndex,
                                NATIVE_TYPE_NUMERIC,
                                2, 
                                1, 
                                0, 
                                sexpString[2],
                                NULL, 
                                1,    
                                1);   
  *v = getConcordanceIndex( fastFlag,
                            size,
                            time,
                            censoring,
                            predicted,
                            denom,
                            weight);
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
SEXP rfsrcCIndexFenwick(SEXP sexp_traceFlag,
                        SEXP sexp_eventType,
                        SEXP sexp_size,
                        SEXP sexp_time,
                        SEXP sexp_status,
                        SEXP sexp_predicted,
                        SEXP sexp_denom,
                        SEXP sexp_weight) {
  uint traceFlag = INTEGER(sexp_traceFlag)[0];
  setUserTraceFlag(traceFlag);
  setNativeGlobalEnv(&RF_nativeIndex, &RF_stackCount);
  uint eventType = (uint) INTEGER(sexp_eventType)[0];
  uint size      = (uint) INTEGER(sexp_size)[0];
  double *time      = REAL(sexp_time);      time--;
  double *status    = REAL(sexp_status);    status--;
  double *predicted = REAL(sexp_predicted); predicted--;
  double *denom     = REAL(sexp_denom);     denom--;
  double *weight    = REAL(sexp_weight);    weight--;
  double *v;
  char  *sexpString[3] = {
    "",              
    "",              
    "err"            
  };
  ${trace.token}  setTraceFlag(traceFlag, 0);
  ${memor.token}  setTraceFlag(traceFlag, 0);
  ${stack.token}  setMaxMemoryAllocation(0);
  ${stack.token}  setMinMemoryAllocation(0);
  ${trace.token}  if (getTraceFlag(0) & SUMM_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\n\nNative code rfsrc() nominal entry:  @PROJECT_BUILD@ \n");
  ${trace.token}  }
  if (eventType == 0) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Parameter verification failed.");
    RF_nativeError("rfsrcCIndexFenwick: eventType must be >= 1");
    RF_nativeExit();
  }
  if (sexp_weight == R_NilValue) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Parameter verification failed.");
    RF_nativeError("rfsrcCIndexFenwick: weight must be supplied");
  }
  if (getUserTraceFlag()) {
    RF_nativePrint("rfsrcCIndexFenwick: eventType=%u, n=%u (IPCW+Fenwick)\n", eventType, size);
  }
  RF_stackCount = 1;
  initProtect(RF_stackCount);
  stackAuxiliaryInfoList(&RF_snpAuxiliaryInfoList, RF_stackCount);
  v = (double*) stackAndProtect(RF_GROW,
                                &RF_nativeIndex,
                                NATIVE_TYPE_NUMERIC,
                                2, 
                                1, 
                                0, 
                                sexpString[2],
                                NULL, 
                                1,    
                                1);   
  if (size == 0) {
    *v = RF_nativeNaN;
  }
  else {
    *v = getCRConcordanceIndexIPCW_Fenwick(size,
                                           time,
                                           status,
                                           predicted,
                                           denom,
                                           weight,
                                           eventType);
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
SEXP rfsrcTestSEXP(SEXP sexp_size) {
  setNativeGlobalEnv(&RF_nativeIndex, &RF_stackCount);
  ulong size = (ulong) REAL(sexp_size)[0];
  char *v;
  char  *sexpString[3] = {
    "",              
    "",              
    "dummy"          
  };
  RF_stackCount = 1;
  initProtect(RF_stackCount);
  stackAuxiliaryInfoList(&RF_snpAuxiliaryInfoList, RF_stackCount);
  v = (char*) stackAndProtect(RF_GROW, &RF_nativeIndex, NATIVE_TYPE_CHARACTER, 2, size, 0, sexpString[2], NULL, 1, size);
  v --;
  unstackAuxiliaryInfoAndList(FALSE, RF_snpAuxiliaryInfoList, RF_stackCount);
  ${stack.token}  if (getTraceFlag(0) & SUMM_DEF_TRACE) {
  ${stack.token}    memoryCheck();
  ${stack.token}  }
  R_ReleaseObject(RF_sexpVector[RF_OUTP_ID]);
  R_ReleaseObject(RF_sexpVector[RF_STRG_ID]);
  return RF_sexpVector[RF_OUTP_ID];
}
