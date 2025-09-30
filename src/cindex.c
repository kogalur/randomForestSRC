
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "cindex.h"
#include "stackOutput.h"
#include "nativeUtil.h"
#include "survivalE.h"
SEXP rfsrcCIndex(SEXP sexp_traceFlag,
                 SEXP sexp_size,
                 SEXP sexp_time,
                 SEXP sexp_censoring,
                 SEXP sexp_predicted,
                 SEXP sexp_denom) {
  uint    traceFlag   = INTEGER(sexp_traceFlag)[0];
  setUserTraceFlag(traceFlag);
  setNativeGlobalEnv(&RF_nativeIndex, &RF_stackCount);
  uint    size        = (uint) INTEGER(sexp_size)[0];
  double *time        = REAL(sexp_time); time--;
  double *censoring   = REAL(sexp_censoring); censoring--;
  double *predicted   = REAL(sexp_predicted); predicted--;
  double *denom       = REAL(sexp_denom); denom--;
  double *v;
  char  *sexpString[3] = {
    "",              
    "",              
    "err"            
  };
  RF_stackCount = 1;
  initProtect(RF_stackCount);
  stackAuxiliaryInfoList(&RF_snpAuxiliaryInfoList, RF_stackCount);
  v = (double*) stackAndProtect(RF_GROW, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, 2, 1, 0, sexpString[2], NULL, 1, 1);
  *v = getConcordanceIndex( 1,
                            size,
                            time,
                            censoring,
                            predicted,
                            denom);
  unstackAuxiliaryInfoAndList(FALSE, RF_snpAuxiliaryInfoList, RF_stackCount);
  R_ReleaseObject(RF_sexpVector[RF_OUTP_ID]);
  R_ReleaseObject(RF_sexpVector[RF_STRG_ID]);
  return RF_sexpVector[RF_OUTP_ID];
}
SEXP rfsrcCIndexNew(SEXP sexp_traceFlag,
                    SEXP sexp_size,
                    SEXP sexp_time,
                    SEXP sexp_censoring,
                    SEXP sexp_predicted,
                    SEXP sexp_denom) {
  uint    traceFlag   = INTEGER(sexp_traceFlag)[0];
  setUserTraceFlag(traceFlag);
  setNativeGlobalEnv(&RF_nativeIndex, &RF_stackCount);
  uint    size        = (uint) INTEGER(sexp_size)[0];
  double *time        = REAL(sexp_time); time--;
  double *censoring   = REAL(sexp_censoring); censoring--;
  double *predicted   = REAL(sexp_predicted); predicted--;
  double *denom       = REAL(sexp_denom); denom--;
  double *v;
  char  *sexpString[3] = {
    "",              
    "",              
    "err"            
  };
  stackAuxiliaryInfoList(&RF_snpAuxiliaryInfoList, RF_stackCount);
  v = (double*) stackAndProtect(RF_GROW, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, 2, 1, 0, sexpString[2], NULL, 1, 1);
  *v = getConcordanceIndexNew( 1,
                               size,
                               time,
                               censoring,
                               predicted,
                               denom);
  unstackAuxiliaryInfoAndList(FALSE, RF_snpAuxiliaryInfoList, RF_stackCount);
  R_ReleaseObject(RF_sexpVector[RF_OUTP_ID]);
  R_ReleaseObject(RF_sexpVector[RF_STRG_ID]);
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
  R_ReleaseObject(RF_sexpVector[RF_OUTP_ID]);
  R_ReleaseObject(RF_sexpVector[RF_STRG_ID]);
  return RF_sexpVector[RF_OUTP_ID];
}
