
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "nativeUtil.h"
#include "stackOutput.h"
#include "nrutil.h"
#include "error.h"
void setNativeGlobalEnv(uint *nativeIndex, uint *stackCount) {
  *nativeIndex = 0;
  *stackCount  = 0;
}
void *copy1DObject(SEXP arr, char type, uint size, char actual) {
  void   *buffer;
  char   *cbuffer;
  double *dbuffer;
  uint i;
  buffer = NULL;
  if (size > 0) {
    switch (type) {
    case NATIVE_TYPE_CHARACTER:
      cbuffer = cvector(1, size);
      for (i = 1; i <= size; i++) {
        cbuffer[i] = ((char*) CHAR(STRING_ELT(AS_CHARACTER(arr), i-1)))[0];
      }
      buffer = cbuffer;
      break;
    case NATIVE_TYPE_NUMERIC:
      dbuffer = dvector(1, size);
      for (i = 1; i <= size; i++) {
        dbuffer[i] = ((double*) REAL(arr))[i-1];
      }
      buffer = dbuffer;
      break;
    }
  }
  return buffer;
}
void free_1DObject(void *arr, char type, uint size) {
  if (size > 0) {
    switch (type) {
    case NATIVE_TYPE_CHARACTER:
      free_cvector((char *) arr, 1, size);
      break;
    case NATIVE_TYPE_NUMERIC:
      free_dvector((double *) arr, 1, size);
      break;
    }
  }
}
void *copy2DObject(SEXP arr, char type, char flag, uint row, uint col) {
  void *buffer;
  double *darray;
  uint   *iarray;
  uint i;
  buffer = NULL;  
  if (flag > 0) {
    switch (type) {
    case NATIVE_TYPE_NUMERIC:
      darray = REAL(arr);
      buffer = (double **) new_vvector(1, row, NRUTIL_DPTR);
      for (i = 1; i <= row; i++) {
        ((double **) buffer)[i] = (darray + ((i-1) * col) - 1);
      }
      break;
    case NATIVE_TYPE_INTEGER:
      iarray = (uint *) INTEGER(arr);
      buffer = (uint **) new_vvector(1, row, NRUTIL_UPTR);
      for (i = 1; i <= row; i++) {
        ((uint **) buffer)[i] = (iarray + ((i-1) * col) - 1);
      }
      break;
    }
  }
  return buffer;
}
void free_2DObject(void *arr, char type, char flag, uint row, uint col) {
  if (flag > 0) {
    switch (type) {
    case NATIVE_TYPE_NUMERIC:
      free_new_vvector((double **) arr, 1, row, NRUTIL_DPTR);
      break;
    case NATIVE_TYPE_INTEGER:
      free_new_vvector((uint **) arr, 1, row, NRUTIL_UPTR);
      break;
    }
  }
}
void initProtect(uint  stackCount) {
  if (stackCount > 0) {
    PROTECT(RF_sexpVector[RF_OUTP_ID] = allocVector(VECSXP, stackCount));
    PROTECT(RF_sexpVector[RF_STRG_ID] = allocVector(STRSXP, stackCount));
    setAttrib(RF_sexpVector[RF_OUTP_ID], R_NamesSymbol, RF_sexpVector[RF_STRG_ID]);
    R_PreserveObject(RF_sexpVector[RF_OUTP_ID]);
    R_PreserveObject(RF_sexpVector[RF_STRG_ID]);
    UNPROTECT(2);
  }
}
void *stackAndProtect(char   mode,
                      uint  *sexpIndex,
                      char   sexpType,
                      uint   sexpIdentity,
                      ulong  size,
                      double value,
                      char  *sexpString,
                      void  *auxiliaryPtr,
                      uint   auxiliaryDimCnt,
                      ...) {
  void *v;
  SEXP thisVector;
  thisVector = NULL;  
  v          = NULL;  
  if (sizeof(ulong) > sizeof(uint)) {
    if (size > UINT_MAX) {
      if (TRUE) {
        RF_nativePrint("\nRF-SRC:  *** WARNING *** ");
        RF_nativePrint("\nRF-SRC:  S.E.X.P. vector element length exceeds 32-bits:  %20lu", size);
        RF_nativePrint("\nRF-SRC:  S.E.X.P. ALLOC:  %s ", sexpString);
        RF_nativePrint("\nRF-SRC:  Please Reduce Dimensionality If Possible.");
      }
    }
  }
  va_list list;
  va_start(list, auxiliaryDimCnt);
  int *auxiliaryDim = ivector(1, auxiliaryDimCnt);
  for (uint i = 1; i <= auxiliaryDimCnt; i++) {
    auxiliaryDim[i] = va_arg(list, int);
  }
  va_end(list);
  if (!(size > 0)) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  SEXP vector element is of size zero (0) and of aux dimensionality:  %20d", auxiliaryDimCnt);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  switch(sexpType) {
  case NATIVE_TYPE_NUMERIC:
    thisVector = PROTECT(allocVector(REALSXP, size));
    break;
  case NATIVE_TYPE_INTEGER:
    thisVector = PROTECT(allocVector(INTSXP, size));
    break;
  case NATIVE_TYPE_CHARACTER:
    thisVector = PROTECT(allocVector(STRSXP, size));
    break;
  default:
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  SEXP vector element type unknown:  %20d", sexpType);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
    break;
  }
  SET_VECTOR_ELT(RF_sexpVector[RF_OUTP_ID], *sexpIndex, thisVector);
  SET_STRING_ELT(RF_sexpVector[RF_STRG_ID], *sexpIndex, mkChar(sexpString));
  UNPROTECT(1);
  switch(sexpType) {
  case NATIVE_TYPE_NUMERIC:
    v = (double*) REAL(thisVector);
    for (ulong i = 0; i < size; i++) {
      ((double*) v)[i] = value;
    }
    break;
  case NATIVE_TYPE_INTEGER:
    v = (uint*) INTEGER(thisVector);
    for (ulong i = 0; i < size; i++) {
      ((uint*) v)[i] = 0;
    }
    break;
  case NATIVE_TYPE_CHARACTER:
    v = (char*) CHAR(thisVector);
    for (ulong i = 0; i < size; i++) {
      ((char*) v)[i] = 0x00;
    }
    break;
  }
  allocateAuxiliaryInfo((mode == RF_GROW) ? FALSE : TRUE,
                        sexpType,
                        sexpString,
                        RF_snpAuxiliaryInfoList,
                        *sexpIndex,
                        v,
                        auxiliaryPtr,
                        auxiliaryDimCnt,
                        auxiliaryDim);
  free_ivector(auxiliaryDim, 1, auxiliaryDimCnt);
  (*sexpIndex) ++;
  return v;
}
void setUserTraceFlag (uint traceFlag) {
  RF_userTraceFlag = traceFlag;
}
uint getUserTraceFlag (void) {
  return RF_userTraceFlag;
}
