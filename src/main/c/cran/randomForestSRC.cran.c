SEXP rfsrcCIndex(SEXP sexp_traceFlag,
                 SEXP sexp_size,
                 SEXP sexp_time,
                 SEXP sexp_censoring,
                 SEXP sexp_predicted,
                 SEXP sexp_denom) {
  uint    traceFlag   = INTEGER(sexp_traceFlag)[0];
  setUserTraceFlag(traceFlag);
  setNativeGlobalEnv();
  uint    size        = (uint) INTEGER(sexp_size)[0];
  double *time        = REAL(sexp_time); time--;
  double *censoring   = REAL(sexp_censoring); censoring--;
  double *predicted   = REAL(sexp_predicted); predicted--;
  uint   *denom       = (uint*) INTEGER(sexp_denom); denom--;
  double *v;
  char  *sexpString[3] = {
    "",              
    "",              
    "err"            
  };
  RF_stackCount = 1;
  initProtect(RF_stackCount);
  stackAuxiliaryInfoList();
  v = (double*) stackAndProtect(&RF_nativeIndex, NATIVE_TYPE_NUMERIC, 2, 1, sexpString, NULL, 1, 0);
  *v = getConcordanceIndex( 1,
                            size,
                            time,
                            censoring,
                            predicted,
                            denom);
  unstackAuxiliaryInfoAndList();
  UNPROTECT(RF_stackCount + 2);
  return RF_sexpVector[RF_OUTP_ID];
}
SEXP rfsrcTestSEXP(SEXP sexp_size) {
  setNativeGlobalEnv();
  ulong size = (ulong) REAL(sexp_size)[0];
  char  *sexpString[3] = {
    "",              
    "",              
    "dummy"          
  };
  char *v;
  RF_stackCount = 1;
  initProtect(RF_stackCount);
  stackAuxiliaryInfoList();
  v = (char*) stackAndProtect(&RF_nativeIndex, NATIVE_TYPE_CHARACTER, 2, size, sexpString, NULL, 1, 0);
  v --;
  unstackAuxiliaryInfoAndList();
  UNPROTECT(RF_stackCount + 2);
  return RF_sexpVector[RF_OUTP_ID];
}
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
  setNativeGlobalEnv();
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
  if (sizeIJ > 0) {
    rowI        = (uint*) INTEGER(sexp_rowI); rowI--;
    rowJ        = (uint*) INTEGER(sexp_rowJ); rowJ--;
  }
  else {
    sizeIJ = (n * (n-1)) >> 1;
    rowI = uivector(1, sizeIJ);
    rowJ = uivector(1, sizeIJ);
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
  stackAuxiliaryInfoList();
  dist = (double*) stackAndProtect(&RF_nativeIndex, NATIVE_TYPE_NUMERIC, 2, sizeIJ, sexpString, NULL, 1, 0);
  dist --;
  xMatrix = (double **) new_vvector(1, p, NRUTIL_DPTR);
  for (i = 1; i <= p; i++) {
    xMatrix[i] = (x + ((i-1) * n) - 1);
  }
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
  for (k = 1; k <= sizeIJ; k++) {
    dist[k] = euclidean(n, p, rowI[k], rowJ[k], xMatrix);
  }
  unstackAuxiliaryInfoAndList();
  UNPROTECT(RF_stackCount + 2);
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
SEXP rfsrcGrow(SEXP traceFlag,
               SEXP seedPtr,
               SEXP optLow,
               SEXP optHigh,
               SEXP splitRule,
               SEXP nsplit,
               SEXP mtry,
               SEXP ytry,
               SEXP nodeSize,
               SEXP nodeDepth,
               SEXP crWeightSize,
               SEXP crWeight,
               SEXP ntree,
               SEXP observationSize,
               SEXP ySize,
               SEXP rType,
               SEXP rLevels,
               SEXP rData,
               SEXP xSize,
               SEXP xType,
               SEXP xLevels,
               SEXP bootstrapSize,
               SEXP bootstrap,
               SEXP caseWeight,
               SEXP xSplitStatWt,
               SEXP yWeight,
               SEXP xWeight,
               SEXP xData,
               SEXP timeInterestSize,
               SEXP timeInterest,
               SEXP nImpute,
               SEXP numThreads) {
  setUserTraceFlag(INTEGER(traceFlag)[0]);
  setNativeGlobalEnv();
  int seedValue           = INTEGER(seedPtr)[0];
  RF_opt                  = INTEGER(optLow)[0];
  RF_optHigh              = INTEGER(optHigh)[0];
  RF_splitRule            = INTEGER(splitRule)[0];
  RF_nsplit               = INTEGER(nsplit)[0];
  RF_mtry                 = INTEGER(mtry)[0];
  RF_ytry                 = INTEGER(ytry)[0];
  RF_nodeSize             = INTEGER(nodeSize)[0];
  RF_nodeDepth            = INTEGER(nodeDepth)[0];
  RF_crWeightSize         = INTEGER(crWeightSize)[0];
  RF_crWeight             = (double *) copy1DObject(crWeight, NATIVE_TYPE_NUMERIC, RF_crWeightSize);
  RF_ntree                = INTEGER(ntree)[0];
  RF_observationSize      = INTEGER(observationSize)[0];
  RF_ySize                = INTEGER(ySize)[0];
  RF_rType                = (char *) copy1DObject(rType, NATIVE_TYPE_CHARACTER, RF_ySize);
  RF_rLevels              = INTEGER(rLevels); RF_rLevels--;
  RF_responseIn           = (double **) copy2DObject(rData, NATIVE_TYPE_NUMERIC, RF_ySize > 0, RF_ySize, RF_observationSize);
  RF_xSize                = INTEGER(xSize)[0];
  RF_xType                = (char *) copy1DObject(xType, NATIVE_TYPE_CHARACTER, RF_xSize);
  RF_xLevels              = INTEGER(xLevels); RF_xLevels--;
  RF_bootstrapSize        = INTEGER(bootstrapSize)[0];
  RF_bootstrapIn          = (uint **) copy2DObject(bootstrap, NATIVE_TYPE_INTEGER, (RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2), RF_ntree, RF_observationSize);
  RF_caseWeight           = REAL(caseWeight);  RF_caseWeight--;
  RF_xSplitStatWt          = REAL(xSplitStatWt);  RF_xSplitStatWt--;
  RF_yWeight               = REAL(yWeight);  RF_yWeight--;
  RF_xWeight              = REAL(xWeight);  RF_xWeight--;
  RF_observationIn        = (double **) copy2DObject(xData, NATIVE_TYPE_NUMERIC, TRUE, RF_xSize, RF_observationSize);
  RF_timeInterestSize     = INTEGER(timeInterestSize)[0];
  RF_timeInterest         = REAL(timeInterest);  RF_timeInterest--;
  RF_nImpute              = INTEGER(nImpute)[0];
  RF_numThreads           = INTEGER(numThreads)[0];
  processDefaultGrow();
  stackAuxiliaryInfoList();
  rfsrc(RF_GROW, seedValue);
  free_1DObject(RF_crWeight, NATIVE_TYPE_NUMERIC, RF_crWeightSize);
  free_1DObject(RF_rType, NATIVE_TYPE_CHARACTER, RF_ySize);
  free_1DObject(RF_xType, NATIVE_TYPE_CHARACTER, RF_xSize);
  free_2DObject(RF_responseIn, NATIVE_TYPE_NUMERIC, RF_ySize > 0, RF_ySize, RF_observationSize);
  free_2DObject(RF_bootstrapIn, NATIVE_TYPE_INTEGER, (RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2), RF_ntree, RF_observationSize);
  free_2DObject(RF_observationIn, NATIVE_TYPE_NUMERIC, TRUE, RF_xSize, RF_observationSize);  
  unstackAuxiliaryInfoAndList();
  if (RF_nativeIndex != RF_stackCount) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Stack imbalance in PROTECT/UNPROTECT:  %10d + 1 versus %10d  ", RF_nativeIndex, RF_stackCount);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  UNPROTECT(RF_stackCount + 2);
  return RF_sexpVector[RF_OUTP_ID];
}
SEXP rfsrcPredict(SEXP traceFlag,
                  SEXP seedPtr,
                  SEXP optLow,
                  SEXP optHigh,
                  SEXP ntree,
                  SEXP observationSize,
                  SEXP ySize,
                  SEXP rType,
                  SEXP rTarget,
                  SEXP rTargetCount,
                  SEXP rLevels,
                  SEXP rData,
                  SEXP xSize,
                  SEXP xType,
                  SEXP xLevels,
                  SEXP xData,
                  SEXP bootstrapSize,
                  SEXP bootstrap,
                  SEXP caseWeight,
                  SEXP timeInterestSize,
                  SEXP timeInterest,
                  SEXP treeID,
                  SEXP nodeID,
                  SEXP parmID,
                  SEXP contPT,
                  SEXP mwcpSZ,
                  SEXP mwcpPT,
                  SEXP tnRMBR,
                  SEXP tnAMBR,
                  SEXP tnRCNT,
                  SEXP tnACNT,
                  SEXP totalNodeCount,
                  SEXP seed,
                  SEXP numThreads,
                  SEXP ptnCount,
                  SEXP intrPredictorSize,
                  SEXP intrPredictor,
                  SEXP sobservationSize,
                  SEXP sobservationIndv,
                                  
                  SEXP partialType,
                  SEXP partialXvar,
                  SEXP partialLength,
                  SEXP partialValue,
                  SEXP partialLength2,
                  SEXP partialXvar2,
                  SEXP partialValue2,
                  SEXP fobservationSize,
                  SEXP frSize,
                  SEXP frData,
                  SEXP fxData,
                  SEXP tnSURV,
                  SEXP tnMORT,
                  SEXP tnNLSN,
                  SEXP tnCSHZ,
                  SEXP tnCIFN,
                  SEXP tnREGR,
                  SEXP tnCLAS) {
  setUserTraceFlag(INTEGER(traceFlag)[0]);
  setNativeGlobalEnv();
  int seedValue           = INTEGER(seedPtr)[0];
  RF_opt                  = INTEGER(optLow)[0];
  RF_optHigh              = INTEGER(optHigh)[0];
  RF_ntree                = INTEGER(ntree)[0];
  RF_observationSize      = INTEGER(observationSize)[0];
  RF_ySize                = INTEGER(ySize)[0];
  RF_rType                = (char *) copy1DObject(rType, NATIVE_TYPE_CHARACTER, RF_ySize);
  RF_rTarget              = (uint *) INTEGER(rTarget); RF_rTarget --;
  RF_rTargetCount         = INTEGER(rTargetCount)[0];
  RF_rLevels              = INTEGER(rLevels); RF_rLevels--;
  RF_responseIn           = (double **) copy2DObject(rData, NATIVE_TYPE_NUMERIC, RF_ySize > 0, RF_ySize, RF_observationSize);
  RF_xSize                = INTEGER(xSize)[0];
  RF_xType                = (char *) copy1DObject(xType, NATIVE_TYPE_CHARACTER, RF_xSize);
  RF_xLevels              = INTEGER(xLevels); RF_xLevels--;
  RF_observationIn        = (double **) copy2DObject(xData, NATIVE_TYPE_NUMERIC, TRUE, RF_xSize, RF_observationSize);
  RF_bootstrapSize        = INTEGER(bootstrapSize)[0];
  RF_bootstrapIn          = (uint **) copy2DObject(bootstrap, NATIVE_TYPE_INTEGER, (RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2), RF_ntree, RF_observationSize);
  RF_caseWeight           = REAL(caseWeight);  RF_caseWeight--;
  RF_timeInterestSize     = INTEGER(timeInterestSize)[0];
  RF_timeInterest         = REAL(timeInterest);  RF_timeInterest --;
  RF_treeID_              = (uint *) INTEGER(treeID);  RF_treeID_ --;
  RF_nodeID_              = (uint *) INTEGER(nodeID);  RF_nodeID_ --;
  RF_parmID_              = (uint *) INTEGER(parmID);  RF_parmID_ --;
  RF_contPT_              = REAL(contPT);  RF_contPT_ --;
  RF_mwcpSZ_              = (uint *) INTEGER(mwcpSZ);  RF_mwcpSZ_ --;
  RF_mwcpPT_              = (uint *) INTEGER(mwcpPT);  RF_mwcpPT_ --;
  RF_RMBR_ID_             = (uint *) INTEGER(tnRMBR);
  RF_AMBR_ID_             = (uint *) INTEGER(tnAMBR);
  RF_TN_RCNT_             = (uint *) INTEGER(tnRCNT);
  RF_TN_ACNT_             = (uint *) INTEGER(tnACNT);
  RF_totalNodeCount       = INTEGER(totalNodeCount)[0];
  RF_seed_                = INTEGER(seed); RF_seed_ --;
  RF_numThreads           = INTEGER(numThreads)[0];
  RF_ptnCount             = INTEGER(ptnCount)[0];
  RF_intrPredictorSize    = INTEGER(intrPredictorSize)[0];
  RF_intrPredictor        = (uint *) INTEGER(intrPredictor);  RF_intrPredictor --;
  RF_sobservationSize     = INTEGER(sobservationSize)[0];
  RF_sobservationIndv     = (uint *) INTEGER(sobservationIndv);  RF_sobservationIndv --;
   
  RF_partialType          = INTEGER(partialType)[0];
  RF_partialXvar          = INTEGER(partialXvar)[0];
  RF_partialLength        = INTEGER(partialLength)[0];
  RF_partialValue         = REAL(partialValue); RF_partialValue --;
  RF_partialLength2       = INTEGER(partialLength2)[0];
  RF_partialXvar2         = (uint *) INTEGER(partialXvar2); RF_partialXvar2 --;
  RF_partialValue2        = REAL(partialValue2); RF_partialValue2 --;
  RF_fobservationSize     = INTEGER(fobservationSize)[0];
  RF_frSize               = INTEGER(frSize)[0];
  RF_fresponseIn          = (double **) copy2DObject(frData, NATIVE_TYPE_NUMERIC, RF_frSize > 0, RF_frSize, RF_fobservationSize);
  RF_fobservationIn       = (double **) copy2DObject(fxData, NATIVE_TYPE_NUMERIC, RF_fobservationSize > 0, RF_xSize, RF_fobservationSize);
  RF_TN_SURV_ = REAL(tnSURV);
  RF_TN_MORT_ = REAL(tnMORT);
  RF_TN_NLSN_ = REAL(tnNLSN) ;
  RF_TN_CSHZ_ = REAL(tnCSHZ);
  RF_TN_CIFN_ = REAL(tnCIFN);
  RF_TN_REGR_ = REAL(tnREGR);
  RF_TN_CLAS_ = (uint *) INTEGER(tnCLAS);
  processDefaultPredict();
  stackAuxiliaryInfoList();
  rfsrc((RF_fobservationSize > 0)? RF_PRED : RF_REST, seedValue);
  free_1DObject(RF_rType, NATIVE_TYPE_CHARACTER, RF_ySize);
  free_1DObject(RF_xType, NATIVE_TYPE_CHARACTER, RF_xSize);
  free_2DObject(RF_responseIn, NATIVE_TYPE_NUMERIC, RF_ySize > 0, RF_ySize, RF_observationSize);
  free_2DObject(RF_observationIn, NATIVE_TYPE_NUMERIC, TRUE, RF_xSize, RF_observationSize);
  free_2DObject(RF_bootstrapIn, NATIVE_TYPE_INTEGER, (RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2), RF_ntree, RF_observationSize);
  free_2DObject(RF_fresponseIn, NATIVE_TYPE_NUMERIC, RF_frSize > 0, RF_frSize, RF_fobservationSize);
  free_2DObject(RF_fobservationIn, NATIVE_TYPE_NUMERIC, RF_fobservationSize > 0 , RF_xSize, RF_fobservationSize);
  unstackAuxiliaryInfoAndList();
  if (RF_nativeIndex != RF_stackCount) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Stack imbalance in PROTECT/UNPROTECT:  %10d + 1 versus %10d  ", RF_nativeIndex, RF_stackCount);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  UNPROTECT(RF_stackCount + 2);
  return RF_sexpVector[RF_OUTP_ID];
}
void exit2R() {
  error("\nRF-SRC:  The application will now exit.\n");
}
void printR(char *format, ...) {
  char *buffer;
  va_list aptr;
  buffer = (char *) malloc(sizeof(char) * 1023);
  va_start(aptr, format);
  vsprintf(buffer, format, aptr);
  va_end(aptr);
  Rprintf(buffer);
  free((char *) buffer);
}
void setNativeGlobalEnv() {
  RF_nativeIndex = RF_stackCount = 0;
}
void *copy1DObject(SEXP arr, char type, uint size) {
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
  PROTECT(RF_sexpVector[RF_OUTP_ID] = allocVector(VECSXP, stackCount));
  PROTECT(RF_sexpVector[RF_STRG_ID] = allocVector(STRSXP, stackCount));
  setAttrib(RF_sexpVector[RF_OUTP_ID], R_NamesSymbol, RF_sexpVector[RF_STRG_ID]);
}
void *stackAndProtect(uint  *sexpIndex,
                      char   sexpType,
                      uint   sexpIdentity,
                      ulong  size,
                      char **sexpString,
                      void  *auxiliaryPtr,
                      uint   auxiliaryDimSize,
                      ...) {
  void *v;
  SEXP thisVector;
  thisVector = NULL;  
  if (((*sexpIndex) >> 6) > 0) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  S.E.X.P. vector list limit exceeded:  %20d", *sexpIndex);
          RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
          RF_nativeExit();
  }
  if (sizeof(ulong) > sizeof(uint)) {
    if (size > UINT_MAX) {
      if (TRUE) {
        RF_nativePrint("\nRF-SRC:  *** WARNING *** ");
        RF_nativePrint("\nRF-SRC:  S.E.X.P. vector element length exceeds 32-bits:  %20lu", size);
        RF_nativePrint("\nRF-SRC:  S.E.X.P. ALLOC:  %s ", sexpString[sexpIdentity]);
        RF_nativePrint("\nRF-SRC:  Please Reduce Dimensionality If Possible.");
      }
    }
  }
  va_list list;
  va_start(list, auxiliaryDimSize);
  uint *auxiliaryDim = uivector(1, auxiliaryDimSize);
  for (int i = 1; i <= auxiliaryDimSize; i++) {
    auxiliaryDim[i] = va_arg(list, unsigned int);
  }
  va_end(list);
  switch(sexpType) {
  case NATIVE_TYPE_NUMERIC:
    PROTECT(thisVector = NEW_NUMERIC(size));
    break;
  case NATIVE_TYPE_INTEGER:
    PROTECT(thisVector = NEW_INTEGER(size));
    break;
  case NATIVE_TYPE_CHARACTER:
    PROTECT(thisVector = NEW_CHARACTER(size));
    break;
  default:
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  SEXP vector element type unknown:  %20d", sexpType);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
    break;
  }
  SET_VECTOR_ELT(RF_sexpVector[RF_OUTP_ID], *sexpIndex, thisVector);
  SET_STRING_ELT(RF_sexpVector[RF_STRG_ID], *sexpIndex, mkChar(sexpString[sexpIdentity]));
  (*sexpIndex) ++;
  switch(sexpType) {
  case NATIVE_TYPE_NUMERIC:
    v = (double*) NUMERIC_POINTER(thisVector);
    break;
  case NATIVE_TYPE_INTEGER:
    v = (uint*) INTEGER_POINTER(thisVector);
    break;
  case NATIVE_TYPE_CHARACTER:
    v = (char*) CHARACTER_POINTER(thisVector);
    break;
  default:
    v = NULL;
    break;
  }
  allocateAuxiliaryInfo(sexpType,
                        sexpIdentity,
                        size,
                        v,
                        auxiliaryPtr,
                        auxiliaryDimSize,
                        auxiliaryDim);
  return v;
}
void setUserTraceFlag (uint traceFlag) {
  RF_userTraceFlag = traceFlag;
}
uint getUserTraceFlag () {
  return RF_userTraceFlag;
}
