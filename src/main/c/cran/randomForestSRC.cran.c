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
  memoryCheck();
  R_ReleaseObject(RF_sexpVector[RF_OUTP_ID]);
  R_ReleaseObject(RF_sexpVector[RF_STRG_ID]);
  return RF_sexpVector[RF_OUTP_ID];
}
SEXP rfsrcTestSEXP(SEXP sexp_size) {
  setNativeGlobalEnv();
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
  memoryCheck();
  R_ReleaseObject(RF_sexpVector[RF_OUTP_ID]);
  R_ReleaseObject(RF_sexpVector[RF_STRG_ID]);
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
  memoryCheck();
  R_ReleaseObject(RF_sexpVector[RF_OUTP_ID]);
  R_ReleaseObject(RF_sexpVector[RF_STRG_ID]);
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
               SEXP lot,
               SEXP baseLearn,
               SEXP vtry,
               SEXP vtryArray,
               SEXP vtryExperimental,
               SEXP ytry,
               SEXP nodeSize,
               SEXP nodeDepth,
               SEXP crWeightSize,
               SEXP crWeight,
               SEXP vimpThreshold,
               SEXP ntree,
               SEXP observationSize,
               SEXP yInfo,
               SEXP yLevels,
               SEXP yData,
               SEXP xInfo,
               SEXP xLevels,
               SEXP xData,
               SEXP sampleInfo,
               SEXP xWeightStat,
               SEXP yWeight,
               SEXP xWeight,
               SEXP timeInterest,
               SEXP nImpute,
               SEXP perfBlock,
               SEXP quantileInfo,
               SEXP xPreSort,
               SEXP numThreads) {
  clock_t cpuTimeStart = clock();
  setUserTraceFlag(INTEGER(traceFlag)[0]);
  setNativeGlobalEnv();
  int seedValue           = INTEGER(seedPtr)[0];
  RF_opt                  = INTEGER(optLow)[0];
  RF_optHigh              = INTEGER(optHigh)[0];
  RF_splitRule            = INTEGER(splitRule)[0];
  RF_nsplit               = INTEGER(nsplit)[0];
  RF_mtry                 = INTEGER(mtry)[0];
  RF_hdim = 0;
  RF_lotSize = 0;
  RF_lotLag = 0;
  RF_lotStrikeout = 0;
  if (VECTOR_ELT(lot, 0) != R_NilValue) {
    RF_hdim = INTEGER(VECTOR_ELT(lot, 0))[0];
  }
  if (RF_hdim > 0) {
    if (VECTOR_ELT(lot, 1) != R_NilValue) {
      RF_lotSize = INTEGER(VECTOR_ELT(lot, 1))[0];
    }
    if (VECTOR_ELT(lot, 2) != R_NilValue) {
      RF_lotLag          = INTEGER(VECTOR_ELT(lot, 2))[0];
    }
    if (VECTOR_ELT(lot, 3) != R_NilValue) {
      RF_lotStrikeout    = INTEGER(VECTOR_ELT(lot, 3))[0];
    }
  }
  RF_baseLearnDepthINTR = 0;
  RF_baseLearnRuleINTR  = AUGT_INTR_NONE;
  RF_baseLearnDepthSYTH = 0;
  RF_baseLearnDimReduce = FALSE;
  if (VECTOR_ELT(baseLearn, 0) != R_NilValue) {
    RF_baseLearnDepthINTR = INTEGER(VECTOR_ELT(baseLearn, 0))[0];
  }
  if (RF_baseLearnDepthINTR > 1) {
    if (VECTOR_ELT(baseLearn, 1) != R_NilValue) {
      RF_baseLearnRuleINTR = INTEGER(VECTOR_ELT(baseLearn, 1))[0];
    }
  }
  if (VECTOR_ELT(baseLearn, 2) != R_NilValue) {
    RF_baseLearnDepthSYTH = INTEGER(VECTOR_ELT(baseLearn, 2))[0];
  }
  if (VECTOR_ELT(baseLearn, 3) != R_NilValue) {
    RF_baseLearnDimReduce = (INTEGER(VECTOR_ELT(baseLearn, 3))[0] ? TRUE : FALSE);
  }
  RF_ytry                 = INTEGER(ytry)[0];
  RF_nodeSize             = INTEGER(nodeSize)[0];
  RF_nodeDepth            = INTEGER(nodeDepth)[0];
  RF_crWeightSize         = INTEGER(crWeightSize)[0];
  RF_crWeight             = REAL(crWeight); RF_crWeight--;
  RF_vimpThreshold        = REAL(vimpThreshold)[0];
  RF_ntree                = INTEGER(ntree)[0];
  RF_observationSize      = INTEGER(observationSize)[0];
  RF_ySize                = INTEGER(VECTOR_ELT(yInfo, 0))[0];
  if(VECTOR_ELT(yInfo, 1) != R_NilValue) {
    RF_rType = (char *) copy1DObject(VECTOR_ELT(yInfo, 1), NATIVE_TYPE_CHARACTER, RF_ySize, TRUE);
  }
  else {
    RF_rType = NULL;
  }
  if(VECTOR_ELT(yInfo, 2) != R_NilValue) {
    RF_rLevelsMax         = (uint *) INTEGER(VECTOR_ELT(yInfo, 2)); RF_rLevelsMax--;
  }
  else {
    RF_rLevelsMax         = NULL;
  }
  if(VECTOR_ELT(yInfo, 3) != R_NilValue) {
    RF_rLevelsCnt         = (uint *) INTEGER(VECTOR_ELT(yInfo, 3)); RF_rLevelsCnt--;
  }
  else {
    RF_rLevelsCnt         = NULL;
  }
  if(VECTOR_ELT(yInfo, 4) != R_NilValue) {
    RF_subjIn             =  (uint *) INTEGER(VECTOR_ELT(yInfo, 4));  RF_subjIn --;
  }
  else {
    RF_subjIn             = NULL;
  }
  if(VECTOR_ELT(yInfo, 5) != R_NilValue) {
    RF_eventTypeSize      =  INTEGER(VECTOR_ELT(yInfo, 5))[0];
  }
  else {
    RF_eventTypeSize      = 0;
  }
  if(VECTOR_ELT(yInfo, 6) != R_NilValue) {
    if (RF_eventTypeSize > 0) {
      RF_eventType        =  (uint *) INTEGER(VECTOR_ELT(yInfo, 6));  RF_eventType --;
    }
    else {
      RF_eventType        = NULL;
    }
  }
  else {
    RF_eventType          = NULL;
  }
  RF_rLevelsSEXP = yLevels;
  RF_rLevels = NULL;
  if(RF_ySize > 0) {
    RF_responseIn         = (double **) copy2DObject(yData, NATIVE_TYPE_NUMERIC, TRUE, RF_ySize, RF_observationSize);
  }
  else {
    RF_responseIn = NULL;
  }
  RF_xSize                 = INTEGER(VECTOR_ELT(xInfo, 0))[0];
  if(VECTOR_ELT(xInfo, 1) != R_NilValue) {  
    RF_xType                = (char *) copy1DObject(VECTOR_ELT(xInfo, 1), NATIVE_TYPE_CHARACTER, RF_xSize, TRUE);
  }
  else {
    RF_xType = NULL;
  }
  if(VECTOR_ELT(xInfo, 2) != R_NilValue) {  
    RF_xLevelsMax           = (uint *) INTEGER(VECTOR_ELT(xInfo, 2)); RF_xLevelsMax --;
  }
  else {
    RF_xLevelsMax = NULL;
  }
  if(VECTOR_ELT(xInfo, 3) != R_NilValue) {
    RF_xLevelsCnt         = (uint *) INTEGER(VECTOR_ELT(xInfo, 3)); RF_xLevelsCnt --;
  }
  else {
    RF_xLevelsCnt         = NULL;
  }
  if(VECTOR_ELT(xInfo, 4) != R_NilValue) {
    RF_xtType             = (uint *) INTEGER(VECTOR_ELT(xInfo, 4)); RF_xtType --;
  }
  else {
    RF_xtType             = NULL;
  }
  if(VECTOR_ELT(xInfo, 5) != R_NilValue) {
    RF_stType             = (uint *) INTEGER(VECTOR_ELT(xInfo, 5)); RF_stType --;
  }
  else {
    RF_stType             = NULL;
  }
  RF_xLevelsSEXP = xLevels;
  RF_xLevels = NULL;
  if (RF_xSize > 0) {
    RF_observationIn      = (double **) copy2DObject(xData, NATIVE_TYPE_NUMERIC, TRUE, RF_xSize, RF_observationSize);
  }
  else {
    RF_observationIn = NULL;
  }
  RF_subjSize             = RF_observationSize;
  RF_subjWeight           = NULL;
  RF_bootstrapSize        = RF_observationSize;
  RF_bootstrapIn          = NULL;
  if(VECTOR_ELT(sampleInfo, 0) != R_NilValue) {
    RF_subjSize = INTEGER(VECTOR_ELT(sampleInfo, 0))[0];
    if(VECTOR_ELT(sampleInfo, 1) != R_NilValue) {
      RF_subjWeight = REAL(VECTOR_ELT(sampleInfo, 1)); RF_subjWeight--;
    }
    if(VECTOR_ELT(sampleInfo, 2) != R_NilValue) {  
      RF_bootstrapSize        = INTEGER(VECTOR_ELT(sampleInfo, 2))[0];
      if(VECTOR_ELT(sampleInfo, 3) != R_NilValue) {
        RF_bootstrapIn = (uint **) copy2DObject(VECTOR_ELT(sampleInfo, 3), NATIVE_TYPE_INTEGER, (RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2), RF_ntree, RF_subjSize);
      }
    }
    else {
      RF_bootstrapSize        = RF_subjSize;
    }
  }
  RF_xWeightStat          = REAL(xWeightStat);  RF_xWeightStat--;
  RF_yWeight              = REAL(yWeight);  RF_yWeight--;
  RF_xWeight              = REAL(xWeight);  RF_xWeight--;
  RF_timeInterestSize = INTEGER(VECTOR_ELT(timeInterest, 0))[0];
  if (VECTOR_ELT(timeInterest, 1) != R_NilValue) {
    RF_timeInterest         = (double *) REAL(VECTOR_ELT(timeInterest, 1));
    RF_timeInterest --;
  }
  else {
    RF_timeInterest = NULL;
  }
  RF_nImpute              = INTEGER(nImpute)[0];
  RF_perfBlock            = INTEGER(perfBlock)[0];
  RF_quantileSize = INTEGER(VECTOR_ELT(quantileInfo, 0))[0];
  if (VECTOR_ELT(quantileInfo, 1) != R_NilValue) {
    RF_quantile = (double *) REAL(VECTOR_ELT(quantileInfo, 1));
    RF_quantile --;
  }
  else {
    RF_quantile = NULL;
  }
  RF_qEpsilon = REAL(VECTOR_ELT(quantileInfo, 2))[0];
  RF_xPreSort            = REAL(xPreSort)[0];
  RF_vtry                = INTEGER(vtry)[0];
  RF_vtryArray           = (uint **) copy2DObject(vtryArray, NATIVE_TYPE_INTEGER, RF_vtry > 0, RF_ntree, RF_xSize);
  RF_vtryMode            = RF_VTRY_NULL;
  RF_vtryBlockSize       = 0;
  if (vtryExperimental != R_NilValue) {
    RF_vtry = INTEGER(VECTOR_ELT(vtryExperimental, 0))[0];
    if (RF_vtry > 0) {
      if (VECTOR_ELT(vtryExperimental, 1) != R_NilValue) {
        RF_vtryArray = (uint **) copy2DObject(VECTOR_ELT(vtryExperimental, 1), NATIVE_TYPE_INTEGER, RF_vtry > 0, RF_ntree, RF_xSize);
        RF_vtryBlockSize = INTEGER(VECTOR_ELT(vtryExperimental, 2))[0]; 
        RF_vtryMode  = INTEGER(VECTOR_ELT(vtryExperimental, 3))[0];
      }
      else {
        RF_vtry = 0;
      }
    }
  }
  RF_numThreads           = INTEGER(numThreads)[0];
  processDefaultGrow();
  rfsrc(RF_GROW, seedValue);
  free_1DObject(RF_rType, NATIVE_TYPE_CHARACTER, RF_ySize);
  free_1DObject(RF_xType, NATIVE_TYPE_CHARACTER, RF_xSize);
  free_2DObject(RF_responseIn, NATIVE_TYPE_NUMERIC, RF_ySize > 0, RF_ySize, RF_observationSize);
  free_2DObject(RF_bootstrapIn, NATIVE_TYPE_INTEGER, (RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2), RF_ntree, RF_subjSize);
  free_2DObject(RF_observationIn, NATIVE_TYPE_NUMERIC, TRUE, RF_xSize, RF_observationSize);  
  free_2DObject(RF_vtryArray, NATIVE_TYPE_INTEGER, RF_vtry > 0, RF_ntree, RF_xSize);  
  memoryCheck();
  RF_cpuTime_[1] = (double) (clock() - cpuTimeStart) / CLOCKS_PER_SEC;
  R_ReleaseObject(RF_sexpVector[RF_OUTP_ID]);
  R_ReleaseObject(RF_sexpVector[RF_STRG_ID]);  
  return RF_sexpVector[RF_OUTP_ID];
}
SEXP rfsrcPredict(SEXP traceFlag,
                  SEXP seedPtr,
                  SEXP optLow,
                  SEXP optHigh,
                  SEXP vimpThreshold,
                  SEXP ntree,
                  SEXP observationSize,
                  SEXP yInfo,
                  SEXP yLevels,
                  SEXP yData,
                  SEXP xInfo,
                  SEXP xLevels,
                  SEXP xData,
                  SEXP sampleInfo,
                  SEXP timeInterest,
                  SEXP totalNodeCount,
                  SEXP tLeafCount,
                  SEXP seedInfo,
                  SEXP hdim,
                  SEXP baseLearn,
                  SEXP treeID,
                  SEXP nodeID,
                  SEXP brnodeID,
                  SEXP hc_zero,
                  SEXP hc_oneAugIntr,
                  SEXP hc_oneAugSyth,
                  SEXP hc_one,
                  SEXP hc_parmID,
                  SEXP hc_contPT,
                  SEXP hc_contPTR,
                  SEXP hc_mwcpSZ,
                  SEXP hc_fsrecID,
                  SEXP hc_mwcpPT,
                  SEXP hc_augmXone,
                  SEXP hc_augmXtwo,
                  SEXP hc_augmXS,
                  SEXP hc_augmSythTop,
                  SEXP tnRMBR,
                  SEXP tnAMBR,
                  SEXP tnRCNT,
                  SEXP tnACNT,
                  SEXP tnSURV,
                  SEXP tnMORT,
                  SEXP tnNLSN,
                  SEXP tnCSHZ,
                  SEXP tnCIFN,
                  SEXP tnREGR,
                  SEXP tnCLAS,
                  SEXP yTarget,
                  SEXP ptnCount,
                  SEXP xMarginalInfo,
                  SEXP intrPredictorInfo,
                  SEXP partial,
                  SEXP fobservationSize,
                  SEXP frSize,
                  SEXP frData,
                  SEXP fxData,
                  SEXP perfBlock,
                  SEXP quantile,
                  SEXP getTree,
                  SEXP numThreads) {
  char mode;
  uint i;
clock_t cpuTimeStart = clock();
  setUserTraceFlag(INTEGER(traceFlag)[0]);
  setNativeGlobalEnv();
  int seedValue           = INTEGER(seedPtr)[0];
  RF_opt                  = INTEGER(optLow)[0];
  RF_optHigh              = INTEGER(optHigh)[0];
  RF_vimpThreshold        = REAL(vimpThreshold)[0];
  RF_ntree                = INTEGER(ntree)[0];
  RF_observationSize      = INTEGER(observationSize)[0];
  RF_ySize                = INTEGER(VECTOR_ELT(yInfo, 0))[0];
  if(VECTOR_ELT(yInfo, 1) != R_NilValue) {
    RF_rType = (char *) copy1DObject(VECTOR_ELT(yInfo, 1), NATIVE_TYPE_CHARACTER, RF_ySize, TRUE);
  }
  else {
    RF_rType = NULL;
  }
  if(VECTOR_ELT(yInfo, 2) != R_NilValue) {
    RF_rLevelsMax         = (uint *) INTEGER(VECTOR_ELT(yInfo, 2)); RF_rLevelsMax--;
  }
  else {
    RF_rLevelsMax         = NULL;
  }
  if(VECTOR_ELT(yInfo, 3) != R_NilValue) {
    RF_rLevelsCnt         = (uint *) INTEGER(VECTOR_ELT(yInfo, 3)); RF_rLevelsCnt--;
  }
  else {
    RF_rLevelsCnt = NULL;
  }
  if(VECTOR_ELT(yInfo, 4) != R_NilValue) {
    RF_subjIn             =  (uint *) INTEGER(VECTOR_ELT(yInfo, 4));  RF_subjIn --;
  }
  else {
    RF_subjIn             = NULL;
  }
  if(VECTOR_ELT(yInfo, 5) != R_NilValue) {
    RF_eventTypeSize      =  INTEGER(VECTOR_ELT(yInfo, 5))[0];
  }
  else {
    RF_eventTypeSize      = 0;
  }
  if(VECTOR_ELT(yInfo, 6) != R_NilValue) {
    if (RF_eventTypeSize > 0) {
      RF_eventType        =  (uint *) INTEGER(VECTOR_ELT(yInfo, 6));  RF_eventType --;
    }
    else {
      RF_eventType        = NULL;
    }
  }
  else {
    RF_eventType          = NULL;
  }
  RF_rLevelsSEXP = yLevels;
  RF_rLevels = NULL;
  if(RF_ySize > 0) {
    RF_responseIn           = (double **) copy2DObject(yData, NATIVE_TYPE_NUMERIC, TRUE, RF_ySize, RF_observationSize);
  }
  else {
    RF_responseIn = NULL;
  }
  RF_xSize                 = INTEGER(VECTOR_ELT(xInfo, 0))[0];
  if(VECTOR_ELT(xInfo, 1) != R_NilValue) {
    RF_xType                = (char *) copy1DObject(VECTOR_ELT(xInfo, 1), NATIVE_TYPE_CHARACTER, RF_xSize, TRUE);
  }
  else {
    RF_xType = NULL;
  }
  if(VECTOR_ELT(xInfo, 2) != R_NilValue) {  
    RF_xLevelsMax           = (uint *) INTEGER(VECTOR_ELT(xInfo, 2)); RF_xLevelsMax--;
  }
  else {
    RF_xLevelsMax = NULL;
  }
  if(VECTOR_ELT(xInfo, 3) != R_NilValue) {
    RF_xLevelsCnt           = (uint *) INTEGER(VECTOR_ELT(xInfo, 3)); RF_xLevelsCnt --;
  }
  else {
    RF_xLevelsCnt = NULL;
  }
  RF_xtType = NULL;
  RF_stType = NULL;
  RF_xLevelsSEXP = xLevels;
  RF_xLevels = NULL;
  if (RF_xSize > 0) {
    RF_observationIn      = (double **) copy2DObject(xData, NATIVE_TYPE_NUMERIC, TRUE, RF_xSize, RF_observationSize);
  }
  else {
    RF_observationIn = NULL;
  }
  RF_subjSize             = 0;
  RF_subjWeight           = NULL;
  RF_bootstrapSize        = 0;
  RF_bootstrapIn          = NULL;
  if(VECTOR_ELT(sampleInfo, 0) != R_NilValue) {
    RF_subjSize = INTEGER(VECTOR_ELT(sampleInfo, 0))[0];
    if(VECTOR_ELT(sampleInfo, 1) != R_NilValue) {
      RF_subjWeight = REAL(VECTOR_ELT(sampleInfo, 1)); RF_subjWeight--;
    }
    if(VECTOR_ELT(sampleInfo, 2) != R_NilValue) {  
      RF_bootstrapSize        = INTEGER(VECTOR_ELT(sampleInfo, 2))[0];
      if(VECTOR_ELT(sampleInfo, 3) != R_NilValue) {
        RF_bootstrapIn = (uint **) copy2DObject(VECTOR_ELT(sampleInfo, 3), NATIVE_TYPE_INTEGER, (RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2), RF_ntree, RF_subjSize);
      }
    }
    else {
      RF_bootstrapSize        = RF_subjSize;
    }
  }
  RF_timeInterestSize = INTEGER(VECTOR_ELT(timeInterest, 0))[0];
  if (VECTOR_ELT(timeInterest, 1) != R_NilValue) {
    RF_timeInterest         = (double *) REAL(VECTOR_ELT(timeInterest, 1));
    RF_timeInterest --;
  }
  else {
    RF_timeInterest = NULL;
  }
  RF_totalNodeCount       = INTEGER(totalNodeCount)[0];
  RF_tLeafCount           = (uint *) INTEGER(tLeafCount); RF_tLeafCount --;
  RF_seed_                = (int *) INTEGER(VECTOR_ELT(seedInfo, 0)); RF_seed_ --;
  if (VECTOR_ELT(seedInfo, 1) != R_NilValue) {
    RF_seedVimp_          = (int *) INTEGER(VECTOR_ELT(seedInfo, 1)); RF_seedVimp_ --;
  }
  else {
    RF_seedVimp_          = NULL;
  }
  RF_optLoGrow            = (uint) INTEGER(VECTOR_ELT(seedInfo, 2))[0]; 
  RF_hdim                 = INTEGER(hdim)[0];
  RF_baseLearnDepthINTR = 0;
  RF_baseLearnRuleINTR  = AUGT_INTR_NONE;
  RF_baseLearnDepthSYTH = 0;
  if (VECTOR_ELT(baseLearn, 0) != R_NilValue) {
    RF_baseLearnDepthINTR = INTEGER(VECTOR_ELT(baseLearn, 0))[0];
  }
  if (RF_baseLearnDepthINTR > 1) {
    if (VECTOR_ELT(baseLearn, 1) != R_NilValue) {
      RF_baseLearnRuleINTR = INTEGER(VECTOR_ELT(baseLearn, 1))[0];
    }
  }
  if (VECTOR_ELT(baseLearn, 2) != R_NilValue) {
    RF_baseLearnDepthSYTH = INTEGER(VECTOR_ELT(baseLearn, 2))[0];
  }
  RF_treeID_              = (uint *) INTEGER(treeID);   RF_treeID_ --;
  RF_nodeID_              = (uint *) INTEGER(nodeID);   RF_nodeID_ --;
  RF_brnodeID_            = (uint *) INTEGER(brnodeID); RF_brnodeID_ --;
  RF_RMBR_ID_             = (uint *) INTEGER(tnRMBR);
  RF_AMBR_ID_             = (uint *) INTEGER(tnAMBR);
  RF_TN_RCNT_             = (uint *) INTEGER(tnRCNT);
  RF_TN_ACNT_             = (uint *) INTEGER(tnACNT);
  RF_perfBlock            = INTEGER(perfBlock)[0];
  RF_quantileSize = INTEGER(VECTOR_ELT(quantile, 0))[0];
  if (VECTOR_ELT(quantile, 1) != R_NilValue) {
    RF_quantile = (double *) REAL(VECTOR_ELT(quantile, 1));
    RF_quantile --;
  }
  else {
    RF_quantile = NULL;
  }
  RF_qEpsilon = REAL(VECTOR_ELT(quantile, 2))[0];
  RF_numThreads           = INTEGER(numThreads)[0];
  RF_ptnCount             = INTEGER(ptnCount)[0];
  RF_rTargetCount         = INTEGER(VECTOR_ELT(yTarget, 0))[0];
  if (VECTOR_ELT(yTarget, 1) != R_NilValue) {
    RF_rTarget         = (uint *) INTEGER(VECTOR_ELT(yTarget, 1));
    RF_rTarget --;
  }
  else {
    RF_rTarget = NULL;
  }
  RF_intrPredictorSize = INTEGER(VECTOR_ELT(intrPredictorInfo, 0))[0];
  if (VECTOR_ELT(intrPredictorInfo, 1) != R_NilValue) {
    RF_intrPredictor = (uint *) INTEGER(VECTOR_ELT(intrPredictorInfo, 1));
    RF_intrPredictor --;
  }
  else {
    RF_intrPredictor = NULL;
  }
  if (xMarginalInfo != R_NilValue) {
    RF_xMarginalSize     = INTEGER(VECTOR_ELT(xMarginalInfo, 0))[0];
    if (VECTOR_ELT(xMarginalInfo, 1) != R_NilValue) {
      RF_xMarginal         = (uint *) INTEGER(VECTOR_ELT(xMarginalInfo, 1));
      RF_xMarginal --;
    }
    else {
      RF_xMarginal = NULL;
    }
  }
  else {
    RF_xMarginalSize = 0;
    RF_xMarginal = NULL;
  }
  RF_partialType          = INTEGER(VECTOR_ELT(partial, 0))[0];
  RF_partialXvar          = INTEGER(VECTOR_ELT(partial, 1))[0];
  RF_partialLength        = INTEGER(VECTOR_ELT(partial, 2))[0];
  if (VECTOR_ELT(partial, 3) != R_NilValue) {
    RF_partialValue         = REAL(VECTOR_ELT(partial, 3)); RF_partialValue --;
  }
  else {
    RF_partialValue = NULL;
  }
  RF_partialLength2       = INTEGER(VECTOR_ELT(partial, 4))[0];
  if (VECTOR_ELT(partial, 5) != R_NilValue) {
    RF_partialXvar2         = (uint *) INTEGER(VECTOR_ELT(partial, 5)); RF_partialXvar2 --;
  }
  else {
    RF_partialXvar2 = NULL;
  }
  if (VECTOR_ELT(partial, 6) != R_NilValue) {
    RF_partialValue2        = REAL(VECTOR_ELT(partial, 6)); RF_partialValue2 --;
  }
  else {
    RF_partialValue2 = NULL;
  }
  RF_fobservationSize     = INTEGER(fobservationSize)[0];
  RF_frSize               = INTEGER(frSize)[0];
  RF_fresponseIn          = (double **) copy2DObject(frData, NATIVE_TYPE_NUMERIC, RF_frSize > 0, RF_frSize, RF_fobservationSize);
  RF_fobservationIn       = (double **) copy2DObject(fxData, NATIVE_TYPE_NUMERIC, RF_fobservationSize > 0, RF_xSize, RF_fobservationSize);
  RF_getTree = (uint *) INTEGER(getTree);  RF_getTree --;
  RF_TN_SURV_ = REAL(tnSURV);
  RF_TN_MORT_ = REAL(tnMORT);
  RF_TN_NLSN_ = REAL(tnNLSN);
  RF_TN_CSHZ_ = REAL(tnCSHZ);
  RF_TN_CIFN_ = REAL(tnCIFN);
  RF_TN_REGR_ = REAL(tnREGR);  
  RF_TN_CLAS_ = (uint *) INTEGER(tnCLAS);
  processDefaultPredict();
  mode = (RF_fobservationSize > 0)? RF_PRED : RF_REST;  
  stackForestObjectsAuxOnly(mode);
  RF_parmID_[1]              = (int *) INTEGER(VECTOR_ELT(hc_zero, 0));   RF_parmID_[1]  --;
  RF_contPT_[1]              =             REAL(VECTOR_ELT(hc_zero, 1));  RF_contPT_[1]  --;
  RF_mwcpSZ_[1]              = (uint *) INTEGER(VECTOR_ELT(hc_zero, 2));  RF_mwcpSZ_[1]  --;
  RF_fsrecID_[1]             = (uint *) INTEGER(VECTOR_ELT(hc_zero, 3));  RF_fsrecID_[1] --;
  if (VECTOR_ELT(hc_zero, 2) != R_NilValue) {
    RF_mwcpPT_[1]            = (uint *) INTEGER(VECTOR_ELT(hc_zero, 4));  RF_mwcpPT_[1]  --;
  }
  else {
    RF_mwcpPT_[1] = NULL;
  }
  if (RF_baseLearnDepthINTR > 1) {
    RF_pairCT_              = (uint *) INTEGER(VECTOR_ELT(hc_oneAugIntr, 0));  RF_pairCT_ --;
    RF_augmX1_[1]           = (int *)  INTEGER(VECTOR_ELT(hc_oneAugIntr, 1));  RF_augmX1_[1]  --;
    RF_augmX2_[1]           = (int *)  INTEGER(VECTOR_ELT(hc_oneAugIntr, 2));  RF_augmX2_[1]  --;
  }
  else {
    RF_pairCT_ = NULL;
    RF_augmX1_ = NULL;
    RF_augmX2_ = NULL;
  }
  if (RF_baseLearnDepthSYTH > 1) {
    RF_lotsSZ_              = (uint *) INTEGER(VECTOR_ELT(hc_oneAugSyth, 0));  RF_lotsSZ_     --;
    RF_augmXS_[1]           = (int *)  INTEGER(VECTOR_ELT(hc_oneAugSyth, 1));  RF_augmXS_[1]  --;
  }
  else {
    RF_lotsSZ_ = NULL;
    RF_augmXS_ = NULL;
  }
  if (RF_hdim > 0) {
    RF_hcDim_                = (uint *) INTEGER(VECTOR_ELT(hc_one, 0));  RF_hcDim_  --;
    RF_contPTR_[1]           =          REAL(VECTOR_ELT(hc_one, 1));     RF_contPTR_[1] --;
  }
  else {
    RF_hcDim_      = NULL;
    RF_contPTR_    = NULL;
  }
  if (RF_hdim > 1) {
    for (i = 2; i <= RF_hdim; i++) {
      RF_parmID_[i]            = (int *)  INTEGER(VECTOR_ELT(hc_parmID, i-2));  RF_parmID_[i] --;
      RF_contPT_[i]            =          REAL(VECTOR_ELT(hc_contPT, i-2));     RF_contPT_[i] --;
      RF_contPTR_[i]           =          REAL(VECTOR_ELT(hc_contPTR, i-2));    RF_contPTR_[i] --;
      RF_mwcpSZ_[i]            = (uint *) INTEGER(VECTOR_ELT(hc_mwcpSZ, i-2));  RF_mwcpSZ_[i] --;
      RF_fsrecID_[i]           = (uint *) INTEGER(VECTOR_ELT(hc_fsrecID, i-2)); RF_fsrecID_[i] --;
      RF_mwcpPT_[i]            = (uint *) INTEGER(VECTOR_ELT(hc_mwcpPT, i-2));  RF_mwcpPT_[i] --;
      if (RF_baseLearnDepthINTR > 1) {      
        RF_augmX1_[i]            = (int *) INTEGER(VECTOR_ELT(hc_augmXone, i-2));  RF_augmX1_[i] --;
        RF_augmX2_[i]            = (int *) INTEGER(VECTOR_ELT(hc_augmXtwo, i-2));  RF_augmX2_[i] --;
      }
      if (RF_baseLearnDepthSYTH > 1) {
        RF_augmXS_[i]            = (int *) INTEGER(VECTOR_ELT(hc_augmXS, i-2));  RF_augmXS_[i] --;
      }
    }
  }
  if (RF_hdim > 0) {
    RF_nodeCountSyth_ = NULL;
    if (RF_baseLearnDepthSYTH > 1) {
      if (hc_augmSythTop != R_NilValue) {
        RF_nodeCountSyth_   = (uint *)  INTEGER(VECTOR_ELT(hc_augmSythTop, 0));  RF_nodeCountSyth_   --;
        RF_syth_treeID_     = (uint *)  INTEGER(VECTOR_ELT(hc_augmSythTop, 1));  RF_syth_treeID_     --;
        RF_syth_nodeID_     = (uint *)  INTEGER(VECTOR_ELT(hc_augmSythTop, 2));  RF_syth_nodeID_     --;
        RF_syth_hcDim_      = (uint *)  INTEGER(VECTOR_ELT(hc_augmSythTop, 3));  RF_syth_hcDim_      --;
        RF_syth_parmID_[1]  =  (int *)  INTEGER(VECTOR_ELT(hc_augmSythTop, 4));  RF_syth_parmID_[1]  --;
        RF_syth_contPT_[1]  =           REAL(VECTOR_ELT(hc_augmSythTop, 5));     RF_syth_contPT_[1]  --;
        RF_syth_contPTR_[1] =           REAL(VECTOR_ELT(hc_augmSythTop, 6));     RF_syth_contPTR_[1] --;
        RF_syth_mwcpSZ_[1]  = (uint *)  INTEGER(VECTOR_ELT(hc_augmSythTop, 7));  RF_syth_mwcpSZ_[1]  --;
        if (VECTOR_ELT(hc_augmSythTop, 8) != R_NilValue) {
          RF_syth_mwcpPT_[1] = (uint *) INTEGER(VECTOR_ELT(hc_augmSythTop, 8)); RF_syth_mwcpPT_[1]  --;
        }
        else {
          RF_syth_mwcpPT_[1] = NULL;
        }
      }
    }
  }
  rfsrc(mode, seedValue);
  unstackForestObjectsAuxOnly(mode);
  free_1DObject(RF_rType, NATIVE_TYPE_CHARACTER, RF_ySize);
  free_1DObject(RF_xType, NATIVE_TYPE_CHARACTER, RF_xSize);
  if (RF_responseIn != NULL) free_2DObject(RF_responseIn, NATIVE_TYPE_NUMERIC, RF_ySize > 0, RF_ySize, RF_observationSize);
  if (RF_observationIn != NULL) free_2DObject(RF_observationIn, NATIVE_TYPE_NUMERIC, TRUE, RF_xSize, RF_observationSize);
  free_2DObject(RF_bootstrapIn, NATIVE_TYPE_INTEGER, (RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2), RF_ntree, RF_subjSize);
  free_2DObject(RF_fresponseIn, NATIVE_TYPE_NUMERIC, RF_frSize > 0, RF_frSize, RF_fobservationSize);
  free_2DObject(RF_fobservationIn, NATIVE_TYPE_NUMERIC, RF_fobservationSize > 0 , RF_xSize, RF_fobservationSize);
  memoryCheck();
  RF_cpuTime_[1] = (double) (clock() - cpuTimeStart) / CLOCKS_PER_SEC;
  R_ReleaseObject(RF_sexpVector[RF_OUTP_ID]);
  R_ReleaseObject(RF_sexpVector[RF_STRG_ID]);
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
                      uint   auxiliaryDimSize,
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
  va_start(list, auxiliaryDimSize);
  int *auxiliaryDim = ivector(1, auxiliaryDimSize);
  for (uint i = 1; i <= auxiliaryDimSize; i++) {
    auxiliaryDim[i] = va_arg(list, int);
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
  SET_STRING_ELT(RF_sexpVector[RF_STRG_ID], *sexpIndex, mkChar(sexpString));
  UNPROTECT(1);
  switch(sexpType) {
  case NATIVE_TYPE_NUMERIC:
    v = (double*) NUMERIC_POINTER(thisVector);
    for (ulong i = 0; i < size; i++) {
      ((double*) v)[i] = value;
    }
    break;
  case NATIVE_TYPE_INTEGER:
    v = (uint*) INTEGER_POINTER(thisVector);
    for (ulong i = 0; i < size; i++) {
      ((uint*) v)[i] = 0;
    }
    break;
  case NATIVE_TYPE_CHARACTER:
    v = (char*) CHARACTER_POINTER(thisVector);
    for (ulong i = 0; i < size; i++) {
      ((char*) v)[i] = 0;
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
                        auxiliaryDimSize,
                        auxiliaryDim);
  free_ivector(auxiliaryDim, 1, auxiliaryDimSize);
  (*sexpIndex) ++;
  return v;
}
void setUserTraceFlag (uint traceFlag) {
  RF_userTraceFlag = traceFlag;
}
uint getUserTraceFlag () {
  return RF_userTraceFlag;
}
