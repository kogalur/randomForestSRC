
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "entry.h"
#include "entryGeneric.h"
#include "nativeUtil.h"
#include "rfsrc.h"
#include "stackOutput.h"
SEXP rfsrcGrow(SEXP traceFlag,
               SEXP seedPtr,
               SEXP optLow,
               SEXP optHigh,
               SEXP optSup,
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
               SEXP qStarPlus,
               SEXP xPreSort,
               SEXP numThreads) {
  clock_t cpuTimeStart = clock();
  setUserTraceFlag(INTEGER(traceFlag)[0]);
  setNativeGlobalEnv(&RF_nativeIndex, &RF_stackCount);
  int seedValue           = INTEGER(seedPtr)[0];
  RF_opt                  = INTEGER(optLow)[0];
  RF_optHigh              = INTEGER(optHigh)[0];
  RF_optSup               = INTEGER(optSup)[0];
  RF_splitRule            = INTEGER(splitRule)[0];
  RF_nsplit               = INTEGER(nsplit)[0];
  RF_mtry                 = INTEGER(mtry)[0];
  RF_ytry                 = INTEGER(ytry)[0];
  RF_nodeSize             = INTEGER(nodeSize)[0];
  RF_nodeDepth            = INTEGER(nodeDepth)[0];
  RF_crWeight             = NULL;
  RF_crWeightSize         = INTEGER(crWeightSize)[0];
  if (RF_crWeightSize > 0) {
    RF_crWeight             = REAL(crWeight); RF_crWeight--;
  }
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
  RF_qStarPlus = NULL;
  if(RF_ySize > 0) {
    if (qStarPlus != R_NilValue) {
      RF_qStarPlus         = (double **) copy2DObject(qStarPlus, NATIVE_TYPE_NUMERIC, TRUE, RF_ySize, RF_ySize);
    }
  }
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
  free_2DObject(RF_qStarPlus, NATIVE_TYPE_NUMERIC, RF_qStarPlus != NULL, RF_ySize, RF_ySize);
  free_2DObject(RF_vtryArray, NATIVE_TYPE_INTEGER, RF_vtry > 0, RF_ntree, RF_xSize);  
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
                  SEXP nodeSZ,
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
  clock_t cpuTimeStart = clock();
  setUserTraceFlag(INTEGER(traceFlag)[0]);
  setNativeGlobalEnv(&RF_nativeIndex, &RF_stackCount);
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
  if (xData != R_NilValue) {
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
  RF_treeID_              = (uint *) INTEGER(treeID);   RF_treeID_ --;
  RF_nodeID_              = (uint *) INTEGER(nodeID);   RF_nodeID_ --;
  RF_nodeSZ_              = (uint *) INTEGER(nodeSZ);   RF_nodeSZ_ --;
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
  if (RF_partialLength > 0) {
    RF_partialValue         = REAL(VECTOR_ELT(partial, 3)); RF_partialValue --;
  }
  else {
    RF_partialValue = NULL;
  }
  RF_partialLength2       = INTEGER(VECTOR_ELT(partial, 4))[0];
  if (RF_partialLength2 > 0) {
    RF_partialXvar2         = (uint *) INTEGER(VECTOR_ELT(partial, 5)); RF_partialXvar2 --;
  }
  else {
    RF_partialXvar2 = NULL;
  }
  if (RF_partialLength2 > 0) {
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
  if (VECTOR_ELT(hc_zero, 4) != R_NilValue) {
    RF_mwcpPT_[1]            = (uint *) INTEGER(VECTOR_ELT(hc_zero, 4));  RF_mwcpPT_[1]  --;
  }
  else {
    RF_mwcpPT_[1] = NULL;
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
  RF_cpuTime_[1] = (double) (clock() - cpuTimeStart) / CLOCKS_PER_SEC;
  R_ReleaseObject(RF_sexpVector[RF_OUTP_ID]);
  R_ReleaseObject(RF_sexpVector[RF_STRG_ID]);
  return RF_sexpVector[RF_OUTP_ID];
}
