
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "stackOutput.h"
#include "importance.h"
#include "splitGreedy.h"
#include "sexpOutgoing.h"
#include "nativeUtil.h"
#include "nrutil.h"
#include "error.h"
void stackDefinedOutputObjects(char      mode,
                               char    **RF_sexpString,
                               Node   ***pRF_root,
                               uint    **pRF_tLeafCount,
                               double  **pRF_proximity,
                               double  **pRF_distance,
                               double  **pRF_weight,
                               double  **p_imputation,
                               double ***pRF_sImputeResponsePtr,
                               double ***pRF_sImputePredictorPtr,
                               uint    **pRF_varUsed,
                               uint   ***pRF_varUsedPtr,
                               double  **p_splitDepth) {
  uint sexpIdentity;
  ulong localSize;
  uint xVimpSize;
  uint  obsSize;
  uint  mRecordSize;
  uint *mRecordIndex;
  double **responsePtr;
  double **predictorPtr;
  uint     rspSize;
  uint     dpthDimOne;
  double **ensembleDen;
  double **ensembleSRG;
  double **ensembleMRT;
  double **ensembleCIF;
  double **ensembleSRV;
  double **ensembleCLS;
  double **ensembleRGR;
  double **ensembleQNT;
  double ****ensembleSRGptr;
  double  ***ensembleMRTptr;
  double  ***ensembleSRVptr;
  double ****ensembleCIFptr;
  double ****ensembleCLSptr;
  double  ***ensembleRGRptr;
  double      ****ensembleQNTptr;
  uint         ***quantileStreamSize;
  LookUpInfo  ****quantileSearchTree;
  QuantileObj ****quantileHead;
  QuantileObj ****quantileTail;
  uint         ***quantileLinkLength;
  double ****ensembleSRGnum;
  double  ***ensembleMRTnum;
  double  ***ensembleSRVnum;
  double ****ensembleCIFnum;
  double ****ensembleCLSnum;
  double  ***ensembleRGRnum;
  char oobFlag, fullFlag;
  uint dimThree;
  uint i, j, k, m;
  xVimpSize      = 0;  
  dpthDimOne     = 0;  
  obsSize        = 0;  
  mRecordSize    = 0;  
  rspSize        = 0;  
  sexpIdentity   = 0;  
  responsePtr    = NULL;  
  predictorPtr   = NULL;  
  mRecordIndex   = NULL;  
  RF_identityMembershipIndexSize = (RF_bootstrapSize > RF_observationSize ) ? RF_bootstrapSize : RF_observationSize;
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    RF_identityMembershipIndexSize = (RF_bootstrapSize > RF_observationSize ) ? RF_bootstrapSize : RF_observationSize;
  }
  else {
    RF_identityMembershipIndexSize = (RF_bootstrapSize > RF_observationSize ) ? RF_bootstrapSize : RF_observationSize;
  }
  RF_identityMembershipIndex = uivector(1, RF_identityMembershipIndexSize);
  for (i = 1; i <= RF_identityMembershipIndexSize; i++) {
    RF_identityMembershipIndex[i] = i;
  }
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    mRecordSize = RF_fmRecordSize;
    rspSize = RF_frSize;
    responsePtr  = RF_fresponseIn;
    predictorPtr = RF_fobservationIn;
    mRecordIndex = RF_fmRecordIndex;
    RF_stackCount = 0;
    RF_stackCount ++;
    RF_stackCount ++;
    if (RF_opt & OPT_FENS) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        RF_stackCount += 3;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          RF_stackCount += 1;
        }
        if (RF_rTargetNonFactorCount > 0) {
          RF_stackCount += 1;
          if (RF_opt & OPT_QUANTLE) {
            RF_stackCount += 1;
          }
        }
      }
    }
    if (RF_opt & OPT_PERF) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        RF_stackCount += 1;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          RF_stackCount += 1;
          if (RF_optHigh & OPT_CSE) {
            RF_stackCount += 1;
          }
        }
        if (RF_rTargetNonFactorCount > 0) {
          RF_stackCount += 1;
          if (RF_optHigh & OPT_CSE) {
            RF_stackCount += 1;
          }
        }
        if ((RF_rTargetFactorCount > 0) || (RF_rTargetNonFactorCount > 0)) {
          if (RF_optHigh & OPT_CSE) {
            RF_stackCount += 1;
          }
        }
      }
    }
    if (RF_opt & OPT_PROX) {
      RF_stackCount += 1;
    }
    if (RF_optHigh & OPT_DIST) {
      RF_stackCount += 1;
    }
    if (RF_optHigh & OPT_WGHT) {
      RF_stackCount += 1;
    }
    if (RF_opt & OPT_MISS_OUT) {
      RF_stackCount += 1;
    }
    if (RF_opt & (OPT_SPLDPTH_1 | OPT_SPLDPTH_2)) {
      RF_stackCount += 1;
    }
    if (RF_opt & OPT_CASE_DPTH) {
      RF_stackCount += 1;
    }
    if (RF_opt & OPT_VIMP) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        RF_stackCount += 2;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          RF_stackCount += 2;
          if (RF_optHigh & OPT_CSV) {
            RF_stackCount += 1;
          }
        }
        if (RF_rTargetNonFactorCount > 0) {
          RF_stackCount += 2;
          if (RF_optHigh & OPT_CSV) {
            RF_stackCount += 1;
          }          
        }
        if ((RF_rTargetFactorCount > 0) || (RF_rTargetNonFactorCount > 0)) {
          if (RF_optHigh & OPT_CSV) {
            RF_stackCount += 1;
          }
        }
      }
    }
    if (RF_optHigh & OPT_MEMB_PRUN) {
      RF_stackCount += 1;
    }
    if (RF_optHigh & OPT_MEMB_USER) {
      RF_stackCount += 2;
    }
    break;
  default:
    obsSize = RF_observationSize;
    mRecordSize = RF_mRecordSize;
    rspSize = RF_ySize;
    responsePtr  = RF_responseIn;
    predictorPtr = RF_observationIn;
    mRecordIndex = RF_mRecordIndex;
    RF_stackCount = 0;
    RF_stackCount ++;
    if (RF_opt & OPT_LEAF) {
      RF_stackCount += 1;
    }
    if (RF_opt & OPT_FENS) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        RF_stackCount += 3;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          RF_stackCount += 1;
        }
        if (RF_rTargetNonFactorCount > 0) {
          RF_stackCount += 1;
          if (RF_opt & OPT_QUANTLE) {
            RF_stackCount += 1;
          }
        }
      }
    }
    if (RF_opt & OPT_OENS) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        RF_stackCount += 3;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          RF_stackCount += 1;
        }
        if (RF_rTargetNonFactorCount > 0) {
          RF_stackCount += 1;
          if (RF_opt & OPT_QUANTLE) {
            RF_stackCount += 1;
          }
        }
      }
    }
    if (RF_opt & OPT_PERF) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        RF_stackCount += 1;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          RF_stackCount += 1;
          if (RF_optHigh & OPT_CSE) {
            RF_stackCount += 1;
          }
        }
        if (RF_rTargetNonFactorCount > 0) {
          RF_stackCount += 1;
          if (RF_optHigh & OPT_CSE) {
            RF_stackCount += 1;
          }
        }
        if ((RF_rTargetFactorCount > 0) || (RF_rTargetNonFactorCount > 0)) {
          if (RF_optHigh & OPT_CSE) {
            RF_stackCount += 1;
          }
        }
      }
    }
    if (RF_opt & OPT_EMPR_RISK) {
      RF_stackCount += 3;
    }
    if (RF_opt & OPT_PROX) {
      RF_stackCount += 1;
    }
    if (RF_optHigh & OPT_DIST) {
      RF_stackCount += 1;
    }
    if (RF_optHigh & OPT_WGHT) {
      RF_stackCount += 1;
    }
    if (RF_opt & OPT_SEED) {
      if (RF_opt & OPT_TREE) {
        RF_stackCount += 1;
        uint bnpSize = getVimpRecoverySeedDimension(mode, RF_opt);
        if (bnpSize > 0) {
          RF_stackCount += 1;
        }
        RF_stackCount += 1;
        RF_stackCount += 4;
        RF_stackCount += 5;
        RF_stackCount += 1;
        RF_stackCount += 1;
      }
      else {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  S.E.X.P. TREE output request inconsistent.");
        RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
        RF_nativeExit();
      }
    }
    if (RF_opt & OPT_MISS_OUT) {
      RF_stackCount += 1;
    }
    if (RF_opt & (OPT_VARUSED_F | OPT_VARUSED_T)) {
      RF_stackCount += 1;
    }
    if (RF_opt & (OPT_SPLDPTH_1 | OPT_SPLDPTH_2)) {
      RF_stackCount += 1;
    }
    if (RF_opt & OPT_CASE_DPTH) {
      RF_stackCount += 1;
    }
    if (RF_opt & OPT_VIMP) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        RF_stackCount += 2;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          RF_stackCount += 2;
          if (RF_optHigh & OPT_CSV) {
            RF_stackCount += 1;
          }
        }
        if (RF_rTargetNonFactorCount > 0) {
          RF_stackCount += 2;
          if (RF_optHigh & OPT_CSV) {
            RF_stackCount += 1;
          }          
        }
        if ((RF_rTargetFactorCount > 0) || (RF_rTargetNonFactorCount > 0)) {
          if (RF_optHigh & OPT_CSV) {
            RF_stackCount += 1;
          }
        }
      }
    }
    if ((RF_vtry > 0) && (RF_vtryMode != RF_VTRY_NULL)) {
      RF_stackCount += 1;
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        RF_stackCount += 1;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          RF_stackCount += 1;
        }
        if (RF_rTargetNonFactorCount > 0) {
          RF_stackCount += 1;
        }
      }
    }
    if (RF_optHigh & OPT_MEMB_PRUN) {
      RF_stackCount += 1;
    }
    if (RF_optHigh & OPT_MEMB_USER) {
      RF_stackCount += 2;
    }
    if (RF_optHigh & OPT_MEMB_OUTG) {
      RF_stackCount += 10;
    }
    if (RF_optHigh & OPT_TERM_OUTG) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        RF_stackCount += 1;
        if (!(RF_opt & OPT_COMP_RISK)) {
          RF_stackCount += 1;
          RF_stackCount += 1;
        }
        else {
          RF_stackCount += 2;
        }
      }
      else {
        if (RF_rTargetNonFactorCount > 0) {
          RF_stackCount += 1;
        }
        if (RF_rTargetFactorCount > 0) {
          RF_stackCount += 1;
        }
      }
    }
    if (RF_optHigh & OPT_PART_PLOT) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        RF_stackCount += 1;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          RF_stackCount += 1;
        }
        if (RF_rTargetNonFactorCount > 0) {
          RF_stackCount += 1;
        }
      }
    }
    break;
  }  
  if (RF_xMarginalSize > 0) {
    if (RF_xMarginalSize <= RF_xSize) {
      for (uint i = 1; i <= RF_xMarginalSize; i++) {
        if ((RF_xMarginal[i] < 1) || (RF_xMarginal[i] > RF_xSize)) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  Parameter verification failed.");
          RF_nativeError("\nRF-SRC:  Marginal predictor must be greater than zero and less than or equal to %10d:  %10d \n", RF_xSize, RF_xMarginal[i]);
          RF_nativeExit();
        }
      }
    }
    else {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Parameter verification failed.");
      RF_nativeError("\nRF-SRC:  Number of marginal predictors must be greater than zero and less than or equal to %10d:  %10d \n", RF_xSize, RF_xMarginalSize);
      RF_nativeExit();
    }
    RF_utTermMembership      =  (uint ***) new_vvector(1, RF_ntree, NRUTIL_UPTR2);
    RF_utTermMembershipCount =  (uint **) new_vvector(1, RF_ntree, NRUTIL_UPTR);
    RF_utTermMembershipAlloc =  (uint **) new_vvector(1, RF_ntree, NRUTIL_UPTR);
    RF_xMarginalFlag         =  uivector(1, RF_xSize);
    for (i = 1; i <= RF_xSize; i++) {
      RF_xMarginalFlag[i] = FALSE;
    }
    for (i = 1; i <= RF_xMarginalSize; i++) {
      RF_xMarginalFlag[RF_xMarginal[i]] = TRUE;
    }
  }
  initProtect(RF_stackCount);
  stackAuxiliaryInfoList(&RF_snpAuxiliaryInfoList, RF_stackCount);
  oobFlag = fullFlag = FALSE;
  if ((RF_opt & OPT_FENS) || (RF_opt & OPT_OENS)) {
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    if (RF_opt & OPT_OENS) {
      oobFlag = TRUE;
    }
    while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
      ensembleDen    = NULL;
      ensembleSRG    = NULL;
      ensembleSRGptr = NULL;
      ensembleSRGnum = NULL;
      ensembleMRT    = NULL;
      ensembleMRTptr = NULL;
      ensembleMRTnum = NULL;
      ensembleSRV    = NULL;
      ensembleSRVptr = NULL;
      ensembleSRVnum = NULL;
      ensembleCIF    = NULL;
      ensembleCIFptr = NULL;
      ensembleCIFnum = NULL;
      ensembleCLS    = NULL;
      ensembleCLSptr = NULL;
      ensembleCLSnum = NULL;
      ensembleRGR    = NULL;
      ensembleRGRptr = NULL;
      ensembleRGRnum = NULL;
      ensembleQNT        = NULL;
      ensembleQNTptr     = NULL;
      quantileStreamSize = NULL;
      quantileSearchTree = NULL;
      quantileHead       = NULL;
      quantileTail       = NULL;
      quantileLinkLength = NULL;
      if (oobFlag == TRUE) {
        ensembleDen    = &RF_oobEnsembleDen;
        ensembleSRG    = &RF_oobEnsembleSRG_;
        ensembleSRGptr = &RF_oobEnsembleSRGptr;
        ensembleSRGnum = &RF_oobEnsembleSRGnum;
        ensembleMRT    = &RF_oobEnsembleMRT_;
        ensembleMRTptr = &RF_oobEnsembleMRTptr;
        ensembleMRTnum = &RF_oobEnsembleMRTnum;
        ensembleSRV    = &RF_oobEnsembleSRV_;
        ensembleSRVptr = &RF_oobEnsembleSRVptr;
        ensembleSRVnum = &RF_oobEnsembleSRVnum;
        ensembleCIF    = &RF_oobEnsembleCIF_;
        ensembleCIFptr = &RF_oobEnsembleCIFptr;
        ensembleCIFnum = &RF_oobEnsembleCIFnum;
        ensembleCLS    = &RF_oobEnsembleCLS_;
        ensembleCLSptr = &RF_oobEnsembleCLSptr;
        ensembleCLSnum = &RF_oobEnsembleCLSnum;
        ensembleRGR    = &RF_oobEnsembleRGR_;
        ensembleRGRptr = &RF_oobEnsembleRGRptr;
        ensembleRGRnum = &RF_oobEnsembleRGRnum;
        ensembleQNT         = &RF_oobEnsembleQNT_;
        ensembleQNTptr      = &RF_oobEnsembleQNTptr;        
        quantileStreamSize  = &RF_oobQuantileStreamSize;
        quantileSearchTree  = &RF_oobQuantileSearchTree;
        quantileHead        = &RF_oobQuantileHead;
        quantileTail        = &RF_oobQuantileTail;
        quantileLinkLength  = &RF_oobQuantileLinkLength;
      }
      else {
        ensembleDen    = &RF_fullEnsembleDen;
        ensembleSRG    = &RF_fullEnsembleSRG_;
        ensembleSRGptr = &RF_fullEnsembleSRGptr;
        ensembleSRGnum = &RF_fullEnsembleSRGnum;
        ensembleMRT    = &RF_fullEnsembleMRT_;
        ensembleMRTptr = &RF_fullEnsembleMRTptr;
        ensembleMRTnum = &RF_fullEnsembleMRTnum;        
        ensembleSRV    = &RF_fullEnsembleSRV_;
        ensembleSRVptr = &RF_fullEnsembleSRVptr;
        ensembleSRVnum = &RF_fullEnsembleSRVnum;
        ensembleCIF    = &RF_fullEnsembleCIF_;
        ensembleCIFptr = &RF_fullEnsembleCIFptr;
        ensembleCIFnum = &RF_fullEnsembleCIFnum;
        ensembleCLS    = &RF_fullEnsembleCLS_;
        ensembleCLSptr = &RF_fullEnsembleCLSptr;
        ensembleCLSnum = &RF_fullEnsembleCLSnum;
        ensembleRGR    = &RF_fullEnsembleRGR_;
        ensembleRGRptr = &RF_fullEnsembleRGRptr;
        ensembleRGRnum = &RF_fullEnsembleRGRnum;
        ensembleQNT         = &RF_fullEnsembleQNT_;
        ensembleQNTptr      = &RF_fullEnsembleQNTptr;        
        quantileStreamSize  = &RF_fullQuantileStreamSize;
        quantileSearchTree  = &RF_fullQuantileSearchTree;
        quantileHead        = &RF_fullQuantileHead;
        quantileTail        = &RF_fullQuantileTail;
        quantileLinkLength  = &RF_fullQuantileLinkLength;
      }
      *ensembleDen = dvector(1, obsSize);
      for (i = 1; i <= obsSize; i++) {
        (*ensembleDen)[i] = 0;
      }
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        (oobFlag == TRUE) ? (sexpIdentity = RF_OSRG_ID) : ((fullFlag == TRUE) ? sexpIdentity = RF_ASRG_ID : TRUE);
        localSize = (ulong) RF_eventTypeSize * RF_sortedTimeInterestSize * obsSize;
        *ensembleSRG = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, sexpIdentity, localSize, 0, RF_sexpString[sexpIdentity], ensembleSRGptr, 3, RF_eventTypeSize, RF_sortedTimeInterestSize, obsSize);
        *ensembleSRGnum = (double ***) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR2);
        for (j = 1; j <= RF_eventTypeSize; j++) {
          (*ensembleSRGnum)[j] = (double **) new_vvector(1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
          for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
            (*ensembleSRGnum)[j][k]  = dvector(1, obsSize);
            for (i = 1; i <= obsSize; i++) {
              (*ensembleSRGnum)[j][k][i] = 0.0;
            }
          }
        }
        (oobFlag == TRUE) ? (sexpIdentity = RF_OMRT_ID) : ((fullFlag == TRUE) ? sexpIdentity = RF_AMRT_ID: TRUE);
        localSize = (ulong) RF_eventTypeSize * obsSize;
        *ensembleMRT = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, sexpIdentity, localSize, 0, RF_sexpString[sexpIdentity], ensembleMRTptr, 2, RF_eventTypeSize, obsSize);
        *ensembleMRTnum = (double **) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR);
        for (j = 1; j <= RF_eventTypeSize; j++) {
          (*ensembleMRTnum)[j] = dvector(1, obsSize);
          for (i = 1; i <= obsSize; i++) {
            (*ensembleMRTnum)[j][i] = 0.0;
          }
        }
        if (!(RF_opt & OPT_COMP_RISK)) {
          (oobFlag == TRUE) ? (sexpIdentity = RF_OSRV_ID) : ((fullFlag == TRUE) ? sexpIdentity = RF_ASRV_ID: TRUE);
          localSize = (ulong) RF_sortedTimeInterestSize * obsSize;
          *ensembleSRV = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, sexpIdentity, localSize, 0, RF_sexpString[sexpIdentity], ensembleSRVptr, 2, RF_sortedTimeInterestSize, obsSize);
          *ensembleSRVnum = (double **) new_vvector(1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
          for (j = 1; j <= RF_sortedTimeInterestSize; j++) {
            (*ensembleSRVnum)[j]  = dvector(1, obsSize);
            for (i = 1; i <= obsSize; i++) {
              (*ensembleSRVnum)[j][i]  = 0.0;
            }
          }
        }  
        else {
          (oobFlag == TRUE) ? (sexpIdentity = RF_OCIF_ID) : ((fullFlag == TRUE) ? sexpIdentity = RF_ACIF_ID: TRUE);
          localSize = (ulong) RF_eventTypeSize * RF_sortedTimeInterestSize * obsSize;
          *ensembleCIF = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, sexpIdentity, localSize, 0, RF_sexpString[sexpIdentity], ensembleCIFptr, 3, RF_eventTypeSize, RF_sortedTimeInterestSize, obsSize);
          *ensembleCIFnum = (double ***) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR2);
          for (j = 1; j <= RF_eventTypeSize; j++) {
            (*ensembleCIFnum)[j] = (double **) new_vvector(1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
            for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
              (*ensembleCIFnum)[j][k]  = dvector(1, obsSize);
              for (i = 1; i <= obsSize; i++) {
                (*ensembleCIFnum)[j][k][i] = 0.0;
              }
            }
          }
        }  
      }  
      else {
        if (RF_rTargetFactorCount > 0) {
          (oobFlag == TRUE) ? (sexpIdentity = RF_OCLS_ID) : ((fullFlag == TRUE) ? sexpIdentity = RF_ACLS_ID: TRUE);
          localSize = 0;
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              localSize += (ulong) obsSize;
            }
          }
          *ensembleCLS = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, sexpIdentity, localSize, 0, RF_sexpString[sexpIdentity], ensembleCLSptr, 3, RF_rTargetFactorCount, 0, obsSize);
          *ensembleCLSnum = (double ***) new_vvector(1, RF_rTargetFactorCount, NRUTIL_DPTR2);
          localSize = 0;
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            (*ensembleCLSnum)[j] = (double **) new_vvector(1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]], NRUTIL_DPTR);
            for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              (*ensembleCLSnum)[j][k]  = dvector(1, obsSize);
              localSize += (ulong) obsSize;
              for (i = 1; i <= obsSize; i++) {
                (*ensembleCLSnum)[j][k][i] = 0.0;
              }
            }
          }
        }
        if (RF_rTargetNonFactorCount > 0) {
          (oobFlag == TRUE) ? (sexpIdentity = RF_ORGR_ID) : ((fullFlag == TRUE) ? sexpIdentity = RF_ARGR_ID: TRUE);
          localSize = (ulong) RF_rTargetNonFactorCount * obsSize;
          *ensembleRGR = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, sexpIdentity, localSize, 0, RF_sexpString[sexpIdentity], ensembleRGRptr, 2, RF_rTargetNonFactorCount, obsSize);
          (*ensembleRGRnum) = (double **) new_vvector(1, RF_rTargetNonFactorCount, NRUTIL_DPTR);
          for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
            (*ensembleRGRnum)[j] = dvector(1, obsSize);
            for (i = 1; i <= obsSize; i++) {
              (*ensembleRGRnum)[j][i] = 0.0;
            }
          }
          if (RF_opt & OPT_QUANTLE) {
            if ((RF_qEpsilon <= 0.0) || (RF_qEpsilon >= 0.50)) {
              RF_nativeError("\nRF-SRC:  *** ERROR *** ");
              RF_nativeError("\nRF-SRC:  Parameter verification failed.");
              RF_nativeError("\nRF-SRC:  Epsilon-approximate quantile threshold must be in the range (0, 1/2):  %10.4f \n", RF_qEpsilon);
              RF_nativeExit();
            }
            else {
              RF_inv_2qEpsilon = 1 / (2 * RF_qEpsilon);
            }
            (oobFlag == TRUE) ? (sexpIdentity = RF_OQNT_ID) : ((fullFlag == TRUE) ? sexpIdentity = RF_AQNT_ID: TRUE);
            localSize = (ulong) RF_rTargetNonFactorCount * RF_quantileSize * obsSize;
            *ensembleQNT = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, sexpIdentity, localSize, 0, RF_sexpString[sexpIdentity], ensembleQNTptr, 3, RF_rTargetNonFactorCount, RF_quantileSize, obsSize);
            *quantileStreamSize = (uint **) new_vvector(1, RF_rTargetNonFactorCount, NRUTIL_UPTR);
            for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
              (*quantileStreamSize)[j] = uivector(1, obsSize);
              for (i = 1; i <= obsSize; i++) {
                (*quantileStreamSize)[j][i] = 0;
              }
            }
            *quantileSearchTree = (LookUpInfo ***) new_vvector(1, RF_rTargetNonFactorCount, NRUTIL_SPTR2);
            for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
              (*quantileSearchTree)[j] = (LookUpInfo **) new_vvector(1, obsSize, NRUTIL_SPTR);
              for (i = 1; i <= obsSize; i++) {
                (*quantileSearchTree)[j][i] = NULL;
              }
            }
            *quantileHead = (QuantileObj ***) new_vvector(1, RF_rTargetNonFactorCount, NRUTIL_QPTR2);
            *quantileTail = (QuantileObj ***) new_vvector(1, RF_rTargetNonFactorCount, NRUTIL_QPTR2);
            *quantileLinkLength = (uint **) new_vvector(1, RF_rTargetNonFactorCount, NRUTIL_UPTR2);
            for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
              (*quantileHead)[j] = (QuantileObj **) new_vvector(1, obsSize, NRUTIL_QPTR);
              (*quantileTail)[j] = (QuantileObj **) new_vvector(1, obsSize, NRUTIL_QPTR);
              (*quantileLinkLength)[j] = uivector(1, obsSize);
              for (i = 1; i <= obsSize; i++) {
                (*quantileHead)[j][i] = NULL;
                (*quantileTail)[j][i] = NULL;
                (*quantileLinkLength)[j][i] = 0;
              }
            }
          }
        }
      }
      if (oobFlag == TRUE) {
        oobFlag = FALSE;
      }
      else {
        fullFlag = FALSE;
      }
    }  
  }
  if (RF_opt & OPT_PERF) {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      {
        localSize = (ulong) RF_ntree * RF_eventTypeSize; 
        RF_perfMRT_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_ER_SURV, localSize, RF_nativeNaN, RF_sexpString[RF_ER_SURV], &RF_perfMRTptr, 2, RF_ntree, RF_eventTypeSize);
      }
    }  
    else {
      if (RF_rTargetFactorCount > 0) {
        localSize = 0;
        for (j = 1; j <= RF_rTargetFactorCount; j++) {
          for (k = 1; k <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
            localSize += (ulong) RF_ntree;
          }
        }
        RF_perfCLS_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_ER_CLAS, localSize, RF_nativeNaN, RF_sexpString[RF_ER_CLAS], &RF_perfCLSptr, 3, RF_ntree, RF_rTargetFactorCount, -1);
        if (RF_optHigh & OPT_CSE) {
          localSize = (ulong) RF_rTargetFactorCount * obsSize;
          RF_cseNumCLS_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_CSE_CLS, localSize, 0, RF_sexpString[RF_CSE_CLS], &RF_cseNumCLSptr, 2, RF_rTargetFactorCount, obsSize);
        }
      }
      if (RF_rTargetNonFactorCount > 0) {
        localSize = (ulong) RF_ntree * RF_rTargetNonFactorCount;
        RF_perfRGR_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_ER_REGR, localSize, RF_nativeNaN, RF_sexpString[RF_ER_REGR], &RF_perfRGRptr, 2, RF_ntree, RF_rTargetNonFactorCount);
        if (RF_optHigh & OPT_CSE) {
          localSize = (ulong) RF_rTargetNonFactorCount * obsSize;
          RF_cseNumRGR_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_CSE_RGR, localSize, 0, RF_sexpString[RF_CSE_RGR], &RF_cseNumRGRptr, 2, RF_rTargetNonFactorCount, obsSize);
        } 
      }
      if ((RF_rTargetFactorCount > 0) || (RF_rTargetNonFactorCount > 0)) {
        if (RF_optHigh & OPT_CSE) {
          RF_cseDen_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_CSE_DEN, obsSize, 0, RF_sexpString[RF_CSE_DEN], &RF_cseDENptr, 1, obsSize);
        }
      }
    }
  }  
  if (RF_opt & OPT_VIMP) {
    RF_vimpMRTstd = NULL;
    RF_vimpCLSstd = NULL;
    RF_vimpRGRstd = NULL;
    RF_vimpEnsembleDen = NULL;
    RF_blkEnsembleMRTnum = NULL;
    RF_blkEnsembleCLSnum = NULL;
    RF_blkEnsembleRGRnum = NULL;
    RF_blkEnsembleDen    = NULL;
    RF_vimpMRTblk = NULL;
    RF_vimpCLSblk = NULL;
    RF_vimpRGRblk = NULL;
    RF_perfMRTblk = NULL;
    RF_perfCLSblk = NULL;
    RF_perfRGRblk = NULL;
    RF_vimpMRTptr = NULL;
    RF_vimpCLSptr = NULL;
    RF_vimpRGRptr = NULL;
  }  
  if ((RF_vtry > 0) && (RF_vtryMode != RF_VTRY_NULL)) {
    RF_perfMRTblk = NULL;
    RF_perfCLSblk = NULL;
    RF_perfRGRblk = NULL;
    RF_holdMRTptr = NULL;
    RF_holdCLSptr = NULL;
    RF_holdRGRptr = NULL;
    RF_holdoutMap = NULL;
    RF_runningHoldoutCount = NULL;
    RF_blockSerialTreeIndex = NULL; 
  }
  if (RF_opt & OPT_VIMP) {
    if (RF_opt & OPT_VIMP_JOIN) {
      xVimpSize = 1;
    }
    else {
      xVimpSize = RF_intrPredictorSize;
    }
    RF_vimpMembership = (Terminal ****) new_vvector(1, xVimpSize, NRUTIL_NPTR3);
    for (k = 1; k <= xVimpSize; k++) {
      RF_vimpMembership[k] = (Terminal ***) new_vvector(1,  RF_ntree, NRUTIL_NPTR2);
    }
    for (k = 1; k <= xVimpSize; k++) {
      for (i = 1; i <= RF_ntree; i++) {
        RF_vimpMembership[k][i] = NULL;
      }
    }
    RF_vimpEnsembleDen  = (double **) new_vvector(1, xVimpSize, NRUTIL_DPTR);
    for (j = 1; j <= xVimpSize; j++) {
      RF_vimpEnsembleDen[j] = dvector(1, obsSize);
      for (i = 1; i <= obsSize; i++) {
        RF_vimpEnsembleDen[j][i] = 0.0;
      }
    }
    RF_blkEnsembleDen = dvector(1, obsSize);
    for (i = 1; i <= obsSize; i++) {
      RF_blkEnsembleDen[i] = 0.0;
    }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      {
        localSize = (ulong) xVimpSize * RF_eventTypeSize;
        RF_vimpMRT_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_VMP_SRG, localSize, RF_nativeNaN, RF_sexpString[RF_VMP_SRG], &RF_vimpMRTptr, 2, xVimpSize, RF_eventTypeSize);
        RF_vimpMRTstd = (double ***) new_vvector(1, xVimpSize, NRUTIL_DPTR2);
        for (j = 1; j <= xVimpSize; j++) {
          RF_vimpMRTstd[j]  = (double **) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR);
          for (k = 1; k <= RF_eventTypeSize; k++) {
            RF_vimpMRTstd[j][k] = dvector(1, obsSize);
            for (m = 1; m <= obsSize; m++) {
              RF_vimpMRTstd[j][k][m] = 0.0;
            }
          }
        }
        RF_blkEnsembleMRTnum = (double **) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR);
        for (j = 1; j <= RF_eventTypeSize; j++) {
          RF_blkEnsembleMRTnum[j] = dvector(1, obsSize);
          for (k = 1; k <= obsSize; k++) {
            RF_blkEnsembleMRTnum[j][k] = 0.0;
          }
        }
        RF_vimpMRTblk = (double ***) new_vvector(1, RF_perfBlockCount, NRUTIL_DPTR2);
        for (i = 1; i <= RF_perfBlockCount; i++) {
          RF_vimpMRTblk[i]  = (double **) new_vvector(1, xVimpSize, NRUTIL_DPTR);
          for (j = 1; j <= xVimpSize; j++) {
            RF_vimpMRTblk[i][j] = dvector(1, RF_eventTypeSize);
            for (k = 1; k <= RF_eventTypeSize; k++) {
              RF_vimpMRTblk[i][j][k] = RF_nativeNaN;
            }
          }
        }
        localSize = localSize / xVimpSize * RF_perfBlockCount;
        RF_perfBlockMRT_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_BLK_SRG, localSize, RF_nativeNaN, RF_sexpString[RF_BLK_SRG], &RF_perfMRTblk, 2, RF_perfBlockCount, RF_eventTypeSize);
      }
    }  
    else {
      if (RF_rTargetFactorCount > 0) {
        localSize = 0;
        for (j = 1; j <= RF_rTargetFactorCount; j++) {
          localSize += (ulong) 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]];
        }
        localSize = (ulong) xVimpSize * localSize;
        RF_vimpCLS_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_VMP_CLS, localSize, RF_nativeNaN, RF_sexpString[RF_VMP_CLS], &RF_vimpCLSptr, 3, xVimpSize, RF_rTargetFactorCount, -1);
        RF_vimpCLSstd = (double ****) new_vvector(1, xVimpSize, NRUTIL_DPTR3);        
        for (i = 1; i <= xVimpSize; i++) {
          RF_vimpCLSstd[i] = (double ***) new_vvector(1, RF_rTargetFactorCount, NRUTIL_DPTR2);
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            RF_vimpCLSstd[i][j]  = (double **) new_vvector(1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]], NRUTIL_DPTR);
            for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              RF_vimpCLSstd[i][j][k] = dvector(1, obsSize);
              for (m = 1; m <= obsSize; m++) {
                RF_vimpCLSstd[i][j][k][m] = 0.0;
              }
            }
          }
        }
        RF_blkEnsembleCLSnum = (double ***) new_vvector(1, RF_rTargetFactorCount, NRUTIL_DPTR2);
        for (j = 1; j <= RF_rTargetFactorCount; j++) {
          RF_blkEnsembleCLSnum[j] = (double **) new_vvector(1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]], NRUTIL_DPTR);
          for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
            RF_blkEnsembleCLSnum[j][k]  = dvector(1, obsSize);
            for (m = 1; m <= obsSize; m++) {
              RF_blkEnsembleCLSnum[j][k][m]  = 0.0;
            }
          }
        }
        RF_vimpCLSblk = (double ****) new_vvector(1, RF_perfBlockCount, NRUTIL_DPTR3);
        for (i = 1; i <= RF_perfBlockCount; i++) {
          RF_vimpCLSblk[i] = (double ***) new_vvector(1, xVimpSize, NRUTIL_DPTR2);
          for (j = 1; j <= xVimpSize; j++) {            
            RF_vimpCLSblk[i][j] = (double **) new_vvector(1, RF_rTargetFactorCount, NRUTIL_DPTR);
            for (k = 1; k <= RF_rTargetFactorCount; k++) {
              RF_vimpCLSblk[i][j][k]  = dvector(1, 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[k]]]);
              for (m = 1; m <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[k]]]; m++) {
                RF_vimpCLSblk[i][j][k][m]  = RF_nativeNaN;
              }
            }
          }
        }
        localSize = localSize / xVimpSize * RF_perfBlockCount;
        RF_perfBlockCLS_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_BLK_CLS, localSize, RF_nativeNaN, RF_sexpString[RF_BLK_CLS], &RF_perfCLSblk, 3, RF_perfBlockCount, RF_rTargetFactorCount, -1);
        if (RF_optHigh & OPT_CSV) {
          localSize = (ulong) xVimpSize * RF_rTargetFactorCount * obsSize;
          RF_csvNumCLS_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_CSV_CLS, localSize, 0, RF_sexpString[RF_CSV_CLS], &RF_csvNumCLSptr, 3, xVimpSize, RF_rTargetFactorCount, obsSize);
        }
      }
      if (RF_rTargetNonFactorCount > 0) {
        localSize = (ulong) xVimpSize * RF_rTargetNonFactorCount;
        RF_vimpRGR_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_VMP_RGR, localSize, RF_nativeNaN, RF_sexpString[RF_VMP_RGR], &RF_vimpRGRptr, 2, xVimpSize, RF_rTargetNonFactorCount);
        RF_vimpRGRstd = (double ***) new_vvector(1, xVimpSize, NRUTIL_DPTR2);
        for (j = 1; j <= xVimpSize; j++) {
          RF_vimpRGRstd[j]  = (double **) new_vvector(1, RF_rTargetNonFactorCount, NRUTIL_DPTR);
          for (k = 1; k <= RF_rTargetNonFactorCount; k++) {
            RF_vimpRGRstd[j][k] = dvector(1, obsSize);
            for (m = 1; m <= obsSize; m++) {
              RF_vimpRGRstd[j][k][m] = 0.0;
            }
          }
        }
        RF_blkEnsembleRGRnum = (double **) new_vvector(1, RF_rTargetNonFactorCount, NRUTIL_DPTR);
        for (j = 1; j <= RF_rTargetNonFactorCount; j++) {            
          RF_blkEnsembleRGRnum[j] = dvector(1, obsSize);   
          for (k= 1; k <= obsSize; k++) {
            RF_blkEnsembleRGRnum[j][k] = 0.0;
          }
        }
        RF_vimpRGRblk = (double ***) new_vvector(1, RF_perfBlockCount, NRUTIL_DPTR2);
        for (i = 1; i <= RF_perfBlockCount; i++) {
          RF_vimpRGRblk[i] = (double **) new_vvector(1, xVimpSize, NRUTIL_DPTR);
          for (j = 1; j <= xVimpSize; j++) {            
            RF_vimpRGRblk[i][j] = dvector(1, RF_rTargetNonFactorCount);             
            for (k= 1; k <= RF_rTargetNonFactorCount; k++) {
              RF_vimpRGRblk[i][j][k] = RF_nativeNaN;
            }
          }
        }
        localSize = localSize / xVimpSize * RF_perfBlockCount;
        RF_perfBlockRGR_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_BLK_RGR, localSize, RF_nativeNaN, RF_sexpString[RF_BLK_RGR], &RF_perfRGRblk, 2, RF_perfBlockCount, RF_rTargetNonFactorCount);
        if (RF_optHigh & OPT_CSV) {
          localSize = (ulong) xVimpSize * RF_rTargetNonFactorCount * obsSize;
          RF_csvNumRGR_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_CSV_RGR, localSize, 0, RF_sexpString[RF_CSV_RGR], &RF_csvNumRGRptr, 3, xVimpSize, RF_rTargetNonFactorCount, obsSize);
        }
      }
      if ((RF_rTargetFactorCount > 0) || (RF_rTargetNonFactorCount > 0)) {
        if (RF_optHigh & OPT_CSV) {
          localSize = (ulong) xVimpSize * obsSize;          
          RF_csvDen_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_CSV_DEN, localSize, 0, RF_sexpString[RF_CSV_DEN], &RF_csvDENptr, 2, xVimpSize, obsSize);
        }
      }
    }
  }  
  if ((RF_vtry > 0) && (RF_vtryMode != RF_VTRY_NULL)) {
    localSize = (ulong) RF_xSize;
    RF_holdBLK_ = (uint*) stackAndProtect(mode, &RF_nativeIndex,
                                          NATIVE_TYPE_INTEGER,
                                          RF_HLDOUT_BLK,
                                          localSize,
                                          0,
                                          RF_sexpString[RF_HLDOUT_BLK],
                                          &RF_holdBLKptr,
                                          1,
                                          RF_xSize);
    RF_holdoutMap = (uint **) new_vvector(1, RF_xSize, NRUTIL_UPTR);
    uint *localBlockID = uivector(1, RF_xSize);
    uint *iter = uivector(1, RF_xSize);
    RF_blockSerialTreeIndex = (uint ***) new_vvector(1, RF_xSize, NRUTIL_UPTR2);
    RF_runningHoldoutCount = (uint **) new_vvector(1, RF_xSize, NRUTIL_UPTR);
    for (j = 1; j <= RF_xSize; j++) {
      localBlockID[j]     = 0;
      iter[j]             = 0;
      RF_blockSerialTreeIndex[j] = NULL;
      RF_runningHoldoutCount[j]  = NULL;
      RF_holdoutMap[j] = uivector(1, RF_ntree);
    }
    for (j = 1; j <= RF_xSize; j++) {
      for (i = 1; i <= RF_ntree; i++) {
        if (RF_vtryArray[i][j] > 0) {
          if (iter[j] == 0) {
            localBlockID[j] ++;
          }
          iter[j] ++;
          RF_holdoutMap[j][i] = localBlockID[j];            
          if (iter[j] < RF_vtryBlockSize) {
          }
          else {
            iter[j] = 0;
          }
        }
        else {
          RF_holdoutMap[j][i] = 0;
        }
      }
    }
    for (j = 1; j <= RF_xSize; j++) {
      if ((iter[j] > 0) && (iter[j] < RF_vtryBlockSize)) {
        RF_holdBLKptr[j] = localBlockID[j] - 1;
      }
      else {
        RF_holdBLKptr[j] = localBlockID[j];
      }
    }
    free_uivector(localBlockID, 1, RF_xSize);
    free_uivector(iter, 1, RF_xSize);
    for (j = 1; j <= RF_xSize; j++) {
      if (RF_holdBLKptr[j] > 0) {
        RF_blockSerialTreeIndex[j] = new_vvector(1, RF_holdBLKptr[j], NRUTIL_UPTR);
        RF_runningHoldoutCount[j]  = uivector(1, RF_holdBLKptr[j]);
        for (k = 1; k <= RF_holdBLKptr[j]; k++) {
          RF_blockSerialTreeIndex[j][k] = uivector(1, RF_vtryBlockSize);
          for (m = 1; m <= RF_vtryBlockSize; m++) {
            RF_blockSerialTreeIndex[j][k][m] = 0; 
          }
          RF_runningHoldoutCount[j][k] = 0;
        }
      }
    }
    for (j = 1; j <= RF_xSize; j++) {
      for (i = 1; i <= RF_ntree; i++) {
        if ((RF_vtryArray[i][j] > 0) && (RF_holdoutMap[j][i] <= RF_holdBLKptr[j])) {
          RF_runningHoldoutCount[j][RF_holdoutMap[j][i]] ++;
          RF_blockSerialTreeIndex[j][RF_holdoutMap[j][i]][RF_runningHoldoutCount[j][RF_holdoutMap[j][i]]] = i;
        }
      }
    }
    for (j = 1; j <= RF_xSize; j++) {
      if (RF_holdBLKptr[j] > 0) {
        for (k = 1; k <= RF_holdBLKptr[j]; k++) {
          RF_runningHoldoutCount[j][k] = 0;
        }
      }
    }
    xVimpSize = RF_xSize;
    RF_holdEnsembleDen  = (double ***) new_vvector(1, xVimpSize, NRUTIL_DPTR2);
    for (j = 1; j <= xVimpSize; j++) {
      if (RF_holdBLKptr[j] > 0) {
        RF_holdEnsembleDen[j] = (double **) new_vvector(1, RF_holdBLKptr[j], NRUTIL_DPTR);
        for (k = 1; k <= RF_holdBLKptr[j]; k++) {
          RF_holdEnsembleDen[j][k] = NULL;
        }
      }
      else {
        RF_holdEnsembleDen[j] = NULL;
      }
    }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      {
        localSize = 0;
        for (j = 1; j <= RF_xSize; j++) {
          if (RF_holdBLKptr[j] > 0) {
            localSize += (ulong) RF_holdBLKptr[j] * RF_eventTypeSize;
          }
        }
        RF_holdMRT_ = (double*) stackAndProtect(mode,
                                                &RF_nativeIndex,
                                                NATIVE_TYPE_NUMERIC,
                                                RF_HLDOUT_SRG,
                                                localSize,
                                                RF_nativeNaN,
                                                RF_sexpString[RF_HLDOUT_SRG],
                                                &RF_holdMRTptr,
                                                3,
                                                xVimpSize,
                                                -3,  
                                                RF_eventTypeSize);
        RF_holdMRTstd = (double ****) new_vvector(1, xVimpSize, NRUTIL_DPTR3);
        for (j = 1; j <= xVimpSize; j++) {
          if (RF_holdBLKptr[j] > 0) {
            RF_holdMRTstd[j] = (double ***) new_vvector(1, RF_holdBLKptr[j], NRUTIL_DPTR2);
            for (k = 1; k <= RF_holdBLKptr[j]; k++) {
              RF_holdMRTstd[j][k] = NULL;
            }
          }
          else {
            RF_holdMRTstd[j]  = NULL;
          }
        }
      }
    }  
    else {
      if (RF_rTargetFactorCount > 0) {
        ulong localSize2 = 0;
        for (i = 1; i <= RF_rTargetFactorCount; i++) {
          localSize2 += (ulong) 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[i]]];
        }
        localSize = 0;
        for (j = 1; j <= RF_xSize; j++) {
          if (RF_holdBLKptr[j] > 0) {
            localSize += (ulong) RF_holdBLKptr[j] * localSize2;
          }
        }
        RF_holdCLS_ = (double*) stackAndProtect(mode,
                                                &RF_nativeIndex,
                                                NATIVE_TYPE_NUMERIC,
                                                RF_HLDOUT_CLS,
                                                localSize,
                                                RF_nativeNaN,
                                                RF_sexpString[RF_HLDOUT_CLS],
                                                &RF_holdCLSptr,
                                                4,
                                                xVimpSize,
                                                -3,   
                                                RF_rTargetFactorCount,
                                                -1);  
        RF_holdCLSstd = (double *****) new_vvector(1, xVimpSize, NRUTIL_DPTR4);        
        for (j = 1; j <= xVimpSize; j++) {
          if (RF_holdBLKptr[j] > 0) {
            RF_holdCLSstd[j] = (double ****) new_vvector(1, RF_holdBLKptr[j], NRUTIL_DPTR3);
            for (k = 1; k <= RF_holdBLKptr[j]; k++) {
              RF_holdCLSstd[j][k] = NULL;
            }
          }
          else {
            RF_holdCLSstd[j]  = NULL;
          }
        }
      }
      if (RF_rTargetNonFactorCount > 0) {
        localSize = 0;
        for (j = 1; j <= RF_xSize; j++) {
          if (RF_holdBLKptr[j] > 0) {
            localSize += (ulong) RF_holdBLKptr[j] * RF_rTargetNonFactorCount;
          }
        }
        RF_holdRGR_ = (double*) stackAndProtect(mode,
                                                &RF_nativeIndex,
                                                NATIVE_TYPE_NUMERIC,
                                                RF_HLDOUT_RGR,
                                                localSize,
                                                RF_nativeNaN,
                                                RF_sexpString[RF_HLDOUT_RGR],
                                                &RF_holdRGRptr,
                                                3,
                                                xVimpSize,
                                                -3,  
                                                RF_rTargetNonFactorCount);
        RF_holdRGRstd = (double ****) new_vvector(1, xVimpSize, NRUTIL_DPTR3);
        for (j = 1; j <= xVimpSize; j++) {
          if (RF_holdBLKptr[j] > 0) {
            RF_holdRGRstd[j] = (double ***) new_vvector(1, RF_holdBLKptr[j], NRUTIL_DPTR2);
            for (k = 1; k <= RF_holdBLKptr[j]; k++) {
              RF_holdRGRstd[j][k] = NULL;
            }
          }
          else {
            RF_holdRGRstd[j]  = NULL;
          }
        }
      }
    }
  }
  if (RF_opt & OPT_PROX) {
    localSize = ((ulong) (obsSize + 1) * obsSize) >> 1;
    *pRF_proximity = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_PROX_ID, localSize, 0, RF_sexpString[RF_PROX_ID], NULL, 1, localSize);
    RF_proximityDen = dvector(1, localSize);
    (*pRF_proximity) --;
    RF_proximityPtr = (double **) new_vvector(1, obsSize, NRUTIL_DPTR);
    RF_proximityDenPtr = (double **) new_vvector(1, obsSize, NRUTIL_DPTR);
    RF_proximityPtr[1] = *pRF_proximity;
    RF_proximityDenPtr[1] = RF_proximityDen;    
    RF_proximityPtr[1][1] = RF_proximityDenPtr[1][1] = 0.0;
    for (i = 2; i <= obsSize; i++) {
      RF_proximityPtr[i] = RF_proximityPtr[i-1] + i - 1;
      RF_proximityDenPtr[i] = RF_proximityDenPtr[i-1] + i - 1;
      for (j = 1; j <= i; j++) {
        RF_proximityPtr[i][j] = 0.0;
        RF_proximityDenPtr[i][j] = 0.0;
      }
    }
  }
  if (RF_optHigh & OPT_DIST) {
    localSize = ((ulong) (obsSize + 1) * obsSize) >> 1;
    *pRF_distance = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_DIST_ID, localSize, 0, RF_sexpString[RF_DIST_ID], NULL, 1, localSize);
    RF_distanceDen = dvector(1, localSize);
    (*pRF_distance) --;
    RF_distancePtr = (double **) new_vvector(1, obsSize, NRUTIL_DPTR);
    RF_distanceDenPtr = (double **) new_vvector(1, obsSize, NRUTIL_DPTR);
    RF_distancePtr[1] = *pRF_distance;
    RF_distanceDenPtr[1] = RF_distanceDen;    
    RF_distancePtr[1][1] = RF_distanceDenPtr[1][1] = 0.0;
    for (i = 2; i <= obsSize; i++) {
      RF_distancePtr[i] = RF_distancePtr[i-1] + i - 1;
      RF_distanceDenPtr[i] = RF_distanceDenPtr[i-1] + i - 1;
      for (j = 1; j <= i; j++) {
        RF_distancePtr[i][j] = 0.0;
        RF_distanceDenPtr[i][j] = 0.0;
      }
    }
  }
  if (RF_optHigh & OPT_WGHT) {
    localSize = (ulong) obsSize * RF_observationSize;
    *pRF_weight = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_WGHT_ID, localSize, 0, RF_sexpString[RF_WGHT_ID], &RF_weightPtr, 2, obsSize, RF_observationSize);
    RF_weightDenom = uivector(1, obsSize);
    for (k = 1; k <= obsSize; k++) {
      RF_weightDenom[k] = 0.0;
    }
  }
  if (RF_opt & OPT_LEAF) {
    localSize = (ulong) RF_ntree;
    *pRF_tLeafCount = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_LEAF_CT, localSize, 0, RF_sexpString[RF_LEAF_CT], NULL, 1, localSize);
    (*pRF_tLeafCount) --;
    if (mode == RF_GROW) {
      RF_tLeafCount = *pRF_tLeafCount;
    }
    else {
      for (i = 1; i <= RF_ntree; i++) {
        (*pRF_tLeafCount)[i] = RF_tLeafCount[i];
      }
    }
  }
  if (RF_opt & OPT_SEED) {
    localSize = (ulong) RF_ntree;
    RF_seed_ = (int*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_SEED_ID, localSize, 0, RF_sexpString[RF_SEED_ID], NULL, 1, localSize);
    RF_seed_ --;
    for (i = 1; i <= RF_ntree; i++) {
      RF_seed_[i] = -1;
    }
    uint bnpSize = getVimpRecoverySeedDimension(mode, RF_opt);
    if (bnpSize > 0) {
      localSize = (ulong) bnpSize;
      RF_seedVimp_ = (int*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_SEED_VM, localSize, 0, RF_sexpString[RF_SEED_VM], NULL, 1, localSize);
      RF_seedVimp_ --;
      for (i = 1; i <= bnpSize; i++) {
        RF_seedVimp_[i] = -1;
      }
    }
    RF_optLoGrow_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_OPT_LO_GROW, 1, 0, RF_sexpString[RF_OPT_LO_GROW], NULL, 1, 1);
    RF_optLoGrow_ --;
    RF_optLoGrow_[1] = RF_optLoGrow = RF_opt;
  }
  if (RF_opt & OPT_MISS_OUT) {
    localSize = (ulong) (1 + rspSize + RF_xSize) * mRecordSize;
    *p_imputation = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_MISS_ID, localSize, RF_nativeNaN, RF_sexpString[RF_MISS_ID], &RF_sImputeDataPtr, 2, 1 + rspSize + RF_xSize, mRecordSize);
    if (rspSize > 0) {
      *pRF_sImputeResponsePtr = (double **) new_vvector(1, rspSize, NRUTIL_DPTR);
      for (i = 1; i <= rspSize; i++) {
        (*pRF_sImputeResponsePtr)[i]  = (*p_imputation)  + (i * mRecordSize) - 1;
      }
    }
    *pRF_sImputePredictorPtr = (double **) new_vvector(1, RF_xSize, NRUTIL_DPTR);
    for (i = 1; i <= RF_xSize; i++) {
      (*pRF_sImputePredictorPtr)[i]  = (*p_imputation)  + ((rspSize + i) * mRecordSize) - 1;
    }
    for (i = 1; i <= mRecordSize; i++) {
      (*p_imputation)[i-1] = (double) mRecordIndex[i];
      if (rspSize > 0) {
        for (j = 1; j <= rspSize; j++) {
          (*pRF_sImputeResponsePtr)[j][i] = responsePtr[j][mRecordIndex[i]];
        }
      }
      for (j = 1; j <= RF_xSize; j++) {
        (*pRF_sImputePredictorPtr)[j][i] = predictorPtr[j][mRecordIndex[i]];
      }
    }
  }
  if (RF_opt & (OPT_VARUSED_F | OPT_VARUSED_T)) {
    if (RF_opt & OPT_VARUSED_T) {
      localSize = (ulong) RF_ntree * RF_xSize;
      *pRF_varUsed = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_VUSE_ID, localSize, 0, RF_sexpString[RF_VUSE_ID], pRF_varUsedPtr, 2, RF_ntree, RF_xSize);
    }
    else {
      localSize = (ulong) 1 * RF_xSize;
      *pRF_varUsed = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_VUSE_ID, localSize, 0, RF_sexpString[RF_VUSE_ID], NULL, 1, localSize);
      *pRF_varUsedPtr = uimatrix(1, RF_ntree, 1, RF_xSize);
      for (i = 1; i <= RF_ntree; i++) {
        for (j = 1; j <= RF_xSize; j++) {
          (*pRF_varUsedPtr)[i][j] = 0;
        }
      }
    }
    (*pRF_varUsed) --;
  }
  if (RF_opt & (OPT_SPLDPTH_1 | OPT_SPLDPTH_2)) {
    if (RF_opt & OPT_SPLDPTH_1) {
      dpthDimOne = 1;
    }
    else {
      dpthDimOne = RF_ntree;
    }
    localSize = (ulong) dpthDimOne * RF_xSize * RF_observationSize;
    *p_splitDepth = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_DPTH_ID, localSize, 0, RF_sexpString[RF_DPTH_ID], &RF_splitDepthPtr, 3, dpthDimOne, RF_xSize, RF_observationSize);
  }
  if (RF_opt & OPT_CASE_DPTH) {
    localSize = (ulong) RF_ntree * obsSize;
    RF_CASE_DPTH_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_CASE_DEPTH, localSize, 0, RF_sexpString[RF_CASE_DEPTH], &RF_CASE_DPTH_ptr, 2, RF_ntree, obsSize);
  }
  if (RF_optHigh & OPT_MEMB_PRUN) {
    localSize = (ulong) RF_ntree * obsSize;
    RF_PRUN_ID_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_PRUN_ID, localSize, 0, RF_sexpString[RF_PRUN_ID], &RF_PRUN_ID_ptr, 2, RF_ntree, obsSize);
  }
  if (RF_optHigh & OPT_MEMB_USER) {
    localSize = (ulong) RF_ntree * obsSize;
    RF_MEMB_ID_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_MEMB_ID, localSize, 0, RF_sexpString[RF_MEMB_ID], &RF_MEMB_ID_ptr, 2, RF_ntree, obsSize);
    localSize = (ulong) RF_ntree * RF_observationSize;
    RF_BOOT_CT_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_BOOT_CT, localSize, 0, RF_sexpString[RF_BOOT_CT], &RF_BOOT_CT_ptr, 2, RF_ntree, RF_observationSize);
  }
  if (RF_optHigh & OPT_PART_PLOT) {
    RF_partSURVptr = NULL;
    RF_partCLASptr = NULL;
    RF_partREGRptr = NULL;
    if ((RF_partialXvar < 1) || (RF_partialXvar > RF_xSize)) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Partial x-variable is out of range:  %10d ", RF_partialXvar);
      RF_nativeExit();
    }
    if (RF_partialLength2 > 0) {
      for (i = 1; i <= RF_partialLength2; i++) {
        if ((RF_partialXvar2[i] < 1) || (RF_partialXvar2[i] > RF_xSize)) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  Second order partial x-variable is out of range:  idx = %10d, val =  %10d ", i, RF_partialXvar2[i]);
          RF_nativeExit();
        }
        if (RF_partialXvar2[i] == RF_partialXvar) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  First and Second order partial x-variables are identical:  idx = %10d, val =  %10d ", i, RF_partialXvar2[i]);
          RF_nativeExit();
        }
        for (j = i + 1; j <= RF_partialLength2; j++) {
          if (RF_partialXvar2[i] == RF_partialXvar2[j]) {
            RF_nativeError("\nRF-SRC:  *** ERROR *** ");
            RF_nativeError("\nRF-SRC:  Second order partial x-variables are not unique:  idx = %10d, idx =  %10d, val = %10d", i, j, RF_partialXvar2[i]);
            RF_nativeExit();
          }
        }
      }
    }
    RF_partMembership = (Terminal ***) new_vvector(1, RF_ntree, NRUTIL_TPTR2);
    for (i = 1; i <= RF_ntree; i++) {
      RF_partMembership[i] = NULL;
    }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      {
        RF_partialTimeLength = RF_sortedTimeInterestSize;
        RF_partialTime = RF_timeInterest;
        if ((RF_partialTimeLength < 1) || (RF_partialTimeLength > RF_timeInterestSize)) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  Partial time length is out of range:  %10d ", RF_partialTimeLength);
          RF_nativeExit();
        }
        if (RF_eventTypeSize > 1) {
          if ((RF_partialType != RF_PART_YRLS) &&
              (RF_partialType != RF_PART_CIFN) &&
              (RF_partialType != RF_PART_CHFN)) {
            RF_nativeError("\nRF-SRC:  *** ERROR *** ");
            RF_nativeError("\nRF-SRC:  Partial type is out of range:  %10d ", RF_partialType);
            RF_nativeExit();
          }
        }
        else {
          if ((RF_partialType != RF_PART_MORT) &&
              (RF_partialType != RF_PART_NLSN) &&
              (RF_partialType != RF_PART_SURV)) {
            RF_nativeError("\nRF-SRC:  *** ERROR *** ");
            RF_nativeError("\nRF-SRC:  Partial type is out of range:  %10d ", RF_partialType);
            RF_nativeExit();
          }
        }
        dimThree = RF_partialTimeLength;
        if (!(RF_opt & OPT_COMP_RISK)) {
          if (RF_partialType == RF_PART_MORT) {
            dimThree = 1;
          }
        }
        else {
          if (RF_partialType == RF_PART_YRLS) {
            dimThree = 1;
          }
        }
        localSize = (ulong) RF_partialLength * RF_eventTypeSize * dimThree * RF_observationSize;
        RF_partial_SURV_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_PART_SR, localSize, 0, RF_sexpString[RF_PART_SR], &RF_partSURVptr, 4, RF_partialLength, RF_eventTypeSize, dimThree, RF_observationSize);
      }
    }  
    else {
      if (RF_rTargetFactorCount > 0) {
        localSize = 0;
        for (k = 1; k <= RF_rTargetFactorCount; k++) {
          localSize += (ulong) 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[k]]];
        }
        localSize = (ulong) RF_partialLength * localSize * RF_observationSize;
        RF_partial_CLAS_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_PART_CL, localSize, 0, RF_sexpString[RF_PART_CL], &RF_partCLASptr, 4, RF_partialLength, RF_rTargetFactorCount, -1, RF_observationSize);
      }
      if (RF_rTargetNonFactorCount > 0) {
        localSize = (ulong) RF_partialLength * RF_rTargetNonFactorCount * RF_observationSize;
        RF_partial_REGR_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_PART_RG, localSize, 0, RF_sexpString[RF_PART_RG], &RF_partREGRptr, 3, RF_partialLength, RF_rTargetNonFactorCount, RF_observationSize);
      }
    }
  }  
  RF_cpuTime_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_CPU_TIME, 1, 0, RF_sexpString[RF_CPU_TIME], NULL, 1, 1);
  RF_cpuTime_ --;
}
void unstackDefinedOutputObjects(char mode) {
  uint obsSize;
  uint xVimpSize;
  ulong localSize;
  char oobFlag, fullFlag;
  uint rspSize;
  double        **ensembleDen;
  uint         ***quantileStreamSize;
  LookUpInfo  ****quantileSearchTree;
  QuantileObj ****quantileHead;
  QuantileObj ****quantileTail;
  uint         ***quantileLinkLength;
  double ****ensembleSRGnum;
  double  ***ensembleMRTnum;
  double  ***ensembleSRVnum;
  double ****ensembleCIFnum;
  double ****ensembleCLSnum;
  double  ***ensembleRGRnum;
  uint i, j, k;
  obsSize        = 0;  
  xVimpSize      = 0;  
  rspSize        = 0;  
  free_uivector(RF_identityMembershipIndex, 1, RF_identityMembershipIndexSize);
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    rspSize = RF_frSize;
    break;
  default:
    obsSize = RF_observationSize;
    rspSize = RF_ySize;
    break;
  }
  if (RF_xMarginalSize > 0) {
    free_new_vvector(RF_utTermMembership,      1, RF_ntree, NRUTIL_UPTR2);
    free_new_vvector(RF_utTermMembershipCount, 1, RF_ntree, NRUTIL_UPTR);
    free_new_vvector(RF_utTermMembershipAlloc, 1, RF_ntree, NRUTIL_UPTR);
    free_uivector(RF_xMarginalFlag, 1, RF_xSize);
  }
  oobFlag = fullFlag = FALSE;
  if ((RF_opt & OPT_FENS) || (RF_opt & OPT_OENS)) {
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    if (RF_opt & OPT_OENS) {
      oobFlag = TRUE;
    }
    while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
      if (oobFlag == TRUE) {
        ensembleDen    = &RF_oobEnsembleDen;
        ensembleSRGnum = &RF_oobEnsembleSRGnum;
        ensembleMRTnum = &RF_oobEnsembleMRTnum;
        ensembleSRVnum = &RF_oobEnsembleSRVnum;
        ensembleCIFnum = &RF_oobEnsembleCIFnum;
        ensembleCLSnum = &RF_oobEnsembleCLSnum;
        ensembleRGRnum = &RF_oobEnsembleRGRnum;
        quantileStreamSize = &RF_oobQuantileStreamSize;
        quantileSearchTree = &RF_oobQuantileSearchTree;
        quantileHead       = &RF_oobQuantileHead;
        quantileTail       = &RF_oobQuantileTail;
        quantileLinkLength = &RF_oobQuantileLinkLength;
      }
      else {
        ensembleDen    = &RF_fullEnsembleDen;
        ensembleSRGnum = &RF_fullEnsembleSRGnum;
        ensembleMRTnum = &RF_fullEnsembleMRTnum;
        ensembleSRVnum = &RF_fullEnsembleSRVnum;
        ensembleCIFnum = &RF_fullEnsembleCIFnum;
        ensembleCLSnum = &RF_fullEnsembleCLSnum;
        ensembleRGRnum = &RF_fullEnsembleRGRnum;
        quantileStreamSize = &RF_fullQuantileStreamSize;
        quantileSearchTree = &RF_fullQuantileSearchTree;
        quantileHead       = &RF_fullQuantileHead;
        quantileTail       = &RF_fullQuantileTail;
        quantileLinkLength = &RF_fullQuantileLinkLength;
      }
      free_dvector(*ensembleDen, 1, obsSize);
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        for (j = 1; j <= RF_eventTypeSize; j++) {
          for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
            free_dvector((*ensembleSRGnum)[j][k], 1, obsSize);
          }
          free_new_vvector((*ensembleSRGnum)[j], 1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
        }
        free_new_vvector(*ensembleSRGnum, 1, RF_eventTypeSize, NRUTIL_DPTR2);
        for (j = 1; j <= RF_eventTypeSize; j++) {
          free_dvector((*ensembleMRTnum)[j], 1, obsSize);
        }
        free_new_vvector(*ensembleMRTnum, 1, RF_eventTypeSize, NRUTIL_DPTR);
        if (!(RF_opt & OPT_COMP_RISK)) {
          for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
            free_dvector((*ensembleSRVnum)[k], 1, obsSize);
          }
          free_new_vvector(*ensembleSRVnum, 1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
        }
        else {
          for (j = 1; j <= RF_eventTypeSize; j++) {
            for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
              free_dvector((*ensembleCIFnum)[j][k], 1, obsSize);
            }
            free_new_vvector((*ensembleCIFnum)[j], 1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
          }
          free_new_vvector(*ensembleCIFnum, 1, RF_eventTypeSize, NRUTIL_DPTR2);            
        }  
      }  
      else {
        if (RF_rTargetFactorCount > 0) {
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              free_dvector((*ensembleCLSnum)[j][k], 1, obsSize);
            }
            free_new_vvector((*ensembleCLSnum)[j], 1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]], NRUTIL_DPTR);
          }
          free_new_vvector(*ensembleCLSnum, 1, RF_rTargetFactorCount, NRUTIL_DPTR2);
        }
        if (RF_rTargetNonFactorCount > 0) {
          for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
            free_dvector((*ensembleRGRnum)[j], 1, obsSize);
          }
          free_new_vvector((*ensembleRGRnum), 1, RF_rTargetNonFactorCount, NRUTIL_DPTR);        
          if (RF_opt & OPT_QUANTLE) {
            for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
              free_uivector((*quantileStreamSize)[j], 1, obsSize);
            }
            free_new_vvector(*quantileStreamSize, 1, RF_rTargetNonFactorCount, NRUTIL_UPTR);            
            for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
              for (i = 1; i <= obsSize; i++) {
                freeLookUpTree((*quantileSearchTree)[j][i]);
              }
              free_new_vvector((*quantileSearchTree)[j], 1, obsSize, NRUTIL_SPTR);
            }
            free_new_vvector(*quantileSearchTree, 1, RF_rTargetNonFactorCount, NRUTIL_SPTR2);
            for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
              for (i = 1; i <= obsSize; i++) {
                freeQuantileObjList((*quantileHead)[j][i]);
              }
              free_new_vvector((*quantileHead)[j], 1, obsSize, NRUTIL_QPTR);
              free_new_vvector((*quantileTail)[j], 1, obsSize, NRUTIL_QPTR);
              free_uivector((*quantileLinkLength)[j], 1, obsSize);
            }
            free_new_vvector(*quantileHead, 1, RF_rTargetNonFactorCount, NRUTIL_QPTR2);
            free_new_vvector(*quantileTail, 1, RF_rTargetNonFactorCount, NRUTIL_QPTR2);
            free_new_vvector(*quantileLinkLength, 1, RF_rTargetNonFactorCount, NRUTIL_UPTR2);
          }
        }
      }
      if (oobFlag == TRUE) {
        oobFlag = FALSE;
      }
      else {
        fullFlag = FALSE;
      }
    }  
  }
  if (RF_opt & OPT_PERF) {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      {
      }
    }  
    else {
      if (RF_rTargetFactorCount > 0) {
      }
      if (RF_rTargetNonFactorCount > 0) {
      }
    }
  }  
  if (RF_opt & OPT_VIMP) {
    if (RF_opt & OPT_VIMP_JOIN) {
      xVimpSize = 1;
    }
    else {
      xVimpSize = RF_intrPredictorSize;
    }
    for (k = 1; k <= xVimpSize; k++) {
      free_new_vvector(RF_vimpMembership[k], 1,  RF_ntree, NRUTIL_NPTR2);
    }
    free_new_vvector(RF_vimpMembership, 1, xVimpSize, NRUTIL_NPTR3);
    for (j = 1; j <= xVimpSize; j++) {
      free_dvector(RF_vimpEnsembleDen[j], 1, obsSize);
    }
    free_new_vvector(RF_vimpEnsembleDen, 1, xVimpSize, NRUTIL_DPTR);
    free_dvector(RF_blkEnsembleDen, 1, obsSize);
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      {
        for (j = 1; j <= xVimpSize; j++) {
          for (k = 1; k <= RF_eventTypeSize; k++) {
            free_dvector(RF_vimpMRTstd[j][k], 1, obsSize);
          }
          free_new_vvector(RF_vimpMRTstd[j], 1, RF_eventTypeSize, NRUTIL_DPTR);
        }
        free_new_vvector(RF_vimpMRTstd, 1, xVimpSize, NRUTIL_DPTR2);
        for (j = 1; j <= RF_eventTypeSize; j++) {
          free_dvector(RF_blkEnsembleMRTnum[j], 1, obsSize);
        }
        free_new_vvector(RF_blkEnsembleMRTnum, 1, RF_eventTypeSize, NRUTIL_DPTR);
        for (i = 1; i <= RF_perfBlockCount; i++) {
          for (j = 1; j <= xVimpSize; j++) {
            free_dvector(RF_vimpMRTblk[i][j], 1, RF_eventTypeSize);
          }
          free_new_vvector(RF_vimpMRTblk[i], 1, xVimpSize, NRUTIL_DPTR);
        }
        free_new_vvector(RF_vimpMRTblk, 1, RF_perfBlockCount, NRUTIL_DPTR2);
      }
    }  
    else {
      if (RF_rTargetFactorCount > 0) {
        for (i = 1; i <= xVimpSize; i++) {
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              free_dvector(RF_vimpCLSstd[i][j][k], 1, obsSize);
            }
            free_new_vvector(RF_vimpCLSstd[i][j], 1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]], NRUTIL_DPTR);
          }
          free_new_vvector(RF_vimpCLSstd[i], 1, RF_rTargetFactorCount, NRUTIL_DPTR2);
        }
        free_new_vvector(RF_vimpCLSstd, 1, xVimpSize, NRUTIL_DPTR3);        
        for (j = 1; j <= RF_rTargetFactorCount; j++) {
          for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
            free_dvector(RF_blkEnsembleCLSnum[j][k], 1, obsSize);
          }
          free_new_vvector(RF_blkEnsembleCLSnum[j], 1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]], NRUTIL_DPTR);
        }
        free_new_vvector(RF_blkEnsembleCLSnum, 1, RF_rTargetFactorCount, NRUTIL_DPTR2);
        for (i = 1; i <= RF_perfBlockCount; i++) {
          for (j = 1; j <= xVimpSize; j++) {            
            for (k = 1; k <= RF_rTargetFactorCount; k++) {
              free_dvector(RF_vimpCLSblk[i][j][k], 1, 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[k]]]);
            }
            free_new_vvector(RF_vimpCLSblk[i][j], 1, RF_rTargetFactorCount, NRUTIL_DPTR);
          }
          free_new_vvector(RF_vimpCLSblk[i], 1, xVimpSize, NRUTIL_DPTR2);
        }
        free_new_vvector(RF_vimpCLSblk, 1, RF_perfBlockCount, NRUTIL_DPTR3);
      }
      if (RF_rTargetNonFactorCount > 0) {
        for (j = 1; j <= xVimpSize; j++) {
          for (k = 1; k <= RF_rTargetNonFactorCount; k++) {
            free_dvector(RF_vimpRGRstd[j][k], 1, obsSize);
          }
          free_new_vvector(RF_vimpRGRstd[j], 1, RF_rTargetNonFactorCount, NRUTIL_DPTR);
        }
        free_new_vvector(RF_vimpRGRstd, 1, xVimpSize, NRUTIL_DPTR2);
        for (j = 1; j <= RF_rTargetNonFactorCount; j++) {            
          free_dvector(RF_blkEnsembleRGRnum[j], 1, obsSize);   
        }
        free_new_vvector(RF_blkEnsembleRGRnum, 1, RF_rTargetNonFactorCount, NRUTIL_DPTR);
        for (i = 1; i <= RF_perfBlockCount; i++) {
          for (j = 1; j <= xVimpSize; j++) {            
            free_dvector(RF_vimpRGRblk[i][j], 1, RF_rTargetNonFactorCount);             
          }
          free_new_vvector(RF_vimpRGRblk[i], 1, xVimpSize, NRUTIL_DPTR);
        }
        free_new_vvector(RF_vimpRGRblk, 1, RF_perfBlockCount, NRUTIL_DPTR2);
      }
    }
  }  
  if ((RF_vtry > 0) && (RF_vtryMode != RF_VTRY_NULL)) {
    xVimpSize = RF_xSize;
    for (j = 1; j <= xVimpSize; j++) {
      if (RF_holdBLKptr[j] > 0) {
        free_new_vvector(RF_holdEnsembleDen[j], 1, RF_holdBLKptr[j], NRUTIL_DPTR);
      }
    }
    free_new_vvector(RF_holdEnsembleDen, 1, xVimpSize, NRUTIL_DPTR2);
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      {
        for (j = 1; j <= xVimpSize; j++) {
          if (RF_holdBLKptr[j] > 0) {
            free_new_vvector(RF_holdMRTstd[j], 1, RF_holdBLKptr[j], NRUTIL_DPTR2);
          }
        }
        free_new_vvector(RF_holdMRTstd, 1, xVimpSize, NRUTIL_DPTR3);
      }
    }  
    else {
      if (RF_rTargetFactorCount > 0) {
        for (j = 1; j <= xVimpSize; j++) {
          if (RF_holdBLKptr[j] > 0) {
            free_new_vvector(RF_holdCLSstd[j], 1, RF_holdBLKptr[j], NRUTIL_DPTR3);
          }
        }
        free_new_vvector(RF_holdCLSstd, 1, xVimpSize, NRUTIL_DPTR4);        
      }
      if (RF_rTargetNonFactorCount > 0) {
        for (j = 1; j <= xVimpSize; j++) {
          if (RF_holdBLKptr[j] > 0) {
            free_new_vvector(RF_holdRGRstd[j], 1, RF_holdBLKptr[j], NRUTIL_DPTR2);
          }
        }
        free_new_vvector(RF_holdRGRstd, 1, xVimpSize, NRUTIL_DPTR3);
      }
    }
      for (j = 1; j <= RF_xSize; j++) {
        if (RF_holdBLKptr[j] > 0) {
          free_uivector(RF_runningHoldoutCount[j], 1, RF_holdBLKptr[j]);
        }
      }
      free_new_vvector(RF_runningHoldoutCount, 1, RF_xSize, NRUTIL_UPTR);
      for (j = 1; j <= RF_xSize; j++) {
        if (RF_holdBLKptr[j] > 0) {
          for (k = 1; k <= RF_holdBLKptr[j]; k++) {
            free_uivector(RF_blockSerialTreeIndex[j][k], 1, RF_vtryBlockSize);
          }
          free_new_vvector(RF_blockSerialTreeIndex[j], 1, RF_holdBLKptr[j], NRUTIL_UPTR);
        }
      }
      free_new_vvector(RF_blockSerialTreeIndex, 1, RF_xSize, NRUTIL_UPTR2);      
      for (j = 1; j <= RF_xSize; j++) {
        free_uivector(RF_holdoutMap[j], 1, RF_ntree);
      }
      free_new_vvector(RF_holdoutMap, 1, RF_xSize, NRUTIL_UPTR);
  }
  if (RF_opt & OPT_PROX) {
    localSize = ((obsSize + 1)  * obsSize) >> 1;
    free_dvector(RF_proximityDen, 1, localSize);
    free_new_vvector(RF_proximityPtr, 1, obsSize, NRUTIL_DPTR);
    free_new_vvector(RF_proximityDenPtr, 1, obsSize, NRUTIL_DPTR);
  }
  if (RF_optHigh & OPT_DIST) {
    localSize = ((obsSize + 1)  * obsSize) >> 1;
    free_dvector(RF_distanceDen, 1, localSize);
    free_new_vvector(RF_distancePtr, 1, obsSize, NRUTIL_DPTR);
    free_new_vvector(RF_distanceDenPtr, 1, obsSize, NRUTIL_DPTR);
  }
  if (RF_optHigh & OPT_WGHT) {
    free_uivector(RF_weightDenom, 1, obsSize);
  }
  if (RF_opt & OPT_MISS_OUT) {
    if (rspSize > 0) {
      free_new_vvector(RF_sImputeResponsePtr, 1, rspSize, NRUTIL_DPTR);
    }
    free_new_vvector(RF_sImputePredictorPtr, 1, RF_xSize, NRUTIL_DPTR);
  }
  if (RF_opt & (OPT_VARUSED_F | OPT_VARUSED_T)) {
    if (RF_opt & OPT_VARUSED_T) {
    }
    else {
      free_uimatrix(RF_varUsedPtr, 1, RF_ntree, 1, RF_xSize);
    }
  }
  if (RF_opt & (OPT_SPLDPTH_1 | OPT_SPLDPTH_2)) {
  }
  if (RF_optHigh & OPT_MEMB_PRUN) {
  }
  if (RF_opt & OPT_CASE_DPTH) {
  }
  if (RF_optHigh & OPT_PART_PLOT) {
    free_new_vvector(RF_partMembership, 1, RF_ntree, NRUTIL_TPTR2);
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      {
      }
    }
    else {
      if (RF_rTargetFactorCount > 0) {
      }
      if (RF_rTargetNonFactorCount > 0) {
      }
    }
  }
}
void stackForestObjectsPtrOnly(char mode) {
  uint i;
  if (RF_opt & OPT_TREE) {
    RF_totalNodeCount = 0;
    RF_totalMWCPCount = 0;
    if (mode == RF_GROW) {
      RF_treeID_ptr     = (uint **) new_vvector(1, RF_ntree, NRUTIL_UPTR);
      RF_nodeID_ptr     = (uint **) new_vvector(1, RF_ntree, NRUTIL_UPTR);
      RF_nodeSZ_ptr     = (uint **) new_vvector(1, RF_ntree, NRUTIL_UPTR);
      RF_blnodeID_ptr   = (uint **) new_vvector(1, RF_ntree, NRUTIL_UPTR);
      RF_brnodeID_ptr   = (uint **) new_vvector(1, RF_ntree, NRUTIL_UPTR);
      RF_parmID_ptr     = (int ***)    new_vvector(1, RF_ntree, NRUTIL_IPTR2);
      RF_contPT_ptr     = (double ***) new_vvector(1, RF_ntree, NRUTIL_DPTR2);
      RF_mwcpSZ_ptr     = (uint ***)   new_vvector(1, RF_ntree, NRUTIL_UPTR2);
      RF_fsrecID_ptr    = (uint ***)   new_vvector(1, RF_ntree, NRUTIL_UPTR2);
      RF_mwcpCT_ptr     = (uint **)    new_vvector(1, RF_ntree, NRUTIL_UPTR2);
      RF_mwcpPT_ptr     = (uint ***)   new_vvector(1, RF_ntree, NRUTIL_UPTR2);
      for (i = 1; i <= RF_ntree; i++) {
        RF_parmID_ptr[i]  = (int **)    new_vvector(1, 1, NRUTIL_IPTR);
        RF_contPT_ptr[i]  = (double **) new_vvector(1, 1, NRUTIL_DPTR);
        RF_mwcpSZ_ptr[i]  = (uint **)   new_vvector(1, 1, NRUTIL_UPTR);
        RF_fsrecID_ptr[i] = (uint **)   new_vvector(1, 1, NRUTIL_UPTR);
        RF_mwcpCT_ptr[i]  =  uivector(1, 1);
        RF_mwcpPT_ptr[i]  = (uint **)   new_vvector(1, 1, NRUTIL_UPTR);
      }
      for (i = 1; i <= RF_ntree; i++) {
        RF_treeID_ptr[i]     = NULL;
        RF_nodeID_ptr[i]     = NULL;
        RF_nodeSZ_ptr[i]     = NULL;
        RF_blnodeID_ptr[i]   = NULL;
        RF_brnodeID_ptr[i]   = NULL;
        RF_parmID_ptr[i][1]  = NULL;
        RF_contPT_ptr[i][1]  = NULL;
        RF_mwcpSZ_ptr[i][1]  = NULL;
        RF_fsrecID_ptr[i][1] = NULL;
        RF_mwcpCT_ptr[i][1]  = 0;
        RF_mwcpPT_ptr[i][1]  = NULL;
      }
    }  
  }  
}
void unstackForestObjectsPtrOnly(char mode) {
  uint treeID;
  if (RF_opt & OPT_TREE) {
    if (mode == RF_GROW) {
      for (treeID = 1; treeID <= RF_ntree; treeID++) {
        unstackTreeObjectsPtrOnly(treeID);
      }
      free_new_vvector(RF_treeID_ptr,     1, RF_ntree, NRUTIL_UPTR);
      free_new_vvector(RF_nodeID_ptr,     1, RF_ntree, NRUTIL_UPTR);
      free_new_vvector(RF_nodeSZ_ptr,     1, RF_ntree, NRUTIL_UPTR);
      free_new_vvector(RF_blnodeID_ptr,   1, RF_ntree, NRUTIL_UPTR);
      free_new_vvector(RF_brnodeID_ptr,   1, RF_ntree, NRUTIL_UPTR);
      free_new_vvector(RF_parmID_ptr,  1,  RF_ntree, NRUTIL_IPTR2);
      free_new_vvector(RF_contPT_ptr,  1,  RF_ntree, NRUTIL_DPTR2);
      free_new_vvector(RF_mwcpSZ_ptr,  1,  RF_ntree, NRUTIL_UPTR2);
      free_new_vvector(RF_fsrecID_ptr, 1,  RF_ntree, NRUTIL_UPTR2);
      free_new_vvector(RF_mwcpCT_ptr, 1,  RF_ntree, NRUTIL_UPTR2);
      free_new_vvector(RF_mwcpPT_ptr, 1,  RF_ntree, NRUTIL_UPTR2);
    }  
  }  
}
void stackTreeObjectsPtrOnly(char mode, uint treeID) {
  uint treeNodeCount;
  uint mwcpSize;
  uint treeMWCPCount;
  if (RF_opt & OPT_TREE) {
    if (mode == RF_GROW) {
      treeNodeCount = RF_nodeCount[treeID];
      if (RF_xFactorCount > 0) {
        mwcpSize = (RF_xMaxFactorLevel >> (3 + ulog2(sizeof(uint)))) + ((RF_xMaxFactorLevel & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
      }
      else {
        mwcpSize = 0;
      }
      treeMWCPCount = mwcpSize * treeNodeCount;
      RF_treeID_ptr[treeID]     = uivector(1, treeNodeCount);
      RF_nodeID_ptr[treeID]     = uivector(1, treeNodeCount);
      RF_nodeSZ_ptr[treeID]     = uivector(1, treeNodeCount);
      RF_blnodeID_ptr[treeID]   = uivector(1, treeNodeCount);
      RF_brnodeID_ptr[treeID]   = uivector(1, treeNodeCount);
      RF_parmID_ptr[treeID][1]  = ivector(1, treeNodeCount);
      RF_contPT_ptr[treeID][1]  = dvector(1, treeNodeCount);
      RF_mwcpSZ_ptr[treeID][1]  = uivector(1, treeNodeCount);
      RF_fsrecID_ptr[treeID][1] = uivector(1, treeNodeCount);
      if (treeMWCPCount > 0) {
        RF_mwcpPT_ptr[treeID][1]  = uivector(1, treeMWCPCount);
      }
    }  
  }  
}
void unstackTreeObjectsPtrOnly(uint treeID) {
  uint treeNodeCount;
  uint mwcpSize;
  uint treeMWCPCount;
  treeNodeCount = RF_nodeCount[treeID];
  if (RF_xFactorCount > 0) {
    mwcpSize = (RF_xMaxFactorLevel >> (3 + ulog2(sizeof(uint)))) + ((RF_xMaxFactorLevel & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
  }
  else {
    mwcpSize = 0;
  }
  treeMWCPCount = mwcpSize * treeNodeCount;
  free_uivector(RF_treeID_ptr[treeID],     1, treeNodeCount);
  free_uivector(RF_nodeID_ptr[treeID],     1, treeNodeCount);
  free_uivector(RF_nodeSZ_ptr[treeID],     1, treeNodeCount);
  free_uivector(RF_blnodeID_ptr[treeID],   1, treeNodeCount);
  free_uivector(RF_brnodeID_ptr[treeID],   1, treeNodeCount);
  free_ivector(RF_parmID_ptr[treeID][1], 1, treeNodeCount);
  free_dvector(RF_contPT_ptr[treeID][1], 1, treeNodeCount);
  free_uivector(RF_mwcpSZ_ptr[treeID][1], 1, treeNodeCount);
  free_uivector(RF_fsrecID_ptr[treeID][1], 1, treeNodeCount);
  if (treeMWCPCount > 0) {
    free_uivector(RF_mwcpPT_ptr[treeID][1], 1, treeMWCPCount);
  }
  free_new_vvector(RF_parmID_ptr[treeID], 1, 1, NRUTIL_IPTR);
  free_new_vvector(RF_contPT_ptr[treeID], 1, 1, NRUTIL_DPTR);
  free_new_vvector(RF_mwcpSZ_ptr[treeID], 1, 1, NRUTIL_UPTR);
  free_new_vvector(RF_fsrecID_ptr[treeID], 1, 1, NRUTIL_UPTR);
  free_new_vvector(RF_mwcpPT_ptr[treeID], 1, 1, NRUTIL_UPTR);
  free_uivector(RF_mwcpCT_ptr[treeID], 1, 1);
}
void stackForestObjectsOutput(char mode) {
  ulong  totalNodeCount;
  ulong *totalMWCPCount;
  char *resultStr;
  char *adjStr;
  uint asciiLengthOfHexPortion;
  uint i, j;
  if (RF_opt & OPT_TREE) {
    j = 0;
    asciiLengthOfHexPortion = 0;
    while (j > 0) {
      asciiLengthOfHexPortion++;
      j = j >> 4;
    }
    resultStr = cvector(0, RF_SEXP_ASCII_SIZE + asciiLengthOfHexPortion + 1);
    adjStr =    cvector(0, asciiLengthOfHexPortion + 1);
    if (mode == RF_GROW) {
        totalMWCPCount     = ulvector(1, 1);
        RF_parmID_  = (int **)    new_vvector(1, 1, NRUTIL_IPTR);
        RF_contPT_  = (double **) new_vvector(1, 1, NRUTIL_DPTR);
        RF_mwcpSZ_  = (uint **)   new_vvector(1, 1, NRUTIL_UPTR);
        RF_fsrecID_ = (uint **)   new_vvector(1, 1, NRUTIL_UPTR);
        RF_mwcpPT_  = (uint **)   new_vvector(1, 1, NRUTIL_UPTR);
        RF_mwcpCT_  = (uint **)   new_vvector(1, 1, NRUTIL_UPTR);
        totalMWCPCount[1] = 0;
        for (i = 1; i <= RF_ntree; i++) {
          totalMWCPCount[1] += RF_mwcpCT_ptr[i][1];
        }
      RF_totalNodeCount = 0;
      for (i = 1; i <= RF_ntree; i++) {      
        RF_totalNodeCount += RF_nodeCount[i];
      }
      totalNodeCount = RF_totalNodeCount;
      RF_treeID_ = (uint*)   stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_TREE_ID, totalNodeCount, 0, RF_sexpString[RF_TREE_ID], NULL, 1, totalNodeCount);
      RF_nodeID_ = (uint*)   stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_NODE_ID, totalNodeCount, 0, RF_sexpString[RF_NODE_ID], NULL, 1, totalNodeCount);
      RF_nodeSZ_ = (uint*)   stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_NODE_SZ, totalNodeCount, 0, RF_sexpString[RF_NODE_SZ], NULL, 1, totalNodeCount);
      RF_blnodeID_ = (uint*)   stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_BL_NODE_ID, totalNodeCount, 0, RF_sexpString[RF_BL_NODE_ID], NULL, 1, totalNodeCount);
      RF_brnodeID_ = (uint*)   stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_BR_NODE_ID, totalNodeCount, 0, RF_sexpString[RF_BR_NODE_ID], NULL, 1, totalNodeCount);
      RF_treeID_ --;
      RF_nodeID_ --;
      RF_nodeSZ_ --;
      RF_blnodeID_ --;
      RF_brnodeID_ --;
      RF_parmID_[1]  = (int*)    stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_PARM_ID, totalNodeCount, 0, RF_sexpString[RF_PARM_ID], NULL, 1, totalNodeCount);
      RF_contPT_[1]  = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_CONT_PT, totalNodeCount, 0, RF_sexpString[RF_CONT_PT], NULL, 1, totalNodeCount);
      RF_mwcpSZ_[1]  = (uint*)   stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_MWCP_SZ, totalNodeCount, 0, RF_sexpString[RF_MWCP_SZ], NULL, 1, totalNodeCount);
      RF_fsrecID_[1] = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_FS_REC_ID, totalNodeCount, 0, RF_sexpString[RF_FS_REC_ID], NULL, 1, totalNodeCount);
      RF_parmID_[1]  --;
      RF_contPT_[1]  --;
      RF_mwcpSZ_[1]  --;
      RF_fsrecID_[1] --;
      if (totalMWCPCount[1] > 0) {
        RF_mwcpPT_[1] = (uint*)   stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_MWCP_PT, totalMWCPCount[1], 0, RF_sexpString[RF_MWCP_PT], NULL, 1, totalMWCPCount[1]);
        RF_mwcpPT_[1] --;
      }
      else {
        RF_mwcpPT_[1] = (uint*)   stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_MWCP_PT, 1, 0, RF_sexpString[RF_MWCP_PT], NULL, 1, 1);
      }
      RF_mwcpCT_[1] = (uint*)   stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_MWCP_CT, (ulong) RF_ntree, 0, RF_sexpString[RF_MWCP_CT], NULL, 1, RF_ntree);
      RF_mwcpCT_[1] --;
      free_ulvector(totalMWCPCount,     1, 1);
    }  
    free_cvector(resultStr, 0, RF_SEXP_ASCII_SIZE + asciiLengthOfHexPortion + 1);
    free_cvector(adjStr, 0, asciiLengthOfHexPortion + 1);
  }  
}
void writeForestObjectsOutput(char mode) {
  ulong *totalMWCPCount;
  ulong offset;
  uint treeID;
  uint j;
  if (RF_opt & OPT_TREE) {
    if (mode == RF_GROW) {
      offset = 0;
      for (treeID = 1; treeID <= RF_ntree; treeID++) {
        RF_mwcpCT_[1][treeID] = RF_mwcpCT_ptr[treeID][1];
        for (j = 1; j <= RF_nodeCount[treeID]; j++) {
          offset ++;
          RF_treeID_[offset] = RF_treeID_ptr[treeID][j];
          RF_nodeID_[offset] = RF_nodeID_ptr[treeID][j];
          RF_nodeSZ_[offset] = RF_nodeSZ_ptr[treeID][j];
          RF_blnodeID_[offset] = RF_blnodeID_ptr[treeID][j];
          RF_brnodeID_[offset] = RF_brnodeID_ptr[treeID][j];
          RF_parmID_[1][offset]  = RF_parmID_ptr[treeID][1][j];
          RF_contPT_[1][offset]  = RF_contPT_ptr[treeID][1][j];
          RF_mwcpSZ_[1][offset]  = RF_mwcpSZ_ptr[treeID][1][j];
          RF_fsrecID_[1][offset]  = RF_fsrecID_ptr[treeID][1][j];
        }
      }
      if (offset != RF_totalNodeCount) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Inconsistent total node count encountered during writing of forest topology.");
        RF_nativeError("\nRF-SRC:  Offset versus total was:  %20lu, %20lu", offset, RF_totalNodeCount);
        RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
        RF_nativeExit();
      }
      totalMWCPCount = ulvector(1, 1);
      totalMWCPCount[1] = 0;
      for (treeID = 1; treeID <= RF_ntree; treeID++) {
        totalMWCPCount[1] += RF_mwcpCT_[1][treeID];
      }
      offset = 0;
      for (treeID = 1; treeID <= RF_ntree; treeID++) {
        for (j = 1; j <= RF_mwcpCT_[1][treeID]; j++) {
          offset ++;
          RF_mwcpPT_[1][offset] = RF_mwcpPT_ptr[treeID][1][j];
        }
      }
      if (offset != totalMWCPCount[1]) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Inconsistent MWCP count encountered during writing of forest topology.");
        RF_nativeError("\nRF-SRC:  Offset versus total MWCP count was:  %20lu, %20lu", offset, totalMWCPCount[1]);
        RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
        RF_nativeExit();
      }
      free_ulvector(totalMWCPCount, 1, 1);
    }  
  }  
}
void stackForestObjectsAuxOnly(char mode) {
  uint i;
  if (mode != RF_GROW) {
    RF_parmID_   = (int **)    new_vvector(1, 1, NRUTIL_UPTR);
    RF_contPT_   = (double **) new_vvector(1, 1, NRUTIL_DPTR);
    RF_mwcpSZ_   = (uint **)   new_vvector(1, 1, NRUTIL_UPTR);
    RF_fsrecID_  = (uint **)  new_vvector(1, 1, NRUTIL_UPTR);
    RF_mwcpPT_   = (uint **)   new_vvector(1, 1, NRUTIL_UPTR);
    RF_mwcpCT_   = (uint **)   new_vvector(1, 1, NRUTIL_UPTR);
    RF_restoreTreeID = uivector(1, RF_ntree);
    RF_restoreTreeOffset = ulvector(1, RF_ntree);
    for (i = 1; i <= RF_ntree; i++) {
      RF_restoreTreeID[i] = 0;
      RF_restoreTreeOffset[i] = 0;
    }
    RF_restoreMWCPoffset = new_vvector(1, 1, NRUTIL_LPTR);
    RF_restoreMWCPoffset[1] = ulvector(1, RF_ntree);
    for (i = 1; i <= RF_ntree; i++) {
      RF_restoreMWCPoffset[1][i] = 0;
    }
    RF_mwcpCT_[1] = uivector(1, RF_ntree);
    for (i = 1; i <= RF_ntree; i++) {
      RF_mwcpCT_[1][i] = 0;
    }
  }
}
void unstackForestObjectsAuxOnly(char mode) {
  if (mode == RF_GROW) {
    if (RF_opt & OPT_TREE) {
      free_new_vvector(RF_parmID_,  1, 1, NRUTIL_IPTR);
      free_new_vvector(RF_contPT_,  1, 1, NRUTIL_DPTR);
      free_new_vvector(RF_mwcpSZ_,  1, 1, NRUTIL_UPTR);
      free_new_vvector(RF_fsrecID_, 1, 1, NRUTIL_UPTR);
      free_new_vvector(RF_mwcpPT_,  1, 1, NRUTIL_UPTR);
      free_new_vvector(RF_mwcpCT_,  1, 1, NRUTIL_UPTR);
    }
  }
  if (mode != RF_GROW) {
    free_new_vvector(RF_parmID_,  1, 1, NRUTIL_IPTR);
    free_new_vvector(RF_contPT_,  1, 1, NRUTIL_DPTR);
    free_new_vvector(RF_mwcpSZ_,  1, 1, NRUTIL_UPTR);
    free_new_vvector(RF_fsrecID_, 1, 1, NRUTIL_UPTR);
    free_new_vvector(RF_mwcpPT_,  1, 1, NRUTIL_UPTR);
    free_uivector(RF_restoreTreeID, 1, RF_ntree);
    free_ulvector(RF_restoreTreeOffset, 1, RF_ntree);
    free_ulvector(RF_restoreMWCPoffset[1], 1, RF_ntree);
    free_new_vvector(RF_restoreMWCPoffset, 1, 1, NRUTIL_LPTR);
    free_uivector(RF_mwcpCT_[1], 1, RF_ntree);
    free_new_vvector(RF_mwcpCT_, 1, 1, NRUTIL_UPTR);
  }
}
void verifyAndRegisterCustomSplitRules(void) {
  uint familyConstant;
  uint i, j;
  familyConstant = 0;  
  if (RF_splitRule == CUST_SPLIT) {
    RF_splitCustomIdx = (RF_optHigh & OPT_SPLT_CUST) >> 8;
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 16; j++) {
        customFunctionArray[i][j] = NULL;
      }
    }
    registerCustomFunctions();
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      {
        if (!(RF_opt & OPT_COMP_RISK)) {
          familyConstant = SURV_FAM;
        }
        else {
          familyConstant = CRSK_FAM;
        }
        if (customFunctionArray[familyConstant][RF_splitCustomIdx] == NULL) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  Custom split rule not registered:  %10d", RF_splitCustomIdx + 1);
          RF_nativeError("\nRF-SRC:  Please register the rule and recompile the package.");
          RF_nativeExit();
        }
      }
    }
    else {
      if (RF_rTargetFactorCount > 0) {
        if (customFunctionArray[CLAS_FAM][RF_splitCustomIdx] == NULL) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  Custom split rule not registered:  %10d", RF_splitCustomIdx + 1);
          RF_nativeError("\nRF-SRC:  Please register the rule and recompile the package.");
          RF_nativeExit();
        }
      }
      if (RF_rTargetNonFactorCount > 0) {
        if (customFunctionArray[REGR_FAM][RF_splitCustomIdx] == NULL) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  Custom split rule not registered:  %10d", RF_splitCustomIdx + 1);
          RF_nativeError("\nRF-SRC:  Please register the rule and recompile the package.");
          RF_nativeExit();
        }
      }
    }
  }
}
void stackAuxiliaryInfoList(SNPAuxiliaryInfo ***list, uint count) {
   *list = new_vvector(0, count, NRUTIL_XPTR);
   for (uint i = 0; i <= count; i++) {
     (*list)[i] = NULL;
  }
 }
void allocateAuxiliaryInfo(char   targetFlag,
                           char   type,
                           char  *stringIdentifier,
                           SNPAuxiliaryInfo **list,
                           uint   slot,
                           void  *snpPtr,
                           void  *auxiliaryArrayPtr,
                           uint   dimSize,
                           int   *dim) {
  uint dim1, dim2, dim3, dim4;
  ulong offset;
  uint stringLength;
  SNPAuxiliaryInfo *auxInfoPtr = (SNPAuxiliaryInfo*) gblock((size_t) sizeof(SNPAuxiliaryInfo));
  list[slot] = auxInfoPtr;
  auxInfoPtr -> slot = slot;
  auxInfoPtr -> type = type;
  stringLength = strlen(stringIdentifier) + 1;
  auxInfoPtr -> identity = cvector(1, stringLength);
  strcpy(auxInfoPtr -> identity, stringIdentifier);
  auxInfoPtr -> snpPtr = snpPtr;
  auxInfoPtr -> auxiliaryArrayPtr = auxiliaryArrayPtr;
  auxInfoPtr -> dimSize = dimSize;
  (auxInfoPtr -> dim) = ivector(1, dimSize);
  for (uint i = 1; i <= dimSize; i++) {
    (auxInfoPtr -> dim)[i] = dim[i];
  }
  switch(type) {
  case NATIVE_TYPE_NUMERIC:
    if (auxiliaryArrayPtr == NULL) {
    }
    else if (dimSize == 4) {
      offset = 0;
      dim1 = getAuxDim(targetFlag, dim, 0 , 1);
      *((double *****) auxiliaryArrayPtr) = (double ****) new_vvector(1, dim1, NRUTIL_DPTR3);
      for (uint i = 1; i <= dim1; i++) {
        dim2 = getAuxDim(targetFlag, dim, i , 2);
        if (dim2 > 0) {
          (*((double *****) auxiliaryArrayPtr))[i] = (double ***) new_vvector(1, dim2, NRUTIL_DPTR2);
          for (uint j = 1; j <= dim2; j++) {
            dim3 = getAuxDim(targetFlag, dim, j , 3);
            (*((double *****) auxiliaryArrayPtr))[i][j] = (double **) new_vvector(1, dim3, NRUTIL_DPTR);
            for (uint k = 1; k <= dim3; k++) {
              dim4 = getAuxDim(targetFlag, dim, k , 4);
              (*((double *****) auxiliaryArrayPtr))[i][j][k] = (double *) snpPtr + offset - 1;
              offset += dim4;
            }
          }
        }
      }
    }
    else if (dimSize == 3) {
      offset = 0;      
      dim1 = getAuxDim(targetFlag, dim, 0 , 1);
      *((double ****) auxiliaryArrayPtr) = (double ***) new_vvector(1, dim1, NRUTIL_DPTR2);
      for (uint i = 1; i <= dim1; i++) {
        dim2 = getAuxDim(targetFlag, dim, i , 2);
        if (dim2 > 0) {
          (*((double ****) auxiliaryArrayPtr))[i] = (double **) new_vvector(1, dim2, NRUTIL_DPTR);
          for (uint j = 1; j <= dim2; j++) {
            dim3 = getAuxDim(targetFlag, dim, j , 3);
            (*((double ****) auxiliaryArrayPtr))[i][j] = (double *) snpPtr + offset - 1;
            offset += dim3;
          }
        }
      }
    }
    else if (dimSize == 2) {
      offset = 0;
      dim1 = getAuxDim(targetFlag, dim, 0 , 1);
      *((double ***) auxiliaryArrayPtr) = (double **) new_vvector(1, dim1, NRUTIL_DPTR);
      for (uint i = 1; i <= dim1; i++) {
        dim2 = getAuxDim(targetFlag, dim, i , 2);
        (*((double ***) auxiliaryArrayPtr))[i] = (double *) snpPtr + offset - 1;
          offset += dim2;
      }
    }
    else if (dimSize == 1) {
      *((double **) auxiliaryArrayPtr) = (double *) snpPtr - 1;
    }
    else {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Invalid ( > 4 ) dimension in stackAndProtect() auxiliary arrays:  %10d", dimSize);
      RF_nativeExit();
    }
    break;
  case NATIVE_TYPE_INTEGER:
    if (auxiliaryArrayPtr == NULL) {
    }
    else if (dimSize == 4) {
      offset = 0;
      dim1 = getAuxDim(targetFlag, dim, 0 , 1);
      *((uint *****) auxiliaryArrayPtr) = (uint ****) new_vvector(1, dim1, NRUTIL_UPTR3);
      for (uint i = 1; i <= dim1; i++) {
        dim2 = getAuxDim(targetFlag, dim, i , 2);
        (*((uint *****) auxiliaryArrayPtr))[i] = (uint ***) new_vvector(1, dim2, NRUTIL_UPTR2);
        for (uint j = 1; j <= dim2; j++) {
            dim3 = getAuxDim(targetFlag, dim, j , 3);
          (*((uint *****) auxiliaryArrayPtr))[i][j] = (uint **) new_vvector(1, dim3, NRUTIL_UPTR);
          for (uint k = 1; k <= dim3; k++) {
            dim4 = getAuxDim(targetFlag, dim, k , 4);
            (*((uint *****) auxiliaryArrayPtr))[i][j][k] = (uint *) snpPtr + offset - 1;
            offset += dim4;
          }
        }
      }
    }
    else if (dimSize == 3) {
      offset = 0;
      dim1 = getAuxDim(targetFlag, dim, 0 , 1);
      *((uint ****) auxiliaryArrayPtr) = (uint ***) new_vvector(1, dim1, NRUTIL_UPTR2);
      for (uint i = 1; i <= dim1; i++) {
        dim2 = getAuxDim(targetFlag, dim, i , 2);
        (*((uint ****) auxiliaryArrayPtr))[i] = (uint **) new_vvector(1, dim2, NRUTIL_UPTR);
        for (uint j = 1; j <= dim2; j++) {
          dim3 = getAuxDim(targetFlag, dim, j , 3);
          (*((uint ****) auxiliaryArrayPtr))[i][j] = (uint *) snpPtr + offset - 1;
            offset += dim3;
        }
      }
    }
    else if (dimSize == 2) {
      offset = 0;
      dim1 = getAuxDim(targetFlag, dim, 0 , 1);
      *((uint ***) auxiliaryArrayPtr) = (uint **) new_vvector(1, dim1, NRUTIL_UPTR);
      for (uint i = 1; i <= dim1; i++) {
        dim2 = getAuxDim(targetFlag, dim, i , 2);
        (*((uint ***) auxiliaryArrayPtr))[i] = (uint *) snpPtr + offset - 1;
          offset += dim2;
      }
    }
    else if (dimSize == 1) {
      *((uint **) auxiliaryArrayPtr) = (uint *) snpPtr - 1;
    }
    else {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Invalid ( > 4 ) dimension in stackAndProtect() auxiliary arrays:  %10d", dimSize);
      RF_nativeExit();
    }
    break;
  }
}
uint getAuxDim(char flag, int *dim, uint iterIndex, uint slot) {
  uint result = 0;
  if (slot == 1) {
    result = dim[slot];
  }
  else if (dim[slot] >= 1) {
    result = dim[slot];
  }
  else if (dim[slot] == 0) {
    if (flag) {
      result = RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[iterIndex]]];
    }
    else {
      result = RF_rFactorSize[iterIndex];
    }
  }
  else if (dim[slot] == -1) {
    if (flag) {
      result = 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[iterIndex]]];
    }
    else {
      result = 1 + RF_rFactorSize[iterIndex];
    }
  }
  else if (dim[slot] == -2) {
    result = RF_tLeafCount[iterIndex];
  }
  else if (dim[slot] == -3) {
    result = RF_holdBLKptr[iterIndex];
  }
  else {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Inconsistent internal dimension of auxiliary array in getAuxDim():  %10d", dim[slot]);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
  }
  return result;
}
void unstackAuxiliaryInfoAndList(char targetFlag, SNPAuxiliaryInfo **list, uint count) {
  SNPAuxiliaryInfo *auxInfoPtr;
  int  *dim;
  uint  dimSize;
  uint  dim1, dim2, dim3;
  uint stringLength;
  for (uint ii = 0; ii < count; ii++) {
     auxInfoPtr = list[ii];
     if (auxInfoPtr != NULL) {
       dim = auxInfoPtr -> dim;
       dimSize = auxInfoPtr -> dimSize;
       stringLength = strlen(auxInfoPtr -> identity) + 1;
       free_cvector(auxInfoPtr -> identity, 1, stringLength);
       switch(auxInfoPtr -> type) {
       case NATIVE_TYPE_NUMERIC:
         if ((auxInfoPtr -> auxiliaryArrayPtr) == NULL) {
         }
         else if (dimSize == 4) {
           dim1 = getAuxDim(targetFlag, dim, 0 , 1);
           for (uint i = 1; i <= dim1; i++) {
             dim2 = getAuxDim(targetFlag, dim, i , 2);
             if (dim2 > 0) {
               for (uint j = 1; j <= dim2; j++) {
                 dim3 = getAuxDim(targetFlag, dim, j , 3);
                 free_new_vvector((*((double *****) (auxInfoPtr -> auxiliaryArrayPtr)))[i][j], 1, dim3, NRUTIL_DPTR);
               }
               free_new_vvector((*((double *****) (auxInfoPtr -> auxiliaryArrayPtr)))[i], 1, dim2, NRUTIL_DPTR2);
             }
           }
           free_new_vvector((*((double *****) (auxInfoPtr -> auxiliaryArrayPtr))), 1, dim1, NRUTIL_DPTR3);
         }
         else if (dimSize == 3) {
           dim1 = getAuxDim(targetFlag, dim, 0 , 1);
           for (uint i = 1; i <= dim1; i++) {
             dim2 = getAuxDim(targetFlag, dim, i , 2);
             if (dim2 > 0) {
               free_new_vvector((*((double ****) (auxInfoPtr -> auxiliaryArrayPtr)))[i], 1, dim2, NRUTIL_DPTR);
             }
           }
           free_new_vvector((*((double ****) (auxInfoPtr -> auxiliaryArrayPtr))), 1, dim1, NRUTIL_DPTR2);
         }
         else if (dimSize == 2) {
           dim1 = getAuxDim(targetFlag, dim, 0 , 1);
           free_new_vvector((*((double ***) (auxInfoPtr -> auxiliaryArrayPtr))), 1, dim1, NRUTIL_DPTR);
         }
         else if (dimSize == 1) {
         }
         break;
       case NATIVE_TYPE_INTEGER:
         if ((auxInfoPtr -> auxiliaryArrayPtr) == NULL) {
         }
         else if (dimSize == 4) {
           dim1 = getAuxDim(targetFlag, dim, 0 , 1);
           for (uint i = 1; i <= dim1; i++) {
             dim2 = getAuxDim(targetFlag, dim, i , 2);
             for (uint j = 1; j <= dim2; j++) {
               dim3 = getAuxDim(targetFlag, dim, j , 3);
               free_new_vvector((*((uint *****) (auxInfoPtr -> auxiliaryArrayPtr)))[i][j], 1, dim3, NRUTIL_UPTR);
             }
             free_new_vvector((*((uint *****) (auxInfoPtr -> auxiliaryArrayPtr)))[i], 1, dim2, NRUTIL_UPTR2);
           }
           free_new_vvector((*((uint *****) (auxInfoPtr -> auxiliaryArrayPtr))), 1,  dim1, NRUTIL_UPTR3);
         }
         else if (dimSize == 3) {
           dim1 = getAuxDim(targetFlag, dim, 0 , 1);
           for (uint i = 1; i <= dim1; i++) {
             dim2 = getAuxDim(targetFlag, dim, i , 2);             
             free_new_vvector((*((uint ****) (auxInfoPtr -> auxiliaryArrayPtr)))[i], 1, dim2, NRUTIL_UPTR);
           }
           free_new_vvector((*((uint ****) (auxInfoPtr -> auxiliaryArrayPtr))), 1, dim1, NRUTIL_UPTR2);
         }
         else if (dimSize == 2) {
           dim1 = getAuxDim(targetFlag, dim, 0 , 1);
           free_new_vvector((*((uint ***) (auxInfoPtr -> auxiliaryArrayPtr))), 1, dim1, NRUTIL_UPTR);
         }
         else if (dimSize == 1) {
         }
         break;
       }
       free_ivector(auxInfoPtr -> dim, 1, auxInfoPtr -> dimSize);
       free_gblock(auxInfoPtr, sizeof(SNPAuxiliaryInfo));
     }
   }
  free_new_vvector(list, 0, count, NRUTIL_XPTR);
}
