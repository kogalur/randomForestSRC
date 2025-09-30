
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "rfsrcUtil.h"
#include "survival.h"
#include "survivalE.h"
#include "classification.h"
#include "regression.h"
#include "impute.h"
#include "termOps.h"
#include "nrutil.h"
void updateTerminalNodeOutcomes(char       mode,
                                uint       treeID,
                                Terminal  *parent,
                                uint      *repMembrIndx,
                                uint       repMembrSize,
                                uint      *genMembrIndx,
                                uint       genMembrSize,
                                uint      *rmbrIterator,
                                uint      *ambrIterator) {
  uint clasIterator, regrIterator;
  uint i;
  if (!(RF_optHigh & OPT_JIT_TOP)) {
    if (RF_optHigh & OPT_MEMB_INCG) {
      for (i = 1; i <= genMembrSize; i++) {
        ++(*ambrIterator);
        RF_tTermMembership[treeID][RF_AMBR_ID_ptr[treeID][(*ambrIterator)]] = parent;
      }
    }
    else if (RF_optHigh & OPT_MEMB_OUTG) {
      for (i = 1; i <= genMembrSize; i++) {
        RF_tTermMembership[treeID][genMembrIndx[i]] = parent;
        RF_AMBR_ID_ptr[treeID][++(*ambrIterator)] = genMembrIndx[i];
      }
    }
    else {
      for (i = 1; i <= genMembrSize; i++) {
        RF_tTermMembership[treeID][genMembrIndx[i]] = parent;
      }
    }
  }  
  if ((RF_opt & OPT_PERF) ||
      (RF_opt & OPT_OENS) ||
      (RF_opt & OPT_FENS)) {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      getAtRiskAndEventCount(treeID, parent, repMembrIndx, repMembrSize, genMembrIndx, genMembrSize, rmbrIterator);
      if (!(RF_optHigh & OPT_TERM_INCG)) {
        getLocalRatio(treeID, parent);
        getLocalSurvival(treeID, parent);
        if (!(RF_opt & OPT_COMP_RISK)) {
          getLocalNelsonAalen(treeID, parent);
        }
        else {
          getLocalCSH(treeID, parent);
          getLocalCIF(treeID, parent);
        }
        unstackAtRiskAndEventCount(parent);
      }  
      if (!(RF_opt & OPT_COMP_RISK)) {
        getSurvival(treeID, parent);
        getNelsonAalen(treeID, parent);
      }
      else {
        getCSH(treeID, parent);
        getCIF(treeID, parent);
      }
      getMortality(treeID, parent);
      freeTerminalNodeLocalSurvivalStructures(parent);
    }
    else {
      clasIterator = regrIterator = *rmbrIterator;
      if (RF_rFactorCount > 0) {
        getMultiClassProb(treeID, parent, repMembrIndx, repMembrSize, genMembrIndx, genMembrSize, & clasIterator);
        *rmbrIterator = clasIterator;
      }
      if (RF_rNonFactorCount > 0) {
        getMeanResponse(treeID, parent, repMembrIndx, repMembrSize, genMembrIndx, genMembrSize, & regrIterator);
        *rmbrIterator = regrIterator;
      }
    }
  }
  else {
    getMembrCountOnly(treeID, parent, repMembrIndx, repMembrSize, genMembrIndx, genMembrSize);
  }
}
void getMembrCountOnly (uint       treeID,
                        Terminal  *parent,
                        uint      *repMembrIndx,
                        uint       repMembrSize,
                        uint      *genMembrIndx,
                        uint       genMembrSize) {
  if ( !(RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2) ) {
    parent -> membrCount = genMembrSize;
  }
  else {
    parent -> membrCount = repMembrSize;
    if (RF_optHigh & OPT_MEMB_OUTG) {
    }
    if (RF_optHigh & OPT_MEMB_INCG) {
      parent -> membrCount = RF_TN_RCNT_ptr[treeID][parent -> nodeID];
    }
  }
  if ((parent -> membrCount) == 0) {
    if (!(RF_opt & OPT_OUTC_TYPE)) {
    }
  }
}
void updateEnsemble (char mode, uint b) {
  LeafLinkedObj *leafLinkedPtr;
  char      potentiallyMixedMultivariate;
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      updateEnsembleSurvival(mode, b, FALSE);
    }
    else {
      potentiallyMixedMultivariate = FALSE;
      if (RF_rTargetFactorCount > 0) {
        updateEnsembleMultiClass(mode, b, FALSE, potentiallyMixedMultivariate);
        potentiallyMixedMultivariate = TRUE;
      }
      if (RF_rTargetNonFactorCount > 0) {
        updateEnsembleMean(mode, b, FALSE, potentiallyMixedMultivariate);
        if (RF_opt & OPT_QUANTLE) {        
          updateQuantileStream(mode, b);
        }
        potentiallyMixedMultivariate = TRUE;
      }
    }
    switch (mode) {
    case RF_GROW:
      if (!(RF_optHigh & OPT_TERM_OUTG)) {
        if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
          leafLinkedPtr = RF_leafLinkedObjHead[b] -> fwdLink;
          while (leafLinkedPtr != NULL) {
            freeTerminalNodeSurvivalStructuresIntermediate(leafLinkedPtr -> termPtr);
            leafLinkedPtr = leafLinkedPtr -> fwdLink;
          }
          if (!(RF_opt & OPT_VIMP)) {
            leafLinkedPtr = RF_leafLinkedObjHead[b] -> fwdLink;
            while (leafLinkedPtr != NULL) {
              freeTerminalNodeSurvivalStructuresFinal(leafLinkedPtr -> termPtr);
              leafLinkedPtr = leafLinkedPtr -> fwdLink;
            }
          }
        }
        else {
          if (!(RF_opt & OPT_VIMP)) {
            leafLinkedPtr = RF_leafLinkedObjHead[b] -> fwdLink;
            while (leafLinkedPtr != NULL) {
              freeTerminalNodeNonSurvivalStructures(leafLinkedPtr -> termPtr);
              leafLinkedPtr = leafLinkedPtr -> fwdLink;
            }
          }
        }
      }
      break;
    default:
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        if (!(RF_optHigh & OPT_PART_PLOT)) {
          leafLinkedPtr = RF_leafLinkedObjHead[b] -> fwdLink;
          while (leafLinkedPtr != NULL) {
            freeTerminalNodeSurvivalStructuresIntermediate(leafLinkedPtr -> termPtr);
            leafLinkedPtr = leafLinkedPtr -> fwdLink;
          }
        }
        if (!(RF_opt & OPT_VIMP) && !(RF_optHigh & OPT_PART_PLOT)) {
          leafLinkedPtr = RF_leafLinkedObjHead[b] -> fwdLink;
          while (leafLinkedPtr != NULL) {
            freeTerminalNodeSurvivalStructuresFinal(leafLinkedPtr -> termPtr);
            leafLinkedPtr = leafLinkedPtr -> fwdLink;
          }
        }
      }
      else {
        if (!(RF_opt & OPT_VIMP) && !(RF_optHigh & OPT_PART_PLOT)) {
          leafLinkedPtr = RF_leafLinkedObjHead[b] -> fwdLink;
          while (leafLinkedPtr != NULL) {
            freeTerminalNodeNonSurvivalStructures(leafLinkedPtr -> termPtr);
            leafLinkedPtr = leafLinkedPtr -> fwdLink;
          }
        }
      }
      break;
    }
}
void summarizeFaithfulBlockPerformance (char        mode,
                                        uint        b,
                                        uint        blockID,
                                        double    **blkEnsembleMRTnum,
                                        double   ***blkEnsembleCLSnum,
                                        double    **blkEnsembleRGRnum,
                                        double     *blkEnsembleDen,
                                        double    **responsePtr,
                                        double    **perfMRTblk,
                                        double   ***perfCLSblk,
                                        double    **perfRGRblk) {
  uint      obsSize;
  obsSize = (mode == RF_PRED) ?  RF_fobservationSize : RF_observationSize;
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    getPerformance(b,
                   mode,
                   obsSize,
                   responsePtr,
                   blkEnsembleDen,
                   blkEnsembleMRTnum,
                   NULL,
                   NULL,
                   perfMRTblk[blockID],
                   NULL,
                   NULL);
  }
  else {
    if (RF_rTargetFactorCount > 0) {
      getPerformance(b,
                     mode,
                     obsSize,
                     responsePtr,
                     blkEnsembleDen,
                     NULL,
                     blkEnsembleCLSnum,
                     NULL,
                     NULL,
                     perfCLSblk[blockID],
                     NULL);
    }
    if (RF_rTargetNonFactorCount > 0) {
      getPerformance(b,
                     mode,
                     obsSize,
                     responsePtr,
                     blkEnsembleDen,
                     NULL,
                     NULL,
                     blkEnsembleRGRnum,
                     NULL,
                     NULL,
                     perfRGRblk[blockID]);
    }
  }
}
void summarizeHoldoutBlockPerformance (char        mode,
                                       uint        b,
                                       uint        xVarIdx,
                                       uint        blockID,
                                       double    **responsePtr,
                                       double    **holdMRTstd,
                                       double   ***holdCLSstd,
                                       double    **holdRGRstd,
                                       double     *holdEnsembleDen,
                                       double     *holdMRTptr,
                                       double    **holdCLSptr,
                                       double     *holdRGRptr) {
  uint      obsSize;
  obsSize = RF_observationSize;
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    getPerformance(b,
                   mode,
                   obsSize,
                   responsePtr,
                   holdEnsembleDen,
                   holdMRTstd,
                   NULL,
                   NULL,
                   holdMRTptr,
                   NULL,
                   NULL);
  }
  else {
    if (RF_rTargetFactorCount > 0) {
      if (holdCLSstd != NULL) {
        getPerformance(b,
                       mode,
                       obsSize,
                       responsePtr,
                       holdEnsembleDen,
                       NULL,
                       holdCLSstd,
                       NULL,
                       NULL,
                       holdCLSptr,
                       NULL);
      }
    }
    if (RF_rTargetNonFactorCount > 0) {
      if (holdRGRstd != NULL) {
        getPerformance(b,
                       mode,
                       obsSize,
                       responsePtr,
                       holdEnsembleDen,
                       NULL,
                       NULL,
                       holdRGRstd,
                       NULL,
                       NULL,
                       holdRGRptr);
      }
    }
  }
}
char stackAndImputePerfResponse(char      mode,
                                char      multImpFlag,
                                uint      b,
                                uint      loSerialTreeID,
                                uint      hiSerialTreeID,
                                uint     *serialTreePtr,
                                double ***responsePtr) {
  double **mResponsePtr;
  uint     obsSize;
  char     imputeFlag;
  uint     i, p;
  imputeFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    if (b == 0) {
      *responsePtr = RF_fresponseIn;
    }
    else {
      *responsePtr = RF_fresponse[b];
    }
    if (RF_fmRecordSize > 0) {
      if(RF_fmResponseFlag == TRUE) {
        imputeFlag = TRUE;
      }
    }
    break;
  default:
    obsSize  = RF_observationSize;
    if (b == 0) {
      *responsePtr = RF_responseIn;
    }
    else {
      *responsePtr = RF_response[b];
    }
    if (multImpFlag == FALSE) {
      if (RF_mRecordSize > 0) {
        if(RF_mResponseFlag == TRUE) {
          imputeFlag = TRUE;
        }
      }
    }
    break;
  }
  if (imputeFlag) {
    mResponsePtr   = dmatrix(1, RF_ySize, 1, obsSize);
    for (p = 1; p <= RF_ySize; p++) {
      for (i = 1; i <= obsSize; i++) {
        mResponsePtr[p][i] = (*responsePtr)[p][i];
      }
    }
    if (b != 0) {
      imputeResponse(mode, loSerialTreeID, hiSerialTreeID, serialTreePtr, mResponsePtr);
      *responsePtr = mResponsePtr;
    }
    else {
      imputeUpdateShadow(mode, mResponsePtr, NULL);
      *responsePtr = mResponsePtr;
    }      
  }
  return imputeFlag;
}
void unstackPerfResponse(char mode, char flag, double **mResponsePtr) {
  uint obsSize;
  if (flag == TRUE) {
    obsSize = (mode == RF_PRED) ?  RF_fobservationSize : RF_observationSize;
    free_dmatrix(mResponsePtr, 1, RF_ySize, 1, obsSize);
  }
}
void getPerformance(uint      serialTreeID,
                    char      mode,
                    uint      obsSize,
                    double  **responsePtr,
                    double    *denomPtr,
                    double   **outcomeMRT,
                    double  ***outcomeCLS,
                    double   **outcomeRGR,
                    double   *perfMRTptr,
                    double  **perfCLSptr,
                    double   *perfRGRptr) {
  uint j, k;
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (!(RF_opt & OPT_COMP_RISK)) {
      perfMRTptr[1] = getConcordanceIndex(1,
                                          obsSize,
                                          responsePtr[RF_timeIndex],
                                          responsePtr[RF_statusIndex],
                                          outcomeMRT[1],
                                          denomPtr);
    }
    else {
      double *cpv = dvector(1, RF_eventTypeSize);
      getCRPerformance(mode,
                       obsSize,
                       responsePtr,
                       outcomeMRT,
                       denomPtr,
                       cpv);
      for (j=1; j <= RF_eventTypeSize; j++) {
        perfMRTptr[j] = cpv[j];
      }
      free_dvector(cpv, 1, RF_eventTypeSize);
    }
  }
  else {
    if (perfCLSptr != NULL) {
      for (j = 1; j <= RF_rTargetFactorCount; j++) {
        if (RF_opt & OPT_PERF_CALB) {
          double *cpv = dvector(1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]);
          perfCLSptr[j][1] = getBrierScore(obsSize,
                                           RF_rTargetFactor[j],                                                            
                                           responsePtr[RF_rTargetFactor[j]],
                                           outcomeCLS[j],
                                           denomPtr,
                                           cpv);
          for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
            perfCLSptr[j][1+k] = cpv[k];
          }
          free_dvector(cpv, 1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]);
        }
        else if ((RF_opt & OPT_PERF_GMN2) && (RF_rFactorMinorityFlag[RF_rFactorMap[RF_rTargetFactor[j]]] == TRUE)) {
          double *maxVote = dvector(1, obsSize);
          getMaxVote(obsSize,
                     RF_rTargetFactor[j],
                     outcomeCLS[j],
                     denomPtr,
                     maxVote);
          perfCLSptr[j][1] = getGMeanIndex(obsSize,
                                           RF_rTargetFactor[j],
                                           responsePtr[RF_rTargetFactor[j]],
                                           denomPtr,
                                           maxVote);
          free_dvector(maxVote, 1, obsSize);
        }
        else {
          double *cpv = dvector(1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]);
          double *maxVote = dvector(1, obsSize);
          getMaxVote(obsSize,
                     RF_rTargetFactor[j],
                     outcomeCLS[j],
                     denomPtr,
                     maxVote);
          perfCLSptr[j][1] = getClassificationIndex(obsSize,
                                                    RF_rTargetFactor[j],
                                                    responsePtr[RF_rTargetFactor[j]],
                                                    denomPtr,
                                                    maxVote);
          getConditionalClassificationIndex(obsSize,
                                            RF_rTargetFactor[j],
                                            responsePtr[RF_rTargetFactor[j]],
                                            outcomeCLS[j],
                                            maxVote,
                                            denomPtr,
                                            cpv);
          for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
            perfCLSptr[j][1+k] = cpv[k];
          }
          free_dvector(cpv, 1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]);
          free_dvector(maxVote, 1, obsSize);
        }
      }
    }
    if (perfRGRptr != NULL) {
      for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
        perfRGRptr[j] = getMeanSquareError(obsSize,
                                           responsePtr[RF_rTargetNonFactor[j]],
                                           outcomeRGR[j],
                                           denomPtr);
      }
    }
  }
}
void normalizeEnsembleEstimates(char mode, char final) {
  char oobFlag, fullFlag;
  uint      obsSize;
  double ***ensembleSRGptr;
  double  **ensembleMRTptr;
  double  **ensembleSRVptr;
  double ***ensembleCIFptr;
  double ***ensembleCLSptr;
  double  **ensembleRGRptr;
  double ***ensembleSRGnum;
  double  **ensembleMRTnum;
  double  **ensembleSRVnum;
  double ***ensembleCIFnum;
  double ***ensembleCLSnum;
  double  **ensembleRGRnum;
  double      ***ensembleQNTptr;
  QuantileObj ***quantileHead;
  uint         **quantileStreamSize;
  double   *ensembleDen;
  uint i, j, k;
  oobFlag = fullFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    oobFlag = FALSE;
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    break;
  default:
    obsSize = RF_observationSize;
    if (RF_opt & OPT_OENS) {
      oobFlag = TRUE;
    }
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    break;
  }
  while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
    if (oobFlag == TRUE) {
      ensembleDen    = RF_oobEnsembleDen;
      ensembleSRGptr = RF_oobEnsembleSRGptr;
      ensembleMRTptr = RF_oobEnsembleMRTptr;
      ensembleSRVptr = RF_oobEnsembleSRVptr;
      ensembleCIFptr = RF_oobEnsembleCIFptr;
      ensembleCLSptr = RF_oobEnsembleCLSptr;
      ensembleRGRptr = RF_oobEnsembleRGRptr;
      ensembleSRGnum = RF_oobEnsembleSRGnum;
      ensembleMRTnum = RF_oobEnsembleMRTnum;
      ensembleSRVnum = RF_oobEnsembleSRVnum;
      ensembleCIFnum = RF_oobEnsembleCIFnum;
      ensembleCLSnum = RF_oobEnsembleCLSnum;
      ensembleRGRnum = RF_oobEnsembleRGRnum;
      ensembleQNTptr     = RF_oobEnsembleQNTptr;
      quantileHead       = RF_oobQuantileHead;
      quantileStreamSize = RF_oobQuantileStreamSize;
    }
    else {
      ensembleDen    = RF_fullEnsembleDen;
      ensembleSRGptr = RF_fullEnsembleSRGptr;
      ensembleMRTptr = RF_fullEnsembleMRTptr;
      ensembleSRVptr = RF_fullEnsembleSRVptr;
      ensembleCIFptr = RF_fullEnsembleCIFptr;
      ensembleCLSptr = RF_fullEnsembleCLSptr;
      ensembleRGRptr = RF_fullEnsembleRGRptr;
      ensembleSRGnum = RF_fullEnsembleSRGnum;
      ensembleMRTnum = RF_fullEnsembleMRTnum;
      ensembleSRVnum = RF_fullEnsembleSRVnum;
      ensembleCIFnum = RF_fullEnsembleCIFnum;
      ensembleCLSnum = RF_fullEnsembleCLSnum;
      ensembleRGRnum = RF_fullEnsembleRGRnum;
      ensembleQNTptr     = RF_fullEnsembleQNTptr;
      quantileHead       = RF_fullQuantileHead;
      quantileStreamSize = RF_fullQuantileStreamSize;
    }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      for (i = 1; i <= obsSize; i++) {
        if (ensembleDen[i] != 0) {
          if (!(RF_opt & OPT_COMP_RISK)) {
            if (final) {
              for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
                ensembleSRGptr[1][k][i] = ensembleSRGnum[1][k][i] / ensembleDen[i];
                ensembleSRVptr[k][i]    = ensembleSRVnum[k][i] / ensembleDen[i];
              }
            }
            ensembleMRTptr[1][i] = ensembleMRTnum[1][i] / ensembleDen[i];
          }
          else {
            for(j = 1; j <= RF_eventTypeSize; j ++) {
              if (final) {              
                for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
                  ensembleSRGptr[j][k][i] = ensembleSRGnum[j][k][i] / ensembleDen[i];
                  ensembleCIFptr[j][k][i] = ensembleCIFnum[j][k][i] / ensembleDen[i];
                }
              }
              ensembleMRTptr[j][i] = ensembleMRTnum[j][i] / ensembleDen[i];
            }
          }
        }
        else {
          if (!(RF_opt & OPT_COMP_RISK)) {
            if (final) {
              for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
                ensembleSRGptr[1][k][i] = RF_nativeNaN;
                ensembleSRVptr[k][i]    = RF_nativeNaN;
              }
            }
            ensembleMRTptr[1][i] = RF_nativeNaN;
          }
          else {
            for(j = 1; j <= RF_eventTypeSize; j ++) {
              if (final) {              
                for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
                  ensembleSRGptr[j][k][i] = RF_nativeNaN;
                  ensembleCIFptr[j][k][i] = RF_nativeNaN;
                }
              }
              ensembleMRTptr[j][i] = RF_nativeNaN;
            }
          }          
        }
      }
    }  
    else {
      if (RF_rTargetFactorCount > 0) {
        for (i = 1; i <= obsSize; i++) {
          if (ensembleDen[i] != 0) {
            for (j = 1; j <= RF_rTargetFactorCount; j++) {
              for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
                ensembleCLSptr[j][k][i] = ensembleCLSnum[j][k][i] / ensembleDen[i];
              }
            }
          }
          else {
            for (j = 1; j <= RF_rTargetFactorCount; j++) {
              for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
                ensembleCLSptr[j][k][i] = RF_nativeNaN;
              }
            }
          }
        }
      }
      if (RF_rTargetNonFactorCount > 0) {      
        for (i = 1; i <= obsSize; i++) {
          if (ensembleDen[i] != 0) {
            for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
              ensembleRGRptr[j][i] = ensembleRGRnum[j][i] / ensembleDen[i];
            }
          }
          else {
            for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
              ensembleRGRptr[j][i] = RF_nativeNaN;
            }
          }
        }
        if (final) {        
          if (RF_opt & OPT_QUANTLE) {
            for (i = 1; i <= obsSize; i++) {
              if (ensembleDen[i] != 0) {
                for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
                  for (k = 1; k <= RF_quantileSize; k++) {
                    ensembleQNTptr[j][k][i] = getApproxQuantile(quantileHead[j][i], RF_quantile[k], quantileStreamSize[j][i]);
                  }
                }
              }
              else {
                for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
                  for (k = 1; k <= RF_quantileSize; k++) {
                    ensembleQNTptr[j][k][i] = RF_nativeNaN;
                  }
                }
              }
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
char getPerfFlag (char mode, uint serialTreeID) {
  char result;
  if (RF_opt & OPT_PERF) {
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  if (result) {
    if (serialTreeID % RF_perfBlock == 0){
    }
    else {
      if (serialTreeID == RF_ntree) {
      }
      else {
        result = FALSE;
      }
    }
  }
  return result;
}
void getVariablesUsed(uint treeID, Node *parent, uint *varUsedVector) {
  SplitInfo *info;
  if (RF_tLeafCount[treeID] > 0) {
    if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
      info = parent -> splitInfo;
      varUsedVector[info -> randomVar[1]] ++;
      getVariablesUsed(treeID, parent ->  left, varUsedVector);
      getVariablesUsed(treeID, parent -> right, varUsedVector);
    }
  }
  return;
}
