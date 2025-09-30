
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "processEnsemble.h"
#include "rfsrcUtil.h"
#include "tree.h"
#include "treeUtil.h"
#include "marginal.h"
#include "importance.h"
#include "partial.h"
#include "factorOps.h"
#include "nrutil.h"
#include "nativeUtil.h"
#include "error.h"
void stackFactorInSitu(uint treeID) {
  uint j;
  if (RF_rFactorCount + RF_xFactorCount > 0) {
    RF_factorList[treeID] = (Factor **) new_vvector(1, RF_maxFactorLevel, NRUTIL_FPTR);
    for (j = 1; j <= RF_maxFactorLevel; j++) {
      RF_factorList[treeID][j] = NULL;
    }
    for (j = 1; j <= RF_xFactorCount; j++) {
      if (RF_factorList[treeID][RF_xFactorSize[j]] == NULL) {
        RF_factorList[treeID][RF_xFactorSize[j]] = makeFactor(RF_xFactorSize[j], FALSE);
      }
    }
    for (j = 1; j <= RF_rFactorCount; j++) {
      if (RF_factorList[treeID][RF_rFactorSize[j]] == NULL) {
        RF_factorList[treeID][RF_rFactorSize[j]] = makeFactor(RF_rFactorSize[j], FALSE);
      }
    }
  }
}
void unstackFactorInSitu(uint treeID) {
  uint j;
  if (RF_rFactorCount + RF_xFactorCount > 0) {
    if (RF_factorList[treeID] != NULL) {
      for (j = 1; j <= RF_maxFactorLevel; j++) {
        if (RF_factorList[treeID][j] != NULL) {
          freeFactor(RF_factorList[treeID][j]);
        }
      }
      free_new_vvector(RF_factorList[treeID], 1, RF_maxFactorLevel, NRUTIL_FPTR);
      RF_factorList[treeID] = NULL;
    }
  }
}
void processEnsembleInSitu(char mode, char multImpFlag, uint b) {
  char perfFlag;
  double **responsePtr;
  char     rImputeFlag;
  uint lo, hi;
    if ((RF_opt & OPT_PERF) ||
        (RF_opt & OPT_OENS) ||
        (RF_opt & OPT_FENS)) {
#ifdef _OPENMP
      omp_set_lock(&RF_lockPerf);
#endif
      RF_serialTreeIndex[++RF_serialTreeID] = b;
      perfFlag = getPerfFlag(mode, RF_serialTreeID);
      if (!perfFlag) {
#ifdef _OPENMP
        omp_unset_lock(&RF_lockPerf);
#endif
      }
#ifdef _OPENMP
      omp_set_lock(&RF_lockEnsbUpdtCount);
#endif
      RF_ensbUpdtCount++;
#ifdef _OPENMP
      omp_unset_lock(&RF_lockEnsbUpdtCount);
#endif
      updateEnsemble(mode, b);
      if (RF_opt & OPT_VIMP) {
        if (RF_opt & OPT_VIMP_JOIN) {
          stackVimpMembership(mode, & RF_vimpMembership[1][b]);
          getVimpMembership(mode, b, RF_vimpMembership[1][b], 0);
          updateEnsembleVimp(mode, b, RF_vimpMembership[1][b], 1);
          unstackVimpMembership(mode, RF_vimpMembership[1][b]);
        }
        else {
          for (uint p = 1; p <= RF_intrPredictorSize; p++) {
            uint pp = RF_intrPredictor[p];            
            stackVimpMembership(mode, & RF_vimpMembership[p][b]);
            getVimpMembership(mode, b, RF_vimpMembership[p][b], pp);
            updateEnsembleVimp(mode, b, RF_vimpMembership[p][b], p);
            unstackVimpMembership(mode, RF_vimpMembership[p][b]);
          }
        }
      }  
#ifdef _OPENMP
      omp_set_lock(&RF_lockEnsbUpdtCount);
#endif
      RF_ensbUpdtCount--;
#ifdef _OPENMP
      omp_unset_lock(&RF_lockEnsbUpdtCount);
#endif
      if (perfFlag) {
        char ensbUpdtCountFlag = FALSE;
        while (!ensbUpdtCountFlag) {
#ifdef _OPENMP
          omp_set_lock(&RF_lockEnsbUpdtCount);
#endif
          if (RF_ensbUpdtCount == 0) {
            ensbUpdtCountFlag = TRUE;
          }
#ifdef _OPENMP
          omp_unset_lock(&RF_lockEnsbUpdtCount);
#endif
        }
        normalizeEnsembleEstimates(mode, FALSE);
        rImputeFlag = stackAndImputePerfResponse(mode,
                                                 multImpFlag,
                                                 b,
                                                 1,
                                                 RF_serialTreeID,
                                                 RF_serialTreeIndex,
                                                 &responsePtr);
        summarizeFaithfulBlockPerformance(mode,
                                          b,
                                          RF_serialTreeID,
                                          (mode == RF_PRED) ? RF_fullEnsembleMRTptr : RF_oobEnsembleMRTptr,
                                          (mode == RF_PRED) ? RF_fullEnsembleCLSptr : RF_oobEnsembleCLSptr,
                                          (mode == RF_PRED) ? RF_fullEnsembleRGRptr : RF_oobEnsembleRGRptr,
                                          (mode == RF_PRED) ? RF_fullEnsembleDen    : RF_oobEnsembleDen,
                                          responsePtr,
                                          RF_perfMRTptr,
                                          RF_perfCLSptr,
                                          RF_perfRGRptr);
        unstackPerfResponse(mode, rImputeFlag, responsePtr);
        if (RF_opt & OPT_VIMP) {
          RF_serialBlockID ++;
          normalizeBlockedEnsembleEstimates(mode,
                                            RF_blkEnsembleMRTnum,
                                            RF_blkEnsembleCLSnum,
                                            RF_blkEnsembleRGRnum,
                                            RF_blkEnsembleDen);
          lo = ((RF_serialBlockID - 1) * RF_perfBlock) + 1;
          hi = RF_serialBlockID * RF_perfBlock;
          if (hi <= RF_ntree) {
            rImputeFlag = stackAndImputePerfResponse(mode, multImpFlag, b, lo, hi, RF_serialTreeIndex, &responsePtr);
            summarizeFaithfulBlockPerformance(mode,
                                              b,
                                              RF_serialBlockID,
                                              RF_blkEnsembleMRTnum,
                                              RF_blkEnsembleCLSnum,
                                              RF_blkEnsembleRGRnum,
                                              RF_blkEnsembleDen,
                                              responsePtr,
                                              RF_perfMRTblk,
                                              RF_perfCLSblk,
                                              RF_perfRGRblk);
            if (RF_opt & OPT_VIMP_JOIN) {
              summarizePerturbedPerformance(mode, b, RF_serialBlockID, 1, responsePtr);
            }
            else {
              for (uint p = 1; p <= RF_intrPredictorSize; p++) {
                summarizePerturbedPerformance(mode, b, RF_serialBlockID, p, responsePtr);
              }
            }
            unstackPerfResponse(mode, rImputeFlag, responsePtr);     
            resetBlockedEnsembleEstimates(mode);
          }
        }
#ifdef _OPENMP
        omp_unset_lock(&RF_lockPerf);
#endif
      }  
    }  
    if (RF_opt & (OPT_SPLDPTH_1 | OPT_SPLDPTH_2)) {
      updateSplitDepth(b, RF_root[b], RF_maxDepth[b]);
    }
    if (RF_opt & OPT_CASE_DPTH) {
      updateCaseDepth(mode, b);
    }
    if (RF_opt & (OPT_VARUSED_F | OPT_VARUSED_T)) {
      getVariablesUsed(b, RF_root[b], RF_varUsedPtr[b]);
    }
    if (RF_optHigh & OPT_PART_PLOT) {
      getAndUpdatePartialMembership(b, RF_root[b]);
    }
      if (RF_xMarginalSize > 0) {
        getMarginalMembership(mode, b);
      }
    if (RF_optHigh & OPT_WGHT) {
      updateWeight(mode, b);
    }
    if (RF_optHigh & OPT_DIST) {
      updateDistance(mode, b);
    }
    if (RF_opt & OPT_PROX) {
      updateProximity(mode, b);
    }
      if (RF_xMarginalSize > 0) {
        releaseMarginalMembership(mode, b);
      }
}
void processEnsemblePost(char mode) {
  char perfFlag;
  double **responsePtr;
  uint lo, hi;
  uint bb;
  if ((RF_opt & OPT_PERF) ||
      (RF_opt & OPT_OENS) ||
      (RF_opt & OPT_FENS)) {
    if (getUserTraceFlag()) {
      RF_userTimeStart = time(NULL);
      RF_userTimeSplit = RF_userTreeID = 0;
      RF_nativePrint("\n");
    }
    for (bb = 1; bb <= RF_getTreeCount; bb++) {
      uint b = RF_getTreeIndex[bb]; 
      if (RF_tLeafCount[b] > 0) {
        RF_serialTreeIndex[++RF_serialTreeID] = b;
        perfFlag = getPerfFlag(mode, RF_serialTreeID);
        updateEnsemble(mode, b);
        if (RF_opt & OPT_VIMP) {
          if (RF_opt & OPT_VIMP_JOIN) {
            stackVimpMembership(mode, & RF_vimpMembership[1][b]);              
            getVimpMembership(mode, b, RF_vimpMembership[1][b], 0);
            updateEnsembleVimp(mode, b, RF_vimpMembership[1][b], 1);
            unstackVimpMembership(mode, RF_vimpMembership[1][b]);
          }
          else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
            for (uint p = 1; p <= RF_intrPredictorSize; p++) {
              uint pp = RF_intrPredictor[p];
              stackVimpMembership(mode, & RF_vimpMembership[p][b]);
              getVimpMembership(mode, b, RF_vimpMembership[p][b], pp);
              updateEnsembleVimp(mode, b, RF_vimpMembership[p][b], p);
              unstackVimpMembership(mode, RF_vimpMembership[p][b]);
            }
          }
        }  
        if (perfFlag) {
          normalizeEnsembleEstimates(mode, FALSE);
          stackAndImputePerfResponse(mode,
                                     FALSE, 
                                     b,
                                     1,
                                     RF_serialTreeID,
                                     RF_serialTreeIndex,
                                     &responsePtr);
          summarizeFaithfulBlockPerformance(mode,
                                            b,
                                            RF_serialTreeID,
                                            (mode == RF_PRED) ? RF_fullEnsembleMRTptr : RF_oobEnsembleMRTptr,
                                            (mode == RF_PRED) ? RF_fullEnsembleCLSptr : RF_oobEnsembleCLSptr,
                                            (mode == RF_PRED) ? RF_fullEnsembleRGRptr : RF_oobEnsembleRGRptr,
                                            (mode == RF_PRED) ? RF_fullEnsembleDen    : RF_oobEnsembleDen,
                                            responsePtr,
                                            RF_perfMRTptr,
                                            RF_perfCLSptr,
                                            RF_perfRGRptr);
          unstackPerfResponse(mode, FALSE, responsePtr);
          if (RF_opt & OPT_VIMP) {
            RF_serialBlockID ++;
            normalizeBlockedEnsembleEstimates(mode,
                                              RF_blkEnsembleMRTnum,
                                              RF_blkEnsembleCLSnum,
                                              RF_blkEnsembleRGRnum,
                                              RF_blkEnsembleDen);
            lo = ((RF_serialBlockID - 1) * RF_perfBlock) + 1;
            hi = RF_serialBlockID * RF_perfBlock;
            if (hi <= RF_ntree) {
              stackAndImputePerfResponse(mode,
                                         FALSE,
                                         b,
                                         lo,
                                         hi,
                                         RF_serialTreeIndex,
                                         &responsePtr);
              summarizeFaithfulBlockPerformance(mode,
                                                b,
                                                RF_serialBlockID,
                                                RF_blkEnsembleMRTnum,
                                                RF_blkEnsembleCLSnum,
                                                RF_blkEnsembleRGRnum,
                                                RF_blkEnsembleDen,
                                                responsePtr,
                                                RF_perfMRTblk,
                                                RF_perfCLSblk,
                                                RF_perfRGRblk);
              if (RF_opt & OPT_VIMP_JOIN) {
                summarizePerturbedPerformance(mode, b, RF_serialBlockID, 1, responsePtr);
              }
              else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
                for (uint p = 1; p <= RF_intrPredictorSize; p++) {
                  summarizePerturbedPerformance(mode, b, RF_serialBlockID, p, responsePtr);
                }
              }
              unstackPerfResponse(mode, FALSE, responsePtr);     
              resetBlockedEnsembleEstimates(mode);
            }  
          }  
        }  
      }  
      if (getUserTraceFlag()) {
#ifdef _OPENMP
#pragma omp critical (_update_timer)
#endif
        { 
          double userTimeElapsedFromStart;
          double userTimeElapsedFromSplit;
          double userTimeRemaining;
          time_t current;
          RF_userTreeID++;
          current = time(NULL);
          userTimeElapsedFromSplit = (double) (current - RF_userTimeSplit);
          if ((userTimeElapsedFromSplit) > (double) getUserTraceFlag()) {
            userTimeElapsedFromStart = (double) (current - RF_userTimeStart);
            userTimeRemaining = (userTimeElapsedFromStart / RF_userTreeID * RF_ntree) - userTimeElapsedFromStart;
            RF_nativePrint("Ensembles Processed:  %6d,    Time Remaining (sec):  %6.0f \n", RF_userTreeID, ceil(userTimeRemaining));
            RF_userTimeSplit = current;
          }
        }  
      }
    }  
  }  
  if (RF_opt & (OPT_SPLDPTH_1 | OPT_SPLDPTH_2)) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
    for (bb = 1; bb <= RF_getTreeCount; bb++) {
      uint b = RF_getTreeIndex[bb];
      if (RF_tLeafCount[b] > 0) {
        updateSplitDepth(b, RF_root[b], RF_maxDepth[b]);
      }
    }
  }
  if (RF_opt & OPT_CASE_DPTH) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
    for (bb = 1; bb <= RF_getTreeCount; bb++) {
      uint b = RF_getTreeIndex[bb];
      if (RF_tLeafCount[b] > 0) {
        updateCaseDepth(mode, b);
      }
    }
  }
  if (RF_opt & (OPT_VARUSED_F | OPT_VARUSED_T)) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
    for (bb = 1; bb <= RF_getTreeCount; bb++) {
      uint b = RF_getTreeIndex[bb];      
      if (RF_tLeafCount[b] > 0) {
        getVariablesUsed(b, RF_root[b], RF_varUsedPtr[b]);
      }
    }
  }
  if (RF_optHigh & OPT_PART_PLOT) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
    for (bb = 1; bb <= RF_getTreeCount; bb++) {
      uint b = RF_getTreeIndex[bb];
      if (RF_tLeafCount[b] > 0) {
        getAndUpdatePartialMembership(b, RF_root[b]);
      }
    }
  }
  for (bb = 1; bb <= RF_getTreeCount; bb++) {
    uint b = RF_getTreeIndex[bb];
    if (RF_tLeafCount[b] > 0) {
        if (RF_xMarginalSize > 0) {
          getMarginalMembership(mode, b);
        }
        if (RF_optHigh & OPT_WGHT) {
          updateWeight(mode, b);
        }
      if (RF_optHigh & OPT_DIST) {
        updateDistance(mode, b);
      }
      if (RF_opt & OPT_PROX) {
        updateProximity(mode, b);
      }
        if (RF_xMarginalSize > 0) {
          releaseMarginalMembership(mode, b);
        }
    }  
  }  
}
void processEnsembleHoldout(uint xVarIdx, uint b) {
  Terminal *terminalNode;
  double   *denomPtr;  
  uint  *membershipIndex;
  uint   membershipSize;
  uint   obsSize;
  uint   blockID;
  uint i, j, k, m, ii;
  obsSize = RF_observationSize;
    blockID = RF_holdoutMap[xVarIdx][b];
    if (blockID > 0) {
      if (blockID <= RF_holdBLKptr[xVarIdx]) {
#ifdef _OPENMP
        omp_set_lock(&(RF_lockVimpHoldout[xVarIdx][blockID]));
#endif
        if (RF_holdEnsembleDen[xVarIdx][blockID] == NULL) {
          RF_holdEnsembleDen[xVarIdx][blockID] = dvector(1, obsSize);
          for (m = 1; m <= obsSize; m++) {
            RF_holdEnsembleDen[xVarIdx][blockID][m] = 0;
          }
        }
        if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
          {
            if (RF_holdMRTstd[xVarIdx][blockID] == NULL) {
              RF_holdMRTstd[xVarIdx][blockID] = (double **) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR);
              for (k = 1; k <= RF_eventTypeSize; k++) {
                RF_holdMRTstd[xVarIdx][blockID][k] = dvector(1, obsSize);
                for (m = 1; m <= obsSize; m++) {
                  RF_holdMRTstd[xVarIdx][blockID][k][m] = 0.0;
                }
              }
            }
          }
        }
        else {
          if (RF_rTargetFactorCount > 0) {
            if (RF_holdCLSstd[xVarIdx][blockID] == NULL) {
              RF_holdCLSstd[xVarIdx][blockID] = (double ***) new_vvector(1, RF_rTargetFactorCount, NRUTIL_DPTR2);
              for (j = 1; j <= RF_rTargetFactorCount; j++) {
                RF_holdCLSstd[xVarIdx][blockID][j]  = (double **) new_vvector(1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]], NRUTIL_DPTR);
                for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
                  RF_holdCLSstd[xVarIdx][blockID][j][k] = dvector(1, obsSize);
                  for (m = 1; m <= obsSize; m++) {
                    RF_holdCLSstd[xVarIdx][blockID][j][k][m] = 0.0;
                  }
                }
              }
            }
          }
          if (RF_rTargetNonFactorCount > 0) {
            if (RF_holdRGRstd[xVarIdx][blockID] == NULL) {
              RF_holdRGRstd[xVarIdx][blockID] = (double **) new_vvector(1, RF_rTargetNonFactorCount, NRUTIL_DPTR);
              for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
                RF_holdRGRstd[xVarIdx][blockID][j] = dvector(1, obsSize);
                for (m = 1; m <= obsSize; m++) {
                  RF_holdRGRstd[xVarIdx][blockID][j][m] = 0.0;
                }
              }
            }
          }
        }
        membershipSize  = RF_oobSize[b];
        membershipIndex = RF_oobMembershipIndex[b];
        denomPtr = RF_holdEnsembleDen[xVarIdx][blockID];
        for (i = 1; i <= membershipSize; i++) {
          ii = membershipIndex[i];
          terminalNode = RF_tTermMembership[b][ii];
          if ((terminalNode -> membrCount) > 0) {
            denomPtr[ii] ++;
            if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
              {
                for (k = 1; k <= RF_eventTypeSize; k++) {
                  RF_holdMRTstd[xVarIdx][blockID][k][ii] += terminalNode -> mortality[k];
                }
              }
            }
            else {
              if (RF_rTargetFactorCount > 0) {
                for (j = 1; j <= RF_rTargetFactorCount; j++) {
                  for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
                    RF_holdCLSstd[xVarIdx][blockID][j][k][ii] += (double) (terminalNode -> multiClassProb)[RF_rFactorMap[RF_rTargetFactor[j]]][k] / (double) (terminalNode -> membrCount);
                  }
                }
              }
              if (RF_rTargetNonFactorCount > 0) {
                for (j=1; j <= RF_rTargetNonFactorCount; j++) {
                  RF_holdRGRstd[xVarIdx][blockID][j][ii] += (terminalNode -> meanResponse)[RF_rNonFactorMap[RF_rTargetNonFactor[j]]];
                }
              }
            }
          }
          else {
              RF_nativePrint("\nRF-SRC:  *** ERROR *** ");
              RF_nativePrint("\nRF-SRC:  NA encountered for VIMP outcome in terminal node:  %10d", terminalNode -> nodeID);
              RF_nativePrint("\nRF-SRC:  Please Contact Technical Support.");
              RF_nativeExit();
          }
        }  
        RF_runningHoldoutCount[xVarIdx][blockID] ++;
        if (RF_runningHoldoutCount[xVarIdx][blockID] == RF_vtryBlockSize) {
          double **responsePtr;
          char     rImputeFlag;
          rImputeFlag = stackAndImputePerfResponse(RF_GROW,
                                                   FALSE, 
                                                   b,
                                                   1,
                                                   RF_vtryBlockSize,
                                                   RF_blockSerialTreeIndex[xVarIdx][blockID],
                                                   &responsePtr);
          if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
            {
              normalizeBlockedEnsembleEstimates(RF_GROW,
                                                RF_holdMRTstd[xVarIdx][blockID],
                                                NULL,
                                                NULL,
                                                RF_holdEnsembleDen[xVarIdx][blockID]);
              summarizeHoldoutBlockPerformance(RF_GROW,
                                               b,
                                               xVarIdx,
                                               blockID,
                                               responsePtr,
                                               RF_holdMRTstd[xVarIdx][blockID],
                                               NULL,
                                               NULL,
                                               RF_holdEnsembleDen[xVarIdx][blockID],
                                               RF_holdMRTptr[xVarIdx][blockID],
                                               NULL,
                                               NULL);
              if (RF_holdMRTstd[xVarIdx][blockID] != NULL) {
                for (k = 1; k <= RF_eventTypeSize; k++) {
                  free_dvector(RF_holdMRTstd[xVarIdx][blockID][k], 1, obsSize);
                }
                free_new_vvector(RF_holdMRTstd[xVarIdx][blockID], 1, RF_eventTypeSize, NRUTIL_DPTR);
                RF_holdMRTstd[xVarIdx][blockID] = NULL;
              }
            }
          }
          else {
            if (RF_rTargetFactorCount > 0) {
              normalizeBlockedEnsembleEstimates(RF_GROW,
                                                NULL,
                                                RF_holdCLSstd[xVarIdx][blockID],
                                                NULL,
                                                RF_holdEnsembleDen[xVarIdx][blockID]);
              summarizeHoldoutBlockPerformance(RF_GROW,
                                               b,
                                               xVarIdx,
                                               blockID,
                                               responsePtr,
                                               NULL,
                                               RF_holdCLSstd[xVarIdx][blockID],
                                               NULL,
                                               RF_holdEnsembleDen[xVarIdx][blockID],
                                               NULL,
                                               RF_holdCLSptr[xVarIdx][blockID],
                                               NULL);
              if (RF_holdCLSstd[xVarIdx][blockID] != NULL) {
                for (j = 1; j <= RF_rTargetFactorCount; j++) {
                  for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
                    free_dvector(RF_holdCLSstd[xVarIdx][blockID][j][k], 1, obsSize);
                  }
                  free_new_vvector(RF_holdCLSstd[xVarIdx][blockID][j], 1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]], NRUTIL_DPTR);
                }
                free_new_vvector(RF_holdCLSstd[xVarIdx][blockID], 1, RF_rTargetFactorCount, NRUTIL_DPTR2);
                RF_holdCLSstd[xVarIdx][blockID] = NULL;
              }
            }
            if (RF_rTargetNonFactorCount > 0) {
              normalizeBlockedEnsembleEstimates(RF_GROW,
                                                NULL,
                                                NULL,
                                                RF_holdRGRstd[xVarIdx][blockID],
                                                RF_holdEnsembleDen[xVarIdx][blockID]);
              summarizeHoldoutBlockPerformance(RF_GROW,
                                               b,
                                               xVarIdx,
                                               blockID,
                                               responsePtr,
                                               NULL,
                                               NULL,
                                               RF_holdRGRstd[xVarIdx][blockID],
                                               RF_holdEnsembleDen[xVarIdx][blockID],
                                               NULL,
                                               NULL,
                                               RF_holdRGRptr[xVarIdx][blockID]);
              if (RF_holdRGRstd[xVarIdx][blockID] != NULL) {
                for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
                  free_dvector(RF_holdRGRstd[xVarIdx][blockID][j], 1, obsSize);
                }
                free_new_vvector(RF_holdRGRstd[xVarIdx][blockID], 1, RF_rTargetNonFactorCount, NRUTIL_DPTR);
                RF_holdRGRstd[xVarIdx][blockID] = NULL;
              }
            }
          }
          if (RF_holdEnsembleDen[xVarIdx][blockID] != NULL) {
            free_dvector(RF_holdEnsembleDen[xVarIdx][blockID], 1, obsSize);
            RF_holdEnsembleDen[xVarIdx][blockID] = NULL;
          }
          unstackPerfResponse(RF_GROW, rImputeFlag, responsePtr);                
        } 
#ifdef _OPENMP
        omp_unset_lock(&(RF_lockVimpHoldout[xVarIdx][blockID]));
#endif
      }  
    }  
}
void processEnsembleHoldoutPost(uint bb) {
  if (RF_tLeafCount[RF_getTreeIndex[bb]] > 0) {
    for (uint p = 1; p <= RF_xSize; p++) {
      processEnsembleHoldout(p, RF_getTreeIndex[bb]);
    }
  }
}
