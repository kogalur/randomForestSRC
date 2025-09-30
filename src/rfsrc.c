
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "rfsrc.h"
#include "stackPreDefined.h"
#include "splitUtil.h"
#include "splitGreedy.h"
#include "splitMult.h"
#include "stackOutput.h"
#include "stackParallel.h"
#include "tree.h"
#include "rfsrcUtil.h"
#include "treeUtil.h"
#include "treeJIT.h"
#include "impute.h"
#include "sexpOutgoing.h"
#include "importance.h"
#include "importanceRand.h"
#include "importanceAnti.h"
#include "importancePerm.h"
#include "partial.h"
#include "splitRegr.h"
#include "splitClas.h"
#include "stack.h"
#include "stackOutputQQ.h"
#include "classification.h"
#include "nodeOps.h"
#include "random.h"
#include "nrutil.h"
#include "nativeUtil.h"
#include "error.h"
void rfsrc(char mode, int seedValue) {
  uint   adj;
  ulong *mwcpOffset;
  uint previousTreeID;
  uint r;
  uint b, p;
  uint seedValueLC;
  seedValueLC    = 0; 
  if (seedValue >= 0) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Parameter verification failed.");
    RF_nativeError("\nRF-SRC:  Random seed must be less than zero.  \n");
    RF_nativeExit();
  }
  if (RF_nImpute < 1) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Parameter verification failed.");
    RF_nativeError("\nRF-SRC:  Number imputations must be greater than zero:  %10d \n", RF_ntree);
    RF_nativeExit();
  }
  if (RF_ntree < 1) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Parameter verification failed.");
    RF_nativeError("\nRF-SRC:  Number of bootstrap iterations must be greater than zero:  %10d \n", RF_ntree);
    RF_nativeExit();
  }
  if (RF_observationSize < 1) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Parameter verification failed.");
    RF_nativeError("\nRF-SRC:  Number of individuals must be greater than one:  %10d \n", RF_observationSize);
    RF_nativeExit();
  }
  if (RF_xSize < 1) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Parameter verification failed.");
    RF_nativeError("\nRF-SRC:  Number of parameters must be greater than zero:  %10d \n", RF_xSize);
    RF_nativeExit();
  }
  if ((RF_perfBlock < 1) || (RF_perfBlock > RF_ntree)) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Parameter verification failed.");
    RF_nativeError("\nRF-SRC:  Invalid value specified for error block count:  %10d \n", RF_perfBlock);
    RF_nativeExit();
  }
#ifdef _OPENMP
  if (RF_numThreads < 0) {
    RF_numThreads = omp_get_max_threads();
  }
  else if (RF_numThreads == 0) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Parameter verification failed.");
    RF_nativeError("\nRF-SRC:  Number of threads must not be zero:  %10d \n", RF_numThreads);
    RF_nativeExit();
  }
  else {
    RF_numThreads = (RF_numThreads < omp_get_max_threads()) ? (RF_numThreads) : (omp_get_max_threads());
  }
#endif
  stackIncomingArrays(mode);
  stackPreDefinedCommonArrays(mode,
                              &RF_nodeMembership,  
                              &RF_tTermMembership,
                              &RF_tTermList,
                              &RF_root);
  switch (mode) {
  case RF_PRED:
    stackPreDefinedPredictArrays();
    break;
  case RF_REST:
    stackPreDefinedRestoreArrays();
    break;
  default:
    stackPreDefinedGrowthArrays();
    break;
  }
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    stackAndInitializeTimeAndSubjectArrays(mode);
  }
  stackFactorArrays(mode);
  initializeFactorArrays(mode);
  stackMissingArraysPhase1(mode);
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      stackCompetingArrays(mode);
  }
  if (RF_rFactorCount > 0) {
    stackClassificationArrays(mode);
  }
  RF_perfBlockCount = (uint) floor(((double) RF_ntree) / RF_perfBlock);
  freeNode = & freeNodeGeneric;
  switch (mode) {
  case RF_GROW:
    stackForestObjectsPtrOnly(mode);
    getPreSplitResult = & getPreSplitResultGeneric;
    stackRandomCovariates   = & stackRandomCovariatesGeneric;
    unstackRandomCovariates = & unstackRandomCovariatesGeneric;
    selectRandomCovariates  = & selectRandomCovariatesGeneric;
    stackRandomResponses   = & stackRandomResponsesGeneric;
    unstackRandomResponses = & unstackRandomResponsesGeneric;
    selectRandomResponses  = & selectRandomResponsesGenericVector;
    virtuallySplitNode = & virtuallySplitNodeGeneric;
    updateMaximumSplit = & updateMaximumSplitGeneric;
    forkAndUpdate = & forkAndUpdateGeneric;
    freeNode = & freeNodeGeneric;
    unstackSplitVector = & unstackSplitVectorGeneric;
    growTree = & growTreeRecursive;
    regressionXwghtSplit = & regressionXwghtSplitCur;
    classificationXwghtSplit = & classificationXwghtSplitCur;
    multivariateSplit        = & multivariateSplitNew;
    unsupervisedSplit        = & unsupervisedSplitMiss;
    randomSplit = & randomSplitGeneric;
    if (RF_mRecordSize == 0) {
      getPreSplitResult = & getPreSplitResultNoMiss;
      multivariateSplit        = & multivariateSplitNew3;
      unsupervisedSplit        = & unsupervisedSplitNew;
    }
    if (RF_yWeightType == RF_WGHT_UNIFORM){
      stackRandomResponses   = & stackRandomResponsesSimple;
      unstackRandomResponses = & unstackRandomResponsesSimple;
      selectRandomResponses  = & selectRandomResponsesSimpleVector;
    }
    RF_optHigh = RF_optHigh & (~OPT_JIT_TOP);
    RF_opt = RF_opt & (~OPT_ANON);
    acquireTree = &acquireTreeGeneric;
    antiMembership = &antiMembershipGeneric;
    randomMembership = &randomMembershipGeneric;
    getMembership = &getMembershipGeneric;
    partialMembership = &partialMembershipGeneric;
    RF_splitRuleObj = makeSplitRuleObj(RF_splitRule);
    break;
  default:
    acquireTree = &acquireTreeGeneric;
    antiMembership = &antiMembershipGeneric;
    randomMembership = &randomMembershipGeneric;
    getMembership = &getMembershipGeneric;
    partialMembership = &partialMembershipGeneric;      
    if (RF_opt & OPT_ANON) {
      if (mode == RF_REST) {
        RF_optHigh = RF_optHigh | (~OPT_JIT_TOP);
      }
    }
    if (RF_optHigh & OPT_JIT_TOP) {
      if ((RF_mRecordSize == 0) &&
          (RF_optHigh & OPT_MEMB_INCG) &&
          (RF_optHigh & OPT_TERM_INCG) &&
          !(RF_optHigh & OPT_MEMB_PRUN) &&
          !(RF_optHigh & OPT_PART_PLOT) &&
          !(RF_xMarginalSize > 0)) {
        acquireTree = &acquireTreeJIT;
        antiMembership = &antiMembershipJIT;
        randomMembership = &randomMembershipJIT;
        getMembership = &getMembershipJIT;
        partialMembership = &partialMembershipJIT;
        if (RF_fmResponseFlag) {
          RF_opt = RF_opt & (~OPT_PERF);
        }
      }
      else {
        RF_optHigh = RF_optHigh & (~OPT_JIT_TOP);
      }
    }
    else {
    }
    if (mode == RF_REST) {
      if (!(RF_optHigh & OPT_JIT_TOP)) {
        if (RF_mRecordSize > 0) {
          RF_optHigh = RF_optHigh & (~OPT_MEMB_INCG);
          RF_optHigh = RF_optHigh & (~OPT_TERM_INCG);
        }
      }
    }
    if (mode == RF_PRED) {
      if (!(RF_optHigh & OPT_JIT_TOP)) {
        if (RF_fmRecordSize > 0) {
          RF_optHigh = RF_optHigh & (~OPT_MEMB_INCG);
          RF_optHigh = RF_optHigh & (~OPT_TERM_INCG);
        }
      }
    }
    adj = 1;
    mwcpOffset = ulvector(1, adj);
    for (uint j = 1; j <= adj; j++) {
      mwcpOffset[j] = 0;
    }
    previousTreeID = b = 0;
    for (ulong ui = 1; ui <= RF_totalNodeCount; ui++) {
      if ((RF_treeID_[ui] > 0) && (RF_treeID_[ui] <= RF_ntree)) {
        if (RF_treeID_[ui] != previousTreeID) {
          previousTreeID = RF_restoreTreeID[++b] = RF_treeID_[ui];
          RF_restoreTreeOffset[RF_treeID_[ui]] = ui;
        }
        RF_nodeCount[RF_treeID_[ui]] ++;
        for (uint j = 1; j <= adj; j++) {
          RF_mwcpCT_[j][RF_treeID_[ui]] += RF_mwcpSZ_[j][ui];
        }
      }
      else {
        RF_nativeError("\nRF-SRC:  Diagnostic Trace of Tree Record:  \n");
        RF_nativeError("\nRF-SRC:      treeID     nodeID ");
        RF_nativeError("\nRF-SRC:  %10d %10d \n", RF_treeID_[ui], RF_nodeID_[ui]);
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Invalid forest input record at line:  %20lu", ui);
        RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
        RF_nativeExit();
      }
    }
    for (b = 1; b <= RF_ntree; b++) {
      for (uint j = 1; j <= adj; j++) {
        if (RF_mwcpCT_[j][RF_restoreTreeID[b]] > 0) {
          RF_restoreMWCPoffset[j][RF_restoreTreeID[b]] = mwcpOffset[j];
          mwcpOffset[j] = mwcpOffset[j] + RF_mwcpCT_[j][RF_restoreTreeID[b]];
        }
        else {
          RF_restoreMWCPoffset[j][RF_restoreTreeID[b]] = 0;
        }
      }
    }
    free_ulvector(mwcpOffset, 1, adj);
    RF_totalTerminalCount = 0;
    for (b = 1; b <= RF_ntree; b++) {
      RF_totalTerminalCount += (ulong) RF_tLeafCount[b];
    }
    break;
  }
  getConditionalClassificationIndex = &getConditionalClassificationIndexGrow;
  getGMeanIndex = &getGMeanIndexGrow;
  if (mode == RF_PRED) {
    getConditionalClassificationIndex = &getConditionalClassificationIndexPred;
    getGMeanIndex = &getGMeanIndexPred;
  }
  stackMissingArraysPhase2(mode);
  stackDefinedOutputObjects(mode,
                            RF_sexpString,
                            & RF_root,
                            & RF_tLeafCount_,
                            & RF_proximity_,
                            & RF_distance_,
                            & RF_weight_,
                            & RF_imputation_,
                            & RF_sImputeResponsePtr,
                            & RF_sImputePredictorPtr,
                            & RF_varUsed_,
                            & RF_varUsedPtr,
                            & RF_splitDepth_);
  verifyAndRegisterCustomSplitRules();
  if ((RF_optHigh & OPT_MEMB_INCG) || (RF_optHigh & OPT_TERM_INCG)) {
    RF_incStackCount = 0;
    stackAuxiliaryInfoList(&RF_incAuxiliaryInfoList, 8);
  }
  else {
    RF_incStackCount = 0;
    RF_incAuxiliaryInfoList = NULL;
  }
  stackTNQualitativeObjectsKnown(mode,
                                 & RF_RMBR_ID_,
                                 & RF_AMBR_ID_,
                                 & RF_TN_RCNT_,
                                 & RF_TN_ACNT_,
                                 & RF_OOB_SZ_,
                                 & RF_IBG_SZ_);
  stackTNQuantitativeForestObjectsPtrOnly(mode);
  uint bnpSize;
  bnpSize = getVimpRecoverySeedDimension(mode, RF_opt);
  ran1A = &randomChainParallel;
  ran1B = &randomUChainParallel;
  ran1D = &randomChainParallelVimp;
  randomSetChain = &randomSetChainParallel;
  randomSetUChain = &randomSetUChainParallel;
  randomSetChainVimp = &randomSetChainParallelVimp;
  randomGetChain = &randomGetChainParallel;
  randomGetUChain = &randomGetUChainParallel;
  randomGetChainVimp = &randomGetChainParallelVimp;
  randomStack(RF_ntree, bnpSize);
  if (mode == RF_GROW) {
    seedValueLC = abs(seedValue);
    lcgenerator(&seedValueLC, TRUE);
    for (b = 1; b <= RF_ntree; b++) {
      lcgenerator(&seedValueLC, FALSE);
      lcgenerator(&seedValueLC, FALSE);
      while(seedValueLC == 0) {
        lcgenerator(&seedValueLC, FALSE);
      }
      randomSetChain(b, -seedValueLC);
    }
    for (b = 1; b <= RF_ntree; b++) {
      lcgenerator(&seedValueLC, FALSE);
      lcgenerator(&seedValueLC, FALSE);
      while(seedValueLC == 0) {
        lcgenerator(&seedValueLC, FALSE);
      }
      randomSetUChain(b, -seedValueLC);
    }
    if (bnpSize > 0) {
      for (p = 1; p <= bnpSize; p++) {
        lcgenerator(&seedValueLC, FALSE);
        lcgenerator(&seedValueLC, FALSE);
        while(seedValueLC == 0) {
          lcgenerator(&seedValueLC, FALSE);
        }
        randomSetChainVimp(p, -seedValueLC);
      }
    }
  }  
  else {
    for (b = 1; b <= RF_ntree; b++) {
      randomSetChain(b , RF_seed_[b]);
    }
    seedValueLC = abs(seedValue);
    lcgenerator(&seedValueLC, TRUE);
    for (b = 1; b <= RF_ntree; b++) {
      lcgenerator(&seedValueLC, FALSE);
      lcgenerator(&seedValueLC, FALSE);
      while(seedValueLC == 0) {
        lcgenerator(&seedValueLC, FALSE);
      }
      randomSetUChain(b, -seedValueLC);
    }
    if (RF_opt & OPT_VIMP) {
      if ( (RF_opt & (OPT_VIMP | OPT_VIMP_JOIN | OPT_VIMP_TYP1 | OPT_VIMP_TYP2)) ==
           (RF_optLoGrow & (OPT_VIMP | OPT_VIMP_JOIN | OPT_VIMP_TYP1 | OPT_VIMP_TYP2)) )  {
        for (p = 1; p <= bnpSize; p++) {
          randomSetChainVimp(p , RF_seedVimp_[p]);
        }
      }
      else {
        for (p = 1; p <= bnpSize; p++) {
          lcgenerator(&seedValueLC, FALSE);
          lcgenerator(&seedValueLC, FALSE);
          while(seedValueLC == 0) {
            lcgenerator(&seedValueLC, FALSE);
          }
          randomSetChainVimp(p, -seedValueLC);
        }
      }
    }
    else {
    }
  }  
#ifdef _OPENMP
  stackLocksOpenMP(mode);
#endif
  for (r = 1; r <= RF_nImpute; r++) {
    if (getUserTraceFlag()) {
      if (RF_nImpute == 1) {
      }
      else {
        RF_nativePrint("\nImpute Iteration:  %6d", r);
      }
    }
    if (r == RF_nImpute) {
      if (mode == RF_GROW) {
        if (RF_opt & OPT_SEED) {
          for (b = 1; b <= RF_ntree; b++) {
            if (r > 1) {
              lcgenerator(&seedValueLC, FALSE);
              lcgenerator(&seedValueLC, FALSE);
              while(seedValueLC == 0) {
                lcgenerator(&seedValueLC, FALSE);
              }
              randomSetChain(b, -seedValueLC);
            }
            RF_seed_[b] = randomGetChain(b);
          }
          if (bnpSize > 0) {
            for (p = 1; p <= bnpSize; p++) {  
              if (r > 1) {
                lcgenerator(&seedValueLC, FALSE);
                lcgenerator(&seedValueLC, FALSE);
                while(seedValueLC == 0) {
                  lcgenerator(&seedValueLC, FALSE);
                }
                randomSetChainVimp(p, -seedValueLC);
              }
              RF_seedVimp_[p] = randomGetChainVimp(p);
            }
          }
        }
      }
    }  
    for(b = 1; b <= RF_ntree; b++) {
      RF_serialTreeIndex[b] = 0;
    }
    RF_serialTreeID = 0;
    RF_ensbUpdtCount = 0;
    RF_serialBlockID = 0;
    if (getUserTraceFlag()) {
      RF_userTimeStart = time(NULL);
      RF_userTimeSplit = RF_userTreeID = 0;
      RF_nativePrint("\n");
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
    for (b = 1; b <= RF_getTreeCount; b++) {
      acquireTree(mode, r, RF_getTreeIndex[RF_getTreeCount - b + 1]);
    }
    if (r == RF_nImpute) {
      RF_rejectedTreeCount = RF_validTreeCount = RF_stumpedTreeCount = 0;
      RF_totalTerminalCount = 0;
      for (b = 1; b <= RF_ntree; b++) {
        RF_totalTerminalCount += (ulong) RF_tLeafCount[b];
        if (RF_tLeafCount[b] == 0) {
          RF_rejectedTreeCount ++;
        }
        else {
          RF_validTreeCount ++;
          if (RF_tLeafCount[b] == 1) {
            RF_stumpedTreeCount ++;
          }
        }
      }
      if (RF_opt & OPT_PROX) {
        finalizeProximity(mode);
      }
      if (RF_optHigh & OPT_DIST) {
        finalizeDistance(mode);
      }
      if (RF_optHigh & OPT_WGHT) {
        finalizeWeight(mode);
      }  
      stackForestObjectsOutput(mode);
      writeForestObjectsOutput(mode);
      stackTNQuantitativeForestObjectsOutput(mode);
      writeTNQuantitativeForestObjectsOutput(mode);
      stackTNQualitativeObjectsUnknown(mode,
                                       & RF_TN_RCNT_,
                                       & RF_TN_ACNT_,
                                       & RF_TN_OCNT_,
                                       & RF_TN_ICNT_);
      if (RF_optHigh & OPT_MEMB_OUTG) {
        for (b = 1; b <= RF_ntree; b++) {
          RF_OOB_SZ_[b] = RF_oobSize[b];
          RF_IBG_SZ_[b] = RF_ibgSize[b];
        }
      }
      stackTNQualitativeObjectsUnknownMembership(mode, & RF_OMBR_ID_, & RF_IMBR_ID_);
    }  
    if (RF_opt & OPT_MISS_OUT) {
      switch (mode) {
      case RF_PRED:
        imputeSummary(RF_PRED, ACTIVE);
        break;
      default:
        if (r == 1) {
          if ( (!(RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
               ( (RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ) {
            imputeSummary(RF_GROW, FALSE);
            if (RF_timeIndex > 0) {
              if (RF_mTimeFlag == TRUE) {
                imputeMultipleTime(FALSE);
              }
            }
          }
          else {
            imputeSummary(RF_GROW, ACTIVE);
            if (RF_timeIndex > 0) {
              if (RF_mTimeFlag == TRUE) {
                imputeMultipleTime(ACTIVE);
              }
            }
          }
        }  
        else {
          if (r < RF_nImpute) {
            if ( (!(RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
                 ( (RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ) {
              imputeSummary(RF_GROW, FALSE);
              if (RF_timeIndex > 0) {
                if (RF_mTimeFlag == TRUE) {
                  imputeMultipleTime(FALSE);
                }
              }
            }
            else {
              imputeSummary(RF_GROW, ACTIVE);
              if (RF_timeIndex > 0) {
                if (RF_mTimeFlag == TRUE) {
                  imputeMultipleTime(ACTIVE);
                }
              }
            }
          }
          else {
            if (RF_opt & OPT_IMPU_ONLY) {
              if ( (!(RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
                   ( (RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ) {
                imputeSummary(RF_GROW, TRUE);
                if (RF_timeIndex > 0) {
                  if (RF_mTimeFlag == TRUE) {
                    imputeMultipleTime(TRUE);
                  }
                }
              }
              else {
                imputeSummary(RF_GROW, ACTIVE);
                if (RF_timeIndex > 0) {
                  if (RF_mTimeFlag == TRUE) {
                    imputeMultipleTime(ACTIVE);
                  }
                }
              }
            }
            else {
            }
          }
        }
        break;
      }
    }  
    if (r < RF_nImpute) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
      for (uint bb = 1; bb <= RF_getTreeCount; bb++) {
        freeLeafLinkedObjList(RF_leafLinkedObjHead[RF_getTreeIndex[bb]]);
        if (RF_tLeafCount[RF_getTreeIndex[bb]] > 0) {
          free_new_vvector(RF_tTermList[RF_getTreeIndex[bb]], 1, RF_tLeafCount[RF_getTreeIndex[bb]], NRUTIL_TPTR);
        }
        free_new_vvector(RF_tTermMembership[RF_getTreeIndex[bb]], 1, RF_observationSize, NRUTIL_TPTR);
        if (mode == RF_PRED) {
          free_new_vvector(RF_ftTermMembership[RF_getTreeIndex[bb]], 1, RF_fobservationSize, NRUTIL_TPTR);
        }
      }
    }
    if (getUserTraceFlag()) {
      RF_nativePrint("\n\n");
    }
  }  
  if (RF_validTreeCount > 0) {
    if (RF_opt & OPT_VIMP) {
      finalizeVimpPerformance(mode);
    }
    if (RF_optHigh & OPT_PART_PLOT) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
      for (p = 1; p <= RF_partialLength; p++) {
        summarizePartialCalculations(0, p);
      }
    }
    normalizeEnsembleEstimates(mode, TRUE);
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
    for (uint bb = 1; bb <= RF_getTreeCount; bb++) {
      freeLeafLinkedObjList(RF_leafLinkedObjHead[RF_getTreeIndex[bb]]);
      if (RF_tLeafCount[RF_getTreeIndex[bb]] > 0) {
        free_new_vvector(RF_tTermList[RF_getTreeIndex[bb]], 1, RF_tLeafCount[RF_getTreeIndex[bb]], NRUTIL_TPTR);
      }
      free_new_vvector(RF_tTermMembership[RF_getTreeIndex[bb]], 1, RF_observationSize, NRUTIL_TPTR);
      if (mode == RF_PRED) {
        free_new_vvector(RF_ftTermMembership[RF_getTreeIndex[bb]], 1, RF_fobservationSize, NRUTIL_TPTR);
      }
    }
    if (RF_opt & (OPT_VARUSED_F | OPT_VARUSED_T)) {
      if (RF_opt & OPT_VARUSED_F) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
        for (uint jj = 1; jj <= RF_xSize; jj++) {
          RF_varUsed_[jj] = 0;
          for (uint bb = 1; bb <= RF_ntree; bb++) {
            RF_varUsed_[jj] += RF_varUsedPtr[bb][jj];
          }
        }
      }
      else {
      }
    }
    if (RF_opt & (OPT_SPLDPTH_1 | OPT_SPLDPTH_2)) {
      if (RF_opt & OPT_SPLDPTH_1) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
        for (uint jj = 1; jj <= RF_xSize; jj++) {
          for (uint ii = 1; ii <= RF_observationSize; ii++) {
            RF_splitDepthPtr[1][jj][ii] = RF_splitDepthPtr[1][jj][ii] / (RF_validTreeCount);
          }
        }
      }
      else {
      }
    }
  }  
  else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
    for (uint bb = 1; bb <= RF_getTreeCount; bb++) {
      freeLeafLinkedObjList(RF_leafLinkedObjHead[RF_getTreeIndex[bb]]);
      free_new_vvector(RF_tTermMembership[RF_getTreeIndex[bb]], 1, RF_observationSize, NRUTIL_TPTR);
      if (mode == RF_PRED) {
        free_new_vvector(RF_ftTermMembership[RF_getTreeIndex[bb]], 1, RF_fobservationSize, NRUTIL_TPTR);
      }
    }
    RF_nativePrint("\nRF-SRC:  *** WARNING *** ");
    RF_nativePrint("\nRF-SRC:  Insufficient trees for analysis.  \n");
  }
  switch (mode) {
  case RF_GROW:
    freeSplitRuleObj(RF_splitRuleObj);
    unstackForestObjectsPtrOnly(mode);
    unstackTNQuantitativeForestObjectsPtrOnly(mode);
    unstackForestObjectsAuxOnly(mode);
    break;
  default:
    break;
  }
  unstackAuxiliaryInfoAndList((mode == RF_GROW) ? FALSE : TRUE, RF_snpAuxiliaryInfoList, RF_stackCount);
  if ((RF_optHigh & OPT_MEMB_INCG) || (RF_optHigh & OPT_TERM_INCG)) {
    unstackAuxiliaryInfoAndList(FALSE, RF_incAuxiliaryInfoList, 8);
  }
  unstackDefinedOutputObjects(mode);
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    unstackCompetingArrays(mode);
  }
  if (RF_rFactorCount > 0) {
    unstackClassificationArrays(mode);
  }
  unstackMissingArrays(mode);
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    unstackTimeAndSubjectArrays(mode);
  }
  switch (mode) {
  case RF_PRED:
    unstackPreDefinedPredictArrays();
    break;
  case RF_REST:
    unstackPreDefinedRestoreArrays();
    break;
  default:
    unstackPreDefinedGrowthArrays();
    break;
  }
  unstackPreDefinedCommonArrays(mode,
                                RF_nodeMembership,
                                RF_tTermMembership,
                                RF_tTermList,
                                RF_root);
  unstackIncomingArrays(mode);
  randomUnstack(RF_ntree, bnpSize);
#ifdef _OPENMP
  unstackLocksOpenMP(mode);
#endif
  unstackFactorArrays(mode);
}
