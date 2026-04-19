
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "tree.h"
#include "impute.h"
#include "stackOutput.h"
#include "treeUtil.h"
#include "treeJIT.h"
#include "processEnsemble.h"
#include "splitGreedy.h"
#include "regression.h"
#include "importance.h"
#include "stackOutputQQ.h"
#include "nodeOps.h"
#include "node.h"
#include "nativeUtil.h"
#include "nrutil.h"
#include "error.h"
void acquireTreeGeneric(char mode, uint r, uint b) {
  char multImpFlag;
  uint *fallMembrIndx;
  uint bootMembrIndxIter;
  uint rmbrIterator;
  uint ambrIterator;
  char result;
  uint offset;
  uint i;
  ${trace.token}  uint ibgSampleSize;
  ${trace.token}  if (getTraceFlag(b) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nTree ran1A():  %20d", randomGetChain(b));
  ${trace.token}    RF_nativePrint("\nTree ran1B():  %20d", randomGetUChain(b));
  ${trace.token}  }
#ifdef _OPENMP
  ${trace.token}  if (getTraceFlag(b) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\n Hello from thread %d of %d", omp_get_thread_num(), omp_get_num_threads());
  ${trace.token}    RF_nativePrint("\n Parallel treeID:  %10d ", b);
  ${trace.token}  }
#endif
  ${trace.token}    if (getTraceFlag(b) & SUMM_LOW_TRACE) {
  ${trace.token}      RF_nativePrint("\nAllocating for treeID:  %10d", b);
  ${trace.token}    }
  RF_root[b] = makeNode((mode == RF_GROW) ? RF_xSize : 0);
  for (i = 1; i <= RF_root[b] -> xSize; i++) {
    RF_root[b] -> permissible[i] = TRUE;
  }
  stackAuxiliary(mode, b);
  if (mode == RF_PRED) {
    fallMembrIndx = uivector(1, RF_fobservationSize);
  }
  else {
    fallMembrIndx = NULL;
  }
  RF_root[b] -> allMembrSizeAlloc = RF_root[b] -> allMembrSize = RF_observationSize;
  RF_root[b] -> allMembrIndx = uivector(1, RF_root[b] -> allMembrSizeAlloc);
  stackShadow(mode, b);
  stackFactorInSitu(b);
  RF_tTermMembership[b] = (Terminal **) new_vvector(1, RF_observationSize, NRUTIL_TPTR);
  RF_leafLinkedObjHead[b] = RF_leafLinkedObjTail[b] = makeLeafLinkedObj();
  if (mode == RF_PRED) {
    RF_ftTermMembership[b] = (Terminal **) new_vvector(1, RF_fobservationSize, NRUTIL_TPTR);
  }
  if (mode == RF_GROW) {
    if (RF_nImpute > 1) {
      if (r > 1) {
        if (RF_mRecordSize > 0) {
          imputeUpdateShadow(RF_GROW,
                             RF_response[b],
                             RF_observation[b]);
        }
        if (RF_timeIndex > 0) {
          if (RF_mTimeFlag == TRUE) {
            updateTimeIndexArray(0,
                                 NULL,
                                 RF_observationSize,
                                 RF_time[b],
                                 FALSE,
                                 FALSE,
                                 RF_masterTimeIndex[b]);
          }
        }
      }  
    }  
  }  
  RF_root[b] -> parent = NULL;
  RF_root[b] -> nodeID = 1;
  RF_maxDepth[b] = 0;
  bootMembrIndxIter = 0;
  ${trace.token}  if (getTraceFlag(b) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nBootstrap sample:  %10d ", b);
  ${trace.token}  }
  for (i = 1; i <= RF_observationSize; i++) {
    RF_root[b] -> allMembrIndx[i] = i;
    RF_nodeMembership[b][i] = RF_root[b];
  }
  if (RF_optHigh & OPT_MEMB_PRUN) {
    Node     ***gNodeMembership;
    uint        obsSize;
    switch (mode) {
    case RF_PRED:
      obsSize = RF_fobservationSize;
      gNodeMembership = RF_fnodeMembership;
      break;
    default:
      obsSize = RF_observationSize;
      gNodeMembership = RF_nodeMembership;
      break;
    }
    for (i = 1; i <= obsSize; i++) {
      RF_pNodeMembership[b][i] = gNodeMembership[b][i];
    }
  }
  multImpFlag = FALSE;
  if (mode == RF_GROW) {
    if (r > 1) {
      multImpFlag = TRUE;
    }
  }
  if ((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP)) {
    getVariance = & getVarianceClassicNoMiss;
  }
  else {
    getVariance = & getVarianceClassic;
  }
  switch (mode) {
  case RF_GROW:
    ${trace.token}      if (getTraceFlag(b) & SUMM_LOW_TRACE) {
    ${trace.token}        RF_nativePrint("\nStart of growTree():  (MI iter %10d, tree %10d)", r, b);
    ${trace.token}      }
    RF_tLeafCount[b] = 0;
    rmbrIterator = ambrIterator = 0;
    result = growTree (r,
                       TRUE,
                       multImpFlag,
                       b,
                       RF_root[b],
                       & bootMembrIndxIter,
                       & rmbrIterator,
                       & ambrIterator);
    if (result) {
      RF_nodeCount[b]  = (RF_tLeafCount[b] << 1) - 1;
    }
    else {
      RF_nodeCount[b] = 1;
    }
    ${trace.token}    if (getTraceFlag(b) & SUMM_LOW_TRACE) {
    ${trace.token}      RF_nativePrint("\nEnd of grow tree:  (MII r, treeID) = (%10d, %10d)", r, b);
    ${trace.token}      RF_nativePrint("\nTotal node count:  %10d", RF_nodeCount[b]);
    ${trace.token}      RF_nativePrint("\nTotal leaf count:  %10d", RF_tLeafCount[b]);
    ${trace.token}    }
    ${trace.token}    if (getTraceFlag(b) & SUMM_MED_TRACE) {
    ${trace.token}      RF_nativePrint("\nTerminal Object Pointers:  %10d", b);
    ${trace.token}      RF_nativePrint("\n       index              address       nodeID       termID\n");
    ${trace.token}      LeafLinkedObj *leafLinkedPtr;
    ${trace.token}      leafLinkedPtr = RF_leafLinkedObjHead[b] -> fwdLink;
    ${trace.token}      i = 0;
    ${trace.token}      while (leafLinkedPtr != NULL) {
    ${trace.token}        i++;
    ${trace.token}        RF_nativePrint("%12d %20x %12d %12d \n", i, leafLinkedPtr, (leafLinkedPtr -> nodePtr) -> nodeID, (leafLinkedPtr -> termPtr) -> nodeID);
    ${trace.token}        leafLinkedPtr = leafLinkedPtr -> fwdLink;
    ${trace.token}      }
    ${trace.token}    }
    break;
  default:
    if (mode == RF_PRED) {
      for (i = 1; i <= RF_fobservationSize; i++) {
        fallMembrIndx[i] = i;
        RF_fnodeMembership[b][i] = RF_root[b];
      }
    }
    ${trace.token}     if (getTraceFlag(b) & SUMM_LOW_TRACE) {
    ${trace.token}       RF_nativePrint("\nStart of restore process:  (%10d, %10d)", r, b);
    ${trace.token}     }
    ${trace.token}     if (getTraceFlag(b) & SPLT_DEF_TRACE) {
    ${trace.token}        RF_nativePrint("\nIncoming Tree Record (partial):  ");
    ${trace.token}        RF_nativePrint("\n    offset     treeID     nodeID   brnodeID");
    ${trace.token}        RF_nativePrint("     parmID     mwcpSZ \n");
    ${trace.token}      }
    ${trace.token}      if (getTraceFlag(b) & SUMM_LOW_TRACE) {
    ${trace.token}        RF_nativePrint("\nStart of restore tree:  %10d", b);
    ${trace.token}      }
    restoreTree(mode, b, RF_root[b]);
    ${trace.token}      if (getTraceFlag(b) & SUMM_LOW_TRACE) {
    ${trace.token}        RF_nativePrint("\nEnd of restoreTree():  %10d", b);
    ${trace.token}        RF_nativePrint("\nTotal node count:  %10d", RF_nodeCount[b]);
    ${trace.token}        RF_nativePrint("\nTotal leaf count:  %10d", RF_tLeafCount[b]);
    ${trace.token}      }
    ${trace.token}      if (getTraceFlag(b) & SUMM_LOW_TRACE) {
    ${trace.token}        RF_nativePrint("\nBeginning of test (and/or) train Membership Restoration:  (%10d, %10d)", r, b);
    ${trace.token}      }
    rmbrIterator = ambrIterator = 0;
    ${trace.token}      if (getTraceFlag(b) & SUMM_LOW_TRACE) {
    ${trace.token}        RF_nativePrint("\nStart of test (and/or) train Node Membership Restoration:  (%10d, %10d)", r, b);
    ${trace.token}      }
    result = restoreNodeMembership(mode,
                                   TRUE,
                                   b,
                                   RF_root[b],
                                   NULL,
                                   0,
                                   RF_root[b] -> allMembrIndx,
                                   RF_root[b] -> allMembrSize,
                                   fallMembrIndx,
                                   RF_fobservationSize,
                                   & bootMembrIndxIter,
                                   & rmbrIterator,
                                   & ambrIterator);
    ${trace.token}      if (getTraceFlag(b) & SUMM_LOW_TRACE) {
    ${trace.token}        RF_nativePrint("\nEnd of restoreNodeMembership():  %10d", b);
    ${trace.token}        RF_nativePrint("\nTotal node count:  %10d", RF_nodeCount[b]);
    ${trace.token}        RF_nativePrint("\nTotal leaf count:  %10d", RF_tLeafCount[b]);
    ${trace.token}      }
    ${trace.token}      if (getTraceFlag(b) & SUMM_LOW_TRACE) {
    ${trace.token}        RF_nativePrint("\nEnd of test (and/or) train Node Membership Restoration:  (%10d, %10d)", r, b);
    ${trace.token}      }
    ${trace.token}    if (getTraceFlag(b) & SUMM_MED_TRACE) {
    ${trace.token}      RF_nativePrint("\nTerminal Object Pointers:  %10d", b);
    ${trace.token}      RF_nativePrint("\n       index        leafLinkedPtr        nodeID       termID\n");
    ${trace.token}      LeafLinkedObj *leafLinkedPtr;
    ${trace.token}      leafLinkedPtr = RF_leafLinkedObjHead[b] -> fwdLink;
    ${trace.token}      i = 0;
    ${trace.token}      while (leafLinkedPtr != NULL) {
    ${trace.token}        i++;
    ${trace.token}        RF_nativePrint("%12d %20x %12d %12d \n", i, leafLinkedPtr, (leafLinkedPtr -> nodePtr) -> nodeID, (leafLinkedPtr -> termPtr) -> nodeID);
    ${trace.token}        leafLinkedPtr = leafLinkedPtr -> fwdLink;
    ${trace.token}      }
    ${trace.token}    }
    break;
  }
  if (result) {
    RF_tTermList[b] = (Terminal **) new_vvector(1, RF_tLeafCount[b], NRUTIL_TPTR);
    LeafLinkedObj *leafLinkedPtr;
    leafLinkedPtr = RF_leafLinkedObjHead[b] -> fwdLink;
    while (leafLinkedPtr != NULL) {
      RF_tTermList[b][(leafLinkedPtr -> termPtr) -> nodeID] = leafLinkedPtr -> termPtr;
      leafLinkedPtr = leafLinkedPtr -> fwdLink;
    }
    ${trace.token}  ibgSampleSize = RF_observationSize - RF_oobSize[b];
    if (mode != RF_PRED) {
      if (RF_mRecordSize > 0) {
        for (i = 1; i <= RF_mRecordSize; i++) {
          if (RF_bootMembershipFlag[b][RF_mRecordIndex[i]] == TRUE) {
            RF_dmRecordBootFlag[b][i] = TRUE;
          }
          else {
            RF_dmRecordBootFlag[b][i] = FALSE;
          }
        }
      }  
    }  
    if (r == RF_nImpute) {      
      if ( (RF_vtry == 0) || ((RF_vtry > 0) && (RF_vtryMode == RF_VTRY_NULL)) ) {
        if (RF_ptnCount > 0) {
          updatePruning(mode, b);
        }
        processEnsembleInSitu(mode, multImpFlag, b);
      }
      else {
        for (uint p = 1; p <= RF_xSize; p++) {
          processEnsembleHoldout(p, b);
        }
      }
    }
    ${trace.token}    if (getTraceFlag(b) & SUMM_MED_TRACE) {
    ${trace.token}      RF_nativePrint("\n\nFinal Bootstrap Information:  ");
    ${trace.token}      RF_nativePrint("\n  orgIndex       flag \n");
    ${trace.token}      RF_nativePrint("\n  orgIndex   bootFlag    oobFlag  bootCount \n");
    ${trace.token}      for (i = 1; i <= RF_observationSize; i++) {
    ${trace.token}        RF_nativePrint("%10d %10d %10d %10d \n", i, RF_bootMembershipFlag[b][i], RF_oobMembershipFlag[b][i], RF_bootMembershipCount[b][i]);
    ${trace.token}      }
    ${trace.token}      RF_nativePrint("\nIBG Size:  %10d ",   ibgSampleSize);
    ${trace.token}      RF_nativePrint("\nOOB Size:  %10d \n", RF_oobSize[b]);
    ${trace.token}      RF_nativePrint("\nFinal (Terminal Object) Membership (in-bag data):  %10d", b);
    ${trace.token}      RF_nativePrint("\n        index         leaf");
    ${trace.token}      for (i = 1; i <= RF_ibgSize[b]; i++) {
    ${trace.token}        RF_nativePrint("\n %12d %12d", i, RF_tTermMembership[b][RF_ibgMembershipIndex[b][i]] -> nodeID);
    ${trace.token}      }
    ${trace.token}      RF_nativePrint("\nFinal (Terminal Object) Membership (oob data):  %10d", b);
    ${trace.token}      RF_nativePrint("\n        index         leaf");
    ${trace.token}      for (i = 1; i <= RF_oobSize[b]; i++) {
    ${trace.token}        RF_nativePrint("\n %12d %12d", i, RF_tTermMembership[b][RF_oobMembershipIndex[b][i]]-> nodeID);
    ${trace.token}      }
    ${trace.token}      if (RF_optHigh & OPT_FENS) {
    ${trace.token}        RF_nativePrint("\nFinal (Terminal Object) Membership (all data):  %10d", b);
    ${trace.token}        RF_nativePrint("\n        index         leaf");
    ${trace.token}        for (i=1; i <= RF_observationSize; i++) {
    ${trace.token}          RF_nativePrint("\n %12d %12d", i, RF_tTermMembership[b][i] -> nodeID);
    ${trace.token}        }
    ${trace.token}      }
    ${trace.token}      if (mode == RF_PRED) {
    ${trace.token}        RF_nativePrint("\nFinal PRED Membership (all data):  %10d", b);
    ${trace.token}        RF_nativePrint("\n     index       leaf\n");
    ${trace.token}        for (i=1; i <=  RF_fobservationSize; i++) {
    ${trace.token}          RF_nativePrint("%10d %10d \n", i, RF_fnodeMembership[b][i] -> nodeID);
    ${trace.token}        }
    ${trace.token}      }
    ${trace.token}    }
  }  
  if (r == RF_nImpute) {
    if (mode == RF_GROW) {
      if (RF_opt & OPT_TREE) {
        offset = 0;
        stackTreeObjectsPtrOnly(mode, b);
        saveTree(b, RF_root[b], & offset);
        if (RF_tLeafCount[b] > 0) {
          stackTNQuantitativeTreeObjectsPtrOnly(b);
          saveTNQuantitativeTreeObjects(b);
        }          
      }  
    }  
  }  
  ${trace.token}    if (getTraceFlag(b) & SUMM_MED_TRACE) {
  ${trace.token}      RF_nativePrint("\nDe-allocating for treeID:  %10d", b);
  ${trace.token}    }
  unstackShadow(mode, b);
  unstackFactorInSitu(b);
  if (mode == RF_PRED) {
    free_uivector(fallMembrIndx, 1, RF_fobservationSize);
  }
  unstackAuxiliary(mode, b);
  freeTree(b, RF_root[b]);
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
        RF_nativePrint("Trees Grown:  %6d,    Time Remaining (sec):  %6.0f \n", RF_userTreeID, ceil(userTimeRemaining));
        RF_userTimeSplit = current;
      }
    }  
  }
  ${trace.token}  if (getTraceFlag(b) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nEnd of acquisition:  (%10d, %10d)", r, b);
  ${trace.token}  }
}
void finalizeWeight(char mode) {
  uint    obsSize;
  ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nfinalizeWeight() ENTRY ...\n");
  ${trace.token}  }
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    break;
  default:
    obsSize = RF_observationSize;
    break;
  }
  ${trace.token}  if((RF_optHigh & OPT_WGHT_IBG) && (RF_optHigh & OPT_WGHT_OOB)) {
  ${trace.token}    if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}      RF_nativePrint("\n  Weight:  FULL  \n");
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}  else {
  ${trace.token}    if((RF_optHigh & OPT_WGHT_IBG)  && !(RF_optHigh & OPT_WGHT_OOB)) {
  ${trace.token}      if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}        RF_nativePrint("\n  Weight:  IBG  \n");
  ${trace.token}      }
  ${trace.token}    }
  ${trace.token}    else if(!(RF_optHigh & OPT_WGHT_IBG)  && (RF_optHigh & OPT_WGHT_OOB)) {
  ${trace.token}      if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}        RF_nativePrint("\n  Weight:  OOB  \n");
  ${trace.token}      }
  ${trace.token}    }
  ${trace.token}    else {
  ${trace.token}      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
  ${trace.token}      RF_nativeError("\nRF-SRC:  Illegal finalizeWeight() call.");
  ${trace.token}      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
  ${trace.token}      RF_nativeExit();
  ${trace.token}    }
  ${trace.token}  }
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
  for (uint i = 1; i <= obsSize; i++) {
    if (RF_weightDenom[i] > 0) {    
      for (uint j = 1; j <= RF_observationSize; j++) {
        RF_weightPtr[i][j] = RF_weightPtr[i][j] /  RF_weightDenom[i];
      }
    }
    else {        
      for (uint j = 1; j <= RF_observationSize; j++) {
        RF_weightPtr[i][j] = RF_nativeNaN;
      }
    }
  }
  ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nfinalizeWeight() EXIT ...\n");
  ${trace.token}  }
}
void updateWeight(char mode, uint b) {
  uint     **utTermMembership;
  uint      *utTermMembershipCount;
  Terminal **itTermMembership, **gtTermMembership;
  uint  *gMembershipIndex, *iMembershipIndex;
  uint   gMembershipSize,   iMembershipSize;
  uint  *bootMembershipCount;
  uint  mtnmFlag;
  ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdateWeight() ENTRY ...\n");
  ${trace.token}  }
  gMembershipIndex = NULL;
  gMembershipSize  = 0;
  gtTermMembership = NULL;
  if((RF_optHigh & OPT_WGHT_IBG) && (RF_optHigh & OPT_WGHT_OOB)) {
    ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
    ${trace.token}    RF_nativePrint("\n  Weight:  FULL  \n");
    ${trace.token}  }
    switch (mode) {
    case RF_PRED:
      gMembershipSize = RF_fobservationSize;
      gMembershipIndex = RF_fidentityMembershipIndex;
      gtTermMembership = RF_ftTermMembership[b];
      break;
    default:
      gMembershipSize = RF_observationSize;
      gMembershipIndex = RF_identityMembershipIndex;
      gtTermMembership = RF_tTermMembership[b];
      break;
    }
  }
  else {
    if((RF_optHigh & OPT_WGHT_IBG)  && !(RF_optHigh & OPT_WGHT_OOB)) {
      ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
      ${trace.token}    RF_nativePrint("\n  Weight:  IBG  \n");
      ${trace.token}  }
      gMembershipSize  = RF_observationSize;
      gMembershipIndex = RF_identityMembershipIndex;
      gtTermMembership = RF_tTermMembership[b];
    }
    else if(!(RF_optHigh & OPT_WGHT_IBG)  && (RF_optHigh & OPT_WGHT_OOB)) {
      ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
      ${trace.token}    RF_nativePrint("\n  Weight:  OOB  \n");
      ${trace.token}  }
      gMembershipIndex = RF_oobMembershipIndex[b];
      gMembershipSize  = RF_oobSize[b];
      gtTermMembership = RF_tTermMembership[b];
    }
    else {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Illegal updateWeight() call.");
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  iMembershipIndex = RF_ibgMembershipIndex[b];
  iMembershipSize  = RF_ibgSize[b];
  bootMembershipCount = RF_bootMembershipCount[b];
  itTermMembership = RF_tTermMembership[b];
  if (RF_xMarginalSize > 0) {
    mtnmFlag = TRUE;
  }
  else {
    mtnmFlag = FALSE;
  }
  if (!mtnmFlag) {
    for (uint i = 1; i <= gMembershipSize; i++) {
      uint ii, jj;
      ii = gMembershipIndex[i];
#ifdef _OPENMP
      omp_set_lock(&(RF_lockWeightRow[ii]));
#endif
      RF_weightDenom[ii] ++;
      for (uint j = 1; j <= iMembershipSize; j++) {
        jj = iMembershipIndex[j];
        if ( gtTermMembership[ii] == itTermMembership[jj] ) {
          RF_weightPtr[ii][jj] +=  (double) bootMembershipCount[jj] / (double) (gtTermMembership[ii] -> membrCount);
        }
      }
#ifdef _OPENMP
      omp_unset_lock(&(RF_lockWeightRow[ii]));
#endif
    }
  }
  else {
    utTermMembershipCount = RF_utTermMembershipCount[b];
    utTermMembership      = RF_utTermMembership[b];
    for (uint i = 1; i <= gMembershipSize; i++) {
      uint rowNodeID, colNodeID;
      uint ii, jj;
      uint xMembrCount;
      ii = gMembershipIndex[i];
#ifdef _OPENMP
      omp_set_lock(&(RF_lockWeightRow[ii]));
#endif
      RF_weightDenom[ii] ++;
      for (uint j = 1; j <= iMembershipSize; j++) {
        jj = iMembershipIndex[j];
        colNodeID = itTermMembership[jj] -> nodeID;
        xMembrCount = 0;
        for (uint k = 1; k <= utTermMembershipCount[ii]; k++) {
          xMembrCount += RF_tTermList[b][utTermMembership[ii][k]] -> membrCount;
          rowNodeID = utTermMembership[ii][k];
          if ( colNodeID == rowNodeID ) {
            for (uint kk = k+1; kk <= utTermMembershipCount[ii]; kk++) {
              xMembrCount += RF_tTermList[b][utTermMembership[ii][kk]] -> membrCount;
            }
            RF_weightPtr[ii][jj] +=  (double) bootMembershipCount[jj] / (double) xMembrCount;
            goto wghtMarginal;
          }
        }
      wghtMarginal:
        continue;
      }
#ifdef _OPENMP
      omp_unset_lock(&(RF_lockWeightRow[ii]));
#endif
    }
  }
  ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdateWeight() EXIT ...\n");
  ${trace.token}  }
}
void finalizeProximity(char mode) {
  uint  obsSize;
  ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nfinalizeProximity() ENTRY ...\n");
  ${trace.token}  }
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    break;
  default:
    obsSize = RF_observationSize;
    break;
  }
  ${trace.token}  if((RF_opt & OPT_PROX_IBG) && (RF_opt & OPT_PROX_OOB)) {
  ${trace.token}    if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}      RF_nativePrint("\n  Proximity:  FULL  \n");
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}  else {
  ${trace.token}    if((RF_opt & OPT_PROX_IBG)  && !(RF_opt & OPT_PROX_OOB)) {
  ${trace.token}      if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}        RF_nativePrint("\n  Proximity:  IBG  \n");
  ${trace.token}      }
  ${trace.token}    }
  ${trace.token}    else if(!(RF_opt & OPT_PROX_IBG)  && (RF_opt & OPT_PROX_OOB)) {
  ${trace.token}      if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}        RF_nativePrint("\n  Proximity:  OOB  \n");
  ${trace.token}      }
  ${trace.token}    }
  ${trace.token}    else {
  ${trace.token}      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
  ${trace.token}      RF_nativeError("\nRF-SRC:  Illegal finalizeProximity() call.");
  ${trace.token}      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
  ${trace.token}      RF_nativeExit();
  ${trace.token}    }
  ${trace.token}  }
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
  for (uint i = 1; i <= obsSize; i++) {
    for (uint j = 1; j <= i; j++) {
      if (RF_proximityDenPtr[i][j] > 0) {
        RF_proximityPtr[i][j] = RF_proximityPtr[i][j] /  RF_proximityDenPtr[i][j];
      }
      else {
        RF_proximityPtr[i][j] = RF_nativeNaN;
      }
    }
  }
  ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nfinalizeProximity() EXIT ...\n");
  ${trace.token}  }
}
void updateProximity(char mode, uint b) {
  uint     **utTermMembership;
  uint      *utTermMembershipCount;
  Terminal **tTermMembership;
  uint  *membershipIndex;
  uint   membershipSize;
  uint  mtnmFlag;
  ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdateProximity() ENTRY ...\n");
  ${trace.token}  }
  membershipSize  = 0;  
  membershipIndex = NULL;  
  if((RF_opt & OPT_PROX_IBG) && (RF_opt & OPT_PROX_OOB)) {
    ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
    ${trace.token}    RF_nativePrint("\n  Proximity:  FULL  \n");
    ${trace.token}  }
    switch (mode) {
    case RF_PRED:
      membershipSize = RF_fobservationSize;
      membershipIndex = RF_fidentityMembershipIndex;
      tTermMembership = RF_ftTermMembership[b];
      break;
    default:
      membershipSize = RF_observationSize;
      membershipIndex = RF_identityMembershipIndex;
      tTermMembership = RF_tTermMembership[b];
      break;
    }
  }
  else {
    if((RF_opt & OPT_PROX_IBG)  && !(RF_opt & OPT_PROX_OOB)) {
      ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
      ${trace.token}    RF_nativePrint("\n  Proximity:  IBG  \n");
      ${trace.token}  }
      membershipIndex = RF_ibgMembershipIndex[b];
      membershipSize  = RF_ibgSize[b];
    }
    else if(!(RF_opt & OPT_PROX_IBG)  && (RF_opt & OPT_PROX_OOB)) {
      ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
      ${trace.token}    RF_nativePrint("\n  Proximity:  OOB  \n");
      ${trace.token}  }
      membershipIndex = RF_oobMembershipIndex[b];
      membershipSize  = RF_oobSize[b];
    }
    else {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Illegal updateProximity() call.");
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
    tTermMembership = RF_tTermMembership[b];
  }
  if (RF_xMarginalSize > 0) {
    mtnmFlag = TRUE;
  }
  else {
    mtnmFlag = FALSE;
  }
  if (!mtnmFlag) {
    for (uint i = 1; i <= membershipSize; i++) {
      uint ii, jj;
      ii = membershipIndex[i];
      for (uint j = 1; j <= i; j++) {
        jj = membershipIndex[j];
        rfsrc_omp_atomic_update(&RF_proximityDenPtr[ii][jj], 1.0);
        if ( tTermMembership[ii] == tTermMembership[jj] ) {
          rfsrc_omp_atomic_update(&RF_proximityPtr[ii][jj], 1.0);
        }
      }
    }
  }
  else {
    utTermMembership =  RF_utTermMembership[b];
    utTermMembershipCount =  RF_utTermMembershipCount[b];
    for (uint i = 1; i <= membershipSize; i++) {
      uint ii, jj;
      ii = membershipIndex[i];
      for (uint j = 1; j <= i; j++) {
        jj = membershipIndex[j];
        rfsrc_omp_atomic_update(&RF_proximityDenPtr[ii][jj], 1.0);
        for (uint ki = 1; ki <= utTermMembershipCount[ii]; ki++) {
          for (uint kj = 1; kj <= utTermMembershipCount[jj]; kj++) {
            if ( utTermMembership[ii][ki] == utTermMembership[jj][kj] ) {
              rfsrc_omp_atomic_update(&RF_proximityPtr[ii][jj], 1.0);
              goto proxMarginal;
            }
          }
        }
      proxMarginal:
        continue;
      }
    }
  }
  ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdateProximity() EXIT ...\n");
  ${trace.token}  }
}
void finalizeDistance(char mode) {
  uint  obsSize;
  ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nfinalizeDistance() ENTRY ...\n");
  ${trace.token}  }
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    break;
  default:
    obsSize = RF_observationSize;
    break;
  }
  ${trace.token}  if((RF_optHigh & OPT_DIST_IBG) && (RF_optHigh & OPT_DIST_OOB)) {
  ${trace.token}    if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}      RF_nativePrint("\n  Distance:  FULL  \n");
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}  else {
  ${trace.token}    if((RF_optHigh & OPT_DIST_IBG)  && !(RF_optHigh & OPT_DIST_OOB)) {
  ${trace.token}      if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}        RF_nativePrint("\n  Distance:  IBG  \n");
  ${trace.token}      }
  ${trace.token}    }
  ${trace.token}    else if(!(RF_optHigh & OPT_DIST_IBG)  && (RF_optHigh & OPT_DIST_OOB)) {
  ${trace.token}      if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}        RF_nativePrint("\n  Distance:  OOB  \n");
  ${trace.token}      }
  ${trace.token}    }
  ${trace.token}    else {
  ${trace.token}      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
  ${trace.token}      RF_nativeError("\nRF-SRC:  Illegal finalizeDistance() call.");
  ${trace.token}      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
  ${trace.token}      RF_nativeExit();
  ${trace.token}    }
  ${trace.token}  }
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
  for (uint i = 1; i <= obsSize; i++) {
    for (uint j = 1; j <= i; j++) {
      if (RF_distanceDenPtr[i][j] > 0) {
        RF_distancePtr[i][j] = RF_distancePtr[i][j] /  RF_distanceDenPtr[i][j];
      }
      else {
        RF_distancePtr[i][j] = RF_nativeNaN;
      }
    }
  }
  ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nfinalizeDistance() EXIT ...\n");
  ${trace.token}  }
}
void updateDistance(char mode, uint b) {
  uint     **utTermMembership;
  uint      *utTermMembershipCount;
  Terminal **tTermMembership;
  uint  *membershipIndex;
  uint   membershipSize;
  uint  mtnmFlag;
  ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdateDistance() ENTRY ...\n");
  ${trace.token}  }
  membershipSize  = 0;  
  membershipIndex = NULL;  
  if((RF_optHigh & OPT_DIST_IBG) && (RF_optHigh & OPT_DIST_OOB)) {
    ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
    ${trace.token}    RF_nativePrint("\n  Distance:  FULL  \n");
    ${trace.token}  }
    switch (mode) {
    case RF_PRED:
      membershipSize = RF_fobservationSize;
      membershipIndex = RF_fidentityMembershipIndex;
      tTermMembership = RF_ftTermMembership[b];
      break;
    default:
      membershipSize = RF_observationSize;
      membershipIndex = RF_identityMembershipIndex;
      tTermMembership = RF_tTermMembership[b];
      break;
    }
  }
  else {
    if((RF_optHigh & OPT_DIST_IBG)  && !(RF_optHigh & OPT_DIST_OOB)) {
      ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
      ${trace.token}    RF_nativePrint("\n  Distance:  IBG  \n");
      ${trace.token}  }
      membershipIndex = RF_ibgMembershipIndex[b];
      membershipSize  = RF_ibgSize[b];
    }
    else if(!(RF_optHigh & OPT_DIST_IBG)  && (RF_optHigh & OPT_DIST_OOB)) {
      ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
      ${trace.token}    RF_nativePrint("\n  Distance:  OOB  \n");
      ${trace.token}  }
      membershipIndex = RF_oobMembershipIndex[b];
      membershipSize  = RF_oobSize[b];
    }
    else {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Illegal updateDistance() call.");
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
    tTermMembership = RF_tTermMembership[b];
  }
  if (RF_xMarginalSize > 0) {
    mtnmFlag = TRUE;
  }
  else {
    mtnmFlag = FALSE;
  }
  if (!mtnmFlag) {
    for (uint i = 1; i <= membershipSize; i++) {
      Node *iNodeMembership, *jNodeMembership;
      Node *deepNodeMembership, *shallowNodeMembership;
      uint  iEdgeCount, jEdgeCount;
      uint  iDepthCount, jDepthCount;
      uint  *deepEdgeCount;
      double realDistance;
      uint ii, jj;
      ii = membershipIndex[i];
      iNodeMembership = tTermMembership[ii] -> mate;
      iDepthCount = iNodeMembership -> depth;
      for (uint j = 1; j <= i; j++) {
        jj = membershipIndex[j];
        rfsrc_omp_atomic_update(&RF_distanceDenPtr[ii][jj], 1.0);
        jNodeMembership = tTermMembership[jj] -> mate;
        jDepthCount = jNodeMembership -> depth;
        iEdgeCount = jEdgeCount = 0;
        if ((iNodeMembership -> depth) > (jNodeMembership -> depth)) {
          deepNodeMembership = iNodeMembership;
          shallowNodeMembership = jNodeMembership;
          deepEdgeCount = & iEdgeCount;
        }
        else {
          deepNodeMembership = jNodeMembership;
          shallowNodeMembership = iNodeMembership;
          deepEdgeCount = & jEdgeCount;
        }
        ${trace.token}    if (getTraceFlag(0) & ENSB_LOW_TRACE) {
        ${trace.token}      RF_nativePrint("\n  Distance : CaseID (i, j) = (%10d, %10d)", ii, jj);
        ${trace.token}      RF_nativePrint("\n  Distance : NodeID (i, j) = (%10d, %10d)", iNodeMembership -> nodeID, jNodeMembership -> nodeID);
        ${trace.token}      RF_nativePrint("\n  Distance : Depth  (i, j) = (%10d, %10d)", iNodeMembership -> depth, jNodeMembership -> depth);
        ${trace.token}    }
        while ((deepNodeMembership -> depth) > (shallowNodeMembership -> depth)) {
          deepNodeMembership = deepNodeMembership -> parent;
          (*deepEdgeCount) ++;
          ${trace.token}    if (getTraceFlag(0) & ENSB_LOW_TRACE) {
          ${trace.token}      RF_nativePrint("\n  Distance : Parsing deeper backwards : depth = %10d", deepNodeMembership -> depth);
          ${trace.token}    }
        }
        while (deepNodeMembership != shallowNodeMembership) {
          deepNodeMembership = deepNodeMembership -> parent;
          shallowNodeMembership = shallowNodeMembership -> parent;
          iEdgeCount ++;
          jEdgeCount ++;
          ${trace.token}    if (getTraceFlag(0) & ENSB_LOW_TRACE) {
          ${trace.token}      RF_nativePrint("\n  Distance : Parsing both backwards : depth = (%10d, %10d)", deepNodeMembership -> depth, shallowNodeMembership -> depth);
          ${trace.token}    }
        }
        if (iDepthCount > 0) { 
          ${trace.token}    if (getTraceFlag(0) & ENSB_LOW_TRACE) {
          ${trace.token}      if (iEdgeCount == 0) {
          ${trace.token}        if (jEdgeCount == 0) { 
          ${trace.token}          RF_nativePrint("\n  Distance : Concurrent node membership : pair = (%10d, %10d)", ii, jj);
          ${trace.token}        }
          ${trace.token}        else {
          ${trace.token}          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          ${trace.token}          RF_nativeError("\nRF-SRC:  Inconsistent edge counts in updateDistance().  Both must be zero.");
          ${trace.token}          RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
          ${trace.token}          RF_nativeExit();
          ${trace.token}        }
          ${trace.token}      }
          ${trace.token}    }
          realDistance = (double) (iEdgeCount + jEdgeCount) / (iDepthCount + jDepthCount);
        }
        else {
          realDistance = 0;
        }
        ${trace.token}    if (getTraceFlag(0) & ENSB_LOW_TRACE) {
        ${trace.token}      RF_nativePrint("\n  Distance : Real Increment: %10f \n", realDistance);
        ${trace.token}    }
        rfsrc_omp_atomic_update(&RF_distancePtr[ii][jj], realDistance);
      }
    }
  }
  else {
    utTermMembership =  RF_utTermMembership[b];
    utTermMembershipCount =  RF_utTermMembershipCount[b];
    for (uint i = 1; i <= membershipSize; i++) {
      Node *iNodeMembership, *jNodeMembership;
      Node *deepNodeMembership, *shallowNodeMembership;
      uint  iEdgeCount, jEdgeCount;
      uint  iDepthCount, jDepthCount;
      uint  *deepEdgeCount;
      double realDistance;
      double minDistance;
      uint ii, jj;
      ii = membershipIndex[i];
      for (uint j = 1; j <= i; j++) {
        jj = membershipIndex[j];
        rfsrc_omp_atomic_update(&RF_distanceDenPtr[ii][jj], 1.0);
        minDistance = 1.0;
        for (uint ki = 1; ki <= utTermMembershipCount[ii]; ki++) {
          iNodeMembership = RF_tTermList[b][utTermMembership[ii][ki]] -> mate;
          iDepthCount = iNodeMembership -> depth;
          for (uint kj = 1; kj <= utTermMembershipCount[jj]; kj++) {
            jNodeMembership = RF_tTermList[b][utTermMembership[jj][kj]] -> mate;
            jDepthCount = jNodeMembership -> depth;
            iEdgeCount = jEdgeCount = 0;
            if ((iNodeMembership -> depth) > (jNodeMembership -> depth)) {
              deepNodeMembership = iNodeMembership;
              shallowNodeMembership = jNodeMembership;
              deepEdgeCount = & iEdgeCount;
            }
            else {
              deepNodeMembership = jNodeMembership;
              shallowNodeMembership = iNodeMembership;
              deepEdgeCount = & jEdgeCount;
            }
            ${trace.token}    if (getTraceFlag(0) & ENSB_LOW_TRACE) {
            ${trace.token}      RF_nativePrint("\n  Distance : CaseID (i, j) = (%10d, %10d)", ii, jj);
            ${trace.token}      RF_nativePrint("\n  Distance : NodeID (i, j) = (%10d, %10d)", iNodeMembership -> nodeID, jNodeMembership -> nodeID);
            ${trace.token}      RF_nativePrint("\n  Distance : Depth  (i, j) = (%10d, %10d)", iNodeMembership -> depth, jNodeMembership -> depth);
            ${trace.token}    }
            while ((deepNodeMembership -> depth) > (shallowNodeMembership -> depth)) {
              deepNodeMembership = deepNodeMembership -> parent;
              (*deepEdgeCount) ++;
              ${trace.token}    if (getTraceFlag(0) & ENSB_LOW_TRACE) {
              ${trace.token}      RF_nativePrint("\n  Distance : Parsing deeper backwards : depth = %10d", deepNodeMembership -> depth);
              ${trace.token}    }
            }
            while (deepNodeMembership != shallowNodeMembership) {
              deepNodeMembership = deepNodeMembership -> parent;
              shallowNodeMembership = shallowNodeMembership -> parent;
              iEdgeCount ++;
              jEdgeCount ++;
              ${trace.token}    if (getTraceFlag(0) & ENSB_LOW_TRACE) {
              ${trace.token}      RF_nativePrint("\n  Distance : Parsing both backwards : depth = (%10d, %10d)", deepNodeMembership -> depth, shallowNodeMembership -> depth);
              ${trace.token}    }
            }
            if (iDepthCount > 0) { 
              ${trace.token}    if (getTraceFlag(0) & ENSB_LOW_TRACE) {
              ${trace.token}      if (iEdgeCount == 0) {
              ${trace.token}        if (jEdgeCount == 0) { 
              ${trace.token}          RF_nativePrint("\n  Distance : Concurrent node membership : pair = (%10d, %10d)", ii, jj);
              ${trace.token}        }
              ${trace.token}        else {
              ${trace.token}          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
              ${trace.token}          RF_nativeError("\nRF-SRC:  Inconsistent edge counts in updateDistance().  Both must be zero.");
              ${trace.token}          RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
              ${trace.token}          RF_nativeExit();
              ${trace.token}        }
              ${trace.token}      }
              ${trace.token}    }
              realDistance = (double) (iEdgeCount + jEdgeCount) / (iDepthCount + jDepthCount);
            }
            else {
              realDistance = 0;
            }
            ${trace.token}    if (getTraceFlag(0) & ENSB_LOW_TRACE) {
            ${trace.token}      RF_nativePrint("\n  Distance : Real Current : %10f \n", realDistance);
            ${trace.token}    }
            if (realDistance < minDistance) {
              minDistance = realDistance;
              ${trace.token}    if (getTraceFlag(0) & ENSB_LOW_TRACE) {
              ${trace.token}      RF_nativePrint("\n  Distance : Real Minimum Update: %10f \n", minDistance);
              ${trace.token}    }
              if (minDistance == 0) {
                goto distMarginal;
              }
            }
          }
        }
      distMarginal:
        continue;
        ${trace.token}    if (getTraceFlag(0) & ENSB_LOW_TRACE) {
        ${trace.token}      RF_nativePrint("\n  Distance : Minimum Distance: %10f \n", minDistance);
        ${trace.token}    }
        rfsrc_omp_atomic_update(&RF_distancePtr[ii][jj], minDistance);
      }
    }
  }
  ${trace.token}  if (getTraceFlag(0) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdateDistance() EXIT ...\n");
  ${trace.token}  }
}
void updateSplitDepth(uint treeID, Node *rootPtr, uint maxDepth) {
  Node  *parent;
  double *localSplitDepth;
  uint index;
  uint i, j, k;
  if (RF_tLeafCount[treeID] > 0) {
    index = 0;  
    if (RF_opt & (OPT_SPLDPTH_1 | OPT_SPLDPTH_2)) {
      if (RF_opt & OPT_SPLDPTH_1) {
        index = 1;
      }
      else {
        index = treeID;
      }
    }
    else {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Illegal updateSplitDepth() call.");
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
    localSplitDepth = dvector(1, RF_xSize);
    ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
    ${trace.token}    RF_nativePrint("\nTerminal Node Predictor Split Depths for Tree:  %10d", treeID);
    ${trace.token}    RF_nativePrint("\n      leaf      count     splits --> \n");
    ${trace.token}    LeafLinkedObj *leafLinkedPtr;
    ${trace.token}    leafLinkedPtr = RF_leafLinkedObjHead[treeID] -> fwdLink;
    ${trace.token}    while (leafLinkedPtr != NULL) {
    ${trace.token}      k = 0;
    ${trace.token}      parent = leafLinkedPtr -> nodePtr;
    ${trace.token}      if (parent != NULL) {
    ${trace.token}        k++;    
    ${trace.token}        RF_nativePrint("%10d %10d ", parent -> nodeID, parent -> depth);
    ${trace.token}        for (j = 1; j <= parent -> depth; j++) {
    ${trace.token}          RF_nativePrint("%10d ", (parent -> splitDepth)[j]);
    ${trace.token}        }
    ${trace.token}      }
    ${trace.token}      else {
    ${trace.token}        RF_nativePrint("%10d ", k);
    ${trace.token}      }
    ${trace.token}      RF_nativePrint("\n");
    ${trace.token}      leafLinkedPtr = leafLinkedPtr -> fwdLink;
    ${trace.token}    }
    ${trace.token}  }
    for (i = 1; i <= RF_observationSize; i++) {
      for (j = 1; j <= RF_xSize; j++) {
        localSplitDepth[j] = RF_nativeNaN;
      }
      parent = RF_tTermList[treeID][RF_tTermMembership[treeID][i] -> nodeID] -> mate;
      for (k = 1; k <= parent -> depth; k++) {
        if (RF_nativeIsNaN(localSplitDepth[(parent -> splitDepth)[k]])) {
          localSplitDepth[(parent -> splitDepth)[k]] = (double) k;
        }
      }
      for (j = 1; j <= RF_xSize; j++) {
        if (RF_nativeIsNaN(localSplitDepth[j])) {
          localSplitDepth[j] = (double) maxDepth + 1;
        }
      }
      if (RF_opt & OPT_SPLDPTH_1) {
#ifdef _OPENMP
#pragma omp critical (_update_splitdepth)
#endif
        {  
          for (j = 1; j <= RF_xSize; j++) {
            RF_splitDepthPtr[index][j][i] += localSplitDepth[j];
          }
        }
      }
      else {
        for (j = 1; j <= RF_xSize; j++) {
          RF_splitDepthPtr[index][j][i] += localSplitDepth[j];
        }
      }
    }
    free_dvector(localSplitDepth, 1, RF_xSize);
    freeSplitPath(treeID);
  }
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nRF-SRC:  SplitDepth update complete.");
  ${trace.token}  }
}
char pruneBranch(uint obsSize,
                 uint treeID,
                 Node **nodesAtDepth,
                 uint nadCount,
                 uint ptnTarget,
                 uint ptnCurrent) {
  char pruneFlag;
  uint i, j;
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\npruneBranch() ENTRY ...\n");
  ${trace.token}  }
  pruneFlag = TRUE;
  double *varianceAtDepth =  dvector(1, nadCount);
  uint   *vadSortedIndex  = uivector(1, nadCount);
  for (i = 1; i <= nadCount; i++) {
    varianceAtDepth[i] = nodesAtDepth[i] -> variance;
  }
  indexx(nadCount, varianceAtDepth, vadSortedIndex);
  j = nadCount;
  while ((j >= 1) && pruneFlag) {
    ${trace.token}  if (getTraceFlag(treeID) & SUMM_LOW_TRACE) {
    ${trace.token}    RF_nativePrint("\n  Node pseudo-terminated at index:  %10d", vadSortedIndex[j]);
    ${trace.token}  }
    nodesAtDepth[vadSortedIndex[j]] -> pseudoTerminal = TRUE;
    (nodesAtDepth[vadSortedIndex[j]] -> left)  -> pseudoTerminal = FALSE;
    (nodesAtDepth[vadSortedIndex[j]] -> right) -> pseudoTerminal = FALSE;
    for (i = 1; i <= obsSize; i++) {
      if ( (RF_pNodeMembership[treeID][i] == nodesAtDepth[vadSortedIndex[j]] -> left) ||
           (RF_pNodeMembership[treeID][i] == nodesAtDepth[vadSortedIndex[j]] -> right)) {
        RF_pNodeMembership[treeID][i] = nodesAtDepth[vadSortedIndex[j]];
      }
    }
    j --;
    ptnCurrent --;
    if (ptnCurrent <= ptnTarget) {
      pruneFlag = FALSE;
    }
  }
  free_dvector(varianceAtDepth, 1, nadCount);
  free_uivector(vadSortedIndex, 1, nadCount);
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\npruneBranch() EXIT ...\n");
  ${trace.token}  }
  return pruneFlag;
}
uint pruneTree(uint obsSize, uint treeID, uint ptnTarget) {
  Node **nodesAtDepth;
  uint   ptnCurrent;
  uint   nadCount;
  uint   tagDepth;
  char   pruneFlag;
  uint   i;
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\npruneTree() ENTRY ...\n");
  ${trace.token}  }
  if (ptnTarget < 1) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Illegal target PTN count in pruneTree():  %10d", ptnTarget);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  if (RF_tLeafCount[treeID] == 0) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Illegal call to pruneTree() on a rejected tree:  %10d", treeID);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  nodesAtDepth = (Node **) new_vvector(1, RF_tLeafCount[treeID], NRUTIL_NPTR);
  ptnCurrent = RF_tLeafCount[treeID];
  tagDepth = getMaximumDepth(RF_root[treeID]) - 1;
  pruneFlag = (ptnCurrent > ptnTarget) && (tagDepth > 0);
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\n  Prune flag:   %10d", pruneFlag);
  ${trace.token}    RF_nativePrint("\n  PTN Current:  %10d", ptnCurrent);
  ${trace.token}    RF_nativePrint("\n  PTN Target:   %10d", ptnTarget);
  ${trace.token}    RF_nativePrint("\n  Maximum depth for tree %10d:  %10d", treeID, tagDepth);
  ${trace.token}  }
  while (pruneFlag) {
    for (i = 1; i <= RF_tLeafCount[treeID]; i++) {
      nodesAtDepth[i] = NULL;
    }
    nadCount = 0;
    getNodesAtDepth(RF_root[treeID], tagDepth, nodesAtDepth, &nadCount);
    ${trace.token}  if (getTraceFlag(treeID) & SUMM_LOW_TRACE) {
    ${trace.token}    RF_nativePrint("\n  Nodes at depth %10d:  %10d", tagDepth, nadCount);
    ${trace.token}  }
    ${trace.token}  if (getTraceFlag(treeID) & SUMM_LOW_TRACE) {
    ${trace.token}    RF_nativePrint("\n  PTN Current:  %10d", ptnCurrent);
    ${trace.token}  }
    ${trace.token}  if (getTraceFlag(treeID) & SUMM_LOW_TRACE) {
    ${trace.token}    RF_nativePrint("\nDiagnostic Tree Output before pruning:  ");
    ${trace.token}    RF_nativePrint("\n    treeID     nodeID     parmID        depth     ptnFlg \n");
    ${trace.token}    printTreeInfo(treeID, RF_root[treeID]);
    ${trace.token}  }
    pruneFlag = pruneBranch(obsSize, treeID, nodesAtDepth, nadCount, ptnTarget, ptnCurrent);
    ${trace.token}  if (getTraceFlag(treeID) & SUMM_LOW_TRACE) {
    ${trace.token}    RF_nativePrint("\nDiagnostic Tree Output after pruning:  ");
    ${trace.token}    RF_nativePrint("\n    treeID     nodeID     parmID        depth     ptnFlg \n");
    ${trace.token}    printTreeInfo(treeID, RF_root[treeID]);
    ${trace.token}  }
    if(pruneFlag) {
      ptnCurrent -= nadCount;
      tagDepth --;
    }
    else {
      ptnCurrent = ptnTarget;
    }
  }
  free_new_vvector(nodesAtDepth, 1, RF_tLeafCount[treeID], NRUTIL_NPTR);
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\npruneTree() EXIT ...\n");
  ${trace.token}  }
  return ptnCurrent;
}
void stackAuxiliary(char mode, uint b) {
  uint obsSize;
  ${trace.token}  if (getTraceFlag(b) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAuxiliary() ENTRY ...\n");
  ${trace.token}  }
  RF_nodeMembership[b] = (Node **) new_vvector(1, RF_observationSize, NRUTIL_NPTR);
  RF_bootMembershipFlag[b] = cvector(1, RF_observationSize);
  RF_oobMembershipFlag[b] = cvector(1, RF_observationSize);
  RF_bootMembershipCount[b] = uivector(1, RF_observationSize);
  RF_ibgMembershipIndex[b] = uivector(1, RF_observationSize);
  RF_oobMembershipIndex[b] = uivector(1, RF_observationSize);
  RF_bootMembershipIndex[b] = uivector(1, RF_identityMembershipIndexSize);
  if (mode == RF_PRED) {
    RF_fnodeMembership[b] = (Node **) new_vvector(1, RF_fobservationSize, NRUTIL_NPTR);
  }
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    break;
  default:
    obsSize = RF_observationSize;
    break;
  }
  if (RF_optHigh & OPT_MEMB_PRUN) {
    RF_pNodeMembership[b] = (Node **) new_vvector(1, obsSize, NRUTIL_NPTR);
  }
  ${trace.token}  if (getTraceFlag(b) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAuxiliary() Exit ...\n");
  ${trace.token}  }
}
void unstackAuxiliary(char mode, uint b) {
  uint obsSize;
  ${trace.token}  if (getTraceFlag(b) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nunstackAuxiliary() ENTRY ...\n");
  ${trace.token}  }
  free_new_vvector(RF_nodeMembership[b], 1, RF_observationSize, NRUTIL_NPTR);
  free_cvector(RF_bootMembershipFlag[b], 1, RF_observationSize);
  free_cvector(RF_oobMembershipFlag[b], 1, RF_observationSize);
  free_uivector(RF_bootMembershipCount[b], 1, RF_observationSize);
  free_uivector(RF_ibgMembershipIndex[b], 1, RF_observationSize);
  free_uivector(RF_oobMembershipIndex[b], 1, RF_observationSize);
  free_uivector(RF_bootMembershipIndex[b], 1, RF_identityMembershipIndexSize);
  if (mode == RF_PRED) {
    free_new_vvector(RF_fnodeMembership[b],  1, RF_fobservationSize, NRUTIL_NPTR);
  }
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    break;
  default:
    obsSize = RF_observationSize;
    break;
  }
  if (RF_optHigh & OPT_MEMB_PRUN) {
    free_new_vvector(RF_pNodeMembership[b], 1, obsSize, NRUTIL_NPTR);
  }
  ${trace.token}  if (getTraceFlag(b) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nunstackAuxiliary() Exit ...\n");
  ${trace.token}  }
}
void printPseudoTNInfo(char mode, uint b) {
  uint i;
  RF_pNodeList[b] = (Node **) new_vvector(1, RF_pLeafCount[b] + 1, NRUTIL_NPTR);
  i = 0;
  getPTNodeList(RF_root[b], RF_pNodeList[b], &i);
  ${trace.token}     if (getTraceFlag(b) & SUMM_MED_TRACE) {
  ${trace.token}       uint obsSize = (mode == RF_PRED) ? RF_fobservationSize : RF_observationSize;
  ${trace.token}       RF_nativePrint("\nPseudo-Terminal Node Pointers:  %10d", b);
  ${trace.token}       RF_nativePrint("\n       index         leaf       nodeID\n");
  ${trace.token}       for (i = 1; i <=  RF_pLeafCount[b]; i++) {
  ${trace.token}         RF_nativePrint("%12d %12x %12d \n", i, RF_pNodeList[b][i], RF_pNodeList[b][i] -> nodeID);
  ${trace.token}       }
  ${trace.token}   
  ${trace.token}       RF_nativePrint("\nFinal Pseudo-Membership (all data):  %10d", b);
  ${trace.token}       RF_nativePrint("\n       index         node         leaf\n");
  ${trace.token}       for (i = 1; i <=  obsSize; i++) {
  ${trace.token}         RF_nativePrint("%12d %12x %12d \n", i, RF_pNodeMembership[b][i], RF_pNodeMembership[b][i] -> nodeID);
  ${trace.token}       }
  ${trace.token}     }
  free_new_vvector(RF_pNodeList[b], 1, RF_pLeafCount[b] + 1, NRUTIL_NPTR);
}      
void getPTNodeList(Node    *parent,
                   Node   **list,
                   uint    *offset) {
  if (!(parent -> pseudoTerminal)) {
    getPTNodeList(parent ->  left, list, offset);
    getPTNodeList(parent -> right, list, offset);
  }
  else {
    (*offset) ++;
    list[*offset] = parent;
  }
}
void getSplitPath(uint treeID, Node *parent) {
  Node *reversePtr;
  uint i;
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetSplitPath() ENTRY ...\n");
  ${trace.token}  }
  if (!(RF_opt & (OPT_SPLDPTH_1 | OPT_SPLDPTH_2))) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Call to calculate split depth without the option being active.");
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  if (parent -> depth > 0) {
    RF_maxDepth[treeID] = (parent -> depth > RF_maxDepth[treeID]) ? parent -> depth : RF_maxDepth[treeID];
    stackSplitDepth(parent, parent -> depth);
    reversePtr = parent;
    for (i = 1; i <= parent -> depth; i++) {
      if ((reversePtr -> parent) == NULL) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Reverse parsing of tree failed in restoreTree().");
        RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
        RF_nativeExit();
      }
      (parent -> splitDepth)[(parent -> depth) - i + 1] = ((reversePtr -> parent) -> splitInfo) -> randomVar[1];
      reversePtr = reversePtr -> parent;
    }
  }
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetSplitPath() EXIT ...\n");
  ${trace.token}  }
}
void freeSplitPath(uint treeID) {
  LeafLinkedObj *leafLinkedPtr;
  leafLinkedPtr = RF_leafLinkedObjHead[treeID] -> fwdLink;
  while (leafLinkedPtr != NULL) {
    unstackSplitDepth(leafLinkedPtr -> nodePtr);
    leafLinkedPtr = leafLinkedPtr -> fwdLink;
  }
}
uint getMaximumDepth(Node *parent) {
  uint result, rLeft, rRight;
  result = parent -> depth;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    rLeft = getMaximumDepth(parent ->  left);
    rRight = getMaximumDepth(parent -> right);
    result = (rLeft > rRight) ? rLeft : rRight;
  }
  return result;
}
void getNodesAtDepth(Node *parent, uint tagDepth, Node **nodesAtDepth, uint *nadCount) {
  char recurseFlag;
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetNodesAtDepth() ENTRY ...\n");
  ${trace.token}    }
  ${trace.token}  }
  recurseFlag = TRUE;
  if (tagDepth == parent -> depth) {
    if (parent -> splitInfo != NULL) {
      (*nadCount) ++;
      nodesAtDepth[*nadCount] = parent;
    }
    recurseFlag = FALSE;
  }
  else {
    if (((parent -> left) == NULL) && ((parent -> right) == NULL)) {
      recurseFlag = FALSE;
    }
  }
  if (recurseFlag) {
    getNodesAtDepth(parent ->  left, tagDepth, nodesAtDepth, nadCount);
    getNodesAtDepth(parent -> right, tagDepth, nodesAtDepth, nadCount);
  }
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetNodesAtDepth() EXIT ...\n");
  ${trace.token}    }
  ${trace.token}  }
}
