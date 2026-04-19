
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "treeUtil.h"
#include "impute.h"
#include "splitGreedy.h"
#include "bootstrap.h"
#include "rfsrcUtil.h"
#include "tree.h"
#include "nrutil.h"
#include "error.h"
char growTreeRecursive (uint     r,
                        char     rootFlag,
                        char     multImpFlag,
                        uint     treeID,
                        Node    *parent,
                        uint    *bootMembrIndxIter,
                        uint    *rmbrIterator,
                        uint    *ambrIterator) {
  char  bootResult;
  char  splitResult;
  char  forkResult;
  char leftResult, rghtResult;
  char terminalFlag;
  uint *bootMembrIndx;
  uint bootMembrSize;
  SplitInfoMax *splitInfoMax;
  SplitInfo *splitInfo;
  uint *repMembrIndx, *allMembrIndx;
  uint  repMembrSize,  allMembrSize;
  uint i, p;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\ngrowTreeRecursive(%10d) ENTRY ...\n", treeID);
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\ngrowTreeRecursive (%10d) leaf:  %10d at depth %10d \n", treeID, parent -> nodeID, parent -> depth);
  ${trace.token}    RF_nativePrint("\n  called with   rep size:  %10d", parent -> repMembrSize);
  ${trace.token}    RF_nativePrint("\n  called with   all size:  %10d", parent -> allMembrSize);
  ${trace.token}  }
  bootResult = TRUE;
  terminalFlag = TRUE;
  splitInfo = NULL;
  allMembrIndx = parent -> allMembrIndx;
  allMembrSize = parent -> allMembrSize;
  repMembrIndx = parent -> repMembrIndx;
  repMembrSize = parent -> repMembrSize;
  if (rootFlag) {
    if ( (!(RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
         ( (RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ) {
      bootMembrIndx  = uivector(1, RF_bootstrapSize);
      bootMembrSize = RF_bootstrapSize;
    }
    else {
      bootMembrIndx  = uivector(1, allMembrSize);
      bootMembrSize = allMembrSize;
    }
    bootResult = bootstrap (RF_GROW,
                            treeID,
                            parent,
                            allMembrIndx,
                            allMembrSize,
                            bootMembrIndx,
                            bootMembrSize);
    repMembrIndx = parent -> repMembrIndx = bootMembrIndx;
    repMembrSize = parent -> repMembrSize = parent -> repMembrSizeAlloc = bootMembrSize;
    if (rootFlag & bootResult) {
      if (RF_vtry > 0) {
        if (RF_vtryMode == RF_VTRY_NULL) {
          for (p = 1; p <= RF_xSize; p++) {
            if (RF_vtryArray[treeID][p] > 0) {
              (parent -> permissible)[p] = FALSE;
            }
          }
        }
        else {
          if (RF_vtryMode == RF_VTRY_HOLD) {
            for (p = 1; p <= RF_xSize; p++) {
              if (RF_vtryArray[treeID][p] > 0) {
                (parent -> permissible)[p] = FALSE;
              }
            }
          }
          else {
          }
        }
      }
      if (RF_mRecordSize > 0) {
        for (p = 1; p <= RF_mpIndexSize; p++) {
          if (RF_mpIndex[p] > 0) {
            if (parent -> mpSign[p] == -1) {
              (parent -> permissible)[RF_mpIndex[p]] = FALSE;
            }
          }
        }
      }  
      parent -> permissibleIndxSize = 0;
      for (p = 1; p <= RF_xSize; p++) {
        if ((parent -> permissible)[p] == TRUE) {
          parent -> permissibleIndx[++ (parent -> permissibleIndxSize)] = p;
        }
      }
      parent -> permissibleReIndxFlag = FALSE;
      parent -> permissibleOwnershipFlag = TRUE;
    }  
  }  
  else {
    bootMembrIndx = repMembrIndx;
    bootMembrSize = repMembrSize;
    parent -> mpSign = (parent -> parent) -> mpSign;
  }
  if (bootResult) {
    if (multImpFlag == FALSE) {
      if (RF_mRecordSize > 0) {
        imputeNode(RF_GROW,
                   FALSE,  
                   TRUE,   
                   treeID,
                   parent,
                   repMembrIndx,
                   repMembrSize,
                   allMembrIndx,
                   allMembrSize);
        if (RF_timeIndex > 0) {
          if (RF_mTimeFlag == TRUE) {
            updateTimeIndexArray(treeID,
                                 allMembrIndx,
                                 allMembrSize,
                                 RF_time[treeID],
                                 FALSE,
                                 FALSE,
                                 RF_masterTimeIndex[treeID]);
          }
        }
      }
    }  
  }  
  if (bootResult) {
    if (rootFlag) {
      RF_tLeafCount[treeID] = 1;
    }
    splitInfoMax = makeSplitInfoMax(0);
    splitResult = getBestSplit(treeID,
                               parent,
                               RF_splitRule,
                               splitInfoMax,
                               multImpFlag);
    if (splitResult == TRUE) {
      splitInfo = makeSplitInfo(0);
      if (splitInfoMax -> size > 0) {
        if (splitInfoMax -> indicator != NULL) {
          splitInfo -> size = splitInfoMax -> size;
          splitInfo -> indicator = splitInfoMax -> indicator;
          splitInfoMax -> indicator = NULL;
          splitInfoMax -> size      = 0;
        }
      }
      else {
        splitInfo -> size = 0;
        splitInfo -> indicator = NULL;
      }
      splitInfo -> mwcpSizeAbs = uivector(1, 1);
      splitInfo -> randomVar   = ivector(1, 1);
      splitInfo -> randomPts   = new_vvector(1, 1, NRUTIL_VPTR);
      (splitInfo -> mwcpSizeAbs)[1] = splitInfoMax -> splitValueMaxFactSize;
      (splitInfo -> randomVar)[1] = splitInfoMax -> splitParameterMax;
      if ((splitInfo -> mwcpSizeAbs)[1] > 0) {
        (splitInfo -> randomPts)[1] = uivector(1, (splitInfo -> mwcpSizeAbs)[1]);
        for (i = 1; i <= (splitInfo -> mwcpSizeAbs)[1]; i++) {
          ((uint *) (splitInfo -> randomPts)[1])[i] = (splitInfoMax -> splitValueMaxFactPtr)[i];
        }
      }
      else {
        (splitInfo -> randomPts)[1] = dvector(1, 1);
        ((double *) (splitInfo -> randomPts)[1])[1] = splitInfoMax -> splitValueMaxCont;
      }
      freeSplitInfoMax(splitInfoMax);
      terminalFlag = FALSE;
      forkResult = forkAndUpdate(treeID,
                                 parent,
                                 repMembrIndx,
                                 repMembrSize,
                                 allMembrIndx,
                                 allMembrSize,
                                 multImpFlag,
                                 splitInfo,
                                 &RF_tLeafCount[treeID],
                                 RF_nodeMembership[treeID]);
      if (forkResult == TRUE) {
        ${trace.token}          if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
        ${trace.token}          if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
        ${trace.token}            RF_nativePrint("\ngNode Membership (%10d):  \n", treeID);
        ${trace.token}            for (i=1; i <=  RF_observationSize; i++) {
        ${trace.token}              RF_nativePrint("%10d %10d \n", i, RF_nodeMembership[treeID][i] -> nodeID);
        ${trace.token}            }
        ${trace.token}          }
        ${trace.token}          }
        ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
        ${trace.token}          RF_nativePrint("\ngrowTreeRecursive(%10d) LEFT:  \n", treeID);
        ${trace.token}        }
        leftResult = growTreeRecursive (r,
                               FALSE,
                               multImpFlag,
                               treeID,
                               parent -> left,
                               bootMembrIndxIter,
                               rmbrIterator,
                               ambrIterator);
        if(!leftResult) {
          ${trace.token}        if (getTraceFlag(treeID) & BOOT_MED_TRACE) {
          ${trace.token}          RF_nativePrint("\ngrowTreeRecursive(%10d) LEFT call failed:  ", treeID);
          ${trace.token}        }
        }
        ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
        ${trace.token}          RF_nativePrint("\ngrowTreeRecursive(%10d) RIGHT:  \n", treeID);
        ${trace.token}        }
        rghtResult = growTreeRecursive (r,
                               FALSE,
                               multImpFlag,
                               treeID,
                               parent -> right,
                               bootMembrIndxIter,
                               rmbrIterator,
                               ambrIterator);
        if(!rghtResult) {
          ${trace.token}        if (getTraceFlag(treeID) & BOOT_MED_TRACE) {
          ${trace.token}          RF_nativePrint("\ngrowTreeRecursive(%10d) RGHT call failed:  ", treeID);
          ${trace.token}        }
        }
        free_uivector((parent -> left)  -> allMembrIndx, 1, (parent -> left)  -> allMembrSizeAlloc);
        free_uivector((parent -> right) -> allMembrIndx, 1, (parent -> right) -> allMembrSizeAlloc);
        (parent -> left) -> allMembrIndx = (parent -> right) -> allMembrIndx = NULL;
        (parent -> left) -> allMembrSize = (parent -> right) -> allMembrSize = 0;
        (parent -> left) -> allMembrSizeAlloc = (parent -> right) -> allMembrSizeAlloc = 0;
        free_uivector((parent -> left)  -> repMembrIndx, 1, (parent -> left)  -> repMembrSizeAlloc);
        free_uivector((parent -> right) -> repMembrIndx, 1, (parent -> right) -> repMembrSizeAlloc);
        (parent -> left) -> repMembrIndx = (parent -> right) -> repMembrIndx = NULL;
        (parent -> left) -> repMembrSizeAlloc = (parent -> right) -> repMembrSizeAlloc = 0;
      }
      else {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  forkAndUpdate(%10d) failed.", treeID);
        RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
        RF_nativeExit();
      }
    }  
    else {
      parent -> splitFlag = FALSE;
      free_cvector(parent -> permissible, 1, parent -> xSize);
      free_uivector(parent -> permissibleIndx, 1, parent -> xSize);
      parent -> permissible = NULL;
      parent -> permissibleIndx = NULL;
      parent -> permissibleIndxSize = 0;
      parent -> splitInfo = NULL;
      freeSplitInfoMax(splitInfoMax);
      ${trace.token}      if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
      ${trace.token}        RF_nativePrint("\ngetBestSplit() FAILED ...\n");
      ${trace.token}      }
    }
  }  
  else {
    if (rootFlag) {
      if (!bootResult) {
        terminalFlag = FALSE;
        ${trace.token}      if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
        ${trace.token}        RF_nativePrint("\nbootstrap(%10d) FAILED at root node.  Tree rejected.  \n", treeID);
        ${trace.token}      }
      }
    }
  }  
  if (terminalFlag) {
    ${trace.token}  if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
    ${trace.token}    RF_nativePrint("\n Terminal node encountered.  Saving pointer for leaf:  %10d %20x ",  parent -> nodeID, parent);
    ${trace.token}  }
    RF_leafLinkedObjTail[treeID] = makeAndSpliceLeafLinkedObj(RF_leafLinkedObjTail[treeID],
                                                              parent,
                                                              repMembrSize,
                                                              allMembrSize);
    parent -> pseudoTerminal = TRUE;
    if (RF_opt & OPT_MISS_OUT) {
      imputeNodeAndSummarize(r,
                             RF_GROW,
                             treeID,
                             parent,
                             repMembrIndx,
                             repMembrSize,
                             allMembrIndx,
                             allMembrSize,
                             NULL,
                             0);
    }
    if (r == RF_nImpute) {
      if (RF_optHigh & OPT_MEMB_USER) {
        for (i = 1; i <= allMembrSize; i++) {
          RF_MEMB_ID_ptr[treeID][allMembrIndx[i]] = parent -> nodeID;
        }
      }
      if (RF_optHigh & OPT_MEMB_OUTG) {
        Terminal *termPtr;
        uint iter;
        termPtr = RF_leafLinkedObjTail[treeID] -> termPtr;
        termPtr -> oobMembrSizeAlloc = allMembrSize;
        termPtr -> oobMembrIndx = uivector(1, termPtr -> oobMembrSizeAlloc);      
        iter = 0;
        for (i = 1; i <= allMembrSize; i++) {
          if (RF_bootMembershipFlag[treeID][allMembrIndx[i]] == FALSE) {
            termPtr -> oobMembrIndx[++iter] = allMembrIndx[i];
          }
        }
        termPtr -> oobMembrSize = iter;
        termPtr -> ibgMembrSizeAlloc = allMembrSize;
        termPtr -> ibgMembrIndx = uivector(1, termPtr -> ibgMembrSizeAlloc);      
        iter = 0;
        for (i = 1; i <= allMembrSize; i++) {
          if (RF_bootMembershipFlag[treeID][allMembrIndx[i]] == TRUE) {
            termPtr -> ibgMembrIndx[++iter] = allMembrIndx[i];
          }
        }
        termPtr -> ibgMembrSize = iter;
      }
      updateTerminalNodeOutcomes(RF_GROW,
                                 treeID,
                                 RF_leafLinkedObjTail[treeID] -> termPtr,
                                 repMembrIndx,
                                 repMembrSize,
                                 allMembrIndx,
                                 allMembrSize,
                                 rmbrIterator,
                                 ambrIterator);
      if (RF_opt & (OPT_SPLDPTH_1 | OPT_SPLDPTH_2)) {
        getSplitPath(treeID, parent);
      }
    }
    else {
      initTerminalNodeMembership(treeID,
                                 RF_leafLinkedObjTail[treeID] -> termPtr,
                                 allMembrIndx,
                                 allMembrSize);
    }
  }  
  if (rootFlag) {
    if ( (!(RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
         ( (RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ) {
      free_uivector(bootMembrIndx, 1, RF_bootstrapSize);
    }
    else {
      free_uivector(bootMembrIndx, 1, allMembrSize);
    }
    parent -> repMembrIndx = NULL;
    parent -> repMembrSizeAlloc = 0;
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\ngrowTreeRecursive(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
  return bootResult;
}
void freeTree(uint treeID, Node *parent) {
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfreeTree() ENTRY ...\n");
  ${trace.token}  }
  if (parent != NULL) {
    if ((parent -> left) != NULL) {
      freeTree(treeID, parent -> left);
    }
    if ((parent -> right) != NULL) {
      freeTree(treeID, parent -> right);
    }
    freeNode(parent);
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfreeTree() EXIT ...\n");
  ${trace.token}  }
}
void saveStatistics(char     mode,
                    uint     b,
                    Node    *parent,
                    uint    *offset,
                    double  *spltST,
                    uint    *dpthST) {
  ${trace.token}  if (getTraceFlag(b) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nsaveStatistics() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(b) & SPLT_HGH_TRACE) {
  ${trace.token}    getNodeInfo(parent);
  ${trace.token}  }
  spltST[++(*offset)] = parent -> splitStatistic;
  dpthST[(*offset)]   = parent -> depth;
  ${trace.token}  if (getTraceFlag(b) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\n   topology      depth               spltST \n");
  ${trace.token}    RF_nativePrint(" %10d %10d %20.8f ", *offset, dpthST[*offset], spltST[*offset]);
  ${trace.token}  }
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    saveStatistics(mode, b, parent ->  left, offset, spltST, dpthST);
    saveStatistics(mode, b, parent -> right, offset, spltST, dpthST);
  }
  ${trace.token}  if (getTraceFlag(b) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nsaveStatistics() EXIT ...\n");
  ${trace.token}  }
}
void initTerminalNodeMembership(uint       treeID,
                                Terminal  *parent,
                                uint      *allMembrIndx,
                                uint       allMembrSize) {
  uint i;
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\ninitTerminalNodeMembership() Entry ...\n");
  ${trace.token}  }
  for (i = 1; i <= allMembrSize; i++) {
    RF_tTermMembership[treeID][allMembrIndx[i]] = parent;
  }
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\ninitTerminalNodeMembership() Exit ...\n");
  ${trace.token}  }
}
void updatePruning(char mode, uint treeID) {
  Terminal ***gTermMembership;
  uint        obsSize;
  uint i;
  ${trace.token}  if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdatePruning() ENTRY ...\n");
  ${trace.token}  }
  if (RF_optHigh & OPT_MEMB_PRUN) {
    switch (mode) {
    case RF_PRED:
      obsSize = RF_fobservationSize;
      gTermMembership = RF_ftTermMembership;
      break;
    default:
      obsSize = RF_observationSize;
      gTermMembership = RF_tTermMembership;
      break;
    }
    ${trace.token}      if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
    ${trace.token}        RF_nativePrint("\nUpdate of PTNM with final TNM: ");
    ${trace.token}        RF_nativePrint("\n      index              address     nodeID\n");
    ${trace.token}      }
    for (i = 1; i <= obsSize; i++) {
      RF_pNodeMembership[treeID][i] = gTermMembership[treeID][i] -> mate;
      ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
      ${trace.token}    RF_nativePrint("\n %10d %20x %10d", i, RF_pNodeMembership[treeID][i], RF_pNodeMembership[treeID][i] -> nodeID);
      ${trace.token}  }
    }
    RF_pLeafCount[treeID] = pruneTree(obsSize, treeID, RF_ptnCount);
    for (i=1; i <= obsSize; i++) {
      RF_PRUN_ID_ptr[treeID][i] = RF_pNodeMembership[treeID][i] -> nodeID;
    }
    ${trace.token}    if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
    ${trace.token}      printPseudoTNInfo(mode, treeID);
    ${trace.token}    }
  }
  ${trace.token}  if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdatePruning() EXIT ...\n");
  ${trace.token}  }
}
void updateCaseDepth(char mode, uint treeID) {
  Terminal ***gTermMembership;
  uint        obsSize;
  uint i;
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdateCaseDepth() ENTRY ...\n");
  ${trace.token}  }
  if (RF_opt & OPT_CASE_DPTH) {
    switch (mode) {
    case RF_PRED:
      obsSize = RF_fobservationSize;
      gTermMembership = RF_ftTermMembership;
      break;
    default:
      obsSize = RF_observationSize;
      gTermMembership = RF_tTermMembership;
      break;
    }
    ${trace.token}      if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
    ${trace.token}        RF_nativePrint("\nUpdate of CASE_DPTH for treeID %10d: ", treeID);
    ${trace.token}        RF_nativePrint("\n      index              address     nodeID     depth\n");
    ${trace.token}      }
    for (i = 1; i <= obsSize; i++) {
      RF_CASE_DPTH_ptr[treeID][i] = gTermMembership[treeID][i] -> mate -> depth;
      ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
      ${trace.token}    RF_nativePrint("\n %10d %20x %10d %10d", i, gTermMembership[treeID][i], gTermMembership[treeID][i] -> nodeID, gTermMembership[treeID][i] -> mate -> depth);
      ${trace.token}  }
    }
  }
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdateCaseDepth() EXIT ...\n");
  ${trace.token}  }
}
