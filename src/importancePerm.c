
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "importancePerm.h"
#include "polarity.h"
#include "splitInfo.h"
#include "splitGreedy.h"
#include "rfsrcUtil.h"
#include "nodeOps.h"
#include "nrutil.h"
${trace.token} #include "error.h"
void getPermuteMembership (char       mode,
                           uint       treeID,
                           Terminal **vimpMembership,
                           uint       p) {
  Node    *rootPtr;
  uint     obsSize;
  uint    *membershipIndex;
  uint     membershipSize;
  double **predictorPtr;
  ${trace.token}  if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetPermuteMembership() ENTRY.");
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nType is VIMP_PERM.");
  ${trace.token}    }
  ${trace.token}  }
  rootPtr = RF_root[treeID];
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    membershipSize = RF_fobservationSize;
    membershipIndex = RF_fidentityMembershipIndex;
    predictorPtr = RF_fobservation[treeID];
    break;
  default:
    obsSize = RF_observationSize;
    membershipSize  = RF_oobSize[treeID];
    membershipIndex = RF_oobMembershipIndex[treeID];
    predictorPtr = RF_observation[treeID];
    break;
  }
  double **shadowVIMP = (double **) new_vvector(1, RF_xSize, NRUTIL_DPTR);
  for (uint j = 1; j <= RF_xSize; j++) {
    shadowVIMP[j] = predictorPtr[j];
  }
  if (RF_opt & OPT_VIMP_JOIN) {
    for (uint pp = 1; pp <= RF_intrPredictorSize; pp++) {
      uint *permuteVIMP = uivector(1, membershipSize + 1);
      uint targetCov = RF_intrPredictor[pp];
      shadowVIMP[targetCov] = dvector(1, obsSize);
      ${trace.token}  if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
      ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
      ${trace.token}      RF_nativePrint("\nRandom seed for ran1D():  %20d", randomGetChainVimp(treeID));
      ${trace.token}    }
      ${trace.token}  }
      permute(4, treeID, membershipSize, permuteVIMP);
      for (uint i = 1; i <= membershipSize; i++) {
        uint ii = membershipIndex[i];
        shadowVIMP[targetCov][ii] = predictorPtr[targetCov][membershipIndex[permuteVIMP[i]]];
      }
      free_uivector(permuteVIMP, 1, membershipSize + 1);
    }
    for (uint i = 1; i <= membershipSize; i++) {
      uint ii = membershipIndex[i];
      vimpMembership[ii] = getMembership(treeID, rootPtr, ii, shadowVIMP) -> mate;
    }
    for (uint pp = 1; pp <= RF_intrPredictorSize; pp++) {
      uint targetCov = RF_intrPredictor[pp];
      free_dvector(shadowVIMP[targetCov], 1, obsSize);
    }
  }
  else {
    uint *permuteVIMP = uivector(1, membershipSize + 1);
    uint targetCov = p;
    shadowVIMP[targetCov] = dvector(1, obsSize);
    ${trace.token}  if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
    ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
    ${trace.token}      RF_nativePrint("\nRandom seed for ran1D():  %20d", randomGetChainVimp(treeID));
    ${trace.token}    }
    ${trace.token}  }
    permute(4, treeID, membershipSize, permuteVIMP);
    for (uint i = 1; i <= membershipSize; i++) {
      uint ii = membershipIndex[i];
      shadowVIMP[targetCov][ii] = predictorPtr[targetCov][membershipIndex[permuteVIMP[i]]];
    }
    for (uint i = 1; i <= membershipSize; i++) {
      uint ii = membershipIndex[i];
      vimpMembership[ii] = getMembership(treeID, rootPtr, ii, shadowVIMP) -> mate;
    }
    free_uivector(permuteVIMP, 1, membershipSize + 1);
    free_dvector(shadowVIMP[targetCov], 1, obsSize);
  }
  free_new_vvector(shadowVIMP, 1, RF_xSize, NRUTIL_DPTR);
  ${trace.token}  if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetPermuteMembership() EXIT.");
  ${trace.token}    }
  ${trace.token}  }
}
Node *getMembershipGeneric(uint     treeID,
                           Node    *parent,
                           uint     individual,
                           double **xArray) {
  char daughterFlag;
  Node *result = parent;
  SplitInfo *info;
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetMembershipGeneric() ENTRY:  idx= %10d: ", individual);
  ${trace.token}    }
  ${trace.token}  }
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    info = parent -> splitInfo;
    daughterFlag = getDaughterPolarity(0, info, individual, xArray);
    ${trace.token}    if (getTraceFlag(0) & VIMP_LOW_TRACE) {
    ${trace.token}      if (getTraceFlag(0) & !TURN_OFF_TRACE) {
    ${trace.token}          if(daughterFlag == LEFT) {
    ${trace.token}            RF_nativePrint("\nNon-Greedy Daughter LEFT :  %10d ", individual);
    ${trace.token}          }
    ${trace.token}          else {
    ${trace.token}            RF_nativePrint("\nNon-Greedy Daughter RGHT :  %10d ", individual);
    ${trace.token}          }
    ${trace.token}        }
    ${trace.token}      }
    if (daughterFlag == LEFT) {
      result = getMembershipGeneric(treeID, parent ->  left, individual, xArray);
    }
    else {
      result = getMembershipGeneric(treeID, parent -> right, individual, xArray);
    }
  }
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetMembershipGeneric() EXIT:   idx= %10d: ", individual);
  ${trace.token}    }
  ${trace.token}  }
  return result;
}
Node *getMembershipJIT(uint     treeID,
                         Node    *root,
                         uint     individual,
                         double **xArray) {
  ulong rootIndex, nodeAbsIndex;
  uint  rmbrAbsOffset, ambrAbsOffset;
  uint rmbrDummyIter, ambrDummyIter;
  ulong *rootMWCPoffset, *nodeAbsMWCPoffset;
  SplitInfo *info;
  char parseFlag;
  char daughterFlag;
  uint adj, i, k;
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetMembershipJIT() ENTRY:  idx= %10d: ", individual);
  ${trace.token}    }
  ${trace.token}  }
  rmbrDummyIter = ambrDummyIter = 0;
  nodeAbsIndex = rootIndex = RF_restoreTreeOffset[treeID];
  rmbrAbsOffset = ambrAbsOffset = 0;
  rootMWCPoffset    = ulvector(1, 1);
  nodeAbsMWCPoffset = ulvector(1, 1);
  rootMWCPoffset[1] = RF_restoreMWCPoffset[1][treeID];
  parseFlag = TRUE;
  Node *parent = root;
  while (parseFlag) {
    ${trace.token}      if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
    ${trace.token}        RF_nativePrint("\n JIT parent and absolute offset:  %20x %10d", parent, nodeAbsIndex); 
    ${trace.token}      }
    if (parent -> nodeID == 0) {
      ${trace.token}      if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
      ${trace.token}        RF_nativePrint("\n JIT parent is not initialized:  (parent = %20x) -> nodeID %10d", parent, parent -> nodeID); 
      ${trace.token}      }
      parent -> nodeID = RF_nodeID_[nodeAbsIndex];
      parent -> bnodeID = nodeAbsIndex - rootIndex + 1;
      if (RF_parmID_[1][nodeAbsIndex] != 0) {
        parent -> blnodeID = parent -> bnodeID + 1;
        parent -> brnodeID = RF_brnodeID_[nodeAbsIndex];
        info = parent -> splitInfo = makeSplitInfo(0);
        adj = 1;
        info -> mwcpSizeAbs = uivector(1, adj);
        info -> randomVar   = ivector(1, adj);
        info -> randomPts   = new_vvector(1, adj, NRUTIL_VPTR);
        for (k = 1; k <= adj; k++) {
          ${trace.token}      if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
          ${trace.token}        RF_nativePrint("\n  getMembershipJIT() (h-idx, (parent -> splitInfo) -> randomVar[k])   = (%10d, %10d)", k, RF_parmID_[k][nodeAbsIndex]); 
          ${trace.token}      }
          info -> randomVar[k] = RF_parmID_[k][nodeAbsIndex];
          info -> mwcpSizeAbs[k] = RF_mwcpSZ_[k][nodeAbsIndex];
          ${trace.token}      if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
          ${trace.token}        RF_nativePrint("\n  getMembershipJIT() (h-idx, (parent -> splitInfo) -> mwcpSizeAbs[k]) = (%10d, %10d)", k, RF_mwcpSZ_[k][nodeAbsIndex]);
          ${trace.token}      }
          if (RF_mwcpSZ_[k][nodeAbsIndex] > 0) {
            nodeAbsMWCPoffset[k] = rootMWCPoffset[k] + RF_fsrecID_[k][nodeAbsIndex];
            ${trace.token}      if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
            ${trace.token}        RF_nativePrint("\n (h-idx, absolute mwcp offset) -> (%10d %10d)", k, nodeAbsMWCPoffset[k]);
            ${trace.token}      }
            ${trace.token}      if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
            ${trace.token}        RF_nativePrint("\n  mwcpPT (reversed):  [%20x] ", RF_mwcpPT_[k]);
            ${trace.token}      }
            info -> randomPts[k] = uivector(1, RF_mwcpSZ_[k][nodeAbsIndex]);
            for (i = 1; i <= RF_mwcpSZ_[k][nodeAbsIndex]; i++) {
              ((uint *) info -> randomPts[k])[i] = RF_mwcpPT_[k][nodeAbsMWCPoffset[k]];
              ${trace.token}        if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
              ${trace.token}          RF_nativePrint("%8x ", RF_mwcpPT_[k][nodeAbsMWCPoffset[k]]);
              ${trace.token}        }
              nodeAbsMWCPoffset[k] ++;
            }
          }
          else {
            info -> randomPts[k] = dvector(1, 1);
            ((double *) info -> randomPts[k])[1] =  RF_contPT_[k][nodeAbsIndex];
          }
        }
        ${trace.token}  if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
        ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
        ${trace.token}      RF_nativePrint("\n getMembershipJIT() parent -> splitInfo:  %20x", info);
        ${trace.token}      getSplitObjectInfo(parent -> splitInfo);
        ${trace.token}    }
        ${trace.token}  }
      }
      else {
        info = parent -> splitInfo = NULL;
      }
    }  
    else {
      info = parent -> splitInfo;
      ${trace.token}      if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
      ${trace.token}        RF_nativePrint("\n JIT parent is initialized:  (parent = %20x) -> nodeID %10d", parent, parent -> nodeID); 
      ${trace.token}      }
    }
    if (info != NULL) {
      daughterFlag = getDaughterPolarity(0, info, individual, xArray);
      ${trace.token}      if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
      ${trace.token}        if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
      ${trace.token}          if(daughterFlag == LEFT) {
      ${trace.token}            RF_nativePrint("\nJIT Membership Daughter LEFT :  indv = %10d", individual);
      ${trace.token}          }
      ${trace.token}          else {
      ${trace.token}            RF_nativePrint("\nJIT Membership Daughter RGHT :  indv = %10d", individual);
      ${trace.token}          }
      ${trace.token}        }
      ${trace.token}      }
      if (daughterFlag == LEFT) {
        nodeAbsIndex = nodeAbsIndex + 1;
        if (parent -> left == NULL) {
          parent -> left  = makeNode(0);
        }
        setParent(parent -> left, parent);
        parent = parent -> left;
      }
      else {
        nodeAbsIndex = rootIndex + (parent -> brnodeID) - 1;
        i = parent -> nodeID;
        while(i < RF_nodeID_[nodeAbsIndex]) {
          rmbrAbsOffset = rmbrAbsOffset + RF_TN_RCNT_ptr[treeID][i];
          ambrAbsOffset = ambrAbsOffset + RF_TN_ACNT_ptr[treeID][i];
          i++;
        }
        if (parent -> right == NULL) {
          parent -> right = makeNode(0);
        }
        setParent(parent -> right, parent);
        parent = parent -> right;
      }
      ${trace.token}            if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
      ${trace.token}              RF_nativePrint("\nJIT running RMBR and AMBR offsets:  %10d and %10d", rmbrAbsOffset, ambrAbsOffset);
      ${trace.token}            }
    }
    else {
      if (RF_tTermList[treeID][parent -> nodeID] != NULL) {
        ${trace.token}            if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
        ${trace.token}              RF_nativePrint("\nJIT terminal node exists:  %20x -> %10d", RF_tTermList[treeID][parent -> nodeID], RF_tTermList[treeID][parent -> nodeID] -> nodeID);
        ${trace.token}            }
      }
      else {
        RF_leafLinkedObjTail[treeID] = makeAndSpliceLeafLinkedObj(RF_leafLinkedObjTail[treeID],
                                                                  parent,
                                                                  RF_TN_RCNT_ptr[treeID][parent -> nodeID],
                                                                  RF_TN_ACNT_ptr[treeID][parent -> nodeID]);
        RF_tTermList[treeID][parent -> nodeID] = RF_leafLinkedObjTail[treeID] -> termPtr;
        ${trace.token}            if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
        ${trace.token}              RF_nativePrint("\nJIT terminal node does not exist (creating it):  %20x -> %10d", RF_tTermList[treeID][parent -> nodeID], RF_tTermList[treeID][parent -> nodeID] -> nodeID);
        ${trace.token}            }
        updateTerminalNodeOutcomes(RF_PRED,
                                   treeID,
                                   parent -> mate,
                                   & RF_RMBR_ID_ptr[treeID][rmbrAbsOffset],
                                   RF_TN_RCNT_ptr[treeID][parent -> nodeID],
                                   & RF_AMBR_ID_ptr[treeID][ambrAbsOffset],
                                   RF_TN_ACNT_ptr[treeID][parent -> nodeID],
                                   & rmbrDummyIter,
                                   & ambrDummyIter);
      }
      ${trace.token}            if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
      ${trace.token}              RF_nativePrint("\nJIT final RMBR and AMBR offsets:  %10d and %10d", rmbrAbsOffset, ambrAbsOffset);
      ${trace.token}              RF_nativePrint("\nJIT final RCNT and ACNT sizes:    %10d and %10d", RF_TN_RCNT_ptr[treeID][parent -> nodeID], RF_TN_ACNT_ptr[treeID][parent -> nodeID]);
      ${trace.token}            }
      parseFlag = FALSE;
    }
  }  
  free_ulvector(rootMWCPoffset, 1, 1);
  free_ulvector(nodeAbsMWCPoffset, 1, 1);
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetMembershipJIT() EXIT:   idx= %10d: ", individual);
  ${trace.token}    }
  ${trace.token}  }
    return parent;
}
void permute(uint ranGenID, uint parallelID, uint n, uint *indx) {
  uint i,j,k;
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\npermute() ENTRY... \n");
  ${trace.token}    }
  ${trace.token}  }
  for (i=1; i<= n; i++) {
    indx[i] = 0;
  }
  for (i=n; i > 0; i--) {
    k = (uint) ceil(ran1D(parallelID)*(i*1.0));
    for (j = 1; k > 0; j++) {
      if (indx[j] == 0) {
        k--;
      }
    }
    indx[j-1] = i;
  }
  ${trace.token}  if (getTraceFlag(0) & OUTP_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nRequest for Permutation:  \n");
  ${trace.token}    RF_nativePrint("     index    permute\n");
  ${trace.token}    for (i=1; i <= n; i++) {
  ${trace.token}      RF_nativePrint("%10d %10d \n", i, indx[i]);
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\npermute() EXIT... \n");
  ${trace.token}    }
  ${trace.token}  }
}
