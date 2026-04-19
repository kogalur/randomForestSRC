
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "importanceAnti.h"
#include "polarity.h"
#include "splitInfo.h"
#include "splitGreedy.h"
#include "rfsrcUtil.h"
#include "nodeOps.h"
#include "nrutil.h"
${trace.token} #include "error.h"
void getAntiMembership (char       mode,
                        uint       treeID,
                        Terminal **vimpMembership,
                        uint       p) {
  Node    *rootPtr;
  uint    *membershipIndex;
  uint     membershipSize;
  double **xArray;
  ${trace.token}  if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetAntiMembership() ENTRY.");
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nType is VIMP_ANTI.");
  ${trace.token}    }
  ${trace.token}  }
  rootPtr = RF_root[treeID];
  switch (mode) {
  case RF_PRED:
    membershipSize = RF_fobservationSize;
    membershipIndex = RF_fidentityMembershipIndex;
    xArray = RF_fobservation[treeID];
    break;
  default:
    membershipSize  = RF_oobSize[treeID];
    membershipIndex = RF_oobMembershipIndex[treeID];
    xArray = RF_observation[treeID];
    break;
  }
  for (uint i = 1; i <= membershipSize; i++) {
    uint ii;
    ii = membershipIndex[i];
    vimpMembership[ii] = antiMembership(treeID, rootPtr, ii, p, xArray) -> mate;
  }
  ${trace.token}  if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetAntiMembership() EXIT.");
  ${trace.token}    }
  ${trace.token}  }
}
Node *antiMembershipGeneric(uint     treeID,
                            Node    *parent,
                            uint     individual,                            
                            uint     vimpX,
                            double **xArray) {
  char daughterFlag;
  char antiSplitFlag;
  Node *result;
  SplitInfo *info;
  double alpha;
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nantiMembershipGeneric() ENTRY... \n");
  ${trace.token}    }
  ${trace.token}  }
  result = parent;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    info = parent -> splitInfo;
    antiSplitFlag = FALSE;
    if (vimpX > 0) {
      if ((uint) info -> randomVar[1] == vimpX) {
        antiSplitFlag = TRUE;
      }
    }
    else {
      if(RF_importanceFlag[info -> randomVar[1]] == TRUE) {
        antiSplitFlag = TRUE;
      }
    }
    daughterFlag = getDaughterPolarity(0, info, individual, xArray);
    if(antiSplitFlag == TRUE) {
      alpha = ran1D(treeID);
      if (alpha <= RF_vimpThreshold) {
        if (daughterFlag == LEFT) {
          daughterFlag = RIGHT;
        }
        else {
          daughterFlag = LEFT;
        }
      }
      else {
      }
    }  
    else {
      ${trace.token}      if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
      ${trace.token}        if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
      ${trace.token}          if(daughterFlag == LEFT) {
      ${trace.token}            RF_nativePrint("\nAnti Faithful Daughter LEFT :  indv = %10d, vimpX =  %10d, nodeID = %10d", individual, vimpX, (parent -> left) -> nodeID);
      ${trace.token}          }
      ${trace.token}          else {
      ${trace.token}            RF_nativePrint("\nAnti Faithful Daughter RGHT :  indv = %10d, vimpX =  %10d, nodeID = %10d", individual, vimpX, (parent -> right) -> nodeID);
      ${trace.token}          }
      ${trace.token}        }
      ${trace.token}      }
    }  
    if (daughterFlag == LEFT) {
      result = antiMembershipGeneric(treeID, parent ->  left, individual, vimpX, xArray);
    }
    else {
      result = antiMembershipGeneric(treeID, parent -> right, individual, vimpX, xArray);
    }
  }
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nantiMembershipGeneric() EXIT... \n");
  ${trace.token}    }
  ${trace.token}  }
  return result;
}
Node *antiMembershipJIT(uint     treeID,
                        Node    *root,
                        uint     individual,
                        uint     vimpX,
                        double **xArray) {
  ulong rootIndex, nodeAbsIndex;
  uint  rmbrAbsOffset, ambrAbsOffset;
  uint rmbrDummyIter, ambrDummyIter;
  ulong *rootMWCPoffset, *nodeAbsMWCPoffset;
  SplitInfo *info;
  char parseFlag;
  char daughterFlag;
  char antiSplitFlag;
  double alpha;
  uint adj, i, k;
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nantiMembershipJIT() ENTRY... \n");
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
          ${trace.token}        RF_nativePrint("\n  antiMembershipJIT() (h-idx, (parent -> splitInfo) -> randomVar[k])   = (%10d, %10d)", k, RF_parmID_[k][nodeAbsIndex]); 
          ${trace.token}      }
          info -> randomVar[k] = RF_parmID_[k][nodeAbsIndex];
          info -> mwcpSizeAbs[k] = RF_mwcpSZ_[k][nodeAbsIndex];
          ${trace.token}      if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
          ${trace.token}        RF_nativePrint("\n  antiMembershipJIT() (h-idx, (parent -> splitInfo) -> mwcpSizeAbs[k]) = (%10d, %10d)", k, RF_mwcpSZ_[k][nodeAbsIndex]);
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
        ${trace.token}      RF_nativePrint("\n antiMembershipJIT() parent -> splitInfo:  %20x", info);
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
      antiSplitFlag = FALSE;
      if (vimpX > 0) {
        if ((uint) info -> randomVar[1] == vimpX) {
          antiSplitFlag = TRUE;
        }
      }
      else {
        if(RF_importanceFlag[info -> randomVar[1]] == TRUE) {
          antiSplitFlag = TRUE;
        }
      }
      daughterFlag = getDaughterPolarity(0, info, individual, xArray);
      if(antiSplitFlag == TRUE) {
        alpha = ran1D(treeID);
        if (alpha <= RF_vimpThreshold) {
          if (daughterFlag == LEFT) {
            daughterFlag = RIGHT;
          }
          else {
            daughterFlag = LEFT;
          }
          ${trace.token}      if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
          ${trace.token}        if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
          ${trace.token}          if(daughterFlag == LEFT) {
          ${trace.token}            RF_nativePrint("\nAnti Reversed Daughter LEFT :  indv = %10d, vimpX =  %10d", individual, vimpX);
          ${trace.token}          }
          ${trace.token}          else {
          ${trace.token}            RF_nativePrint("\nAnti Reversed Daughter RGHT :  indv = %10d, vimpX =  %10d", individual, vimpX);
          ${trace.token}          }
          ${trace.token}        }
          ${trace.token}      }
        }
        else {
          ${trace.token}      if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
          ${trace.token}        if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
          ${trace.token}          if(daughterFlag == LEFT) {
          ${trace.token}            RF_nativePrint("\nAnti (Alpha) Faithful Daughter LEFT :  indv = %10d, vimpX =  %10d", individual, vimpX);
          ${trace.token}          }
          ${trace.token}          else {
          ${trace.token}            RF_nativePrint("\nAnti (Alpha) Faithful Daughter RGHT :  indv = %10d, vimpX =  %10d", individual, vimpX);
          ${trace.token}          }
          ${trace.token}        }
          ${trace.token}      }
        }
      }  
      else {
        ${trace.token}      if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
        ${trace.token}        if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
        ${trace.token}          if(daughterFlag == LEFT) {
        ${trace.token}            RF_nativePrint("\nAnti Faithful Daughter LEFT :  indv = %10d, vimpX =  %10d", individual, vimpX);
        ${trace.token}          }
        ${trace.token}          else {
        ${trace.token}            RF_nativePrint("\nAnti Faithful Daughter RGHT :  indv = %10d, vimpX =  %10d", individual, vimpX);
        ${trace.token}          }
        ${trace.token}        }
        ${trace.token}      }
      }  
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
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nantiMembershipJIT() EXIT... \n");
  ${trace.token}    }
  ${trace.token}  }
  return parent;
}
