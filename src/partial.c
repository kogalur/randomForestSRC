
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "partial.h"
#include "splitGreedy.h"
#include "factorOps.h"
#include "rfsrcUtil.h"
#include "nodeOps.h"
#include "nrutil.h"
${trace.token} #include "error.h"
void getAndUpdatePartialMembership(uint treeID, Node *root) {
  uint i, j;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetAndUpdatePartialMembership() ENTRY ...\n");
  ${trace.token}  }
  Terminal **membership =  (Terminal **) new_vvector(1, RF_observationSize, NRUTIL_TPTR);
  if (!(RF_optHigh & OPT_JIT_TOP)) {
    for (i = 1; i <= RF_partialLength; i++) {
      partialMembershipGeneric(treeID,
                               root,
                               i,
                               RF_identityMembershipIndex,
                               RF_observationSize,
                               RF_observation[treeID],
                               membership);
      ${trace.token}    if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
      ${trace.token}      RF_nativePrint("\nFinal PARTIAL Membership (all data):  (treeID = %10d, partialID = %10d)", treeID, i);
      ${trace.token}      RF_nativePrint("\n       index         leaf ->\n");
      ${trace.token}      for (uint m=1; m <= RF_observationSize; m++) {
      ${trace.token}        RF_nativePrint("\n %12d %12d", m, membership[m] -> nodeID);
      ${trace.token}      }
      ${trace.token}    }
      updatePartialCalculations(treeID, i, membership);
    }
  }
  else {
    for (i = 1; i <= RF_partialLength; i++) {
      for (j = 1; j <= RF_observationSize; j++) {
        partialMembershipJIT(treeID,
                             root,
                             i,
                             NULL,
                             RF_identityMembershipIndex[j], 
                             RF_observation[treeID],
                             membership);
      }
      ${trace.token}    if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
      ${trace.token}      RF_nativePrint("\nFinal PARTIAL Membership (all data):  (treeID = %10d, partialID = %10d)", treeID, i);
      ${trace.token}      RF_nativePrint("\n       index         leaf ->\n");
      ${trace.token}      for (uint m=1; m <= RF_observationSize; m++) {
      ${trace.token}        RF_nativePrint("\n %12d %12d", m, membership[m] -> nodeID);
      ${trace.token}      }
      ${trace.token}    }
      updatePartialCalculations(treeID, i, membership);
    }
  }
  free_new_vvector(membership, 1, RF_observationSize, NRUTIL_TPTR);
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetAndUpdatePartialMembership() EXIT ...\n");
  ${trace.token}  }
}
void partialMembershipGeneric(uint       treeID,
                              Node      *parent,
                              uint       partialIndex,
                              uint      *allMembrIndx,
                              uint       allMembrSize,
                              double   **xArray,
                              Terminal **membership) {
  char terminalFlag;
  uint *leftAllMembrIndx;
  uint *rghtAllMembrIndx;
  uint leftAllMembrSize;
  uint rghtAllMembrSize;
  uint jLeft;
  uint jRght;
  uint obsSize;
  uint primaryPartialIndex, secondaryPartialIndex;
  uint factorValue;
  double continuousValue;
  char daughterFlag;
  SplitInfo *info;
  uint i, k;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\npartialMembershipGeneric() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\npartialMembershipGeneric (%10d) leaf:  %10d \n", treeID, parent -> nodeID);
  ${trace.token}    RF_nativePrint("\n  called with (size, partIdx, partVal):  (%10d %10d %10.4f)",
  ${trace.token}      allMembrSize, RF_partialXvar, RF_partialValue[partialIndex]);
  ${trace.token}  }
  terminalFlag = TRUE;
  obsSize = RF_observationSize;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    info = parent -> splitInfo;
    terminalFlag = FALSE;
    ${trace.token}    if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
    ${trace.token}      RF_nativePrint("\nForking On:  ");
    ${trace.token}      getNodeInfo(parent);
    ${trace.token}    }
    leftAllMembrIndx = rghtAllMembrIndx = NULL;
    leftAllMembrSize = rghtAllMembrSize = 0;
    uint *indicator = uivector(1, obsSize);
    leftAllMembrSize = rghtAllMembrSize = 0;
    primaryPartialIndex = secondaryPartialIndex = k = 0;
    if ((uint) info -> randomVar[1] ==  RF_partialXvar) {
      primaryPartialIndex = RF_partialXvar;
    }
    else {
      for (k = 1; k <= RF_partialLength2; k++) {
        if ((uint) info -> randomVar[1] ==  RF_partialXvar2[k]) {
          secondaryPartialIndex = k;
        }
      }
    }
    for (i = 1; i <= allMembrSize; i++) {
      ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
      ${trace.token}            RF_nativePrint("\npPartial Node Membership:  %10d %10d ", i, allMembrIndx[i]);
      ${trace.token}        }
      if (info -> mwcpSizeAbs[1] > 0) {
        if (primaryPartialIndex > 0) { 
          factorValue = (uint) RF_partialValue[partialIndex];
          ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
          ${trace.token}          RF_nativePrint("  --> primary partial split on %10d", factorValue);
          ${trace.token}        }
        }
        else if (secondaryPartialIndex > 0) {
          factorValue = (uint) RF_partialValue2[secondaryPartialIndex];            
          ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
          ${trace.token}          RF_nativePrint("  --> secondary partial split on %10d", factorValue);
          ${trace.token}        }
        }
        else {
          factorValue = (uint) xArray[info -> randomVar[1]][allMembrIndx[i]];
          ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
          ${trace.token}          RF_nativePrint("  --> indigenous split on %10d", factorValue); 
          ${trace.token}        }
        }
        ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
        ${trace.token}        if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
        ${trace.token}          RF_nativePrint("\nNon-Greedy Daughter Comparison (value, const):  (%10d, ", factorValue);
        ${trace.token}          for (uint m = 1; m <= info -> mwcpSizeAbs[1]; m++) {
        ${trace.token}            RF_nativePrint(" %10x", ((uint*) info -> randomPts[1])[m]);
        ${trace.token}          }
        ${trace.token}          RF_nativePrint(")");
        ${trace.token}        }
        ${trace.token}      }
        daughterFlag = splitOnFactor(factorValue, (uint*) info -> randomPts[1]);
      }
      else {
        if (primaryPartialIndex > 0) { 
          continuousValue = RF_partialValue[partialIndex];
          ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
          ${trace.token}          RF_nativePrint("  --> primary partial split on %10.4f", continuousValue);
          ${trace.token}        }
        }
        else if (secondaryPartialIndex > 0) {
          continuousValue = RF_partialValue2[secondaryPartialIndex];            
          ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
          ${trace.token}          RF_nativePrint("  --> secondary partial split on %10.4f", continuousValue);
          ${trace.token}        }
        }
        else {
          continuousValue = xArray[info -> randomVar[1]][allMembrIndx[i]];
          ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
          ${trace.token}          RF_nativePrint("  --> indigenous split on %10.4f", continuousValue); 
          ${trace.token}        }
        }
        ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
        ${trace.token}        if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
        ${trace.token}          RF_nativePrint("\nNon-Greedy Daughter Comparison (value, const):  (%10.4f, %10.4f)", continuousValue, ((double*) info -> randomPts[1])[1]);
        ${trace.token}        }
        ${trace.token}      }
        daughterFlag =  (( ((double*) info -> randomPts[1])[1] - continuousValue) >= 0.0) ? LEFT : RIGHT;
      }
      ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
      ${trace.token}        if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
      ${trace.token}          if(daughterFlag == LEFT) {
      ${trace.token}            RF_nativePrint("\nNon-Greedy Daughter LEFT :  %10d %10d", i, allMembrIndx[i]);
      ${trace.token}          }
      ${trace.token}          else {
      ${trace.token}            RF_nativePrint("\nNon-Greedy Daughter RGHT :  %10d %10d", i, allMembrIndx[i]);
      ${trace.token}          }
      ${trace.token}        }
      ${trace.token}      }
      indicator[allMembrIndx[i]] = daughterFlag;
      if (daughterFlag == LEFT) {
        leftAllMembrSize ++;
        ${trace.token}          if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
        ${trace.token}            RF_nativePrint(" --> LEFT ");
        ${trace.token}          }
      }
      else {
        rghtAllMembrSize ++;
        ${trace.token}          if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
        ${trace.token}            RF_nativePrint(" --> RGHT ");
        ${trace.token}          }
      }
    }  
    ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
    ${trace.token}          RF_nativePrint("\n");
    ${trace.token}        }
    leftAllMembrIndx  = uivector(1, leftAllMembrSize + 1);
    rghtAllMembrIndx  = uivector(1, rghtAllMembrSize + 1);
    jLeft = jRght = 0;
    for (i = 1; i <= allMembrSize; i++) {
      if (indicator[allMembrIndx[i]] == LEFT) {
        leftAllMembrIndx[++jLeft] = allMembrIndx[i];
      }
      else {
        rghtAllMembrIndx[++jRght] = allMembrIndx[i];
      }
    }
    free_uivector(indicator, 1, obsSize);
    ${trace.token}        if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
    ${trace.token}          RF_nativePrint("\npartialMembershipGeneric(%10d) LEFT:  \n", treeID);
    ${trace.token}        }
    partialMembershipGeneric(treeID,
                          parent -> left,                             
                          partialIndex,
                          leftAllMembrIndx,
                          leftAllMembrSize,
                          xArray,
                          membership);
    ${trace.token}        if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
    ${trace.token}          RF_nativePrint("\npartialMembershipGeneric(%10d) RGHT:  \n", treeID);
    ${trace.token}        }
    partialMembershipGeneric(treeID,
                          parent -> right,
                          partialIndex,
                          rghtAllMembrIndx,
                          rghtAllMembrSize,
                          xArray,
                          membership);
    free_uivector(leftAllMembrIndx, 1, leftAllMembrSize + 1);
    free_uivector(rghtAllMembrIndx, 1, rghtAllMembrSize + 1);
  }  
  else {
    ${trace.token}  if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
    ${trace.token}    RF_nativePrint("\n  Terminal node encountered for leaf:  %10d %20x ",  parent -> nodeID, parent);
    ${trace.token}  }
  }
  if (terminalFlag) {
    for (i = 1; i <= allMembrSize; i++) {
      membership[allMembrIndx[i]] = parent -> mate;
    }
    ${trace.token}            if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
    ${trace.token}              RF_nativePrint("\nPartial Split Terminal Node Membership: ");
    ${trace.token}              RF_nativePrint("\n      index       indv membership ");
    ${trace.token}              for (i = 1; i <= allMembrSize; i++) {
    ${trace.token}                RF_nativePrint("\n %10d %10d ", i, allMembrIndx[i]) ;
    ${trace.token}              }
    ${trace.token}            }
  }  
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\npartialMembershipGeneric() EXIT ...\n");
  ${trace.token}  }
}
void partialMembershipJIT(uint       treeID,
                          Node      *root,
                          uint       partialIndex,
                          uint       *nullMembrIndx,
                          uint       individual,
                          double   **xArray,
                          Terminal **membership) {
  ulong rootIndex, nodeAbsIndex;
  uint  rmbrAbsOffset, ambrAbsOffset;
  uint rmbrDummyIter, ambrDummyIter;
  ulong *rootMWCPoffset, *nodeAbsMWCPoffset;
  SplitInfo *info;
  char parseFlag;
  char daughterFlag;
  uint adj, i, k;
  uint primaryPartialIndex, secondaryPartialIndex;
  uint factorValue;
  double continuousValue;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\npartialMembershipJIT() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\npartialMembershipJIT (%10d) \n", treeID);
  ${trace.token}    RF_nativePrint("\n  called with (indv, partIdx, partVal):  (%10d %10d %10.4f)",
  ${trace.token}      individual, RF_partialXvar, RF_partialValue[partialIndex]);
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
      primaryPartialIndex = secondaryPartialIndex = k = 0;
      if ((uint) info -> randomVar[1] ==  RF_partialXvar) {
        primaryPartialIndex = RF_partialXvar;
      }
      else {
        for (k = 1; k <= RF_partialLength2; k++) {
          if ((uint) info -> randomVar[1] ==  RF_partialXvar2[k]) {
            secondaryPartialIndex = k;
          }
        }
      }
      ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
      ${trace.token}            RF_nativePrint("\nPartial Node Membership:  %10d ", individual);
      ${trace.token}        }
      if (info -> mwcpSizeAbs[1] > 0) {
        if (primaryPartialIndex > 0) { 
          factorValue = (uint) RF_partialValue[partialIndex];
          ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
          ${trace.token}          RF_nativePrint("  --> primary partial split on %10d", factorValue);
          ${trace.token}        }
        }
        else if (secondaryPartialIndex > 0) {
          factorValue = (uint) RF_partialValue2[secondaryPartialIndex];            
          ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
          ${trace.token}          RF_nativePrint("  --> secondary partial split on %10d", factorValue);
          ${trace.token}        }
        }
        else {
          factorValue = (uint) xArray[info -> randomVar[1]][individual];
          ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
          ${trace.token}          RF_nativePrint("  --> indigenous split on %10d", factorValue); 
          ${trace.token}        }
        }
        ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
        ${trace.token}        if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
        ${trace.token}          RF_nativePrint("\nNon-Greedy Daughter Comparison (value, const):  (%10d, ", factorValue);
        ${trace.token}          for (uint m = 1; m <= info -> mwcpSizeAbs[1]; m++) {
        ${trace.token}            RF_nativePrint(" %10x", ((uint*) info -> randomPts[1])[m]);
        ${trace.token}          }
        ${trace.token}          RF_nativePrint(")");
        ${trace.token}        }
        ${trace.token}      }
        daughterFlag = splitOnFactor(factorValue, (uint*) info -> randomPts[1]);
      }
      else {
        if (primaryPartialIndex > 0) { 
          continuousValue = RF_partialValue[partialIndex];
          ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
          ${trace.token}          RF_nativePrint("  --> primary partial split on %10.4f", continuousValue);
          ${trace.token}        }
        }
        else if (secondaryPartialIndex > 0) {
          continuousValue = RF_partialValue2[secondaryPartialIndex];            
          ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
          ${trace.token}          RF_nativePrint("  --> secondary partial split on %10.4f", continuousValue);
          ${trace.token}        }
        }
        else {
          continuousValue = xArray[info -> randomVar[1]][individual];
          ${trace.token}        if (getTraceFlag(treeID) & FORK_DEF_TRACE) {
          ${trace.token}          RF_nativePrint("  --> indigenous split on %10.4f", continuousValue); 
          ${trace.token}        }
        }
        ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
        ${trace.token}        if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
        ${trace.token}          RF_nativePrint("\nNon-Greedy Daughter Comparison (value, const):  (%10.4f, %10.4f)", continuousValue, ((double*) info -> randomPts[1])[1]);
        ${trace.token}        }
        ${trace.token}      }
        daughterFlag =  (( ((double*) info -> randomPts[1])[1] - continuousValue) >= 0.0) ? LEFT : RIGHT;
      }
      ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
      ${trace.token}        if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
      ${trace.token}          if(daughterFlag == LEFT) {
      ${trace.token}            RF_nativePrint(" --> LEFT");
      ${trace.token}          }
      ${trace.token}          else {
      ${trace.token}            RF_nativePrint(" --> RGHT");
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
        ${trace.token}            if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
        ${trace.token}              RF_nativePrint("\nJIT running RMBR and AMBR offsets:  %10d and %10d", rmbrAbsOffset, ambrAbsOffset);
        ${trace.token}            }
      }
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
  membership[individual] = parent -> mate;
  free_ulvector(rootMWCPoffset, 1, 1);
  free_ulvector(nodeAbsMWCPoffset, 1, 1);
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\npartialMembershipJIT() EXIT ...\n");
  ${trace.token}  }
}
void updatePartialCalculations (uint       treeID,
                                uint       pVarIdx,
                                Terminal **partialMembership) {
  Terminal *terminalNode;
  uint  *membershipIndex;
  uint   membershipSize;
  uint   i, j, k;
  uint   ii;
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdatePartialCalculations() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & OUTP_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\n  pVarIdx:  %10d ", pVarIdx);
  ${trace.token}  }
  if (RF_tLeafCount[treeID] > 0) {
  if (RF_opt & OPT_OENS) {
    membershipSize  = RF_oobSize[treeID];
    membershipIndex = RF_oobMembershipIndex[treeID];
  }
  else {
    membershipSize  = RF_observationSize;
    membershipIndex = RF_identityMembershipIndex;
  }
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (RF_eventTypeSize > 1) {
      if (RF_partialType == RF_PART_YRLS) {
        for (i = 1; i <= membershipSize; i++) {
          ii = membershipIndex[i];
          terminalNode = partialMembership[ii];
#ifdef _OPENMP
          omp_set_lock(&(RF_lockPartial[ii]));
#endif
          for (j = 1; j <= RF_eventTypeSize; j++) {
            RF_partSURVptr[pVarIdx][j][1][ii] += terminalNode -> mortality[j];
          }
#ifdef _OPENMP
          omp_unset_lock(&(RF_lockPartial[ii]));
#endif
        }
      }
      else if (RF_partialType == RF_PART_CIFN) {
        for (i = 1; i <= membershipSize; i++) {
          ii = membershipIndex[i];
          terminalNode = partialMembership[ii];
#ifdef _OPENMP
          omp_set_lock(&(RF_lockPartial[ii]));
#endif
          for (j = 1; j <= RF_eventTypeSize; j++) {
            for (k = 1; k <= RF_partialTimeLength; k++) {
              RF_partSURVptr[pVarIdx][j][k][ii] += terminalNode -> CIF[j][k];
            }
          }
#ifdef _OPENMP
          omp_unset_lock(&(RF_lockPartial[ii]));
#endif
        }
      }
      else if (RF_partialType == RF_PART_CHFN) {
        for (i = 1; i <= membershipSize; i++) {
          ii = membershipIndex[i];
          terminalNode = partialMembership[ii];
#ifdef _OPENMP
          omp_set_lock(&(RF_lockPartial[ii]));
#endif
          for (j = 1; j <= RF_eventTypeSize; j++) {
            for (k = 1; k <= RF_partialTimeLength; k++) {
              RF_partSURVptr[pVarIdx][j][k][ii] += terminalNode -> CSH[j][k];
            }
          }
#ifdef _OPENMP
          omp_unset_lock(&(RF_lockPartial[ii]));
#endif
        }
      }
    }   
    else {
      if (RF_partialType == RF_PART_MORT) {
        for (i = 1; i <= membershipSize; i++) {
          ii = membershipIndex[i];
          terminalNode = partialMembership[ii];
#ifdef _OPENMP
          omp_set_lock(&(RF_lockPartial[ii]));
#endif
            RF_partSURVptr[pVarIdx][1][1][ii] += terminalNode -> mortality[1];
#ifdef _OPENMP
            omp_unset_lock(&(RF_lockPartial[ii]));
#endif
        }
      }
      else if (RF_partialType == RF_PART_NLSN) {
        for (i = 1; i <= membershipSize; i++) {
          ii = membershipIndex[i];
          terminalNode = partialMembership[ii];
#ifdef _OPENMP
          omp_set_lock(&(RF_lockPartial[ii]));
#endif
          for (k = 1; k <= RF_partialTimeLength; k++) {
            RF_partSURVptr[pVarIdx][1][k][ii] += terminalNode -> nelsonAalen[k];
          }
#ifdef _OPENMP
          omp_unset_lock(&(RF_lockPartial[ii]));
#endif
        }
      }
      else if (RF_partialType == RF_PART_SURV) {
        for (i = 1; i <= membershipSize; i++) {
          ii = membershipIndex[i];
          terminalNode = partialMembership[ii];
#ifdef _OPENMP
          omp_set_lock(&(RF_lockPartial[ii]));
#endif
          for (k = 1; k <= RF_partialTimeLength; k++) {
            RF_partSURVptr[pVarIdx][1][k][ii] += terminalNode -> survival[k];
          }
#ifdef _OPENMP
          omp_unset_lock(&(RF_lockPartial[ii]));
#endif
        }
      }
    }
  }
  else {
    if (RF_rTargetFactorCount > 0) {
      for (i = 1; i <= membershipSize; i++) {
        ii = membershipIndex[i];
        terminalNode = partialMembership[ii];
#ifdef _OPENMP
        omp_set_lock(&(RF_lockPartial[ii]));
#endif
        for (j = 1; j <= RF_rTargetFactorCount; j++) {
          for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
            RF_partCLASptr[pVarIdx][j][k+1][ii] += (double) (terminalNode -> multiClassProb)[RF_rFactorMap[RF_rTargetFactor[j]]][k] / (double) (terminalNode -> membrCount);
          }
        }
#ifdef _OPENMP
        omp_unset_lock(&(RF_lockPartial[ii]));
#endif
      }
    }
    if (RF_rTargetNonFactorCount > 0) {
      for (i = 1; i <= membershipSize; i++) {
        ii = membershipIndex[i];
        terminalNode = partialMembership[ii];
#ifdef _OPENMP
        omp_set_lock(&(RF_lockPartial[ii]));
#endif
        for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
          RF_partREGRptr[pVarIdx][j][ii] += (terminalNode -> meanResponse)[RF_rNonFactorMap[RF_rTargetNonFactor[j]]];
        }
#ifdef _OPENMP
        omp_unset_lock(&(RF_lockPartial[ii]));
#endif
      }
    }
  }
  ${trace.token}    if (getTraceFlag(treeID) & OUTP_DEF_TRACE) {
  ${trace.token}      RF_nativePrint("\nPartial outcome calculation:  treeID %10d \n", treeID);
  ${trace.token}      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
  ${trace.token}        if (RF_eventTypeSize > 1) {
  ${trace.token}          if (RF_partialType == RF_PART_YRLS) {
  ${trace.token}            RF_nativePrint("\nRF_PART_YRLS:  \n");
  ${trace.token}            for (j=1; j <= RF_eventTypeSize; j++) {
  ${trace.token}              RF_nativePrint("\n  [pVarIdx][event] = [%10d][%10d] \n", pVarIdx, j);
  ${trace.token}              RF_nativePrint("          ");
  ${trace.token}              for (i = 1; i <= RF_observationSize; i++) {
  ${trace.token}                RF_nativePrint("%10d", i);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}              for (i=1; i <= RF_observationSize; i++) {
  ${trace.token}                RF_nativePrint("%10.4f", RF_partSURVptr[pVarIdx][j][1][i]);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}            }
  ${trace.token}          }
  ${trace.token}          else if (RF_partialType == RF_PART_CIFN) {
  ${trace.token}            RF_nativePrint("\nRF_PART_CIFN:  \n");
  ${trace.token}            for (i = 1; i <= RF_observationSize; i++) {
  ${trace.token}              RF_nativePrint("Indv %10d: \n", i);
  ${trace.token}              RF_nativePrint("                 time ");
  ${trace.token}              for (j = 1; j <= RF_eventTypeSize; j++) {
  ${trace.token}                RF_nativePrint("%10d ", j);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}              for (k=1; k <= RF_partialTimeLength; k++) {
  ${trace.token}                RF_nativePrint("%10d %10.4f ", k, RF_partialTime[k]);
  ${trace.token}                for (j = 1; j <= RF_eventTypeSize; j++) {
  ${trace.token}                  RF_nativePrint("%10.4f ", RF_partSURVptr[pVarIdx][j][k][i]);
  ${trace.token}                }
  ${trace.token}                RF_nativePrint("\n");
  ${trace.token}              }
  ${trace.token}            }
  ${trace.token}          }
  ${trace.token}          else if (RF_partialType == RF_PART_CHFN) {
  ${trace.token}            RF_nativePrint("\nRF_PART_CHFN:  \n");
  ${trace.token}            for (i = 1; i <= RF_observationSize; i++) {
  ${trace.token}              RF_nativePrint("Indv %10d: \n", i);
  ${trace.token}              RF_nativePrint("                 time ");
  ${trace.token}              for (j = 1; j <= RF_eventTypeSize; j++) {
  ${trace.token}                RF_nativePrint("%10d ", j);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}              for (k=1; k <= RF_partialTimeLength; k++) {
  ${trace.token}                RF_nativePrint("%10d %10.4f ", k, RF_partialTime[k]);
  ${trace.token}                for (j = 1; j <= RF_eventTypeSize; j++) {
  ${trace.token}                  RF_nativePrint("%10.4f ", RF_partSURVptr[pVarIdx][j][k][i]);
  ${trace.token}                }
  ${trace.token}                RF_nativePrint("\n");
  ${trace.token}              }
  ${trace.token}            }
  ${trace.token}          }
  ${trace.token}        }
  ${trace.token}        else {
  ${trace.token}          if (RF_partialType == RF_PART_MORT) {
  ${trace.token}            RF_nativePrint("\nRF_PART_MORT:  \n");
  ${trace.token}            for (i = 1; i <= RF_observationSize; i++) {
  ${trace.token}              RF_nativePrint("%10d %10.4f", i, RF_partSURVptr[pVarIdx][1][1][i]);
  ${trace.token}            }
  ${trace.token}            RF_nativePrint("\n");
  ${trace.token}           
  ${trace.token}          }
  ${trace.token}          else if (RF_partialType == RF_PART_NLSN) {
  ${trace.token}            RF_nativePrint("\nRF_PART_NLSN:  \n");
  ${trace.token}            RF_nativePrint("      index       time        obs -> ");  
  ${trace.token}            for (k = 1; k <= RF_partialTimeLength; k++) {
  ${trace.token}              RF_nativePrint("\n %10d %10.4f ", k, RF_partialTime[k]);
  ${trace.token}              for (i = 1; i <= RF_observationSize; i++) {
  ${trace.token}                RF_nativePrint(" %10.4f", RF_partSURVptr[pVarIdx][1][k][i]);
  ${trace.token}              }
  ${trace.token}            }
  ${trace.token}          }
  ${trace.token}          else if (RF_partialType == RF_PART_SURV) {
  ${trace.token}            RF_nativePrint("\nRF_PART_SURV:  \n");
  ${trace.token}            RF_nativePrint("      index       time        obs -> ");  
  ${trace.token}            for (k = 1; k <= RF_partialTimeLength; k++) {
  ${trace.token}              RF_nativePrint("\n %10d %10.4f ", k, RF_partialTime[k]);
  ${trace.token}              for (i = 1; i <= RF_observationSize; i++) {
  ${trace.token}                RF_nativePrint(" %10.4f", RF_partSURVptr[pVarIdx][1][k][i]);
  ${trace.token}              }
  ${trace.token}            }
  ${trace.token}          }
  ${trace.token}        }
  ${trace.token}      }
  ${trace.token}      else {
  ${trace.token}        if (RF_rTargetFactorCount > 0) {
  ${trace.token}          for (j = 1; j <= RF_rTargetFactorCount; j++) {
  ${trace.token}            for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
  ${trace.token}              RF_nativePrint("\n  [pVarIdx][rTarget][level] = [%10d][%10d][%10d] \n", pVarIdx, j, k);
  ${trace.token}              RF_nativePrint("          ");
  ${trace.token}              for (i=1; i <= RF_observationSize; i++) {
  ${trace.token}                RF_nativePrint("%10d", i);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}              RF_nativePrint("          ");
  ${trace.token}              for (i=1; i <= RF_observationSize; i++) {
  ${trace.token}                RF_nativePrint("%10.4f", RF_partCLASptr[pVarIdx][j][k+1][i]);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}            }
  ${trace.token}          }
  ${trace.token}        }
  ${trace.token}        if (RF_rTargetNonFactorCount > 0) {
  ${trace.token}          for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
  ${trace.token}            RF_nativePrint("\n  [pVarIdx][rTarget] = [%10d][%10d] \n", pVarIdx, j);
  ${trace.token}            RF_nativePrint("          ");
  ${trace.token}            for (i=1; i <= RF_observationSize; i++) {
  ${trace.token}              RF_nativePrint("%10d", i);
  ${trace.token}            }
  ${trace.token}            RF_nativePrint("\n");
  ${trace.token}            RF_nativePrint("          ");
  ${trace.token}            for (i=1; i <= RF_observationSize; i++) {
  ${trace.token}              RF_nativePrint("%10.4f", RF_partREGRptr[pVarIdx][j][i]); 
  ${trace.token}            }
  ${trace.token}            RF_nativePrint("\n");
  ${trace.token}          }
  ${trace.token}        }
  ${trace.token}      }
  ${trace.token}      RF_nativePrint("\n");
  ${trace.token}    }
  }
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}      RF_nativePrint("\nupdatePartialCalculations() EXIT ...\n");
  ${trace.token}    }
}
void summarizePartialCalculations(uint       treeID,
                                  uint       pVarIdx) {
  double *ensembleDen;
  uint    membershipSize;
  uint i, j, k;
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nsummarizePartialCalculations() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & OUTP_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\n  pVarIdx:  %10d ", pVarIdx);
  ${trace.token}  }
  membershipSize  = RF_observationSize;
  ensembleDen = RF_oobEnsembleDen;
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (RF_eventTypeSize > 1) {
      if (RF_partialType == RF_PART_YRLS) {
        for (i = 1; i <= membershipSize; i++) {
          if (ensembleDen[i] > 0) {
            for (j = 1; j <= RF_eventTypeSize; j++) {
              RF_partSURVptr[pVarIdx][j][1][i] = RF_partSURVptr[pVarIdx][j][1][i] / ensembleDen[i];
            }
          }
        }
      }
      else if (RF_partialType == RF_PART_CIFN) {
        for (i = 1; i <= membershipSize; i++) {
          if (ensembleDen[i] > 0) {
            for (j = 1; j <= RF_eventTypeSize; j++) {
              for (k = 1; k <= RF_partialTimeLength; k++) {
                RF_partSURVptr[pVarIdx][j][k][i] = RF_partSURVptr[pVarIdx][j][k][i] / ensembleDen[i];
              }
            }
          }
        }
      }
      else if (RF_partialType == RF_PART_CHFN) {
        for (i = 1; i <= membershipSize; i++) {
          if (ensembleDen[i] > 0) {
            for (j = 1; j <= RF_eventTypeSize; j++) {
              for (k = 1; k <= RF_partialTimeLength; k++) {
                RF_partSURVptr[pVarIdx][j][k][i] = RF_partSURVptr[pVarIdx][j][k][i] / ensembleDen[i];
              }
            }
          }
        }
      }
    }   
    else {
      if (RF_partialType == RF_PART_MORT) {
        for (i = 1; i <= membershipSize; i++) {
          if (ensembleDen[i] > 0) {
            RF_partSURVptr[pVarIdx][1][1][i] = RF_partSURVptr[pVarIdx][1][1][i] / ensembleDen[i];
          }
        }
      }
      else if (RF_partialType == RF_PART_NLSN) {
        for (i = 1; i <= membershipSize; i++) {
          if (ensembleDen[i] > 0) {
            for (k = 1; k <= RF_partialTimeLength; k++) {
              RF_partSURVptr[pVarIdx][1][k][i] = RF_partSURVptr[pVarIdx][1][k][i] / ensembleDen[i];
            }
          }
        }
      }
      else if (RF_partialType == RF_PART_SURV) {
        for (i = 1; i <= membershipSize; i++) {
          if (ensembleDen[i] > 0) {
            for (k = 1; k <= RF_partialTimeLength; k++) {
              RF_partSURVptr[pVarIdx][1][k][i] =  RF_partSURVptr[pVarIdx][1][k][i] / ensembleDen[i];
            }
          }
        }
      }
    }
  }
  else {
    if (RF_rTargetFactorCount > 0) {
      for (i = 1; i <= membershipSize; i++) {
        if (ensembleDen[i] > 0) {
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              RF_partCLASptr[pVarIdx][j][k+1][i] = RF_partCLASptr[pVarIdx][j][k+1][i] / ensembleDen[i];
            }
            RF_partCLASptr[pVarIdx][j][1][i] = RF_nativeNaN;
          }
        }
        else {
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {          
              RF_partCLASptr[pVarIdx][j][k+1][i] = RF_nativeNaN;
            }
            RF_partCLASptr[pVarIdx][j][1][i] = RF_nativeNaN;
          }
        }
      }
    }
    if (RF_rTargetNonFactorCount > 0) {
      for (i = 1; i <= membershipSize; i++) {
        if (ensembleDen[i] > 0) {
          for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
            RF_partREGRptr[pVarIdx][j][i] = RF_partREGRptr[pVarIdx][j][i]  / ensembleDen[i];
          }
        }
      }
    }
  }
  ${trace.token}    if (getTraceFlag(treeID) & OUTP_DEF_TRACE) {
  ${trace.token}      RF_nativePrint("\nPartial outcomes after normalization:  pVarIdx %10d \n", pVarIdx);
  ${trace.token}      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
  ${trace.token}        if (RF_eventTypeSize > 1) {
  ${trace.token}          if (RF_partialType == RF_PART_YRLS) {
  ${trace.token}            RF_nativePrint("\nRF_PART_YRLS:  \n");
  ${trace.token}            for (j=1; j <= RF_eventTypeSize; j++) {
  ${trace.token}              RF_nativePrint("\n  [pVarIdx][event] = [%10d][%10d] \n", pVarIdx, j);
  ${trace.token}              RF_nativePrint("          ");
  ${trace.token}              for (i=1; i <= RF_observationSize; i++) {
  ${trace.token}                RF_nativePrint("%10d", i);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}              for (i=1; i <= RF_observationSize; i++) {
  ${trace.token}                RF_nativePrint("%10.4f", RF_partSURVptr[pVarIdx][j][1][i]);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}            }
  ${trace.token}          }
  ${trace.token}          else if (RF_partialType == RF_PART_CIFN) {
  ${trace.token}            RF_nativePrint("\nRF_PART_CIFN:  \n");
  ${trace.token}            for (i = 1; i <= RF_observationSize; i++) {
  ${trace.token}              RF_nativePrint("Indv %10d: \n", i);
  ${trace.token}              RF_nativePrint("                 time ");
  ${trace.token}              for (j = 1; j <= RF_eventTypeSize; j++) {
  ${trace.token}                RF_nativePrint("%10d ", j);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}              for (k=1; k <= RF_partialTimeLength; k++) {
  ${trace.token}                RF_nativePrint("%10d %10.4f ", k, RF_partialTime[k]);
  ${trace.token}                for (j = 1; j <= RF_eventTypeSize; j++) {
  ${trace.token}                  RF_nativePrint("%10.4f ", RF_partSURVptr[pVarIdx][j][k][i]);
  ${trace.token}                }
  ${trace.token}                RF_nativePrint("\n");
  ${trace.token}              }
  ${trace.token}            }
  ${trace.token}          }
  ${trace.token}          else if (RF_partialType == RF_PART_CHFN) {
  ${trace.token}            RF_nativePrint("\nRF_PART_CHFN:  \n");
  ${trace.token}            for (i = 1; i <= RF_observationSize; i++) {
  ${trace.token}              RF_nativePrint("Indv %10d: \n", i);
  ${trace.token}              RF_nativePrint("                 time ");
  ${trace.token}              for (j = 1; j <= RF_eventTypeSize; j++) {
  ${trace.token}                RF_nativePrint("%10d ", j);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}              for (k=1; k <= RF_partialTimeLength; k++) {
  ${trace.token}                RF_nativePrint("%10d %10.4f ", k, RF_partialTime[k]);
  ${trace.token}                for (j = 1; j <= RF_eventTypeSize; j++) {
  ${trace.token}                  RF_nativePrint("%10.4f ", RF_partSURVptr[pVarIdx][j][k][i]);
  ${trace.token}                }
  ${trace.token}                RF_nativePrint("\n");
  ${trace.token}              }
  ${trace.token}            }
  ${trace.token}          }
  ${trace.token}        }
  ${trace.token}        else {
  ${trace.token}          if (RF_partialType == RF_PART_MORT) {
  ${trace.token}            RF_nativePrint("\nRF_PART_MORT:  \n");
  ${trace.token}            for (i = 1; i <= RF_observationSize; i++) {
  ${trace.token}              RF_nativePrint("%10d %10.4f", i, RF_partSURVptr[pVarIdx][1][1][i]);
  ${trace.token}            }
  ${trace.token}            RF_nativePrint("\n");
  ${trace.token}           
  ${trace.token}          }
  ${trace.token}          else if (RF_partialType == RF_PART_NLSN) {
  ${trace.token}            RF_nativePrint("\nRF_PART_NLSN:  \n");
  ${trace.token}            RF_nativePrint("      index       time        obs -> ");  
  ${trace.token}            for (k = 1; k <= RF_partialTimeLength; k++) {
  ${trace.token}              RF_nativePrint("\n %10d %10.4f ", k, RF_partialTime[k]);
  ${trace.token}              for (i = 1; i <= RF_observationSize; i++) {
  ${trace.token}                RF_nativePrint(" %10.4f", RF_partSURVptr[pVarIdx][1][k][i]);
  ${trace.token}              }
  ${trace.token}            }
  ${trace.token}          }
  ${trace.token}          else if (RF_partialType == RF_PART_SURV) {
  ${trace.token}            RF_nativePrint("\nRF_PART_SURV:  \n");
  ${trace.token}            RF_nativePrint("      index       time        obs -> ");  
  ${trace.token}            for (k = 1; k <= RF_partialTimeLength; k++) {
  ${trace.token}              RF_nativePrint("\n %10d %10.4f ", k, RF_partialTime[k]);
  ${trace.token}              for (i = 1; i <= RF_observationSize; i++) {
  ${trace.token}                RF_nativePrint(" %10.4f", RF_partSURVptr[pVarIdx][1][k][i]);
  ${trace.token}              }
  ${trace.token}            }
  ${trace.token}          }
  ${trace.token}        }
  ${trace.token}      }
  ${trace.token}      else {
  ${trace.token}        if (RF_rTargetFactorCount > 0) {
  ${trace.token}          for (j = 1; j <= RF_rTargetFactorCount; j++) {
  ${trace.token}            for (k = 1; k <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
  ${trace.token}              RF_nativePrint("\n  [pVarIdx][rTarget][level] = [%10d][%10d][%10d] \n", pVarIdx, j, k-1);
  ${trace.token}              RF_nativePrint("          ");
  ${trace.token}              for (i=1; i <= RF_observationSize; i++) {
  ${trace.token}                RF_nativePrint("%10d", i);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}              RF_nativePrint("          ");
  ${trace.token}              for (i=1; i <= RF_observationSize; i++) {
  ${trace.token}                RF_nativePrint("%10.4f", RF_partCLASptr[pVarIdx][j][k][i]);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}            }
  ${trace.token}          }
  ${trace.token}        }
  ${trace.token}        if (RF_rTargetNonFactorCount > 0) {
  ${trace.token}          for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
  ${trace.token}            RF_nativePrint("\n  [pVarIdx][rTarget] = [%10d][%10d] \n", pVarIdx, j);
  ${trace.token}            RF_nativePrint("          ");
  ${trace.token}            for (i=1; i <= RF_observationSize; i++) {
  ${trace.token}              RF_nativePrint("%10d", i);
  ${trace.token}            }
  ${trace.token}            RF_nativePrint("\n");
  ${trace.token}            RF_nativePrint("          ");
  ${trace.token}            for (i=1; i <= RF_observationSize; i++) {
  ${trace.token}              RF_nativePrint("%10.4f", RF_partREGRptr[pVarIdx][j][i]); 
  ${trace.token}            }
  ${trace.token}            RF_nativePrint("\n");
  ${trace.token}          }
  ${trace.token}        }
  ${trace.token}      }
  ${trace.token}      RF_nativePrint("\n");
  ${trace.token}    }
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nsummarizePartialCalculations() EXIT ...\n");
  ${trace.token}  }
}
