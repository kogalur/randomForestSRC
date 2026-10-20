
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "importanceRand.h"
#include "polarity.h"
#include "splitInfo.h"
#include "splitGreedy.h"
#include "rfsrcUtil.h"
#include "nodeOps.h"
#include "nrutil.h"
void getRandomMembership (char       mode,
                          uint       treeID,
                          Terminal **vimpMembership,
                          uint       p) {
  Node    *rootPtr;
  uint    *membershipIndex;
  uint     membershipSize;
  double **xArray;
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
    vimpMembership[ii] = randomMembership(treeID, rootPtr, ii, p, xArray) -> mate;
  }
}
Node *randomMembershipGeneric(uint     treeID,
                              Node    *parent,
                              uint     individual,                            
                              uint     vimpX,
                              double **xArray) {
  char daughterFlag;
  char randomSplitFlag;
  Node *result;
  SplitInfo *info;
  uint repMembrSize, leftRepMembrSize;
  double alpha;
  result = parent;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    info = parent -> splitInfo;
    randomSplitFlag = FALSE;
    if (vimpX > 0) {
      if ((uint) info -> randomVar[1] == vimpX) {
        randomSplitFlag = TRUE;
      }
    }
    else {
      if (RF_importanceFlag[info -> randomVar[1]] == TRUE) {
        randomSplitFlag = TRUE;
      }
    }
    if (randomSplitFlag == TRUE) {
      repMembrSize = parent -> repMembrSize;
      leftRepMembrSize = (parent -> left) -> repMembrSize;
      alpha = ran1D(treeID);
      if (alpha <= RF_vimpThreshold) {
        if (alpha <= (double) leftRepMembrSize / (repMembrSize * RF_vimpThreshold)) {
          daughterFlag = LEFT;
        }
        else {
          daughterFlag = RIGHT;
        }
      }
      else {
        daughterFlag = getDaughterPolarity(0, info, individual, xArray);
      }
    }  
    else {
      daughterFlag = getDaughterPolarity(0, info, individual, xArray);
    }  
    if (daughterFlag == LEFT) {
      result = randomMembershipGeneric(treeID, parent ->  left, individual, vimpX, xArray);
    }
    else {
      result = randomMembershipGeneric(treeID, parent -> right, individual, vimpX, xArray);
    }
  }
  return result;
}
Node *randomMembershipJIT(uint     treeID,
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
  char randomSplitFlag;
  uint repMembrSize, leftRepMembrSize;
  double alpha;
  uint adj, i, k;
  rmbrDummyIter = ambrDummyIter = 0;
  nodeAbsIndex = rootIndex = RF_restoreTreeOffset[treeID];
  rmbrAbsOffset = ambrAbsOffset = 0;
  rootMWCPoffset    = ulvector(1, 1);
  nodeAbsMWCPoffset = ulvector(1, 1);
  rootMWCPoffset[1] = RF_restoreMWCPoffset[1][treeID];
  parseFlag = TRUE;
  Node *parent = root;
  while (parseFlag) {
    if (parent -> nodeID == 0) {
      parent -> nodeID = RF_nodeID_[nodeAbsIndex];
      parent -> repMembrSize = RF_nodeSZ_[nodeAbsIndex];
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
          info -> randomVar[k] = RF_parmID_[k][nodeAbsIndex];
          info -> mwcpSizeAbs[k] = RF_mwcpSZ_[k][nodeAbsIndex];
          if (RF_mwcpSZ_[k][nodeAbsIndex] > 0) {
            nodeAbsMWCPoffset[k] = rootMWCPoffset[k] + RF_fsrecID_[k][nodeAbsIndex];
            info -> randomPts[k] = uivector(1, RF_mwcpSZ_[k][nodeAbsIndex]);
            for (i = 1; i <= RF_mwcpSZ_[k][nodeAbsIndex]; i++) {
              ((uint *) info -> randomPts[k])[i] = RF_mwcpPT_[k][nodeAbsMWCPoffset[k]];
              nodeAbsMWCPoffset[k] ++;
            }
          }
          else {
            info -> randomPts[k] = dvector(1, 1);
            ((double *) info -> randomPts[k])[1] =  RF_contPT_[k][nodeAbsIndex];
          }
        }
      }
      else {
        info = parent -> splitInfo = NULL;
      }
    }  
    else {
      info = parent -> splitInfo;
    }
    if (info != NULL) {
      randomSplitFlag = FALSE;
      if (vimpX > 0) {
        if ((uint) info -> randomVar[1] == vimpX) {
          randomSplitFlag = TRUE;
        }
      }
      else {
        if(RF_importanceFlag[info -> randomVar[1]] == TRUE) {
          randomSplitFlag = TRUE;
        }
      }
      if(randomSplitFlag == TRUE) {
        repMembrSize = parent -> repMembrSize;
        if (parent -> left != NULL) {
          leftRepMembrSize = (parent -> left) -> repMembrSize;
        }
        else {
          leftRepMembrSize = RF_nodeSZ_[nodeAbsIndex + 1];
        }
        alpha = ran1D(treeID);
        if (alpha <= RF_vimpThreshold) {
          if (alpha <= (double) leftRepMembrSize / (repMembrSize * RF_vimpThreshold)) {
            daughterFlag = LEFT;
          }
          else {
            daughterFlag = RIGHT;
          }
        }
        else {
          daughterFlag = getDaughterPolarity(0, info, individual, xArray);
        }
      }  
      else {
        daughterFlag = getDaughterPolarity(0, info, individual, xArray);
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
    }
    else {
      if (RF_tTermList[treeID][parent -> nodeID] != NULL) {
      }
      else {
        RF_leafLinkedObjTail[treeID] = makeAndSpliceLeafLinkedObj(RF_leafLinkedObjTail[treeID],
                                                                  parent,
                                                                  RF_TN_RCNT_ptr[treeID][parent -> nodeID],
                                                                  RF_TN_ACNT_ptr[treeID][parent -> nodeID]);
        RF_tTermList[treeID][parent -> nodeID] = RF_leafLinkedObjTail[treeID] -> termPtr;
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
      parseFlag = FALSE;
    }
  }  
  free_ulvector(rootMWCPoffset, 1, 1);
  free_ulvector(nodeAbsMWCPoffset, 1, 1);
  return parent;
}
