
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
void getAndUpdatePartialMembership(uint treeID, Node *root) {
  uint i, j;
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
      updatePartialCalculations(treeID, i, membership);
    }
  }
  free_new_vvector(membership, 1, RF_observationSize, NRUTIL_TPTR);
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
  terminalFlag = TRUE;
  obsSize = RF_observationSize;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    info = parent -> splitInfo;
    terminalFlag = FALSE;
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
      if (info -> mwcpSizeAbs[1] > 0) {
        if (primaryPartialIndex > 0) { 
          factorValue = (uint) RF_partialValue[partialIndex];
        }
        else if (secondaryPartialIndex > 0) {
          factorValue = (uint) RF_partialValue2[secondaryPartialIndex];            
        }
        else {
          factorValue = (uint) xArray[info -> randomVar[1]][allMembrIndx[i]];
        }
        daughterFlag = splitOnFactor(factorValue, (uint*) info -> randomPts[1]);
      }
      else {
        if (primaryPartialIndex > 0) { 
          continuousValue = RF_partialValue[partialIndex];
        }
        else if (secondaryPartialIndex > 0) {
          continuousValue = RF_partialValue2[secondaryPartialIndex];            
        }
        else {
          continuousValue = xArray[info -> randomVar[1]][allMembrIndx[i]];
        }
        daughterFlag =  (( ((double*) info -> randomPts[1])[1] - continuousValue) >= 0.0) ? LEFT : RIGHT;
      }
      indicator[allMembrIndx[i]] = daughterFlag;
      if (daughterFlag == LEFT) {
        leftAllMembrSize ++;
      }
      else {
        rghtAllMembrSize ++;
      }
    }  
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
    partialMembershipGeneric(treeID,
                          parent -> left,                             
                          partialIndex,
                          leftAllMembrIndx,
                          leftAllMembrSize,
                          xArray,
                          membership);
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
  }
  if (terminalFlag) {
    for (i = 1; i <= allMembrSize; i++) {
      membership[allMembrIndx[i]] = parent -> mate;
    }
  }  
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
      if (info -> mwcpSizeAbs[1] > 0) {
        if (primaryPartialIndex > 0) { 
          factorValue = (uint) RF_partialValue[partialIndex];
        }
        else if (secondaryPartialIndex > 0) {
          factorValue = (uint) RF_partialValue2[secondaryPartialIndex];            
        }
        else {
          factorValue = (uint) xArray[info -> randomVar[1]][individual];
        }
        daughterFlag = splitOnFactor(factorValue, (uint*) info -> randomPts[1]);
      }
      else {
        if (primaryPartialIndex > 0) { 
          continuousValue = RF_partialValue[partialIndex];
        }
        else if (secondaryPartialIndex > 0) {
          continuousValue = RF_partialValue2[secondaryPartialIndex];            
        }
        else {
          continuousValue = xArray[info -> randomVar[1]][individual];
        }
        daughterFlag =  (( ((double*) info -> randomPts[1])[1] - continuousValue) >= 0.0) ? LEFT : RIGHT;
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
  membership[individual] = parent -> mate;
  free_ulvector(rootMWCPoffset, 1, 1);
  free_ulvector(nodeAbsMWCPoffset, 1, 1);
}
void updatePartialCalculations (uint       treeID,
                                uint       pVarIdx,
                                Terminal **partialMembership) {
  Terminal *terminalNode;
  uint  *membershipIndex;
  uint   membershipSize;
  uint   i, j, k;
  uint   ii;
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
  }
}
void summarizePartialCalculations(uint       treeID,
                                  uint       pVarIdx) {
  double *ensembleDen;
  uint    membershipSize;
  uint i, j, k;
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
}
