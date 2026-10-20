
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "classification.h"
#include "termOps.h"
#include "nrutil.h"
#include "error.h"
void getMultiClassProb (uint       treeID,
                        Terminal  *parent,
                        uint      *repMembrIndx,
                        uint       repMembrSize,
                        uint      *allMembrIndx,
                        uint       allMembrSize,
                        uint      *rmbrIterator) {
  uint *membershipIndex;
  uint  membershipSize;
  double maxValue, maxClass;
  uint i, j, k;
  if ( !(RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2) ) {
    membershipIndex = allMembrIndx;
    membershipSize = parent -> membrCount = allMembrSize;
    if (RF_optHigh & OPT_MEMB_INCG) {
      membershipIndex = RF_AMBR_ID_ptr[treeID];
    }
  }
  else {
    membershipIndex = repMembrIndx;
    membershipSize = parent -> membrCount = repMembrSize;
    if (RF_optHigh & OPT_MEMB_INCG) {
      membershipIndex = RF_RMBR_ID_ptr[treeID];
    }
  }
  if (membershipSize == 0) {
    if (!(RF_opt & OPT_OUTC_TYPE)) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Zero node count encountered in (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  if (!(RF_optHigh & OPT_TERM_INCG)) {
    stackMultiClassProb(parent, RF_rFactorCount, RF_rFactorSize);
    for (j=1; j <= RF_rFactorCount; j++) {
      for (k=1; k <= RF_rFactorSize[j]; k++) {
        (parent -> multiClassProb)[j][k] = 0;
      }
    }
    if (RF_optHigh & OPT_MEMB_OUTG) {
      for (i = 1; i <= membershipSize; i++) {
        RF_RMBR_ID_ptr[treeID][++(*rmbrIterator)] = membershipIndex[i];
        for (j=1; j <= RF_rFactorCount; j++) {
          (parent -> multiClassProb)[j][(uint) RF_response[treeID][RF_rFactorIndex[j]][membershipIndex[i]]] ++;
        }
      }
    }
    else if (RF_optHigh & OPT_MEMB_INCG) {
      for (i = 1; i <= membershipSize; i++) {
        ++(*rmbrIterator);
        for (j=1; j <= RF_rFactorCount; j++) {
          (parent -> multiClassProb)[j][(uint) RF_response[treeID][RF_rFactorIndex[j]][ membershipIndex[*rmbrIterator] ]] ++;
        }
      }
    }
    else {
      for (i = 1; i <= membershipSize; i++) {
        for (j=1; j <= RF_rFactorCount; j++) {
          (parent -> multiClassProb)[j][(uint) RF_response[treeID][RF_rFactorIndex[j]][membershipIndex[i]]] ++;
        }
      }
    }
    for (j = 1; j <= RF_rFactorCount; j++) {
      maxValue = 0;
      maxClass = 0;
      for (k=1; k <= RF_rFactorSize[j]; k++) {
        if (maxValue < (double) (parent -> multiClassProb[j][k])) {
          maxValue = (double) parent -> multiClassProb[j][k];
          maxClass = (double) k;
        }
      }
      (parent -> maxClass)[j] = maxClass;
    }
  }
  else {
    stackMultiClassProb(parent, RF_rFactorCount, RF_rFactorSize);
    for (j = 1; j <= RF_rFactorCount; j++) {
      for (k=1; k <= RF_rFactorSize[j]; k++) {      
        (parent -> multiClassProb)[j][k] = RF_TN_CLAS_ptr[treeID][parent -> nodeID][j][k];      
      }
    }
    for (j = 1; j <= RF_rFactorCount; j++) {
      maxValue = 0;
      maxClass = 0;
      for (k=1; k <= RF_rFactorSize[j]; k++) {
        if (maxValue < (double) (parent -> multiClassProb[j][k])) {
          maxValue = (double) parent -> multiClassProb[j][k];
          maxClass = (double) k;
        }
      }
      (parent -> maxClass)[j] = maxClass;
    }
  }
}
void updateEnsembleMultiClass(char      mode,
                              uint      treeID,
                              char      normalizationFlag,
                              char      omitDenominator) {
  char oobFlag, fullFlag, outcomeFlag;
  Terminal ***termMembershipPtr;
  uint    *membershipIndex;
  uint     membershipSize;
  double   ***ensembleCLSptr;
  double   ***ensembleCLSnum;
  double     *ensembleDen;
#ifdef _OPENMP
  omp_lock_t   *lockDENptr;
#endif
  ensembleCLSnum = NULL;  
  ensembleDen    = NULL;  
  oobFlag = fullFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    termMembershipPtr = RF_ftTermMembership;
    break;
  default:
    if (RF_opt & OPT_OENS) {
      if (RF_oobSize[treeID] > 0) {
        oobFlag = TRUE;
      }
    }
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    termMembershipPtr = RF_tTermMembership;
    break;
  }
  outcomeFlag = TRUE;
  while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
    if (oobFlag == TRUE) {
      ensembleCLSptr = RF_oobEnsembleCLSptr;
      ensembleCLSnum = RF_oobEnsembleCLSnum;
      ensembleDen    = RF_oobEnsembleDen;
      membershipSize  = RF_oobSize[treeID];
      membershipIndex = RF_oobMembershipIndex[treeID];
#ifdef _OPENMP
      lockDENptr      = RF_lockDENoens;
#endif
    }
    else {
      ensembleCLSptr = RF_fullEnsembleCLSptr;
      ensembleCLSnum = RF_fullEnsembleCLSnum;
      ensembleDen    = RF_fullEnsembleDen;
      switch (mode) {
      case RF_PRED:
        membershipSize = RF_fobservationSize;
        membershipIndex = RF_fidentityMembershipIndex;
        break;
      default:
        membershipSize  = RF_observationSize;
        membershipIndex = RF_identityMembershipIndex;
        break;
      }
#ifdef _OPENMP
      lockDENptr      = RF_lockDENfens;
#endif
    }
    for (uint i = 1; i <= membershipSize; i++) {
      Terminal *parent;
      char selectionFlag;
      uint j, k, ii;
      ii = membershipIndex[i];
      parent = termMembershipPtr[treeID][ii];
      selectionFlag = TRUE;
      if (RF_opt & OPT_OUTC_TYPE) {
        if ((parent -> membrCount) > 0) {
        }
        else {
          selectionFlag = FALSE;
        }
      }
      if (selectionFlag) {
#ifdef _OPENMP        
        omp_set_lock(&(lockDENptr[ii]));
#endif
        if(!omitDenominator) {
          ensembleDen[ii] ++;          
          if (outcomeFlag == TRUE) {
            if (RF_optHigh & OPT_CSE) {              
              RF_cseDENptr[ii] ++;
            }
            if (RF_opt & OPT_VIMP) {
              RF_blkEnsembleDen[ii] ++;
            }
          }
        }
        for (j = 1; j <= RF_rTargetFactorCount; j++) {
          for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
            ensembleCLSnum[j][k][ii] += (double) (parent -> multiClassProb)[RF_rFactorMap[RF_rTargetFactor[j]]][k] / (double) (parent -> membrCount);
          }
        }
        if (outcomeFlag == TRUE) {
          if (RF_optHigh & OPT_CSE) {              
            for (j = 1; j <= RF_rTargetFactorCount; j++) {
              RF_cseNumCLSptr[j][ii] +=
                ( (parent -> maxClass)[RF_rFactorMap[RF_rTargetFactor[j]]] ==
                  (uint) RF_response[treeID][RF_rFactorMap[RF_rTargetFactor[j]]][ii] ) ? 1 : 0;
            }
          }
          if (RF_opt & OPT_VIMP) {
            for (j = 1; j <= RF_rTargetFactorCount; j++) {
              for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
                RF_blkEnsembleCLSnum[j][k][ii] += (double) (parent -> multiClassProb)[RF_rFactorMap[RF_rTargetFactor[j]]][k] / (double) (parent -> membrCount);
              }
            }
          }
        }
        if (outcomeFlag && normalizationFlag) {
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              ensembleCLSptr[j][k][ii] = ensembleCLSnum[j][k][ii] / ensembleDen[ii];
            }
          }
        }
#ifdef _OPENMP
        omp_unset_lock(&(lockDENptr[ii]));
#endif
      }
    }  
    if (outcomeFlag == TRUE) {
      outcomeFlag = FALSE;
    }
    if (oobFlag == TRUE) {
      oobFlag = FALSE;
    }
    else {
      fullFlag = FALSE;
    }
  }  
}
double getBrierScore(uint     obsSize,
                     uint     rTarget,
                     double  *responsePtr,
                     double **outcomeCLS,
                     double  *denomCount,
                     double  *cpv) {
  uint k;
  uint against;
  uint *oaaResponse;
  uint cumDenomCount;
  double result;
  oaaResponse = uivector(1, obsSize);
  result = 0.0;
  cumDenomCount = 0;
  for (k = 1; k <= obsSize; k ++) {
    if (denomCount[k] != 0) {
      cumDenomCount += 1;
    }
  }
  for (against = 1; against <= RF_rFactorSize[RF_rFactorMap[rTarget]]; against++) {
    for (k = 1; k <= obsSize; k ++) {
      if ((uint) responsePtr[k] == against) {
        oaaResponse[k] = 1;
      }
      else {
        oaaResponse[k] = 0;
      }
    }
    cpv[against] = 0.0;
    for (k = 1; k <= obsSize; k ++) {
      if (denomCount[k] != 0) {
        cpv[against] += pow(((double) oaaResponse[k] - outcomeCLS[against][k]), 2.0);
      }
    }
    if (cumDenomCount == 0) {
      cpv[against] = RF_nativeNaN;
    }
    else {
      cpv[against] = cpv[against] / (double) cumDenomCount;
      result += cpv[against];
    }
  }
  if (cumDenomCount == 0) {
    result = RF_nativeNaN;
  }
  else {
    result = result  * RF_rFactorSize[RF_rFactorMap[rTarget]] / (RF_rFactorSize[RF_rFactorMap[rTarget]] - 1);
  }
  for (against = 1; against <= RF_rFactorSize[RF_rFactorMap[rTarget]]; against++) {
    cpv[against] = cpv[against] * RF_rFactorSize[RF_rFactorMap[rTarget]] / (RF_rFactorSize[RF_rFactorMap[rTarget]] - 1);
  }
  free_uivector(oaaResponse, 1, obsSize);
  return result;
}
void getConditionalClassificationIndexGrow(uint     size,
                                       uint     rTarget,
                                       double  *responsePtr,
                                       double **outcomeCLS,
                                       double  *maxVote,
                                       double  *denomCount,
                                       double  *cpv) {
  uint i, k;
  uint cumDenomCount;
  uint *condClassificationCount;
  cumDenomCount = 0;
  condClassificationCount = uivector(1, RF_rFactorSize[RF_rFactorMap[rTarget]]);
  for (k=1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
    cpv[k] = condClassificationCount[k] = 0;
  }
  for (i = 1; i <= size; i++) {
    condClassificationCount[(uint) responsePtr[i]] ++;
    if (denomCount[i] != 0) {
      cumDenomCount += 1;
      if (responsePtr[i] == maxVote[i]) {
        cpv[(uint) responsePtr[i]] += 1.0;
      }
    }
  }  
  if (cumDenomCount == 0) {
    for (k=1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
      cpv[k] = RF_nativeNaN;
    }
  }
  else {
    for (k=1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
      if (condClassificationCount[k] != 0) {
        cpv[k] = 1.0 - cpv[k] / (double) condClassificationCount[k];
      }
      else {
        cpv[k] = RF_nativeNaN;
      }
    }
  }
  free_uivector(condClassificationCount, 1, RF_rFactorSize[RF_rFactorMap[rTarget]]);
  return;
}
void getConditionalClassificationIndexPred(uint     size,
                                           uint     rTarget,
                                           double  *responsePtr,
                                           double **outcomeCLS,
                                           double  *maxVote,
                                           double  *denomCount,
                                           double  *cpv) {
  uint i, k;
  uint cumDenomCount;
  uint *condClassificationCount;
  cumDenomCount = 0;
  condClassificationCount = uivector(1, RF_rFactorSize[RF_rFactorMap[rTarget]]);
  for (k=1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
    cpv[k] = condClassificationCount[k] = 0;
  }
  for (i = 1; i <= size; i++) {
    if ( (uint) responsePtr[i] <= RF_rFactorSize[RF_rFactorMap[rTarget]]) {
      condClassificationCount[(uint) responsePtr[i]] ++;
      if (denomCount[i] != 0) {
        cumDenomCount += 1;
        if (responsePtr[i] == maxVote[i]) {
          cpv[(uint) responsePtr[i]] += 1.0;
        }
      }
    }
  }  
  if (cumDenomCount == 0) {
    for (k=1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
      cpv[k] = RF_nativeNaN;
    }
  }
  else {
    for (k=1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
      if (condClassificationCount[k] != 0) {
        cpv[k] = 1.0 - cpv[k] / (double) condClassificationCount[k];
      }
      else {
        cpv[k] = RF_nativeNaN;
      }
    }
  }
  free_uivector(condClassificationCount, 1, RF_rFactorSize[RF_rFactorMap[rTarget]]);
  return;
}
double getClassificationIndex(uint     size,
                              uint     rTarget,
                              double  *responsePtr,
                              double  *denomCount,
                              double  *maxVote) {
  uint i;
  uint cumDenomCount;
  double result;
  cumDenomCount = 0;
  result = 0.0;
  for (i=1; i <= size; i++) {
    if (denomCount[i] > 0) {
      cumDenomCount += 1;
      if (responsePtr[i] == maxVote[i]) {
        result += 1.0;
      }
    }
    else {
      maxVote[i] = RF_nativeNaN;
    }
  }  
  if (cumDenomCount == 0) {
    result = RF_nativeNaN;
  }
  else {
    result = 1.0 - result / (double) cumDenomCount;
  }
  return result;
}
double getGMeanIndexGrow(uint    size,
                         uint    rTarget,
                         double *responsePtr,
                         double *denomCount,
                         double *maxVote) {
  uint i, k;
  uint cumDenomCount;
  double *trueRate, *falseRate;
  double denom, result;
  cumDenomCount = 0;
  result = 1.0;
  trueRate  = dvector(1, RF_rFactorSize[RF_rFactorMap[rTarget]]);
  falseRate = dvector(1, RF_rFactorSize[RF_rFactorMap[rTarget]]);
  for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
    trueRate[k] = falseRate[k] = 0;
  }
  for (i = 1; i <= size; i++) {
    if (denomCount[i] > 0) {
      cumDenomCount += 1;
      if (responsePtr[i] == maxVote[i]) {
        trueRate[(uint) responsePtr[i]] += 1.0;
      }
      else {
        falseRate[(uint) responsePtr[i]] += 1.0;
      }
    }
  }  
  if (cumDenomCount == 0) {
    result = RF_nativeNaN;
  }
  else {
    for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
      denom = trueRate[k] + falseRate[k];
      if (denom > 0) {
        result = result * trueRate[k] / denom; 
      }
      else {
        result = result * (1 + trueRate[k]) / (1 + denom); 
      }
    }
    result = 1.0 - sqrt(result);
  }
  free_dvector(trueRate, 1, RF_rFactorSize[RF_rFactorMap[rTarget]]);
  free_dvector(falseRate, 1, RF_rFactorSize[RF_rFactorMap[rTarget]]);
  return result;
}
double getGMeanIndexPred(uint    size,
                         uint    rTarget,
                         double *responsePtr,
                         double *denomCount,
                         double *maxVote) {
  uint i, k;
  uint cumDenomCount;
  double *trueRate, *falseRate;
  double denom, result;
  cumDenomCount = 0;
  result = 1.0;
  trueRate  = dvector(1, RF_rFactorSizeTest[RF_rFactorMap[rTarget]]);
  falseRate = dvector(1, RF_rFactorSizeTest[RF_rFactorMap[rTarget]]);
  for (k = 1; k <= RF_rFactorSizeTest[RF_rFactorMap[rTarget]]; k++) {
    trueRate[k] = falseRate[k] = 0;
  }
  for (i = 1; i <= size; i++) {
    if (denomCount[i] > 0) {
      cumDenomCount += 1;
      if (responsePtr[i] == maxVote[i]) {
        trueRate[(uint) responsePtr[i]] += 1.0;
      }
      else {
        falseRate[(uint) responsePtr[i]] += 1.0;
      }
    }
  }  
  if (cumDenomCount == 0) {
    result = RF_nativeNaN;
  }
  else {
    for (k = 1; k <= RF_rFactorSizeTest[RF_rFactorMap[rTarget]]; k++) {
      denom = trueRate[k] + falseRate[k];
      if (denom > 0) {
        result = result * trueRate[k] / denom; 
      }
      else {
        result = result * (1 + trueRate[k]) / (1 + denom); 
      }
    }
    result = 1.0 - sqrt(result);
  }
  free_dvector(trueRate, 1, RF_rFactorSizeTest[RF_rFactorMap[rTarget]]);
  free_dvector(falseRate, 1, RF_rFactorSizeTest[RF_rFactorMap[rTarget]]);
  return result;
}
void restoreMultiClassProb(uint treeID) {
  LeafLinkedObj *leafLinkedPtr;
  Terminal *parent;
  uint leaf;
  uint j, k;
  leafLinkedPtr = RF_leafLinkedObjHead[treeID] -> fwdLink;
  while (leafLinkedPtr != NULL) {
    parent = leafLinkedPtr -> termPtr;
    leaf = parent -> nodeID;
    if ((parent -> membrCount) > 0) {
      for (j = 1; j <= RF_rFactorCount; j++) {
        for (k = 1; k <= RF_rFactorSize[j]; k++) {
          (parent -> multiClassProb)[j][k] = RF_TN_CLAS_ptr[treeID][leaf][j][k];
        }
      }
    }
    else {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Zero node count encountered in restoreMultiClassProb() in (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
    leafLinkedPtr = leafLinkedPtr -> fwdLink;
  }
}
void getMaxVote(uint     size,
                uint     rTarget,
                double **outcomeCLS,
                double  *denomCount,
                double  *maxVote) {
  uint i,k;
  if ((RF_opt & OPT_CLAS_RFQ) && RF_rFactorMinorityFlag[RF_rFactorMap[rTarget]]) {
    uint minorityClass = RF_rFactorMinority[RF_rFactorMap[rTarget]];
    uint majorityClass = RF_rFactorMajority[RF_rFactorMap[rTarget]];
    double threshold   = RF_rFactorThreshold[RF_rFactorMap[rTarget]];
    for (i = 1; i <= size; i++) {
      if (denomCount[i] > 0) {
        if (outcomeCLS[minorityClass][i] >= threshold) {
          maxVote[i] = (double) minorityClass;
        }
        else {
          maxVote[i] = (double) majorityClass;
        }
      }
      else {
        maxVote[i] = RF_nativeNaN;
      }
    }  
  }
  else {
    double maxValue, maxClass;
    for (i = 1; i <= size; i++) {
      if (denomCount[i] > 0) {
        maxValue = 0.0;
        maxClass = 0.0;
        for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
          if (maxValue <= outcomeCLS[k][i]) {
            maxValue = outcomeCLS[k][i];
            maxClass = (double) k;
          }
        }
        maxVote[i] = maxClass;
      }
      else {
        maxVote[i] = RF_nativeNaN;
      }
    }  
 } 
}
