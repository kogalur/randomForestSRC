
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "regression.h"
#include "termOps.h"
#include "error.h"
void getMeanResponse(uint       treeID,
                     Terminal  *parent,
                     uint      *repMembrIndx,
                     uint       repMembrSize,
                     uint      *allMembrIndx,
                     uint       allMembrSize,
                     uint      *rmbrIterator) {
  uint *membershipIndex;
  uint  membershipSize;
  uint i, j;
  membershipIndex = repMembrIndx;
  membershipSize = parent -> membrCount = repMembrSize;
  if (RF_optHigh & OPT_MEMB_INCG) {
    membershipIndex = RF_RMBR_ID_ptr[treeID];
  }
  if (membershipSize == 0) {
    if (!(RF_opt & OPT_OUTC_TYPE)) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Zero node count encountered in (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  if (RF_opt & OPT_QUANTLE) {
    if (membershipSize > 0) {
      stackMemberStream(parent, membershipSize);
    }
  }
  if (!(RF_optHigh & OPT_TERM_INCG)) {
    stackMeanResponse(parent, RF_rNonFactorCount);
    for (j=1; j <= RF_rNonFactorCount; j++) {
      (parent -> meanResponse)[j] = 0.0;
    }
    if (RF_optHigh & OPT_MEMB_OUTG) {
      for (i = 1; i <= membershipSize; i++) {
        RF_RMBR_ID_ptr[treeID][++(*rmbrIterator)] = membershipIndex[i];
        if (RF_opt & OPT_QUANTLE) {
          parent -> membrStream[i] = membershipIndex[i];
        }
        for (j = 1; j <= RF_rNonFactorCount; j++) {
          (parent -> meanResponse)[j] += RF_response[treeID][RF_rNonFactorIndex[j]][membershipIndex[i]];
        }
      }
    }
    else if (RF_optHigh & OPT_MEMB_INCG) {
      for (i = 1; i <= membershipSize; i++) {
        ++(*rmbrIterator);
        if (RF_opt & OPT_QUANTLE) {
          parent -> membrStream[i] = membershipIndex[*rmbrIterator];
        }
        for (j = 1; j <= RF_rNonFactorCount; j++) {
          (parent -> meanResponse)[j] += RF_response[treeID][RF_rNonFactorIndex[j]][ membershipIndex[*rmbrIterator] ];
        }
      }
    }
    else {
      for (i = 1; i <= membershipSize; i++) {
        if (RF_opt & OPT_QUANTLE) {
          parent -> membrStream[i] = membershipIndex[i];
        }
        for (j=1; j <= RF_rNonFactorCount; j++) {
          (parent -> meanResponse)[j] += RF_response[treeID][RF_rNonFactorIndex[j]][membershipIndex[i]];
        }
      }
    }
    if (membershipSize > 0) {
      for (j = 1; j <= RF_rNonFactorCount; j++) {
        (parent -> meanResponse)[j] = (parent -> meanResponse)[j] / (double) membershipSize;
      }
    }
  }
  else {
    stackMeanResponse(parent, RF_rNonFactorCount);
    for (j = 1; j <= RF_rNonFactorCount; j++) {
      (parent -> meanResponse)[j] = RF_TN_REGR_ptr[treeID][parent -> nodeID][j];
    }
  }
}
void updateEnsembleMean(char     mode,
                        uint     treeID,
                        char     normalizationFlag,
                        char     omitDenominator) {
  char oobFlag, fullFlag, selectionFlag, outcomeFlag;
  Terminal ***termMembershipPtr;
  uint    *membershipIndex;
  uint     membershipSize;
  double    **ensembleRGRptr;
  double    **ensembleRGRnum;
  double     *ensembleDen;
#ifdef _OPENMP
  omp_lock_t   *lockDENptr;
#endif
  ensembleRGRnum = NULL;  
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
      ensembleRGRptr = RF_oobEnsembleRGRptr;
      ensembleRGRnum = RF_oobEnsembleRGRnum;
      ensembleDen    = RF_oobEnsembleDen;
      membershipSize  = RF_oobSize[treeID];
      membershipIndex = RF_oobMembershipIndex[treeID];
#ifdef _OPENMP
      lockDENptr      = RF_lockDENoens;
#endif
    }
    else {
      ensembleRGRptr = RF_fullEnsembleRGRptr;
      ensembleRGRnum = RF_fullEnsembleRGRnum;
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
      uint j, ii;
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
        for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
          ensembleRGRnum[j][ii] += (parent -> meanResponse)[RF_rNonFactorMap[RF_rTargetNonFactor[j]]];
        }
        if (outcomeFlag == TRUE) {
          if (RF_optHigh & OPT_CSE) {              
            for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
              RF_cseNumRGRptr[j][ii] += 
                pow((parent -> meanResponse)[RF_rNonFactorMap[RF_rTargetNonFactor[j]]] -
                    RF_response[treeID][RF_rNonFactorMap[RF_rTargetNonFactor[j]]][ii], 2.0);
            }
          }
          if (RF_opt & OPT_VIMP) {
            for (j = 1; j <= RF_rTargetNonFactorCount; j++) {          
              RF_blkEnsembleRGRnum[j][ii] += (parent -> meanResponse)[RF_rNonFactorMap[RF_rTargetNonFactor[j]]];
            }
          }
        }
        if (outcomeFlag && normalizationFlag) {
          for (j=1; j <= RF_rTargetNonFactorCount; j++) {
            ensembleRGRptr[j][ii] = ensembleRGRnum[j][ii] / ensembleDen[ii];
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
double getMeanSquareError(uint    size,
                          double *responsePtr,
                          double *predictedOutcome,
                          double *denomCount) {
  uint i;
  uint cumDenomCount;
  double result;
  cumDenomCount = 0;
  result = 0.0;
  for (i = 1; i <= size; i++) {
    if (denomCount[i] != 0) {
      cumDenomCount += 1;
      result += pow (responsePtr[i] - predictedOutcome[i], 2.0);
    }
  }  
  if (cumDenomCount == 0) {
    result = RF_nativeNaN;
  }
  else {
    result = result / (double) cumDenomCount;
  }
  return result;
}
char getVarianceClassic(uint    repMembrSize,
                        uint   *repMembrIndx,
                        uint    nonMissMembrSize,
                        uint   *nonMissMembrIndx,
                        double *targetResponse,
                        double *mean,
                        double *variance) {
  uint i;
  double meanResult, varResult;
  char result;
  uint *genIndx;
  uint  genSize;
  if (nonMissMembrIndx == NULL) {
    genIndx = RF_identityMembershipIndex;
    genSize = repMembrSize;
  }
  else {
    genIndx = nonMissMembrIndx;
    genSize = nonMissMembrSize;
  }
  meanResult = 0.0;
  for (i = 1; i <= genSize; i++) {
      meanResult += targetResponse[repMembrIndx[genIndx[i]]];
  }
  if (genSize > 0) {
    meanResult = meanResult / (double) genSize;  
  }
  else {
    meanResult = RF_nativeNaN;
  }
  varResult = 0.0;
  if(!RF_nativeIsNaN(meanResult)) {
    for (i = 1; i <= genSize; i++) {
      varResult += pow(meanResult - targetResponse[repMembrIndx[genIndx[i]]], 2.0);
    }
    varResult = varResult / (double) genSize;
    result = ((varResult <= EPSILON) ? FALSE : TRUE);
  }
  else {
    varResult = RF_nativeNaN;
    result = FALSE;
  }    
  if (mean != NULL) *mean = meanResult;
  if (variance != NULL) *variance = varResult;
  return(result);
}
char getVarianceClassicNoMiss(uint    repMembrSize,
                              uint   *repMembrIndx,
                              uint    nonMissMembrSize,
                              uint   *nonMissMembrIndx,
                              double *targetResponse,
                              double *mean,
                              double *variance) {
  uint i;
  double meanResult, varResult;
  char result;
  meanResult = 0.0;
  for (i = 1; i <= repMembrSize; i++) {
      meanResult += targetResponse[repMembrIndx[i]];
  }
  if (repMembrSize > 0) {
    meanResult = meanResult / (double) repMembrSize;  
  }
  else {
    meanResult = RF_nativeNaN;
  }
  varResult = 0.0;
  if(!RF_nativeIsNaN(meanResult)) {
    for (i = 1; i <= repMembrSize; i++) {
      varResult += pow(meanResult - targetResponse[repMembrIndx[i]], 2.0);
    }
    varResult = varResult / (double) repMembrSize;
    result = ((varResult <= EPSILON) ? FALSE : TRUE);
  }
  else {
    varResult = RF_nativeNaN;
    result = FALSE;
  }    
  if (mean != NULL) *mean = meanResult;
  if (variance != NULL) *variance = varResult;
  return(result);
}
char getVarianceDoublePass(uint    repMembrSize,
                           uint   *repMembrIndx,
                           uint    nonMissMembrSize,
                           uint   *nonMissMembrIndx,
                           double *targetResponse,
                           double *mean,
                           double *variance) {
  uint i;
  uint denom;
  double meanResult, varResult;
  char result;
  uint *genIndx;
  uint  genSize;
  if (nonMissMembrIndx == NULL) {
    genIndx = RF_identityMembershipIndex;
    genSize = repMembrSize;
  }
  else {
    genIndx = nonMissMembrIndx;
    genSize = nonMissMembrSize;
  }
  denom = 0;
  meanResult = 0.0;
  for (i = 1; i <= genSize; i++) {
    if(!RF_nativeIsNaN(targetResponse[repMembrIndx[genIndx[i]]])) {
      denom ++;
      meanResult += targetResponse[repMembrIndx[genIndx[i]]];
    }
  }
  if (denom > 0) {
    meanResult = meanResult / (double) denom;
  }
  else {
    meanResult = RF_nativeNaN;
  }
  if (mean != NULL) *mean = meanResult;
  varResult = 0.0;
  if(!RF_nativeIsNaN(meanResult)) {
    for (i = 1; i <= genSize; i++) {
      if(!RF_nativeIsNaN(targetResponse[repMembrIndx[genIndx[i]]])) {
        varResult += pow(meanResult - targetResponse[repMembrIndx[genIndx[i]]], 2.0);
      }
    }
    varResult = varResult / (double) denom;
    result = ((varResult <= EPSILON) ? FALSE : TRUE);
  }
  else {
    varResult = RF_nativeNaN;
    result = FALSE;
  }
  if (variance != NULL)  *variance = varResult;
  return(result);
}
char getVarianceSinglePass(uint    repMembrSize,
                           uint   *repMembrIndx,
                           uint    nonMissMembrSize,
                           uint   *nonMissMembrIndx,
                           double *targetResponse,
                           double *mean,
                           double *variance) {
  uint i;
  double meanResultOld, meanResult, varResultOld, varResult;
  char result;
  result = TRUE;
  meanResultOld = targetResponse[repMembrIndx[1]];
  varResultOld  = 0.0;
  for (i = 2; i <= repMembrSize; i++) {
    meanResult = meanResultOld + ((targetResponse[repMembrIndx[i]] - meanResultOld) / i);
    varResult  = varResultOld + ((targetResponse[repMembrIndx[i]] - meanResultOld) * (targetResponse[repMembrIndx[i]] - meanResult)); 
    meanResultOld = meanResult;
    varResultOld = varResult;
  }
  if (repMembrSize > 1) {
    varResultOld = varResultOld / (double) (repMembrSize - 1);
    result = ((varResultOld <= EPSILON) ? FALSE : TRUE);
  }
  *mean = meanResultOld;
  if (variance != NULL) *variance = varResultOld;
  return(result);
}
void restoreMeanResponse(uint treeID) {
  LeafLinkedObj *leafLinkedPtr;
  Terminal *parent;
  uint leaf;
  uint j;
  leafLinkedPtr = RF_leafLinkedObjHead[treeID] -> fwdLink;
  while (leafLinkedPtr != NULL) {
    parent = leafLinkedPtr -> termPtr;
    leaf = parent -> nodeID;
    if ((parent -> membrCount) > 0) {
      for (j = 1; j <= RF_rNonFactorCount; j++) {
        (parent -> meanResponse)[j] = RF_TN_REGR_ptr[treeID][leaf][j];
      }
    }
    else {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Zero node count encountered in restoreMeanResponse() in (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
    leafLinkedPtr = leafLinkedPtr -> fwdLink;
  }
}
void updateQuantileStream(char     mode,
                          uint     treeID) {
  char oobFlag, fullFlag, selectionFlag;
  Terminal ***termMembershipPtr;
  uint    *membershipIndex;
  uint     membershipSize;
  uint         **quantileStreamSize;
  LookUpInfo  ***quantileSearchTree;
  QuantileObj ***quantileHead;
  QuantileObj ***quantileTail;
  uint         **quantileLinkLength;
#ifdef _OPENMP
  omp_lock_t   *lockQNTptr;
#endif
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
  while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
    if (oobFlag == TRUE) {
      quantileStreamSize  = RF_oobQuantileStreamSize;
      quantileSearchTree  = RF_oobQuantileSearchTree;
      quantileHead        = RF_oobQuantileHead;
      quantileTail        = RF_oobQuantileTail;
      quantileLinkLength  = RF_oobQuantileLinkLength;
      membershipSize  = RF_oobSize[treeID];
      membershipIndex = RF_oobMembershipIndex[treeID];
#ifdef _OPENMP
      lockQNTptr      = RF_lockQNToens;
#endif
    }
    else {
      quantileStreamSize  = RF_fullQuantileStreamSize;
      quantileSearchTree  = RF_fullQuantileSearchTree;
      quantileHead        = RF_fullQuantileHead;
      quantileTail        = RF_fullQuantileTail;
      quantileLinkLength  = RF_fullQuantileLinkLength;
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
      lockQNTptr      = RF_lockQNTfens;
#endif
    }
    for (uint i = 1; i <= membershipSize; i++) {
      Terminal *parent;
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
        omp_set_lock(&(lockQNTptr[ii]));
#endif
        for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
          for (k = 1; k <= parent -> membrCount; k++) { 
            insertQuantileObj(&quantileStreamSize[j][ii],
                              &quantileHead[j][ii],
                              &quantileTail[j][ii],
                              &quantileLinkLength[j][ii],
                              RF_response[treeID][RF_rTargetNonFactor[j]][parent -> membrStream[k]],
                              &quantileSearchTree[j][ii]);
          }
        }
#ifdef _OPENMP
        omp_unset_lock(&(lockQNTptr[ii]));
#endif
      }
    }  
    if (oobFlag == TRUE) {
      oobFlag = FALSE;
    }
    else {
      fullFlag = FALSE;
    }
  }  
}
