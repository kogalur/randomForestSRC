
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
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetMultiClassProb() ENTRY ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & ENSB_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nClass proportion vector for (tree, leaf):  (%10d, %10d) \n", treeID, parent -> nodeID);
  ${trace.token}    for (j=1; j <= RF_rFactorCount; j++) {
  ${trace.token}      RF_nativePrint("\nFactor Index:  %10d", RF_rFactorIndex[j]);
  ${trace.token}      RF_nativePrint("\nMember Count:  %10d", parent -> membrCount);
  ${trace.token}      RF_nativePrint("\nClass ->  ");
  ${trace.token}      for (k=1; k <= RF_rFactorSize[j]; k++) {
  ${trace.token}        RF_nativePrint("%10d", k);
  ${trace.token}      }
  ${trace.token}      RF_nativePrint("\n");
  ${trace.token}      RF_nativePrint("%10d", j);
  ${trace.token}      for (k=1; k <= RF_rFactorSize[j]; k++) {
  ${trace.token}        RF_nativePrint("%10d", (parent -> multiClassProb)[j][k]);
  ${trace.token}      }
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetMultiClassProb() EXIT ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdateEnsembleMultiClass() ENTRY ...\n");
  ${trace.token}  }
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
          ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
          ${trace.token}    RF_nativePrint("\nCLAS OUTC_TYPE case no predicted value:  %10d \n", ii);
          ${trace.token}  }
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
    ${trace.token}  if (getTraceFlag(treeID) & SUMM_HGH_TRACE) {
    ${trace.token}    uint obsSize = (mode == RF_PRED) ? RF_fobservationSize : RF_observationSize;
    ${trace.token}    if (oobFlag == TRUE) {
    ${trace.token}      RF_nativePrint("\nOOB Ensemble calculations follow: \n");
    ${trace.token}    }
    ${trace.token}    else {
    ${trace.token}      if (fullFlag == TRUE) {
    ${trace.token}        RF_nativePrint("\nFULL Ensemble calculations follow: \n");
    ${trace.token}      }
    ${trace.token}    }
    ${trace.token}    RF_nativePrint("\nCLAS Numerator calculation: \n");
    ${trace.token}    RF_nativePrint("          ");
    ${trace.token}    for (uint i = 1; i <= obsSize; i++) {
    ${trace.token}      RF_nativePrint("%10d", i);
    ${trace.token}    }
    ${trace.token}    for (uint j = 1; j <= RF_rTargetFactorCount; j++) {
    ${trace.token}      RF_nativePrint("\n  Target Factor (idx, rsp): (%10d, %10d)", j, RF_rTargetFactor[j]);
    ${trace.token}      for (uint k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
    ${trace.token}        RF_nativePrint("\n%10d", k);
    ${trace.token}        for (uint i  = 1; i <= obsSize; i++) {
    ${trace.token}          RF_nativePrint("%10.4f", ensembleCLSnum[j][k][i]);
    ${trace.token}        }
    ${trace.token}      }
    ${trace.token}      RF_nativePrint("\n");
    ${trace.token}    }
    ${trace.token}    RF_nativePrint("\nCLAS Denominator calculation: \n");
    ${trace.token}    RF_nativePrint("          ");
    ${trace.token}    for (uint i = 1; i <= obsSize; i++) {
    ${trace.token}      RF_nativePrint("%10d", i);
    ${trace.token}    }
    ${trace.token}    RF_nativePrint("\n          ");
    ${trace.token}    for (uint i = 1; i <= obsSize; i++) {
    ${trace.token}      RF_nativePrint("%10d", (uint) ensembleDen[i]);
    ${trace.token}    }
    ${trace.token}    RF_nativePrint("\n");
    ${trace.token}    if (outcomeFlag == TRUE) {
    ${trace.token}      if (RF_optHigh & OPT_CSE) {
    ${trace.token}        RF_nativePrint("\nCLAS CSE Denominator calculation: \n");
    ${trace.token}        RF_nativePrint("                    ");
    ${trace.token}        for (uint i = 1; i <= obsSize; i++) {
    ${trace.token}          RF_nativePrint("%10d", i);
    ${trace.token}        }
    ${trace.token}        RF_nativePrint("\n                    ");
    ${trace.token}        for (uint i = 1; i <= obsSize; i++) {
    ${trace.token}          RF_nativePrint("%10d", RF_cseDENptr[i]);
    ${trace.token}        }
    ${trace.token}        RF_nativePrint("\nCLAS CSE Numerator calculation: \n");
    ${trace.token}        for (uint j = 1; j <= RF_rTargetFactorCount; j++) {
    ${trace.token}          RF_nativePrint("\n");
    ${trace.token}          RF_nativePrint("%20d", j);
    ${trace.token}          for (uint i = 1; i <= obsSize; i++) {
    ${trace.token}            RF_nativePrint("%10.4f", RF_cseNumCLSptr[j][i]);
    ${trace.token}          }
    ${trace.token}        }
    ${trace.token}        RF_nativePrint("\n");
    ${trace.token}      }
    ${trace.token}    }
    ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdateEnsembleMultiClass() EXIT ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetBrierScore() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}  uint j;
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
    ${trace.token}  if (getTraceFlag(0) & ENSB_LOW_TRACE) {
    ${trace.token}    
    ${trace.token}    RF_nativePrint("\nBrier Score Conditional Classification complete for level:  %10d ", against);
    ${trace.token}    RF_nativePrint("\nLevel     ");
    ${trace.token}    for (j=1; j <= RF_rFactorSize[RF_rFactorMap[rTarget]]; j++) {
    ${trace.token}      RF_nativePrint(" %10d", j);
    ${trace.token}    }
    ${trace.token}    RF_nativePrint("\n          ");
    ${trace.token}    for (j=1; j <= RF_rFactorSize[RF_rFactorMap[rTarget]]; j++) {
    ${trace.token}      RF_nativePrint(" %10.4f", cpv[j]);
    ${trace.token}    }
    ${trace.token}    RF_nativePrint("\n");
    ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(0) & ENSB_LOW_TRACE) {
  ${trace.token}    
  ${trace.token}    RF_nativePrint("\nBrier Score calculation complete:");
  ${trace.token}    RF_nativePrint("\nResult:                          %20.4f", result);
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}  }
  free_uivector(oaaResponse, 1, obsSize);
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetBrierScore() EXIT() ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetConditionalClassificationIndexGrow() ENTRY ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(0) & ENSB_LOW_TRACE) {
  ${trace.token}    
  ${trace.token}    RF_nativePrint("\nConditional Classification and error update complete:  ");
  ${trace.token}    RF_nativePrint("\nTarget:  %10d ", rTarget);
  ${trace.token}    RF_nativePrint("\nLevel     ");
  ${trace.token}    for (k=1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
  ${trace.token}      RF_nativePrint(" %10d", k);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n          ");
  ${trace.token}    for (k=1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
  ${trace.token}      RF_nativePrint(" %10.4f", cpv[k]);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetConditionalClassificationIndexGrow() EXIT() ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetConditionalClassificationIndexGrow() ENTRY ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(0) & ENSB_LOW_TRACE) {
  ${trace.token}    
  ${trace.token}    RF_nativePrint("\nConditional Classification and error update complete:  ");
  ${trace.token}    RF_nativePrint("\nTarget:  %10d ", rTarget);
  ${trace.token}    RF_nativePrint("\nLevel     ");
  ${trace.token}    for (k=1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
  ${trace.token}      RF_nativePrint(" %10d", k);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n          ");
  ${trace.token}    for (k=1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
  ${trace.token}      RF_nativePrint(" %10.4f", cpv[k]);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetConditionalClassificationIndexGrow() EXIT() ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetClassificationIndex() ENTRY ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(0) & ENSB_LOW_TRACE) {
  ${trace.token}    
  ${trace.token}    RF_nativePrint("\nPredicted Outcome used in Classification Index Calculations:  ");
  ${trace.token}    RF_nativePrint("\n        count     OOBcount     Response      Outcome");
  ${trace.token}    for (i=1; i <= size; i++) {
  ${trace.token}      RF_nativePrint("\n %12d %12d %12.4f %12.4f", i, denomCount[i], responsePtr[i], maxVote[i]);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}    
  ${trace.token}    RF_nativePrint("\nClassification and error update complete:  ");
  ${trace.token}    RF_nativePrint("\nCount of cumulative OOB count:         %20d", cumDenomCount);
  ${trace.token}    RF_nativePrint("\nResult:                                %20.4f", result);
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetClassificationIndex() EXIT() ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetGMeanIndexGrow() ENTRY ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(0) & ENSB_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nConfusion Matrix:");
  ${trace.token}    RF_nativePrint("\n                              Predicted L1         Predicted L2");
  ${trace.token}    RF_nativePrint("\n              True L1 %20d %20d", (uint) trueRate[1],  (uint) falseRate[1]);
  ${trace.token}    RF_nativePrint("\n              True L2 %20d %20d", (uint) falseRate[2], (uint) trueRate[2]);  
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}  }
  free_dvector(trueRate, 1, RF_rFactorSize[RF_rFactorMap[rTarget]]);
  free_dvector(falseRate, 1, RF_rFactorSize[RF_rFactorMap[rTarget]]);
  ${trace.token}  if (getTraceFlag(0) & ENSB_LOW_TRACE) {
  ${trace.token}    
  ${trace.token}    RF_nativePrint("\nPredicted Outcome used in Classification Index Calculations:  ");
  ${trace.token}    RF_nativePrint("\n        count     OOBcount     Response      Outcome");
  ${trace.token}    for (i=1; i <= size; i++) {
  ${trace.token}      RF_nativePrint("\n %12d %12.4f %12.4f %12.4f", i, denomCount[i], responsePtr[i], maxVote[i]);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}    
  ${trace.token}    RF_nativePrint("\nClassification and error update complete:  ");
  ${trace.token}    RF_nativePrint("\nCount of cumulative OOB count:         %20d", cumDenomCount);
  ${trace.token}    RF_nativePrint("\nResult:                                %20.4f", result);
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetGMeanIndexGrow() EXIT() ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetGMeanIndexGrow() ENTRY ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(0) & ENSB_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nConfusion Matrix:");
  ${trace.token}    RF_nativePrint("\n          ");
  ${trace.token}    for (i=1; i <= RF_rFactorSizeTest[RF_rFactorMap[rTarget]]; i++) {
  ${trace.token}      RF_nativePrint(" %10d", i);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\ntrue:     ");
  ${trace.token}    for (i=1; i <= RF_rFactorSizeTest[RF_rFactorMap[rTarget]]; i++) {
  ${trace.token}      RF_nativePrint(" %10d", (uint) trueRate[i]);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\nfalse:    ");
  ${trace.token}    for (i=1; i <= RF_rFactorSizeTest[RF_rFactorMap[rTarget]]; i++) {
  ${trace.token}      RF_nativePrint(" %10d", (uint) falseRate[i]);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}  }
  free_dvector(trueRate, 1, RF_rFactorSizeTest[RF_rFactorMap[rTarget]]);
  free_dvector(falseRate, 1, RF_rFactorSizeTest[RF_rFactorMap[rTarget]]);
  ${trace.token}  if (getTraceFlag(0) & ENSB_LOW_TRACE) {
  ${trace.token}    
  ${trace.token}    RF_nativePrint("\nPredicted Outcome used in Classification Index Calculations:  ");
  ${trace.token}    RF_nativePrint("\n        count     OOBcount     Response      Outcome");
  ${trace.token}    for (i=1; i <= size; i++) {
  ${trace.token}      RF_nativePrint("\n %12d %12.4f %12.4f %12.4f", i, denomCount[i], responsePtr[i], maxVote[i]);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}    
  ${trace.token}    RF_nativePrint("\nClassification and error update complete:  ");
  ${trace.token}    RF_nativePrint("\nCount of cumulative OOB count:         %20d", cumDenomCount);
  ${trace.token}    RF_nativePrint("\nResult:                                %20.4f", result);
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetGMeanIndexGrow() EXIT() ...\n");
  ${trace.token}  }
  return result;
}
void restoreMultiClassProb(uint treeID) {
  LeafLinkedObj *leafLinkedPtr;
  Terminal *parent;
  uint leaf;
  uint j, k;
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nrestoreMultiClassProb() ENTRY ...\n");
  ${trace.token}  }
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
    ${trace.token}  if (getTraceFlag(treeID) & SUMM_HGH_TRACE) {
    ${trace.token}    RF_nativePrint("\nClass proportion vector for (tree, leaf):  (%10d, %10d) \n", treeID, leaf);
    ${trace.token}    for (j=1; j <= RF_rFactorCount; j++) {
    ${trace.token}      RF_nativePrint("\nFactor Index:  %10d ", RF_rFactorIndex[j]);
    ${trace.token}      RF_nativePrint("\nMember Count:  %10d", parent -> membrCount);
    ${trace.token}      RF_nativePrint("\nClass ->  ");
    ${trace.token}      for (k=1; k <= RF_rFactorSize[j]; k++) {
    ${trace.token}        RF_nativePrint("%10d", k);
    ${trace.token}      }
    ${trace.token}      RF_nativePrint("\n");
    ${trace.token}      RF_nativePrint("%10d", j);
    ${trace.token}      for (k=1; k <= RF_rFactorSize[j]; k++) {
    ${trace.token}        RF_nativePrint("%10d", (parent -> multiClassProb)[j][k]);
    ${trace.token}      }
    ${trace.token}    }
    ${trace.token}  }
  }
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nrestoreMultiClassProb() EXIT ...\n");
  ${trace.token}  }
}
void getMaxVote(uint     size,
                uint     rTarget,
                double **outcomeCLS,
                double  *denomCount,
                double  *maxVote) {
  uint i,k;
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetMaxVote() ENTRY() ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(0) & SUMM_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\n  Target Factor: %10d", rTarget);
  ${trace.token}    for (int k=1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
  ${trace.token}      RF_nativePrint("\n%10d", k);
  ${trace.token}      for (i = 1; i <= size; i++) {
  ${trace.token}        RF_nativePrint("%10.4f", outcomeCLS[k][i]);
  ${trace.token}      }
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}    RF_nativePrint("\nMax Vote calculation: \n");
  ${trace.token}    RF_nativePrint("          ");
  ${trace.token}    for (i = 1; i <= size; i++) {
  ${trace.token}      RF_nativePrint("%10d", i);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}    RF_nativePrint("          ");
  ${trace.token}    for (i = 1; i <= size; i++) {
  ${trace.token}      RF_nativePrint("%10.4f", maxVote[i]);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetMaxVote() EXIT() ...\n");
  ${trace.token}  }
}
