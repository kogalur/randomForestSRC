
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "survivalE.h"
#include "impute.h"
#include "nrutil.h"
#include "error.h"
void updateEnsembleSurvival(char mode,
                            uint treeID,
                            char normalizationFlag) {
  char oobFlag, fullFlag, outcomeFlag;
  Terminal ***termMembershipPtr;
  uint    *membershipIndex;
  uint     membershipSize;
  double  **ensembleMRTptr;
  double ***ensembleSRGnum;
  double ***ensembleCIFnum;
  double  **ensembleSRVnum;
  double  **ensembleMRTnum;
  double   *ensembleDen;
#ifdef _OPENMP
  omp_lock_t   *lockDENptr;
#endif
  ensembleSRGnum = NULL;  
  ensembleCIFnum = NULL;  
  ensembleSRVnum = NULL;  
  ensembleMRTnum = NULL;  
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
      ensembleMRTptr = RF_oobEnsembleMRTptr;        
      ensembleSRGnum = RF_oobEnsembleSRGnum;
      ensembleMRTnum = RF_oobEnsembleMRTnum;
      ensembleSRVnum = RF_oobEnsembleSRVnum;
      ensembleCIFnum = RF_oobEnsembleCIFnum;
      ensembleDen    = RF_oobEnsembleDen;
      membershipSize  = RF_oobSize[treeID];
      membershipIndex = RF_oobMembershipIndex[treeID];
#ifdef _OPENMP
      lockDENptr      = RF_lockDENoens;
#endif
    }
    else {
      ensembleMRTptr = RF_fullEnsembleMRTptr;        
      ensembleSRGnum = RF_fullEnsembleSRGnum;
      ensembleMRTnum = RF_fullEnsembleMRTnum;        
      ensembleSRVnum = RF_fullEnsembleSRVnum;
      ensembleCIFnum = RF_fullEnsembleCIFnum;
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
        ensembleDen[ii] ++;
        if (outcomeFlag == TRUE) {
          if (RF_opt & OPT_VIMP) {
            RF_blkEnsembleDen[ii] ++;
          }
        }
        if (!(RF_opt & OPT_COMP_RISK)) {
          for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
            ensembleSRGnum[1][k][ii] += parent -> nelsonAalen[k];
            ensembleSRVnum[k][ii] += parent -> survival[k];
          }
          ensembleMRTnum[1][ii] += parent -> mortality[1];
          if (outcomeFlag == TRUE) {
            if (RF_opt & OPT_VIMP) {
              RF_blkEnsembleMRTnum[1][ii] += parent -> mortality[1];
            }
          }
          if (outcomeFlag && normalizationFlag) {
            ensembleMRTptr[1][ii] = ensembleMRTnum[1][ii] / ensembleDen[ii];
          }
        }
        else {
          for (j = 1; j <= RF_eventTypeSize; j++) {
            for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
              ensembleSRGnum[j][k][ii] += parent -> CSH[j][k];
              ensembleCIFnum[j][k][ii] += parent -> CIF[j][k];
            }
            ensembleMRTnum[j][ii] += parent -> mortality[j];
            if (outcomeFlag == TRUE) {
              if (RF_opt & OPT_VIMP) {
                RF_blkEnsembleMRTnum[j][ii] += parent -> mortality[j];
              }
            }
            if (outcomeFlag && normalizationFlag) {
              ensembleMRTptr[j][ii] = ensembleMRTnum[j][ii] / ensembleDen[ii];
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
void getEnsembleMortalityCR(char      mode,
                            uint      treeID,
                            uint      obsSize,
                            double  **ensembleMRTptr,
                            double   *ensembleDen,
                            double  **cMortality) {
  uint i, j;
  for (i = 1; i <= obsSize; i++) {
    if (ensembleDen[i] != 0) {
      for (j = 1; j <= RF_eventTypeSize; j ++) {
        cMortality[j][i] = ensembleMRTptr[j][i] / ensembleDen[i];
      }
    }
    else {
      for (j = 1; j <= RF_eventTypeSize; j ++) {
        cMortality[j][i] = RF_nativeNaN;
      }
    }
  }
}
void getEnsembleMortality(char      mode,
                          uint      treeID,
                          uint      obsSize,
                          double  **ensembleMRTptr,
                          double   *ensembleDen,
                          double   *mortality) {
  uint i;
  for (i = 1; i <= obsSize; i++) {
    if (ensembleDen[i] != 0) {
      mortality[i] = ensembleMRTptr[1][i] / ensembleDen[i];
    }
    else {
      mortality[i] = RF_nativeNaN;
    }
  }
}
void getConditionalConcordanceArrays(uint     j,
                                     double  *timePtr,
                                     double  *statusPtr,
                                     double  *mortalityPtr,
                                     double  *genericEnsembleDenPtr,
                                     uint    *meIndividualSize,
                                     uint   **eIndividual,
                                     double  *subsettedTime,
                                     double  *subsettedStatus,
                                     double  *subsettedMortality,
                                     double  *subsettedEnsembleDen) {
  uint i;
  if (!(RF_opt & OPT_COMP_RISK)) {
    RF_nativePrint("\nRF-SRC:  *** ERROR *** ");
    RF_nativePrint("\nRF-SRC:  Attempt to update event type subsets in a non-CR analysis.");
    RF_nativePrint("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  for (i = 1; i <= meIndividualSize[j]; i++) {
    subsettedTime[i]        = timePtr[eIndividual[j][i]];
    subsettedStatus[i]      = statusPtr[eIndividual[j][i]];
    subsettedMortality[i]   = mortalityPtr[eIndividual[j][i]];
    subsettedEnsembleDen[i] = genericEnsembleDenPtr[eIndividual[j][i]];
  }
}
double getConcordanceIndex(int     polarity,
                           uint    size,
                           double *timePtr,
                           double *statusPtr,
                           double *predictedOutcome,
                           double *denCount) {
  uint i,j;
  long long concordancePairSize;
  long long concordanceWorseCount;
  double result;
  concordancePairSize = concordanceWorseCount = 0;
  for (i=1; i < size; i++) {
    for (j=i+1; j <= size; j++) {
      if (denCount[i] != 0  && denCount[j] != 0) {
        if ( ((timePtr[i] - timePtr[j] > EPSILON) && (statusPtr[j] > 0)) ||
             ((fabs(timePtr[i] - timePtr[j]) <= EPSILON) && (statusPtr[j] > 0) && (statusPtr[i] == 0)) ) {
          concordancePairSize += 2;
          if (predictedOutcome[j] - predictedOutcome[i] > EPSILON) {
            concordanceWorseCount += 2;
          }
          else if (fabs(predictedOutcome[j] - predictedOutcome[i]) <= EPSILON) {
            concordanceWorseCount += 1;
          }
        }
        else if ( ((timePtr[j] - timePtr[i]) > EPSILON  && (statusPtr[i] > 0)) ||
                  ((fabs(timePtr[j] - timePtr[i]) <= EPSILON)  && (statusPtr[i] > 0) && (statusPtr[j] == 0)) ) {
          concordancePairSize += 2;
          if ( predictedOutcome[i] - predictedOutcome[j] > EPSILON ) {
            concordanceWorseCount += 2;
          }
          else if (fabs(predictedOutcome[i] - predictedOutcome[j]) <= EPSILON) {
            concordanceWorseCount += 1;
          }
        }
        else if ( (fabs(timePtr[i]- timePtr[j]) <= EPSILON) && (statusPtr[i] > 0) && (statusPtr[j] > 0) ) {
          concordancePairSize += 2;
          if (fabs(predictedOutcome[i] - predictedOutcome[j]) < EPSILON) {
            concordanceWorseCount += 2;
          }
          else {
            concordanceWorseCount += 1;
          }
        }
      }  
    }  
  }  
  if (concordancePairSize == 0) {
    result = RF_nativeNaN;
  }
  else {
    result = 1.0 - ((double) concordanceWorseCount / (double) concordancePairSize);
  }
  return result;
}
double getConcordanceIndexNew(int     polarity,
                              uint    size,
                              double *timePtr,
                              double *statusPtr,
                              double *predicted,
                              double *denCount) {
  uint i,j;
  long long concordancePairSize;
  long long concordanceBetterCount;
  double result;
  uint *timePtrIndxx = uivector(1, size);
  indexx(size, timePtr, timePtrIndxx);
  uint *riskSetSize = uivector(1, size);
  double *predictedOrdByTime = dvector(1, size);
  double *statusOrdByTime    = dvector(1, size);
  double *denCountOrdByTime  = dvector(1, size);
  for (i = 1; i <= size; i++) {
    riskSetSize[i] = size - i;
    predictedOrdByTime[i] = predicted[timePtrIndxx[i]];
    statusOrdByTime[i]    = statusPtr[timePtrIndxx[i]];
    denCountOrdByTime[i]  = denCount[timePtrIndxx[i]];
  }
  uint *predictedIndxx = uivector(1, size);
  indexx(size, predictedOrdByTime, predictedIndxx);
  uint *predictedRank = uivector(1, size);
  for (i = 1; i <= size; i++) {
    predictedRank[predictedIndxx[i]] = size - i;
  }
  concordancePairSize = concordanceBetterCount = 0;
  for (i = 1; i <= size; i++) {
    if (statusOrdByTime[i] > 0) {
      for (j = i+1; j <= size; j++) {
        if (denCountOrdByTime[i] != 0  && denCountOrdByTime[j] != 0) {
          if (predictedRank[i] > predictedRank[j]) {
            concordanceBetterCount += 1;
          }
          else {
          }
        }
        else {
        }
      }
      concordancePairSize += riskSetSize[i];
    }
  }
  if (concordancePairSize == 0) {
    result = RF_nativeNaN;
  }
  else {
    result = ((double) concordanceBetterCount / (double) concordancePairSize);
  }
  free_uivector(predictedRank, 1, size);
  free_uivector(predictedIndxx, 1, size);
  free_uivector(timePtrIndxx, 1, size);
  free_dvector(predictedOrdByTime, 1, size);
  free_dvector(statusOrdByTime,    1, size);
  free_dvector(denCountOrdByTime,  1, size);
  free_uivector(riskSetSize, 1, size);
  return result;
}
void getCRPerformance (char     mode,
                       uint     obsSize,
                       double **responsePtr,
                       double **yearsLost,
                       double  *denom,
                       double  *performanceVector) {
  uint   mRecordSize;
  int  **mpSign;
  uint  *mRecordIndex;
  uint  *meIndividualSize;
  uint **eIndividual;
  double concordanceIndex;
  uint j;
  if (!(RF_opt & OPT_COMP_RISK)) {
    RF_nativePrint("\nRF-SRC:  *** ERROR *** ");
    RF_nativePrint("\nRF-SRC:  Attempt at conditional performance updates in a non-CR analysis.");
    RF_nativePrint("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  if (RF_mStatusSize > 0) {
    switch (mode) {
    case RF_PRED:
      mRecordSize = RF_fmRecordSize;
      mpSign = RF_fmpSign;
      mRecordIndex = RF_fmRecordIndex;
      break;
    default:
      mRecordSize = RF_mRecordSize;
      mpSign = RF_mpSign;
      mRecordIndex = RF_mRecordIndex;
      break;
    }
    meIndividualSize  = uivector(1, RF_eventTypeSize);
    eIndividual = (uint **) new_vvector(1, RF_eventTypeSize, NRUTIL_UPTR);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      eIndividual[j] = uivector(1, RF_eIndividualSize[j] + RF_mStatusSize + 1);
    }
    updateEventTypeSubsets(responsePtr[RF_statusIndex], mRecordSize, mpSign, mRecordIndex, meIndividualSize, eIndividual);
  }
  else {
    meIndividualSize  = RF_eIndividualSize;
    eIndividual = RF_eIndividualIn;
  }
  double *subsettedTime      = dvector(1, obsSize);
  double *subsettedStatus    = dvector(1, obsSize);
  double *subsettedMortality = dvector(1, obsSize);
  double *subsettedEnsembleDen = dvector(1, obsSize);
  for (j = 1; j <= RF_eventTypeSize; j++) {
    getConditionalConcordanceArrays(j,
                                    responsePtr[RF_timeIndex],
                                    responsePtr[RF_statusIndex],
                                    yearsLost[j],
                                    denom,
                                    meIndividualSize,
                                    eIndividual,
                                    subsettedTime,
                                    subsettedStatus,
                                    subsettedMortality,
                                    subsettedEnsembleDen);
    concordanceIndex = getConcordanceIndex(1,
                                           meIndividualSize[j],
                                           subsettedTime,
                                           subsettedStatus,
                                           subsettedMortality,
                                           subsettedEnsembleDen);
    if (RF_nativeIsNaN(concordanceIndex)) {
      performanceVector[j] = RF_nativeNaN;
    }
    else {
      performanceVector[j] = concordanceIndex;
    }
  }
  if (RF_mStatusSize > 0) {
    free_uivector(meIndividualSize, 1, RF_eventTypeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      free_uivector(eIndividual[j], 1, RF_eIndividualSize[j] + RF_mStatusSize + 1);
    }
    free_new_vvector(eIndividual, 1, RF_eventTypeSize, NRUTIL_UPTR);
  }
  free_dvector(subsettedTime, 1, obsSize);
  free_dvector(subsettedStatus, 1, obsSize);
  free_dvector(subsettedMortality, 1, obsSize);
  free_dvector(subsettedEnsembleDen, 1, obsSize);
}
uint getTimeInterestIndex(double *array, uint length, double value) {
  uint low, high, mid, result;
  if (value <= array[1]) {
    result = 1;
  }
  else if (value > array[length]) {
    result = length + 1;
  }
  else {
    low  = 1;
    high = length;;
    while (low < high) {
      mid  = (low + high) >> 1;
      if (value > array[mid]) {
        if (low == mid) {
          low = high;
        }
        else {
          low = mid;
        }
      }
      else {
        if (low == mid) {
          low = high;
        }
        else {
          high = mid;
        }
      }
    }
    result = high;
  }
  return result;
}
