
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
                                     double  *weight,
                                     uint    *meIndividualSize,
                                     uint   **eIndividual,
                                     double  *subsettedTime,
                                     double  *subsettedStatus,
                                     double  *subsettedMortality,
                                     double  *subsettedEnsembleDen,
                                     double  *subsettedWeight) {
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
  if (weight != NULL) {
    for (i = 1; i <= meIndividualSize[j]; i++) {
      subsettedWeight[i] = weight[eIndividual[j][i]];
    }
  }
}
double getConcordanceIndex(int     fastAction,
                           uint    size, 
                           double *timePtr, 
                           double *statusPtr, 
                           double *predicted,
                           double *denom,
                           double *weight) {
  double  (*getConcordanceIndexActual) (uint, double*, double*, double*, double*);
  char    fastFlag; 
  double *genericWeight;
  double  result;
  uint    i, j;
  getConcordanceIndexActual = NULL;  
  if (fastAction == 1) {
    fastFlag = TRUE;
  }
  else if (fastAction == 0) {
    fastFlag = FALSE;
  }
  else {
    fastFlag = FALSE;
    if (size > 500) {
      fastFlag = TRUE;
    }
    else {
      j = 0;
      for (i = 1; i <= size; i++) {
        if (statusPtr[i] > 0 && denom[i] > 0) {
          j++;
        }
      }
      if (j > 250) {
        fastFlag = TRUE;
      }
    }
  }
  if (weight == NULL) {
    genericWeight = denom;
  }
  else {
    genericWeight = dvector(1, size);
    for (i = 1; i <= size; i++) {
      if (denom[i] > 0) {
        genericWeight[i] = weight[i];
      }
      else {
        genericWeight[i] = 0;
      }
    }
  }
  if ((weight != NULL) &&  (fastFlag == TRUE)) {
    getConcordanceIndexActual = & getConcordanceIndexUnoFenwick;
  }
  else if ((weight == NULL) &&  (fastFlag == TRUE)) {
    getConcordanceIndexActual = & getConcordanceIndexFenwick;    
  }
  else if ((weight != NULL) &&  (fastFlag == FALSE)) {
    getConcordanceIndexActual = & getConcordanceIndexUno;    
  }
  else if ((weight == NULL) &&  (fastFlag == FALSE)) {
    getConcordanceIndexActual = & getConcordanceIndexOriginal;    
  }
  result = getConcordanceIndexActual(size,
                                     timePtr,
                                     statusPtr,
                                     predicted,
                                     genericWeight);
  if (weight == NULL) {
  }
  else {
    free_dvector(genericWeight, 1, size);
  }
  return result;
}
double getConcordanceIndexOriginal(uint    size,
                                   double *timePtr,
                                   double *statusPtr,
                                   double *predictedOutcome,
                                   double *denom) {
  uint i,j;
  double concordancePairCount;
  double concordanceCount;
  double result;
  concordancePairCount = concordanceCount = 0;
  for (i=1; i < size; i++) {
    for (j=i+1; j <= size; j++) {
      if (denom[i] != 0  && denom[j] != 0) {
        if ( ((timePtr[i] - timePtr[j] > EPSILON) && (statusPtr[j] > 0)) ||
             ((fabs(timePtr[i] - timePtr[j]) <= EPSILON) && (statusPtr[j] > 0) && (statusPtr[i] == 0)) ) {
          concordancePairCount += 2;
          if (predictedOutcome[j] - predictedOutcome[i] > EPSILON) {
            concordanceCount += 2;
          }
          else if (fabs(predictedOutcome[j] - predictedOutcome[i]) <= EPSILON) {
            concordanceCount += 1;
          }
        }
        else if ( ((timePtr[j] - timePtr[i]) > EPSILON  && (statusPtr[i] > 0)) ||
                  ((fabs(timePtr[j] - timePtr[i]) <= EPSILON)  && (statusPtr[i] > 0) && (statusPtr[j] == 0)) ) {
          concordancePairCount += 2;
          if ( predictedOutcome[i] - predictedOutcome[j] > EPSILON ) {
            concordanceCount += 2;
          }
          else if (fabs(predictedOutcome[i] - predictedOutcome[j]) <= EPSILON) {
            concordanceCount += 1;
          }
        }
        else if ( (fabs(timePtr[i]- timePtr[j]) <= EPSILON) && (statusPtr[i] > 0) && (statusPtr[j] > 0) ) {
          concordancePairCount += 2;
          if (fabs(predictedOutcome[i] - predictedOutcome[j]) < EPSILON) {
            concordanceCount += 2;
          }
          else {
            concordanceCount += 1;
          }
        }
      }  
    }  
  }  
  if (concordancePairCount == 0) {
    result = RF_nativeNaN;
  }
  else {
    result = 1.0 - (concordanceCount / concordancePairCount);
  }
  return result;
}
double getConcordanceIndexUno(uint    size,
                              double *timePtr,
                              double *statusPtr,
                              double *predictedOutcome,
                              double *weight) {
  uint i, j;
  double concordancePairWeight;
  double concordanceWeight;
  double w;
  double result;
  concordancePairWeight = concordanceWeight = 0.0;
  for (i = 1; i < size; i++) {
    for (j = i + 1; j <= size; j++) {
      if (weight[i] != 0 && weight[j] != 0) {
        if ((timePtr[i] - timePtr[j] > EPSILON) && (statusPtr[j] > 0) && (weight[j] > 0.0)) {
          w = weight[j] * 2;
          concordancePairWeight += w;
          if (predictedOutcome[j] - predictedOutcome[i] > EPSILON) {
            concordanceWeight += w;
          }
          else if (fabs(predictedOutcome[j] - predictedOutcome[i]) <= EPSILON) {
            concordanceWeight += 0.5 * w;
          }
        }
        else if ((timePtr[j] - timePtr[i] > EPSILON) && (statusPtr[i] > 0) && (weight[i] > 0.0)) {
          w = weight[i] * 2;
          concordancePairWeight += w;
          if (predictedOutcome[i] - predictedOutcome[j] > EPSILON) {
            concordanceWeight += w;
          }
          else if (fabs(predictedOutcome[i] - predictedOutcome[j]) <= EPSILON) {
            concordanceWeight += 0.5 * w;
          }
        }
        else {
          if ((fabs(timePtr[i] - timePtr[j]) <= EPSILON) && (statusPtr[j] > 0) && (statusPtr[i] == 0)) {
            w = weight[j] * 2;
            concordancePairWeight += w;
            if (predictedOutcome[j] - predictedOutcome[i] > EPSILON) {
              concordanceWeight += w;
            }
            else if (fabs(predictedOutcome[j] - predictedOutcome[i]) <= EPSILON) {
              concordanceWeight += 0.5 * w;
            }
          }
          else if ((fabs(timePtr[j] - timePtr[i]) <= EPSILON)  && (statusPtr[i] > 0) && (statusPtr[j] == 0)) {
            w = weight[i] * 2;
            concordancePairWeight += w;
            if (predictedOutcome[i] - predictedOutcome[j] > EPSILON) {
              concordanceWeight += w;
            }
            else if (fabs(predictedOutcome[i] - predictedOutcome[j]) <= EPSILON) {
              concordanceWeight += 0.5 * w;
            }
          }
          else if ((fabs(timePtr[j] - timePtr[i]) <= EPSILON)  && (statusPtr[i] > 0) && (statusPtr[j] > 0)) {
            w = (weight[i] + weight[j]);
            concordancePairWeight += w;
            if (fabs(predictedOutcome[i] - predictedOutcome[j]) < EPSILON) {
              concordanceWeight += w;
            }
            else { 
              concordanceWeight += 0.5 * w;
            }
          }
        }
      }
    }
  }
  if (concordancePairWeight <= 0.0) {
    result = RF_nativeNaN;
  }
  else {
    result = 1.0 - (concordanceWeight / concordancePairWeight);
  }
  return result;
}
static inline void rfsrc_bitAdd(uint *bit, uint m, uint idx, uint delta) {
  while (idx <= m) {
    bit[idx] += delta;
    idx += idx & (uint)(-((int) idx));
  }
}
static inline uint rfsrc_bitSum(const uint *bit, uint idx) {
  uint s = 0;
  while (idx > 0) {
    s += bit[idx];
    idx -= idx & (uint)(-((int) idx));
  }
  return s;
}
static inline void rfsrc_bitAddD(double *bit, uint m, uint idx, double delta) {
  while (idx <= m) {
    bit[idx] += delta;
    idx += idx & (uint)(-((int) idx));
  }
}
static inline double rfsrc_bitSumD(const double *bit, uint idx) {
  double s = 0.0;
  while (idx > 0) {
    s += bit[idx];
    idx -= idx & (uint)(-((int) idx));
  }
  return s;
}
double getConcordanceIndexFenwick(uint    size,
                                  double *timePtr,
                                  double *statusPtr,
                                  double *predictedOutcome,
                                  double *denom) {
  uint i;
  double result;
  uint n = 0;
  for (i = 1; i <= size; i++) {
    if (denom[i] != 0) {
      n++;
    }
  }
  if (n < 2) {
    return RF_nativeNaN;
  }
  double *t = dvector(1, n);
  double *s = dvector(1, n);
  double *p = dvector(1, n);
  uint kk = 0;
  for (i = 1; i <= size; i++) {
    if (denom[i] != 0) {
      kk++;
      t[kk] = timePtr[i];
      s[kk] = statusPtr[i];
      p[kk] = predictedOutcome[i];
    }
  }
  uint *pIndxx = uivector(1, n);
  indexx(n, p, pIndxx);
  uint *rank = uivector(1, n);
  uint m = 1;
  double last = p[pIndxx[1]];
  rank[pIndxx[1]] = m;
  for (i = 2; i <= n; i++) {
    double cur = p[pIndxx[i]];
    if ((cur - last) > EPSILON) {
      m++;
      last = cur;
    }
    rank[pIndxx[i]] = m;
  }
  uint *tIndxx = uivector(1, n);
  indexx(n, t, tIndxx);
  uint *bit = uivector(1, m);
  for (i = 1; i <= m; i++) {
    bit[i] = 0;
  }
  uint bitTotal = 0;
  uint *rankCount = uivector(1, m);
  for (i = 1; i <= m; i++) {
    rankCount[i] = 0;
  }
  uint *touched = uivector(1, n);
  uint *eventPos = uivector(1, n);
  double concordancePairCount = 0.0;
  double concordanceCount     = 0.0;
  int pos = (int) n;
  while (pos >= 1) {
    double curTime = t[tIndxx[pos]];
    int start = pos;
    while (start > 1 && (fabs(t[tIndxx[start - 1]] - curTime) <= EPSILON)) {
      start--;
    }
    for (int q = start; q <= pos; q++) {
      uint j = tIndxx[q];
      if (s[j] == 0) {
        rfsrc_bitAdd(bit, m, rank[j], 1);
        bitTotal++;
      }
    }
    uint d = 0;
    for (int q = start; q <= pos; q++) {
      uint j = tIndxx[q];
      if (s[j] > 0) {
        eventPos[++d] = j;
      }
    }
    if (d > 0) {
      uint touchedSize = 0;
      for (uint e = 1; e <= d; e++) {
        uint r = rank[eventPos[e]];
        if (rankCount[r] == 0) {
          touched[++touchedSize] = r;
        }
        rankCount[r]++;
      }
      double tiePairs = 0.0;
      for (uint u = 1; u <= touchedSize; u++) {
        uint r = touched[u];
        uint c = rankCount[r];
        if (c >= 2) {
          tiePairs += ((double) c * (double) (c - 1)) / 2.0;
        }
        rankCount[r] = 0;
      }
      if (d >= 2) {
        concordancePairCount += (double) d * (double) (d - 1);        
        concordanceCount     += 0.5 * (double) d * (double) (d - 1)   
          + tiePairs;                           
      }
    }
    if (bitTotal > 0) {
      for (uint e = 1; e <= d; e++) {
        uint j = eventPos[e];
        uint r = rank[j];
        uint less = (r > 1) ? rfsrc_bitSum(bit, r - 1) : 0;
        uint leq  = rfsrc_bitSum(bit, r);
        uint eq   = leq - less;
        concordancePairCount += 2.0 * (double) bitTotal;
        concordanceCount     += 2.0 * (double) less + 1.0 * (double) eq;
      }
    }
    for (uint e = 1; e <= d; e++) {
      uint j = eventPos[e];
      rfsrc_bitAdd(bit, m, rank[j], 1);
      bitTotal++;
    }
    pos = start - 1;
  }
  free_uivector(eventPos, 1, n);
  free_uivector(touched, 1, n);
  free_uivector(rankCount, 1, m);
  free_uivector(bit, 1, m);
  free_uivector(tIndxx, 1, n);
  free_uivector(rank, 1, n);
  free_uivector(pIndxx, 1, n);
  free_dvector(p, 1, n);
  free_dvector(s, 1, n);
  free_dvector(t, 1, n);
  if (concordancePairCount == 0.0) {
    result = RF_nativeNaN;
  }
  else {
    result = 1.0 - (concordanceCount / concordancePairCount);
  }
  return result;
}
double getConcordanceIndexUnoFenwick(uint    size,
                                     double *timePtr,
                                     double *statusPtr,
                                     double *predictedOutcome,
                                     double *weight) {
  uint i;
  double result;
  uint n = 0;
  for (i = 1; i <= size; i++) {
    if (weight[i] != 0.0) {
      n++;
    }
  }
  if (n < 2) {
    return RF_nativeNaN;
  }
  double *t = dvector(1, n);
  double *s = dvector(1, n);
  double *p = dvector(1, n);
  double *w = dvector(1, n);
  uint kk = 0;
  for (i = 1; i <= size; i++) {
    if (weight[i] != 0.0) {
      kk++;
      t[kk] = timePtr[i];
      s[kk] = statusPtr[i];
      p[kk] = predictedOutcome[i];
      w[kk] = weight[i];
    }
  }
  uint *pIndxx = uivector(1, n);
  indexx(n, p, pIndxx);
  uint *rank = uivector(1, n);
  uint m = 1;
  double last = p[pIndxx[1]];
  rank[pIndxx[1]] = m;
  for (i = 2; i <= n; i++) {
    double cur = p[pIndxx[i]];
    if ((cur - last) > EPSILON) {
      m++;
      last = cur;
    }
    rank[pIndxx[i]] = m;
  }
  uint *tIndxx = uivector(1, n);
  indexx(n, t, tIndxx);
  uint *bit = uivector(1, m);
  for (i = 1; i <= m; i++) {
    bit[i] = 0;
  }
  uint bitTotal = 0;
  uint   *rankCount     = uivector(1, m);
  double *rankWeightSum = dvector(1, m);
  for (i = 1; i <= m; i++) {
    rankCount[i] = 0;
    rankWeightSum[i] = 0.0;
  }
  uint *touched  = uivector(1, n);
  uint *eventPos = uivector(1, n);
  double concordancePairWeight = 0.0;  
  double concordanceWeight     = 0.0;  
  int pos = (int) n;
  while (pos >= 1) {
    double curTime = t[tIndxx[pos]];
    int start = pos;
    while (start > 1 && (fabs(t[tIndxx[start - 1]] - curTime) <= EPSILON)) {
      start--;
    }
    for (int q = start; q <= pos; q++) {
      uint j = tIndxx[q];
      if (s[j] == 0) {
        rfsrc_bitAdd(bit, m, rank[j], 1);
        bitTotal++;
      }
    }
    uint d = 0;
    for (int q = start; q <= pos; q++) {
      uint j = tIndxx[q];
      if (s[j] > 0) {
        eventPos[++d] = j;
      }
    }
    if (d > 0) {
      double sumW = 0.0;
      uint touchedSize = 0;
      for (uint e = 1; e <= d; e++) {
        uint j = eventPos[e];
        sumW += w[j];
        uint r = rank[j];
        if (rankCount[r] == 0) {
          touched[++touchedSize] = r;
        }
        rankCount[r]++;
        rankWeightSum[r] += w[j];
      }
      if (d >= 2) {
        double denomTie = (double) (d - 1) * sumW;
        double tieMass  = 0.0;
        for (uint u = 1; u <= touchedSize; u++) {
          uint r = touched[u];
          uint c = rankCount[r];
          if (c >= 2) {
            tieMass += (double) (c - 1) * rankWeightSum[r]; 
          }
          rankCount[r] = 0;
          rankWeightSum[r] = 0.0;
        }
        concordancePairWeight += denomTie;
        concordanceWeight     += 0.5 * denomTie + 0.5 * tieMass;
      }
      else {
        for (uint u = 1; u <= touchedSize; u++) {
          uint r = touched[u];
          rankCount[r] = 0;
          rankWeightSum[r] = 0.0;
        }
      }
    }
    if (bitTotal > 0) {
      for (uint e = 1; e <= d; e++) {
        uint j = eventPos[e];
        if (w[j] > 0.0) {
          uint r = rank[j];
          uint less = (r > 1) ? rfsrc_bitSum(bit, r - 1) : 0;
          uint leq  = rfsrc_bitSum(bit, r);
          uint eq   = leq - less;
          double wi = w[j];
          concordancePairWeight += (2.0 * wi) * (double) bitTotal;
          concordanceWeight     += (2.0 * wi) * (double) less + (wi) * (double) eq;
        }
      }
    }
    for (uint e = 1; e <= d; e++) {
      uint j = eventPos[e];
      rfsrc_bitAdd(bit, m, rank[j], 1);
      bitTotal++;
    }
    pos = start - 1;
  }
  free_uivector(eventPos, 1, n);
  free_uivector(touched, 1, n);
  free_dvector(rankWeightSum, 1, m);
  free_uivector(rankCount, 1, m);
  free_uivector(bit, 1, m);
  free_uivector(tIndxx, 1, n);
  free_uivector(rank, 1, n);
  free_uivector(pIndxx, 1, n);
  free_dvector(w, 1, n);
  free_dvector(p, 1, n);
  free_dvector(s, 1, n);
  free_dvector(t, 1, n);
  if (concordancePairWeight <= 0.0) {
    result = RF_nativeNaN;
  }
  else {
    result = 1.0 - (concordanceWeight / concordancePairWeight);
  }
  return result;
}
double getCRConcordanceIndexIPCW_Fenwick(uint    size,
                                         double *timePtr,
                                         double *statusPtr,
                                         double *predictedOutcome,
                                         double *denom,
                                         double *weight,
                                         uint    eventType) {
  uint i;
  double result;
  uint n = 0;
  for (i = 1; i <= size; i++) {
    if ((denom[i] != 0.0) && (weight[i] != 0.0) &&
        !RF_nativeIsNaN(timePtr[i]) && !RF_nativeIsNaN(statusPtr[i]) && !RF_nativeIsNaN(predictedOutcome[i])) {
      n++;
    }
  }
  if (n < 2) {
    return RF_nativeNaN;
  }
  double  *t   = dvector(1, n);
  double  *p   = dvector(1, n);
  double *w2   = dvector(1, n);  
  double *w1   = dvector(1, n);  
  uint   *st   = uivector(1, n);
  uint *pIndxx = uivector(1, n);
  uint *rank   = uivector(1, n);
  uint kk = 0;
  for (i = 1; i <= size; i++) {
    if ((denom[i] != 0.0) && (weight[i] != 0.0) &&
        !RF_nativeIsNaN(timePtr[i]) && !RF_nativeIsNaN(statusPtr[i]) && !RF_nativeIsNaN(predictedOutcome[i])) {
      kk++;
      t[kk]  = timePtr[i];
      st[kk] = (uint) statusPtr[i];
      p[kk]  = predictedOutcome[i];
      w2[kk] = weight[i];
      w1[kk] = sqrt(weight[i]);
    }
  }
  indexx(n, p, pIndxx);
  uint m = 1;
  double last = p[pIndxx[1]];
  rank[pIndxx[1]] = m;
  for (i = 2; i <= n; i++) {
    double cur = p[pIndxx[i]];
    if ((cur - last) > EPSILON) {
      m++;
      last = cur;
    }
    rank[pIndxx[i]] = m;
  }
  uint *tIndxx = uivector(1, n);
  indexx(n, t, tIndxx);
  double denomW = 0.0;
  double numerW = 0.0;
  uint *bitCount = uivector(1, m);
  for (i = 1; i <= m; i++) {
    bitCount[i] = 0;
  }
  uint bitTotal = 0;
  uint   *rankCount     = uivector(1, m);
  double *rankWeightSum = dvector(1, m);
  for (i = 1; i <= m; i++) {
    rankCount[i] = 0;
    rankWeightSum[i] = 0.0;
  }
  uint *touched  = uivector(1, n);
  uint *casePos  = uivector(1, n);
  int pos = (int) n;
  while (pos >= 1) {
    double curTime = t[tIndxx[pos]];
    int start = pos;
    while (start > 1 && (fabs(t[tIndxx[start - 1]] - curTime) <= EPSILON)) {
      start--;
    }
    for (int q = start; q <= pos; q++) {
      uint j = tIndxx[q];
      if (st[j] == 0) {
        rfsrc_bitAdd(bitCount, m, rank[j], 1);
        bitTotal++;
      }
    }
    uint d = 0;
    for (int q = start; q <= pos; q++) {
      uint j = tIndxx[q];
      if (st[j] == eventType) {
        casePos[++d] = j;
      }
    }
    if (d > 0) {
      double sumW = 0.0;
      uint touchedSize = 0;
      for (uint e = 1; e <= d; e++) {
        uint j = casePos[e];
        sumW += w2[j];
        uint r = rank[j];
        if (rankCount[r] == 0) {
          touched[++touchedSize] = r;
        }
        rankCount[r]++;
        rankWeightSum[r] += w2[j];
      }
      if (d >= 2) {
        double denomTie = (double) (d - 1) * sumW;
        double tieMass  = 0.0;
        for (uint u = 1; u <= touchedSize; u++) {
          uint r = touched[u];
          uint c = rankCount[r];
          if (c >= 2) {
            tieMass += (double) (c - 1) * rankWeightSum[r];
          }
          rankCount[r] = 0;
          rankWeightSum[r] = 0.0;
        }
        denomW += denomTie;
        numerW += 0.5 * denomTie + 0.5 * tieMass;
      }
      else {
        for (uint u = 1; u <= touchedSize; u++) {
          uint r = touched[u];
          rankCount[r] = 0;
          rankWeightSum[r] = 0.0;
        }
      }
    }
    if (bitTotal > 0) {
      for (uint e = 1; e <= d; e++) {
        uint j = casePos[e];
        uint r = rank[j];
        uint less = (r > 1) ? rfsrc_bitSum(bitCount, r - 1) : 0;
        uint leq  = rfsrc_bitSum(bitCount, r);
        uint eq   = leq - less;
        double wi2 = w2[j];
        denomW += (2.0 * wi2) * (double) bitTotal;
        numerW += (2.0 * wi2) * (double) less + (wi2) * (double) eq;
      }
    }
    for (int q = start; q <= pos; q++) {
      uint j = tIndxx[q];
      if (st[j] > 0) {
        rfsrc_bitAdd(bitCount, m, rank[j], 1);
        bitTotal++;
      }
    }
    pos = start - 1;
  }
  free_uivector(casePos, 1, n);
  free_uivector(touched, 1, n);
  free_dvector(rankWeightSum, 1, m);
  free_uivector(rankCount, 1, m);
  free_uivector(bitCount, 1, m);
  double *bitW = dvector(1, m);
  for (i = 1; i <= m; i++) {
    bitW[i] = 0.0;
  }
  double bitWTotal = 0.0;
  int qpos = 1;
  while (qpos <= (int) n) {
    double curTime = t[tIndxx[qpos]];
    int end = qpos;
    while (end < (int) n && (fabs(t[tIndxx[end + 1]] - curTime) <= EPSILON)) {
      end++;
    }
    for (int q = qpos; q <= end; q++) {
      uint j = tIndxx[q];
      if ((st[j] > 0) && (st[j] != eventType)) {
        double wj1 = w1[j];
        if (wj1 > 0.0) {
          rfsrc_bitAddD(bitW, m, rank[j], wj1);
          bitWTotal += wj1;
        }
      }
    }
    for (int q = qpos; q <= end; q++) {
      uint j = tIndxx[q];
      if (st[j] == eventType) {
        if (bitWTotal > 0.0) {
          uint r = rank[j];
          double lessW = (r > 1) ? rfsrc_bitSumD(bitW, r - 1) : 0.0;
          double leqW  = rfsrc_bitSumD(bitW, r);
          double eqW   = leqW - lessW;
          double wi1 = w1[j];
          denomW += (2.0 * wi1) * bitWTotal;
          numerW += (2.0 * wi1) * lessW + (wi1) * eqW;
        }
      }
    }
    qpos = end + 1;
  }
  free_dvector(bitW, 1, m);
  free_uivector(tIndxx, 1, n);
  free_uivector(rank, 1, n);
  free_uivector(pIndxx, 1, n);
  free_uivector(st, 1, n);
  free_dvector(w1, 1, n);
  free_dvector(w2, 1, n);
  free_dvector(p, 1, n);
  free_dvector(t, 1, n);
  if (denomW <= 0.0) {
    result = RF_nativeNaN;
  }
  else {
    result = 1.0 - (numerW / denomW);
  }
  return result;
}
void getCRPerformance (char     mode,
                       uint     obsSize,
                       double **responsePtr,
                       double **yearsLost,
                       double  *denom,
                       double  *weight,
                       double  *performanceVector) {
  uint   mRecordSize;
  int  **mpSign;
  uint  *mRecordIndex;
  uint  *meIndividualSize;
  uint **eIndividual;
  double *subsettedWeight;
  double concordanceIndex;
  uint j;
  if (!(RF_opt & OPT_COMP_RISK)) {
    RF_nativePrint("\nRF-SRC:  *** ERROR *** ");
    RF_nativePrint("\nRF-SRC:  Attempt at conditional performance updates in a non-CR analysis.");
    RF_nativePrint("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  if (weight != NULL) {
    for (j = 1; j <= RF_eventTypeSize; j++) {
      concordanceIndex = getCRConcordanceIndexIPCW_Fenwick(obsSize,
                                                           responsePtr[RF_timeIndex],
                                                           responsePtr[RF_statusIndex],
                                                           yearsLost[j],
                                                           denom,
                                                           weight,
                                                           RF_eventType[j]);
      if (RF_nativeIsNaN(concordanceIndex)) {
        performanceVector[j] = RF_nativeNaN;
      }
      else {
        performanceVector[j] = concordanceIndex;
      }
    }
  }
  else {
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
    if (weight != NULL) {
      subsettedWeight = dvector(1, obsSize);
    }
    else {
      subsettedWeight = NULL;
    }
    for (j = 1; j <= RF_eventTypeSize; j++) {
      getConditionalConcordanceArrays(j,
                                      responsePtr[RF_timeIndex],
                                      responsePtr[RF_statusIndex],
                                      yearsLost[j],
                                      denom,
                                      weight,
                                      meIndividualSize,
                                      eIndividual,
                                      subsettedTime,
                                      subsettedStatus,
                                      subsettedMortality,
                                      subsettedEnsembleDen,
                                      subsettedWeight);
      concordanceIndex = getConcordanceIndex(-1,
                                             meIndividualSize[j],
                                             subsettedTime,
                                             subsettedStatus,
                                             subsettedMortality,
                                             subsettedEnsembleDen,
                                             subsettedWeight);
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
    if (weight != NULL) {
      free_dvector(subsettedWeight, 1, obsSize);
    }
  }
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
