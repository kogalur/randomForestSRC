
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "splitUtil.h"
#include "factor.h"
#include "factorOps.h"
#include "nrutil.h"
#include "error.h"
char getPreSplitResultGeneric (uint      treeID,
                               Node     *parent,
                               char      multImpFlag,
                               char      multVarFlag) {
  uint i, r;
  char mResponseFlag;
  char result;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetPreSplitResultGeneric(%10d) ENTRY ...\n", treeID);
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\n  called with   rep size:  %10d", parent -> repMembrSize);
  ${trace.token}  }
  result = TRUE;
  if (result) {
    if (parent -> repMembrSize >= (2 * RF_nodeSize)) {
      result = TRUE;
    }
    else {
      result = FALSE;
      ${trace.token}    if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
      ${trace.token}      RF_nativePrint("\nLess than twice the minimum number of replicates encountered.  ");
      ${trace.token}      RF_nativePrint("\nNode will not be split.  \n");
      ${trace.token}    }
    }
  }
  if (result) {
    if (RF_nodeDepth < 0) {
      result = TRUE;
    }
    else {
      if (parent -> depth < (uint) RF_nodeDepth) {
        result = TRUE;
      }
      else {
        result = FALSE;
        ${trace.token}    if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
        ${trace.token}      RF_nativePrint("\nMaximum node depth encountered.  ");
        ${trace.token}      RF_nativePrint("\nNode will not be split.  \n");
        ${trace.token}    }
      }
    }
  }
  if (result) {
    if ((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP) || multVarFlag) {
      parent -> nonMissMembrSizeStatic = parent -> repMembrSize;
      parent -> nonMissMembrIndxStatic = RF_identityMembershipIndex;
      ${trace.token}    if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
      ${trace.token}      RF_nativePrint("\nNon-miss membership index assuming identity mapping.  ");
      ${trace.token}    }
    }
    else {
      parent -> nonMissMembrIndxStatic = uivector(1, parent -> repMembrSize);
      parent -> nonMissMembrSizeStatic = 0;
      for (i = 1; i <= parent -> repMembrSize; i++) {
        mResponseFlag = FALSE;
        if (RF_mRecordMap[parent -> repMembrIndx[i]] > 0) {
          for (r = 1; r <= RF_ySize; r++) {
            if (RF_mpSign[r][RF_mRecordMap[parent -> repMembrIndx[i]]] == 1) {
              mResponseFlag = TRUE;
            }
          }
        }
        if (!mResponseFlag) {
          (parent -> nonMissMembrSizeStatic) ++;
          (parent -> nonMissMembrIndxStatic)[parent -> nonMissMembrSizeStatic] = i;
        }
      }  
    }  
    if (!multVarFlag) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        uint q,k,m;
        uint *evntProp = uivector(1, RF_eventTypeSize + 1);
        for (q = 1; q <= RF_eventTypeSize + 1; q++) {
          evntProp[q] = 0;
        }
        for (i = 1; i <= parent -> nonMissMembrSizeStatic; i++) {
          m = (uint) RF_status[treeID][(parent -> repMembrIndx)[parent -> nonMissMembrIndxStatic[i]]];
          if (m > 0) {
            evntProp[RF_eventTypeIndex[m]] ++;
          }
          else {
            evntProp[RF_eventTypeSize + 1] ++;
          }
        }
        k = 0;
        for (q = 1; q <= RF_eventTypeSize + 1; q++) {
          if(evntProp[q] > 0) {
            k ++;
          }
        }
        if (k == 0) {
          result = FALSE;
        }
        else {
          if (k == 1) {
            if (evntProp[RF_eventTypeSize + 1] > 0) {
              result = FALSE;
            }
            else {
              result = getVariance(parent -> repMembrSize,
                                   parent -> repMembrIndx,
                                   parent -> nonMissMembrSizeStatic,
                                   parent -> nonMissMembrIndxStatic,
                                   RF_time[treeID],
                                   & (parent -> mean),
                                   NULL);
            }
          }
        }
        ${trace.token}    if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
        ${trace.token}      if (!result) {
        ${trace.token}        RF_nativePrint("\nNode purity encountered.  ");
        ${trace.token}        RF_nativePrint("\nNode will not be split.  \n");
        ${trace.token}      }
        ${trace.token}    }
        free_uivector(evntProp, 1, RF_eventTypeSize + 1);
      }
      else {
        result = getVariance(parent -> repMembrSize,
                             parent -> repMembrIndx,
                             parent -> nonMissMembrSizeStatic,
                             parent -> nonMissMembrIndxStatic,
                             RF_response[treeID][1],
                             & (parent -> mean),
                             NULL);
      }
    }
    if (!result) {
      parent -> nonMissMembrSizeStatic = 0;
      if (!((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP) || multVarFlag)) {
        free_uivector(parent -> nonMissMembrIndxStatic, 1, parent -> repMembrSize);
        parent -> nonMissMembrIndxStatic = NULL;
        parent -> nonMissMembrIndx       = NULL;
      }
    }
  }
  parent -> nonMissMembrIndx = parent -> nonMissMembrIndxStatic;
  parent -> nonMissMembrSize  = parent -> nonMissMembrSizeStatic;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetPreSplitResultGeneric(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
  return result;
}
char getPreSplitResultNoMiss (uint      treeID,
                              Node     *parent,
                              char      multImpFlag,
                              char      multVarFlag) {
  uint i;
  char result;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetPreSplitResultNoMiss(%10d) ENTRY ...\n", treeID);
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\n  called with   rep size:  %10d", parent -> repMembrSize);
  ${trace.token}  }
  result = TRUE;
  if (result) {
    if (parent -> repMembrSize >= (2 * RF_nodeSize)) {
      result = TRUE;
    }
    else {
      result = FALSE;
      ${trace.token}    if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
      ${trace.token}      RF_nativePrint("\nLess than twice the minimum number of replicates encountered.  ");
      ${trace.token}      RF_nativePrint("\nNode will not be split.  \n");
      ${trace.token}    }
    }
  }
  if (result) {
    if (RF_nodeDepth < 0) {
      result = TRUE;
    }
    else {
      if (parent -> depth < (uint) RF_nodeDepth) {
        result = TRUE;
      }
      else {
        result = FALSE;
        ${trace.token}    if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
        ${trace.token}      RF_nativePrint("\nMaximum node depth encountered.  ");
        ${trace.token}      RF_nativePrint("\nNode will not be split.  \n");
        ${trace.token}    }
      }
    }
  }
  if (result) {
    parent -> nonMissMembrSizeStatic = parent -> repMembrSize;
    parent -> nonMissMembrIndxStatic = RF_identityMembershipIndex;
    if (!multVarFlag) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        uint q,k,m;
        uint *evntProp = uivector(1, RF_eventTypeSize + 1);
        for (q = 1; q <= RF_eventTypeSize + 1; q++) {
          evntProp[q] = 0;
        }
        for (i = 1; i <= parent -> nonMissMembrSizeStatic; i++) {
          m = (uint) RF_status[treeID][(parent -> repMembrIndx)[parent -> nonMissMembrIndxStatic[i]]];
          if (m > 0) {
            evntProp[RF_eventTypeIndex[m]] ++;
          }
          else {
            evntProp[RF_eventTypeSize + 1] ++;
          }
        }
        k = 0;
        for (q = 1; q <= RF_eventTypeSize + 1; q++) {
          if(evntProp[q] > 0) {
            k ++;
          }
        }
        if (k == 0) {
          result = FALSE;
        }
        else {
          if (k == 1) {
            if (evntProp[RF_eventTypeSize + 1] > 0) {
              result = FALSE;
            }
            else {
              result = getVariance(parent -> repMembrSize,
                                   parent -> repMembrIndx,
                                   parent -> nonMissMembrSizeStatic,
                                   parent -> nonMissMembrIndxStatic,
                                   RF_time[treeID],
                                   & (parent -> mean),
                                   NULL);
            }
          }
        }
        ${trace.token}    if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
        ${trace.token}      if (!result) {
        ${trace.token}        RF_nativePrint("\nNode purity encountered.  ");
        ${trace.token}        RF_nativePrint("\nNode will not be split.  \n");
        ${trace.token}      }
        ${trace.token}    }
        free_uivector(evntProp, 1, RF_eventTypeSize + 1);
      }
      else {
        result = getVariance(parent -> repMembrSize,
                             parent -> repMembrIndx,
                             parent -> nonMissMembrSizeStatic,
                             parent -> nonMissMembrIndxStatic,
                             RF_response[treeID][1],
                             & (parent -> mean),
                             NULL);
      }
    }
    if (!result) {
      parent -> nonMissMembrSizeStatic = 0;
    }
  }
  parent -> nonMissMembrIndx = parent -> nonMissMembrIndxStatic;
  parent -> nonMissMembrSize  = parent -> nonMissMembrSizeStatic;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetPreSplitResultNoMiss(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
  return result;
}
void unstackPreSplit (char      preliminaryResult,
                      Node     *parent,
                      char      multImpFlag,
                      char      multVarFlag) {
  ${trace.token}  if (getTraceFlag(0) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nunstackPreSplit() ENTRY ...\n");
  ${trace.token}  }
  if (preliminaryResult) {
    if (!((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP) || multVarFlag)) {
      free_uivector(parent -> nonMissMembrIndxStatic, 1, parent -> repMembrSize);
    }
  }
  else {
  }
  ${trace.token}  if (getTraceFlag(0) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nunstackPreSplit() EXIT ...\n");
  ${trace.token}  }
}
void stackSplitPreliminary(uint     nodeSize,
                           char   **localSplitIndicator,
                           double **splitVector) {
  ${trace.token}  if (getTraceFlag(0) & SPLT_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nstackSplitPreliminary() ENTRY ...\n");
  ${trace.token}    }
  ${trace.token}  }
  *localSplitIndicator = cvector(1, nodeSize);
  *splitVector = dvector(1, nodeSize);
  ${trace.token}  if (getTraceFlag(0) & SPLT_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nstackSplitPreliminary() EXIT ...\n");
  ${trace.token}    }
  ${trace.token}  }
}
void unstackSplitPreliminary(uint    nodeSize,
                             char   *localSplitIndicator,
                             double *splitVector) {
  ${trace.token}  if (getTraceFlag(0) & SPLT_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nunstackSplitPreliminary() ENTRY ...\n");
  ${trace.token}    }
  ${trace.token}  }
  free_cvector(localSplitIndicator, 1, nodeSize);
  free_dvector(splitVector, 1, nodeSize);
  ${trace.token}  if (getTraceFlag(0) & SPLT_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nunstackSplitPreliminary() EXIT ...\n");
  ${trace.token}    }
  ${trace.token}  }
}
DistributionObj *stackRandomCovariatesGeneric(uint treeID, Node *parent) {
  uint actualWeightType;
  uint *augmentationSize;
  char *permissible;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackRandomCovariatesGeneric(%10d) ENTRY ...\n", treeID);
  ${trace.token}  }
  DistributionObj *obj = makeDistributionObjRaw();
  actualWeightType = RF_xWeightType;
  augmentationSize    = NULL;
  permissible    = parent -> permissible;
  obj -> permissibleIndex    = NULL;
  obj -> permissible         = permissible;
  obj -> permissibleSize     = parent -> xSize;
  obj -> augmentationSize    = augmentationSize;
  obj -> weightType          = actualWeightType;
  obj -> weight              = RF_xWeight;
  obj -> weightSorted        = RF_xWeightSorted;
  obj -> densityAllocSize    = RF_xWeightDensitySize;
  initializeCDFNew(treeID, obj);
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackRandomCovariatesGeneric(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
  return obj;
}
void unstackRandomCovariatesGeneric(uint treeID, DistributionObj *obj) {
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nunstackRandomCovariatesGeneric(%10d) ENTRY ...\n", treeID);
  ${trace.token}  }
  if (obj -> augmentationSize != NULL) {
    free_uivector(obj -> augmentationSize, 1, 2);
    obj -> augmentationSize = NULL;
  }
  discardCDFNew(treeID, obj);
  freeDistributionObjRaw(obj);
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nunstackRandomCovariatesGeneric(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
}
char selectRandomCovariatesGeneric(uint     treeID,
                                   Node     *parent,
                                   DistributionObj *distributionObj,
                                   char     *factorFlag,
                                   uint     *covariate,
                                   uint     *covariateCount) {
  char xVarFound;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nselectRandomCovariatesGeneric(%10d) ENTRY ...\n", treeID);
  ${trace.token}  }
  (*covariate) = UINT_MAX;
  xVarFound = FALSE;
  *factorFlag = FALSE;
  while ( ((*covariateCount) < RF_mtry) && ((*covariate) != 0) && (xVarFound == FALSE)) {
    (*covariateCount) ++;
    ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
    ${trace.token}        RF_nativePrint("\nCovariate counter:  %10d \n", (*covariateCount));
    ${trace.token}      }
    *covariate = sampleFromCDFNew(ran1B, treeID, distributionObj);
    if (*covariate != 0) {
      updateCDFNew(treeID, distributionObj);
      ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
      ${trace.token}        RF_nativePrint("\nCovariate:  %10d \n", *covariate);
      ${trace.token}      }
      xVarFound = TRUE;
        if (RF_xType[*covariate] == 'C') {
          *factorFlag = TRUE;
          ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
          ${trace.token}        RF_nativePrint("\nCandidate covariate %10d is a factor.", *covariate);
          ${trace.token}      }
        }
    }  
  }  
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nselectRandomCovariatesGeneric(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
  return xVarFound;
}
uint stackAndConstructSplitVectorGenericPhase1 (uint     treeID,
                                                Node    *parent,
                                                uint     covariate,
                                                ...) {
  uint offset;
  double *nonMissSplit;
  char mPredictorFlag;
  uint vectorSize;
  uint i, ii;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAndConstructSplitVectorGenericPhase1(%10d) ENTRY ...\n", treeID);
  ${trace.token}  }
  uint  *repMembrIndx = parent -> repMembrIndx;
  uint   repMembrSize = parent -> repMembrSize;
  uint  *nonMissMembrIndxStatic = parent -> nonMissMembrIndxStatic;
  uint   nonMissMembrSizeStatic = parent -> nonMissMembrSizeStatic;
  uint **nonMissMembrIndx = & (parent -> nonMissMembrIndx);
  uint  *nonMissMembrSize = & (parent -> nonMissMembrSize);
  va_list list;
  va_start(list, covariate);
  double *splitVector = va_arg(list, double*);
  uint **indxx = (uint**) va_arg(list, uint**);
  char   multImpFlag = (char) va_arg(list, uint);
  nonMissSplit = dvector(1, repMembrSize);
  if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
    *nonMissMembrSize = nonMissMembrSizeStatic;
    *nonMissMembrIndx = nonMissMembrIndxStatic;
  }
  else {
    *nonMissMembrSize = 0;
    *nonMissMembrIndx = uivector(1, nonMissMembrSizeStatic);
  }
  vectorSize = 0;
  if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      for (i = 1; i <= repMembrSize; i++) {
        nonMissSplit[i] = RF_observation[treeID][covariate][repMembrIndx[i]];
      }
      ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
      ${trace.token}          RF_nativePrint("\nNote that all responses are inculded below, and not just the ytry targets:");
      ${trace.token}        RF_nativePrint("\nReplicate:           indx    targ indx        covar       responses ->");
      ${trace.token}        for (i = 1; i <= repMembrSize; i++) {
      ${trace.token}          RF_nativePrint("\n               %10d %12d %12.4f",
      ${trace.token}                    i,
      ${trace.token}                    repMembrIndx[i],
      ${trace.token}                    RF_observation[treeID][covariate][repMembrIndx[i]]);
      ${trace.token}          for (uint r = 1; r<= RF_ySize; r++) {
      ${trace.token}            RF_nativePrint(" %12.4f", RF_response[treeID][r][repMembrIndx[i]]);
      ${trace.token}          }
      ${trace.token}        }
      ${trace.token}      }
  }
  else {
    offset = RF_ySize + covariate;
    (*nonMissMembrSize) = 0;
    for (i = 1; i <= nonMissMembrSizeStatic; i++) {
      ii = nonMissMembrIndxStatic[i];
      mPredictorFlag = FALSE;
      if (RF_mRecordMap[repMembrIndx[ii]] > 0) {
        if (RF_mpSign[offset][RF_mRecordMap[repMembrIndx[ii]]] == 1) {
          mPredictorFlag = TRUE;
        }
      }
      if (!mPredictorFlag) {
        (*nonMissMembrSize) ++;
        (*nonMissMembrIndx)[*nonMissMembrSize] = ii;
        nonMissSplit[*nonMissMembrSize] = RF_observation[treeID][covariate][repMembrIndx[(*nonMissMembrIndx)[*nonMissMembrSize]]];
      }
    }  
    ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
    ${trace.token}          RF_nativePrint("\nNote that all responses are inculded below, and not just the ytry targets:");
    ${trace.token}          RF_nativePrint("\nReplicate:           indx    targ indx        covar         responses ->");
    ${trace.token}          for (i = 1; i <= (*nonMissMembrSize); i++) {
    ${trace.token}            RF_nativePrint("\n               %10d %12d %12.4f",
    ${trace.token}                    i,
    ${trace.token}                    repMembrIndx[(*nonMissMembrIndx)[i]],
    ${trace.token}                    RF_observation[treeID][covariate][repMembrIndx[(*nonMissMembrIndx)[i]]]);
    ${trace.token}            for (uint r = 1; r<= RF_ySize; r++) {
    ${trace.token}              RF_nativePrint(" %12.4f", RF_response[treeID][r][repMembrIndx[(*nonMissMembrIndx)[i]]]);
    ${trace.token}            }
    ${trace.token}          }
    ${trace.token}        }
  }  
  if ((*nonMissMembrSize) >= 2) {
    (*indxx) = uivector(1, repMembrSize);
    indexx((*nonMissMembrSize),
           nonMissSplit,
           (*indxx));
    splitVector[1] = nonMissSplit[(*indxx)[1]];
    vectorSize = 1;
    for (i = 2; i <= (*nonMissMembrSize); i++) {
      if (nonMissSplit[(*indxx)[i]] > splitVector[vectorSize]) {
        vectorSize ++;
        splitVector[vectorSize] = nonMissSplit[(*indxx)[i]];
      }
    }
    ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
    ${trace.token}        RF_nativePrint("\n\nNon-miss member size:        %10d ", (*nonMissMembrSize));
    ${trace.token}        RF_nativePrint("\nPermissible split size:      %10d ", vectorSize);
    ${trace.token}      }
    if(vectorSize >= 2) {
      ${trace.token}        if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
      ${trace.token}          RF_nativePrint("\nCovariate accepted:          %10d \n", covariate);
      ${trace.token}          RF_nativePrint("\nsplit points:       index        value");
      ${trace.token}          for (i = 1; i <= vectorSize; i++) {
      ${trace.token}            RF_nativePrint("\n               %10d %12.4f", i, splitVector[i]);
      ${trace.token}          }
      ${trace.token}          RF_nativePrint("\n");
      ${trace.token}          RF_nativePrint("\nNon-miss Repl:       indx        indxx  unord value    ord value");
      ${trace.token}          for (i = 1; i <= (*nonMissMembrSize); i++) {
      ${trace.token}            RF_nativePrint("\n               %10d %12d %12.4f %12.4f ",
      ${trace.token}                    i,
      ${trace.token}                    (*indxx)[i],
      ${trace.token}                    nonMissSplit[i],
      ${trace.token}                    nonMissSplit[(*indxx)[i]]
      ${trace.token}                   );
      ${trace.token}          }
      ${trace.token}        }
    }
    else {
      vectorSize = 0;
      if (covariate <= RF_xSize) {
        (parent -> permissible)[covariate] = FALSE;
        parent -> permissibleReIndxFlag = TRUE;
      }
      free_uivector(*indxx, 1, repMembrSize);
      if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
        *nonMissMembrSize = 0;
        *nonMissMembrIndx = NULL;
      }
      else {
        free_uivector(*nonMissMembrIndx, 1, nonMissMembrSizeStatic);
        *nonMissMembrSize = 0;
        *nonMissMembrIndx = NULL;
      }
    ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
    ${trace.token}        RF_nativePrint("\nCovariate rejected due to non-miss split being less than two (2):  %10d ", covariate);
    ${trace.token}      }
    }
  }
  else {
    vectorSize = 0;
    if (covariate <= RF_xSize) {
      (parent -> permissible)[covariate] = FALSE;
      parent -> permissibleReIndxFlag = TRUE;
    }
    ${trace.token}        if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
    ${trace.token}          RF_nativePrint("\nCovariate rejected:  %10d \n", covariate);
    ${trace.token}        }
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      *nonMissMembrSize = 0;
      *nonMissMembrIndx = NULL;
    }
    else {
      free_uivector(*nonMissMembrIndx, 1, nonMissMembrSizeStatic);
      *nonMissMembrSize = 0;
      *nonMissMembrIndx = NULL;
    }
  }
  free_dvector(nonMissSplit, 1, repMembrSize);
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAndConstructSplitVectorGenericPhase1(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
  return vectorSize;
}
uint stackAndConstructSplitVectorGenericPhase2 (uint     treeID,
                                                Node    *parent,
                                                uint     covariate,
                                                double  *splitVector,
                                                uint     vectorSize,
                                                char    *factorFlag,
                                                char    *deterministicSplitFlag,
                                                uint    *mwcpSizeAbsolute,
                                                void   **splitVectorPtr) {
  uint repMembrSize;
  uint  sworIndex;
  uint *sworVector;
  uint  sworVectorSize;
  uint j, j2, k2;
  uint factorSizeAbsolute;
  uint offset;
  uint splitLength;
  uint relativePair;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAndConstructSplitVectorGenericPhase2(%10d) ENTRY ...\n", treeID);
  ${trace.token}  }
  repMembrSize = parent -> repMembrSize;
  splitLength = 0;  
  (*splitVectorPtr) = NULL;  
  if (vectorSize < 2) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Split vector is of insufficient size in Stack Phase II allocation:  %10d", vectorSize);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\n(TreeID %10d):  Constructing split vector for (parameter, of size):  (%10d, %10d) \n", treeID, covariate, vectorSize);
  ${trace.token}  }
  if (*factorFlag) {
    if(RF_factorList[treeID][vectorSize] == NULL) {
      RF_factorList[treeID][vectorSize] = makeFactor(vectorSize, FALSE);
    }
    factorSizeAbsolute = RF_xFactorSize[RF_xFactorMap[covariate]];
    *mwcpSizeAbsolute = RF_factorList[treeID][factorSizeAbsolute] -> mwcpSize;
    ${trace.token}    if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
    ${trace.token}    RF_nativePrint("\n(%10d):  (Absolute Factor Size, Absolute MWCP Size):  (%10d, %10d)", treeID, factorSizeAbsolute, *mwcpSizeAbsolute);
    ${trace.token}    }
    if (RF_splitRule == RAND_SPLIT) {
      splitLength = 1 + 1;
      *deterministicSplitFlag = FALSE;
    }
    else {
      if (RF_nsplit == 0) {
        *deterministicSplitFlag = TRUE;
        if ((RF_factorList[treeID][vectorSize] -> r) > MAX_EXACT_LEVEL) {
          *deterministicSplitFlag = FALSE;
        }
        else {
          if ( *((uint *) RF_factorList[treeID][vectorSize] -> complementaryPairCount) >= repMembrSize ) {
            *deterministicSplitFlag = FALSE;
          }
        }
        if (*deterministicSplitFlag == FALSE) {
          splitLength = repMembrSize + 1;
          ${trace.token}          if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
          ${trace.token}            if (vectorSize <= MAX_EXACT_LEVEL) {
          ${trace.token}              RF_nativePrint("\n(%10d):  Factor override to random (pSplit, nsplit, ndSize):  (%10d, %10d, %10d) \n",
          ${trace.token}                      treeID,
          ${trace.token}                      *((uint*) RF_factorList[treeID][vectorSize] -> complementaryPairCount),
          ${trace.token}                      RF_nsplit,
          ${trace.token}                      repMembrSize);
          ${trace.token}            }
          ${trace.token}            else {
          ${trace.token}              RF_nativePrint("\n(%10d):  Factor override to random (pSplit, nsplit, ndSize):  (%24.0f, %10d, %10d) \n",
          ${trace.token}                      treeID,
          ${trace.token}                      *((double*) RF_factorList[treeID][vectorSize] -> complementaryPairCount),
          ${trace.token}                      RF_nsplit,
          ${trace.token}                      repMembrSize);
          ${trace.token}            }
          ${trace.token}          }
        }
        else {
          splitLength = *((uint*) RF_factorList[treeID][vectorSize] -> complementaryPairCount) + 1;
        }
      }
      else {
        *deterministicSplitFlag = FALSE;
        if ((RF_factorList[treeID][vectorSize] -> r) <= MAX_EXACT_LEVEL) {
          if (*((uint*) RF_factorList[treeID][vectorSize] -> complementaryPairCount) <= ((RF_nsplit <= repMembrSize) ? RF_nsplit : repMembrSize)) {
            splitLength = *((uint*) RF_factorList[treeID][vectorSize] -> complementaryPairCount) + 1;
            *deterministicSplitFlag = TRUE;
            ${trace.token}            if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
            ${trace.token}              RF_nativePrint("\nFactor override to determ (pSplit, nsplit, ndSize):  (%10d, %10d, %10d) \n",
            ${trace.token}                      *((uint*) RF_factorList[treeID][vectorSize] -> complementaryPairCount),
            ${trace.token}                      RF_nsplit,
            ${trace.token}                      repMembrSize);
            ${trace.token}            }
          }
        }
        if (*deterministicSplitFlag == FALSE) {
          splitLength = 1 + ((RF_nsplit <= repMembrSize) ? RF_nsplit : repMembrSize);
        }
      }  
    }  
    (*splitVectorPtr) = uivector(1, splitLength * (*mwcpSizeAbsolute));
    for (offset = 1; offset <= *mwcpSizeAbsolute; offset++) {
      ((uint*) (*splitVectorPtr) + ((splitLength - 1) * (*mwcpSizeAbsolute)))[offset] = 0;
    }
    if (*deterministicSplitFlag) {
      bookFactor(RF_factorList[treeID][vectorSize]);
      j2 = 0;
      for (j = 1; j <= RF_factorList[treeID][vectorSize] -> cardinalGroupCount; j++) {
        for (k2 = 1; k2 <= ((uint*) RF_factorList[treeID][vectorSize] -> cardinalGroupSize)[j]; k2++) {
          ++j2;
          relativePair = (RF_factorList[treeID][vectorSize] -> cardinalGroupBinary)[j][k2];
          convertRelToAbsBinaryPair(treeID,
                                    vectorSize,
                                    factorSizeAbsolute,
                                    relativePair,
                                    splitVector,
                                    (uint*) (*splitVectorPtr) + ((j2 - 1) * (*mwcpSizeAbsolute)));
        }
      }
    }  
    else {
      for (j = 1; j < splitLength; j++) {
        getRandomPair(treeID, vectorSize, factorSizeAbsolute, splitVector, (uint*) (*splitVectorPtr) + ((j - 1) * (*mwcpSizeAbsolute)));
      }
    }
  }  
  else {
    if (RF_splitRule == RAND_SPLIT) {
      splitLength = 1 + 1;
      *deterministicSplitFlag = FALSE;
    }
    else {
      if(RF_nsplit == 0) {
        splitLength = vectorSize;
        (*splitVectorPtr) = splitVector;
        *deterministicSplitFlag = TRUE;
      }
      else {
        if (vectorSize <= RF_nsplit + 1) {
          splitLength = vectorSize;
          (*splitVectorPtr) = splitVector;
          *deterministicSplitFlag = TRUE;
          ${trace.token}          if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
          ${trace.token}            RF_nativePrint("\nContinuous override to determ (pSplit, nsplit, ndSize):  (%10d, %10d, %10d) \n",
          ${trace.token}                    vectorSize,
          ${trace.token}                    RF_nsplit,
          ${trace.token}                    repMembrSize);
          ${trace.token}          }
        }
        else {
          splitLength = RF_nsplit + 1;
          *deterministicSplitFlag = FALSE;
        }
      }  
    }  
    if (*deterministicSplitFlag == FALSE) {
      (*splitVectorPtr) = dvector(1, splitLength);
      ((double*) (*splitVectorPtr))[splitLength] = 0;
      if (RF_splitRule == RAND_SPLIT) {
        ((double*) (*splitVectorPtr))[1]  = splitVector[(uint) ceil(ran1B(treeID) * ((vectorSize - 1) * 1.0))];
      }
      else {
        ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
        ${trace.token}    RF_nativePrint("\nConstructing SWOR of (populated) length:  %10d - 1 \n", splitLength);
        ${trace.token}  }
        sworVector = uivector(1, vectorSize);
        sworVectorSize = vectorSize - 1;
        for (j = 1; j <= sworVectorSize; j++) {
          sworVector[j] = j;
        }
        for (j = 1; j < splitLength; j++) {
          sworIndex = (uint) ceil(ran1B(treeID) * (sworVectorSize * 1.0));
          ((double*) (*splitVectorPtr))[j]  = splitVector[sworVector[sworIndex]];
          sworVector[sworIndex] = sworVector[sworVectorSize];
          sworVectorSize --;
        }
        free_uivector (sworVector, 1, vectorSize);
        qksort(((double*) (*splitVectorPtr)), splitLength-1);
      }
    }
  }  
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAndConstructSplitVectorGenericPhase2(%10d) length:  %10d", treeID, splitLength);
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAndConstructSplitVectorGenericPhase2(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
  return splitLength;
}
void unstackSplitVectorGeneric(uint   treeID,
                               Node  *parent,
                               uint   splitLength,
                               char   factorFlag,
                               uint   splitVectorSize,
                               uint   mwcpSizeAbsolute,
                               char   deterministicSplitFlag,
                               void  *splitVectorPtr,
                               char   multImpFlag,
                               uint  *indxx) {
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nunstackSplitVectorGeneric(%10d) ENTRY ...\n", treeID);
  ${trace.token}  }
  if (splitLength > 0) {
    if (factorFlag == TRUE) {
      free_uivector(splitVectorPtr, 1, splitLength * mwcpSizeAbsolute);
      if (deterministicSplitFlag == FALSE) {
        if (splitVectorSize > SAFE_FACTOR_SIZE) {
          unbookFactor(RF_factorList[treeID][splitVectorSize]);
        }
      }
    }
    else {
      if (deterministicSplitFlag == FALSE) {
        free_dvector(splitVectorPtr, 1, splitLength);
      }
    }
    if (indxx != NULL) {
      free_uivector((indxx), 1, parent -> repMembrSize);
    }
    if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
      free_uivector(parent -> nonMissMembrIndx, 1, parent -> nonMissMembrSizeStatic);
    }
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nunstackSplitVectorGeneric(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
}
uint virtuallySplitNodeGeneric(uint  treeID,
                               Node *parent,
                               char  factorFlag,
                               uint  mwcpSizeAbsolute,
                               double *observation,
                               uint *indxx,
                               void *splitVectorPtr,
                               uint  offset,
                               char *localSplitIndicator,
                               uint *leftSize,
                               uint  priorMembrIter,
                               uint *currentMembrIter) {
  uint *repMembrIndx, *nonMissMembrIndx;
  uint nonMissMembrSize;
  char daughterFlag;
  char iterFlag;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}  RF_nativePrint("\nvirtuallySplitNodeGeneric(%10d) ENTRY ...\n", treeID);
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    if (factorFlag == TRUE) {
  ${trace.token}      RF_nativePrint("\nVirtually splitting on vector element with factor (index, binary words (hex)):  ");
  ${trace.token}      RF_nativePrint("( %10d, ", offset);
  ${trace.token}      for (uint k = mwcpSizeAbsolute; k >= 1; k--) {
  ${trace.token}        RF_nativePrint("%8x", ((uint*) splitVectorPtr + ((offset - 1) * mwcpSizeAbsolute))[k]);
  ${trace.token}      }
  ${trace.token}      RF_nativePrint(")");
  ${trace.token}    }
  ${trace.token}    else {
  ${trace.token}      RF_nativePrint("\nVirtually splitting on vector element with real (index, value):  ( %10d, %12.4f) \n", offset, ((double*) splitVectorPtr)[offset]);
  ${trace.token}    }
  ${trace.token}  }
  iterFlag = TRUE;
  repMembrIndx = parent -> repMembrIndx;
  nonMissMembrIndx = parent -> nonMissMembrIndx;
  nonMissMembrSize = parent -> nonMissMembrSize;
  *currentMembrIter = priorMembrIter;
  while (iterFlag) {
    (*currentMembrIter) ++;
    ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
    ${trace.token}        if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
    ${trace.token}          RF_nativePrint("\nCurrent   (indx):  %10d ", *currentMembrIter);
    ${trace.token}          if (factorFlag != TRUE) {
    ${trace.token}            RF_nativePrint("\nCurrent  (indxx):  %10d ", indxx[*currentMembrIter]);
    ${trace.token}            RF_nativePrint("\nCurrent (nmmIdx):  %10d ", nonMissMembrIndx[indxx[*currentMembrIter]]);
    ${trace.token}          }
    ${trace.token}          else {    
    ${trace.token}            RF_nativePrint("\nCurrent (nmmIdx):  %10d ", nonMissMembrIndx[*currentMembrIter]);
    ${trace.token}          }
    ${trace.token}        }
    ${trace.token}      }
    if (factorFlag == TRUE) {
      daughterFlag = splitOnFactor((uint)  observation[    repMembrIndx[nonMissMembrIndx[*currentMembrIter]]     ],
                                   (uint*) splitVectorPtr + ((offset - 1) * mwcpSizeAbsolute));
      localSplitIndicator[     nonMissMembrIndx[*currentMembrIter]   ] = daughterFlag;
      if ((*currentMembrIter) == nonMissMembrSize) {
        iterFlag = FALSE;
      }
      ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
      ${trace.token}        if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
      ${trace.token}          if (daughterFlag == LEFT) {
      ${trace.token}            RF_nativePrint("\nMember of LEFT Daughter (index):  %10d %10d %10d %12.4f", *currentMembrIter, nonMissMembrIndx[*currentMembrIter], repMembrIndx[nonMissMembrIndx[*currentMembrIter]], observation[ repMembrIndx[nonMissMembrIndx[*currentMembrIter]] ]);
      ${trace.token}          }
      ${trace.token}          else {
      ${trace.token}            RF_nativePrint("\nMember of RGHT Daughter (index):  %10d %10d %10d %12.4f", *currentMembrIter, nonMissMembrIndx[*currentMembrIter], repMembrIndx[nonMissMembrIndx[*currentMembrIter]], observation[ repMembrIndx[nonMissMembrIndx[*currentMembrIter]] ]);
      ${trace.token}          }
      ${trace.token}        }
      ${trace.token}      }
    }
    else {
      if ((((double*) splitVectorPtr)[offset] - observation[   repMembrIndx[nonMissMembrIndx[indxx[*currentMembrIter]]]    ]) >= 0.0) {
        daughterFlag = LEFT;
      }
      else {
        daughterFlag = RIGHT;
        iterFlag = FALSE;
      }
      localSplitIndicator[     nonMissMembrIndx[indxx[*currentMembrIter]]   ] = daughterFlag;
      ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
      ${trace.token}        if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
      ${trace.token}          if (daughterFlag == LEFT) {
      ${trace.token}            RF_nativePrint("\nMember of LEFT Daughter (index):  %10d %10d %10d %10d %12.4f", *currentMembrIter, indxx[*currentMembrIter], nonMissMembrIndx[indxx[*currentMembrIter]], repMembrIndx[nonMissMembrIndx[indxx[*currentMembrIter]]], observation[ repMembrIndx[nonMissMembrIndx[indxx[*currentMembrIter]]] ]);
      ${trace.token}          }
      ${trace.token}          else {
      ${trace.token}            RF_nativePrint("\nMember of RGHT Daughter (index):  %10d %10d %10d %10d %12.4f", *currentMembrIter, indxx[*currentMembrIter], nonMissMembrIndx[indxx[*currentMembrIter]], repMembrIndx[nonMissMembrIndx[indxx[*currentMembrIter]]], observation[ repMembrIndx[nonMissMembrIndx[indxx[*currentMembrIter]]] ]);
      ${trace.token}          }
      ${trace.token}        }
      ${trace.token}      }
    }
    if (daughterFlag == LEFT) {
      (*leftSize) ++;
    }
  }  
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nLeft Size:  %10d, Right Size:  %10d", *leftSize, nonMissMembrSize - *leftSize);
  ${trace.token}  }
  if ((*leftSize == 0) || (*leftSize == nonMissMembrSize)) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Left or Right Daughter of size zero:  (%10d, %10d)", *leftSize, nonMissMembrSize);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nvirtuallySplitNodeGeneric(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
  return (*leftSize);
}
char summarizeSplitResult(SplitInfoMax *splitInfoMax) {
  char result;
  ${trace.token}  uint k;
  ${trace.token}  if (getTraceFlag(0) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nsummarizeSplitResult() ENTRY ...\n");
  ${trace.token}  }
  if (!RF_nativeIsNaN(splitInfoMax -> deltaMax)) {
    splitInfoMax -> splitStatistic = splitInfoMax -> deltaMax;
    result = TRUE;
    ${trace.token}    if (getTraceFlag(0) & SPLT_MED_TRACE) {
    ${trace.token}      RF_nativePrint("\nBest Split Statistics: \n");
    ${trace.token}      RF_nativePrint("  SplitParm        Delta \n");
    ${trace.token}      RF_nativePrint(" %10d %12.4f \n", splitInfoMax -> splitParameterMax, splitInfoMax -> deltaMax);
    ${trace.token}      if (RF_xType[splitInfoMax -> splitParameterMax] == 'C') {
    ${trace.token}        RF_nativePrint(" at MWCPsize= %2d, mwcp= ", splitInfoMax -> splitValueMaxFactSize);
    ${trace.token}        for (k = splitInfoMax -> splitValueMaxFactSize; k >= 1; k--) {
    ${trace.token}          RF_nativePrint("%8x ", splitInfoMax -> splitValueMaxFactPtr[k]);
    ${trace.token}        }
    ${trace.token}        RF_nativePrint("\n");
    ${trace.token}      }
    ${trace.token}      else {
    ${trace.token}        RF_nativePrint(" at %12.4f \n", splitInfoMax -> splitValueMaxCont);
    ${trace.token}      }
    ${trace.token}    }
  }
  else {
    result = FALSE;
  }
  ${trace.token}  if (getTraceFlag(0) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nsummarizeSplitResult(%1d) EXIT ...\n", result);
  ${trace.token}  }
  return result;
}
char updateMaximumSplitGeneric(uint    treeID,
                               Node   *parent,
                               double  delta,
                               uint    covariate,
                               uint    index,
                               char    factorFlag,
                               uint    mwcpSizeAbsolute,
                               uint    repMembrSize,
                               char  **polarity,
                               void   *splitVectorPtr,
                               SplitInfoMax *splitInfoMax) {
  char flag;
  uint k;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdateMaximumSplitGeneric() ENTRY ...\n");
  ${trace.token}  }
  if(RF_nativeIsNaN(delta)) {
    flag = FALSE;
  }
  else {
    delta = delta * RF_xWeightStat[covariate];
    if(RF_nativeIsNaN(splitInfoMax -> deltaMax)) {
      flag = TRUE;
    }
    else {
      if ((delta - (splitInfoMax -> deltaMax)) > EPSILON) {
        flag = TRUE;
      }
      else {
        flag = FALSE;
      }
    }
  }
  if (flag) {
    splitInfoMax -> deltaMax = delta;
    splitInfoMax -> splitParameterMax = covariate;
    ${trace.token}    if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
    ${trace.token}      RF_nativePrint("\n\nRunning Split Statistics Updated: \n");
    ${trace.token}      RF_nativePrint("  SplitParm  SplitValIdx        Delta \n");
    ${trace.token}      RF_nativePrint(" %10d %12d %12.4f \n", covariate, index, splitInfoMax -> deltaMax);
    ${trace.token}    }
      ${trace.token}      if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
      ${trace.token}        RF_nativePrint(" with nominal x-var");
      ${trace.token}      }
    if (factorFlag == TRUE) {
      if (splitInfoMax -> splitValueMaxFactSize > 0) {
        if (splitInfoMax -> splitValueMaxFactSize != mwcpSizeAbsolute) {
          free_uivector(splitInfoMax -> splitValueMaxFactPtr, 1, splitInfoMax -> splitValueMaxFactSize);
          splitInfoMax -> splitValueMaxFactSize = mwcpSizeAbsolute;
          splitInfoMax -> splitValueMaxFactPtr = uivector(1, splitInfoMax -> splitValueMaxFactSize);
        }
      }
      else {
        splitInfoMax -> splitValueMaxFactSize = mwcpSizeAbsolute;
        splitInfoMax -> splitValueMaxFactPtr = uivector(1, splitInfoMax -> splitValueMaxFactSize);
      }
      splitInfoMax -> splitValueMaxCont = RF_nativeNaN;
      for (k=1; k <= splitInfoMax -> splitValueMaxFactSize; k++) {
        (splitInfoMax -> splitValueMaxFactPtr)[k] =
          ((uint*) splitVectorPtr + ((index - 1) * (splitInfoMax -> splitValueMaxFactSize)))[k];
      }
      ${trace.token}      if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
      ${trace.token}        RF_nativePrint(" at MWCPsize= %2d, mwcp= ", splitInfoMax -> splitValueMaxFactSize);
      ${trace.token}        for (k = splitInfoMax -> splitValueMaxFactSize; k >= 1; k--) {
      ${trace.token}          RF_nativePrint("%8x ", (splitInfoMax -> splitValueMaxFactPtr)[k]);
      ${trace.token}        }
      ${trace.token}        RF_nativePrint("\n");
      ${trace.token}      }
    }
    else {
      if (splitInfoMax -> splitValueMaxFactSize > 0) {
        free_uivector(splitInfoMax -> splitValueMaxFactPtr, 1, splitInfoMax -> splitValueMaxFactSize);
        splitInfoMax -> splitValueMaxFactSize = 0;
        splitInfoMax -> splitValueMaxFactPtr = NULL;
      }
      else {
      }
      splitInfoMax -> splitValueMaxCont = ((double*) splitVectorPtr)[index];
      ${trace.token}      if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
      ${trace.token}        RF_nativePrint(" at continuous point %12.4f \n", splitInfoMax -> splitValueMaxCont);
      ${trace.token}      }
    }
  }
  else {
    ${trace.token}    if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
    ${trace.token}      RF_nativePrint("\n\nRunning Split Statistics NOT Updated: \n");
    ${trace.token}      RF_nativePrint("  SplitParm  SplitValIdx        Delta     DeltaMax\n");
    ${trace.token}      RF_nativePrint(" %10d %12d %12.4f %12.4f\n", covariate, index, delta, splitInfoMax -> deltaMax);
    ${trace.token}    }
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdateMaximumSplitGeneric() EXIT ...\n");
  ${trace.token}  }
  return flag;
}
void getReweightedRandomPair (uint    treeID,
                              uint    relativeFactorSize,
                              uint    absoluteFactorSize,
                              double *absoluteLevel,
                              uint   *result) {
  uint randomGroupIndex;
  ${trace.token}  if (getTraceFlag(treeID) & FACT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetReweightedRandomPair() ENTRY ...\n");
  ${trace.token}  }
  if(RF_factorList[treeID][relativeFactorSize] == NULL) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Factor not allocated for size:  %10d", relativeFactorSize);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  randomGroupIndex = (uint) ceil(ran1B(treeID) * ((RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount) * 1.0));
  ${trace.token}  if (getTraceFlag(treeID) & FACT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nRandomly Selected Group Index:  %10d", randomGroupIndex);
  ${trace.token}  }
  createRandomBinaryPair(treeID, relativeFactorSize, absoluteFactorSize, randomGroupIndex, absoluteLevel, result);
  ${trace.token}  if (getTraceFlag(treeID) & FACT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetReweightedRandomPair() EXIT ...\n");
  ${trace.token}  }
}
void getRandomPair (uint treeID, uint relativeFactorSize, uint absoluteFactorSize, double *absoluteLevel, uint *result) {
  uint randomGroupIndex;
  double randomValue;
  uint k;
  ${trace.token}  if (getTraceFlag(treeID) & FACT_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetRandomPair() ENTRY ...\n");
  ${trace.token}  }
  if(RF_factorList[treeID][relativeFactorSize] == NULL) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Factor not allocated for size:  %10d", relativeFactorSize);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  double *cdf = dvector(1, RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount);
  if (relativeFactorSize <= MAX_EXACT_LEVEL) {
    for (k=1; k <= RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount; k++) {
      cdf[k] = (double) ((uint*) RF_factorList[treeID][relativeFactorSize] -> cardinalGroupSize)[k];
    }
  }
  else {
    for (k=1; k <= RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount; k++) {
      cdf[k] = ((double*) RF_factorList[treeID][relativeFactorSize] -> cardinalGroupSize)[k];
    }
  }
  ${trace.token}  if (getTraceFlag(treeID) & FACT_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nFactor (relativeFactorSize, cardinalGroupCount):  (%10d, %10d) \n", relativeFactorSize, RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount);
  ${trace.token}  }
  for (k=2; k <= RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount; k++) {
    cdf[k] += cdf[k-1];
  }
  ${trace.token}  if (getTraceFlag(treeID) & FACT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nUpdated CDF based on cardinal group size:  ");
  ${trace.token}    RF_nativePrint("\n     index          cdf");
  ${trace.token}    for (k=1; k <= RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount; k++) {
  ${trace.token}      RF_nativePrint("\n%10d  %12f", k, cdf[k]);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}  }
  randomValue = ceil((ran1B(treeID) * cdf[RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount]));
  randomGroupIndex = 1;
  while (randomValue > cdf[randomGroupIndex]) {
    randomGroupIndex ++;
  }
  free_dvector(cdf, 1, RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount);
  ${trace.token}  if (getTraceFlag(treeID) & FACT_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nRandomly Selected Group Index:  %10d", randomGroupIndex);
  ${trace.token}  }
  createRandomBinaryPair(treeID, relativeFactorSize, absoluteFactorSize, randomGroupIndex, absoluteLevel, result);
  ${trace.token}  if (getTraceFlag(treeID) & FACT_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetRandomPair() EXIT ...\n");
  ${trace.token}  }
}
void createRandomBinaryPair(uint    treeID,
                            uint    relativeFactorSize,
                            uint    absoluteFactorSize,
                            uint    groupIndex,
                            double *absoluteLevel,
                            uint   *pair) {
  uint mwcpLevelIdentifier;
  uint mwcpSizeAbsolute;
  uint offset, levelSize, levelIndex;
  uint k;
  ${trace.token}  if (getTraceFlag(treeID) & FACT_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\ncreateRandomBinaryPair() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & FACT_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nrelativeSize absoluteSize   groupIndex ");
  ${trace.token}    RF_nativePrint("\n%12d %12d %12d", relativeFactorSize, absoluteFactorSize, groupIndex);
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}  }
  levelIndex = 0;  
  mwcpSizeAbsolute = RF_factorList[treeID][absoluteFactorSize] -> mwcpSize;
  uint *levelVector = uivector(1, relativeFactorSize);
  uint *randomLevel = uivector(1, groupIndex);
  for (k = 1; k <= relativeFactorSize; k++) {
    levelVector[k] = k;
  }
  levelSize = relativeFactorSize;
  for (k = 1; k <= groupIndex; k++) {
    randomLevel[k] = sampleUniformlyFromVector(treeID,
                                               levelVector,
                                               levelSize,
                                               &levelIndex);
    levelVector[levelIndex] = levelVector[levelSize];
    levelSize --;
  }
  ${trace.token}  if (getTraceFlag(treeID) & FACT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nAbsolute Levels:  ");
  ${trace.token}    RF_nativePrint("\n     index      level");
  ${trace.token}    for (k=1; k <= relativeFactorSize; k++) {
  ${trace.token}      RF_nativePrint("\n%10d  %10d", k, (uint) absoluteLevel[k]);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}    RF_nativePrint("\nRandomly Selected Levels Prior to Remapping:  ");
  ${trace.token}    RF_nativePrint("\n     index      level");
  ${trace.token}    for (k=1; k <= groupIndex; k++) {
  ${trace.token}      RF_nativePrint("\n%10d  %10d", k, randomLevel[k]);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}  }
  for (k = 1; k <= groupIndex; k++) {
    randomLevel[k] = (uint) absoluteLevel[randomLevel[k]];
  }
  ${trace.token}  if (getTraceFlag(treeID) & FACT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nRandomly Selected Levels After Remapping:  ");
  ${trace.token}    RF_nativePrint("\n     index      level");
  ${trace.token}    for (k=1; k <= groupIndex; k++) {
  ${trace.token}      RF_nativePrint("\n%10d  %10d", k, randomLevel[k]);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}  }
  for (offset = 1; offset <= mwcpSizeAbsolute; offset++) {
    pair[offset] = 0;
  }
  for (k = 1; k <= groupIndex; k++) {
    mwcpLevelIdentifier = (randomLevel[k] >> (3 + ulog2(sizeof(uint)))) + ((randomLevel[k] & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
    ${trace.token}    if (getTraceFlag(treeID) & FACT_HGH_TRACE) {
    ${trace.token}      RF_nativePrint("\n MWCP Level Identifier:   %10d ", mwcpLevelIdentifier);
    ${trace.token}      RF_nativePrint("\n upower() bit:  %10d ", randomLevel[k] - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1 );
    ${trace.token}    }
    pair[mwcpLevelIdentifier] += upower(2, randomLevel[k] - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1 );
  }
  free_uivector(levelVector, 1, relativeFactorSize);
  free_uivector(randomLevel, 1, groupIndex);
  ${trace.token}  if (getTraceFlag(treeID) & FACT_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\ncreateRandomBinaryPair() EXIT ...\n");
  ${trace.token}  }
}
void convertRelToAbsBinaryPair(uint    treeID,
                               uint    relativeFactorSize,
                               uint    absoluteFactorSize,
                               uint    relativePair,
                               double *absoluteLevel,
                               uint   *pair) {
  uint mwcpLevelIdentifier;
  uint mwcpSizeAbsolute;
  uint coercedAbsoluteLevel;
  uint k, offset;
  ${trace.token}  if (getTraceFlag(0) & FACT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nconvertRelToAbsBinaryPair() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(0) & FACT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nrelativeSize absoluteSize      relPair");
  ${trace.token}    RF_nativePrint("\n%12d %12d %12x", relativeFactorSize, absoluteFactorSize, relativePair);
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}  }
  mwcpSizeAbsolute = RF_factorList[treeID][absoluteFactorSize] -> mwcpSize;
  ${trace.token}  if (getTraceFlag(0) & FACT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nAbsolute Levels:  ");
  ${trace.token}    RF_nativePrint("\n     index      level");
  ${trace.token}    for (k=1; k <= relativeFactorSize; k++) {
  ${trace.token}      RF_nativePrint("\n%10d  %10d", k, (uint) absoluteLevel[k]);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\n");
  ${trace.token}    RF_nativePrint("\nRelative Pair:  %8x \n", relativePair);
  ${trace.token}  }
  for (offset = 1; offset <= mwcpSizeAbsolute; offset++) {
    pair[offset] = 0;
  }
  for (k = 1; k <= relativeFactorSize; k++) {
    if (relativePair & ((uint) 0x01)) {
      coercedAbsoluteLevel = (uint) absoluteLevel[k];
      mwcpLevelIdentifier = (coercedAbsoluteLevel >> (3 + ulog2(sizeof(uint)))) + ((coercedAbsoluteLevel & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
      pair[mwcpLevelIdentifier] += upower(2, coercedAbsoluteLevel - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1 );
      ${trace.token}      if (getTraceFlag(0) & FACT_HGH_TRACE) {
      ${trace.token}        RF_nativePrint("\n MWCP Level Identifier:   %10d ", mwcpLevelIdentifier);
      ${trace.token}        RF_nativePrint("\n upower() bit:  %10d ", coercedAbsoluteLevel - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1);
      ${trace.token}      }
    }
    relativePair = relativePair >> 1;
  }
  ${trace.token}  if (getTraceFlag(0) & FACT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nconvertRelToAbsBinaryPair() EXIT ...\n");
  ${trace.token}  }
}
