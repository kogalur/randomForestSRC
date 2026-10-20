
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
  result = TRUE;
  if (result) {
    if (parent -> repMembrSize >= (2 * RF_nodeSize)) {
      result = TRUE;
    }
    else {
      result = FALSE;
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
      }
    }
  }
  if (result) {
    if ((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP) || multVarFlag) {
      parent -> nonMissMembrSizeStatic = parent -> repMembrSize;
      parent -> nonMissMembrIndxStatic = RF_identityMembershipIndex;
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
  return result;
}
char getPreSplitResultNoMiss (uint      treeID,
                              Node     *parent,
                              char      multImpFlag,
                              char      multVarFlag) {
  uint i;
  char result;
  result = TRUE;
  if (result) {
    if (parent -> repMembrSize >= (2 * RF_nodeSize)) {
      result = TRUE;
    }
    else {
      result = FALSE;
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
  return result;
}
void unstackPreSplit (char      preliminaryResult,
                      Node     *parent,
                      char      multImpFlag,
                      char      multVarFlag) {
  if (preliminaryResult) {
    if (!((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP) || multVarFlag)) {
      free_uivector(parent -> nonMissMembrIndxStatic, 1, parent -> repMembrSize);
    }
  }
  else {
  }
}
void stackSplitPreliminary(uint     nodeSize,
                           char   **localSplitIndicator,
                           double **splitVector) {
  *localSplitIndicator = cvector(1, nodeSize);
  *splitVector = dvector(1, nodeSize);
}
void unstackSplitPreliminary(uint    nodeSize,
                             char   *localSplitIndicator,
                             double *splitVector) {
  free_cvector(localSplitIndicator, 1, nodeSize);
  free_dvector(splitVector, 1, nodeSize);
}
DistributionObj *stackRandomCovariatesGeneric(uint treeID, Node *parent) {
  uint actualWeightType;
  uint *augmentationSize;
  char *permissible;
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
  return obj;
}
void unstackRandomCovariatesGeneric(uint treeID, DistributionObj *obj) {
  if (obj -> augmentationSize != NULL) {
    free_uivector(obj -> augmentationSize, 1, 2);
    obj -> augmentationSize = NULL;
  }
  discardCDFNew(treeID, obj);
  freeDistributionObjRaw(obj);
}
char selectRandomCovariatesGeneric(uint     treeID,
                                   Node     *parent,
                                   DistributionObj *distributionObj,
                                   char     *factorFlag,
                                   uint     *covariate,
                                   uint     *covariateCount) {
  char xVarFound;
  (*covariate) = UINT_MAX;
  xVarFound = FALSE;
  *factorFlag = FALSE;
  while ( ((*covariateCount) < RF_mtry) && ((*covariate) != 0) && (xVarFound == FALSE)) {
    (*covariateCount) ++;
    *covariate = sampleFromCDFNew(ran1B, treeID, distributionObj);
    if (*covariate != 0) {
      updateCDFNew(treeID, distributionObj);
      xVarFound = TRUE;
        if (RF_xType[*covariate] == 'C') {
          *factorFlag = TRUE;
        }
    }  
  }  
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
    if(vectorSize >= 2) {
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
    }
  }
  else {
    vectorSize = 0;
    if (covariate <= RF_xSize) {
      (parent -> permissible)[covariate] = FALSE;
      parent -> permissibleReIndxFlag = TRUE;
    }
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
  repMembrSize = parent -> repMembrSize;
  splitLength = 0;  
  (*splitVectorPtr) = NULL;  
  if (vectorSize < 2) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Split vector is of insufficient size in Stack Phase II allocation:  %10d", vectorSize);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  if (*factorFlag) {
    if(RF_factorList[treeID][vectorSize] == NULL) {
      RF_factorList[treeID][vectorSize] = makeFactor(vectorSize, FALSE);
    }
    factorSizeAbsolute = RF_xFactorSize[RF_xFactorMap[covariate]];
    *mwcpSizeAbsolute = RF_factorList[treeID][factorSizeAbsolute] -> mwcpSize;
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
  iterFlag = TRUE;
  repMembrIndx = parent -> repMembrIndx;
  nonMissMembrIndx = parent -> nonMissMembrIndx;
  nonMissMembrSize = parent -> nonMissMembrSize;
  *currentMembrIter = priorMembrIter;
  while (iterFlag) {
    (*currentMembrIter) ++;
    if (factorFlag == TRUE) {
      daughterFlag = splitOnFactor((uint)  observation[    repMembrIndx[nonMissMembrIndx[*currentMembrIter]]     ],
                                   (uint*) splitVectorPtr + ((offset - 1) * mwcpSizeAbsolute));
      localSplitIndicator[     nonMissMembrIndx[*currentMembrIter]   ] = daughterFlag;
      if ((*currentMembrIter) == nonMissMembrSize) {
        iterFlag = FALSE;
      }
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
    }
    if (daughterFlag == LEFT) {
      (*leftSize) ++;
    }
  }  
  if ((*leftSize == 0) || (*leftSize == nonMissMembrSize)) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Left or Right Daughter of size zero:  (%10d, %10d)", *leftSize, nonMissMembrSize);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  return (*leftSize);
}
char summarizeSplitResult(SplitInfoMax *splitInfoMax) {
  char result;
  if (!RF_nativeIsNaN(splitInfoMax -> deltaMax)) {
    splitInfoMax -> splitStatistic = splitInfoMax -> deltaMax;
    result = TRUE;
  }
  else {
    result = FALSE;
  }
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
    }
  }
  else {
  }
  return flag;
}
void getReweightedRandomPair (uint    treeID,
                              uint    relativeFactorSize,
                              uint    absoluteFactorSize,
                              double *absoluteLevel,
                              uint   *result) {
  uint randomGroupIndex;
  if(RF_factorList[treeID][relativeFactorSize] == NULL) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Factor not allocated for size:  %10d", relativeFactorSize);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  randomGroupIndex = (uint) ceil(ran1B(treeID) * ((RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount) * 1.0));
  createRandomBinaryPair(treeID, relativeFactorSize, absoluteFactorSize, randomGroupIndex, absoluteLevel, result);
}
void getRandomPair (uint treeID, uint relativeFactorSize, uint absoluteFactorSize, double *absoluteLevel, uint *result) {
  uint randomGroupIndex;
  double randomValue;
  uint k;
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
  for (k=2; k <= RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount; k++) {
    cdf[k] += cdf[k-1];
  }
  randomValue = ceil((ran1B(treeID) * cdf[RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount]));
  randomGroupIndex = 1;
  while (randomValue > cdf[randomGroupIndex]) {
    randomGroupIndex ++;
  }
  free_dvector(cdf, 1, RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount);
  createRandomBinaryPair(treeID, relativeFactorSize, absoluteFactorSize, randomGroupIndex, absoluteLevel, result);
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
  for (k = 1; k <= groupIndex; k++) {
    randomLevel[k] = (uint) absoluteLevel[randomLevel[k]];
  }
  for (offset = 1; offset <= mwcpSizeAbsolute; offset++) {
    pair[offset] = 0;
  }
  for (k = 1; k <= groupIndex; k++) {
    mwcpLevelIdentifier = (randomLevel[k] >> (3 + ulog2(sizeof(uint)))) + ((randomLevel[k] & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
    pair[mwcpLevelIdentifier] += upower(2, randomLevel[k] - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1 );
  }
  free_uivector(levelVector, 1, relativeFactorSize);
  free_uivector(randomLevel, 1, groupIndex);
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
  mwcpSizeAbsolute = RF_factorList[treeID][absoluteFactorSize] -> mwcpSize;
  for (offset = 1; offset <= mwcpSizeAbsolute; offset++) {
    pair[offset] = 0;
  }
  for (k = 1; k <= relativeFactorSize; k++) {
    if (relativePair & ((uint) 0x01)) {
      coercedAbsoluteLevel = (uint) absoluteLevel[k];
      mwcpLevelIdentifier = (coercedAbsoluteLevel >> (3 + ulog2(sizeof(uint)))) + ((coercedAbsoluteLevel & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
      pair[mwcpLevelIdentifier] += upower(2, coercedAbsoluteLevel - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1 );
    }
    relativePair = relativePair >> 1;
  }
}
