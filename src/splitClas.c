
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "splitClas.h"
#include "splitUtil.h"
#include "splitInfo.h"
#include "node.h"
#include "nrutil.h"
char classificationXwghtSplitCur (uint       treeID,
                                  Node      *parent,
                                  SplitInfoMax *splitInfoMax,
                                  GreedyObj    *greedyMembr,
                                  char       multImpFlag) {
  uint     covariate;
  uint     covariateCount;
  double  *splitVector;
  uint     splitVectorSize;
  uint   *indxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  uint splitLength;
  void *splitVectorPtr;
  double *observation;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char preliminaryResult, result;
  double delta;
  uint j, k, p;
  mwcpSizeAbsolute       = 0;     
  indxx = NULL;
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        multImpFlag,
                                        FALSE);
  if (preliminaryResult) {
    uint  repMembrSize = parent -> repMembrSize;
    uint *repMembrIndx = parent -> repMembrIndx;
    uint  nonMissMembrSize;
    uint *nonMissMembrIndx;
    stackSplitPreliminary(repMembrSize,
                          & localSplitIndicator,
                          & splitVector);
    DistributionObj *distributionObj = stackRandomCovariates(treeID, parent);
    uint responseClassCount = RF_classLevelSize[1];
    uint *parentClassProp = uivector(1, responseClassCount);
    uint *leftClassProp   = uivector(1, responseClassCount);
    uint *rghtClassProp   = uivector(1, responseClassCount);
    double sumLeft, sumRght;
    double leftTemp1, rghtTemp1;
    delta = 0.0;  
    double deltaMax;
    uint   indexMax;
    if ((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP)) {
      for (p=1; p <= responseClassCount; p++) {
        parentClassProp[p] = 0;
      }
      for (j = 1; j <= repMembrSize; j++) {
        parentClassProp[RF_classLevelIndex[1][ (uint) RF_response[treeID][1][ repMembrIndx[j] ]]] ++;
      }
    }
    covariateCount = 0;
    while (selectRandomCovariates(treeID,
                                  parent,
                                  distributionObj,
                                  & factorFlag,
                                  & covariate,
                                  & covariateCount)) {
      splitVectorSize = stackAndConstructSplitVectorGenericPhase1(treeID,
                                                                    parent,
                                                                    covariate,
                                                                    splitVector,
                                                                  & indxx,
                                                                  multImpFlag);
      if (splitVectorSize >= 2) {
        splitLength = stackAndConstructSplitVectorGenericPhase2(treeID,
                                                                parent,
                                                                covariate,
                                                                splitVector,
                                                                splitVectorSize,
                                                                & factorFlag,
                                                                & deterministicSplitFlag,
                                                                & mwcpSizeAbsolute,
                                                                & splitVectorPtr);
        nonMissMembrIndx = parent -> nonMissMembrIndx;
        nonMissMembrSize = parent -> nonMissMembrSize;
        observation = RF_observation[treeID][covariate];
        leftSize = 0;
        priorMembrIter = 0;
        if (factorFlag == FALSE) {
          if ((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP)) {
            for (p = 1; p <= responseClassCount; p++) {
              rghtClassProp[p] = parentClassProp[p];
              leftClassProp[p] = 0;
            }
            for (j = 1; j <= repMembrSize; j++) {
              localSplitIndicator[j] = RIGHT;
            }    
          }
          else {
            for (p=1; p <= responseClassCount; p++) {
              rghtClassProp[p] = 0;
              leftClassProp[p] = 0;
            }
            for (j = 1; j <= nonMissMembrSize; j++) {
              rghtClassProp[RF_classLevelIndex[1][ (uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[j]] ]]] ++;
            }
            for (j = 1; j <= nonMissMembrSize; j++) {
              localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
            }
          }
        }
        deltaMax = RF_nativeNaN;
        indexMax =  0;
        for (j = 1; j < splitLength; j++) {
          if (factorFlag == TRUE) {
            priorMembrIter = 0;
            leftSize = 0;
          }
          virtuallySplitNode(treeID,
                             parent,
                             factorFlag,
                             mwcpSizeAbsolute,
                             observation,
                             indxx,
                             splitVectorPtr,
                             j,
                             localSplitIndicator,
                             & leftSize,
                             priorMembrIter,
                             & currentMembrIter);
          rghtSize = nonMissMembrSize - leftSize;
          if ((leftSize != 0) && (rghtSize != 0)) {
            if (factorFlag == TRUE) {
              for (p=1; p <= responseClassCount; p++) {
                leftClassProp[p] = 0;
                rghtClassProp[p] = 0;
              }
              for (k = 1; k <= nonMissMembrSize; k++) {
                if (localSplitIndicator[ nonMissMembrIndx[k] ] == LEFT)  {
                  leftClassProp[RF_classLevelIndex[1][ (uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ]]] ++;
                }
                else {
                  rghtClassProp[RF_classLevelIndex[1][ (uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ]]] ++;
                }
              }
            }
            else {
              for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                leftClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] ++;
                rghtClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] --;
              }
            }
            sumLeft = sumRght = 0.0;
            switch(RF_splitRule) {
            case CLAS_WT_OFF:
              for (p=1; p <= responseClassCount; p++) {
                sumLeft += pow ((double) leftClassProp[p] / (double) leftSize, 2.0);
                sumRght += pow ((double) (rghtClassProp[p]) / (double) rghtSize, 2.0);
              }
              delta = sumLeft + sumRght;
              break;
            case CLAS_WT_HVY:
              for (p=1; p <= responseClassCount; p++) {
                sumLeft += (double) upower(leftClassProp[p], 2);
                sumRght += (double) upower(rghtClassProp[p], 2);
              }
              delta =
                (sumLeft / (double) (upower(nonMissMembrSize, 2))) +
                (sumRght / (double) (upower(nonMissMembrSize, 2))) -
                pow ((double) leftSize / nonMissMembrSize, 2.0) -
                pow ((double) rghtSize / nonMissMembrSize, 2.0) + 2.0;
              break;
            default:
              for (p=1; p <= responseClassCount; p++) {
                sumLeft += (double) upower(leftClassProp[p], 2);
                sumRght += (double) upower(rghtClassProp[p], 2);
              }
              leftTemp1 = sumLeft / leftSize;
              rghtTemp1 = sumRght / rghtSize;
              delta = (leftTemp1 + rghtTemp1) / nonMissMembrSize;
              break;
            }
          }
          else {
            delta = RF_nativeNaN;
          }
          if (!RF_nativeIsNaN(delta)) {
            if(RF_nativeIsNaN(deltaMax)) {
              deltaMax = delta;
              indexMax = j;
            }
            else {
              if ((delta - deltaMax) > EPSILON) {
                deltaMax = delta;
                indexMax = j;
              }
            }
          }
          if (factorFlag == FALSE) {
            priorMembrIter = currentMembrIter - 1;
          }
        }  
        updateMaximumSplit(treeID,
                             parent,
                             deltaMax,
                             covariate,
                             indexMax,
                             factorFlag,
                             mwcpSizeAbsolute,
                             repMembrSize,
                             & localSplitIndicator,
                             splitVectorPtr,
                             splitInfoMax);
        unstackSplitVector(treeID,
                           parent,
                           splitLength,
                           factorFlag,            
                           splitVectorSize,
                           mwcpSizeAbsolute,
                           deterministicSplitFlag,
                           splitVectorPtr,
                           multImpFlag,
                           indxx);
      }
    }  
    unstackRandomCovariates(treeID, distributionObj);
    free_uivector (parentClassProp, 1, responseClassCount);
    free_uivector (leftClassProp,   1, responseClassCount);
    free_uivector (rghtClassProp,   1, responseClassCount);
    unstackSplitPreliminary(repMembrSize,
                            localSplitIndicator,
                            splitVector);
  }  
  unstackPreSplit(preliminaryResult,
                  parent,
                  multImpFlag,
                  FALSE);  
  result = summarizeSplitResult(splitInfoMax);
  return result;
}
char classificationAreaUnderROCSplit (uint       treeID,
                                      Node      *parent,
                                      SplitInfoMax *splitInfoMax,
                                      GreedyObj    *greedyMembr,
                                      char       multImpFlag) {
  uint     covariate;
  uint     covariateCount;
  double  *splitVector;
  uint     splitVectorSize;
  uint   *indxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  uint splitLength;
  void *splitVectorPtr;
  double *observation;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char preliminaryResult, result;
  double delta;
  uint j, k, p;
  mwcpSizeAbsolute       = 0;     
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        multImpFlag,
                                        FALSE);
  if (preliminaryResult) {
    uint  repMembrSize = parent -> repMembrSize;
    uint *repMembrIndx = parent -> repMembrIndx;
    uint  nonMissMembrSize;
    uint *nonMissMembrIndx;
    stackSplitPreliminary(repMembrSize,
                          & localSplitIndicator,
                          & splitVector);
    DistributionObj *distributionObj = stackRandomCovariates(treeID, parent);
    uint responseClassCount = RF_classLevelSize[1];
    uint *parentClassProp = uivector(1, responseClassCount);
    uint *leftClassProp   = uivector(1, responseClassCount);
    uint *rghtClassProp   = uivector(1, responseClassCount);
    double alpha, beta;
    uint sLen;
    delta = 0;  
    double deltaMax;
    uint   indexMax;
    if (responseClassCount == 2) {
      sLen = 0;  
    }
    else {
      sLen = responseClassCount * (responseClassCount - 1) / 2;
    }
    if ((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP)) {
      for (p=1; p <= responseClassCount; p++) {
        parentClassProp[p] = 0;
      }
      for (j = 1; j <= repMembrSize; j++) {
        parentClassProp[RF_classLevelIndex[1][ (uint) RF_response[treeID][1][ repMembrIndx[j] ]]] ++;
      }
    }
    covariateCount = 0;
    while (selectRandomCovariates(treeID,
                                  parent,
                                  distributionObj,
                                  & factorFlag,
                                  & covariate,
                                  & covariateCount)) {
      splitVectorSize = stackAndConstructSplitVectorGenericPhase1(treeID,
                                                                  parent,
                                                                  covariate,
                                                                  splitVector,
                                                                  & indxx,
                                                                  multImpFlag);
      if (splitVectorSize >= 2) {
        splitLength = stackAndConstructSplitVectorGenericPhase2(treeID,
                                                                parent,
                                                                covariate,
                                                                splitVector,
                                                                splitVectorSize,
                                                                & factorFlag,
                                                                & deterministicSplitFlag,
                                                                & mwcpSizeAbsolute,
                                                                & splitVectorPtr);
        nonMissMembrIndx = parent -> nonMissMembrIndx;
        nonMissMembrSize = parent -> nonMissMembrSize;
        observation = RF_observation[treeID][covariate];
        leftSize = 0;
        priorMembrIter = 0;
        if (factorFlag == FALSE) {
          if ((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP)) {
            for (p = 1; p <= responseClassCount; p++) {
              rghtClassProp[p] = parentClassProp[p];
              leftClassProp[p] = 0;
            }
            for (j = 1; j <= nonMissMembrSize; j++) {
              localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
            }
          }
          else {
            for (p=1; p <= responseClassCount; p++) {
              rghtClassProp[p] = 0;
              leftClassProp[p] = 0;
            }
            for (j = 1; j <= nonMissMembrSize; j++) {
              rghtClassProp[RF_classLevelIndex[1][ (uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[j]] ]]] ++;
            }
            for (j = 1; j <= nonMissMembrSize; j++) {
              localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
            }
          }
        }
        deltaMax = RF_nativeNaN;
        indexMax =  0;
        for (j = 1; j < splitLength; j++) {
          if (factorFlag == TRUE) {
            priorMembrIter = 0;
            leftSize = 0;
          }
          virtuallySplitNode(treeID,
                             parent,
                             factorFlag,
                             mwcpSizeAbsolute,
                             observation,
                             indxx,
                             splitVectorPtr,
                             j,
                             localSplitIndicator,
                             & leftSize,
                             priorMembrIter,
                             & currentMembrIter);
          rghtSize = nonMissMembrSize - leftSize;
          if ((leftSize != 0) && (rghtSize != 0)) {
            if (factorFlag == TRUE) {
              for (p=1; p <= responseClassCount; p++) {
                leftClassProp[p] = 0;
                rghtClassProp[p] = 0;
              }
              for (k = 1; k <= nonMissMembrSize; k++) {
                if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT)  {
                  leftClassProp[RF_classLevelIndex[1][ (uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] ++;
                }
                else {
                  rghtClassProp[RF_classLevelIndex[1][ (uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ]]] ++;
                }
              }
            }
            else {
              for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                leftClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] ++;
                rghtClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] --;
              }
            }
            if (responseClassCount == 2) {
              alpha = (double) rghtClassProp[1] / (double) (leftClassProp[1] + rghtClassProp[1]);
              beta  = (double) rghtClassProp[2] / (double) (leftClassProp[2] + rghtClassProp[2]);
              delta = 1.0 - alpha + beta;
              if (delta < (1.0 - beta + alpha)) {
                delta =  1.0 - beta + alpha;
              }
            }
            else {
              delta = 0.0;
              for (uint k1 = 1; k1 <= responseClassCount - 1; k1++) {
                for (uint k2 = k1+1; k2 <= responseClassCount; k2++) {
                  alpha  = (double) rghtClassProp[k1] / (double) (leftClassProp[k1] + rghtClassProp[k1]);
                  beta  = (double) rghtClassProp[k2] / (double) (leftClassProp[k2] + rghtClassProp[k2]);
                  if ((1.0 - alpha + beta) > (1.0 - beta + alpha)) {
                    delta +=  1.0 - alpha + beta;
                  }
                  else {
                    delta +=  1.0 - beta + alpha;
                  }
                }
              }
              delta = delta / sLen;
            }
          }
          else {
            delta = RF_nativeNaN;
          }
          if (!RF_nativeIsNaN(delta)) {
            if(RF_nativeIsNaN(deltaMax)) {
              deltaMax = delta;
              indexMax = j;
            }
            else {
              if ((delta - deltaMax) > EPSILON) {
                deltaMax = delta;
                indexMax = j;
              }
            }
          }
          if (factorFlag == FALSE) {
            priorMembrIter = currentMembrIter - 1;
          }
        }  
        updateMaximumSplit(treeID,
                           parent,
                           deltaMax,
                           covariate,
                           indexMax,
                           factorFlag,
                           mwcpSizeAbsolute,
                           repMembrSize,
                           & localSplitIndicator,
                           splitVectorPtr,
                           splitInfoMax);
        unstackSplitVector(treeID,
                              parent,
                              splitLength,
                              factorFlag,            
                              splitVectorSize,
                              mwcpSizeAbsolute,
                              deterministicSplitFlag,
                              splitVectorPtr,
                              multImpFlag,
                              indxx);
      }
    }  
    unstackRandomCovariates(treeID, distributionObj);
    free_uivector (parentClassProp, 1, responseClassCount);
    free_uivector (leftClassProp,   1, responseClassCount);
    free_uivector (rghtClassProp,   1, responseClassCount);
    unstackSplitPreliminary(repMembrSize,
                            localSplitIndicator,
                            splitVector);
  }  
  unstackPreSplit(preliminaryResult,
                  parent,
                  multImpFlag,
                  FALSE);  
  result = summarizeSplitResult(splitInfoMax);
  return result;
}
char classificationEntropySplit (uint       treeID,
                                 Node      *parent,
                                 SplitInfoMax *splitInfoMax,
                                 GreedyObj    *greedyMembr,
                                 char       multImpFlag) {
  uint     covariate;
  uint     covariateCount;
  double  *splitVector;
  uint     splitVectorSize;
  uint   *indxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  uint splitLength;
  void *splitVectorPtr;
  double *observation;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char preliminaryResult, result;
  double delta, deltaLeft, deltaRght;
  uint j, k, p;
  mwcpSizeAbsolute       = 0;     
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        multImpFlag,
                                        FALSE);
  if (preliminaryResult) {
    uint  repMembrSize = parent -> repMembrSize;
    uint *repMembrIndx = parent -> repMembrIndx;
    uint  nonMissMembrSize;
    uint *nonMissMembrIndx;
    stackSplitPreliminary(repMembrSize,
                          & localSplitIndicator,
                          & splitVector);
    DistributionObj *distributionObj = stackRandomCovariates(treeID, parent);
    uint responseClassCount = RF_classLevelSize[1];
    uint *parentClassProp = uivector(1, responseClassCount);
    uint *leftClassProp   = uivector(1, responseClassCount);
    uint *rghtClassProp   = uivector(1, responseClassCount);
    uint *leftClassRatio  = uivector(1, responseClassCount);
    uint *rghtClassRatio  = uivector(1, responseClassCount);
    delta = deltaLeft = deltaRght = 0;  
    double deltaMax;
    uint   indexMax;
    if ((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP)) {
      for (p=1; p <= responseClassCount; p++) {
        parentClassProp[p] = 0;
      }
      for (j = 1; j <= repMembrSize; j++) {
        parentClassProp[RF_classLevelIndex[1][ (uint) RF_response[treeID][1][ repMembrIndx[j] ]]] ++;
      }
    }
    covariateCount = 0;
    while (selectRandomCovariates(treeID,
                                  parent,
                                  distributionObj,
                                  & factorFlag,
                                  & covariate,
                                  & covariateCount)) {
      splitVectorSize = stackAndConstructSplitVectorGenericPhase1(treeID,
                                                                  parent,
                                                                  covariate,
                                                                  splitVector,
                                                                  & indxx,
                                                                  multImpFlag);
      if (splitVectorSize >= 2) {
        splitLength = stackAndConstructSplitVectorGenericPhase2(treeID,
                                                                parent,
                                                                covariate,
                                                                splitVector,
                                                                splitVectorSize,
                                                                & factorFlag,
                                                                & deterministicSplitFlag,
                                                                & mwcpSizeAbsolute,
                                                                & splitVectorPtr);
        nonMissMembrIndx = parent -> nonMissMembrIndx;
        nonMissMembrSize = parent -> nonMissMembrSize;
        observation = RF_observation[treeID][covariate];
        leftSize = 0;
        priorMembrIter = 0;
        if (factorFlag == FALSE) {
          if ((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP)) {
            for (p = 1; p <= responseClassCount; p++) {
              rghtClassProp[p] = parentClassProp[p];
              leftClassProp[p] = 0;
            }
            for (j = 1; j <= nonMissMembrSize; j++) {
              localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
            }
          }
          else {
            for (p=1; p <= responseClassCount; p++) {
              rghtClassProp[p] = 0;
              leftClassProp[p] = 0;
            }
            for (j = 1; j <= nonMissMembrSize; j++) {
              rghtClassProp[RF_classLevelIndex[1][ (uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[j]] ]]] ++;
            }
            for (j = 1; j <= nonMissMembrSize; j++) {
              localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
            }
          }
        }
        deltaMax = RF_nativeNaN;
        indexMax =  0;
        for (j = 1; j < splitLength; j++) {
          if (factorFlag == TRUE) {
            priorMembrIter = 0;
            leftSize = 0;
          }
          virtuallySplitNode(treeID,
                             parent,
                             factorFlag,
                             mwcpSizeAbsolute,
                             observation,
                             indxx,
                             splitVectorPtr,
                             j,
                             localSplitIndicator,
                             & leftSize,
                             priorMembrIter,
                             & currentMembrIter);
          rghtSize = nonMissMembrSize - leftSize;
          if ((leftSize != 0) && (rghtSize != 0)) {
            if (factorFlag == TRUE) {
              for (p=1; p <= responseClassCount; p++) {
                leftClassProp[p] = 0;
                rghtClassProp[p] = 0;
              }
              for (k = 1; k <= nonMissMembrSize; k++) {
                if (localSplitIndicator[ nonMissMembrIndx[k] ] == LEFT)  {
                  leftClassProp[RF_classLevelIndex[1][ (uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ]]] ++;
                }
                else {
                  rghtClassProp[RF_classLevelIndex[1][ (uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ]]] ++;
                }
              }
            }
            else {
              for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                leftClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] ++;
                rghtClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] --;
              }
            }
            deltaLeft = deltaRght = 0.0;
            for (uint k1 = 1; k1 <= responseClassCount; k1++) {
              leftClassRatio[k1] = leftClassProp[k1] / (double) leftSize;
              rghtClassRatio[k1] = rghtClassProp[k1] / (double) rghtSize;
              if (leftClassRatio[k1] > 0) {
                deltaLeft += leftClassRatio[k1] / log(leftClassRatio[k1]);
              }
              if (rghtClassRatio[k1] > 0) {
                deltaRght += rghtClassRatio[k1] / log(rghtClassRatio[k1]);
              }
            }
            delta = (deltaLeft * leftSize / (double) repMembrSize) + (deltaRght * rghtSize / (double) repMembrSize);
          }          
          else {
            delta = RF_nativeNaN;
          }
          if (!RF_nativeIsNaN(delta)) {
            if(RF_nativeIsNaN(deltaMax)) {
              deltaMax = delta;
              indexMax = j;
            }
            else {
              if ((delta - deltaMax) > EPSILON) {
                deltaMax = delta;
                indexMax = j;
              }
            }
          }
          if (factorFlag == FALSE) {
            priorMembrIter = currentMembrIter - 1;
          }
        }  
        updateMaximumSplit(treeID,
                           parent,
                           deltaMax,
                           covariate,
                           indexMax,
                           factorFlag,
                           mwcpSizeAbsolute,
                           repMembrSize,
                           & localSplitIndicator,
                           splitVectorPtr,
                           splitInfoMax);
        unstackSplitVector(treeID,
                              parent,
                              splitLength,
                              factorFlag,            
                              splitVectorSize,
                              mwcpSizeAbsolute,
                              deterministicSplitFlag,
                              splitVectorPtr,
                              multImpFlag,
                              indxx);
      }  
    }  
    unstackRandomCovariates(treeID, distributionObj);
    free_uivector (parentClassProp, 1, responseClassCount);
    free_uivector (leftClassProp,   1, responseClassCount);
    free_uivector (rghtClassProp,   1, responseClassCount);
    free_uivector (leftClassRatio,  1, responseClassCount);
    free_uivector (rghtClassRatio,  1, responseClassCount);
    unstackSplitPreliminary(repMembrSize,
                            localSplitIndicator,
                            splitVector);
  }  
  unstackPreSplit(preliminaryResult,
                  parent,
                  multImpFlag,
                  FALSE);  
  result = summarizeSplitResult(splitInfoMax);
  return result;
}
