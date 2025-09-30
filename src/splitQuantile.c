
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "splitQuantile.h"
#include "splitUtil.h"
#include "nrutil.h"
char locallyAdaptiveQuantileRegrSplit (uint       treeID,
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
  uint j, jj, k, p;
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
    uint responseClassCount = RF_quantileSize + 1;
    uint *pseudoResponse = uivector(1, repMembrSize);
    double *sortedResponse = dvector(1, repMembrSize);
    double *quantileValue = dvector(1, RF_quantileSize);
    uint *leftClassProp   = uivector(1, responseClassCount);
    uint *rghtClassProp   = uivector(1, responseClassCount);
    double sumLeft, sumRght, sumLeftSqr, sumRghtSqr;
    double sumLeftMean, sumRghtMean;
    double meanLeft, meanRght;
    double centeredResponse, adaptiveResponse;
    char centerFlag;
    sumLeftMean = sumRghtMean = 0.0;  
    delta = 0.0;  
    centerFlag = TRUE;
    for (j = 1; j <= repMembrSize; j++) {
      sortedResponse[j] = RF_response[treeID][1][ repMembrIndx[j]];
    }
    hpsort(sortedResponse, repMembrSize);
    if (centerFlag) {
      centeredResponse = 0.0;
      for (j = 1; j <= repMembrSize; j++) {
        centeredResponse += sortedResponse[j];
      }
      centeredResponse = centeredResponse / repMembrSize;
      for (j = 1; j <= repMembrSize; j++) {
        sortedResponse[j] = sortedResponse[j] - centeredResponse;
      }
    }
    for (k = 1; k <= RF_quantileSize; k++) {
      quantileValue[k] = quantile7(sortedResponse, repMembrSize, RF_quantile[k]);
    }
    double deltaMax;
    uint   indexMax;
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
          sumLeftMean = sumRghtMean = 0.0;
          for (j = 1; j <= nonMissMembrSize; j++) {
            sumRghtMean += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[j]] ];
          }
          for (j = 1; j <= repMembrSize; j++) {
            localSplitIndicator[j] = RIGHT;
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
              sumLeftMean = sumRghtMean = 0.0;
              for (k = 1; k <= nonMissMembrSize; k++) {
                if (localSplitIndicator[ nonMissMembrIndx[k] ] == LEFT) {
                  sumLeftMean += RF_response[treeID][1][ repMembrIndx[ nonMissMembrIndx[k]] ];
                }
                else {
                  sumRghtMean += RF_response[treeID][1][ repMembrIndx[ nonMissMembrIndx[k]] ];
                }
              }
            }
            else {
              for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                sumLeftMean += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                sumRghtMean -= RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
              }
            }
            meanLeft = sumLeftMean / leftSize;
            meanRght = sumRghtMean / rghtSize;
            for (jj = 1; jj <= nonMissMembrSize; jj++) {
              if (localSplitIndicator[ nonMissMembrIndx[jj] ] == LEFT)  {
                adaptiveResponse = RF_response[treeID][1][ repMembrIndx[ nonMissMembrIndx[jj]] ] - meanLeft;
              }
              else {
                adaptiveResponse = RF_response[treeID][1][ repMembrIndx[ nonMissMembrIndx[jj]] ] - meanRght;
              }
              for (k = 1; k <= RF_quantileSize; k++) {
                if (adaptiveResponse <= quantileValue[k]) {
                  pseudoResponse[ nonMissMembrIndx[jj] ] = k;
                  k = RF_quantileSize;
                }
                else {
                  if (k == RF_quantileSize) {
                    pseudoResponse[ nonMissMembrIndx[jj] ] = k + 1;
                  }
                }
              }
            }
            for (p = 1; p <= responseClassCount; p++) {
              leftClassProp[p] = rghtClassProp[p] = 0;
            }
            for (k = 1; k <= nonMissMembrSize; k++) {
              if (localSplitIndicator[ nonMissMembrIndx[k] ] == LEFT)  {
                leftClassProp[ pseudoResponse[ nonMissMembrIndx[k] ]] ++;
              }
              else {
                rghtClassProp[ pseudoResponse[ nonMissMembrIndx[k] ]] ++;
              }
            }
            sumLeft = sumRght = 0.0;
            for (p=1; p <= responseClassCount; p++) {
              sumLeft += (double) upower(leftClassProp[p], 2);
              sumRght += (double) upower(rghtClassProp[p], 2);
            }
            sumLeftSqr = sumLeft / leftSize;
            sumRghtSqr  = sumRght / rghtSize;
            delta = (sumLeftSqr + sumRghtSqr) / nonMissMembrSize;
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
    free_uivector (leftClassProp,   1, responseClassCount);
    free_uivector (rghtClassProp,   1, responseClassCount);
    free_uivector(pseudoResponse, 1, repMembrSize);
    free_dvector(sortedResponse, 1, repMembrSize);
    free_dvector(quantileValue, 1, RF_quantileSize);
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
char quantileRegrSplit (uint       treeID,
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
    uint responseClassCount = RF_quantileSize + 1;
    uint *pseudoResponse = uivector(1, repMembrSize);
    double *sortedResponse = dvector(1, repMembrSize);
    double *quantileValue = dvector(1, RF_quantileSize);
    uint *parentClassProp = uivector(1, responseClassCount);
    uint *leftClassProp   = uivector(1, responseClassCount);
    uint *rghtClassProp   = uivector(1, responseClassCount);
    double sumLeft, sumRght, sumLeftSqr, sumRghtSqr;
    delta = 0;  
    for (j = 1; j <= repMembrSize; j++) {
      sortedResponse[j] = RF_response[treeID][1][ repMembrIndx[j]];
    }
    hpsort(sortedResponse, repMembrSize);
    for (k = 1; k <= RF_quantileSize; k++) {
      quantileValue[k] = quantile7(sortedResponse, repMembrSize, RF_quantile[k]);
    }
    for (j = 1; j <= repMembrSize; j++) {
      for (k = 1; k <= RF_quantileSize; k++) {
        if (RF_response[treeID][1][ repMembrIndx[j]] <= quantileValue[k]) {
          pseudoResponse[j] = k;
          k = RF_quantileSize;
        }
        else {
          if (k == RF_quantileSize) {
            pseudoResponse[j] = k + 1;
          }
        }
      }
    }
    double deltaMax;
    uint   indexMax;
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
        for (p=1; p <= responseClassCount; p++) {
          parentClassProp[p] = 0;
        }
        for (j = 1; j <= nonMissMembrSize; j++) {
          parentClassProp[ pseudoResponse[ nonMissMembrIndx[j] ]] ++;
        }
        leftSize = 0;
        priorMembrIter = 0;
        if (factorFlag == FALSE) {
          for (j = 1; j <= nonMissMembrSize; j++) {
            localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
          }
          for (p = 1; p <= responseClassCount; p++) {
            rghtClassProp[p] = parentClassProp[p];
            leftClassProp[p] = 0;
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
              }
              for (k = 1; k <= nonMissMembrSize; k++) {
                if (localSplitIndicator[ nonMissMembrIndx[k] ] == LEFT)  {
                  leftClassProp[ pseudoResponse[ nonMissMembrIndx[k] ]] ++;
                }
              }
              for (p=1; p <= responseClassCount; p++) {
                rghtClassProp[p] = parentClassProp[p] - leftClassProp[p];
              }
            }
            else {
              for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                leftClassProp[ pseudoResponse[ nonMissMembrIndx[indxx[k]] ]] ++;
                rghtClassProp[ pseudoResponse[ nonMissMembrIndx[indxx[k]] ]] --;
              }
            }
            sumLeft = sumRght = 0.0;
            for (p=1; p <= responseClassCount; p++) {
              sumLeft += (double) upower(leftClassProp[p], 2);
              sumRght += (double) upower(rghtClassProp[p], 2);
            }
            sumLeftSqr = sumLeft / leftSize;
            sumRghtSqr  = sumRght / rghtSize;
            delta = (sumLeftSqr + sumRghtSqr) / nonMissMembrSize;
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
    free_uivector(pseudoResponse, 1, repMembrSize);
    free_dvector(sortedResponse, 1, repMembrSize);
    free_dvector(quantileValue, 1, RF_quantileSize);
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
double quantile7 (double *r, uint s, double p) {
  double result;
  double delta;
  uint i;
  i = floor(1 + ((s-1) * p));
  delta = 1.0 + ((s-1) * p) - i;
  result = ((1.0 - delta) * r[i]) + (delta * r[i+1]);
  return result;
}
