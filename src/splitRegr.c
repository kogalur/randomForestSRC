
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "splitRegr.h"
#include "splitUtil.h"
#include "nrutil.h"
char regressionXwghtSplitCur (uint       treeID,
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
  uint j, k;
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
    double sumLeft, sumRght, sumLeftSqr, sumRghtSqr;
    double sumRghtSave, sumRghtSqrSave;
    double leftTemp1, rghtTemp1, leftTemp2, rghtTemp2;
    sumLeft = sumRght = sumLeftSqr = sumRghtSqr = 0.0;  
    sumRghtSave = sumRghtSqrSave = 0.0;                 
    delta = 0.0;                                        
    double deltaMax;
    uint   indexMax;
    if ((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP)) {
      sumRghtSave = (parent -> mean) * repMembrSize;
      switch(RF_splitRule) {
      case REGR_WT_OFF:
      case REGR_WT_HVY:
        sumRghtSqrSave = 0.0;
        for (j = 1; j <= repMembrSize; j++) {
          sumRghtSqrSave += pow (RF_response[treeID][1][ repMembrIndx[j] ], 2.0);
        }
        break;
      default:
        break;
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
          sumLeft = sumLeftSqr = 0.0;
          if ((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP)) {
            sumRght = sumRghtSave;
            sumRghtSqr = sumRghtSqrSave;
            for (j = 1; j <= repMembrSize; j++) {
              localSplitIndicator[j] = RIGHT;
            }
          }
          else {
            switch(RF_splitRule) {
            case REGR_WT_OFF:
            case REGR_WT_HVY:
              sumRght = sumRghtSqr = 0.0;
              for (j = 1; j <= nonMissMembrSize; j++) {
                sumRght += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[j]] ];
                sumRghtSqr += pow (RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[j]] ], 2.0);
              }
              break;
            default:
              sumRght = 0.0;
              for (j = 1; j <= nonMissMembrSize; j++) {
                sumRght += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[j]] ];
              }
              break;
            }
            for (j = 1; j <= nonMissMembrSize; j++) {
              localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
            }
          }          
        }
        deltaMax =  RF_nativeNaN;
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
              switch (RF_splitRule) {
              case REGR_WT_OFF:
              case REGR_WT_HVY:
                sumLeft = sumRght = 0.0;
                sumLeftSqr = sumRghtSqr = 0.0;
                for (k = 1; k <= nonMissMembrSize; k++) {
                  if (localSplitIndicator[ nonMissMembrIndx[k] ] == LEFT) {
                    sumLeft += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ];
                    sumLeftSqr += pow (RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ], 2.0);
                  }
                  else {
                    sumRght += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ];
                    sumRghtSqr += pow (RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ], 2.0);
                  }
                } 
                break;
              default:
                sumLeft = sumRght = 0.0;
                for (k = 1; k <= nonMissMembrSize; k++) {
                  if (localSplitIndicator[ nonMissMembrIndx[k] ] == LEFT) {
                    sumLeft += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ];
                  }
                  else {
                    sumRght += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ];
                  }
                } 
                break;
              }
            }
            else {
              switch(RF_splitRule) {
              case REGR_WT_OFF:
              case REGR_WT_HVY:
                for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                  sumLeft    += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                  sumLeftSqr += pow (RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ], 2.0);
                  sumRght    -= RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                  sumRghtSqr -= pow (RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ], 2.0);
                }
                break;
              default:
                for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                  sumLeft += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                  sumRght -= RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                }
                break;
              }
            }
            switch(RF_splitRule) {
            case REGR_WT_OFF:
              leftTemp1 = pow (sumLeft, 2.0) / (double) upower (leftSize, 2); 
              rghtTemp1 = pow (sumRght, 2.0) / (double) upower (rghtSize, 2); 
              leftTemp2 = sumLeftSqr / leftSize; 
              rghtTemp2 = sumRghtSqr / rghtSize; 
              delta = leftTemp1 + rghtTemp1 - leftTemp2 - rghtTemp2;
              break;
            case REGR_WT_HVY:
              leftTemp1 = pow (sumLeft, 2.0) / (double) upower (nonMissMembrSize, 2);
              rghtTemp1 = pow (sumRght, 2.0) / (double) upower (nonMissMembrSize, 2);
              leftTemp2 = sumLeftSqr * leftSize / (double) upower (nonMissMembrSize, 2);
              rghtTemp2 = sumRghtSqr * rghtSize / (double) upower (nonMissMembrSize, 2);
              delta = leftTemp1 + rghtTemp1 - leftTemp2 - rghtTemp2;
              break;
            default:
              leftTemp1 = pow (sumLeft, 2.0) / leftSize;
              rghtTemp1 = pow (sumRght, 2.0) / rghtSize;
              delta = leftTemp1 + rghtTemp1;
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
              else {
              }
            }
          }
          else {
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
