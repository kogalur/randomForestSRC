
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "splitRegr.h"
#include "splitUtil.h"
#include "nrutil.h"
${trace.token} #include "error.h"
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
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nregressionXwghtSplitCur(%10d) ENTRY ...\n", treeID);
  ${trace.token}  }
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
          ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
          ${trace.token}    RF_nativePrint("\nNon-miss Node Size:  %10d, Non-miss Left Size:  %10d, Non-miss Right Size:  %10d", nonMissMembrSize, leftSize, rghtSize);
          ${trace.token}  }
          ${trace.token}        if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
          ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
          ${trace.token}            RF_nativePrint("\n PriorIter:     %10d  CurrentIter:   %10d", priorMembrIter, currentMembrIter);
          ${trace.token}          }
          ${trace.token}        }
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
            ${trace.token}          if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
            ${trace.token}            if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
            ${trace.token}              RF_nativePrint("\nVirtual (non-miss) Membership:  ");
            ${trace.token}              if (factorFlag) {
            ${trace.token}                for (k = 1; k <= nonMissMembrSize; k++) {
            ${trace.token}                  if (localSplitIndicator[   nonMissMembrIndx[k]  ] == LEFT) {
            ${trace.token}                    RF_nativePrint("\n %10d %10d %10d %12.4f %12.4f --> LEFT ", k, nonMissMembrIndx[k],  repMembrIndx[nonMissMembrIndx[k]], observation[ repMembrIndx[nonMissMembrIndx[k]] ], RF_response[treeID][1][  repMembrIndx[nonMissMembrIndx[k]]  ]);
            ${trace.token}                  }
            ${trace.token}                  else {
            ${trace.token}                    RF_nativePrint("\n %10d %10d %10d %12.4f %12.4f --> RGHT ", k, nonMissMembrIndx[k],  repMembrIndx[nonMissMembrIndx[k]], observation[ repMembrIndx[nonMissMembrIndx[k]] ], RF_response[treeID][1][  repMembrIndx[nonMissMembrIndx[k]]  ]);
            ${trace.token}                  }
            ${trace.token}                }
            ${trace.token}              }
            ${trace.token}              else {
            ${trace.token}                for (k = 1; k <= nonMissMembrSize; k++) {
            ${trace.token}                  if (localSplitIndicator[   nonMissMembrIndx[indxx[k]]  ] == LEFT) {
            ${trace.token}                    RF_nativePrint("\n %10d %10d %10d %12.4f %12.4f --> LEFT ", k, nonMissMembrIndx[indxx[k]],  repMembrIndx[nonMissMembrIndx[indxx[k]]], observation[ repMembrIndx[nonMissMembrIndx[indxx[k]]] ], RF_response[treeID][1][  repMembrIndx[nonMissMembrIndx[indxx[k]]]  ]);
            ${trace.token}                  }
            ${trace.token}                  else {
            ${trace.token}                    RF_nativePrint("\n %10d %10d %10d %12.4f %12.4f --> RGHT ", k, nonMissMembrIndx[indxx[k]],  repMembrIndx[nonMissMembrIndx[indxx[k]]], observation[ repMembrIndx[nonMissMembrIndx[indxx[k]]] ], RF_response[treeID][1][  repMembrIndx[nonMissMembrIndx[indxx[k]]]  ]);
            ${trace.token}                  }
            ${trace.token}                }
            ${trace.token}              }
            ${trace.token}            }
            ${trace.token}          }
          }
          else {
            delta = RF_nativeNaN;
          }
          if (!RF_nativeIsNaN(delta)) {
            if(RF_nativeIsNaN(deltaMax)) {
              deltaMax = delta;
              indexMax = j;
              ${trace.token}    if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
              ${trace.token}      RF_nativePrint("\n\nVirtual Split Statistics Updated: \n");
              ${trace.token}      RF_nativePrint("  SplitParm  SplitValIdx        Delta \n");
              ${trace.token}      RF_nativePrint(" %10d %12d %12.4f \n", covariate, indexMax, deltaMax);
              ${trace.token}    }
            }
            else {
              if ((delta - deltaMax) > EPSILON) {
                deltaMax = delta;
                indexMax = j;
                ${trace.token}    if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
                ${trace.token}      RF_nativePrint("\n\nVirtual Split Statistics Updated: \n");
                ${trace.token}      RF_nativePrint("  SplitParm  SplitValIdx        Delta \n");
                ${trace.token}      RF_nativePrint(" %10d %12d %12.4f \n", covariate, indexMax, deltaMax);
                ${trace.token}    }
              }
              else {
                ${trace.token}    if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
                ${trace.token}      RF_nativePrint("\n\nVirtual Split Statistics Not Updated: \n");
                ${trace.token}      RF_nativePrint("  SplitParm  SplitValIdx        Delta \n");
                ${trace.token}      RF_nativePrint(" %10d %12d %12.4f \n", covariate, j, delta);
                ${trace.token}    }
              }
            }
          }
          else {
            ${trace.token}    if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
            ${trace.token}      RF_nativePrint("\n\nVirtual Split Statistics Not Updated: \n");
            ${trace.token}      RF_nativePrint("  SplitParm  SplitValIdx        Delta \n");
            ${trace.token}      RF_nativePrint(" %10d %12d %12.4f \n", covariate, j, delta);
            ${trace.token}    }
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
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nregressionXwghtSplitCur(%10d) result:  %10d", treeID, result);
  ${trace.token}    RF_nativePrint("\nregressionXwghtSplitCur(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
  return result;
}
