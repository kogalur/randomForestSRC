
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "splitMahalanobis.h"
#include "splitUtil.h"
#include "regression.h"
#include "svdUtil.h"
#include "nrutil.h"
${trace.token} #include "error.h"
char mahalanobis (uint       treeID,
                  Node      *parent,
                  SplitInfoMax *splitInfoMax,
                  GreedyObj    *greedyMembr,
                  char       multImpFlag) {
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nmahalanobis(%10d) ENTRY ...\n", treeID);
  ${trace.token}  }
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
  double deltaPartial;
  double partialLeft;
  double partialRght;
  uint i, j, r, r1, r2;
  localSplitIndicator    = NULL;  
  splitVector            = NULL;  
  splitVectorSize        = 0;     
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        multImpFlag,
                                        TRUE);
  if (preliminaryResult) {
    uint  repMembrSize = parent -> repMembrSize;
    uint *repMembrIndx = parent -> repMembrIndx;
    char   *impurity   = cvector(1, RF_ySize);
    double *mean       = dvector(1, RF_ySize);
    double *variance   = dvector(1, RF_ySize);
    char impuritySummary;
    impuritySummary = FALSE;
    for (r = 1; r <= RF_ySize; r++)  {
      impurity[r] = getVarianceDoublePass(repMembrSize,
                                          repMembrIndx,
                                          0,
                                          NULL,
                                          RF_response[treeID][r],
                                          &mean[r],
                                          &variance[r]);
      impuritySummary = impuritySummary | impurity[r];
    }
    if (impuritySummary) {
      stackSplitPreliminary(repMembrSize,
                            & localSplitIndicator,
                            & splitVector);
      DistributionObj *distributionObj = stackRandomCovariates(treeID, parent);
      uint *impureIdx;
      uint  impureIdxCount;
      double **elStar;
      double **elStarTranspose;
      double **qStar;
      double **qStarPlus;
      double **u, *w, **v;      
      double *leftMean;
      double *rghtMean;
      double **leftCentered;
      double **rghtCentered;
      double **leftCenteredT;
      double **rghtCenteredT;
      double **tempResult;
      double **tempResult2;
      double *sumLeft         = dvector(1, RF_ySize);
      double *sumRght         = dvector(1, RF_ySize);
      impureIdx = uivector(1, RF_ySize);
      impureIdxCount = 0;
      for (r = 1; r <= RF_ySize; r++) {
        if (impurity[r]) {
          impureIdx[++impureIdxCount] = r;
        }  
      }  
      leftMean      = dvector(1, impureIdxCount);
      rghtMean      = dvector(1, impureIdxCount);
      leftCentered  = dmatrix(1, impureIdxCount, 1, 1);
      rghtCentered  = dmatrix(1, impureIdxCount, 1, 1);
      leftCenteredT = dmatrix(1, 1, 1, impureIdxCount);
      rghtCenteredT = dmatrix(1, 1, 1, impureIdxCount);
      elStar = dmatrix(1, repMembrSize, 1, impureIdxCount);
      ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
      ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
      ${trace.token}              RF_nativePrint("\nMean of responses:       index   response");
      ${trace.token}              for (uint rr = 1; rr <= impureIdxCount; rr++) {
      ${trace.token}                RF_nativePrint("\n                    %10d %10.4f", impureIdx[rr], mean[impureIdx[rr]]);
      ${trace.token}              }
      ${trace.token}            }
      ${trace.token}          }
      for (i = 1; i <= repMembrSize; i++) {
        for (uint rr = 1; rr <= impureIdxCount; rr++) {
          r = impureIdx[rr];
          elStar[i][rr] = RF_response[treeID][r][ repMembrIndx[i] ] - mean[r];
        }
      }
      qStarPlus = NULL;
      if (RF_qStarPlus != NULL) {
        qStarPlus = dmatrix(1, impureIdxCount, 1, impureIdxCount);
        for (uint rr1 = 1; rr1 <= impureIdxCount; rr1++) {
          r1 = impureIdx[rr1];
          for (uint rr2 = 1; rr2 <= impureIdxCount; rr2++) {
            r2 = impureIdx[rr2];
            qStarPlus[rr1][rr2] = RF_qStarPlus[r1][r2];
          }
        }
        elStarTranspose = NULL;
        qStar = NULL;
      }
      else {
        ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
        ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
        ${trace.token}              RF_nativePrint("\n elStar (of dim n x r'):  ");
        ${trace.token}              matrixPrint(elStar, repMembrSize, impureIdxCount);
        ${trace.token}            }
        ${trace.token}          }
        elStarTranspose = matrixTrans(elStar, repMembrSize, impureIdxCount);
        ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
        ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
        ${trace.token}              RF_nativePrint("\n elStarTranspose (of dim r' x n):  ");
        ${trace.token}              matrixPrint(elStarTranspose, impureIdxCount, repMembrSize);
        ${trace.token}            }
        ${trace.token}          }
        qStar = matrixMult(elStarTranspose, elStar, impureIdxCount, repMembrSize, impureIdxCount);
        ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
        ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
        ${trace.token}              RF_nativePrint("\n qStar (of dim r' x r'):  ");
        ${trace.token}              matrixPrint(qStar, impureIdxCount, impureIdxCount);
        ${trace.token}            }
        ${trace.token}          }
        svdcmp(qStar, impureIdxCount, impureIdxCount, &u, &w, &v);
        ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
        ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
        ${trace.token}              RF_nativePrint("\n SVD Check of qStar (qStar =? u x w x v):  ");
        ${trace.token}              svdchk(qStar, impureIdxCount, impureIdxCount, u, w, v);
        ${trace.token}            }
        ${trace.token}          }
        qStarPlus = svdinv(u, w, v, impureIdxCount, impureIdxCount, impureIdxCount);
        ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
        ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
        ${trace.token}            if (qStarPlus != NULL) {
        ${trace.token}              RF_nativePrint("\n qStar x qStarPlus (of dim r' x r'):  ");
        ${trace.token}              double **qStarChk = matrixMult(qStar, qStarPlus, impureIdxCount, impureIdxCount, impureIdxCount);
        ${trace.token}              matrixPrint(qStarChk, impureIdxCount, impureIdxCount);
        ${trace.token}              RF_nativePrint("\n qStarPlus x qStar (of dim r' x r'):  ");
        ${trace.token}              double **qStarChk2 = matrixMult(qStarPlus, qStar, impureIdxCount, impureIdxCount, impureIdxCount);
        ${trace.token}              matrixPrint(qStarChk2, impureIdxCount, impureIdxCount);
        ${trace.token}              free_dmatrix(qStarChk, 1, impureIdxCount, 1, impureIdxCount);
        ${trace.token}              free_dmatrix(qStarChk2, 1, impureIdxCount, 1, impureIdxCount);
        ${trace.token}            }
        ${trace.token}            }
        ${trace.token}          }
      }
      ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
      ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
      ${trace.token}            if (qStarPlus != NULL) {
      ${trace.token}              RF_nativePrint("\n qStarPlus (of dim r' x r'):  ");
      ${trace.token}              matrixPrint(qStarPlus, impureIdxCount, impureIdxCount);
      ${trace.token}            }
      ${trace.token}            }
      ${trace.token}          }
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
          observation = RF_observation[treeID][covariate];
          leftSize = 0;
          priorMembrIter = 0;
          if (factorFlag == FALSE) {
            for (j = 1; j <= repMembrSize; j++) {
              localSplitIndicator[j] = RIGHT;
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
            rghtSize = repMembrSize - leftSize;
            ${trace.token}   if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
            ${trace.token}     RF_nativePrint("\nNon-miss Node Size:  %10d, Non-miss Left Size:  %10d, Non-miss Right Size:  %10d", repMembrSize, leftSize, rghtSize);
            ${trace.token}   }
            ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
            ${trace.token}            if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
            ${trace.token}              RF_nativePrint("\n PriorIter:     %10d  CurrentIter:   %10d", priorMembrIter, currentMembrIter);
            ${trace.token}            }
            ${trace.token}          }
            if ((leftSize != 0) && (rghtSize != 0)) {
              delta        = 0.0;
              deltaPartial = 0.0;
              if (qStarPlus != NULL) {
                for (uint rr = 1; rr <= impureIdxCount; rr++) {
                  leftMean[rr] = rghtMean[rr] = 0.0;
                }
                for (i = 1; i <= repMembrSize; i++) {
                  if (localSplitIndicator[i] == LEFT) {
                    for (uint rr = 1; rr <= impureIdxCount; rr++) {
                      r = impureIdx[rr];
                      leftMean[rr] += elStar[i][rr];
                    }
                  }
                  else {
                    for (uint rr = 1; rr <= impureIdxCount; rr++) {
                      r = impureIdx[rr];
                      rghtMean[rr] += elStar[i][rr];
                    }
                  }
                }
                ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
                ${trace.token}              RF_nativePrint("\nLeft and Right Response Means:       index   response       left      right");
                ${trace.token}            }
                ${trace.token}          }
                for (uint rr = 1; rr <= impureIdxCount; rr++) {
                  leftMean[rr] = leftMean[rr] / leftSize;
                  rghtMean[rr] = rghtMean[rr] / rghtSize;
                }
                ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
                ${trace.token}              for (uint rr = 1; rr <= impureIdxCount; rr++) {
                ${trace.token}                RF_nativePrint("\n                                %10d %10d %10.4f %10.4f",
                ${trace.token}                rr,
                ${trace.token}                impureIdx[rr],
                ${trace.token}                leftMean[rr],
                ${trace.token}                rghtMean[rr]);
                ${trace.token}              }
                ${trace.token}            }
                ${trace.token}          }
                partialLeft = partialRght = 0.0;
                for (i = 1; i <= repMembrSize; i++) {
                  if (localSplitIndicator[i] == LEFT) {
                    for (uint rr = 1; rr <= impureIdxCount; rr++) {
                      r = impureIdx[rr];
                      leftCentered[rr][1]  = elStar[i][rr] - leftMean[rr];  
                      leftCenteredT[1][rr] = leftCentered[rr][1];
                    }
                    tempResult = matrixMult(leftCenteredT, qStarPlus, 1, impureIdxCount, impureIdxCount);
                    tempResult2 = matrixMult(tempResult, leftCentered, 1, impureIdxCount, 1);
                    partialLeft += tempResult2[1][1];
                  }
                  else {
                    for (uint rr = 1; rr <= impureIdxCount; rr++) {
                      r = impureIdx[rr];
                      rghtCentered[rr][1]  = elStar[i][rr] - rghtMean[rr];  
                      rghtCenteredT[1][rr] = rghtCentered[rr][1];
                    }
                    tempResult = matrixMult(rghtCenteredT, qStarPlus, 1, impureIdxCount, impureIdxCount);
                    tempResult2 = matrixMult(tempResult, rghtCentered, 1, impureIdxCount, 1);
                    partialRght += tempResult2[1][1];
                  }
                  free_dmatrix(tempResult,  1, 1, 1, impureIdxCount);
                  free_dmatrix(tempResult2, 1, 1, 1, 1);
                }  
                deltaPartial = ( ((double) leftSize / repMembrSize) * partialLeft) + ( ((double) rghtSize / repMembrSize) * partialRght);
                ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
                ${trace.token}              RF_nativePrint("\nMahalanobis Internals:                left                right                  sum");
                ${trace.token}              RF_nativePrint("\n                      %20.8f %20.8f %20.8f", ( ((double) leftSize / repMembrSize) * partialLeft), ( ((double) rghtSize / repMembrSize) * partialRght), deltaPartial);
                ${trace.token}            }
                ${trace.token}          }
                delta = 1.0 - (deltaPartial / impureIdxCount);
              }
              else {
                for (uint rr = 1; rr <= impureIdxCount; rr++) {
                  r = impureIdx[rr];
                  sumLeft[r] = sumRght[r] = 0.0;
                }
                for (i = 1; i <= repMembrSize; i++) {
                  if (localSplitIndicator[i] == LEFT) {
                    for (uint rr = 1; rr <= impureIdxCount; rr++) {
                      r = impureIdx[rr];
                      sumLeft[r] += RF_response[treeID][r][ repMembrIndx[i] ] - mean[r];
                    }
                  }
                  else {
                    for (uint rr = 1; rr <= impureIdxCount; rr++) {
                      r = impureIdx[rr];
                      sumRght[r] += RF_response[treeID][r][ repMembrIndx[i] ] - mean[r];
                    }
                  }
                }
                for (uint rr = 1; rr <= impureIdxCount; rr++) {
                  r = impureIdx[rr];
                  partialLeft = pow (sumLeft[r], 2.0) / (leftSize * variance[r]);
                  partialRght = pow (sumRght[r], 2.0) / (rghtSize * variance[r]);
                  ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                  ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
                  ${trace.token}              RF_nativePrint("\nPartial Delta (Regr): leftSize   rghtSize       mean   variance     sumLFT     sumRGT partialLFT partialRGT        sum");
                  ${trace.token}              RF_nativePrint("\n                    %10d %10d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f",
                  ${trace.token}                       leftSize,
                  ${trace.token}                       rghtSize,
                  ${trace.token}                       mean[r],
                  ${trace.token}                       variance[r],
                  ${trace.token}                       sumLeft[r],
                  ${trace.token}                       sumRght[r],
                  ${trace.token}                       partialLeft,
                  ${trace.token}                       partialRght,
                  ${trace.token}                       partialLeft+partialRght);
                  ${trace.token}            }
                  ${trace.token}          }
                  delta = partialLeft + partialRght;
                }
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
            ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
            ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
            ${trace.token}              RF_nativePrint("\nVirtual (non-miss) Membership:  ");
            ${trace.token}              for (uint k = 1; k <= repMembrSize; k++) {
            ${trace.token}                if (localSplitIndicator[ indxx[k] ] == LEFT) {
            ${trace.token}                  RF_nativePrint("\n %10d %10d %10d %12.4f %12.4f --> LEFT ", k, indxx[k], repMembrIndx[indxx[k]], observation[ repMembrIndx[indxx[k]] ], RF_response[treeID][1][  repMembrIndx[indxx[k]]  ]);
            ${trace.token}                }
            ${trace.token}                else {
            ${trace.token}                  RF_nativePrint("\n %10d %10d %10d %12.4f %12.4f --> RGHT ", k, indxx[k], repMembrIndx[indxx[k]], observation[ repMembrIndx[indxx[k]] ], RF_response[treeID][1][  repMembrIndx[indxx[k]]  ]);
            ${trace.token}                }
            ${trace.token}              }
            ${trace.token}            }
            ${trace.token}          }
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
      if (qStarPlus != NULL)  {
        free_dmatrix(qStarPlus, 1, impureIdxCount, 1, impureIdxCount);
      }
      if (RF_qStarPlus == NULL) {
        free_svdcmp(qStar, impureIdxCount, impureIdxCount, u, w, v);
        free_dmatrix(elStarTranspose, 1, impureIdxCount, 1, repMembrSize);
      }
      free_dmatrix(elStar, 1, repMembrSize, 1, impureIdxCount);      
      free_dvector(leftMean, 1, impureIdxCount);
      free_dvector(rghtMean, 1, impureIdxCount);
      free_dmatrix(leftCentered, 1, impureIdxCount, 1, 1);
      free_dmatrix(rghtCentered, 1, impureIdxCount, 1, 1);
      free_dmatrix(leftCenteredT, 1, 1, 1, impureIdxCount);
      free_dmatrix(rghtCenteredT, 1, 1, 1, impureIdxCount);
      free_uivector(impureIdx, 1, RF_ySize);
      free_dvector(sumLeft,     1, RF_ySize);
      free_dvector(sumRght,     1, RF_ySize);
      unstackRandomCovariates(treeID, distributionObj);
      unstackSplitPreliminary(repMembrSize,
                              localSplitIndicator,
                              splitVector);
    }  
    free_cvector(impurity,   1, RF_ySize);
    free_dvector(mean,       1, RF_ySize);
    free_dvector(variance,   1, RF_ySize);
  }  
  unstackPreSplit(preliminaryResult,
                  parent,
                  multImpFlag,
                  TRUE);  
  result = summarizeSplitResult(splitInfoMax);
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nmahalanobis(%10d) result:  %10d", treeID, result);
  ${trace.token}    RF_nativePrint("\nmahalanobis(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
  return result;
}
