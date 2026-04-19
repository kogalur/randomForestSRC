
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "splitCustomDriver.h"
#include "splitUtil.h"
#include "splitUtilSurv.h"
#include "nrutil.h"
${trace.token} #include "error.h"
char customMultivariateSplit (uint       treeID,
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
  double delta, deltaPartial;
  uint   deltaNorm;
  uint j, k, m, r, rr;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\ncustomMultivariateSplit() ENTRY ...\n");
  ${trace.token}  }
  localSplitIndicator    = NULL;  
  splitVector            = NULL;  
  splitVectorSize        = 0;     
  mwcpSizeAbsolute       = 0;     
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        multImpFlag,
                                        TRUE);
  if (preliminaryResult) {
    uint  repMembrSize = parent -> repMembrSize;
    uint *repMembrIndx = parent -> repMembrIndx;
    uint  nonMissMembrSize;
    uint *nonMissMembrIndx;
    char   *impurity   = cvector(1, RF_ySize);
    double *mean       = dvector(1, RF_ySize);
    double *variance   = dvector(1, RF_ySize);
    char impuritySummary;
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      impuritySummary = FALSE;
      for (r = 1; r <= RF_ySize; r++)  {
        if (RF_yWeight[r] == 0) {
          impurity[r] = FALSE;
        }
        else {
          impurity[r] = getVariance(repMembrSize,
                                    repMembrIndx,
                                    0,
                                    NULL,
                                    RF_response[treeID][r],
                                    &mean[r],
                                    &variance[r]);
        }
        impuritySummary = impuritySummary | impurity[r];
      }
    }
    else {
      impuritySummary = TRUE;
    }
    if (impuritySummary) {
      stackSplitPreliminary(repMembrSize,
                            & localSplitIndicator,
                            & splitVector);
      DistributionObj *distributionObj = stackRandomCovariates(treeID, parent);
      char **secondNonMissMembrFlag = (char **) new_vvector(1, RF_ySize, NRUTIL_CPTR);
      uint  *secondNonMissMembrSize =           uivector(1, RF_ySize);
      uint  *secondNonMissMembrLeftSize =       uivector(1, RF_ySize);
      uint  *secondNonMissMembrRghtSize =       uivector(1, RF_ySize);
      char  *tempNonMissMembrFlag = 0;
      uint  *tempNonMissMembrIndx;
      char   mResponseFlag;
      char   nonMissImpuritySummary;
      covariateCount = 0;
      double deltaMax;
      uint   indexMax;
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
          if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
            tempNonMissMembrFlag = cvector(1, nonMissMembrSize);
            for (k = 1; k <= nonMissMembrSize; k++) {
              tempNonMissMembrFlag[k] = TRUE;
            }
            for (r = 1; r <= RF_ySize; r++) {
              secondNonMissMembrFlag[r] = tempNonMissMembrFlag;
              secondNonMissMembrSize[r] = nonMissMembrSize;
            }
            nonMissImpuritySummary = TRUE;
          }
          else {
            tempNonMissMembrIndx = uivector(1, nonMissMembrSize);
            nonMissImpuritySummary = FALSE;
            for (r = 1; r <= RF_ySize; r++)  {
              secondNonMissMembrFlag[r] = cvector(1, nonMissMembrSize);
              j = 0;
              for (k = 1; k <= nonMissMembrSize; k++) {
                mResponseFlag = FALSE;
                if (RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[k]] ] > 0) {
                  if (RF_mpSign[r][RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[k]] ]] == 1) {
                    mResponseFlag = TRUE;
                  }
                }
                if (!mResponseFlag) {
                  j ++;
                  tempNonMissMembrIndx[j] = nonMissMembrIndx[k];
                  secondNonMissMembrFlag[r][k] = TRUE;
                }
                else {
                  secondNonMissMembrFlag[r][k] = FALSE;
                }
              }  
              secondNonMissMembrSize[r] = j;
              if (RF_yWeight[r] == 0) {
                impurity[r] = FALSE;
              }
              else {
                impurity[r] = getVariance(repMembrSize,
                                          repMembrIndx,
                                          secondNonMissMembrSize[r],
                                          tempNonMissMembrIndx,
                                          RF_response[treeID][r],
                                          &mean[r],
                                          &variance[r]);
              }
              nonMissImpuritySummary = nonMissImpuritySummary | impurity[r];
              secondNonMissMembrLeftSize[r] = secondNonMissMembrRghtSize[r] = 0;
            }  
            ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
            ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
            ${trace.token}              RF_nativePrint("\nImpurity:");
            ${trace.token}              for (r=1; r <= RF_ySize; r++) {
            ${trace.token}                RF_nativePrint(" %10d", r);
            ${trace.token}              }
            ${trace.token}              RF_nativePrint("\n          ");
            ${trace.token}              for (r=1; r <= RF_ySize; r++) {
            ${trace.token}                RF_nativePrint( "%10d,", impurity[r]);
            ${trace.token}              }
            ${trace.token}            }
            ${trace.token}          }
            free_uivector(tempNonMissMembrIndx, 1, nonMissMembrSize);
          }  
          if (nonMissImpuritySummary) {
            leftSize = 0;
            priorMembrIter = 0;
            if (factorFlag == FALSE) {
              for (j = 1; j <= nonMissMembrSize; j++) {
                localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
              }
              for (r = 1; r <= RF_ySize; r++) {
                if (impurity[r]) {
                  secondNonMissMembrLeftSize[r] = 0;
                  secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r];
                }
              }
            }
            double *userResponse = dvector(1, nonMissMembrSize);
            char   *userSplitIndicator = cvector(1, nonMissMembrSize);
            double **userFeature = NULL;
            if (RF_yIndexZeroSize > 0) {
              userFeature = dmatrix(1, RF_yIndexZeroSize, 1, nonMissMembrSize);
            }
            deltaMax =  RF_nativeNaN;
            indexMax =  0;
            for (j = 1; j < splitLength; j++) {
              if (factorFlag == TRUE) {
                priorMembrIter = 0;
                leftSize = 0;
                for (r = 1; r <= RF_ySize; r++) {
                  secondNonMissMembrLeftSize[r] = 0;
                  secondNonMissMembrRghtSize[r] = 0;
                }
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
              ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
              ${trace.token}            if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
              ${trace.token}              RF_nativePrint("\n PriorIter:     %10d  CurrentIter:   %10d", priorMembrIter, currentMembrIter);
              ${trace.token}            }
              ${trace.token}          }
              if ((leftSize != 0) && (rghtSize != 0)) {
                delta        = 0.0;
                deltaPartial = 0.0;
                deltaNorm    = 0;
                for (r = 1; r <= RF_ySize; r++) {
                  if (impurity[r]) {
                    if (factorFlag == TRUE) {
                      for (k = 1; k <= nonMissMembrSize; k++) {
                        if (secondNonMissMembrFlag[r][indxx[k]] == TRUE) {
                          if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                            secondNonMissMembrLeftSize[r] ++;
                          }
                          else {
                          }
                        }
                      }
                      secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r] - secondNonMissMembrLeftSize[r];
                    }
                    else {
                      for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                        if (secondNonMissMembrFlag[r][indxx[k]] == TRUE) {
                          secondNonMissMembrLeftSize[r] ++;
                          secondNonMissMembrRghtSize[r] --;
                        }
                      }
                    }  
                    if ((secondNonMissMembrLeftSize[r] > 0) && (secondNonMissMembrRghtSize[r] > 0)) {
                      m = 0;
                      ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                      ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
                      ${trace.token}               RF_nativePrint("\n          k     indexx     nonMissMembrIndx         repMembrIndx  2ndNonMissMembrFlag             response");
                      ${trace.token}            }
                      ${trace.token}          }
                      for (k = 1; k <= nonMissMembrSize; k++) {
                        ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                        ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
                        ${trace.token}              RF_nativePrint("\n %10d %10d %20d %20d %20d %20.4f", k, indxx[k], nonMissMembrIndx[indxx[k]], repMembrIndx[nonMissMembrIndx[indxx[k]]], secondNonMissMembrFlag[r][k], RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]);
                        ${trace.token}            }
                        ${trace.token}          }
                        if (secondNonMissMembrFlag[r][indxx[k]] == TRUE) {
                          userResponse[++m] = RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                          userSplitIndicator[m] = localSplitIndicator[ nonMissMembrIndx[indxx[k]] ];
                          for (rr = 1; rr <= RF_yIndexZeroSize; rr++) {
                            userFeature[rr][m] = RF_response[treeID][RF_yIndexZero[rr]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                          }
                        }
                      }
                      if ((RF_rType[r] == 'B') ||
                          (RF_rType[r] == 'I') ||
                          (RF_rType[r] == 'C')) {
                        ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                        ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
                        ${trace.token}               RF_nativePrint("\nGetting custom family %10d at %10d with %20x", CLAS_FAM, RF_splitCustomIdx, customFunctionArray[CLAS_FAM][RF_splitCustomIdx]);
                        ${trace.token}            }
                        ${trace.token}          }
                        deltaPartial = customFunctionArray[CLAS_FAM][RF_splitCustomIdx](m,
                                                                                        userSplitIndicator,
                                                                                        NULL,
                                                                                        NULL,
                                                                                        0,
                                                                                        0,
                                                                                        NULL,
                                                                                        userResponse,
                                                                                        mean[r],
                                                                                        variance[r],
                                                                                        RF_rFactorSize[RF_rFactorMap[r]],
                                                                                        userFeature,
                                                                                        RF_yIndexZeroSize);
                        ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                        ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
                        ${trace.token}              RF_nativePrint("\nPartial Delta (Clas):  2ndLeft    2ndRght       mean   variance  rFactSize      delta");
                        ${trace.token}              RF_nativePrint("\n                    %10d %10d %10.4f %10.4f %10d %10.4f",
                        ${trace.token}                       secondNonMissMembrLeftSize[r],
                        ${trace.token}                       secondNonMissMembrRghtSize[r],
                        ${trace.token}                       mean[r],
                        ${trace.token}                       variance[r],
                        ${trace.token}                       RF_rFactorSize[RF_rFactorMap[r]],
                        ${trace.token}                       deltaPartial);
                        ${trace.token}            }
                        ${trace.token}          }
                      }
                      else {
                        ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                        ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
                        ${trace.token}               RF_nativePrint("\nGetting custom family %10d at %10d with %20x", REGR_FAM, RF_splitCustomIdx, customFunctionArray[REGR_FAM][RF_splitCustomIdx]);
                        ${trace.token}            }
                        ${trace.token}          }
                        deltaPartial = (customFunctionArray[REGR_FAM][RF_splitCustomIdx])(m,
                                                                                          userSplitIndicator,
                                                                                          NULL,
                                                                                          NULL,
                                                                                          0,
                                                                                          0,
                                                                                          NULL,
                                                                                          userResponse,
                                                                                          mean[r],
                                                                                          variance[r],
                                                                                          0,
                                                                                          userFeature,
                                                                                          RF_yIndexZeroSize);
                        ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                        ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
                        ${trace.token}              RF_nativePrint("\nPartial Delta (Regr):  2ndLeft    2ndRght       mean   variance  rFactSize      delta");
                        ${trace.token}              RF_nativePrint("\n                    %10d %10d %10.4f %10.4f %10d %10.4f",
                        ${trace.token}                       secondNonMissMembrLeftSize[r],
                        ${trace.token}                       secondNonMissMembrRghtSize[r],
                        ${trace.token}                       mean[r],
                        ${trace.token}                       variance[r],
                        ${trace.token}                       0,
                        ${trace.token}                       deltaPartial);
                        ${trace.token}            }
                        ${trace.token}          }
                      }
                      if (!RF_nativeIsNaN(deltaPartial)) {
                        deltaNorm ++;
                        delta += deltaPartial;
                      }
                    }
                  }  
                }  
                ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
                ${trace.token}              RF_nativePrint("\nDelta Normalization:   %10d", deltaNorm);
                ${trace.token}            }
                ${trace.token}          }
                if (deltaNorm > 0) {
                  delta = delta / (double) deltaNorm;
                }
                else {
                  delta = RF_nativeNaN;
                }
              }
              else {
                delta = RF_nativeNaN;
              }
              ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
              ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
              ${trace.token}              RF_nativePrint("\nVirtual (non-miss) Membership:  ");
              ${trace.token}              for (k = 1; k <= nonMissMembrSize; k++) {
              ${trace.token}                if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
              ${trace.token}                  RF_nativePrint("\n %10d %10d %12.4f %12.4f %12.4f --> LEFT ", k, nonMissMembrIndx[indxx[k]], repMembrIndx[nonMissMembrIndx[indxx[k]]], RF_observation[treeID][covariate][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ], RF_response[treeID][1][  repMembrIndx[nonMissMembrIndx[indxx[k]]]  ]);
              ${trace.token}                }
              ${trace.token}                else {
              ${trace.token}                  RF_nativePrint("\n %10d %10d %12.4f %12.4f %12.4f --> RGHT ", k, nonMissMembrIndx[indxx[k]], repMembrIndx[nonMissMembrIndx[indxx[k]]], RF_observation[treeID][covariate][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ], RF_response[treeID][1][  repMembrIndx[nonMissMembrIndx[indxx[k]]]  ]);
              ${trace.token}                }
              ${trace.token}              }
              ${trace.token}            }
              ${trace.token}          }
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
            if (RF_yIndexZeroSize > 0) {
              free_dmatrix (userFeature, 1, RF_yIndexZeroSize, 1, nonMissMembrSize);
            }
            free_dvector (userResponse, 1, nonMissMembrSize);
            free_cvector (userSplitIndicator, 1, nonMissMembrSize);
          }  
          if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
            free_cvector(tempNonMissMembrFlag, 1, nonMissMembrSize);
          }
          else {
            for (r = 1; r <= RF_ySize; r++)  {
              free_cvector(secondNonMissMembrFlag[r], 1, nonMissMembrSize);
            }
          }
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
      free_new_vvector(secondNonMissMembrFlag,  1, RF_ySize, NRUTIL_CPTR);
      free_uivector(secondNonMissMembrSize,     1, RF_ySize);
      free_uivector(secondNonMissMembrLeftSize, 1, RF_ySize);
      free_uivector(secondNonMissMembrRghtSize, 1, RF_ySize);
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
  ${trace.token}    RF_nativePrint("\ncustomMultivariateSplit(%10d) result:  %10d", treeID, result);
  ${trace.token}    RF_nativePrint("\ncustomMultivariateSplit(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
  return result;
}
char customSurvivalSplit (uint       treeID,
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
  uint j, k, m, rr;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\ncustomSplitSurvival() ENTRY ...\n");
  ${trace.token}  }
  localSplitIndicator    = NULL;  
  splitVector            = NULL;  
  splitVectorSize        = 0;     
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
    uint *localEventTimeCount, *localEventTimeIndex;
    uint  localEventTimeSize;
    uint *nodeParentEvent,  *nodeLeftEvent,  *nodeRightEvent;
    uint *nodeParentAtRisk, *nodeLeftAtRisk, *nodeRightAtRisk;
    localEventTimeSize = 0;  
    delta = 0.0;             
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      stackAndGetSplitSurv(treeID,
                           parent,
                           TRUE,
                           & localEventTimeCount,
                           & localEventTimeIndex,
                           & localEventTimeSize,
                           & nodeParentEvent,
                           & nodeParentAtRisk,
                           & nodeLeftEvent,
                           & nodeLeftAtRisk,
                           & nodeRightEvent,
                           & nodeRightAtRisk);
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
        if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
          stackAndGetSplitSurv(treeID,
                               parent,
                               TRUE,
                               & localEventTimeCount,
                               & localEventTimeIndex,
                               & localEventTimeSize,
                               & nodeParentEvent,
                               & nodeParentAtRisk,
                               & nodeLeftEvent,
                               & nodeLeftAtRisk,
                               & nodeRightEvent,
                               & nodeRightAtRisk);
        }
        if (localEventTimeSize > 0) {
          leftSize = 0;
          priorMembrIter = 0;
          if (factorFlag == FALSE) {
            for (j = 1; j <= nonMissMembrSize; j++) {
              localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
            }
          }
          double *userTime  = dvector(1, nonMissMembrSize);
          double *userEvent = dvector(1, nonMissMembrSize);
          double *userEventTime = dvector(1, localEventTimeSize);
          char   *userSplitIndicator = cvector(1, nonMissMembrSize);
          uint   *userSort  = uivector(1, nonMissMembrSize);
          double *tempTime  = dvector(1, nonMissMembrSize);
          double **userFeature = NULL;
          if (RF_yIndexZeroSize > 0) {
            userFeature = dmatrix(1, RF_yIndexZeroSize, 1, nonMissMembrSize);
          }
          for (k = 1; k <= nonMissMembrSize; k++) {
            tempTime[k]  = RF_time[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
          }
          indexx(nonMissMembrSize, tempTime, userSort);
          for (k = 1; k <= nonMissMembrSize; k++) {
            userTime[k]  = RF_time[treeID][ repMembrIndx[nonMissMembrIndx[indxx[userSort[k]]]] ];
            userEvent[k] = RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[userSort[k]]]] ];
            for (rr = 1; rr <= RF_yIndexZeroSize; rr++) {
              userFeature[rr][k] = RF_response[treeID][RF_yIndexZero[rr]][ repMembrIndx[nonMissMembrIndx[indxx[userSort[k]]]] ];
            }
          }
          for (m = 1; m <= localEventTimeSize; m++) {
            userEventTime[m] = RF_masterTime[localEventTimeIndex[m]];
          }
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
            ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
            ${trace.token}            if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
            ${trace.token}              RF_nativePrint("\n PriorIter:     %10d  CurrentIter:   %10d", priorMembrIter, currentMembrIter);
            ${trace.token}            }
            ${trace.token}          }
            for (k = 1; k <= nonMissMembrSize; k++) {
              userSplitIndicator[k] = localSplitIndicator[ nonMissMembrIndx[indxx[userSort[k]]] ];
            }
            ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
            ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
            ${trace.token}               RF_nativePrint("\nGetting custom family %10d at %10d with %20x", SURV_FAM, RF_splitCustomIdx, customFunctionArray[SURV_FAM][RF_splitCustomIdx]);
            ${trace.token}            }
            ${trace.token}          }
            ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
            ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
            ${trace.token}      RF_nativePrint("\n Outgoing custom split event type count (unused): %10d", RF_eventType[RF_eventTypeSize]);
            ${trace.token}      RF_nativePrint("\n Outgoing custom split event times:");
            ${trace.token}      for (k = 1; k <= localEventTimeSize; k++) {
            ${trace.token}        RF_nativePrint("\n %10d %10.4f", k, userEventTime[k]);
            ${trace.token}      }
            ${trace.token}      RF_nativePrint("\n");
            ${trace.token}      RF_nativePrint("\n Outgoing custom split membership:");
            ${trace.token}      RF_nativePrint("\n      index  membrship       time     status");
            ${trace.token}      for (k = 1; k <= nonMissMembrSize; k++) {
            ${trace.token}        RF_nativePrint("\n %10d %10d %10.4f %10.4f", k, userSplitIndicator[k], userTime[k], userEvent[k]);
            ${trace.token}      }
            ${trace.token}    }  
            ${trace.token}  }  
            if ((leftSize != 0) && (rghtSize != 0)) {
              delta = customFunctionArray[SURV_FAM][RF_splitCustomIdx](nonMissMembrSize,
                                                                       userSplitIndicator,
                                                                       userTime,
                                                                       userEvent,
                                                                       0,
                                                                       localEventTimeSize,
                                                                       userEventTime,
                                                                       NULL,
                                                                       0,
                                                                       0,
                                                                       0,
                                                                       userFeature,
                                                                       RF_yIndexZeroSize);
            }
            else {
              delta = RF_nativeNaN;
            }
            updateMaximumSplit(treeID,
                               parent,
                               delta,
                               covariate,
                               j,
                               factorFlag,
                               mwcpSizeAbsolute,
                               repMembrSize,
                               & localSplitIndicator,
                               splitVectorPtr,
                               splitInfoMax);
            if (factorFlag == FALSE) {
              priorMembrIter = currentMembrIter - 1;
            }
          }  
          if (RF_yIndexZeroSize > 0) {
            free_dmatrix (userFeature, 1, RF_yIndexZeroSize, 1, nonMissMembrSize);
          }
          free_uivector(userSort, 1, nonMissMembrSize);
          free_dvector(tempTime, 1, nonMissMembrSize);
          free_dvector (userTime, 1, nonMissMembrSize);
          free_dvector (userEvent, 1, nonMissMembrSize);
          free_dvector (userEventTime, 1, localEventTimeSize);
          free_cvector (userSplitIndicator, 1, nonMissMembrSize);
        }  
        else {
          ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
          ${trace.token}  if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
          ${trace.token}    RF_nativePrint("\nCovariate ignored due to zero localEventTimeSize:  %10d", covariate);
          ${trace.token}  }
          ${trace.token}  }
        }
        if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
          unstackSplitSurv(treeID,
                           parent,
                           localEventTimeCount,
                           localEventTimeIndex,
                           localEventTimeSize,
                           nodeParentEvent,
                           nodeParentAtRisk,
                           nodeLeftEvent,
                           nodeLeftAtRisk,
                           nodeRightEvent,
                           nodeRightAtRisk);
        }
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
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      unstackSplitSurv(treeID,
                       parent,
                       localEventTimeCount,
                       localEventTimeIndex,
                       localEventTimeSize,
                       nodeParentEvent,
                       nodeParentAtRisk,
                       nodeLeftEvent,
                       nodeLeftAtRisk,
                       nodeRightEvent,
                       nodeRightAtRisk);
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
  ${trace.token}    RF_nativePrint("\ncustomSplitSurvival(%10d) result:  %10d", treeID, result);
  ${trace.token}    RF_nativePrint("\ncustomSplitSurvival(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
  return result;
}
char customCompetingRiskSplit (uint       treeID,
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
  uint j, k, m, rr;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\ncustomSplitCompetingRisk() ENTRY ...\n");
  ${trace.token}  }
  localSplitIndicator    = NULL;  
  splitVector            = NULL;  
  splitVectorSize        = 0;     
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
    uint *localEventTimeCount, *localEventTimeIndex;
    uint  localEventTimeSize;
    uint *nodeParentEvent,  *nodeLeftEvent,  *nodeRightEvent;
    uint *nodeParentAtRisk, *nodeLeftAtRisk, *nodeRightAtRisk;
    localEventTimeSize = 0;  
    delta = 0.0;  
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      stackAndGetSplitSurv(treeID,
                           parent,
                           TRUE,
                           & localEventTimeCount,
                           & localEventTimeIndex,
                           & localEventTimeSize,
                           & nodeParentEvent,
                           & nodeParentAtRisk,
                           & nodeLeftEvent,
                           & nodeLeftAtRisk,
                           & nodeRightEvent,
                           & nodeRightAtRisk);
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
        if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
          stackAndGetSplitSurv(treeID,
                               parent,
                               TRUE,
                               & localEventTimeCount,
                               & localEventTimeIndex,
                               & localEventTimeSize,
                               & nodeParentEvent,
                               & nodeParentAtRisk,
                               & nodeLeftEvent,
                               & nodeLeftAtRisk,
                               & nodeRightEvent,
                               & nodeRightAtRisk);
        }
        if (localEventTimeSize > 0) {
          leftSize = 0;
          priorMembrIter = 0;
          if (factorFlag == FALSE) {
            for (j = 1; j <= nonMissMembrSize; j++) {
              localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
            }
          }
          double *userTime  = dvector(1, nonMissMembrSize);
          double *userEvent = dvector(1, nonMissMembrSize);
          double *userEventTime = dvector(1, localEventTimeSize);
          char   *userSplitIndicator = cvector(1, nonMissMembrSize);
          uint   *userSort  = uivector(1, nonMissMembrSize);
          double *tempTime  = dvector(1, nonMissMembrSize);
          double **userFeature = NULL;
          if (RF_yIndexZeroSize > 0) {
            userFeature = dmatrix(1, RF_yIndexZeroSize, 1, nonMissMembrSize);
          }
          for (k = 1; k <= nonMissMembrSize; k++) {
            tempTime[k]  = RF_time[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
          }
          indexx(nonMissMembrSize, tempTime, userSort);
          for (k = 1; k <= nonMissMembrSize; k++) {
            userTime[k]  = RF_time[treeID][ repMembrIndx[nonMissMembrIndx[indxx[userSort[k]]]] ];
            userEvent[k] = RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[userSort[k]]]] ];
            for (rr = 1; rr <= RF_yIndexZeroSize; rr++) {
              userFeature[rr][k] = RF_response[treeID][RF_yIndexZero[rr]][ repMembrIndx[nonMissMembrIndx[indxx[userSort[k]]]] ];
            }
          }
          for (m = 1; m <= localEventTimeSize; m++) {
            userEventTime[m] = RF_masterTime[localEventTimeIndex[m]];
          }
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
            ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
            ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
            ${trace.token}              RF_nativePrint("\n PriorIter:     %10d  CurrentIter:   %10d", priorMembrIter, currentMembrIter);
            ${trace.token}            }
            ${trace.token}          }
            for (k = 1; k <= nonMissMembrSize; k++) {
              userSplitIndicator[k] = localSplitIndicator[ nonMissMembrIndx[indxx[userSort[k]]] ];
            }
            ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
            ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
            ${trace.token}      RF_nativePrint("\nGetting custom family %10d at %10d with %20x", CRSK_FAM, RF_splitCustomIdx, customFunctionArray[CRSK_FAM][RF_splitCustomIdx]);
            ${trace.token}    }
            ${trace.token}  }
            ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
            ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
            ${trace.token}      RF_nativePrint("\n Outgoing custom split event type count: %10d", RF_eventType[RF_eventTypeSize]);
            ${trace.token}      RF_nativePrint("\n Outgoing custom split event times:");
            ${trace.token}      for (k = 1; k <= localEventTimeSize; k++) {
            ${trace.token}        RF_nativePrint("\n %10d %10.4f", k, userEventTime[k]);
            ${trace.token}      }
            ${trace.token}      RF_nativePrint("\n");
            ${trace.token}      RF_nativePrint("\n Outgoing custom split membership:");
            ${trace.token}      RF_nativePrint("\n      index  membrship       time     status");
            ${trace.token}      for (k = 1; k <= nonMissMembrSize; k++) {
            ${trace.token}        RF_nativePrint("\n %10d %10d %10.4f %10.4f", k, userSplitIndicator[k], userTime[k], userEvent[k]);
            ${trace.token}      }
            ${trace.token}    }  
            ${trace.token}  }  
            if ((leftSize != 0) && (rghtSize != 0)) {
              delta = customFunctionArray[CRSK_FAM][RF_splitCustomIdx](nonMissMembrSize,
                                                                       userSplitIndicator,
                                                                       userTime,
                                                                       userEvent,
                                                                       RF_eventType[RF_eventTypeSize],
                                                                       localEventTimeSize,
                                                                       userEventTime,
                                                                       NULL,
                                                                       0,
                                                                       0,
                                                                       0,
                                                                       userFeature,
                                                                       RF_yIndexZeroSize);
            }
            else {
              delta = RF_nativeNaN;
            }
            updateMaximumSplit(treeID,
                               parent,
                               delta,
                               covariate,
                               j,
                               factorFlag,
                               mwcpSizeAbsolute,
                               repMembrSize,
                               & localSplitIndicator,
                               splitVectorPtr,
                               splitInfoMax);
            if (factorFlag == FALSE) {
              priorMembrIter = currentMembrIter - 1;
            }
          }  
          free_uivector(userSort, 1, nonMissMembrSize);
          free_dvector(tempTime, 1, nonMissMembrSize);
          free_dvector (userTime, 1, nonMissMembrSize);
          free_dvector (userEvent, 1, nonMissMembrSize);
          free_dvector (userEventTime, 1, localEventTimeSize);
          free_cvector (userSplitIndicator, 1, nonMissMembrSize);
        }  
        else {
          ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
          ${trace.token}  if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
          ${trace.token}    RF_nativePrint("\nCovariate ignored due to zero localEventTimeSize:  %10d", covariate);
          ${trace.token}  }
          ${trace.token}  }
        }
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
        if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
          unstackSplitSurv(treeID,
                           parent,
                           localEventTimeCount,
                           localEventTimeIndex,
                           localEventTimeSize,
                           nodeParentEvent,
                           nodeParentAtRisk,
                           nodeLeftEvent,
                           nodeLeftAtRisk,
                           nodeRightEvent,
                           nodeRightAtRisk);
        }
      }
    }  
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      unstackSplitSurv(treeID,
                       parent,
                       localEventTimeCount,
                       localEventTimeIndex,
                       localEventTimeSize,
                       nodeParentEvent,
                       nodeParentAtRisk,
                       nodeLeftEvent,
                       nodeLeftAtRisk,
                       nodeRightEvent,
                       nodeRightAtRisk);
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
  ${trace.token}    RF_nativePrint("\ncustomSplitCompetingRisk(%10d) result:  %10d", treeID, result);
  ${trace.token}    RF_nativePrint("\ncustomSplitCompetingRisk(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
  return result;
}
