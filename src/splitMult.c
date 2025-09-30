
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "splitMult.h"
#include "splitUtil.h"
#include "nrutil.h"
char unsupervisedSplitMiss (uint       treeID,
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
  uint   deltaNorm;
  uint i, j, k, p, r;
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
    stackSplitPreliminary(repMembrSize,
                          & localSplitIndicator,
                          & splitVector);
    DistributionObj *distributionObj = stackRandomCovariates(treeID, parent);
    char   *impurity   = cvector(1, RF_ytry);
    double *mean       = dvector(1, RF_ytry);
    double *variance   = dvector(1, RF_ytry);
    uint  **parentClassProp = (uint **) new_vvector(1, RF_ytry, NRUTIL_UPTR);
    uint  **leftClassProp   = (uint **) new_vvector(1, RF_ytry, NRUTIL_UPTR);
    uint  **rghtClassProp   = (uint **) new_vvector(1, RF_ytry, NRUTIL_UPTR);
    double *sumLeft      = dvector(1, RF_ytry);
    double *sumRght      = dvector(1, RF_ytry);
    double *sumRghtSave  = dvector(1, RF_ytry);
    uint *pseudoResponseClassSize = uivector(1, RF_ytry);
    uint *pseudoResponse = uivector(1, RF_ytry);
    char **secondNonMissMembrFlag = (char **) new_vvector(1, RF_ytry, NRUTIL_CPTR);
    uint  *secondNonMissMembrSize =              uivector(1, RF_ytry);
    uint  *secondNonMissMembrLeftSize =          uivector(1, RF_ytry);
    uint  *secondNonMissMembrRghtSize =          uivector(1, RF_ytry);
    char  *tempNonMissMembrFlag = 0;
    uint  *tempNonMissMembrIndx;
    char   mResponseFlag;
    uint   localIndex = 0; 
    uint   localSize;
    char    nonMissImpuritySummary;
    double sumLeftSqr, sumRghtSqr;
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
        uint *pseudoResponseIndex = uivector(1, RF_xSize);
        for (i = 1; i <= RF_xSize; i++) {
          pseudoResponseIndex[i] = i;
        }
        pseudoResponseIndex[covariate] = pseudoResponseIndex[RF_xSize];
        localSize = RF_xSize - 1;
        for (r = 1; r <= RF_ytry; r++) {
          pseudoResponse[r] = sampleUniformlyFromVector(treeID, pseudoResponseIndex, RF_xSize, & localIndex);
          pseudoResponseIndex[localIndex] = pseudoResponseIndex[localSize];
          localSize --;
        }
        free_uivector(pseudoResponseIndex, 1, RF_xSize);
        if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
          tempNonMissMembrFlag = cvector(1, nonMissMembrSize);
          tempNonMissMembrIndx = uivector(1, nonMissMembrSize);
          for (k = 1; k <= nonMissMembrSize; k++) {
            tempNonMissMembrFlag[k] = TRUE;
            tempNonMissMembrIndx[k] = k;
          }
          for (r = 1; r <= RF_ytry; r++) {
            secondNonMissMembrFlag[r] = tempNonMissMembrFlag;
            secondNonMissMembrSize[r] = nonMissMembrSize;
          }
          nonMissImpuritySummary = FALSE;
          for (r = 1; r <= RF_ytry; r++)  {
            impurity[r] = getVariance(repMembrSize,
                                      repMembrIndx,
                                      secondNonMissMembrSize[r],
                                      tempNonMissMembrIndx,
                                      RF_observation[treeID][pseudoResponse[r]],
                                      &mean[r],
                                      &variance[r]);
            nonMissImpuritySummary = nonMissImpuritySummary | impurity[r];
            secondNonMissMembrLeftSize[r] = secondNonMissMembrRghtSize[r] = 0;
          }
          free_uivector(tempNonMissMembrIndx, 1, nonMissMembrSize);
        }
        else {
          tempNonMissMembrIndx = uivector(1, nonMissMembrSize);
          nonMissImpuritySummary = FALSE;
          for (r = 1; r <= RF_ytry; r++)  {
            secondNonMissMembrFlag[r] = cvector(1, nonMissMembrSize);
            j = 0;
            for (k = 1; k <= nonMissMembrSize; k++) {
              mResponseFlag = FALSE;
              if (RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[k]] ] > 0) {
                if (RF_mpSign[pseudoResponse[r]][RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[k]] ]] == 1) {
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
            impurity[r] = getVariance(repMembrSize,
                                      repMembrIndx,
                                      secondNonMissMembrSize[r],
                                      tempNonMissMembrIndx,
                                      RF_observation[treeID][pseudoResponse[r]],
                                      &mean[r],
                                      &variance[r]);
            nonMissImpuritySummary = nonMissImpuritySummary | impurity[r];
            secondNonMissMembrLeftSize[r] = secondNonMissMembrRghtSize[r] = 0;
          }  
          free_uivector(tempNonMissMembrIndx, 1, nonMissMembrSize);
        }  
        if (nonMissImpuritySummary) {
          for (r = 1; r <= RF_ytry; r++) {
            pseudoResponseClassSize[r] = 0;
            parentClassProp[r] = leftClassProp[r] = rghtClassProp[r] = NULL;
            sumLeft[r] = sumRght[r] = sumRghtSave[r] = 0.0;
          }
          for (r = 1; r <= RF_ytry; r++) {
            if (impurity[r]) {
              if ((RF_xType[pseudoResponse[r]] == 'B') ||
                  (RF_xType[pseudoResponse[r]] == 'C')) {
                pseudoResponseClassSize[r] = RF_xFactorSize[RF_xFactorMap[pseudoResponse[r]]];
                parentClassProp[r] = uivector(1, pseudoResponseClassSize[r]);
                leftClassProp[r]   = uivector(1, pseudoResponseClassSize[r]);
                rghtClassProp[r]   = uivector(1, pseudoResponseClassSize[r]);
                for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                  parentClassProp[r][p] = 0;
                }
                for (j = 1; j <= nonMissMembrSize; j++) {
                  if (secondNonMissMembrFlag[r][j] == TRUE) {
                    parentClassProp[r][ (uint) RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[j]] ]] ++;
                  }
                }
              }
              else {
                sumRghtSave[r] = 0.0;
                for (j = 1; j <= nonMissMembrSize; j++) {
                  if (secondNonMissMembrFlag[r][j] == TRUE) {
                    sumRghtSave[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[j]] ] - mean[r];
                  }
                }
              }
            }  
          }  
          leftSize = 0;
          priorMembrIter = 0;
          if (factorFlag == FALSE) {
            for (j = 1; j <= nonMissMembrSize; j++) {
              localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
            }
            for (r = 1; r <= RF_ytry; r++) {
              if (impurity[r]) {
                if ((RF_xType[pseudoResponse[r]] == 'B') ||
                    (RF_xType[pseudoResponse[r]] == 'C')) {
                  for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                    rghtClassProp[r][p] = parentClassProp[r][p];
                    leftClassProp[r][p] = 0;
                  }
                }
                else {
                  sumRght[r] = sumRghtSave[r];
                  sumLeft[r] = 0.0;
                  secondNonMissMembrLeftSize[r] = 0;
                  secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r];
                }
              }
            }
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
            if ((leftSize != 0) && (rghtSize != 0)) {
              delta     = 0.0;
              deltaNorm = 0;
              for (r = 1; r <= RF_ytry; r++) {
                if (impurity[r]) {
                  if (factorFlag == TRUE) {
                    if ((RF_xType[pseudoResponse[r]] == 'B') ||
                        (RF_xType[pseudoResponse[r]] == 'C')) {
                      for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                        leftClassProp[r][p] = 0;
                      }
                      for (k = 1; k <= nonMissMembrSize; k++) {
                        if (secondNonMissMembrFlag[r][indxx[k]] == TRUE) {
                          if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                            leftClassProp[r][ (uint) RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]  ++;
                            secondNonMissMembrLeftSize[r] ++;
                          }
                          else {
                          }
                        }
                      }
                      for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                        rghtClassProp[r][p] = parentClassProp[r][p] - leftClassProp[r][p];
                      }
                    }
                    else {
                      sumLeft[r] = sumRght[r] = 0.0;
                      for (k = 1; k <= nonMissMembrSize; k++) {
                        if (secondNonMissMembrFlag[r][indxx[k]] == TRUE) {
                          if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                            sumLeft[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                            secondNonMissMembrLeftSize[r] ++;
                          }
                          else {
                            sumRght[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                          }
                        }
                      }
                    }
                    secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r] - secondNonMissMembrLeftSize[r];
                  }
                  else {
                    for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                      if (secondNonMissMembrFlag[r][indxx[k]] == TRUE) {
                        if ((RF_xType[pseudoResponse[r]] == 'B') ||
                            (RF_xType[pseudoResponse[r]] == 'C')) {
                          leftClassProp[r][(uint) RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]] ++;
                          rghtClassProp[r][(uint) RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]] --;
                        }
                        else {
                          sumLeft[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                          sumRght[r] -= RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                        }
                        secondNonMissMembrLeftSize[r] ++;
                        secondNonMissMembrRghtSize[r] --;
                      }
                    }
                  }  
                  if ((secondNonMissMembrLeftSize[r] > 0) && (secondNonMissMembrRghtSize[r] > 0)) {
                    deltaNorm ++;
                    if ((RF_xType[pseudoResponse[r]] == 'B') ||
                        (RF_xType[pseudoResponse[r]] == 'C')) {
                      sumLeft[1] = sumRght[1] = 0;
                      for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                        sumLeft[1] += (double) upower(leftClassProp[r][p], 2);
                        sumRght[1] += (double) upower(rghtClassProp[r][p], 2);
                      }
                      sumLeftSqr = sumLeft[1] / secondNonMissMembrLeftSize[r];
                      sumRghtSqr  = sumRght[1] / secondNonMissMembrRghtSize[r];
                    }
                    else {
                      sumLeftSqr = pow (sumLeft[r], 2.0) / (secondNonMissMembrLeftSize[r] * variance[r]);
                      sumRghtSqr = pow (sumRght[r], 2.0) / (secondNonMissMembrRghtSize[r] * variance[r]);
                    }
                    delta += sumLeftSqr + sumRghtSqr;
                  }
                }  
              }  
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
          for (r = 1; r <= RF_ytry; r++) {
            if (impurity[r]) {
              if ((RF_xType[pseudoResponse[r]] == 'B') ||
                  (RF_xType[pseudoResponse[r]] == 'C')) {
                free_uivector (parentClassProp[r], 1, pseudoResponseClassSize[r]);
                free_uivector (leftClassProp[r],   1, pseudoResponseClassSize[r]);
                free_uivector (rghtClassProp[r],   1, pseudoResponseClassSize[r]);
              }
              else {
              }
            }
          }  
        }  
        if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
          free_cvector(tempNonMissMembrFlag, 1, nonMissMembrSize);
        }
        else {
          for (r = 1; r <= RF_ytry; r++)  {
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
    free_new_vvector(parentClassProp, 1, RF_ytry, NRUTIL_UPTR);
    free_new_vvector(leftClassProp,   1, RF_ytry, NRUTIL_UPTR);
    free_new_vvector(rghtClassProp,   1, RF_ytry, NRUTIL_UPTR);
    free_dvector(sumLeft,     1, RF_ytry);
    free_dvector(sumRght,     1, RF_ytry);
    free_dvector(sumRghtSave, 1, RF_ytry);
    free_uivector(pseudoResponseClassSize, 1, RF_ytry);
    free_uivector(pseudoResponse, 1, RF_ytry);
    free_new_vvector(secondNonMissMembrFlag,  1, RF_ytry, NRUTIL_CPTR);
    free_uivector(secondNonMissMembrSize,     1, RF_ytry);
    free_uivector(secondNonMissMembrLeftSize, 1, RF_ytry);
    free_uivector(secondNonMissMembrRghtSize, 1, RF_ytry);
    unstackRandomCovariates(treeID, distributionObj);
    unstackSplitPreliminary(repMembrSize,
                            localSplitIndicator,
                            splitVector);
    free_cvector(impurity,   1, RF_ytry);
    free_dvector(mean,     1, RF_ytry);
    free_dvector(variance, 1, RF_ytry);
  }  
  unstackPreSplit(preliminaryResult,
                  parent,
                  multImpFlag,
                  TRUE);  
  result = summarizeSplitResult(splitInfoMax);
  return result;
}
char unsupervisedSplitNew (uint       treeID,
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
  uint   deltaNorm;
  uint i, j, k, p, r;
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
    stackSplitPreliminary(repMembrSize,
                          & localSplitIndicator,
                          & splitVector);
    DistributionObj *distributionObj = stackRandomCovariates(treeID, parent);
    char   *impurity   = cvector(1, RF_ytry);
    double *mean       = dvector(1, RF_ytry);
    double *variance   = dvector(1, RF_ytry);
    uint  **parentClassProp = (uint **) new_vvector(1, RF_ytry, NRUTIL_UPTR);
    uint  **leftClassProp   = (uint **) new_vvector(1, RF_ytry, NRUTIL_UPTR);
    uint  **rghtClassProp   = (uint **) new_vvector(1, RF_ytry, NRUTIL_UPTR);
    double *sumLeft      = dvector(1, RF_ytry);
    double *sumRght      = dvector(1, RF_ytry);
    double *sumRghtSave  = dvector(1, RF_ytry);
    double *sumLeftSqr      = dvector(1, RF_ytry);
    double *sumRghtSqr      = dvector(1, RF_ytry);
    double *sumRghtSqrSave  = dvector(1, RF_ytry);
    uint *pseudoResponseClassSize = uivector(1, RF_ytry);
    uint *pseudoResponse = uivector(1, RF_ytry);
    uint   localIndex = 0; 
    uint   localSize;
    char    impuritySummary;
    double partialLeft, partialRght;
    double partialLeftSqr, partialRghtSqr;
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
        uint *pseudoResponseIndex = uivector(1, RF_xSize);
        for (i = 1; i <= RF_xSize; i++) {
          pseudoResponseIndex[i] = i;
        }
        pseudoResponseIndex[covariate] = pseudoResponseIndex[RF_xSize];
        localSize = RF_xSize - 1;
        for (r = 1; r <= RF_ytry; r++) {
          pseudoResponse[r] = sampleUniformlyFromVector(treeID, pseudoResponseIndex, RF_xSize, & localIndex);
          pseudoResponseIndex[localIndex] = pseudoResponseIndex[localSize];
          localSize --;
        }
        free_uivector(pseudoResponseIndex, 1, RF_xSize);
        impuritySummary = FALSE;
        for (r = 1; r <= RF_ytry; r++)  {
          impurity[r] = getVariance(repMembrSize,
                                    repMembrIndx,
                                    0,
                                    NULL,
                                    RF_observation[treeID][pseudoResponse[r]],
                                    &mean[r],
                                    &variance[r]);
          impuritySummary = impuritySummary | impurity[r];
        }
        if (impuritySummary) {
          for (r = 1; r <= RF_ytry; r++) {
            pseudoResponseClassSize[r] = 0;
            parentClassProp[r] = leftClassProp[r] = rghtClassProp[r] = NULL;
            if ((RF_xType[pseudoResponse[r]] == 'B') ||
                (RF_xType[pseudoResponse[r]] == 'C')) {
              pseudoResponseClassSize[r] = RF_xFactorSize[RF_xFactorMap[pseudoResponse[r]]];
              parentClassProp[r] = uivector(1, pseudoResponseClassSize[r]);
              leftClassProp[r]   = uivector(1, pseudoResponseClassSize[r]);
              rghtClassProp[r]   = uivector(1, pseudoResponseClassSize[r]);
            }
            else {
            }
          }  
          for (r = 1; r <= RF_ytry; r++) {
            if (impurity[r]) {
              if ((RF_xType[pseudoResponse[r]] == 'B') ||
                  (RF_xType[pseudoResponse[r]] == 'C')) {
                for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                  parentClassProp[r][p] = 0;
                }
                for (j = 1; j <= repMembrSize; j++) {
                  parentClassProp[r][ (uint) RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[j] ]] ++;
                }
              }
              else {
                sumRghtSave[r] = 0.0;
                sumRghtSqrSave[r] = 0.0;
                switch(RF_splitRule) {
                case USPV_WT_OFF:
                case USPV_WT_HVY:
                  for (j = 1; j <= repMembrSize; j++) {
                    sumRghtSave[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[j] ];
                    sumRghtSqrSave[r] += pow(RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[j] ], 2.0);
                  }
                  break;
                default:
                  for (j = 1; j <= repMembrSize; j++) {
                    sumRghtSave[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[j] ] - mean[r];
                  }
                  break;
                }      
              }
            }  
          }  
          leftSize = 0;
          priorMembrIter = 0;
          if (factorFlag == FALSE) {
            for (j = 1; j <= repMembrSize; j++) {
              localSplitIndicator[ j ] = RIGHT;
            }
            for (r = 1; r <= RF_ytry; r++) {
              if (impurity[r]) {
                if ((RF_xType[pseudoResponse[r]] == 'B') ||
                    (RF_xType[pseudoResponse[r]] == 'C')) {
                  for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                    rghtClassProp[r][p] = parentClassProp[r][p];
                    leftClassProp[r][p] = 0;
                  }
                }
                else {
                  sumRght[r] = sumRghtSave[r];
                  sumRghtSqr[r] = sumRghtSqrSave[r];
                  sumLeft[r] = 0.0;
                  sumLeftSqr[r] = 0.0;
                }
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
            rghtSize = repMembrSize - leftSize;
            if ((leftSize != 0) && (rghtSize != 0)) {
              delta     = 0.0;
              deltaNorm = 0;
              for (r = 1; r <= RF_ytry; r++) {
                if (impurity[r]) {
                  if (factorFlag == TRUE) {
                    if ((RF_xType[pseudoResponse[r]] == 'B') ||
                        (RF_xType[pseudoResponse[r]] == 'C')) {
                      for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                        leftClassProp[r][p] = 0;
                      }
                      for (k = 1; k <= repMembrSize; k++) {
                        if (localSplitIndicator[ indxx[k] ] == LEFT) {
                          leftClassProp[r][ (uint) RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[ indxx[k]] ]]  ++;
                        }
                        else {
                        }
                      }
                      for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                        rghtClassProp[r][p] = parentClassProp[r][p] - leftClassProp[r][p];
                      }
                    }
                    else {
                      switch (RF_splitRule) {
                      case USPV_WT_OFF:
                      case USPV_WT_HVY:
                        sumLeft[r] = sumRght[r] = 0.0;
                        sumLeftSqr[r] = sumRghtSqr[r] = 0.0;
                        for (k = 1; k <= repMembrSize; k++) {
                          if (localSplitIndicator[ indxx[k] ] == LEFT) {
                            sumLeft[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[ indxx[k]] ];
                            sumLeftSqr[r] += pow(RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[ indxx[k] ]], 2.0);
                          }
                          else {
                            sumRght[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[ indxx[k]] ];
                            sumRghtSqr[r] += pow(RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[ indxx[k] ]], 2.0);
                          }
                        }
                        break;
                      default:
                        sumLeft[r] = sumRght[r] = 0.0;
                        for (k = 1; k <= repMembrSize; k++) {
                          if (localSplitIndicator[ indxx[k] ] == LEFT) {
                            sumLeft[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[ indxx[k] ]] - mean[r];
                          }
                          else {
                            sumRght[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[ indxx[k] ]] - mean[r];
                          }
                        }
                        break;
                      }
                    }
                  }
                  else {
                    for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                      if ((RF_xType[pseudoResponse[r]] == 'B') ||
                          (RF_xType[pseudoResponse[r]] == 'C')) {
                        leftClassProp[r][ (uint) RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[ indxx[k]] ]] ++;
                        rghtClassProp[r][ (uint) RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[ indxx[k]] ]]  --;
                      }
                      else {
                        switch(RF_splitRule) {
                        case USPV_WT_OFF:
                        case USPV_WT_HVY:
                          sumLeft[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[indxx[k]] ];
                          sumLeftSqr[r] += pow(RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[indxx[k]] ], 2.0);
                          sumRght[r] -= RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[indxx[k]] ];
                          sumRghtSqr[r] -= pow(RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[indxx[k]] ], 2.0);                              
                          break;
                        default:
                          sumLeft[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[indxx[k]] ] - mean[r];
                          sumRght[r] -= RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[indxx[k]] ] - mean[r];
                          break;
                        }
                      }
                    }
                  }  
                  if ((leftSize > 0) && (rghtSize > 0)) {
                    deltaNorm ++;
                    if ((RF_xType[pseudoResponse[r]] == 'B') ||
                        (RF_xType[pseudoResponse[r]] == 'C')) {
                      partialLeft = partialRght = 0.0;
                      switch(RF_splitRule) {
                      case USPV_WT_OFF:
                        for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                          partialLeft += pow((double) leftClassProp[r][p] / leftSize, 2.0);
                          partialRght += pow((double) (rghtClassProp[r][p]) / rghtSize, 2.0);
                        }
                        delta += partialLeft + partialRght;
                        break;
                      case USPV_WT_HVY:
                        for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                          partialLeft += (double) upower(leftClassProp[r][p], 2);
                          partialRght += (double) upower(rghtClassProp[r][p], 2);
                        }
                        delta +=
                          (partialLeft / (double) (upower(repMembrSize, 2))) +
                          (partialRght / (double) (upower(repMembrSize, 2))) -
                          pow((double) leftSize / repMembrSize, 2.0) -
                          pow((double) rghtSize / repMembrSize, 2.0) + 2.0;
                        break;
                      default:
                        for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                          partialLeft += (double) upower(leftClassProp[r][p], 2);
                          partialRght += (double) upower(rghtClassProp[r][p], 2);
                        }
                        delta += ((partialLeft / leftSize) + (partialRght / rghtSize)) / repMembrSize;
                        break;
                      }
                    }
                    else {
                      switch(RF_splitRule) {
                      case USPV_WT_OFF:
                        partialLeft = pow (sumLeft[r], 2.0) / (double) upower (leftSize, 2);
                        partialRght = pow (sumRght[r], 2.0) / (double) upower (rghtSize, 2);
                        partialLeftSqr = sumLeftSqr[r] / leftSize;
                        partialRghtSqr = sumRghtSqr[r] / rghtSize;
                        delta += partialLeft + partialRght - partialLeftSqr - partialRghtSqr;
                        break;
                      case USPV_WT_HVY:
                        partialLeft = pow (sumLeft[r], 2.0) / (double) upower (repMembrSize, 2);
                        partialRght = pow (sumRght[r], 2.0) / (double) upower (repMembrSize, 2);
                        partialLeftSqr = sumLeftSqr[r] * leftSize / (double) upower (repMembrSize, 2);
                        partialRghtSqr = sumRghtSqr[r] * rghtSize / (double) upower (repMembrSize, 2);
                        delta += partialLeft + partialRght - partialLeftSqr - partialRghtSqr;
                        break;
                      default:
                        partialLeft = pow (sumLeft[r], 2.0) / (leftSize * variance[r]);
                        partialRght = pow (sumRght[r], 2.0) / (rghtSize * variance[r]);
                        delta += partialLeft + partialRght;
                        break;
                      }
                    }
                  }
                }  
              }  
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
          for (r = 1; r <= RF_ytry; r++) {
            if ((RF_xType[pseudoResponse[r]] == 'B') ||
                (RF_xType[pseudoResponse[r]] == 'C')) {
              free_uivector (parentClassProp[r], 1, pseudoResponseClassSize[r]);
              free_uivector (leftClassProp[r],   1, pseudoResponseClassSize[r]);
              free_uivector (rghtClassProp[r],   1, pseudoResponseClassSize[r]);
            }
            else {
            }
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
    free_new_vvector(parentClassProp, 1, RF_ytry, NRUTIL_UPTR);
    free_new_vvector(leftClassProp,   1, RF_ytry, NRUTIL_UPTR);
    free_new_vvector(rghtClassProp,   1, RF_ytry, NRUTIL_UPTR);
    free_dvector(sumLeft,     1, RF_ytry);
    free_dvector(sumRght,     1, RF_ytry);
    free_dvector(sumRghtSave, 1, RF_ytry);
    free_dvector(sumLeftSqr,     1, RF_ytry);
    free_dvector(sumRghtSqr,     1, RF_ytry);
    free_dvector(sumRghtSqrSave, 1, RF_ytry);
    free_uivector(pseudoResponseClassSize, 1, RF_ytry);
    free_uivector(pseudoResponse, 1, RF_ytry);
    unstackRandomCovariates(treeID, distributionObj);
    unstackSplitPreliminary(repMembrSize,
                            localSplitIndicator,
                            splitVector);
    free_cvector(impurity,   1, RF_ytry);
    free_dvector(mean,       1, RF_ytry);
    free_dvector(variance,   1, RF_ytry);
  }  
  unstackPreSplit(preliminaryResult,
                  parent,
                  multImpFlag,
                  TRUE);  
  result = summarizeSplitResult(splitInfoMax);
  return result;
}
char multivariateSplitOld (uint       treeID,
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
  uint   deltaNorm;
  uint j, k, p, r;
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
    uint  nonMissMembrSize;
    uint *nonMissMembrIndx;
    char   *impurity   = cvector(1, RF_ySize);
    double *mean       = dvector(1, RF_ySize);
    double *variance   = dvector(1, RF_ySize);
    char impuritySummary;
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      impuritySummary = FALSE;
      for (r = 1; r <= RF_ySize; r++)  {
        impurity[r] = getVariance(repMembrSize,
                                  repMembrIndx,
                                  0,
                                  NULL,
                                  RF_response[treeID][r],
                                  &mean[r],
                                  &variance[r]);
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
      uint **parentClassProp = (uint **) new_vvector(1, RF_ySize, NRUTIL_UPTR);
      uint **leftClassProp   = (uint **) new_vvector(1, RF_ySize, NRUTIL_UPTR);
      uint **rghtClassProp   = (uint **) new_vvector(1, RF_ySize, NRUTIL_UPTR);
      double *sumLeft         = dvector(1, RF_ySize);
      double *sumRght         = dvector(1, RF_ySize);
      double *sumRghtSave     = dvector(1, RF_ySize);
      char **secondNonMissMembrFlag = (char **) new_vvector(1, RF_ySize, NRUTIL_CPTR);
      uint  *secondNonMissMembrSize =           uivector(1, RF_ySize);
      uint  *secondNonMissMembrLeftSize =       uivector(1, RF_ySize);
      uint  *secondNonMissMembrRghtSize =       uivector(1, RF_ySize);
      for (r = 1; r <= RF_ySize; r++) {
        parentClassProp[r] = leftClassProp[r] = rghtClassProp[r] = NULL;
        if ((RF_rType[r] == 'B') ||
            (RF_rType[r] == 'I') ||
            (RF_rType[r] == 'C')) {
          parentClassProp[r] = uivector(1, RF_classLevelSize[RF_rFactorMap[r]]);
          leftClassProp[r]   = uivector(1, RF_classLevelSize[RF_rFactorMap[r]]);
          rghtClassProp[r]   = uivector(1, RF_classLevelSize[RF_rFactorMap[r]]);
        }
        else {
        }
      }  
      char  *tempNonMissMembrFlag = 0;
      uint  *tempNonMissMembrIndx;
      char   mResponseFlag;
      char   nonMissImpuritySummary;
      double partialLeft, partialRght;
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
              impurity[r] = getVariance(repMembrSize,
                                        repMembrIndx,
                                        secondNonMissMembrSize[r],
                                        tempNonMissMembrIndx,
                                        RF_response[treeID][r],
                                        &mean[r],
                                        &variance[r]);
              nonMissImpuritySummary = nonMissImpuritySummary | impurity[r];
              secondNonMissMembrLeftSize[r] = secondNonMissMembrRghtSize[r] = 0;
            }  
            free_uivector(tempNonMissMembrIndx, 1, nonMissMembrSize);
          }  
          if (nonMissImpuritySummary) {
            for (r = 1; r <= RF_ySize; r++) {
              if (impurity[r]) {
                if ((RF_rType[r] == 'B') ||
                    (RF_rType[r] == 'I') ||
                    (RF_rType[r] == 'C')) {
                  for (p=1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                    parentClassProp[r][p] = 0;
                  }
                  for (j = 1; j <= nonMissMembrSize; j++) {
                    if (secondNonMissMembrFlag[r][j] == TRUE) {
                      parentClassProp[r][RF_classLevelIndex[RF_rFactorMap[r]][(uint) RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[j]] ]]] ++;
                    }
                  }
                }
                else {
                  sumRghtSave[r] = 0.0;
                  for (j = 1; j <= nonMissMembrSize; j++) {
                    if (secondNonMissMembrFlag[r][j] == TRUE) {
                      sumRghtSave[r] += RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[j]] ] - mean[r];
                    }
                  }
                }
              }  
            }  
            leftSize = 0;
            priorMembrIter = 0;
            if (factorFlag == FALSE) {
              for (j = 1; j <= nonMissMembrSize; j++) {
                localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
              }
              for (r = 1; r <= RF_ySize; r++) {
                if (impurity[r]) {
                  if ((RF_rType[r] == 'B') ||
                      (RF_rType[r] == 'I') ||
                      (RF_rType[r] == 'C')) {
                    for (p=1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                      rghtClassProp[r][p] = parentClassProp[r][p];
                      leftClassProp[r][p] = 0;
                    }
                  }
                  else {
                    sumRght[r] = sumRghtSave[r];
                    sumLeft[r] = 0.0;
                  }
                  secondNonMissMembrLeftSize[r] = 0;
                  secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r];
                }
              }
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
              if ((leftSize != 0) && (rghtSize != 0)) {
                delta     = 0.0;
                deltaNorm = 0;
                for (r = 1; r <= RF_ySize; r++) {
                  if (impurity[r]) {
                    if (factorFlag == TRUE) {
                      if ((RF_rType[r] == 'B') ||
                          (RF_rType[r] == 'I') ||
                          (RF_rType[r] == 'C')) {
                        for (p=1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                          leftClassProp[r][p] = 0;
                        }
                        for (k = 1; k <= nonMissMembrSize; k++) {
                          if (secondNonMissMembrFlag[r][indxx[k]] == TRUE) {
                            if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                              leftClassProp[r][RF_classLevelIndex[RF_rFactorMap[r]][(uint) RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] ++;
                              secondNonMissMembrLeftSize[r] ++;
                            }
                            else {
                            }
                          }
                        }
                        for (p=1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                          rghtClassProp[r][p] = parentClassProp[r][p] - leftClassProp[r][p];
                        }
                      }
                      else {
                        sumLeft[r] = sumRght[r] = 0.0;
                        for (k = 1; k <= nonMissMembrSize; k++) {
                          if (secondNonMissMembrFlag[r][indxx[k]] == TRUE) {
                            if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                              sumLeft[r] += RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                              secondNonMissMembrLeftSize[r] ++;
                            }
                            else {
                              sumRght[r] += RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                            }
                          }
                        }
                      }
                      secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r] - secondNonMissMembrLeftSize[r];
                    }
                    else {
                      for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                        if (secondNonMissMembrFlag[r][indxx[k]] == TRUE) {
                          if ((RF_rType[r] == 'B') ||
                              (RF_rType[r] == 'I') ||
                              (RF_rType[r] == 'C')) {
                            leftClassProp[r][RF_classLevelIndex[RF_rFactorMap[r]][(uint) RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] ++;
                            rghtClassProp[r][RF_classLevelIndex[RF_rFactorMap[r]][(uint) RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] --;
                          }
                          else {
                            sumLeft[r] += RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                            sumRght[r] -= RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                          }
                          secondNonMissMembrLeftSize[r] ++;
                          secondNonMissMembrRghtSize[r] --;
                        }
                      }
                    }  
                    if ((secondNonMissMembrLeftSize[r] > 0) && (secondNonMissMembrRghtSize[r] > 0)) {
                      deltaNorm ++;
                      if ((RF_rType[r] == 'B') ||
                          (RF_rType[r] == 'I') ||
                          (RF_rType[r] == 'C')) {
                        partialLeft = partialRght = 0;
                        for (p = 1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                          partialLeft += (double) upower(leftClassProp[r][p], 2);
                          partialRght += (double) upower(rghtClassProp[r][p], 2);
                        }
                        partialLeft = partialLeft / secondNonMissMembrLeftSize[r];
                        partialRght = partialRght / secondNonMissMembrRghtSize[r];
                      }
                      else {
                        partialLeft = pow (sumLeft[r], 2.0) / (secondNonMissMembrLeftSize[r] * variance[r]);
                        partialRght = pow (sumRght[r], 2.0) / (secondNonMissMembrRghtSize[r] * variance[r]);
                      }
                      delta += partialLeft + partialRght;
                    }
                  }  
                }  
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
      for (r = 1; r <= RF_ySize; r++) {
        if ((RF_rType[r] == 'B') ||
            (RF_rType[r] == 'I') ||
            (RF_rType[r] == 'C')) {
          free_uivector (parentClassProp[r], 1, RF_classLevelSize[RF_rFactorMap[r]]);
          free_uivector (leftClassProp[r], 1, RF_classLevelSize[RF_rFactorMap[r]]);
          free_uivector (rghtClassProp[r], 1, RF_classLevelSize[RF_rFactorMap[r]]);
        }
        else {
        }
      }
      free_new_vvector(parentClassProp, 1, RF_ySize, NRUTIL_UPTR);
      free_new_vvector(leftClassProp,   1, RF_ySize, NRUTIL_UPTR);
      free_new_vvector(rghtClassProp,   1, RF_ySize, NRUTIL_UPTR);
      free_dvector(sumLeft,     1, RF_ySize);
      free_dvector(sumRght,     1, RF_ySize);
      free_dvector(sumRghtSave, 1, RF_ySize);
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
  return result;
}
char multivariateSplitNew (uint       treeID,
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
  uint   deltaNorm;
  uint j, k, p, r;
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
    uint  nonMissMembrSize;
    uint *nonMissMembrIndx;
    uint ySize;
    if (RF_ytry == 0) {
      ySize = RF_ySize;
    }
    else {
      ySize = RF_ytry;
    }
    char   *impurity   = cvector(1, ySize);
    double *mean       = dvector(1, ySize);
    double *variance   = dvector(1, ySize);
    DistributionObj *distributionObj = stackRandomResponses(treeID, parent);
    uint *response = uivector(1, ySize);
    uint  responseCount = 0;
    selectRandomResponses(treeID,
                          parent,
                          distributionObj,
                          response,
                          & responseCount);
    unstackRandomResponses(treeID, distributionObj);
    char impuritySummary;
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      impuritySummary = FALSE;
      for (r = 1; r <= ySize; r++)  {
        impurity[r] = getVariance(repMembrSize,
                                  repMembrIndx,
                                  0,
                                  NULL,
                                  RF_response[treeID][response[r]],
                                  &mean[r],
                                  &variance[r]);
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
      uint **parentClassProp = (uint **) new_vvector(1, ySize, NRUTIL_UPTR);
      uint **leftClassProp   = (uint **) new_vvector(1, ySize, NRUTIL_UPTR);
      uint **rghtClassProp   = (uint **) new_vvector(1, ySize, NRUTIL_UPTR);
      double *sumLeft         = dvector(1, ySize);
      double *sumRght         = dvector(1, ySize);
      double *sumRghtSave     = dvector(1, ySize);
      char **secondNonMissMembrFlag = (char **) new_vvector(1, ySize, NRUTIL_CPTR);
      uint  *secondNonMissMembrSize =           uivector(1, ySize);
      uint  *secondNonMissMembrLeftSize =       uivector(1, ySize);
      uint  *secondNonMissMembrRghtSize =       uivector(1, ySize);
      for (r = 1; r <= ySize; r++) {
        parentClassProp[r] = leftClassProp[r] = rghtClassProp[r] = NULL;
        if ((RF_rType[response[r]] == 'B') ||
            (RF_rType[response[r]] == 'I') ||
            (RF_rType[response[r]] == 'C')) {
          parentClassProp[r] = uivector(1, RF_classLevelSize[RF_rFactorMap[response[r]]]);
          leftClassProp[r]   = uivector(1, RF_classLevelSize[RF_rFactorMap[response[r]]]);
          rghtClassProp[r]   = uivector(1, RF_classLevelSize[RF_rFactorMap[response[r]]]);
        }
        else {
        }
      }  
      char  *tempNonMissMembrFlag = 0;
      uint  *tempNonMissMembrIndx;
      char   mResponseFlag;
      char   nonMissImpuritySummary;
      double partialLeft, partialRght;
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
          if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
            tempNonMissMembrFlag = cvector(1, nonMissMembrSize);
            for (k = 1; k <= nonMissMembrSize; k++) {
              tempNonMissMembrFlag[k] = TRUE;
            }
            for (r = 1; r <= ySize; r++) {
              secondNonMissMembrFlag[r] = tempNonMissMembrFlag;
              secondNonMissMembrSize[r] = nonMissMembrSize;
            }
            nonMissImpuritySummary = TRUE;
          }
          else {
            tempNonMissMembrIndx = uivector(1, nonMissMembrSize);
            nonMissImpuritySummary = FALSE;
            for (r = 1; r <= ySize; r++)  {
              secondNonMissMembrFlag[r] = cvector(1, nonMissMembrSize);
              j = 0;
              for (k = 1; k <= nonMissMembrSize; k++) {
                mResponseFlag = FALSE;
                if (RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[k]] ] > 0) {
                  if (RF_mpSign[response[r]][RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[k]] ]] == 1) {
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
              impurity[r] = getVariance(repMembrSize,
                                        repMembrIndx,
                                        secondNonMissMembrSize[r],
                                        tempNonMissMembrIndx,
                                        RF_response[treeID][response[r]],
                                        &mean[r],
                                        &variance[r]);
              nonMissImpuritySummary = nonMissImpuritySummary | impurity[r];
              secondNonMissMembrLeftSize[r] = secondNonMissMembrRghtSize[r] = 0;
            }  
            free_uivector(tempNonMissMembrIndx, 1, nonMissMembrSize);
          }  
          if (nonMissImpuritySummary) {
            for (r = 1; r <= ySize; r++) {
              if (impurity[r]) {
                if ((RF_rType[response[r]] == 'B') ||
                    (RF_rType[response[r]] == 'I') ||
                    (RF_rType[response[r]] == 'C')) {
                  for (p = 1; p <= RF_classLevelSize[RF_rFactorMap[response[r]]]; p++) {
                    parentClassProp[r][p] = 0;
                  }
                  for (j = 1; j <= nonMissMembrSize; j++) {
                    if (secondNonMissMembrFlag[r][j] == TRUE) {
                      parentClassProp[r][RF_classLevelIndex[RF_rFactorMap[response[r]]][(uint) RF_response[treeID][response[r]][ repMembrIndx[nonMissMembrIndx[j]] ]]] ++;
                    }
                  }
                }
                else {
                  sumRghtSave[r] = 0.0;
                  for (j = 1; j <= nonMissMembrSize; j++) {
                    if (secondNonMissMembrFlag[r][j] == TRUE) {
                      sumRghtSave[r] += RF_response[treeID][response[r]][ repMembrIndx[nonMissMembrIndx[j]] ] - mean[r];
                    }
                  }
                }
              }  
            }  
            leftSize = 0;
            priorMembrIter = 0;
            if (factorFlag == FALSE) {
              for (j = 1; j <= nonMissMembrSize; j++) {
                localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
              }
              for (r = 1; r <= ySize; r++) {
                if (impurity[r]) {
                  if ((RF_rType[response[r]] == 'B') ||
                      (RF_rType[response[r]] == 'I') ||
                      (RF_rType[response[r]] == 'C')) {
                    for (p = 1; p <= RF_classLevelSize[RF_rFactorMap[response[r]]]; p++) {
                      rghtClassProp[r][p] = parentClassProp[r][p];
                      leftClassProp[r][p] = 0;
                    }
                  }
                  else {
                    sumRght[r] = sumRghtSave[r];
                    sumLeft[r] = 0.0;
                  }
                  secondNonMissMembrLeftSize[r] = 0;
                  secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r];
                }
              }
            }
            deltaMax =  RF_nativeNaN;
            indexMax =  0;
            for (j = 1; j < splitLength; j++) {
              if (factorFlag == TRUE) {
                priorMembrIter = 0;
                leftSize = 0;
                for (r = 1; r <= ySize; r++) {
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
              if ((leftSize != 0) && (rghtSize != 0)) {
                delta     = 0.0;
                deltaNorm = 0;
                for (r = 1; r <= ySize; r++) {
                  if (impurity[r]) {
                    if (factorFlag == TRUE) {
                      if ((RF_rType[response[r]] == 'B') ||
                          (RF_rType[response[r]] == 'I') ||
                          (RF_rType[response[r]] == 'C')) {
                        for (p=1; p <= RF_classLevelSize[RF_rFactorMap[response[r]]]; p++) {
                          leftClassProp[r][p] = 0;
                        }
                        for (k = 1; k <= nonMissMembrSize; k++) {
                          if (secondNonMissMembrFlag[r][indxx[k]] == TRUE) {
                            if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                              leftClassProp[r][RF_classLevelIndex[RF_rFactorMap[response[r]]][(uint) RF_response[treeID][response[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] ++;
                              secondNonMissMembrLeftSize[r] ++;
                            }
                            else {
                            }
                          }
                        }
                        for (p = 1; p <= RF_classLevelSize[RF_rFactorMap[response[r]]]; p++) {
                          rghtClassProp[r][p] = parentClassProp[r][p] - leftClassProp[r][p];
                        }
                      }
                      else {
                        sumLeft[r] = sumRght[r] = 0.0;
                        for (k = 1; k <= nonMissMembrSize; k++) {
                          if (secondNonMissMembrFlag[r][indxx[k]] == TRUE) {
                            if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                              sumLeft[r] += RF_response[treeID][response[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                              secondNonMissMembrLeftSize[r] ++;
                            }
                            else {
                              sumRght[r] += RF_response[treeID][response[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                            }
                          }
                        }
                      }
                      secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r] - secondNonMissMembrLeftSize[r];
                    }
                    else {
                      for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                        if (secondNonMissMembrFlag[r][indxx[k]] == TRUE) {
                          if ((RF_rType[response[r]] == 'B') ||
                              (RF_rType[response[r]] == 'I') ||
                              (RF_rType[response[r]] == 'C')) {
                            leftClassProp[r][RF_classLevelIndex[RF_rFactorMap[response[r]]][(uint) RF_response[treeID][response[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] ++;
                            rghtClassProp[r][RF_classLevelIndex[RF_rFactorMap[response[r]]][(uint) RF_response[treeID][response[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] --;
                          }
                          else {
                            sumLeft[r] += RF_response[treeID][response[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                            sumRght[r] -= RF_response[treeID][response[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                          }
                          secondNonMissMembrLeftSize[r] ++;
                          secondNonMissMembrRghtSize[r] --;
                        }
                      }
                    }  
                    if ((secondNonMissMembrLeftSize[r] > 0) && (secondNonMissMembrRghtSize[r] > 0)) {
                      deltaNorm ++;
                      if ((RF_rType[response[r]] == 'B') ||
                          (RF_rType[response[r]] == 'I') ||
                          (RF_rType[response[r]] == 'C')) {
                        partialLeft = partialRght = 0;
                        for (p = 1; p <= RF_classLevelSize[RF_rFactorMap[response[r]]]; p++) {
                          partialLeft += (double) upower(leftClassProp[r][p], 2);
                          partialRght += (double) upower(rghtClassProp[r][p], 2);
                        }
                        partialLeft = partialLeft / secondNonMissMembrLeftSize[r];
                        partialRght = partialRght / secondNonMissMembrRghtSize[r];
                        delta += (partialLeft + partialRght) / secondNonMissMembrSize[r];
                      }
                      else {
                        partialLeft = pow (sumLeft[r], 2.0) / (secondNonMissMembrLeftSize[r] * variance[r]);
                        partialRght = pow (sumRght[r], 2.0) / (secondNonMissMembrRghtSize[r] * variance[r]);
                        delta += partialLeft + partialRght;
                      }
                    }
                  }  
                }  
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
          }  
          if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
            free_cvector(tempNonMissMembrFlag, 1, nonMissMembrSize);
          }
          else {
            for (r = 1; r <= ySize; r++)  {
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
      for (r = 1; r <= ySize; r++) {
        if ((RF_rType[response[r]] == 'B') ||
            (RF_rType[response[r]] == 'I') ||
            (RF_rType[response[r]] == 'C')) {
          free_uivector (parentClassProp[r], 1, RF_classLevelSize[RF_rFactorMap[response[r]]]);
          free_uivector (leftClassProp[r], 1, RF_classLevelSize[RF_rFactorMap[response[r]]]);
          free_uivector (rghtClassProp[r], 1, RF_classLevelSize[RF_rFactorMap[response[r]]]);
        }
        else {
        }
      }
      free_new_vvector(parentClassProp, 1, ySize, NRUTIL_UPTR);
      free_new_vvector(leftClassProp,   1, ySize, NRUTIL_UPTR);
      free_new_vvector(rghtClassProp,   1, ySize, NRUTIL_UPTR);
      free_dvector(sumLeft,     1, ySize);
      free_dvector(sumRght,     1, ySize);
      free_dvector(sumRghtSave, 1, ySize);
      free_new_vvector(secondNonMissMembrFlag,  1, ySize, NRUTIL_CPTR);
      free_uivector(secondNonMissMembrSize,     1, ySize);
      free_uivector(secondNonMissMembrLeftSize, 1, ySize);
      free_uivector(secondNonMissMembrRghtSize, 1, ySize);
      unstackRandomCovariates(treeID, distributionObj);
      unstackSplitPreliminary(repMembrSize,
                              localSplitIndicator,
                              splitVector);
    }  
    free_cvector(impurity,   1, ySize);
    free_dvector(mean,       1, ySize);
    free_dvector(variance,   1, ySize);
    free_uivector(response, 1, ySize);
  }  
  unstackPreSplit(preliminaryResult,
                  parent,
                  multImpFlag,
                  TRUE);  
  result = summarizeSplitResult(splitInfoMax);
  return result;
}
char multivariateSplitNew3 (uint       treeID,
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
  uint   deltaNorm;
  uint j, k, p, r;
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
    uint ySize;
    if (RF_ytry == 0) {
      ySize = RF_ySize;
    }
    else {
      ySize = RF_ytry;
    }
    char   *impurity   = cvector(1, ySize);
    double *mean       = dvector(1, ySize);
    double *variance   = dvector(1, ySize);
    DistributionObj *distributionObj = stackRandomResponses(treeID, parent);
    uint *response = uivector(1, ySize);
    uint  responseCount = 0;
    selectRandomResponses(treeID,
                          parent,
                          distributionObj,
                          response,
                          & responseCount);
    unstackRandomResponses(treeID, distributionObj);
    char impuritySummary;
    impuritySummary = FALSE;
    for (r = 1; r <= ySize; r++)  {
      impurity[r] = getVariance(repMembrSize,
                                repMembrIndx,
                                0,
                                NULL,
                                RF_response[treeID][response[r]],
                                &mean[r],
                                &variance[r]);
      impuritySummary = impuritySummary | impurity[r];
    }
    if (impuritySummary) {
      stackSplitPreliminary(repMembrSize,
                            & localSplitIndicator,
                            & splitVector);
      DistributionObj *distributionObj = stackRandomCovariates(treeID, parent);
      uint **parentClassProp = (uint **) new_vvector(1, ySize, NRUTIL_UPTR);
      uint **leftClassProp   = (uint **) new_vvector(1, ySize, NRUTIL_UPTR);
      uint **rghtClassProp   = (uint **) new_vvector(1, ySize, NRUTIL_UPTR);
      double *sumLeft         = dvector(1, ySize);
      double *sumRght         = dvector(1, ySize);
      double *sumRghtSave     = dvector(1, ySize);
      double *sumLeftSqr      = dvector(1, ySize);
      double *sumRghtSqr      = dvector(1, ySize);
      double *sumRghtSqrSave  = dvector(1, ySize);
      for (r = 1; r <= ySize; r++) {
        parentClassProp[r] = leftClassProp[r] = rghtClassProp[r] = NULL;
        if ((RF_rType[response[r]] == 'B') ||
            (RF_rType[response[r]] == 'I') ||
            (RF_rType[response[r]] == 'C')) {
          parentClassProp[r] = uivector(1, RF_classLevelSize[RF_rFactorMap[response[r]]]);
          leftClassProp[r]   = uivector(1, RF_classLevelSize[RF_rFactorMap[response[r]]]);
          rghtClassProp[r]   = uivector(1, RF_classLevelSize[RF_rFactorMap[response[r]]]);
        }
        else {
        }
      }  
      for (r = 1; r <= ySize; r++) {
        if (impurity[r]) {
          if ((RF_rType[response[r]] == 'B') ||
              (RF_rType[response[r]] == 'I') ||
              (RF_rType[response[r]] == 'C')) {
            for (p = 1; p <= RF_classLevelSize[RF_rFactorMap[response[r]]]; p++) {
              parentClassProp[r][p] = 0;
            }
            for (j = 1; j <= repMembrSize; j++) {
              parentClassProp[r][RF_classLevelIndex[RF_rFactorMap[response[r]]][(uint) RF_response[treeID][response[r]][ repMembrIndx[j] ]]] ++;
            }
          }
          else {
            sumRghtSave[r] = 0.0;
            sumRghtSqrSave[r] = 0.0;
            switch(RF_splitRule) {
            case MV_WT_OFF:
            case MV_WT_HVY:
              for (j = 1; j <= repMembrSize; j++) {
                sumRghtSave[r] += RF_response[treeID][response[r]][ repMembrIndx[j] ];
                sumRghtSqrSave[r] += pow(RF_response[treeID][response[r]][ repMembrIndx[j] ], 2.0);
              }
              break;
            default:
              for (j = 1; j <= repMembrSize; j++) {
                sumRghtSave[r] += RF_response[treeID][response[r]][ repMembrIndx[j] ] - mean[r];
              }
              break;
            }      
          }
        }  
      }  
      double partialLeft, partialRght;
      double partialLeftSqr, partialRghtSqr;
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
          if (TRUE) {
            leftSize = 0;
            priorMembrIter = 0;
            if (factorFlag == FALSE) {
              for (j = 1; j <= repMembrSize; j++) {
                localSplitIndicator[ j ] = RIGHT;
              }
              for (r = 1; r <= ySize; r++) {
                if (impurity[r]) {
                  if ((RF_rType[response[r]] == 'B') ||
                      (RF_rType[response[r]] == 'I') ||
                      (RF_rType[response[r]] == 'C')) {
                    for (p = 1; p <= RF_classLevelSize[RF_rFactorMap[response[r]]]; p++) {
                      rghtClassProp[r][p] = parentClassProp[r][p];
                      leftClassProp[r][p] = 0;
                    }
                  }
                  else {
                    sumRght[r] = sumRghtSave[r];
                    sumRghtSqr[r] = sumRghtSqrSave[r];
                    sumLeft[r] = 0.0;
                    sumLeftSqr[r] = 0.0;
                  }
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
              rghtSize = repMembrSize - leftSize;
              if ((leftSize != 0) && (rghtSize != 0)) {
                delta     = 0.0;
                deltaNorm = 0;
                for (r = 1; r <= ySize; r++) {
                  if (impurity[r]) {
                    if (factorFlag == TRUE) {
                      if ((RF_rType[response[r]] == 'B') ||
                          (RF_rType[response[r]] == 'I') ||
                          (RF_rType[response[r]] == 'C')) {
                        for (p=1; p <= RF_classLevelSize[RF_rFactorMap[response[r]]]; p++) {
                          leftClassProp[r][p] = 0;
                        }
                        for (k = 1; k <= repMembrSize; k++) {
                          if (localSplitIndicator[ indxx[k] ] == LEFT) {
                            leftClassProp[r][RF_classLevelIndex[RF_rFactorMap[response[r]]][(uint) RF_response[treeID][response[r]][ repMembrIndx[ indxx[k] ] ]]] ++;
                          }
                          else {
                          }
                        }
                        for (p = 1; p <= RF_classLevelSize[RF_rFactorMap[response[r]]]; p++) {
                          rghtClassProp[r][p] = parentClassProp[r][p] - leftClassProp[r][p];
                        }
                      }
                      else {
                        switch (RF_splitRule) {
                        case MV_WT_OFF:
                        case MV_WT_HVY:
                          sumLeft[r] = sumRght[r] = 0.0;
                          sumLeftSqr[r] = sumRghtSqr[r] = 0.0;
                          for (k = 1; k <= repMembrSize; k++) {
                            if (localSplitIndicator[ indxx[k] ] == LEFT) {
                              sumLeft[r] += RF_response[treeID][response[r]][ repMembrIndx[ indxx[k]] ];
                              sumLeftSqr[r] += pow(RF_response[treeID][response[r]][ repMembrIndx[ indxx[k] ]], 2.0);
                            }
                            else {
                              sumRght[r] += RF_response[treeID][response[r]][ repMembrIndx[ indxx[k]] ];
                              sumRghtSqr[r] += pow(RF_response[treeID][response[r]][ repMembrIndx[ indxx[k] ]], 2.0);
                            }
                          }
                          break;
                        default:
                          sumLeft[r] = sumRght[r] = 0.0;
                          for (k = 1; k <= repMembrSize; k++) {
                            if (localSplitIndicator[ indxx[k] ] == LEFT) {
                              sumLeft[r] += RF_response[treeID][response[r]][ repMembrIndx[ indxx[k] ]] - mean[r];
                            }
                            else {
                              sumRght[r] += RF_response[treeID][response[r]][ repMembrIndx[ indxx[k] ]] - mean[r];
                            }
                          }
                          break;
                        }
                      }
                    }
                    else {
                      for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                        if ((RF_rType[response[r]] == 'B') ||
                            (RF_rType[response[r]] == 'I') ||
                            (RF_rType[response[r]] == 'C')) {
                          leftClassProp[r][RF_classLevelIndex[RF_rFactorMap[response[r]]][(uint) RF_response[treeID][response[r]][ repMembrIndx[indxx[k]] ]]] ++;
                          rghtClassProp[r][RF_classLevelIndex[RF_rFactorMap[response[r]]][(uint) RF_response[treeID][response[r]][ repMembrIndx[indxx[k]] ]]] --;
                        }
                        else {
                          switch(RF_splitRule) {
                          case MV_WT_OFF:
                          case MV_WT_HVY:
                            sumLeft[r] += RF_response[treeID][response[r]][ repMembrIndx[indxx[k]] ];
                            sumLeftSqr[r] += pow(RF_response[treeID][response[r]][ repMembrIndx[indxx[k]] ], 2.0);
                            sumRght[r] -= RF_response[treeID][response[r]][ repMembrIndx[indxx[k]] ];
                            sumRghtSqr[r] -= pow(RF_response[treeID][response[r]][ repMembrIndx[indxx[k]] ], 2.0);                              
                            break;
                          default:
                            sumLeft[r] += RF_response[treeID][response[r]][ repMembrIndx[indxx[k]] ] - mean[r];
                            sumRght[r] -= RF_response[treeID][response[r]][ repMembrIndx[indxx[k]] ] - mean[r];
                            break;
                          }
                        }
                      }
                    }  
                    if ((leftSize > 0) && (rghtSize > 0)) {
                      deltaNorm ++;
                      if ((RF_rType[response[r]] == 'B') ||
                          (RF_rType[response[r]] == 'I') ||
                          (RF_rType[response[r]] == 'C')) {
                        partialLeft = partialRght = 0.0;
                        switch(RF_splitRule) {
                        case MV_WT_OFF:
                          for (p = 1; p <= RF_classLevelSize[RF_rFactorMap[response[r]]]; p++) {
                            partialLeft += pow((double) leftClassProp[r][p] / leftSize, 2.0);
                            partialRght += pow((double) (rghtClassProp[r][p]) / rghtSize, 2.0);
                          }
                          delta += partialLeft + partialRght;
                          break;
                        case MV_WT_HVY:
                          for (p = 1; p <= RF_classLevelSize[RF_rFactorMap[response[r]]]; p++) {
                            partialLeft += (double) upower(leftClassProp[r][p], 2);
                            partialRght += (double) upower(rghtClassProp[r][p], 2);
                          }
                          delta +=
                            (partialLeft / (double) (upower(repMembrSize, 2))) +
                            (partialRght / (double) (upower(repMembrSize, 2))) -
                            pow((double) leftSize / repMembrSize, 2.0) -
                            pow((double) rghtSize / repMembrSize, 2.0) + 2.0;
                          break;
                        default:
                          for (p = 1; p <= RF_classLevelSize[RF_rFactorMap[response[r]]]; p++) {
                            partialLeft += (double) upower(leftClassProp[r][p], 2);
                            partialRght += (double) upower(rghtClassProp[r][p], 2);
                          }
                          delta += ((partialLeft / leftSize) + (partialRght / rghtSize)) / repMembrSize;
                          break;
                        }
                      }
                      else {
                        switch(RF_splitRule) {
                        case MV_WT_OFF:
                          partialLeft = pow (sumLeft[r], 2.0) / (double) upower (leftSize, 2);
                          partialRght = pow (sumRght[r], 2.0) / (double) upower (rghtSize, 2);
                          partialLeftSqr = sumLeftSqr[r] / leftSize;
                          partialRghtSqr = sumRghtSqr[r] / rghtSize;
                          delta += partialLeft + partialRght - partialLeftSqr - partialRghtSqr;
                          break;
                        case MV_WT_HVY:
                          partialLeft = pow (sumLeft[r], 2.0) / (double) upower (repMembrSize, 2);
                          partialRght = pow (sumRght[r], 2.0) / (double) upower (repMembrSize, 2);
                          partialLeftSqr = sumLeftSqr[r] * leftSize / (double) upower (repMembrSize, 2);
                          partialRghtSqr = sumRghtSqr[r] * rghtSize / (double) upower (repMembrSize, 2);
                          delta += partialLeft + partialRght - partialLeftSqr - partialRghtSqr;
                          break;
                        default:
                          partialLeft = pow (sumLeft[r], 2.0) / (leftSize * variance[r]);
                          partialRght = pow (sumRght[r], 2.0) / (rghtSize * variance[r]);
                          delta += partialLeft + partialRght;
                          break;
                        }
                      }
                    }  
                  }  
                }  
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
      for (r = 1; r <= ySize; r++) {
        if ((RF_rType[response[r]] == 'B') ||
            (RF_rType[response[r]] == 'I') ||
            (RF_rType[response[r]] == 'C')) {
          free_uivector (parentClassProp[r], 1, RF_classLevelSize[RF_rFactorMap[response[r]]]);
          free_uivector (leftClassProp[r], 1, RF_classLevelSize[RF_rFactorMap[response[r]]]);
          free_uivector (rghtClassProp[r], 1, RF_classLevelSize[RF_rFactorMap[response[r]]]);
        }
        else {
        }
      }
      free_new_vvector(parentClassProp, 1, ySize, NRUTIL_UPTR);
      free_new_vvector(leftClassProp,   1, ySize, NRUTIL_UPTR);
      free_new_vvector(rghtClassProp,   1, ySize, NRUTIL_UPTR);
      free_dvector(sumLeft,     1, ySize);
      free_dvector(sumRght,     1, ySize);
      free_dvector(sumRghtSave, 1, ySize);
      free_dvector(sumLeftSqr,     1, ySize);
      free_dvector(sumRghtSqr,     1, ySize);
      free_dvector(sumRghtSqrSave, 1, ySize);
      unstackRandomCovariates(treeID, distributionObj);
      unstackSplitPreliminary(repMembrSize,
                              localSplitIndicator,
                              splitVector);
    }  
    free_cvector(impurity,   1, ySize);
    free_dvector(mean,       1, ySize);
    free_dvector(variance,   1, ySize);
    free_uivector(response, 1, ySize);
  }  
  unstackPreSplit(preliminaryResult,
                  parent,
                  multImpFlag,
                  TRUE);  
  result = summarizeSplitResult(splitInfoMax);
  return result;
}
DistributionObj *stackRandomResponsesSimple(uint treeID, Node *parent) {
  DistributionObj *obj = makeDistributionObjRaw();
  obj -> indexSize = RF_ySize;
  obj -> index = NULL;
  obj -> uIndexAllocSize = 0;
  if (RF_ytry > 1) {
    if (RF_ytry < obj -> indexSize) {
      obj -> uIndexAllocSize = obj -> indexSize;
      obj -> index = uivector(1, obj -> uIndexAllocSize);
      for (uint p = 1; p <= obj -> uIndexAllocSize; p++) {
        (obj -> index)[p] = p;
      }
    }
  }
  return obj;
}
void unstackRandomResponsesSimple(uint treeID, DistributionObj *obj) {
  if (obj -> uIndexAllocSize > 0) {
    free_uivector(obj -> index, 1, obj -> uIndexAllocSize);
  }
  freeDistributionObjRaw(obj);
}
char selectRandomResponsesSimpleVector(uint  treeID,
                                        Node *parent,
                                       DistributionObj *distributionObj,
                                       uint *response,
                                       uint *responseCount) {
  uint r;
  char yVarSearch;
  yVarSearch = TRUE;
  *responseCount = 0;
  while (yVarSearch) {
    if (distributionObj -> indexSize > 0) {
      if (RF_ytry == 1) {
        (*responseCount) ++;
        distributionObj -> slot = (uint) ceil(ran1B(treeID) * ((distributionObj -> indexSize) * 1.0));
        response[*responseCount] = distributionObj -> slot;
        yVarSearch = FALSE;
      }
      else {
        if (RF_ytry >= RF_ySize) {
          for (r = 1; r <= RF_ySize; r++) {
            (*responseCount) ++;
            response[*responseCount] = r;
          }
          yVarSearch = FALSE;
        }
        else {
          (*responseCount) ++;
          distributionObj -> slot = (uint) ceil(ran1B(treeID) * ((distributionObj -> indexSize) * 1.0));
          response[*responseCount] = (distributionObj -> index)[distributionObj -> slot];
          (distributionObj -> index)[distributionObj -> slot] = (distributionObj -> index)[distributionObj -> indexSize];
          (distributionObj -> indexSize) --;
        }
      }
      if ((*responseCount) >= RF_ytry) {
        yVarSearch = FALSE;
      }
    }
    else {
      yVarSearch = FALSE;
    }
  }  
  return TRUE;
}
DistributionObj *stackRandomResponsesGeneric(uint treeID, Node *parent) {
  uint actualWeightType;
  DistributionObj *obj = makeDistributionObjRaw();
  actualWeightType = RF_yWeightType;
  obj -> permissibleIndex    = NULL;
  obj -> permissibleSize     = RF_ySize;
  obj -> permissible         = cvector(1, obj -> permissibleSize);
  for (uint r = 1; r <= RF_ySize; r++) {
    obj -> permissible[r] = TRUE;
  }
  obj -> weightType          = actualWeightType;
  obj -> weight              = RF_yWeight;
  obj -> weightSorted        = RF_yWeightSorted;
  obj -> densityAllocSize    = RF_yWeightDensitySize;
  initializeCDFNew(treeID, obj);
  return obj;
}
void unstackRandomResponsesGeneric(uint treeID, DistributionObj *obj) {
  free_cvector(obj -> permissible, 1, obj -> permissibleSize);
  discardCDFNew(treeID, obj);
  freeDistributionObjRaw(obj);
}
char selectRandomResponsesGenericVector(uint     treeID,
                                        Node     *parent,
                                        DistributionObj *distributionObj,
                                        uint     *response,
                                        uint     *responseCount) {
  uint r;
  uint selectedValue;
  char yVarSearch;
  *responseCount = 0;
  if (RF_ytry == 1) {
    selectedValue = sampleFromCDFNew(ran1B, treeID, distributionObj);
    if (selectedValue != 0) {
      (*responseCount) ++;
      response[*responseCount] = selectedValue;
    }    
  }
  else if (RF_ytry >= RF_ySize) {
    for (r = 1; r <= RF_ySize; r++) {
      (*responseCount) ++;
      response[*responseCount] = r;
    }
  }
  else {
    yVarSearch = TRUE;
    while (yVarSearch) {
      selectedValue = sampleFromCDFNew(ran1B, treeID, distributionObj);
      if (selectedValue != 0) {
        (*responseCount) ++;
        response[*responseCount] = selectedValue;
        updateCDFNew(treeID, distributionObj);
        if ((*responseCount) >= RF_ytry) {
          yVarSearch = FALSE;
        }
      }
      else {
        yVarSearch = FALSE;
      }
    }  
  }
  return TRUE;
}
