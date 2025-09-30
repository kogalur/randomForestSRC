
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "split.h"
#include "splitCustom.h"
#include "nrutil.h"
#include "splitSurv.h"
#include "splitClas.h"
#include "splitQuantile.h"
#include "splitCustomDriver.h"
#include "splitMahalanobis.h"
#include "splitGreedy.h"
#include "splitUtil.h"
#include "error.h"
SplitRuleObj *makeSplitRuleObj(uint rule) {
  SplitRuleObj *obj = (SplitRuleObj*) gblock((size_t) sizeof(SplitRuleObj));
  switch(rule) {
  case SURV_LGRNK:
    obj -> function = &logRankNCR;
    break;
  case SURV_LRSCR:
    obj -> function = &logRankNCR;
    break;
  case SURV_CR_LAU:
    obj -> function = &logRankCR;
    break;
  case SURV_CR_GEN:
    obj -> function = &logRankCR;
    break;
  case SURV_BSG1:
    obj -> function = &brierScoreGradient1;
    break;
  case RAND_SPLIT:
    obj -> function = randomSplit;
    break;
  case REGR_NRM:
    obj -> function = regressionXwghtSplit;
    break;
  case REGR_WT_OFF:
    obj -> function = regressionXwghtSplit;
    break;
  case REGR_WT_HVY:
    obj -> function = regressionXwghtSplit;
    break;
  case REGR_QUANT:
    obj -> function = &quantileRegrSplit;
    break;
  case LARG_QUANT:
    obj -> function = &locallyAdaptiveQuantileRegrSplit;
    break;
  case CLAS_NRM:
    obj -> function = classificationXwghtSplit;
    break;
  case CLAS_WT_OFF:
    obj -> function = classificationXwghtSplit;
    break;
  case CLAS_WT_HVY:
    obj -> function = classificationXwghtSplit;
    break;
  case CLAS_AU_ROC:
    obj -> function = &classificationAreaUnderROCSplit;
    break;
  case CLAS_ENTROP:
    obj -> function = &classificationEntropySplit;
    break;
  case MV_NRM:
    obj -> function = multivariateSplit;
    break;
  case MV_WT_OFF:
    obj -> function = multivariateSplit;
    break;
  case MV_WT_HVY:
    obj -> function = multivariateSplit;
    break;
  case USPV_NRM:
    obj -> function = unsupervisedSplit;
    break;
  case USPV_WT_OFF:
    obj -> function = unsupervisedSplit;
    break;
  case USPV_WT_HVY:
    obj -> function = unsupervisedSplit;
    break;
  case CUST_SPLIT:
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      if (!(RF_opt & OPT_COMP_RISK)) {
        obj -> function = &customSurvivalSplit;
      }
      else {
        obj -> function = &customCompetingRiskSplit;
      }
    }
    else {
      obj -> function = &customMultivariateSplit;
    }
    break;
  case MAHALANOBIS:
    obj -> function = &mahalanobis;
    break;
  default:
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Split rule not found:  %10d", rule);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
    break;
  }
  return obj;
}
void freeSplitRuleObj(SplitRuleObj *obj) {
  free_gblock(obj, (size_t) sizeof(SplitRuleObj));
}
char getBestSplit(uint       treeID,
                  Node      *parent,
                  uint       splitRule,
                  SplitInfoMax *splitInfoMax,
                  char       multImpFlag) {
  char  result;
  result = RF_splitRuleObj -> function(treeID,
                                       parent,
                                       splitInfoMax,
                                       NULL,  
                                       multImpFlag);
  return result;
}
char randomSplitGeneric(uint       treeID,
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
  char multVarFlag;
  double delta;
  uint j;
  mwcpSizeAbsolute       = 0;     
  multVarFlag = TRUE;
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    multVarFlag = FALSE;
  }
  else {
    if (((RF_rFactorCount == 0) && (RF_rNonFactorCount == 1)) ||
        ((RF_rFactorCount == 1) && (RF_rNonFactorCount == 0))) {
      multVarFlag = FALSE;
    }
  }
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        multImpFlag,
                                        multVarFlag);
  if(preliminaryResult) {
    uint  repMembrSize = parent -> repMembrSize;
    uint  nonMissMembrSize;
    uint *nonMissMembrIndx;
    stackSplitPreliminary(repMembrSize,
                          & localSplitIndicator,
                          & splitVector);
    DistributionObj *distributionObj = stackRandomCovariates(treeID, parent);
    covariateCount = 0;
    while ( (RF_nativeIsNaN(splitInfoMax -> deltaMax)) &&
            selectRandomCovariates(treeID,
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
          for (j = 1; j <= nonMissMembrSize; j++) {
            localSplitIndicator[ nonMissMembrIndx[indxx[j]] ] = RIGHT;
          }
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
          if ((leftSize != 0) && (rghtSize != 0)) {
            delta = 0;
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
          j = splitLength;
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
    unstackRandomCovariates(treeID, distributionObj);
    unstackSplitPreliminary(repMembrSize,
                            localSplitIndicator,
                            splitVector);
  }  
  unstackPreSplit(preliminaryResult,
                  parent,
                  multImpFlag,
                  multVarFlag);
  result = summarizeSplitResult(splitInfoMax);
  return result;
}
char randomSplitSimple(uint       treeID,
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
  char multVarFlag;
  double delta;
  uint j;
  mwcpSizeAbsolute       = 0;     
  multVarFlag = TRUE;
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    multVarFlag = FALSE;
  }
  else {
    if (((RF_rFactorCount == 0) && (RF_rNonFactorCount == 1)) ||
        ((RF_rFactorCount == 1) && (RF_rNonFactorCount == 0))) {
      multVarFlag = FALSE;
    }
  }
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        multImpFlag,
                                        multVarFlag);
  if(preliminaryResult) {
    uint  repMembrSize = parent -> repMembrSize;
    uint  nonMissMembrSize;
    stackSplitPreliminary(repMembrSize,
                          & localSplitIndicator,
                          & splitVector);
    DistributionObj *distributionObj = stackRandomCovariates(treeID, parent);
    double deltaMax;
    uint   indexMax;
    covariateCount = 0;
    while ( (RF_nativeIsNaN(splitInfoMax -> deltaMax)) &&
            selectRandomCovariates(treeID,
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
        nonMissMembrSize = parent -> nonMissMembrSize;
        observation = RF_observation[treeID][covariate];
        deltaMax = RF_nativeNaN;
        indexMax =  0;
        for (j = 1; j < splitLength; j++) {
          priorMembrIter = 0;
          leftSize = 0;
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
            delta = 0;
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
          j = splitLength;
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
                  multVarFlag);
  result = summarizeSplitResult(splitInfoMax);
  return result;
}
void registerThis (customFunction func, unsigned int family, unsigned int slot) {
  if ((slot >= 1) && (slot <= 16)) {
    customFunctionArray[family][slot-1] = func;
  }
  else {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Invalid slot for custom split rule:  %10d", slot);
    RF_nativeError("\nRF-SRC:  The slot must be an integer within [1, 16].");
    RF_nativeExit();
  }    
}
