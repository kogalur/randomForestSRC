
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "splitSurv.h"
#include "splitUtil.h"
#include "splitUtilSurv.h"
#include "nrutil.h"
#include "error.h"
char logRankNCR (uint       treeID,
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
  uint j, k, m;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nlogRankNCR() ENTRY ...\n");
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
    uint   *survivalTimeIndexRank;
    double *survivalRank;
    double  meanSurvRank, varSurvRank;
    double deltaNum, deltaNumAdj, deltaDen;
    uint   tIndx;
    meanSurvRank = varSurvRank = 0;  
    survivalTimeIndexRank = NULL;    
    survivalRank = NULL;             
    localEventTimeSize = 0;          
    delta = deltaNum = 0.0;          
    switch(RF_splitRule) {
    case SURV_LGRNK:
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
      break;
    case SURV_LRSCR:
      survivalTimeIndexRank = uivector(1, repMembrSize);
      survivalRank = dvector(1, repMembrSize);
      localEventTimeSize = 1;
      break;
    default:
      break;
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
        switch(RF_splitRule) {
        case SURV_LGRNK:
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
          break;
        case SURV_LRSCR:
          localEventTimeSize = 1;
          break;
        default:
          break;
        }
        if (localEventTimeSize > 0) {
          leftSize = 0;
          priorMembrIter = 0;
          if (factorFlag == FALSE) {
            for (j = 1; j <= nonMissMembrSize; j++) {
              localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
            }
            switch(RF_splitRule) {
            case SURV_LGRNK:
              for (m = 1; m <= localEventTimeSize; m++) {
                nodeLeftEvent[m] = nodeLeftAtRisk[m] = 0;
              }
              break;
            case SURV_LRSCR:
              deltaNum =  0.0;
              break;
            default:
              break;
            }
          }
          switch(RF_splitRule) {
          case SURV_LGRNK:
            break;
          case SURV_LRSCR:
            for (k = 1; k <= nonMissMembrSize; k++) {
              survivalTimeIndexRank[k] = 0;
              for (j = 1; j <= nonMissMembrSize; j++) {
                if ( RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[j]] ]  <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] ) {
                  survivalTimeIndexRank[k] ++;
                }
              }
            }
            ${trace.token}      if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
            ${trace.token}        RF_nativePrint("\nLocal Membership Information for Parent Node: \n");
            ${trace.token}        RF_nativePrint("       RANK    NMISidx       REPLidx    SPLTval    TIMEidx -> SORTidx \n");
            ${trace.token}        for (k = 1; k <=  nonMissMembrSize; k++) {
            ${trace.token}          RF_nativePrint(" %10d %10d %10d %10.4f %10d %10d\n", k,
            ${trace.token}                  nonMissMembrIndx[indxx[k]], repMembrIndx[nonMissMembrIndx[indxx[k]]], RF_observation[treeID][covariate][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ],
            ${trace.token}                  RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ], survivalTimeIndexRank[k]);
            ${trace.token}        }
            ${trace.token}      }
            meanSurvRank = varSurvRank = 0;
            for (k = 1; k <= nonMissMembrSize; k++) {
              survivalRank[k] = 0;
              for (j = 1; j <= survivalTimeIndexRank[k]; j++) {
                survivalRank[k] = survivalRank[k] + (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[j]]] ] / (nonMissMembrSize - survivalTimeIndexRank[j] + 1) );
              }
              survivalRank[k] = RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - survivalRank[k];
              meanSurvRank = meanSurvRank + survivalRank[k];
              varSurvRank = varSurvRank +  pow(survivalRank[k], 2.0);
            }
            varSurvRank = ( varSurvRank - (pow(meanSurvRank, 2.0) / nonMissMembrSize) ) / (nonMissMembrSize - 1);
            meanSurvRank = meanSurvRank / nonMissMembrSize;
            break;
          default:
            break;
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
            ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
            ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
            ${trace.token}              RF_nativePrint("\n PriorIter:     %10d  CurrentIter:   %10d", priorMembrIter, currentMembrIter);
            ${trace.token}            }
            ${trace.token}          }
            if ((leftSize != 0) && (rghtSize != 0)) {
              if (factorFlag == TRUE) {
                switch(RF_splitRule) {
                case SURV_LGRNK:
                  for (m = 1; m <= localEventTimeSize; m++) {
                    nodeLeftEvent[m] = nodeLeftAtRisk[m] = 0;
                  }
                  for (k = 1; k <= nonMissMembrSize; k++) {
                    if (localSplitIndicator[  nonMissMembrIndx[k]  ] == LEFT) {
                      tIndx = 0;  
                      for (m = 1; m <= localEventTimeSize; m++) {
                        if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[k]] ]) {
                          ${trace.token}            if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                          ${trace.token}              if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
                          ${trace.token}                RF_nativePrint("\nNLAR:  %10d %10d %10d ", m, localEventTimeIndex[m], RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[k]] ]);
                          ${trace.token}              }
                          ${trace.token}            }
                          tIndx = m;
                          nodeLeftAtRisk[tIndx] ++;
                        }
                        else {
                          m = localEventTimeSize;
                        }
                      }
                      if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[k]] ] > 0) {
                        nodeLeftEvent[tIndx] ++;
                      }
                    }
                    else {
                    }
                  } 
                  break;
                case SURV_LRSCR:
                  deltaNum = 0.0;
                  for (k = 1; k <= nonMissMembrSize; k++) {
                    if (localSplitIndicator[ nonMissMembrIndx[k] ] == LEFT) {
                      deltaNum = deltaNum + survivalRank[k];
                    }
                  }
                  break;
                default:
                  break;
                }
              }
              else {
                switch(RF_splitRule) {
                case SURV_LGRNK:
                  for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                    tIndx = 0;  
                    for (m = 1; m <= localEventTimeSize; m++) {
                      if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]) {
                        ${trace.token}            if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                        ${trace.token}              if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
                        ${trace.token}                RF_nativePrint("\nNLAR:  %10d %10d %10d ", m, localEventTimeIndex[m], RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]);
                        ${trace.token}              }
                        ${trace.token}            }
                        tIndx = m;
                        nodeLeftAtRisk[tIndx] ++;
                      }
                      else {
                        m = localEventTimeSize;
                      }
                    }
                    if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
                      nodeLeftEvent[tIndx] ++;
                    }
                  }
                  break;
                case SURV_LRSCR:
                  for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                    deltaNum = deltaNum + survivalRank[ indxx[k] ];
                  }
                  break;
                default:
                  break;
                }
              }
              switch(RF_splitRule) {
              case SURV_LGRNK:
                ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                ${trace.token}  if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
                ${trace.token}    if (localEventTimeSize > 0) {
                ${trace.token}      RF_nativePrint("\nRunning Split Risk Counts: \n");
                ${trace.token}      RF_nativePrint("     timIdx    PARrisk    LFTrisk   PARevent    LFTevent \n");
                ${trace.token}      for (k=1; k <=  localEventTimeSize; k++) {
                ${trace.token}        nodeRightAtRisk[k] = nodeParentAtRisk[k] - nodeLeftAtRisk[k],
                ${trace.token}        RF_nativePrint(" %10d %10d %10d %10d %10d\n", k,
                ${trace.token}                nodeParentAtRisk[k], nodeLeftAtRisk[k],
                ${trace.token}                nodeParentEvent[k], nodeLeftEvent[k]);
                ${trace.token}      }
                ${trace.token}      uint totalEventCount = 0;
                ${trace.token}      uint leftEventCount  = 0;
                ${trace.token}      uint rightEventCount = 0;
                ${trace.token}      for (k=1; k <=  localEventTimeSize; k++) {
                ${trace.token}        nodeRightEvent[k] = nodeParentEvent[k] - nodeLeftEvent[k],
                ${trace.token}        totalEventCount += nodeParentEvent[k];
                ${trace.token}        leftEventCount  += nodeLeftEvent[k];
                ${trace.token}        rightEventCount += nodeRightEvent[k];
                ${trace.token}      }
                ${trace.token}      RF_nativePrint("\nRunning Split Total LFT & RGT Events: \n");
                ${trace.token}      RF_nativePrint("                                             %10d %10d %10d \n", totalEventCount, leftEventCount, rightEventCount);
                ${trace.token}    }
                ${trace.token}  }
                ${trace.token}  }
                delta = deltaNum = deltaDen =  0.0;
                for (k = 1; k <= localEventTimeSize; k++) {
                  deltaNum = deltaNum + ((double) nodeLeftEvent[k] - ((double) ( nodeLeftAtRisk[k] * nodeParentEvent[k]) / nodeParentAtRisk[k]));
                  ${trace.token}            if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
                  ${trace.token}              if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
                  ${trace.token}                RF_nativePrint("\nPartial Sum deltaNum:  %10d %10.4f", k, deltaNum);
                  ${trace.token}              }
                  ${trace.token}            }
                  if (nodeParentAtRisk[k] >= 2) {
                    deltaDen = deltaDen + (
                                           ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k]) *
                                           (1.0 - ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k])) *
                                           ((double) (nodeParentAtRisk[k] - nodeParentEvent[k]) / (nodeParentAtRisk[k] - 1)) * nodeParentEvent[k]
                                           );
                    ${trace.token}            if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
                    ${trace.token}              if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
                    ${trace.token}                RF_nativePrint("\nPartial Sum deltaDen:  %10d %10.4f", k, deltaDen);
                    ${trace.token}              }
                    ${trace.token}            }
                  }
                }
                deltaNum = fabs(deltaNum);
                deltaDen = sqrt(deltaDen);
                if (deltaDen <= EPSILON) {
                  if (deltaNum <= EPSILON) {
                    delta = 0.0;
                  }
                  else {
                    delta = deltaNum / deltaDen;
                  }
                }
                else {
                  delta = deltaNum / deltaDen;
                }
                break;
              case SURV_LRSCR:
                deltaNumAdj  = deltaNum - (leftSize * meanSurvRank);
                deltaDen     = leftSize * (1.0 - (leftSize / nonMissMembrSize)) * varSurvRank;
                deltaNumAdj = fabs(deltaNumAdj);
                deltaDen = sqrt(deltaDen);
                if (deltaDen <= EPSILON) {
                  if (deltaNumAdj <= EPSILON) {
                    delta = 0.0;
                  }
                  else {
                    delta = deltaNumAdj / deltaDen;
                  }
                }
                else {
                  delta = deltaNumAdj / deltaDen;
                }
                break;
              default:
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
        }  
        else {
          ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
          ${trace.token}  if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
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
        switch(RF_splitRule) {
        case SURV_LGRNK:
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
          break;
        case SURV_LRSCR:
          break;
        default:
          break;
        }
      }
    }  
    switch(RF_splitRule) {
    case SURV_LGRNK:
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
    break;
    case SURV_LRSCR:
      free_uivector(survivalTimeIndexRank, 1, repMembrSize);
      free_dvector(survivalRank, 1, repMembrSize);
      break;
    default:
      break;
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
  ${trace.token}    RF_nativePrint("\nlogRankNCR(%10d) result:  %10d", treeID, result);
  ${trace.token}    RF_nativePrint("\nlogRankNCR(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
  return result;
}
char logRankCR (uint       treeID,
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
  uint j, k, m;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nlogRank() ENTRY ...\n");
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
    uint **nodeParentEventCR, **nodeLeftEventCR;
    uint **nodeParentInclusiveAtRisk, **nodeLeftInclusiveAtRisk;
    nodeParentInclusiveAtRisk = nodeLeftInclusiveAtRisk = NULL;  
    nodeParentEventCR = nodeLeftEventCR = NULL;  
    double deltaNum, deltaSubNum, deltaDen, deltaSubDen;
    uint   tIndx;
    uint   q, s, r;
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
      nodeParentEventCR = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
      nodeLeftEventCR = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
      if (RF_splitRule == SURV_CR_LAU) {
        nodeParentInclusiveAtRisk = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
        nodeLeftInclusiveAtRisk = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
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
          if (localEventTimeSize > 0) {
            nodeParentEventCR = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
            nodeLeftEventCR = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
          }
          if (RF_splitRule == SURV_CR_LAU) {
            if (localEventTimeSize > 0) {
              nodeParentInclusiveAtRisk = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
              nodeLeftInclusiveAtRisk = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
            }
          }
        }
        if (localEventTimeSize > 0) {
          leftSize = 0;
          priorMembrIter = 0;
          for (m = 1; m <= localEventTimeSize; m++) {
            for (q = 1; q <= RF_eventTypeSize; q++) {
              nodeParentEventCR[q][m] = 0;
            }
          }
          for (k = 1; k <= nonMissMembrSize; k++) {
            if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[k]] ] > 0) {
              tIndx = 0;  
              for (m = 1; m <= localEventTimeSize; m++) {
                if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[k]] ]) {
                  ${trace.token}            if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                  ${trace.token}              if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
                  ${trace.token}                RF_nativePrint("\nNLAR:  %10d %10d %10d ", m, localEventTimeIndex[m], RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[k]] ]);
                  ${trace.token}              }
                  ${trace.token}            }
                  tIndx = m;
                }
                else {
                  m = localEventTimeSize;
                }
              }
              nodeParentEventCR[RF_eventTypeIndex[(uint) RF_status[treeID][ repMembrIndx[nonMissMembrIndx[k]] ]]][tIndx] ++;
            }
          }
          if (RF_splitRule == SURV_CR_LAU) {
            for (m = 1; m <= localEventTimeSize; m++) {
              for (q = 1; q <= RF_eventTypeSize; q++) {
                if (RF_crWeight[q] > 0) {
                  nodeParentInclusiveAtRisk[q][m] = nodeParentAtRisk[m];
                  for (s = 1; s < m; s++) {
                    for (r = 1; r <= RF_eventTypeSize; r++) {
                      if (q != r) {
                        nodeParentInclusiveAtRisk[q][m]  += nodeParentEventCR[r][s];
                      }
                    }
                  }
                }
              }
            }
          }
          if (factorFlag == FALSE) {
            for (j = 1; j <= nonMissMembrSize; j++) {
              localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
            }
            for (m = 1; m <= localEventTimeSize; m++) {
              nodeLeftAtRisk[m] = 0;
              for (q = 1; q <= RF_eventTypeSize; q++) {
                nodeLeftEventCR[q][m] = 0;
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
            ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
            ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
            ${trace.token}              RF_nativePrint("\n PriorIter:     %10d  CurrentIter:   %10d", priorMembrIter, currentMembrIter);
            ${trace.token}            }
            ${trace.token}          }
            if ((leftSize != 0) && (rghtSize != 0)) {
              if (factorFlag == TRUE) {
                for (m = 1; m <= localEventTimeSize; m++) {
                  nodeLeftAtRisk[m] = 0;
                  for (q = 1; q <= RF_eventTypeSize; q++) {
                    nodeLeftEventCR[q][m] = 0;
                  }
                }
                for (k = 1; k <= nonMissMembrSize; k++) {
                  if (localSplitIndicator[  nonMissMembrIndx[k]  ] == LEFT) {
                    tIndx = 0;  
                    for (m = 1; m <= localEventTimeSize; m++) {
                      if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[k]] ]) {
                        ${trace.token}            if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                        ${trace.token}              if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
                        ${trace.token}                RF_nativePrint("\nNLAR:  %10d %10d %10d ", m, localEventTimeIndex[m], RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[k]] ]);
                        ${trace.token}              }
                        ${trace.token}            }
                        tIndx = m;
                        nodeLeftAtRisk[tIndx] ++;
                      }
                      else {
                        m = localEventTimeSize;
                      }
                    }
                    if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[k]] ] > 0) {
                      nodeLeftEventCR[RF_eventTypeIndex[(uint) RF_status[treeID][ repMembrIndx[nonMissMembrIndx[k]] ]]][tIndx] ++;
                    }
                  }
                }
              }
              else {
                for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                  tIndx = 0;  
                  for (m = 1; m <= localEventTimeSize; m++) {
                    if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]) {
                      ${trace.token}            if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
                      ${trace.token}              if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
                      ${trace.token}                RF_nativePrint("\nNLAR:  %10d %10d %10d ", m, localEventTimeIndex[m], RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]);
                      ${trace.token}              }
                      ${trace.token}            }
                      tIndx = m;
                      nodeLeftAtRisk[tIndx] ++;
                    }
                    else {
                      m = localEventTimeSize;
                    }
                  }
                  if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
                    nodeLeftEventCR[RF_eventTypeIndex[(uint) RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]][tIndx] ++;
                  }
                }
              }
              if (RF_splitRule == SURV_CR_LAU) {
                for (m=1; m <= localEventTimeSize; m++) {
                  for (q = 1; q <= RF_eventTypeSize; q++) {
                    if (RF_crWeight[q] > 0) {
                      nodeLeftInclusiveAtRisk[q][m] = nodeLeftAtRisk[m];
                      for (s = 1; s < m; s++) {
                        for (r = 1; r <= RF_eventTypeSize; r++) {
                          if (q != r) {
                            nodeLeftInclusiveAtRisk[q][m] += nodeLeftEventCR[r][s];
                          }
                        }
                      }
                    }
                  }
                }
              }
              ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
              ${trace.token}            if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
              ${trace.token}              RF_nativePrint("\nLogRankCR Weights:  \n");
              ${trace.token}              for (q=1; q <= RF_eventTypeSize; q++) {
              ${trace.token}                RF_nativePrint(" %10d", RF_eventType[q]);
              ${trace.token}              }
              ${trace.token}              RF_nativePrint("\n");
              ${trace.token}              for (q=1; q <= RF_eventTypeSize; q++) {
              ${trace.token}                RF_nativePrint(" %10.4f", RF_crWeight[q]);
              ${trace.token}              }
              ${trace.token}              RF_nativePrint("\n");
              ${trace.token}            }
              ${trace.token}          }
              delta = deltaNum = deltaDen =  0.0;
              if (RF_splitRule == SURV_CR_LAU) {
                for (q = 1; q <= RF_eventTypeSize; q++) {
                  if (RF_crWeight[q] > 0) {
                    deltaSubNum = 0;
                    for (m = 1; m <= localEventTimeSize; m++) {
                      deltaSubNum = deltaSubNum + (nodeLeftEventCR[q][m] - (nodeParentEventCR[q][m] * ((double) nodeLeftInclusiveAtRisk[q][m] / nodeParentInclusiveAtRisk[q][m])));
                    }
                    deltaNum = deltaNum + (RF_crWeight[q] * deltaSubNum);
                    deltaSubDen = 0;
                    for (m = 1; m <= localEventTimeSize; m++) {
                      if (nodeParentAtRisk[m] >= 2) {
                        deltaSubDen = deltaSubDen  + (
                                                      (nodeParentEventCR[q][m] * ((double) nodeLeftInclusiveAtRisk[q][m] / nodeParentInclusiveAtRisk[q][m])) *
                                                      (1.0 - ((double) nodeLeftInclusiveAtRisk[q][m] / nodeParentInclusiveAtRisk[q][m])) *
                                                      ((double) (nodeParentInclusiveAtRisk[q][m] - nodeParentEventCR[q][m]) / (nodeParentInclusiveAtRisk[q][m] - 1))
                                                      );
                      }
                    }
                    deltaDen = deltaDen + (RF_crWeight[q] * RF_crWeight[q] * deltaSubDen);
                  }
                }
              }
              else {
                for (q = 1; q <= RF_eventTypeSize; q++) {
                  if (RF_crWeight[q] > 0) {
                    deltaSubNum = 0;
                    for (m=1; m <= localEventTimeSize; m++) {
                      deltaSubNum = deltaSubNum + (nodeLeftEventCR[q][m] - (nodeParentEventCR[q][m] * ((double) nodeLeftAtRisk[m] / nodeParentAtRisk[m])));
                    }
                    deltaNum = deltaNum + (RF_crWeight[q] * deltaSubNum);
                    deltaSubDen = 0;
                    for (m = 1; m <= localEventTimeSize; m++) {
                      if (nodeParentAtRisk[m] >= 2) {
                        deltaSubDen = deltaSubDen  + (
                                                      (nodeParentEventCR[q][m] * ((double) nodeLeftAtRisk[m] / nodeParentAtRisk[m])) *
                                                      (1.0 - ((double) nodeLeftAtRisk[m] / nodeParentAtRisk[m])) *
                                                      ((double) (nodeParentAtRisk[m] - nodeParentEventCR[q][m]) / (nodeParentAtRisk[m] - 1))
                                                      );
                      }
                    }
                    deltaDen = deltaDen + (RF_crWeight[q] * RF_crWeight[q] * deltaSubDen);
                  }
                }
              }
              ${trace.token}          if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
              ${trace.token}            RF_nativePrint("\nSum deltaNum:  %16.4f", deltaNum);
              ${trace.token}            RF_nativePrint("\nSum deltaDen:  %16.4f", deltaDen);
              ${trace.token}            RF_nativePrint("\n");
              ${trace.token}          }
              deltaNum = fabs(deltaNum);
              deltaDen = sqrt(deltaDen);
              if (deltaDen <= EPSILON) {
                if (deltaNum <= EPSILON) {
                  delta = 0.0;
                }
                else {
                  delta = deltaNum / deltaDen;
                }
              }
              else {
                delta = deltaNum / deltaDen;
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
          if (localEventTimeSize > 0) {
            free_uimatrix(nodeParentEventCR, 1, RF_eventTypeSize, 1, localEventTimeSize);
            free_uimatrix(nodeLeftEventCR, 1, RF_eventTypeSize, 1, localEventTimeSize);
          }
          if (RF_splitRule == SURV_CR_LAU) {          
            if (localEventTimeSize > 0) {
              free_uimatrix(nodeParentInclusiveAtRisk, 1, RF_eventTypeSize, 1, localEventTimeSize);
              free_uimatrix(nodeLeftInclusiveAtRisk, 1, RF_eventTypeSize, 1, localEventTimeSize);
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
      free_uimatrix(nodeParentEventCR, 1, RF_eventTypeSize, 1, localEventTimeSize);
      free_uimatrix(nodeLeftEventCR, 1, RF_eventTypeSize, 1, localEventTimeSize);
      if (RF_splitRule == SURV_CR_LAU) {      
        free_uimatrix(nodeParentInclusiveAtRisk, 1, RF_eventTypeSize, 1, localEventTimeSize);
        free_uimatrix(nodeLeftInclusiveAtRisk, 1, RF_eventTypeSize, 1, localEventTimeSize);
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
  ${trace.token}    RF_nativePrint("\nlogRankCR(%1d) EXIT ...\n", result);
  ${trace.token}  }
  return result;
}
char brierScoreGradient1 (uint       treeID,
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
  ${trace.token}    RF_nativePrint("\nbrierScoreGradient1() ENTRY ...\n");
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
    uint  nonMissMembrSizeStatic = parent -> nonMissMembrSizeStatic;
    uint *nonMissMembrIndxStatic = parent -> nonMissMembrIndxStatic;
    uint  repMembrSize = parent -> repMembrSize;
    uint *repMembrIndx = parent -> repMembrIndx;
    uint  nonMissMembrSize;
    uint *nonMissMembrIndx;
    stackSplitPreliminary(repMembrSize,
                          & localSplitIndicator,
                          & splitVector);
    DistributionObj *distributionObj = stackRandomCovariates(treeID, parent);
    uint  eventTimeSize;
    uint *eventTimeCount, *eventTimeIndex;
    uint *parentEvent,  *leftEvent,  *rightEvent; 
    uint *parentAtRisk, *leftAtRisk, *rightAtRisk;
    uint  revEventTimeSize;
    uint *revEventTimeCount, *revEventTimeIndex;
    uint *revParentEvent,  *revLeftEvent,  *revRightEvent; 
    uint *revParentAtRisk, *revLeftAtRisk, *revRightAtRisk;
    double *parentSurvival, *revParentSurvival;
    double  *gHat;
    double **w_ktm;
    uint     tIndx;
    double **gamma_ktm;
    uint    *qeTimeIndex;
    uint     qeTimeSize;
    double *leftGammaSum, *rightGammaSum;
    double *leftGammaBar, *rightGammaBar;
    double  sumLeft, sumRght;
    char adHocFlag;
    leftGammaSum = rightGammaSum = leftGammaBar = rightGammaBar = NULL;
    sumRght = 0.0;      
    eventTimeSize = 0;  
    delta = 0;          
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      stackAndGetSplitSurv(treeID,
                           parent,
                           TRUE, 
                           & eventTimeCount,
                           & eventTimeIndex,
                           & eventTimeSize,
                           & parentEvent,
                           & parentAtRisk,
                           & leftEvent,
                           & leftAtRisk,
                           & rightEvent,
                           & rightAtRisk);
      stackAndGetSplitSurv(treeID,
                           parent,
                           FALSE, 
                           & revEventTimeCount,
                           & revEventTimeIndex,
                           & revEventTimeSize,
                           & revParentEvent,
                           & revParentAtRisk,
                           & revLeftEvent,
                           & revLeftAtRisk,
                           & revRightEvent,
                           & revRightAtRisk);
      stackAndGetSplitSurv2(treeID,
                            parent,
                            eventTimeSize,
                            parentEvent,
                            parentAtRisk,
                            & parentSurvival);
      stackAndGetSplitSurv2(treeID,
                            parent,
                            revEventTimeSize,
                            revParentEvent,
                            revParentAtRisk,
                            & revParentSurvival);
      stackAndGetQETime(treeID,
                        parent,
                        eventTimeIndex,
                        eventTimeSize,
                        parentSurvival,
                        & qeTimeIndex,
                        & qeTimeSize);
      stackAndGetLocalGamma(treeID,
                            parent,
                            repMembrIndx,
                            repMembrSize,
                            nonMissMembrIndxStatic,
                            nonMissMembrSizeStatic,
                            eventTimeIndex,
                            eventTimeSize,
                            revEventTimeIndex,
                            revEventTimeSize,
                            revParentSurvival,
                              qeTimeIndex,
                              qeTimeSize,
                            & gHat,
                            & w_ktm,
                            & gamma_ktm);
      leftGammaSum  = dvector(1, qeTimeSize + 1);
      rightGammaSum = dvector(1, qeTimeSize + 1);
      leftGammaBar  = dvector(1, qeTimeSize + 1);
      rightGammaBar = dvector(1, qeTimeSize + 1);
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
        if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
          stackAndGetSplitSurv(treeID,
                               parent,
                               TRUE,
                               & eventTimeCount,
                               & eventTimeIndex,
                               & eventTimeSize,
                               & parentEvent,
                               & parentAtRisk,
                               & leftEvent,
                               & leftAtRisk,        
                               & rightEvent,
                               & rightAtRisk);
          stackAndGetSplitSurv(treeID,
                               parent,
                               FALSE,
                               & revEventTimeCount,
                               & revEventTimeIndex,
                               & revEventTimeSize,
                               & revParentEvent,
                               & revParentAtRisk,
                               & revLeftEvent,      
                               & revLeftAtRisk,     
                               & revRightEvent,     
                               & revRightAtRisk);   
          stackAndGetSplitSurv2(treeID,
                                parent,
                                eventTimeSize,
                                parentEvent,
                                parentAtRisk,
                                & parentSurvival);
          stackAndGetSplitSurv2(treeID,
                                parent,
                                revEventTimeSize,
                                revParentEvent,
                                revParentAtRisk,
                                & revParentSurvival);
          stackAndGetQETime(treeID,
                            parent,
                            eventTimeIndex,
                            eventTimeSize,
                            parentSurvival,
                            & qeTimeIndex,
                            & qeTimeSize);
          stackAndGetLocalGamma(treeID,
                                parent,
                                repMembrIndx,
                                repMembrSize,
                                nonMissMembrIndx,
                                nonMissMembrSize,
                                eventTimeIndex,
                                eventTimeSize,
                                revEventTimeIndex,
                                revEventTimeSize,
                                revParentSurvival,
                                qeTimeIndex,
                                qeTimeSize,
                                & gHat,
                                & w_ktm,
                                & gamma_ktm);
          leftGammaSum  = dvector(1, qeTimeSize + 1);
          rightGammaSum = dvector(1, qeTimeSize + 1);
          leftGammaBar  = dvector(1, qeTimeSize + 1);
          rightGammaBar = dvector(1, qeTimeSize + 1);
        }
        if ((eventTimeSize > 0) && (qeTimeSize > 0)) {
          leftSize = 0;
          priorMembrIter = 0;
          if (factorFlag == FALSE) {
            for (j = 1; j <= nonMissMembrSize; j++) {
              localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
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
            ${trace.token}          if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
            ${trace.token}            if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
            ${trace.token}              RF_nativePrint("\n PriorIter:     %10d  CurrentIter:   %10d", priorMembrIter, currentMembrIter);
            ${trace.token}            }
            ${trace.token}          }
            if ((leftSize != 0) && (rghtSize != 0)) {
              sumLeft = sumRght = 0.0;
              tIndx = 1;
              adHocFlag = FALSE;
              while ((tIndx <= qeTimeSize) && !adHocFlag) {
                if (qeTimeIndex[tIndx] > 0) {
                  leftGammaSum[tIndx] = rightGammaSum[tIndx] = 0.0;
                  for (k = 1; k <= nonMissMembrSize; k++) {
                    if (RF_nativeIsNaN(gamma_ktm[qeTimeIndex[tIndx]][k])) {
                      adHocFlag = TRUE;
                      k = nonMissMembrSize;
                      tIndx = qeTimeSize;
                    }
                    else {
                      if (localSplitIndicator[  nonMissMembrIndx[k]  ] == LEFT) {
                        leftGammaSum[tIndx]  += gamma_ktm[qeTimeIndex[tIndx]][k];
                      }
                      else {
                        rightGammaSum[tIndx] += gamma_ktm[qeTimeIndex[tIndx]][k];
                      }
                    }
                  }
                  if (!adHocFlag) {
                    leftGammaBar[tIndx] = leftGammaSum[tIndx] / leftSize;
                    rightGammaBar[tIndx] = rightGammaSum[tIndx] / rghtSize;
                    sumLeft  += pow(leftGammaBar[tIndx],  2);
                    sumRght += pow(rightGammaBar[tIndx], 2);
                  }
                }
                tIndx++;
              }  
              if (!adHocFlag) {
                delta = ( (((double) leftSize / nonMissMembrSize) * sumLeft) + (((double) rghtSize / nonMissMembrSize) * sumRght));
              }
              else {
                delta = RF_nativeNaN;
              }
              ${trace.token}    if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
              ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
              ${trace.token}      RF_nativePrint("\nRKM-BS Final:  (sumLeft, sumRght, statistic) = (%10.4f, %10.4f, %10.4f) ", sumLeft, sumRght, delta);
              ${trace.token}    }
              ${trace.token}    }
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
        else {
          ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
          ${trace.token}  if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
          ${trace.token}    RF_nativePrint("\nCovariate ignored due to zero eventTimeSize:  %10d", covariate);
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
          free_dvector(leftGammaSum,  1, qeTimeSize + 1);
          free_dvector(rightGammaSum, 1, qeTimeSize + 1);
          free_dvector(leftGammaBar,  1, qeTimeSize + 1);
          free_dvector(rightGammaBar, 1, qeTimeSize + 1);
          unstackLocalGamma(treeID,
                            nonMissMembrSize,
                            eventTimeIndex,
                            eventTimeSize,
                            qeTimeIndex,
                            qeTimeSize,
                            gamma_ktm);
          unstackQETime(treeID,
                        eventTimeSize,
                        qeTimeIndex);
          unstackSplitSurv(treeID,
                           parent,
                           eventTimeCount,
                           eventTimeIndex,
                           eventTimeSize,
                           parentEvent,
                           parentAtRisk,
                           leftEvent,
                           leftAtRisk,
                           rightEvent,
                           rightAtRisk);
          unstackSplitSurv(treeID,
                           parent,
                           revEventTimeCount,
                           revEventTimeIndex,
                           revEventTimeSize,
                           revParentEvent,
                           revParentAtRisk,
                           revLeftEvent,
                           revLeftAtRisk,
                           revRightEvent,
                           revRightAtRisk);
          unstackAndGetSplitSurv2(treeID,
                                  parent,
                                  eventTimeSize,
                                  parentSurvival);
          unstackAndGetSplitSurv2(treeID,
                                  parent,
                                  revEventTimeSize,
                                  revParentSurvival);
        }
      }      
    }  
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      free_dvector(leftGammaSum,  1, qeTimeSize + 1);
      free_dvector(rightGammaSum, 1, qeTimeSize + 1);
      free_dvector(leftGammaBar,  1, qeTimeSize + 1);
      free_dvector(rightGammaBar, 1, qeTimeSize + 1);
      unstackLocalGamma(treeID,
                        nonMissMembrSizeStatic,
                        eventTimeIndex,
                        eventTimeSize,
                        qeTimeIndex,
                        qeTimeSize,
                        gamma_ktm);
      unstackQETime(treeID,
                    eventTimeSize,
                    qeTimeIndex);
      unstackSplitSurv(treeID,
                       parent,
                       eventTimeCount,
                       eventTimeIndex,
                       eventTimeSize,
                       parentEvent,
                       parentAtRisk,
                       leftEvent,
                       leftAtRisk,
                       rightEvent,
                       rightAtRisk);
      unstackSplitSurv(treeID,
                       parent,
                       revEventTimeCount,
                       revEventTimeIndex,
                       revEventTimeSize,
                       revParentEvent,
                       revParentAtRisk,
                       revLeftEvent,
                       revLeftAtRisk,
                       revRightEvent,
                       revRightAtRisk);
      unstackAndGetSplitSurv2(treeID,
                              parent,
                              eventTimeSize,
                              parentSurvival);
      unstackAndGetSplitSurv2(treeID,
                              parent,
                              revEventTimeSize,
                              revParentSurvival);
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
  ${trace.token}    RF_nativePrint("\nbrierScoreGradient1(%10d) result:  %10d", treeID, result);
  ${trace.token}    RF_nativePrint("\nbrierScoreGradient1(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
  return result;
}
