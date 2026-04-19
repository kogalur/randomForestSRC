
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "splitUtilSurv.h"
#include "factorOps.h"
#include "nrutil.h"
#include "error.h"
void stackAndGetSplitSurv(uint    treeID,
                          Node   *parent,
                          char    eventType,
                          uint  **eventTimeCount,
                          uint  **eventTimeIndex,
                          uint   *eventTimeSize,
                          uint  **parentEvent,
                          uint  **parentAtRisk,
                          uint  **leftEvent,
                          uint  **leftAtRisk,
                          uint  **rightEvent,
                          uint  **rightAtRisk) {
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAndGetSplitSurv() ENTRY ...\n");
  ${trace.token}  }
  uint *repMembrIndx = parent -> repMembrIndx;
  uint  repMembrSize = parent -> repMembrSize;
  uint *nonMissMembrIndx = parent -> nonMissMembrIndx;
  uint  nonMissMembrSize = parent -> nonMissMembrSize;
  *eventTimeCount = uivector(1, RF_masterTimeSize);
  *eventTimeIndex = uivector(1, RF_masterTimeSize);
  *eventTimeSize = getEventTime(treeID,
                                parent,
                                repMembrIndx,
                                repMembrSize,
                                nonMissMembrIndx,
                                nonMissMembrSize,
                                eventType,
                                *eventTimeCount,
                                *eventTimeIndex);
  stackSplitEventAndRisk(treeID,
                         parent,
                         *eventTimeSize,
                          parentEvent,
                          parentAtRisk,
                          leftEvent,
                          leftAtRisk,
                          rightEvent,
                          rightAtRisk);
  getSplitEventAndRisk( treeID,
                        parent,
                        repMembrIndx,
                        repMembrSize,
                        nonMissMembrIndx,
                        nonMissMembrSize,
                        *eventTimeCount,
                        *eventTimeIndex,
                        *eventTimeSize,
                        *parentEvent,
                        *parentAtRisk);
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAndGetSplitSurv() EXIT ...\n");
  ${trace.token}  }
}
void unstackSplitSurv(uint    treeID,
                      Node   *parent,
                      uint *eventTimeCount,
                      uint *eventTimeIndex,
                      uint  eventTimeSize,
                      uint *parentEvent,
                      uint *parentAtRisk,
                      uint *leftEvent,
                      uint *leftAtRisk,
                      uint *rightEvent,
                      uint *rightAtRisk) {
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nunstackSplitSurv() ENTRY ...\n");
  ${trace.token}    }
  ${trace.token}  }
  free_uivector(eventTimeCount, 1, RF_masterTimeSize);
  free_uivector(eventTimeIndex, 1, RF_masterTimeSize);
  unstackSplitEventAndRisk(treeID,
                           parent,
                           eventTimeSize,
                           parentEvent,
                           parentAtRisk,
                           leftEvent,
                           leftAtRisk,
                           rightEvent,
                           rightAtRisk);
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nunstackSplitSurv() EXIT ...\n");
  ${trace.token}    }
  ${trace.token}  }
}
void stackSplitSurv3(uint    treeID,
                     Node   *parent,
                     uint   eventTimeSize,
                     double **leftLocalRatio,
                     double **rightLocalRatio,
                     double **leftLocalSurvival,
                     double **rightLocalSurvival,
                     uint   revEventTimeSize,
                     double **leftRevLocalRatio,
                     double **rightRevLocalRatio,
                     double **leftRevLocalSurvival,
                     double **rightRevLocalSurvival,
                     double **leftBS,
                     double **rightBS) {
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nstackSplitSurv3() ENTRY ...\n");
  ${trace.token}    }
  ${trace.token}  }
  if (eventTimeSize > 0) {
    *leftLocalRatio     = dvector(1, eventTimeSize);
    *rightLocalRatio    = dvector(1, eventTimeSize);
    *leftLocalSurvival  = dvector(1, eventTimeSize);
    *rightLocalSurvival  = dvector(1, eventTimeSize);
    *leftBS              = dvector(1, eventTimeSize);
    *rightBS             = dvector(1, eventTimeSize);
  }
  else {
    *leftLocalRatio = *rightLocalRatio = *leftLocalSurvival = *rightLocalSurvival = *leftBS = *rightBS = NULL;
  }
  if (revEventTimeSize > 0) {
    *leftRevLocalRatio     = dvector(1, revEventTimeSize);
    *rightRevLocalRatio    = dvector(1, revEventTimeSize);
    *leftRevLocalSurvival  = dvector(1, revEventTimeSize);
    *rightRevLocalSurvival = dvector(1, revEventTimeSize);
  }
  else {
    *leftRevLocalRatio = *rightRevLocalRatio = *leftRevLocalSurvival = *rightRevLocalSurvival = NULL;
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nstackSplitSurv3() EXIT ...\n");
  ${trace.token}    }
  ${trace.token}  }
}
void unstackSplitSurv3(uint    treeID,
                       Node   *parent,
                       uint   eventTimeSize,
                       double *leftLocalRatio,
                       double *rightLocalRatio,
                       double *leftLocalSurvival,
                       double *rightLocalSurvival,
                       uint   revEventTimeSize,
                       double *leftRevLocalRatio,
                       double *rightRevLocalRatio,
                       double *leftRevLocalSurvival,
                       double *rightRevLocalSurvival,
                       double *leftBS,
                       double *rightBS) {
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nunstackSplitSurv3() ENTRY ...\n");
  ${trace.token}    }
  ${trace.token}  }
  if (eventTimeSize > 0) {
    free_dvector(leftLocalRatio, 1, eventTimeSize);
    free_dvector(rightLocalRatio, 1, eventTimeSize);
    free_dvector(leftLocalSurvival, 1, eventTimeSize);
    free_dvector(rightLocalSurvival, 1, eventTimeSize);
    free_dvector(leftBS, 1, eventTimeSize);
    free_dvector(rightBS, 1, eventTimeSize);
  }
  if (revEventTimeSize > 0) {
    free_dvector(leftRevLocalRatio, 1, revEventTimeSize);
    free_dvector(rightRevLocalRatio, 1, revEventTimeSize);
    free_dvector(leftRevLocalSurvival, 1, revEventTimeSize);
    free_dvector(rightRevLocalSurvival, 1, revEventTimeSize);
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nunstackSplitSurv3() EXIT ...\n");
  ${trace.token}    }
  ${trace.token}  }
}
uint getEventTime(uint   treeID,
                  Node  *parent,
                  uint   *repMembrIndx,
                  uint    repMembrSize,
                  uint   *nonMissMembrIndx,
                  uint    nonMissMembrSize,
                  char    eventType,
                  uint   *eventTimeCount,
                  uint   *eventTimeIndex) {
  uint i;
  uint eventTimeSize;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetEventTime() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}  uint parentEventCount = 0;
  eventTimeSize = 0;
  for (i=1; i <= RF_masterTimeSize; i++) {
    eventTimeCount[i] = 0;
  }
  if (eventType) {
    for (i = 1; i <= nonMissMembrSize; i++) {
      if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[i]] ] > 0) {
        eventTimeCount[RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[i]] ]] ++;
        ${trace.token} parentEventCount ++;
      }
    }
  }
  else {
    for (i = 1; i <= nonMissMembrSize; i++) {
      if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[i]] ] == 0) {
        eventTimeCount[RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[i]] ]] ++;
        ${trace.token} parentEventCount ++;
      }
    }
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nNode specific analysis for:  (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
  ${trace.token}    if (eventType) {
  ${trace.token}      RF_nativePrint("\nParent Raw Event Count:  %10d \n", parentEventCount);
  ${trace.token}    }
  ${trace.token}    else {
  ${trace.token}      RF_nativePrint("\nParent Raw Reverse Event Count:  %10d \n", parentEventCount);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\nLocal Membership Index for Parent Node: \n");
  ${trace.token}    RF_nativePrint("\n       idx   replcIdx mastTimIdx     status         time \n");
  ${trace.token}    for (i = 1; i <= nonMissMembrSize; i++) {
  ${trace.token}      RF_nativePrint("%10d %10d %10d %10d %12.4f \n", i, repMembrIndx[nonMissMembrIndx[i]], RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[i]] ], (uint) RF_status[treeID][ repMembrIndx[nonMissMembrIndx[i]] ], RF_time[treeID][ repMembrIndx[nonMissMembrIndx[i]] ]);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\nLocal Time Counts:  \n");
  ${trace.token}    for (i=1; i <= RF_masterTimeSize; i++) {
  ${trace.token}      RF_nativePrint("%10d %10d %12.4f \n", i, eventTimeCount[i], RF_masterTime[i]);
  ${trace.token}    }
  ${trace.token}  }
  for (i=1; i <= RF_masterTimeSize; i++) {
    if (eventTimeCount[i] > 0) {
      eventTimeIndex[++eventTimeSize] = i;
    }
  }
  ${trace.token}    if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}      RF_nativePrint("\nNode specific Event Times for:  (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);  
  ${trace.token}      RF_nativePrint("\n  Event Times (i, RF_masterTimeIndex): \n");
  ${trace.token}      for (i=1; i <= eventTimeSize; i++) {
  ${trace.token}        RF_nativePrint("%10d %10d \n", i, eventTimeIndex[i]);
  ${trace.token}      }
  ${trace.token}    }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetEventTime() EXIT ...\n");
  ${trace.token}  }
  return (eventTimeSize);
}
void stackSplitEventAndRisk(uint    treeID,
                            Node   *parent,
                            uint    genEventTimeSize,
                            uint  **genParentEvent,
                            uint  **genParentAtRisk,
                            uint  **genLeftEvent,
                            uint  **genLeftAtRisk,
                            uint  **genRightEvent,
                            uint  **genRightAtRisk) {
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nstackSplitEventAndRisk() ENTRY ...\n");
  ${trace.token}    }
  ${trace.token}  }
  if (genEventTimeSize > 0) {
    *genParentEvent  = uivector(1, genEventTimeSize);
    *genParentAtRisk = uivector(1, genEventTimeSize);
    *genLeftEvent    = uivector(1, genEventTimeSize);
    *genLeftAtRisk   = uivector(1, genEventTimeSize);
    *genRightEvent   = uivector(1, genEventTimeSize);
    *genRightAtRisk  = uivector(1, genEventTimeSize);
  }
  else {
    *genParentEvent  = *genParentAtRisk = *genLeftEvent  = *genLeftAtRisk = *genRightEvent  = *genRightAtRisk = NULL;
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nstackSplitEventAndRisk() EXIT ...\n");
  ${trace.token}    }
  ${trace.token}  }
}
void unstackSplitEventAndRisk(uint    treeID,
                              Node   *parent,
                              uint    genEventTimeSize,
                              uint   *genParentEvent,
                              uint   *genParentAtRisk,
                              uint   *genLeftEvent,
                              uint   *genLeftAtRisk,
                              uint   *genRightEvent,
                              uint   *genRightAtRisk) {
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nunstackSplitEventAndRisk() ENTRY ...\n");
  ${trace.token}    }
  ${trace.token}  }
  if (genEventTimeSize > 0) {
    free_uivector(genParentEvent, 1, genEventTimeSize);
    free_uivector(genParentAtRisk, 1, genEventTimeSize);
    free_uivector(genLeftEvent, 1, genEventTimeSize);
    free_uivector(genLeftAtRisk, 1, genEventTimeSize);
    free_uivector(genRightEvent, 1, genEventTimeSize);
    free_uivector(genRightAtRisk, 1, genEventTimeSize);
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nunstackSplitEventAndRisk() EXIT ...\n");
  ${trace.token}    }
  ${trace.token}  }
}
void getSplitEventAndRisk(uint    treeID,
                          Node   *parent,
                          uint   *repMembrIndx,
                          uint    repMembrSize,
                          uint   *nonMissMembrIndx,
                          uint    nonMissMembrSize,
                          uint   *eventTimeCount,
                          uint   *eventTimeIndex,
                          uint    eventTimeSize,
                          uint   *parentEvent,
                          uint   *parentAtRisk) {
  uint i, j;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nNode specific Event Counts for:  (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
  ${trace.token}    RF_nativePrint("\nLocal Event Counts (timIdx, events): \n");
  ${trace.token}  }
  for (i = 1; i <= eventTimeSize; i++) {
    parentAtRisk[i] = 0;
    parentEvent[i]  = eventTimeCount[eventTimeIndex[i]];
    ${trace.token}    if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
    ${trace.token}      RF_nativePrint("%10d %10d \n", i, parentEvent[i]);
    ${trace.token}    }
    for (j = 1; j <= nonMissMembrSize; j++) {
      if (eventTimeIndex[i] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[j]] ]) {
        parentAtRisk[i] ++;
      }
    }
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nNode specific At Risk Counts for:  (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
  ${trace.token}    RF_nativePrint("\nLocal At Risk Counts (timIdx, at risk): \n");
  ${trace.token}    for (i=1; i <= eventTimeSize; i++) {
  ${trace.token}      RF_nativePrint("%10d %10d \n", i, parentAtRisk[i]);
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetSplitEventAndRisk() EXIT ...\n");
  ${trace.token}  }
}
void stackAndGetSplitSurv2(uint     treeID,
                           Node    *parent,
                           uint     eventTimeSize,
                           uint    *parentEvent,
                           uint    *parentAtRisk,
                           double **localSurvival) {
  double *localRatio;
  uint q;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAndGetSplitSurv2() ENTRY ...\n");
  ${trace.token}  }
  localRatio = dvector(1, eventTimeSize + 1);
  *localSurvival = dvector(1, eventTimeSize + 1);
  for (q = 1; q <= eventTimeSize; q++) {
    if (parentEvent[q] > 0) {
      if (parentAtRisk[q] >= 1) {
        localRatio[q] = ((double) parentEvent[q] / parentAtRisk[q]);
      }
      else {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Zero At Risk Count encountered in local ratio calculation for (tree, leaf) = (%10d, %10d)", treeID, parent -> nodeID);
        RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
        RF_nativeExit();
      }
    }
    else {
      localRatio[q] = 0.0;
    }
    (*localSurvival)[q] = 1.0 - localRatio[q];
  }
  for (q = 2; q <= eventTimeSize; q++) {
    (*localSurvival)[q] *= (*localSurvival)[q-1];
  }
  free_dvector(localRatio, 1, eventTimeSize + 1);
  ${trace.token}    if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}      if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}        RF_nativePrint("\nNode specific local survival function for:  (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
  ${trace.token}        RF_nativePrint("    mTimIdx   survival \n");
  ${trace.token}        for (q = 1; q <= eventTimeSize; q++) {
  ${trace.token}          RF_nativePrint(" %10d %10.4f \n", q, (*localSurvival)[q]);
  ${trace.token}        }
  ${trace.token}      }
  ${trace.token}    }      
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAndGetSplitSurv2() EXIT ...\n");
  ${trace.token}  }
}
void unstackAndGetSplitSurv2(uint     treeID,
                             Node    *parent,
                             uint     eventTimeSize,
                             double  *localSurvival) {
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nunstackAndGetSplitSurv2() ENTRY ...\n");
  ${trace.token}    }
  ${trace.token}  }
  free_dvector(localSurvival, 1, eventTimeSize + 1);
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nunstackAndGetSplitSurv2() EXIT ...\n");
  ${trace.token}    }
  ${trace.token}  }
}
void stackAndGetFZhat(uint  treeID,
                      Node *parent,
                      uint *repMembrIndx,
                      uint  repMembrSize,
                      uint *nonMissMembrIndx,
                      uint  nonMissMembrSize,
                      uint *eventTimeIndex,
                      uint  eventTimeSize,
                      uint *revEventTimeIndex,
                      uint  revEventTimeSize,
                      double *revParentSurvival,
                      double **fZHat) {
  double gHatPrevious, gHatCurrent;
  double w_it, w_jt, y_it;
  double denom;
  char escapeFlag;
  uint tIndx;
  uint i, j, t;
  char adHocFlag, adHocFlagSum;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAndGetFZhat() ENTRY ...\n");
  ${trace.token}  }
  *fZHat        = dvector(1, eventTimeSize);
  escapeFlag = FALSE;
  gHatPrevious = gHatCurrent = 1.0;
  tIndx = 1;
  while (!escapeFlag) {
    if(tIndx <= revEventTimeSize) {
      if(revEventTimeIndex[tIndx] < eventTimeIndex[1]) {
        gHatPrevious = gHatCurrent = revParentSurvival[tIndx];
        tIndx++;
      }
      else {
        escapeFlag = TRUE;
      }
    }
    else {
      escapeFlag = TRUE;
    }
  }
  tIndx = 1;
  for (t = 1; t <= eventTimeSize; t++) {
    escapeFlag = FALSE;
    gHatPrevious = gHatCurrent;
    if (revEventTimeSize > 0) {
      while (!escapeFlag) {
        if(tIndx <= revEventTimeSize) {
          if(revEventTimeIndex[tIndx] < eventTimeIndex[t]) {
            gHatCurrent = revParentSurvival[tIndx];
            tIndx++;
          }
          else {
            escapeFlag = TRUE;
          }
        }
        else {
          escapeFlag = TRUE;
        }
      }
    }
    denom = 0.0;
    adHocFlagSum = FALSE;
    for (j = 1; j <= repMembrSize; j++) {
      if (RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[j]] ] > eventTimeIndex[t]) {
        if (gHatCurrent > 0.0) {        
          w_jt = 1.0 / gHatCurrent;
        }
        else {
          adHocFlagSum = TRUE;
        }
      }
      else {
        if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[j]] ] > 0) {
          if (gHatPrevious > 0.0) {
            w_jt = 1.0 / gHatPrevious;
          }
          else {
            adHocFlagSum = TRUE;
          }
        }
        else {
          w_jt = 0.0;
        }
      }
      if (!adHocFlagSum) {
        denom += w_jt;
      }
      else {
        denom = RF_nativeNaN;
      }
    }
    (*fZHat)[t] = 0.0;
    for (i = 1; i <= repMembrSize; i++) {
      adHocFlag = FALSE;
      if (RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[i]] ] > eventTimeIndex[t] ) {      
        y_it = 1.0;
        if (gHatCurrent > 0.0) {
          w_it = 1.0 / gHatCurrent;
        }
        else {
          adHocFlag = TRUE;
        }
        if (!adHocFlagSum && !adHocFlag) {
          (*fZHat)[t] += (w_it / denom) * y_it;
        }
        else {
          (*fZHat)[t] = RF_nativeNaN;
        }
      }
      else {
        y_it = 0.0;
        (*fZHat)[t] += 0.0;
      }
    }
  }
  ${trace.token}    if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nFZhat for virtual split in (treeID, nodeID) = (%10d, %10d) \n", treeID, parent -> nodeID);
  ${trace.token}      RF_nativePrint("          k    timeIdx       time      fZHat \n");
  ${trace.token}      for (uint k = 1; k <= eventTimeSize; k++) {
  ${trace.token}        RF_nativePrint(" %10d %10d %10.4f %10.4f \n",
  ${trace.token}            k, eventTimeIndex[k], RF_masterTime[eventTimeIndex[k]], (*fZHat)[k]);
  ${trace.token}      }
  ${trace.token}    }
  ${trace.token}    } 
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAndGetFZhat() EXIT ...\n");
  ${trace.token}  }
}
void unstackFZhat(uint  treeID,
                  Node *parent,
                  uint  eventTimeSize,
                  double *fZHat) {
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nunstackFZhat() ENTRY ...\n");
  ${trace.token}  }
  free_dvector(fZHat, 1, eventTimeSize);
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nunstackFZhat() EXIT ...\n");
  ${trace.token}  }
}
double getW_kt(uint  treeID,
               Node *parent,
               uint  indv,
               uint  tIndx,
               uint *eventTimeIndex,
               uint *revEventTimeIndex,
               uint  revEventTimeSize,
               double *revParentSurvival,
               double *gHatPreviousX,
               double *gHatCurrentX) {
  double w_kt;
  double gHatPrevious, gHatCurrent;
  char escapeFlag;
  uint t;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetW_kt() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}    if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngHat for virtual split in (treeID, nodeID) = (%10d, %10d) \n", treeID, parent -> nodeID);
  ${trace.token}      RF_nativePrint("       indv      tIndx   evtmIndx  gHatPrevInc  gHatCurrInc  gHatPrevOut  gHatCurrOut         W_kt \n");
  ${trace.token}        RF_nativePrint(" %10d %10d %10d %12.4f %12.4f",
  ${trace.token}            indv, tIndx, eventTimeIndex[tIndx], *gHatPreviousX, *gHatCurrentX);
  ${trace.token}    }
  ${trace.token}    } 
  escapeFlag = FALSE;
  gHatPrevious = gHatCurrent = 1.0;
  t = 1;
  w_kt = 0.0;  
  while (!escapeFlag) {
    if(t <= revEventTimeSize) {
      if(revEventTimeIndex[t] < eventTimeIndex[tIndx]) {
        gHatCurrent = revParentSurvival[t];
        t++;
      }
      else {
        escapeFlag = TRUE;
      }
    }
    else {
      escapeFlag = TRUE;
    }
  }
  if (t > 1) {
    gHatPrevious = revParentSurvival[t-1];
  }
  if (RF_masterTimeIndex[treeID][ indv ] > eventTimeIndex[tIndx] ) {      
    if (gHatCurrent > 0.0) {
      w_kt = 1.0 / gHatCurrent;
    }
    else {
      w_kt = RF_nativeNaN;
    }
  }
  else {
    if (RF_status[treeID][ indv ] > 0) {
      if (gHatPrevious > 0.0) {      
        w_kt = 1.0 / gHatPrevious;
      }
      else {
        w_kt = RF_nativeNaN;
      }
    }
    else {
      w_kt = 0.0;
    }
  }
  *gHatPreviousX = gHatPrevious;
  *gHatCurrentX  = gHatCurrent;
  ${trace.token}    if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}        RF_nativePrint(" %12.4f %12.4f %12.4f \n", *gHatPreviousX, *gHatCurrentX, w_kt);
  ${trace.token}    }
  ${trace.token}    } 
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetW_kt() EXIT ...\n");
  ${trace.token}  }
  return w_kt;
}
void stackAndGetLocalGamma(uint  treeID,
                           Node *parent,
                           uint *repMembrIndx,
                           uint  repMembrSize,
                           uint *nonMissMembrIndx,
                           uint  nonMissMembrSize,
                           uint *eventTimeIndex,
                           uint  eventTimeSize,
                           uint *revEventTimeIndex,
                           uint  revEventTimeSize,
                           double *revParentSurvival,
                           uint      *qeTimeIndex,
                           uint       qeTimeSize,
                           double   **gHat,
                           double  ***w_ktm,
                           double    ***gamma_ktm) {
  uint tIndx;
  double gHatCurrent;
  double *fHatDenom;
  double *fHat;
  double *y_kt;
  uint tauTimeIdx;
  char escapeFlag;
  char adHocFlagSum;
  uint k, t;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAndGetLocalGamma() ENTRY ...\n");
  ${trace.token}  }
  if ((eventTimeSize > 0) && (qeTimeSize > 0)) {
    if (qeTimeIndex[qeTimeSize] > 0) {
      tauTimeIdx = eventTimeIndex[qeTimeIndex[qeTimeSize]];
      ${trace.token}    if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
      ${trace.token}      if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
      ${trace.token}        RF_nativePrint("\n\nLocal WIBS Tau Time Specified: (tauTimeIdx, time absolute) = (%10d, %10.4f) ", tauTimeIdx, RF_masterTime[tauTimeIdx]);
      ${trace.token}      }
      ${trace.token}    }
    }
    else {
      tauTimeIdx = 0;
      ${trace.token}    if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
      ${trace.token}      if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
      ${trace.token}        RF_nativePrint("\n\nLocal WIBS Tau Time Specified: (tauTimeIdx, time absolute) = (%10d, Not Applicable) ", tauTimeIdx);
      ${trace.token}      }
      ${trace.token}    }
    }
      *gamma_ktm = (double **) new_vvector(1, eventTimeSize, NRUTIL_DPTR);
      *gHat = dvector(0, eventTimeSize);
      y_kt = dvector(1, nonMissMembrSize);
      *w_ktm = (double **) new_vvector(1, eventTimeSize, NRUTIL_DPTR);
      fHat = dvector(1, eventTimeSize);
      fHatDenom = dvector(1, eventTimeSize);
      escapeFlag = FALSE;
      (*gHat)[0] = 1.0;
      tIndx = 1;
      while (!escapeFlag) {
        if(tIndx <= revEventTimeSize) {
          if(revEventTimeIndex[tIndx] < eventTimeIndex[1]) {
            (*gHat)[0] = revParentSurvival[tIndx];
            tIndx++;
          }
          else {
            escapeFlag = TRUE;
          }
        }
        else {
          escapeFlag = TRUE;
        }
      }
      gHatCurrent = (*gHat)[0];
      tIndx = 1;
      for (t = 1; t <= eventTimeSize; t++) {
        (*w_ktm)[t] = dvector(1, nonMissMembrSize);
        fHatDenom[t] = 0.0;
        adHocFlagSum = FALSE;
        ${trace.token}    if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
        ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
        ${trace.token}      RF_nativePrint("\nw_kt[][] for virtual split in (treeID, nodeID) = (%10d, %10d) \n", treeID, parent -> nodeID);
        ${trace.token}      RF_nativePrint("      tIndx   evtmIndx          k       indv         y_kt         w_kt \n");  
        ${trace.token}    }
        ${trace.token}    } 
        escapeFlag = FALSE;
        if (revEventTimeSize > 0) {
          while (!escapeFlag) {
            if(tIndx <= revEventTimeSize) {
              if(revEventTimeIndex[tIndx] < eventTimeIndex[t]) {
                gHatCurrent= revParentSurvival[tIndx];
                tIndx++;
              }
              else {
                escapeFlag = TRUE;
              }
            }
            else {
              escapeFlag = TRUE;
            }
          }
        }
        (*gHat)[t] = gHatCurrent;
        ${trace.token}    if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
        ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
        ${trace.token}      RF_nativePrint("\ngHatCurrent for virtual split in (treeID, nodeID) = (%10d, %10d) \n", treeID, parent -> nodeID);
        ${trace.token}      RF_nativePrint("      tIndx   evtmIndx    eventTime         gHat \n");  
        ${trace.token}      for (uint tt = 1; tt <= eventTimeSize; tt++) {
        ${trace.token}        RF_nativePrint(" %10d %10d %12.4f %12.4f \n", tt, eventTimeIndex[tt], RF_masterTime[eventTimeIndex[tt]], (*gHat)[tt]);
        ${trace.token}      }
        ${trace.token}    }
        ${trace.token}    } 
        for (k = 1; k <= nonMissMembrSize; k++) {
          if (RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[k]] ] > eventTimeIndex[t]) {
            y_kt[k] = 1.0;
            if ((*gHat)[t] > 0.0) {        
              (*w_ktm)[t][k] = 1.0 / (*gHat)[t];
            }
            else {
              (*w_ktm)[t][k] = RF_nativeNaN;
              adHocFlagSum = TRUE;
            }
          }
          else {
            y_kt[k] = 0.0;
            if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[k]] ] > 0) {
              if ((*gHat)[t-1] > 0.0) {
                (*w_ktm)[t][k] = 1.0 / (*gHat)[t-1];
              }
              else {
                (*w_ktm)[t][k] = RF_nativeNaN;
                adHocFlagSum = TRUE;
              }
            }
            else {
              (*w_ktm)[t][k] = 0.0;
            }
          }
          if (!adHocFlagSum) {
            fHatDenom[t] += (*w_ktm)[t][k];
          }
          ${trace.token}    if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
          ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
          ${trace.token}      RF_nativePrint(" %10d %10d %10d %10d %12.4f %12.4f \n", t, eventTimeIndex[t], k, repMembrIndx[nonMissMembrIndx[k]], y_kt[k], (*w_ktm)[t][k]);
          ${trace.token}    }
          ${trace.token}    } 
        }
        if (adHocFlagSum) {
          fHatDenom[t] = RF_nativeNaN;
        }
        fHat[t] = 0.0;
        for (k = 1; k <= nonMissMembrSize; k++) {
          if (y_kt[k] != 0) { 
            if ( !RF_nativeIsNaN((*w_ktm)[t][k]) && !RF_nativeIsNaN(fHatDenom[t])) {
              fHat[t] += (*w_ktm)[t][k] / fHatDenom[t];
            }
            else {
              fHat[t] = RF_nativeNaN;
              k = nonMissMembrSize;
            }
          }
          else {
          }
        }
        ${trace.token}    if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
        ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
        ${trace.token}      RF_nativePrint("\nfHat and fHatDenom for virtual split in (treeID, nodeID) = (%10d, %10d) \n", treeID, parent -> nodeID);
        ${trace.token}      RF_nativePrint("      tIndx   evtmIndx    fHatDenom         fHat \n");  
        ${trace.token}      RF_nativePrint(" %10d %10d %12.4f %12.4f \n", t, eventTimeIndex[t], fHatDenom[t], fHat[t]);
        ${trace.token}    }
        ${trace.token}    } 
        if (eventTimeIndex[t] <= tauTimeIdx) {
          (*gamma_ktm)[t] = dvector(1, nonMissMembrSize);
          ${trace.token}    if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
          ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
          ${trace.token}      RF_nativePrint("\ngamma_kt[][] for virtual split in (treeID, nodeID) = (%10d, %10d) \n", treeID, parent -> nodeID);
          ${trace.token}      RF_nativePrint("   evtmIndx  eventTime          k     gamma_kt \n");  
          ${trace.token}    }
          ${trace.token}    } 
          for (k = 1; k <= nonMissMembrSize; k++) {
            if (RF_nativeIsNaN((*w_ktm)[t][k]) || RF_nativeIsNaN(fHat[t])) {
              (*gamma_ktm)[t][k] = RF_nativeNaN;
              k = nonMissMembrSize;
            }
            else {
              (*gamma_ktm)[t][k] = - 2.0 * (*w_ktm)[t][k] * (y_kt[k] - fHat[t]);
            }
          }
        }
        else {
          (*gamma_ktm)[t] = NULL;
          ${trace.token}    if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
          ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
          ${trace.token}      RF_nativePrint("\ngamma_kt[][] for virtual split in (treeID, nodeID) = (%10d, %10d) \n", treeID, parent -> nodeID);
          ${trace.token}      RF_nativePrint(" %10d %10d not in bounds of tauTimeIdx \n", t, eventTimeIndex[t]);
          ${trace.token}    }
          ${trace.token}    } 
        }
      }
      free_dvector(*gHat, 0, eventTimeSize);
      free_dvector(y_kt, 1, nonMissMembrSize);
      for (t = 1; t <= eventTimeSize; t++) {
        if ((*w_ktm)[t] != NULL) {
          free_dvector((*w_ktm)[t], 1, nonMissMembrSize);
        }
      }
      free_new_vvector(*w_ktm, 1, eventTimeSize, NRUTIL_DPTR);  
      free_dvector(fHat, 1, eventTimeSize);
      free_dvector(fHatDenom, 1, eventTimeSize);
  }
  else {
    ${trace.token}    if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
    ${trace.token}      if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
    ${trace.token}        RF_nativePrint("\nZero event time size or QE time size in stackAndGetLocalGamma(), gamma calculations omitted.");
    ${trace.token}      }
    ${trace.token}    }
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAndGetLocalGamma() EXIT ...\n");
  ${trace.token}  }
}
void  unstackLocalGamma(uint    treeID,
                        uint    nonMissMembrSize,
                        uint   *eventTimeIndex,
                        uint    eventTimeSize,
                           uint      *qeTimeIndex,
                        uint       qeTimeSize,
                        double **gamma_ktm) {
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nunstackLocalGamma() ENTRY ...");
  ${trace.token}  }
  if ((eventTimeSize > 0) && (qeTimeSize > 0)) {
    for (uint t = 1; t <= eventTimeSize; t++) {
      if (gamma_ktm[t] != NULL) {
        free_dvector(gamma_ktm[t], 1, nonMissMembrSize);
        gamma_ktm[t] = NULL;
      }
    }  
    free_new_vvector(gamma_ktm, 1, eventTimeSize, NRUTIL_DPTR);
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nunstackLocalGamma() EXIT ...");
  ${trace.token}  }
}
void stackAndGetQTime(uint  treeID,
                      Node *parent,
                      uint  eventTimeSize,
                      double *survival,
                      uint  **quantileTime) {
  uint itr;
  uint k;
  char found;
  *quantileTime = uivector(1, RF_quantileSize);
  itr = 1;
  for (k = 1; k <= RF_quantileSize; k++) {
    found = FALSE;
    while (found == FALSE) {
      (*quantileTime)[k] = itr;
      if (itr > eventTimeSize) {
        found = TRUE;
      }
      else {
        if (survival[itr] > (1.0 - RF_quantile[k])) {
          itr ++;
        }
        else {
          found = TRUE;
        }
      }
    }
    (*quantileTime)[k] --;
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}  if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}    RF_nativePrint("\nQS Split Rule quantile times for split statistic (treeID, nodeID) = (%10d, %10d) \n", treeID, parent -> nodeID);
  ${trace.token}    RF_nativePrint("          k      quant  1 - quant       surv      qTime \n");
  ${trace.token}    for (k = 1; k <= RF_quantileSize; k++) {
  ${trace.token}      RF_nativePrint(" %10d %10.4f %10.4f %10.4f %10d \n",
  ${trace.token}         k, RF_quantile[k], 1 - RF_quantile[k], ((*quantileTime)[k] == 0) ? 1.0 : survival[(*quantileTime)[k]], (*quantileTime)[k]);
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}  }
}
void unstackQTime(uint  *quantileTime) {
  free_uivector(quantileTime, 1, RF_quantileSize);
}
void stackAndGetQETime(uint    treeID,
                       Node   *parent,
                       uint   *eventTimeIndex,                       
                       uint    eventTimeSize,
                       double *survival,
                       uint  **qeTimeIndex,
                       uint   *qeTimeSize) {
  uint k;
  uint uLimit;
  uint itr;
  char found;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAndGetQETime() ENTRY ...\n");
  ${trace.token}  }
  if (RF_splitRule == SURV_BSG1) {
    ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
    ${trace.token}  if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
    ${trace.token}    RF_nativePrint("\nQuantile value mapping using BSG1");
    ${trace.token}  }
    ${trace.token}  }
    *qeTimeIndex = uivector(1, RF_quantileSize);
    itr = 1;
    for (k = 1; k <= RF_quantileSize; k++) {
      found = FALSE;
      while (found == FALSE) {
        (*qeTimeIndex)[k] = itr;
        if (itr > eventTimeSize) {
          found = TRUE;
        }
        else {
          if (survival[itr] > (1.0 - RF_quantile[k])) {
            itr ++;
          }
          else {
            found = TRUE;
          }
        }
      }
      (*qeTimeIndex)[k] --;
    }
    *qeTimeSize = RF_quantileSize;
    ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
    ${trace.token}  if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
    ${trace.token}    RF_nativePrint("\nQuantile value mapping for split statistic (treeID, nodeID) = (%10d, %10d) \n", treeID, parent -> nodeID);
    ${trace.token}    RF_nativePrint("          k      quant  1 - quant       surv      qTime \n");
    ${trace.token}    for (k = 1; k <= RF_quantileSize; k++) {
    ${trace.token}      RF_nativePrint(" %10d %10.4f %10.4f %10.4f %10d \n",
    ${trace.token}         k, RF_quantile[k], 1 - RF_quantile[k], ((*qeTimeIndex)[k] == 0) ? 1.0 : survival[(*qeTimeIndex)[k]], (*qeTimeIndex)[k]);
    ${trace.token}    }
    ${trace.token}  }
    ${trace.token}  }
  }
  else {
    ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
    ${trace.token}  if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
    ${trace.token}    RF_nativePrint("\nQuantile value mapping using BSG1 or BSG2");
    ${trace.token}  }
    ${trace.token}  }
    *qeTimeIndex = uivector(1, eventTimeSize + 1);
    *qeTimeSize = 0;
    if (RF_splitRule == SURV_BSG1) {
      uLimit = (uint) ceil((double) RF_masterTimeSize * RF_quantile[1]);
      for (k = 1; k <= eventTimeSize; k++) {      
        if (eventTimeIndex[k] <= uLimit) {
          (*qeTimeIndex)[k] = k;          
          (*qeTimeSize) ++;
        }
      }
    }
    else {
      uLimit = (uint) ceil((double) eventTimeSize * RF_quantile[1]);
      for (k = 1; k <= uLimit; k++) {
        (*qeTimeIndex)[k] = k;
      }
      *qeTimeSize = uLimit;
    }
    ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
    ${trace.token}  if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
    ${trace.token}    RF_nativePrint("\nQuantile value mapping for split statistic (treeID, nodeID) = (%10d, %10d) \n", treeID, parent -> nodeID);
    ${trace.token}    RF_nativePrint("          k      qTime \n");
    ${trace.token}    for (k = 1; k <= *qeTimeSize; k++) {
    ${trace.token}      RF_nativePrint(" %10d %10d \n",
    ${trace.token}         k, (*qeTimeIndex)[k]);
    ${trace.token}    }
    ${trace.token}  }
    ${trace.token}  }
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nstackAndGetQETime() EXIT ...\n");
  ${trace.token}  }
}
void unstackQETime(uint treeID, uint eventTimeSize, uint  *qeTimeIndex) {
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nunstackQETime() ENTRY ...\n");
  ${trace.token}  }
  if (RF_splitRule == SURV_BSG1) {
    free_uivector(qeTimeIndex, 1, RF_quantileSize);
  }
  else {
    free_uivector(qeTimeIndex, 1, eventTimeSize + 1);
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\nunstackQETime() EXIT ...\n");
  ${trace.token}  }
}
