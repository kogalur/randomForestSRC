
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
  eventTimeSize = 0;
  for (i=1; i <= RF_masterTimeSize; i++) {
    eventTimeCount[i] = 0;
  }
  if (eventType) {
    for (i = 1; i <= nonMissMembrSize; i++) {
      if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[i]] ] > 0) {
        eventTimeCount[RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[i]] ]] ++;
      }
    }
  }
  else {
    for (i = 1; i <= nonMissMembrSize; i++) {
      if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[i]] ] == 0) {
        eventTimeCount[RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[i]] ]] ++;
      }
    }
  }
  for (i=1; i <= RF_masterTimeSize; i++) {
    if (eventTimeCount[i] > 0) {
      eventTimeIndex[++eventTimeSize] = i;
    }
  }
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
  if (genEventTimeSize > 0) {
    free_uivector(genParentEvent, 1, genEventTimeSize);
    free_uivector(genParentAtRisk, 1, genEventTimeSize);
    free_uivector(genLeftEvent, 1, genEventTimeSize);
    free_uivector(genLeftAtRisk, 1, genEventTimeSize);
    free_uivector(genRightEvent, 1, genEventTimeSize);
    free_uivector(genRightAtRisk, 1, genEventTimeSize);
  }
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
  for (i = 1; i <= eventTimeSize; i++) {
    parentAtRisk[i] = 0;
    parentEvent[i]  = eventTimeCount[eventTimeIndex[i]];
    for (j = 1; j <= nonMissMembrSize; j++) {
      if (eventTimeIndex[i] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[j]] ]) {
        parentAtRisk[i] ++;
      }
    }
  }
}
void stackAndGetSplitSurv2(uint     treeID,
                           Node    *parent,
                           uint     eventTimeSize,
                           uint    *parentEvent,
                           uint    *parentAtRisk,
                           double **localSurvival) {
  double *localRatio;
  uint q;
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
}
void unstackAndGetSplitSurv2(uint     treeID,
                             Node    *parent,
                             uint     eventTimeSize,
                             double  *localSurvival) {
  free_dvector(localSurvival, 1, eventTimeSize + 1);
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
}
void unstackFZhat(uint  treeID,
                  Node *parent,
                  uint  eventTimeSize,
                  double *fZHat) {
  free_dvector(fZHat, 1, eventTimeSize);
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
  if ((eventTimeSize > 0) && (qeTimeSize > 0)) {
    if (qeTimeIndex[qeTimeSize] > 0) {
      tauTimeIdx = eventTimeIndex[qeTimeIndex[qeTimeSize]];
    }
    else {
      tauTimeIdx = 0;
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
        if (eventTimeIndex[t] <= tauTimeIdx) {
          (*gamma_ktm)[t] = dvector(1, nonMissMembrSize);
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
  }
}
void  unstackLocalGamma(uint    treeID,
                        uint    nonMissMembrSize,
                        uint   *eventTimeIndex,
                        uint    eventTimeSize,
                           uint      *qeTimeIndex,
                        uint       qeTimeSize,
                        double **gamma_ktm) {
  if ((eventTimeSize > 0) && (qeTimeSize > 0)) {
    for (uint t = 1; t <= eventTimeSize; t++) {
      if (gamma_ktm[t] != NULL) {
        free_dvector(gamma_ktm[t], 1, nonMissMembrSize);
        gamma_ktm[t] = NULL;
      }
    }  
    free_new_vvector(gamma_ktm, 1, eventTimeSize, NRUTIL_DPTR);
  }
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
  if (RF_splitRule == SURV_BSG1) {
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
  }
  else {
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
  }
}
void unstackQETime(uint treeID, uint eventTimeSize, uint  *qeTimeIndex) {
  if (RF_splitRule == SURV_BSG1) {
    free_uivector(qeTimeIndex, 1, RF_quantileSize);
  }
  else {
    free_uivector(qeTimeIndex, 1, eventTimeSize + 1);
  }
}
