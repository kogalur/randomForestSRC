
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "survival.h"
#include "rfsrcUtil.h"
#include "termOps.h"
#include "nrutil.h"
#include "error.h"
void getAtRiskAndEventCount(uint       treeID,
                             Terminal  *parent,
                             uint      *repMembrIndx,
                             uint       repMembrSize,
                             uint      *allMembrIndx,
                             uint       allMembrSize,
                             uint      *rmbrIterator) {
  uint *membershipIndex;
  uint  membershipSize;
  uint i, j, k;
  uint ii;
  char eventFlag;
  if ( !(RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2) ) {
    membershipIndex = allMembrIndx;
    membershipSize = parent -> membrCount = allMembrSize;
    if (RF_optHigh & OPT_MEMB_INCG) {
      membershipIndex = RF_AMBR_ID_ptr[treeID];
    }
  }
  else {
    membershipIndex = repMembrIndx;
    membershipSize = parent -> membrCount = repMembrSize;
    if (RF_optHigh & OPT_MEMB_INCG) {
      membershipIndex = RF_RMBR_ID_ptr[treeID];
    }
  }
  if (membershipSize == 0) {
    if (!(RF_opt & OPT_OUTC_TYPE)) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Zero node count encountered in (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  if (!(RF_optHigh & OPT_TERM_INCG)) {
    stackAtRiskAndEventCount(parent, RF_eventTypeSize, RF_masterTimeSize);
    for (j = 1; j <= RF_masterTimeSize; j++) {
      (parent -> atRiskCount)[j] = 0;
      for (k = 1; k <= RF_eventTypeSize; k++) {
        (parent -> eventCount)[k][j] = 0;
      }
    }
    if (RF_optHigh & OPT_MEMB_OUTG) {
      for (i = 1; i <= membershipSize; i++) {
        ii = membershipIndex[i];
        RF_RMBR_ID_ptr[treeID][++(*rmbrIterator)] = ii;
        for (j = 1; j <= RF_masterTimeIndex[treeID][ii]; j++) {
          (parent -> atRiskCount)[j] ++;
        }
        if (RF_status[treeID][ii] > 0) {
          if (RF_eventTypeSize > 1) {
            k = RF_eventTypeIndex[(uint) RF_status[treeID][ii]];
          }
          else {
            k = 1;
          }
          (parent -> eventCount)[k][RF_masterTimeIndex[treeID][ii]] ++;
        }
      }
    }
    else if (RF_optHigh & OPT_MEMB_INCG) {
      for (i = 1; i <= membershipSize; i++) {
        ii = membershipIndex[++(*rmbrIterator)];
        for (j = 1; j <= RF_masterTimeIndex[treeID][ii]; j++) {
          (parent -> atRiskCount)[j] ++;
        }
        if (RF_status[treeID][ii] > 0) {
          if (RF_eventTypeSize > 1) {
            k = RF_eventTypeIndex[(uint) RF_status[treeID][ii]];
          }
          else {
            k = 1;
          }
          (parent -> eventCount)[k][RF_masterTimeIndex[treeID][ii]] ++;
        }
      }
    }
    else {
      for (i = 1; i <= membershipSize; i++) {
        ii = membershipIndex[i];
        for (j = 1; j <= RF_masterTimeIndex[treeID][ii]; j++) {
          (parent -> atRiskCount)[j] ++;
        }
        if (RF_status[treeID][ii] > 0) {
          if (RF_eventTypeSize > 1) {
            k = RF_eventTypeIndex[(uint) RF_status[treeID][ii]];
          }
          else {
            k = 1;
          }
          (parent -> eventCount)[k][RF_masterTimeIndex[treeID][ii]] ++;
        }
      }
    }
    uint *tempEventTimeIndex = uivector(1, RF_masterTimeSize);
    parent -> eTimeSize = 0;
    i = 0;    
    for (j = 1; j <= RF_masterTimeSize; j++) {
      eventFlag = FALSE;
      for (k = 1; k <= RF_eventTypeSize; k++) {
        if ((parent -> eventCount)[k][j] > 0) {
          eventFlag = TRUE;
          k = RF_eventTypeSize;
        }
      }
      if (eventFlag == TRUE) {
        tempEventTimeIndex[++i] = j;        
        (parent -> eTimeSize)++;
      }
    }
    stackEventTimeIndex(parent, parent -> eTimeSize);
    for (j = 1; j <= parent -> eTimeSize; j++) {
      (parent -> eventTimeIndex)[j] = tempEventTimeIndex[j];
    }
    free_uivector(tempEventTimeIndex, 1, RF_masterTimeSize);
  }
  else {
  }
}
void getLocalRatio(uint treeID, Terminal *parent) {
  uint j, q;
  if (parent -> membrCount > 0) {
    if(parent -> eTimeSize > 0) {
      stackLocalRatio(parent, RF_eventTypeSize, parent -> eTimeSize);
      for (j = 1; j <= RF_eventTypeSize; j++) {
        for (q = 1; q <= parent -> eTimeSize; q++) {
          if ((parent -> eventCount)[j][(parent -> eventTimeIndex)[q]] > 0) {
            if ((parent -> atRiskCount)[(parent -> eventTimeIndex)[q]] >= 1) {
              (parent -> localRatio)[j][q] = (double) ((parent -> eventCount)[j][(parent -> eventTimeIndex)[q]]) / (parent -> atRiskCount)[(parent -> eventTimeIndex)[q]];
            }
            else {
              RF_nativeError("\nRF-SRC:  *** ERROR *** ");
              RF_nativeError("\nRF-SRC:  Zero At Risk Count encountered in local ratio calculation for (tree, leaf) = (%10d, %10d)", treeID, parent -> nodeID);
              RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
              RF_nativeExit();
            }
          }
          else {
            (parent -> localRatio)[j][q] = 0.0;
          }
        }
      }
    }
  }
}
void getLocalSurvival(uint treeID, Terminal *parent) {
  uint j, q;
  if(parent -> eTimeSize > 0) {
    stackLocalSurvival(parent, parent -> eTimeSize);
    for (q = 1; q <= parent -> eTimeSize; q++) {
      (parent -> localSurvival)[q] = 0.0;
      for (j = 1; j <= RF_eventTypeSize; j++) {
        (parent -> localSurvival)[q] += (parent -> localRatio)[j][q];
      }
      (parent -> localSurvival)[q] = 1.0 - (parent -> localSurvival)[q];
    }  
    for (q = 2; q <= parent -> eTimeSize; q++) {
      (parent -> localSurvival)[q] *= (parent -> localSurvival)[q-1];
    }
  }
}
void getLocalNelsonAalen(uint treeID, Terminal *parent) {
  uint q;
  if (parent -> eTimeSize > 0) {
    stackLocalNelsonAalen(parent, parent -> eTimeSize);
    for (q = 1; q <= parent -> eTimeSize; q++) {
      (parent -> localNelsonAalen)[q] = (parent -> localRatio)[1][q];
    }
    for (q = 2; q <= parent -> eTimeSize; q++) {
      (parent -> localNelsonAalen)[q] += (parent -> localNelsonAalen)[q-1];
    }
  }
}
void getLocalCSH(uint treeID, Terminal *parent) {
  uint j, q;
    if (parent -> eTimeSize > 0) {
      stackLocalCSH(parent, RF_eventTypeSize, parent -> eTimeSize);
      for (j = 1; j <= RF_eventTypeSize; j++) {
        for (q = 1; q <= parent -> eTimeSize; q++) {
          (parent -> localCSH)[j][q] = (parent -> localRatio)[j][q];
        }
        for (q = 2; q <= parent -> eTimeSize; q++) {
          (parent -> localCSH)[j][q] += (parent -> localCSH)[j][q-1];
        }
      }
    }
}
void getLocalCIF(uint treeID, Terminal *parent) {
  uint j, q;
  if(parent -> eTimeSize > 0) {
    stackLocalCIF(parent, RF_eventTypeSize, parent -> eTimeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      (parent -> localCIF)[j][1] = (parent -> localRatio)[j][1];
      for (q = 2; q <= parent -> eTimeSize; q++) {
        (parent -> localCIF)[j][q] = (parent -> localSurvival)[q-1] * (parent -> localRatio)[j][q];
      }
      for (q = 2; q <= parent -> eTimeSize; q++) {
        (parent -> localCIF)[j][q] += (parent -> localCIF)[j][q-1];
      }
    }
  }
}
void mapLocalToTimeInterest(uint      treeID,
                            Terminal *parent,
                            void     *genericLocal,
                            void     *genericGlobal) {
  uint itIndex, etIndex, lookAheadIndex;
  char mapFlag, transitFlag;
  uint j;
  if (!(RF_opt & OPT_COMP_RISK)) {
    if ((parent -> eTimeSize) > 0) {
      itIndex = 1;
      etIndex = 1;
      mapFlag = TRUE;
      while(mapFlag) {
        if (RF_timeInterest[itIndex] < RF_masterTime[(parent -> eventTimeIndex)[etIndex]] ) {
          if (itIndex > 1) {
            ((double *) genericGlobal)[itIndex] = ((double *) genericGlobal)[itIndex-1];
          }
          itIndex++;
        }
        else {
          lookAheadIndex = etIndex;
          transitFlag = TRUE;
          while (transitFlag) {
            if (RF_timeInterest[itIndex] >= RF_masterTime[(parent -> eventTimeIndex)[lookAheadIndex]] ) {
              ((double *) genericGlobal)[itIndex] = ((double *) genericLocal)[lookAheadIndex];
              lookAheadIndex++;
              if (lookAheadIndex > (parent -> eTimeSize)) {
                transitFlag = FALSE;
              }
            }
            else {
              transitFlag = FALSE;
            }
          }
          itIndex++;
          etIndex = lookAheadIndex;
        }
        if(etIndex > (parent -> eTimeSize)) {
          while(itIndex <= RF_sortedTimeInterestSize) {
            ((double *) genericGlobal)[itIndex] = ((double *) genericGlobal)[itIndex-1];
            itIndex++;
          }
        }
        if(itIndex > RF_sortedTimeInterestSize) {
          mapFlag = FALSE;
        }
      }
    }
  }  
  else {
    if ((parent -> eTimeSize) > 0) {
      itIndex = 1;
      etIndex = 1;
      mapFlag = TRUE;
      while(mapFlag) {
        if (RF_timeInterest[itIndex] < RF_masterTime[(parent -> eventTimeIndex)[etIndex]] ) {
          if (itIndex > 1) {
            for (j = 1; j <= RF_eventTypeSize; j++) {
              ((double **) genericGlobal)[j][itIndex] = ((double **) genericGlobal)[j][itIndex-1];
            }
          }
          itIndex++;
        }
        else {
          lookAheadIndex = etIndex;
          transitFlag = TRUE;
          while (transitFlag) {
            if (RF_timeInterest[itIndex] >= RF_masterTime[(parent -> eventTimeIndex)[lookAheadIndex]] ) {
              for (j = 1; j <= RF_eventTypeSize; j++) {
                ((double **) genericGlobal)[j][itIndex] = ((double **) genericLocal)[j][lookAheadIndex];
              }
              lookAheadIndex++;
              if (lookAheadIndex > (parent -> eTimeSize)) {
                transitFlag = FALSE;
              }
            }
            else {
              transitFlag = FALSE;
            }
          }
          itIndex++;
          etIndex = lookAheadIndex;
        }
        if(etIndex > (parent -> eTimeSize)) {
          while(itIndex <= RF_sortedTimeInterestSize) {
              for (j = 1; j <= RF_eventTypeSize; j++) {
                ((double **) genericGlobal)[j][itIndex] = ((double **) genericGlobal)[j][itIndex-1];
              }
              itIndex++;
          }
        }
        if(itIndex > RF_sortedTimeInterestSize) {
          mapFlag = FALSE;
        }
      }
    }    
  }  
}
void getSurvival(uint treeID, Terminal *parent) {
  uint k;
  if (!(RF_optHigh & OPT_TERM_INCG)) {
    stackSurvival(parent, RF_sortedTimeInterestSize);
    for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
      (parent -> survival)[k] = 1.0;
    }
    mapLocalToTimeInterest(treeID,
                           parent,
                           parent -> localSurvival,
                           parent -> survival);
  }
  else {
    stackSurvival(parent, RF_sortedTimeInterestSize);
    for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
      (parent -> survival)[k] = RF_TN_SURV_ptr[treeID][parent -> nodeID][k];
    }
  }
}
void getNelsonAalen(uint treeID, Terminal *parent) {
  uint k;
  if (!(RF_optHigh & OPT_TERM_INCG)) {
    stackNelsonAalen(parent, RF_sortedTimeInterestSize);
    for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
      (parent -> nelsonAalen)[k] = 0.0;
    }
    mapLocalToTimeInterest(treeID,
                           parent,
                           parent -> localNelsonAalen,
                           parent -> nelsonAalen);
  }
  else {
    stackNelsonAalen(parent, RF_sortedTimeInterestSize);
    for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
      (parent -> nelsonAalen)[k] = RF_TN_NLSN_ptr[treeID][parent -> nodeID][k];
    }
  }
}
void getCSH(uint treeID, Terminal *parent) {
  uint j, k;
  if (!(RF_optHigh & OPT_TERM_INCG)) {
    stackCSH(parent, RF_eventTypeSize, RF_sortedTimeInterestSize);
    for (j=1; j <= RF_eventTypeSize; j++) {
      for (k=1; k <= RF_sortedTimeInterestSize; k++) {
        (parent -> CSH)[j][k] = 0.0;
      }
    }
    mapLocalToTimeInterest(treeID,
                           parent,
                           parent -> localCSH,
                           parent -> CSH);
  }
  else {
    stackCSH(parent, RF_eventTypeSize, RF_sortedTimeInterestSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
        (parent -> CSH)[j][k] = RF_TN_CSHZ_ptr[treeID][parent -> nodeID][j][k];
      }
    }
  }
}
void getCIF(uint treeID, Terminal *parent) {
  uint j, k;
  if (!(RF_optHigh & OPT_TERM_INCG)) {
    stackCIF(parent, RF_eventTypeSize, RF_sortedTimeInterestSize);
    for (j=1; j <= RF_eventTypeSize; j++) {
      for (k=1; k <= RF_sortedTimeInterestSize; k++) {
        (parent -> CIF)[j][k] = 0.0;
      }
    }
    mapLocalToTimeInterest(treeID,
                           parent,
                           parent -> localCIF,
                           parent -> CIF);
  }
  else {
    stackCIF(parent, RF_eventTypeSize, RF_sortedTimeInterestSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
        (parent -> CIF)[j][k] = RF_TN_CIFN_ptr[treeID][parent -> nodeID][j][k];
      }
    }
  }
}
void getMortality(uint treeID, Terminal *parent) {
  uint j, q;
  if (!(RF_optHigh & OPT_TERM_INCG)) {
    stackMortality(parent, RF_eventTypeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      (parent -> mortality)[j] = 0.0;
    }
    if (!(RF_opt & OPT_COMP_RISK)) {
      for (q = 1; q <= RF_sortedTimeInterestSize; q++) {
        (parent -> mortality)[1] += (parent -> nelsonAalen)[q];
      }
    }
    else {
      for (j = 1; j <= RF_eventTypeSize; j ++) {
        for (q = 1; q <= RF_sortedTimeInterestSize - 1; q++) {
          (parent -> mortality)[j] += (parent -> CIF)[j][q] * (RF_timeInterest[q+1] - RF_timeInterest[q]);
        }
      }
    }
  }
  else {
    stackMortality(parent, RF_eventTypeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      (parent -> mortality)[j] = RF_TN_MORT_ptr[treeID][parent -> nodeID][j];
    }
  }
}
