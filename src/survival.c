
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
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetAtRiskAndEventCount() ENTRY ...\n");
  ${trace.token}  }
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
    ${trace.token}        if (getTraceFlag(treeID) & ENSB_HGH_TRACE) {
    ${trace.token}          if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
    ${trace.token}            RF_nativePrint("\nNode Specific At Risk and Event Counts for (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
    ${trace.token}            RF_nativePrint("Mstr Time ");
    ${trace.token}            for (i=1; i <= RF_masterTimeSize; i++) {
    ${trace.token}              RF_nativePrint("%10d", i);
    ${trace.token}            }
    ${trace.token}            RF_nativePrint("\n");
    ${trace.token}            RF_nativePrint("At Risk   ");
    ${trace.token}            for (i=1; i <= RF_masterTimeSize; i++) {
    ${trace.token}              RF_nativePrint("%10d", (parent -> atRiskCount)[i]);
    ${trace.token}            }
    ${trace.token}            RF_nativePrint("\n");
    ${trace.token}            for (j=1; j <= RF_eventTypeSize; j++) {
    ${trace.token}              RF_nativePrint("Ev %7d", j);
    ${trace.token}              for (i=1; i <= RF_masterTimeSize; i++) {
    ${trace.token}                RF_nativePrint("%10d", (parent -> eventCount)[j][i]);
    ${trace.token}              }
    ${trace.token}              RF_nativePrint("\n");
    ${trace.token}            }
    ${trace.token}            RF_nativePrint("\nEv Tm Idx ");
    ${trace.token}            for (i=1; i <= parent -> eTimeSize; i++) {
    ${trace.token}              RF_nativePrint("%10d", i);
    ${trace.token}            }
    ${trace.token}            RF_nativePrint("\n          ");
    ${trace.token}            for (i=1; i <= parent -> eTimeSize; i++) {
    ${trace.token}              RF_nativePrint("%10d", (parent -> eventTimeIndex)[i]);
    ${trace.token}            }
    ${trace.token}            RF_nativePrint("\n");
    ${trace.token}          }
    ${trace.token}        }
  }
  else {
  }
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetAtRiskAndEventCount() EXIT ...\n");
  ${trace.token}  }
}
void getLocalRatio(uint treeID, Terminal *parent) {
  uint j, q;
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetLocalRatio() ENTRY ...\n");
  ${trace.token}  }
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
    ${trace.token}      if (getTraceFlag(treeID) & ENSB_HGH_TRACE) {
    ${trace.token}        RF_nativePrint("\nNode specific local ratios:  [RF_eventTypeSize] x [nodeEventTimeSize] for:  (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
    ${trace.token}        RF_nativePrint("              mTimIdx       time ");
    ${trace.token}        for (j=1; j <= RF_eventTypeSize; j++) {
    ${trace.token}          RF_nativePrint("%10d ", j);
    ${trace.token}        }
    ${trace.token}        RF_nativePrint("\n");
    ${trace.token}        for (q=1; q <= parent -> eTimeSize; q++) {
    ${trace.token}          RF_nativePrint("%10d %10d %10.4f ", q, (parent -> eventTimeIndex)[q], RF_masterTime[(parent -> eventTimeIndex)[q]]);
    ${trace.token}          for (j=1; j <= RF_eventTypeSize; j++) {
    ${trace.token}            RF_nativePrint("%10.4f ", (parent -> localRatio)[j][q]);
    ${trace.token}          }
    ${trace.token}          RF_nativePrint("\n");
    ${trace.token}        }
    ${trace.token}      }
  }
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetLocalRatio() EXIT ...\n");
  ${trace.token}  }
}
void getLocalSurvival(uint treeID, Terminal *parent) {
  uint j, q;
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetLocalSurvival() ENTRY ...\n");
  ${trace.token}  }
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
    ${trace.token}      if (getTraceFlag(treeID) & ENSB_HGH_TRACE) {
    ${trace.token}        RF_nativePrint("\nNode specific local survival function of length [parent -> eTimeSize] for:  (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
    ${trace.token}        RF_nativePrint("              mTimIdx       time   survival \n");
    ${trace.token}        for (q = 1; q <= parent -> eTimeSize; q++) {
    ${trace.token}          RF_nativePrint("%10d %10d %10.4f %10.4f \n", q, (parent -> eventTimeIndex)[q], RF_masterTime[(parent -> eventTimeIndex)[q]], (parent -> localSurvival)[q]);
    ${trace.token}        }
    ${trace.token}      }
  }
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetLocalSurvival() EXIT ...\n");
  ${trace.token}  }
}
void getLocalNelsonAalen(uint treeID, Terminal *parent) {
  uint q;
  ${trace.token}  uint j;
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetLocalNelsonAalen() ENTRY ...\n");
  ${trace.token}  }
  if (parent -> eTimeSize > 0) {
    stackLocalNelsonAalen(parent, parent -> eTimeSize);
    for (q = 1; q <= parent -> eTimeSize; q++) {
      (parent -> localNelsonAalen)[q] = (parent -> localRatio)[1][q];
    }
    for (q = 2; q <= parent -> eTimeSize; q++) {
      (parent -> localNelsonAalen)[q] += (parent -> localNelsonAalen)[q-1];
    }
  }
  ${trace.token}  if (getTraceFlag(treeID) & ENSB_HGH_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
  ${trace.token}        RF_nativePrint("\nLocal Nelson-Aalen estimator for (tree, leaf):  (%10d, %10d) \n", treeID, parent -> nodeID);
  ${trace.token}        for (j=1; j <= parent -> eTimeSize; j++) {
  ${trace.token}          RF_nativePrint("%10d", j);
  ${trace.token}        }
  ${trace.token}        RF_nativePrint("\n");
  ${trace.token}        for (j=1; j <= parent -> eTimeSize; j++) {
  ${trace.token}          RF_nativePrint("%10.4f", parent -> localNelsonAalen[j]);
  ${trace.token}        }
  ${trace.token}        RF_nativePrint("\n");
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetLocalNelsonAalen() EXIT ...\n");
  ${trace.token}  }
}
void getLocalCSH(uint treeID, Terminal *parent) {
  uint j, q;
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetLocalCSH() ENTRY ...\n");
  ${trace.token}  }
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
      ${trace.token}      if (getTraceFlag(treeID) & ENSB_HGH_TRACE) {
      ${trace.token}        RF_nativePrint("\nNode specific local CSH:  [RF_eventTypeSize] x [nodeEventTimeSize] for:  (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
      ${trace.token}        RF_nativePrint("              mTimIdx       time ");
      ${trace.token}        for (j=1; j <= RF_eventTypeSize; j++) {
      ${trace.token}          RF_nativePrint("%10d ", j);
      ${trace.token}        }
      ${trace.token}        RF_nativePrint("\n");
      ${trace.token}        for (q=1; q <= parent -> eTimeSize; q++) {
      ${trace.token}          RF_nativePrint("%10d %10d %10.4f ", q, (parent -> eventTimeIndex)[q], RF_masterTime[(parent -> eventTimeIndex)[q]]);
      ${trace.token}          for (j=1; j <= RF_eventTypeSize; j++) {
      ${trace.token}            RF_nativePrint("%10.4f ", (parent -> localCSH)[j][q]);
      ${trace.token}          }
      ${trace.token}          RF_nativePrint("\n");
      ${trace.token}        }
      ${trace.token}      }
    }
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetLocalCSH() EXIT ...\n");
  ${trace.token}  }
}
void getLocalCIF(uint treeID, Terminal *parent) {
  uint j, q;
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetLocalCIF() ENTRY ...\n");
  ${trace.token}  }
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
    ${trace.token}      if (getTraceFlag(treeID) & ENSB_HGH_TRACE) {
    ${trace.token}        RF_nativePrint("\nNode specific local CIF:  [RF_eventTypeSize] x [nodeEventTimeSize] for:  (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
    ${trace.token}        RF_nativePrint("              mTimIdx       time ");
    ${trace.token}        for (j=1; j <= RF_eventTypeSize; j++) {
    ${trace.token}          RF_nativePrint("%10d ", j);
    ${trace.token}        }
    ${trace.token}        RF_nativePrint("\n");
    ${trace.token}        for (q=1; q <= parent -> eTimeSize; q++) {
    ${trace.token}          RF_nativePrint("%10d %10d %10.4f ", q, (parent -> eventTimeIndex)[q], RF_masterTime[(parent -> eventTimeIndex)[q]]);
    ${trace.token}          for (j=1; j <= RF_eventTypeSize; j++) {
    ${trace.token}            RF_nativePrint("%10.4f ", (parent -> localCIF)[j][q]);
    ${trace.token}          }
    ${trace.token}          RF_nativePrint("\n");
    ${trace.token}        }
    ${trace.token}      }
  }
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetLocalCIF() EXIT ...\n");
  ${trace.token}  }
}
void mapLocalToTimeInterest(uint      treeID,
                            Terminal *parent,
                            void     *genericLocal,
                            void     *genericGlobal) {
  uint itIndex, etIndex, lookAheadIndex;
  char mapFlag, transitFlag;
  uint j;
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nmapTimeInterest() ENTRY ...\n");
  ${trace.token}  }
  if (!(RF_opt & OPT_COMP_RISK)) {
    if ((parent -> eTimeSize) > 0) {
      itIndex = 1;
      etIndex = 1;
      mapFlag = TRUE;
      ${trace.token}  if (getTraceFlag(treeID) & ENSB_HGH_TRACE) {
      ${trace.token}    RF_nativePrint("\n Map:  (treeID, nodeID) = (%10d, %10d)", treeID, parent -> nodeID);
      ${trace.token}    RF_nativePrint("\n Size: (itSize, etSize) = (%10d, %10d)", RF_sortedTimeInterestSize, parent -> eTimeSize);
      ${trace.token}  }
      while(mapFlag) {
        if (RF_timeInterest[itIndex] < RF_masterTime[(parent -> eventTimeIndex)[etIndex]] ) {
          if (itIndex > 1) {
            ${trace.token}  if (getTraceFlag(treeID) & ENSB_HGH_TRACE) {
            ${trace.token}  RF_nativePrint("\n Flat-Line index:   (%10d, %10d) (%10.4f, %10.4f)", itIndex, etIndex, RF_timeInterest[itIndex], RF_masterTime[(parent -> eventTimeIndex)[etIndex]]);
            ${trace.token}  }
            ((double *) genericGlobal)[itIndex] = ((double *) genericGlobal)[itIndex-1];
          }
          itIndex++;
        }
        else {
          lookAheadIndex = etIndex;
          transitFlag = TRUE;
          while (transitFlag) {
            if (RF_timeInterest[itIndex] >= RF_masterTime[(parent -> eventTimeIndex)[lookAheadIndex]] ) {
              ${trace.token}  if (getTraceFlag(treeID) & ENSB_HGH_TRACE) {
              ${trace.token}  RF_nativePrint("\n Transition index:  (%10d, %10d) (%10.4f, %10.4f)", itIndex, lookAheadIndex, RF_timeInterest[itIndex], RF_masterTime[(parent -> eventTimeIndex)[lookAheadIndex]]);
              ${trace.token}  }
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
            ${trace.token}  if (getTraceFlag(treeID) & ENSB_HGH_TRACE) {
            ${trace.token}  RF_nativePrint("\n Tail index:        (%10d,           ) (%10.4f,           )", itIndex, RF_timeInterest[itIndex]);
            ${trace.token}  }
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
      ${trace.token}  if (getTraceFlag(treeID) & ENSB_HGH_TRACE) {
      ${trace.token}    RF_nativePrint("\n Map:  (treeID, nodeID) = (%10d, %10d)", treeID, parent -> nodeID);
      ${trace.token}    RF_nativePrint("\n Size: (itSize, etSize) = (%10d, %10d)", RF_sortedTimeInterestSize, parent -> eTimeSize);
      ${trace.token}  }
      while(mapFlag) {
        if (RF_timeInterest[itIndex] < RF_masterTime[(parent -> eventTimeIndex)[etIndex]] ) {
          if (itIndex > 1) {
            ${trace.token}  if (getTraceFlag(treeID) & ENSB_HGH_TRACE) {
            ${trace.token}  RF_nativePrint("\n Flat-Line index:   (%10d, %10d) (%10.4f, %10.4f)", itIndex, etIndex, RF_timeInterest[itIndex], RF_masterTime[(parent -> eventTimeIndex)[etIndex]]);
            ${trace.token}  }
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
              ${trace.token}  if (getTraceFlag(treeID) & ENSB_HGH_TRACE) {
              ${trace.token}  RF_nativePrint("\n Transition index:  (%10d, %10d) (%10.4f, %10.4f)", itIndex, lookAheadIndex, RF_timeInterest[itIndex], RF_masterTime[(parent -> eventTimeIndex)[lookAheadIndex]]);
              ${trace.token}  }
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
            ${trace.token}  if (getTraceFlag(treeID) & ENSB_HGH_TRACE) {
            ${trace.token}  RF_nativePrint("\n Tail index:        (%10d,           ) (%10.4f,           )", itIndex, RF_timeInterest[itIndex]);
            ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nmapTimeInterest() EXIT ...\n");
  ${trace.token}  }
}
void getSurvival(uint treeID, Terminal *parent) {
  uint k;
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetSurvival() ENTRY ...\n");
  ${trace.token}  }
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
  ${trace.token}      if (getTraceFlag(treeID) & ENSB_LOW_TRACE) {
  ${trace.token}      if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}        RF_nativePrint("\nNode specific non-local survival function [RF_sortedTimeInterestSize] x [RF_tLeafCount[treeID]] for:  (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
  ${trace.token}        RF_nativePrint("                 time ");
  ${trace.token}        RF_nativePrint("\n");
  ${trace.token}        for (k=1; k <= RF_sortedTimeInterestSize; k++) {
  ${trace.token}          RF_nativePrint("%10d %10.4f %10.4f", k, RF_timeInterest[k], (parent -> survival)[k]);
  ${trace.token}          RF_nativePrint("\n");
  ${trace.token}        }
  ${trace.token}      }
  ${trace.token}      }
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetSurvival() EXIT ...\n");
  ${trace.token}  }
}
void getNelsonAalen(uint treeID, Terminal *parent) {
  uint k;
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetNelsonAalen() ENTRY ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(treeID) & ENSB_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nNelson-Aalen estimator matrix for (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
  ${trace.token}      for (k=1; k <= RF_sortedTimeInterestSize; k++) {
  ${trace.token}        RF_nativePrint("%10d", k);
  ${trace.token}      }
  ${trace.token}      RF_nativePrint("\n");  
  ${trace.token}      for (k=1; k <= RF_sortedTimeInterestSize; k++) {
  ${trace.token}          RF_nativePrint("%10.4f", parent -> nelsonAalen[k]);
  ${trace.token}      }
  ${trace.token}      RF_nativePrint("\n");
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetNelsonAalen() EXIT ...\n");
  ${trace.token}  }
}
void getCSH(uint treeID, Terminal *parent) {
  uint j, k;
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetCSH() ENTRY ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(treeID) & ENSB_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}        RF_nativePrint("\nNode specific non-local CSH:  [RF_eventTypeSize] x [nodeEventTimeSize] for:  (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
  ${trace.token}        RF_nativePrint("                 time ");
  ${trace.token}        for (j=1; j <= RF_eventTypeSize; j++) {
  ${trace.token}          RF_nativePrint("%10d ", j);
  ${trace.token}        }
  ${trace.token}        RF_nativePrint("\n");
  ${trace.token}        for (k=1; k <= RF_sortedTimeInterestSize; k++) {
  ${trace.token}          RF_nativePrint("%10d %10.4f ", k, RF_timeInterest[k]);
  ${trace.token}          for (j=1; j <= RF_eventTypeSize; j++) {
  ${trace.token}            RF_nativePrint("%10.4f ", (parent -> CSH)[j][k]);
  ${trace.token}          }
  ${trace.token}          RF_nativePrint("\n");
  ${trace.token}        }
  ${trace.token}      }
  ${trace.token}      }
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetCSH() EXIT ...\n");
  ${trace.token}  }
}
void getCIF(uint treeID, Terminal *parent) {
  uint j, k;
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetCIF() ENTRY ...\n");
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(treeID) & ENSB_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}        RF_nativePrint("\nNode specific non-local CIF:  [RF_eventTypeSize] x [nodeEventTimeSize] for:  (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
  ${trace.token}        RF_nativePrint("                 time ");
  ${trace.token}        for (j=1; j <= RF_eventTypeSize; j++) {
  ${trace.token}          RF_nativePrint("%10d ", j);
  ${trace.token}        }
  ${trace.token}        RF_nativePrint("\n");
  ${trace.token}        for (k=1; k <= RF_sortedTimeInterestSize; k++) {
  ${trace.token}          RF_nativePrint("%10d %10.4f ", k, RF_timeInterest[k]);
  ${trace.token}          for (j=1; j <= RF_eventTypeSize; j++) {
  ${trace.token}            RF_nativePrint("%10.4f ", (parent -> CIF)[j][k]);
  ${trace.token}          }
  ${trace.token}          RF_nativePrint("\n");
  ${trace.token}        }
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetCIF() EXIT ...\n");
  ${trace.token}  }
}
void getMortality(uint treeID, Terminal *parent) {
  uint j, q;
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetMortality() ENTRY ...\n");
  ${trace.token}  }
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
  ${trace.token}      if (getTraceFlag(treeID) & ENSB_LOW_TRACE) {
  ${trace.token}      if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}        RF_nativePrint("\nNode specific mortality:  [RF_eventTypeSize] for:  (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
  ${trace.token}        for (j=1; j <= RF_eventTypeSize; j++) {
  ${trace.token}          RF_nativePrint("%10d ", j);
  ${trace.token}        }
  ${trace.token}        RF_nativePrint("\n");
  ${trace.token}        for (j=1; j <= RF_eventTypeSize; j++) {
  ${trace.token}          RF_nativePrint("%10.4f ", (parent -> mortality)[j]);
  ${trace.token}        }
  ${trace.token}        RF_nativePrint("\n");
  ${trace.token}      }
  ${trace.token}      }
  ${trace.token}  if (getTraceFlag(treeID) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetMortality() EXIT ...\n");
  ${trace.token}  }
}
