
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "termOps.h"
#include "nrutil.h"
#include "error.h"
Terminal *makeTerminal(void) {
  Terminal *parent = (Terminal*) gblock((size_t) sizeof(Terminal));
  parent -> lmiIndex      = NULL;
  parent -> lmiValue      = NULL;
  parent -> lmiSize       = 0;
  parent -> lmiAllocSize  = 0;
  parent -> nodeID        = 0;
  parent -> mate          = NULL;
  parent -> eTypeSize            = 0;
  parent -> mTimeSize            = 0;
  parent -> eTimeSize            = 0;
  parent -> sTimeSize            = 0;
  parent -> atRiskCount          = NULL;
  parent -> eventCount           = NULL;
  parent -> eventTimeIndex       = NULL;
  parent -> localRatio           = NULL;
  parent -> localCSH             = NULL;
  parent -> localCIF             = NULL;
  parent -> localSurvival        = NULL;
  parent -> localNelsonAalen     = NULL;
  parent -> CSH                  = NULL;
  parent -> CIF                  = NULL;
  parent -> survival             = NULL;
  parent -> nelsonAalen          = NULL;
  parent -> rfCount              = 0;
  parent -> rfSize               = NULL;
  parent -> multiClassProb       = NULL;
  parent -> maxClass             = NULL;
  parent -> rnfCount             = 0;
  parent -> meanResponse         = NULL;
  parent -> weight               = 0.0;
  parent -> inbagProxy           = 0;
  parent -> membrStream       = NULL;
  parent -> oobMembrSizeAlloc = 0;
  parent -> oobMembrSize      = 0;
  parent -> oobMembrIndx      = NULL;
  parent -> ibgMembrSizeAlloc = 0;
  parent -> ibgMembrSize      = 0;
  parent -> ibgMembrIndx      = NULL;
  return parent;
}
void freeTerminal(Terminal        *parent) {
  unstackTermLMIIndex(parent);
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    freeTerminalNodeSurvivalStructuresIntermediate(parent);
    freeTerminalNodeSurvivalStructuresFinal(parent);
  }
  else {
    freeTerminalNodeNonSurvivalStructures(parent);
  }
  if (parent -> oobMembrIndx != NULL) {
    if (parent -> oobMembrSizeAlloc > 0) {
      free_uivector(parent -> oobMembrIndx, 1, parent -> oobMembrSizeAlloc);
      parent -> oobMembrSize = parent -> oobMembrSizeAlloc = 0;
    }
  }
  if (parent -> ibgMembrIndx != NULL) {
    if (parent -> ibgMembrSizeAlloc > 0) {
      free_uivector(parent -> ibgMembrIndx, 1, parent -> ibgMembrSizeAlloc);
      parent -> ibgMembrSize = parent -> ibgMembrSizeAlloc = 0;
    }
  }
  free_gblock(parent, (size_t) sizeof(Terminal));
}
void freeTerminalNodeLocalSurvivalStructures(Terminal *tTerm) {
  unstackLocalRatio(tTerm);
  unstackLocalSurvival(tTerm);
  unstackLocalNelsonAalen(tTerm);
  if (tTerm -> eTypeSize > 1) {
    unstackLocalCSH(tTerm);
    unstackLocalCIF(tTerm);
  }
  unstackEventTimeIndex(tTerm);
}
void freeTerminalNodeSurvivalStructuresIntermediate(Terminal *tTerm) {
  unstackSurvival(tTerm);
  unstackNelsonAalen(tTerm);
  unstackCSH(tTerm);
  unstackCIF(tTerm);
}
void freeTerminalNodeSurvivalStructuresFinal(Terminal *tTerm) {
  unstackMortality(tTerm);
}
void freeTerminalNodeNonSurvivalStructures(Terminal *tTerm) {
  unstackMultiClassProb(tTerm);
  unstackMeanResponse(tTerm);
  unstackMemberStream(tTerm);
}
void stackAtRiskAndEventCount(Terminal *tTerm, unsigned int eTypeSize, unsigned int mTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> mTimeSize > 0) {
    if (tTerm -> mTimeSize != mTimeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  mTimeSize has been previously defined:  %10d vs %10d", tTerm -> mTimeSize, mTimeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> mTimeSize = mTimeSize;
  }
  tTerm -> atRiskCount     = uivector(1, mTimeSize);
  tTerm -> eventCount      = uimatrix(1, eTypeSize, 1, mTimeSize);
}
void unstackAtRiskAndEventCount(Terminal *tTerm) {
  if (tTerm -> atRiskCount != NULL) {
    free_uivector(tTerm -> atRiskCount, 1, tTerm -> mTimeSize);
    tTerm -> atRiskCount = NULL;
  }
  if (tTerm -> eventCount != NULL) {
    free_uimatrix(tTerm -> eventCount, 1, tTerm -> eTypeSize, 1, tTerm -> mTimeSize);
    tTerm -> eventCount = NULL;
  }
}
void stackEventTimeIndex(Terminal *tTerm, unsigned int eTimeSize) {
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> eventTimeIndex  = uivector(1, eTimeSize + 1);
}
void unstackEventTimeIndex(Terminal *tTerm) {
  if (tTerm -> eventTimeIndex != NULL) {
    free_uivector(tTerm -> eventTimeIndex, 1, tTerm -> eTimeSize + 1);
    tTerm -> eventTimeIndex = NULL;
  }
}
void stackLocalRatio(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> localRatio = dmatrix(1, eTypeSize, 1, tTerm -> eTimeSize);
}
void unstackLocalRatio(Terminal *tTerm) {
  if(tTerm -> eTimeSize > 0) {
    if (tTerm -> localRatio != NULL) {
      free_dmatrix(tTerm -> localRatio, 1, tTerm -> eTypeSize, 1, tTerm -> eTimeSize);
      tTerm -> localRatio = NULL;
    }
  }
}
void stackLocalSurvival(Terminal *tTerm, unsigned int eTimeSize) {
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> localSurvival = dvector(1, tTerm -> eTimeSize);
}
void unstackLocalSurvival(Terminal *tTerm) {
  if(tTerm -> eTimeSize > 0) {
    if (tTerm -> localSurvival != NULL) {
      free_dvector(tTerm -> localSurvival, 1, tTerm -> eTimeSize);
      tTerm -> localSurvival = NULL;
    }
  }
}
void stackLocalNelsonAalen(Terminal *tTerm, unsigned int eTimeSize) {
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> localNelsonAalen = dvector(1, tTerm -> eTimeSize);
}
void unstackLocalNelsonAalen(Terminal *tTerm) {
  if(tTerm -> eTimeSize > 0) {
    if (tTerm -> localNelsonAalen != NULL) {
      free_dvector(tTerm -> localNelsonAalen, 1, tTerm -> eTimeSize);
      tTerm -> localNelsonAalen = NULL;
    }
  }
}
void stackLocalCSH(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> localCSH = dmatrix(1, eTypeSize, 1, tTerm -> eTimeSize);
}
void unstackLocalCSH(Terminal *tTerm) {
  if(tTerm -> eTimeSize > 0) {
    if (tTerm -> localCSH != NULL) {
      free_dmatrix(tTerm -> localCSH, 1, tTerm -> eTypeSize, 1, tTerm -> eTimeSize);
      tTerm -> localCSH = NULL;
    }
  }
}
void stackLocalCIF(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> localCIF = dmatrix(1, eTypeSize, 1, tTerm -> eTimeSize);
}
void unstackLocalCIF(Terminal *tTerm) {
  if(tTerm -> eTimeSize > 0) {
    if (tTerm -> localCIF != NULL) {
      free_dmatrix(tTerm -> localCIF, 1, tTerm -> eTypeSize, 1, tTerm -> eTimeSize);
      tTerm -> localCIF = NULL;
    }
  }
}
void stackNelsonAalen(Terminal *tTerm, unsigned int sTimeSize) {
  if (tTerm -> sTimeSize > 0) {
    if (tTerm -> sTimeSize != sTimeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  sTimeSize has been previously defined:  %10d vs %10d", tTerm -> sTimeSize, sTimeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> sTimeSize = sTimeSize;
  }
  tTerm -> nelsonAalen = dvector(1, tTerm -> sTimeSize);
}
void unstackNelsonAalen(Terminal *tTerm) {
  if(tTerm -> sTimeSize > 0) {
    if (tTerm -> nelsonAalen != NULL) {
      free_dvector(tTerm -> nelsonAalen, 1, tTerm -> sTimeSize);
      tTerm -> nelsonAalen = NULL;
    }
  }
}
void stackSurvival(Terminal *tTerm, unsigned int sTimeSize) {
  if (tTerm -> sTimeSize > 0) {
    if (tTerm -> sTimeSize != sTimeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  sTimeSize has been previously defined:  %10d vs %10d", tTerm -> sTimeSize, sTimeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> sTimeSize = sTimeSize;
  }
  tTerm -> survival = dvector(1, tTerm -> sTimeSize);
}
void unstackSurvival(Terminal *tTerm) {
  if(tTerm -> sTimeSize > 0) {
    if (tTerm -> survival != NULL) {
      free_dvector(tTerm -> survival, 1, tTerm -> sTimeSize);
      tTerm -> survival = NULL;
    }
  }
}
void stackCSH(Terminal *tTerm, unsigned int eTypeSize, unsigned int sTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> sTimeSize > 0) {
    if (tTerm -> sTimeSize != sTimeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  sTimeSize has been previously defined:  %10d vs %10d", tTerm -> sTimeSize, sTimeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> sTimeSize = sTimeSize;
  }
  tTerm -> CSH = dmatrix(1, eTypeSize, 1, tTerm -> sTimeSize);
}
void unstackCSH(Terminal *tTerm) {
  if(tTerm -> eTypeSize > 0) {
    if(tTerm -> sTimeSize > 0) {
      if (tTerm -> CSH != NULL) {
        free_dmatrix(tTerm -> CSH, 1, tTerm -> eTypeSize, 1, tTerm -> sTimeSize);
        tTerm -> CSH = NULL;
      }
    }
  }
}
void stackCIF(Terminal *tTerm, unsigned int eTypeSize, unsigned int sTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> sTimeSize > 0) {
    if (tTerm -> sTimeSize != sTimeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  sTimeSize has been previously defined:  %10d vs %10d", tTerm -> sTimeSize, sTimeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> sTimeSize = sTimeSize;
  }
  tTerm -> CIF = dmatrix(1, eTypeSize, 1, tTerm -> sTimeSize);
}
void unstackCIF(Terminal *tTerm) {
  if(tTerm -> eTypeSize > 0) {
    if(tTerm -> sTimeSize > 0) {
      if (tTerm -> CIF != NULL) {
        free_dmatrix(tTerm -> CIF, 1, tTerm -> eTypeSize, 1, tTerm -> sTimeSize);
        tTerm -> CIF = NULL;
      }
    }
  }
}
void stackMortality(Terminal *tTerm, unsigned int eTypeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  tTerm -> mortality = dvector(1, eTypeSize);
}
void unstackMortality(Terminal *tTerm) {
  if(tTerm -> eTypeSize > 0) {
    if (tTerm -> mortality != NULL) {
      free_dvector(tTerm -> mortality, 1, tTerm -> eTypeSize);
      tTerm -> mortality = NULL;
    }
  }
}
void stackMultiClassProb(Terminal *tTerm, unsigned int rfCount, unsigned int *rfSize) {
  unsigned int j;
  if (tTerm -> rfCount > 0) {
    if (tTerm -> rfCount != rfCount) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  rfCount has been previously defined:  %10d vs %10d", tTerm -> rfCount, rfCount);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> rfCount = rfCount;
  }
  tTerm -> rfSize = uivector(1, tTerm -> rfCount);
  tTerm -> multiClassProb = (unsigned int **) new_vvector(1, tTerm -> rfCount, NRUTIL_UPTR);
  for (j = 1; j <= tTerm -> rfCount; j++) {
    (tTerm -> rfSize)[j] = rfSize[j];
    (tTerm -> multiClassProb)[j] = uivector(1, (tTerm -> rfSize)[j]);
  }
  tTerm -> maxClass = dvector(1, tTerm -> rfCount);
}
void stackMultiClassProbPartial(Terminal *tTerm, unsigned int rfCount) {
  if (tTerm -> rfCount > 0) {
    if (tTerm -> rfCount != rfCount) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  rfCount has been previously defined:  %10d vs %10d", tTerm -> rfCount, rfCount);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> rfCount = rfCount;
  }
  tTerm -> maxClass = dvector(1, tTerm -> rfCount);
}
void unstackMultiClassProb(Terminal *tTerm) {
  unsigned int j;
  if (tTerm -> rfCount > 0) {
    if (tTerm -> rfSize != NULL) {
      if (tTerm -> multiClassProb != NULL) {
        for (j = 1; j <= tTerm -> rfCount; j++) {
          if (tTerm -> multiClassProb[j] != NULL) {
            free_uivector(tTerm -> multiClassProb[j], 1, tTerm -> rfSize[j]);
            tTerm -> multiClassProb[j] = NULL;
          }
        }
        free_new_vvector(tTerm -> multiClassProb, 1, tTerm -> rfCount, NRUTIL_UPTR);
        tTerm -> multiClassProb = NULL;
      }
      free_uivector(tTerm -> rfSize, 1, tTerm -> rfCount);
      tTerm -> rfSize = NULL;
    }
  }
  if (tTerm -> rfCount > 0) {
    if (tTerm -> maxClass != NULL) {
      free_dvector(tTerm -> maxClass, 1, tTerm -> rfCount);
      tTerm -> maxClass = NULL;
    }
  }
}
void stackMeanResponse(Terminal *tTerm, unsigned int rnfCount) {
  if (tTerm -> rnfCount > 0) {
    if (tTerm -> rnfCount != rnfCount) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  rnfCount has been previously defined:  %10d vs %10d", tTerm -> rnfCount, rnfCount);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> rnfCount = rnfCount;
  }
  tTerm -> meanResponse = dvector(1, tTerm -> rnfCount);
}
void unstackMeanResponse(Terminal *tTerm) {
  if (tTerm -> rnfCount > 0) {
    if (tTerm -> meanResponse != NULL) {
      free_dvector(tTerm -> meanResponse, 1, tTerm -> rnfCount);
      tTerm -> meanResponse = NULL;
    }
  }
}
void stackMemberStream(Terminal *tTerm, unsigned int membrCount) {
  if (tTerm -> membrCount > 0) {
    if (tTerm -> membrCount != membrCount) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  membrSize has been previously defined:  %10d vs %10d", tTerm -> membrCount, membrCount);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tTerm -> membrCount = membrCount;
  }
  tTerm -> membrStream = uivector(1, tTerm -> membrCount);
}
void unstackMemberStream(Terminal *tTerm) {
  if (tTerm -> membrCount > 0) {
    if (tTerm -> membrStream != NULL) {
      free_uivector(tTerm -> membrStream, 1, tTerm -> membrCount);
      tTerm -> membrStream = NULL;
    }
  }
}
void stackTermLMIIndex(Terminal *tTerm, unsigned int size) {
  if (tTerm -> lmiAllocSize > 0) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  lmiIndex has been previously defined:  %10d vs %10d", tTerm -> lmiAllocSize, size);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  else {
    tTerm -> lmiAllocSize = size;
    tTerm -> lmiSize = size;
  }
  tTerm -> lmiIndex = uivector(1, tTerm -> lmiAllocSize);
  tTerm -> lmiValue = dvector(1, tTerm -> lmiAllocSize);
}
void unstackTermLMIIndex(Terminal *tTerm) {
  if(tTerm -> lmiAllocSize > 0) {
    if (tTerm -> lmiIndex != NULL) {
      free_uivector(tTerm -> lmiIndex, 1, tTerm -> lmiAllocSize);
      free_dvector(tTerm -> lmiValue, 1, tTerm -> lmiAllocSize);
      tTerm -> lmiIndex = NULL;
      tTerm -> lmiValue = NULL;
      tTerm ->lmiAllocSize = 0;
      tTerm ->lmiSize = 0;
    }
  }
}
