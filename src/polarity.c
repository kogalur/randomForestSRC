
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "polarity.h"
#include "factorOps.h"
#include "nrutil.h"
${trace.token} #include "error.h"
char getDaughterPolarity(uint treeID, SplitInfo *info, uint indv, void *value, ...) {
  char (*getDaughterPolarityGeneric) (uint       treeID,
                                      SplitInfo *info,
                                      uint       indv,
                                      void      *value,
                                      ...);
  void *obsLocal;
  char daughterFlag;
    obsLocal = ((double **) value)[info -> randomVar[1]];
    if (info -> mwcpSizeAbs[1] > 0) {
      getDaughterPolarityGeneric = &getDaughterPolaritySimpleFactor;
    }
    else {
      getDaughterPolarityGeneric = &getDaughterPolaritySimpleNonFactor;
    }
  daughterFlag = getDaughterPolarityGeneric(0, info, indv, obsLocal);
  return daughterFlag;
}
char getDaughterPolaritySimpleFactor(uint treeID, SplitInfo *info, uint indv, void *value, ...) {
  char daughterFlag;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetDaughterPolaritySimpleFactor(%10d) ENTRY ...\n", treeID);
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}        if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}          RF_nativePrint("\nNon-Greedy Daughter Comparison (value, const):  (%10d, ", (uint) ((double *) value)[indv]);
  ${trace.token}          for (uint m = 1; m <= info -> mwcpSizeAbs[1]; m++) {
  ${trace.token}            RF_nativePrint(" %10x", ((uint*) info -> randomPts[1])[m]);
  ${trace.token}          }
  ${trace.token}          RF_nativePrint(")");
  ${trace.token}        }
  ${trace.token}      }
  daughterFlag = splitOnFactor((uint) ((double *) value)[indv], (uint*) info -> randomPts[1]);
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetDaughterPolaritySimpleFactor(%10d) EXIT ...\n", treeID);
  ${trace.token}    }
  ${trace.token}  }
  return daughterFlag;
}
char getDaughterPolaritySimpleNonFactor(uint treeID, SplitInfo *info, uint indv, void *value, ...) {
  char daughterFlag;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetDaughterPolaritySimpleNonFactor(%10d) ENTRY ...\n", treeID);
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}        if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}          RF_nativePrint("\nNon-Greedy Daughter Comparison (value, const):  (%10.4f, %10.4f)", ((double *) value)[indv], ((double*) info -> randomPts[1])[1]);
  ${trace.token}        }
  ${trace.token}      }
  daughterFlag =  (( ((double*) info -> randomPts[1])[1] - ((double *) value)[indv]) >= 0.0) ? LEFT : RIGHT;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetDaughterPolaritySimpleNonFactor(%10d) EXIT ...\n", treeID);
  ${trace.token}    }
  ${trace.token}  }
  return daughterFlag;
}
char getDaughterPolaritySimpleFactorSingle(uint treeID, SplitInfo *info, uint indv, void *value, ...) {
  char daughterFlag;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetDaughterPolaritySimpleFactorSingle(%10d) ENTRY ...\n", treeID);
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}        if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
  ${trace.token}          RF_nativePrint("\nNon-Greedy Daughter Comparison (value, const):  (%10d, ", *((uint *) ((double *) value)));
  ${trace.token}          for (uint m = 1; m <= info -> mwcpSizeAbs[1]; m++) {
  ${trace.token}            RF_nativePrint(" %10x", ((uint*) info -> randomPts[1])[m]);
  ${trace.token}          }
  ${trace.token}          RF_nativePrint(")");
  ${trace.token}        }
  ${trace.token}      }
  daughterFlag = splitOnFactor(*((uint *) ((double *) value)), (uint*) info -> randomPts[1]);
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetDaughterPolaritySimpleFactorSingle(%10d) EXIT ...\n", treeID);
  ${trace.token}    }
  ${trace.token}  }
  return daughterFlag;
}
char getDaughterPolaritySimpleNonFactorSingle(uint treeID, SplitInfo *info, uint indv, void *value, ...) {
  char daughterFlag;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetDaughterPolaritySimpleNonFactorSingle(%10d) ENTRY ...\n", treeID);
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}      if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}        if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
  ${trace.token}          RF_nativePrint("\nNon-Greedy Daughter Comparison (value, const):  (%10.4f, %10.4f)", *((double *) value), ((double*) info -> randomPts[1])[1]);
  ${trace.token}        }
  ${trace.token}      }
  daughterFlag =  (( ((double*) info -> randomPts[1])[1] - (*((double *) value))) >= 0.0) ? LEFT : RIGHT;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    if (getTraceFlag(treeID) & TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetDaughterPolaritySimpleNonFactorSingle(%10d) EXIT ...\n", treeID);
  ${trace.token}    }
  ${trace.token}  }
  return daughterFlag;
}
