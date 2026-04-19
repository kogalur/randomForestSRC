
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "importance.h"
#include "importanceAnti.h"
#include "importancePerm.h"
#include "importanceRand.h"
#include "survivalE.h"
#include "rfsrcUtil.h"
#include "polarity.h"
#include "nrutil.h"
#include "error.h"
Node *identifyExtrapolatedMembership (Node     *parent,
                                      double **yShadow,
                                      double **xShadow) {
  void *obsLocal;
  char daughterFlag;
  Node *result = parent;
  SplitInfo *info;
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\nidentifyExtrapolatedMembership() ENTRY:  (nodeID, depth)  = (%10d, %10d)", parent -> nodeID, parent -> depth);
  ${trace.token}    }
  ${trace.token}  }
  daughterFlag = NEITHER;  
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    info = parent -> splitInfo;
    if (info -> randomVar[1] > 0) {
      obsLocal = xShadow;
    }
    else {
      obsLocal = yShadow;
    }
    daughterFlag = getDaughterPolarity(0, info, 1, obsLocal);    
    ${trace.token}      if (getTraceFlag(0) & VIMP_LOW_TRACE) {
    ${trace.token}        if (getTraceFlag(0) & !TURN_OFF_TRACE) {
    ${trace.token}          if(daughterFlag == LEFT) {
    ${trace.token}            RF_nativePrint("\nNon-Greedy Daughter LEFT :  %10d ", 0);
    ${trace.token}          }
    ${trace.token}          else {
    ${trace.token}            RF_nativePrint("\nNon-Greedy Daughter RGHT :  %10d ", 0);
    ${trace.token}          }
    ${trace.token}        }
    ${trace.token}      }
    if (daughterFlag == LEFT) {
      result = identifyExtrapolatedMembership(parent ->  left, yShadow, xShadow);
    }
    else {
      result = identifyExtrapolatedMembership(parent -> right, yShadow, xShadow);
    }
  }
  ${trace.token}  if (getTraceFlag(0) & TURN_OFF_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}      RF_nativePrint("\nidentifyExtrapolatedMembership() EXIT:   (nodeID, depth)  = (%10d, %10d)", parent -> nodeID, parent -> depth);
  ${trace.token}    }
  ${trace.token}  }
  return result;
}
void getVimpMembership (char       mode,
                        uint       treeID,
                        Terminal **vimpMembership,
                        uint p) {
  char result;
  ${trace.token}  uint  obsSize;
  ${trace.token}  char *membershipFlag;
  ${trace.token}  uint  j;
  ${trace.token}  if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetVimpMembership() ENTRY ...\n");
  ${trace.token}    }
  ${trace.token}  }
  if (!(RF_opt & OPT_VIMP)) {
    RF_nativePrint("\nRF-SRC:  *** ERROR *** ");
    RF_nativePrint("\nRF-SRC:  Attempt to compute variable importance though not requested.");
    RF_nativePrint("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  result = FALSE;
  switch (mode) {
  case RF_PRED:
    result = TRUE;
    break;
  default:
    if (RF_oobSize[treeID] > 0) {
      result = TRUE;
    }
    break;
  }
  ${trace.token}  switch (mode) {
  ${trace.token}  case RF_PRED:
  ${trace.token}  obsSize = RF_fobservationSize;
  ${trace.token}    membershipFlag = RF_testMembershipFlag;
  ${trace.token}    break;
  ${trace.token}  default:
  ${trace.token}  obsSize = RF_observationSize;
  ${trace.token}    membershipFlag = RF_bootMembershipFlag[treeID];
  ${trace.token}    break;
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
  ${trace.token}  if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}    RF_nativePrint("\n\nVariable Importance boot flag for tree:  %10d", treeID);
  ${trace.token}    RF_nativePrint("\n     index       flag \n");
  ${trace.token}    for (j=1; j <=  obsSize; j++) {
  ${trace.token}      RF_nativePrint("%10d %10d \n", j, membershipFlag[j]);
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}  }
  if (result == TRUE) {
    if (!(RF_opt & OPT_VIMP_TYP1) && !(RF_opt & OPT_VIMP_TYP2)) {
      getAntiMembership(mode, treeID, vimpMembership, p);
    }
    else if ((RF_opt & OPT_VIMP_TYP1) && !(RF_opt & OPT_VIMP_TYP2)) { 
      getPermuteMembership(mode, treeID, vimpMembership, p);
    }
    else if (!(RF_opt & OPT_VIMP_TYP1) && (RF_opt & OPT_VIMP_TYP2)) { 
      getRandomMembership(mode, treeID, vimpMembership, p);
    }
    else {
      RF_nativePrint("\nRF-SRC:  *** ERROR *** ");
      RF_nativePrint("\nRF-SRC:  Unknown VIMP type encountered:  %10d", RF_opt);
      RF_nativePrint("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    ${trace.token}  if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
    ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
    ${trace.token}      RF_nativePrint("\nVIMP omitted since OOB sample size is zero. \n");
    ${trace.token}    }
    ${trace.token}  }
  }
  ${trace.token}  if (getTraceFlag(treeID) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}      RF_nativePrint("\ngetVimpMembership() EXIT ...\n");
  ${trace.token}    }
  ${trace.token}  }
}
void updateEnsembleVimp (char       mode,
                         uint       treeID,
                         Terminal **noiseMembership,
                         uint       xVarIdx) {
  uint  *membershipIndex;
  uint   membershipSize;
  double *denomPtr;
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdateEnsembleVimp() ENTRY ...\n");
  ${trace.token}  }
  switch (mode) {
  case RF_PRED:
    membershipSize = RF_fobservationSize;
    membershipIndex = RF_fidentityMembershipIndex;
    break;
  default:
    membershipSize  = RF_oobSize[treeID];
    membershipIndex = RF_oobMembershipIndex[treeID];
    break;
  }
  denomPtr = RF_vimpEnsembleDen[xVarIdx];
  for (uint i = 1; i <= membershipSize; i++) {
    uint ii = membershipIndex[i];
    Terminal *terminalNode = noiseMembership[ii];
    if ((terminalNode -> membrCount) > 0) {
      rfsrc_omp_atomic_update(&denomPtr[ii], 1.0);
#ifdef _OPENMP
      omp_set_lock(&(RF_lockVimp[xVarIdx][ii]));
#endif
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        for (uint k = 1; k <= RF_eventTypeSize; k++) {
          RF_vimpMRTstd[xVarIdx][k][ii] += terminalNode -> mortality[k];
        }
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          for (uint j = 1; j <= RF_rTargetFactorCount; j++) {
            for (uint k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              RF_vimpCLSstd[xVarIdx][j][k][ii] += (double) (terminalNode -> multiClassProb)[RF_rFactorMap[RF_rTargetFactor[j]]][k] / (double) (terminalNode -> membrCount);
            }
          }
          if (RF_optHigh & OPT_CSV) {          
            for (uint j = 1; j <= RF_rTargetFactorCount; j++) {
              RF_csvNumCLSptr[xVarIdx][j][ii] += 
                ( (terminalNode -> maxClass)[RF_rFactorMap[RF_rTargetFactor[j]]] ==
                  (uint) RF_response[treeID][RF_rFactorMap[RF_rTargetFactor[j]]][ii] ) ? 1 : 0;
            }
          }
        }
        if (RF_rTargetNonFactorCount > 0) {
          for (uint j = 1; j <= RF_rTargetNonFactorCount; j++) {
            RF_vimpRGRstd[xVarIdx][j][ii] += (terminalNode -> meanResponse)[RF_rNonFactorMap[RF_rTargetNonFactor[j]]];
          }
          if (RF_optHigh & OPT_CSV) {          
            for (uint j = 1; j <= RF_rTargetNonFactorCount; j++) {
              RF_csvNumRGRptr[xVarIdx][j][ii] += 
                pow((terminalNode -> meanResponse)[RF_rNonFactorMap[RF_rTargetNonFactor[j]]] -
                    RF_response[treeID][RF_rNonFactorMap[RF_rTargetNonFactor[j]]][ii], 2.0);
            }
          }
        }
        if ((RF_rTargetFactorCount > 0) || (RF_rTargetNonFactorCount > 0)) {
          if (RF_optHigh & OPT_CSV) {
            RF_csvDENptr[xVarIdx][ii] ++;
          }
        }
      }
#ifdef _OPENMP
      omp_unset_lock(&(RF_lockVimp[xVarIdx][ii]));
#endif
    }
    else {
      if (!(RF_opt & OPT_OUTC_TYPE)) {
        RF_nativePrint("\nRF-SRC:  *** ERROR *** ");
        RF_nativePrint("\nRF-SRC:  NA encountered for VIMP outcome in terminal node:  %10d", terminalNode -> nodeID);
        RF_nativePrint("\nRF-SRC:  Please Contact Technical Support.");
        RF_nativeExit();
      }
    }
  }  
  ${trace.token}    if (getTraceFlag(treeID) & OUTP_DEF_TRACE) {
  ${trace.token}      RF_nativePrint("\nGeneric variable importance outcome calculation:  treeID %10d \n", treeID);
  ${trace.token}      uint obsSize;
  ${trace.token}      if (mode == RF_PRED) {
  ${trace.token}        obsSize = RF_fobservationSize;
  ${trace.token}      }
  ${trace.token}      else {
  ${trace.token}        obsSize = RF_observationSize;
  ${trace.token}      }
  ${trace.token}      RF_nativePrint("\n  VIMP Ensemble Denom Count: \n");
  ${trace.token}      RF_nativePrint("          ");
  ${trace.token}      for (uint i=1; i <= obsSize; i++) {
  ${trace.token}        RF_nativePrint("%10d", i);
  ${trace.token}      }
  ${trace.token}      RF_nativePrint("\n");
  ${trace.token}      RF_nativePrint("          ");
  ${trace.token}      for (uint i=1; i <= obsSize; i++) {
  ${trace.token}        RF_nativePrint("%10.0f", denomPtr[i]);
  ${trace.token}      }
  ${trace.token}      RF_nativePrint("\n");
  ${trace.token}      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
  ${trace.token}        for (uint k=1; k <= RF_eventTypeSize; k++) {
  ${trace.token}          RF_nativePrint("\n  [xVarIdx][event] = [%10d][%10d] \n", xVarIdx, k);
  ${trace.token}          RF_nativePrint("          ");
  ${trace.token}          for (uint i=1; i <= obsSize; i++) {
  ${trace.token}            RF_nativePrint("%10d", i);
  ${trace.token}          }
  ${trace.token}          RF_nativePrint("\n");
  ${trace.token}          for (uint i=1; i <= obsSize; i++) {
  ${trace.token}            RF_nativePrint("%10.4f", RF_vimpMRTstd[xVarIdx][k][i]);
  ${trace.token}          }
  ${trace.token}          RF_nativePrint("\n");
  ${trace.token}        }
  ${trace.token}      }
  ${trace.token}      else {
  ${trace.token}        if (RF_rTargetFactorCount > 0) {
  ${trace.token}            for (uint j=1; j <= RF_rTargetFactorCount; j++) {
  ${trace.token}              for (uint k=1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
  ${trace.token}                RF_nativePrint("\n  [xVarIdx][rTarget][level] = [%10d][%10d][%10d] \n", xVarIdx, j, k);
  ${trace.token}                RF_nativePrint("          ");
  ${trace.token}                for (uint i=1; i <= obsSize; i++) {
  ${trace.token}                  RF_nativePrint("%10d", i);
  ${trace.token}                }
  ${trace.token}                RF_nativePrint("\n");
  ${trace.token}                RF_nativePrint("          ");
  ${trace.token}                for (uint i=1; i <= obsSize; i++) {
  ${trace.token}                  RF_nativePrint("%10.4f", RF_vimpCLSstd[xVarIdx][j][k][i]);
  ${trace.token}                }
  ${trace.token}                RF_nativePrint("\n");
  ${trace.token}              }
  ${trace.token}            }
  ${trace.token}          if (RF_optHigh & OPT_CSV) {
  ${trace.token}            RF_nativePrint("\nCLAS CSV Denominator calculation: \n");
  ${trace.token}            RF_nativePrint("          ");
  ${trace.token}            for (uint i=1; i <= obsSize; i++) {
  ${trace.token}              RF_nativePrint("%10d", i);
  ${trace.token}            }
  ${trace.token}            RF_nativePrint("\n");
  ${trace.token}            RF_nativePrint("          ");
  ${trace.token}            for (uint i=1; i <= obsSize; i++) {
  ${trace.token}              RF_nativePrint("%10d", RF_csvDENptr[xVarIdx][i]);
  ${trace.token}            }
  ${trace.token}            RF_nativePrint("\n");
  ${trace.token}            RF_nativePrint("\nCLAS CSV Numerator calculation: \n");
  ${trace.token}            for (uint j=1; j <= RF_rTargetFactorCount; j++) {
  ${trace.token}              RF_nativePrint("\n  [xVarIdx][rTarget] = [%10d][%10d] \n", xVarIdx, j);
  ${trace.token}              RF_nativePrint("          ");
  ${trace.token}              for (uint i=1; i <= obsSize; i++) {
  ${trace.token}                RF_nativePrint("%10d", i);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}              RF_nativePrint("          ");
  ${trace.token}              for (uint i=1; i <= obsSize; i++) {
  ${trace.token}                RF_nativePrint("%10.4f", RF_csvNumCLSptr[xVarIdx][j][i]);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}            }
  ${trace.token}          }
  ${trace.token}        }
  ${trace.token}        if (RF_rTargetNonFactorCount > 0) {
  ${trace.token}          for (uint j=1; j <= RF_rTargetNonFactorCount; j++) {
  ${trace.token}            RF_nativePrint("\n  [xVarIdx][rTarget] = [%10d][%10d] \n", xVarIdx, j);
  ${trace.token}            RF_nativePrint("          ");
  ${trace.token}            for (uint i=1; i <= obsSize; i++) {
  ${trace.token}              RF_nativePrint("%10d", i);
  ${trace.token}            }
  ${trace.token}            RF_nativePrint("\n");
  ${trace.token}            RF_nativePrint("          ");
  ${trace.token}            for (uint i=1; i <= obsSize; i++) {
  ${trace.token}              RF_nativePrint("%10.4f", RF_vimpRGRstd[xVarIdx][j][i]);
  ${trace.token}            }
  ${trace.token}            RF_nativePrint("\n");
  ${trace.token}          }
  ${trace.token}          if (RF_optHigh & OPT_CSV) {
  ${trace.token}            RF_nativePrint("\nREGR CSV Denominator calculation: \n");
  ${trace.token}            RF_nativePrint("          ");
  ${trace.token}            for (uint i=1; i <= obsSize; i++) {
  ${trace.token}              RF_nativePrint("%10d", i);
  ${trace.token}            }
  ${trace.token}            RF_nativePrint("\n");
  ${trace.token}            RF_nativePrint("          ");
  ${trace.token}            for (uint i=1; i <= obsSize; i++) {
  ${trace.token}              RF_nativePrint("%10d", RF_csvDENptr[xVarIdx][i]);
  ${trace.token}            }
  ${trace.token}            RF_nativePrint("\n");
  ${trace.token}            RF_nativePrint("\nREGR CSV Numerator calculation: \n");
  ${trace.token}            for (uint j=1; j <= RF_rTargetNonFactorCount; j++) {
  ${trace.token}              RF_nativePrint("\n  [xVarIdx][rTarget] = [%10d][%10d] \n", xVarIdx, j);
  ${trace.token}              RF_nativePrint("          ");
  ${trace.token}              for (uint i=1; i <= obsSize; i++) {
  ${trace.token}                RF_nativePrint("%10d", i);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}              RF_nativePrint("          ");
  ${trace.token}              for (uint i=1; i <= obsSize; i++) {
  ${trace.token}                RF_nativePrint("%10.4f", RF_csvNumRGRptr[xVarIdx][j][i]);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}            }
  ${trace.token}          }
  ${trace.token}        }
  ${trace.token}      }
  ${trace.token}      RF_nativePrint("\n");
  ${trace.token}    }
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdateEnsembleVimp() EXIT ...\n");
  ${trace.token}  }
}
void summarizePerturbedPerformance(char       mode,
                                   uint       treeID,
                                   uint       bb,
                                   uint       p,
                                   double   **responsePtr) {
  uint      obsSize;
  double   *vimpDenom;
  double   **ensembleMRT;
  double  ***ensembleCLS;
  double   **ensembleRGR;
  double    *vimpMRTptr;
  double   **vimpCLSptr;
  double    *vimpRGRptr;
  uint i, j, k;
  ${trace.token}  uint n;
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nsummarizePerturbedPerformance(%10d, %10d) ENTRY ...\n", bb, p);
  ${trace.token}  }
  obsSize = (mode == RF_PRED) ?  RF_fobservationSize : RF_observationSize;  
  vimpDenom = RF_vimpEnsembleDen[p];
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    ${trace.token}  if (getTraceFlag(0) & OUTP_DEF_TRACE) {
    ${trace.token}    RF_nativePrint("\nUpdating Ensemble Mortality for VIMP (covariate):  %10d:  \n", p);
    ${trace.token}  }
    if (!(RF_opt & OPT_COMP_RISK)) {
      getEnsembleMortality(mode, treeID, obsSize, RF_vimpMRTstd[p], vimpDenom, RF_vimpMRTstd[p][1]);
    }
    else {
      getEnsembleMortalityCR(mode, treeID, obsSize, RF_vimpMRTstd[p], vimpDenom, RF_vimpMRTstd[p]);
    }
  }  
  else {
    if (RF_rTargetFactorCount > 0) {
      for (i = 1; i <= obsSize; i++) {
        if(vimpDenom[i] > 0) {
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              RF_vimpCLSstd[p][j][k][i] = RF_vimpCLSstd[p][j][k][i] / vimpDenom[i];
            }
          }
        }
      }
    }
    if (RF_rTargetNonFactorCount > 0) {
      for (i = 1; i <= obsSize; i++) {
        if(vimpDenom[i] > 0) {
          for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
            RF_vimpRGRstd[p][j][i] = RF_vimpRGRstd[p][j][i] / vimpDenom[i];
          }
        }
      }
    }
  }
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}      RF_nativePrint("\nPerformance Calculation for VIMP (covariate):  %10d", p);
  ${trace.token}    }
  ensembleMRT = NULL;
  ensembleCLS = NULL;
  ensembleRGR = NULL;
  vimpMRTptr = NULL;
  vimpCLSptr = NULL;
  vimpRGRptr = NULL;
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    ensembleMRT =  RF_vimpMRTstd[p];
  }
  else {
    if (RF_rTargetFactorCount > 0) {
      ensembleCLS = RF_vimpCLSstd[p];
    }
    if (RF_rTargetNonFactorCount > 0) {
      ensembleRGR = RF_vimpRGRstd[p];
    }
  }
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    vimpMRTptr = RF_vimpMRTblk[bb][p];
  }
  else {
    if (RF_rTargetFactorCount > 0) {
      vimpCLSptr = RF_vimpCLSblk[bb][p];
    }
    if (RF_rTargetNonFactorCount > 0) {
      vimpRGRptr = RF_vimpRGRblk[bb][p];
    }    
  }
  ${trace.token}  if (getTraceFlag(0) & OUTP_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nEnsemble input variable importance calculation:  \n");
  ${trace.token}    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
  ${trace.token}      RF_nativePrint("\nMRT variable importance calculation for [covariate] : [p] =   [%10d] \n", p);
  ${trace.token}      for (k = 1; k <= RF_eventTypeSize; k++) {
  ${trace.token}        RF_nativePrint("\n  [event] = [%10d] \n", k);
  ${trace.token}        RF_nativePrint("          ");
  ${trace.token}        for (n=1; n <= obsSize; n++) {
  ${trace.token}          RF_nativePrint("%10d", n);
  ${trace.token}        }
  ${trace.token}        RF_nativePrint("\n");
  ${trace.token}        for (n=1; n <= obsSize; n++) {
  ${trace.token}          RF_nativePrint("%10.4f", ensembleMRT[k][n]);
  ${trace.token}        }
  ${trace.token}        RF_nativePrint("\n");
  ${trace.token}      }
  ${trace.token}    }
  ${trace.token}    else {
  ${trace.token}      if (RF_rTargetFactorCount > 0) {
  ${trace.token}          for (j=1; j <= RF_rTargetFactorCount; j++) {
  ${trace.token}          RF_nativePrint("\nCLS variable importance calculation for [covariate][target] : [p][j] =   [%10d][%10d] \n", p, j);
  ${trace.token}            for (k=1; k <=  RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
  ${trace.token}              RF_nativePrint("\n  [level] = [%10d] \n", k);
  ${trace.token}              RF_nativePrint("          ");
  ${trace.token}              for (n=1; n <= obsSize; n++) {
  ${trace.token}                RF_nativePrint("%10d", n);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}              RF_nativePrint("          ");
  ${trace.token}              for (n=1; n <= obsSize; n++) {
  ${trace.token}                RF_nativePrint("%10.4f", ensembleCLS[j][k][n]);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}            }
  ${trace.token}          }
  ${trace.token}      }
  ${trace.token}      if (RF_rTargetNonFactorCount > 0) {
  ${trace.token}        for (j=1; j <= RF_rTargetNonFactorCount; j++) {
  ${trace.token}          RF_nativePrint("\nRGR variable importance calculation for [covariate][target] : [p][j] =   [%10d][%10d] \n", p, j);
  ${trace.token}          RF_nativePrint("          ");
  ${trace.token}          for (n=1; n <= obsSize; n++) {
  ${trace.token}            RF_nativePrint("%10d", n);
  ${trace.token}          }
  ${trace.token}          RF_nativePrint("\n");
  ${trace.token}          RF_nativePrint("          ");
  ${trace.token}          for (n=1; n <= obsSize; n++) {
  ${trace.token}            RF_nativePrint("%10.4f", ensembleRGR[j][n]);
  ${trace.token}          }
  ${trace.token}          RF_nativePrint("\n");
  ${trace.token}        }
  ${trace.token}      }
  ${trace.token}    }
  ${trace.token}  }
  getPerformance(treeID,            
                 mode,
                 obsSize,
                 responsePtr,
                 vimpDenom,
                 ensembleMRT,       
                 ensembleCLS,       
                 ensembleRGR,       
                 vimpMRTptr,        
                 vimpCLSptr,        
                 vimpRGRptr);       
  ${trace.token}  if (getTraceFlag(0) & OUTP_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nEnsemble output variable importance calculation:  \n");
  ${trace.token}    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
  ${trace.token}      RF_nativePrint("\nMRT perturbed performance for [covariate] : [p] = [%10d] \n", p);
  ${trace.token}      RF_nativePrint("Event:    ");
  ${trace.token}      for (k = 1; k <= RF_eventTypeSize; k++) {
  ${trace.token}        RF_nativePrint("%10d", k);
  ${trace.token}      }
  ${trace.token}      RF_nativePrint("\n");
  ${trace.token}      RF_nativePrint("          ");
  ${trace.token}      for (k = 1; k <= RF_eventTypeSize; k++) {
  ${trace.token}        RF_nativePrint(" %10.4f ", vimpMRTptr[k]);
  ${trace.token}      }
  ${trace.token}    }
  ${trace.token}    else {
  ${trace.token}      if (RF_rTargetFactorCount > 0) {
  ${trace.token}        for (j=1; j <= RF_rTargetFactorCount; j++) {
  ${trace.token}          RF_nativePrint("\nCLS perturbed performance for [covariate][target] : [p][j] =  [%10d][%10d] \n", p, j);
  ${trace.token}          RF_nativePrint("Level:    ");
  ${trace.token}          for (k=1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
  ${trace.token}            RF_nativePrint("%10d", k);
  ${trace.token}          }
  ${trace.token}          RF_nativePrint("\n");
  ${trace.token}          RF_nativePrint("          ");
  ${trace.token}          for (k=1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
  ${trace.token}            RF_nativePrint("%10.4f", vimpCLSptr[j][k]);
  ${trace.token}          }
  ${trace.token}          RF_nativePrint("\n");
  ${trace.token}        }
  ${trace.token}      }
  ${trace.token}      if (RF_rTargetNonFactorCount > 0) {
  ${trace.token}          RF_nativePrint("\nRGR perturbed performance for [covariate] : [p] = [%10d] \n", p);
  ${trace.token}        RF_nativePrint("Resp:     ");
  ${trace.token}        for (j=1; j <= RF_rTargetNonFactorCount; j++) {
  ${trace.token}          RF_nativePrint("%10d", j);
  ${trace.token}        }
  ${trace.token}        RF_nativePrint("\n");
  ${trace.token}        RF_nativePrint("          ");
  ${trace.token}        for (j=1; j <= RF_rTargetNonFactorCount; j++) {
  ${trace.token}          RF_nativePrint("%10.4f", vimpRGRptr[j]);
  ${trace.token}        }
  ${trace.token}        RF_nativePrint("\n");
  ${trace.token}      }
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nsummarizePerturbedPerformance(%10d, %10d) EXIT ...\n", bb, p);
  ${trace.token}  }
}
void finalizeVimpPerformance(char       mode) {
  double  ***vimpMRTptr;
  double ****vimpCLSptr;
  double  ***vimpRGRptr;
  double   **perfMRTptr;
  double  ***perfCLSptr;
  double   **perfRGRptr;
  uint xVimpCount;
  double result;
  uint cumDenomCount;
  uint genericBlockCount;
  uint i, j, k, p;
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nfinalizeVimpPerformance() ENTRY ...\n");
  ${trace.token}  }
  if (RF_opt & OPT_VIMP_JOIN) {
    xVimpCount = 1;
  }
  else {
    xVimpCount = RF_intrPredictorSize;
  }
    vimpMRTptr = RF_vimpMRTblk;
    vimpCLSptr = RF_vimpCLSblk;
    vimpRGRptr = RF_vimpRGRblk;
    perfMRTptr = RF_perfMRTblk;
    perfCLSptr = RF_perfCLSblk;
    perfRGRptr = RF_perfRGRblk;
    genericBlockCount = RF_perfBlockCount;
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nVIMP BEFORE Finalization:  \n");
  ${trace.token}    for (p=1; p <= xVimpCount; p++) {
  ${trace.token}      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
  ${trace.token}        RF_nativePrint("\nMRT values for covariate [p]:  %10d \n", p);
  ${trace.token}        RF_nativePrint("Event:              ");
  ${trace.token}        for (k = 1; k <= RF_eventTypeSize; k++) {
  ${trace.token}          RF_nativePrint(" %10d", k);
  ${trace.token}        }
  ${trace.token}        RF_nativePrint("\n");
  ${trace.token}        for (i=1; i <= RF_perfBlockCount; i++) {
  ${trace.token}          RF_nativePrint("Perturbed:%10d", i);
  ${trace.token}          for (k = 1; k <= RF_eventTypeSize; k++) {
  ${trace.token}            RF_nativePrint(" %10.4f ", RF_vimpMRTblk[i][p][k]);
  ${trace.token}          }
  ${trace.token}        }
  ${trace.token}      }
  ${trace.token}      else {
  ${trace.token}        if (RF_rTargetFactorCount > 0) {
  ${trace.token}          for (j=1; j <= RF_rTargetFactorCount; j++) {
  ${trace.token}            RF_nativePrint("\nCLS values for [covariate][target] : [p][j] =  [%10d][%10d] \n", p, j);
  ${trace.token}            RF_nativePrint("Level:              ");
  ${trace.token}            for (k=1; k <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
  ${trace.token}              RF_nativePrint(" %10d", k);
  ${trace.token}            }
  ${trace.token}            RF_nativePrint("\n");
  ${trace.token}            for (i=1; i <= RF_perfBlockCount; i++) {
  ${trace.token}              RF_nativePrint("Perturbed:%10d", i);
  ${trace.token}              for (k=1; k <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
  ${trace.token}                RF_nativePrint(" %10.4f", RF_vimpCLSblk[i][p][j][k]);
  ${trace.token}              }
  ${trace.token}              RF_nativePrint("\n");
  ${trace.token}            }
  ${trace.token}          }
  ${trace.token}        }
  ${trace.token}        if (RF_rTargetNonFactorCount > 0) {
  ${trace.token}          RF_nativePrint("\nRGR values for [covariate][target] : [p] =  [%10d] \n", p);
  ${trace.token}          RF_nativePrint("Resp:               ");
  ${trace.token}          for (j=1; j <= RF_rTargetNonFactorCount; j++) {
  ${trace.token}            RF_nativePrint(" %10d", j);
  ${trace.token}          }
  ${trace.token}          for (i=1; i <= RF_perfBlockCount; i++) {
  ${trace.token}          RF_nativePrint("\n");
  ${trace.token}          RF_nativePrint("Faithful: %10d", i);
  ${trace.token}          for (j=1; j <= RF_rTargetNonFactorCount; j++) {
  ${trace.token}            RF_nativePrint(" %10.4f", RF_perfRGRblk[i][j]);
  ${trace.token}          }
  ${trace.token}          RF_nativePrint("\n");
  ${trace.token}          }
  ${trace.token}          for (i=1; i <= RF_perfBlockCount; i++) {  
  ${trace.token}          RF_nativePrint("\n");
  ${trace.token}          RF_nativePrint("Perturbed:%10d", i);
  ${trace.token}          for (j=1; j <= RF_rTargetNonFactorCount; j++) {
  ${trace.token}            RF_nativePrint(" %10.4f", RF_vimpRGRblk[i][p][j]);
  ${trace.token}          }
  ${trace.token}          RF_nativePrint("\n");
  ${trace.token}          }
  ${trace.token}        }
  ${trace.token}      }
  ${trace.token}    }
  ${trace.token}  }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
      ${trace.token}    RF_nativePrint("\nVIMP DURING Finalization:  \n");
      ${trace.token}  }
      for(p = 1; p <= xVimpCount; p++) {
        for (k = 1; k <= RF_eventTypeSize; k++) {
          cumDenomCount = 0;
          result = 0.0;
          for (i = 1; i <= genericBlockCount; i++) {
            if(!RF_nativeIsNaN(vimpMRTptr[i][p][k])) {
              if(!RF_nativeIsNaN(perfMRTptr[i][k])) {
                ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
                ${trace.token}    RF_nativePrint("\nMRT variable importance value for [block][covariate][event] : [i][p][k] = [%10d][%10d][%10d] = %10.4f versus %10.4f \n", i, p, k, vimpMRTptr[i][p][k], perfMRTptr[i][k]);
                ${trace.token}  }
                result += vimpMRTptr[i][p][k] - perfMRTptr[i][k];
                cumDenomCount ++;
              }
            }
          }
          if (cumDenomCount != 0) {
            RF_vimpMRTptr[p][k] = result / (double) cumDenomCount;
          }
          else {
            RF_vimpMRTptr[p][k] = RF_nativeNaN;
          }
          ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
          ${trace.token}    RF_nativePrint("\nMRT variable importance denominator count for [covariate] : [p] = [%10d] = %10d \n", p, cumDenomCount);
          ${trace.token}  }
        }
      }
    }
    else {
      if (RF_rTargetFactorCount > 0) {
        ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
        ${trace.token}    RF_nativePrint("\nVIMP DURING Finalization:  \n");
        ${trace.token}  }
        for(p = 1; p <= xVimpCount; p++) {
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            for (k = 1; k <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              cumDenomCount = 0;
              result = 0.0;
              for (i = 1; i <= genericBlockCount; i++) {
                if(!RF_nativeIsNaN(vimpCLSptr[i][p][j][k])) {
                  if(!RF_nativeIsNaN(perfCLSptr[i][j][k])) {
                    ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
                    ${trace.token}    RF_nativePrint("\nCLS variable importance value for [block][covariate][target][level] : [i][p][j][k] = [%10d][%10d][%10d][%10d] = %10.4f versus %10.4f \n", i, p, j, k, vimpCLSptr[i][p][j][k], perfCLSptr[i][j][k]);
                    ${trace.token}  }
                    result += vimpCLSptr[i][p][j][k] - perfCLSptr[i][j][k];
                    cumDenomCount ++;
                  }
                }
              }
              if (cumDenomCount != 0) {
                RF_vimpCLSptr[p][j][k] = result / (double) cumDenomCount;
              }
              else {
                RF_vimpCLSptr[p][j][k] = RF_nativeNaN;
              }
              ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
              ${trace.token}    RF_nativePrint("\nCLS variable importance denominator count for [covariate][target][level] : [p][j][k] = [%10d][%10d][%10d] = %10d \n", p, j, k, cumDenomCount);
              ${trace.token}  }
            }
          }
        }
      }
      if (RF_rTargetNonFactorCount > 0) {
        for(p = 1; p <= xVimpCount; p++) {
          for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
            cumDenomCount = 0;
            result = 0.0;
            for (i = 1; i <= genericBlockCount; i++) {
              if(!RF_nativeIsNaN(vimpRGRptr[i][p][j])) {
                if(!RF_nativeIsNaN(perfRGRptr[i][j])) {
                  result += vimpRGRptr[i][p][j] - perfRGRptr[i][j];
                  cumDenomCount ++;
                }
              }
            }
            if (cumDenomCount != 0) {
              RF_vimpRGRptr[p][j] = result / (double) cumDenomCount;
            }
            else {
              RF_vimpRGRptr[p][j] = RF_nativeNaN;
            }
            ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
            ${trace.token}    RF_nativePrint("\nRGR variable importance denominator count for [covariate][target] : [p][j] = [%10d][%10d] = %10d \n", p, j, cumDenomCount);
            ${trace.token}  }
          }
        }
      }    
    }
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nVIMP AFTER Finalization:  \n");
  ${trace.token}    for (p=1; p <= xVimpCount; p++) {
  ${trace.token}      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
  ${trace.token}        RF_nativePrint("\nMRT variable importance value for [covariate] : [p] = [%10d] \n", p);
  ${trace.token}        RF_nativePrint("Event:    ");
  ${trace.token}        for (k = 1; k <= RF_eventTypeSize; k++) {
  ${trace.token}          RF_nativePrint(" %10d", k);
  ${trace.token}        }
  ${trace.token}        RF_nativePrint("\n");
  ${trace.token}        RF_nativePrint("          ");
  ${trace.token}        for (k = 1; k <= RF_eventTypeSize; k++) {
  ${trace.token}          RF_nativePrint(" %10.4f ", RF_vimpMRTptr[p][k]);
  ${trace.token}        }
  ${trace.token}      }
  ${trace.token}      else {
  ${trace.token}        if (RF_rTargetFactorCount > 0) {
  ${trace.token}          for (j=1; j <= RF_rTargetFactorCount; j++) {
  ${trace.token}            RF_nativePrint("\nCLS variable importance value for [covariate][target] : [p][j] = [%10d][%10d] \n", p, j);
  ${trace.token}            RF_nativePrint("Level:    ");
  ${trace.token}            for (k=1; k <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
  ${trace.token}              RF_nativePrint(" %10d", k);
  ${trace.token}            }
  ${trace.token}            RF_nativePrint("\n");
  ${trace.token}            RF_nativePrint("          ");
  ${trace.token}            for (k=1; k <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
  ${trace.token}              RF_nativePrint(" %10.4f", RF_vimpCLSptr[p][j][k]);
  ${trace.token}            }
  ${trace.token}            RF_nativePrint("\n");
  ${trace.token}          }
  ${trace.token}        }
  ${trace.token}        if (RF_rTargetNonFactorCount > 0) {
  ${trace.token}          RF_nativePrint("\nRGR variable importance value for [covariate] : [p] = [%10d] \n", p);
  ${trace.token}          RF_nativePrint("Resp:     ");
  ${trace.token}          for (j=1; j <= RF_rTargetNonFactorCount; j++) {
  ${trace.token}            RF_nativePrint(" %10d", j);
  ${trace.token}          }
  ${trace.token}          RF_nativePrint("\n");
  ${trace.token}          RF_nativePrint("          ");
  ${trace.token}          for (j=1; j <= RF_rTargetNonFactorCount; j++) {
  ${trace.token}            RF_nativePrint(" %10.4f", RF_vimpRGRptr[p][j]);
  ${trace.token}          }
  ${trace.token}          RF_nativePrint("\n");
  ${trace.token}        }
  ${trace.token}      }
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nfinalizeVimpPerformance() EXIT ...\n");
  ${trace.token}  }
}
void stackVimpMembership(char mode, Terminal ***membership) {
  uint obsSize;
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {  
  ${trace.token}      RF_nativePrint("\nstackVimpMembership() ENTRY ...\n");
  ${trace.token}    }
  ${trace.token}  }
  (*membership) = NULL;
  if (RF_opt & OPT_VIMP) {
    switch (mode) {
    case RF_PRED:
      obsSize = RF_fobservationSize;
      break;
    default:
      obsSize = RF_observationSize;
      break;
    }
    *membership = (Terminal **) new_vvector(1, obsSize, NRUTIL_TPTR);
  }
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {  
  ${trace.token}      RF_nativePrint("\nstackVimpMembership() EXIT ...\n");
  ${trace.token}    }
  ${trace.token}  }
}
void unstackVimpMembership(char mode, Terminal **membership) {
  uint obsSize;
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {  
  ${trace.token}      RF_nativePrint("\nunstackVimpMembership() ENTRY ...\n");
  ${trace.token}    }
  ${trace.token}  }
  if (RF_opt & OPT_VIMP) {
    switch (mode) {
    case RF_PRED:
      obsSize = RF_fobservationSize;
      break;
    default:
      obsSize = RF_observationSize;
      break;
    }
    free_new_vvector(membership, 1, obsSize, NRUTIL_TPTR);
  }
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {  
  ${trace.token}      RF_nativePrint("\nunstackVimpMembership() EXIT ...\n");
  ${trace.token}    }
  ${trace.token}  }
}
void normalizeBlockedEnsembleEstimates(char      mode,
                                       double  **ensembleMRTptr,
                                       double ***ensembleCLSptr,
                                       double  **ensembleRGRptr,
                                       double   *ensembleDen) {
  uint      obsSize;
  uint i, j, k;
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nnormalizeBlockedEnsembleEstimates() ENTRY ...\n");
  ${trace.token}  }
  obsSize = (mode == RF_PRED) ? RF_fobservationSize : RF_observationSize;
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      for (i = 1; i <= obsSize; i++) {
        if (ensembleDen[i] != 0) {
          if (!(RF_opt & OPT_COMP_RISK)) {
            ensembleMRTptr[1][i] = ensembleMRTptr[1][i] / ensembleDen[i];
          }
          else {
            for(j = 1; j <= RF_eventTypeSize; j ++) {
              ensembleMRTptr[j][i] = ensembleMRTptr[j][i] / ensembleDen[i];
            }
          }
        }
        else {
        }
      }
    }  
    else {
      if (RF_rTargetFactorCount > 0) {
        if (ensembleCLSptr != NULL) {
        for (i = 1; i <= obsSize; i++) {
          if (ensembleDen[i] != 0) {
            for (j = 1; j <= RF_rTargetFactorCount; j++) {
              for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
                ensembleCLSptr[j][k][i] = ensembleCLSptr[j][k][i] / ensembleDen[i];
              }
            }
          }
          else {
          }
        }
        }
      }
      if (RF_rTargetNonFactorCount > 0) {      
        if (ensembleRGRptr != NULL) {
        for (i = 1; i <= obsSize; i++) {
          if (ensembleDen[i] != 0) {
            for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
              ensembleRGRptr[j][i] = ensembleRGRptr[j][i] / ensembleDen[i];
            }
          }
          else {
          }
        }
        }
      }
    }
    ${trace.token}    if (getTraceFlag(0) & ENSB_LOW_TRACE) {
    ${trace.token}      RF_nativePrint("\nNormalized Blocked Ensembles Follow:  \n");
    ${trace.token}      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    ${trace.token}        if (!(RF_opt & OPT_COMP_RISK)) {
    ${trace.token}            RF_nativePrint("\nEnsemble MRT:  \n");
    ${trace.token}            RF_nativePrint("          ");
    ${trace.token}            for (i=1; i <= obsSize; i++) {
    ${trace.token}              RF_nativePrint("%10d", i);
    ${trace.token}            }
    ${trace.token}            RF_nativePrint("\n");
    ${trace.token}            RF_nativePrint("          ");
    ${trace.token}            for (i=1; i <= obsSize; i++) {
    ${trace.token}              RF_nativePrint("%10.4f", ensembleMRTptr[1][i]);
    ${trace.token}            }
    ${trace.token}            RF_nativePrint("\n");
    ${trace.token}        }
    ${trace.token}        else {
    ${trace.token}          for (j = 1; j <= RF_eventTypeSize; j++) {
    ${trace.token}            RF_nativePrint("\nEnsemble MRT:  Event %10d \n", j);
    ${trace.token}            RF_nativePrint("          ");
    ${trace.token}            for (i=1; i <= obsSize; i++) {
    ${trace.token}              RF_nativePrint("%10d", i);
    ${trace.token}            }
    ${trace.token}            RF_nativePrint("\n");
    ${trace.token}            RF_nativePrint("          ");
    ${trace.token}            for (i=1; i <= obsSize; i++) {
    ${trace.token}              RF_nativePrint("%10.4f", ensembleMRTptr[j][i]);
    ${trace.token}            }
    ${trace.token}            RF_nativePrint("\n");
    ${trace.token}          }
    ${trace.token}        }
    ${trace.token}      }
    ${trace.token}      else {
    ${trace.token}        if (RF_rTargetFactorCount > 0) {
    ${trace.token}          RF_nativePrint("\nEnsemble:  \n");
    ${trace.token}          RF_nativePrint("          ");
    ${trace.token}          for (i=1; i <= obsSize; i++) {
    ${trace.token}            RF_nativePrint("%10d", i);
    ${trace.token}          }
    ${trace.token}          for (j = 1; j <= RF_rTargetFactorCount; j++) {
    ${trace.token}            RF_nativePrint("\nTarget:   %10d \n", RF_rTargetFactor[j]);
    ${trace.token}            for (k=1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
    ${trace.token}              RF_nativePrint("%10d", k);
    ${trace.token}              for (i=1; i <= obsSize; i++) {
    ${trace.token}                RF_nativePrint("%10.4f", ensembleCLSptr[j][k][i]);
    ${trace.token}              }
    ${trace.token}              RF_nativePrint("\n");
    ${trace.token}            }
    ${trace.token}          }
    ${trace.token}        }
    ${trace.token}        if (RF_rTargetNonFactorCount > 0) {
    ${trace.token}          RF_nativePrint("\nEnsemble:  \n");
    ${trace.token}          RF_nativePrint("          ");
    ${trace.token}          for (i=1; i <= obsSize; i++) {
    ${trace.token}            RF_nativePrint("%10d", i);
    ${trace.token}          }
    ${trace.token}          for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
    ${trace.token}            RF_nativePrint("\nTarget:   %10d \n", RF_rTargetNonFactor[j]);
    ${trace.token}            RF_nativePrint("          ");
    ${trace.token}            for (i=1; i <= obsSize; i++) {
    ${trace.token}              RF_nativePrint("%10.4f", ensembleRGRptr[j][i]);
    ${trace.token}            }
    ${trace.token}            RF_nativePrint("\n");
    ${trace.token}          }
    ${trace.token}        }
    ${trace.token}      }
    ${trace.token}      RF_nativePrint("\nRF-SRC:  Ensemble outputs finalized. \n");
    ${trace.token}    }
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nnormalizeBlockedEnsembleEstimates() EXIT ...\n");
  ${trace.token}  }
}
void resetBlockedEnsembleEstimates(char mode) {
  uint      obsSize, xVimpCount;
  uint i, j, k, p;
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nresetBlockedEnsembleEstimates() ENTRY ...\n");
  ${trace.token}  }
  obsSize = (mode == RF_PRED) ? RF_fobservationSize : RF_observationSize;
  if (RF_opt & OPT_VIMP_JOIN) {
    xVimpCount = 1;
  }
  else {
    xVimpCount = RF_intrPredictorSize;
  }
  for (i = 1; i <= obsSize; i++) {
    RF_blkEnsembleDen[i]  = 0.0;
    for (p = 1; p <= xVimpCount; p++) {
      RF_vimpEnsembleDen[p][i] = 0.0;
    }
  }
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    for(j = 1; j <= RF_eventTypeSize; j ++) {
      for (i = 1; i <= obsSize; i++) {
        RF_blkEnsembleMRTnum[j][i] = 0.0;
        for (p = 1; p <= xVimpCount; p++) {
          RF_vimpMRTstd[p][j][i] = 0.0;
        }
      }
    }
  }  
  else {
    if (RF_rTargetFactorCount > 0) {
      for (j = 1; j <= RF_rTargetFactorCount; j++) {
        for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
          for (i = 1; i <= obsSize; i++) {
            RF_blkEnsembleCLSnum[j][k][i] = 0.0;
            for (p = 1; p <= xVimpCount; p++) {
              RF_vimpCLSstd[p][j][k][i] = 0.0;
            }
          }
        }
      }
    }
    if (RF_rTargetNonFactorCount > 0) {      
      for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
        for (i = 1; i <= obsSize; i++) {
          RF_blkEnsembleRGRnum[j][i] = 0.0;
          for (p = 1; p <= xVimpCount; p++) {
            RF_vimpRGRstd[p][j][i] = 0.0;
          }
        }
      }
    }
  }
  ${trace.token}  if (getTraceFlag(0) & VIMP_LOW_TRACE) {
  ${trace.token}    RF_nativePrint("\nresetBlockedEnsembleEstimates() EXIT ...\n");
  ${trace.token}  }
}
#ifdef _OPENMP
void rfsrc_omp_atomic_update(double *addr, double incr) {
#pragma omp atomic update
  (*addr) += incr;
}
#else
void rfsrc_omp_atomic_update(double *addr, double incr) { 
  (*addr) += incr;
}
#endif
uint getVimpRecoverySeedDimension(char mode, uint opt) {
  uint bnpSize;
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetVimpRecoverySeedDimension() ENTRY ...\n");
  ${trace.token}  }
  bnpSize = 0;
  if (opt & OPT_VIMP) {
      bnpSize = RF_ntree;
  }
  ${trace.token}    if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}      RF_nativePrint("\n Vimp Recovery Seed Dimension:  (mode, opt, dim) = (%10d, %10d, %10d)", mode, opt, bnpSize);
  ${trace.token}    }
  ${trace.token}  if (getTraceFlag(0) & SUMM_MED_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetVimpRecoverySeedDimension() EXIT ...\n");
  ${trace.token}  }
  return bnpSize;
}
