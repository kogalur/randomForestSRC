
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
    if (daughterFlag == LEFT) {
      result = identifyExtrapolatedMembership(parent ->  left, yShadow, xShadow);
    }
    else {
      result = identifyExtrapolatedMembership(parent -> right, yShadow, xShadow);
    }
  }
  return result;
}
void getVimpMembership (char       mode,
                        uint       treeID,
                        Terminal **vimpMembership,
                        uint p) {
  char result;
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
  }
}
void updateEnsembleVimp (char       mode,
                         uint       treeID,
                         Terminal **noiseMembership,
                         uint       xVarIdx) {
  uint  *membershipIndex;
  uint   membershipSize;
  double *denomPtr;
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
  obsSize = (mode == RF_PRED) ?  RF_fobservationSize : RF_observationSize;  
  vimpDenom = RF_vimpEnsembleDen[p];
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
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
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      for(p = 1; p <= xVimpCount; p++) {
        for (k = 1; k <= RF_eventTypeSize; k++) {
          cumDenomCount = 0;
          result = 0.0;
          for (i = 1; i <= genericBlockCount; i++) {
            if(!RF_nativeIsNaN(vimpMRTptr[i][p][k])) {
              if(!RF_nativeIsNaN(perfMRTptr[i][k])) {
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
        }
      }
    }
    else {
      if (RF_rTargetFactorCount > 0) {
        for(p = 1; p <= xVimpCount; p++) {
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            for (k = 1; k <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              cumDenomCount = 0;
              result = 0.0;
              for (i = 1; i <= genericBlockCount; i++) {
                if(!RF_nativeIsNaN(vimpCLSptr[i][p][j][k])) {
                  if(!RF_nativeIsNaN(perfCLSptr[i][j][k])) {
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
          }
        }
      }    
    }
}
void stackVimpMembership(char mode, Terminal ***membership) {
  uint obsSize;
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
}
void unstackVimpMembership(char mode, Terminal **membership) {
  uint obsSize;
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
}
void normalizeBlockedEnsembleEstimates(char      mode,
                                       double  **ensembleMRTptr,
                                       double ***ensembleCLSptr,
                                       double  **ensembleRGRptr,
                                       double   *ensembleDen) {
  uint      obsSize;
  uint i, j, k;
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
}
void resetBlockedEnsembleEstimates(char mode) {
  uint      obsSize, xVimpCount;
  uint i, j, k, p;
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
  bnpSize = 0;
  if (opt & OPT_VIMP) {
      bnpSize = RF_ntree;
  }
  return bnpSize;
}
