
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "stackParallel.h"
#include "nrutil.h"
#ifdef _OPENMP
void stackLocksOpenMP(char mode) {
  uint i, j;
  omp_init_lock(&RF_lockEnsbUpdtCount);
  omp_init_lock(&RF_lockPerf);
  if (RF_optHigh & OPT_PART_PLOT) {
    RF_lockPartial = ompvector(1, RF_observationSize);
    for (i = 1; i <= RF_observationSize; i++) {
      omp_init_lock(&(RF_lockPartial[i]));
    }
  }
  if (RF_optHigh & OPT_WGHT) {
    uint gMembershipSize;
    if((RF_optHigh & OPT_WGHT_IBG) && (RF_optHigh & OPT_WGHT_OOB)) {
      switch (mode) {
      case RF_PRED:
        gMembershipSize = RF_fobservationSize;
        break;
      default:
        gMembershipSize = RF_observationSize;
        break;
      }
    }
    else {
      gMembershipSize  = RF_observationSize;
    }
    RF_lockWeight = (omp_lock_t **) new_vvector(1, gMembershipSize, NRUTIL_OMPLPTR);
    for (i = 1; i <= gMembershipSize; i++) {
      RF_lockWeight[i] = ompvector(1, RF_observationSize);
    }
    for (i = 1; i <= gMembershipSize; i++) {
      for (j = 1; j <= RF_observationSize; j++) {
        omp_init_lock(&(RF_lockWeight[i][j]));
      }
    }
    RF_lockWeightRow = ompvector(1, gMembershipSize);
    for (i = 1; i <= gMembershipSize; i++) {
      omp_init_lock(&(RF_lockWeightRow[i]));
    }
  }
  if (RF_opt & OPT_VIMP) {
      uint obsSize, xVimpSize;
      if (RF_opt & OPT_VIMP_JOIN) {
        xVimpSize = 1;
      }
      else {
        xVimpSize = RF_intrPredictorSize;
      }
      switch (mode) {
      case RF_PRED:
        obsSize = RF_fobservationSize;
        break;
      default:
        obsSize  = RF_observationSize;
        break;
      }
      RF_lockVimp = (omp_lock_t **) new_vvector(1, xVimpSize, NRUTIL_OMPLPTR);
      for (i = 1; i <= xVimpSize; i++) {
        RF_lockVimp[i] = ompvector(1, obsSize);
      }
      for (i = 1; i <= xVimpSize; i++) {
        for (j = 1; j <= obsSize; j++) {
          omp_init_lock(&(RF_lockVimp[i][j]));
        }
      }
      RF_lockVimpRow = ompvector(1, obsSize);
      for (i = 1; i <= obsSize; i++) {
        omp_init_lock(&(RF_lockVimpRow[i]));
      }
      RF_lockVimpCol = ompvector(1, xVimpSize);
      for (i = 1; i <= xVimpSize; i++) {
        omp_init_lock(&(RF_lockVimpCol[i]));
      }
  }
  if ((RF_vtry > 0) && (RF_vtryMode != RF_VTRY_NULL)) {
    uint xVimpSize;
    xVimpSize = RF_xSize;
    RF_lockVimpHoldout = (omp_lock_t **) new_vvector(1, xVimpSize, NRUTIL_OMPLPTR);
    for (i = 1; i <= xVimpSize; i++) {
      if (RF_holdBLKptr[i] > 0) {
        RF_lockVimpHoldout[i] = ompvector(1, RF_holdBLKptr[i]);
        for (j = 1; j <= RF_holdBLKptr[i]; j++) {
          omp_init_lock(&(RF_lockVimpHoldout[i][j]));
        }
      }
    }
  }
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    omp_lock_t   **lockDENptr;
    uint obsSize;  
    char oobFlag, fullFlag;
    oobFlag = fullFlag = FALSE;
    switch (mode) {
    case RF_PRED:
      if (RF_opt & OPT_FENS) {
        fullFlag = TRUE;
      }
      break;
    default:
      if (RF_opt & OPT_OENS) {
        oobFlag = TRUE;
      }
      if (RF_opt & OPT_FENS) {
        fullFlag = TRUE;
      }
      break;
    }
    while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
      if (oobFlag == TRUE) {
        lockDENptr = &RF_lockDENoens;
        obsSize = RF_observationSize;
      }
      else {
        lockDENptr = &RF_lockDENfens;        
        obsSize = (mode == RF_PRED) ? RF_fobservationSize : RF_observationSize;
      }
      *lockDENptr = ompvector(1, obsSize);
      for (i = 1; i <= obsSize; i++) {
        omp_init_lock(&((*lockDENptr)[i]));
      }
      if (oobFlag == TRUE) {
        oobFlag = FALSE;
      }
      else {
        fullFlag = FALSE;
      }
    }
  }
  else {
    char  potentiallyMixedMultivariate = FALSE;
    if (RF_rTargetFactorCount > 0) {
      omp_lock_t   **lockDENptr;
      uint obsSize;  
      char oobFlag, fullFlag;
      oobFlag = fullFlag = FALSE;
      switch (mode) {
      case RF_PRED:
        if (RF_opt & OPT_FENS) {
          fullFlag = TRUE;
        }
        break;
      default:
        if (RF_opt & OPT_OENS) {
          oobFlag = TRUE;
        }
        if (RF_opt & OPT_FENS) {
          fullFlag = TRUE;
        }
        break;
      }
      while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
        if (oobFlag == TRUE) {
          lockDENptr = &RF_lockDENoens;
          obsSize = RF_observationSize;
        }
        else {
          lockDENptr = &RF_lockDENfens;
          obsSize = (mode == RF_PRED) ? RF_fobservationSize : RF_observationSize;
        }
        if (!potentiallyMixedMultivariate) {
          *lockDENptr = ompvector(1, obsSize);
          for (i = 1; i <= obsSize; i++) {
            omp_init_lock(&((*lockDENptr)[i]));
          }
        }
        if (oobFlag == TRUE) {
          oobFlag = FALSE;
        }
        else {
          fullFlag = FALSE;
        }
      }
      potentiallyMixedMultivariate = TRUE;
    }
    if (RF_rTargetNonFactorCount > 0) {
      omp_lock_t   **lockDENptr;
      omp_lock_t   **lockQNTptr;
      uint obsSize;  
      char oobFlag, fullFlag;
      oobFlag = fullFlag = FALSE;
      switch (mode) {
      case RF_PRED:
        if (RF_opt & OPT_FENS) {
          fullFlag = TRUE;
        }
        break;
      default:
        if (RF_opt & OPT_OENS) {
          oobFlag = TRUE;
        }
        if (RF_opt & OPT_FENS) {
          fullFlag = TRUE;
        }
        break;
      }
      while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
        if (oobFlag == TRUE) {
          lockDENptr = &RF_lockDENoens;
          lockQNTptr = &RF_lockQNToens;
          obsSize = RF_observationSize;
        }
        else {
          lockDENptr = &RF_lockDENfens;
          lockQNTptr = &RF_lockQNTfens;
          obsSize = (mode == RF_PRED) ? RF_fobservationSize : RF_observationSize;
        }
        if (!potentiallyMixedMultivariate) {
          *lockDENptr = ompvector(1, obsSize);
          for (i = 1; i <= obsSize; i++) {
            omp_init_lock(&((*lockDENptr)[i]));
          }
        }
        *lockQNTptr = ompvector(1, obsSize);        
        for (i = 1; i <= obsSize; i++) {
          omp_init_lock(&((*lockQNTptr)[i]));
        }
        if (oobFlag == TRUE) {
          oobFlag = FALSE;
        }
        else {
          fullFlag = FALSE;
        }
      }
      potentiallyMixedMultivariate = TRUE;
    }
  }
}
#else
void stackLocksOpenMP(char mode) { }
#endif
#ifdef _OPENMP
void unstackLocksOpenMP(char mode) {
  uint i, j;
  omp_destroy_lock(&RF_lockEnsbUpdtCount);
  omp_destroy_lock(&RF_lockPerf);
  if (RF_optHigh & OPT_PART_PLOT) {
    for (i = 1; i <= RF_observationSize; i++) {
      omp_destroy_lock(&(RF_lockPartial[i]));
    }
    free_ompvector(RF_lockPartial, 1, RF_observationSize);
  }
  if (RF_optHigh & OPT_WGHT) {
    uint gMembershipSize;
    if((RF_optHigh & OPT_WGHT_IBG) && (RF_optHigh & OPT_WGHT_OOB)) {
      switch (mode) {
      case RF_PRED:
        gMembershipSize = RF_fobservationSize;
        break;
      default:
        gMembershipSize = RF_observationSize;
        break;
      }
    }
    else {
      gMembershipSize  = RF_observationSize;
    }
    for (i = 1; i <= gMembershipSize; i++) {
      for (j = 1; j <= RF_observationSize; j++) {
        omp_destroy_lock(&(RF_lockWeight[i][j]));
      }
    }
    for (i = 1; i <= gMembershipSize; i++) {    
      free_ompvector(RF_lockWeight[i], 1, RF_observationSize);
    }
    free_new_vvector(RF_lockWeight, 1, gMembershipSize, NRUTIL_OMPLPTR);
    for (i = 1; i <= gMembershipSize; i++) {
      omp_destroy_lock(&(RF_lockWeightRow[i]));
    }
    free_ompvector(RF_lockWeightRow, 1, gMembershipSize);
  }
  if (RF_opt & OPT_VIMP) {
      uint obsSize, xVimpSize;
      if (RF_opt & OPT_VIMP_JOIN) {
        xVimpSize = 1;
      }
      else {
        xVimpSize = RF_intrPredictorSize;
      }
      switch (mode) {
      case RF_PRED:
        obsSize = RF_fobservationSize;
        break;
      default:
        obsSize  = RF_observationSize;
        break;
      }
      for (i = 1; i <= xVimpSize; i++) {
        for (j = 1; j <= obsSize; j++) {
          omp_destroy_lock(&(RF_lockVimp[i][j]));
        }
      }
      for (i = 1; i <= xVimpSize; i++) {
        free_ompvector(RF_lockVimp[i], 1, obsSize);
      }
      free_new_vvector(RF_lockVimp, 1, xVimpSize, NRUTIL_OMPLPTR);
      for (i = 1; i <= obsSize; i++) {
        omp_destroy_lock(&(RF_lockVimpRow[i]));
      }
      free_ompvector(RF_lockVimpRow, 1, obsSize);
      for (i = 1; i <= xVimpSize; i++) {
        omp_destroy_lock(&(RF_lockVimpCol[i]));
      }
      free_ompvector(RF_lockVimpCol, 1, xVimpSize);
  }
  if ((RF_vtry > 0) && (RF_vtryMode != RF_VTRY_NULL)) {
    uint xVimpSize;
    xVimpSize = RF_xSize;
    for (i = 1; i <= xVimpSize; i++) {
      if (RF_holdBLKptr[i] > 0) {
        for (j = 1; j <= RF_holdBLKptr[i]; j++) {
          omp_destroy_lock(&(RF_lockVimpHoldout[i][j]));
        }
        free_ompvector(RF_lockVimpHoldout[i], 1, RF_holdBLKptr[i]);
      }
    }
    free_new_vvector(RF_lockVimpHoldout, 1, xVimpSize, NRUTIL_OMPLPTR);    
  }
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    omp_lock_t   *lockDENptr;
    uint obsSize;  
    char oobFlag, fullFlag;
    oobFlag = fullFlag = FALSE;
    switch (mode) {
    case RF_PRED:
      if (RF_opt & OPT_FENS) {
        fullFlag = TRUE;
      }
      break;
    default:
      if (RF_opt & OPT_OENS) {
        oobFlag = TRUE;
      }
      if (RF_opt & OPT_FENS) {
        fullFlag = TRUE;
      }
      break;
    }
    while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
      if (oobFlag == TRUE) {
        lockDENptr = RF_lockDENoens;
        obsSize = RF_observationSize;
      }
      else {
        lockDENptr = RF_lockDENfens;
        obsSize = (mode == RF_PRED) ? RF_fobservationSize : RF_observationSize;
      }
      for (i = 1; i <= obsSize; i++) {
        omp_destroy_lock(&(lockDENptr[i]));
      }
      free_ompvector(lockDENptr, 1, obsSize);
      if (oobFlag == TRUE) {
        oobFlag = FALSE;
      }
      else {
        fullFlag = FALSE;
      }
    }
  }
  else {
    char  potentiallyMixedMultivariate = FALSE;
    if (RF_rTargetFactorCount > 0) {
      omp_lock_t   *lockDENptr;
      uint obsSize;  
      char oobFlag, fullFlag;
      oobFlag = fullFlag = FALSE;
      switch (mode) {
      case RF_PRED:
        if (RF_opt & OPT_FENS) {
          fullFlag = TRUE;
        }
        break;
      default:
        if (RF_opt & OPT_OENS) {
          oobFlag = TRUE;
        }
        if (RF_opt & OPT_FENS) {
          fullFlag = TRUE;
        }
        break;
      }
      while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
        if (oobFlag == TRUE) {
          lockDENptr = RF_lockDENoens;
          obsSize = RF_observationSize;
        }
        else {
          lockDENptr = RF_lockDENfens;
          obsSize = (mode == RF_PRED) ? RF_fobservationSize : RF_observationSize;
        }
        if (!potentiallyMixedMultivariate) {
          for (i = 1; i <= obsSize; i++) {
            omp_destroy_lock(&(lockDENptr[i]));
          }
          free_ompvector(lockDENptr, 1, obsSize);
        }
        if (oobFlag == TRUE) {
          oobFlag = FALSE;
        }
        else {
          fullFlag = FALSE;
        }
      }
      potentiallyMixedMultivariate = TRUE;
    }
    if (RF_rTargetNonFactorCount > 0) {
      omp_lock_t   *lockDENptr;
      omp_lock_t   *lockQNTptr;
      uint obsSize;  
      char oobFlag, fullFlag;
      oobFlag = fullFlag = FALSE;
      switch (mode) {
      case RF_PRED:
        if (RF_opt & OPT_FENS) {
          fullFlag = TRUE;
        }
        break;
      default:
        if (RF_opt & OPT_OENS) {
          oobFlag = TRUE;
        }
        if (RF_opt & OPT_FENS) {
          fullFlag = TRUE;
        }
        break;
      }
      while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
        if (oobFlag == TRUE) {
          lockDENptr = RF_lockDENoens;
          lockQNTptr = RF_lockQNToens;
          obsSize = RF_observationSize;
        }
        else {
          lockDENptr = RF_lockDENfens;
          lockQNTptr = RF_lockQNTfens;
          obsSize = (mode == RF_PRED) ? RF_fobservationSize : RF_observationSize;
        }
        if (!potentiallyMixedMultivariate) {
          for (i = 1; i <= obsSize; i++) {
            omp_destroy_lock(&(lockDENptr[i]));
          }
          free_ompvector(lockDENptr, 1, obsSize);
        }
          for (i = 1; i <= obsSize; i++) {
            omp_destroy_lock(&(lockQNTptr[i]));
          }
          free_ompvector(lockQNTptr, 1, obsSize);
          if (oobFlag == TRUE) {
          oobFlag = FALSE;
        }
        else {
          fullFlag = FALSE;
        }
      }
      potentiallyMixedMultivariate = TRUE;
    }
  }
}
#else
void unstackLocksOpenMP(char mode) { }
#endif
