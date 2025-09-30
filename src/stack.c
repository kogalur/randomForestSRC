
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "stack.h"
#include "impute.h"
#include "nrutil.h"
#include "error.h"
/*
#include "nativeUtil.h"
*/
void stackAndInitializeTimeAndSubjectArrays(char mode) {
  uint i, j;
  uint leadingIndex;
  uint adjObsSize;
  if (!(RF_opt & OPT_ANON)) {
    RF_masterTime  = dvector(1, RF_observationSize);
    RF_masterTimeIndexIn  = uivector(1, RF_observationSize);
    RF_masterTimeSize = 0;
    for (j = 1; j <= RF_observationSize; j++) {
      if (!RF_nativeIsNaN(RF_responseIn[RF_timeIndex][j])) {
        RF_masterTimeSize ++;
        RF_masterTime[RF_masterTimeSize] = RF_responseIn[RF_timeIndex][j];
      }
    }
    adjObsSize = RF_observationSize;
    qksort(RF_masterTime, RF_masterTimeSize);
    leadingIndex = 1;
    for (i=2; i <= RF_masterTimeSize; i++) {
      if (RF_masterTime[i] > RF_masterTime[leadingIndex]) {
        leadingIndex++;
        RF_masterTime[leadingIndex] = RF_masterTime[i];
      }
    }
    RF_masterTimeSize = leadingIndex;
    for (i= RF_masterTimeSize + 1; i <= adjObsSize; i++) {
      RF_masterTime[i] = 0;
    }
  }
  if (!(RF_opt & OPT_IMPU_ONLY)) {
    qksort(RF_timeInterest, RF_timeInterestSize);
    RF_sortedTimeInterestSize = 1;
    for (i=2; i <= RF_timeInterestSize; i++) {
      if (RF_timeInterest[i] > RF_timeInterest[RF_sortedTimeInterestSize]) {
        (RF_sortedTimeInterestSize) ++;
        RF_timeInterest[RF_sortedTimeInterestSize] = RF_timeInterest[i];
      }
    }
    if (RF_sortedTimeInterestSize != RF_timeInterestSize) {
      RF_nativePrint("\nRF-SRC:  *** WARNING *** ");
      RF_nativePrint("\nRF-SRC:  Time points of interest are not unique.");
      RF_nativePrint("\nRF-SRC:  Any ensemble matricies will be");
      RF_nativePrint("\nRF-SRC:  resized as [N'] x [n], where N' is the");
      RF_nativePrint("\nRF-SRC:  unique time points of interest and n is");
      RF_nativePrint("\nRF-SRC:  number of observations in the data.");
    }
    for (i = RF_sortedTimeInterestSize + 1; i <= RF_timeInterestSize; i++) {
      RF_timeInterest[i] = 0;
    }
  }  
}
void unstackTimeAndSubjectArrays(char mode) {
  if (!(RF_opt & OPT_ANON)) {
    free_dvector(RF_masterTime, 1, RF_observationSize);
    free_uivector(RF_masterTimeIndexIn, 1, RF_observationSize);
  }
}
void stackFactorArrays(char mode) {
  uint i, k;
  stackFactorGeneric(TRUE,
                     RF_ySize,
                     RF_rType,
                     &RF_rFactorMap,
                     &RF_rFactorCount,
                     &RF_rFactorIndex,
                     &RF_rFactorSize,
                     &RF_rNonFactorMap,
                     &RF_rNonFactorCount,
                     &RF_rNonFactorIndex);
  stackFactorGeneric(FALSE,
                     RF_xSize,
                     RF_xType,
                     &RF_xFactorMap,
                     &RF_xFactorCount,
                     &RF_xFactorIndex,
                     &RF_xFactorSize,
                     &RF_xNonFactorMap,
                     &RF_xNonFactorCount,
                     &RF_xNonFactorIndex);
  if (RF_xFactorCount > 0) {
    
    RF_xLevels = (uint **) new_vvector(1, RF_xFactorCount, NRUTIL_UPTR);  
    for (k = 1; k <= RF_xFactorCount; k++) {
      if (RF_xLevelsCnt[k] > 0) {
        RF_xLevels[k] = (uint *) INTEGER(VECTOR_ELT(RF_xLevelsSEXP, k-1));
        RF_xLevels[k] --;
      }
      else {
        RF_nativeError("\nRF-SRC: *** ERROR *** ");
        RF_nativeError("\nRF-SRC: Inconsistent zero-level count in factor:  compressed-index = %10d, x-index = %10d", k, RF_xFactorIndex[k]);
        RF_nativeError("\nRF-SRC: Please Contact Technical Support.");
        RF_nativeExit();
      }
    }
    
    /*
    RF_xLevels = (uint **) copy2DObject(RF_xLevelsJNIE, NATIVE_TYPE_INTEGER, &RF_nat2DInfoListSize);
    for (k = 1; k <= RF_xFactorCount; k++) {    
      if (RF_xLevelsCnt[k] > 0) {
      }
      else {
        RF_nativeError("\nRF-SRC: *** ERROR *** ");
        RF_nativeError("\nRF-SRC: Inconsistent zero-level count in factor:  compressed-index = %10d, x-index = %10d", k, RF_xFactorIndex[k]);
        RF_nativeError("\nRF-SRC: Please Contact Technical Support.");
        RF_nativeExit();
      }
    }
    */    
    /*
    RF_xLevels = (uint **) copyXDObject(RF_xLevelsPYTH, &RF_natXDInfoListSize);
    for (k = 1; k <= RF_xFactorCount; k++) {
      if (RF_xLevelsCnt[k] > 0) {
      }
      else {
        RF_nativeError("\nRF-SRC: *** ERROR *** ");
        RF_nativeError("\nRF-SRC: Inconsistent zero-level count in factor:  compressed-index = %10d, x-index = %10d", k, RF_xFactorIndex[k]);
        RF_nativeError("\nRF-SRC: Please Contact Technical Support.");
        RF_nativeExit();
      }
    }
    */    
  }
  if (RF_ySize == 0) {
    RF_rTarget = NULL;
    RF_rTargetCount = 0;
  }
  else {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      RF_rTarget = NULL;
      RF_rTargetCount = 0;
    }
    else {
      if (mode == RF_GROW) {
        RF_rTargetCount = RF_ySize;
        RF_rTarget = uivector(1 , RF_rTargetCount);
        for (i = 1; i <= RF_rTargetCount; i++) {
          RF_rTarget[i] = i;
        }
      }
      else {
      }
      RF_rTargetFactor    = uivector(1, RF_rTargetCount);
      RF_rTargetNonFactor = uivector(1, RF_rTargetCount);
      RF_rTargetFactorCount = RF_rTargetNonFactorCount = 0;
      for (i = 1; i <= RF_rTargetCount; i++) {
        if ((RF_rTarget[i] < 1) || (RF_rTarget[i] > RF_ySize)) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  Target response is out of range for [C+], [R+], [M+]:  %10d %10d ", i, RF_rTarget[i]);
          RF_nativeExit();
        }
        if ((RF_rType[RF_rTarget[i]] == 'B') ||
            (RF_rType[RF_rTarget[i]] == 'I') ||
            (RF_rType[RF_rTarget[i]] == 'C')) {
          RF_rTargetFactor[++RF_rTargetFactorCount] = RF_rTarget[i];
        }
        else {
          RF_rTargetNonFactor[++RF_rTargetNonFactorCount] = RF_rTarget[i];
        }
      }
    }  
  }  
}
void stackFactorGeneric(char    respFlag,
                        uint    size,
                        char   *type,
                        uint  **p_factorMap,
                        uint   *factorCount,
                        uint  **p_factorIndex,
                        uint  **p_factorSize,
                        uint  **p_nonfactorMap,
                        uint   *nonfactorCount,
                        uint  **p_nonfactorIndex) {
  uint i, j;
  if (size > 0) {
    *p_factorMap    = uivector(1, size);
    *p_nonfactorMap = uivector(1, size);
    *factorCount    = 0;
    *nonfactorCount = 0;
    for (i = 1; i <= size; i++) {
      (*p_factorMap)[i]    = 0;
      (*p_nonfactorMap)[i] = 0;
      if ((type[i] == 'B') ||
          ((type[i] == 'I') && respFlag) ||
          (type[i] == 'C')) {
        (*factorCount) ++;
        (*p_factorMap)[i] = *factorCount;
      }
      else {
        (*nonfactorCount) ++;
        (*p_nonfactorMap)[i] = *nonfactorCount;
      }
    }
    if (*factorCount > 0) {
      *p_factorIndex = uivector(1, *factorCount);
      j = 0;
      for (i = 1; i <= size; i++) {
        if ((*p_factorMap)[i] > 0) {
          (*p_factorIndex)[++j] = i;
        }
      }
      *p_factorSize = uivector(1, *factorCount);
    }
    if (*nonfactorCount > 0) {
      *p_nonfactorIndex = uivector(1, *nonfactorCount);
      j = 0;
      for (i = 1; i <= size; i++) {
        if ((*p_nonfactorMap)[i] > 0) {
          (*p_nonfactorIndex)[++j] = i;
        }
      }
    }
  }
  else {
    *factorCount    = 0;
    *nonfactorCount = 0;
  }
}
void unstackFactorArrays(char mode) {
  if (RF_ySize > 0) {
    free_uivector(RF_rFactorMap, 1, RF_ySize);
    if (RF_rFactorCount > 0) {
      free_uivector(RF_rFactorIndex, 1, RF_rFactorCount);
      free_uivector(RF_rFactorSize, 1, RF_rFactorCount);
    }
    free_uivector(RF_rNonFactorMap, 1, RF_ySize);
    if (RF_rNonFactorCount > 0) {
      free_uivector(RF_rNonFactorIndex, 1, RF_rNonFactorCount);
    }
  }
  free_uivector(RF_xFactorMap, 1, RF_xSize);
  if (RF_xFactorCount > 0) {
    free_uivector(RF_xFactorIndex, 1, RF_xFactorCount);
    free_uivector(RF_xFactorSize, 1, RF_xFactorCount);
    free_new_vvector(RF_xLevels, 1, RF_xFactorCount, NRUTIL_UPTR);  
  }
  free_uivector(RF_xNonFactorMap, 1, RF_xSize);
  if (RF_xNonFactorCount > 0) {
    free_uivector(RF_xNonFactorIndex, 1, RF_xNonFactorCount);
  }
  if ((RF_rFactorCount + RF_xFactorCount) > 0) {
    free_new_vvector(RF_factorList, 1, RF_ntree, NRUTIL_FPTR2);
  }
  if (RF_ySize == 0) {
  }
  else {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    }
    else {
      free_uivector(RF_rTargetFactor, 1, RF_rTargetCount);
      free_uivector(RF_rTargetNonFactor, 1, RF_rTargetCount);
      if (mode == RF_GROW) {
        free_uivector(RF_rTarget, 1 , RF_rTargetCount);
      }
      else {
      }
    }
  }
}
char stackMissingArraysPhase1(char mode) {
  char result;
  char mFlag;
  uint i, j;
  result = TRUE;
  if (!(RF_opt & OPT_ANON)) {
    if (!(RF_optHigh & OPT_DATA_PASG)) {
      for (j = 1; j <= RF_ySize; j++) {
        for (i = 1; i <= RF_observationSize; i++) {
          if (!RF_nativeIsNaN(RF_responseIn[j][i])) {
            if (!isfinite(RF_responseIn[j][i])) {
              result = FALSE;
              RF_nativePrint("\nRF-SRC:  train response elements must not be plus or minus infinity:  [%10d, %10d] = %12.4f \n", j, i, RF_responseIn[j][i]);
            }
          }
        }
        if (result == FALSE) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  Plus or Minus Infinity detected.");
          RF_nativeExit();
        }
      }  
      for (j = 1; j <= RF_xSize; j++) {
        for (i = 1; i <= RF_observationSize; i++) {
          if (!RF_nativeIsNaN(RF_observationIn[j][i])) {
            if (!isfinite(RF_observationIn[j][i])) {
              result = FALSE;
              RF_nativeError("\nRF-SRC:  train x-variable elements must not be plus or minus infinity:  [%10d, %10d] = %12.4f \n", j, i, RF_observationIn[j][i]);
            }
          }
        }
        if (result == FALSE) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  Plus or Minus Infinity detected.");
          RF_nativeExit();
        }
      }
    }
    if (!(RF_optHigh & OPT_DATA_PASG)) {
      for (j = 1; j <= RF_ySize; j++) {
        if (j == RF_timeIndex) {
          for (i = 1; i <= RF_observationSize; i++) {
            if (!RF_nativeIsNaN(RF_responseIn[RF_timeIndex][i])) {
              if (RF_responseIn[RF_timeIndex][i] < 0) {
                result = FALSE;
                RF_nativePrint("\nRF-SRC:  train time elements must be greater than or equal to zero or NA:  [%10d] = %12.4f \n", i, RF_responseIn[RF_timeIndex][i]);
              }
            }
          }
        }
        if (j == RF_statusIndex) {
          for (i = 1; i <= RF_observationSize; i++) {
            if (!RF_nativeIsNaN(RF_responseIn[RF_statusIndex][i])) {
              if (RF_responseIn[RF_statusIndex][i] < 0) {
                result = FALSE;
                RF_nativePrint("\nRF-SRC:  train status elements must be greater than or equal to zero or NA:  [%10d] = %12.4f \n", i, RF_responseIn[RF_statusIndex][i]);
              }
            }
          }
        }
        if (j == RF_statusIndex) {
          mFlag = FALSE;
          for (i = 1; i <= RF_observationSize; i++) {
            if (!RF_nativeIsNaN(RF_responseIn[RF_statusIndex][i])) {
              if (RF_responseIn[RF_statusIndex][i] >= 0) {
                mFlag = TRUE;
                i = RF_observationSize;
              }
            }
          }
          if (mFlag == FALSE) {
            RF_nativePrint("\nRF-SRC:  All train status elements are censored or missing. \n");
            result = FALSE;
          }
        }
        if ((RF_statusIndex == 0) && (RF_timeIndex == 0)) {
          mFlag = FALSE;
          for (i = 1; i <= RF_observationSize; i++) {
            if (!RF_nativeIsNaN(RF_responseIn[j][i])) {
              mFlag = TRUE;
              i = RF_observationSize;
            }
          }
          if (mFlag == FALSE) {
            RF_nativePrint("\nRF-SRC:  All train outcome/response elements are missing for:  %10d \n", j);
            result = FALSE;
          }
        }
        if (result == FALSE) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  Missingness verification failed.");
          RF_nativeExit();
        }
      }  
    }
    RF_response = (double ***) new_vvector(1, RF_ntree, NRUTIL_DPTR2);
    if (RF_ySize > 0) {
      for (i = 1 ; i <= RF_ntree; i++) {
        RF_response[i] = RF_responseIn;
      }
      RF_time = NULL;
      RF_masterTimeIndex = NULL;
      if (RF_statusIndex > 0) {
        RF_time = (double **) new_vvector(1, RF_ntree, NRUTIL_DPTR);
        RF_masterTimeIndex = (uint **) new_vvector(1, RF_ntree, NRUTIL_UPTR);
        for (i = 1 ; i <= RF_ntree; i++) {
          RF_time[i] = RF_responseIn[RF_timeIndex];
          RF_masterTimeIndex[i] = RF_masterTimeIndexIn;
        }
        updateTimeIndexArray(0,
                             NULL,
                             RF_observationSize,
                             RF_responseIn[RF_timeIndex],
                             TRUE,
                             FALSE,
                             RF_masterTimeIndexIn);
        RF_status = (double **) new_vvector(1, RF_ntree, NRUTIL_DPTR);
        for (i = 1 ; i <= RF_ntree; i++) {
          RF_status[i] = RF_responseIn[RF_statusIndex];
        }
      }
    }
    else {
      for (i = 1 ; i <= RF_ntree; i++) {
        RF_response[i] = NULL;
      }
    }
    RF_observation = (double ***) new_vvector(1, RF_ntree, NRUTIL_DPTR2);
    for (i = 1 ; i <= RF_ntree; i++) {
      RF_observation[i] = RF_observationIn;
    }
    if (RF_optHigh & OPT_DATA_PASG) {
      RF_mStatusFlag = RF_mTimeFlag = RF_mResponseFlag = RF_mPredictorFlag = FALSE;
      RF_mRecordSize = 0;
      RF_mRecordMap = NULL;
    }
    else {
      RF_mRecordMap = uivector(1, RF_observationSize);
      RF_mRecordSize = getRecordMap(RF_mRecordMap,
                                    RF_observationSize,
                                    RF_responseIn,
                                    RF_observationIn);
      if (RF_mRecordSize == 0) {
        RF_mStatusFlag = RF_mTimeFlag = RF_mResponseFlag = RF_mPredictorFlag = FALSE;
      }
      else {
        RF_optHigh = RF_optHigh & (~OPT_MEMB_INCG);
        RF_optHigh = RF_optHigh & (~OPT_TERM_INCG);
        stackMissingSignatures(RF_observationSize,
                               RF_ySize,
                               RF_responseIn,
                               RF_observationIn,
                               RF_mRecordMap,
                               RF_mRecordSize,
                               & RF_mRecordIndex,
                               & RF_mpIndexSize,
                               & RF_mpSign,
                               & RF_mpIndex,
                               & RF_mrFactorSize,
                               & RF_mrFactorIndex,
                               & RF_mxFactorSize,
                               & RF_mxFactorIndex,
                               & RF_mTimeFlag,
                               & RF_mStatusFlag,
                               & RF_mResponseFlag,
                               & RF_mPredictorFlag);
        if (RF_mResponseFlag == TRUE) {
          for (i = 1 ; i <= RF_ntree; i++) {
            RF_response[i] = NULL;
            if (RF_timeIndex > 0) {
              RF_time[i] = NULL;
              RF_masterTimeIndex[i] = NULL;
            }
            if (RF_statusIndex > 0) {
              RF_status[i] = NULL;
            }
          }
        }
        if (RF_mPredictorFlag == TRUE) {
          for (i = 1 ; i <= RF_ntree; i++) {
            RF_observation[i] = NULL;
          }
        }
      }  
    }  
  }  
  else {
    RF_mStatusFlag = RF_mTimeFlag = RF_mResponseFlag = RF_mPredictorFlag = FALSE;
    RF_mRecordSize = 0;
    RF_mRecordMap = NULL;
  }
  if (mode == RF_PRED) {
    if (!(RF_optHigh & OPT_DATA_PASP)) {
      if (RF_frSize > 0) {
        for (j = 1; j <= RF_ySize; j++) {      
          for (i = 1 ; i <= RF_fobservationSize; i++) {
            if (!RF_nativeIsNaN(RF_fresponseIn[j][i])) {
              if (!isfinite(RF_fresponseIn[j][i])) {
                result = FALSE;
                RF_nativeError("\nRF-SRC:  test response elements must not be plus or minus infinity:  [%10d, %10d] = %12.4f \n", j, i, RF_fresponseIn[j][i]);
              }
            }
          }
          if (result == FALSE) {
            RF_nativeError("\nRF-SRC:  *** ERROR *** ");
            RF_nativeError("\nRF-SRC:  Plus or Minus Infinity detected.");
            RF_nativeExit();
          }
        }
      }
      for (j = 1; j <= RF_xSize; j++) {      
        for (i = 1; i <= RF_fobservationSize; i++) {
          if (!RF_nativeIsNaN(RF_fobservationIn[j][i])) {
            if (!isfinite(RF_fobservationIn[j][i])) {
              result = FALSE;
              RF_nativeError("\nRF-SRC:  test x-variable elements must not be plus or minus infinity:  [%10d, %10d] = %12.4f \n", j, i, RF_fobservationIn[j][i]);
            }
          }
        }
        if (result == FALSE) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  Plus or Minus Infinity detected.");
          RF_nativeExit();
        }
      }
    }
    if (!(RF_optHigh & OPT_DATA_PASP)) {    
      if (RF_frSize > 0) {
        if (RF_timeIndex > 0) {
          for (i = 1 ; i <= RF_fobservationSize; i++) {
            if (!RF_nativeIsNaN(RF_fresponseIn[RF_timeIndex][i])) {
              if (RF_fresponseIn[RF_timeIndex][i] < 0) {
                result = FALSE;
                RF_nativeError("\nRF-SRC:  PRED time elements must be greater than or equal to zero or NA:  [%10d] = %12.4f \n", i, RF_fresponseIn[RF_timeIndex][i]);
              }
            }
          }
        }
      }
      if (RF_frSize > 0) {
        if (RF_statusIndex > 0) {
          for (i = 1 ; i <= RF_fobservationSize; i++) {
            if (!RF_nativeIsNaN(RF_fresponseIn[RF_statusIndex][i])) {
              if (RF_fresponseIn[RF_statusIndex][i] < 0) {
                result = FALSE;
                RF_nativeError("\nRF-SRC:  PRED status elements must be greater than or equal to zero or NA:  [%10d] = %12.4f \n", i, RF_fresponseIn[RF_statusIndex][i]);
              }
            }
          }
        }
      }
      if (result == FALSE) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Missingness verification failed.");
        RF_nativeExit();
      }
    }
    RF_fobservation = (double ***) new_vvector(1, RF_ntree, NRUTIL_DPTR2);
    for (i = 1 ; i <= RF_ntree; i++) {
      RF_fobservation[i] = RF_fobservationIn;
    }
    RF_fresponse = (double ***) new_vvector(1, RF_ntree, NRUTIL_DPTR2);
    if (RF_frSize > 0) {
      for (i = 1 ; i <= RF_ntree; i++) {
        RF_fresponse[i] = RF_fresponseIn;
      }
      if (RF_timeIndex > 0) {
        RF_ftime = (double **) new_vvector(1, RF_ntree, NRUTIL_DPTR);
        for (i = 1 ; i <= RF_ntree; i++) {
          RF_ftime[i] = RF_fresponseIn[RF_timeIndex];
        }
      }
      if (RF_statusIndex > 0) {
        RF_fstatus = (double **) new_vvector(1, RF_ntree, NRUTIL_DPTR);
        for (i = 1 ; i <= RF_ntree; i++) {
          RF_fstatus[i] = RF_fresponseIn[RF_statusIndex];
        }
      }
    }
    else {
      for (i = 1 ; i <= RF_ntree; i++) {
        RF_fresponse[i] = NULL;
      }
    }
    if (RF_optHigh & OPT_DATA_PASP) {
      RF_fmStatusFlag = RF_fmTimeFlag = RF_fmResponseFlag = RF_fmPredictorFlag = FALSE;
      RF_fmRecordSize = 0;
      RF_fmRecordMap = NULL;
    }
    else {
    RF_fmRecordMap = uivector(1, RF_fobservationSize);
    RF_fmRecordSize = getRecordMap(RF_fmRecordMap,
                                 RF_fobservationSize,
                                 RF_fresponseIn,
                                 RF_fobservationIn);
    if (RF_fmRecordSize == 0) {
      RF_fmStatusFlag = RF_fmTimeFlag = RF_fmResponseFlag = RF_fmPredictorFlag = FALSE;
    }  
    else {
      if (RF_opt & OPT_ANON) {
      }
      else {
      }
      stackMissingSignatures(RF_fobservationSize,
                             RF_frSize,
                             RF_fresponseIn,
                             RF_fobservationIn,
                             RF_fmRecordMap,
                             RF_fmRecordSize,
                             & RF_fmRecordIndex,
                             & RF_fmpIndexSize,
                             & RF_fmpSign,
                             & RF_fmpIndex,
                             & RF_fmrFactorSize,
                             & RF_fmrFactorIndex,
                             & RF_fmxFactorSize,
                             & RF_fmxFactorIndex,
                             & RF_fmTimeFlag,
                             & RF_fmStatusFlag,
                             & RF_fmResponseFlag,
                             & RF_fmPredictorFlag);
      if (FALSE) {
        if (RF_frSize > 0) {
          if (RF_fmResponseFlag == TRUE) {
            for (i = 1 ; i <= RF_ntree; i++) {
              RF_fresponse[i] = NULL;
              if (RF_timeIndex > 0) {
                RF_ftime[i] = NULL;
              }
              if (RF_statusIndex > 0) {
                RF_fstatus[i] = NULL;
              }
            }
          }
        }
        if (RF_fmPredictorFlag == TRUE) {
          for (i = 1 ; i <= RF_ntree; i++) {
            RF_fobservation[i] = NULL;
          }
        }
      }
    }  
    }  
  }  
  return result;
}
char stackMissingArraysPhase2(char mode) {
  char result;
  char mFlag;
  char dualUseFlag;
  uint recordSize;
  uint i, j;
  result = TRUE;
  if (RF_opt & OPT_ANON) {
    result = FALSE;
    if (RF_fmResponseFlag == TRUE) {
      RF_opt = RF_opt & (~OPT_PERF);
    }
    if (RF_fmPredictorFlag == TRUE) {
      if (RF_optHigh & OPT_JIT_TOP) {
      }
      else {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  An anonymous forest with missingness in the test data requires the JITT flag to be asserted");
        RF_nativeError("\nRF-SRC:  Please adjust your script accordingly.");
        RF_nativeExit();
      }      
    }
  }
  if (RF_optHigh & OPT_JIT_TOP) {
    result = FALSE;
  }
  if (result) {
    dualUseFlag = FALSE;
    switch (mode) {
    case RF_PRED:
      if (RF_fmRecordSize > 0) {
        recordSize = RF_fmRecordSize;
        dualUseFlag = TRUE;
        mFlag = ACTIVE;
      }
      else {
        RF_opt = RF_opt & (~OPT_MISS_OUT);
      }
      break;
    default:
      RF_fmRecordSize = 0;
      if (RF_mRecordSize > 0) {
        recordSize = RF_mRecordSize;
        dualUseFlag = TRUE;
        mFlag = FALSE;
      }
      else {
        RF_opt = RF_opt & (~OPT_MISS_OUT);
        RF_nImpute = 1;
      }
      break;
    }  
    if (dualUseFlag == TRUE) {
      RF_dmRecordBootFlag = cmatrix(1, RF_ntree, 1, recordSize);
      for (j = 1; j <= RF_ntree; j++) {
        for (i = 1; i <= recordSize; i++) {
          RF_dmRecordBootFlag[j][i] = mFlag;
        }
      }
    }
  }  
  else {
    RF_opt = RF_opt & (~OPT_MISS_OUT);    
  }
  return result;
}
void unstackMissingArrays(char mode) {
  char dualUseFlag;
  uint recordSize;
  if (!(RF_opt & OPT_ANON)) {  
    free_new_vvector(RF_response, 1, RF_ntree, NRUTIL_DPTR2);
    if (RF_ySize > 0) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        free_new_vvector(RF_time, 1, RF_ntree, NRUTIL_DPTR);
        free_new_vvector(RF_masterTimeIndex, 1, RF_ntree, NRUTIL_UPTR);
        free_new_vvector(RF_status, 1, RF_ntree, NRUTIL_DPTR);
      }
    }
    free_new_vvector(RF_observation, 1, RF_ntree, NRUTIL_DPTR2);
    if (RF_optHigh & OPT_DATA_PASG) {
    }
    else {
      free_uivector(RF_mRecordMap, 1, RF_observationSize);
    }
    if (RF_mRecordSize == 0) {
    }
    else {
      unstackMissingSignatures(RF_ySize,
                               RF_mRecordSize,
                               RF_mRecordIndex,
                               RF_mpIndexSize,
                               RF_mpSign,
                               RF_mpIndex,
                               RF_mrFactorSize,
                               RF_mrFactorIndex,
                               RF_mxFactorSize,
                               RF_mxFactorIndex);
    }
  }  
  else {
  }  
  if (mode == RF_PRED) {
    free_new_vvector(RF_fobservation, 1, RF_ntree, NRUTIL_DPTR2);
    free_new_vvector(RF_fresponse, 1, RF_ntree, NRUTIL_DPTR2);
    if (RF_frSize > 0) {
      if (RF_timeIndex > 0) {
        free_new_vvector(RF_ftime, 1, RF_ntree, NRUTIL_DPTR);
      }
      if (RF_statusIndex > 0) {
        free_new_vvector(RF_fstatus, 1, RF_ntree, NRUTIL_DPTR);
      }
    }
    if (RF_optHigh & OPT_DATA_PASP) {
    }
    else {
      free_uivector(RF_fmRecordMap, 1, RF_fobservationSize);
    }
    if (RF_fmRecordSize == 0) {
    }
    else {
      unstackMissingSignatures(RF_frSize,
                               RF_fmRecordSize,
                               RF_fmRecordIndex,
                               RF_fmpIndexSize,
                               RF_fmpSign,
                               RF_fmpIndex,
                               RF_fmrFactorSize,
                               RF_fmrFactorIndex,
                               RF_fmxFactorSize,
                               RF_fmxFactorIndex);
    }
  }
  if (RF_opt & OPT_MISS_OUT) {
    dualUseFlag = FALSE;
    switch (mode) {
    case RF_PRED:
      if (RF_fmRecordSize > 0) {
        dualUseFlag = TRUE;
        recordSize = RF_fmRecordSize;
      }
      break;
    default:
      if (RF_mRecordSize > 0) {
        dualUseFlag = TRUE;
        recordSize = RF_mRecordSize;
      }
      break;
    }  
    if (dualUseFlag == TRUE) {
      free_cmatrix(RF_dmRecordBootFlag, 1, RF_ntree, 1, recordSize);
    }
  }
}
void stackMissingSignatures(uint     obsSize,
                            uint     rspSize,
                            double **responsePtr,
                            double **predictorPtr,
                            uint    *recordMap,
                            uint     recordSize,
                            uint   **p_recordIndex,
                            uint    *p_pIndexSize,
                            int   ***p_pSign,
                            int    **p_pIndex,
                            uint    *pRF_mrFactorSize,
                            uint   **pRF_mrFactorIndex,
                            uint    *pRF_mxFactorSize,
                            uint   **pRF_mxFactorIndex,
                            char    *pRF_mTimeFlag,
                            char    *pRF_mStatusFlag,
                            char    *pRF_mResponseFlag,
                            char    *pRF_mPredictorFlag) {
  uint i, j, p;
  if (recordSize < 1) {
    RF_nativePrint("\nRF-SRC:  *** ERROR *** ");
    RF_nativePrint("\nRF-SRC:  Attempt to allocate for missingness in its absence.");
    RF_nativePrint("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  *p_recordIndex = uivector(1, recordSize);
  i = 0;
  for (j = 1; j <= obsSize; j++) {
    if (recordMap[j] > 0) {
      i++;
      (*p_recordIndex)[i] = j;
    }
  }
  *p_pSign = imatrix(1, rspSize + RF_xSize, 1, recordSize);
  for (j = 1; j <= recordSize; j++) {
    for (i = 1; i <= rspSize + RF_xSize; i++) {
      (*p_pSign)[i][j] = 0;
    }
  }
  for (j = 1; j <= recordSize; j++) {
    for (i = 1; i <= rspSize; i++) {
      if (RF_nativeIsNaN(responsePtr[i][(*p_recordIndex)[j]])) {
        (*p_pSign)[i][j] = 1;
      }
    }
    for (i = 1; i <= RF_xSize; i++) {
      if (RF_nativeIsNaN(predictorPtr[i][(*p_recordIndex)[j]])) {
        (*p_pSign)[rspSize + i][j] = 1;
      }
    }
  }
  *pRF_mStatusFlag = *pRF_mTimeFlag = *pRF_mResponseFlag = *pRF_mPredictorFlag = FALSE;
  *p_pIndex = ivector(1, rspSize + RF_xSize);
  *p_pIndexSize = 0;
  for (i = 1; i <= rspSize; i++) {
    (*p_pIndex)[i] = 0;
    for (j = 1; j <= recordSize; j++) {
      if ((*p_pSign)[i][j] == 1) {
        (*p_pIndexSize) ++;
        (*p_pIndex)[*p_pIndexSize] = - i;
        *pRF_mResponseFlag = TRUE;
        if (i == RF_timeIndex) {
          *pRF_mTimeFlag = TRUE;
        }
        else if (i == RF_statusIndex) {
          *pRF_mStatusFlag = TRUE;
        }
        j = recordSize;
      }
    }
  }  
  for (i = rspSize + 1; i <= rspSize + RF_xSize; i++) {
    (*p_pIndex)[i] = 0;
    for (j = 1; j <= recordSize; j++) {
      if ((*p_pSign)[i][j] == 1) {
        (*p_pIndexSize) ++;
        (*p_pIndex)[*p_pIndexSize] =  i - rspSize;
        *pRF_mPredictorFlag = TRUE;
        j = recordSize;
      }
    }
  }  
  if (rspSize > 0) {
    *pRF_mrFactorIndex = uivector(1, rspSize);
    for (p = 1; p <= rspSize; p++) {
      (*pRF_mrFactorIndex)[p] = 0;
    }
  }
  *pRF_mxFactorIndex = uivector(1, RF_xSize);
  for (p = 1; p <= RF_xSize; p++) {
    (*pRF_mxFactorIndex)[p] = 0;
  }
  *pRF_mrFactorSize = *pRF_mxFactorSize = 0;
  for (p = 1; p <= *p_pIndexSize; p++) {
    if ((*p_pIndex)[p] < 0) {
      if ((RF_rType[(uint) abs((*p_pIndex)[p])] == 'B') ||
          (RF_rType[(uint) abs((*p_pIndex)[p])] == 'I') ||
          (RF_rType[(uint) abs((*p_pIndex)[p])] == 'C')) {
        (*pRF_mrFactorSize) ++;
        (*pRF_mrFactorIndex)[*pRF_mrFactorSize] = (uint) abs((*p_pIndex)[p]);
      }
    }
    else {
      if ((RF_xType[(*p_pIndex)[p]] == 'B') ||
          (RF_xType[(*p_pIndex)[p]] == 'C')) {
        (*pRF_mxFactorSize) ++;
        (*pRF_mxFactorIndex)[*pRF_mxFactorSize] = (*p_pIndex)[p];
      }
    }
  }
}
void unstackMissingSignatures(uint      rspSize,
                              uint      recordSize,
                              uint     *recordIndex,
                              uint      vIndexSize,
                              int     **vSign,
                              int      *vIndex,
                              uint      mrFactorSize,
                              uint     *mrFactorIndex,
                              uint      mxFactorSize,
                              uint     *mxFactorIndex) {
  if (recordSize == 0) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Attempt to deallocate for missingness in its absence.");
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  free_uivector(recordIndex, 1, recordSize);
  free_imatrix(vSign, 1, rspSize + RF_xSize, 1, recordSize);
  free_ivector(vIndex, 1, rspSize + RF_xSize);
  if (rspSize > 0) {
    free_uivector(mrFactorIndex, 1, rspSize);
  }
  free_uivector(mxFactorIndex, 1, RF_xSize);
}
void initializeFactorArrays(char mode) {
  uint j;
  if (RF_rFactorCount + RF_xFactorCount > 0) {
    RF_rMaxFactorLevel = 0;
    for (j = 1; j <= RF_rFactorCount; j++) {
      RF_rFactorSize[j] = RF_rLevelsMax[RF_rFactorIndex[j]];
      if (RF_rMaxFactorLevel < RF_rFactorSize[j]) {
        RF_rMaxFactorLevel = RF_rFactorSize[j];
      }
    }
    RF_xMaxFactorLevel = 0;
    for (j = 1; j <= RF_xFactorCount; j++) {
      RF_xFactorSize[j] = RF_xLevelsMax[RF_xFactorIndex[j]];
      if (RF_xMaxFactorLevel < RF_xFactorSize[j]) {
        RF_xMaxFactorLevel = RF_xFactorSize[j];
      }
    }
    RF_maxFactorLevel = (RF_xMaxFactorLevel > RF_rMaxFactorLevel) ? RF_xMaxFactorLevel : RF_rMaxFactorLevel;
    RF_factorList = (Factor ***) new_vvector(1, RF_ntree, NRUTIL_FPTR2);
    for (j = 1; j <= RF_ntree; j++) {
      RF_factorList[j] = NULL;
    }
  }
}
char stackCompetingArrays(char mode) {
  uint obsSize;
  double  *statusPtr;
  uint    *mRecordMap;
  int    **mpSign;
  uint     mRecordSize;
  char eventSubsetFlag;
  char statusFlag;
  uint *eventCounter;
  uint i, j;
  if (RF_statusIndex == 0) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Attempt to stack competing risk structures in the absence of SURV data.");
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  switch (mode) {
  case RF_GROW:
    if ((RF_splitRule == SURV_CR_LAU) || (RF_splitRule == SURV_CR_GEN)) {
      if (RF_eventTypeSize > 1) {
        RF_opt = RF_opt | OPT_COMP_RISK;
      }
      else {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Parameter verification failed.");
        RF_nativeError("\nRF-SRC:  Competing Risk analysis has been requested.");
        RF_nativeError("\nRF-SRC:  The train data set does not contain competing risks.");
        RF_nativeExit();
      }
    }
    else {
      if (RF_splitRule == CUST_SPLIT) {
        if (RF_eventTypeSize > 1) {
          RF_opt = RF_opt | OPT_COMP_RISK;
        }
        else {
          RF_opt = RF_opt & (~OPT_COMP_RISK);
        }
      }
      else {
        RF_opt = RF_opt & (~OPT_COMP_RISK);
      }
    }
    break;
  default:
    break;
  }
  if (RF_eventTypeSize == 0) {
    if ((RF_opt & OPT_OUTC_TYPE) && !(RF_opt & OPT_PERF) && !(RF_opt & OPT_VIMP)) {
      RF_opt                  = RF_opt & (~OPT_OENS);
      RF_opt                  = RF_opt & (~OPT_FENS);
    }
    else {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Parameter verification failed.");
      RF_nativeError("\nRF-SRC:  Performance or vimp has been requested.");
      RF_nativeError("\nRF-SRC:  The train or pseudo-train data set does not contain any events.");
      RF_nativeExit();
    }
  }
  else {
    hpsortui(RF_eventType, RF_eventTypeSize);
    RF_eventTypeIndex  = uivector(1, RF_eventType[RF_eventTypeSize]);
    for (j = 1; j <= RF_eventType[RF_eventTypeSize]; j++) {
      RF_eventTypeIndex[j] = 0;
    }
    for (j = 1; j <= RF_eventTypeSize; j++) {
      RF_eventTypeIndex[RF_eventType[j]] = j;
    }
  }
  switch (mode) {
  case RF_GROW:
    if (RF_splitRule == RAND_SPLIT) {
      if (RF_eventTypeSize == 1) {
      }
      else {
        RF_opt = RF_opt | OPT_COMP_RISK;
      }
    }
    if (RF_splitRule == SURV_CR_LAU) {
      if (RF_eventTypeSize == 1) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Split rule specified is for Competing Risk scenarios only.");
        RF_nativeError("\nRF-SRC:  The data set does not contain multiple events.");
        RF_nativeExit();
      }
      if(RF_crWeightSize != RF_eventTypeSize) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Parameter verification failed.");
        RF_nativeError("\nRF-SRC:  Competing risk weight vector must be of size equal to number of event types:  %12d != %12d \n", RF_crWeightSize, RF_eventTypeSize);
        RF_nativeExit();
      }
      i = 0;
      for (j = 1; j <= RF_eventTypeSize; j++) {
        if(RF_crWeight[j] < 0) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  Parameter verification failed.");
          RF_nativeError("\nRF-SRC:  Competing risk weight elements must be greater than or equal to zero:  %12.4f \n", RF_crWeight[j]);
          RF_nativeExit();
        }
        else {
          if(RF_crWeight[j] == 0) {
            i ++;
          }
        }
      }
      if (i == RF_eventTypeSize) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Parameter verification failed.");
        RF_nativeError("\nRF-SRC:  Competing risk weight elements are all zero. \n");
        RF_nativeExit();
      }
    }
    break;
  default:
    if (RF_opt & OPT_COMP_RISK) {
      if (RF_eventTypeSize == 1) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  CR analysis has been specified in !GROW mode.");
        RF_nativeError("\nRF-SRC:  However, the GROW data set does not contain multiple events.");
        RF_nativeExit();
      }
    }
    break;
  }
  switch (mode) {
  case RF_PRED:
    if (RF_frSize > 0) {
      getEventInfo(mode);
    }
    else {
      RF_feventTypeSize = RF_mStatusSize = 0;
    }
    break;
  default:
    getEventInfo(mode);
    break;
  } 
  if (RF_eventTypeSize > 1) {
    if (mode == RF_PRED) {
      if (RF_feventTypeSize > 0) {
        eventSubsetFlag = TRUE;
      }
      else {
        eventSubsetFlag = FALSE;
      }
    }
    else if (mode == RF_REST) {
      if (!(RF_opt & OPT_ANON)) {
        eventSubsetFlag = TRUE;
      }
      else {
        eventSubsetFlag = FALSE;
      }
    }
    else {
      eventSubsetFlag = TRUE;        
    }
  }
  else {
    eventSubsetFlag = FALSE;
  }
  if (eventSubsetFlag == TRUE) {
    switch (mode) {
    case RF_PRED:
      obsSize    = RF_fobservationSize;
      statusPtr  = RF_fresponseIn[RF_statusIndex];
      mpSign     = RF_fmpSign;
      mRecordMap = RF_fmRecordMap;
      mRecordSize = RF_fmRecordSize;
      break;
    default:
      obsSize    = RF_observationSize;
      statusPtr  = RF_responseIn[RF_statusIndex];
      mpSign     = RF_mpSign;
      mRecordMap = RF_mRecordMap;
      mRecordSize = RF_mRecordSize;
      break;
    }
    RF_eIndividualSize = uivector(1, RF_eventTypeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      RF_eIndividualSize[j] = 0;
    }
    if (mRecordSize > 0) {
      for (i = 1; i <= obsSize; i++) {
        statusFlag = FALSE;
        if (mRecordMap[i] == 0) {
          statusFlag = TRUE;
        }
        else {
          if (mpSign[RF_statusIndex][mRecordMap[i]] == 0) {
            statusFlag = TRUE;
          }
        }
        if (statusFlag == TRUE) {
          if ((uint) statusPtr[i] > 0) {
            RF_eIndividualSize[RF_eventTypeIndex[(uint) statusPtr[i]]] ++;
          }
          else {
            for (j=1; j <= RF_eventTypeSize; j++) {
              RF_eIndividualSize[j] ++;
            }
          }
        }
      } 
    }
    else {
      for (i = 1; i <= obsSize; i++) {
        if ((uint) statusPtr[i] > 0) {
          RF_eIndividualSize[RF_eventTypeIndex[(uint) statusPtr[i]]] ++;
        }
        else {
          for (j=1; j <= RF_eventTypeSize; j++) {
            RF_eIndividualSize[j] ++;
          }
        }
      } 
    }
    RF_eIndividualIn = (uint **) new_vvector(1, RF_eventTypeSize, NRUTIL_UPTR);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      RF_eIndividualIn[j] = uivector(1, RF_eIndividualSize[j] + RF_mStatusSize + 1);
    }
    eventCounter = uivector(1, RF_eventTypeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      eventCounter[j] = 0;
    }
    if (mRecordSize > 0) {    
      for (i = 1; i <= obsSize; i++) {
        statusFlag = FALSE;
        if (mRecordMap[i] == 0) {
          statusFlag = TRUE;
        }
        else {
          if (mpSign[RF_statusIndex][mRecordMap[i]] == 0) {
            statusFlag = TRUE;
          }
        }
        if (statusFlag == TRUE) {
          if ((uint) statusPtr[i] > 0) {
            j = RF_eventTypeIndex[(uint) statusPtr[i]];
            eventCounter[j] ++;
            RF_eIndividualIn[j][eventCounter[j]] = i;
          }
          else {
            for (j=1; j <= RF_eventTypeSize; j++) {
              eventCounter[j] ++;
              RF_eIndividualIn[j][eventCounter[j]] = i;
            }
          }
        }
      }
    }
    else {
      for (i = 1; i <= obsSize; i++) {
        if ((uint) statusPtr[i] > 0) {
          j = RF_eventTypeIndex[(uint) statusPtr[i]];
          eventCounter[j] ++;
          RF_eIndividualIn[j][eventCounter[j]] = i;
        }
        else {
          for (j=1; j <= RF_eventTypeSize; j++) {
            eventCounter[j] ++;
            RF_eIndividualIn[j][eventCounter[j]] = i;
          }
        }
      }
    }
    free_uivector(eventCounter, 1, RF_eventTypeSize);
  }  
  return TRUE;
}
void getEventInfo(char mode) {
  uint    obsSize;
  double *status;
  uint   *mRecordMap;
  int   **mpSign;
  uint    mRecordSize;
  uint statusFlag;
  uint leadingIndex;
  uint i, j;
  uint jgrow;
  if (RF_statusIndex == 0) {
    RF_nativeError("\nRF-SRC: *** ERROR *** ");
    RF_nativeError("\nRF-SRC: Attempt to stack competing risk structures in the absence of SURV data.");
    RF_nativeError("\nRF-SRC: Please Contact Technical Support.");
    RF_nativeExit();
  }
  if (mode == RF_PRED) {
    obsSize    = RF_fobservationSize;
    status     = RF_fresponseIn[RF_statusIndex];
    mRecordMap = RF_fmRecordMap;
    mpSign     = RF_fmpSign;
    mRecordSize = RF_fmRecordSize;
  }
  else {
    obsSize    = RF_observationSize;
    status     = RF_responseIn[RF_statusIndex];
    mRecordMap = RF_mRecordMap;
    mpSign     = RF_mpSign;
    mRecordSize = RF_mRecordSize;    
  }
  RF_mStatusSize = 0;
  uint *eventTypeLocal = uivector(1, obsSize);
  uint eventTypeSizeLocal = 0;
  if (mRecordSize > 0) {
    for (i = 1; i <= obsSize; i++) {
      eventTypeLocal[i] = 0;
      statusFlag = FALSE;
      if (mRecordMap[i] == 0) {
        statusFlag = TRUE;
      }
      else {
        if (mpSign[RF_statusIndex][mRecordMap[i]] == 0) {
          statusFlag = TRUE;
        }
      }
      if (statusFlag == TRUE) {
        if ((uint) status[i] > 0) {
          eventTypeSizeLocal ++;
          eventTypeLocal[eventTypeSizeLocal] = (uint) status[i];
        } 
        else {
        }
      }
      else {
        RF_mStatusSize ++;
      }
    }  
  }
  else {
    for (i = 1; i <= obsSize; i++) {
      eventTypeLocal[i] = 0;
      if ((uint) status[i] > 0) {
        eventTypeSizeLocal ++;
        eventTypeLocal[eventTypeSizeLocal] = (uint) status[i];
      } 
      else {
      }
    }
  }
  if (mode == RF_PRED) {
    if(eventTypeSizeLocal > 0) {
      hpsortui(eventTypeLocal, eventTypeSizeLocal);
      leadingIndex = 1;
      for (i = 2; i <= eventTypeSizeLocal; i++) {
        if (eventTypeLocal[i] > eventTypeLocal[leadingIndex]) {
          leadingIndex++;
          eventTypeLocal[leadingIndex] = eventTypeLocal[i];
        }
      }
      eventTypeSizeLocal = leadingIndex;
    }
    if (eventTypeSizeLocal > 0) {
      RF_feventTypeSize = eventTypeSizeLocal;
    }
    else {
      RF_feventTypeSize = 0;
    }
    if (RF_feventTypeSize == 0) {
      if (!(RF_opt & OPT_PERF) && !(RF_opt & OPT_VIMP)) {
      }
      else {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Parameter verification failed.");
        RF_nativeError("\nRF-SRC:  Performance or vimp has been requested.");
        RF_nativeError("\nRF-SRC:  The test or pseudo-train data set does not contain any events.");
        RF_nativeExit();
      }
    }
    else {
      char consistencyFlag = TRUE;
      if (RF_eventTypeSize > 1) {
        for (j = 1; j <= RF_feventTypeSize; j++) {
          for (jgrow = 1; jgrow <= RF_eventTypeSize; jgrow++) {
            if (eventTypeLocal[j] != RF_eventType[jgrow]) {
              if (jgrow == RF_eventTypeSize) {
                consistencyFlag = FALSE;
              }
            }
            else {
              jgrow = RF_eventTypeSize;
            }
          }
        }
      }
      if (consistencyFlag == FALSE) {
        RF_nativeError("\nRF-SRC: *** ERROR *** ");
        RF_nativeError("\nRF-SRC: Unknown event type encountered in PRED mode. ");
        RF_nativeError("\nRF-SRC: Please Contact Technical Support.");
        RF_nativeExit();
      }
    }
  }
  free_uivector(eventTypeLocal, 1, obsSize);
}
void unstackCompetingArrays(char mode) {
  char eventSubsetFlag;
  uint j;
    if (RF_statusIndex == 0) {
      RF_nativeError("\nRF-SRC: *** ERROR *** ");
      RF_nativeError("\nRF-SRC: Attempt to unstack competing risk structures in the absence of SURV data.");
      RF_nativeError("\nRF-SRC: Please Contact Technical Support.");
      RF_nativeExit();
    }
    if (RF_eventTypeSize == 0) {
    }
    else {
      free_uivector(RF_eventTypeIndex, 1, RF_eventType[RF_eventTypeSize]);
    }
    if (RF_eventTypeSize > 1) {
      if (mode == RF_PRED) {
        if (RF_feventTypeSize > 0) {
          eventSubsetFlag = TRUE;
        }
        else {
          eventSubsetFlag = FALSE;
        }
      }
      else if (mode == RF_REST) {
        if (!(RF_opt & OPT_ANON)) {
          eventSubsetFlag = TRUE;
        }
        else {
          eventSubsetFlag = FALSE;
        }
      }
      else {
        eventSubsetFlag = TRUE;        
      }
    }
    else {
      eventSubsetFlag = FALSE;
    }
    if (eventSubsetFlag == TRUE) {
      for (j = 1; j <= RF_eventTypeSize; j++) {
        free_uivector(RF_eIndividualIn[j], 1, RF_eIndividualSize[j] + RF_mStatusSize + 1);
      }
      free_new_vvector(RF_eIndividualIn, 1, RF_eventTypeSize, NRUTIL_UPTR);
      free_uivector(RF_eIndividualSize, 1, RF_eventTypeSize);
    }  
}
char stackClassificationArrays(char mode) {
  uint  minorityClassID, minorityClassCnt;
  uint  majorityClassID, majorityClassCnt;
  uint i, j, k;
  if (RF_rFactorCount == 0) {
    RF_nativeError("\nRF-SRC: *** ERROR *** ");
    RF_nativeError("\nRF-SRC: Attempt to stack classification structures in the absence of CLAS data.");
    RF_nativeError("\nRF-SRC: Please Contact Technical Support.");
    RF_nativeExit();
  }
  RF_classLevelSize = uivector(1, RF_rFactorCount);
  RF_classLevel = (uint **) new_vvector(1, RF_rFactorCount, NRUTIL_UPTR);
  RF_rFactorMinorityFlag = cvector(1, RF_rFactorCount);
  
  RF_rLevels = (uint **) new_vvector(1, RF_rFactorCount, NRUTIL_UPTR);  
  for (k = 1; k <= RF_rFactorCount; k++) {
    if (RF_rLevelsCnt[k] > 0) {
      RF_classLevelSize[k] = RF_rLevelsCnt[k];
      RF_rLevels[k] = (uint *) INTEGER(VECTOR_ELT(RF_rLevelsSEXP, k-1));
      RF_rLevels[k] --;
      RF_classLevel[k] = RF_rLevels[k];
    }
    else {
      RF_nativeError("\nRF-SRC: *** ERROR *** ");
      RF_nativeError("\nRF-SRC: Inconsistent zero-level count in factor:  compressed-index = %10d, y-index = %10d", k, RF_rFactorIndex[k]);
      RF_nativeError("\nRF-SRC: Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  
  /*
  RF_rLevels = (uint **) copy2DObject(RF_rLevelsJNIE, NATIVE_TYPE_INTEGER, &RF_nat2DInfoListSize);
  for (k = 1; k <= RF_rFactorCount; k++) {
    if (RF_rLevelsCnt[k] > 0) {
      RF_classLevelSize[k] = RF_rLevelsCnt[k];
      RF_classLevel[k] = RF_rLevels[k];
    }
    else {
      RF_nativeError("\nRF-SRC: *** ERROR *** ");
      RF_nativeError("\nRF-SRC: Inconsistent zero-level count in factor:  compressed-index = %10d, y-index = %10d", k, RF_rFactorIndex[k]);
      RF_nativeError("\nRF-SRC: Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  */
  /*
  RF_rLevels = (uint **) copyXDObject(RF_xLevelsPYTH, &RF_natXDInfoListSize);
  for (k = 1; k <= RF_rFactorCount; k++) {
    if (RF_rLevelsCnt[k] > 0) {
      RF_classLevelSize[k] = RF_rLevelsCnt[k];
      RF_classLevel[k] = RF_rLevels[k];
    }
    else {
      RF_nativeError("\nRF-SRC: *** ERROR *** ");
      RF_nativeError("\nRF-SRC: Inconsistent zero-level count in factor:  compressed-index = %10d, y-index = %10d", k, RF_rFactorIndex[k]);
      RF_nativeError("\nRF-SRC: Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  */
  RF_classLevelIndex = (uint **) new_vvector(1, RF_rFactorCount, NRUTIL_UPTR);
  for (k = 1; k <= RF_rFactorCount; k++) {
    RF_rFactorMinorityFlag[k] = FALSE;
    RF_classLevelIndex[k] = uivector(1, RF_rFactorSize[k]);
    for (j = 1; j <= RF_rFactorSize[k]; j++) {
      RF_classLevelIndex[k][j] = 0;
    }
    for (j = 1; j <= RF_classLevelSize[k]; j++) {
      RF_classLevelIndex[k][RF_classLevel[k][j]] = j;
    }
  }  
  if (RF_opt & OPT_PERF) {
    if (RF_opt & OPT_CLAS_RFQ) {
      RF_rFactorMinority = uivector(1, RF_rFactorCount);
      RF_rFactorMajority = uivector(1, RF_rFactorCount);
      RF_rFactorThreshold = dvector(1, RF_rFactorCount);
      uint totalCount;
      for (j = 1; j <= RF_rFactorCount; j++) {
        uint *levelCount = uivector(1, RF_rFactorSize[j]);
        totalCount = 0;
        for (k = 1; k <= RF_rFactorSize[j]; k++) {
          levelCount[k] = 0;
        }
        for (i = 1; i <= RF_observationSize; i++) {
          if (!RF_nativeIsNaN(RF_responseIn[RF_rFactorIndex[j]][i])) {
            levelCount[(uint) RF_responseIn[RF_rFactorIndex[j]][i]] ++;
            totalCount ++;
          }
        }
        minorityClassCnt = levelCount[1];
        minorityClassID = 1;
        for (k = 1; k <= RF_rFactorSize[j]; k++) {
          if (levelCount[k] > 0) {
            if (levelCount[k] < minorityClassCnt) {
              minorityClassCnt = levelCount[k];
              minorityClassID = k;
            }
          }
        }
        RF_rFactorMinority[j] = minorityClassID;
        majorityClassCnt = levelCount[1];
        majorityClassID = 1;
        for (k = 1; k <= RF_rFactorSize[j]; k++) {
          if (levelCount[k] >= majorityClassCnt) {
            majorityClassCnt = levelCount[k];
            majorityClassID = k;
          }
        }
        RF_rFactorMajority[j] = majorityClassID;
        RF_rFactorThreshold[j] = (double) levelCount[RF_rFactorMinority[j]] / totalCount;
        free_uivector(levelCount, 1, RF_rFactorSize[j]);
      }
    }
    for (j = 1; j <= RF_rFactorCount; j++) {
      if (RF_rFactorSize[j] == 2) {
        RF_rFactorMinorityFlag[j] = TRUE;
      }
    }
  }
  if (mode == RF_PRED) {
    RF_rFactorSizeTest = uivector(1, RF_rFactorCount);
    if (RF_frSize > 0) {
      for (k = 1; k <= RF_rFactorCount; k++) {
        RF_rFactorSizeTest[k] = RF_rFactorSize[k];
        for (i = 1; i <= RF_fobservationSize; i++) {
          if (!RF_nativeIsNaN(RF_fresponseIn[RF_rFactorIndex[k]][i])) {
            if ((uint) RF_fresponseIn[RF_rFactorIndex[k]][i] > RF_rFactorSize[k]) {
              RF_rFactorSizeTest[k] = (uint) RF_fresponseIn[RF_rFactorIndex[k]][i];
            }
          }
        }
      }
    }
  }
  return TRUE;
}
void unstackClassificationArrays(char mode) {
  uint k;
  if (RF_rFactorCount == 0) {
    RF_nativeError("\nRF-SRC: *** ERROR *** ");
    RF_nativeError("\nRF-SRC: Attempt to unstack classification structures in the absence of CLAS data.");
    RF_nativeError("\nRF-SRC: Please Contact Technical Support.");
    RF_nativeExit();
  }
  for (k = 1; k <= RF_rFactorCount; k++) {
    free_uivector(RF_classLevelIndex[k], 1, RF_rFactorSize[k]);
  }
  free_new_vvector(RF_classLevelIndex, 1, RF_rFactorCount, NRUTIL_UPTR);
  free_uivector(RF_classLevelSize, 1, RF_rFactorCount);
  free_new_vvector(RF_classLevel, 1, RF_rFactorCount, NRUTIL_UPTR);
  free_cvector(RF_rFactorMinorityFlag, 1, RF_rFactorCount);
  
  free_new_vvector(RF_rLevels, 1, RF_rFactorCount, NRUTIL_UPTR);  
  
  if (RF_opt & OPT_PERF) {
    if (RF_opt & OPT_CLAS_RFQ) {
      free_dvector(RF_rFactorThreshold, 1, RF_rFactorCount);
      free_uivector(RF_rFactorMinority, 1, RF_rFactorCount);
      free_uivector(RF_rFactorMajority, 1, RF_rFactorCount);
    }
  }
  if (mode == RF_PRED) {
    free_uivector(RF_rFactorSizeTest, 1, RF_rFactorCount);
  }
}
