
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "stackPreDefined.h"
#include "impute.h"
#include "nrutil.h"
#include "error.h"
void stackIncomingResponseArrays(char mode) {
  uint i, j;
  RF_timeIndex = RF_statusIndex = 0;
  RF_masterTime = NULL;
  RF_masterTimeIndexIn = NULL;
  if (RF_ySize > 0) {
    RF_yIndex = uivector(1, RF_ySize);
    RF_yIndexZero = uivector(1, RF_ySize);
    j = 0;
    for (i = 1; i <= RF_ySize; i++) {
      if ((RF_rType[i] != 'B') &&
          (RF_rType[i] != 'R') &&
          (RF_rType[i] != 'I') &&
          (RF_rType[i] != 'C') &&
          (RF_rType[i] != 't') &&
          (RF_rType[i] != 'T') &&
          (RF_rType[i] != 'S')) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Invalid type:  [%10d] = %2c", i, RF_rType[i]);
        RF_nativeError("\nRF-SRC:  Variables must be [B], [R], [I], [C], [t], [T], [S].");
        RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
        RF_nativeExit();
      }
      RF_yIndex[i] = RF_yIndexZero[i] = 0;
      if (RF_rType[i] == 'T') {
        RF_timeIndex = i;
      }
      else if (RF_rType[i] == 'S') {
        RF_statusIndex = i;
      }
      else {
        RF_yIndex[++j] = i;
      }
    }
    if (mode == RF_PRED) {
      if (RF_frSize > 0) {
        if (RF_ySize != RF_frSize) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  train and test outcome/response matrices must be of the same dimension.  ");
          RF_nativeError("\nRF-SRC:  train vs test:  %10d vs %10d  ", RF_ySize, RF_frSize);
          RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
          RF_nativeExit();
        }
      }
      else {
        if ((RF_opt & OPT_PERF) | (RF_opt & OPT_VIMP)) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  test outcome/response matrix must be present when PERF or VIMP is requested.  ");
          RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
          RF_nativeExit();
        }
      }
    }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      RF_ptnCount = 0;
    }
    RF_ySizeProxy = RF_ySize - ((RF_timeIndex == 0) ? 0:1) - ((RF_statusIndex == 0) ? 0:1);
    RF_yIndexZeroSize = 0;
  }
  else {
    RF_rType      = NULL;
    RF_responseIn = NULL;
    RF_ySizeProxy = 0;
    RF_yIndexZeroSize = 0;
  }
  if (RF_opt & OPT_ANON) {
    if (mode != RF_PRED) {
      RF_opt = RF_opt & (~OPT_PERF);
      RF_opt = RF_opt & (~OPT_VIMP);
    }
  }
}
void unstackIncomingResponseArrays(char mode) {
  if (RF_ySize > 0) {
    free_uivector(RF_yIndex, 1, RF_ySize);
    free_uivector(RF_yIndexZero, 1, RF_ySize);
  }
}
void stackIncomingCovariateArrays(char mode) {
  uint i;
  for (i = 1; i <= RF_xSize; i++) {
    if ((RF_xType[i] != 'B') &&
        (RF_xType[i] != 'R') &&
        (RF_xType[i] != 'I') &&
        (RF_xType[i] != 'C')) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Invalid type:  [%10d] = %2c", i, RF_xType[i]);
      RF_nativeError("\nRF-SRC:  Variables must be [B], [R], [I] or [C].");
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
}
void unstackIncomingCovariateArrays(char mode) {
}
void stackIncomingArrays(char mode) {
  stackIncomingResponseArrays(mode);
  stackIncomingCovariateArrays(mode);
  if (mode == RF_GROW) {
    if (RF_nodeSize < 1) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Parameter verification failed.");
      RF_nativeError("\nRF-SRC:  Minimum node size must be greater than zero:  %10d \n", RF_nodeSize);
      RF_nativeExit();
    }
    if (RF_bootstrapSize < 1) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Parameter verification failed.");
      RF_nativeError("\nRF-SRC:  Bootstrap size must be greater than zero:  %12d \n", RF_bootstrapSize);
      RF_nativeExit();
    }
    if ( RF_splitRule > MAXM_SPLIT) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Parameter verification failed.");
      RF_nativeError("\nRF-SRC:  Invalid split rule:  %10d \n", RF_splitRule);
      RF_nativeExit();
    }
    if ((RF_splitRule == USPV_NRM) || (RF_splitRule == USPV_WT_OFF) || (RF_splitRule == USPV_WT_HVY)) {
      if ( RF_xSize < 2) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Parameter verification failed.");
        RF_nativeError("\nRF-SRC:  Number of covariates must be greater than or equal to two (2) with specified split rule:  %10d \n", RF_xSize);
        RF_nativeExit();
      }
      if ( ((int) (RF_xSize - RF_ytry) < 1) || (RF_mtry > RF_xSize) ) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Parameter verification failed.");
        RF_nativeError("\nRF-SRC:  ytry and mtry must be within range:  %10d %10d \n", RF_ytry,  RF_mtry);
        RF_nativeExit();
      }
    }
    else {
      if (RF_ySize == 0) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Parameter verification failed.");
        RF_nativeError("\nRF-SRC:  Number of response variables must be greater than zero:  %10d \n", RF_ySize);
        RF_nativeExit();
      }
      if ( ((RF_mtry < 1) || (RF_mtry > RF_xSize)) ) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Parameter verification failed.");
        RF_nativeError("\nRF-SRC:  Number of random covariate parameters must be greater");
        RF_nativeError("\nRF-SRC:  than zero and less than or equal to the total number of covariates:  %10d \n", RF_mtry);
        RF_nativeExit();
      }
    }
    if ((RF_splitRule != USPV_NRM) && (RF_splitRule != USPV_WT_OFF) && (RF_splitRule != USPV_WT_HVY) && (RF_splitRule != RAND_SPLIT)) {
      if ((RF_timeIndex != 0) && (RF_statusIndex != 0)) {
      }
      else {
        if (RF_ySizeProxy == 0) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  No non-[S] and non-[C] responses found.");
          RF_nativeExit();
        }
        if (RF_ytry > RF_ySizeProxy) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  Parameter verification failed.");
          RF_nativeError("\nRF-SRC:  ytry must be within range:  %10d \n", RF_ytry);
          RF_nativeExit();
        }
      }
    }
    for (uint i = 1; i <= RF_xSize; i++) {
      if(RF_xWeightStat[i] < 0) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Parameter verification failed.");
        RF_nativeError("\nRF-SRC:  Split statistical weight elements must be greater than or equal to zero:  %12.4f \n", RF_xWeightStat[i]);
        RF_nativeExit();
      }
    }
    if(RF_ySize > 0) {
      for (uint i = 1; i <= RF_ySize; i++) {
        if(RF_yWeight[i] < 0) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  Parameter verification failed.");
          RF_nativeError("\nRF-SRC:  Y-weight elements must be greater than or equal to zero:  %12.4f \n", RF_yWeight[i]);
          RF_nativeExit();
        }
      }
    }
    for (uint i = 1; i <= RF_xSize; i++) {
      if(RF_xWeight[i] < 0) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Parameter verification failed.");
        RF_nativeError("\nRF-SRC:  X-weight elements must be greater than or equal to zero:  %12.4f \n", RF_xWeight[i]);
        RF_nativeExit();
      }
    }
    if ((RF_timeIndex == 0) && (RF_statusIndex == 0)) {
      if ((RF_splitRule != RAND_SPLIT)  &&
          (RF_splitRule != REGR_NRM)    &&
          (RF_splitRule != REGR_WT_OFF) &&
          (RF_splitRule != REGR_WT_HVY) &&
          (RF_splitRule != REGR_SGS)    &&
          (RF_splitRule != REGR_QUANT)  &&
          (RF_splitRule != LARG_QUANT)  &&
          (RF_splitRule != CLAS_NRM)    &&
          (RF_splitRule != CLAS_WT_OFF) &&
          (RF_splitRule != CLAS_WT_HVY) &&
          (RF_splitRule != CLAS_SGS)    &&
          (RF_splitRule != CLAS_AU_ROC) &&
          (RF_splitRule != CLAS_ENTROP) &&
          (RF_splitRule != MV_NRM)      &&
          (RF_splitRule != MV_WT_OFF)   &&
          (RF_splitRule != MV_WT_HVY)   &&
          (RF_splitRule != USPV_NRM)    &&
          (RF_splitRule != USPV_WT_OFF) &&
          (RF_splitRule != USPV_WT_HVY) &&
          (RF_splitRule != CUST_SPLIT)  &&
          (RF_splitRule != MAHALANOBIS)) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  !SURV data and split rule specified are incompatible.");
        RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
        RF_nativeExit();
      }
      if  ((RF_splitRule == REGR_QUANT) || (RF_splitRule == LARG_QUANT)) {
        if (RF_quantileSize > 0) {
        }
        else {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  Quantile regression split rules require the presence of a probability vector.");
          RF_nativeExit();
        }
      }
      if  (RF_splitRule == MAHALANOBIS) {
      }
    }
    else if ((RF_timeIndex != 0) && (RF_statusIndex != 0)) {
      if ((RF_splitRule != SURV_LGRNK)  &&
          (RF_splitRule != SURV_LRSCR)  &&
          (RF_splitRule != SURV_BSG1)   &&
          (RF_splitRule != SURV_CR_LAU) &&
          (RF_splitRule != SURV_CR_GEN) &&
          (RF_splitRule != RAND_SPLIT)  &&
          (RF_splitRule != CUST_SPLIT)) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  SURV data and split rule specified are incompatible.");
        RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
        RF_nativeExit();
      }
    }
    else {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Data set contains mixed outcomes with no comatible split rule.");
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  if (RF_quantileSize > 0) {
    for (uint i = 1; i <= RF_quantileSize; i++) {
      if ((0 < RF_quantile[i]) && (RF_quantile[i] <= 1.0)) {
      }
      else {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Parameter verification failed.");
        RF_nativeError("\nRF-SRC:  Quantile value is out of range (0, 1):  %.10e ", RF_quantile[i]);
        RF_nativeExit();
      }
    }
  }
}
void unstackIncomingArrays(char mode) {
  unstackIncomingResponseArrays(mode);
  unstackIncomingCovariateArrays(mode);
}
void stackPreDefinedCommonArrays(char          mode,
                                 Node      ****nodeMembership,
                                 Terminal  ****tTermMembership,
                                 Terminal  ****tTermList,
                                 Node       ***root) {
  uint i, j, k;
  *nodeMembership = (Node ***)     new_vvector(1, RF_ntree, NRUTIL_NPTR2);
  *tTermMembership = (Terminal ***) new_vvector(1, RF_ntree, NRUTIL_TPTR2);
  *tTermList = (Terminal ***) new_vvector(1, RF_ntree, NRUTIL_NPTR2);
  RF_nodeCount = uivector(1, RF_ntree);
  for (i = 1; i <= RF_ntree; i++) {
    RF_nodeCount[i] = 0;    
  }
  RF_leafLinkedObjHead = (LeafLinkedObj **) new_vvector(1, RF_ntree, NRUTIL_LEAFPTR);
  RF_leafLinkedObjTail = (LeafLinkedObj **) new_vvector(1, RF_ntree, NRUTIL_LEAFPTR);
  RF_bootMembershipIndex = (uint **) new_vvector(1, RF_ntree, NRUTIL_UPTR);
  if ( (RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2) ) {
    for (i = 1; i <= RF_ntree; i++) {
      k = 0;
      for (j = 1; j <= RF_subjSize; j++) {
        k += RF_bootstrapIn[i][j];
      }
      if(k != RF_bootstrapSize) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Parameter verification failed.");
        RF_nativeError("\nRF-SRC:  Bootstrap size implied by samp matrix inconsistent:  %12d found vs. %12d specified \n", k, RF_bootstrapSize);
        RF_nativeExit();
      }
    }
  }
  RF_bootMembershipFlag = (char **) new_vvector(1, RF_ntree, NRUTIL_CPTR);
  RF_bootMembershipCount = (uint **) new_vvector(1, RF_ntree, NRUTIL_UPTR);
  RF_oobMembershipFlag = (char **) new_vvector(1, RF_ntree, NRUTIL_CPTR);
  RF_ibgMembershipIndex = (uint **) new_vvector(1, RF_ntree, NRUTIL_UPTR);
  RF_oobMembershipIndex = (uint **) new_vvector(1, RF_ntree, NRUTIL_UPTR);
  RF_oobSize = uivector(1, RF_ntree);
  RF_ibgSize = uivector(1, RF_ntree);
  RF_maxDepth = uivector(1, RF_ntree);
  RF_orderedTreeIndex = uivector(1, RF_ntree);
  for (i = 1; i <= RF_ntree; i++) {
    RF_orderedTreeIndex[i] = i;
  }
  RF_serialTreeIndex = uivector(1, RF_ntree);
  *root = (Node **) new_vvector(1, RF_ntree, NRUTIL_NPTR);
  for (i = 1; i <= RF_ntree; i++) {
    (*root)[i] = NULL;
  }
  if (RF_ptnCount > 0) {
    RF_pNodeMembership = (Node ***)     new_vvector(1, RF_ntree, NRUTIL_NPTR2);
    RF_pTermMembership = (Terminal ***) new_vvector(1, RF_ntree, NRUTIL_NPTR2);
    RF_pNodeList = (Node ***)     new_vvector(1, RF_ntree, NRUTIL_NPTR2);
    RF_pTermList = (Terminal ***) new_vvector(1, RF_ntree, NRUTIL_NPTR2);
    RF_pLeafCount = uivector(1, RF_ntree);
  }
  if ( ( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
       (!(RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ) {
    for (i = 1; i <= RF_subjSize; i++) {
      if(RF_subjWeight[i] < 0) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Parameter verification failed.");
        RF_nativeError("\nRF-SRC:  Subject-weight elements must be greater than or equal to zero:  %12.4f \n", RF_subjWeight[i]);
        RF_nativeExit();
      }
    }
    stackWeights(RF_subjWeight,
                 RF_subjSize,
                 &RF_subjWeightType,
                 &RF_subjWeightSorted,
                 &RF_subjWeightDensitySize); 
  }
  RF_getTreeIndex = uivector(1, RF_ntree);
  if (mode == RF_GROW) {
    for (i = 1; i <= RF_ntree; i++) {
      RF_getTreeIndex[i] = i;
    }
    RF_getTreeCount = RF_ntree;
  }
  else {
    RF_getTreeCount = 0;    
    for (i = 1; i <= RF_ntree; i++) {
      if (RF_getTree[i] != 0) {
        RF_getTreeIndex[++RF_getTreeCount] = i;
      }
    }
  }
}
void unstackPreDefinedCommonArrays(char         mode,
                                   Node      ***nodeMembership,
                                   Terminal  ***tTermMembership,
                                   Terminal  ***tTermList,
                                   Node       **root) {
  free_new_vvector(nodeMembership, 1, RF_ntree, NRUTIL_NPTR2);
  free_new_vvector(tTermMembership, 1, RF_ntree, NRUTIL_TPTR2);
  free_new_vvector(tTermList, 1, RF_ntree, NRUTIL_TPTR2);
  free_new_vvector(RF_leafLinkedObjHead, 1, RF_ntree, NRUTIL_LEAFPTR);
  free_new_vvector(RF_leafLinkedObjTail, 1, RF_ntree, NRUTIL_LEAFPTR);
  free_uivector(RF_nodeCount, 1, RF_ntree);
  free_new_vvector(RF_bootMembershipIndex, 1, RF_ntree, NRUTIL_UPTR);
  if ( (RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2) ) {
  }
  free_new_vvector(RF_bootMembershipFlag, 1, RF_ntree, NRUTIL_CPTR);
  free_new_vvector(RF_bootMembershipCount, 1, RF_ntree, NRUTIL_UPTR);
  free_new_vvector(RF_oobMembershipFlag, 1, RF_ntree, NRUTIL_CPTR);
  free_new_vvector(RF_ibgMembershipIndex, 1, RF_ntree, NRUTIL_UPTR);
  free_new_vvector(RF_oobMembershipIndex, 1, RF_ntree, NRUTIL_UPTR);
  free_uivector(RF_oobSize, 1, RF_ntree);
  free_uivector(RF_ibgSize, 1, RF_ntree);
  free_uivector(RF_maxDepth, 1, RF_ntree);
  free_uivector(RF_orderedTreeIndex, 1, RF_ntree);
  free_uivector(RF_serialTreeIndex, 1, RF_ntree);
  free_new_vvector(root, 1, RF_ntree, NRUTIL_NPTR);
  if (RF_ptnCount > 0) {
    free_new_vvector(RF_pNodeMembership, 1, RF_ntree, NRUTIL_NPTR2);
    free_new_vvector(RF_pTermMembership, 1, RF_ntree, NRUTIL_NPTR2);
    free_new_vvector(RF_pNodeList, 1, RF_ntree, NRUTIL_NPTR2);
    free_new_vvector(RF_pTermList, 1, RF_ntree, NRUTIL_NPTR2);
    free_uivector(RF_pLeafCount, 1, RF_ntree);
  }
  if ( (!(RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
       ( (RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ) {
    unstackWeights(RF_subjWeightType, RF_subjSize, RF_subjWeightSorted); 
  }
  free_uivector(RF_getTreeIndex, 1, RF_ntree);
}
void stackPreDefinedGrowthArrays(void) {
  uint i;
  if (RF_opt & OPT_VIMP) {
    RF_intrPredictorSize = RF_xSize;
    RF_intrPredictor = uivector(1, RF_intrPredictorSize);
    for (i = 1; i <= RF_intrPredictorSize; i++) {
      RF_intrPredictor[i] = i;
    }
    RF_importanceFlag = cvector(1, RF_xSize);
    for (i = 1; i <= RF_xSize; i++) {
      RF_importanceFlag[i] = TRUE;
    }
  }
  else {
    RF_intrPredictorSize = 0;    
  }
  stackWeights(RF_xWeight,
               RF_xSize,
               &RF_xWeightType,
               &RF_xWeightSorted,
               &RF_xWeightDensitySize); 
  if(RF_ySize > 0) {
    stackWeights(RF_yWeight,
                 RF_ySize,
                 &RF_yWeightType,
                 &RF_yWeightSorted,
                 &RF_yWeightDensitySize); 
    RF_yIndexZeroSize = 0;
    for (i = 1; i <= RF_ySizeProxy; i++) {
      if (RF_yWeight[RF_yIndex[i]] == 0) {
        RF_yIndexZero[++RF_yIndexZeroSize] = RF_yIndex[i];
      }
    }
  }
}
void unstackPreDefinedGrowthArrays(void) {
  if (RF_opt & OPT_VIMP) {
    free_uivector(RF_intrPredictor, 1, RF_intrPredictorSize);
    free_cvector(RF_importanceFlag, 1, RF_xSize);
  }
  unstackWeights(RF_xWeightType,
                 RF_xSize,
                 RF_xWeightSorted); 
  if(RF_ySize > 0) {
    unstackWeights(RF_yWeightType,
                   RF_ySize,
                   RF_yWeightSorted); 
  }
}
void stackPreDefinedRestoreArrays(void) {
  uint i;
  if (RF_opt & OPT_VIMP) {
    checkInteraction();
    RF_importanceFlag = cvector(1, RF_xSize);
    for (i = 1; i <= RF_xSize; i++) {
      RF_importanceFlag[i] = FALSE;
    }
    for (i = 1; i <= RF_intrPredictorSize; i++) {
      RF_importanceFlag[RF_intrPredictor[i]] = TRUE;
    }
  }
}
void unstackPreDefinedRestoreArrays(void) {
  if (RF_opt & OPT_VIMP) {
    free_cvector(RF_importanceFlag, 1, RF_xSize);
  }
}
void stackPreDefinedPredictArrays(void) {
  uint i;
  RF_fnodeMembership = (Node ***)     new_vvector(1, RF_ntree, NRUTIL_NPTR2);
  RF_ftTermMembership = (Terminal ***) new_vvector(1, RF_ntree, NRUTIL_TPTR2);
  RF_fidentityMembershipIndex = uivector(1, RF_fobservationSize);
  for (i = 1; i <= RF_fobservationSize; i++) {
    RF_fidentityMembershipIndex[i] = i;
  }
  RF_testMembershipFlag = cvector(1, RF_fobservationSize);
  for (i = 1; i <= RF_fobservationSize; i++) {
    RF_testMembershipFlag[i] = ACTIVE;
  }
  if (RF_opt & OPT_VIMP) {
    checkInteraction();
    RF_importanceFlag = cvector(1, RF_xSize);
    for (i = 1; i <= RF_xSize; i++) {
      RF_importanceFlag[i] = FALSE;
    }
    for (i = 1; i <= RF_intrPredictorSize; i++) {
      RF_importanceFlag[RF_intrPredictor[i]] = TRUE;
    }
  }
}
void unstackPreDefinedPredictArrays(void) {
  free_new_vvector(RF_fnodeMembership, 1, RF_ntree, NRUTIL_NPTR2);
  free_new_vvector(RF_ftTermMembership, 1, RF_ntree, NRUTIL_TPTR2);
  free_uivector(RF_fidentityMembershipIndex, 1, RF_fobservationSize);
  free_cvector(RF_testMembershipFlag, 1, RF_fobservationSize);
  if (RF_opt & OPT_VIMP) {
    free_cvector(RF_importanceFlag, 1, RF_xSize);
  }
}
void checkInteraction(void) {
  uint leadingIndex, i;
  if((RF_intrPredictorSize <= 0) || (RF_intrPredictorSize > RF_xSize)) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Parameter verification failed.");
    RF_nativeError("\nRF-SRC:  Number of predictors to be perturbed must be greater than zero and less than or equal to %10d:  %10d \n", RF_xSize, RF_intrPredictorSize);
    RF_nativeExit();
  }
  uint *intrPredictorCopy = uivector(1, RF_intrPredictorSize);
  for (i=1; i <= RF_intrPredictorSize; i++) {
    intrPredictorCopy[i] = RF_intrPredictor[i];
  }
  hpsortui(intrPredictorCopy, RF_intrPredictorSize);
  leadingIndex = 1;
  for (i=2; i <= RF_intrPredictorSize; i++) {
    if (intrPredictorCopy[i] > intrPredictorCopy[leadingIndex]) {
      leadingIndex++;
    }
  }
  if (RF_intrPredictorSize != leadingIndex) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Parameter verification failed.");
    RF_nativeError("\nRF-SRC:  Interaction terms are not unique.");
    RF_nativeError("\nRF-SRC:  Only %10d of %10d are unique.", leadingIndex, RF_intrPredictorSize);
    RF_nativeExit();
  }
  free_uivector(intrPredictorCopy, 1, RF_intrPredictorSize);
  for (i=1; i <= RF_intrPredictorSize; i++) {
    if (RF_intrPredictor[i] > RF_xSize) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Parameter verification failed.");
      RF_nativeError("\nRF-SRC:  Interaction terms are not coherent.");
      RF_nativeError("\nRF-SRC:  Predictor encountered is %10d, maximum allowable is %10d.", RF_intrPredictor[i], RF_xSize);
      RF_nativeExit();
    }
  }
}
void stackWeights(double *weight,
                  uint    size,
                  uint   *weightType,
                  uint  **weightSorted,
                  uint   *weightDensitySize) {
  char uniformFlag, integerFlag;
  double meanWeight;
  uint i;
  *weightSorted      = NULL;
  *weightDensitySize = 0;
  meanWeight = getMeanValue(weight, size);
  uniformFlag = TRUE;
  i = 0;
  while (uniformFlag && (i < size)) {
    ++i;
    if (fabs(weight[i] - meanWeight) > 0) {
      uniformFlag = FALSE;
    }
  }
  if (uniformFlag) {
    *weightType = RF_WGHT_UNIFORM;
  } 
  else {
    integerFlag = TRUE;
    i = 0;
    while (integerFlag && (i < size)) {
      i++;
      if (fabs(round(weight[i]) - weight[i]) > 0.0) {
        integerFlag = FALSE;
      }
    }
    if(integerFlag) {
      *weightType = RF_WGHT_INTEGER;
    }
    else {
      *weightType = RF_WGHT_GENERIC;
    }
  }
  switch (*weightType) {
  case RF_WGHT_UNIFORM:
    break;
  case RF_WGHT_INTEGER:
    *weightSorted = uivector(1, size);
    indexx(size, weight, *weightSorted);
    *weightDensitySize = 0;
    for (i = 1; i <= size; i++) {
      *weightDensitySize += (uint) weight[i];
    }
    break;
  case RF_WGHT_GENERIC:
    *weightSorted = uivector(1, size);
    indexx(size, weight, *weightSorted);
    break;
  }
}
void unstackWeights(uint    weightType,
                    uint    size,
                    uint   *weightSorted) {
  switch (weightType) {
  case RF_WGHT_UNIFORM:
    break;
  case RF_WGHT_INTEGER:
    free_uivector(weightSorted, 1, size);
    break;
  case RF_WGHT_GENERIC:
    free_uivector(weightSorted, 1, size);
    break;
  }
}
