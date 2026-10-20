
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "stackOutputQQ.h"
#include "stackOutput.h"
#include "nativeUtil.h"
#include "sexpOutgoing.h"
#include "nrutil.h"
void stackTNQualitativeObjectsKnown(char     mode,
                                    uint   **pRF_RMBR_ID_,
                                    uint   **pRF_AMBR_ID_,
                                    uint   **pRF_TN_RCNT_,
                                    uint   **pRF_TN_ACNT_,
                                    uint   **pRF_OOB_SZ_,
                                    uint   **pRF_IBG_SZ_) {
  ulong localSize;
  if (RF_optHigh & OPT_MEMB_OUTG) {
    localSize = (ulong) RF_ntree * RF_bootstrapSize;
    *pRF_RMBR_ID_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_RMBR_ID, localSize, 0, RF_sexpString[RF_RMBR_ID], &RF_RMBR_ID_ptr, 2, RF_ntree, RF_bootstrapSize);
    localSize = (ulong) RF_ntree * RF_observationSize;
    *pRF_AMBR_ID_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_AMBR_ID, localSize, 0, RF_sexpString[RF_AMBR_ID], &RF_AMBR_ID_ptr, 2, RF_ntree, RF_observationSize);
    *pRF_OOB_SZ_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_OOB_SZ, RF_ntree, 0, RF_sexpString[RF_OOB_SZ], NULL, 1, RF_ntree);
    (*pRF_OOB_SZ_) --;
    *pRF_IBG_SZ_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_IBG_SZ, RF_ntree, 0, RF_sexpString[RF_IBG_SZ], NULL, 1, RF_ntree);
    (*pRF_IBG_SZ_) --;
  }
  else if (RF_optHigh & OPT_MEMB_INCG) {
    int *dim = ivector(1, 2);
    dim[1] = RF_ntree;
    dim[2] = RF_bootstrapSize;
    allocateAuxiliaryInfo(FALSE,
                          NATIVE_TYPE_INTEGER,
                          RF_sexpString[RF_RMBR_ID],
                          RF_incAuxiliaryInfoList,
                          RF_incStackCount,
                          *pRF_RMBR_ID_,
                          &RF_RMBR_ID_ptr,
                          2,
                          dim);
    RF_incStackCount ++;
    dim[1] = RF_ntree;
    dim[2] = RF_observationSize;
    allocateAuxiliaryInfo(FALSE,
                          NATIVE_TYPE_INTEGER,
                          RF_sexpString[RF_AMBR_ID],
                          RF_incAuxiliaryInfoList,
                          RF_incStackCount,                          
                          *pRF_AMBR_ID_,
                          &RF_AMBR_ID_ptr,
                          2,
                          dim);
    RF_incStackCount ++;
    dim[1] = RF_ntree;
    dim[2] = -2;
    allocateAuxiliaryInfo(FALSE,
                          NATIVE_TYPE_INTEGER,
                          RF_sexpString[RF_TN_RCNT],
                          RF_incAuxiliaryInfoList,
                          RF_incStackCount,
                          *pRF_TN_RCNT_,
                          &RF_TN_RCNT_ptr,
                          2,
                          dim);
    RF_incStackCount ++;
    allocateAuxiliaryInfo(FALSE,
                          NATIVE_TYPE_INTEGER,
                          RF_sexpString[RF_TN_ACNT],
                          RF_incAuxiliaryInfoList,
                          RF_incStackCount,
                          *pRF_TN_ACNT_,
                          &RF_TN_ACNT_ptr,
                          2,
                          dim);
    RF_incStackCount ++;
    free_ivector(dim, 1, 2);
  }
}
void stackTNQualitativeObjectsUnknown(char     mode,
                                      uint   **pRF_TN_RCNT_,
                                      uint   **pRF_TN_ACNT_,
                                      uint   **pRF_TN_OCNT_,
                                      uint   **pRF_TN_ICNT_) {
  LeafLinkedObj *leafLinkedPtr;
  ulong localSize;
  uint i;
  if (RF_optHigh & OPT_MEMB_OUTG) {
    localSize = RF_totalTerminalCount;
    *pRF_TN_RCNT_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_TN_RCNT, localSize, 0, RF_sexpString[RF_TN_RCNT], &RF_TN_RCNT_ptr, 2, RF_ntree, -2);
    *pRF_TN_ACNT_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_TN_ACNT, localSize, 0, RF_sexpString[RF_TN_ACNT], &RF_TN_ACNT_ptr, 2, RF_ntree, -2);
    *pRF_TN_OCNT_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_TN_OCNT, localSize, 0, RF_sexpString[RF_TN_OCNT], &RF_TN_OCNT_ptr, 2, RF_ntree, -2);
    *pRF_TN_ICNT_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_TN_ICNT, localSize, 0, RF_sexpString[RF_TN_ICNT], &RF_TN_ICNT_ptr, 2, RF_ntree, -2);
    for (i = 1; i <= RF_ntree; i++) {
      leafLinkedPtr = RF_leafLinkedObjHead[i] -> fwdLink;
      while (leafLinkedPtr != NULL) {
        RF_TN_RCNT_ptr[i][(leafLinkedPtr -> termPtr) -> nodeID] = leafLinkedPtr -> ibgMembrCount;
        RF_TN_ACNT_ptr[i][(leafLinkedPtr -> termPtr) -> nodeID] = leafLinkedPtr -> allMembrCount;
        RF_TN_OCNT_ptr[i][(leafLinkedPtr -> termPtr) -> nodeID] = leafLinkedPtr -> oobMembrCount = leafLinkedPtr -> termPtr -> oobMembrSize;
        RF_TN_ICNT_ptr[i][(leafLinkedPtr -> termPtr) -> nodeID] = leafLinkedPtr -> ibgMembrCount = leafLinkedPtr -> termPtr -> ibgMembrSize;        
        leafLinkedPtr = leafLinkedPtr -> fwdLink;
      }
    }
  }
}
void stackTNQuantitativeForestObjectsPtrOnly(char mode) {
  if (RF_optHigh & OPT_TERM_OUTG) {
    RF_totalTerminalCount = 0;
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      RF_TN_MORT_ptr = (double ***) new_vvector(1, RF_ntree, NRUTIL_DPTR2);
      if (!(RF_opt & OPT_COMP_RISK)) {
        RF_TN_SURV_ptr = (double ***) new_vvector(1, RF_ntree, NRUTIL_DPTR2);
        RF_TN_NLSN_ptr = (double ***) new_vvector(1, RF_ntree, NRUTIL_DPTR2);
      }
      else {
        RF_TN_CSHZ_ptr = (double ****) new_vvector(1, RF_ntree, NRUTIL_DPTR3);
        RF_TN_CIFN_ptr = (double ****) new_vvector(1, RF_ntree, NRUTIL_DPTR3);
      }
    }
    else {
      if (RF_rNonFactorCount > 0) {
        RF_TN_REGR_ptr = (double ***) new_vvector(1, RF_ntree, NRUTIL_DPTR2);
      }
      if (RF_rFactorCount > 0) {
        RF_TN_CLAS_ptr = (uint ****) new_vvector(1, RF_ntree, NRUTIL_UPTR3);
      }
    }
  }
  else if (RF_optHigh & OPT_TERM_INCG) {
    int *dim = ivector(1, 4);
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      dim[1] = RF_ntree;
      dim[2] = -2;
      dim[3] = RF_eventTypeSize;
      allocateAuxiliaryInfo(FALSE,
                            NATIVE_TYPE_NUMERIC,
                            RF_sexpString[RF_TN_MORT],
                            RF_incAuxiliaryInfoList,
                            RF_incStackCount,
                            RF_TN_MORT_,
                            &RF_TN_MORT_ptr,
                            3,
                            dim);
      RF_incStackCount ++;
      if (!(RF_opt & OPT_COMP_RISK)) {
        dim[1] = RF_ntree;
        dim[2] = -2;
        dim[3] = RF_sortedTimeInterestSize;
        allocateAuxiliaryInfo(FALSE,
                              NATIVE_TYPE_NUMERIC,
                              RF_sexpString[RF_TN_SURV],
                              RF_incAuxiliaryInfoList,
                              RF_incStackCount,
                              RF_TN_SURV_,
                              &RF_TN_SURV_ptr,
                              3,
                              dim);
        RF_incStackCount ++;
        allocateAuxiliaryInfo(FALSE,
                              NATIVE_TYPE_NUMERIC,
                              RF_sexpString[RF_TN_NLSN],
                              RF_incAuxiliaryInfoList,
                              RF_incStackCount,
                              RF_TN_NLSN_,
                              &RF_TN_NLSN_ptr,
                              3,
                              dim);
        RF_incStackCount ++;
      }
      else {
        dim[1] = RF_ntree;
        dim[2] = -2;
        dim[3] = RF_eventTypeSize;
        dim[4] = RF_sortedTimeInterestSize;
        allocateAuxiliaryInfo(FALSE,
                              NATIVE_TYPE_NUMERIC,
                              RF_sexpString[RF_TN_CSHZ],
                              RF_incAuxiliaryInfoList,
                              RF_incStackCount,
                              RF_TN_CSHZ_,
                              &RF_TN_CSHZ_ptr,
                              4,
                              dim);
        RF_incStackCount ++;
        allocateAuxiliaryInfo(FALSE,
                              NATIVE_TYPE_NUMERIC,
                              RF_sexpString[RF_TN_CIFN],
                              RF_incAuxiliaryInfoList,
                              RF_incStackCount,
                              RF_TN_CIFN_,
                              &RF_TN_CIFN_ptr,
                              4,
                              dim);
        RF_incStackCount ++;
      }
    }
    else {
      if (RF_rNonFactorCount > 0) {
        dim[1] = RF_ntree;
        dim[2] = -2;
        dim[3] = RF_rNonFactorCount;
        allocateAuxiliaryInfo(FALSE,
                              NATIVE_TYPE_NUMERIC,
                              RF_sexpString[RF_TN_REGR],
                              RF_incAuxiliaryInfoList,
                              RF_incStackCount,
                              RF_TN_REGR_,
                              &RF_TN_REGR_ptr,
                              3,
                              dim);
        RF_incStackCount ++;
      }
      if (RF_rFactorCount > 0) {
        dim[1] = RF_ntree;
        dim[2] = -2;
        dim[3] = RF_rFactorCount;
        dim[4] = 0;
        allocateAuxiliaryInfo(FALSE,
                              NATIVE_TYPE_INTEGER,
                              RF_sexpString[RF_TN_CLAS],
                              RF_incAuxiliaryInfoList,
                              RF_incStackCount,
                              RF_TN_CLAS_,
                              &RF_TN_CLAS_ptr,
                              4,
                              dim);
        RF_incStackCount ++;
      }
    }
    free_ivector(dim, 1, 4);
  }
}
void unstackTNQuantitativeForestObjectsPtrOnly(char mode) {
  uint i;
  if (RF_optHigh & OPT_TERM_OUTG) {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      for (i = 1; i <= RF_ntree; i++) {
        unstackTNQuantitativeTreeObjectsPtrOnly(i);
      }
      free_new_vvector(RF_TN_MORT_ptr, 1, RF_ntree, NRUTIL_DPTR2);
      if (!(RF_opt & OPT_COMP_RISK)) {
        free_new_vvector(RF_TN_SURV_ptr, 1, RF_ntree, NRUTIL_DPTR2);
        free_new_vvector(RF_TN_NLSN_ptr, 1, RF_ntree, NRUTIL_DPTR2);          
      }
      else {
        free_new_vvector(RF_TN_CSHZ_ptr, 1, RF_ntree, NRUTIL_DPTR3);
        free_new_vvector(RF_TN_CIFN_ptr, 1, RF_ntree, NRUTIL_DPTR3);
      }
    }
    else {
      if ((RF_rNonFactorCount > 0) || (RF_rFactorCount > 0)) {
        for (i = 1; i <= RF_ntree; i++) {
          unstackTNQuantitativeTreeObjectsPtrOnly(i);
        }
        if (RF_rNonFactorCount > 0) {
          free_new_vvector(RF_TN_REGR_ptr, 1, RF_ntree, NRUTIL_DPTR2);
        }
        if (RF_rFactorCount > 0) {
          free_new_vvector(RF_TN_CLAS_ptr, 1, RF_ntree, NRUTIL_UPTR2);
        }
      }
    }
  }
  else if (RF_optHigh & OPT_TERM_INCG) {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    }
    else {
      if (RF_rNonFactorCount > 0) {
      }
      if (RF_rFactorCount > 0) {
      }
    }
  }
}
void stackTNQuantitativeTreeObjectsPtrOnly(uint treeID) {
  uint i, j;
  if (RF_optHigh & OPT_TERM_OUTG) {
    if (RF_tLeafCount[treeID] > 0 ) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        RF_TN_MORT_ptr[treeID] = (double **) new_vvector(1, RF_tLeafCount[treeID], NRUTIL_DPTR);
        for (i = 1; i <= RF_tLeafCount[treeID]; i++) {
          RF_TN_MORT_ptr[treeID][i] = dvector(1, RF_eventTypeSize);
        }
        if (!(RF_opt & OPT_COMP_RISK)) {
          RF_TN_SURV_ptr[treeID] = (double **) new_vvector(1, RF_tLeafCount[treeID], NRUTIL_DPTR);
          RF_TN_NLSN_ptr[treeID] = (double **) new_vvector(1, RF_tLeafCount[treeID], NRUTIL_DPTR);
          for (i = 1; i <= RF_tLeafCount[treeID]; i++) {
            RF_TN_SURV_ptr[treeID][i] = dvector(1, RF_sortedTimeInterestSize);
            RF_TN_NLSN_ptr[treeID][i] = dvector(1, RF_sortedTimeInterestSize);
          }
        }
        else {
          RF_TN_CSHZ_ptr[treeID] = (double ***) new_vvector(1, RF_tLeafCount[treeID], NRUTIL_DPTR2);
          RF_TN_CIFN_ptr[treeID] = (double ***) new_vvector(1, RF_tLeafCount[treeID], NRUTIL_DPTR2);
          for (i = 1; i <= RF_tLeafCount[treeID]; i++) {
            RF_TN_CSHZ_ptr[treeID][i] = (double **) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR);
            RF_TN_CIFN_ptr[treeID][i] = (double **) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR);
            for (j = 1; j <= RF_eventTypeSize; j++) {
              RF_TN_CSHZ_ptr[treeID][i][j] = dvector(1, RF_sortedTimeInterestSize);
              RF_TN_CIFN_ptr[treeID][i][j] = dvector(1, RF_sortedTimeInterestSize);
            }
          }
        }
      }
      else {
        if (RF_rNonFactorCount > 0) {
          RF_TN_REGR_ptr[treeID] = (double **) new_vvector(1, RF_tLeafCount[treeID], NRUTIL_DPTR);
          for (i = 1; i <= RF_tLeafCount[treeID]; i++) {
            RF_TN_REGR_ptr[treeID][i] = dvector(1, RF_rNonFactorCount);
          }
        }
        if (RF_rFactorCount > 0) {
          RF_TN_CLAS_ptr[treeID] = (uint ***) new_vvector(1, RF_tLeafCount[treeID], NRUTIL_UPTR2);
          for (i = 1; i <= RF_tLeafCount[treeID]; i++) {
            RF_TN_CLAS_ptr[treeID][i] = (uint **) new_vvector(1, RF_rFactorCount, NRUTIL_UPTR);
            for (j = 1; j <= RF_rFactorCount; j++) {
              RF_TN_CLAS_ptr[treeID][i][j] = uivector(1, RF_rFactorSize[j]);
            }
          }
        }
      }
    }
  }
}
void unstackTNQuantitativeTreeObjectsPtrOnly(uint treeID) {
  uint i, j;
  if (RF_optHigh & OPT_TERM_OUTG) {
    if (RF_tLeafCount[treeID] > 0 ) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        for (i = 1; i <= RF_tLeafCount[treeID]; i++) {
          free_dvector(RF_TN_MORT_ptr[treeID][i], 1, RF_eventTypeSize);
        }
        free_new_vvector(RF_TN_MORT_ptr[treeID], 1, RF_tLeafCount[treeID], NRUTIL_DPTR);
        if (!(RF_opt & OPT_COMP_RISK)) {
          for (i = 1; i <= RF_tLeafCount[treeID]; i++) {
            free_dvector(RF_TN_SURV_ptr[treeID][i], 1, RF_sortedTimeInterestSize);
            free_dvector(RF_TN_NLSN_ptr[treeID][i], 1, RF_sortedTimeInterestSize);
          }
          free_new_vvector(RF_TN_SURV_ptr[treeID], 1, RF_tLeafCount[treeID], NRUTIL_DPTR);
          free_new_vvector(RF_TN_NLSN_ptr[treeID], 1, RF_tLeafCount[treeID], NRUTIL_DPTR);
        }
        else {
          for (i = 1; i <= RF_tLeafCount[treeID]; i++) {
            for (j = 1; j <= RF_eventTypeSize; j++) {
              free_dvector(RF_TN_CSHZ_ptr[treeID][i][j], 1, RF_sortedTimeInterestSize);
              free_dvector(RF_TN_CIFN_ptr[treeID][i][j], 1, RF_sortedTimeInterestSize);
            }
            free_new_vvector(RF_TN_CSHZ_ptr[treeID][i], 1, RF_eventTypeSize, NRUTIL_DPTR);
            free_new_vvector(RF_TN_CIFN_ptr[treeID][i], 1, RF_eventTypeSize, NRUTIL_DPTR);
          }
          free_new_vvector(RF_TN_CSHZ_ptr[treeID], 1, RF_tLeafCount[treeID], NRUTIL_DPTR2);
          free_new_vvector(RF_TN_CIFN_ptr[treeID], 1, RF_tLeafCount[treeID], NRUTIL_DPTR2);
        }
      }
      else {
        if (RF_rNonFactorCount > 0) {
          for (i = 1; i <= RF_tLeafCount[treeID]; i++) {
            free_dvector(RF_TN_REGR_ptr[treeID][i], 1, RF_rNonFactorCount);
          }
          free_new_vvector(RF_TN_REGR_ptr[treeID], 1, RF_tLeafCount[treeID], NRUTIL_DPTR);
        }
        if (RF_rFactorCount > 0) {
          for (i = 1; i <= RF_tLeafCount[treeID]; i++) {
            for (j = 1; j <= RF_rFactorCount; j++) {
              free_uivector(RF_TN_CLAS_ptr[treeID][i][j], 1, RF_rFactorSize[j]);
            }
            free_new_vvector(RF_TN_CLAS_ptr[treeID][i], 1, RF_rFactorCount, NRUTIL_UPTR);
          }
          free_new_vvector(RF_TN_CLAS_ptr[treeID], 1, RF_tLeafCount[treeID], NRUTIL_UPTR2);
        }
      }
    }
  }
}
void saveTNQuantitativeTreeObjects(uint treeID) {
  LeafLinkedObj *leafLinkedPtr;
  Terminal *parent;
  uint leaf;
  uint j, k;
  if (RF_optHigh & OPT_TERM_OUTG) {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      leafLinkedPtr = RF_leafLinkedObjHead[treeID] -> fwdLink;
      while (leafLinkedPtr != NULL) {
        parent = leafLinkedPtr -> termPtr;
        leaf = parent -> nodeID;
        for (j = 1; j <= RF_eventTypeSize; j++) {
          RF_TN_MORT_ptr[treeID][leaf][j] = parent -> mortality[j];
        }
        leafLinkedPtr = leafLinkedPtr -> fwdLink;
      }
      if (!(RF_opt & OPT_COMP_RISK)) {
        leafLinkedPtr = RF_leafLinkedObjHead[treeID] -> fwdLink;
        while (leafLinkedPtr != NULL) {
          parent = leafLinkedPtr -> termPtr;
          leaf = parent -> nodeID;
          for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
            RF_TN_SURV_ptr[treeID][leaf][k] = parent -> survival[k];
            RF_TN_NLSN_ptr[treeID][leaf][k] = parent -> nelsonAalen[k];
          }
          leafLinkedPtr = leafLinkedPtr -> fwdLink;
        }
      }
      else {
        leafLinkedPtr = RF_leafLinkedObjHead[treeID] -> fwdLink;
        while (leafLinkedPtr != NULL) {
          parent = leafLinkedPtr -> termPtr;
          leaf = parent -> nodeID;
          for (j = 1; j <= RF_eventTypeSize; j++) {
            for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
              RF_TN_CSHZ_ptr[treeID][leaf][j][k] = parent -> CSH[j][k];
              RF_TN_CIFN_ptr[treeID][leaf][j][k] = parent -> CIF[j][k];
            }
          }
          leafLinkedPtr = leafLinkedPtr -> fwdLink;
        }
      }
    }
    else {
      if (RF_rNonFactorCount > 0) {
        leafLinkedPtr = RF_leafLinkedObjHead[treeID] -> fwdLink;
        while (leafLinkedPtr != NULL) {
          parent = leafLinkedPtr -> termPtr;
          leaf = parent -> nodeID;
          for (j = 1; j <= RF_rNonFactorCount; j++) {
            RF_TN_REGR_ptr[treeID][leaf][j] = (parent -> meanResponse)[j];
          }
          leafLinkedPtr = leafLinkedPtr -> fwdLink;
        }
      }
      if (RF_rFactorCount > 0) {
        leafLinkedPtr = RF_leafLinkedObjHead[treeID] -> fwdLink;
        while (leafLinkedPtr != NULL) {
          parent = leafLinkedPtr -> termPtr;
          leaf = parent -> nodeID;
          for (j = 1; j <= RF_rFactorCount; j++) {
            for (k = 1; k <= RF_rFactorSize[j]; k++) {
              RF_TN_CLAS_ptr[treeID][leaf][j][k] = (parent -> multiClassProb)[j][k];
            }
          }
          leafLinkedPtr = leafLinkedPtr -> fwdLink;
        }
      }
    }
  }
}
void stackTNQuantitativeForestObjectsOutput(char mode) {
  ulong localSize;
  uint tnDimOne, tnDimTwo;
  uint j;
  if (RF_optHigh & OPT_TERM_OUTG) {
    tnDimOne = RF_totalTerminalCount;
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      localSize = (ulong) tnDimOne * RF_eventTypeSize;
      RF_TN_MORT_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_TN_MORT, localSize, RF_nativeNaN, RF_sexpString[RF_TN_MORT], NULL, 1, localSize);
      RF_TN_MORT_ --;
      if (!(RF_opt & OPT_COMP_RISK)) {
        localSize = (ulong) tnDimOne * RF_sortedTimeInterestSize;
        RF_TN_SURV_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_TN_SURV, localSize, RF_nativeNaN, RF_sexpString[RF_TN_SURV], NULL, 1, localSize);
        RF_TN_NLSN_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_TN_NLSN, localSize, RF_nativeNaN, RF_sexpString[RF_TN_NLSN], NULL, 1, localSize);
        RF_TN_SURV_ --;
        RF_TN_NLSN_ --;
      }
      else {
        localSize = (ulong) tnDimOne * RF_eventTypeSize * RF_sortedTimeInterestSize;
        RF_TN_CSHZ_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_TN_CSHZ, localSize, RF_nativeNaN, RF_sexpString[RF_TN_CSHZ], NULL, 1, localSize);
        RF_TN_CIFN_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_TN_CIFN, localSize, RF_nativeNaN, RF_sexpString[RF_TN_CIFN], NULL, 1, localSize);
        RF_TN_CSHZ_ --;
        RF_TN_CIFN_ --;
      }
    }
    else {
      if (RF_rNonFactorCount > 0) {
        localSize = (ulong) tnDimOne * RF_rNonFactorCount;
        RF_TN_REGR_ = (double*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_NUMERIC, RF_TN_REGR, localSize, RF_nativeNaN, RF_sexpString[RF_TN_REGR], NULL, 1, localSize);
        RF_TN_REGR_ --;
      }
      if (RF_rFactorCount > 0) {
        tnDimTwo = 0;
        for (j = 1; j <= RF_rFactorCount; j++) {
          tnDimTwo += RF_rFactorSize[j];
        }
        localSize = (ulong) tnDimOne * tnDimTwo;
        RF_TN_CLAS_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_TN_CLAS, localSize, RF_nativeNaN, RF_sexpString[RF_TN_CLAS], NULL, 1, localSize);
        RF_TN_CLAS_ --;
      }
    }
  }
}
void writeTNQuantitativeForestObjectsOutput(char mode) {
  uint i, j, k, m;
  ulong iter;
  if (RF_optHigh & OPT_TERM_OUTG) {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      iter = 0;
      for (i = 1; i <= RF_ntree; i++) {        
        for (j = 1; j <= RF_tLeafCount[i]; j++) {        
          for (k = 1; k <= RF_eventTypeSize; k++) {
            RF_TN_MORT_[++iter] = RF_TN_MORT_ptr[i][j][k];
          }
        }
      }
      if (!(RF_opt & OPT_COMP_RISK)) {
        iter = 0;          
        for (i = 1; i <= RF_ntree; i++) {
          for (j = 1; j <= RF_tLeafCount[i]; j++) {
            for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
              RF_TN_SURV_[++iter] = RF_TN_SURV_ptr[i][j][k];
              RF_TN_NLSN_[iter]   = RF_TN_NLSN_ptr[i][j][k];
            }
          }
        }
      }
      else {
        iter = 0;          
        for (i = 1; i <= RF_ntree; i++) {
          for (j = 1; j <= RF_tLeafCount[i]; j++) {
            for (k = 1; k <= RF_eventTypeSize; k++) {
              for (m = 1; m <= RF_sortedTimeInterestSize; m++) {
                RF_TN_CSHZ_[++iter] = RF_TN_CSHZ_ptr[i][j][k][m];
                RF_TN_CIFN_[iter]   = RF_TN_CIFN_ptr[i][j][k][m];
              }
            }
          }
        }
      }
    }
    else {
      if (RF_rNonFactorCount > 0) {
        iter = 0;
        for (i = 1; i <= RF_ntree; i++) {
          for (j = 1; j <= RF_tLeafCount[i]; j++) {
            for (k = 1; k <= RF_rNonFactorCount; k++) {
              RF_TN_REGR_[++iter] = RF_TN_REGR_ptr[i][j][k];
            }
          }
        }
      }
      if (RF_rFactorCount > 0) {
        iter = 0;
        for (i = 1; i <= RF_ntree; i++) {
          for (j = 1; j <= RF_tLeafCount[i]; j++) {
            for (k = 1; k <= RF_rFactorCount; k++) {
              for (m = 1; m <= RF_rFactorSize[k]; m++) {
                RF_TN_CLAS_[++iter] = RF_TN_CLAS_ptr[i][j][k][m];
              }
            }
          }
        }
      }
    }
  }
}
void stackTNQualitativeObjectsUnknownMembership(char   mode, uint **pRF_OMBR_ID_, uint **pRF_IMBR_ID_) {
  Terminal *termPtr;
  ulong localSize;
  uint  treeID;
  uint  i, j;
  uint iter;
  if (RF_optHigh & OPT_MEMB_OUTG) {
    localSize = 0;
    for (uint treeID = 1; treeID <= RF_ntree; treeID++) {
      localSize += RF_oobSize[treeID];
    }
    if (localSize > 0) {
      RF_OMBR_ID_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_OMBR_ID, localSize, 0, RF_sexpString[RF_OMBR_ID], NULL, 1, localSize);
      RF_OMBR_ID_ --;
      iter = 0;
      for (treeID = 1; treeID <= RF_ntree; treeID++) {
        for (i = 1; i <= RF_tLeafCount[treeID]; i++) {
          termPtr = RF_tTermList[treeID][i];
          for (j = 1; j <= termPtr -> oobMembrSize; j ++) {
            RF_OMBR_ID_[++iter] = termPtr -> oobMembrIndx[j];
          }
        }
      }
    }
    else {
      localSize = 0;
      for (uint treeID = 1; treeID <= RF_ntree; treeID++) {
        localSize += RF_tLeafCount[treeID];
      }
      RF_OMBR_ID_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_OMBR_ID, localSize, 0, RF_sexpString[RF_OMBR_ID], NULL, 1, localSize);      
    }
    localSize = 0;
    for (uint treeID = 1; treeID <= RF_ntree; treeID++) {
      localSize += RF_ibgSize[treeID];
    }
    RF_IMBR_ID_ = (uint*) stackAndProtect(mode, &RF_nativeIndex, NATIVE_TYPE_INTEGER, RF_IMBR_ID, localSize, 0, RF_sexpString[RF_IMBR_ID], NULL, 1, localSize);
    RF_IMBR_ID_ --;
    iter = 0;
    for (treeID = 1; treeID <= RF_ntree; treeID++) {
      for (i = 1; i <= RF_tLeafCount[treeID]; i++) {
        termPtr = RF_tTermList[treeID][i];
        for (j = 1; j <= termPtr -> ibgMembrSize; j ++) {
          RF_IMBR_ID_[++iter] = termPtr -> ibgMembrIndx[j];
        }
      }
    }
  }  
}
