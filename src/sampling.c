
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "sampling.h"
#include "nrutil.h"
${trace.token} #include "error.h"
DistributionObj *makeDistributionObjRaw(void) {
  DistributionObj *obj = (DistributionObj*) gblock((size_t) sizeof(DistributionObj));
  return obj;
}
DistributionObj *makeDistributionObjFull(void) {
  DistributionObj *obj = (DistributionObj*) gblock((size_t) sizeof(DistributionObj));
  obj -> permissibleIndex  = NULL;
  obj -> permissible       = NULL;
  obj -> permissibleSize   = 0;
  obj -> augmentationSize    = NULL;
  obj -> weightType       = 0;
  obj -> weight           = NULL;
  obj -> weightSorted     = NULL;
  obj -> cdf     = NULL;
  obj -> cdfSize = 0;
  obj -> cdfSort = NULL;
  obj -> density          = NULL;
  obj -> densityAllocSize = 0;
  obj -> densitySize      = 0;
  obj -> densitySwap      = NULL;
  obj -> index           = NULL;
  obj -> indexSize       = 0;
  obj -> uIndexAllocSize = 0;
  obj -> slot            = 0;
  return obj;
}
void freeDistributionObjRaw(DistributionObj *obj) {
  free_gblock(obj, (size_t) sizeof(DistributionObj));
}
void initializeCDFNew(uint treeID, DistributionObj *obj) {
  char validElement;
  uint i, j, k, kk;
  ${trace.token}  if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\ninitializeCDFNew() ENTRY ...\n");
  ${trace.token}  }
  switch (obj -> weightType) {
  case RF_WGHT_UNIFORM:
    if (obj -> permissible != NULL) {
      if (obj -> augmentationSize != NULL) {
        obj -> uIndexAllocSize = obj -> permissibleSize + 
          obj -> augmentationSize[1] +
          obj -> augmentationSize[2] +
          (RF_xSize * (obj -> augmentationSize[2])) +
          ((obj -> augmentationSize[1]) * (obj -> augmentationSize[2]));
      }
      else {
        obj -> uIndexAllocSize = obj -> permissibleSize;
      }
      obj -> index = uivector(1, obj -> uIndexAllocSize);
      obj -> indexSize = 0;
      for (k = 1; k <= obj -> permissibleSize; k++) {
        if (obj -> permissible[k]) {
          obj -> index[++(obj -> indexSize)] = k;
        }
      }
    }
    else {
      obj -> index = uivector(1, obj -> permissibleSize);
      obj -> indexSize = obj -> uIndexAllocSize = obj -> permissibleSize;
      for (k=1; k <= obj -> permissibleSize; k++) {
        obj -> index[k] = obj -> permissibleIndex[k];
      }
    }
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
    ${trace.token}        RF_nativePrint("\nUniform of Size: %10d", obj -> uIndexAllocSize);
    ${trace.token}      }
    ${trace.token}      }
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
    ${trace.token}        RF_nativePrint("\nVector Weights:  Uniform");
    ${trace.token}        RF_nativePrint("\n     index      absolute");
    ${trace.token}        for (k=1; k <= obj -> uIndexAllocSize; k++) {
    ${trace.token}          RF_nativePrint("\n%10d  %10d", k, obj -> index[k]);
    ${trace.token}        }
    ${trace.token}        RF_nativePrint("\n");
    ${trace.token}      }
    ${trace.token}      }
    break;
  case RF_WGHT_INTEGER:
    obj -> density = uivector(1, obj -> densityAllocSize);
    obj -> densitySize = 0;
    obj -> densitySwap = (uint **) new_vvector(1, obj -> permissibleSize, NRUTIL_UPTR);
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}        RF_nativePrint("\nDensity Swap Vector:  ");
    ${trace.token}      }
    for (k = obj -> permissibleSize; k >= 1; k--) {
      kk = obj -> weightSorted[k];
      validElement = TRUE;
      if (obj -> permissible != NULL) {
        if (obj -> permissible[kk] == FALSE) {
          validElement = FALSE;
        }
      }
      if (validElement) {
        j = (uint) (obj -> weight)[kk];
        if (j > 0) {
          (obj -> densitySwap)[kk] = uivector(1, j);
          for (i = 1; i <= j; i++) {
            (obj -> density)[++(obj -> densitySize)] = kk;
            (obj -> densitySwap)[kk][i] = obj -> densitySize;
          }
          ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
          ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
          ${trace.token}        RF_nativePrint("\nSwap Vector for (Weight Index, of Length):  (%10d, %10d)", kk, j);
          ${trace.token}        RF_nativePrint("\n     index          swap");
          ${trace.token}        for (i=1; i <= j; i++) {
          ${trace.token}          RF_nativePrint("\n%10d  %12d", i, (obj -> densitySwap)[kk][i]);
          ${trace.token}        }
          ${trace.token}        RF_nativePrint("\n");
          ${trace.token}      }
          ${trace.token}      }
        }
        else {
          (obj -> densitySwap)[kk] = NULL;
          ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
          ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
          ${trace.token}        RF_nativePrint("\nSwap Vector for (Weight Index, of Length):  (%10d, %10d)", kk, 0);
          ${trace.token}      }
          ${trace.token}      }
        }
      }
      else {
        (obj -> densitySwap)[kk] = NULL;
          ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
          ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
          ${trace.token}        RF_nativePrint("\nSwap Vector for (Weight Index, of Length):  (%10d, %10d)", kk, 0);
          ${trace.token}      }
          ${trace.token}      }
      }
    }
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}        RF_nativePrint("\nVector Weights Integer Size: %10d", obj -> densitySize);
    ${trace.token}      }
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
    ${trace.token}        RF_nativePrint("\nVector Weights:  Integer");
    ${trace.token}        RF_nativePrint("\n     index          CWDV");
    ${trace.token}        for (i=1; i <= (obj -> densitySize); i++) {
    ${trace.token}          RF_nativePrint("\n%10d  %12d", i, (obj -> density)[i]);
    ${trace.token}        }
    ${trace.token}        RF_nativePrint("\n");
    ${trace.token}      }
    ${trace.token}      }
    break;
  case RF_WGHT_GENERIC:
    obj -> index = uivector(1, obj -> permissibleSize);
    obj -> cdf     =  dvector(1, obj -> permissibleSize);
    obj -> cdfSize = 0;
    i = 0;
    for (k = 1; k <= obj -> permissibleSize; k++) {
      kk = obj -> weightSorted[k];
      validElement = TRUE;
      if (obj -> permissible != NULL) {
        if (obj -> permissible[kk] == FALSE) {
          validElement = FALSE;
        }
      }
      if (validElement) {
        if (obj -> weight[kk] > 0) {
          (obj -> index)[++i] = kk;
          (obj -> cdfSize) ++;
          (obj -> cdf)[obj -> cdfSize] = obj -> weight[kk];
        }
      }
    }
    for (k = 2; k <= obj -> cdfSize; k++) {
      (obj -> cdf)[k] += (obj -> cdf)[k-1];
    }
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}        RF_nativePrint("\nIncoming Weights Generic Size:  %10d", obj -> cdfSize);
    ${trace.token}      }
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
    ${trace.token}        RF_nativePrint("\nIncoming Weights:  Generic");
    ${trace.token}        RF_nativePrint("\n  (note that sorting is ignored)");
    ${trace.token}        RF_nativePrint("\n       slot      index                            CDF");
    ${trace.token}        for (k=1; k <= obj -> cdfSize; k++) {
    ${trace.token}          RF_nativePrint("\n %10d %10d %30.24e", k, obj -> index[k], obj -> cdf[k]);
    ${trace.token}        }
    ${trace.token}        RF_nativePrint("\n");
    ${trace.token}      }
    ${trace.token}      }
    break;
  }
  ${trace.token}  if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\ninitializeCDFNew() EXIT ...\n");
  ${trace.token}  }
}
uint sampleFromCDFNew (float (*genericGenerator) (uint), uint treeID, DistributionObj *obj) {
  double randomValue;
  double midValue;
  char flag;
  uint low, mid, high, value;
  uint p;
  ${trace.token}  if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nsampleFromCDFNew() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token} p = 0;
  value = 0;  
  switch (obj -> weightType) {
  case RF_WGHT_UNIFORM:
    if (obj -> indexSize > 0) {
      obj -> slot = (uint) ceil(genericGenerator(treeID) * (obj -> indexSize * 1.0));
      value = obj -> index[obj -> slot];
    }
    else {
      value = obj -> slot = 0;
    }
    break;
  case RF_WGHT_INTEGER:
    if (obj -> densitySize > 0) {
      p = (uint) ceil(genericGenerator(treeID) * (obj -> densitySize * 1.0));
      value = obj -> slot = obj -> density[p];
    }
    else {
      value = obj -> slot = 0;
    }
    break;
  case RF_WGHT_GENERIC:
    if (obj -> cdf[obj -> cdfSize] > 0) {
      randomValue = genericGenerator(treeID) * (obj -> cdf)[obj -> cdfSize];
      ${trace.token}  if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
      ${trace.token}    RF_nativePrint("\nRandom Value:  %30.24e", randomValue);
      ${trace.token}  }
      low  = mid = 1;
      high = obj -> cdfSize;
      while (low < high) {
        mid  = (low + high) >> 1;
        ${trace.token}  if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
        ${trace.token}    RF_nativePrint("\nBinary Search (low, mid, high) = (%10d, %10d, %10d)", low, mid, high);
        ${trace.token}  }
        if (randomValue > obj -> cdf[mid]) {
          ${trace.token}  if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
          ${trace.token}    RF_nativePrint("\nBinary Search:  Upper with (random, cdf[mid] = (%30.24e %30.24e)", randomValue, obj -> cdf[mid]);
          ${trace.token}  }
          if (low == mid) {
            low = high;
            mid = high;
            if (obj -> cdf[mid] == 0) {
              mid ++;
              ${trace.token}  if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
              ${trace.token}    RF_nativePrint("\nZero level CDF adjustment:  %10d", mid);
              ${trace.token}  }
            }
          }
          else {
            low = mid;
          }
        }
        else {
          ${trace.token}  if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
          ${trace.token}    RF_nativePrint("\nBinary Search:  Lower with (random, cdf[mid]) = (%30.24e %30.24e)", randomValue, obj -> cdf[mid]);
          ${trace.token}  }
          if (low == mid) {
            low = high;
            if (obj -> cdf[mid] == 0) {
              mid ++;
              ${trace.token}  if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
              ${trace.token}    RF_nativePrint("\nZero level CDF adjustment:  %10d", mid);
              ${trace.token}  }
            }
          }
          else {
            high = mid;
          }
        }
      }
      midValue = obj -> cdf[mid];
      flag = TRUE;
      while (flag) {
        if (mid > 1) {
          if (midValue == obj -> cdf[mid-1]) {
            mid --;
            ${trace.token}  if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
            ${trace.token}    RF_nativePrint("\nFlat-line mid adjustment:  %10d", mid);
            ${trace.token}  }
          }
          else {
            flag = FALSE;
          }
        }
        else {
          flag = FALSE;
        }
      }
      obj -> slot = mid;
      value = obj -> index[obj -> slot];
    }
    else {
      value = obj -> slot = 0;
    }
    break;
  }
  ${trace.token}  if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nSelected element:  %10d", value);
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nsampleFromCDFNew() EXIT ...\n");
  ${trace.token}  }
  return value;
}
void updateCDFNew(uint    treeID, DistributionObj *obj) {
  uint sourcePt;
  uint stepIndex;
  uint currCov, nextCov;
  double oldStepValue, newStepValue;
  char flag;
  uint   i, j, k;
  ${trace.token}  if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdateCDF() ENTRY ...\n");
  ${trace.token}  }
  switch (obj -> weightType) {
  case RF_WGHT_UNIFORM:
    obj -> index[obj -> slot] = obj -> index[obj -> indexSize];
    (obj -> indexSize) --;
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
    ${trace.token}        RF_nativePrint("\nUpdated Incoming Weights Uniform Size:  %10d ", obj -> indexSize);
    ${trace.token}      }
    ${trace.token}      }
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
    ${trace.token}        RF_nativePrint("\nUpdated Incoming Weights:  Uniform");
    ${trace.token}        RF_nativePrint("\n     index      absolute");
    ${trace.token}        for (k=1; k <= obj -> indexSize; k++) {
    ${trace.token}          RF_nativePrint("\n%10d  %10d", k, obj -> index[k]);
    ${trace.token}        }
    ${trace.token}        RF_nativePrint("\n");
    ${trace.token}      }
    ${trace.token}      }
    break;
  case RF_WGHT_INTEGER:
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}        RF_nativePrint("\nIncoming Weights Integer Size (before):  %10d", obj -> densitySize);
    ${trace.token}      }
    currCov = nextCov = obj -> density[obj -> densitySize];
    i = 0;
    j = (uint) (obj -> weight)[currCov];
    k = (uint) (obj -> weight)[obj -> slot];
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}        if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
    ${trace.token}          RF_nativePrint("\nLength of end element weight:       %10d", j);
    ${trace.token}          RF_nativePrint("\nLength of selected element weight:  %10d", k);
    ${trace.token}        }
    ${trace.token}      }
    while(i < k) {
      if (obj -> density[obj -> densitySize] == obj -> slot) {
        ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
        ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
        ${trace.token}        RF_nativePrint("\nWeight iterator while found:  %10d", i);
        ${trace.token}        RF_nativePrint("\nTo be deleted element found at end of density vector at end position:    %10d", obj -> densitySize);
        ${trace.token}        RF_nativePrint("\nTo be deleted element found at end of density vector at swap index:      %10d", k);
        ${trace.token}      }
        ${trace.token}      }
        obj -> density[obj -> densitySize] = 0;
        (obj -> densitySize) --;
        (obj -> densitySwap)[obj -> slot][k] = 0;
        k--;
        if (obj -> densitySize > 0) {
          currCov = nextCov = obj -> density[obj -> densitySize];
          j = (uint) (obj -> weight)[currCov];
        }
      }
      else {
        i++;
        sourcePt = obj -> densitySwap[obj -> slot][i];
        obj -> density[sourcePt] = obj -> density[obj -> densitySize];
        ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
        ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
        ${trace.token}        RF_nativePrint("\nWeight iterator while not found:  %10d", i);
        ${trace.token}        RF_nativePrint("\nTo be deleted element found in density vector at position:      %10d", sourcePt);
        ${trace.token}        RF_nativePrint("\nTo be deleted element found in density vector at swap index:    %10d", i);
        ${trace.token}        RF_nativePrint("\nTo be swapped element found in density vector at end position:  %10d", obj -> densitySize);
        ${trace.token}        RF_nativePrint("\nTo be swapped element found in density vector at swap index:    %10d", j);
        ${trace.token}      }
        ${trace.token}      }
        obj -> density[obj -> densitySize] = 0;
        (obj -> densitySize) --;
        obj -> densitySwap[currCov][j] = obj -> densitySwap[obj -> slot][i];
        obj -> densitySwap[obj -> slot][i] = 0;
        nextCov = obj -> density[obj -> densitySize];
        if (nextCov == currCov) {
          j--;
        }
        else {
          ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
          ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
          ${trace.token}        RF_nativePrint("\nDone swapping element:  ");
          ${trace.token}        RF_nativePrint("\nOld x-var:   %10d", currCov);
          ${trace.token}        RF_nativePrint("\nOld weight:  %10d", (uint) (obj -> weight)[currCov]);
          ${trace.token}        RF_nativePrint("\nNew x-var:   %10d", nextCov);
          ${trace.token}        RF_nativePrint("\nNew weight:  %10d", (uint) (obj -> weight)[nextCov]);
          ${trace.token}      }
          ${trace.token}      }
          hpsortui(obj -> densitySwap[currCov], (uint) (obj -> weight)[currCov]);
          currCov = nextCov = obj -> density[obj -> densitySize];
          j = (uint) (obj -> weight)[currCov];
        }
      }
    }
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
    ${trace.token}        RF_nativePrint("\nWeight Density Vector Integer Size (after):   %10d", obj -> densitySize);
    ${trace.token}      }
    ${trace.token}      }
    if (obj -> densitySize > 0) {
      if (nextCov == currCov) {
        hpsortui(obj -> densitySwap[currCov], (uint) (obj -> weight)[currCov]);
      }
    }
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
    ${trace.token}        RF_nativePrint("\nUpdated Incoming Weights:  Integer");
    ${trace.token}        RF_nativePrint("\n     index          CWDV");
    ${trace.token}        for (i=1; i <= obj -> densitySize; i++) {
    ${trace.token}          RF_nativePrint("\n%10d  %12d", i, obj -> density[i]);
    ${trace.token}        }
    ${trace.token}        RF_nativePrint("\n");
    ${trace.token}      }
    ${trace.token}      }
    break;
  case RF_WGHT_GENERIC:
    stepIndex = obj -> slot;
    if (stepIndex == 1) {
      newStepValue = 0;
    }
    else {
      newStepValue = obj -> cdf[stepIndex - 1];
    }
    oldStepValue = obj -> cdf[stepIndex];
    k = stepIndex;
    flag = TRUE;
    while (flag) {
      if ((obj -> cdf)[k] == oldStepValue) {
        (obj -> cdf)[k] = newStepValue;
        k ++;
      }
      else {
        flag = FALSE;
      }
      if (k > obj -> cdfSize) {
        flag = FALSE;
      }
    }
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
    ${trace.token}        RF_nativePrint("\nUpdated Incoming Weights Generic Step:  %10d", stepIndex);
    ${trace.token}      }
    ${trace.token}      }
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}        RF_nativePrint("\nUpdated Outgoing Weights:   ");
    ${trace.token}        RF_nativePrint("\n     index                            CDF");
    ${trace.token}        for (k=1; k <= obj -> cdfSize; k++) {
    ${trace.token}          RF_nativePrint("\n%10d %30.24e", k, obj -> cdf[k]);
    ${trace.token}        }
    ${trace.token}        RF_nativePrint("\n");
    ${trace.token}      }
    break;
  }
  ${trace.token}  if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nupdateCDF() EXIT ...\n");
  ${trace.token}  }
}
void discardCDFNew(uint treeID, DistributionObj *obj) {
  uint k;
  ${trace.token}  if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\ndiscardCDFNew() ENTRY ...\n");
  ${trace.token}  }
  switch (obj -> weightType) {
  case RF_WGHT_UNIFORM:
    free_uivector(obj -> index, 1, obj -> uIndexAllocSize);
    break;
  case RF_WGHT_INTEGER:
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}        RF_nativePrint("\nVector Weights Integer Size: %10d", obj -> permissibleSize);
    ${trace.token}      }
    ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
    ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
    ${trace.token}        RF_nativePrint("\nVector Weights:  Integer");
    ${trace.token}        RF_nativePrint("\n     index          CWDV");
    ${trace.token}        for (uint i=1; i <= obj -> densityAllocSize; i++) {
    ${trace.token}          RF_nativePrint("\n%10d  %12d", i, obj -> density[i]);
    ${trace.token}        }
    ${trace.token}        RF_nativePrint("\n");
    ${trace.token}      }
    ${trace.token}      }
    free_uivector(obj -> density, 1, obj -> densityAllocSize);
    for (k = 1; k <= obj -> permissibleSize; k++) {
      if (obj -> densitySwap[k] != NULL) {
        ${trace.token}      if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
        ${trace.token}      if (getTraceFlag(treeID) | TURN_OFF_TRACE) {
        ${trace.token}        RF_nativePrint("\nSwap Vector for (Weight Index, of Length):  (%10d, %10d)", k, (uint) (obj -> weight)[k]);
        ${trace.token}        RF_nativePrint("\n     index          swap");
        ${trace.token}        for (int i=1; i <= (uint) (obj -> weight)[k]; i++) {
        ${trace.token}          RF_nativePrint("\n%10d  %12d", i, (obj -> densitySwap)[k][i]);
        ${trace.token}        }
        ${trace.token}        RF_nativePrint("\n");
        ${trace.token}      }
        ${trace.token}      }
        free_uivector(obj -> densitySwap[k], 1, (uint) (obj -> weight)[k]);
        obj -> densitySwap[k] = NULL;
      }
    }
    free_new_vvector(obj -> densitySwap, 1, obj -> permissibleSize, NRUTIL_UPTR);
    break;
  case RF_WGHT_GENERIC:
    free_uivector(obj -> index, 1, obj -> permissibleSize);
    free_dvector(obj -> cdf, 1, obj -> permissibleSize);
    break;
  }
  ${trace.token}  if (getTraceFlag(treeID) & SAMP_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\ndiscardCDFNew() EXIT ...\n");
  ${trace.token}  }
}
uint sampleUniformlyFromVector (uint    treeID,
                                uint   *index,
                                uint    size,
                                uint   *sampleSlot) {
  uint result;
  if (size > 0) {
    (*sampleSlot) = (uint) ceil(ran1B(treeID) * (size * 1.0));
    result = index[*sampleSlot];
  }
  else {
    result = 0;
  }
  return result;
}
