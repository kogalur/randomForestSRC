
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "sampling.h"
#include "nrutil.h"
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
    break;
  case RF_WGHT_INTEGER:
    obj -> density = uivector(1, obj -> densityAllocSize);
    obj -> densitySize = 0;
    obj -> densitySwap = (uint **) new_vvector(1, obj -> permissibleSize, NRUTIL_UPTR);
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
        }
        else {
          (obj -> densitySwap)[kk] = NULL;
        }
      }
      else {
        (obj -> densitySwap)[kk] = NULL;
      }
    }
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
    break;
  }
}
uint sampleFromCDFNew (float (*genericGenerator) (uint), uint treeID, DistributionObj *obj) {
  double randomValue;
  double midValue;
  char flag;
  uint low, mid, high, value;
  uint p;
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
      low  = mid = 1;
      high = obj -> cdfSize;
      while (low < high) {
        mid  = (low + high) >> 1;
        if (randomValue > obj -> cdf[mid]) {
          if (low == mid) {
            low = high;
            mid = high;
            if (obj -> cdf[mid] == 0) {
              mid ++;
            }
          }
          else {
            low = mid;
          }
        }
        else {
          if (low == mid) {
            low = high;
            if (obj -> cdf[mid] == 0) {
              mid ++;
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
  return value;
}
void updateCDFNew(uint    treeID, DistributionObj *obj) {
  uint sourcePt;
  uint stepIndex;
  uint currCov, nextCov;
  double oldStepValue, newStepValue;
  char flag;
  uint   i, j, k;
  switch (obj -> weightType) {
  case RF_WGHT_UNIFORM:
    obj -> index[obj -> slot] = obj -> index[obj -> indexSize];
    (obj -> indexSize) --;
    break;
  case RF_WGHT_INTEGER:
    currCov = nextCov = obj -> density[obj -> densitySize];
    i = 0;
    j = (uint) (obj -> weight)[currCov];
    k = (uint) (obj -> weight)[obj -> slot];
    while(i < k) {
      if (obj -> density[obj -> densitySize] == obj -> slot) {
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
        obj -> density[obj -> densitySize] = 0;
        (obj -> densitySize) --;
        obj -> densitySwap[currCov][j] = obj -> densitySwap[obj -> slot][i];
        obj -> densitySwap[obj -> slot][i] = 0;
        nextCov = obj -> density[obj -> densitySize];
        if (nextCov == currCov) {
          j--;
        }
        else {
          hpsortui(obj -> densitySwap[currCov], (uint) (obj -> weight)[currCov]);
          currCov = nextCov = obj -> density[obj -> densitySize];
          j = (uint) (obj -> weight)[currCov];
        }
      }
    }
    if (obj -> densitySize > 0) {
      if (nextCov == currCov) {
        hpsortui(obj -> densitySwap[currCov], (uint) (obj -> weight)[currCov]);
      }
    }
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
    break;
  }
}
void discardCDFNew(uint treeID, DistributionObj *obj) {
  uint k;
  switch (obj -> weightType) {
  case RF_WGHT_UNIFORM:
    free_uivector(obj -> index, 1, obj -> uIndexAllocSize);
    break;
  case RF_WGHT_INTEGER:
    free_uivector(obj -> density, 1, obj -> densityAllocSize);
    for (k = 1; k <= obj -> permissibleSize; k++) {
      if (obj -> densitySwap[k] != NULL) {
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
