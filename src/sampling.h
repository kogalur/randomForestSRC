#ifndef RF_SAMPLING_H
#define RF_SAMPLING_H
typedef struct distributionObj DistributionObj;
struct distributionObj {
  uint *permissibleIndex;
  char *permissible;
  uint permissibleSize;
  uint *augmentationSize;
  uint weightType;
  double *weight;
  uint *weightSorted;
  uint densityAllocSize;
  double *cdf;
  uint    cdfSize;
  uint   *cdfSort;
  uint   *density;
  uint    densitySize;
  uint  **densitySwap;
  uint *index;
  uint  indexSize;
  uint  uIndexAllocSize;
  uint  slot;
};
DistributionObj *makeDistributionObjRaw(void);
DistributionObj *makeDistributionObjFull(void);
void freeDistributionObjRaw(DistributionObj *obj);
void initializeCDFNew(uint treeID, DistributionObj *obj);
uint sampleFromCDFNew (float (*genericGenerator) (uint), uint treeID, DistributionObj *obj);
void updateCDFNew(uint    treeID, DistributionObj *obj);
void discardCDFNew(uint treeID, DistributionObj *obj);
uint sampleUniformlyFromVector (uint    treeID,
                                uint   *index,
                                uint    size,
                                uint   *sampleSlot);
#endif
