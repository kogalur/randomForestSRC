#ifndef RF_POLARITY_H
#define RF_POLARITY_H
#include "splitInfo.h"
char getDaughterPolaritySimpleFactor   (uint treeID, SplitInfo *info, uint index, void *value, ...);
char getDaughterPolaritySimpleNonFactor(uint treeID, SplitInfo *info, uint index, void *value, ...);
char getDaughterPolarity               (uint treeID, SplitInfo *info, uint index, void *value, ...);
char getDaughterPolaritySimpleFactorSingle(uint treeID, SplitInfo *info, uint index, void *value, ...);
char getDaughterPolaritySimpleNonFactorSingle(uint treeID, SplitInfo *info, uint index, void *value, ...);
#endif
