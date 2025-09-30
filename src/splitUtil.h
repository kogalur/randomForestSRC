#ifndef RF_SPLIT_UTIL_H
#define RF_SPLIT_UTIL_H
#include "node.h"
#include "splitInfo.h"
#include "sampling.h"
char getPreSplitResultGeneric (uint      treeID,
                               Node     *parent,
                               char      multImpFlag,
                               char      multVarFlag);
char getPreSplitResultNoMiss (uint      treeID,
                              Node     *parent,
                              char      multImpFlag,
                              char      multVarFlag);
void unstackPreSplit (char      preliminaryResult,
                      Node     *parent,
                      char      multImpFlag,
                      char      multVarFlag);
void stackSplitPreliminary(uint     nodeSize,
                           char   **localSplitIndicator,
                           double **splitVector);
void unstackSplitPreliminary(uint    nodeSize,
                             char   *localSplitIndicator,
                             double *splitVector);
DistributionObj *stackRandomCovariatesGeneric(uint treeID, Node *parent);
void unstackRandomCovariatesGeneric(uint treeID, DistributionObj *obj);
char selectRandomCovariatesGeneric(uint     treeID,
                                   Node     *parent,
                                   DistributionObj *distributionObj,
                                   char     *factorFlag,
                                   uint     *covariate,
                                   uint     *covariateCount);
uint stackAndConstructSplitVectorGenericPhase1 (uint     treeID,
                                                Node    *parent,
                                                uint     covariate,
                                                ...);
uint stackAndConstructSplitVectorGenericPhase2 (uint     treeID,
                                                Node    *parent,
                                                uint     covariate,
                                                double  *splitVector,
                                                uint     vectorSize,
                                                char    *factorFlag,
                                                char    *deterministicSplitFlag,
                                                uint    *mwcpSizeAbsolute,
                                                void   **splitVectorPtr);
void unstackSplitVectorGeneric(uint   treeID,
                               Node  *parent,
                               uint   splitLength,
                               char   factorFlag,
                               uint   splitVectorSize,
                               uint   mwcpSizeAbsolute,
                               char   deterministicSplitFlag,
                               void  *splitVectorPtr,
                               char   multImpFlag,
                               uint  *indxx);
uint virtuallySplitNodeGeneric(uint  treeID,
                               Node *parent,
                               char  factorFlag,
                               uint  mwcpSizeAbsolute,
                               double *observation,
                               uint *indxx,
                               void *splitVectorPtr,
                               uint  offset,
                               char *localSplitIndicator,
                               uint *leftSize,
                               uint  priorMembrIter,
                               uint *currentMembrIter);
char summarizeSplitResult(SplitInfoMax *splitInfoMax);
char updateMaximumSplitGeneric(uint    treeID,
                               Node   *parent,
                               double  delta,
                               uint    covariate,
                               uint    index,
                               char    factorFlag,
                               uint    mwcpSizeAbsolute,
                               uint    repMembrSize,
                               char  **polarity,
                               void   *splitVectorPtr,
                               SplitInfoMax *splitInfoMax);
void getReweightedRandomPair(uint    treeID,
                             uint    relativefactorSize,
                             uint    absoluteFactorSize,
                             double *absoluteLevel,
                             uint   *result);
void getRandomPair(uint treeID, uint relativeFactorSize, uint absoluteFactorSize, double *absoluteLevel, uint *result);
void createRandomBinaryPair(uint    treeID,
                            uint    relativeFactorSize,
                            uint    absoluteFactorSize,
                            uint    groupSize,
                            double *absolutelevel,
                            uint   *pair);
void convertRelToAbsBinaryPair(uint    treeID,
                               uint    relativeFactorSize,
                               uint    absoluteFactorSize,
                               uint    relativePair,
                               double *absoluteLevel,
                               uint   *pair);
#endif
