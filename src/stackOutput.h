#ifndef RF_STACK_OUTPUT_H
#define RF_STACK_OUTPUT_H
#include "node.h"
#include "snpAuxiliaryInfo.h"
void stackDefinedOutputObjects(char      mode,
                               char    **sexpString,
                               Node   ***pRF_root,
                               uint    **pRF_tLeafCount_,
                               double  **pRF_proximity_,
                               double  **pRF_distance_,
                               double  **pRF_weight_,
                               double  **p_imputation_,
                               double ***pRF_sImputeResponsePtr,
                               double ***pRF_sImputePredictorPtr,
                               uint    **pRF_varUsed_,
                               uint   ***pRF_varUsedPtr,
                               double  **p_splitDepth_);
void unstackDefinedOutputObjects(char      mode);
void stackForestObjectsPtrOnly(char mode);
void stackTreeObjectsPtrOnly(char mode, uint treeID);
void stackForestObjectsOutput(char mode);
void writeForestObjectsOutput(char mode);
void unstackForestObjectsPtrOnly(char mode);
void unstackTreeObjectsPtrOnly(uint treeID);
void stackForestObjectsAuxOnly(char mode);
void unstackForestObjectsAuxOnly(char mode);
void unstackAuxStatisticalStructures(char mode);
void restackTermListAndQualitativeObjectsUnknown(uint treeID, uint length);
void verifyAndRegisterCustomSplitRules(void);
extern void registerCustomFunctions(void);
void stackAuxiliaryInfoList(SNPAuxiliaryInfo ***list, uint count);
void allocateAuxiliaryInfo(char   targetFlag,
                           char   type,
                           char  *stringIdentifier,
                           SNPAuxiliaryInfo **list,
                           uint   slot,
                           void  *snpPtr,
                           void  *auxiliaryArrayPtr,
                           uint   dimSize,
                           int   *dim);
uint getAuxDim(char flag, int *dim, uint preIndex, uint postIndex);
void unstackAuxiliaryInfoAndList(char targetFlag, SNPAuxiliaryInfo **list, uint count);
#endif
