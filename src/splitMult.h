#ifndef RF_SPLIT_MULT_H
#define RF_SPLIT_MULT_H
#include "node.h"
#include "splitInfo.h"
char unsupervisedSplitMiss(uint       treeID,
                           Node      *parent,
                           SplitInfoMax *splitInfoMax,
                           GreedyObj    *greedyMembr,
                           char       multImpFlag);
char unsupervisedSplitNew(uint       treeID,
                          Node      *parent,
                          SplitInfoMax *splitInfoMax,
                          GreedyObj    *greedyMembr,
                          char       multImpFlag);
char multivariateSplitOld (uint       treeID,
                           Node      *parent,
                           SplitInfoMax *splitInfoMax,
                           GreedyObj    *greedyMembr,
                           char       multImpFlag);
char multivariateSplitNew (uint       treeID,
                           Node      *parent,
                           SplitInfoMax *splitInfoMax,
                           GreedyObj    *greedyMembr,
                           char       multImpFlag);
char multivariateSplitNew3 (uint       treeID,
                           Node      *parent,
                           SplitInfoMax *splitInfoMax,
                           GreedyObj    *greedyMembr,
                           char       multImpFlag);
DistributionObj *stackRandomResponsesSimple(uint treeID, Node *parent);
void unstackRandomResponsesSimple(uint treeID, DistributionObj *obj);
char selectRandomResponsesSimpleVector(uint  treeID,
                                       Node *parent,
                                       DistributionObj *distributionObj,
                                       uint *response,
                                       uint *responseCount);
DistributionObj *stackRandomResponsesGeneric(uint treeID, Node *parent);
void unstackRandomResponsesGeneric(uint treeID, DistributionObj *obj);
char selectRandomResponsesGenericVector(uint     treeID,
                                        Node     *parent,
                                        DistributionObj *distributionObj,
                                        uint     *covariate,
                                        uint     *covariateCount);
#endif
