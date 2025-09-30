#ifndef RF_SPLIT_MAHALANOBIS_H
#define RF_SPLIT_MAHALANOBIS_H
#include "node.h"
#include "splitInfo.h"
char mahalanobis (uint       treeID,
                  Node      *parent,
                  SplitInfoMax *splitInfoMax,
                  GreedyObj    *greedyMembr,
                  char       multImpFlag);
#endif
