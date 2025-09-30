#ifndef RF_SPLIT_REGR_H
#define RF_SPLIT_REGR_H
#include "node.h"
#include "splitInfo.h"
char regressionXwghtSplitCur (uint treeID, Node *parent, SplitInfoMax *splitInfoMax, GreedyObj *greedyMembr, char multImpFlag);
#endif
