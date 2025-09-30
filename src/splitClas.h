#ifndef RF_SPLIT_CLAS_H
#define RF_SPLIT_CLAS_H
char classificationXwghtSplitCur (uint treeID, Node *parent, SplitInfoMax *splitInfoMax, GreedyObj *greedyMembr, char multImpFlag);
char classificationAreaUnderROCSplit (uint treeID, Node *parent, SplitInfoMax *splitInfoMax, GreedyObj *greedyMembr, char multImpFlag);
char classificationEntropySplit      (uint treeID, Node *parent, SplitInfoMax *splitInfoMax, GreedyObj *greedyMembr, char multImpFlag);
#endif
