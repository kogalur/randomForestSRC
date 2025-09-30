#ifndef RF_SPLIT_GREEDY_H
#define RF_SPLIT_GREEDY_H
char summarizeSplitResultGreedy(SplitInfo *info);
SplitInfo *makeSplitInfo(uint indicatorSize);
void freeSplitInfo(SplitInfo *info);
SplitInfoMax *makeSplitInfoMax(uint size);
void freeSplitInfoMax(SplitInfoMax *info);
char forkAndUpdateGeneric(uint       treeID,
                          Node      *parent,
                          uint      *repMembrIndx,
                          uint       repMembrSize,
                          uint      *allMembrIndx,
                          uint       allMembrSize,
                          char       multImpFlag,
                          SplitInfo *info,
                          uint      *leafCount,
                          Node     **nodeMembership);
char forkNode(Node      *parent,
              SplitInfo *info);
void saveTree(uint b, Node *parent, uint *offset);
void restoreTree(char mode, uint b, Node *parent);
void integerToHexString(uint n, char *s);
uint numHexDigits(unsigned n);
double standardVector(uint treeID,
                      char standardFlag,
                      GreedyObj *greedyMembr,
                      double    *rawVector,
                      uint      *repMembrIndx,
                      uint      repMembrSize);
double getL2Loss(uint    treeID,
                 double *response,
                 uint   *repMembrIndx,
                 uint    repMembrSize,
                 uint   *allMembrIndx,
                 uint    allMembrSize,
                 char   *membershipFlag,
                 char    selectFlag);
double getNegLogLikelihood(uint    treeID,
                           uint    maxLevel,
                           double *response,
                           uint   *repMembrIndx,
                           uint    repMembrSize,
                           uint   *allMembrIndx,
                           uint    allMembrSize,
                           char   *membershipFlag,
                           char    selectFlag);
GreedyObj *makeGreedyObj(Node *parent, GreedyObj *head);
void freeGreedyObj(GreedyObj *gObj);
void freeGreedyObjList(GreedyObj *gObj);
GreedyObj *findGreedyObj(GreedyObj *head, Node *parent);
#endif
