#ifndef RF_TREE_UTIL_H
#define RF_TREE_UTIL_H
#include "node.h"
#include "terminal.h"
char growTreeRecursive(uint     r,
                       char     rootFlag,
                       char     multImpFlag,
                       uint     b,
                       Node    *parent,
                       uint    *bootMembrIndxIter,
                       uint    *rmbrIterator,
                       uint    *ambrIterator);
void freeTree(uint treeID, Node *parent);
void saveStatistics(char     mode,
                    uint     b,
                    Node    *parent,
                    uint    *offset,
                    double  *spltST,
                    uint    *dpthST);
void initTerminalNodeMembership(uint       treeID,
                                Terminal  *parent,
                                uint      *allMembrIndx,
                                uint       allMembrSize);
void updatePruning(char mode, uint treeID);
void updateCaseDepth(char mode, uint treeID);
#endif
