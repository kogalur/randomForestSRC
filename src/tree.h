#ifndef RF_TREE_H
#define RF_TREE_H
#include "node.h"
void acquireTreeGeneric(char mode, uint r, uint b);
void updateWeight(char mode, uint b);
void finalizeWeight(char mode);
void updateDistance(char mode, uint b);
void finalizeDistance(char mode);
void updateProximity(char mode, uint b);
void finalizeProximity(char mode);
void updateSplitDepth(uint treeID, Node *rootPtr, uint maxDepth);
char pruneBranch(uint obsSize, uint treeID, Node **nodesAtDepth, uint nadCount, uint ptnTarget, uint ptnCurrent);
uint pruneTree(uint obsSize, uint treeID, uint ptnCount);
void stackAuxiliary(char mode, uint b);
void unstackAuxiliary(char mode, uint b);
void printPseudoTNInfo(char mode, uint b);
void getPTNodeList(Node    *parent,
                   Node   **list,
                   uint    *offset);
void getSplitPath(uint treeID, Node *parent);
void freeSplitPath(uint treeID);
uint getMaximumDepth(Node *parent);
void getNodesAtDepth(Node *parent, uint tagDepth, Node **nodesAtDepth, uint *nadCount);
#endif
