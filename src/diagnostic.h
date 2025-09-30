#ifndef  RF_DIAGNOSTIC_H
#define  RF_DIAGNOSTIC_H
#include "splitInfo.h"
#include "node.h"
#include "terminal.h"
void getSplitObjectInfo(SplitInfo *info);
void getNodeInfo(Node *leaf);
void getTerminalInfo(Terminal *termPtr);
Node *getTerminalNode(uint treeID, uint leaf);
void getRawNodeSize(uint  type,
                    uint  treeID,
                    Node *parent,
                    uint *repMembrIndx,
                    uint *repMembrSize,
                    uint *allMembrIndx,
                    uint *allMembrSize);
void printTreeInfo(uint treeID, Node *parent);
void initTimer(void);
void printTimer(void);
void printParameters(char mode);
#endif
