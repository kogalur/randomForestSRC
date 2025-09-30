#ifndef RF_TREE_JIT_H
#define RF_TREE_JIT_H
#include "node.h"
#include "terminal.h"
void acquireTreeJIT(char mode, uint r, uint treeID);
void restoreTerminalNodeJIT(uint treeID, Node *root, uint indv, double **xArray, Terminal **termMembership);
#endif
