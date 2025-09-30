#ifndef RF_MARGINAL_H
#define RF_MARGINAL_H
#include "node.h"
void getMarginalMembership(char mode, uint treeID);
void releaseMarginalMembership(char mode, uint treeID);
void marginalMembership(uint     treeID,
                        Node    *parent,
                        uint    *gAllMembrIndx,
                        uint     gAllMembrSize,
                        uint     obsSize,
                        double **xArray);
#endif
