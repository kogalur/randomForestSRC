#ifndef RF_IMPORTANCE_ANTI_H
#define RF_IMPORTANCE_ANTI_H
#include "terminal.h"
#include "node.h"
void getAntiMembership(char       mode,
                       uint       treeID,
                       Terminal **vimpMembership,
                       uint       p);
Node *antiMembershipGeneric(uint     treeID,
                            Node    *parent,
                            uint     individual,
                            uint     vimpX,
                            double **xArray);
Node *antiMembershipJIT(uint     treeID,
                        Node    *parent,
                        uint     individual,
                        uint     vimpX,
                        double **xArray);
#endif
