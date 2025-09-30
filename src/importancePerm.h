#ifndef RF_IMPORTANCE_PERM_H
#define RF_IMPORTANCE_PERM_H
#include "terminal.h"
#include "node.h"
void getPermuteMembership(char       mode,
                          uint       treeID,
                          Terminal **vimpMembership,
                          uint       p);
Node *permuteMembershipGeneric(uint     treeID,
                               Node    *parent,
                               uint     individual,
                               uint     vimpX,
                               double **xArray);
Node *permuteMembershipJIT(uint     treeID,
                           Node    *parent,
                           uint     individual,
                           uint     vimpX,
                           double **xArray);
Node *getMembershipGeneric(uint     treeID,
                           Node    *parent,
                           uint     individual,
                           double **xArray);
Node *getMembershipJIT(uint     treeID,
                       Node    *parent,
                       uint     individual,
                       double **xArray);
void permute(uint ranGenID, uint parallelID, uint n, uint *indx);
#endif
