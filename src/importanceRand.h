#ifndef RF_IMPORTANCE_RAND_H
#define RF_IMPORTANCE_RAND_H
#include "importanceRand.h"
#include "terminal.h"
#include "node.h"
void getRandomMembership (char       mode,
                          uint       treeID,
                          Terminal **vimpMembership,
                          uint       p);
Node *randomMembershipGeneric(uint     treeID,
                              Node    *parent,
                              uint     individual,                            
                              uint     vimpX,
                              double **xArray);
Node *randomMembershipJIT(uint     treeID,
                          Node    *parent,
                          uint     individual,
                          uint     vimpX,
                          double **xArray);
#endif
