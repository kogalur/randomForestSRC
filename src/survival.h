#ifndef RF_SURVIVAL_H
#define RF_SURVIVAL_H
#include "terminal.h"
void getAtRiskAndEventCount (uint       treeID,
                              Terminal  *parent,
                              uint      *repMembrIndx,
                              uint       repMembrSize,
                              uint      *allMembrIndx,
                              uint       allMembrSize,
                              uint      *rmbrIterator);
void getLocalRatio (uint treeID, Terminal *parent);
void getRevLocalRatio(uint treeID, Terminal *parent);
void getLocalCSH (uint treeID, Terminal *parent);
void getLocalCIF (uint treeID, Terminal *parent);
void mapLocalToTimeInterest(uint      treeID,
                            Terminal *parent,
                            void     *genericLocal,
                            void     *genericGlobal);
void getLocalSurvival (uint treeID, Terminal *parent);
void getLocalNelsonAalen (uint treeID, Terminal *parent);
void getSurvival (uint treeID, Terminal *parent);
void getMortality (uint treeID, Terminal *parent);
void getNelsonAalen (uint treeID, Terminal *parent);
void getCSH (uint treeID, Terminal *parent);
void getCIF (uint treeID, Terminal *parent);
#endif
