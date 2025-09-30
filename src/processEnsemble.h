#ifndef RF_PROCESS_ENSEMBLE_H
#define RF_PROCESS_ENSEMBLE_H
void stackFactorInSitu(uint treeID);
void unstackFactorInSitu(uint treeID);
void processEnsembleInSitu(char mode, char multImpFlag, uint b);
void processEnsemblePost(char mode);
void processEnsembleHoldout(uint xVarIdx, uint b);
void processEnsembleHoldoutPost(uint bb);
#endif
