#ifndef RF_IMPORTANCE_H
#define RF_IMPORTANCE_H
#include "terminal.h"
Node *identifyExtrapolatedMembership (Node      *parent,
                                      double  **yShadow,
                                      double  **xShadow);
void getVimpMembership(char       mode,
                       uint       treeID,
                       Terminal **vimpMembership,
                       uint       p);
void updateEnsembleVimp (char       mode,
                         uint       treeID,
                         Terminal **vimpMembership,
                         uint       xVarIdx);
void summarizePerturbedPerformance(char       mode,
                                   uint       treeID,
                                   uint       bb,
                                   uint       p,
                                   double   **responsePtr);
void finalizeVimpPerformance(char mode);
void  stackVimpMembership(char mode, Terminal ***membership);
void  unstackVimpMembership(char mode, Terminal **membership);
void normalizeBlockedEnsembleEstimates(char      mode,
                                       double  **blkEnsembleMRTnum,
                                       double ***blkEnsembleCLSnum,
                                       double  **blkEnsembleRGRnum,
                                       double   *blkEnsembleDen);
void resetBlockedEnsembleEstimates(char mode);
void rfsrc_omp_atomic_update(double *addr, double incr);
uint getVimpRecoverySeedDimension(char mode, uint opt);
#endif
