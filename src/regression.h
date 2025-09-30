#ifndef RF_REGRESSION_H
#define RF_REGRESSION_H
#include "terminal.h"
void getMeanResponse(uint       treeID,
                     Terminal  *parent,
                     uint      *repMembrIndx,
                     uint       repMembrSize,
                     uint      *allMembrIndx,
                     uint       allMembrSize,
                     uint      *rmbrIterator);
void updateEnsembleMean(char     mode,
                        uint     treeID,
                        char     perfFlag,
                        char     omitDenominator);
double getMeanSquareError(uint    size,
                          double *responsePtr,
                          double *predictedOutcome,
                          double *denomCount);
char getVarianceClassic(uint    repMembrSize,
                        uint   *repMembrIndx,
                        uint    nonMissMembrSize,
                        uint   *nonMissMembrIndx,
                        double *targetResponse,
                        double *mean,
                        double *variance);
char getVarianceClassicNoMiss(uint    repMembrSize,
                              uint   *repMembrIndx,
                              uint    nonMissMembrSize,
                              uint   *nonMissMembrIndx,
                              double *targetResponse,
                              double *mean,
                              double *variance);
char getVarianceSinglePass(uint    repMembrSize,
                           uint   *repMembrIndx,
                           uint    nonMissMembrSize,
                           uint   *nonMissMembrIndx,
                           double *targetResponse,
                           double *mean,
                           double *variance);
char getVarianceDoublePass(uint    repMembrSize,
                           uint   *repMembrIndx,
                           uint    nonMissMembrSize,
                           uint   *nonMissMembrIndx,
                           double *targetResponse,
                           double *mean,
                           double *variance);
void updateQuantileStream(char     mode,
                          uint     treeID);
#endif
