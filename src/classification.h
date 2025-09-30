#ifndef  RF_CLASSIFICATION_H
#define  RF_CLASSIFICATION_H
#include "terminal.h"
void getMultiClassProb (uint       treeID,
                        Terminal  *parent,
                        uint      *repMembrIndx,
                        uint       repMembrSize,
                        uint      *allMembrIndx,
                        uint       allMembrSize,
                        uint      *rmbrIterator);
void updateEnsembleMultiClass(char     mode,
                              uint     treeID,
                              char     perfFlag,
                              char     omitDenominator);
double getBrierScore(uint     obsSize,
                     uint     rTarget,
                     double  *responsePtr,
                     double **outcomeCLS,
                     double  *denomCount,
                     double  *cpv);
void getConditionalClassificationIndexGrow(uint     size,
                                           uint     rTarget,
                                           double  *responsePtr,
                                           double **outcomeCLS,
                                           double  *maxVote,
                                           double  *denomCount,
                                           double  *cpv);
void getConditionalClassificationIndexPred(uint     size,
                                           uint     rTarget,
                                           double  *responsePtr,
                                           double **outcomeCLS,
                                           double  *maxVote,
                                           double  *denomCount,
                                           double  *cpv);
double getClassificationIndex(uint     size,
                              uint     rTarget,
                              double  *responsePtr,
                              double  *denomCount,
                              double  *maxVote);
double getGMeanIndexGrow(uint    size,
                         uint    rTarget,
                         double *responsePtr,
                         double *denomCount,
                         double *maxVote);
double getGMeanIndexPred(uint    size,
                         uint    rTarget,
                         double *responsePtr,
                         double *denomCount,
                         double *maxVote);
void getMaxVote(uint     size,
                uint     rTarget,
                double **outcomeCLS,
                double  *denomCount,
                double  *maxVote);
#endif
