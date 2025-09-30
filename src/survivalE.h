#ifndef RF_SURVIVAL_E_H
#define RF_SURVIVAL_E_H
void updateEnsembleSurvival(char mode,
                            uint treeID,
                            char perfFlag);
void getEnsembleMortality(char      mode, 
                          uint      treeID,
                          uint      obsSize,
                          double  **ensembleMRTptr,
                          double   *ensembleDen,
                          double   *mortality);
void getEnsembleMortalityCR(char      mode, 
                            uint      treeID,
                            uint      obsSize,
                            double  **ensembleMRTptr,
                            double   *ensembleDen,
                            double  **cMortality);
void getConditionalConcordanceArrays(uint     j, 
                                     double  *timePtr, 
                                     double  *statusPtr, 
                                     double  *mortalityPtr, 
                                     double  *genericEnsembleDenPtr,
                                     uint    *meIndividualSize,
                                     uint   **eIndividual,
                                     double  *subsettedTime,
                                     double  *subsettedStatus,
                                     double  *subsettedMortality,
                                     double   *subsettedEnsembleDen);
double getConcordanceIndex(int     polarity,
                           uint    size, 
                           double *timePtr, 
                           double *statusPtr, 
                           double *predictedOutcome,
                           double *oobCount);
double getConcordanceIndexNew(int     polarity,
                              uint    size, 
                              double *timePtr, 
                              double *statusPtr, 
                              double *predicted,
                              double *oobCount);
void getCRPerformance (char     mode,
                       uint     obsSize,
                       double **responsePtr,
                       double **yearsLost,
                       double  *denom,
                       double  *performanceVector);
uint getTimeInterestIndex(double *array, uint length, double value);
#endif
