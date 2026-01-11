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
                                     double  *weight,
                                     uint    *meIndividualSize,
                                     uint   **eIndividual,
                                     double  *subsettedTime,
                                     double  *subsettedStatus,
                                     double  *subsettedMortality,
                                     double  *subsettedEnsembleDen,
                                     double  *subsettedWeight);
double getConcordanceIndex(int     fastAction,
                           uint    size, 
                           double *timePtr, 
                           double *statusPtr, 
                           double *predicted,
                           double *denom,
                           double *weight);
double getConcordanceIndexOriginal(uint    size, 
                                   double *timePtr, 
                                   double *statusPtr, 
                                   double *predicted,
                                   double *denom);
double getConcordanceIndexUno(uint    size, 
                              double *timePtr, 
                              double *statusPtr, 
                              double *predicted,
                              double *weight);
double getConcordanceIndexFenwick(uint    size, 
                                  double *timePtr, 
                                  double *statusPtr, 
                                  double *predicted,
                                  double *denom);
double getConcordanceIndexUnoFenwick(uint    size, 
                                     double *timePtr, 
                                     double *statusPtr, 
                                     double *predicted,
                                     double *weight);
double getCRConcordanceIndexIPCW_Fenwick(uint    size,
                                         double *timePtr,
                                         double *statusPtr,
                                         double *predictedOutcome,
                                         double *denom,
                                         double *weight,
                                         uint    eventType);
void getCRPerformance (char     mode,
                       uint     obsSize,
                       double **responsePtr,
                       double **yearsLost,
                       double  *denom,
                       double  *weight,
                       double  *performanceVector);
uint getTimeInterestIndex(double *array, uint length, double value);
#endif
