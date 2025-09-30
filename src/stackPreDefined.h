#ifndef RF_STACK_PRE_DEFINED_H
#define RF_STACK_PRE_DEFINED_H
#include "node.h"
#include "terminal.h"
void stackIncomingResponseArrays(char mode);
void unstackIncomingResponseArrays(char mode);
void unstackIncomingCovariateArrays(char mode);
void unstackIncomingCovariateArrays(char mode);
void stackIncomingArrays(char mode);
void unstackIncomingArrays(char mode);
void checkInteraction(void);
void stackPreDefinedCommonArrays(char          mode,
                                 Node      ****nodeMembership,
                                 Terminal  ****tTermMembership,
                                 Terminal  ****tTermList,
                                 Node       ***root);
void unstackPreDefinedCommonArrays(char          mode,
                                   Node      ***nodeMembership,
                                   Terminal  ***tTermMembership,
                                   Terminal  ***tTermList,
                                   Node       **root);
void stackPreDefinedGrowthArrays(void);
void unstackPreDefinedGrowthArrays(void);
void stackPreDefinedRestoreArrays(void);
void unstackPreDefinedRestoreArrays(void);
void stackPreDefinedPredictArrays(void);
void unstackPreDefinedPredictArrays(void);
void stackWeights(double *weight,
                  uint    size,
                  uint   *weightType,
                  uint  **weightSorted,
                  uint   *weightDensitySize);
void unstackWeights(uint    weightType,
                    uint    size,
                    uint   *weightSorted);
#endif
