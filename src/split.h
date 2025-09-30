#ifndef RF_SPLIT_H
#define RF_SPLIT_H
#include "splitInfo.h"
typedef struct splitRuleObj SplitRuleObj;
struct splitRuleObj {
  char (*function) (uint,
                    Node*,
                    SplitInfoMax*,
                    GreedyObj*,
                    char);
};
SplitRuleObj *makeSplitRuleObj(uint rule);
void freeSplitRuleObj(SplitRuleObj *obj);
char getBestSplit(uint       treeID,
                  Node      *parent,
                  uint       splitRule,
                  SplitInfoMax *splitInfoMax,
                  char       multImpFlag);
char randomSplitGeneric(uint       treeID,
                        Node      *parent,
                        SplitInfoMax *splitInfoMax,
                        GreedyObj    *greedyMembr,
                        char       multImpFlag);
char randomSplitSimple(uint       treeID,
                       Node      *parent,
                       SplitInfoMax *splitInfoMax,
                       GreedyObj    *greedyMembr,
                       char       multImpFlag);
typedef double (*customFunction) (unsigned int n,
                                  char        *membership,
                                  double      *time,
                                  double      *event,
                                  unsigned int eventTypeSize,
                                  unsigned int eventTimeSize,
                                  double      *eventTime,
                                  double      *response,
                                  double       mean,
                                  double       variance,
                                  unsigned int maxLevel,
                                  double     **feature,
                                  unsigned int featureCount);
void regCustomFunctionClassification (customFunction func, uint i);
void regCustomFunctionRegression (customFunction func, uint i);
void regCustomFunctionSurvival (customFunction func, uint i);
void regCustomFunctionCompetingRisk (customFunction func, uint i);
void registerThis (customFunction func, unsigned int family, unsigned int slot);
#endif
