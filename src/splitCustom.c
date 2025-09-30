
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include  <stdlib.h>
#include  <math.h>
#include  "splitCustom.h"
void registerCustomFunctions(void) {
  registerThis (&getCustomSplitStatisticMultivariateClassification, CLAS_FAM, 1);
  registerThis (&getCustomSplitStatisticMultivariateRegression, REGR_FAM, 1);
  registerThis (&getCustomSplitStatisticSurvival, SURV_FAM, 1);
  registerThis (&getCustomSplitStatisticCompetingRisk, CRSK_FAM, 1);
}
double getCustomSplitStatisticMultivariateRegression (unsigned int n,
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
                                                      unsigned int featureCount)
{
  double sumLeft, sumRght;
  double sumLeftSqr, sumRghtSqr;
  double delta;
  unsigned int leftSize, rghtSize;
  unsigned int i;
  sumLeft = sumRght = 0.0;
  leftSize = rghtSize = 0;
  delta = 0.0;
  for (i = 1; i <= n; i++) {
    if (membership[i] == LEFT) {
      sumLeft += response[i] - mean;
      leftSize ++;
    }
    else {
      sumRght += response[i] - mean;
      rghtSize ++;
    }
  }
  sumLeftSqr = pow(sumLeft, 2.0) / ((double) leftSize * variance);
  sumRghtSqr = pow(sumRght, 2.0) / ((double) rghtSize * variance);
  delta = sumLeftSqr + sumRghtSqr;
  return delta;
}
double getCustomSplitStatisticMultivariateClassification (unsigned int n,
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
                                                          unsigned int featureCount)
{
  unsigned int   *leftClassProp, *rghtClassProp;
  double sumLeftSqr, sumRghtSqr;
  double delta;
  double sumLeft, sumRght;
  unsigned int leftSize, rghtSize;
  unsigned int i, p;
  sumLeft = sumRght = 0.0;
  leftSize = rghtSize = 0;
  delta = 0.0;
  leftClassProp = alloc_uivector(maxLevel);
  rghtClassProp = alloc_uivector(maxLevel);
  for (p = 1; p <= maxLevel; p++) {
    leftClassProp[p] = rghtClassProp[p] = 0;
  }
  for (i = 1; i <= n; i++) {
    if (membership[i] == LEFT) {
      leftClassProp[(unsigned int) response[i]] ++;
      leftSize ++;
    }
    else {
      rghtClassProp[(unsigned int) response[i]] ++;
      rghtSize ++;
    }
  }
  for (p = 1; p <= maxLevel; p++) {
    sumLeft += pow((double) leftClassProp[p], 2.0);
    sumRght += pow((double) rghtClassProp[p], 2.0);
  }
  dealloc_uivector(leftClassProp, maxLevel);
  dealloc_uivector(rghtClassProp, maxLevel);
  sumLeftSqr = sumLeft / leftSize;
  sumRghtSqr = sumRght / rghtSize;
  delta = sumLeftSqr + sumRghtSqr;
  return delta;
}
double getCustomSplitStatisticSurvival (unsigned int n,
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
                                        unsigned int featureCount)
{
  unsigned int   *nodeLeftEvent,  *nodeParentEvent;
  unsigned int   *nodeLeftAtRisk, *nodeParentAtRisk;
  double delta, deltaNum, deltaDen;
  unsigned int i, k;
  deltaNum = deltaDen = 0.0;
  nodeLeftEvent   = alloc_uivector(eventTimeSize);
  nodeParentEvent = alloc_uivector(eventTimeSize);
  nodeLeftAtRisk    = alloc_uivector(eventTimeSize);
  nodeParentAtRisk  = alloc_uivector(eventTimeSize);
  for (k = 1; k <= eventTimeSize; k++) {
    nodeParentEvent[k]  = 0;
    nodeParentAtRisk[k] = 0;
    nodeLeftEvent[k]    = 0;
    nodeLeftAtRisk[k]   = 0;
  }
  k = eventTimeSize;
  i = n;
  while ((i > 0) && (k > 0)) {
    if (eventTime[k] <= time[i]) {
      nodeParentAtRisk[k] ++;
      if (membership[i] == LEFT) {
        nodeLeftAtRisk[k] ++;
      }
      if (eventTime[k] == time[i]) {
        if (event[i] > 0) {
          nodeParentEvent[k] ++;
          if (membership[i] == LEFT) {
            nodeLeftEvent[k] ++;
          }
        }
      }
      i--;
    }
    else {
      k--;
    }
  }
  for (k = eventTimeSize; k > 1; k--) {
    nodeParentAtRisk[k-1] = nodeParentAtRisk[k] + nodeParentAtRisk[k-1];
    nodeLeftAtRisk[k-1] = nodeLeftAtRisk[k] + nodeLeftAtRisk[k-1];
  }
  for (k = 1; k <= eventTimeSize; k++) {
    deltaNum = deltaNum + ((double) nodeLeftEvent[k] - ((double) ( nodeLeftAtRisk[k] * nodeParentEvent[k]) / nodeParentAtRisk[k]));
    if (nodeParentAtRisk[k] >= 2) {
      deltaDen = deltaDen + (
                             ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k]) *
                             (1.0 - ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k])) *
                             ((double) (nodeParentAtRisk[k] - nodeParentEvent[k]) / (nodeParentAtRisk[k] - 1)) * nodeParentEvent[k]
                             );
    }
  }
  dealloc_uivector(nodeLeftEvent, eventTimeSize);
  dealloc_uivector(nodeParentEvent, eventTimeSize);
  dealloc_uivector(nodeLeftAtRisk, eventTimeSize);
  dealloc_uivector(nodeParentAtRisk, eventTimeSize);
  deltaNum = fabs(deltaNum);
  deltaDen = sqrt(deltaDen);
  if (deltaDen <= 1.0e-9) {
    if (deltaNum <= 1.0e-9) {
      delta = 0.0;
    }
    else {
      delta = deltaNum / deltaDen;
    }
  }
  else {
    delta = deltaNum / deltaDen;
  }
  return delta;
}
double getCustomSplitStatisticCompetingRisk (unsigned int n,
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
                                             unsigned int featureCount)
{
  unsigned int   *nodeLeftEvent,  *nodeParentEvent;
  unsigned int   *nodeLeftAtRisk, *nodeParentAtRisk;
  unsigned int   **nodeLeftEventCR,  **nodeParentEventCR;
  unsigned int   **nodeLeftInclusiveAtRisk, **nodeParentInclusiveAtRisk;
  double delta, deltaNum, deltaSubNum, deltaDen, deltaSubDen;
  unsigned int i, j, k, r, s;
  deltaNum = deltaDen = 0.0;
  nodeLeftEvent   = alloc_uivector(eventTimeSize);
  nodeParentEvent = alloc_uivector(eventTimeSize);
  nodeLeftAtRisk    = alloc_uivector(eventTimeSize);
  nodeParentAtRisk  = alloc_uivector(eventTimeSize);
  nodeParentEventCR = alloc_uimatrix(eventTypeSize, eventTimeSize);
  nodeLeftEventCR = alloc_uimatrix(eventTypeSize, eventTimeSize);
  nodeParentInclusiveAtRisk = alloc_uimatrix(eventTypeSize, eventTimeSize);
  nodeLeftInclusiveAtRisk = alloc_uimatrix(eventTypeSize, eventTimeSize);
  for (k = 1; k <= eventTimeSize; k++) {
    nodeParentEvent[k]  = 0;
    nodeParentAtRisk[k] = 0;
    nodeLeftEvent[k]    = 0;
    nodeLeftAtRisk[k]   = 0;
    for (j = 1; j <= eventTypeSize; j++) {
      nodeParentEventCR[j][k]         = 0;
      nodeLeftEventCR[j][k]           = 0;
      nodeParentInclusiveAtRisk[j][k] = 0;
      nodeLeftInclusiveAtRisk[j][k]   = 0;
    }
  }
  k = eventTimeSize;
  i = n;
  while ((i > 0) && (k > 0)) {
    if (eventTime[k] <= time[i]) {
      nodeParentAtRisk[k] ++;
      if (membership[i] == LEFT) {
        nodeLeftAtRisk[k] ++;
      }
      if (eventTime[k] == time[i]) {
        if (event[i] > 0) {
          nodeParentEventCR[(unsigned int) event[i]][k] ++;
          nodeParentEvent[k] ++;
          if (membership[i] == LEFT) {
            nodeLeftEventCR[(unsigned int) event[i]][k] ++;
          }
        }
      }
      i--;
    }
    else {
      k--;
    }
  }
  for (k = eventTimeSize; k > 1; k--) {
    nodeParentAtRisk[k-1] = nodeParentAtRisk[k] + nodeParentAtRisk[k-1];
    nodeLeftAtRisk[k-1] = nodeLeftAtRisk[k] + nodeLeftAtRisk[k-1];
  }
  for (k = 1; k <= eventTimeSize; k++) {
    for (j = 1; j <= eventTypeSize; j++) {
      nodeParentInclusiveAtRisk[j][k] = nodeParentAtRisk[k];
      nodeLeftInclusiveAtRisk[j][k] = nodeLeftAtRisk[k];
      for (s = 1; s < k; s++) {
        for (r = 1; r <= eventTypeSize; r++) {
          if (j != r) {
            nodeParentInclusiveAtRisk[j][k]  += nodeParentEventCR[r][s];
            nodeLeftInclusiveAtRisk[j][k]  += nodeLeftEventCR[r][s];
          }
        }
      }
    }
  }
  for (j = 1; j <= eventTypeSize; j++) {
    deltaSubNum = 0;
    for (k = 1; k <= eventTimeSize; k++) {
      deltaSubNum = deltaSubNum + (nodeLeftEventCR[j][k] - (nodeParentEventCR[j][k] * ((double) nodeLeftInclusiveAtRisk[j][k] / nodeParentInclusiveAtRisk[j][k])));
    }
    deltaNum = deltaNum + deltaSubNum;
    deltaSubDen = 0;
    for (k = 1; k <= eventTimeSize; k++) {
      if (nodeParentAtRisk[k] >= 2) {
        deltaSubDen = deltaSubDen  + (
                                      (nodeParentEventCR[j][k] * ((double) nodeLeftInclusiveAtRisk[j][k] / nodeParentInclusiveAtRisk[j][k])) *
                                      (1.0 - ((double) nodeLeftInclusiveAtRisk[j][k] / nodeParentInclusiveAtRisk[j][k])) *
                                      ((double) (nodeParentInclusiveAtRisk[j][k] - nodeParentEventCR[j][k]) / (nodeParentInclusiveAtRisk[j][k] - 1))
                                      );
      }
    }
    deltaDen = deltaDen + deltaSubDen;
  }
  dealloc_uivector(nodeLeftEvent, eventTimeSize);
  dealloc_uivector(nodeParentEvent, eventTimeSize);
  dealloc_uivector(nodeLeftAtRisk, eventTimeSize);
  dealloc_uivector(nodeParentAtRisk, eventTimeSize);
  dealloc_uimatrix(nodeParentEventCR, eventTypeSize, eventTimeSize);
  dealloc_uimatrix(nodeLeftEventCR, eventTypeSize, eventTimeSize);
  dealloc_uimatrix(nodeParentInclusiveAtRisk, eventTypeSize, eventTimeSize);
  dealloc_uimatrix(nodeLeftInclusiveAtRisk, eventTypeSize, eventTimeSize);
  deltaNum = fabs(deltaNum);
  deltaDen = sqrt(deltaDen);
  if (deltaDen <= 1.0e-9) {
    if (deltaNum <= 1.0e-9) {
      delta = 0.0;
    }
    else {
      delta = deltaNum / deltaDen;
    }
  }
  else {
    delta = deltaNum / deltaDen;
  }
  return delta;
}
unsigned int *alloc_uivector(unsigned int nh)
{
  return (unsigned int *) malloc((size_t) ((nh+1) * (sizeof(unsigned int))));
}
void dealloc_uivector(unsigned int *v, unsigned int nh)
{
  free((char *) v);
}
double *alloc_dvector(double *v, unsigned int nh)
{
  return (double *) malloc((size_t) ((nh+1) * (sizeof(double))));
}
void dealloc_dvector(double *v, unsigned int nh)
{
  free((char *) v);
}
unsigned int **alloc_uimatrix(unsigned int n2h, unsigned int nh)
{
  unsigned int **v = (unsigned int **) malloc((size_t) ((n2h+1) * (sizeof(unsigned int *))));
  for (unsigned int i = 1; i <= n2h; i++) {
    v[i] = alloc_uivector(nh);
  }
  return v;
}
void dealloc_uimatrix(unsigned int **v, unsigned int n2h, unsigned int nh)
{
  for (unsigned int i = 1; i <= n2h; i++) {
    dealloc_uivector(v[i], nh);
  }
  free((char *) v);
}
