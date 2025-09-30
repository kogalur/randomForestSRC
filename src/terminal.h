#ifndef RF_TERMINAL_H
#define RF_TERMINAL_H
typedef struct terminal Terminal;
struct terminal {
  unsigned int nodeID;
  struct node *mate;
  unsigned int *lmiIndex;
  unsigned int  lmiAllocSize, lmiSize;
  double       *lmiValue;
  unsigned int eTypeSize;
  unsigned int mTimeSize;
  unsigned int eTimeSize;
  unsigned int sTimeSize;
  unsigned int *atRiskCount;
  unsigned int **eventCount;
  unsigned int *eventTimeIndex;
  double **localRatio;
  double **localCSH;
  double **localCIF;
  double *localSurvival;
  double *localNelsonAalen;
  double **CSH;
  double **CIF;
  double *survival;
  double *nelsonAalen;
  double *mortality;
  unsigned int   rnfCount;
  double        *meanResponse;
  unsigned int   rfCount;
  unsigned int  *rfSize;
  unsigned int **multiClassProb;
  double        *maxClass;
  double weight;
  unsigned int membrCount;
  unsigned int *membrStream;
  unsigned int inbagProxy;
  uint repMembrSizeAlloc, oobMembrSizeAlloc, ibgMembrSizeAlloc;
  uint repMembrSize, oobMembrSize, ibgMembrSize;
  uint *repMembrIndx, *oobMembrIndx, *ibgMembrIndx;
};
#endif
