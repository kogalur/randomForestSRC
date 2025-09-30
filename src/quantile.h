#ifndef RF_QUANTILE_H
#define RF_QUANTILE_H
typedef struct quantileObj QuantileObj;
struct quantileObj {
  double v;
  uint g;
  uint dlt;
  QuantileObj *fwdLink;
  QuantileObj *bakLink;
};
typedef struct lookUpInfo LookUpInfo;
struct lookUpInfo {
  QuantileObj *qPtr;
  LookUpInfo *rootPtr;
  LookUpInfo *leftPtr;
  LookUpInfo *rghtPtr;
};
QuantileObj *makeQuantileObj(double value);
void freeQuantileObj(QuantileObj *obj);
void freeQuantileObjList(QuantileObj *obj);
QuantileObj *insertQuantileObj(uint *qStreamSize, QuantileObj **head, QuantileObj **tail, uint *quantileLinkLength, double value, LookUpInfo **tree);
QuantileObj *findInsertionPoint(QuantileObj *head, double value, LookUpInfo *tree);
double getApproxQuantile(QuantileObj *head, double phi, uint streamSize);
void populateBand(uint p, uint *band);
void makeLookUpTree(LookUpInfo *infoObj, QuantileObj *qObj, uint size, uint depth);
void findApproximateInsertionPoint(QuantileObj *head, LookUpInfo *tree, double value, QuantileObj **insertPtr);
LookUpInfo *makeLookUpInfo(void);
void freeLookUpInfo(LookUpInfo *obj);
void freeLookUpTree(LookUpInfo *obj);
void testQuantile(uint treeID);
#endif
