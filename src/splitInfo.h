#ifndef RF_SPLIT_INFO_H
#define RF_SPLIT_INFO_H
typedef struct node Node;
typedef struct splitInfo SplitInfo;
struct splitInfo {
  uint     size;
  char    *indicator;
  int    *randomVar;
  uint    *mwcpSizeAbs;
  void   **randomPts;
};
typedef struct greedyObj GreedyObj;
struct greedyObj {
  Node *parent;
  GreedyObj *fwdLink;
  GreedyObj *bakLink;
  GreedyObj *head;
  SplitInfo *splitInfo;
  uint inbagProxy;
  uint nodeID;
  uint depth;
  char leafFlag;
  double *standardResponse;  
  uint *membershipComplement;
  double G_nR_h_l;
  double G_nR_h_r;
  double sgStat;
  double eRisk;
  double oobEmprRisk;
};
typedef struct splitInfoMax SplitInfoMax;
struct splitInfoMax {
  uint   size;
char  *indicator;
  double deltaMax;
  int    splitParameterMax;
  double splitValueMaxCont;
  uint   splitValueMaxFactSize;
  uint  *splitValueMaxFactPtr;
  double splitStatistic;
};
#endif
