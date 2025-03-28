typedef struct augmentationObj AugmentationObj;
struct augmentationObj {
  double **observationIntr;
  double **fobservationIntr;
  double **observationSyth;
  double **fobservationSyth;
  double **observationXS;
  double **fobservationXS;
  double **observationIS;
  double **fobservationIS;
  unsigned int     pairCount;
  unsigned int     sythCount;
  struct node     *lotsRoot;
  unsigned int     lotsSize;
  struct leafLinkedObj *leafLinkedObjHead;
  int    *pairOneX;
  int    *pairTwoX;
  int    *sythX;
  unsigned int obsSize;
  unsigned int fobsSize;
  char *permissible;
};
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef NOT_HAVE_PTHREAD
#include <pthread.h>
#endif

#include <R_ext/Print.h>
#include <Rdefines.h>
#define RF_nativePrint printR
#define RF_nativeError printR
#define RF_nativeExit  exit2R
#define RF_nativeNaN NA_REAL
#define RF_nativeIsNaN ISNA

 
#ifndef NULL
#define NULL 0
#endif
#ifndef TRUE
#define TRUE      0x01
#endif
#ifndef FALSE
#define FALSE     0x00
#endif
#ifndef uint
typedef unsigned int  uint;
#endif
#ifndef ulong
typedef unsigned long ulong;
#endif
#define RF_OUTP_ID   0  
#define RF_STRG_ID   1  
#define RF_ASRG_ID   2  
#define RF_OSRG_ID   3  
#define RF_ACIF_ID   4  
#define RF_OCIF_ID   5  
#define RF_ASRV_ID   6  
#define RF_OSRV_ID   7  
#define RF_AMRT_ID   8  
#define RF_OMRT_ID   9  
#define RF_ACLS_ID  10  
#define RF_OCLS_ID  11  
#define RF_ARGR_ID  12  
#define RF_ORGR_ID  13  
#define RF_ER_SURV  14  
#define RF_ER_CLAS  15  
#define RF_ER_REGR  16  
#define RF_PROX_ID  17  
#define RF_LEAF_CT  18  
#define RF_SEED_ID  19  
#define RF_SEED_VM  20  
#define RF_AQNT_ID  21  
#define RF_OQNT_ID  22  
#define RF_BLK_SRG  23  
#define RF_BLK_CLS  24  
#define RF_BLK_RGR  25  
#define RF_VMP_SRG  26  
#define RF_VMP_CLS  27  
#define RF_VMP_RGR  28  
#define RF_MISS_ID  29  
#define RF_XXXX_30  30  
#define RF_VUSE_ID  31  
#define RF_DPTH_ID  32  
#define RF_MEMB_ID  33  
#define RF_PRUN_ID  34  
#define RF_BOOT_CT  35  
#define RF_RMBR_ID  36  
#define RF_AMBR_ID  37  
#define RF_TN_RCNT  38  
#define RF_TN_ACNT  39  
#define RF_SPLT_ST  40  
#define RF_DPTH_ST  41  
#define RF_WGHT_ID  42  
#define RF_TN_SURV  43  
#define RF_TN_MORT  44  
#define RF_TN_NLSN  45  
#define RF_TN_CSHZ  46  
#define RF_TN_CIFN  47  
#define RF_TN_REGR  48  
#define RF_TN_CLAS  49  
#define RF_TN_KHZF  50  
#define RF_CSE_RGR  51  
#define RF_CSV_RGR  52  
#define RF_CSE_CLS  53  
#define RF_CSV_CLS  54  
#define RF_XXXX_55  55  
#define RF_PART_SR   56  
#define RF_PART_CL   57  
#define RF_PART_RG   58  
#define RF_DIST_ID   59  
#define RF_TREE_ID   60  
#define RF_NODE_ID   61  
#define RF_PARM_ID   62  
#define RF_CONT_PT   63  
#define RF_MWCP_SZ   64  
#define RF_MWCP_PT   65  
#define RF_MWCP_CT   66  
#define RF_HC_DIM    67  
#define RF_CONT_PTR  68  
#define RF_PAIR_CT   69  
#define RF_AUGM_X1   70  
#define RF_AUGM_X2   71  
#define RF_AUGM_XS   72  
#define RF_EMP_RSK     73  
#define RF_OEMP_RSK    74  
#define RF_STAT_LOT    75  
#define RF_HLDOUT_BLK  76  
#define RF_HLDOUT_SRG  77  
#define RF_HLDOUT_CLS  78  
#define RF_HLDOUT_RGR  79  
#define RF_AKHZ_ID       80  
#define RF_OKHZ_ID       81  
#define RF_SYTH_SZ       82  
#define RF_SYTH_TREE_ID  83  
#define RF_SYTH_NODE_ID  84  
#define RF_SYTH_PARM_ID  85  
#define RF_SYTH_CONT_PT  86  
#define RF_SYTH_MWCP_SZ  87  
#define RF_SYTH_MWCP_PT  88  
#define RF_SYTH_MWCP_CT  89  
#define RF_SYTH_HC_DIM   90  
#define RF_SYTH_CONT_PTR 91  
#define RF_SYTH_NODE_CT  92  
#define RF_TDC_MEMB_ID   93  
#define RF_CSE_DEN       94  
#define RF_CSV_DEN       95  
#define RF_BL_NODE_ID     96  
#define RF_BR_NODE_ID     97  
#define RF_FS_REC_ID      98  
#define RF_NODE_SZ        99  
#define RF_CASE_DEPTH    100  
#define RF_OPT_LO_GROW   101  
#define RF_OPT_HI_GROW   102  
#define RF_CPU_TIME      103  
#define RF_TN_OCNT       104  
#define RF_TN_ICNT       105  
#define RF_OMBR_ID       106 
#define RF_IMBR_ID       107 
#define RF_OOB_SZ        108 
#define RF_IBG_SZ        109 
#define RF_SEXP_CNT      110  
#define RF_SEXP_ASCII_SIZE 16
#define OPT_FENS      0x00000001 
#define OPT_OENS      0x00000002 
#define OPT_PERF      0x00000004 
#define OPT_PERF_CALB 0x00000008 
#define OPT_LEAF      0x00000010 
#define OPT_TREE      0x00000020 
#define OPT_SEED      0x00000040 
#define OPT_MISS_OUT  0x00000080 
#define OPT_VIMP_TYP1 0x00000100 
#define OPT_VIMP_TYP2 0x00000200 
#define OPT_VIMP_JOIN 0x00000400 
#define OPT_CASE_DPTH 0x00000800 
#define OPT_VARUSED_F 0x00001000 
#define OPT_VARUSED_T 0x00002000 
#define OPT_PERF_GMN2 0x00004000 
#define OPT_CLAS_RFQ  0x00008000 
#define OPT_IMPU_ONLY 0x00010000 
#define OPT_OUTC_TYPE 0x00020000 
#define OPT_EMPR_RISK 0x00040000 
#define OPT_BOOT_TYP1 0x00080000 
#define OPT_BOOT_TYP2 0x00100000 
#define OPT_COMP_RISK 0x00200000 
#define OPT_SPLDPTH_1 0x00400000 
#define OPT_SPLDPTH_2 0x00800000 
#define OPT_QUANTLE   0x01000000 
#define OPT_VIMP      0x02000000 
#define OPT_ANON      0x04000000 
#define OPT_NODE_STAT 0x08000000 
#define OPT_PROX      0x10000000 
#define OPT_PROX_IBG  0x20000000 
#define OPT_PROX_OOB  0x40000000 
#define OPT_PROX_FUL  0x60000000 
#define OPT_WGHT      0x00000001 
#define OPT_WGHT_IBG  0x00000002 
#define OPT_WGHT_OOB  0x00000004 
#define OPT_WGHT_FUL  0x00000006 
#define OPT_MISS_SKIP 0x00000010 
#define OPT_MEMB_PRUN 0x00000020 
#define OPT_MEMB_USER 0x00000040 
#define OPT_SPLT_CUST 0x00000F00 
#define OPT_BOOT_SWOR 0x00001000 
#define OPT_TREE_ERR  0x00002000 
#define OPT_PART_PLOT 0x00004000 
#define OPT_DATA_PASG 0x00008000 
#define OPT_MEMB_OUTG 0x00010000 
#define OPT_MEMB_INCG 0x00020000 
#define OPT_TERM_OUTG 0x00040000 
#define OPT_TERM_INCG 0x00080000 
#define OPT_DIST      0x00100000 
#define OPT_DIST_IBG  0x00200000 
#define OPT_DIST_OOB  0x00400000 
#define OPT_DIST_FUL  0x00600000 
#define OPT_JIT_TOP   0x00800000 
#define OPT_TDC_BIT1  0x01000000 
#define OPT_TDC_BIT2  0x02000000 
#define OPT_TDC_BIT3  0x04000000 
#define OPT_DATA_PASP 0x08000000 
#define OPT_CSE       0x10000000 
#define OPT_CSV       0x20000000 
#define OPT_EXPERMNT1 0x00000001 
#define OPT_EXPERMNT2 0x00000002 
#define OPT_EXPERMNT3 0x00000004 
#define RF_PART_MORT 1
#define RF_PART_NLSN 2
#define RF_PART_SURV 3
#define RF_PART_YRLS 1
#define RF_PART_CIFN 2 
#define RF_PART_CHFN 3
#define ACTIVE    0x02
#define LEFT      0x01
#define RIGHT     0x02
#define NEITHER   0x00
#define BOTH      0x03
#define EPSILON 1.0e-9
#define SVD_EPS 1.0e-9
#define RF_GROW   0x01
#define RF_PRED   0x02
#define RF_REST   0x04
#define RF_PART   0x08
#define SURV_LGRNK   1
#define SURV_LRSCR   2
#define SURV_CR_LAU  3
#define RAND_SPLIT   4
#define REGR_NRM     5 
#define REGR_WT_OFF  6 
#define REGR_WT_HVY  7 
#define CLAS_NRM     8 
#define CLAS_WT_OFF  9 
#define CLAS_WT_HVY 10 
#define USPV_NRM    11 
#define USPV_WT_OFF 12 
#define USPV_WT_HVY 13 
#define MV_NRM      14 
#define MV_WT_OFF   15 
#define MV_WT_HVY   16 
#define CUST_SPLIT  17
#define REGR_QUANT  18
#define LARG_QUANT  19
#define SURV_BSG1   20
#define CLAS_AU_ROC 21
#define CLAS_ENTROP 22
#define REGR_SGS    23
#define CLAS_SGS    24
#define SURV_SGS    25
#define SURV_TDC    26
#define MAHALANOBIS 27
#define SURV_CR_GEN 28
#define MAXM_SPLIT  29 
#define AUGT_INTR_NONE      0
#define AUGT_INTR_MULT      1
#define AUGT_INTR_DIVS      2
#define AUGT_INTR_ADDT      3
#define AUGT_INTR_SUBT      4
#define AUGT_INTR      1
#define AUGT_SYTH      2
#define RF_VTRY_NULL 0
#define RF_VTRY_HOLD 2
#define CLAS_FAM     0
#define REGR_FAM     1
#define SURV_FAM     2
#define CRSK_FAM     3
#define APROX 0
#define EXACT 1
#define MAX_EXACT_LEVEL sizeof(uint) * 8
#define SAFE_FACTOR_SIZE 16
#define MAX_HYPER_SIZE 15
#define MARGINAL_SIZE 8
#define RF_WGHT_UNIFORM 1
#define RF_WGHT_INTEGER 2
#define RF_WGHT_GENERIC 3
#define RF_DISTANCE_EUCLIDEAN 1
#define NATIVE_TYPE_CHARACTER 0
#define NATIVE_TYPE_INTEGER   1
#define NATIVE_TYPE_NUMERIC   2
#define NATIVE_TYPE_LIST      3
typedef struct node Node;
struct node {
  unsigned int nodeID;
  unsigned int bnodeID;
  unsigned int blnodeID;
  unsigned int brnodeID;
  unsigned int fsrecID;
  struct node *parent;
  struct node *left;
  struct node *right;
  struct terminal *mate;
  unsigned int xSize;
  char *permissible;
  uint *permissibleIndx;
  uint  permissibleIndxSize;
  char  permissibleReIndxFlag;
  char  permissibleOwnershipFlag;
  char splitFlag;
  double splitStatistic;
  double mean;
  double variance;
  unsigned int depth;
  unsigned int *splitDepth;
  char pseudoTerminal;
  unsigned int mpIndexSize;
  unsigned int fmpIndexSize;
  int *mpSign;
  int *fmpSign;
  char imputed;
  unsigned int *lmpIndex;
  unsigned int  lmpIndexAllocSize, lmpIndexActualSize;
  double *lmpValue;
  unsigned int *flmpIndex;
  unsigned int  flmpIndexAllocSize, flmpIndexActualSize;
  double *flmpValue;
  struct splitInfo *splitInfo;
  unsigned int *repMembrIndx;
  unsigned int *allMembrIndx;
  unsigned int  repMembrSizeAlloc;
  unsigned int  allMembrSizeAlloc;
  unsigned int  repMembrSize;
  unsigned int  allMembrSize;
  unsigned int  oobMembrSizeAlloc;
  unsigned int  oobMembrSize;
  unsigned int *oobMembrIndx;
  unsigned int *nonMissMembrIndxStatic;
  unsigned int  nonMissMembrSizeStatic;
  unsigned int *nonMissMembrIndx;
  unsigned int  nonMissMembrSize;
  double timeCutLeft;
  double timeCutRight;
  char   xtdcSplitFlag;
  char   ttdcSplitFlag;
  AugmentationObj *augmentationObj;
  struct node *lotsRoot;
  unsigned int lotsSize;
  uint *minRank;
  uint *maxRank;
  double sumRght;
};
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
  double *atRiskTime;
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
  double *empiricalHazard;
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
  double timeCutLeft;
  double timeCutRight;
  double *localEmpiricalHazard;
  uint repMembrSizeAlloc, oobMembrSizeAlloc, ibgMembrSizeAlloc;
  uint repMembrSize, oobMembrSize, ibgMembrSize;
  uint *repMembrIndx, *oobMembrIndx, *ibgMembrIndx;
};
typedef struct leafLinkedObj LeafLinkedObj;
struct leafLinkedObj {
  struct leafLinkedObj *fwdLink;
  struct leafLinkedObj *bakLink;
  struct node     *nodePtr;
  struct terminal *termPtr;
  struct terminal *termPtrAux;
  uint nodeID;
  uint ibgMembrCount;
  uint allMembrCount;
  uint oobMembrCount;
};
typedef struct leafLinkedObjSimple LeafLinkedObjSimple;
struct leafLinkedObjSimple {
  struct leafLinkedObjSimple *fwdLink;
  struct leafLinkedObjSimple *bakLink;
  struct node     *nodePtr;
};
LeafLinkedObj *makeLeafLinkedObj(void);
LeafLinkedObjSimple *makeLeafLinkedObjSimple(void);
LeafLinkedObj *makeAndSpliceLeafLinkedObj(LeafLinkedObj *tail,
                                          Node *nodePtr,
                                          uint ibgCount,
                                          uint allCount);
LeafLinkedObjSimple *makeAndSpliceLeafLinkedObjSimple(LeafLinkedObjSimple *tail,
                                                      Node *nodePtr);
LeafLinkedObj *makeAndSpliceLeafLinkedObjAux(LeafLinkedObj *tail,
                                             Terminal *termPtrAux);
void freeLeafLinkedObj(LeafLinkedObj *obj);
void freeLeafLinkedObjSimple(LeafLinkedObjSimple *obj);
void freeLeafLinkedObjList(LeafLinkedObj *obj);
void freeLeafLinkedObjListRev(LeafLinkedObj *obj);
typedef struct splitInfo SplitInfo;
struct splitInfo {
  uint     size;
  char    *indicator;
  uint     hcDim;
  int    *randomVar;
  uint    *mwcpSizeAbs;
  void   **randomPts;
  void   **randomPtsRight;
  uint  pairCT;  
  int  *augmX1;  
  int  *augmX2;  
  int  *augmXS;  
  char  sythFlag;
  uint *sythIndx;
  double timeCutLeft;
  double timeCutRight;
  uint splitRank;
};
typedef struct hcDimension HCDimension;
struct hcDimension {
  uint size;
  char    *splitFlag;
  uint hcDim;
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
typedef struct lotObj LatOptTreeObj;
struct lotObj {
  uint firstIn;
  uint lastIn;
  uint size;
  uint strikeout;
  double *risk;
  double firstOD;
  uint treeSize;
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
  uint   splitAugmMaxPairOne;
  uint   splitAugmMaxPairTwo;
  uint   splitAugmMaxSyth;
  uint splitRank;
  double splitStatistic;
};
AugmentationObj *getAugmentationObj(uint treeID, char multImpFlag, Node *parent);
char growAugmentationLOT(char           multImpFlag,
                         uint           treeID,
                         char           augmentationType,
                         Node          *subRoot,
                         uint          *leafCount,
                         Node         **nodeMembership,
                         int           *xSplitList,
                         uint          *xSplitCount,
                         LeafLinkedObj *leafLinkedObjTail);
AugmentationObj *makeAugmentationObj(uint obsSize, uint fobsSize);
void populateAugmentationObjIntr(AugmentationObj *obj,
                                 uint     pairCount,
                                 int     *pairOneX,
                                 int     *pairTwoX,
                                 double **observationIntr);
void populateAugmentationObjSyth(AugmentationObj *obj,
                                 uint             lotsSize,
                                 Node            *lotsRoot,
                                 uint             sythCount,
                                 int             *sythX,
                                 double         **observationSyth);
void populateAugmentationObjXS(AugmentationObj *obj,
                               double         **observationXS);
void populateAugmentationObjIS(AugmentationObj *obj,
                               double         **observationIS);
void restoreAugmentationObj(char       mode,
                            char       type,
                            uint       treeID,
                            uint      *repMembrIndx,
                            uint       repMembrSize,
                            uint      *allMembrIndx,
                            uint       allMembrSize,
                            uint      *ngAllMembrIndx,
                            uint       ngAllMembrSize,
                            SplitInfo *info,
                            Node      *parent);
void freeAugmentationObj(char mode, AugmentationObj *obj);
void freeAugmentationObjList(GreedyObj *gObj);
char restoreNodeMembershipSyth(char  type,
                               uint  treeID,
                               Node *parent,
                               uint *repMembrIndx,
                               uint  repMembrSize,
                               uint *allMembrIndx,
                               uint  allMembrSize,
                               uint *ngAllMembrIndx,
                               uint  ngAllMembrSize,
                               double **rtobservationSyth,
                               LeafLinkedObj **leafLinkedObjTail);
char bootstrap (char     mode,
                uint     treeID,
                Node    *nodePtr,
                uint    *subIndex,
                uint     subsetSize,
                uint    *index,
                uint     indexSize);
char getNodeSign (char mode, uint treeID, Node *nodePtr, uint *bmIndex, uint repMembrSize);
char bootstrapSubject (char     mode,
                       uint     treeID,
                       Node    *nodePtr,
                       uint   **index,
                       uint    *indexSize);
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
void getSplitObjectInfo(SplitInfo *info);
void getNodeInfo(Node *leaf);
void getTerminalInfo(Terminal *termPtr);
Node *getTerminalNode(uint treeID, uint leaf);
void getRawNodeSize(uint  type,
                    uint  treeID,
                    Node *parent,
                    uint *repMembrIndx,
                    uint *repMembrSize,
                    uint *allMembrIndx,
                    uint *allMembrSize);
void printTreeInfo(uint treeID, Node *parent);
void initTimer(void);
void printTimer(void);
void printParameters(char mode);
void processDefaultGrow(void);
void processDefaultPredict(void);
typedef struct factor Factor;
struct factor {
  unsigned int r; 
  unsigned int cardinalGroupCount; 
  void *complementaryPairCount;
  void *cardinalGroupSize; 
  unsigned int **cardinalGroupBinary;
  unsigned int mwcpSize;
};
Factor *makeFactor(uint r, char bookFlag);
void freeFactor(Factor *f);
char bookFactor(Factor *f);
char unbookFactor(Factor *f);
void bookPair (uint    levelCount,
               uint    groupIndex,
               uint    levelIndex,
               uint   *row,
               uint   *level,
               Factor *f);
void nChooseK (uint n, uint r, char type, void *result);
char reduceFraction(uint *numerator, uint *denominator);
char splitOnFactor(uint level, uint *mwcp);
Node *identifyExtrapolatedMembership (Node      *parent,
                                      double  **yShadow,
                                      double  **xShadow);
void getVimpMembership(char       mode,
                       uint       treeID,
                       Terminal **vimpMembership,
                       uint       p);
void updateEnsembleVimp (char       mode,
                         uint       treeID,
                         Terminal **vimpMembership,
                         uint       xVarIdx);
void summarizePerturbedPerformance(char       mode,
                                   uint       treeID,
                                   uint       bb,
                                   uint       p,
                                   double   **responsePtr);
void finalizeVimpPerformance(char mode);
void  stackVimpMembership(char mode, Terminal ***membership);
void  unstackVimpMembership(char mode, Terminal **membership);
void normalizeBlockedEnsembleEstimates(char      mode,
                                       double  **blkEnsembleMRTnum,
                                       double ***blkEnsembleCLSnum,
                                       double  **blkEnsembleRGRnum,
                                       double   *blkEnsembleDen);
void resetBlockedEnsembleEstimates(char mode);
void rfsrc_omp_atomic_update(double *addr, double incr);
uint getVimpRecoverySeedDimension(char mode, uint opt);
void getAntiMembership(char       mode,
                       uint       treeID,
                       Terminal **vimpMembership,
                       uint       p);
Node *antiMembershipGeneric(uint     treeID,
                            Node    *parent,
                            uint     individual,
                            uint     vimpX,
                            double **xArray);
Node *antiMembershipJIT(uint     treeID,
                        Node    *parent,
                        uint     individual,
                        uint     vimpX,
                        double **xArray);
void getPermuteMembership(char       mode,
                          uint       treeID,
                          Terminal **vimpMembership,
                          uint       p);
Node *permuteMembershipGeneric(uint     treeID,
                               Node    *parent,
                               uint     individual,
                               uint     vimpX,
                               double **xArray);
Node *permuteMembershipJIT(uint     treeID,
                           Node    *parent,
                           uint     individual,
                           uint     vimpX,
                           double **xArray);
Node *getMembershipGeneric(uint     treeID,
                           Node    *parent,
                           uint     individual,
                           double **xArray);
Node *getMembershipJIT(uint     treeID,
                       Node    *parent,
                       uint     individual,
                       double **xArray);
void permute(uint ranGenID, uint parallelID, uint n, uint *indx);
void getRandomMembership (char       mode,
                          uint       treeID,
                          Terminal **vimpMembership,
                          uint       p);
Node *randomMembershipGeneric(uint     treeID,
                              Node    *parent,
                              uint     individual,                            
                              uint     vimpX,
                              double **xArray);
Node *randomMembershipJIT(uint     treeID,
                          Node    *parent,
                          uint     individual,
                          uint     vimpX,
                          double **xArray);
char imputeNode (char     type,
                 char     termFlag,
                 char     chainFlag,
                 uint     treeID, 
                 Node    *nodePtr,
                 uint    *repAbsIdx,
                 uint     repNodeSize,
                 uint    *iAbsIdx,
                 uint     iNodeSize);
char restoreNodeMembership(char  mode, 
                           char  rootFlag,
                           uint  treeID, 
                           Node *parent, 
                           uint *repMembrIndx,
                           uint  repMembrSize,
                           uint *allMembrIndx,
                           uint  allMembrSize,
                           uint *ngAllMembrIndx,
                           uint  ngAllMembrSize,
                           uint *bootMembrIndxIter,
                           uint *rmbrIterator,
                           uint *ambrIterator);
void imputeUpdateShadow (char      mode, 
                         double  **shadowResponse, 
                         double  **shadowPredictor);
void imputeNodeAndSummarize(uint     r,
                            char     mode,
                            uint     treeID,
                            Node    *parent,
                            uint    *repMembrIndx,
                            uint     repMembrSize,
                            uint    *allMembrIndx,
                            uint     allMembrSize,
                            uint    *ngAllMembrIndx,
                            uint     ngAllMembrSize);
void imputeSummary(char      mode,
                   char      selectionFlag);
void imputeResponse(char      mode,
                    uint      loSerialTreeID,
                    uint      hiSerialTreeID,
                    uint     *serialTreePtr,
                    double  **tempResponse);
void imputeCommon(char      mode,
                  uint      loSerialTreeID,
                  uint      hiSerialTreeID,
                  uint     *serialTreePtr,
                  char      selectionFlag,
                  char      predictorFlag);
void imputeMultipleTime (char selectionFlag);
double getNearestMasterTime (double   meanvalue,
                             char     chainFlag,
                             uint     treeID);
double getMaximalValue(double *value, uint size, char chainFlag, uint treeID);
double getMedianValue(double *value, uint size);
double getMeanValue(double *value, uint size);
double getSampleValue(double *value, uint size, char chainFlag, uint treeID);
uint getRecordMap(uint     *map, 
                  uint      size, 
                  double  **resp, 
                  double  **data);
void updateTimeIndexArray(uint    treeID,
                          uint   *allMemberIndx,
                          uint    allMembrSize,
                          double *startTime,
                          double *time, 
                          char    naflag,
                          char    idFlag,
                          uint   *startMasterTimeIndex,
                          uint   *masterTimeIndex);
void updateEventTypeSubsets(double *summaryStatus, 
                            uint    mRecordSize,
                            int   **mpSign,
                            uint   *mRecordIndex,
                            uint   *meIndividualSize,
                            uint  **eIndividual);
void stackShadow (char mode, uint treeID);
void unstackShadow (char mode, uint treeID);
char xferMissingness(char type, Node *source, Terminal *destination);
void getMarginalMembership(char mode, uint treeID);
void releaseMarginalMembership(char mode, uint treeID);
void marginalMembership(uint     treeID,
                        Node    *parent,
                        uint    *gAllMembrIndx,
                        uint     gAllMembrSize,
                        uint     obsSize,
                        double **xArray);
Node *makeNode(unsigned int xSize);
void freeNodeGeneric(Node *parent);
void freeNodeNew(Node *parent);
void setParent(
  Node *daughter,
  Node *parent
);
void setLeftDaughter(
   Node *daughter,
   Node *parent
);
void setRightDaughter(
  Node *daughter,
  Node *parent
);
void stackMPSign(Node *node, unsigned int mpIndexSize);
void unstackMPSign(Node *node);
void stackFMPSign(Node *node, unsigned int fmpIndexSize);
void unstackFMPSign(Node *node);
void stackNodeLMPIndex(Node *node, unsigned int size);
void unstackNodeLMPIndex(Node *node);
void stackNodeFLMPIndex(Node *node, unsigned int size);
void unstackNodeFLMPIndex(Node *node);
void stackSplitDepth(Node *tNode, unsigned int depth);
void unstackSplitDepth(Node *tNode);
#define FREE_ARG char*
#define NR_END 2
enum alloc_type{
  NRUTIL_DPTR,   
  NRUTIL_UPTR,   
  NRUTIL_IPTR,   
  NRUTIL_CPTR,   
  NRUTIL_NPTR,   
  NRUTIL_TPTR,   
  NRUTIL_FPTR,   
  NRUTIL_LPTR,   
  NRUTIL_DPTR2,  
  NRUTIL_UPTR2,  
  NRUTIL_IPTR2,  
  NRUTIL_CPTR2,  
  NRUTIL_NPTR2,  
  NRUTIL_TPTR2,  
  NRUTIL_FPTR2,  
  NRUTIL_DPTR3,  
  NRUTIL_UPTR3,  
  NRUTIL_NPTR3,  
  NRUTIL_DPTR4,  
  NRUTIL_UPTR4,  
  NRUTIL_XPTR,   
  NRUTIL_QPTR,   
  NRUTIL_QPTR2,  
  NRUTIL_SPTR,   
  NRUTIL_SPTR2,  
  NRUTIL_VPTR,   
  NRUTIL_GPTR,   
  NRUTIL_OMPLPTR,  
  NRUTIL_OMPLPTR2, 
  NRUTIL_PTMLPTR,  
  NRUTIL_PTMLPTR2, 
  NRUTIL_LEAFPTR,  
  NRUTIL_LEAFPTR2, 
  NRUTIL_SRTLNKPTR, 
  NRUTIL_TARPTR,   
};
unsigned int upower (unsigned int x, unsigned int n);
unsigned int upower2 (unsigned int n);
unsigned int ulog2 (unsigned int n);
void hpsort(double *ra, unsigned int n);
void hpsortui(unsigned int *ra, unsigned int n);
void hpsorti(int *ra, unsigned int n);
void qksort(double *arr, unsigned int n);
void indexx(unsigned int n, double *arr, unsigned int *indx);
void nrerror(char error_text[]);
void *gblock(size_t size);
void free_gblock(void *v, size_t size);
void *gvector(unsigned long long nl, unsigned long long nh, size_t size);
void free_gvector(void *v, unsigned long long nl, unsigned long long nh, size_t size);
char *cvector(unsigned long long nl, unsigned long long nh);
void free_cvector(char *v, unsigned long long nl, unsigned long long nh);
char **cmatrix(unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch);
void free_cmatrix(char **v, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch);
int *ivector(unsigned long long nl, unsigned long long nh);
void free_ivector(int *v, unsigned long long nl, unsigned long long nh);
int **imatrix(unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch);
void free_imatrix(int **v, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch);
unsigned int *uivector(unsigned long long nl, unsigned long long nh);
void free_uivector(unsigned int *v, unsigned long long nl, unsigned long long nh);
unsigned int **uimatrix(unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch);
void free_uimatrix(unsigned int **v, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch);
unsigned long *ulvector(unsigned long long nl, unsigned long long nh);
void free_ulvector(unsigned long *v, unsigned long long nl, unsigned long long nh);
double *dvector(unsigned long long nl, unsigned long long nh);
void free_dvector(double *v, unsigned long long nl, unsigned long long nh);
double **dmatrix(unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch);
void free_dmatrix(double **v, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch);
double ***dmatrix3(unsigned long long n3l, unsigned long long n3h, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch);
void free_dmatrix3(double ***v, unsigned long long n3l, unsigned long long n3h, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch);
double ****dmatrix4(unsigned long long n4l, unsigned long long n4h, unsigned long long n3l, unsigned long long n3h, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch);
void free_dmatrix4(double ****v, unsigned long long n4l, unsigned long long n4h, unsigned long long n3l, unsigned long long n3h, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch);
#ifdef _OPENMP
omp_lock_t *ompvector(unsigned long long nl, unsigned long long nh);
void free_ompvector(omp_lock_t *v, unsigned long long nl, unsigned long long nh);
#endif
#ifdef NOT_HAVE_PTHREAD
pthread_mutex_t *ptmvector(unsigned long long nl, unsigned long long nh);
void free_ptmvector(pthread_mutex_t *v, unsigned long long nl, unsigned long long nh);
#endif
#ifdef NOT_HAVE_PTHREAD
pthread_t *pthreadvector(unsigned long long nl, unsigned long long nh);
void free_pthreadvector(pthread_t *v, unsigned long long nl, unsigned long long nh);
#endif
void *new_vvector(unsigned long long nl, unsigned long long nh, enum alloc_type type);
void free_new_vvector(void *v, unsigned long long nl, unsigned long long nh, enum alloc_type type);
void nrCopyMatrix(
  unsigned int **new,
  unsigned int **old,
  unsigned int nrow,
  unsigned int ncol
);
void nrCopyVector(
  char *new,
  char *old,
  unsigned int ncol
);
void testEndianness(void);
void acquireTreePosix (void *arg);
typedef struct _acquireTreeArg acquireTreeArg;
struct _acquireTreeArg {
  char mode;
  uint r;
  uint b;
};
void execPreSortPosix (void *arg);
typedef struct _execPreSortArg execPreSortArg;
struct _execPreSortArg {
  uint  b;
  uint  p;
  uint *index;
  uint  indexSize;
};
void processEnsembleHoldoutPostPosix (void *arg);
typedef struct _processEnsembleHoldoutPostArg processEnsembleHoldoutPostArg;
struct _processEnsembleHoldoutPostArg {
  uint bb;
};
void getAndUpdatePartialMembership(uint treeID, Node *root);
void partialMembershipGeneric(uint       treeID,
                              Node      *parent,
                              uint       partialIndex,
                              uint      *allMembrIndx,
                              uint       allMembrSize,
                              double   **xArray,
                              Terminal **membership);
void partialMembershipJIT(uint       treeID,
                          Node      *root,
                          uint       partialIndex,
                          uint      *nullMembrIndx,
                          uint       individual,
                          double   **xArray,
                          Terminal **membership);
void updatePartialCalculations (uint       treeID,
                                uint       pVarIdx,
                                Terminal **partialMembership);
void summarizePartialCalculations(uint       treeID,
                                  uint       pVarIdx);
char getDaughterPolaritySimpleFactor   (uint treeID, SplitInfo *info, uint index, void *value, ...);
char getDaughterPolaritySimpleNonFactor(uint treeID, SplitInfo *info, uint index, void *value, ...);
char getDaughterPolaritySimpleTime     (uint treeID, SplitInfo *info, uint index, void *value, ...);
char getDaughterPolarityComplex        (uint treeID, SplitInfo *info, uint index, void *value, ...);
char getDaughterPolarity               (uint treeID, SplitInfo *info, uint index, void *value, ...);
char getDaughterPolaritySimpleFactorSingle(uint treeID, SplitInfo *info, uint index, void *value, ...);
char getDaughterPolaritySimpleNonFactor(uint treeID, SplitInfo *info, uint index, void *value, ...);
void stackFactorInSitu(uint treeID);
void unstackFactorInSitu(uint treeID);
void processEnsembleInSitu(char mode, char multImpFlag, uint b);
void processEnsemblePost(char mode);
void processEnsembleHoldout(uint xVarIdx, uint b);
void processEnsembleHoldoutPost(uint bb);
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
void randomStack(uint bSize, uint bnpSize);
void randomUnstack(uint bSize, uint bnpSize);
void randomSetChainParallel(uint b, int value);
void randomSetUChainParallel(uint b, int value);
void randomSetUChainParallelVimp(uint b, int value);
void randomSetChainParallelVimp(uint p, int value);
void randomSetChainSerial(uint b, int value);
void randomSetUChainSerial(uint b, int value);
void randomSetUChainSerialVimp(uint b, int value);
void randomSetChainSerialVimp(uint p, int value);
int randomGetChainParallel(uint b);
int randomGetUChainParallel(uint b);
int randomGetUChainParallelVimp(uint b);
int randomGetChainParallelVimp(uint p);
int randomGetChainSerial(uint b);
int randomGetUChainSerial(uint b);
int randomGetUChainSerialVimp(uint b);
int randomGetChainSerialVimp(uint p);
float randomChainParallel(uint b);
float randomUChainParallel(uint b);
float randomUChainParallelVimp(uint b);
float randomChainParallelVimp(uint p);
float randomChainSerial(uint b);
float randomUChainSerial(uint b);
float randomUChainSerialVimp(uint b);
float randomChainSerialVimp(uint p);
float ran1_generic(int *iy, int *iv, int *idum);
void lcgenerator(unsigned int *seed, unsigned char reset);
float ran1_original(int *idum);
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
void rfsrc(char mode, int seedValue);
void updateTerminalNodeOutcomes(char       mode,
                                uint       treeID,
                                Terminal  *parent,
                                uint      *repMembrIndx,
                                uint       repMembrSize,
                                uint      *allMembrIndx,
                                uint       allMembrSize,
                                uint      *rbmrIterator,
                                uint      *ambrIterator);
void getMembrCountOnly (uint       treeID,
                        Terminal  *parent,
                        uint      *repMembrIndx,
                        uint       repMembrSize,
                        uint      *allMembrIndx,
                        uint       allMembrSize);
void updateEnsemble (char mode, uint b);
void summarizeFaithfulBlockPerformance (char        mode,
                                        uint        b,
                                        uint        blockID,
                                        double    **blkEnsembleMRTnum,
                                        double   ***blkEnsembleCLSnum,
                                        double    **blkEnsembleRGRnum,
                                        double     *blkEnsembleDen,
                                        double    **responsePtr,
                                        double    **perfMRTblk,
                                        double   ***perfCLSblk,
                                        double    **perfRGRblk);
void summarizeHoldoutBlockPerformance (char        mode,
                                       uint        b,
                                       uint        xVarIdx,
                                       uint        blockID,
                                       double    **responsePtr,
                                       double    **holdMRTstd,
                                       double   ***holdCLSstd,
                                       double    **holdRGRstd,
                                       double     *holdEnsembleDen,
                                       double     *holdMRTptr,
                                       double    **holdCLSptr,
                                       double     *holdRGRptr);
char stackAndImputePerfResponse(char      mode,
                                char      multImpFlag,
                                uint      b,
                                uint      loSerialTreeID,
                                uint      hiSerialTreeID,
                                uint     *serialTreePtr,
                                double ***responsePtr);
void unstackPerfResponse(char mode, char flag, double **mResponsePtr);
void getPerformance(uint      serialTreeID,
                    char      mode,
                    uint      obsSize,
                    double  **responsePtr,
                    double    *denomPtr,
                    double   **outcomeMRT,
                    double  ***outcomeCLS,
                    double   **outcomeRGR,
                    double   *perfMRTptr,
                    double  **perfCLSptr,
                    double   *perfRGRptr);
void normalizeEnsembleEstimates(char mode, char final);
char getPerfFlag (char mode, uint serialTreeID);
void getVariablesUsed(uint treeID, Node *rootPtr, uint *varUsedVector);
typedef struct distributionObj DistributionObj;
struct distributionObj {
  uint *permissibleIndex;
  char *permissible;
  uint permissibleSize;
  uint *augmentationSize;
  uint weightType;
  double *weight;
  uint *weightSorted;
  uint densityAllocSize;
  double *cdf;
  uint    cdfSize;
  uint   *cdfSort;
  uint   *density;
  uint    densitySize;
  uint  **densitySwap;
  uint *index;
  uint  indexSize;
  uint  uIndexAllocSize;
  uint  slot;
};
DistributionObj *makeDistributionObjRaw(void);
DistributionObj *makeDistributionObjFull(void);
void freeDistributionObjRaw(DistributionObj *obj);
void initializeCDFNew(uint treeID, DistributionObj *obj);
uint sampleFromCDFNew (float (*genericGenerator) (uint), uint treeID, DistributionObj *obj);
void updateCDFNew(uint    treeID, DistributionObj *obj);
void discardCDFNew(uint treeID, DistributionObj *obj);
uint sampleUniformlyFromVector (uint    treeID,
                                uint   *index,
                                uint    size,
                                uint   *sampleSlot);
typedef struct snpAuxiliaryInfo SNPAuxiliaryInfo;
struct snpAuxiliaryInfo {
  char type;
  char *identity;
  uint slot;
  ulong linearSize;
  void *snpPtr;
  void *auxiliaryArrayPtr;
  uint dimSize;
  int *dim;
};
typedef struct sortedLinkedObj SortedLinkedObj;
struct sortedLinkedObj {
  struct sortedLinkedObj *fwdLink;
  struct sortedLinkedObj *bakLink;
  uint rank;
  uint indx;
};
SortedLinkedObj *makeSortedLinkedObj(void);
void makeAndSpliceSortedLinkedObj(uint treeID,
                                  SortedLinkedObj **headPtr,
                                  SortedLinkedObj **tailPtr,
                                  uint *listLength,
                                  uint rank, uint indx);
void freeSortedLinkedObjList(SortedLinkedObj *obj);
void freeSortedLinkedObj(SortedLinkedObj *obj);
typedef struct splitRuleObj SplitRuleObj;
struct splitRuleObj {
  char (*function) (uint,
                    Node*,
                    SplitInfoMax*,
                    GreedyObj*,
                    char);
};
SplitRuleObj *makeSplitRuleObj(uint rule);
SplitRuleObj *makeSplitRuleObjGreedy(uint rule);
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
char classificationXwghtSplitCur (uint treeID, Node *parent, SplitInfoMax *splitInfoMax, GreedyObj *greedyMembr, char multImpFlag);
char classificationAreaUnderROCSplit (uint treeID, Node *parent, SplitInfoMax *splitInfoMax, GreedyObj *greedyMembr, char multImpFlag);
char classificationEntropySplit      (uint treeID, Node *parent, SplitInfoMax *splitInfoMax, GreedyObj *greedyMembr, char multImpFlag);
char customMultivariateSplit (uint       treeID,
                              Node      *parent,
                              SplitInfoMax *splitInfoMax,
                              GreedyObj    *greedyMembr,
                              char       multImpFlag);
char customSurvivalSplit (uint       treeID,
                          Node      *parent,
                          SplitInfoMax *splitInfoMax,
                          GreedyObj    *greedyMembr,
                          char       multImpFlag);
char customCompetingRiskSplit (uint       treeID,
                               Node      *parent,
                               SplitInfoMax *splitInfoMax,
                               GreedyObj    *greedyMembr,
                               char       multImpFlag);
char getBestSplitLOT(uint         treeID,
                     char         multImpFlag,
                     GreedyObj   *greedyMembr,
                     char         augmFlag,
                     uint         hdimProxy);
char summarizeSplitResultGreedy(SplitInfo *info);
HCDimension *makeHCDimension(uint size);
void freeHCDimension(HCDimension *obj);
SplitInfo *makeSplitInfo(uint indicatorSize);
void freeSplitInfo(SplitInfo *info);
SplitInfoMax *makeSplitInfoMax(uint size);
void freeSplitInfoMax(SplitInfoMax *info);
char forkAndUpdateGeneric(uint       treeID,
                          Node      *parent,
                          uint      *repMembrIndx,
                          uint       repMembrSize,
                          uint      *allMembrIndx,
                          uint       allMembrSize,
                          char       multImpFlag,
                          SplitInfo *info,
                          uint      *leafCount,
                          Node     **nodeMembership);
char forkNode(Node      *parent,
              SplitInfo *info);
void saveTree(uint b, Node *parent, uint *offset, uint *offsetSyth);
void saveTreeSyth(uint b, Node *parent, uint *offset);
void restoreTree(char mode, uint b, Node *parent);
void restoreTreeSyth(uint    b, Node   *parent);
void integerToHexString(uint n, char *s);
uint numHexDigits(unsigned n);
char growHyperCube(char       rootFlag,
                   char       multImpFlag,
                   uint       treeID,
                   uint       hdim,
                   uint      *leafCount,
                   Node     **nodeMembership,
                   GreedyObj *greedyMembr);
double standardVector(uint treeID,
                      char standardFlag,
                      GreedyObj *greedyMembr,
                      double    *rawVector,
                      uint      *repMembrIndx,
                      uint      repMembrSize);
double getL2Loss(uint    treeID,
                 double *response,
                 uint   *repMembrIndx,
                 uint    repMembrSize,
                 uint   *allMembrIndx,
                 uint    allMembrSize,
                 char   *membershipFlag,
                 char    selectFlag);
double getNegLogLikelihood(uint    treeID,
                           uint    maxLevel,
                           double *response,
                           uint   *repMembrIndx,
                           uint    repMembrSize,
                           uint   *allMembrIndx,
                           uint    allMembrSize,
                           char   *membershipFlag,
                           char    selectFlag);
void defineHyperCubeDimension(uint  treeID,
                              Node *parent,
                              uint  proxy,
                              uint  depth,
                              HCDimension *obj);
void defineHyperCube(uint  treeID,
                     Node *parent,
                     uint  proxy,
                     uint  depth,
                     uint *hcCount,
                     uint *hcMapping,
                     SplitInfo *splitInfo);
GreedyObj *makeGreedyObj(Node *parent, GreedyObj *head);
void freeGreedyObj(GreedyObj *gObj);
void freeGreedyObjList(GreedyObj *gObj);
GreedyObj *findGreedyObj(GreedyObj *head, Node *parent);
char getBestSplitHyperCube(uint       treeID,
                           Node      *parent,
                           SplitInfoMax *splitInfoMax,
                           GreedyObj    *greedyMembr,
                           char       multImpFlag);
char regressionSGS (uint       treeID,
                    Node      *parent,
                    SplitInfoMax *splitInfoMax,
                    GreedyObj    *greedyMembr,
                    char       multImpFlag);
char classificationSGS (uint       treeID,
                        Node      *parent,
                        SplitInfoMax *splitInfoMax,
                        GreedyObj    *greedyMembr,
                        char       multImpFlag);
char randomSGS (uint       treeID,
                Node      *parent,
                SplitInfoMax *splitInfoMax,
                GreedyObj    *greedyMembr,
                char       multImpFlag);
LatOptTreeObj *makeLatOptTreeObj(void);
void freeLatOptTreeObj(LatOptTreeObj *lotObj);
void insertRisk(uint treeID, LatOptTreeObj *obj, double value);
char mahalanobis (uint       treeID,
                  Node      *parent,
                  SplitInfoMax *splitInfoMax,
                  GreedyObj    *greedyMembr,
                  char       multImpFlag);
char unsupervisedSplitMiss(uint       treeID,
                           Node      *parent,
                           SplitInfoMax *splitInfoMax,
                           GreedyObj    *greedyMembr,
                           char       multImpFlag);
char unsupervisedSplitNew(uint       treeID,
                          Node      *parent,
                          SplitInfoMax *splitInfoMax,
                          GreedyObj    *greedyMembr,
                          char       multImpFlag);
char multivariateSplitOld (uint       treeID,
                           Node      *parent,
                           SplitInfoMax *splitInfoMax,
                           GreedyObj    *greedyMembr,
                           char       multImpFlag);
char multivariateSplitNew (uint       treeID,
                           Node      *parent,
                           SplitInfoMax *splitInfoMax,
                           GreedyObj    *greedyMembr,
                           char       multImpFlag);
char multivariateSplitNew3 (uint       treeID,
                           Node      *parent,
                           SplitInfoMax *splitInfoMax,
                           GreedyObj    *greedyMembr,
                           char       multImpFlag);
DistributionObj *stackRandomResponsesSimple(uint treeID, Node *parent);
void unstackRandomResponsesSimple(uint treeID, DistributionObj *obj);
char selectRandomResponsesSimpleVector(uint  treeID,
                                       Node *parent,
                                       DistributionObj *distributionObj,
                                       uint *response,
                                       uint *responseCount);
DistributionObj *stackRandomResponsesGeneric(uint treeID, Node *parent);
void unstackRandomResponsesGeneric(uint treeID, DistributionObj *obj);
char selectRandomResponsesGenericVector(uint     treeID,
                                        Node     *parent,
                                        DistributionObj *distributionObj,
                                        uint     *covariate,
                                        uint     *covariateCount);
char locallyAdaptiveQuantileRegrSplit (uint       treeID,
                                       Node      *parent,
                                       SplitInfoMax *splitInfoMax,
                                       GreedyObj    *greedyMembr,
                                       char       multImpFlag);
char quantileRegrSplit (uint       treeID,
                        Node      *parent,
                        SplitInfoMax *splitInfoMax,
                        GreedyObj    *greedyMembr,
                        char       multImpFlag);
double quantile7 (double *r, uint s, double p);
char regressionXwghtSplitCur (uint treeID, Node *parent, SplitInfoMax *splitInfoMax, GreedyObj *greedyMembr, char multImpFlag);
char regressionXwghtSplitOpt (uint treeID, Node *parent, SplitInfoMax *splitInfoMax, GreedyObj *greedyMembr, char multImpFlag);
char logRankNCR(uint       treeID,
                Node      *parent,
                SplitInfoMax *splitInfoMax,
                GreedyObj    *greedyMembr,
                char       multImpFlag);
char logRankCR(uint       treeID,
               Node      *parent,
               SplitInfoMax *splitInfoMax,
               GreedyObj    *greedyMembr,
               char       multImpFlag);
char wiBrierScore (uint       treeID,
                   Node      *parent,
                   SplitInfoMax *splitInfoMax,
                   GreedyObj    *greedyMembr,
                   char       multImpFlag);
char brierScoreGradient1 (uint       treeID,
                          Node      *parent,
                          SplitInfoMax *splitInfoMax,
                          GreedyObj    *greedyMembr,
                          char       multImpFlag);
char tdcGradient(uint       treeID,
                 Node      *parent,
                 SplitInfoMax *splitInfoMax,
                 GreedyObj    *greedyMembr,
                 char       multImpFlag);
char getPreSplitResultGeneric (uint      treeID,
                               Node     *parent,
                               char      multImpFlag,
                               char      multVarFlag);
char getPreSplitResultNoMiss (uint      treeID,
                              Node     *parent,
                              char      multImpFlag,
                              char      multVarFlag);
void unstackPreSplit (char      preliminaryResult,
                      Node     *parent,
                      char      multImpFlag,
                      char      multVarFlag);
void stackSplitPreliminary(uint     nodeSize,
                           char   **localSplitIndicator,
                           double **splitVector);
void unstackSplitPreliminary(uint    nodeSize,
                             char   *localSplitIndicator,
                             double *splitVector);
DistributionObj *stackRandomCovariatesGeneric(uint treeID, Node *parent);
void unstackRandomCovariatesGeneric(uint treeID, DistributionObj *obj);
char selectRandomCovariatesGeneric(uint     treeID,
                                   Node     *parent,
                                   DistributionObj *distributionObj,
                                   char     *factorFlag,
                                   uint     *covariate,
                                   uint     *covariateCount);
uint stackAndConstructSplitVectorGenericPhase1 (uint     treeID,
                                                Node    *parent,
                                                uint     covariate,
                                                ...);
uint stackAndConstructSplitVectorGenericPhase2 (uint     treeID,
                                                Node    *parent,
                                                uint     covariate,
                                                double  *splitVector,
                                                uint     vectorSize,
                                                char    *factorFlag,
                                                char    *deterministicSplitFlag,
                                                uint    *mwcpSizeAbsolute,
                                                void   **splitVectorPtr);
void unstackSplitVectorGeneric(uint   treeID,
                               Node  *parent,
                               uint   splitLength,
                               char   factorFlag,
                               uint   splitVectorSize,
                               uint   mwcpSizeAbsolute,
                               char   deterministicSplitFlag,
                               void  *splitVectorPtr,
                               char   multImpFlag,
                               uint  *indxx);
uint virtuallySplitNodeGeneric(uint  treeID,
                               Node *parent,
                               char  factorFlag,
                               uint  mwcpSizeAbsolute,
                               double *observation,
                               uint *indxx,
                               void *splitVectorPtr,
                               uint  offset,
                               char *localSplitIndicator,
                               uint *leftSize,
                               uint  priorMembrIter,
                               uint *currentMembrIter);
char summarizeSplitResult(SplitInfoMax *splitInfoMax);
char updateMaximumSplitGeneric(uint    treeID,
                               Node   *parent,
                               double  delta,
                               uint    covariate,
                               uint    index,
                               char    factorFlag,
                               uint    mwcpSizeAbsolute,
                               uint    repMembrSize,
                               char  **polarity,
                               void   *splitVectorPtr,
                               SplitInfoMax *splitInfoMax);
void getReweightedRandomPair(uint    treeID,
                             uint    relativefactorSize,
                             uint    absoluteFactorSize,
                             double *absoluteLevel,
                             uint   *result);
void getRandomPair(uint treeID, uint relativeFactorSize, uint absoluteFactorSize, double *absoluteLevel, uint *result);
void createRandomBinaryPair(uint    treeID,
                            uint    relativeFactorSize,
                            uint    absoluteFactorSize,
                            uint    groupSize,
                            double *absolutelevel,
                            uint   *pair);
void convertRelToAbsBinaryPair(uint    treeID,
                               uint    relativeFactorSize,
                               uint    absoluteFactorSize,
                               uint    relativePair,
                               double *absoluteLevel,
                               uint   *pair);
void initPreSortExtra(void);
void initPreSortIntra(uint treeID);
void execPreSort(uint treeID, uint xvar, uint *membrIndx, uint membrSize);
void freePreSort(uint treeID, uint xvar);
void freePreSortIntra(uint treeID);
void freePreSortExtra(void);
DistributionObj *stackRandomCovariatesSimple(uint treeID, Node *parent);
void unstackRandomCovariatesSimple(uint treeID, DistributionObj *obj);
char selectRandomCovariatesSimpleSingle(uint  treeID,
                                        Node *parent,
                                        DistributionObj *distributionObj,
                                        char *factorFlag,
                                        uint *covariate,
                                        uint *covariateCount);
char selectRandomCovariatesSimpleVector(uint  treeID,
                                        Node *parent,
                                        DistributionObj *distributionObj,
                                        char *factorFlag,
                                        uint *covariate,
                                        uint *covariateCount);
uint stackAndConstructSplitVectorSimple (uint     treeID,
                                         Node    *parent,
                                         uint     covariate,
                                         ...);
char forkNodeSimple(Node *parent, SplitInfo *info);
char forkAndUpdateSimple(uint       treeID,
                         Node      *parent,
                         uint      *repMembrIndx,
                         uint       repMembrSize,
                         uint      *allMembrIndx,
                         uint       allMembrSize,
                         char       multImpFlag,
                         SplitInfo *info,
                         uint      *leafCount,
                         Node     **nodeMembership);
char updateMaximumSplitSimple(uint    treeID,
                              Node   *parent,
                              double  delta,
                              uint    covariate,
                              uint    index,
                              char    factorFlag,
                              uint    mwcpSizeAbsolute,
                              uint    repMembrSize,
                              char  **polarity,
                              void   *splitVectorPtr,
                              SplitInfoMax *splitInfoMax);
char growTreeRecursiveSimple (uint     r,
                              char     rootFlag,
                              char     multImpFlag,
                              uint     treeID,
                              Node    *parent,
                              uint    *bootMembrIndxIter,
                              uint    *rmbrIterator,
                              uint    *ombrIterator);
char growTreeNonRecursiveSimple (uint     r,
                                 char     rootFlag,
                                 char     multImpFlag,
                                 uint     treeID,
                                 Node    *root,
                                 uint    *bootMembrIndxIter,
                                 uint    *rmbrIterator,
                                 uint    *ombrIterator);
void stackAndGetSplitSurv(uint    treeID,
                          Node    *parent,
                          char    eventType,
                          uint  **eventTimeCount,
                          uint  **eventTimeIndex,
                          uint   *eventTimeSize,
                          uint  **parentEvent,
                          uint  **parentAtRisk,
                          uint  **leftEvent,
                          uint  **leftAtRisk,
                          uint  **rightEvent,
                          uint  **rightAtRisk);
void unstackSplitSurv(uint    treeID,
                      Node    *parent,
                      uint *eventTimeCount,
                      uint *eventTimeIndex,
                      uint  eventTimeSize,
                      uint *parentEvent,
                      uint *parentAtRisk,
                      uint *leftEvent,
                      uint *leftAtRisk,
                      uint *rightEvent,
                      uint *rightAtRisk);
void stackSplitSurv3(uint    treeID,
                     Node    *parent,
                     uint   eventTimeSize,
                     double **leftLocalRatio,
                     double **rightLocalRatio,
                     double **leftLocalSurvival,
                     double **rightLocalSurvival,
                     uint   revEventTimeSize,
                     double **leftRevLocalRatio,
                     double **rightRevLocalRatio,
                     double **leftRevLocalSurvival,
                     double **rightRevLocalSurvival,
                     double **leftBS,
                     double **rightBS);
void unstackSplitSurv3(uint    treeID,
                       Node    *parent,
                       uint   eventTimeSize,
                       double *leftLocalRatio,
                       double *rightLocalRatio,
                       double *leftLocalSurvival,
                       double *rightLocalSurvival,
                       uint   revEventTimeSize,
                       double *leftRevLocalRatio,
                       double *rightRevLocalRatio,
                       double *leftRevLocalSurvival,
                       double *rightRevLocalSurvival,
                       double *leftBS,
                       double *rightBS);
uint getEventTime(uint   treeID,
                  Node   *parent,                  
                  uint   *repMembrIndx,
                  uint    repMembrSize,
                  uint   *nonMissMembrIndx,
                  uint    nonMissMembrSize,
                  char    eventType,
                  uint   *eventTimeCount,
                  uint   *eventTimeIndex);
void stackSplitEventAndRisk(uint    treeID,
                            Node    *parent,
                            uint    genEventTimeSize,
                            uint  **genParentEvent,
                            uint  **genParentAtRisk,
                            uint  **genLeftEvent,
                            uint  **genLeftAtRisk,
                            uint  **genRightEvent,
                            uint  **genRightAtRisk);
void unstackSplitEventAndRisk(uint    treeID,
                              Node    *parent,
                              uint    eventTimeSize,
                              uint   *genParentEvent,
                              uint   *genParentAtRisk,
                              uint   *genLeftEvent,
                              uint   *genLeftAtRisk,
                              uint   *genRightEvent,
                              uint   *genRightAtRisk);
void getSplitEventAndRisk(uint    treeID,
                          Node    *parent,
                          uint   *repMembrIndx,
                          uint    repMembrSize,
                          uint   *nonMissMembrIndx,
                          uint    nonMissMembrSize,
                          uint   *eventTimeCount,
                          uint   *eventTimeIndex,
                          uint    eventTimeSize,
                          uint   *parentEvent,
                          uint   *parentAtRisk);
void stackAndGetSplitSurv2(uint     treeID,
                           Node    *parent,
                           uint     eventTimeSize,
                           uint    *nodeParentEvent,
                           uint    *nodeParentAtRisk,
                           double **localSurvival);
void unstackAndGetSplitSurv2(uint     treeID,
                             Node    *parent,
                             uint     eventTimeSize,
                             double  *localSurvival);
void stackAndGetFZhat(uint  treeID,
                      Node *parent,
                      uint *repMembrIndx,
                      uint  repMembrSize,
                      uint *nonMissMembrIndx,
                      uint  nonMissMembrSize,
                      uint *eventTimeIndex,
                      uint  eventTimeSize,
                      uint *revEventTimeIndex,
                      uint  revEventTimeSize,
                      double *revParentSurvival,
                      double **fZHat);
void unstackFZhat(uint  treeID,
                  Node *parent,
                  uint  eventTimeSize,
                  double *fZHat);
void stackAndGetFZhatNew(uint  treeID,
                         Node *parent,
                         uint *repMembrIndx,
                         uint  repMembrSize,
                         uint *nonMissMembrIndx,
                         uint  nonMissMembrSize,
                         uint *eventTimeIndex,
                         uint  eventTimeSize,
                         uint *revEventTimeIndex,
                         uint  revEventTimeSize,
                         double *revParentSurvival,
                         double **gHat,
                         double **fZHat);
double getW_kt(uint  treeID,
               Node *parent,
               uint  indv,
               uint  tIndx,
               uint *eventTimeIndex,
               uint *revEventTimeIndex,
               uint  revEventTimeSize,
               double *revParentSurvival,
               double *gHatPrevious,
               double *gHatCurrent);
void stackAndGetQTime(uint  treeID,
                      Node *parent,
                      uint  eventTimeSize,
                      double *survival,
                      uint  **quantileTime);
void unstackQTime(uint  *quantileTime);
void stackAndGetQETime(uint   treeID,
                       Node  *parent,
                       uint  *eventTimeIndex,                       
                       uint   eventTimeSize,
                       double *survival,
                       uint  **qeTimeIndex,
                       uint   *qeTimeSize);
void unstackQETime(uint treeID,
                   uint eventTimeSize,
                   uint  *qeTimeIndex);
void stackAndGetLocalGamma(uint  treeID,
                             Node *parent,
                             uint *repMembrIndx,
                             uint  repMembrSize,
                             uint *nonMissMembrIndx,
                             uint  nonMissMembrSize,
                             uint *eventTimeIndex,
                             uint  eventTimeSize,
                             uint *revEventTimeIndex,
                             uint  revEventTimeSize,
                             double *revParentSurvival,
                           uint      *qeTimeIndex,
                           uint       qeTimeSize,
                           double   **gHat,
                             double  ***w_ktm,
                             double  ***gamma_ktm);
void  unstackLocalGamma(uint    treeID,
                        uint    nonMissMembrSize,
                        uint   *eventTimeIndex,
                        uint    eventTimeSize,
                           uint      *qeTimeIndex,
                        uint       qeTimeSize,
                        double **gamma_ktm);
uint stackAndConstructSplitVectorTDCPhase1 (uint     treeID,
                                            Node    *parent,
                                            uint     covariate,
                                            double  *splitVector,
                                            uint   **indxx);
void unstackSplitVectorTDC(uint   treeID,
                           Node  *parent,
                           uint   covariate,
                           uint   splitLength,
                           char   factorFlag,
                           uint   splitVectorSize,
                           uint   mwcpSizeAbsolute,
                           char   deterministicSplitFlag,
                           void  *splitVectorPtr,
                           char   multImpFlag,
                           uint  *indxx);
uint virtuallySplitNodeTDC(uint  treeID,
                           Node  *parent,
                           char  factorFlag,
                           uint  mwcpSizeAbsolute,
                           double *observation,
                           double **response,
                           uint *indxx,
                           void *splitVectorPtr,
                           uint  offset,
                           char *localSplitIndicator,
                           uint *leftSize,
                           uint *rghtSize,
                           uint  priorMembrIter,
                           uint *currentMembrIter);
void stackAndInitializeTimeAndSubjectArrays(char mode);
void unstackTimeAndSubjectArrays(char mode);
void stackFactorArrays(char mode);
void stackFactorGeneric(char    respFlag,
                        uint    size,
                        char   *type,
                        uint  **p_factorMap,
                        uint   *factorCount,
                        uint  **p_factorIndex,
                        uint  **p_factorSize,
                        uint  **p_nonfactorMap,
                        uint   *nonfactorCount,
                        uint  **p_nonfactorIndex);
void unstackFactorArrays(char mode);
void initializeFactorArrays(char mode);
char stackMissingArraysPhase1(char mode);
char stackMissingArraysPhase2(char mode);
void unstackMissingArrays(char mode);
void stackMissingSignatures(uint     obsSize,
                            uint     rspSize,
                            double **responsePtr,
                            double **predictorPtr,
                            uint    *recordMap,
                            uint     recordSize,
                            uint   **p_recordIndex,
                            uint    *p_vSize,
                            int   ***p_vSign,
                            int    **p_vIndex,
                            uint    *pRF_mrFactorSize,
                            uint   **pRF_mrFactorIndex,
                            uint    *pRF_mxFactorSize,
                            uint   **pRF_mxFactorIndex,
                            char    *pRF_mTimeFlag,
                            char    *pRF_mStatusFlag,
                            char    *pRF_mResponseFlag,
                            char    *pRF_mPredictorFlag);
void unstackMissingSignatures(uint      rspSize,
                              uint      recordSize,
                              uint     *recordIndex,
                              uint      vSize,
                              int     **vSign,
                              int      *vIndex,
                              uint      mrFactorSize,
                              uint     *mrFactorIndex,
                              uint      mxFactorSize,
                              uint     *mxFactorIndex);
char stackCompetingArrays(char mode);
void unstackCompetingArrays(char mode);
char stackClassificationArrays(char mode);
void unstackClassificationArrays(char mode);
void getEventInfo(char mode);
void stackDefinedOutputObjects(char      mode,
                               char    **sexpString,
                               Node   ***pRF_root,
                               uint    **pRF_tLeafCount_,
                               double  **pRF_proximity_,
                               double  **pRF_distance_,
                               double  **pRF_weight_,
                               double  **p_imputation_,
                               double ***pRF_sImputeResponsePtr,
                               double ***pRF_sImputePredictorPtr,
                               uint    **pRF_varUsed_,
                               uint   ***pRF_varUsedPtr,
                               double  **p_splitDepth_);
void unstackDefinedOutputObjects(char      mode);
void stackForestObjectsPtrOnly(char mode);
void stackTreeObjectsPtrOnly(char mode, uint treeID);
void stackForestObjectsOutput(char mode);
void writeForestObjectsOutput(char mode);
void unstackForestObjectsPtrOnly(char mode);
void unstackTreeObjectsPtrOnly(uint treeID);
void stackForestObjectsAuxOnly(char mode);
void unstackForestObjectsAuxOnly(char mode);
void unstackAuxStatisticalStructures(char mode);
void restackTermListAndQualitativeObjectsUnknown(uint treeID, uint length);
void verifyAndRegisterCustomSplitRules(void);
extern void registerCustomFunctions(void);
void stackAuxiliaryInfoList(SNPAuxiliaryInfo ***list, uint count);
void allocateAuxiliaryInfo(char   targetFlag,
                           char   type,
                           char  *stringIdentifier,
                           SNPAuxiliaryInfo **list,
                           uint   slot,
                           void  *snpPtr,
                           void  *auxiliaryArrayPtr,
                           uint   dimSize,
                           int   *dim);
uint getAuxDim(char flag, int *dim, uint preIndex, uint postIndex);
void unstackAuxiliaryInfoAndList(char targetFlag, SNPAuxiliaryInfo **list, uint count);
void memoryCheck(void);
void stackTNQualitativeObjectsKnown(char     mode,
                                    uint   **pRF_RMBR_ID_,
                                    uint   **pRF_AMBR_ID_,
                                    uint   **pRF_TN_RCNT_,
                                    uint   **pRF_TN_ACNT_,
                                    uint   **pRF_OOB_SZ_,
                                    uint   **pRF_IBG_SZ_);
void stackTNQualitativeObjectsUnknown(char     mode,
                                      uint   **pRF_TN_RCNT_,
                                      uint   **pRF_TN_ACNT_,
                                      uint   **pRF_TN_OCNT_,
                                      uint   **pRF_TN_ICNT_);
void stackTNQuantitativeForestObjectsPtrOnly(char mode);
void unstackTNQuantitativeForestObjectsPtrOnly(char mode);
void stackTNQuantitativeTreeObjectsPtrOnly(uint treeID);
void unstackTNQuantitativeTreeObjectsPtrOnly(uint treeID);
void saveTNQuantitativeTreeObjects(uint treeID);
void stackTNQuantitativeForestObjectsOutput(char mode);
void writeTNQuantitativeForestObjectsOutput(char mode);
void stackTNQualitativeObjectsUnknownMembership(char   mode, uint **pRF_OMBR_ID_, uint **pRF_IMBR_ID_);
void stackLocksOpenMP(char mode);
void unstackLocksOpenMP(char mode);
void stackLocksPosix(char mode);
void unstackLocksPosix(char mode);
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
void getAtRiskAndEventCount (uint       treeID,
                              Terminal  *parent,
                              uint      *repMembrIndx,
                              uint       repMembrSize,
                              uint      *allMembrIndx,
                              uint       allMembrSize,
                              uint      *rmbrIterator);
void getLocalRatio (uint treeID, Terminal *parent);
void getRevLocalRatio(uint treeID, Terminal *parent);
void getLocalCSH (uint treeID, Terminal *parent);
void getLocalCIF (uint treeID, Terminal *parent);
void mapLocalToTimeInterest(uint      treeID,
                            Terminal *parent,
                            void     *genericLocal,
                            void     *genericGlobal);
void getLocalSurvival (uint treeID, Terminal *parent);
void getLocalNelsonAalen (uint treeID, Terminal *parent);
void getSurvival (uint treeID, Terminal *parent);
void getMortality (uint treeID, Terminal *parent);
void getNelsonAalen (uint treeID, Terminal *parent);
void getCSH (uint treeID, Terminal *parent);
void getCIF (uint treeID, Terminal *parent);
void getLocalEmpiricalHazard(uint       treeID,
                             Terminal  *parent,
                             uint      *repMembrIndx,
                             uint       repMembrSize,
                             uint      *allMembrIndx,
                             uint       allMembrSize,
                             uint      *rmbrIterator);
void getLocalRatioTDC(uint treeID, Terminal *parent);
void getEmpiricalHazard(uint treeID, Terminal *parent);
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
void updateEnsembleHazard(char     mode,
                          uint     treeID,
                          char     normalizationFlag);
uint getTimeInterestIndex(double *array, uint length, double value);
void svdcmp(double **a, int m, int n, double ***uptr, double **wptr, double ***vptr);
char svdchk(double **a, uint m, uint n, double **u, double *w, double **v);
double **svdinv(double **u, double *w, double **v, uint m, uint n, uint singularity);
void free_svdcmp(double **a, int m, int n, double **u, double *w, double **v);
void svbksb(double **u, double *w, double **v, uint m, uint n, double *b, double *x);
double **matrixCopy(double **a, uint m, uint n);
double **matrixTrans(double **a, uint m, uint n);
double **matrixMult(double **a, double **b, uint m, uint n, uint p);
void matrixPrint(double **x, uint m, uint n);
double pythag(double a, double b);
void harness(void);
Terminal *makeTerminal(void);
void freeTerminal(Terminal *parent);
void stackTermLMIIndex(Terminal *tTerm, unsigned int size);
void unstackTermLMIIndex(Terminal *tTerm);
void freeTerminalNodeLocalSurvivalStructures(Terminal *tTerm);
void freeTerminalNodeSurvivalStructuresIntermediate(Terminal *tTerm);
void freeTerminalNodeSurvivalStructuresFinal(Terminal *tTerm);
void freeTerminalNodeLocalTDC(Terminal *tTerm);
void freeTerminalNodeTDC(Terminal *tTerm);
void freeTerminalNodeNonSurvivalStructures(Terminal *tTerm);
void stackAtRiskAndEventCount(Terminal *tTerm, unsigned int eTypeSize, unsigned int mTimeSize);
void unstackAtRiskAndEventCount(Terminal *tTerm);
void stackEventTimeIndex(Terminal *tTerm, unsigned int mTimeSize);
void unstackEventTimeIndex(Terminal *tTerm);
void stackEventTimeIndexHazard(Terminal *tTerm, unsigned int mTimeSize);
void unstackEventTimeIndexHazard(Terminal *tTerm);
void stackLocalRatio(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize);
void unstackLocalRatio(Terminal *tTerm);
void stackLocalRatioHazard(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize);
void unstackLocalRatioHazard(Terminal *tTerm);
void stackLocalSurvival(Terminal *tTerm, unsigned int eTimeSize);
void unstackLocalSurvival(Terminal *tTerm);
void stackLocalNelsonAalen(Terminal *tTerm, unsigned int eTimeSize);
void unstackLocalNelsonAalen(Terminal *tTerm);
void stackLocalCSH(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize);
void unstackLocalCSH(Terminal *tTerm);
void stackLocalCIF(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize);
void unstackLocalCIF(Terminal *tTerm);
void stackNelsonAalen(Terminal *tTerm, unsigned int sTimeSize);
void unstackNelsonAalen(Terminal *tTerm);
void stackSurvival(Terminal *tTerm, unsigned int sTimeSize);
void unstackSurvival(Terminal *tTerm);
void stackCSH(Terminal *tTerm, unsigned int eTypeSize, unsigned int sTimeSize);
void unstackCSH(Terminal *tTerm);
void stackCIF(Terminal *tTerm, unsigned int eTypeSize, unsigned int sTimeSize);
void unstackCIF(Terminal *tTerm);
void stackMortality(Terminal *tTerm, unsigned int eTypeSize);
void unstackMortality(Terminal *tTerm);
void stackMultiClassProb(Terminal *tTerm, unsigned int rfCount, unsigned int *rfSize);
void stackMultiClassProbPartial(Terminal *tTerm, unsigned int rfCount);
void unstackMultiClassProb(Terminal *tTerm);
void stackMeanResponse(Terminal *tTerm, unsigned int rnfCount);
void unstackMeanResponse(Terminal *tTerm);
void stackMemberStream(Terminal *tTerm, unsigned int membrSize);
void unstackMemberStream(Terminal *tTerm);
void stackLocalEmpiricalHazard(Terminal *tTerm, unsigned int eTimeSize);
void unstackLocalEmpiricalHazard(Terminal *tTerm);
void stackEmpiricalHazard(Terminal *tTerm, unsigned int sTimeSize);
void unstackEmpiricalHazard(Terminal *tTerm);
void acquireTreeGeneric(char mode, uint r, uint b);
void updateWeight(char mode, uint b);
void finalizeWeight(char mode);
void updateDistance(char mode, uint b);
void finalizeDistance(char mode);
void updateProximity(char mode, uint b);
void finalizeProximity(char mode);
void updateSplitDepth(uint treeID, Node *rootPtr, uint maxDepth);
char pruneBranch(uint obsSize, uint treeID, Node **nodesAtDepth, uint nadCount, uint ptnTarget, uint ptnCurrent);
uint pruneTree(uint obsSize, uint treeID, uint ptnCount);
void stackAuxiliary(char mode, uint b);
void unstackAuxiliary(char mode, uint b);
void printPseudoTNInfo(char mode, uint b);
void getPTNodeList(Node    *parent,
                   Node   **list,
                   uint    *offset);
void getSplitPath(uint treeID, Node *parent);
void freeSplitPath(uint treeID);
uint getMaximumDepth(Node *parent);
void getNodesAtDepth(Node *parent, uint tagDepth, Node **nodesAtDepth, uint *nadCount);
void acquireTreeJIT(char mode, uint r, uint treeID);
void restoreTerminalNodeJIT(uint treeID, Node *root, uint indv, double **xArray, Terminal **termMembership);
void getTerminalNodeJIT(uint treeID, Node *root, uint indv, double **xArray, Terminal **termMembership);
char growTreeRecursive(uint     r,
                       char     rootFlag,
                       char     multImpFlag,
                       uint     b,
                       Node    *parent,
                       uint    *bootMembrIndxIter,
                       uint    *rmbrIterator,
                       uint    *ambrIterator);
char growTreeLOT (uint     r,
                  char     multImpFlag,
                  uint     treeID,
                  Node    *root,
                  uint    *bootMembrIndxIter,                     
                  uint    *rmbrIterator,
                  uint    *ambrIterator);
void freeTree(uint treeID, Node *parent);
void saveStatistics(char     mode,
                    uint     b,
                    Node    *parent,
                    uint    *offset,
                    double  *spltST,
                    uint    *dpthST);
void initTerminalNodeMembership(uint       treeID,
                                Terminal  *parent,
                                uint      *allMembrIndx,
                                uint       allMembrSize);
void updatePruning(char mode, uint treeID);
void updateCaseDepth(char mode, uint treeID);
SEXP rfsrcCIndex(SEXP sexp_traceFlag,
                 SEXP sexp_size,
                 SEXP sexp_time,
                 SEXP sexp_censoring,
                 SEXP sexp_predicted,
                 SEXP sexp_denom);
SEXP rfsrcCIndexNew(SEXP sexp_traceFlag,
                    SEXP sexp_size,
                    SEXP sexp_time,
                    SEXP sexp_censoring,
                    SEXP sexp_predicted,
                    SEXP sexp_denom);
SEXP rfsrcTestSEXP(SEXP sexp_size);
SEXP rfsrcDistance(SEXP sexp_metric,
                   SEXP sexp_n,
                   SEXP sexp_p,
                   SEXP sexp_x,
                   SEXP sexp_sizeIJ,
                   SEXP sexp_rowI,
                   SEXP sexp_rowJ,
                   SEXP sexp_numThreads,
                   SEXP sexp_traceFlag);
double euclidean(uint n, uint p, uint i, uint j, double **x);
SEXP rfsrcGrow(SEXP traceFlag,
               SEXP seedPtr,
               SEXP optLow,
               SEXP optHigh,
               SEXP optSup,
               SEXP splitRule,
               SEXP nsplit,
               SEXP mtry,
               SEXP lot,
               SEXP baseLearn,
               SEXP vtry,
               SEXP vtryArray,
               SEXP vtryExperimental,
               SEXP ytry,
               SEXP nodeSize,
               SEXP nodeDepth,
               SEXP crWeightSize,
               SEXP crWeight,
               SEXP vimpThreshold,
               SEXP ntree,
               SEXP observationSize,
               SEXP yInfo,
               SEXP yLevels,
               SEXP yData,
               SEXP xInfo,
               SEXP xLevels,
               SEXP xData,
               SEXP sampleInfo,
               SEXP xWeightStat,
               SEXP yWeight,
               SEXP xWeight,
               SEXP timeInterest,
               SEXP nImpute,
               SEXP perfBlock,
               SEXP quantile,
               SEXP qStarPlus,
               SEXP xPreSort,
               SEXP numThreads);
SEXP rfsrcPredict(SEXP traceFlag,
                  SEXP seedPtr,
                  SEXP optLow,
                  SEXP optHigh,
                  SEXP vimpThreshold,
                  SEXP ntree,
                  SEXP observationSize,
                  SEXP yInfo,
                  SEXP yLevels,
                  SEXP yData,
                  SEXP xInfo,
                  SEXP xLevels,
                  SEXP xData,
                  SEXP sampleInfo,
                  SEXP timeInterestInfo,
                  SEXP totalNodeCount,
                  SEXP tLeafCount,
                  SEXP seedInfo,
                  SEXP hdim,
                  SEXP baseLearn,
                  SEXP treeID,
                  SEXP nodeID,
                  SEXP nodeSZ,
                  SEXP brnodeID,
                  SEXP hc_zero,
                  SEXP hc_oneAugIntr,
                  SEXP hc_oneAugSyth,
                  SEXP hc_one,
                  SEXP hc_parmID,
                  SEXP hc_contPT,
                  SEXP hc_contPTR,
                  SEXP hc_mwcpSZ,
                  SEXP hc_fsrecID,
                  SEXP hc_mwcpPT,
                  SEXP hc_augmXone,
                  SEXP hc_augmXtwo,
                  SEXP hc_augmXS,
                  SEXP hc_augmSythTop,
                  SEXP tnRMBR,
                  SEXP tnAMBR,
                  SEXP tnRCNT,
                  SEXP tnACNT,
                  SEXP tnSURV,
                  SEXP tnMORT,
                  SEXP tnNLSN,
                  SEXP tnCSHZ,
                  SEXP tnCIFN,
                  SEXP tnREGR,
                  SEXP tnCLAS,
                  SEXP rTargetInfo,
                  SEXP ptnCount,
                  SEXP xMarginalInfo,
                  SEXP intrPredictorInfo,
                  SEXP partial,
                  SEXP fobservationSize,
                  SEXP frSize,
                  SEXP frData,
                  SEXP fxData,
                  SEXP perfBlock,
                  SEXP quantileInfo,
                  SEXP getTree,
                  SEXP numThreads);
void exit2R(void);
void printR(char *format, ...);
void setNativeGlobalEnv(uint *nativeIndex, uint *stackCount);
void *copy1DObject(SEXP arr, char type, uint size, char actual);
void *copy2DObject(SEXP arr, char type, char flag, uint row, uint col);
void free_1DObject(void *arr, char type, uint size);
void free_2DObject(void *arr, char type, char flag, uint row, uint col);
void initProtect(uint  stackCount);
void *stackAndProtect(char   mode,
                      uint  *sexpIndex,
                      char   sexpType,
                      uint   sexpIdentity,
                      ulong  size,
                      double value,
                      char  *sexpString,
                      void  *auxiliaryPtr,
                      uint   auxiliaryDimSize,
                      ...);
void setUserTraceFlag (uint traceFlag);
uint getUserTraceFlag (void);
