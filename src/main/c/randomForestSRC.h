@RF_CRANX_BEG@
#include <R_ext/Print.h>
#include <Rdefines.h>
#define RF_nativePrint printR
#define RF_nativeError printR
#define RF_nativeExit  exit2R
#define RF_nativeNaN NA_REAL
#define RF_nativeIsNaN ISNA
@RF_CRANX_END@
@RF_JAVAX_BEG@
#include <jni.h>
#define RF_nativePrint printJ
#define RF_nativeError errorJ
#define RF_nativeExit  exit2J
#define RF_nativeNaN nan("")
#define RF_nativeIsNaN isnan
@RF_JAVAX_END@
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
#define RF_XXXX_19  19  
#define RF_XXXX_20  20  
#define RF_XXXX_21  21  
#define RF_XXXX_22  22  
#define RF_XXXX_23  23  
#define RF_XXXX_24  24  
#define RF_SEED_ID  25  
#define RF_VMP_SRG  26  
#define RF_VMP_CLS  27  
#define RF_VMP_RGR  28  
#define RF_MISS_ID  29  
#define RF_OMIS_ID  30  
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
#define RF_SPLT_VR  41  
#define RF_WGHT_ID  42  
#define RF_TN_SURV  43  
#define RF_TN_MORT  44  
#define RF_TN_NLSN  45  
#define RF_TN_CSHZ  46  
#define RF_TN_CIFN  47  
#define RF_TN_REGR  48  
#define RF_TN_CLAS  49  
#define RF_USPV_ST  50  
#define RF_MTRY_ID  51  
#define RF_MTRY_ST  52  
#define RF_XXXX_53  53  
#define RF_PART_SR  54  
#define RF_PART_CL  55  
#define RF_PART_RG   56  
#define RF_DIST_ID   57  
#define RF_XXXX_58   58  
#define RF_XXXX_59   59  
#define RF_TREE_ID   60  
#define RF_NODE_ID   61  
#define RF_HC_DIM    62  
#define RF_HC_PDIM   63  
#define RF_HC_PIDX   64  
#define RF_OS_PIDX   65  
#define RF_PARM_ID   66  
#define RF_CONT_PT   67  
#define RF_MWCP_SZ   68  
#define RF_MWCP_PT   69  
#define RF_MWCP_CT   70  
#define RF_PARM_ID2  71  
#define RF_CONT_PT2  72  
#define RF_MWCP_SZ2  73  
#define RF_MWCP_PT2  74  
#define RF_MWCP_CT2  75  
#define RF_PARM_ID3  76  
#define RF_CONT_PT3  77  
#define RF_MWCP_SZ3  78  
#define RF_MWCP_PT3  79  
#define RF_MWCP_CT3  80  
#define RF_PARM_ID4  81  
#define RF_CONT_PT4  82  
#define RF_MWCP_SZ4  83  
#define RF_MWCP_PT4  84  
#define RF_MWCP_CT4  85  
#define RF_SEXP_CNT  86  
#define OPT_FENS      0x00000001 
#define OPT_OENS      0x00000002 
#define OPT_PERF      0x00000004 
#define OPT_PERF_CALB 0x00000008 
#define OPT_LEAF      0x00000010 
#define OPT_TREE      0x00000020 
#define OPT_SEED      0x00000040 
#define OPT_MISS      0x00000080 
#define OPT_VIMP_TYP1 0x00000100 
#define OPT_VIMP_TYP2 0x00000200 
#define OPT_VIMP_JOIN 0x00000400 
#define OPT_USPV_STAT 0x00000800 
#define OPT_VARUSED_F 0x00001000 
#define OPT_VARUSED_T 0x00002000 
#define OPT_PERF_GMN2 0x00004000 
#define OPT_CLAS_RFQ  0x00008000 
#define OPT_IMPU_ONLY 0x00010000 
#define OPT_OUTC_TYPE 0x00020000 
#define OPT_BOOT_TYP1 0x00080000 
#define OPT_BOOT_TYP2 0x00100000 
#define OPT_COMP_RISK 0x00200000 
#define OPT_SPLDPTH_1 0x00400000 
#define OPT_SPLDPTH_2 0x00800000 
#define OPT_VIMP_LEOB 0x01000000 
#define OPT_VIMP      0x02000000 
#define OPT_NODE_STAT 0x08000000 
#define OPT_PROX      0x10000000 
#define OPT_PROX_IBG  0x20000000 
#define OPT_PROX_OOB  0x40000000 
#define OPT_PROX_FUL  0x60000000 
#define OPT_WGHT      0x00000001 
#define OPT_WGHT_IBG  0x00000002 
#define OPT_WGHT_OOB  0x00000004 
#define OPT_WGHT_FUL  0x00000006 
#define OPT_MISS_MIA  0x00000008 
#define OPT_MISS_SKIP 0x00000010 
#define OPT_MEMB_PRUN 0x00000020 
#define OPT_MEMB_USER 0x00000040 
#define OPT_MISS_MIAH 0x00000080 
#define OPT_SPLT_CUST 0x00000F00 
#define OPT_BOOT_SWOR 0x00001000 
#define OPT_TREE_ERR  0x00002000 
#define OPT_PART_PLOT 0x00004000 
#define OPT_MEMB_OUTG 0x00010000 
#define OPT_MEMB_INCG 0x00020000 
#define OPT_TERM_OUTG 0x00040000 
#define OPT_TERM_INCG 0x00080000 
#define OPT_DIST      0x00100000 
#define OPT_DIST_IBG  0x00200000 
#define OPT_DIST_OOB  0x00400000 
#define OPT_DIST_FUL  0x00600000 
#define RF_PART_MORT 1
#define RF_PART_NLSN 2
#define RF_PART_SURV 3
#define RF_PART_YRLS 1
#define RF_PART_CIFN 2 
#define RF_PART_CHFN 3
#define ACTIVE    0x02
#define LEFT      0x01
#define RIGHT     0x00
#define NEITHER   0x03
#define EPSILON 1.0e-9
#define RF_GROW   0x01
#define RF_PRED   0x02
#define RF_REST   0x04
#define RF_PART   0x08
#define SURV_LGRNK   1
#define SURV_LRSCR   2
#define SURV_CR_LAU  3
#define SURV_CR_LOG  4
#define RAND_SPLIT   5
#define REGR_WT_NRM  6
#define REGR_WT_OFF  7
#define REGR_WT_HVY  8
#define CLAS_WT_NRM  9
#define CLAS_WT_OFF 10
#define CLAS_WT_HVY 11
#define USPV_SPLIT  12
#define MVRG_SPLIT  13 
#define MVCL_SPLIT  14 
#define MVMX_SPLIT  15 
#define CUST_SPLIT  16
#define SURV_L2IMP  17
#define CLAS_RNKPS  18
#define MAXM_SPLIT  19 
#define SPLIT_MIA_NONE  0  
#define SPLIT_MIA_LEFT  1  
#define SPLIT_MIA_RGHT  2  
#define SPLIT_MIA_UNIT  3  
#define CLAS_FAM     0
#define REGR_FAM     1
#define SURV_FAM     2
#define CRSK_FAM     3
#define APROX 0
#define EXACT 1
#define MAX_EXACT_LEVEL sizeof(uint) * 8
#define SAFE_FACTOR_SIZE 16
#define MAX_HYPER_SIZE 16
#define MARGINAL_SIZE 8
#define RF_WGHT_UNIFORM 1
#define RF_WGHT_INTEGER 2
#define RF_WGHT_GENERIC 3
#define RF_DISTANCE_EUCLIDEAN 1
#define NATIVE_TYPE_CHARACTER 0
#define NATIVE_TYPE_INTEGER   1
#define NATIVE_TYPE_NUMERIC   2
typedef struct node Node;
struct node {
  struct node *parent;
  unsigned int xSize;
  char *permissibleSplit;
  struct terminal *mate;
  struct node *left;
  struct node *right;
  char splitFlag;
  unsigned int nodeID;
  unsigned int splitParameter;
  double splitValueCont;
  unsigned int splitValueFactSize;
  unsigned int *splitValueFactPtr;
  double splitStatistic;
  double mean;
  double variance;
  unsigned int *urStat;
  unsigned int  urStatSize;
  unsigned int  mtrySize;
  unsigned int *mtryIndx;
  double       *mtryStat;
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
};
typedef struct virtualSplitInfo VirtualSplitInfo;
struct virtualSplitInfo {
  uint size;
  char *indicator;
  double delta;
  uint *mwcpSizeAbs;
  uint *randomVar;
  double **splitArray;
  uint    *splitArrayLength;
  void **randomPts;
  uint randomIdx;
  uint hcDim;
  uint hcPartDim;
  uint hcPartIdx;
  uint osPartIdx;
};
typedef struct splitInfo SplitInfo;
struct splitInfo {
  uint     size;
  char    *indicator;
  double   delta;
  uint     hcDim;
  uint    *mwcpSizeAbs;
  uint    *randomVar;
  void   **randomPts;
  uint     hcPartDim;
  uint     hcPartIdx;
  uint     osPartIdx;
};
typedef struct terminal Terminal;
struct terminal {
  unsigned int nodeID;
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
};
char bootstrap (char     mode,
                uint     treeID,
                Node    *nodePtr,
                uint    *subIndex,
                uint     subsetSize,
                uint    *index,
                uint     indexSize);
char getNodeSign (char mode, uint treeID, Node *nodePtr, uint *bmIndex, uint repMembrSize);
void getMultiClassProb (uint       treeID,
                           Terminal  *parent,
                           uint      *repMembrIndx,
                           uint       repMembrSize,
                           uint      *allMembrIndx,
                           uint       allMembrSize,
                           uint      *rmbrIterator);
void updateEnsembleMultiClass(char     mode,
                              uint     treeID,
                              uint     serialTreeID,
                              char     omitDenominator);
double getBrierScore(uint     obsSize,
                     uint     rTarget,
                     double  *responsePtr,
                     double **outcomeCLS,
                     uint    *denomCount,
                     double  *cpv);
void getConditionalClassificationIndex(uint     size,
                                       uint     rTarget,
                                       double  *responsePtr,
                                       double **outcomeCLS,
                                       double  *maxVote,
                                       uint    *denomCount,
                                       double  *cpv);
double getClassificationIndex(uint     size,
                              uint     rTarget,
                              double  *responsePtr,
                              uint    *denomCount,
                              double  *maxVote);
double getGMeanIndex(uint    size,
                     uint    rTarget,
                     double *responsePtr,
                     uint   *denomCount,
                     double *maxVote);
void getMaxVote(uint     size,
                uint     rTarget,
                double **outcomeCLS,
                uint    *denomCount,
                double  *maxVote);
void processDefaultGrow();
void processDefaultPredict();
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
typedef struct hyperCube HyperCube;
struct hyperCube {
  unsigned int r; 
  unsigned int *cSize;
  unsigned int ***partition;
};
typedef struct orthoSlice OrthoSlice;
struct orthoSlice {
  unsigned int s;
  unsigned int osCount;
  unsigned int *polarity; 
};
HyperCube *makeHyperCube(uint r);
void bookPartition(uint s, uint index, uint *row, uint *slice, HyperCube *h);
void freeHyperCube(HyperCube *h);
void unbookHyperCube(HyperCube *h);
OrthoSlice *makeOrthoSlice(uint s);
void freeOrthoSlice(OrthoSlice *os);
Node *identifyPerturbedMembership(Node    *parent,
                                  double **shadowVIMP,
                                  uint     index);
Node *randomizeMembership(Node    *parent,
                          double **predictor,
                          uint     individual,
                          uint     splitParameter,
                          uint     treeID);
Node *antiMembership(Node    *parent,
                     double **predictor,
                     uint     individual,
                     uint     splitParameter,
                     uint     treeID);
void permute(uint ranGenID, uint p, uint n, uint *indx);
void getAntiMembership(char       mode,
                       uint       treeID,
                       Terminal **vimpMembership,
                       uint       p);
void getRandomMembership(char       mode,
                         uint       treeID,
                         Terminal **vimpMembership,
                         uint       p);
void getPermuteMembership(char       mode,
                          uint       treeID,
                          Terminal **vimpMembership,
                          uint       p);
void getVimpMembership(char      mode,
                       uint      treeID,
                       Terminal    **vimpMembership,
                       uint      p);
void updateGenericVimpEnsemble (char       mode,
                                uint       treeID,
                                uint       xVarIdx,
                                uint      *treeDenom,
                                Terminal **noiseMembership,
                                char       ensembleFlag,
                                double  ***genEnsembleMRT,
                                double ****genEnsembleCLS,
                                double  ***genEnsembleRGR);
void updateTreeEnsemble (char       mode,
                         uint       treeID,
                         uint      *treeDenom,
                         double  ***treeEnsembleMRT,
                         double ****treeEnsembleCLS,
                         double  ***treeEnsembleRGR);
void updateVimpEnsemble (char       mode,
                         uint       treeID,
                         Terminal **vimpMembership,
                         uint       p);
void summarizeVimpPerformance(char       mode,
                              uint       treeID,
                              uint       p);
void finalizeVimpPerformance(char mode, uint rejectedTreeCount);
void  stackVimpMembership(char mode, Terminal ***membership);
void  unstackVimpMembership(char mode, Terminal **membership);
void stackTreeEnsemble(char         mode,
                       uint         treeID,
                       uint       **denomTree,
                       double   ****treeOutcome,
                       double  *****sTreeOutcome,
                       double   ****mcTreeOutcome);
void unstackTreeEnsemble(char       mode,
                         uint       treeID,
                         uint      *denomTree,
                         double  ***treeOutcome,
                         double ****sTreeOutcome,
                         double  ***mcTreeOutcome);
void updateVimpCalculations (char mode, uint b, uint intrIndex, Terminal **vimpMembership);
void summarizeTreePerformance(char mode, uint treeID);
void updatePartialCalculations (uint       treeID,
                                uint       pVarIdx,
                                Terminal **partialMembership);
void summarizePartialCalculations(uint       treeID,
                                  uint       pVarIdx);
char imputeNode (char     type,
                 char     termFlag,
                 char     chainFlag,
                 uint     treeID, 
                 Node    *nodePtr,
                 uint    *repAbsIdx,
                 uint     repNodeSize,
                 uint    *iAbsIdx,
                 uint     iNodeSize);
char restoreNodeMembership(uint  r,
                           char  mode, 
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
void imputeUpdateSummary (char     mode, 
                          uint     treeID);
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
                    uint      treeID,
                    double  **tempResponse);
void imputeCommon(char      mode,
                  uint      treeID,
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
                          double *time, 
                          char    naflag,
                          char    idFlag,
                          uint   *masterTimeIndex);
void updateEventTypeSubsets(double *summaryStatus, 
                            uint    mRecordSize,
                            int   **mpSign,
                            uint   *mRecordIndex,
                            uint   *meIndividualSize,
                            uint  **eIndividual);
void stackShadow (char mode, uint treeID);
void unstackShadow (char mode, uint treeID, char respFlag, char covrFlag);
char xferMissingness(char type, Node *source, Terminal *destination);
char getMarginalNodeMembership(char     mode,
                               char     rootFlag,
                               uint     treeID,
                               Node    *parent,
                               uint    *allMembrIndx,
                               uint     allMembrSize,
                               double **observationPtr);
char getPartialNodeMembership(char       rootFlag,
                              uint       treeID,
                              uint       partialIndex,
                              Node      *parent,
                              uint      *allMembrIndx,
                              uint       allMembrSize,
                              double   **observationPtr,
                              Terminal **membership);
Node *makeNode(unsigned int xSize, unsigned int permissibleSizeAlloc, unsigned int urStatSize, unsigned int mtrySize);
void freeNode(Node *parent);
void getNodeInfo(Node *leaf);
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
  NRUTIL_HPTR,   
  NRUTIL_OPTR,   
  NRUTIL_VPTR    
};
unsigned int upower (unsigned int x, unsigned int n);
unsigned int upower2 (unsigned int n);
unsigned int ulog2 (unsigned int n);
void hpsort(double *ra, unsigned int n);
void hpsortui(unsigned int *ra, unsigned int n);
void sort(double *arr, unsigned int n);
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
void testEndianness();
void randomStack(uint bSize, uint xSize);
void randomUnstack(uint bSize, uint xSize);
void randomSetChainParallel(uint b, int value);
void randomSetUChainParallel(uint b, int value);
void randomSetUChainParallelCov(uint b, int value);
void randomSetChainSerial(uint b, int value);
void randomSetUChainSerial(uint b, int value);
void randomSetUChainSerialCov(uint b, int value);
int randomGetChainParallel(uint b);
int randomGetUChainParallel(uint b);
int randomGetUChainParallelCov(uint b);
int randomGetChainSerial(uint b);
int randomGetUChainSerial(uint b);
int randomGetUChainSerialCov(uint b);
float randomChainParallel(uint b);
float randomUChainParallel(uint b);
float randomUChainParallelCov(uint b);
float randomChainSerial(uint b);
float randomUChainSerial(uint b);
float randomUChainSerialCov(uint b);
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
                        uint     serialTreeID,
                        char     omitDenominator);
double getMeanSquareError(uint    size,
                          double *responsePtr,
                          double *predictedOutcome,
                          uint   *oobCount);
char getVariance(uint    repMembrSize,
                 uint   *repMembrIndx,
                 uint    nonMissMembrSize,
                 uint   *nonMissMembrIndx,
                 double *targetResponse,
                 double *mean,
                 double *variance);
void rfsrc(char mode, int seedValue);
void stackGreedyDataStructures();
void unstackGreedyDataStructures();
void updateTerminalNodeOutcomes(char       mode,
                                uint       treeID,
                                Terminal  *parent,
                                uint      *repMembrIndx,
                                uint       repMembrSize,
                                uint      *allMembrIndx,
                                uint       allMembrSize,
                                uint      *rbmrIterator,
                                uint      *ambrIterator);
void updateEnsembleCalculations (char      multipleImputeFlag,
                                 char      mode,
                                 uint      b);
char stackAndImputePerfResponse(char      mode,
                                char      multipleImputeFlag,
                                uint      treeID,
                                uint      serialID,
                                double ***responsePtr);
double **stackAndImputeGenericResponse(char flag, char mode, uint obsSize, uint treeID, uint serialID, double **responsePtr);
void unstackImputeResponse(char flag, uint obsSize, double **mResponsePtr);
void getPerformance(uint      serialTreeID,
                    char      mode,
                    uint      obsSize,
                    double  **responsePtr,
                    uint      *denomPtr,
                    double   **outcomeMRT,
                    double  ***outcomeCLS,
                    double   **outcomeRGR,
                    double   *perfMRTptr,
                    double  **perfCLSptr,
                    double   *perfRGRptr);
void finalizeEnsembleEstimates(char mode);
char getPerformanceFlag (char mode, uint serialTreeID);
void getVariablesUsed(uint treeID, Node *rootPtr, uint *varUsedVector);
void initializeCDF(uint     treeID,
                   uint    *permissibilityIndex,
                   char    *permissibilityFlag,
                   uint     permissibilitySize,
                   uint     weightType,
                   double  *weight,
                   uint    *weightSorted,
                   uint     maxDensitySize,
                   uint   **index,
                   uint    *sampleSize,
                   double **cdf,
                   uint    *cdfSize,
                   uint   **cdfSort,
                   uint   **density,
                   uint    *densitySize,
                   uint  ***densitySwap);
void updateCDF(uint    treeID,
               uint    weightType,
               double *weight,
               uint   *index,
               uint   *sampleSize,
               uint    sampleSlot,
               double *cdf,
               uint    cdfSize,
               uint   *density,
               uint   *densitySize,
               uint  **densitySwap,
               uint    absoluteSlot);
uint sampleFromCDF (float (*genericGenerator) (uint),
                    uint    treeID,
                    uint    weightType,
                    uint   *sampleIndex,
                    uint    sampleSize,
                    uint   *sampleSlot,
                    double *cdf,
                    uint    cdfSize,
                    uint   *cdfSort,
                    uint   *density,
                    uint    densitySize);
void discardCDF(uint     treeID,
                uint     weightType,
                double  *weight,
                uint     permissibilitySize,
                uint    *index,
                uint     indexSize,
                uint    *density,
                uint     densitySize,
                uint   **densitySwap,
                double  *cdf,
                uint     cdfSize,
                uint    *cdfSort);
void duplicateCDF(uint treeID,
                  uint  type,
                  double  *weight,
                  uint     weightSize,
                  uint    *inc_index,
                  uint     inc_sampleSize,
                  double  *inc_cdf,
                  uint     inc_cdfSize,
                  uint    *inc_cdfSort,
                  uint    *inc_density,
                  uint     inc_densitySize,
                  uint   **inc_densitySwap,
                  uint   **index,
                  uint    *sampleSize,
                  double **cdf,
                  uint    *cdfSize,
                  uint   **cdfSort,
                  uint   **density,
                  uint    *densitySize,
                  uint  ***densitySwap);
uint sampleUniformlyFromVector (uint    treeID,
                                uint   *index,
                                uint    size,
                                uint   *sampleSlot);
typedef struct snpAuxiliaryInfo SNPAuxiliaryInfo;
struct snpAuxiliaryInfo {
  char type;
  uint identity;
  ulong linearSize;
  void *snpPtr;
  void *auxiliaryArrayPtr;
  uint dimSize;
  int *dim;
};
char getBestSplit(uint    treeID,
                  Node   *parent,
                  uint   *repMembrIndx,
                  uint    repMembrSize,
                  uint   *allMembrIndx,
                  uint    allMembrSize,
                  uint   *splitParameterMax,
                  double *splitValueMaxCont,
                  uint   *splitValueMaxFactSize,
                  uint  **splitValueMaxFactPtr,
                  double *splitStatistic,
                  char  **splitIndicator,
                  char   *splitMIA,
                  char    multImpFlag);
char randomSplit(uint    treeID,
                 Node   *parent,
                 uint   *repMembrIndx,
                 uint    repMembrSize,
                 uint   *allMembrIndx,
                 uint    allMembrSize,
                 uint   *splitParameterMax,
                 double *splitValueMaxCont,
                 uint   *splitValueMaxFactSize,
                 uint  **splitValueMaxFactPtr,
                 double *splitStatistic,
                 char  **splitIndicator,
                 char   *splitMIA,
                 char    multImpFlag);
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
char getBestSplitGreedy(uint         treeID,
                        Node        *parent,
                        uint        *repMembrIndx,
                        uint         repMembrSize,
                        uint        *allMembrIndx,
                        uint         allMembrSize,
                        char         multImpFlag,
                        SplitInfo  **splitInfo);
char classificationXwghtSplit(uint    treeID,
                              Node   *parent,
                              uint   *repMembrIndx,
                              uint    repMembrSize,
                              uint   *allMembrIndx,
                              uint    allMembrSize,
                              uint   *splitParameterMax,
                              double *splitValueMaxCont,
                              uint   *splitValueMaxFactSize,
                              uint  **splitValueMaxFactPtr,
                              double *splitStatistic,
                              char  **splitIndicator,
                              char   *splitMIA,
                              char    multImpFlag);
char classificationRankedPS(uint    treeID,
                            Node   *parent,
                            uint   *repMembrIndx,
                            uint    repMembrSize,
                            uint   *allMembrIndx,
                            uint    allMembrSize,
                            uint   *splitParameterMax,
                            double *splitValueMaxCont,
                            uint   *splitValueMaxFactSize,
                            uint  **splitValueMaxFactPtr,
                            double *splitStatistic,
                            char  **splitIndicator,
                            char   *splitMIA,
                            char    multImpFlag);
char regressionXwghtSplit(uint    treeID,
                          Node   *parent,
                          uint   *repMembrIndx,
                          uint    repMembrSize,
                          uint   *allMembrIndx,
                          uint    allMembrSize,
                          uint   *splitParameterMax,
                          double *splitValueMaxCont,
                          uint   *splitValueMaxFactSize,
                          uint  **splitValueMaxFactPtr,
                          double *splitStatistic,
                          char  **splitIndicator,
                          char   *splitMIA,
                          char    multImpFlag);
char regressionXwghtSplitGreedy (uint         treeID,
                                 Node        *parent,
                                 uint        *repMembrIndx,
                                 uint         repMembrSize,
                                 uint        *allMembrIndx,
                                 uint         allMembrSize,
                                 char         multImpFlag,
                                 SplitInfo  **splitInfo);
char selectRandomCovariatesGreedy(uint      treeID,
                                  Node     *parent,
                                  uint     *repMembrIndx,
                                  uint      nonMissMembrSizeStatic,
                                  uint     *nonMissMembrIndxStatic,
                                  char      multImpFlag,
                                  char     *cdfInitFlag,
                                  uint    **covariateIndex,
                                  uint     *uniformSize,
                                  double  **cdf,
                                  uint     *cdfSize,
                                  uint    **cdfSort,
                                  uint    **density,
                                  uint     *densitySize,
                                  uint   ***densitySwap,
                                  uint     *hCubeCount,
                                  VirtualSplitInfo *info,
                                  uint     *nonMissMembrSize,
                                  uint    **nonMissMembrIndx);
void unselectRandomCovariatesGreedy(uint      treeID,
                                    char      multImpFlag,
                                    VirtualSplitInfo *info,
                                    uint      nonMissMembrSizeStatic,
                                    uint     *nonMissMembrIndx);
uint stackAndConstructSplitVectorGreedy (uint     treeID,
                                         VirtualSplitInfo *info);
void unstackSplitVectorGreedy(uint   treeID,
                              VirtualSplitInfo *info);
uint virtuallySplitNodeGreedy(uint   treeID,
                              uint  *repMembrIndx,
                              uint  *nonMissMembrIndx,
                              uint   nonMissMembrSize,
                              uint   hcPartDim,
                              uint   hcPartIdx,
                              uint   osPartIdx,
                              uint   randomIdx, 
                              VirtualSplitInfo *info,
                              uint  *leftSize);
char updateMaximumSplitGreedy(uint    treeID,
                              Node   *parent,
                              VirtualSplitInfo *info,
                              SplitInfo       **maxInfo);
char summarizeSplitResultGreedy(SplitInfo *info);
SplitInfo *makeSplitInfo(unsigned int indicatorSize);
void freeSplitInfo(SplitInfo *info);
SplitInfo *makeSplitInfoPartial(SplitInfo *info, uint hcDim);
void freeSplitInfoPartial(SplitInfo *info);
VirtualSplitInfo *makeVirtualSplitInfo();
void freeVirtualSplitInfo(VirtualSplitInfo *info);
void getSplitObjectInfo(SplitInfo *info);
void getVirtualSplitObjectInfo(VirtualSplitInfo *info);
char forkAndUpdate(uint       treeID,
                   Node      *parent,
                   uint      *repMembrIndx,
                   uint       repMembrSize,
                   uint      *allMembrIndx,
                   uint       allMembrSize,
                   char       multImpFlag,
                   SplitInfo *info,
                   char      *indicator,
                   uint      *leftDaughterSize,
                   uint      *rghtDaughterSize);
char forkNode(Node      *parent,
              SplitInfo *info);
void saveTree(uint    b,
              Node   *parent,
              ulong  *offset);
void restoreTree(char    mode,
                 uint    b,
                 Node   *parent);
char logRankNCR(uint    treeID,
                Node   *parent,
                uint   *repMembrIndx,
                uint    repMembrSize,
                uint   *allMembrIndx,
                uint    allMembrSize,
                uint   *splitParameterMax,
                double *splitValueMaxCont,
                uint   *splitValueMaxFactSize,
                uint  **splitValueMaxFactPtr,
                double *splitStatistic,
                char  **splitIndicator,
                char   *splitMIA,
                char    multImpFlag);
char logRankCR(uint    treeID,
               Node   *parent,
               uint   *repMembrIndx,
               uint    repMembrSize,
               uint   *allMembrIndx,
               uint    allMembrSize,
               uint   *splitParameterMax,
               double *splitValueMaxCont,
               uint   *splitValueMaxFactSize,
               uint  **splitValueMaxFactPtr,
               double *splitStatistic,
               char  **splitIndicator,
               char   *splitMIA,
               char    multImpFlag);
char l2Impute(uint    treeID,
              Node   *parent,
              uint   *repMembrIndx,
              uint    repMembrSize,
              uint   *allMembrIndx,
              uint    allMembrSize,
              uint   *splitParameterMax,
              double *splitValueMaxCont,
              uint   *splitValueMaxFactSize,
              uint  **splitValueMaxFactPtr,
              double *splitStatistic,
              char  **splitIndicator,
              char   *splitMIA,
              char    multImpFlag);
void getMembrCountOnly (uint       treeID,
                        Terminal  *parent,
                        uint      *repMembrIndx,
                        uint       repMembrSize,
                        uint      *allMembrIndx,
                        uint       allMembrSize,
                        uint      *rmbrIterator);
char unsupervisedSplit(uint    treeID,
                       Node   *parent,
                       uint   *repMembrIndx,
                       uint    repMembrSize,
                       uint   *allMembrIndx,
                       uint    allMembrSize,
                       uint   *splitParameterMax,
                       double *splitValueMaxCont,
                       uint   *splitValueMaxFactSize,
                       uint  **splitValueMaxFactPtr,
                       double *splitStatistic,
                       char  **splitIndicator,
                       char   *splitMIA,
                       char    multImpFlag);
char multivariateSplit (uint    treeID,
                        Node   *parent,
                        uint   *repMembrIndx,
                        uint    repMembrSize,
                        uint   *allMembrIndx,
                        uint    allMembrSize,
                        uint   *splitParameterMax,
                        double *splitValueMaxCont,
                        uint   *splitValueMaxFactSize,
                        uint  **splitValueMaxFactPtr,
                        double *splitStatistic,
                        char  **splitIndicator,
                        char   *splitMIA,
                        char    multImpFlag);
char customMultivariateSplit (uint    treeID,
                              Node   *parent,
                              uint   *repMembrIndx,
                              uint    repMembrSize,
                              uint   *allMembrIndx,
                              uint    allMembrSize,
                              uint   *splitParameterMax,
                              double *splitValueMaxCont,
                              uint   *splitValueMaxFactSize,
                              uint  **splitValueMaxFactPtr,
                              double *splitStatistic,
                              char  **splitIndicator,
                              char   *splitMIA,
                              char    multImpFlag);
char customSurvivalSplit (uint    treeID,
                         Node   *parent,
                         uint   *repMembrIndx,
                         uint    repMembrSize,
                         uint   *allMembrIndx,
                         uint    allMembrSize,
                         uint   *splitParameterMax,
                         double *splitValueMaxCont,
                         uint   *splitValueMaxFactSize,
                         uint  **splitValueMaxFactPtr,
                         double *splitStatistic,
                         char  **splitIndicator,
                          char   *splitMIA,
                         char    multImpFlag);
char customCompetingRiskSplit (uint    treeID,
                               Node   *parent,
                               uint   *repMembrIndx,
                               uint    repMembrSize,
                               uint   *allMembrIndx,
                               uint    allMembrSize,
                               uint   *splitParameterMax,
                               double *splitValueMaxCont,
                               uint   *splitValueMaxFactSize,
                               uint  **splitValueMaxFactPtr,
                               double *splitStatistic,
                               char  **splitIndicator,
                               char   *splitMIA,
                               char    multImpFlag);
void stackSplitIndicator(uint     nodeSize,
                         char   **localSplitIndicator,
                         double **splitVector);
void unstackSplitIndicator(uint    nodeSize,
                           char   *localSplitIndicator,
                           double *splitvector);
void stackSplitEventTime(uint **localEventTimeCount,
                         uint **localEventTimeIndex);
void unstackSplitEventTime(uint *localEventTimeCount,
                           uint *localEventTimeIndex);
uint getSplitEventTime(uint   treeID,
                       uint   *repMembrIndx,
                       uint    repMembrSize,
                       uint   *nonMissMembrIndx,
                       uint    nonMissMembrSize,
                       uint   *localEventTimeCount,
                       uint   *localEventTimeIndex);
void stackSplitEventAndRisk(uint   eventTimeSize,
                            uint **nodeParentEvent,
                            uint **nodeParentAtRisk,
                            uint **nodeLeftEvent,
                            uint **nodeLeftAtRisk,
                            uint **nodeRightEvent,
                            uint **nodeRightAtRisk);
void unstackSplitEventAndRisk(uint  eventTimeSize,
                              uint *nodeParentEvent,
                              uint *nodeParentAtRisk,
                              uint *nodeLeftEvent,
                              uint *nodeLeftAtRisk,
                              uint *nodeRightEvent,
                              uint *nodeRightAtRisk);
void getSplitEventAndRisk(uint    treeID,
                          uint   *repMembrIndx,
                          uint    repMembrSize,
                          uint   *nonMissMembrIndx,
                          uint    nonMissMembrSize,
                          uint   *localEventTimeCount,
                          uint   *localEventTimeIndex,
                          uint    localEventTimeSize,
                          uint   *nodeParentEvent,
                          uint   *nodeParentAtRisk);
void stackAndGetSplitSurv(uint    treeID,
                          uint   *repMembrIndx,
                          uint    repMembrSize,
                          uint   *nonMissMembrIndx,
                          uint    nonMissMembrSize,
                          uint  **localEventTimeCount,
                          uint  **localEventTimeIndex,
                          uint   *localEventTimeSize,
                          uint  **nodeParentEvent,
                          uint  **nodeParentAtRisk,
                          uint  **nodeLeftEvent,
                          uint  **nodeLeftAtRisk,
                          uint  **nodeRightEvent,
                          uint  **nodeRightAtRisk);
void unstackSplitSurv(uint *localEventTimeCount,
                      uint *localEventTimeIndex,
                      uint  eventTimeSize,
                      uint *nodeParentEvent,
                      uint *nodeParentAtRisk,
                      uint *nodeLeftEvent,
                      uint *nodeLeftAtRisk,
                      uint *nodeRightEvent,
                      uint *nodeRightAtRisk);
void stackAndGetSplitSurvL2(uint     treeID,
                            Node    *parent,
                            uint     localEventTimeSize,
                            uint    *localEventTimeIndex,
                            uint    *nodeParentEvent,
                            uint    *nodeParentAtRisk,
                            double **localRatio,
                            double **localSurvival);
void unstackAndGetSplitSurvL2(uint     localEventTimeSize,
                              double  *localRatio,
                              double  *localSurvival);
void stackAndGetL2Impute(uint     treeID,
                         Node    *parent,
                         uint    *repMembrIndx,
                         uint     repMembrSize,
                         uint    *nonMissMembrIndx,
                         uint     nonMissMembrSize,
                         uint     localEventTimeSize,
                         double  *localSurvival,
                         double **l2Impute);
uint stackAndConstructSplitVector(uint     treeID,
                                  uint     localMembershipSize,
                                  uint     randomCovariateIndex,
                                  double  *splitVector,
                                  uint     splitVectorSize,
                                  char    *factorFlag,
                                  char    *deterministicSplitFlag,
                                  uint    *mwcpSizeAbsolute,
                                  void   **splitVectorPtr);
void unstackSplitVector(uint   treeID,
                        uint   splitVectorSize,
                        uint   splitLength,
                        char   factorFlag,
                        char   deterministicSplitFlag,
                        uint   mwcpSizeAbsolute,
                        void  *splitVectorPtr);
void stackRandomCovariates(uint      treeID,
                           Node     *parent,
                           uint      repMembrSize,
                           char      multImpFlag,
                           uint    **covariateIndex,
                           uint     *covariateSize,
                           double  **cdf,
                           uint     *cdfSize,
                           uint    **cdfSort,
                           uint    **density,
                           uint     *densitySize,
                           uint   ***densitySwap);
void unstackRandomCovariates(uint     treeID,
                             Node     *parent,
                             uint    *covariateIndex,
                             double  *cdf,
                             uint    *cdfSort,
                             uint    *density,
                             uint   **densitySwap);
char selectRandomCovariates(uint     treeID,
                            Node     *parent,
                            uint     *repMembrIndx,
                            uint      repMembrSize,
                            uint     *covariateIndex,
                            uint     *uniformSize,
                            uint     *uniformSelectedSlot,
                            double   *cdf,
                            uint     *cdfSize,
                            uint     *cdfSort,
                            uint     *density,
                            uint     *densitySize,
                            uint    **densitySwap,
                            uint     *covariate,
                            uint     *actualCovariateCount,
                            uint     *candidateCovariateCount,
                            double   *splitVector,
                            uint     *splitVectorSize,
                            uint    **indxx,
                            uint      nonMissMembrSizeStatic,
                            uint     *nonMissMembrIndxStatic,
                            uint     *nonMissMembrSize,
                            uint    **nonMissMembrIndx,
                            char      multImpFlag);
void unselectRandomCovariates(uint      treeID,
                              Node     *parent,
                              uint      repMembrSize,
                              uint     *indxx,
                              uint      nonMissMembrSizeStatic,
                              uint     *nonMissMembrIndx,
                              char      multImpFlag);
uint virtuallySplitNode(uint  treeID,
                           char  factorFlag,
                           uint  mwcpSizeAbsolute,
                           uint  randomCovariate,
                           uint *repMembrIndx,
                           uint  repMembrSize,
                           uint *nonMissMembrIndx,
                           uint  nonMissMembrSize,
                           uint *indxx,
                           void *splitVectorPtr,
                           uint  offset,
                           char *localSplitIndicator,
                           uint *leftSize,
                           uint  priorMembrIter,
                           uint *currentMembrIter);
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
char summarizeSplitResult(uint    splitParameterMax,
                          double  splitValueMaxCont,
                          uint    splitValueMaxFactSize,
                          uint   *splitValueMaxFactPtr,
                          double *splitStatistic,
                          double  deltaMax);
char getPreSplitResult (uint      treeID,
                        Node     *parent,
                        uint      repMembrSize,
                        uint     *repMembrIndx,
                        uint     *nonMissMembrSize,
                        uint    **nonMissMembrIndx,
                        double   *preSplitMean,
                        char      multImpFlag,
                        char      multVarFlag);
void unstackPreSplit (char      preliminaryResult,
                      char      multImpFlag,
                      char      multVarFlag,
                      uint      repMembrSize,
                      uint     *nonMissMembrIndx);
char updateMaximumSplit(uint    treeID,
                        Node   *parent,
                        double  delta,
                        uint    candidateCovariateCount,
                        uint    covariate,
                        uint    index,
                        char    factorFlag,
                        uint    mwcpSizeAbsolute,
                        uint    repMembrSize,
                        char   *localSplitIndicator,
                        double *deltaMax,
                        uint   *splitParameterMax,
                        double *splitValueMaxCont,
                        uint   *splitValueMaxFactSize,
                        uint  **splitValueMaxFactPtr,
                        void   *splitVectorPtr,
                        char  **splitIndicator);
void updateNodeStatistics(Node *parent, double delta, uint candidateCovariateCount, uint covariate);
void initializeTimeArrays(char mode);
void stackFactorArrays();
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
void unstackFactorArrays();
void initializeFactorArrays(char mode);
char stackMissingArrays(char mode);
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
void getEventTypeSize(uint     obsSize,
                      double  *status,
                      uint    *mRecordMap,
                      int    **mpSign,
                      uint    *eventTypeSize,
                      uint    *msize,
                      uint    *eventType);
void getClassLevelSize(uint      obsSize,
                       double  **response,
                       uint     *mRecordMap,
                       int     **mpSign,
                       uint     *classLevelSize,
                       uint    **classLevel);
void stackDefinedOutputObjects(char      mode,
                               char    **sexpString,
                               Node   ***pRF_root,
                               uint    **pRF_tLeafCount_,
                               double  **pRF_proximity_,
                                
                               double  **pRF_weight_,
                               int     **pRF_seed_,
                               double  **p_imputation_,
                               double ***pRF_sImputeResponsePtr,
                               double ***pRF_sImputePredictorPtr,
                               uint    **pRF_varUsed_,
                               uint   ***pRF_varUsedPtr,
                               double  **p_splitDepth_);
void unstackDefinedOutputObjects(char      mode,
                                 Node    **root);
void stackForestOutputObjects(char mode);
void stackAuxForestObjects(char mode);
void unstackAuxForestObjects(char mode);
void stackStatisticalOutputObjects(char     mode,
                                   double **pRF_spltST_,
                                   double **pRF_spltVR_,
                                   uint   **pRF_uspvST_,
                                   uint   **pRF_mtryID_,
                                   double **pRF_mtryST_);
void unstackAuxStatisticalStructures(char mode);
void stackTNQualitativeOutputObjects(char     mode,
                                     uint   **pRF_RMBR_ID_,
                                     uint   **pRF_AMBR_ID_,
                                     uint   **pRF_TN_RCNT_,
                                     uint   **pRF_TN_ACNT_);
void stackTNQuantitativeOutputObjects(char     mode,
                                      double **pRF_TN_SURV_,
                                      double **pRF_TN_MORT_,
                                      double **pRF_TN_NLSN_,
                                      double **pRF_TN_CSHZ_,
                                      double **pRF_TN_CIFN_,
                                      double **pRF_TN_REGR_,
                                      uint   **pRF_TN_CLAS_);
void verifyAndRegisterCustomSplitRules();
extern void registerCustomFunctions();
void stackAuxiliaryInfoList();
void allocateAuxiliaryInfo(char   type,
                           uint   identity,
                           void  *snpPtr,
                           void  *auxiliaryArrayPtr,
                           uint   dimSize,
                           int   *dim);
uint getAuxDim(int *dim, uint preIndex, uint postIndex);
void unstackAuxiliaryInfoAndList();
void memoryCheck();
void stackIncomingResponseArrays(char mode);
void unstackIncomingResponseArrays(char mode);
void unstackIncomingCovariateArrays(char mode);
void unstackIncomingCovariateArrays(char mode);
void stackIncomingArrays(char mode);
void unstackIncomingArrays(char mode);
void checkInteraction();
void stackPreDefinedCommonArrays();
void unstackPreDefinedCommonArrays();
void stackPreDefinedGrowthArrays();
void unstackPreDefinedGrowthArrays();
void stackPreDefinedRestoreArrays();
void unstackPreDefinedRestoreArrays();
void stackPreDefinedPredictArrays();
void unstackPreDefinedPredictArrays();
void stackWeights(double *weight,
                  uint    size,
                  uint   *weightType,
                  uint  **weightSorted,
                  uint   *weightDensitySize);
void unstackWeights(uint    weightType,
                    uint    size,
                    uint   *weightSorted);
void getAtRiskAndEventCounts (uint       treeID,
                              Terminal  *parent,
                              uint      *repMembrIndx,
                              uint       repMembrSize,
                              uint      *allMembrIndx,
                              uint       allMembrSize,
                              uint      *rmbrIterator);
void getLocalRatios (uint treeID, Terminal *parent);
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
void updateEnsembleSurvival(char mode,
                            uint treeID,
                            uint serialTreeID);
void getEnsembleMortality(char      mode, 
                          uint      treeID,
                          uint      obsSize,
                          double  **ensembleMRTptr,
                          uint     *ensembleDen,
                          double   *mortality);
void getEnsembleMortalityCR(char      mode, 
                            uint      treeID,
                            uint      obsSize,
                            double  **ensembleMRTptr,
                            uint     *ensembleDen,
                            double  **cMortality);
void getConditionalConcordanceArrays(uint     j, 
                                     double  *timePtr, 
                                     double  *statusPtr, 
                                     double  *mortalityPtr, 
                                     uint    *genericEnsembleDenPtr,
                                     uint    *meIndividualSize,
                                     uint   **eIndividual,
                                     double  *subsettedTime,
                                     double  *subsettedStatus,
                                     double  *subsettedMortality,
                                     uint    *subsettedEnsembleDen);
double getConcordanceIndex(int     polarity,
                           uint    size, 
                           double *timePtr, 
                           double *statusPtr, 
                           double *predictedOutcome,
                           uint   *oobCount);
void getCRPerformance (char     mode,
                       uint     obsSize,
                       double **responsePtr,
                       double **yearsLost,
                       uint    *denom,
                       double  *performanceVector);
Terminal *makeTerminal();
void freeTerminal(Terminal *parent);
void stackTermLMIIndex(Terminal *tTerm, unsigned int size);
void unstackTermLMIIndex(Terminal *tTerm);
void freeTerminalNodeLocalSurvivalStructures(Terminal *tTerm);
void freeTerminalNodeSurvivalStructuresIntermediate(Terminal *tTerm);
void freeTerminalNodeSurvivalStructuresFinal(Terminal *tTerm);
void freeTerminalNodeNonSurvivalStructures(Terminal *tTerm);
void stackAtRiskAndEventCounts(Terminal *tTerm, unsigned int eTypeSize, unsigned int mTimeSize);
void stackEventTimeIndex(Terminal *tTerm, unsigned int mTimeSize);
void unstackAtRiskAndEventCounts(Terminal *tTerm);
void unstackEventTimeIndex(Terminal *tTerm);
void unstackAtRisk(Terminal *tTerm);
void stackLocalRatio(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize);
void unstackLocalRatio(Terminal *tTerm);
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
void getTerminalInfo(Terminal *termPtr);
void acquireTree(char mode, uint r, uint b);
void updateWeight(char mode, uint b);
void finalizeWeight(char mode);
 
void updateProximity(char mode, uint b);
void finalizeProximity(char mode);
void updateSplitDepth(uint treeID, Node *rootPtr, uint maxDepth);
char pruneBranch(uint obsSize, uint treeID, Node **nodesAtDepth, uint nadCount, uint ptnTarget, uint ptnCurrent);
uint pruneTree(uint obsSize, uint treeID, uint ptnCount);
void unstackAuxiliary2(char mode, uint b);
void unstackAuxiliary(char mode, uint b);
void stackNodeAndTermList(uint treeID, uint length);
void unstackNodeList(uint treeID);
void unstackTermList(uint treeID);
void printPseudoTNInfo(char mode, uint b);
Node *getTerminalNode(uint treeID, uint leaf);
void getRawNodeSize(uint  type,
                    uint  treeID,
                    Node *parent,
                    uint *repMembrIndx,
                    uint *repMembrSize,
                    uint *allMembrIndx,
                    uint *allMembrSize);
char growTree(uint     r,
              char     rootFlag,
              char     multImpFlag,
              uint     b,
              Node    *parent,
              uint    *repMembrIndx,
              uint     repMembrSize,
              uint    *allMembrIndx,
              uint     allMembrSize,
              uint    *bootMembrIndxIter,
              uint    *rmbrIterator,
              uint    *ambrIterator);
void freeTree(uint treeID, Node *parent, char rootFlag);
void getSplitDepth(uint treeID, Node *parent);
void freeSplitDepth(uint treeID);
void saveStatistics(char     mode,
                    uint     b,
                    Node    *parent,
                    ulong   *offset,
                    double  *spltST,
                    double  *spltVR,
                    uint   **uspvST,
                    uint   **mtryID,
                    double **mtryST);
uint getMaximumDepth(Node *parent);
void getNodesAtDepth(Node *parent, uint tagDepth, Node **nodesAtDepth, uint *nadCount);
void getTreeInfo(uint treeID, Node *parent);
void getPTNodeList(Node    *parent,
                   Node   **list,
                   uint    *offset);
void initTerminalNodeMembership(uint       treeID,
                                Terminal  *parent,
                                uint      *allMembrIndx,
                                uint       allMembrSize);
