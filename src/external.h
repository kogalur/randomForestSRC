// *** THIS FILE IS AUTO GENERATED. DO NOT EDIT IT ***
#ifndef RF_EXTERNAL_H
#define RF_EXTERNAL_H

#include "global.h"
/*
#include "nativeInfo.h"
extern NATXDInfo **RF_natXDInfoList;
extern NativeEnsembleInfo **RF_nativeEnsembleInfoList;
extern uint       RF_natXDInfoListSize;
extern uint       RF_nativeEnsembleInfoListSize;
extern PS_arrayXD *RF_rLevelsPYTH;
extern PS_arrayXD *RF_xLevelsPYTH;
*/
/*
extern JNIEnv    *RF_java_env;
extern jobject    RF_java_obj;
extern jclass     RF_java_cls;
extern jclass     RF_java_except_cls;
extern jclass     RF_java_hshmap_cls;
extern jmethodID  RF_java_hshmap_constr;
extern jmethodID  RF_java_hshmap_put;
extern jobject    RF_java_hshmap_obj;
extern jclass     RF_java_ens_cls;
extern jmethodID  RF_java_ens_mid;
extern jmethodID  RF_java_mid_log;
extern jmethodID  RF_java_mid_logError;
extern jmethodID  RF_java_mid_logExit;
extern NAT1DInfo **RF_nat1DInfoList;
extern NAT2DInfo **RF_nat2DInfoList;
extern NativeEnsembleInfo **RF_nativeEnsembleInfoList;
extern uint       RF_nat1DInfoListSize;
extern uint       RF_nat2DInfoListSize;
extern uint       RF_nativeEnsembleInfoListSize;
extern jobject    RF_rLevelsJNIE;
extern jobject    RF_xLevelsJNIE;
*/

#include <Rinternals.h>
#include "nativeInfo.h"
extern SEXP RF_sexpVector[2];
extern SEXP      RF_rLevelsSEXP;
extern SEXP      RF_xLevelsSEXP;

#include "factor.h"
#include "leafLink.h"
#include "node.h"
#include "terminal.h"
#include "snpAuxiliaryInfo.h"
#include "split.h"
#include "quantile.h"
#include "sortedLink.h"
#include "sampling.h"
extern SNPAuxiliaryInfo **RF_snpAuxiliaryInfoList;
extern SNPAuxiliaryInfo **RF_incAuxiliaryInfoList;
extern uint RF_stackCount;
extern uint RF_nativeIndex;
extern uint RF_incStackCount;
extern double   *RF_cpuTime_;
extern uint     *RF_treeID_;
extern uint     *RF_nodeID_;
extern uint     *RF_nodeSZ_;
extern uint     *RF_blnodeID_;
extern uint     *RF_brnodeID_;
extern int      **RF_parmID_;
extern double   **RF_contPT_;
extern uint     **RF_mwcpSZ_;
extern uint     **RF_fsrecID_;
extern uint     **RF_mwcpPT_;
extern uint   **RF_mwcpCT_;
extern ulong     RF_totalNodeCount;
extern ulong     RF_totalTerminalCount;
extern ulong     RF_totalMWCPCount;
extern ulong    *RF_restoreTreeOffset;
extern ulong   **RF_restoreMWCPoffset;
extern uint     *RF_orderedTreeIndex;
extern uint     *RF_serialTreeIndex;
extern uint      RF_serialTreeID;
extern uint      RF_userTreeID;
extern uint     *RF_restoreTreeID;
extern uint      RF_perfBlockCount;
extern uint      RF_serialBlockID;
extern double   *RF_TN_SURV_;
extern double   *RF_TN_MORT_;
extern double   *RF_TN_NLSN_;
extern double   *RF_TN_CSHZ_;
extern double   *RF_TN_CIFN_;
extern double   *RF_TN_REGR_;
extern uint     *RF_TN_CLAS_;
extern int      *RF_seed_;
extern int      *RF_seedVimp_;
extern uint     *RF_optLoGrow_;
extern uint      RF_optLoGrow;
extern uint     *RF_tLeafCount_;
extern uint     *RF_tLeafCount;
extern double   *RF_imputation_;
extern uint     *RF_varUsed_;
extern double   *RF_splitDepth_;
extern uint     *RF_MEMB_ID_;
extern uint     *RF_BOOT_CT_;
extern uint     *RF_PRUN_ID_;
extern uint     *RF_CASE_DPTH_;
extern uint     *RF_RMBR_ID_;
extern uint     *RF_AMBR_ID_;
extern uint     *RF_OMBR_ID_;
extern uint     *RF_IMBR_ID_;
extern uint     *RF_TN_RCNT_;
extern uint     *RF_TN_ACNT_;
extern uint     *RF_TN_OCNT_;
extern uint     *RF_TN_ICNT_;
extern uint     *RF_OOB_SZ_;
extern uint     *RF_IBG_SZ_;
extern double   *RF_perfMRT_;
extern double   *RF_perfCLS_;
extern double   *RF_perfRGR_;
extern double   *RF_holdMRT_;
extern double   *RF_holdCLS_;
extern double   *RF_holdRGR_;
extern uint     *RF_holdBLK_;
extern double   *RF_splitStatLOT_;
extern double   *RF_emprRSK_;
extern double   *RF_oobEmprRSK_;
extern double   *RF_vimpMRT_;
extern double   *RF_vimpCLS_;
extern double   *RF_vimpRGR_;
extern double   *RF_perfBlockMRT_;
extern double   *RF_perfBlockCLS_;
extern double   *RF_perfBlockRGR_;
extern double   *RF_partial_SURV_;
extern double   *RF_partial_CRSK_;
extern double   *RF_partial_CLAS_;
extern double   *RF_partial_REGR_;
extern double   *RF_oobEnsembleSRG_;
extern double   *RF_fullEnsembleSRG_;
extern double   *RF_oobEnsembleCIF_;
extern double   *RF_fullEnsembleCIF_;
extern double   *RF_oobEnsembleSRV_;
extern double   *RF_fullEnsembleSRV_;
extern double   *RF_oobEnsembleMRT_;
extern double   *RF_fullEnsembleMRT_;
extern double   *RF_oobEnsembleCLS_;
extern double   *RF_fullEnsembleCLS_;
extern double   *RF_cseNumCLS_;
extern double   *RF_csvNumCLS_;
extern double   *RF_oobEnsembleRGR_;
extern double   *RF_fullEnsembleRGR_;
extern double   *RF_cseNumRGR_;
extern double   *RF_csvNumRGR_;
extern uint     *RF_cseDen_;
extern uint     *RF_csvDen_;
extern double   *RF_oobEnsembleQNT_;
extern double   *RF_fullEnsembleQNT_;
extern double     *RF_proximity_;
extern double     *RF_distance_;
extern double     *RF_weight_;
extern uint      RF_opt;
extern uint      RF_optHigh;
extern uint      RF_optSup;
extern uint      RF_splitRule;
extern uint      RF_splitCustomIdx;
extern uint      RF_nsplit;
extern uint      RF_nImpute;
extern uint      RF_ntree;
extern uint      RF_nodeSize;
extern int       RF_nodeDepth;
extern uint      RF_crWeightSize;
extern double   *RF_crWeight;
extern uint      RF_perfBlock;
extern uint      RF_mtry;
extern uint      RF_vtry;
extern uint    **RF_vtryArray;
extern uint      RF_vtryMode;
extern uint      RF_vtryBlockSize;
extern uint      RF_subjSize;
extern double   *RF_subjWeight;
extern uint      RF_bootstrapSize;
extern uint    **RF_bootstrapIn;
extern double   *RF_xWeightStat;
extern double   *RF_yWeight;
extern uint      RF_xMarginalSize;
extern uint     *RF_xMarginal;
extern uint     *RF_xMarginalFlag;
extern double   *RF_xWeight;
extern uint      RF_ptnCount;
extern int       RF_numThreads;
extern uint      RF_observationSize;
extern uint      RF_ySize;
extern uint      RF_rTargetCount;
extern uint      RF_rTargetFactorCount;
extern uint      RF_rTargetNonFactorCount;
extern uint     *RF_rTarget;
extern uint     *RF_rTargetFactor;
extern uint     *RF_rTargetNonFactor;
extern double  **RF_responseIn;
extern double  **RF_observationIn;
extern uint     *RF_subjIn;
extern uint      RF_xSize;
extern char     *RF_xType;
extern uint     *RF_xLevelsMax;
extern uint     *RF_xLevelsCnt;
extern uint    **RF_xLevels;
extern uint     *RF_xtType;
extern uint     *RF_stType;
extern char     *RF_rType;
extern uint     *RF_rLevelsMax;
extern uint     *RF_rLevelsCnt;
extern uint    **RF_rLevels;
extern uint      RF_ytry;
extern double    RF_wibsTau;
extern double   **RF_qStarPlus;
extern double    RF_vimpThreshold;
extern uint      RF_fobservationSize;
extern uint      RF_frSize;
extern double  **RF_fresponseIn;
extern double  **RF_fobservationIn;
extern uint      RF_timeIndex;
extern uint      RF_statusIndex;
extern uint     *RF_yIndex;
extern uint      RF_ySizeProxy;
extern uint     *RF_yIndexZero;
extern uint      RF_yIndexZeroSize;
extern char     *RF_testMembershipFlag;  
extern uint      RF_intrPredictorSize;
extern uint     *RF_intrPredictor;
extern char     *RF_importanceFlag;   
extern uint      RF_partialType;
extern uint      RF_partialXvar;
extern uint      RF_partialLength;
extern double   *RF_partialValue;
extern uint      RF_partialLength2;
extern uint     *RF_partialXvar2;
extern double   *RF_partialValue2;
extern uint      RF_partialTimeLength;
extern double   *RF_partialTime;
extern double   *RF_quantile;
extern uint      RF_quantileSize;
extern double    RF_qEpsilon;
extern uint      RF_inv_2qEpsilon;
extern uint     *RF_getTree;
extern uint      RF_xWeightType;
extern uint     *RF_xWeightSorted;
extern uint      RF_xWeightDensitySize;
extern uint      RF_yWeightType;
extern uint     *RF_yWeightSorted;
extern uint      RF_yWeightDensitySize;
extern uint      RF_subjWeightType;
extern uint     *RF_subjWeightSorted;
extern uint      RF_subjWeightDensitySize;
extern uint      RF_eventTypeSize;
extern uint      RF_feventTypeSize;
extern uint      RF_mStatusSize;
extern uint     *RF_eventType;
extern uint     *RF_eventTypeIndex;
extern uint     *RF_eIndividualSize;
extern uint    **RF_eIndividualIn;
extern uint     *RF_subjSlot;
extern uint     *RF_subjSlotCount;
extern uint    **RF_subjList;
extern uint     *RF_caseMap;
extern uint      RF_subjCount;
extern uint      *RF_classLevelSize;
extern uint     **RF_classLevel;
extern uint     **RF_classLevelIndex;
extern uint    ***RF_cIndividualIn;
extern double   *RF_timeInterest;
extern uint      RF_timeInterestSize;
extern uint      RF_sortedTimeInterestSize;
extern double   *RF_masterTime;
extern uint      RF_masterTimeSize;
extern uint     *RF_masterTimeIndexIn;
extern uint     *RF_masterToInterestTimeMap;
extern double    RF_wibsTauTime;
extern uint      RF_wibsTauTimeIdx;
extern uint      RF_rFactorCount;
extern uint     *RF_rFactorMap;
extern uint     *RF_rFactorIndex;
extern uint     *RF_rFactorSize;
extern uint     *RF_rFactorSizeTest;
extern uint      RF_mrFactorSize;
extern uint      RF_fmrFactorSize;
extern uint     *RF_mrFactorIndex;
extern uint     *RF_fmrFactorIndex;
extern uint      RF_rNonFactorCount;
extern uint     *RF_rNonFactorMap;
extern uint     *RF_rNonFactorIndex;
extern uint      RF_xFactorCount;
extern uint     *RF_xFactorMap;
extern uint     *RF_xFactorIndex;
extern uint     *RF_xFactorSize;
extern uint      RF_mxFactorSize;
extern uint      RF_fmxFactorSize;
extern uint     *RF_mxFactorIndex;
extern uint     *RF_fmxFactorIndex;
extern uint      RF_xNonFactorCount;
extern uint     *RF_xNonFactorMap;
extern uint     *RF_xNonFactorIndex;
extern uint      RF_rMaxFactorLevel;
extern uint      RF_xMaxFactorLevel;
extern uint      RF_maxFactorLevel;
extern char      RF_mStatusFlag;
extern char      RF_mTimeFlag;
extern char      RF_mResponseFlag;
extern char      RF_mPredictorFlag;
extern char      RF_fmStatusFlag;
extern char      RF_fmTimeFlag;
extern char      RF_fmResponseFlag;
extern char      RF_fmPredictorFlag;
extern uint     *RF_mRecordMap;
extern uint     *RF_fmRecordMap;
extern uint      RF_mRecordSize;
extern uint      RF_fmRecordSize;
extern uint     *RF_mRecordIndex;
extern uint     *RF_fmRecordIndex;
extern uint      RF_mpIndexSize;
extern uint      RF_fmpIndexSize;
extern int     **RF_mpSign;
extern int     **RF_fmpSign;
extern int      *RF_mpIndex;
extern int      *RF_fmpIndex;
extern uint     *RF_getTreeIndex;
extern uint      RF_getTreeCount;
extern SplitRuleObj *RF_splitRuleObj;
extern double   **RF_importancePtr;
extern double **RF_sImputeResponsePtr;
extern double **RF_sImputePredictorPtr;
extern double **RF_sImputeDataPtr;
extern uint  **RF_MEMB_ID_ptr;
extern uint  **RF_BOOT_CT_ptr;
extern uint  **RF_PRUN_ID_ptr;
extern uint  **RF_CASE_DPTH_ptr;
extern uint     **RF_treeID_ptr;
extern uint     **RF_nodeID_ptr;
extern uint     **RF_nodeSZ_ptr;
extern uint     **RF_blnodeID_ptr;
extern uint     **RF_brnodeID_ptr;
extern int     ***RF_parmID_ptr;
extern double  ***RF_contPT_ptr;
extern uint    ***RF_mwcpSZ_ptr;
extern uint    ***RF_fsrecID_ptr;
extern uint    ***RF_mwcpPT_ptr;
extern uint     **RF_mwcpCT_ptr;
extern double   **RF_spltST_ptr;
extern uint     **RF_dpthST_ptr;
extern uint  **RF_RMBR_ID_ptr;
extern uint  **RF_AMBR_ID_ptr;
extern uint  **RF_OMBR_ID_ptr;
extern uint  **RF_IMBR_ID_ptr;
extern uint  **RF_TN_RCNT_ptr;
extern uint  **RF_TN_ACNT_ptr;
extern uint  **RF_TN_OCNT_ptr;
extern uint  **RF_TN_ICNT_ptr;
extern LeafLinkedObj **RF_leafLinkedObjHead;
extern LeafLinkedObj **RF_leafLinkedObjTail;
extern double **RF_proximityPtr;
extern double **RF_proximityDenPtr;
extern double  *RF_proximityDen;
extern double **RF_distancePtr;
extern double **RF_distanceDenPtr;
extern double  *RF_distanceDen;
extern uint    RF_rejectedTreeCount;
extern uint    RF_validTreeCount;
extern uint    RF_stumpedTreeCount;
extern double  ***RF_TN_SURV_ptr;
extern double  ***RF_TN_MORT_ptr;
extern double  ***RF_TN_NLSN_ptr;
extern double ****RF_TN_CSHZ_ptr;
extern double ****RF_TN_CIFN_ptr;
extern double  ***RF_TN_REGR_ptr;
extern uint   ****RF_TN_CLAS_ptr;
extern Terminal ****RF_vimpMembership;
extern Terminal ***RF_partMembership;
extern double  ***RF_vimpMRTstd;
extern double ****RF_vimpCLSstd;
extern double  ***RF_vimpRGRstd;
extern double   **RF_vimpEnsembleDen;
extern double   **RF_perfMRTptr;
extern double  ***RF_perfCLSptr;
extern double   **RF_perfRGRptr;
extern double   **RF_cseNumCLSptr;
extern double   **RF_cseNumRGRptr;
extern uint      *RF_cseDENptr;
extern double  ***RF_csvNumCLSptr;
extern double  ***RF_csvNumRGRptr;
extern uint     **RF_csvDENptr;
extern double   **RF_blkEnsembleMRTnum;
extern double  ***RF_blkEnsembleCLSnum;
extern double   **RF_blkEnsembleRGRnum;
extern double    *RF_blkEnsembleDen;
extern double  ***RF_vimpMRTblk;
extern double ****RF_vimpCLSblk;
extern double  ***RF_vimpRGRblk;
extern double   **RF_perfMRTblk;
extern double  ***RF_perfCLSblk;
extern double   **RF_perfRGRblk;
extern double   **RF_vimpMRTptr;
extern double  ***RF_vimpCLSptr;
extern double   **RF_vimpRGRptr;
extern double  ***RF_holdMRTptr;
extern double ****RF_holdCLSptr;
extern double  ***RF_holdRGRptr;
extern double  ****RF_holdMRTstd;
extern double *****RF_holdCLSstd;
extern double  ****RF_holdRGRstd;
extern double   ***RF_holdEnsembleDen;
extern uint **RF_holdoutMap;
extern uint *RF_holdBLKptr;
extern uint **RF_runningHoldoutCount;
extern uint ***RF_blockSerialTreeIndex;
extern double   **RF_splitStatLOTptr;
extern double   **RF_emprRSKptr;
extern double   **RF_oobEmprRSKptr;
extern double  ****RF_partSURVptr;
extern double  ****RF_partCLASptr;
extern double   ***RF_partREGRptr;
extern double ***RF_oobEnsembleSRGptr;
extern double ***RF_fullEnsembleSRGptr;
extern double ***RF_oobEnsembleCIFptr;
extern double ***RF_fullEnsembleCIFptr;
extern double  **RF_oobEnsembleSRVptr;
extern double  **RF_fullEnsembleSRVptr;
extern double  **RF_oobEnsembleMRTptr;
extern double  **RF_fullEnsembleMRTptr;
extern double ***RF_oobEnsembleCLSptr;
extern double ***RF_fullEnsembleCLSptr;
extern double  **RF_oobEnsembleRGRptr;
extern double  **RF_fullEnsembleRGRptr;
extern double  ***RF_oobEnsembleQNTptr;
extern double  ***RF_fullEnsembleQNTptr;
extern uint    **RF_oobQuantileStreamSize;
extern uint    **RF_fullQuantileStreamSize;
extern LookUpInfo *** RF_oobQuantileSearchTree;
extern LookUpInfo *** RF_fullQuantileSearchTree;
extern QuantileObj ***RF_oobQuantileHead;
extern QuantileObj ***RF_oobQuantileTail;
extern uint         **RF_oobQuantileLinkLength;
extern QuantileObj ***RF_fullQuantileHead;
extern QuantileObj ***RF_fullQuantileTail;
extern uint         **RF_fullQuantileLinkLength;
extern double ***RF_oobEnsembleSRGnum;
extern double ***RF_fullEnsembleSRGnum;
extern double ***RF_oobEnsembleCIFnum;
extern double ***RF_fullEnsembleCIFnum;
extern double  **RF_oobEnsembleSRVnum;
extern double  **RF_fullEnsembleSRVnum;
extern double  **RF_oobEnsembleMRTnum;
extern double  **RF_fullEnsembleMRTnum;
extern double ***RF_oobEnsembleCLSnum;
extern double ***RF_fullEnsembleCLSnum;
extern double  **RF_oobEnsembleRGRnum;
extern double  **RF_fullEnsembleRGRnum;
extern double   *RF_oobEnsembleDen;
extern double   *RF_fullEnsembleDen;
extern double ***RF_splitDepthPtr;
extern double **RF_weightPtr;
extern uint    *RF_weightDenom;
extern char    **RF_dmRecordBootFlag;
extern uint   **RF_varUsedPtr;
extern uint    *RF_oobSize;
extern uint    *RF_ibgSize;
extern uint    *RF_nodeCount;
extern uint   **RF_mwcpPtr;
extern uint    *RF_pLeafCount;
extern uint    *RF_maxDepth;
extern Node    **RF_root;
extern Node   ***RF_nodeMembership;  
extern Node   ***RF_fnodeMembership;
extern Node   ***RF_pNodeMembership;
extern Node   ***RF_pNodeList;
extern Terminal   ***RF_tTermMembership;
extern Terminal   ***RF_ftTermMembership;
extern uint       ***RF_utTermMembership;
extern uint        **RF_utTermMembershipCount;
extern uint        **RF_utTermMembershipAlloc;
extern Terminal   ***RF_pTermMembership;
extern Terminal   ***RF_tTermList;
extern Terminal   ***RF_pTermList;
extern LeafLinkedObj ***RF_hTermMembership;
extern uint    **RF_bootMembershipIndex;
extern uint     *RF_identityMembershipIndex;
extern uint      RF_identityMembershipIndexSize;
extern uint     *RF_fidentityMembershipIndex;
extern char    **RF_bootMembershipFlag;
extern uint    **RF_bootMembershipCount;
extern char    **RF_oobMembershipFlag;
extern uint    **RF_ibgMembershipIndex;
extern uint    **RF_oobMembershipIndex;
extern double  **RF_status;
extern double  **RF_time;
extern double ***RF_response;
extern double  **RF_ftime;
extern double  **RF_fstatus;
extern double ***RF_fresponse;
extern double ***RF_observation;
extern double ***RF_fobservation;
extern uint    ***RF_observationRank;
extern uint     **RF_observationUniqueSize;
extern SortedLinkedObj **RF_sortedLinkedHead;
extern SortedLinkedObj **RF_sortedLinkedTail;
extern uint    **RF_masterTimeIndex;
extern Factor ***RF_factorList;
extern double   *RF_rFactorThreshold;
extern uint     *RF_rFactorMinority;
extern uint     *RF_rFactorMajority;
extern char     *RF_rFactorMinorityFlag;
extern void (*acquireTree) (char mode, uint r, uint b);
extern Node *(*antiMembership) (uint     treeID,
                         Node    *parent,
                         uint     individual,
                         uint     vimpX,
                         double **xArray);
extern Node *(*randomMembership) (uint     treeID,
                           Node    *parent,
                           uint     individual,
                           uint     vimpX,
                           double **xArray);
extern Node *(*getMembership) (uint     treeID,
                        Node    *parent,
                        uint     individual,
                        double **xArray);
extern void (*partialMembership) (uint       treeID,
                           Node      *parent,
                           uint       partialIndex,
                           uint      *allMembrIndx,
                           uint       obsSize,
                           double   **xArray,
                           Terminal **membership);
extern char (*getVariance) (uint, uint*, uint, uint*, double*, double*, double*);
extern char (*growTree) (uint, char, char, uint, Node*, uint*, uint*, uint*);
extern char (*getPreSplitResult) (uint, Node*, char, char);
extern DistributionObj *(*stackRandomCovariates) (uint, Node*);
extern void (*unstackRandomCovariates) (uint, DistributionObj*);
extern char (*selectRandomCovariates) (uint, Node*, DistributionObj*, char*, uint*, uint*);
extern DistributionObj *(*stackRandomResponses) (uint, Node*);
extern void (*unstackRandomResponses) (uint, DistributionObj*);
extern char (*selectRandomResponses) (uint, Node*, DistributionObj*, uint*, uint*);
extern uint (*virtuallySplitNode) (uint, Node*, char, uint, double*, uint*, void*, uint, char*, uint*, uint, uint*);
extern uint (*stackAndConstructSplitVector) (uint treeID, Node *parent, uint covariate);
extern char (*updateMaximumSplit) (uint    treeID,
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
extern char (*forkAndUpdate) (uint       treeID,
                       Node      *parent,
                       uint      *repMembrIndx,
                       uint       repMembrSize,
                       uint      *allMembrIndx,
                       uint       allMembrSize,
                       char       multImpFlag,
                       SplitInfo *info,
                       uint      *leafCount,
                       Node     **nodeMembership);
extern void (*freeNode) (Node *parent);
extern void (*unstackSplitVector) (uint   treeID,
                            Node  *parent,
                            uint   splitLength,
                            char   factorFlag,
                            uint   splitVectorSize,
                            uint   mwcpSizeAbsolute,
                            char   deterministicSplitFlag,
                            void  *splitVectorPtr,
                            char   multImpFlag,
                            uint  *indxx);
extern char (*randomSplit) (uint, Node*, SplitInfoMax*, GreedyObj*, char);
extern char (*regressionXwghtSplit) (uint, Node*, SplitInfoMax*, GreedyObj*, char);
extern char (*classificationXwghtSplit) (uint, Node*, SplitInfoMax*, GreedyObj*, char);
extern char (*multivariateSplit) (uint, Node*, SplitInfoMax*, GreedyObj*, char);
extern char (*unsupervisedSplit) (uint, Node*, SplitInfoMax*, GreedyObj*, char);
extern void (*getConditionalClassificationIndex) (uint,
                                           uint,
                                           double*,
                                           double**,
                                           double*,
                                           double*,
                                           double*);
extern double (*getGMeanIndex) (uint,
                         uint,
                         double*,
                         double*,
                         double*);
extern float (*ran1A) (uint);
extern void  (*randomSetChain) (uint, int);
extern int   (*randomGetChain) (uint);
extern float (*ran1B) (uint);
extern void  (*randomSetUChain) (uint, int);
extern int   (*randomGetUChain) (uint);
extern void  (*randomSetUChainVimp) (uint, int);
extern int   (*randomGetUChainVimp) (uint);
extern float (*ran1D) (uint);
extern void  (*randomSetChainVimp) (uint, int);
extern int   (*randomGetChainVimp) (uint);
extern uint RF_forestChunkSize;
extern uint RF_forestChunkCount;
extern uint *RF_forestChunkSizeActual;
extern uint **RF_forestChunk;
#ifdef _OPENMP
extern omp_lock_t   *RF_lockPartial;
extern omp_lock_t  **RF_lockWeight;
extern omp_lock_t   *RF_lockWeightRow;
extern omp_lock_t  **RF_lockVimp;
extern omp_lock_t   *RF_lockVimpRow;
extern omp_lock_t   *RF_lockVimpCol;
extern omp_lock_t  **RF_lockVimpHoldout;
extern omp_lock_t  **RF_lockMRToens;
extern omp_lock_t  **RF_lockMRTfens;
extern omp_lock_t  **RF_lockSRVoens;
extern omp_lock_t  **RF_lockSRVfens;
extern omp_lock_t ***RF_lockSRGoens;
extern omp_lock_t ***RF_lockSRGfens;
extern omp_lock_t ***RF_lockCIFoens;
extern omp_lock_t ***RF_lockCIFfens;
extern omp_lock_t ***RF_lockCLSoens;
extern omp_lock_t ***RF_lockCLSfens;
extern omp_lock_t   *RF_lockDENoens;
extern omp_lock_t   *RF_lockDENfens;
extern omp_lock_t   *RF_lockQNToens;
extern omp_lock_t   *RF_lockQNTfens;
extern omp_lock_t      RF_lockPerf;
extern omp_lock_t      RF_lockEnsbUpdtCount;
#endif
extern uint RF_ensbUpdtCount;
extern customFunction customFunctionArray[4][16];
extern uint   RF_userTraceFlag;
extern time_t RF_userTimeStart;
extern time_t RF_userTimeSplit;  

#endif
