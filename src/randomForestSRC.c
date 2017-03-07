
#include "randomForestSRC.h"
//  **********************************************************************
//  **********************************************************************
//  
//    RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
//  
//    This program is free software; you can redistribute it and/or
//    modify it under the terms of the GNU General Public License
//    as published by the Free Software Foundation; either version 3
//    of the License, or (at your option) any later version.
//  
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//  
//    You should have received a copy of the GNU General Public
//    License along with this program; if not, write to the Free
//    Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
//    Boston, MA  02110-1301, USA.
//  
//    ----------------------------------------------------------------
//    Project Partially Funded By: 
//    ----------------------------------------------------------------
//    Dr. Ishwaran's work was funded in part by DMS grant 1148991 from the
//    National Science Foundation and grant R01 CA163739 from the National
//    Cancer Institute.
//  
//    Dr. Kogalur's work was funded in part by grant R01 CA163739 from the 
//    National Cancer Institute.
//    ----------------------------------------------------------------
//    Written by:
//    ----------------------------------------------------------------
//      Hemant Ishwaran, Ph.D.
//      Director of Statistical Methodology
//      Professor, Division of Biostatistics
//      Clinical Research Building, Room 1058
//      1120 NW 14th Street
//      University of Miami, Miami FL 33136
//  
//      email:  hemant.ishwaran@gmail.com
//      URL:    http://web.ccs.miami.edu/~hishwaran
//      --------------------------------------------------------------
//      Udaya B. Kogalur, Ph.D.
//      Adjunct Staff
//      Department of Quantitative Health Sciences
//      Cleveland Clinic Foundation
//      
//      Kogalur & Company, Inc.
//      5425 Nestleway Drive, Suite L1
//      Clemmons, NC 27012
//  
//      email:  ubk@kogalur.com
//      URL:    https://github.com/kogalur/randomForestSRC
//      --------------------------------------------------------------
//  
//  **********************************************************************
//  **********************************************************************


uint RF_stackCount;
char  *sexpString[RF_SEXP_CNT] = {
  "",              
  "",              
  "allEnsbCHF",    
  "oobEnsbCHF",    
  "allEnsbCIF",    
  "oobEnsbCIF",    
  "allEnsbSRV",    
  "oobEnsbSRV",    
  "allEnsbMRT",    
  "oobEnsbMRT",    
  "allEnsbCLS",    
  "oobEnsbCLS",    
  "allEnsbRGR",    
  "oobEnsbRGR",    
  "perfSurv",      
  "perfClas",      
  "perfRegr",      
  "proximity",     
  "leafCount",     
  "treeID",        
  "nodeID",        
  "parmID",        
  "contPT",        
  "mwcpSZ",        
  "mwcpPT",        
  "seed",          
  "vimpSurv",      
  "vimpClas",      
  "vimpRegr",      
  "imputation",    
  "oobImputation", 
  "varUsed",       
  "splitDepth",    
  "nodeMembership",
  "pstnMembership",
  "bootMembership",
  "rmbrMembership",
  "ambrMembership",
  "tnRCNT",        
  "tnACNT",        
  "spltST",        
  "spltVR",        
  "weight",        
  "tnSURV",        
  "tnMORT",        
  "tnNLSN",        
  "tnCSHZ",        
  "tnCIFN",        
  "tnREGR",        
  "tnCLAS",        
  "uspvST",        
  "mtryID",        
  "mtryST",        
  "mwcpCount",     
  "partialSurv",   
  "partialClas",   
  "partialRegr",   
  "distance"       
};
SEXP sexpVector[RF_SEXP_CNT];
uint     *RF_treeID_;
uint     *RF_nodeID_;
uint     *RF_parmID_;
uint     *RF_mwcpSZ_;
double   *RF_contPT_;
uint     *RF_mwcpPT_;
double   *RF_spltST_;
double   *RF_spltVR_;
uint     *RF_uspvST_;
uint     *RF_mtryID_;
double   *RF_mtryST_;
ulong     RF_totalNodeCount;
ulong     RF_totalNodeCount1;
ulong     RF_totalNodeCount2;
ulong     RF_totalTerminalCount;
uint     *RF_theoreticalMaxtLeafCount;
ulong    *RF_restoreTreeOffset;
ulong    *RF_restoreMWCPOffset;
uint     *RF_serialTreeIndex;
uint      RF_serialTreeCount;
uint     *RF_restoreTreeID;
double   *RF_TN_SURV_;
double   *RF_TN_MORT_;
double   *RF_TN_NLSN_;
double   *RF_TN_CSHZ_;
double   *RF_TN_CIFN_;
double   *RF_TN_REGR_;
uint     *RF_TN_CLAS_;
int      *RF_seed_;
uint     *RF_tLeafCount_;
double   *RF_imputation_;
uint     *RF_varUsed_;
double   *RF_splitDepth_;
uint     *RF_MEMB_ID_;
uint     *RF_BOOT_CT_;
uint     *RF_PRUN_ID_;
uint     *RF_RMBR_ID_;
uint     *RF_AMBR_ID_;
uint     *RF_TN_RCNT_;
uint     *RF_TN_ACNT_;
double   *RF_perfMRT_;
double   *RF_perfCLS_;
double   *RF_perfRGR_;
double   *RF_vimpMRT_;
double   *RF_vimpCLS_;
double   *RF_vimpRGR_;
double   *RF_partial_SURV_;
double   *RF_partial_CRSK_;
double   *RF_partial_CLAS_;
double   *RF_partial_REGR_;
double   *RF_oobEnsembleSRG_;
double   *RF_fullEnsembleSRG_;
double   *RF_oobEnsembleCIF_;
double   *RF_fullEnsembleCIF_;
double   *RF_oobEnsembleSRV_;
double   *RF_fullEnsembleSRV_;
double   *RF_oobEnsembleMRT_;
double   *RF_fullEnsembleMRT_;
double   *RF_oobEnsembleCLS_;
double   *RF_fullEnsembleCLS_;
double   *RF_oobEnsembleRGR_;
double   *RF_fullEnsembleRGR_;
double     *RF_proximity_;
uint      RF_opt;
uint      RF_optHigh;
uint      RF_splitRule;
uint      RF_splitCustomIdx;
uint      RF_splitRandomCount;
uint      RF_nImpute;
uint      RF_forestSize;
uint      RF_minimumNodeSize;
int       RF_maximumNodeDepth;
double   *RF_crWeight;
uint      RF_randomCovariateCount;
uint      RF_bootstrapSize;
int      *RF_bootstrap;
uint    **RF_bootstrapIn;
double   *RF_caseWeight;
double   *RF_xWeight;
uint      RF_ptnCount;
int       RF_numThreads;
uint      RF_observationSize;
uint      RF_rSize;
uint      RF_rTargetCount;
uint      RF_rTargetFactorCount;
uint      RF_rTargetNonFactorCount;
uint     *RF_rTarget;
uint     *RF_rTargetFactor;
uint     *RF_rTargetNonFactor;
double   *RF_rData;
uint      RF_xSize;
double   *RF_xData;
double  **RF_responseIn;
double  **RF_observationIn;
SEXP      RF_sexp_xType;
char    **RF_xType;
int      *RF_xLevels;
SEXP      RF_sexp_rType;
char    **RF_rType;
int      *RF_rLevels;
uint      RF_randomResponseCount;
uint      RF_fobservationSize;
uint      RF_frSize;
double   *RF_frData;
double   *RF_fxData;
double  **RF_fresponseIn;
double  **RF_fobservationIn;
uint      RF_timeIndex;
uint      RF_statusIndex;
uint     *RF_yIndex;
uint      RF_ySize;
char     *RF_testMembershipFlag;  
uint      RF_intrPredictorSize;
uint     *RF_intrPredictor;
uint      RF_sobservationSize;
uint     *RF_sobservationIndv;
char     *RF_importanceFlag;   
uint      RF_partialType;
uint      RF_partialXvar;
uint      RF_partialLength;
double   *RF_partialValue;
uint      RF_partialLength2;
uint     *RF_partialXvar2;
double   *RF_partialValue2;
uint      RF_partialTimeLength;
double   *RF_partialTime;
uint      RF_xWeightType;
uint     *RF_xWeightSorted;
uint     *RF_xWeightDensity;
uint      RF_xWeightDensitySize;
uint      RF_caseWeightType;
uint     *RF_caseWeightSorted;
uint     *RF_caseWeightDensity;
uint      RF_caseWeightDensitySize;
uint      RF_eventTypeSize;
uint      RF_feventTypeSize;
uint      RF_mStatusSize;
uint     *RF_eventType;
uint     *RF_eventTypeIndex;
uint     *RF_eIndividualSize;
uint    **RF_eIndividualIn;
uint      *RF_classLevelSize;
uint     **RF_classLevel;
uint     **RF_classLevelIndex;
uint    ***RF_cIndividualIn;
double   *RF_timeInterest;
uint      RF_timeInterestSize;
uint      RF_sortedTimeInterestSize;
double   *RF_masterTime;
uint      RF_masterTimeSize;
uint     *RF_masterTimeIndexIn;
uint      RF_rFactorCount;
uint     *RF_rFactorMap;
uint     *RF_rFactorIndex;
uint     *RF_rFactorSize;
uint      RF_mrFactorSize;
uint      RF_fmrFactorSize;
uint     *RF_mrFactorIndex;
uint     *RF_fmrFactorIndex;
uint      RF_rNonFactorCount;
uint     *RF_rNonFactorMap;
uint     *RF_rNonFactorIndex;
uint      RF_xFactorCount;
uint     *RF_xFactorMap;
uint     *RF_xFactorIndex;
uint     *RF_xFactorSize;
uint      RF_mxFactorSize;
uint      RF_fmxFactorSize;
uint     *RF_mxFactorIndex;
uint     *RF_fmxFactorIndex;
uint      RF_xNonFactorCount;
uint     *RF_xNonFactorMap;
uint     *RF_xNonFactorIndex;
uint      RF_rMaxFactorLevel;
uint      RF_xMaxFactorLevel;
uint      RF_maxFactorLevel;
char      RF_mStatusFlag;
char      RF_mTimeFlag;
char      RF_mResponseFlag;
char      RF_mPredictorFlag;
char      RF_fmStatusFlag;
char      RF_fmTimeFlag;
char      RF_fmResponseFlag;
char      RF_fmPredictorFlag;
uint     *RF_mRecordMap;
uint     *RF_fmRecordMap;
uint      RF_mRecordSize;
uint      RF_fmRecordSize;
uint     *RF_mRecordIndex;
uint     *RF_fmRecordIndex;
uint      RF_mpIndexSize;
uint      RF_fmpIndexSize;
int     **RF_mpSign;
int     **RF_fmpSign;
int      *RF_mpIndex;
int      *RF_fmpIndex;
double   **RF_importancePtr;
double **RF_sImputeResponsePtr;
double **RF_sImputePredictorPtr;
double **RF_sOOBImputeResponsePtr;
double **RF_sOOBImputePredictorPtr;
uint  **RF_MEMB_ID_ptr;
uint  **RF_BOOT_CT_ptr;
uint  **RF_PRUN_ID_ptr;
uint  **RF_RMBR_ID_ptr;
uint  **RF_AMBR_ID_ptr;
uint  **RF_TN_RCNT_ptr;
uint  **RF_TN_ACNT_ptr;
double **RF_proximityPtr;
double **RF_proximityDenPtr;
double  *RF_proximityDen;
uint    RF_rejectedTreeCount;
uint    RF_validTreeCount;
uint    RF_stumpedTreeCount;
uint     **RF_uspvST_ptr;
uint     **RF_mtryID_ptr;
double   **RF_mtryST_ptr;
double  ***RF_TN_SURV_ptr;
double  ***RF_TN_MORT_ptr;
double  ***RF_TN_NLSN_ptr;
double ****RF_TN_CSHZ_ptr;
double ****RF_TN_CIFN_ptr;
double  ***RF_TN_REGR_ptr;
uint   ****RF_TN_CLAS_ptr;
double  **RF_perfMRTptr;
double ***RF_perfCLSptr;
double  **RF_perfRGRptr;
double  **RF_outcomeRGR;
double  **RF_vimpMRTptr;
double ***RF_vimpCLSptr;
double  **RF_vimpRGRptr;
double  ***RF_vimpEnsembleMRT;
double ****RF_vimpEnsembleCLS;
double  ***RF_vimpEnsembleRGR;
double  **RF_perfMRTleo;
double ***RF_perfCLSleo;
double  **RF_perfRGRleo;
double  ***RF_vimpMRTleo;
double ****RF_vimpCLSleo;
double  ***RF_vimpRGRleo;
double  ****RF_partSURVptr;
double  ****RF_partCLASptr;
double   ***RF_partREGRptr;
double ***RF_oobEnsembleSRGptr;
double ***RF_fullEnsembleSRGptr;
double ***RF_oobEnsembleCIFptr;
double ***RF_fullEnsembleCIFptr;
double  **RF_oobEnsembleSRVptr;
double  **RF_fullEnsembleSRVptr;
double  **RF_oobEnsembleMRTptr;
double  **RF_fullEnsembleMRTptr;
double ***RF_oobEnsembleCLSptr;
double ***RF_fullEnsembleCLSptr;
double  **RF_oobEnsembleRGRptr;
double  **RF_fullEnsembleRGRptr;
double ***RF_oobEnsembleSRGnum;
double ***RF_fullEnsembleSRGnum;
double ***RF_oobEnsembleCIFnum;
double ***RF_fullEnsembleCIFnum;
double  **RF_oobEnsembleSRVnum;
double  **RF_fullEnsembleSRVnum;
double  **RF_oobEnsembleMRTnum;
double  **RF_fullEnsembleMRTnum;
double ***RF_oobEnsembleCLSnum;
double ***RF_fullEnsembleCLSnum;
double  **RF_oobEnsembleRGRnum;
double  **RF_fullEnsembleRGRnum;
uint     *RF_oobEnsembleDen;
uint     *RF_fullEnsembleDen;
uint     **RF_vimpEnsembleDen;
double ***RF_splitDepthPtr;
char    **RF_dmRecordBootFlag;
double **RF_performancePtr;
uint   **RF_varUsedPtr;
uint    *RF_oobSize;
uint    *RF_ibgSize;
uint    *RF_soobSize;
uint    *RF_tLeafCount;
uint    *RF_nodeCount;
uint    *RF_mwcpCount;
uint    *RF_mwcpIterator;
uint   **RF_mwcpPtr;
uint    *RF_pLeafCount;
uint    *RF_maxDepth;
Node    **RF_root;
Node   ***RF_nodeMembership;
Node   ***RF_fnodeMembership;
Node   ***RF_pNodeMembership;
Node   ***RF_tNodeList;
uint     *RF_tNodeListLength;
Node   ***RF_pNodeList;
Terminal   ***RF_tTermMembership;
Terminal   ***RF_ftTermMembership;
uint       ***RF_utTermMembership;
uint        **RF_utTermMembershipCount;
uint        **RF_utTermMembershipAlloc;
Terminal   ***RF_pTermMembership;
Terminal   ***RF_tTermList;
Terminal   ***RF_pTermList;
uint    **RF_bootMembershipIndex;
uint     *RF_identityMembershipIndex;
uint     *RF_fidentityMembershipIndex;
char    **RF_bootMembershipFlag;
uint    **RF_bootMembershipCount;
char    **RF_oobMembershipFlag;
uint    **RF_ibgMembershipIndex;
uint    **RF_oobMembershipIndex;
uint     *RF_orderedLeafCount;
Terminal ****RF_vimpMembership;
Terminal ***RF_partMembership;
double  **RF_status;
double  **RF_time;
double ***RF_response;
double  **RF_ftime;
double  **RF_fstatus;
double ***RF_fresponse;
double ***RF_observation;
double ***RF_fobservation;
uint    **RF_masterTimeIndex;
Factor ***RF_factorList;
float (*ran1A) (uint);
void  (*randomSetChain) (uint, int);
int   (*randomGetChain) (uint);
float (*ran1B) (uint);
void  (*randomSetUChain) (uint, int);
int   (*randomGetUChain) (uint);
float (*ran1C) (uint);
void  (*randomSetUChainCov) (uint, int);
int   (*randomGetUChainCov) (uint);
char (*genericSplit) (uint,
                      Node*,
                      uint*,
                      uint,
                      uint*,
                      uint,
                      uint*,
                      double*,
                      uint*,
                      uint**,
                      double*,
                      char**,
                      char*,
                      char);
customFunction customFunctionArray[4][16];
uint   RF_userTraceFlag;
time_t RF_userTimeStart;
time_t RF_userTimeSplit;  
#define IA      16807
#define IM      2147483647
#define AM      (1.0/IM)
#define IQ      127773
#define IR      2836
#define NTAB    32
#define NDIV    (1+(IM-1)/NTAB)
#define EPS     1.2e-7
#define RNMX    (1.0-EPS)
#define LCG_IM  714025
#define LCG_IA  1366
#define LCG_IC  150889
int  *ran1A_iy;
int **ran1A_iv;
int  *ran1B_iy;
int **ran1B_iv;
int  *ran1C_iy;
int **ran1C_iv;
int      *seed1AValue;
int      *seed1BValue;
int      *seed1CValue;
void randomStack(uint bSize, uint pSize) {
  uint b;
  ran1A_iy = ivector(1, bSize);
  ran1A_iv = imatrix(1, bSize, 1, NTAB);
  ran1B_iy = ivector(1, bSize);
  ran1B_iv = imatrix(1, bSize, 1, NTAB);
  ran1C_iy = ivector(1, bSize);
  ran1C_iv = imatrix(1, bSize, 1, NTAB);
  for (b = 1; b <= bSize; b++) {
    ran1A_iy[b] = 0;
    ran1B_iy[b] = 0;
    ran1C_iy[b] = 0;
  }
  seed1AValue = ivector(1, bSize);
  seed1BValue = ivector(1, bSize);
  seed1CValue = ivector(1, bSize);
}
void randomUnstack(uint bSize, uint pSize) {
  free_ivector(ran1A_iy, 1, bSize);
  free_imatrix(ran1A_iv, 1, bSize, 1, NTAB);
  free_ivector(ran1B_iy, 1, bSize);
  free_imatrix(ran1B_iv, 1, bSize, 1, NTAB);
  free_ivector(ran1C_iy, 1, bSize);
  free_imatrix(ran1C_iv, 1, bSize, 1, NTAB);
  free_ivector(seed1AValue, 1, bSize);
  free_ivector(seed1BValue, 1, bSize);
  free_ivector(seed1CValue, 1, bSize);
}
void randomSetChainParallel(uint b, int value) {
  seed1AValue[b] = value;
}
void randomSetUChainParallel(uint b, int value) {
  seed1BValue[b] = value;
}
void randomSetUChainParallelCov(uint b, int value) {
  seed1CValue[b] = value;
}
void randomSetChainSerial(uint b, int value) {
  seed1AValue[1] = value;
}
void randomSetUChainSerial(uint b, int value) {
  seed1BValue[1] = value;
}
void randomSetUChainSerialCov(uint b, int value) {
  seed1CValue[1] = value;
}
int randomGetChainParallel(uint b) {
  return seed1AValue[b];
}
int randomGetUChainParallel(uint b) {
  return seed1BValue[b];
}
int randomGetUChainParallelCov(uint b) {
  return seed1CValue[b];
}
int randomGetChainSerial(uint b) {
  return seed1AValue[1];
}
int randomGetUChainSerial(uint b) {
  return seed1BValue[1];
}
int randomGetUChainSerialCov(uint b) {
  return seed1CValue[1];
}
float randomChainParallel(uint b) {
  return  ran1_generic(& ran1A_iy[b], ran1A_iv[b], & seed1AValue[b]);
}
float randomUChainParallel(uint b) {
  return  ran1_generic(& ran1B_iy[b], ran1B_iv[b], & seed1BValue[b]);
}
float randomUChainParallelCov(uint b) {
  return  ran1_generic(& ran1C_iy[b], ran1C_iv[b], & seed1CValue[b]);
}
float randomChainSerial(uint b) {
  return  ran1_generic(& ran1A_iy[1], ran1A_iv[1], & seed1AValue[1]);
}
float randomUChainSerial(uint b) {
  return  ran1_generic(& ran1B_iy[1], ran1B_iv[1], & seed1BValue[1]);
}
float randomUChainSerialCov(uint b) {
  return  ran1_generic(& ran1C_iy[1], ran1C_iv[1], & seed1CValue[1]);
}
float ran1_generic(int *iy, int *iv, int *idum) {
  int j, k;
  float temp;
  if (*idum <= 0 || !(*iy)) {
    if (-(*idum) < 1) {
      *idum = 1;
    }
    else {
      *idum = -(*idum);
    }
    for (j = NTAB+7; j >= 0; j--) {
      k = (*idum) / IQ;
      *idum = IA * (*idum - k * IQ) - IR * k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    (*iy) = iv[1];
  }
  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0) *idum += IM;
  j = (*iy) / NDIV;
  (*iy) = iv[j];
  iv[j] = *idum;
  if ((temp = AM * (*iy)) > RNMX) {
    return RNMX;
  }
  else {
    return temp;
  }
}
void lcgenerator(unsigned int *seed, unsigned char reset) {
  if (reset) {
    if (*seed >= LCG_IM) (*seed) %= LCG_IM;
  }
  else {
    *seed = (LCG_IA * (*seed) + LCG_IC) % LCG_IM;
  }
}
float ran1_original(int *idum) {
  int j;
  int k;
  static int iy = 0;
  static int iv[NTAB];
  float temp;
  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) {
      *idum = 1;
    }
    else {
      *idum = -(*idum);
    }
    for (j = NTAB+7; j >= 0; j--) {
      k = (*idum) / IQ;
      *idum = IA * (*idum - k * IQ) - IR * k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }
  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0) *idum += IM;
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if ((temp = AM * iy) > RNMX) {
    return RNMX;
  }
  else {
    return temp;
  }
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
unsigned int upower (unsigned int x, unsigned int n) {
  unsigned int p;
  if ((x >= 2) & (n > (sizeof(unsigned int) * 8) - 1)) {
    nrerror("Overflow in upower(), exponent too large.");
  }
  for (p = 1; n > 0; --n) {
    p = p * x;
  }
  return p;
}
unsigned int upower2 (unsigned int n) {
  unsigned int p;
  if (n > (sizeof(unsigned int) * 8) - 1) {
    nrerror("Overflow in upower2(), exponent too large.");
  }
  p = ((unsigned int) 1) << n;
  return p;
}
unsigned int ulog2 (unsigned int n) {
  unsigned int p;
  p = 0;
  while (n > 1) {
    n = n >> 1;
    p++;
  }
  return p;
}
void hpsort(double *ra, unsigned int n) {
  unsigned int i, ir, j, l;
  double rra;
  if (n < 2) return;
  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra = ra[--l];
    }
    else {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1) {
        ra[1] = rra;
        break;
      }
    }
    i = l;
    j = l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
        ra[i] = ra[j];
        i = j;
        j <<= 1;
      }
      else {
        j = ir+1;
      }
    }
    ra[i] = rra;
  }
}
void hpsortui(unsigned int *ra, unsigned int n) {
  unsigned int i, ir, j, l;
  unsigned int rra;
  if (n < 2) return;
  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra = ra[--l];
    }
    else {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1) {
        ra[1] = rra;
        break;
      }
    }
    i = l;
    j = l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
        ra[i] = ra[j];
        i = j;
        j <<= 1;
      }
      else {
        j = ir+1;
      }
    }
    ra[i] = rra;
  }
}
#ifdef SWAP
#undef SWAP
#endif
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50
void sort(double *arr, unsigned int n) {
  unsigned int i, j, k, l;
  unsigned int ir;
  unsigned int *istack, jstack;
  double a, temp;
  if (n < 1) nrerror("\n n of zero (0) length in indexx().");
  l  = 1;
  ir = n;
  jstack = 0;
  istack = uivector(1, NSTACK);
  for (;;) {
    if (ir-l < M) {
      for (j = l+1; j <= ir; j++) {
        a=arr[j];
        for (i = j-1; i >= l; i--) {
          if (arr[i] <= a) break;
          arr[i+1] = arr[i];
        }
        arr[i+1] = a;
      }
      if (jstack == 0) break;
      ir = istack[jstack--];
      l  = istack[jstack--]; 
    } else {
      k = (l+ir) >> 1; 
      SWAP(arr[k], arr[l+1]);
      if (arr[l] > arr[ir]) {
        SWAP(arr[l], arr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1], arr[ir]);
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l], arr[l+1]);
      }
      i = l+1; 
      j = ir;
      a = arr[l+1];
      for (;;) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        SWAP(arr[i], arr[j]);
      }
      arr[l+1] = arr[j]; 
      arr[j] = a;
      jstack += 2;
      if (jstack > NSTACK) nrerror("NSTACK too small in sort().");
      if (ir-i+1 >= j-l) {
        istack[jstack] = ir;
        istack[jstack-1] = i;
        ir = j-1;
      }
      else {
        istack[jstack] = j-1;
        istack[jstack-1] = l;
        l=i;
      }
    }
  }
  free_uivector(istack,1,NSTACK);
}
#undef SWAP
#undef M
#undef NSTACK
#ifdef SWAP
#undef SWAP
#endif
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50
void indexx(unsigned int n, double *arr, unsigned int *indx) {
  unsigned int i, j, k, l;
  unsigned int indxt, itemp, ir;
  unsigned int *istack, jstack;
  double a;
  if (n < 1) nrerror("\n n of zero (0) length in indexx().");
  l  = 1;
  ir = n;
  jstack = 0;
  istack = uivector(1, NSTACK);
  for (j=1; j<=n; j++) indx[j]=j;
  for (;;) {
    if (ir-l < M) {
      for (j = l+1; j <= ir; j++) {
        indxt = indx[j];
        a = arr[indxt];
        for (i=j-1; i>=l; i--) {
          if (arr[indx[i]] <= a) break;
          indx[i+1] = indx[i];
        }
        indx[i+1] = indxt;
      }
      if (jstack == 0) break;
      ir = istack[jstack--];
      l  = istack[jstack--];
    }
    else {
      k = (l+ir) >> 1;
      SWAP(indx[k], indx[l+1]);
      if (arr[indx[l]] > arr[indx[ir]]) {
        SWAP(indx[l], indx[ir])
      }
      if (arr[indx[l+1]] > arr[indx[ir]]) {
        SWAP(indx[l+1], indx[ir])
      }
      if (arr[indx[l]] > arr[indx[l+1]]) {
        SWAP(indx[l], indx[l+1])
      }
      i = l+1;
      j = ir;
      indxt = indx[l+1];
      a = arr[indxt];
      for (;;) {
        do i++; while (arr[indx[i]] < a);
        do j--; while (arr[indx[j]] > a);
        if (j < i) break;
        SWAP(indx[i], indx[j])
      }
      indx[l+1] = indx[j];
      indx[j] = indxt;
      jstack += 2;
      if (jstack > NSTACK) nrerror("NSTACK too small in indexx().");
      if (ir-i+1 >= j-l) {
        istack[jstack] = ir;
        istack[jstack-1] = i;
        ir = j-1;
      }
      else {
        istack[jstack] = j-1;
        istack[jstack-1] = l;
        l = i;
      }
    }
  }
  free_uivector(istack, 1, NSTACK);
}
#undef SWAP
#undef M
#undef NSTACK
#define FREE_ARG char*
#define NR_END 1
void nrerror(char error_text[]) {
  RFprintf("\nRF-SRC");
  RFprintf("\nRF-SRC:  *** ERROR *** ");
  RFprintf("\nRF-SRC:  Numerical Recipes Run-Time Error:");
  RFprintf("\nRF-SRC:  %s", error_text);
  RFprintf("\nRF-SRC:  Please Contact Technical Support.");
  error("\nRF-SRC:  The application will now exit.\n");
}
void *gblock(size_t size) {
  void *v = (void *) malloc(size);
  if (!v) nrerror("\n  Allocation Failure in gblock().");
  return v;
}
void free_gblock(void *v, size_t size) {
  free((FREE_ARG) (v));
}
void *gvector(unsigned long long nl, unsigned long long nh, size_t size) {
  if (nh < nl) nrerror("\n  Illegal indices in gvector().");
  void *v = gblock((size_t) ((nh-nl+1+NR_END) * size));
  return v;
}
void free_gvector(void *v, unsigned long long nl, unsigned long long nh, size_t size) {
  if (nh < nl) nrerror("\n  Illegal indices in free_gvector().");
  free_gblock(v, (nh-nl+1+NR_END) * size);
}
char *cvector(unsigned long long nl, unsigned long long nh) {
  return ((char *) gvector(nl, nh, sizeof(char)) -nl+NR_END);
}
void free_cvector(char *v, unsigned long long nl, unsigned long long nh) {
  free_gvector(v+nl-NR_END, nl, nh, sizeof(char));
}
char **cmatrix(unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  char **v = (char **) new_vvector(nrl, nrh, NRUTIL_CPTR);
  for(unsigned long long i = nrl; i <= nrh; i++) {
    v[i] = cvector(ncl, nch);
  }
  return v;
}
void free_cmatrix(char **v, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  for(unsigned long long i = nrl; i <= nrh; i++) {
    free_cvector(v[i], ncl, nch);
  }
  free_new_vvector(v, nrl, nrh, NRUTIL_CPTR);
}
int *ivector(unsigned long long nl, unsigned long long nh) {
  return ((int *) gvector(nl, nh, sizeof(int)) -nl+NR_END);
}
void free_ivector(int *v, unsigned long long nl, unsigned long long nh) {
  free_gvector(v+nl-NR_END, nl, nh, sizeof(int));
}
int **imatrix(unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  int **v = (int **) new_vvector(nrl, nrh, NRUTIL_IPTR);
  for(unsigned long long i = nrl; i <= nrh; i++) {
    v[i] = ivector(ncl, nch);
  }
  return v;
}
void free_imatrix(int **v, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  for(unsigned long long i = nrl; i <= nrh; i++) {
    free_ivector(v[i], ncl, nch);
  }
  free_new_vvector(v, nrl, nrh, NRUTIL_IPTR);
}
unsigned int *uivector(unsigned long long nl, unsigned long long nh) {
  return ((unsigned int *) gvector(nl, nh, sizeof(unsigned int)) -nl+NR_END);
}
void free_uivector(unsigned int *v, unsigned long long nl, unsigned long long nh) {
  free_gvector(v+nl-NR_END, nl, nh, sizeof(unsigned int));
}
unsigned int **uimatrix(unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  unsigned int **v = (unsigned int **) new_vvector(nrl, nrh, NRUTIL_UPTR);
  for(unsigned long long i = nrl; i <= nrh; i++) {
    v[i] = uivector(ncl, nch);
  }
  return v;
}
void free_uimatrix(unsigned int **v, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  for(unsigned long long i = nrl; i <= nrh; i++) {
    free_uivector(v[i], ncl, nch);
  }
  free_new_vvector(v, nrl, nrh, NRUTIL_UPTR);
}
unsigned long *ulvector(unsigned long long nl, unsigned long long nh) {
  return ((unsigned long *) gvector(nl, nh, sizeof(unsigned long)) -nl+NR_END);
}
void free_ulvector(unsigned long *v, unsigned long long nl, unsigned long long nh) {
  free_gvector(v+nl-NR_END, nl, nh, sizeof(unsigned long));
}
double *dvector(unsigned long long nl, unsigned long long nh) {
  return ((double *) gvector(nl, nh, sizeof(double)) -nl+NR_END);
}
void free_dvector(double *v, unsigned long long nl, unsigned long long nh) {
  free_gvector(v+nl-NR_END, nl, nh, sizeof(double));
}
double **dmatrix(unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  double **v = (double **) new_vvector(nrl, nrh, NRUTIL_DPTR);
  for(unsigned long long i = nrl; i <= nrh; i++) {
    v[i] = dvector(ncl, nch);
  }
  return v;
}
void free_dmatrix(double **v, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  for(unsigned long long i = nrl; i <= nrh; i++) {
    free_dvector(v[i], ncl, nch);
  }
  free_new_vvector(v, nrl, nrh, NRUTIL_DPTR);
}
double ***dmatrix3(unsigned long long n3l, unsigned long long n3h, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  double ***v = (double ***) new_vvector(n3l, n3h, NRUTIL_DPTR2);
  for(unsigned long long i = n3l; i <= n3h; i++) {
    v[i] = dmatrix(nrl, nrh, ncl, nch);
  }
  return v;
}
void free_dmatrix3(double ***v, unsigned long long n3l, unsigned long long n3h, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  for(unsigned long long i = n3l; i <= n3h; i++) {
    free_dmatrix(v[i], nrl, nrh, ncl, nch);
  }
  free_new_vvector(v, n3l, n3h, NRUTIL_DPTR2);
}
double ****dmatrix4(unsigned long long n4l, unsigned long long n4h, unsigned long long n3l, unsigned long long n3h, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  double ****v = (double ****) new_vvector(n4l, n4h, NRUTIL_DPTR3);
  for(unsigned long long i = n4l; i <= n4h; i++) {
    v[i] = dmatrix3(n3l, n3h, nrl, nrh, ncl, nch);
  }
  return v;
}
void free_dmatrix4(double ****v, unsigned long long n4l, unsigned long long n4h, unsigned long long n3l, unsigned long long n3h, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  for(unsigned long long i = n4l; i <= n4h; i++) {
    free_dmatrix3(v[i], n3l, n3h, nrl, nrh, ncl, nch);
  }
  free_new_vvector(v, n4l, n4h, NRUTIL_DPTR3);
}
void *new_vvector(unsigned long long nl, unsigned long long nh, enum alloc_type type) {
  void *v;
  v = NULL;  
  switch(type){
  case NRUTIL_DPTR:
    v = ((double **) gvector(nl, nh, sizeof(double*)) -nl+NR_END);
    break;
  case NRUTIL_UPTR:
    v = ((unsigned int **) gvector(nl, nh, sizeof(unsigned int*)) -nl+NR_END);
    break;
  case NRUTIL_DPTR2:
    v = ((double ***) gvector(nl, nh, sizeof(double**)) -nl+NR_END);
    break;
  case NRUTIL_NPTR:
    v = ((Node **) gvector(nl, nh, sizeof(Node*)) -nl+NR_END);
    break;
  case NRUTIL_NPTR2:
    v = ((Node ***) gvector(nl, nh, sizeof(Node**)) -nl+NR_END);
    break;
  case NRUTIL_CPTR:
    v = ((char **) gvector(nl, nh, sizeof(char*)) -nl+NR_END);
    break;
  case NRUTIL_DPTR4:
    v = ((double *****) gvector(nl, nh, sizeof(double****)) -nl+NR_END);
    break;
  case NRUTIL_TPTR:
    v = ((Terminal **) gvector(nl, nh, sizeof(Terminal*)) -nl+NR_END);
    break;
  case NRUTIL_TPTR2:
    v = ((Terminal ***) gvector(nl, nh, sizeof(Terminal**)) -nl+NR_END);
    break;
  case NRUTIL_IPTR:
    v = ((int **) gvector(nl, nh, sizeof(int*)) -nl+NR_END);
    break;
  case NRUTIL_NPTR3:
    v = ((Node ****) gvector(nl, nh, sizeof(Node***)) -nl+NR_END);
    break;
  case NRUTIL_FPTR:
    v = ((Factor **) gvector(nl, nh, sizeof(Factor*)) -nl+NR_END);
    break;
  case NRUTIL_FPTR2:
    v = ((Factor ***) gvector(nl, nh, sizeof(Factor**)) -nl+NR_END);
    break;
  case NRUTIL_DPTR3:
    v = ((double ****) gvector(nl, nh, sizeof(double***)) -nl+NR_END);
    break;
  case NRUTIL_UPTR3:
    v = ((unsigned int ****) gvector(nl, nh, sizeof(unsigned int***)) -nl+NR_END);
    break;
  case NRUTIL_UPTR4:
    v = ((unsigned int *****) gvector(nl, nh, sizeof(unsigned int****)) -nl+NR_END);
    break;
  case NRUTIL_UPTR2:
    v = ((unsigned int ***) gvector(nl, nh, sizeof(unsigned int**)) -nl+NR_END);
    break;
  default:
    v = NULL;
    nrerror("\n  Illegal case in new_vvector().");
    break;
  }
  return v;
}
void free_new_vvector(void *v, unsigned long long nl, unsigned long long nh, enum alloc_type type) {
  switch(type){
  case NRUTIL_DPTR:
    free_gvector((double*) v+nl-NR_END, nl, nh, sizeof(double*));
    break;
  case NRUTIL_UPTR:
    free_gvector((unsigned int*) v+nl-NR_END, nl, nh, sizeof(unsigned int*));
    break;
  case NRUTIL_DPTR2:
    free_gvector((double**) v+nl-NR_END, nl, nh, sizeof(double**));
    break;
  case NRUTIL_NPTR:
    free_gvector((Node*) v+nl-NR_END, nl, nh, sizeof(Node*));
    break;
  case NRUTIL_NPTR2:
    free_gvector((Node**) v+nl-NR_END, nl, nh, sizeof(Node**));
    break;
  case NRUTIL_CPTR:
    free_gvector((char*) v+nl-NR_END, nl, nh, sizeof(char*));
    break;
  case NRUTIL_DPTR4:
    free_gvector((double****) v+nl-NR_END, nl, nh, sizeof(double****));
    break;
  case NRUTIL_TPTR:
    free_gvector((Terminal*) v+nl-NR_END, nl, nh, sizeof(Terminal*));
    break;
  case NRUTIL_TPTR2:
    free_gvector((Terminal**) v+nl-NR_END, nl, nh, sizeof(Terminal**));
    break;
  case NRUTIL_IPTR:
    free_gvector((int*) v+nl-NR_END, nl, nh, sizeof(int*));
    break;
  case NRUTIL_NPTR3:
    free_gvector((Node***) v+nl-NR_END, nl, nh, sizeof(Node***));
    break;
  case NRUTIL_FPTR:
    free_gvector((Factor*) v+nl-NR_END, nl, nh, sizeof(Factor*));
    break;
  case NRUTIL_FPTR2:
    free_gvector((Factor**) v+nl-NR_END, nl, nh, sizeof(Factor**));
    break;
  case NRUTIL_DPTR3:
    free_gvector((double***) v+nl-NR_END, nl, nh, sizeof(double***));
    break;
  case NRUTIL_UPTR3:
    free_gvector((unsigned int***) v+nl-NR_END, nl, nh, sizeof(unsigned int***));
    break;
  case NRUTIL_UPTR4:
    free_gvector((unsigned int****) v+nl-NR_END, nl, nh, sizeof(unsigned int****));
    break;
  case NRUTIL_UPTR2:
    free_gvector((unsigned int**) v+nl-NR_END, nl, nh, sizeof(unsigned int**));
    break;
  default:
    nrerror("\n  Illegal case in free_new_vvector().");
    break;
  }
}
#undef FREE_ARG
#undef NR_END
void nrCopyMatrix(unsigned int **new, unsigned int **old, unsigned int nrow, unsigned int ncol) {
  unsigned int i,j;
  for (i = 1; i <= nrow; i++) {
    for (j = 1; j <= ncol; j++) {
      new[i][j] = old[i][j];
    }
  }
}
void nrCopyVector(char *new, char *old, unsigned int ncol) {
  unsigned int j;
  for (j = 1; j <= ncol; j++) {
    new[j] = old[j];
  }
}
void testEndianness() {
  unsigned int     test = 0x12345678;
  unsigned int *testPtr = & test;
  RFprintf("\nTest of Endianness:  ");
  RFprintf("%2x %2x %2x %2x \n",
           *((char *) testPtr),
           *((char *) testPtr + 1),
           *((char *) testPtr + 2),
           *((char *) testPtr + 3));
}
Node *makeNode(unsigned int xSize,
               unsigned int permissibleSizeAlloc,
               unsigned int urStatSize,
               unsigned int mtrySize) {
  unsigned int i;
  Node *parent = (Node*) gblock((size_t) sizeof(Node));
  if (xSize > 0) {
    parent -> xSize = xSize;
    parent -> permissibleSplit = cvector(1, xSize);
    for (i = 1; i <= xSize; i++) {
      (parent -> permissibleSplit)[i] = TRUE;
    }
  }
  else {
    parent -> permissibleSplit = NULL;
    parent -> xSize = 0;
  }
  if (permissibleSizeAlloc > 0) {
    (parent -> permissibleSplitIndex) = uivector(1, permissibleSizeAlloc);
    parent -> permissibleSizeAlloc = permissibleSizeAlloc;
    parent -> permissibleSizeActual = 0;
  }
  else {
    (parent -> permissibleSplitIndex) = NULL;
    (parent -> permissibleSizeAlloc)  = 0;
    (parent -> permissibleSizeActual)  = 0;
  }
  parent -> mate               = NULL;
  parent -> left               = NULL;
  parent -> right              = NULL;
  parent -> splitFlag            = TRUE;
  parent -> splitParameter       = 0;
  parent -> splitValueCont       = NA_REAL;
  parent -> splitValueFactSize   = 0;
  parent -> splitValueFactPtr    = NULL;
  parent -> splitStatistic       = NA_REAL;
  parent -> variance             = NA_REAL;
  parent -> mean                 = NA_REAL;
  parent -> urStatSize           = urStatSize;
  if (urStatSize > 0) {
    parent -> urStat = uivector(1, urStatSize);
    for (i = 1; i <= urStatSize; i++) {
      (parent -> urStat)[i] = 0;
    }
  }
  else {
    parent -> urStat = NULL;
  }
  parent -> mtrySize             = mtrySize;
  if (mtrySize > 0) {
    parent -> mtryIndx = uivector(1, mtrySize);
    parent -> mtryStat = dvector(1, mtrySize);
    for (i = 1; i <= mtrySize; i++) {
      (parent -> mtryIndx)[i] = 0;
      (parent -> mtryStat)[i] = NA_REAL;      
    }
  }
  else {
    parent -> mtryIndx = NULL;
    parent -> mtryStat = NULL;
  }
  parent -> nodeID               = 0;
  parent -> depth                = 0;
  parent -> splitDepth           = NULL;
  parent -> pseudoTerminal       = FALSE;
  parent -> mpIndexSize          = 0;
  parent -> fmpIndexSize         = 0;
  parent -> mpSign               = NULL;
  parent -> fmpSign              = NULL;
  parent -> imputed              = FALSE;
  parent -> lmpIndex             = NULL;
  parent -> flmpIndex            = NULL;
  parent -> lmpValue             = NULL;
  parent -> lmpIndexAllocSize    = 0;
  parent -> flmpIndexAllocSize   = 0;
  parent -> lmpIndexActualSize   = 0;
  parent -> flmpIndexActualSize  = 0;
  return parent;
}
void freeNode(Node         *parent) {
  if (parent -> xSize > 0) {
    free_cvector(parent -> permissibleSplit, 1, parent -> xSize);
    parent -> permissibleSplit = NULL;
    parent -> xSize = 0;
  }
  if (parent -> permissibleSizeAlloc > 0) {
    free_uivector(parent -> permissibleSplitIndex, 1, parent -> permissibleSizeAlloc);
    parent -> permissibleSplitIndex = NULL;
    parent -> permissibleSizeAlloc = 0;
  }
  if ((parent -> splitValueFactSize) > 0) {
    if (parent -> splitValueFactPtr != NULL) {
      free_uivector(parent -> splitValueFactPtr, 1, parent -> splitValueFactSize);
      parent -> splitValueFactPtr = NULL;
    }
  }
  if ((parent -> splitParameter) == 0) {
    if ((parent -> depth) > 0) {
      if (parent -> splitDepth != NULL) {
        free_uivector(parent -> splitDepth, 1, parent -> depth);
        parent -> splitDepth = NULL;
      }
    }
  }
  if ((parent -> urStatSize) > 0) {
    if ((parent -> urStat) != NULL) {
      free_uivector(parent -> urStat, 1, parent -> urStatSize);
      parent -> urStat = NULL;
    }
  }
  if ((parent -> mtrySize) > 0) {
    if ((parent -> mtryIndx) != NULL) {
      free_uivector(parent -> mtryIndx, 1, parent -> mtrySize);
      parent -> mtryIndx = NULL;
    }
    if ((parent -> mtryStat) != NULL) {
      free_dvector(parent -> mtryStat, 1, parent -> mtrySize);
      parent -> mtryStat = NULL;
    }
  }
  unstackMPSign(parent);
  unstackFMPSign(parent);
  unstackNodeLMPIndex(parent);
  unstackNodeFLMPIndex(parent);
  free_gblock(parent, sizeof(Node));
}
void getNodeInfo(Node *nodePtr) {
  unsigned int i;
  Rprintf("\nNodeInfo:  %20x", nodePtr);
  Rprintf("\n   LeafCnt   SpltParm  ");
  Rprintf("\n%10d %10d \n", nodePtr -> nodeID, nodePtr -> splitParameter);
  if (nodePtr -> splitValueFactSize > 0) {
    Rprintf("FactorInfo %20x \n", nodePtr -> splitValueFactPtr);
    Rprintf("0x ");
    for (i = nodePtr -> splitValueFactSize; i >= 1; i--) {
      Rprintf("%8x ", (nodePtr -> splitValueFactPtr)[i]);
    }
  }
  else {
    Rprintf(" %12.4f \n", nodePtr -> splitValueCont);
  }
  Rprintf("\nSplit Statistic \n");
  Rprintf(" %12.4f \n", nodePtr -> splitStatistic);
  Rprintf("\nNode Variance \n");
  Rprintf(" %12.4f \n", nodePtr -> variance);
  Rprintf("\nPermissible Flag Size:          %10d", nodePtr -> xSize);
  Rprintf("\nPermissible Index Size Alloc:   %10d", nodePtr -> permissibleSizeAlloc);
  Rprintf("\nPermissible Index Size Actual:  %10d", nodePtr -> permissibleSizeActual);
  Rprintf("\n mpIndexSize   = %20d", nodePtr -> mpIndexSize);
  Rprintf("\n fmpIndexSize  = %20d", nodePtr -> fmpIndexSize);
  Rprintf("\n");
  Rprintf("\n mpSign       = %20x", nodePtr -> mpSign);
  Rprintf("\n fmpSign      = %20x", nodePtr -> fmpSign);
  Rprintf("\n");
  Rprintf("\n lmpIndexActualSize        = %20d", nodePtr -> lmpIndexActualSize);
  Rprintf("\n flmpIndexActualSize       = %20d", nodePtr -> flmpIndexActualSize);
  Rprintf("\n lmpIndexAllocSize         = %20d", nodePtr -> lmpIndexAllocSize);
  Rprintf("\n flmpIndexAllocSize        = %20d", nodePtr -> flmpIndexAllocSize);
  Rprintf("\n");
  Rprintf("\n lmpIndex            = %20x", nodePtr -> lmpIndex);
  Rprintf("\n flmpIndex           = %20x", nodePtr -> flmpIndex);
  Rprintf("\n");
}
void setParent(Node *daughter, Node *parent) {
  daughter -> parent = parent;
}
void setLeftDaughter(Node *daughter, Node *parent) {
  parent -> left = daughter;
}
void setRightDaughter(Node *daughter, Node *parent) {
  parent -> right = daughter;
}
char forkNode(Node         *parent,
              unsigned int  splitParameter,
              double        splitValueMaxCont,
              unsigned int  splitValueMaxFactSize,
              unsigned int *splitValueMaxFactPtr) {
  unsigned int i;
  if (parent == NULL) {
    RFprintf("\nRF-SRC:  *** WARNING *** ");
    RFprintf("\nRF-SRC:  Inconsistent call to forkNode().  ");
    RFprintf("\nRF-SRC:  The parent node is NULL.");
    return FALSE;
  }
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    RFprintf("\nRF-SRC:  *** WARNING *** ");
    RFprintf("\nRF-SRC:  Inconsistent call to forkNode().  ");
    RFprintf("\nRF-SRC:  The daughter nodes are NON-NULL.");
    return FALSE;
  }
  if (parent -> splitFlag == FALSE) {
    RFprintf("\nRF-SRC:  *** WARNING *** ");
    RFprintf("\nRF-SRC:  Inconsistent call to forkNode().  ");
    RFprintf("\nRF-SRC:  The split flag is FALSE.");
    return FALSE;
  }
  if (parent -> xSize < splitParameter) {
    RFprintf("\nRF-SRC:  *** WARNING *** ");
    RFprintf("\nRF-SRC:  Inconsistent call to forkNode().  ");
    RFprintf("\nRF-SRC:  The split parameter index is out of range [1, xSize].");
    return FALSE;
  }
  Node *left  = makeNode(parent -> xSize, parent -> permissibleSizeActual, parent -> urStatSize, parent -> mtrySize);
  Node *right = makeNode(parent -> xSize, parent -> permissibleSizeActual, parent -> urStatSize, parent -> mtrySize);
  parent -> splitParameter = splitParameter;
  parent -> splitValueCont = splitValueMaxCont;
  parent -> splitValueFactSize = splitValueMaxFactSize;
  parent -> splitValueFactPtr = splitValueMaxFactPtr;
  setParent(left, parent);
  setParent(right, parent);
  setLeftDaughter(left, parent);
  setRightDaughter(right, parent);
  if (parent -> xSize > 0) {
    for (i=1; i <= parent -> xSize; i++) {
      left  -> permissibleSplit[i] = right -> permissibleSplit[i] = parent -> permissibleSplit[i];
    }
    free_cvector(parent -> permissibleSplit, 1, parent -> xSize);
    parent -> permissibleSplit = NULL;
    parent -> xSize = 0;
  }
  if (parent -> permissibleSizeActual > 0) {
    for (i=1; i <= parent -> permissibleSizeActual; i++) {
      left -> permissibleSplitIndex[i] =  right -> permissibleSplitIndex[i] = parent -> permissibleSplitIndex[i];
    }
    free_uivector(parent -> permissibleSplitIndex, 1, parent -> permissibleSizeAlloc);
    parent -> permissibleSplitIndex = NULL;
    parent -> permissibleSizeAlloc  = 0;
  }
  parent -> splitFlag = FALSE;
  return TRUE;
}
void stackMPSign(Node *tNode, unsigned int mpIndexSize) {
  if (tNode -> mpIndexSize > 0) {
    if (tNode -> mpIndexSize != mpIndexSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  mpIndexSize has been previously defined:  %10d vs %10d", tNode -> mpIndexSize, mpIndexSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> mpIndexSize = mpIndexSize;
  }
  tNode -> mpSign = ivector(1, tNode -> mpIndexSize);
}
void unstackMPSign(Node *tNode) {
  if(tNode -> mpIndexSize > 0) {
    if (tNode -> mpSign != NULL) {
      free_ivector(tNode -> mpSign, 1, tNode -> mpIndexSize);
      tNode -> mpSign = NULL;
    }
  }
}
void stackFMPSign(Node *tNode, unsigned int fmpIndexSize) {
  if (tNode -> fmpIndexSize > 0) {
    if (tNode -> fmpIndexSize != fmpIndexSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  fmpIndexSize has been previously defined:  %10d vs %10d", tNode -> fmpIndexSize, fmpIndexSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> fmpIndexSize = fmpIndexSize;
  }
  tNode -> fmpSign = ivector(1, tNode -> fmpIndexSize);
}
void unstackFMPSign(Node *tNode) {
  if(tNode -> fmpIndexSize > 0) {
    if (tNode -> fmpSign != NULL) {
      free_ivector(tNode -> fmpSign, 1, tNode -> fmpIndexSize);
      tNode -> fmpSign = NULL;
    }
  }
}
void stackNodeLMPIndex(Node *tNode, unsigned int size) {
  if (tNode -> lmpIndexAllocSize > 0) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  lmpIndex has been previously defined:  %10d vs %10d", tNode -> lmpIndexAllocSize, size);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tNode -> lmpIndexAllocSize = size;
  }
  tNode -> lmpIndex = uivector(1, tNode -> lmpIndexAllocSize);
  tNode -> lmpValue = dvector(1, tNode -> lmpIndexAllocSize);
}
void unstackNodeLMPIndex(Node *tNode) {
  if(tNode -> lmpIndexAllocSize > 0) {
    if (tNode -> lmpIndex != NULL) {
      free_uivector(tNode -> lmpIndex, 1, tNode -> lmpIndexAllocSize);
      free_dvector(tNode -> lmpValue, 1, tNode -> lmpIndexAllocSize);
      tNode -> lmpIndex = NULL;
      tNode -> lmpValue = NULL;
      tNode -> lmpIndexAllocSize = 0;
    }
  }
}
void stackNodeFLMPIndex(Node *tNode, unsigned int size) {
  if (tNode -> flmpIndexAllocSize > 0) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  flmpIndex has been previously defined:  %10d vs %10d", tNode -> flmpIndexAllocSize, size);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tNode -> flmpIndexAllocSize = size;
  }
  tNode -> flmpIndex = uivector(1, tNode -> flmpIndexAllocSize);
  tNode -> flmpValue = dvector(1, tNode -> flmpIndexAllocSize);
}
void unstackNodeFLMPIndex(Node *tNode) {
  if(tNode -> flmpIndexAllocSize > 0) {
    if (tNode -> flmpIndex != NULL) {
      free_uivector(tNode -> flmpIndex, 1, tNode -> flmpIndexAllocSize);
      free_dvector(tNode -> flmpValue, 1, tNode -> flmpIndexAllocSize);
      tNode -> flmpIndex = NULL;
      tNode -> flmpIndexAllocSize = 0;
    }
  }
}
void stackSplitDepth(Node *tNode, unsigned int depth) {
  if (tNode -> depth > 0) {
    if (tNode -> depth != depth) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  depth has been previously defined:  %10d vs %10d", tNode -> depth, depth);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> depth = depth;
  }
  tNode -> splitDepth = uivector(1, tNode -> depth);
}
void unstackSplitDepth(Node *tNode) {
  if (tNode -> splitDepth != NULL) {
    free_uivector(tNode -> splitDepth, 1, tNode -> depth);
    tNode -> splitDepth = NULL;
  }
}
Terminal *makeTerminal() {
  Terminal *parent = (Terminal*) gblock((size_t) sizeof(Terminal));
  parent -> lmiIndex      = NULL;
  parent -> lmiValue      = NULL;
  parent -> lmiSize       = 0;
  parent -> lmiAllocSize  = 0;
  parent -> nodeID     = 0;
  parent -> eTypeSize            = 0;
  parent -> mTimeSize            = 0;
  parent -> eTimeSize            = 0;
  parent -> sTimeSize            = 0;
  parent -> atRiskCount          = NULL;
  parent -> eventCount           = NULL;
  parent -> eventTimeIndex       = NULL;
  parent -> localRatio           = NULL;
  parent -> localCSH             = NULL;
  parent -> localCIF             = NULL;
  parent -> localSurvival        = NULL;
  parent -> localNelsonAalen     = NULL;
  parent -> CSH                  = NULL;
  parent -> CIF                  = NULL;
  parent -> survival             = NULL;
  parent -> nelsonAalen          = NULL;
  parent -> rfCount              = 0;
  parent -> rfSize               = NULL;
  parent -> multiClassProb       = NULL;
  parent -> maxClass             = NULL;
  parent -> rnfCount             = 0;
  parent -> meanResponse         = NULL;
  parent -> weight               = 0.0;
  parent -> membrCount           = 0;
  return parent;
}
void freeTerminal(Terminal        *parent) {
  unstackTermLMIIndex(parent);
  freeTerminalNodeSurvivalStructuresIntermediate(parent);
  freeTerminalNodeSurvivalStructuresFinal(parent);
  freeTerminalNodeNonSurvivalStructures(parent);
  free_gblock(parent, sizeof(Terminal));
}
void freeTerminalNodeLocalSurvivalStructures(Terminal *tTerm) {
  unstackLocalRatio(tTerm);
  unstackLocalSurvival(tTerm);
  unstackLocalNelsonAalen(tTerm);
  if (tTerm -> eTypeSize > 1) {
    unstackLocalCSH(tTerm);
    unstackLocalCIF(tTerm);
  }
  unstackEventTimeIndex(tTerm);
}
void freeTerminalNodeSurvivalStructuresIntermediate(Terminal *tTerm) {
  unstackSurvival(tTerm);
  unstackNelsonAalen(tTerm);
  unstackCSH(tTerm);
  unstackCIF(tTerm);
}
void freeTerminalNodeSurvivalStructuresFinal(Terminal *tTerm) {
  unstackMortality(tTerm);
}
void freeTerminalNodeNonSurvivalStructures(Terminal *tTerm) {
  unstackMultiClassProb(tTerm);
  unstackMeanResponse(tTerm);
}
void stackAtRiskAndEventCounts(Terminal *tTerm, unsigned int eTypeSize, unsigned int mTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> mTimeSize > 0) {
    if (tTerm -> mTimeSize != mTimeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  mTimeSize has been previously defined:  %10d vs %10d", tTerm -> mTimeSize, mTimeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> mTimeSize = mTimeSize;
  }
  tTerm -> atRiskCount     = uivector(1, mTimeSize);
  tTerm -> eventCount      = uimatrix(1, eTypeSize, 1, mTimeSize);
}
void stackEventTimeIndex(Terminal *tTerm, unsigned int eTimeSize) {
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> eventTimeIndex  = uivector(1, eTimeSize + 1);
}
void unstackAtRiskAndEventCounts(Terminal *tTerm) {
  if (tTerm -> atRiskCount != NULL) {
    free_uivector(tTerm -> atRiskCount, 1, tTerm -> mTimeSize);
    tTerm -> atRiskCount = NULL;
  }
  if (tTerm -> eventCount != NULL) {
    free_uimatrix(tTerm -> eventCount, 1, tTerm -> eTypeSize, 1, tTerm -> mTimeSize);
    tTerm -> eventCount = NULL;
  }
}
void unstackEventTimeIndex(Terminal *tTerm) {
  if (tTerm -> eventTimeIndex != NULL) {
    free_uivector(tTerm -> eventTimeIndex, 1, tTerm -> eTimeSize + 1);
    tTerm -> eventTimeIndex = NULL;
  }
}
void stackLocalRatio(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> localRatio = dmatrix(1, eTypeSize, 1, tTerm -> eTimeSize);
}
void unstackLocalRatio(Terminal *tTerm) {
  if(tTerm -> eTimeSize > 0) {
    if (tTerm -> localRatio != NULL) {
      free_dmatrix(tTerm -> localRatio, 1, tTerm -> eTypeSize, 1, tTerm -> eTimeSize);
      tTerm -> localRatio = NULL;
    }
  }
}
void stackLocalSurvival(Terminal *tTerm, unsigned int eTimeSize) {
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> localSurvival = dvector(1, tTerm -> eTimeSize);
}
void unstackLocalSurvival(Terminal *tTerm) {
  if(tTerm -> eTimeSize > 0) {
    if (tTerm -> localSurvival != NULL) {
      free_dvector(tTerm -> localSurvival, 1, tTerm -> eTimeSize);
      tTerm -> localSurvival = NULL;
    }
  }
}
void stackLocalNelsonAalen(Terminal *tTerm, unsigned int eTimeSize) {
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> localNelsonAalen = dvector(1, tTerm -> eTimeSize);
}
void unstackLocalNelsonAalen(Terminal *tTerm) {
  if(tTerm -> eTimeSize > 0) {
    if (tTerm -> localNelsonAalen != NULL) {
      free_dvector(tTerm -> localNelsonAalen, 1, tTerm -> eTimeSize);
      tTerm -> localNelsonAalen = NULL;
    }
  }
}
void stackLocalCSH(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> localCSH = dmatrix(1, eTypeSize, 1, tTerm -> eTimeSize);
}
void unstackLocalCSH(Terminal *tTerm) {
  if(tTerm -> eTimeSize > 0) {
    if (tTerm -> localCSH != NULL) {
      free_dmatrix(tTerm -> localCSH, 1, tTerm -> eTypeSize, 1, tTerm -> eTimeSize);
      tTerm -> localCSH = NULL;
    }
  }
}
void stackLocalCIF(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> localCIF = dmatrix(1, eTypeSize, 1, tTerm -> eTimeSize);
}
void unstackLocalCIF(Terminal *tTerm) {
  if(tTerm -> eTimeSize > 0) {
    if (tTerm -> localCIF != NULL) {
      free_dmatrix(tTerm -> localCIF, 1, tTerm -> eTypeSize, 1, tTerm -> eTimeSize);
      tTerm -> localCIF = NULL;
    }
  }
}
void stackNelsonAalen(Terminal *tTerm, unsigned int sTimeSize) {
  if (tTerm -> sTimeSize > 0) {
    if (tTerm -> sTimeSize != sTimeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  sTimeSize has been previously defined:  %10d vs %10d", tTerm -> sTimeSize, sTimeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> sTimeSize = sTimeSize;
  }
  tTerm -> nelsonAalen = dvector(1, tTerm -> sTimeSize);
}
void unstackNelsonAalen(Terminal *tTerm) {
  if(tTerm -> sTimeSize > 0) {
    if (tTerm -> nelsonAalen != NULL) {
      free_dvector(tTerm -> nelsonAalen, 1, tTerm -> sTimeSize);
      tTerm -> nelsonAalen = NULL;
    }
  }
}
void stackSurvival(Terminal *tTerm, unsigned int sTimeSize) {
  if (tTerm -> sTimeSize > 0) {
    if (tTerm -> sTimeSize != sTimeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  sTimeSize has been previously defined:  %10d vs %10d", tTerm -> sTimeSize, sTimeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> sTimeSize = sTimeSize;
  }
  tTerm -> survival = dvector(1, tTerm -> sTimeSize);
}
void unstackSurvival(Terminal *tTerm) {
  if(tTerm -> sTimeSize > 0) {
    if (tTerm -> survival != NULL) {
      free_dvector(tTerm -> survival, 1, tTerm -> sTimeSize);
      tTerm -> survival = NULL;
    }
  }
}
void stackCSH(Terminal *tTerm, unsigned int eTypeSize, unsigned int sTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> sTimeSize > 0) {
    if (tTerm -> sTimeSize != sTimeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  sTimeSize has been previously defined:  %10d vs %10d", tTerm -> sTimeSize, sTimeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> sTimeSize = sTimeSize;
  }
  tTerm -> CSH = dmatrix(1, eTypeSize, 1, tTerm -> sTimeSize);
}
void unstackCSH(Terminal *tTerm) {
  if(tTerm -> eTypeSize > 0) {
    if(tTerm -> sTimeSize > 0) {
      if (tTerm -> CSH != NULL) {
        free_dmatrix(tTerm -> CSH, 1, tTerm -> eTypeSize, 1, tTerm -> sTimeSize);
        tTerm -> CSH = NULL;
      }
    }
  }
}
void stackCIF(Terminal *tTerm, unsigned int eTypeSize, unsigned int sTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> sTimeSize > 0) {
    if (tTerm -> sTimeSize != sTimeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  sTimeSize has been previously defined:  %10d vs %10d", tTerm -> sTimeSize, sTimeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> sTimeSize = sTimeSize;
  }
  tTerm -> CIF = dmatrix(1, eTypeSize, 1, tTerm -> sTimeSize);
}
void unstackCIF(Terminal *tTerm) {
  if(tTerm -> eTypeSize > 0) {
    if(tTerm -> sTimeSize > 0) {
      if (tTerm -> CIF != NULL) {
        free_dmatrix(tTerm -> CIF, 1, tTerm -> eTypeSize, 1, tTerm -> sTimeSize);
        tTerm -> CIF = NULL;
      }
    }
  }
}
void stackMortality(Terminal *tTerm, unsigned int eTypeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  tTerm -> mortality = dvector(1, eTypeSize);
}
void unstackMortality(Terminal *tTerm) {
  if(tTerm -> eTypeSize > 0) {
    if (tTerm -> mortality != NULL) {
      free_dvector(tTerm -> mortality, 1, tTerm -> eTypeSize);
      tTerm -> mortality = NULL;
    }
  }
}
void stackMultiClassProb(Terminal *tTerm, unsigned int rfCount, unsigned int *rfSize) {
  unsigned int j;
  if (tTerm -> rfCount > 0) {
    if (tTerm -> rfCount != rfCount) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  rfCount has been previously defined:  %10d vs %10d", tTerm -> rfCount, rfCount);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> rfCount = rfCount;
  }
  tTerm -> rfSize = uivector(1, tTerm -> rfCount);
  tTerm -> multiClassProb = (unsigned int **) new_vvector(1, tTerm -> rfCount, NRUTIL_UPTR);
  for (j = 1; j <= tTerm -> rfCount; j++) {
    (tTerm -> rfSize)[j] = rfSize[j];
    (tTerm -> multiClassProb)[j] = uivector(1, (tTerm -> rfSize)[j]);
  }
  tTerm -> maxClass = dvector(1, tTerm -> rfCount);
}
void stackMultiClassProbPartial(Terminal *tTerm, unsigned int rfCount) {
  if (tTerm -> rfCount > 0) {
    if (tTerm -> rfCount != rfCount) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  rfCount has been previously defined:  %10d vs %10d", tTerm -> rfCount, rfCount);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> rfCount = rfCount;
  }
  tTerm -> maxClass = dvector(1, tTerm -> rfCount);
}
void unstackMultiClassProb(Terminal *tTerm) {
  unsigned int j;
  if (tTerm -> rfCount > 0) {
    if (tTerm -> rfSize != NULL) {
      if (tTerm -> multiClassProb != NULL) {
        for (j = 1; j <= tTerm -> rfCount; j++) {
          if (tTerm -> multiClassProb[j] != NULL) {
            free_uivector(tTerm -> multiClassProb[j], 1, tTerm -> rfSize[j]);
            tTerm -> multiClassProb[j] = NULL;
          }
        }
        free_new_vvector(tTerm -> multiClassProb, 1, tTerm -> rfCount, NRUTIL_UPTR);
        tTerm -> multiClassProb = NULL;
      }
      free_uivector(tTerm -> rfSize, 1, tTerm -> rfCount);
      tTerm -> rfSize = NULL;
    }
  }
  if (tTerm -> rfCount > 0) {
    if (tTerm -> maxClass != NULL) {
      free_dvector(tTerm -> maxClass, 1, tTerm -> rfCount);
      tTerm -> maxClass = NULL;
    }
  }
}
void stackMeanResponse(Terminal *tTerm, unsigned int rnfCount) {
  if (tTerm -> rnfCount > 0) {
    if (tTerm -> rnfCount != rnfCount) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  rnfCount has been previously defined:  %10d vs %10d", tTerm -> rnfCount, rnfCount);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> rnfCount = rnfCount;
  }
  tTerm -> meanResponse = dvector(1, tTerm -> rnfCount);
}
void unstackMeanResponse(Terminal *tTerm) {
  if (tTerm -> rnfCount > 0) {
    if (tTerm -> meanResponse != NULL) {
      free_dvector(tTerm -> meanResponse, 1, tTerm -> rnfCount);
      tTerm -> meanResponse = NULL;
    }
  }
}
void stackTermLMIIndex(Terminal *tTerm, unsigned int size) {
  if (tTerm -> lmiAllocSize > 0) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  lmiIndex has been previously defined:  %10d vs %10d", tTerm -> lmiAllocSize, size);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tTerm -> lmiAllocSize = size;
    tTerm -> lmiSize = size;
  }
  tTerm -> lmiIndex = uivector(1, tTerm -> lmiAllocSize);
  tTerm -> lmiValue = dvector(1, tTerm -> lmiAllocSize);
}
void unstackTermLMIIndex(Terminal *tTerm) {
  if(tTerm -> lmiAllocSize > 0) {
    if (tTerm -> lmiIndex != NULL) {
      free_uivector(tTerm -> lmiIndex, 1, tTerm -> lmiAllocSize);
      free_dvector(tTerm -> lmiValue, 1, tTerm -> lmiAllocSize);
      tTerm -> lmiIndex = NULL;
      tTerm -> lmiValue = NULL;
      tTerm ->lmiAllocSize = 0;
      tTerm ->lmiSize = 0;
    }
  }
}
void getTerminalInfo(Terminal *termPtr) {
  Rprintf("\nTerminalInfo:  %20x", termPtr);
  Rprintf("\n  LeafCnt: %10d", termPtr -> nodeID);
  Rprintf("\n");
  Rprintf("\n lmiIndex            = %20x", termPtr -> lmiIndex);
  Rprintf("\n lmiSize             = %20d", termPtr -> lmiSize);
  Rprintf("\n lmiValue            = %20x", termPtr -> lmiValue);
}
Factor *makeFactor(uint r, char bookFlag) {
  uint i;
  Factor *f = (Factor*) gblock((size_t)sizeof(Factor));
  f -> r = r;
  f -> cardinalGroupCount = (uint) floor(r/2);
  f -> mwcpSize = (r >> (3 + ulog2(sizeof(uint)))) + ((r & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
  if (r > 1) {
    if (r <= MAX_EXACT_LEVEL) {
      f -> cardinalGroupSize = uivector(1, (f -> cardinalGroupCount) + 1);
      f -> complementaryPairCount =  ((uint*) (f -> cardinalGroupSize)) + (f -> cardinalGroupCount) + 1;
      *((uint*) f -> complementaryPairCount) = upower2(r-1) - 1;
    }
    else {
      f -> cardinalGroupSize = dvector(1, (f -> cardinalGroupCount) + 1);
      f -> complementaryPairCount =  ((double*) (f -> cardinalGroupSize)) + (f -> cardinalGroupCount) + 1;
      *((double*) f -> complementaryPairCount) = pow(2, r-1) - 1;
    }
    for (i=1; i <= f -> cardinalGroupCount; i++) {
      if (r <= MAX_EXACT_LEVEL) {
        nChooseK(r, i, EXACT, ((uint*) f -> cardinalGroupSize) + i);
      }
      else {
        nChooseK(r, i, APROX, ((double*) f -> cardinalGroupSize) + i);
      }
      f -> cardinalGroupBinary = NULL;
    }
    if (!((f -> r) & 0x01)) {
      if (r <= MAX_EXACT_LEVEL) {
        ((uint*) f -> cardinalGroupSize)[f -> cardinalGroupCount] = ((uint*) f -> cardinalGroupSize)[f -> cardinalGroupCount] >> 1;
      }
      else {
        ((double*) f -> cardinalGroupSize)[f -> cardinalGroupCount] = ((double*) f -> cardinalGroupSize)[f -> cardinalGroupCount] / 2;
      }
    }
    if (bookFlag && (r <= MAX_EXACT_LEVEL)) {
      bookFactor(f);
    }
  }  
  return f;
}
void free_Factor(Factor *f) {
  if (f -> r > 1) {
    unbookFactor(f);
    if (f -> r <= MAX_EXACT_LEVEL) {
      free_uivector(f -> cardinalGroupSize, 1, (f -> cardinalGroupCount) + 1);
    }
    else {
      free_dvector(f -> cardinalGroupSize, 1, (f -> cardinalGroupCount) + 1);
    }
  }
  free_gblock(f, sizeof(Factor));
}
char bookFactor(Factor *f) {
  uint i, j;
  uint row;
  char result;
  if (((f -> r) < 2) || ((f -> r) > MAX_EXACT_LEVEL)) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Minimum or Maximum number of factor levels violated in bookFactor(). ");
    RFprintf("\nRF-SRC:  Requested %10d, Minimum Allowed %10d, Maximum Allowed %10d ", f -> r, 2, MAX_EXACT_LEVEL);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit. \n");
  }
  if (f -> cardinalGroupBinary == NULL) {
    uint *leftLevel = uivector(1, f -> r);
    f -> cardinalGroupBinary = (uint **) new_vvector(1, f -> cardinalGroupCount, NRUTIL_UPTR);
    for (i=1; i <= f -> cardinalGroupCount; i++) {
      (f -> cardinalGroupBinary)[i] = uivector(1, ((uint*) f -> cardinalGroupSize)[i]);
      row = 0;
      for (j = 1; j <= i; j++) {
        leftLevel[j] = 0;
      }
      bookPair(f -> r , i, 1, &row, leftLevel, f);
    }
    free_uivector(leftLevel, 1, f -> r);
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  return result;
}
char unbookFactor(Factor *f) {
  char result;
  uint i;
  if (f -> cardinalGroupBinary != NULL) {
    for (i = 1; i <= f -> cardinalGroupCount; i++) {
      free_uivector((f -> cardinalGroupBinary)[i], 1, ((uint*) f -> cardinalGroupSize)[i]);
    }
    free_new_vvector(f -> cardinalGroupBinary, 1, f -> cardinalGroupCount, NRUTIL_UPTR);
    f -> cardinalGroupBinary = NULL;
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  return result;
}
void bookPair (uint   levelCount,
               uint    groupIndex,
               uint    setColumn,
               uint   *setRow,
               uint   *daughter,
               Factor *f) {
  uint i;
  daughter[setColumn] ++;
  if (setColumn < groupIndex) {
    setColumn ++;
    daughter[setColumn] ++;
    while (daughter[setColumn] < daughter[setColumn-1]) {
      daughter[setColumn] ++;
    }
    bookPair(levelCount, groupIndex, setColumn, setRow, daughter, f);
    daughter[setColumn] = 0;
    setColumn --;
    if ((*setRow) < ((uint*) (f -> cardinalGroupSize))[groupIndex]) {
      if (daughter[setColumn] < levelCount - (groupIndex - setColumn)) {
        bookPair(levelCount, groupIndex, setColumn, setRow, daughter, f);
      }
    }
  }
  else {
    (*setRow)++;
    (f -> cardinalGroupBinary)[groupIndex][*setRow] = 0;
    for (i=1; i <=groupIndex; i++) {
      (f -> cardinalGroupBinary)[groupIndex][*setRow] += upower(2, daughter[i] - 1);
    }
    if ( (levelCount > 2) && (daughter[setColumn] < levelCount)) {
      bookPair(levelCount, groupIndex, setColumn, setRow, daughter, f);
    }
  }
}
void nChooseK (uint n, uint r, char type, void *result) {
  if (type == EXACT) {
    uint total, multiplier, divisor, newMultiplier, newDivisor, k;
    total = 1;
    divisor = 1;
    multiplier = n;
    k = ((r < (n-r)) ? r : (n-r));
    while(divisor <= k) {
      newMultiplier = multiplier;
      newDivisor = divisor;
      reduceFraction(& total, & newDivisor);
      reduceFraction(& newMultiplier, & newDivisor);
      if (newMultiplier > (UINT_MAX / total)) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Arithmetic Overflow Encountered in nChooseK(n, k). ");
        RFprintf("\nRF-SRC:  Incoming parameters are (%10d, %10d). ", n, r);
        RFprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit. \n");
      }
      total = (total * newMultiplier) / newDivisor;
      multiplier--;
      divisor++;
    }
    *((uint*) result) = total;
  }
  else {
    double total, multiplier, divisor, k;
    total = 1;
    divisor = 1;
    multiplier = (double) n;
    k = (double) ((r < (n-r)) ? r : (n-r));
    while(divisor <= k) {
      total = (total * multiplier) / divisor;
      multiplier--;
      divisor++;
    }
    *((double*) result) = total;
  }
}
char reduceFraction(uint *numerator, uint *denominator) {
  uint numRemain, denRemain;
  char result;
  uint i;
  i = 2;
  result = FALSE;
  while (i <= *denominator) {
    numRemain = *numerator % i;
    if (numRemain == 0) {
      denRemain = *denominator % i;
      if (denRemain == 0) {
        *numerator = *numerator / i;
        *denominator = *denominator / i;
        result = TRUE;
      }
    }
    i++;
  }
  return result;
}
char splitOnFactor(uint level, uint *mwcp) {
  char daughterFlag;
  uint mwcpLevelIdentifier = (level >> (3 + ulog2(sizeof(uint)))) + ((level & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
  uint mwcpLevelWord = upower(2, level - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1 );
  daughterFlag = RIGHT;
  if (mwcpLevelWord & mwcp[mwcpLevelIdentifier]) {
    daughterFlag = LEFT;
  }
  return daughterFlag;
}
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
                   uint  ***densitySwap) {
  char validElement;
  uint i, j, k, kk;
  *sampleSize  = 0;
  *cdf         = NULL;
  *cdfSize     = 0;
  *cdfSort     = NULL;
  *density     = NULL;
  *densitySize = 0;
  *densitySwap = NULL;
  switch (weightType) {
  case RF_WGHT_UNIFORM:
    *index = uivector(1, permissibilitySize);
    if (permissibilityFlag != NULL) {
      *sampleSize = 0;
      for (k=1; k <= permissibilitySize; k++) {
        if (permissibilityFlag[k]) {
          (*index)[++(*sampleSize)] = k;
        }
      }
    }
    if (permissibilityIndex != NULL) {
      *sampleSize = permissibilitySize;
      for (k=1; k <= permissibilitySize; k++) {
        (*index)[k] = permissibilityIndex[k];
      }
    }
    break;
  case RF_WGHT_INTEGER:
    (*density) = uivector(1, maxDensitySize);
    (*densitySize) = 0;
    *densitySwap = (uint **) new_vvector(1, permissibilitySize, NRUTIL_UPTR);
    for (k = permissibilitySize; k >= 1; k--) {
      kk = weightSorted[k];
      validElement = TRUE;
      if (permissibilityFlag != NULL) {
        if (permissibilityFlag[kk] == FALSE) {
          validElement = FALSE;
        }
      }
      if (validElement) {
        j = (uint) weight[kk];
        if (j > 0) {
          (*densitySwap)[kk] = uivector(1, j);
          for (i = 1; i <= j; i++) {
            (*density)[++(*densitySize)] = kk;
            (*densitySwap)[kk][i] = (*densitySize);
          }
        }
        else {
          (*densitySwap)[kk] = NULL;
        }
      }
      else {
        (*densitySwap)[kk] = NULL;
      }
    }
    break;
  case RF_WGHT_GENERIC:
    *index = uivector(1, permissibilitySize);
    i      = 0;
    *cdf     = dvector(1, permissibilitySize);
    *cdfSort = uivector(1, permissibilitySize);
    *cdfSize = 0;
    for (k = 1; k <= permissibilitySize; k++) {
      kk = weightSorted[k];
      validElement = TRUE;
      if (permissibilityFlag != NULL) {
        if (permissibilityFlag[kk] == FALSE) {
          validElement = FALSE;
        }
      }
      if (validElement) {
        if (weightSorted[k] > 0) {
          (*index)[kk] = ++ i;
          (*cdfSize) ++;
          (*cdfSort)[(*cdfSize)] = kk;
          (*cdf)[(*cdfSize)] = weight[kk];
        }
        else {
          (*index)[kk] = 0;
        }
      }
      else {
        (*index)[kk] = 0;
      }
    }
    for (k = 2; k <= (*cdfSize); k++) {
      (*cdf)[k] += (*cdf)[k-1];
    }
    break;
  }
}
void updateCDF(uint    treeID,
               uint    weightType,
               double *weight,
               uint   *index,
               uint   *sampleSize,
               uint    sampleSlot,
               double *cdf,
               uint   *cdfSize,
               uint   *density,
               uint   *densitySize,
               uint  **densitySwap,
               uint    absoluteSlot) {
  double stepValue;
  uint sourcePt;
  uint stepIndex;
  uint currCov, nextCov;
  uint   i, j, k;
  switch (weightType) {
  case RF_WGHT_UNIFORM:
    index[sampleSlot] = index[(*sampleSize)];
    (*sampleSize) --;
    break;
  case RF_WGHT_INTEGER:
    currCov = nextCov = density[*densitySize];
    i = 0;
    j = (uint) weight[currCov];
    k = (uint) weight[absoluteSlot];
    while(i < k) {
      if (density[(*densitySize)] == absoluteSlot) {
        density[(*densitySize)] = 0;
        (*densitySize) --;
        densitySwap[absoluteSlot][k] = 0;
        k--;
        if (*densitySize > 0) {
          currCov = nextCov = density[(*densitySize)];
          j = (uint) weight[currCov];
        }
      }
      else {
        i++;
        sourcePt = densitySwap[absoluteSlot][i];
        density[sourcePt] = density[(*densitySize)];
        density[(*densitySize)] = 0;
        (*densitySize) --;
        densitySwap[currCov][j] = densitySwap[absoluteSlot][i];
        densitySwap[absoluteSlot][i] = 0;
        nextCov = density[(*densitySize)];
        if (nextCov == currCov) {
          j--;
        }
        else {
          hpsortui(densitySwap[currCov], (uint) weight[currCov]);
          currCov = nextCov = density[(*densitySize)];
          j = (uint) weight[currCov];
        }
      }
    }
    if (*densitySize > 0) {
      if (nextCov == currCov) {
        hpsortui(densitySwap[currCov], (uint) weight[currCov]);
      }
    }
    break;
  case RF_WGHT_GENERIC:
    stepIndex = index[absoluteSlot];
    stepValue = cdf[stepIndex];
    if (stepIndex > 1) {
      stepValue -= cdf[stepIndex-1];
    }
    for (k = stepIndex; k <= (*cdfSize); k++) {
      cdf[k] = cdf[k] - stepValue;
    }
    break;
  }
}
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
                    uint    densitySize) {
  double randomValue;
  uint low, mid, high, value;
  uint p;
  value = 0;  
  switch (weightType) {
  case RF_WGHT_UNIFORM:
    if (sampleSize > 0) {
      (*sampleSlot) = (uint) ceil(genericGenerator(treeID) * (sampleSize * 1.0));
      value = sampleIndex[(*sampleSlot)];
    }
    else {
      value = 0;
    }
    break;
  case RF_WGHT_INTEGER:
    if (densitySize > 0) {
      p = (uint) ceil(genericGenerator(treeID) * (densitySize * 1.0));
      value = density[p];
    }
    else {
      value = 0;
    }
    break;
  case RF_WGHT_GENERIC:
    if (cdf[cdfSize] > 0) {
      randomValue = genericGenerator(treeID) * cdf[cdfSize];
      low  = 1;
      high = cdfSize;
      while (low < high) {
        mid  = (low + high) >> 1;
        if (randomValue > cdf[mid]) {
          if (low == mid) {
            low = high;
          }
          else {
            low = mid;
          }
        }
        else {
          if (low == mid) {
            low = high;
          }
          else {
            high = mid;
          }
        }
      }
      value = cdfSort[high];
    }
    else {
      value = 0;
    }
    break;
  }
  return value;
}
void discardCDF(uint     treeID,
                uint     permissibilitySize,
                uint     weightType,
                double  *weight,
                uint     maxDensitySize,
                uint    *index,
                uint    *density,
                uint   **densitySwap,
                double  *cdf,
                uint    *cdfSort) {
  uint k;
  switch (weightType) {
  case RF_WGHT_UNIFORM:
    free_uivector(index, 1, permissibilitySize);
    break;
  case RF_WGHT_INTEGER:
    free_uivector(density, 1, maxDensitySize);
    for (k = 1; k <= permissibilitySize; k++) {
      if (densitySwap[k] != NULL) {
        free_uivector(densitySwap[k], 1, (uint) weight[k]);
        densitySwap[k] = NULL;
      }
    }
    free_new_vvector(densitySwap, 1, permissibilitySize, NRUTIL_UPTR);
    break;
  case RF_WGHT_GENERIC:
    free_uivector(index, 1, permissibilitySize);
    free_dvector(cdf, 1, permissibilitySize);
    free_uivector(cdfSort, 1, permissibilitySize);
    break;
  }
}
uint sampleUniformlyFromVector (uint    treeID,
                                uint   *index,
                                uint    size,
                                uint   *sampleSlot) {
  uint result;
  if (size > 0) {
    (*sampleSlot) = (uint) ceil(ran1B(treeID) * (size * 1.0));
    result = index[*sampleSlot];
  }
  else {
    result = 0;
  }
  return result;
}
SEXP rfsrcCIndex(SEXP sexp_traceFlag,
                 SEXP sexp_size,
                 SEXP sexp_time,
                 SEXP sexp_censoring,
                 SEXP sexp_predicted,
                 SEXP sexp_denom) {
  uint    traceFlag   = INTEGER(sexp_traceFlag)[0];
  int     size        = INTEGER(sexp_size)[0];
  double *time        = REAL(sexp_time); time--;
  double *censoring   = REAL(sexp_censoring); censoring--;
  double *predicted   = REAL(sexp_predicted); predicted--;
  uint   *denom       = (uint*) INTEGER(sexp_denom); denom--;
  double *v;
  char  *sexpString[3] = {
    "",              
    "",              
    "err"            
  };
  SEXP sexpVector[3];
  uint sexpIndex;
  uint stackCount;
  setUserTraceFlag(traceFlag);
  stackCount = 1;
  initProtect(sexpVector, stackCount);
  sexpIndex = 0;
  v = (double*) stackAndProtect(&sexpIndex, SEXP_TYPE_NUMERIC, 2, 1, sexpVector, sexpString);
  *v = getConcordanceIndex( 1,
                            size,
                            time,
                            censoring,
                            predicted,
                            denom);
  UNPROTECT(stackCount + 2);
  return sexpVector[0];
}
SEXP rfsrcTestSEXP(SEXP sexp_size) {
  ulong size = (ulong) REAL(sexp_size)[0];
  char  *sexpString[3] = {
    "",              
    "",              
    "dummy"          
  };
  SEXP sexpVector[3];
  uint sexpIndex;
  char *v;
  initProtect(sexpVector, 1);
  sexpIndex = 0;
  v = (char*) stackAndProtect(&sexpIndex, SEXP_TYPE_CHARACTER, 2, size, sexpVector, sexpString);
  v --;
  UNPROTECT(3);
  return sexpVector[0];
}
SEXP rfsrcGrow(SEXP traceFlag,
               SEXP seedPtr,
               SEXP opt,
               SEXP optHigh,
               SEXP splitRule,
               SEXP splitRandomCount,
               SEXP randomCovariateCount,
               SEXP randomResponseCount,
               SEXP minimumNodeSize,
               SEXP maximumNodeDepth,
               SEXP crWeight,
               SEXP forestSize,
               SEXP observationSize,
               SEXP rSize,
               SEXP rType,
               SEXP rLevels,
               SEXP rData,
               SEXP xSize,
               SEXP xType,
               SEXP xLevels,
               SEXP bootstrapSize,
               SEXP bootstrap,
               SEXP caseWeight,
               SEXP xWeight,
               SEXP xData,
               SEXP timeInterestSize,
               SEXP timeInterest,
               SEXP missTree,
               SEXP nImpute,
               SEXP numThreads) {
  int seedValue           = INTEGER(seedPtr)[0];
  RF_opt                  = INTEGER(opt)[0];
  RF_optHigh              = INTEGER(optHigh)[0];
  RF_splitRule            = INTEGER(splitRule)[0];
  RF_splitRandomCount     = INTEGER(splitRandomCount)[0];
  RF_randomCovariateCount = INTEGER(randomCovariateCount)[0];
  RF_randomResponseCount  = INTEGER(randomResponseCount)[0];
  RF_minimumNodeSize      = INTEGER(minimumNodeSize)[0];
  RF_maximumNodeDepth     = INTEGER(maximumNodeDepth)[0];
  RF_crWeight             = REAL(crWeight);  RF_crWeight--;
  RF_forestSize           = INTEGER(forestSize)[0];
  RF_observationSize      = INTEGER(observationSize)[0];
  RF_rSize                = INTEGER(rSize)[0];
  RF_sexp_rType           = rType;
  RF_rLevels              = INTEGER(rLevels); RF_rLevels--;
  RF_rData                = REAL(rData);
  RF_xSize                = INTEGER(xSize)[0];
  RF_sexp_xType           = xType;
  RF_xLevels              = INTEGER(xLevels); RF_xLevels--;
  RF_bootstrapSize        = INTEGER(bootstrapSize)[0];
  RF_bootstrap            = INTEGER(bootstrap);
  RF_caseWeight           = REAL(caseWeight);  RF_caseWeight--;
  RF_xWeight              = REAL(xWeight);  RF_xWeight--;
  RF_xData                = REAL(xData);
  RF_timeInterestSize     = INTEGER(timeInterestSize)[0];
  RF_timeInterest         = REAL(timeInterest);  RF_timeInterest--;
  RF_nImpute              = INTEGER(nImpute)[0];
  RF_numThreads           = INTEGER(numThreads)[0];
  RF_intrPredictorSize    = RF_xSize;
  RF_optHigh = RF_optHigh & (~OPT_MEMB_PRUN);
  RF_ptnCount             = 0;
  RF_sobservationSize = 0;
  RF_optHigh = RF_optHigh & (~OPT_PART_PLOT);
  RF_partialLength = 0;
  RF_opt = RF_opt & (~OPT_VIMP_JOIN);
  RF_opt = RF_opt & (~OPT_OUTC_TYPE);
  RF_opt = RF_opt & (~OPT_COMP_RISK);
  RF_optHigh = RF_optHigh & (~OPT_TERM_INCG);
  RF_optHigh = RF_optHigh & (~OPT_MEMB_INCG);
  RF_frSize = RF_fobservationSize = 0;
  if (RF_opt & OPT_IMPU_ONLY) {
    RF_opt                  = RF_opt & (OPT_IMPU_ONLY | OPT_BOOT_TYP1 | OPT_BOOT_TYP2);
    RF_optHigh              = RF_optHigh & (OPT_MISS_SKIP | OPT_MISS_MIA | OPT_MISS_MIAH | OPT_BOOT_SWOR);
    RF_opt                  = RF_opt | OPT_MISS;
    RF_opt = RF_opt | OPT_LEAF;
  }
  else {
    RF_opt                  = RF_opt | OPT_MISS;
    RF_opt = RF_opt | OPT_LEAF;
    RF_opt                  = RF_opt | OPT_FENS;
    RF_opt                  = RF_opt | OPT_OENS;
  }
  if ( (!(RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ||
       ( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2))) {
    RF_bootstrapSize = RF_observationSize;
    RF_optHigh = RF_optHigh & (~OPT_BOOT_SWOR);
    RF_opt                  = RF_opt & (~OPT_PERF);
    RF_opt                  = RF_opt & (~OPT_VIMP);
    RF_opt                  = RF_opt & (~OPT_OENS);
    if (RF_opt & OPT_PROX) {
      RF_opt = RF_opt | OPT_PROX_IBG;
      RF_opt = RF_opt | OPT_PROX_OOB;
    }
  }
  else {
  }
  if (RF_splitRule == USPV_SPLIT) {
    RF_opt                  = RF_opt & (~OPT_PERF);
    RF_opt                  = RF_opt & (~OPT_VIMP);
    RF_opt                  = RF_opt & (~OPT_OENS);
    RF_opt                  = RF_opt & (~OPT_FENS);
    if (RF_opt & OPT_NODE_STAT) {
      RF_opt = RF_opt | OPT_USPV_STAT;
    }
    RF_rSize = 0;
  }
  else {
    RF_opt = RF_opt & (~OPT_USPV_STAT);
  }
  if (RF_opt & OPT_PERF) {
  }
  else {
    RF_opt                  = RF_opt & (~OPT_VIMP);
  }
  if (RF_opt & OPT_TREE) {
    RF_opt = RF_opt | OPT_SEED;
  }
  else {
    RF_opt = RF_opt & (~OPT_SEED);
  }
  if ((RF_opt & OPT_OENS) | (RF_opt & OPT_FENS)) {
  }
  else {
    RF_optHigh = RF_optHigh & (~OPT_TERM_OUTG);
  }
  if (seedValue >= 0) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Parameter verification failed.");
    RFprintf("\nRF-SRC:  Random seed must be less than zero.  \n");
    RFprintf("\nRF-SRC:  The application will now exit.\n");
    return R_NilValue;
  }
  return rfsrc(RF_GROW, seedValue, INTEGER(traceFlag)[0]);
}
SEXP rfsrcPredict(SEXP traceFlag,
                  SEXP seedPtr,
                  SEXP opt,
                  SEXP optHigh,
                  SEXP forestSize,
                  SEXP observationSize,
                  SEXP rSize,
                  SEXP rType,
                  SEXP rTarget,
                  SEXP rTargetCount,
                  SEXP rLevels,
                  SEXP rData,
                  SEXP xSize,
                  SEXP xType,
                  SEXP xLevels,
                  SEXP xData,
                  SEXP bootstrapSize,
                  SEXP bootstrap,
                  SEXP caseWeight,
                  SEXP timeInterestSize,
                  SEXP timeInterest,
                  SEXP treeID,
                  SEXP nodeID,
                  SEXP parmID,
                  SEXP contPT,
                  SEXP mwcpSZ,
                  SEXP mwcpPT,
                  SEXP tnRMBR,
                  SEXP tnAMBR,
                  SEXP tnRCNT,
                  SEXP tnACNT,
                  SEXP totalNodeCount,
                  SEXP seed,
                  SEXP numThreads,
                  SEXP ptnCount,
                  SEXP intrPredictorSize,
                  SEXP intrPredictor,
                  SEXP sobservationSize,
                  SEXP sobservationIndv,
                  SEXP partialType,
                  SEXP partialXvar,
                  SEXP partialLength,
                  SEXP partialValue,
                  SEXP partialLength2,
                  SEXP partialXvar2,
                  SEXP partialValue2,
                  SEXP fobservationSize,
                  SEXP frSize,
                  SEXP frData,
                  SEXP fxData,
                  SEXP tnSURV,
                  SEXP tnMORT,
                  SEXP tnNLSN,
                  SEXP tnCSHZ,
                  SEXP tnCIFN,
                  SEXP tnREGR,
                  SEXP tnCLAS) {
  char mode;
  int seedValue           = INTEGER(seedPtr)[0];
  RF_opt                  = INTEGER(opt)[0];
  RF_optHigh              = INTEGER(optHigh)[0];
  RF_forestSize           = INTEGER(forestSize)[0];
  RF_observationSize      = INTEGER(observationSize)[0];
  RF_rSize                = INTEGER(rSize)[0];
  RF_sexp_rType           = rType;
  RF_rTarget              = (uint*) INTEGER(rTarget); RF_rTarget --;
  RF_rTargetCount         = INTEGER(rTargetCount)[0];
  RF_rLevels              = INTEGER(rLevels); RF_rLevels--;
  RF_rData                = REAL(rData);
  RF_xSize                = INTEGER(xSize)[0];
  RF_sexp_xType           = xType;
  RF_xLevels              = INTEGER(xLevels); RF_xLevels--;
  RF_xData                = REAL(xData);
  RF_bootstrapSize        = INTEGER(bootstrapSize)[0];
  RF_bootstrap            = INTEGER(bootstrap);
  RF_caseWeight           = REAL(caseWeight);  RF_caseWeight--;
  RF_timeInterestSize     = INTEGER(timeInterestSize)[0];
  RF_timeInterest         = REAL(timeInterest);  RF_timeInterest --;
  RF_treeID_              = (uint*) INTEGER(treeID);  RF_treeID_ --;
  RF_nodeID_              = (uint*) INTEGER(nodeID);  RF_nodeID_ --;
  RF_parmID_              = (uint*) INTEGER(parmID);  RF_parmID_ --;
  RF_contPT_              = REAL(contPT);  RF_contPT_ --;
  RF_mwcpSZ_              = (uint*) INTEGER(mwcpSZ);  RF_mwcpSZ_ --;
  RF_mwcpPT_              = (uint*) INTEGER(mwcpPT);  RF_mwcpPT_ --;
  RF_RMBR_ID_             = (uint*) INTEGER(tnRMBR);
  RF_AMBR_ID_             = (uint*) INTEGER(tnAMBR);
  RF_TN_RCNT_             = (uint*) INTEGER(tnRCNT);
  RF_TN_ACNT_             = (uint*) INTEGER(tnACNT);
  RF_totalNodeCount       = INTEGER(totalNodeCount)[0];
  RF_seed_                = INTEGER(seed); RF_seed_ --;
  RF_numThreads           = INTEGER(numThreads)[0];
  RF_ptnCount             = INTEGER(ptnCount)[0];
  RF_intrPredictorSize    = INTEGER(intrPredictorSize)[0];
  RF_intrPredictor        = (uint*) INTEGER(intrPredictor);  RF_intrPredictor --;
  RF_sobservationSize     = INTEGER(sobservationSize)[0];
  RF_sobservationIndv     = (uint *) INTEGER(sobservationIndv);  RF_sobservationIndv --;
  RF_partialType          = INTEGER(partialType)[0];
  RF_partialXvar          = INTEGER(partialXvar)[0];
  RF_partialLength        = INTEGER(partialLength)[0];
  RF_partialValue         = REAL(partialValue); RF_partialValue --;
  RF_partialLength2       = INTEGER(partialLength2)[0];
  RF_partialXvar2         = (uint *) INTEGER(partialXvar2); RF_partialXvar2 --;
  RF_partialValue2        = REAL(partialValue2); RF_partialValue2 --;
  RF_fobservationSize     = INTEGER(fobservationSize)[0];
  RF_frSize               = INTEGER(frSize)[0];
  RF_frData               = REAL(frData);
  RF_fxData               = REAL(fxData);
  RF_TN_SURV_ = REAL(tnSURV);
  RF_TN_MORT_ = REAL(tnMORT);
  RF_TN_NLSN_ = REAL(tnNLSN) ;
  RF_TN_CSHZ_ = REAL(tnCSHZ);
  RF_TN_CIFN_ = REAL(tnCIFN);
  RF_TN_REGR_ = REAL(tnREGR);
  RF_TN_CLAS_ = (uint*) INTEGER(tnCLAS);
  RF_opt = RF_opt & (~OPT_IMPU_ONLY);
  RF_opt = RF_opt & (~OPT_USPV_STAT);
  RF_opt = RF_opt & (~OPT_TREE);
  RF_opt = RF_opt & (~OPT_SEED);
  RF_nImpute = 1;
  RF_opt = RF_opt | OPT_LEAF;
  RF_opt  = RF_opt | OPT_MISS;
  RF_optHigh = RF_optHigh & (~OPT_MEMB_OUTG);
  RF_optHigh = RF_optHigh & (~OPT_TERM_OUTG);
  if(RF_fobservationSize > 0) {
    mode = RF_PRED;
  }
  else {
    mode = RF_REST;
  }
  switch (mode) {
  case RF_PRED:
    RF_sobservationSize = 0;
    RF_opt = RF_opt & (~OPT_OUTC_TYPE);
    RF_optHigh = RF_optHigh & (~OPT_PART_PLOT);
    RF_partialLength = RF_partialLength2 = 0;
    RF_opt = RF_opt & (~OPT_OENS);
    RF_opt = RF_opt | OPT_FENS;
    if (RF_rSize == 0) {
      RF_opt                  = RF_opt & (~OPT_PERF);
      RF_opt                  = RF_opt & (~OPT_VIMP);
      RF_opt                  = RF_opt & (~OPT_FENS);
    }
    else {
      if (RF_frSize == 0) {
        RF_opt                  = RF_opt & (~OPT_PERF);
        RF_opt                  = RF_opt & (~OPT_VIMP);
      }
    }
    if (RF_opt & OPT_PROX) {
      RF_opt = RF_opt | OPT_PROX_IBG;
      RF_opt = RF_opt | OPT_PROX_OOB;
    }
    break;
  case RF_REST:
    RF_frSize = RF_fobservationSize = 0;
    if (RF_sobservationSize > 0) {
      RF_opt = RF_opt & (~OPT_OUTC_TYPE);
      RF_optHigh = RF_optHigh & (~OPT_PART_PLOT);
      RF_partialLength = RF_partialLength2 = 0;
      RF_opt = RF_opt | OPT_OENS;
      RF_opt = RF_opt | OPT_FENS;
    }
    else if (RF_opt & OPT_OUTC_TYPE) {
      RF_sobservationSize = 0;
      RF_optHigh = RF_optHigh & (~OPT_PART_PLOT);
      RF_partialLength = RF_partialLength2 = 0;
      RF_optHigh = RF_optHigh & (~OPT_MEMB_INCG);
      RF_optHigh = RF_optHigh & (~OPT_TERM_INCG);
      RF_opt = RF_opt | OPT_OENS;
      RF_opt = RF_opt | OPT_FENS;
    }
    else if (RF_optHigh & OPT_PART_PLOT) {
      RF_sobservationSize = 0;
      RF_opt = RF_opt & (~OPT_OUTC_TYPE);
      RF_opt = RF_opt & (~OPT_PERF);
    }
    else {
      RF_opt = RF_opt | OPT_OENS;
      RF_opt = RF_opt | OPT_FENS;
    }
    if(RF_rSize == 0) {
      RF_opt = RF_opt & (~OPT_PERF);
      RF_opt = RF_opt & (~OPT_VIMP);
      RF_opt = RF_opt & (~OPT_OENS);
      RF_opt = RF_opt & (~OPT_FENS);
    }
    if ( (!(RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ||
         ( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2))) { 
      RF_opt                  = RF_opt & (~OPT_PERF);
      RF_opt                  = RF_opt & (~OPT_VIMP);
      RF_opt                  = RF_opt & (~OPT_OENS);
      if (RF_opt & OPT_PROX) {
        RF_opt = RF_opt | OPT_PROX_IBG;
        RF_opt = RF_opt | OPT_PROX_OOB;
      }
    }
    break;
  }
  if ( (!(RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ||
       ( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2))) { 
    RF_bootstrapSize = RF_observationSize;
    RF_optHigh = RF_optHigh & (~OPT_BOOT_SWOR);
    RF_optHigh = RF_optHigh & (~OPT_MEMB_INCG);
    RF_optHigh = RF_optHigh & (~OPT_TERM_INCG);
  }
  if (RF_ptnCount > 0) {
      RF_optHigh = RF_optHigh | OPT_MEMB_PRUN;
      RF_opt = RF_opt | OPT_NODE_STAT;
      RF_opt     = RF_opt & (~OPT_PERF);
      RF_opt     = RF_opt & (~OPT_VIMP);
      RF_opt     = RF_opt & (~OPT_PROX);
      RF_opt     = RF_opt & (~OPT_OENS);
      RF_opt     = RF_opt & (~OPT_FENS);
  }
  else {
    RF_optHigh = RF_optHigh & (~OPT_MEMB_PRUN);
  }
  if (RF_opt & OPT_PERF) {
  }
  else {
    RF_opt                  = RF_opt & (~OPT_VIMP);
  }
  if (seedValue >= 0) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Parameter verification failed.");
    RFprintf("\nRF-SRC:  User random seed must be less than zero.  \n");
    RFprintf("\nRF-SRC:  The application will now exit.\n");
    return R_NilValue;
  }
  return rfsrc(mode, seedValue, INTEGER(traceFlag)[0]);
}
char bootstrap (char     mode,
                uint     treeID,
                Node    *nodePtr,
                uint    *subsetIndex,
                uint     subsetSize,
                uint    *index,
                uint     indexSize) {
  char   *permissibility;
  uint   *caseIndex;
  uint    caseIndexSize;
  uint    caseIndexSlot;
  double *cdf;
  uint    cdfSize;
  uint   *cdfSort;
  uint   *density;
  uint    densitySize;
  uint  **densitySwap;
  char result;
  uint i, j, k;
  caseIndexSlot = 0;  
  result = TRUE;
  if (!(RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2)) {
    for (i=1; i <= subsetSize; i++) {
      index[i] = subsetIndex[i];
    }
  }
  else {
    if ( (RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2) ) {
      i = 0;
      for (k = 1; k <= RF_observationSize; k++) {
        for (j = 1; j <= RF_bootstrapIn[treeID][k]; j++) {
          index[++i] = k;
        }
      }
    }
    else {
      if ((RF_caseWeightType == RF_WGHT_UNIFORM) && !(RF_optHigh & OPT_BOOT_SWOR)) {
        for (i = 1; i <= indexSize; i++) {
          k = (uint) ceil(ran1A(treeID)*(subsetSize * 1.0));
          index[i] = subsetIndex[k];
        }
      }
      else {
        if (RF_caseWeightType != RF_WGHT_UNIFORM) {
          permissibility = cvector(1, RF_observationSize);
          for (i = 1; i <= RF_observationSize; i++) {
            permissibility[i] = FALSE;
          }
          for (i = 1; i <= subsetSize; i++) {
            permissibility[subsetIndex[i]] = TRUE;
          }
        }
        else {
          permissibility = NULL;
        }
        initializeCDF(treeID,
                      (RF_caseWeightType == RF_WGHT_UNIFORM) ? subsetIndex : NULL,
                      (RF_caseWeightType == RF_WGHT_UNIFORM) ? NULL : permissibility,
                      (RF_caseWeightType == RF_WGHT_UNIFORM) ? subsetSize : RF_observationSize,
                      RF_caseWeightType,
                      RF_caseWeight,
                      RF_caseWeightSorted,
                      RF_caseWeightDensitySize,
                      &caseIndex,
                      &caseIndexSize,
                      &cdf,
                      &cdfSize,
                      &cdfSort,
                      &density,
                      &densitySize,
                      &densitySwap);
        for (i = 1; i <= indexSize; i++) {
          index[i] = sampleFromCDF(ran1A,
                                   treeID,
                                   RF_caseWeightType,
                                   caseIndex,
                                   caseIndexSize,
                                   &caseIndexSlot,
                                   cdf,
                                   cdfSize,
                                   cdfSort,
                                   density,
                                   densitySize);
          if (RF_optHigh & OPT_BOOT_SWOR) {
            if (index[i] != 0) {
              updateCDF(treeID,
                        RF_caseWeightType,
                        RF_caseWeight,
                        caseIndex,
                        &caseIndexSize,
                        caseIndexSlot,
                        cdf,
                        &cdfSize,
                        density,
                        &densitySize,
                        densitySwap,
                        index[i]);
            }
            else {
              RFprintf("\nRF-SRC:  *** ERROR *** ");
              RFprintf("\nRF-SRC:  No cases left to select for bootstrap SWOR of size:  %10d", indexSize);
              RFprintf("\nRF-SRC:  Please Contact Technical Support.");
              error("\nRF-SRC:  The application will now exit.\n");
            }
          }
        }
        discardCDF(treeID,
                   (RF_caseWeightType == RF_WGHT_UNIFORM) ? subsetSize : RF_observationSize,
                   RF_caseWeightType,
                   RF_caseWeight,
                   RF_caseWeightDensitySize,
                   caseIndex,
                   density,
                   densitySwap,
                   cdf,
                   cdfSort);
        if (RF_caseWeightType != RF_WGHT_UNIFORM) {
          free_cvector(permissibility, 1, RF_observationSize);
        }
      }
    }
  }
  result = getNodeSign(mode, treeID, nodePtr, index, indexSize);
  if (result == FALSE) {
  }
  if (result == TRUE) {
    if (mode == RF_PRED) {
    }
  } 
  else {
  }
  return result;
}
char getNodeSign (char mode, uint treeID, Node *nodePtr, uint *bmIndex, uint repMembrSize) {
  int   *mvNSptr;
  int   *fmvNSptr;
  char result;
  uint i,p,q,m;
  result = TRUE;
  switch (mode) {
  case RF_PRED:
    if (RF_mRecordSize > 0) {
      stackMPSign(nodePtr, RF_mpIndexSize);
      mvNSptr = nodePtr -> mpSign;
    }
    else {
      mvNSptr = NULL;
    }
    if (RF_fmRecordSize > 0) {
      stackFMPSign(nodePtr, RF_fmpIndexSize);
      fmvNSptr = nodePtr -> fmpSign;
    }
    else {
      fmvNSptr = NULL;
    }
    break;
  default:
    if (RF_mRecordSize > 0) {
      stackMPSign(nodePtr, RF_mpIndexSize);
      mvNSptr = nodePtr -> mpSign;
    }
    else {
      mvNSptr = NULL;
    }
    fmvNSptr = NULL;
    break;
  }  
  if (mvNSptr != NULL) {
    int **mvBootstrapSign = imatrix(1, RF_mpIndexSize, 1, repMembrSize);
    for (p = 1; p <= RF_mpIndexSize; p++) {
      for (i = 1; i <= repMembrSize; i++) {
        mvBootstrapSign[p][i] = 0;
      }
    }
    for (p = 1; p <= RF_mpIndexSize; p++) {
      mvNSptr[p] = 0;
    }
    for (i=1; i <= repMembrSize; i++) {
      m = bmIndex[i];
      if (RF_mRecordMap[m] != 0) {
        for (p = 1; p <= RF_mpIndexSize; p++) {
          if (RF_mpIndex[p] < 0) {
            mvBootstrapSign[p][i] = RF_mpSign[(uint) abs(RF_mpIndex[p])][RF_mRecordMap[m]];
          }
          else {
            mvBootstrapSign[p][i] = RF_mpSign[RF_rSize + (uint) RF_mpIndex[p]][RF_mRecordMap[m]];
          }
        }
      }
      else {
        for (p = 1; p <= RF_mpIndexSize; p++) {
          mvBootstrapSign[p][i] = 0;
        }
      }
      for (p = 1; p <= RF_mpIndexSize; p++) {
        mvNSptr[p] = mvNSptr[p] + mvBootstrapSign[p][i];
      }
    }
    m = 0;
    for (p = 1; p <= RF_mpIndexSize; p++) {
      if (mvNSptr[p] > 0) {
        if (mvNSptr[p] == repMembrSize) {
          mvNSptr[p] = -1;
        }
        else {
          mvNSptr[p] = 1;
        }
      }
      if(RF_mpIndex[p] < 0) {
        if (mvNSptr[p] == -1) result = FALSE;
      }
      else {
        if (mvNSptr[p] == -1) m ++;
      }
    }  
    if (m == RF_mpIndexSize) {
      result = FALSE;
    }
    free_imatrix(mvBootstrapSign, 1, RF_mpIndexSize, 1, repMembrSize);
  }
  if (fmvNSptr != NULL) {
    for (p = 1; p <= RF_fmpIndexSize; p++) {
      fmvNSptr[p] = 1;
    }
    if (RF_mRecordSize > 0) {
      p = q = 1;
      while ((p <= RF_mpIndexSize) && (q <= RF_fmpIndexSize)) {
        if (RF_mpIndex[p] == RF_fmpIndex[q]) {
          if (mvNSptr[p] == -1) {
            fmvNSptr[q] = -1;
          }
          p++;
          q++;
        }
        else if (RF_fmpIndex[q] < 0) {
          if (RF_mpIndex[p] > 0) {
            q++;
          }
          else {
            if (abs(RF_fmpIndex[q]) < abs(RF_mpIndex[p])) {
              q++;
            }
            else {
              p++;
            }
          }
        }
        else {
          if (RF_fmpIndex[q] < RF_mpIndex[p]) {
            q++;
          }
          else {
            p++;
          }
        }
      }  
    }  
  }  
  return result;
}
void initializeTimeArrays(char mode) {
  uint i, j;
  uint leadingIndex;
  if (RF_timeIndex > 0) {
    RF_masterTimeSize = 0;
    for (j = 1; j <= RF_observationSize; j++) {
      if (!ISNA(RF_responseIn[RF_timeIndex][j])) {
        RF_masterTimeSize ++;
        RF_masterTime[RF_masterTimeSize] = RF_responseIn[RF_timeIndex][j];
      }
    }
    sort(RF_masterTime, RF_masterTimeSize);
    leadingIndex = 1;
    for (i=2; i <= RF_masterTimeSize; i++) {
      if (RF_masterTime[i] > RF_masterTime[leadingIndex]) {
        leadingIndex++;
        RF_masterTime[leadingIndex] = RF_masterTime[i];
      }
    }
    RF_masterTimeSize = leadingIndex;
    for (i= RF_masterTimeSize + 1; i <= RF_observationSize; i++) {
      RF_masterTime[i] = 0;
    }
    if (!(RF_opt & OPT_IMPU_ONLY)) {
      sort(RF_timeInterest, RF_timeInterestSize);
      RF_sortedTimeInterestSize = 1;
      for (i=2; i <= RF_timeInterestSize; i++) {
        if (RF_timeInterest[i] > RF_timeInterest[RF_sortedTimeInterestSize]) {
          (RF_sortedTimeInterestSize) ++;
          RF_timeInterest[RF_sortedTimeInterestSize] = RF_timeInterest[i];
        }
      }
      if (RF_sortedTimeInterestSize != RF_timeInterestSize) {
        RFprintf("\nRF-SRC:  *** WARNING *** ");
        RFprintf("\nRF-SRC:  Time points of interest are not unique.");
        RFprintf("\nRF-SRC:  The ensemble estimate output matrix is being");
        RFprintf("\nRF-SRC:  resized as [N'] x [n], where N' is the");
        RFprintf("\nRF-SRC:  unique time points of interest and n is");
        RFprintf("\nRF-SRC:  number of observations in the data.");
      }
      for (i = RF_sortedTimeInterestSize + 1; i <= RF_timeInterestSize; i++) {
        RF_timeInterest[i] = 0;
      }
    }
  }
}
void stackFactorArrays() {
  stackFactorGeneric(RF_rSize,
                     RF_rType,
                     &RF_rFactorMap,
                     &RF_rFactorCount,
                     &RF_rFactorIndex,
                     &RF_rFactorSize,
                     &RF_rNonFactorMap,
                     &RF_rNonFactorCount,
                     &RF_rNonFactorIndex);
  stackFactorGeneric(RF_xSize,
                     RF_xType,
                     &RF_xFactorMap,
                     &RF_xFactorCount,
                     &RF_xFactorIndex,
                     &RF_xFactorSize,
                     &RF_xNonFactorMap,
                     &RF_xNonFactorCount,
                     &RF_xNonFactorIndex);
}
void stackFactorGeneric(uint    size,
                        char  **type,
                        uint  **p_factorMap,
                        uint   *factorCount,
                        uint  **p_factorIndex,
                        uint  **p_factorSize,
                        uint  **p_nonfactorMap,
                        uint   *nonfactorCount,
                        uint  **p_nonfactorIndex) {
  uint i, j;
  if (size > 0) {
    *p_factorMap    = uivector(1, size);
    *p_nonfactorMap = uivector(1, size);
    *factorCount    = 0;
    *nonfactorCount = 0;
    for (i = 1; i <= size; i++) {
      (*p_factorMap)[i]    = 0;
      (*p_nonfactorMap)[i] = 0;
      if ((strcmp(type[i], "C") == 0) ||
          (strcmp(type[i], "I") == 0)) {
        (*factorCount) ++;
        (*p_factorMap)[i] = *factorCount;
      }
      else {
        (*nonfactorCount) ++;
        (*p_nonfactorMap)[i] = *nonfactorCount;
      }
    }
    if (*factorCount > 0) {
      *p_factorIndex = uivector(1, *factorCount);
      j = 0;
      for (i = 1; i <= size; i++) {
        if ((*p_factorMap)[i] > 0) {
          (*p_factorIndex)[++j] = i;
        }
      }
      *p_factorSize = uivector(1, *factorCount);
    }
    if (*nonfactorCount > 0) {
      *p_nonfactorIndex = uivector(1, *nonfactorCount);
      j = 0;
      for (i = 1; i <= size; i++) {
        if ((*p_nonfactorMap)[i] > 0) {
          (*p_nonfactorIndex)[++j] = i;
        }
      }
    }
  }
  else {
    *factorCount    = 0;
    *nonfactorCount = 0;
  }
}
void unstackFactorArrays() {
  uint j, k;
  if (RF_rSize > 0) {
    free_uivector(RF_rFactorMap, 1, RF_rSize);
    if (RF_rFactorCount > 0) {
      free_uivector(RF_rFactorIndex, 1, RF_rFactorCount);
      free_uivector(RF_rFactorSize, 1, RF_rFactorCount);
    }
    free_uivector(RF_rNonFactorMap, 1, RF_rSize);
    if (RF_rNonFactorCount > 0) {
      free_uivector(RF_rNonFactorIndex, 1, RF_rNonFactorCount);
    }
  }
  free_uivector(RF_xFactorMap, 1, RF_xSize);
  if (RF_xFactorCount > 0) {
    free_uivector(RF_xFactorIndex, 1, RF_xFactorCount);
    free_uivector(RF_xFactorSize, 1, RF_xFactorCount);
  }
  free_uivector(RF_xNonFactorMap, 1, RF_xSize);
  if (RF_xNonFactorCount > 0) {
    free_uivector(RF_xNonFactorIndex, 1, RF_xNonFactorCount);
  }
  if ((RF_rFactorCount + RF_xFactorCount) > 0) {
    for (j = 1; j <= RF_forestSize; j++) {
      if (RF_factorList[j] != NULL) {
        for (k = 1; k <= RF_maxFactorLevel; k++) {
          if (RF_factorList[j][k] != NULL) {
            free_Factor(RF_factorList[j][k]);
          }
        }
        free_new_vvector(RF_factorList[j], 1, RF_maxFactorLevel, NRUTIL_FPTR);
      }
    }
    free_new_vvector(RF_factorList, 1, RF_forestSize, NRUTIL_FPTR2);
  }
}
char stackMissingArrays(char mode) {
  char result;
  char mFlag;
  char dualUseFlag;
  uint recordSize;
  uint i, j;
  result = TRUE;
  for (j = 1; j <= RF_rSize; j++) {
    if (j == RF_timeIndex) {
      for (i = 1; i <= RF_observationSize; i++) {
        if (!ISNA(RF_responseIn[RF_timeIndex][i])) {
          if (RF_responseIn[RF_timeIndex][i] < 0) {
            result = FALSE;
            RFprintf("\nRF-SRC:  TRAINING time elements must be greater than or equal to zero or NA:  [%10d] = %12.4f \n", i, RF_responseIn[RF_timeIndex][i]);
          }
        }
      }
    }
    if (j == RF_statusIndex) {
      for (i = 1; i <= RF_observationSize; i++) {
        if (!ISNA(RF_responseIn[RF_statusIndex][i])) {
          if (RF_responseIn[RF_statusIndex][i] < 0) {
            result = FALSE;
            RFprintf("\nRF-SRC:  TRAINING status elements must be greater than or equal to zero or NA:  [%10d] = %12.4f \n", i, RF_responseIn[RF_statusIndex][i]);
          }
        }
      }
    }
    if (j == RF_statusIndex) {
      mFlag = FALSE;
      for (i = 1; i <= RF_observationSize; i++) {
        if (!ISNA(RF_responseIn[RF_statusIndex][i])) {
          if (RF_responseIn[RF_statusIndex][i] >= 0) {
            mFlag = TRUE;
            i = RF_observationSize;
          }
        }
      }
      if (mFlag == FALSE) {
        RFprintf("\nRF-SRC:  All TRAINING status elements are censored or missing. \n");
        result = FALSE;
      }
    }
    if ((RF_statusIndex == 0) && (RF_timeIndex == 0)) {
      mFlag = FALSE;
      for (i = 1; i <= RF_observationSize; i++) {
        if (!ISNA(RF_responseIn[j][i])) {
          mFlag = TRUE;
          i = RF_observationSize;
        }
      }
      if (mFlag == FALSE) {
        RFprintf("\nRF-SRC:  All TRAINING outcome/response elements are missing for:  %10d \n", j);
        result = FALSE;
      }
    }
    if (result == FALSE) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Missingness verification failed.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }  
  if (mode == RF_PRED) {
    if (RF_frSize > 0) {
      if (RF_timeIndex > 0) {
        for (i = 1 ; i <= RF_fobservationSize; i++) {
          if (!ISNA(RF_fresponseIn[RF_timeIndex][i])) {
            if (RF_fresponseIn[RF_timeIndex][i] < 0) {
              result = FALSE;
              RFprintf("\nRF-SRC:  PRED time elements must be greater than or equal to zero or NA:  [%10d] = %12.4f \n", i, RF_fresponseIn[RF_timeIndex][i]);
            }
          }
        }
      }
    }
    if (RF_frSize > 0) {
      if (RF_statusIndex > 0) {
        for (i = 1 ; i <= RF_fobservationSize; i++) {
          if (!ISNA(RF_fresponseIn[RF_statusIndex][i])) {
            if (RF_fresponseIn[RF_statusIndex][i] < 0) {
              result = FALSE;
              RFprintf("\nRF-SRC:  PRED status elements must be greater than or equal to zero or NA:  [%10d] = %12.4f \n", i, RF_fresponseIn[RF_statusIndex][i]);
            }
          }
        }
      }
    }
    if (result == FALSE) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Missingness verification failed.");
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }  
  RF_response = (double ***) new_vvector(1, RF_forestSize, NRUTIL_DPTR2);
  if (RF_rSize > 0) {
    for (i = 1 ; i <= RF_forestSize; i++) {
      RF_response[i] = RF_responseIn;
    }
    if (RF_timeIndex > 0) {
      RF_time = (double **) new_vvector(1, RF_forestSize, NRUTIL_DPTR);
      RF_masterTimeIndex = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
      for (i = 1 ; i <= RF_forestSize; i++) {
        RF_time[i] = RF_responseIn[RF_timeIndex];
        RF_masterTimeIndex[i] = RF_masterTimeIndexIn;
      }
      updateTimeIndexArray(0,
                           NULL,
                           RF_observationSize,
                           RF_responseIn[RF_timeIndex],
                           TRUE,
                           FALSE,
                           RF_masterTimeIndexIn);
    }
    if (RF_statusIndex > 0) {
      RF_status = (double **) new_vvector(1, RF_forestSize, NRUTIL_DPTR);
      for (i = 1 ; i <= RF_forestSize; i++) {
        RF_status[i] = RF_responseIn[RF_statusIndex];
      }
    }
  }
  else {
    for (i = 1 ; i <= RF_forestSize; i++) {
      RF_response[i] = NULL;
    }
  }
  RF_observation = (double ***) new_vvector(1, RF_forestSize, NRUTIL_DPTR2);
  for (i = 1 ; i <= RF_forestSize; i++) {
    RF_observation[i] = RF_observationIn;
  }
  RF_mRecordMap = uivector(1, RF_observationSize);
  RF_mRecordSize = getRecordMap(RF_mRecordMap,
                                RF_observationSize,
                                RF_responseIn,
                                RF_observationIn);
  if (RF_mRecordSize == 0) {
    RF_mStatusFlag = RF_mTimeFlag = RF_mResponseFlag = RF_mPredictorFlag = FALSE;
  }
  else {
    RF_optHigh = RF_optHigh & (~OPT_MEMB_INCG);
    RF_optHigh = RF_optHigh & (~OPT_TERM_INCG);
    stackMissingSignatures(RF_observationSize,
                           RF_rSize,
                           RF_responseIn,
                           RF_observationIn,
                           RF_mRecordMap,
                           RF_mRecordSize,
                           & RF_mRecordIndex,
                           & RF_mpIndexSize,
                           & RF_mpSign,
                           & RF_mpIndex,
                           & RF_mrFactorSize,
                           & RF_mrFactorIndex,
                           & RF_mxFactorSize,
                           & RF_mxFactorIndex,
                           & RF_mTimeFlag,
                           & RF_mStatusFlag,
                           & RF_mResponseFlag,
                           & RF_mPredictorFlag);
    if (RF_mResponseFlag == TRUE) {
      for (i = 1 ; i <= RF_forestSize; i++) {
        RF_response[i] = NULL;
        if (RF_timeIndex > 0) {
          RF_time[i] = NULL;
          RF_masterTimeIndex[i] = NULL;
        }
        if (RF_statusIndex > 0) {
          RF_status[i] = NULL;
        }
      }
    }
    if (RF_mPredictorFlag == TRUE) {
      for (i = 1 ; i <= RF_forestSize; i++) {
        RF_observation[i] = NULL;
      }
    }
  }  
  if (mode == RF_PRED) {
    RF_fobservation = (double ***) new_vvector(1, RF_forestSize, NRUTIL_DPTR2);
    for (i = 1 ; i <= RF_forestSize; i++) {
      RF_fobservation[i] = RF_fobservationIn;
    }
    RF_fmRecordMap = uivector(1, RF_fobservationSize);
    RF_fresponse = (double ***) new_vvector(1, RF_forestSize, NRUTIL_DPTR2);
    if (RF_frSize > 0) {
      for (i = 1 ; i <= RF_forestSize; i++) {
        RF_fresponse[i] = RF_fresponseIn;
      }
      if (RF_timeIndex > 0) {
        RF_ftime = (double **) new_vvector(1, RF_forestSize, NRUTIL_DPTR);
        for (i = 1 ; i <= RF_forestSize; i++) {
          RF_ftime[i] = RF_fresponseIn[RF_timeIndex];
        }
      }
      if (RF_statusIndex > 0) {
        RF_fstatus = (double **) new_vvector(1, RF_forestSize, NRUTIL_DPTR);
        for (i = 1 ; i <= RF_forestSize; i++) {
          RF_fstatus[i] = RF_fresponseIn[RF_statusIndex];
        }
      }
    }
    else {
      for (i = 1 ; i <= RF_forestSize; i++) {
        RF_fresponse[i] = NULL;
      }
    }
    RF_fmRecordSize = getRecordMap(RF_fmRecordMap,
                                 RF_fobservationSize,
                                 RF_fresponseIn,
                                 RF_fobservationIn);
    if (RF_fmRecordSize == 0) {
      RF_fmStatusFlag = RF_fmTimeFlag = RF_fmResponseFlag = RF_fmPredictorFlag = FALSE;
    }  
    else {
      RF_optHigh = RF_optHigh & (~OPT_MEMB_INCG);
      RF_optHigh = RF_optHigh & (~OPT_TERM_INCG);
      stackMissingSignatures(RF_fobservationSize,
                             RF_frSize,
                             RF_fresponseIn,
                             RF_fobservationIn,
                             RF_fmRecordMap,
                             RF_fmRecordSize,
                             & RF_fmRecordIndex,
                             & RF_fmpIndexSize,
                             & RF_fmpSign,
                             & RF_fmpIndex,
                             & RF_fmrFactorSize,
                             & RF_fmrFactorIndex,
                             & RF_fmxFactorSize,
                             & RF_fmxFactorIndex,
                             & RF_fmTimeFlag,
                             & RF_fmStatusFlag,
                             & RF_fmResponseFlag,
                             & RF_fmPredictorFlag);
      if (RF_frSize > 0) {
        if (RF_fmResponseFlag == TRUE) {
          for (i = 1 ; i <= RF_forestSize; i++) {
            RF_fresponse[i] = NULL;
            if (RF_timeIndex > 0) {
              RF_ftime[i] = NULL;
            }
            if (RF_statusIndex > 0) {
              RF_fstatus[i] = NULL;
            }
          }
        }
      }
      if (RF_fmPredictorFlag == TRUE) {
        for (i = 1 ; i <= RF_forestSize; i++) {
          RF_fobservation[i] = NULL;
        }
      }
    }  
  }  
  dualUseFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    if (RF_fmRecordSize > 0) {
      recordSize = RF_fmRecordSize;
      dualUseFlag = TRUE;
      mFlag = ACTIVE;
    }
    else {
      RF_opt = RF_opt & (~OPT_MISS);
    }
    break;
  default:
    RF_fmRecordSize = 0;
    if (RF_mRecordSize > 0) {
      recordSize = RF_mRecordSize;
      dualUseFlag = TRUE;
      mFlag = FALSE;
    }
    else {
      RF_opt = RF_opt & (~OPT_MISS);
      RF_nImpute = 1;
    }
    break;
  }  
  if (dualUseFlag == TRUE) {
    RF_dmRecordBootFlag = cmatrix(1, RF_forestSize, 1, recordSize);
    for (j = 1; j <= RF_forestSize; j++) {
      for (i = 1; i <= recordSize; i++) {
        RF_dmRecordBootFlag[j][i] = mFlag;
      }
    }
  }
  if (RF_rFactorCount + RF_xFactorCount > 0) {
    initializeFactorArrays(mode);
  }
  return result;
}
void unstackMissingArrays(char mode) {
  char dualUseFlag;
  uint recordSize;
  free_new_vvector(RF_response, 1, RF_forestSize, NRUTIL_DPTR2);
  if (RF_rSize > 0) {
    if (RF_timeIndex > 0) {
      free_new_vvector(RF_time, 1, RF_forestSize, NRUTIL_DPTR);
      free_new_vvector(RF_masterTimeIndex, 1, RF_forestSize, NRUTIL_UPTR);
    }
    if (RF_statusIndex > 0) {
      free_new_vvector(RF_status, 1, RF_forestSize, NRUTIL_DPTR);
    }
  }
  free_new_vvector(RF_observation, 1, RF_forestSize, NRUTIL_DPTR2);
  free_uivector(RF_mRecordMap, 1, RF_observationSize);
  if (RF_mRecordSize == 0) {
  }
  else {
    unstackMissingSignatures(RF_rSize,
                             RF_mRecordSize,
                             RF_mRecordIndex,
                             RF_mpIndexSize,
                             RF_mpSign,
                             RF_mpIndex,
                             RF_mrFactorSize,
                             RF_mrFactorIndex,
                             RF_mxFactorSize,
                             RF_mxFactorIndex);
  }
  if (mode == RF_PRED) {
    free_new_vvector(RF_fobservation, 1, RF_forestSize, NRUTIL_DPTR2);
    free_uivector(RF_fmRecordMap, 1, RF_fobservationSize);
    free_new_vvector(RF_fresponse, 1, RF_forestSize, NRUTIL_DPTR2);
    if (RF_frSize > 0) {
      if (RF_timeIndex > 0) {
        free_new_vvector(RF_ftime, 1, RF_forestSize, NRUTIL_DPTR);
      }
      if (RF_statusIndex > 0) {
        free_new_vvector(RF_fstatus, 1, RF_forestSize, NRUTIL_DPTR);
      }
    }
    if (RF_fmRecordSize == 0) {
    }
    else {
      unstackMissingSignatures(RF_frSize,
                               RF_fmRecordSize,
                               RF_fmRecordIndex,
                               RF_fmpIndexSize,
                               RF_fmpSign,
                               RF_fmpIndex,
                               RF_fmrFactorSize,
                               RF_fmrFactorIndex,
                               RF_fmxFactorSize,
                               RF_fmxFactorIndex);
    }
  }
  dualUseFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    if (RF_fmRecordSize > 0) {
      dualUseFlag = TRUE;
      recordSize = RF_fmRecordSize;
    }
    break;
  default:
    if (RF_mRecordSize > 0) {
      dualUseFlag = TRUE;
      recordSize = RF_mRecordSize;
    }
    break;
  }  
  if (dualUseFlag == TRUE) {
    free_cmatrix(RF_dmRecordBootFlag, 1, RF_forestSize, 1, recordSize);
  }
}
void stackMissingSignatures(uint     obsSize,
                            uint     rspSize,
                            double **responsePtr,
                            double **predictorPtr,
                            uint    *recordMap,
                            uint     recordSize,
                            uint   **p_recordIndex,
                            uint    *p_pIndexSize,
                            int   ***p_pSign,
                            int    **p_pIndex,
                            uint    *pRF_mrFactorSize,
                            uint   **pRF_mrFactorIndex,
                            uint    *pRF_mxFactorSize,
                            uint   **pRF_mxFactorIndex,
                            char    *pRF_mTimeFlag,
                            char    *pRF_mStatusFlag,
                            char    *pRF_mResponseFlag,
                            char    *pRF_mPredictorFlag) {
  uint i, j, p;
  if (recordSize < 1) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Attempt to allocate for missingness in its absence.");
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  *p_recordIndex = uivector(1, recordSize);
  i = 0;
  for (j = 1; j <= obsSize; j++) {
    if (recordMap[j] > 0) {
      i++;
      (*p_recordIndex)[i] = j;
    }
  }
  *p_pSign = imatrix(1, rspSize + RF_xSize, 1, recordSize);
  for (j = 1; j <= recordSize; j++) {
    for (i = 1; i <= rspSize + RF_xSize; i++) {
      (*p_pSign)[i][j] = 0;
    }
  }
  for (j = 1; j <= recordSize; j++) {
    for (i = 1; i <= rspSize; i++) {
      if (ISNA(responsePtr[i][(*p_recordIndex)[j]])) {
        (*p_pSign)[i][j] = 1;
      }
    }
    for (i = 1; i <= RF_xSize; i++) {
      if (ISNA(predictorPtr[i][(*p_recordIndex)[j]])) {
        (*p_pSign)[rspSize + i][j] = 1;
      }
    }
  }
  *pRF_mStatusFlag = *pRF_mTimeFlag = *pRF_mResponseFlag = *pRF_mPredictorFlag = FALSE;
  *p_pIndex = ivector(1, rspSize + RF_xSize);
  *p_pIndexSize = 0;
  for (i = 1; i <= rspSize; i++) {
    (*p_pIndex)[i] = 0;
    for (j = 1; j <= recordSize; j++) {
      if ((*p_pSign)[i][j] == 1) {
        (*p_pIndexSize) ++;
        (*p_pIndex)[*p_pIndexSize] = - i;
        *pRF_mResponseFlag = TRUE;
        if (i == RF_timeIndex) {
          *pRF_mTimeFlag = TRUE;
        }
        else if (i == RF_statusIndex) {
          *pRF_mStatusFlag = TRUE;
        }
        j = recordSize;
      }
    }
  }  
  for (i = rspSize + 1; i <= rspSize + RF_xSize; i++) {
    (*p_pIndex)[i] = 0;
    for (j = 1; j <= recordSize; j++) {
      if ((*p_pSign)[i][j] == 1) {
        (*p_pIndexSize) ++;
        (*p_pIndex)[*p_pIndexSize] =  i - rspSize;
        *pRF_mPredictorFlag = TRUE;
        j = recordSize;
      }
    }
  }  
  if (rspSize > 0) {
    *pRF_mrFactorIndex = uivector(1, rspSize);
    for (p = 1; p <= rspSize; p++) {
      (*pRF_mrFactorIndex)[p] = 0;
    }
  }
  *pRF_mxFactorIndex = uivector(1, RF_xSize);
  for (p = 1; p <= RF_xSize; p++) {
    (*pRF_mxFactorIndex)[p] = 0;
  }
  *pRF_mrFactorSize = *pRF_mxFactorSize = 0;
  for (p = 1; p <= *p_pIndexSize; p++) {
    if ((*p_pIndex)[p] < 0) {
      if ((strcmp(RF_rType[(uint) abs((*p_pIndex)[p])], "C") == 0) ||
          (strcmp(RF_rType[(uint) abs((*p_pIndex)[p])], "I") == 0)) {
        (*pRF_mrFactorSize) ++;
        (*pRF_mrFactorIndex)[*pRF_mrFactorSize] = (uint) abs((*p_pIndex)[p]);
      }
    }
    else {
      if ((strcmp(RF_xType[(*p_pIndex)[p]], "C") == 0) ||
          (strcmp(RF_xType[(*p_pIndex)[p]], "I") == 0)) {
        (*pRF_mxFactorSize) ++;
        (*pRF_mxFactorIndex)[*pRF_mxFactorSize] = (*p_pIndex)[p];
      }
    }
  }
}
void unstackMissingSignatures(uint      rspSize,
                              uint      recordSize,
                              uint     *recordIndex,
                              uint      vIndexSize,
                              int     **vSign,
                              int      *vIndex,
                              uint      mrFactorSize,
                              uint     *mrFactorIndex,
                              uint      mxFactorSize,
                              uint     *mxFactorIndex) {
  if (recordSize == 0) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Attempt to deallocate for missingness in its absence.");
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  free_uivector(recordIndex, 1, recordSize);
  free_imatrix(vSign, 1, rspSize + RF_xSize, 1, recordSize);
  free_ivector(vIndex, 1, rspSize + RF_xSize);
  if (rspSize > 0) {
    free_uivector(mrFactorIndex, 1, rspSize);
  }
  free_uivector(mxFactorIndex, 1, RF_xSize);
}
void initializeFactorArrays(char mode) {
  uint i, j;
  uint factorLevel;
  if (!(RF_rFactorCount + RF_xFactorCount > 0)) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Attempt to initialize factorness in its absence.");
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  RF_rMaxFactorLevel = 0;
  for (j = 1; j <= RF_rFactorCount; j++) {
    RF_rFactorSize[j] = RF_rLevels[RF_rFactorIndex[j]];
    for (i = 1; i <= RF_observationSize; i++) {
      if (!ISNA(RF_responseIn[RF_rFactorIndex[j]][i])) {
        if (RF_responseIn[RF_rFactorIndex[j]][i] >= 1) {
          factorLevel = (uint) RF_responseIn[RF_rFactorIndex[j]][i];
          if (RF_rLevels[RF_rFactorIndex[j]] < factorLevel) {
            RFprintf("\nRF-SRC:  *** ERROR *** ");
            RFprintf("\nRF-SRC:  Factor level in data inconsistent with number of levels indicated:  %10d %10d", factorLevel, RF_rLevels[RF_rFactorIndex[j]]);
            RFprintf("\nRF-SRC:  Please Contact Technical Support.");
            error("\nRF-SRC:  The application will now exit.\n");
          }
        }
        else {
          RFprintf("\nRF-SRC:  *** ERROR *** ");
          RFprintf("\nRF-SRC:  Factor level less than one (1):  %10.4f", RF_responseIn[RF_rFactorIndex[j]][i]);
          RFprintf("\nRF-SRC:  Please Contact Technical Support.");
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
    }
    if (RF_rMaxFactorLevel < RF_rFactorSize[j]) {
      RF_rMaxFactorLevel = RF_rFactorSize[j];
    }
  }
  RF_xMaxFactorLevel = 0;
  for (j = 1; j <= RF_xFactorCount; j++) {
    RF_xFactorSize[j] = RF_xLevels[RF_xFactorIndex[j]];
    for (i = 1; i <= RF_observationSize; i++) {
      if (!ISNA(RF_observationIn[RF_xFactorIndex[j]][i])) {
        if (RF_observationIn[RF_xFactorIndex[j]][i] >= 1) {
          factorLevel = (uint) RF_observationIn[RF_xFactorIndex[j]][i];
          if (RF_xLevels[RF_xFactorIndex[j]] < factorLevel) {
            RFprintf("\nRF-SRC:  *** ERROR *** ");
            RFprintf("\nRF-SRC:  Factor level in data inconsistent with number of levels indicated:  %10d %10d", factorLevel, RF_xLevels[RF_xFactorIndex[j]]);
            RFprintf("\nRF-SRC:  Please Contact Technical Support.");
            error("\nRF-SRC:  The application will now exit.\n");
          }
        }
        else {
          RFprintf("\nRF-SRC:  *** ERROR *** ");
          RFprintf("\nRF-SRC:  Factor level less than one (1):  %10.4f", RF_observationIn[RF_xFactorIndex[j]][i]);
          RFprintf("\nRF-SRC:  Please Contact Technical Support.");
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
    }
    if (RF_xMaxFactorLevel < RF_xFactorSize[j]) {
      RF_xMaxFactorLevel = RF_xFactorSize[j];
    }
  }
  RF_maxFactorLevel = (RF_xMaxFactorLevel > RF_rMaxFactorLevel) ? RF_xMaxFactorLevel : RF_rMaxFactorLevel;
  if (mode == RF_PRED) {
    if (RF_frSize > 0) {
      for (j = 1; j <= RF_rFactorCount; j++) {
        factorLevel = 0;
        for (i = 1; i <= RF_fobservationSize; i++) {
          if (!ISNA(RF_fresponseIn[RF_rFactorIndex[j]][i])) {
            if (RF_fresponseIn[RF_rFactorIndex[j]][i] >= 1) {
              factorLevel = (factorLevel > (uint) RF_fresponseIn[RF_rFactorIndex[j]][i]) ? factorLevel : ((uint) RF_fresponseIn[RF_rFactorIndex[j]][i]);
            }
            else {
              RFprintf("\nRF-SRC:  *** ERROR *** ");
              RFprintf("\nRF-SRC:  Factor level less than one (1):  %10.4f", RF_fobservationIn[RF_rFactorIndex[j]][i]);
              RFprintf("\nRF-SRC:  Please Contact Technical Support.");
              error("\nRF-SRC:  The application will now exit.\n");
            }
          }
        }
        if (factorLevel > RF_rFactorSize[j]) {
          RFprintf("\nRF-SRC:  *** ERROR *** ");
          RFprintf("\nRF-SRC:  !GROW factor level greater than maximum GROW factor level:  %10d vs. %10d", factorLevel, RF_rFactorSize[j]);
          RFprintf("\nRF-SRC:  Please Contact Technical Support.");
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
    }
    for (j = 1; j <= RF_xFactorCount; j++) {
      factorLevel = 0;
      for (i = 1; i <= RF_fobservationSize; i++) {
        if (!ISNA(RF_fobservationIn[RF_xFactorIndex[j]][i])) {
          if (RF_fobservationIn[RF_xFactorIndex[j]][i] >= 1) {
            factorLevel = (factorLevel > (uint) RF_fobservationIn[RF_xFactorIndex[j]][i]) ? factorLevel : ((uint) RF_fobservationIn[RF_xFactorIndex[j]][i]);
          }
          else {
            RFprintf("\nRF-SRC:  *** ERROR *** ");
            RFprintf("\nRF-SRC:  Factor level less than one (1):  %10.4f", RF_fobservationIn[RF_xFactorIndex[j]][i]);
            RFprintf("\nRF-SRC:  Please Contact Technical Support.");
            error("\nRF-SRC:  The application will now exit.\n");
          }
        }
      }
      if (factorLevel > RF_xFactorSize[j]) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  !GROW factor level greater than maximum GROW factor level:  %10d vs. %10d", factorLevel, RF_xFactorSize[j]);
        RFprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
  }
  RF_factorList = (Factor ***) new_vvector(1, RF_forestSize, NRUTIL_FPTR2);
  for (j = 1; j <= RF_forestSize; j++) {
    RF_factorList[j] = NULL;
  }
}
char stackCompetingArrays(char mode) {
  uint obsSize;
  double  *statusPtr;
  uint    *mRecordMap;
  int    **mpSign;
  char eventAnalysisFlag, eventSubsetFlag, consistencyFlag;
  char statusFlag;
  uint *eventCounter;
  uint *feventType;
  uint i, j, jgrow;
  if (RF_statusIndex == 0) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Attempt to stack competing risk structures in the absence of SURV data.");
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  RF_eventType = uivector(1, RF_observationSize);
  getEventTypeSize(RF_observationSize,
                   RF_responseIn[RF_statusIndex],
                   RF_mRecordMap,
                   RF_mpSign,
                   & RF_eventTypeSize,
                   & RF_mStatusSize,
                   RF_eventType);
  switch (mode) {
  case RF_GROW:
    if ((RF_splitRule == SURV_CR_LAU) || (RF_splitRule == SURV_CR_LOG)) {
      if (RF_eventTypeSize > 1) {
        RF_opt = RF_opt | OPT_COMP_RISK;
      }
      else {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Parameter verification failed.");
        RFprintf("\nRF-SRC:  Competing Risk analysis has been requested.");
        error("\nRF-SRC:  The training data set does not contain not contain competing risks.");
      }
    }
    else {
      if (RF_splitRule == CUST_SPLIT) {
        if (RF_eventTypeSize > 1) {
          RF_opt = RF_opt | OPT_COMP_RISK;
        }
        else {
          RF_opt = RF_opt & (~OPT_COMP_RISK);
        }
      }
      else {
        RF_opt = RF_opt & (~OPT_COMP_RISK);
      }
    }
    break;
  default:
    break;
  }
  if (RF_eventTypeSize == 0) {
    if ((RF_opt & OPT_OUTC_TYPE) && !(RF_opt & OPT_PERF) && !(RF_opt & OPT_VIMP)) {
      RF_opt                  = RF_opt & (~OPT_OENS);
      RF_opt                  = RF_opt & (~OPT_FENS);
    }
    else {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Parameter verification failed.");
      RFprintf("\nRF-SRC:  Performance or vimp has been requested.");
      error("\nRF-SRC:  The training or pseudo-training data set does not contain any events.");
    }
  }
  else {
    RF_eventTypeIndex  = uivector(1, RF_eventType[RF_eventTypeSize]);
    for (j = 1; j <= RF_eventType[RF_eventTypeSize]; j++) {
      RF_eventTypeIndex[j] = 0;
    }
    for (j = 1; j <= RF_eventTypeSize; j++) {
      RF_eventTypeIndex[RF_eventType[j]] = j;
    }
  }
  switch (mode) {
  case RF_GROW:
    if (RF_splitRule == RAND_SPLIT) {
      if (RF_eventTypeSize == 1) {
      }
      else {
        RF_opt = RF_opt | OPT_COMP_RISK;
      }
    }
    if ((RF_splitRule == SURV_CR_LAU) || (RF_splitRule == SURV_CR_LOG)) {
      if (RF_eventTypeSize == 1) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Split rule specified is for Competing Risk scenarios only.");
        error("\nRF-SRC:  The data set does not contain multiple events.");
      }
      i = 0;
      for (j = 1; j <= RF_eventTypeSize; j++) {
        if(fabs(RF_crWeight[j]) <= EPSILON) {
          i ++;
        }
        else {
          if(RF_crWeight[j] < 0.0) {
            RFprintf("\nRF-SRC:  *** ERROR *** ");
            RFprintf("\nRF-SRC:  Parameter verification failed.");
            RFprintf("\nRF-SRC:  Competing risk weight elements must be greater than or equal to zero:  %12.4f \n", RF_crWeight[j]);
            error("\nRF-SRC:  The application will now exit.\n");
          }
        }
      }
      if (i == RF_eventTypeSize) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Parameter verification failed.");
        RFprintf("\nRF-SRC:  Competing risk weight elements are all zero. \n");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
    break;
  default:
    if (RF_opt & OPT_COMP_RISK) {
      if (RF_eventTypeSize == 1) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  CR analysis has been specified in !GROW mode.");
        error("\nRF-SRC:  However, the GROW data set does not contain multiple events.");
      }
    }
    break;
  }
  switch (mode) {
  case RF_PRED:
    if (RF_frSize > 0) {
      eventAnalysisFlag = TRUE;
    }
    else {
      eventAnalysisFlag = FALSE;
    }
    break;
  default:
    eventAnalysisFlag = FALSE;
    break;
  } 
  if (eventAnalysisFlag == TRUE) {
    feventType = uivector(1, RF_fobservationSize);
    getEventTypeSize(RF_fobservationSize,
                     RF_fresponseIn[RF_statusIndex],
                     RF_fmRecordMap,
                     RF_fmpSign,
                     & RF_feventTypeSize,
                     & RF_mStatusSize,
                     feventType);
    if (RF_feventTypeSize == 0) {
      if (!(RF_opt & OPT_PERF) && !(RF_opt & OPT_VIMP)) {
      }
      else {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Parameter verification failed.");
        RFprintf("\nRF-SRC:  Performance or vimp has been requested.");
        error("\nRF-SRC:  The training or pseudo-training data set does not contain any events.");
      }
    }
    else {
      consistencyFlag = TRUE;
      if (RF_eventTypeSize > 1) {
        for (j = 1; j <= RF_feventTypeSize; j++) {
          for (jgrow = 1; jgrow <= RF_eventTypeSize; jgrow++) {
            if (feventType[j] != RF_eventType[jgrow]) {
              if (jgrow == RF_eventTypeSize) {
                consistencyFlag = FALSE;
              }
            }
            else {
              jgrow = RF_eventTypeSize;
            }
          }
        }
      }
      if (consistencyFlag == FALSE) {
        RFprintf("\nRF-SRC: *** ERROR *** ");
        RFprintf("\nRF-SRC: Unknown event type encountered in !GROW mode. ");
        error("\nRF-SRC: Please Contact Technical Support.");
      }
    }
    free_uivector(feventType, 1, RF_fobservationSize);
  }  
  if (RF_eventTypeSize > 1) {
    if (mode == RF_PRED) {
      if (RF_feventTypeSize > 0) {
        eventSubsetFlag = TRUE;
      }
      else {
        eventSubsetFlag = FALSE;
      }
    }
    else {
      eventSubsetFlag = TRUE;
    }
  }
  else {
    eventSubsetFlag = FALSE;
  }
  if (eventSubsetFlag == TRUE) {
    switch (mode) {
    case RF_PRED:
      obsSize = RF_fobservationSize;
      statusPtr = RF_fresponseIn[RF_statusIndex];
      mpSign = RF_fmpSign;
      mRecordMap = RF_fmRecordMap;
      break;
    default:
      obsSize = RF_observationSize;
      statusPtr = RF_responseIn[RF_statusIndex];
      mpSign = RF_mpSign;
      mRecordMap = RF_mRecordMap;
      break;
    }
    RF_eIndividualSize = uivector(1, RF_eventTypeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      RF_eIndividualSize[j] = 0;
    }
    for (i = 1; i <= obsSize; i++) {
      statusFlag = FALSE;
      if (mRecordMap[i] == 0) {
        statusFlag = TRUE;
      }
      else {
        if (mpSign[RF_statusIndex][mRecordMap[i]] == 0) {
          statusFlag = TRUE;
        }
      }
      if (statusFlag == TRUE) {
        if ((uint) statusPtr[i] > 0) {
          RF_eIndividualSize[RF_eventTypeIndex[(uint) statusPtr[i]]] ++;
        }
        else {
          for (j=1; j <= RF_eventTypeSize; j++) {
            RF_eIndividualSize[j] ++;
          }
        }
      }
    } 
    RF_eIndividualIn = (uint **) new_vvector(1, RF_eventTypeSize, NRUTIL_UPTR);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      RF_eIndividualIn[j] = uivector(1, RF_eIndividualSize[j] + RF_mStatusSize + 1);
    }
    eventCounter = uivector(1, RF_eventTypeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      eventCounter[j] = 0;
    }
    for (i = 1; i <= obsSize; i++) {
      statusFlag = FALSE;
      if (mRecordMap[i] == 0) {
        statusFlag = TRUE;
      }
      else {
        if (mpSign[RF_statusIndex][mRecordMap[i]] == 0) {
          statusFlag = TRUE;
        }
      }
      if (statusFlag == TRUE) {
        if ((uint) statusPtr[i] > 0) {
          j = RF_eventTypeIndex[(uint) statusPtr[i]];
          eventCounter[j] ++;
          RF_eIndividualIn[j][eventCounter[j]] = i;
        }
        else {
          for (j=1; j <= RF_eventTypeSize; j++) {
            eventCounter[j] ++;
            RF_eIndividualIn[j][eventCounter[j]] = i;
          }
        }
      }
    }
    free_uivector(eventCounter, 1, RF_eventTypeSize);
  }  
  return TRUE;
}
void getEventTypeSize(uint obsSize,
                      double *status,
                      uint *mRecordMap,
                      int **mpSign,
                      uint *eventTypeSize,
                      uint *msize,
                      uint *eventType) {
  uint statusFlag;
  uint leadingIndex;
  uint i;
  if (RF_statusIndex == 0) {
    RFprintf("\nRF-SRC: *** ERROR *** ");
    RFprintf("\nRF-SRC: Attempt to stack competing risk structures in the absence of SURV data.");
    RFprintf("\nRF-SRC: Please Contact Technical Support.");
    error("\nRF-SRC: The application will now exit.\n");
  }
  *eventTypeSize = *msize = 0;
  for (i = 1; i <= obsSize; i++) {
    eventType[i] = 0;
    statusFlag = FALSE;
    if (mRecordMap[i] == 0) {
      statusFlag = TRUE;
    }
    else {
      if (mpSign[RF_statusIndex][mRecordMap[i]] == 0) {
        statusFlag = TRUE;
      }
    }
    if (statusFlag == TRUE) {
      if ((uint) status[i] > 0) {
        (*eventTypeSize) ++;
        eventType[*eventTypeSize] = (uint) status[i];
      } 
      else {
      }
    }
    else {
      (*msize) ++;
    }
  }  
  if(*eventTypeSize > 0) {
    hpsortui(eventType, *eventTypeSize);
    leadingIndex = 1;
    for (i = 2; i <= *eventTypeSize; i++) {
      if (eventType[i] > eventType[leadingIndex]) {
        leadingIndex++;
        eventType[leadingIndex] = eventType[i];
      }
    }
    *eventTypeSize = leadingIndex;
  }
  for (i= *eventTypeSize + 1; i <= obsSize; i++) {
    eventType[i] = 0;
  }
}
void unstackCompetingArrays(char mode) {
  char eventSubsetFlag;
  uint j;
  if (RF_statusIndex == 0) {
    RFprintf("\nRF-SRC: *** ERROR *** ");
    RFprintf("\nRF-SRC: Attempt to unstack competing risk structures in the absence of SURV data.");
    RFprintf("\nRF-SRC: Please Contact Technical Support.");
    error("\nRF-SRC: The application will now exit.\n");
  }
  if (RF_eventTypeSize == 0) {
  }
  else {
    free_uivector(RF_eventTypeIndex, 1, RF_eventType[RF_eventTypeSize]);
  }
  free_uivector(RF_eventType, 1, RF_observationSize);
  if (RF_eventTypeSize > 1) {
    if (mode == RF_PRED) {
      if (RF_feventTypeSize > 0) {
        eventSubsetFlag = TRUE;
      }
      else {
        eventSubsetFlag = FALSE;
      }
    }
    else {
      eventSubsetFlag = TRUE;
    }
  }
  else {
    eventSubsetFlag = FALSE;
  }
  if (eventSubsetFlag == TRUE) {
    for (j = 1; j <= RF_eventTypeSize; j++) {
      free_uivector(RF_eIndividualIn[j], 1, RF_eIndividualSize[j] + RF_mStatusSize + 1);
    }
    free_new_vvector(RF_eIndividualIn, 1, RF_eventTypeSize, NRUTIL_UPTR);
    free_uivector(RF_eIndividualSize, 1, RF_eventTypeSize);
  }  
}
char stackClassificationArrays(char mode) {
  char classAnalysisFlag, consistencyFlag;
  uint j, k, jgrow;
  if (RF_rFactorCount == 0) {
    RFprintf("\nRF-SRC: *** ERROR *** ");
    RFprintf("\nRF-SRC: Attempt to stack classification structures in the absence of CLAS data.");
    RFprintf("\nRF-SRC: Please Contact Technical Support.");
    error("\nRF-SRC: The application will now exit.\n");
  }
  RF_classLevel = (uint **) new_vvector(1, RF_rFactorCount, NRUTIL_UPTR);
  RF_classLevelSize = uivector(1, RF_rFactorCount);
  getClassLevelSize(RF_observationSize,
                    RF_responseIn,
                    RF_mRecordMap,
                    RF_mpSign,
                    RF_classLevelSize,
                    RF_classLevel);
  RF_classLevelIndex = (uint **) new_vvector(1, RF_rFactorCount, NRUTIL_UPTR);
  for (k = 1; k <= RF_rFactorCount; k++) {
    RF_classLevelIndex[k] = uivector(1, RF_classLevel[k][RF_classLevelSize[k]]);
    for (j = 1; j <= RF_classLevel[k][RF_classLevelSize[k]]; j++) {
      RF_classLevelIndex[k][j] = 0;
    }
    for (j = 1; j <= RF_classLevelSize[k]; j++) {
      RF_classLevelIndex[k][RF_classLevel[k][j]] = j;
    }
  }  
  switch (mode) {
  case RF_PRED:
    if ((RF_opt & OPT_PERF) | (RF_opt & OPT_VIMP)) {
      classAnalysisFlag = TRUE;
    }
    else {
      classAnalysisFlag = FALSE;
    }
    break;
  default:
    classAnalysisFlag = FALSE;
    break;
  } 
  if (classAnalysisFlag == TRUE) {
    uint **fclassLevel = (uint **) new_vvector(1, RF_rFactorCount, NRUTIL_UPTR);
    uint *fclassLevelSize = uivector(1, RF_rFactorCount);
    getClassLevelSize(RF_fobservationSize,
                      RF_fresponseIn,
                      RF_fmRecordMap,
                      RF_fmpSign,
                      fclassLevelSize,
                      fclassLevel);
    consistencyFlag = TRUE;
    for (j = 1; j <= RF_rFactorCount; j++) {
      for (k = 1; k <= fclassLevelSize[j]; k++) {
        for (jgrow = 1; jgrow <= RF_classLevelSize[j]; jgrow++) {
          if (fclassLevel[j][k] != RF_classLevel[j][jgrow]) {
            if (jgrow == RF_classLevelSize[j]) {
              consistencyFlag = FALSE;
            }
          }
          else {
            jgrow = RF_classLevelSize[j];
          }
        }
      }
    }
    for (j = 1; j <= RF_rFactorCount; j ++) {
      free_uivector(fclassLevel[j], 1, fclassLevelSize[j]);
    }
    free_new_vvector(fclassLevel, 1, RF_rFactorCount, NRUTIL_UPTR);
    free_uivector(fclassLevelSize, 1, RF_rFactorCount);
    if (consistencyFlag == FALSE) {
    }
  }  
  return TRUE;
}
void getClassLevelSize(uint      obsSize,
                       double  **response,
                       uint     *mRecordMap,
                       int     **mpSign,
                       uint     *classLevelSize,
                       uint    **classLevel) {
  uint *rawClassVector;
  uint classFlag;
  uint leadingIndex;
  uint i, j, k;
  if (RF_rFactorCount == 0) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Attempt to stack classification response structures in the absence of CLAS data.");
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  rawClassVector = uivector(1, obsSize);
  for (j = 1; j <= RF_rFactorCount; j++) {
    classLevelSize[j] = 0;
  }
  for (j = 1; j <= RF_rFactorCount; j++) {
    for (i = 1; i <= obsSize; i++) {
      classFlag = FALSE;
      if (mRecordMap[i] == 0) {
        classFlag = TRUE;
      }
      else {
        if (mpSign[RF_rFactorIndex[j]][mRecordMap[i]] == 0) {
          classFlag = TRUE;
        }
      }
      if (classFlag == TRUE) {
        classLevelSize[j] ++;
        rawClassVector[classLevelSize[j]] = (uint) response[RF_rFactorIndex[j]][i];
      }
      else {
      }
    }  
    hpsortui(rawClassVector, classLevelSize[j]);
    leadingIndex = 1;
    for (k=2; k <= classLevelSize[j]; k++) {
      if (rawClassVector[k] > rawClassVector[leadingIndex]) {
        leadingIndex++;
        rawClassVector[leadingIndex] = rawClassVector[k];
      }
    }
    classLevelSize[j] = leadingIndex;
    classLevel[j] = uivector(1, classLevelSize[j]);
    for (k=1; k <= classLevelSize[j]; k++) {
      classLevel[j][k] = rawClassVector[k];
    }
  } 
  free_uivector(rawClassVector, 1, obsSize);
}
void unstackClassificationArrays(char mode) {
  uint j;
  if (RF_rFactorCount == 0) {
    RFprintf("\nRF-SRC: *** ERROR *** ");
    RFprintf("\nRF-SRC: Attempt to unstack classification structures in the absence of CLAS data.");
    RFprintf("\nRF-SRC: Please Contact Technical Support.");
    error("\nRF-SRC: The application will now exit.\n");
  }
  for (j = 1; j <= RF_rFactorCount; j++) {
    free_uivector(RF_classLevelIndex[j], 1, RF_classLevel[j][RF_classLevelSize[j]]);
  }
  free_new_vvector(RF_classLevelIndex, 1, RF_rFactorCount, NRUTIL_UPTR);
  for (j = 1; j <= RF_rFactorCount; j ++) {
    free_uivector(RF_classLevel[j], 1, RF_classLevelSize[j]);
  }
  free_new_vvector(RF_classLevel, 1, RF_rFactorCount, NRUTIL_UPTR);
  free_uivector(RF_classLevelSize, 1, RF_rFactorCount);
}
void stackIncomingResponseArrays(char mode) {
  uint i, j;
  RF_timeIndex = RF_statusIndex = 0;
  if (RF_rSize > 0) {
    RF_rType = (char **) new_vvector(1, RF_rSize, NRUTIL_CPTR);
    RF_yIndex = uivector(1, RF_rSize);
    j = 0;
    for (i = 1; i <= RF_rSize; i++) {
      RF_rType[i] = (char*) CHAR(STRING_ELT(AS_CHARACTER(RF_sexp_rType), i-1));
      if ((strcmp(RF_rType[i], "C") != 0) &&
          (strcmp(RF_rType[i], "I") != 0) &&
          (strcmp(RF_rType[i], "R") != 0) &&
          (strcmp(RF_rType[i], "T") != 0) &&
          (strcmp(RF_rType[i], "S") != 0)) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Invalid type:  [%10d] = %2s", i, RF_rType[i]);
        RFprintf("\nRF-SRC:  Outcomes must be 'C', 'I', 'R', 'T', or 'S'.");
        RFprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
      RF_yIndex[i] = 0;
      if (strcmp(RF_rType[i], "T") == 0) {
        RF_timeIndex = i;
      }
      else if (strcmp(RF_rType[i], "S") == 0) {
        RF_statusIndex = i;
      }
      else {
        RF_yIndex[++j] = i;
      }
    }
    if (mode == RF_PRED) {
      if (RF_frSize > 0) {
        if (RF_rSize != RF_frSize) {
          RFprintf("\nRF-SRC:  *** ERROR *** ");
          RFprintf("\nRF-SRC:  TRAIN and TEST outcome/response matrices must be of the same dimension.  ");
          RFprintf("\nRF-SRC:  TRAIN vs TEST:  %10d vs %10d  ", RF_rSize, RF_frSize);
          RFprintf("\nRF-SRC:  Please Contact Technical Support.");
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
      else {
        if ((RF_opt & OPT_PERF) | (RF_opt & OPT_VIMP)) {
          RFprintf("\nRF-SRC:  *** ERROR *** ");
          RFprintf("\nRF-SRC:  TEST outcome/response matrix must be present when PERF or VIMP is requested.  ");
          RFprintf("\nRF-SRC:  Please Contact Technical Support.");
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
    }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      RF_ySize = 0;
      RF_ptnCount = 0;
    }
    else {
      RF_ySize = RF_rSize - ((RF_timeIndex == 0) ? 0:1) - ((RF_statusIndex == 0) ? 0:1);
      if (RF_ySize != RF_rSize) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Responses must be [S], [C], [R], [R+], [C+], [M+].  ");
        RFprintf("\nRF-SRC:  The application will now exit.\n");
        error("\nRF-SRC:  The application will now exit.\n");
      }
      if (RF_ySize == 1) {
      }
      else {
      }
    }
    RF_responseIn = (double **) new_vvector(1, RF_rSize, NRUTIL_DPTR);
    for (i=1; i <= RF_rSize; i++) {
      RF_responseIn[i] = (RF_rData + ((i-1) * RF_observationSize) - 1);
    }
    if (mode == RF_PRED) {
      if (RF_frSize > 0) {
        RF_fresponseIn = (double **) new_vvector(1, RF_frSize, NRUTIL_DPTR);
        for (i=1; i <= RF_rSize; i++) {
          RF_fresponseIn[i] = (RF_frData + ((i-1) * RF_fobservationSize) - 1);
        }
      }
      else {
        RF_fresponseIn = NULL;
      }
    }
  }
  else {
    RF_rType      = NULL;
    RF_responseIn = NULL;
    RF_ySize = 0;
  }
}
void unstackIncomingResponseArrays(char mode) {
  if (RF_rSize > 0) {
    free_new_vvector(RF_rType, 1, RF_rSize, NRUTIL_CPTR);
    free_uivector(RF_yIndex, 1, RF_rSize);
    free_new_vvector(RF_responseIn, 1, RF_rSize, NRUTIL_DPTR);
    if (mode == RF_PRED) {
      if (RF_frSize > 0) {
        free_new_vvector(RF_fresponseIn, 1, RF_frSize, NRUTIL_DPTR);
      }
    }
  }
}
void stackIncomingCovariateArrays(char mode) {
  uint i;
  RF_xType = (char **) new_vvector(1, RF_xSize, NRUTIL_CPTR);
  for (i = 1; i <= RF_xSize; i++) {
    RF_xType[i] = (char*) CHAR(STRING_ELT(AS_CHARACTER(RF_sexp_xType), i-1));
    if ((strcmp(RF_xType[i], "C") != 0) &&
        (strcmp(RF_xType[i], "I") != 0) &&
        (strcmp(RF_xType[i], "R") != 0)) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Invalid type:  [%10d] = %2s", i, RF_xType[i]);
      RFprintf("\nRF-SRC:  Predictors must be 'C', 'I', or 'R'.");
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  RF_observationIn = (double **) new_vvector(1, RF_xSize, NRUTIL_DPTR);
  for (i=1; i <= RF_xSize; i++) {
    RF_observationIn[i] = (RF_xData + ((i-1) * RF_observationSize) - 1);
  }
  if (mode == RF_PRED) {
    RF_fobservationIn = (double **) new_vvector(1, RF_xSize, NRUTIL_DPTR);
    for (i=1; i <= RF_xSize; i++) {
      RF_fobservationIn[i] = (RF_fxData + ((i-1) * RF_fobservationSize) - 1);
    }
  }
}
void unstackIncomingCovariateArrays(char mode) {
  free_new_vvector(RF_xType, 1, RF_xSize, NRUTIL_CPTR);
  free_new_vvector(RF_observationIn, 1, RF_xSize, NRUTIL_DPTR);
  if (mode == RF_PRED) {
    free_new_vvector(RF_fobservationIn, 1, RF_xSize, NRUTIL_DPTR);
  }
}
void stackIncomingArrays(char mode) {
  stackIncomingResponseArrays(mode);
  stackIncomingCovariateArrays(mode);
  if (mode == RF_GROW) {
    if (RF_minimumNodeSize < 1) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Parameter verification failed.");
      RFprintf("\nRF-SRC:  Minimum node size must be greater than zero:  %10d \n", RF_minimumNodeSize);
      error("\nRF-SRC:  The application will now exit.\n");
    }
    if (RF_bootstrapSize < 1) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Parameter verification failed.");
      RFprintf("\nRF-SRC:  Bootstrap size must be greater than zero:  %12d \n", RF_bootstrapSize);
      error("\nRF-SRC:  The application will now exit.\n");
    }
    if ( RF_splitRule > MAXM_SPLIT) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Parameter verification failed.");
      RFprintf("\nRF-SRC:  Invalid split rule:  %10d \n", RF_splitRule);
      error("\nRF-SRC:  The application will now exit.\n");
    }
    if (!(RF_splitRule == USPV_SPLIT)) {
      if (RF_rSize == 0) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Parameter verification failed.");
        RFprintf("\nRF-SRC:  Number of response variables must be greater than zero:  %10d \n", RF_rSize);
        error("\nRF-SRC:  The application will now exit.\n");
      }
      if ( ((RF_randomCovariateCount < 1) || (RF_randomCovariateCount > RF_xSize)) ) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Parameter verification failed.");
        RFprintf("\nRF-SRC:  Number of random covariate parameters must be greater");
        RFprintf("\nRF-SRC:  than zero and less than or equal to the total number of covariates:  %10d \n", RF_randomCovariateCount);
        error("\nRF-SRC:  The application will now exit.\n");
      }
      if (RF_ySize > 0) {
        if ( (RF_randomResponseCount < 1) || (RF_randomResponseCount > RF_ySize) ) {
          RFprintf("\nRF-SRC:  *** ERROR *** ");
          RFprintf("\nRF-SRC:  Parameter verification failed.");
          RFprintf("\nRF-SRC:  ytry must be within range:  %10d \n", RF_randomResponseCount);
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
    }
    else {
      if ( RF_xSize < 2) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Parameter verification failed.");
        RFprintf("\nRF-SRC:  Number of covariates must be greater than or equal to two (2) with specified split rule:  %10d \n", RF_xSize);
        error("\nRF-SRC:  The application will now exit.\n");
      }
      if ( ((int) (RF_xSize - RF_randomResponseCount) < 1) || (RF_randomCovariateCount > RF_xSize) ) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Parameter verification failed.");
        RFprintf("\nRF-SRC:  ytry and mtry must be within range:  %10d %10d \n", RF_randomResponseCount,  RF_randomCovariateCount);
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
    for (uint i = 1; i <= RF_xSize; i++) {
      if(RF_xWeight[i] < 0) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Parameter verification failed.");
        RFprintf("\nRF-SRC:  X-weight elements must be greater than or equal to zero:  %12.4f \n", RF_xWeight[i]);
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
    if ((RF_timeIndex == 0) && (RF_statusIndex == 0)) {
      if ((RF_splitRule != RAND_SPLIT)  &&
          (RF_splitRule != REGR_WT_NRM) &&
          (RF_splitRule != REGR_WT_OFF) &&
          (RF_splitRule != REGR_WT_HVY) &&
          (RF_splitRule != CLAS_WT_NRM) &&
          (RF_splitRule != CLAS_WT_OFF) &&
          (RF_splitRule != CLAS_WT_HVY) &&
          (RF_splitRule != MVRG_SPLIT)  &&
          (RF_splitRule != MVCL_SPLIT)  &&
          (RF_splitRule != USPV_SPLIT)  &&
          (RF_splitRule != CUST_SPLIT)) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  !SURV data and split rule specified are incompatible.");
        RFprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
    else if ((RF_timeIndex != 0) && (RF_statusIndex != 0)) {
      if ((RF_splitRule != SURV_LGRNK)  &&
          (RF_splitRule != SURV_LRSCR)  &&
          (RF_splitRule != SURV_L2IMP)  &&
          (RF_splitRule != SURV_CR_LOG) &&
          (RF_splitRule != SURV_CR_LAU) &&
          (RF_splitRule != RAND_SPLIT)  &&
          (RF_splitRule != CUST_SPLIT)) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  SURV data and split rule specified are incompatible.");
        RFprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
    else {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Data set contains mixed outcomes with no comatible split rule.");
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
}
void unstackIncomingArrays(char mode) {
  unstackIncomingResponseArrays(mode);
  unstackIncomingCovariateArrays(mode);
}
void stackPreDefinedCommonArrays() {
  uint maxSize;
  uint i, j, k;
  RF_nodeMembership = (Node ***)     new_vvector(1, RF_forestSize, NRUTIL_NPTR2);
  RF_tTermMembership = (Terminal ***) new_vvector(1, RF_forestSize, NRUTIL_NPTR2);
  RF_tNodeList = (Node ***)     new_vvector(1, RF_forestSize, NRUTIL_NPTR2);
  RF_tNodeListLength = uivector(1, RF_forestSize);
  for (i = 1; i <= RF_forestSize; i++) {
    RF_tNodeListLength[i] = 0;
  }
  RF_tTermList = (Terminal ***) new_vvector(1, RF_forestSize, NRUTIL_NPTR2);
  RF_bootMembershipIndex = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
  if ( (RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2) ) {
    RF_bootstrapIn = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
    for (i = 1; i <= RF_forestSize; i++) {
      RF_bootstrapIn[i] = (uint *) ((RF_bootstrap) + ((i-1) * RF_observationSize) - 1);
      k = 0;
      for (j = 1; j <= RF_observationSize; j++) {
        k += RF_bootstrapIn[i][j];
      }
      if(k != RF_bootstrapSize) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Parameter verification failed.");
        RFprintf("\nRF-SRC:  Bootstrap size must be the size specified:  %12d found vs. %12d specified \n", k, RF_bootstrapSize);
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
  }
  RF_bootMembershipFlag = (char **) new_vvector(1, RF_forestSize, NRUTIL_CPTR);
  RF_bootMembershipCount = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
  RF_oobMembershipFlag = (char **) new_vvector(1, RF_forestSize, NRUTIL_CPTR);
  RF_ibgMembershipIndex = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
  RF_oobMembershipIndex = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
  maxSize = RF_observationSize;
  if (RF_bootstrapSize > RF_observationSize ) {
    maxSize = RF_bootstrapSize;
  }
  RF_identityMembershipIndex = uivector(1, maxSize);
  for (i = 1; i <= maxSize; i++) {
    RF_identityMembershipIndex[i] = i;
  }
  RF_oobSize = uivector(1, RF_forestSize);
  RF_ibgSize = uivector(1, RF_forestSize);
  RF_maxDepth = uivector(1, RF_forestSize);
  RF_serialTreeIndex = uivector(1, RF_forestSize);
  if (RF_timeIndex > 0) {
    RF_masterTime  = dvector(1, RF_observationSize);
    RF_masterTimeIndexIn  = uivector(1, RF_observationSize);
  }
  RF_root = (Node **) new_vvector(1, RF_forestSize, NRUTIL_NPTR);
  for (i = 1; i <= RF_forestSize; i++) {
    RF_root[i] = NULL;
  }
  if (RF_ptnCount > 0) {
    RF_pNodeMembership = (Node ***)     new_vvector(1, RF_forestSize, NRUTIL_NPTR2);
    RF_pTermMembership = (Terminal ***) new_vvector(1, RF_forestSize, NRUTIL_NPTR2);
    RF_pNodeList = (Node ***)     new_vvector(1, RF_forestSize, NRUTIL_NPTR2);
    RF_pTermList = (Terminal ***) new_vvector(1, RF_forestSize, NRUTIL_NPTR2);
    RF_pLeafCount = uivector(1, RF_forestSize);
  }
  RF_orderedLeafCount = uivector(1, RF_forestSize);
  for (i = 1; i <= RF_forestSize; i++) {
    RF_orderedLeafCount[i] = 0;
  }
  if ( ( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
       (!(RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
       ( (RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ) {
    for (i = 1; i <= RF_observationSize; i++) {
      if(RF_caseWeight[i] < 0) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Parameter verification failed.");
        RFprintf("\nRF-SRC:  Case-weight elements must be greater than or equal to zero:  %12.4f \n", RF_caseWeight[i]);
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
    stackWeights(RF_caseWeight,
                 RF_observationSize,
                 &RF_caseWeightType,
                 &RF_caseWeightSorted,
                 &RF_caseWeightDensitySize); 
  }
}
void unstackPreDefinedCommonArrays() {
  uint maxSize;
  free_new_vvector(RF_nodeMembership, 1, RF_forestSize, NRUTIL_NPTR2);
  free_new_vvector(RF_tTermMembership, 1, RF_forestSize, NRUTIL_NPTR2);
  free_new_vvector(RF_tNodeList, 1, RF_forestSize, NRUTIL_NPTR2);
  free_uivector(RF_tNodeListLength, 1, RF_forestSize);
  free_new_vvector(RF_tTermList, 1, RF_forestSize, NRUTIL_NPTR2);
  free_new_vvector(RF_bootMembershipIndex, 1, RF_forestSize, NRUTIL_UPTR);
  if ( (RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2) ) {
    free_new_vvector(RF_bootstrapIn, 1, RF_forestSize, NRUTIL_UPTR);
  }
  free_new_vvector(RF_bootMembershipFlag, 1, RF_forestSize, NRUTIL_CPTR);
  free_new_vvector(RF_bootMembershipCount, 1, RF_forestSize, NRUTIL_UPTR);
  free_new_vvector(RF_oobMembershipFlag, 1, RF_forestSize, NRUTIL_CPTR);
  free_new_vvector(RF_ibgMembershipIndex, 1, RF_forestSize, NRUTIL_UPTR);
  free_new_vvector(RF_oobMembershipIndex, 1, RF_forestSize, NRUTIL_UPTR);
  maxSize = RF_observationSize;
  if (RF_bootstrapSize > RF_observationSize ) {
    maxSize = RF_bootstrapSize;
  }
  free_uivector(RF_identityMembershipIndex, 1, maxSize);
  free_uivector(RF_oobSize, 1, RF_forestSize);
  free_uivector(RF_ibgSize, 1, RF_forestSize);
  free_uivector(RF_maxDepth, 1, RF_forestSize);
  free_uivector(RF_serialTreeIndex, 1, RF_forestSize);
  if (RF_timeIndex > 0) {
    free_dvector(RF_masterTime, 1, RF_observationSize);
    free_uivector(RF_masterTimeIndexIn, 1, RF_observationSize);
  }
  free_new_vvector(RF_root, 1, RF_forestSize, NRUTIL_NPTR);
  if (RF_ptnCount > 0) {
    free_new_vvector(RF_pNodeMembership, 1, RF_forestSize, NRUTIL_NPTR2);
    free_new_vvector(RF_pTermMembership, 1, RF_forestSize, NRUTIL_NPTR2);
    free_new_vvector(RF_pNodeList, 1, RF_forestSize, NRUTIL_NPTR2);
    free_new_vvector(RF_pTermList, 1, RF_forestSize, NRUTIL_NPTR2);
    free_uivector(RF_pLeafCount, 1, RF_forestSize);
  }
  free_uivector(RF_orderedLeafCount, 1, RF_forestSize);
  if ( ( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
       (!(RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
       ( (RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ) {
    unstackWeights(RF_caseWeightType, RF_observationSize, RF_caseWeightSorted); 
  }
}
void stackPreDefinedGrowthArrays() {
  uint i;
  if (RF_opt & OPT_VIMP) {
    RF_intrPredictor = uivector(1, RF_intrPredictorSize);
    for (i = 1; i <= RF_intrPredictorSize; i++) {
      RF_intrPredictor[i] = i;
    }
    RF_importanceFlag = cvector(1, RF_xSize);
    for (i = 1; i <= RF_xSize; i++) {
      RF_importanceFlag[i] = TRUE;
    }
  }
  stackWeights(RF_xWeight,
               RF_xSize,
               &RF_xWeightType,
               &RF_xWeightSorted,
               &RF_xWeightDensitySize); 
}
void unstackPreDefinedGrowthArrays() {
  if (RF_opt & OPT_VIMP) {
    free_uivector(RF_intrPredictor, 1, RF_intrPredictorSize);
    free_cvector(RF_importanceFlag, 1, RF_xSize);
  }
  unstackWeights(RF_xWeightType,
                 RF_xSize,
                 RF_xWeightSorted); 
}
void stackPreDefinedRestoreArrays() {
  uint i;
  RF_nodeCount = uivector(1, RF_forestSize);
  RF_mwcpCount = uivector(1, RF_forestSize);
  RF_restoreTreeID = uivector(1, RF_forestSize);
  RF_restoreTreeOffset = ulvector(1, RF_forestSize);
  RF_restoreMWCPOffset = ulvector(1, RF_forestSize);
  RF_mwcpPtr = (uint**) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
  for (i = 1; i <= RF_forestSize; i++) {
    RF_nodeCount[i] = RF_mwcpCount[i] = RF_restoreTreeID[i] = 0;
    RF_restoreTreeOffset[i] = RF_restoreMWCPOffset[i] = 0;
  }
  if (RF_opt & OPT_VIMP) {
    checkInteraction();
    RF_importanceFlag = cvector(1, RF_xSize);
    for (i = 1; i <= RF_xSize; i++) {
      RF_importanceFlag[i] = FALSE;
    }
    for (i = 1; i <= RF_intrPredictorSize; i++) {
      RF_importanceFlag[RF_intrPredictor[i]] = TRUE;
    }
  }
  if(RF_sobservationSize > 0) {
    hpsortui(RF_sobservationIndv, RF_sobservationSize);
    uint j = 1;
    for (uint i = 2; i <= RF_sobservationSize; i++) {
      if (RF_sobservationIndv[i] > RF_sobservationIndv[j]) {
        j ++;
      }
    }
    if (RF_sobservationSize != j) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Parameter verification failed.");
      RFprintf("\nRF-SRC:  Subsetted individuals are not unique:  %10d of %10d are unique.", j, RF_sobservationSize);
      error("\nRF-SRC:  The application will now exit.\n");
    }
    for (uint i = 1; i <= RF_sobservationSize; i++) {
      if (RF_sobservationIndv[i] > RF_observationSize) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Parameter verification failed.");
        RFprintf("\nRF-SRC:  Subsetted individuals are not coherent.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
    RF_soobSize = uivector(1, RF_forestSize);
  }
}
void unstackPreDefinedRestoreArrays() {
  free_uivector(RF_nodeCount, 1, RF_forestSize);
  free_uivector(RF_mwcpCount, 1, RF_forestSize);
  free_uivector(RF_restoreTreeID, 1, RF_forestSize);
  free_ulvector(RF_restoreTreeOffset, 1, RF_forestSize);
  free_ulvector(RF_restoreMWCPOffset, 1, RF_forestSize);
  free_new_vvector(RF_mwcpPtr, 1, RF_forestSize, NRUTIL_UPTR);
  if (RF_opt & OPT_VIMP) {
    free_cvector(RF_importanceFlag, 1, RF_xSize);
  }
  if(RF_sobservationSize > 0) {
    free_uivector(RF_soobSize, 1, RF_forestSize);
  }
}
void stackPreDefinedPredictArrays() {
  uint i;
  RF_fnodeMembership = (Node ***)     new_vvector(1, RF_forestSize, NRUTIL_NPTR2);
  RF_ftTermMembership = (Terminal ***) new_vvector(1, RF_forestSize, NRUTIL_NPTR2);
  RF_fidentityMembershipIndex = uivector(1, RF_fobservationSize);
  for (i = 1; i <= RF_fobservationSize; i++) {
    RF_fidentityMembershipIndex[i] = i;
  }
  RF_testMembershipFlag = cvector(1, RF_fobservationSize);
  for (i = 1; i <= RF_fobservationSize; i++) {
    RF_testMembershipFlag[i] = ACTIVE;
  }
  RF_nodeCount = uivector(1, RF_forestSize);
  RF_mwcpCount = uivector(1, RF_forestSize);
  RF_restoreTreeID = uivector(1, RF_forestSize);
  RF_restoreTreeOffset = ulvector(1, RF_forestSize);
  RF_restoreMWCPOffset = ulvector(1, RF_forestSize);
  RF_mwcpPtr = (uint**) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
  for (i = 1; i <= RF_forestSize; i++) {
    RF_nodeCount[i] = RF_mwcpCount[i] = RF_restoreTreeID[i] = 0;
    RF_restoreTreeOffset[i] = RF_restoreMWCPOffset[i] = 0;
  }
  if (RF_opt & OPT_VIMP) {
    checkInteraction();
    RF_importanceFlag = cvector(1, RF_xSize);
    for (i = 1; i <= RF_xSize; i++) {
      RF_importanceFlag[i] = FALSE;
    }
    for (i = 1; i <= RF_intrPredictorSize; i++) {
      RF_importanceFlag[RF_intrPredictor[i]] = TRUE;
    }
  }
}
void unstackPreDefinedPredictArrays() {
  free_new_vvector(RF_fnodeMembership, 1, RF_forestSize, NRUTIL_NPTR2);
  free_new_vvector(RF_ftTermMembership, 1, RF_forestSize, NRUTIL_NPTR2);
  free_uivector(RF_fidentityMembershipIndex, 1, RF_fobservationSize);
  free_cvector(RF_testMembershipFlag, 1, RF_fobservationSize);
  free_uivector(RF_nodeCount, 1, RF_forestSize);
  free_uivector(RF_mwcpCount, 1, RF_forestSize);
  free_uivector(RF_restoreTreeID, 1, RF_forestSize);
  free_ulvector(RF_restoreTreeOffset, 1, RF_forestSize);
  free_ulvector(RF_restoreMWCPOffset, 1, RF_forestSize);
  free_new_vvector(RF_mwcpPtr, 1, RF_forestSize, NRUTIL_UPTR);
  if (RF_opt & OPT_VIMP) {
    free_cvector(RF_importanceFlag, 1, RF_xSize);
  }
}
void checkInteraction() {
  uint leadingIndex, i;
  if((RF_intrPredictorSize <= 0) || (RF_intrPredictorSize > RF_xSize)) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Parameter verification failed.");
    Rprintf("\nRF-SRC:  Number of predictors to be perturbed must be greater than zero and less than or equal to %10d:  %10d \n", RF_xSize, RF_intrPredictorSize);
    error("\nRF-SRC:  The application will now exit.\n");
  }
  uint *intrPredictorCopy = uivector(1, RF_intrPredictorSize);
  for (i=1; i <= RF_intrPredictorSize; i++) {
    intrPredictorCopy[i] = RF_intrPredictor[i];
  }
  hpsortui(intrPredictorCopy, RF_intrPredictorSize);
  leadingIndex = 1;
  for (i=2; i <= RF_intrPredictorSize; i++) {
    if (intrPredictorCopy[i] > intrPredictorCopy[leadingIndex]) {
      leadingIndex++;
    }
  }
  if (RF_intrPredictorSize != leadingIndex) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Parameter verification failed.");
    Rprintf("\nRF-SRC:  Interaction terms are not unique.");
    Rprintf("\nRF-SRC:  Only %10d of %10d are unique.", leadingIndex, RF_intrPredictorSize);
    error("\nRF-SRC:  The application will now exit.\n");
  }
  free_uivector(intrPredictorCopy, 1, RF_intrPredictorSize);
  for (i=1; i <= RF_intrPredictorSize; i++) {
    if (RF_intrPredictor[i] > RF_xSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  Parameter verification failed.");
      Rprintf("\nRF-SRC:  Interaction terms are not coherent.");
      Rprintf("\nRF-SRC:  Predictor encountered is %10d, maximum allowable is %10d.", RF_intrPredictor[i], RF_xSize);
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
}
void stackWeights(double *weight,
                  uint    size,
                  uint   *weightType,
                  uint  **weightSorted,
                  uint   *weightDensitySize) {
  char uniformFlag, integerFlag;
  double meanWeight;
  uint i;
  *weightSorted      = NULL;
  *weightDensitySize = 0;
  meanWeight = getMeanValue(weight, size);
  uniformFlag = TRUE;
  i = 0;
  while (uniformFlag && (i < size)) {
    ++i;
    if (fabs(weight[i] - meanWeight) > EPSILON) {
      uniformFlag = FALSE;
    }
  }
  if (uniformFlag) {
    *weightType = RF_WGHT_UNIFORM;
  } 
  else {
    integerFlag = TRUE;
    i = 0;
    while (integerFlag && (i < size)) {
      i++;
      if (fabs(round(weight[i]) - weight[i]) > EPSILON) {
        integerFlag = FALSE;
      }
    }
    if(integerFlag) {
      *weightType = RF_WGHT_INTEGER;
    }
    else {
      *weightType = RF_WGHT_GENERIC;
    }
  }
  switch (*weightType) {
  case RF_WGHT_UNIFORM:
    break;
  case RF_WGHT_INTEGER:
    *weightSorted = uivector(1, size);
    indexx(size, weight, *weightSorted);
    *weightDensitySize = 0;
    for (i = 1; i <= size; i++) {
      *weightDensitySize += (uint) weight[i];
    }
    break;
  case RF_WGHT_GENERIC:
    *weightSorted = uivector(1, size);
    indexx(size, weight, *weightSorted);
    break;
  }
}
void unstackWeights(uint    weightType,
                    uint    size,
                    uint   *weightSorted) {
  switch (weightType) {
  case RF_WGHT_UNIFORM:
    break;
  case RF_WGHT_INTEGER:
    free_uivector(weightSorted, 1, size);
    break;
  case RF_WGHT_GENERIC:
    free_uivector(weightSorted, 1, size);
    break;
  }
}
uint stackDefinedOutputObjects(char      mode,
                               char    **sexpString,
                               Node   ***pRF_root,
                               uint    **pRF_tLeafCount,
                               double  **pRF_proximity,
                               int     **pRF_seed,
                               double  **p_imputation,
                               double ***pRF_sImputeResponsePtr,
                               double ***pRF_sImputePredictorPtr,
                               uint    **pRF_varUsed,
                               uint   ***pRF_varUsedPtr,
                               double  **p_splitDepth,
                               uint     *sexpIndex,
                               uint     *stackCount,
                               SEXP     *sexpVector) {
  uint sexpIdentity;
  ulong localSize;
  uint xVimpSize;
  uint  obsSize;
  uint  mRecordSize;
  uint *mRecordIndex;
  double **responsePtr;
  double **predictorPtr;
  uint     rspSize;
  uint     dpthDimOne;
  uint     vuseDimOne;
  uint   **ensembleDen;
  double **ensembleSRG;
  double **ensembleMRT;
  double **ensembleCIF;
  double **ensembleSRV;
  double **ensembleCLS;
  double **ensembleRGR;
  double ****ensembleSRGptr;
  double  ***ensembleMRTptr;
  double  ***ensembleSRVptr;
  double ****ensembleCIFptr;
  double ****ensembleCLSptr;
  double  ***ensembleRGRptr;
  double ****ensembleSRGnum;
  double  ***ensembleMRTnum;
  double  ***ensembleSRVnum;
  double ****ensembleCIFnum;
  double ****ensembleCLSnum;
  double  ***ensembleRGRnum;
  char oobFlag, fullFlag;
  char maxVoteFlag;
  uint dimThree;
  uint i, j, k, m;
  xVimpSize      = 0;  
  dpthDimOne     = 0;  
  vuseDimOne     = 0;  
  obsSize        = 0;  
  mRecordSize    = 0;  
  rspSize        = 0;  
  sexpIdentity   = 0;  
  responsePtr    = NULL;  
  predictorPtr   = NULL;  
  mRecordIndex   = NULL;  
  RF_rTargetFactor = RF_rTargetNonFactor = NULL;
  RF_rTargetFactorCount = RF_rTargetNonFactorCount = 0;
  if (RF_opt & (OPT_VARUSED_F | OPT_VARUSED_T)) {
    if (RF_opt & OPT_VARUSED_F) {
      vuseDimOne = 1;
    }
    else {
      vuseDimOne = RF_forestSize;
    }
  }
  if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
    if (RF_opt & OPT_SPLDPTH_F) {
      dpthDimOne = 1;
    }
    else {
      dpthDimOne = RF_forestSize;
    }
  }
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    mRecordSize = RF_fmRecordSize;
    rspSize = RF_frSize;
    responsePtr  = RF_fresponseIn;
    predictorPtr = RF_fobservationIn;
    mRecordIndex = RF_fmRecordIndex;
    if (RF_rSize == 0) {
      RF_rTarget = NULL;
      RF_rTargetCount = 0;
    }
    else {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        RF_rTarget = NULL;
        RF_rTargetCount = 0;
      }
      else {
        RF_rTargetFactor    = uivector(1, RF_rTargetCount);
        RF_rTargetNonFactor = uivector(1, RF_rTargetCount);
        RF_rTargetFactorCount = RF_rTargetNonFactorCount = 0;
        for (i = 1; i <= RF_rTargetCount; i++) {
          if ((RF_rTarget[i] < 1) || (RF_rTarget[i] > RF_rSize)) {
            RFprintf("\nRF-SRC:  *** ERROR *** ");
            RFprintf("\nRF-SRC:  Target response is out of range for [C+], [R+], [M+]:  %10d %10d ", i, RF_rTarget[i]);
            error("\nRF-SRC:  The application will now exit.\n");
          }
          if ((strcmp(RF_rType[RF_rTarget[i]], "C") == 0) || 
              (strcmp(RF_rType[RF_rTarget[i]], "I") == 0)) {
            RF_rTargetFactor[++RF_rTargetFactorCount] = RF_rTarget[i];
          }
          else {
            RF_rTargetNonFactor[++RF_rTargetNonFactorCount] = RF_rTarget[i];
          }
        }
      }
    }  
    *stackCount = 1;
    if (RF_opt & OPT_FENS) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        (*stackCount) += 3;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          (*stackCount) += 1;
        }
        if (RF_rTargetNonFactorCount > 0) {
          (*stackCount) += 1;
        }
      }
    }
    if (RF_opt & OPT_PERF) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        (*stackCount) += 1;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          (*stackCount) += 1;
        }
        if (RF_rTargetNonFactorCount > 0) {
          (*stackCount) += 1;
        }
      }
    }
    if (RF_opt & OPT_PROX) {
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_NODE_STAT) {
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_MISS) {
      (*stackCount) += 1;
    }
    if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_VIMP) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        (*stackCount) += 1;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          (*stackCount) += 1;
        }
        if (RF_rTargetNonFactorCount > 0) {
          (*stackCount) += 1;
        }
      }
    }
    if (RF_optHigh & OPT_MEMB_PRUN) {
      (*stackCount) += 1;
    }
    if (RF_optHigh & OPT_MEMB_USER) {
      (*stackCount) += 2;
    }
    break;
  default:
    obsSize = RF_observationSize;
    mRecordSize = RF_mRecordSize;
    rspSize = RF_rSize;
    responsePtr  = RF_responseIn;
    predictorPtr = RF_observationIn;
    mRecordIndex = RF_mRecordIndex;
    if (RF_rSize == 0) {
      RF_rTarget = NULL;
      RF_rTargetCount = 0;
    }
    else {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        RF_rTarget = NULL;
        RF_rTargetCount = 0;
      }
      else {
        if (mode == RF_GROW) {
          RF_rTargetCount = rspSize;
          RF_rTarget = uivector(1 , RF_rTargetCount);
          for (i = 1; i <= RF_rTargetCount; i++) {
            RF_rTarget[i] = i;
          }
        }
        RF_rTargetFactor    = uivector(1, RF_rTargetCount);
        RF_rTargetNonFactor = uivector(1, RF_rTargetCount);
        RF_rTargetFactorCount = RF_rTargetNonFactorCount = 0;
        for (i = 1; i <= RF_rTargetCount; i++) {
          if ((RF_rTarget[i] < 1) || (RF_rTarget[i] > RF_rSize)) {
            RFprintf("\nRF-SRC:  *** ERROR *** ");
            RFprintf("\nRF-SRC:  Target response is out of range for [C+], [R+], [M+]:  %10d %10d ", i, RF_rTarget[i]);
            error("\nRF-SRC:  The application will now exit.\n");
          }
          if ((strcmp(RF_rType[RF_rTarget[i]], "C") == 0) ||
              (strcmp(RF_rType[RF_rTarget[i]], "I") == 0)) {
            RF_rTargetFactor[++RF_rTargetFactorCount] = RF_rTarget[i];
          }
          else {
            RF_rTargetNonFactor[++RF_rTargetNonFactorCount] = RF_rTarget[i];
          }
        }
      }  
    }  
    *stackCount = 0;
    if (RF_opt & OPT_LEAF) {
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_FENS) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        (*stackCount) += 3;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          (*stackCount) += 1;
        }
        if (RF_rTargetNonFactorCount > 0) {
          (*stackCount) += 1;
        }
      }
    }
    if (RF_opt & OPT_OENS) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        (*stackCount) += 3;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          (*stackCount) += 1;
        }
        if (RF_rTargetNonFactorCount > 0) {
          (*stackCount) += 1;
        }
      }
    }
    if (RF_opt & OPT_PERF) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        (*stackCount) += 1;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          (*stackCount) += 1;
        }
        if (RF_rTargetNonFactorCount > 0) {
          (*stackCount) += 1;
        }
      }
    }
    if (RF_opt & OPT_PROX) {
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_SEED) {
      if (RF_opt & OPT_TREE) {
        (*stackCount) += 1;
        (*stackCount) += 7;
      }
      else {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  SEXP TREE output request inconsistent.");
        RFprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
    if (RF_opt & OPT_NODE_STAT) {
      (*stackCount) += 1;
      if (mode == RF_GROW) {
        (*stackCount) += 2;
      }
    }
    if (RF_opt & OPT_USPV_STAT) {
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_MISS) {
      (*stackCount) += 1;
    }
    if (RF_opt & (OPT_VARUSED_F | OPT_VARUSED_T)) {
      (*stackCount) += 1;
    }
    if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_VIMP) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        (*stackCount) += 1;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          (*stackCount) += 1;
        }
        if (RF_rTargetNonFactorCount > 0) {
          (*stackCount) += 1;
        }
      }
    }
    if (RF_optHigh & OPT_MEMB_PRUN) {
      (*stackCount) += 1;
    }
    if (RF_optHigh & OPT_MEMB_USER) {
      (*stackCount) += 2;
    }
    if (RF_optHigh & OPT_MEMB_OUTG) {
      (*stackCount) += 4;
    }
    if (RF_optHigh & OPT_TERM_OUTG) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        (*stackCount) += 1;
        if (!(RF_opt & OPT_COMP_RISK)) {
          (*stackCount) += 1;
          (*stackCount) += 1;
        }
        else {
          (*stackCount) += 2;
        }
      }
      else {
        if (RF_rTargetNonFactorCount > 0) {
          (*stackCount) += 1;
        }
        if (RF_rTargetFactorCount > 0) {
          (*stackCount) += 1;
        }
      }
    }
    if (RF_optHigh & OPT_PART_PLOT) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        (*stackCount) += 1;
      }
      else {
        if (RF_rTargetFactorCount > 0) {
          (*stackCount) += 1;
        }
        if (RF_rTargetNonFactorCount > 0) {
          (*stackCount) += 1;
        }
      }
    }
    break;
  }  
  initProtect(sexpVector, *stackCount);
  oobFlag = fullFlag = FALSE;
  if ((RF_opt & OPT_FENS) | (RF_opt & OPT_OENS)) {
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    if (RF_opt & OPT_OENS) {
      oobFlag = TRUE;
    }
    while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
      ensembleDen    = NULL;
      ensembleSRG    = NULL;
      ensembleSRGptr = NULL;
      ensembleSRGnum = NULL;
      ensembleMRT    = NULL;
      ensembleMRTptr = NULL;
      ensembleMRTnum = NULL;
      ensembleSRV    = NULL;
      ensembleSRVptr = NULL;
      ensembleSRVnum = NULL;
      ensembleCIF    = NULL;
      ensembleCIFptr = NULL;
      ensembleCIFnum = NULL;
      ensembleCLS    = NULL;
      ensembleCLSptr = NULL;
      ensembleCLSnum = NULL;
      ensembleRGR    = NULL;
      ensembleRGRptr = NULL;
      ensembleRGRnum = NULL;
      if (oobFlag == TRUE) {
        ensembleDen    = &RF_oobEnsembleDen;
        ensembleSRG    = &RF_oobEnsembleSRG_;
        ensembleSRGptr = &RF_oobEnsembleSRGptr;
        ensembleSRGnum = &RF_oobEnsembleSRGnum;
        ensembleMRT    = &RF_oobEnsembleMRT_;
        ensembleMRTptr = &RF_oobEnsembleMRTptr;
        ensembleMRTnum = &RF_oobEnsembleMRTnum;
        ensembleSRV    = &RF_oobEnsembleSRV_;
        ensembleSRVptr = &RF_oobEnsembleSRVptr;
        ensembleSRVnum = &RF_oobEnsembleSRVnum;
        ensembleCIF    = &RF_oobEnsembleCIF_;
        ensembleCIFptr = &RF_oobEnsembleCIFptr;
        ensembleCIFnum = &RF_oobEnsembleCIFnum;
        ensembleCLS    = &RF_oobEnsembleCLS_;
        ensembleCLSptr = &RF_oobEnsembleCLSptr;
        ensembleCLSnum = &RF_oobEnsembleCLSnum;
        ensembleRGR    = &RF_oobEnsembleRGR_;
        ensembleRGRptr = &RF_oobEnsembleRGRptr;
        ensembleRGRnum = &RF_oobEnsembleRGRnum;
      }
      else {
        ensembleDen    = &RF_fullEnsembleDen;
        ensembleSRG    = &RF_fullEnsembleSRG_;
        ensembleSRGptr = &RF_fullEnsembleSRGptr;
        ensembleSRGnum = &RF_fullEnsembleSRGnum;
        ensembleMRT    = &RF_fullEnsembleMRT_;
        ensembleMRTptr = &RF_fullEnsembleMRTptr;
        ensembleMRTnum = &RF_fullEnsembleMRTnum;        
        ensembleSRV    = &RF_fullEnsembleSRV_;
        ensembleSRVptr = &RF_fullEnsembleSRVptr;
        ensembleSRVnum = &RF_fullEnsembleSRVnum;
        ensembleCIF    = &RF_fullEnsembleCIF_;
        ensembleCIFptr = &RF_fullEnsembleCIFptr;
        ensembleCIFnum = &RF_fullEnsembleCIFnum;
        ensembleCLS    = &RF_fullEnsembleCLS_;
        ensembleCLSptr = &RF_fullEnsembleCLSptr;
        ensembleCLSnum = &RF_fullEnsembleCLSnum;
        ensembleRGR    = &RF_fullEnsembleRGR_;
        ensembleRGRptr = &RF_fullEnsembleRGRptr;
        ensembleRGRnum = &RF_fullEnsembleRGRnum;
      }
      *ensembleDen = uivector(1, obsSize);
      for (i = 1; i <= obsSize; i++) {
        (*ensembleDen)[i] = 0;
      }
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        {
          (oobFlag == TRUE) ? (sexpIdentity = RF_OSRG_ID) : ((fullFlag == TRUE) ? sexpIdentity = RF_ASRG_ID : TRUE);
          localSize = (ulong) RF_eventTypeSize * RF_sortedTimeInterestSize * obsSize;
          *ensembleSRG = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, sexpIdentity, localSize, sexpVector, sexpString);
          *ensembleSRGptr = (double ***) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR2);
          *ensembleSRGnum = (double ***) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR2);
          for (j = 1; j <= RF_eventTypeSize; j++) {
            (*ensembleSRGptr)[j] = (double **) new_vvector(1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
            (*ensembleSRGnum)[j] = (double **) new_vvector(1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
            for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
              (*ensembleSRGptr)[j][k]  = (*ensembleSRG) + ((j-1) * RF_sortedTimeInterestSize * obsSize) + ((k-1) * obsSize) - 1;
              (*ensembleSRGnum)[j][k]  = dvector(1, obsSize);
              for (i = 1; i <= obsSize; i++) {
                (*ensembleSRGptr)[j][k][i] = 0.0;
                (*ensembleSRGnum)[j][k][i] = 0.0;
              }
            }
          }
          (oobFlag == TRUE) ? (sexpIdentity = RF_OMRT_ID) : ((fullFlag == TRUE) ? sexpIdentity = RF_AMRT_ID: TRUE);
          localSize = (ulong) RF_eventTypeSize * obsSize;
          *ensembleMRT = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, sexpIdentity, localSize, sexpVector, sexpString);
          *ensembleMRTptr = (double **) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR);
          *ensembleMRTnum = (double **) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR);
          for (j = 1; j <= RF_eventTypeSize; j++) {
            (*ensembleMRTptr)[j] = (*ensembleMRT) + ((j-1) * obsSize) - 1;
            (*ensembleMRTnum)[j] = dvector(1, obsSize);
            for (i = 1; i <= obsSize; i++) {
              (*ensembleMRTptr)[j][i] = 0.0;
              (*ensembleMRTnum)[j][i] = 0.0;
            }
          }
          if (!(RF_opt & OPT_COMP_RISK)) {
            (oobFlag == TRUE) ? (sexpIdentity = RF_OSRV_ID) : ((fullFlag == TRUE) ? sexpIdentity = RF_ASRV_ID: TRUE);
            localSize = (ulong) RF_sortedTimeInterestSize * obsSize;
            *ensembleSRV = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, sexpIdentity, localSize, sexpVector, sexpString);
            *ensembleSRVptr = (double **) new_vvector(1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
            *ensembleSRVnum = (double **) new_vvector(1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
            for (j = 1; j <= RF_sortedTimeInterestSize; j++) {
              (*ensembleSRVptr)[j]  = (*ensembleSRV) + ((j-1) * obsSize) - 1;
              (*ensembleSRVnum)[j]  = dvector(1, obsSize);
              for (i = 1; i <= obsSize; i++) {
                (*ensembleSRVptr)[j][i]  = 0.0;
                (*ensembleSRVnum)[j][i]  = 0.0;
              }
            }
          }  
          else {
            (oobFlag == TRUE) ? (sexpIdentity = RF_OCIF_ID) : ((fullFlag == TRUE) ? sexpIdentity = RF_ACIF_ID: TRUE);
            localSize = (ulong) RF_eventTypeSize * RF_sortedTimeInterestSize * obsSize;
            *ensembleCIF = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, sexpIdentity, localSize, sexpVector, sexpString);
            *ensembleCIFptr = (double ***) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR2);
            *ensembleCIFnum = (double ***) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR2);
            for (j = 1; j <= RF_eventTypeSize; j++) {
              (*ensembleCIFptr)[j] = (double **) new_vvector(1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
              (*ensembleCIFnum)[j] = (double **) new_vvector(1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
              for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
                (*ensembleCIFptr)[j][k]  = (*ensembleCIF) + ((j-1) * RF_sortedTimeInterestSize * obsSize) + ((k-1) * obsSize) - 1;
                (*ensembleCIFnum)[j][k]  = dvector(1, obsSize);
                for (i = 1; i <= obsSize; i++) {
                  (*ensembleCIFptr)[j][k][i] = 0.0;
                  (*ensembleCIFnum)[j][k][i] = 0.0;
                }
              }
            }
          }  
        }
      }  
      else {
        if (RF_rTargetFactorCount > 0) {
          (oobFlag == TRUE) ? (sexpIdentity = RF_OCLS_ID) : ((fullFlag == TRUE) ? sexpIdentity = RF_ACLS_ID: TRUE);
          localSize = 0;
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              localSize += (ulong) obsSize;
            }
          }
          *ensembleCLS = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, sexpIdentity, localSize, sexpVector, sexpString);
          *ensembleCLSptr = (double ***) new_vvector(1, RF_rTargetFactorCount, NRUTIL_DPTR2);
          *ensembleCLSnum = (double ***) new_vvector(1, RF_rTargetFactorCount, NRUTIL_DPTR2);
          localSize = 0;
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            (*ensembleCLSptr)[j] = (double **) new_vvector(1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]], NRUTIL_DPTR);
            (*ensembleCLSnum)[j] = (double **) new_vvector(1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]], NRUTIL_DPTR);
            for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              (*ensembleCLSptr)[j][k]  = (*ensembleCLS) + localSize - 1;
              (*ensembleCLSnum)[j][k]  = dvector(1, obsSize);
              localSize += (ulong) obsSize;
              for (i = 1; i <= obsSize; i++) {
                (*ensembleCLSptr)[j][k][i] = 0.0;
                (*ensembleCLSnum)[j][k][i] = 0.0;
              }
            }
          }
        }
        if (RF_rTargetNonFactorCount > 0) {
          (oobFlag == TRUE) ? (sexpIdentity = RF_ORGR_ID) : ((fullFlag == TRUE) ? sexpIdentity = RF_ARGR_ID: TRUE);
          localSize = (ulong) RF_rTargetNonFactorCount * obsSize;
          *ensembleRGR = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, sexpIdentity, localSize, sexpVector, sexpString);
          (*ensembleRGRptr) = (double **) new_vvector(1, RF_rTargetNonFactorCount, NRUTIL_DPTR);
          (*ensembleRGRnum) = (double **) new_vvector(1, RF_rTargetNonFactorCount, NRUTIL_DPTR);
          for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
            (*ensembleRGRptr)[j] = (*ensembleRGR) + ((j-1) * obsSize) - 1;
            (*ensembleRGRnum)[j] = dvector(1, obsSize);
            for (i = 1; i <= obsSize; i++) {
              (*ensembleRGRptr)[j][i] = NA_REAL;
              (*ensembleRGRnum)[j][i] = 0.0;
            }
          }
        }
      }
      if (oobFlag == TRUE) {
        oobFlag = FALSE;
      }
      else {
        fullFlag = FALSE;
      }
    }  
  }
  if (RF_opt & OPT_PERF) {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      {
        localSize = (ulong) RF_forestSize * RF_eventTypeSize; 
        RF_perfMRT_ = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, RF_ER_SURV, localSize, sexpVector, sexpString);
        RF_perfMRTptr = (double **) new_vvector(1, RF_forestSize, NRUTIL_DPTR);
        for (i = 1; i <= RF_forestSize; i++) {
          RF_perfMRTptr[i]  = (RF_perfMRT_)  + ((i-1) * RF_eventTypeSize) - 1;
          for (k = 1; k <= RF_eventTypeSize; k++) {
            RF_perfMRTptr[i][k] = NA_REAL;
          }
        }
      }
    }  
    else {
      if (RF_rTargetFactorCount > 0) {
        localSize = 0;
        for (j = 1; j <= RF_rTargetFactorCount; j++) {
          for (k = 1; k <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
            localSize += (ulong) RF_forestSize;
          }
        }
        RF_perfCLS_ = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, RF_ER_CLAS, localSize, sexpVector, sexpString);        
        RF_perfCLSptr = (double ***) new_vvector(1, RF_forestSize, NRUTIL_DPTR2);
        localSize = 0;
        for (i = 1; i <= RF_forestSize; i++) {
          RF_perfCLSptr[i] = (double **) new_vvector(1, RF_rTargetFactorCount, NRUTIL_DPTR);
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            RF_perfCLSptr[i][j]  = (RF_perfCLS_) + localSize - 1;
            localSize += (ulong) 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]];
            for (k = 1; k <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              RF_perfCLSptr[i][j][k]  = NA_REAL;
            }
          }
        }
      }
      if (RF_rTargetNonFactorCount > 0) {
        localSize = (ulong) RF_forestSize * RF_rTargetNonFactorCount;
        RF_perfRGR_ = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, RF_ER_REGR, localSize, sexpVector, sexpString);
        RF_perfRGRptr = (double **) new_vvector(1, RF_forestSize, NRUTIL_DPTR);
        for (i = 1; i <= RF_forestSize; i++) {
          RF_perfRGRptr[i] = (RF_perfRGR_) + ((i-1) * RF_rTargetNonFactorCount) - 1;
          for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
            RF_perfRGRptr[i][j] = NA_REAL;
          }
        }
      } 
    }
  }  
  if (RF_opt & OPT_VIMP) {
    RF_vimpEnsembleDen = NULL;
    RF_vimpEnsembleMRT = NULL;
    RF_vimpEnsembleCLS = NULL;
    RF_vimpEnsembleRGR = NULL;
    RF_vimpMRTleo = NULL;
    RF_vimpCLSleo = NULL;
    RF_vimpRGRleo = NULL;
    RF_vimpMRTptr = NULL;
    RF_vimpRGRptr = NULL;
    RF_vimpCLSptr = NULL;
    if (RF_opt & OPT_VIMP_JOIN) {
      xVimpSize = 1;
    }
    else {
      xVimpSize = RF_intrPredictorSize;
    }
    RF_vimpMembership = (Terminal ****) new_vvector(1, xVimpSize, NRUTIL_NPTR3);
    for (k = 1; k <= xVimpSize; k++) {
      RF_vimpMembership[k] = (Terminal ***) new_vvector(1,  RF_forestSize, NRUTIL_NPTR2);
    }
    for (k = 1; k <= xVimpSize; k++) {
      for (i = 1; i <= RF_forestSize; i++) {
        RF_vimpMembership[k][i] = NULL;
      }
    }
    RF_vimpEnsembleDen  = (uint **) new_vvector(1, xVimpSize, NRUTIL_UPTR);
    for (j = 1; j <= xVimpSize; j++) {
      RF_vimpEnsembleDen[j] = uivector(1, obsSize);
      for (i = 1; i <= obsSize; i++) {
        RF_vimpEnsembleDen[j][i] = 0;
      }
    }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      {
        localSize = (ulong) xVimpSize * RF_eventTypeSize;
        RF_vimpMRT_ = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, RF_VMP_SRG, localSize, sexpVector, sexpString);        
        RF_vimpMRTptr = (double **) new_vvector(1, xVimpSize, NRUTIL_DPTR);
        for (j = 1; j <= xVimpSize; j++) {
          RF_vimpMRTptr[j]  = RF_vimpMRT_ + ((j-1) * RF_eventTypeSize) - 1;
          for (k = 1; k <= RF_eventTypeSize; k++) {
            RF_vimpMRTptr[j][k] = NA_REAL;
          }
        }
        RF_vimpEnsembleMRT = (double ***) new_vvector(1, xVimpSize, NRUTIL_DPTR2);
        for (j = 1; j <= xVimpSize; j++) {
          RF_vimpEnsembleMRT[j]  = (double **) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR);
          for (k = 1; k <= RF_eventTypeSize; k++) {
            RF_vimpEnsembleMRT[j][k] = dvector(1, obsSize);
            for (m = 1; m <= obsSize; m++) {
              RF_vimpEnsembleMRT[j][k][m] = 0;
            }
          }
        }
        if(RF_opt & OPT_VIMP_LEOB) {
          RF_perfMRTleo = (double **) new_vvector(1, RF_forestSize, NRUTIL_DPTR);
          for (i = 1; i <= RF_forestSize; i++) {
            RF_perfMRTleo[i]  = dvector(1, RF_eventTypeSize);
            for (k = 1; k <= RF_eventTypeSize; k++) {
              RF_perfMRTleo[i][k] = NA_REAL;
            }
          }
          RF_vimpMRTleo = (double ***) new_vvector(1, RF_forestSize, NRUTIL_DPTR2);
          for (i = 1; i <= RF_forestSize; i++) {
            RF_vimpMRTleo[i]  = (double **) new_vvector(1, xVimpSize, NRUTIL_DPTR);
            for (j = 1; j <= xVimpSize; j++) {
              RF_vimpMRTleo[i][j] = dvector(1, RF_eventTypeSize);
              for (k = 1; k <= RF_eventTypeSize; k++) {
                RF_vimpMRTleo[i][j][k] = NA_REAL;
              }
            }
          }
        }
      }
    }  
    else {
      if (RF_rTargetFactorCount > 0) {
        localSize = 0;
        for (j = 1; j <= RF_rTargetFactorCount; j++) {
          localSize += (ulong) 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]];
        }
        localSize = (ulong) xVimpSize * localSize;
        RF_vimpCLS_ = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, RF_VMP_CLS, localSize, sexpVector, sexpString);
        RF_vimpCLSptr = (double ***) new_vvector(1, xVimpSize, NRUTIL_DPTR2);
        localSize = 0;
        for (i = 1; i <= xVimpSize; i++) {
          RF_vimpCLSptr[i] = (double **) new_vvector(1, RF_rTargetFactorCount, NRUTIL_DPTR);
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            RF_vimpCLSptr[i][j]  = (RF_vimpCLS_) + localSize - 1;
            localSize += (ulong) 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]];
            for (k = 1; k <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              RF_vimpCLSptr[i][j][k]  = NA_REAL;
            }
          }
        }
        maxVoteFlag = FALSE;
        if(RF_opt & OPT_VIMP_LEOB) {
          if (!(RF_opt & OPT_PERF_CALB)) {
            maxVoteFlag = TRUE;
          }
        }
        if (maxVoteFlag) {
          RF_vimpEnsembleCLS = (double ****) new_vvector(1, xVimpSize, NRUTIL_DPTR3);        
          for (i = 1; i <= xVimpSize; i++) {
            RF_vimpEnsembleCLS[i] = (double ***) new_vvector(1, RF_rTargetFactorCount, NRUTIL_DPTR2);
            for (j = 1; j <= RF_rTargetFactorCount; j++) {
              RF_vimpEnsembleCLS[i][j]  = (double **) new_vvector(1, 1, NRUTIL_DPTR);
              RF_vimpEnsembleCLS[i][j][1] = dvector(1, obsSize);
              for (m = 1; m <= obsSize; m++) {
                RF_vimpEnsembleCLS[i][j][1][m] = 0;
              }
            }
          }
        }
        else {
          RF_vimpEnsembleCLS = (double ****) new_vvector(1, xVimpSize, NRUTIL_DPTR3);        
          for (i = 1; i <= xVimpSize; i++) {
            RF_vimpEnsembleCLS[i] = (double ***) new_vvector(1, RF_rTargetFactorCount, NRUTIL_DPTR2);
            for (j = 1; j <= RF_rTargetFactorCount; j++) {
              RF_vimpEnsembleCLS[i][j]  = (double **) new_vvector(1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]], NRUTIL_DPTR);
              for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
                RF_vimpEnsembleCLS[i][j][k] = dvector(1, obsSize);
                for (m = 1; m <= obsSize; m++) {
                  RF_vimpEnsembleCLS[i][j][k][m] = 0;
                }
              }
            }
          }
        }
        if(RF_opt & OPT_VIMP_LEOB) {
          RF_perfCLSleo = (double ***) new_vvector(1, RF_forestSize, NRUTIL_DPTR2);
          for (i = 1; i <= RF_forestSize; i++) {
            RF_perfCLSleo[i] = (double **) new_vvector(1, RF_rTargetFactorCount, NRUTIL_DPTR);
            for (j = 1; j <= RF_rTargetFactorCount; j++) {
              RF_perfCLSleo[i][j]  = dvector(1, 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]);
              for (k = 1; k <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
                RF_perfCLSleo[i][j][k]  = NA_REAL;
              }
            }
          }
          RF_vimpCLSleo = (double ****) new_vvector(1, RF_forestSize, NRUTIL_DPTR3);
          for (i = 1; i <= RF_forestSize; i++) {
            RF_vimpCLSleo[i] = (double ***) new_vvector(1, xVimpSize, NRUTIL_DPTR2);
            for (j = 1; j <= xVimpSize; j++) {            
              RF_vimpCLSleo[i][j] = (double **) new_vvector(1, RF_rTargetFactorCount, NRUTIL_DPTR);
              for (k = 1; k <= RF_rTargetFactorCount; k++) {
                RF_vimpCLSleo[i][j][k]  = dvector(1, 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[k]]]);
                for (m = 1; m <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[k]]]; m++) {
                  RF_vimpCLSleo[i][j][k][m]  = NA_REAL;
                }
              }
            }
          }
        }        
      }
      if (RF_rTargetNonFactorCount > 0) {
        localSize = (ulong) xVimpSize * RF_rTargetNonFactorCount;
        RF_vimpRGR_ = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, RF_VMP_RGR, localSize, sexpVector, sexpString);        
        RF_vimpRGRptr = (double **) new_vvector(1, xVimpSize, NRUTIL_DPTR);
        for (i = 1; i <= xVimpSize; i++) {
          RF_vimpRGRptr[i] = (RF_vimpRGR_) + ((i-1) * RF_rTargetNonFactorCount) - 1;
          for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
            RF_vimpRGRptr[i][j] = NA_REAL;
          }
        }
        RF_vimpEnsembleRGR = (double ***) new_vvector(1, xVimpSize, NRUTIL_DPTR2);
        for (j = 1; j <= xVimpSize; j++) {
          RF_vimpEnsembleRGR[j]  = (double **) new_vvector(1, RF_rTargetNonFactorCount, NRUTIL_DPTR);
          for (k = 1; k <= RF_rTargetNonFactorCount; k++) {
            RF_vimpEnsembleRGR[j][k] = dvector(1, obsSize);
            for (m = 1; m <= obsSize; m++) {
              RF_vimpEnsembleRGR[j][k][m] = 0;
            }
          }
        }
        if(RF_opt & OPT_VIMP_LEOB) {
          RF_perfRGRleo = (double **) new_vvector(1, RF_forestSize, NRUTIL_DPTR);
          for (i = 1; i <= RF_forestSize; i++) {
            RF_perfRGRleo[i] = dvector(1, RF_rTargetNonFactorCount); 
            for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
              RF_perfRGRleo[i][j] = NA_REAL;
            }
          }
          RF_vimpRGRleo = (double ***) new_vvector(1, RF_forestSize, NRUTIL_DPTR2);
          for (i = 1; i <= RF_forestSize; i++) {
            RF_vimpRGRleo[i] = (double **) new_vvector(1, xVimpSize, NRUTIL_DPTR);
            for (j = 1; j <= xVimpSize; j++) {            
              RF_vimpRGRleo[i][j] = dvector(1, RF_rTargetNonFactorCount);             
              for (k= 1; k <= RF_rTargetNonFactorCount; k++) {
                RF_vimpRGRleo[i][j][k] = NA_REAL;
              }
            }
          }
        }
      }
    }
  }  
  if (RF_opt & OPT_PROX) {
    localSize = ((ulong) (obsSize + 1) * obsSize) >> 1;
    *pRF_proximity = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, RF_PROX_ID, localSize, sexpVector, sexpString);
    RF_proximityDen = dvector(1, localSize);
    (*pRF_proximity) --;
    RF_proximityPtr = (double **) new_vvector(1, obsSize, NRUTIL_DPTR);
    RF_proximityDenPtr = (double **) new_vvector(1, obsSize, NRUTIL_DPTR);
    RF_proximityPtr[1] = *pRF_proximity;
    RF_proximityDenPtr[1] = RF_proximityDen;    
    RF_proximityPtr[1][1] = RF_proximityDenPtr[1][1] = 0;
    for (i = 2; i <= obsSize; i++) {
      RF_proximityPtr[i] = RF_proximityPtr[i-1] + i - 1;
      RF_proximityDenPtr[i] = RF_proximityDenPtr[i-1] + i - 1;
      for (j = 1; j <= i; j++) {
        RF_proximityPtr[i][j] = 0;
        RF_proximityDenPtr[i][j] = 0;
      }
    }
  }
  if (RF_opt & OPT_LEAF) {
    localSize = (ulong) RF_forestSize;
    *pRF_tLeafCount = (uint*) stackAndProtect(sexpIndex, SEXP_TYPE_INTEGER, RF_LEAF_ID, localSize, sexpVector, sexpString);
    (*pRF_tLeafCount) --;
    for (i = 1; i <= RF_forestSize; i++) {
      (*pRF_tLeafCount)[i] = 0;
    }
    RF_tLeafCount = *pRF_tLeafCount;
  }
  if (RF_opt & OPT_SEED) {
    localSize = (ulong) RF_forestSize;
    *pRF_seed = (int*) stackAndProtect(sexpIndex, SEXP_TYPE_INTEGER, RF_SEED_ID, localSize, sexpVector, sexpString);
    (*pRF_seed) --;
    for (i = 1; i <= RF_forestSize; i++) {
      (*pRF_seed)[i] = -1;
    }
  }
  if (RF_opt & OPT_MISS) {
    localSize = (ulong) (RF_xSize + rspSize + 1) * mRecordSize;
    *p_imputation = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, RF_MISS_ID, localSize, sexpVector, sexpString);
    if (rspSize > 0) {
      *pRF_sImputeResponsePtr = (double **) new_vvector(1, rspSize, NRUTIL_DPTR);
      for (i = 1; i <= rspSize; i++) {
        (*pRF_sImputeResponsePtr)[i]  = (*p_imputation)  + (i * mRecordSize) - 1;
      }
    }
    *pRF_sImputePredictorPtr = (double **) new_vvector(1, RF_xSize, NRUTIL_DPTR);
    for (i = 1; i <= RF_xSize; i++) {
      (*pRF_sImputePredictorPtr)[i]  = (*p_imputation)  + ((rspSize + i) * mRecordSize) - 1;
    }
    for (i = 1; i <= mRecordSize; i++) {
      (*p_imputation)[i-1] = (double) mRecordIndex[i];
      if (rspSize > 0) {
        for (j = 1; j <= rspSize; j++) {
          (*pRF_sImputeResponsePtr)[j][i] = responsePtr[j][mRecordIndex[i]];
        }
      }
      for (j = 1; j <= RF_xSize; j++) {
        (*pRF_sImputePredictorPtr)[j][i] = predictorPtr[j][mRecordIndex[i]];
      }
    }
  }
  if (RF_opt & (OPT_VARUSED_F | OPT_VARUSED_T)) {
    localSize = (ulong) vuseDimOne * RF_xSize;
    *pRF_varUsed = (uint*) stackAndProtect(sexpIndex, SEXP_TYPE_INTEGER, RF_VUSE_ID, localSize, sexpVector, sexpString);
    if (RF_opt & OPT_VARUSED_T) {
      *pRF_varUsedPtr = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
      for (i = 1; i <= RF_forestSize; i++) {
        (*pRF_varUsedPtr)[i] = (*pRF_varUsed) + ((i-1)*(RF_xSize)) - 1;
      }
    }
    else {
      *pRF_varUsedPtr = uimatrix(1, RF_forestSize, 1, RF_xSize);
    }
    for (i = 1; i <= RF_forestSize; i++) {
      for (j = 1; j <= RF_xSize; j++) {
        (*pRF_varUsedPtr)[i][j] = 0;
      }
    }
    (*pRF_varUsed) --;
  }
  if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
    localSize = (ulong) dpthDimOne * RF_xSize * RF_observationSize;
    *p_splitDepth = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, RF_DPTH_ID, localSize, sexpVector, sexpString);
    RF_splitDepthPtr  = (double ***) new_vvector(1, dpthDimOne, NRUTIL_DPTR2);
    for (j = 1; j <= dpthDimOne; j++) {
      RF_splitDepthPtr[j] = (double **) new_vvector(1, RF_xSize, NRUTIL_DPTR);
      for (k = 1; k <= RF_xSize; k++) {
        RF_splitDepthPtr[j][k]  = (*p_splitDepth) +
          ((j-1) * RF_xSize * RF_observationSize) +
          ((k-1) * RF_observationSize) - 1;
      }
    }
    for (i = 1; i <= RF_observationSize; i++) {
      for (j = 1; j <= dpthDimOne; j++) {
        for (k = 1; k <= RF_xSize; k++) {
          RF_splitDepthPtr[j][k][i] = 0;
        }
      }
    }
  }
  if (RF_optHigh & OPT_MEMB_PRUN) {
    localSize = (ulong) RF_forestSize * obsSize;
    RF_PRUN_ID_ = (uint*) stackAndProtect(sexpIndex, SEXP_TYPE_INTEGER, RF_PRUN_ID, localSize, sexpVector, sexpString);
    RF_PRUN_ID_ptr = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
    for (i = 1; i <= RF_forestSize; i++) {
      (RF_PRUN_ID_ptr)[i] = (RF_PRUN_ID_) + ((i-1) * obsSize) - 1;
      for (j = 1; j <= obsSize; j++) {
        (RF_PRUN_ID_ptr)[i][j] = 0;
      }
    }
  }
  if (RF_optHigh & OPT_MEMB_USER) {
    localSize = (ulong) RF_forestSize * obsSize;
    RF_MEMB_ID_ = (uint*) stackAndProtect(sexpIndex, SEXP_TYPE_INTEGER, RF_MEMB_ID, localSize, sexpVector, sexpString);
    RF_MEMB_ID_ptr = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
    for (i = 1; i <= RF_forestSize; i++) {
      (RF_MEMB_ID_ptr)[i] = (RF_MEMB_ID_) + ((i-1) * obsSize) - 1;
      for (j = 1; j <= obsSize; j++) {
        (RF_MEMB_ID_ptr)[i][j] = 0;
      }
    }
    localSize = (ulong) RF_forestSize * RF_observationSize;
    RF_BOOT_CT_ = (uint*) stackAndProtect(sexpIndex, SEXP_TYPE_INTEGER, RF_BOOT_CT, localSize, sexpVector, sexpString);
    RF_BOOT_CT_ptr = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
    for (i = 1; i <= RF_forestSize; i++) {
      (RF_BOOT_CT_ptr)[i] = (RF_BOOT_CT_) + ((i-1) * RF_observationSize) - 1;
      for (j = 1; j <= RF_observationSize; j++) {
        (RF_BOOT_CT_ptr)[i][j] = 0;
      }
    }
  }
  if (RF_optHigh & OPT_PART_PLOT) {
    RF_partSURVptr = NULL;
    RF_partCLASptr = NULL;
    RF_partREGRptr = NULL;
    if ((RF_partialXvar < 1) || (RF_partialXvar > RF_xSize)) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Partial x-variable is out of range:  %10d ", RF_partialXvar);
      error("\nRF-SRC:  The application will now exit.\n");
    }
    if (RF_partialLength2 > 0) {
      for (i = 1; i <= RF_partialLength2; i++) {
        if ((RF_partialXvar2[i] < 1) || (RF_partialXvar2[i] > RF_xSize)) {
          RFprintf("\nRF-SRC:  *** ERROR *** ");
          RFprintf("\nRF-SRC:  Second order partial x-variable is out of range:  idx = %10d, val =  %10d ", i, RF_partialXvar2[i]);
          error("\nRF-SRC:  The application will now exit.\n");
        }
        if (RF_partialXvar2[i] == RF_partialXvar) {
          RFprintf("\nRF-SRC:  *** ERROR *** ");
          RFprintf("\nRF-SRC:  First and Second order partial x-variables are identical:  idx = %10d, val =  %10d ", i, RF_partialXvar2[i]);
          error("\nRF-SRC:  The application will now exit.\n");
        }
        for (j = i + 1; j <= RF_partialLength2; j++) {
          if (RF_partialXvar2[i] == RF_partialXvar2[j]) {
            RFprintf("\nRF-SRC:  *** ERROR *** ");
            RFprintf("\nRF-SRC:  Second order partial x-variables are not unique:  idx = %10d, idx =  %10d, val = %10d", i, j, RF_partialXvar2[i]);
            error("\nRF-SRC:  The application will now exit.\n");
          }
        }
      }
    }
    RF_partMembership = (Terminal ***) new_vvector(1, RF_forestSize, NRUTIL_TPTR2);
    for (i = 1; i <= RF_forestSize; i++) {
      RF_partMembership[i] = NULL;
    }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      {
        RF_partialTimeLength = RF_sortedTimeInterestSize;
        RF_partialTime = RF_timeInterest;
        if ((RF_partialTimeLength < 1) || (RF_partialTimeLength > RF_timeInterestSize)) {
          RFprintf("\nRF-SRC:  *** ERROR *** ");
          RFprintf("\nRF-SRC:  Partial time length is out of range:  %10d ", RF_partialTimeLength);
          error("\nRF-SRC:  The application will now exit.\n");
        }
        if (RF_eventTypeSize > 1) {
          if ((RF_partialType != RF_PART_YRLS) &&
              (RF_partialType != RF_PART_CIFN) &&
              (RF_partialType != RF_PART_CHFN)) {
            RFprintf("\nRF-SRC:  *** ERROR *** ");
            RFprintf("\nRF-SRC:  Partial type is out of range:  %10d ", RF_partialType);
            error("\nRF-SRC:  The application will now exit.\n");
          }
        }
        else {
          if ((RF_partialType != RF_PART_MORT) &&
              (RF_partialType != RF_PART_NLSN) &&
              (RF_partialType != RF_PART_SURV)) {
            RFprintf("\nRF-SRC:  *** ERROR *** ");
            RFprintf("\nRF-SRC:  Partial type is out of range:  %10d ", RF_partialType);
            error("\nRF-SRC:  The application will now exit.\n");
          }
        }
        dimThree = RF_partialTimeLength;
        if (!(RF_opt & OPT_COMP_RISK)) {
          if (RF_partialType == RF_PART_MORT) {
            dimThree = 1;
          }
        }
        else {
          if (RF_partialType == RF_PART_YRLS) {
            dimThree = 1;
          }
        }
        localSize = (ulong) RF_partialLength * RF_eventTypeSize * dimThree * RF_observationSize;
        RF_partial_SURV_ = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, RF_PART_SR, localSize, sexpVector, sexpString);
        RF_partSURVptr = (double ****) new_vvector(1, RF_partialLength, NRUTIL_DPTR3);
        for (j = 1; j <= RF_partialLength; j++) {
          RF_partSURVptr[j] = (double ***) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR2);
          for (k = 1; k <= RF_eventTypeSize; k++) {
            RF_partSURVptr[j][k] = (double **) new_vvector(1, dimThree, NRUTIL_DPTR);
            for (m = 1; m <= dimThree; m++) {
              RF_partSURVptr[j][k][m] = RF_partial_SURV_ +
                ((j-1) * RF_eventTypeSize * dimThree * RF_observationSize) +
                ((k-1) * dimThree * RF_observationSize) +
                ((m-1) * RF_observationSize) - 1;
              for (i = 1; i <= RF_observationSize; i++) {
                RF_partSURVptr[j][k][m][i] = 0;
              }
            }
          }
        }
      }
    }  
    else {
      if (RF_rTargetFactorCount > 0) {
        localSize = 0;
        for (k = 1; k <= RF_rTargetFactorCount; k++) {
          localSize += (ulong) 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[k]]];
        }
        localSize = (ulong) RF_partialLength * localSize * RF_observationSize;
        RF_partial_CLAS_ = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, RF_PART_CL, localSize, sexpVector, sexpString);
        RF_partCLASptr = (double ****) new_vvector(1, RF_partialLength, NRUTIL_DPTR3);
        localSize = 0;
        for (j = 1; j <= RF_partialLength; j++) {
          RF_partCLASptr[j] = (double ***) new_vvector(1, RF_rTargetFactorCount, NRUTIL_DPTR2);
          for (k = 1; k <= RF_rTargetFactorCount; k++) {
            RF_partCLASptr[j][k] = (double **) new_vvector(1, 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[k]]], NRUTIL_DPTR);
            for (m = 1; m <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[k]]]; m++) {
              RF_partCLASptr[j][k][m] = (RF_partial_CLAS_) + localSize - 1;
              for (i = 1; i <= RF_observationSize; i++) {
                RF_partCLASptr[j][k][m][i] = 0;
              }
              localSize += (ulong) obsSize;
            }
          }
        }
      }
      if (RF_rTargetNonFactorCount > 0) {
        localSize = (ulong) RF_partialLength * RF_rTargetNonFactorCount * RF_observationSize;
        RF_partial_REGR_ = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, RF_PART_RG, localSize, sexpVector, sexpString);
        RF_partREGRptr = (double ***) new_vvector(1, RF_partialLength, NRUTIL_DPTR2);
        for (j = 1; j <= RF_partialLength; j++) {
          RF_partREGRptr[j] = (double **) new_vvector(1, RF_rTargetNonFactorCount, NRUTIL_DPTR);
          for (k = 1; k <= RF_rTargetNonFactorCount; k++) {
            RF_partREGRptr[j][k] = (RF_partial_REGR_) +
              ((j-1) * RF_rTargetNonFactorCount * RF_observationSize) +
              ((k-1) * RF_observationSize) - 1;
            for (i = 1; i <= RF_observationSize; i++) {
              RF_partREGRptr[j][k][i] = 0;
            }
          }
        }
      }
    }
  }  
  return (*sexpIndex);
}
void unstackDefinedOutputObjects(char      mode,
                                 Node    **root) {
  uint obsSize;
  uint xVimpSize;
  ulong proximitySize;
  char oobFlag, fullFlag;
  uint rspSize;
  uint dpthDimOne;
  uint     **ensembleDen;
  double ****ensembleSRGptr;
  double  ***ensembleMRTptr;
  double  ***ensembleSRVptr;
  double ****ensembleCIFptr;
  double ****ensembleCLSptr;
  double  ***ensembleRGRptr;
  double ****ensembleSRGnum;
  double  ***ensembleMRTnum;
  double  ***ensembleSRVnum;
  double ****ensembleCIFnum;
  double ****ensembleCLSnum;
  double  ***ensembleRGRnum;
  char maxVoteFlag;
  uint dimThree;
  uint i, j, k;
  obsSize        = 0;  
  xVimpSize      = 0;  
  proximitySize  = 0;  
  rspSize        = 0;  
  dpthDimOne     = 0;  
  if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
    if (RF_opt & OPT_SPLDPTH_F) {
      dpthDimOne = 1;
    }
    else {
      dpthDimOne = RF_forestSize;
    }
  }
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    rspSize = RF_frSize;
    break;
  default:
    obsSize = RF_observationSize;
    rspSize = RF_rSize;
    break;
  }
  if (RF_opt & OPT_PROX) {
    proximitySize = ((obsSize + 1)  * obsSize) >> 1;
  }
  oobFlag = fullFlag = FALSE;
  if ((RF_opt & OPT_FENS) | (RF_opt & OPT_OENS)) {
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    if (RF_opt & OPT_OENS) {
      oobFlag = TRUE;
    }
    while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
      if (oobFlag == TRUE) {
        ensembleDen    = &RF_oobEnsembleDen;
        ensembleSRGptr = &RF_oobEnsembleSRGptr;
        ensembleSRGnum = &RF_oobEnsembleSRGnum;
        ensembleMRTptr = &RF_oobEnsembleMRTptr;
        ensembleMRTnum = &RF_oobEnsembleMRTnum;        
        ensembleSRVptr = &RF_oobEnsembleSRVptr;
        ensembleSRVnum = &RF_oobEnsembleSRVnum;
        ensembleCIFptr = &RF_oobEnsembleCIFptr;
        ensembleCIFnum = &RF_oobEnsembleCIFnum;
        ensembleCLSptr = &RF_oobEnsembleCLSptr;
        ensembleCLSnum = &RF_oobEnsembleCLSnum;
        ensembleRGRptr = &RF_oobEnsembleRGRptr;
        ensembleRGRnum = &RF_oobEnsembleRGRnum;
      }
      else {
        ensembleDen    = &RF_fullEnsembleDen;
        ensembleSRGptr = &RF_fullEnsembleSRGptr;
        ensembleSRGnum = &RF_fullEnsembleSRGnum;
        ensembleMRTptr = &RF_fullEnsembleMRTptr;
        ensembleMRTnum = &RF_fullEnsembleMRTnum;        
        ensembleSRVptr = &RF_fullEnsembleSRVptr;
        ensembleSRVnum = &RF_fullEnsembleSRVnum;
        ensembleCIFptr = &RF_fullEnsembleCIFptr;
        ensembleCIFnum = &RF_fullEnsembleCIFnum;
        ensembleCLSptr = &RF_fullEnsembleCLSptr;
        ensembleCLSnum = &RF_fullEnsembleCLSnum;
        ensembleRGRptr = &RF_fullEnsembleRGRptr;
        ensembleRGRnum = &RF_fullEnsembleRGRnum;
      }
      free_uivector(*ensembleDen, 1, obsSize);
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        {
          for (j = 1; j <= RF_eventTypeSize; j++) {
            for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
              free_dvector((*ensembleSRGnum)[j][k], 1, obsSize);
            }
            free_new_vvector((*ensembleSRGptr)[j], 1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
            free_new_vvector((*ensembleSRGnum)[j], 1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
          }
          free_new_vvector(*ensembleSRGptr, 1, RF_eventTypeSize, NRUTIL_DPTR2);
          free_new_vvector(*ensembleSRGnum, 1, RF_eventTypeSize, NRUTIL_DPTR2);
          for (j = 1; j <= RF_eventTypeSize; j++) {
            free_dvector((*ensembleMRTnum)[j], 1, obsSize);
          }
          free_new_vvector(*ensembleMRTptr, 1, RF_eventTypeSize, NRUTIL_DPTR);
          free_new_vvector(*ensembleMRTnum, 1, RF_eventTypeSize, NRUTIL_DPTR);
          if (!(RF_opt & OPT_COMP_RISK)) {
            for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
              free_dvector((*ensembleSRVnum)[k], 1, obsSize);
            }
            free_new_vvector(*ensembleSRVptr, 1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
            free_new_vvector(*ensembleSRVnum, 1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
          }
          else {
            for (j = 1; j <= RF_eventTypeSize; j++) {
              for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
                free_dvector((*ensembleCIFnum)[j][k], 1, obsSize);
              }
              free_new_vvector((*ensembleCIFptr)[j], 1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
              free_new_vvector((*ensembleCIFnum)[j], 1, RF_sortedTimeInterestSize, NRUTIL_DPTR);
            }
            free_new_vvector(*ensembleCIFptr, 1, RF_eventTypeSize, NRUTIL_DPTR2);
            free_new_vvector(*ensembleCIFnum, 1, RF_eventTypeSize, NRUTIL_DPTR2);            
          }  
        }
      }  
      else {
        if (RF_rTargetFactorCount > 0) {
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              free_dvector((*ensembleCLSnum)[j][k], 1, obsSize);
            }
            free_new_vvector((*ensembleCLSptr)[j], 1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]], NRUTIL_DPTR);
            free_new_vvector((*ensembleCLSnum)[j], 1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]], NRUTIL_DPTR);
          }
          free_new_vvector(*ensembleCLSptr, 1, RF_rTargetFactorCount, NRUTIL_DPTR2);
          free_new_vvector(*ensembleCLSnum, 1, RF_rTargetFactorCount, NRUTIL_DPTR2);
        }
        if (RF_rTargetNonFactorCount > 0) {
          for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
            free_dvector((*ensembleRGRnum)[j], 1, obsSize);
          }
          free_new_vvector((*ensembleRGRptr), 1, RF_rTargetNonFactorCount, NRUTIL_DPTR);
          free_new_vvector((*ensembleRGRnum), 1, RF_rTargetNonFactorCount, NRUTIL_DPTR);        
        }
      }
      if (oobFlag == TRUE) {
        oobFlag = FALSE;
      }
      else {
        fullFlag = FALSE;
      }
    }  
  }
  if (RF_opt & OPT_PERF) {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      {
        free_new_vvector(RF_perfMRTptr, 1, RF_forestSize, NRUTIL_DPTR);
      }
    }  
    else {
      if (RF_rTargetFactorCount > 0) {
        for (i = 1; i <= RF_forestSize; i++) {
          free_new_vvector(RF_perfCLSptr[i], 1, RF_rTargetFactorCount, NRUTIL_DPTR);
        }
        free_new_vvector(RF_perfCLSptr, 1, RF_forestSize, NRUTIL_DPTR2);
      }
      if (RF_rTargetNonFactorCount > 0) {
        free_new_vvector(RF_perfRGRptr, 1, RF_forestSize, NRUTIL_DPTR);
      }
    }
  }  
  if (RF_opt & OPT_VIMP) {
    if (RF_opt & OPT_VIMP_JOIN) {
      xVimpSize = 1;
    }
    else {
      xVimpSize = RF_intrPredictorSize;
    }
    for (k = 1; k <= xVimpSize; k++) {
      free_new_vvector(RF_vimpMembership[k], 1,  RF_forestSize, NRUTIL_NPTR2);
    }
    free_new_vvector(RF_vimpMembership, 1, xVimpSize, NRUTIL_NPTR3);
    for (j = 1; j <= xVimpSize; j++) {
      free_uivector(RF_vimpEnsembleDen[j], 1, obsSize);
    }
    free_new_vvector(RF_vimpEnsembleDen, 1, xVimpSize, NRUTIL_UPTR);
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      {
        free_new_vvector(RF_vimpMRTptr, 1, xVimpSize, NRUTIL_DPTR);
        for (j = 1; j <= xVimpSize; j++) {
          for (k = 1; k <= RF_eventTypeSize; k++) {
            free_dvector(RF_vimpEnsembleMRT[j][k], 1, obsSize);
          }
          free_new_vvector(RF_vimpEnsembleMRT[j], 1, RF_eventTypeSize, NRUTIL_DPTR);
        }
        free_new_vvector(RF_vimpEnsembleMRT, 1, xVimpSize, NRUTIL_DPTR2);
        if(RF_opt & OPT_VIMP_LEOB) {
          for (i = 1; i <= RF_forestSize; i++) {
            free_dvector(RF_perfMRTleo[i], 1, RF_eventTypeSize);
          }
          free_new_vvector(RF_perfMRTleo, 1, RF_forestSize, NRUTIL_DPTR);
          for (i = 1; i <= RF_forestSize; i++) {
            for (j = 1; j <= xVimpSize; j++) {
              free_dvector(RF_vimpMRTleo[i][j], 1, RF_eventTypeSize);
            }
            free_new_vvector(RF_vimpMRTleo[i], 1, xVimpSize, NRUTIL_DPTR);
          }
          free_new_vvector(RF_vimpMRTleo, 1, RF_forestSize, NRUTIL_DPTR2);
        }
      }
    }  
    else {
      if (RF_rTargetFactorCount > 0) {
        for (i = 1; i <= xVimpSize; i++) {
          free_new_vvector(RF_vimpCLSptr[i], 1, RF_rTargetFactorCount, NRUTIL_DPTR);
        }
        free_new_vvector(RF_vimpCLSptr, 1, xVimpSize, NRUTIL_DPTR2);
        maxVoteFlag = FALSE;
        if(RF_opt & OPT_VIMP_LEOB) {
          if (!(RF_opt & OPT_PERF_CALB)) {
            maxVoteFlag = TRUE;
          }
        }
        if (maxVoteFlag) {
          for (i = 1; i <= xVimpSize; i++) {
            for (j = 1; j <= RF_rTargetFactorCount; j++) {
              free_dvector(RF_vimpEnsembleCLS[i][j][1], 1, obsSize);
              free_new_vvector(RF_vimpEnsembleCLS[i][j], 1, 1, NRUTIL_DPTR);
            }
            free_new_vvector(RF_vimpEnsembleCLS[i], 1, RF_rTargetFactorCount, NRUTIL_DPTR2);
          }
          free_new_vvector(RF_vimpEnsembleCLS, 1, xVimpSize, NRUTIL_DPTR3);
        }
        else {
          for (i = 1; i <= xVimpSize; i++) {
            for (j = 1; j <= RF_rTargetFactorCount; j++) {
              for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
                free_dvector(RF_vimpEnsembleCLS[i][j][k], 1, obsSize);
              }
              free_new_vvector(RF_vimpEnsembleCLS[i][j], 1, RF_rFactorSize[RF_rTargetFactor[j]], NRUTIL_DPTR);
            }
            free_new_vvector(RF_vimpEnsembleCLS[i], 1, RF_rTargetFactorCount, NRUTIL_DPTR2);
          }
          free_new_vvector(RF_vimpEnsembleCLS, 1, xVimpSize, NRUTIL_DPTR3);        
        }
        if(RF_opt & OPT_VIMP_LEOB) {
          for (i = 1; i <= RF_forestSize; i++) {
            for (j = 1; j <= RF_rTargetFactorCount; j++) {
              free_dvector(RF_perfCLSleo[i][j], 1, 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]);
            }
            free_new_vvector(RF_perfCLSleo[i], 1, RF_rTargetFactorCount, NRUTIL_DPTR);
          }
          free_new_vvector(RF_perfCLSleo, 1, RF_forestSize, NRUTIL_DPTR2);
          for (i = 1; i <= RF_forestSize; i++) {
            for (j = 1; j <= xVimpSize; j++) {            
              for (k = 1; k <= RF_rTargetFactorCount; k++) {
                free_dvector(RF_vimpCLSleo[i][j][k], 1, 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[k]]]);
              }
              free_new_vvector(RF_vimpCLSleo[i][j], 1, RF_rTargetFactorCount, NRUTIL_DPTR);
            }
            free_new_vvector(RF_vimpCLSleo[i], 1, xVimpSize, NRUTIL_DPTR2);
          }
          free_new_vvector(RF_vimpCLSleo, 1, RF_forestSize, NRUTIL_DPTR3);
        }
      }
      if (RF_rTargetNonFactorCount > 0) {
        free_new_vvector(RF_vimpRGRptr, 1, xVimpSize, NRUTIL_DPTR);
        for (j = 1; j <= xVimpSize; j++) {
          for (k = 1; k <= RF_rTargetNonFactorCount; k++) {
            free_dvector(RF_vimpEnsembleRGR[j][k], 1, obsSize);
          }
          free_new_vvector(RF_vimpEnsembleRGR[j], 1, RF_rTargetNonFactorCount, NRUTIL_DPTR);
        }
        free_new_vvector(RF_vimpEnsembleRGR, 1, xVimpSize, NRUTIL_DPTR2);
        if(RF_opt & OPT_VIMP_LEOB) {
          for (i = 1; i <= RF_forestSize; i++) {
            free_dvector(RF_perfRGRleo[i], 1, RF_rTargetNonFactorCount); 
          }
          free_new_vvector(RF_perfRGRleo, 1, RF_forestSize, NRUTIL_DPTR);
          for (i = 1; i <= RF_forestSize; i++) {
            for (j = 1; j <= xVimpSize; j++) {            
              free_dvector(RF_vimpRGRleo[i][j], 1, RF_rTargetNonFactorCount);             
            }
            free_new_vvector(RF_vimpRGRleo[i], 1, xVimpSize, NRUTIL_DPTR);
          }
          free_new_vvector(RF_vimpRGRleo, 1, RF_forestSize, NRUTIL_DPTR2);
        }
      }
    }
  }  
  if (RF_opt & OPT_PROX) {
    free_dvector(RF_proximityDen, 1, proximitySize);
    free_new_vvector(RF_proximityPtr, 1, obsSize, NRUTIL_DPTR);
    free_new_vvector(RF_proximityDenPtr, 1, obsSize, NRUTIL_DPTR);
  }
  if (RF_opt & OPT_MISS) {
    if (rspSize > 0) {
      free_new_vvector(RF_sImputeResponsePtr, 1, rspSize, NRUTIL_DPTR);
    }
    free_new_vvector(RF_sImputePredictorPtr, 1, RF_xSize, NRUTIL_DPTR);
  }
  if (RF_opt & (OPT_VARUSED_F | OPT_VARUSED_T)) {
    if (RF_opt & OPT_VARUSED_T) {
      free_new_vvector(RF_varUsedPtr, 1, RF_forestSize, NRUTIL_UPTR);
    }
    else {
      free_uimatrix(RF_varUsedPtr, 1, RF_forestSize, 1, RF_xSize);
    }
  }
  if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
    for (j = 1; j <= dpthDimOne; j++) {
      free_new_vvector(RF_splitDepthPtr[j], 1, RF_xSize, NRUTIL_DPTR);
    }
    free_new_vvector(RF_splitDepthPtr, 1, dpthDimOne, NRUTIL_DPTR2);
  }
  if (RF_optHigh & OPT_MEMB_PRUN) {
    free_new_vvector(RF_PRUN_ID_ptr, 1, RF_forestSize, NRUTIL_UPTR);
  }
  if (RF_optHigh & OPT_MEMB_USER) {
    free_new_vvector(RF_MEMB_ID_ptr, 1, RF_forestSize, NRUTIL_UPTR);
    free_new_vvector(RF_BOOT_CT_ptr, 1, RF_forestSize, NRUTIL_UPTR);
  }
  if (RF_optHigh & OPT_PART_PLOT) {
    free_new_vvector(RF_partMembership, 1, RF_forestSize, NRUTIL_TPTR2);
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      {
        dimThree = RF_partialTimeLength;
        if (!(RF_opt & OPT_COMP_RISK)) {
          if (RF_partialType == RF_PART_MORT) {
            dimThree = 1;
          }
        }
        else {
          if (RF_partialType == RF_PART_YRLS) {
            dimThree = 1;
          }
        }
        for (j = 1; j <= RF_partialLength; j++) {
          for (k = 1; k <= RF_eventTypeSize; k++) {
            free_new_vvector(RF_partSURVptr[j][k], 1, dimThree, NRUTIL_DPTR);
          }
          free_new_vvector(RF_partSURVptr[j], 1, RF_eventTypeSize, NRUTIL_DPTR2);
        }
        free_new_vvector(RF_partSURVptr, 1, RF_partialLength, NRUTIL_DPTR3);
      }
    }
    else {
      if (RF_rTargetFactorCount > 0) {
        for (j = 1; j <= RF_partialLength; j++) {
          for (k = 1; k <= RF_rTargetFactorCount; k++) {
            free_new_vvector(RF_partCLASptr[j][k], 1, 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[k]]], NRUTIL_DPTR);
          }
          free_new_vvector(RF_partCLASptr[j], 1, RF_rTargetFactorCount, NRUTIL_DPTR2);
        }
        free_new_vvector(RF_partCLASptr, 1, RF_partialLength, NRUTIL_DPTR3);
      }
      if (RF_rTargetNonFactorCount > 0) {
        for (j = 1; j <= RF_partialLength; j++) {
          free_new_vvector(RF_partREGRptr[j], 1, RF_rTargetNonFactorCount, NRUTIL_DPTR);
        }
        free_new_vvector(RF_partREGRptr, 1, RF_partialLength, NRUTIL_DPTR2);
      }
    }
  }
  switch (mode) {
  case RF_PRED:
    if (RF_rSize == 0) {
    }
    else {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      }
      else {
        free_uivector(RF_rTargetFactor, 1, RF_rTargetCount);
        free_uivector(RF_rTargetNonFactor, 1, RF_rTargetCount);
      }
    }
    break;
  default:
    if (RF_rSize == 0) {
    }
    else {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      }
      else {
        if (mode == RF_GROW) {
          free_uivector(RF_rTarget, 1 , RF_rTargetCount);
        }
        free_uivector(RF_rTargetFactor, 1, RF_rTargetCount);
        free_uivector(RF_rTargetNonFactor, 1, RF_rTargetCount);
      }
    }
    break;
  }
}
uint stackForestOutputObjects(char     mode,
                                uint   **pRF_treeID,
                                uint   **pRF_nodeID,
                                uint   **pRF_parmID,
                                double **pRF_contPT,
                                uint   **pRF_mwcpSZ,
                                uint   **pRF_mwcpPT,
                                uint   **pRF_mwcpCT,
                                uint    *sexpIndex,
                                char   **sexpString,
                                SEXP    *sexpVector) {
  ulong totalNodeCount;
  ulong totalMWCPCount;
  ulong localSize;
  uint  mwcpSize;
  uint i;
  if (mode == RF_GROW) {
    if (RF_opt & OPT_TREE) {
      RF_totalNodeCount1  = 1;
      totalNodeCount = (ulong) ((RF_theoreticalMaxtLeafCount[1] << 1) - 1) * RF_forestSize;
      if (RF_xFactorCount > 0) {
        mwcpSize = (RF_xMaxFactorLevel >> (3 + ulog2(sizeof(uint)))) + ((RF_xMaxFactorLevel & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
      }
      else {
        mwcpSize = 0;
      }
      totalMWCPCount = (ulong) mwcpSize * totalNodeCount;
      *pRF_treeID = (uint*)   stackAndProtect(sexpIndex, SEXP_TYPE_INTEGER, RF_TREE_ID, totalNodeCount, sexpVector, sexpString);
      *pRF_nodeID = (uint*)   stackAndProtect(sexpIndex, SEXP_TYPE_INTEGER, RF_NODE_ID, totalNodeCount, sexpVector, sexpString);
      *pRF_parmID = (uint*)   stackAndProtect(sexpIndex, SEXP_TYPE_INTEGER, RF_PARM_ID, totalNodeCount, sexpVector, sexpString);
      *pRF_contPT = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, RF_CONT_PT, totalNodeCount, sexpVector, sexpString);
      *pRF_mwcpSZ = (uint*)   stackAndProtect(sexpIndex, SEXP_TYPE_INTEGER, RF_MWCP_SZ, totalNodeCount, sexpVector, sexpString);
      *pRF_mwcpPT = (uint*)   stackAndProtect(sexpIndex, SEXP_TYPE_INTEGER, RF_MWCP_PT, totalMWCPCount, sexpVector, sexpString);
      (*pRF_treeID) --;
      (*pRF_nodeID) --;
      (*pRF_parmID) --;
      (*pRF_contPT) --;
      (*pRF_mwcpSZ) --;
      (*pRF_mwcpPT) --;
      RF_mwcpIterator = (*pRF_mwcpPT);
      localSize = (ulong) RF_forestSize;
      *pRF_mwcpCT = (uint*)   stackAndProtect(sexpIndex, SEXP_TYPE_INTEGER, RF_MWCP_CT, localSize, sexpVector, sexpString);
      (*pRF_mwcpCT) --;
      for (i = 1; i <= RF_forestSize; i++) {
        (*pRF_mwcpCT)[i] = 0;
      }
    }
  }
  return (*sexpIndex);
}
uint stackStatisticalOutputObjects(char     mode,
                                   double **pRF_spltST,
                                   double **pRF_spltVR,
                                   uint   **pRF_uspvST,
                                   uint   **pRF_mtryID,
                                   double **pRF_mtryST,
                                   uint    *sexpIndex,
                                   char   **sexpString,
                                   SEXP    *sexpVector) {
  ulong totalNodeCount;
  ulong localSize;
  uint i;
  RF_totalNodeCount2 = 1;
  switch (mode) {
  case RF_GROW:
    totalNodeCount = ((RF_theoreticalMaxtLeafCount[1] << 1) - 1) * RF_forestSize;
    if (RF_opt & OPT_USPV_STAT) {
      localSize = (ulong) totalNodeCount * RF_randomResponseCount;
      *pRF_uspvST = (uint*) stackAndProtect(sexpIndex, SEXP_TYPE_INTEGER, RF_USPV_ST, localSize, sexpVector, sexpString);
      RF_uspvST_ptr = (uint**) new_vvector(1, totalNodeCount, NRUTIL_UPTR);
      for (i = 1; i <= totalNodeCount; i++) {
        RF_uspvST_ptr[i] = (*pRF_uspvST) + (i-1) * RF_randomResponseCount - 1;
      }
    }
    break;
  default:
    totalNodeCount = RF_totalNodeCount;
    break;
  }
  if (RF_opt & OPT_NODE_STAT) {
    localSize = totalNodeCount;
    *pRF_spltST = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, RF_SPLT_ST, localSize, sexpVector, sexpString);
    (*pRF_spltST) --;
    *pRF_spltVR = NULL;
    if (mode == RF_GROW) {
      localSize = (ulong) totalNodeCount * RF_randomCovariateCount;
      *pRF_mtryID = (uint*) stackAndProtect(sexpIndex, SEXP_TYPE_INTEGER, RF_MTRY_ID, localSize, sexpVector, sexpString);
      RF_mtryID_ptr = (uint**) new_vvector(1, totalNodeCount, NRUTIL_UPTR);
      for (ulong ui = 1; ui <= totalNodeCount; ui++) {
        RF_mtryID_ptr[ui] = (*pRF_mtryID) + (ui-1) * RF_randomCovariateCount - 1;
      }
      *pRF_mtryST = (double*) stackAndProtect(sexpIndex, SEXP_TYPE_NUMERIC, RF_MTRY_ST, localSize, sexpVector, sexpString);
      RF_mtryST_ptr = (double**) new_vvector(1, totalNodeCount, NRUTIL_DPTR);
      for (ulong ui = 1; ui <= totalNodeCount; ui++) {
        RF_mtryST_ptr[ui] = (*pRF_mtryST) + (ui-1) * RF_randomCovariateCount - 1;
      }
    }
  }
  return (*sexpIndex);
}
void unstackAuxStatisticalStructures (char mode) {
  ulong totalNodeCount;
  if (mode == RF_GROW) {
    totalNodeCount = (ulong) ((RF_theoreticalMaxtLeafCount[1] << 1) - 1) * RF_forestSize;
    if (RF_opt & OPT_USPV_STAT) {
      free_new_vvector(RF_uspvST_ptr, 1, totalNodeCount, NRUTIL_UPTR);
    }
    if (RF_opt & OPT_NODE_STAT) {
      free_new_vvector(RF_mtryID_ptr, 1, totalNodeCount, NRUTIL_UPTR);    
      free_new_vvector(RF_mtryST_ptr, 1, totalNodeCount, NRUTIL_DPTR);
    }
  }
}
uint stackTNQualitativeOutputObjects(char     mode,
                                     uint   **pRF_RMBR_ID_,
                                     uint   **pRF_AMBR_ID_,
                                     uint   **pRF_TN_RCNT_,
                                     uint   **pRF_TN_ACNT_,
                                     uint     sexpIndex,
                                     char   **sexpString,
                                     SEXP    *sexpVector) {
  ulong localSize;
  if (RF_optHigh & OPT_MEMB_OUTG) {
    localSize = (ulong) RF_forestSize * RF_bootstrapSize;
    *pRF_RMBR_ID_ = (uint*) stackAndProtect(&sexpIndex, SEXP_TYPE_INTEGER, RF_RMBR_ID, localSize, sexpVector, sexpString);
    localSize = (ulong) RF_forestSize * RF_observationSize;
    *pRF_AMBR_ID_ = (uint*) stackAndProtect(&sexpIndex, SEXP_TYPE_INTEGER, RF_AMBR_ID, localSize, sexpVector, sexpString);
    localSize = RF_totalTerminalCount;
    *pRF_TN_RCNT_ = (uint*) stackAndProtect(&sexpIndex, SEXP_TYPE_INTEGER, RF_TN_RCNT, localSize, sexpVector, sexpString);
    *pRF_TN_ACNT_ = (uint*) stackAndProtect(&sexpIndex, SEXP_TYPE_INTEGER, RF_TN_ACNT, localSize, sexpVector, sexpString);
  }
  return (sexpIndex);
}
void stackAuxTNQualitativeStructures(char     mode,
                                     uint    *pRF_RMBR_ID_,
                                     uint    *pRF_AMBR_ID_,
                                     uint    *pRF_TN_RCNT_,
                                     uint    *pRF_TN_ACNT_) {
  uint *tLeafCount;
  ulong offset;
  uint i, j;
  if ((RF_optHigh & OPT_MEMB_INCG) || (RF_optHigh & OPT_MEMB_OUTG)) {
    RF_RMBR_ID_ptr = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
    for (i = 1; i <= RF_forestSize; i++) {
      (RF_RMBR_ID_ptr)[i] = (pRF_RMBR_ID_) + ((i-1) * RF_bootstrapSize) - 1;
    }
    if (RF_optHigh & OPT_MEMB_OUTG) {
      for (i = 1; i <= RF_forestSize; i++) {
        for (j = 1; j <= RF_bootstrapSize; j++) {
          (RF_RMBR_ID_ptr)[i][j] = 0;
        }
      }
    }
    RF_AMBR_ID_ptr = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
    for (i = 1; i <= RF_forestSize; i++) {
      (RF_AMBR_ID_ptr)[i] = (pRF_AMBR_ID_) + ((i-1) * RF_observationSize) - 1;
    }
    if (RF_optHigh & OPT_MEMB_OUTG) {
      for (i = 1; i <= RF_forestSize; i++) {
        for (j = 1; j <= RF_observationSize; j++) {
          (RF_AMBR_ID_ptr)[i][j] = 0;
        }
      }
    }
    if (mode == RF_GROW) {
      tLeafCount = RF_theoreticalMaxtLeafCount;
    }
    else {
      tLeafCount = RF_tLeafCount;
    }
    offset = 0;
    RF_TN_RCNT_ptr = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
    RF_TN_ACNT_ptr = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
    for (i = 1; i <= RF_forestSize; i++) {
      if (tLeafCount[i] > 0) {
        RF_TN_RCNT_ptr[i] = (pRF_TN_RCNT_) + offset - 1;
        RF_TN_ACNT_ptr[i] = (pRF_TN_ACNT_) + offset - 1;
        offset += (ulong) tLeafCount[i];
      }
      else {
        RF_TN_RCNT_ptr[i] = NULL;
        RF_TN_ACNT_ptr[i] = NULL;
      }
    }
  }
}
void unstackAuxTNQualitativeStructures(char mode) {
  if ((RF_optHigh & OPT_MEMB_INCG) || (RF_optHigh & OPT_MEMB_OUTG)) {
    free_new_vvector(RF_RMBR_ID_ptr, 1, RF_forestSize, NRUTIL_UPTR);
    free_new_vvector(RF_AMBR_ID_ptr, 1, RF_forestSize, NRUTIL_UPTR);
    free_new_vvector(RF_TN_RCNT_ptr, 1, RF_forestSize, NRUTIL_UPTR);
    free_new_vvector(RF_TN_ACNT_ptr, 1, RF_forestSize, NRUTIL_UPTR);
  }
}
uint stackTNQuantitativeOutputObjects(char     mode,
                                      double **pRF_TN_SURV_,
                                      double **pRF_TN_MORT_,
                                      double **pRF_TN_NLSN_,
                                      double **pRF_TN_CSHZ_,
                                      double **pRF_TN_CIFN_,
                                      double **pRF_TN_REGR_,
                                      uint   **pRF_TN_CLAS_,
                                      uint     sexpIndex,
                                      char   **sexpString,
                                      SEXP    *sexpVector) {
  ulong tnDimOne, localSize;
  uint  tnDimTwo;
  uint j;
  if (RF_optHigh & OPT_TERM_OUTG) {
    tnDimOne = RF_totalTerminalCount;
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      localSize = (ulong) tnDimOne * RF_eventTypeSize;
      *pRF_TN_MORT_ = (double*) stackAndProtect(&sexpIndex, SEXP_TYPE_NUMERIC, RF_TN_MORT, localSize, sexpVector, sexpString);
      if (!(RF_opt & OPT_COMP_RISK)) {
        localSize = (ulong) tnDimOne * RF_sortedTimeInterestSize;
        *pRF_TN_SURV_ = (double*) stackAndProtect(&sexpIndex, SEXP_TYPE_NUMERIC, RF_TN_SURV, localSize, sexpVector, sexpString);
        localSize = (ulong) tnDimOne * RF_sortedTimeInterestSize;
        *pRF_TN_NLSN_ = (double*) stackAndProtect(&sexpIndex, SEXP_TYPE_NUMERIC, RF_TN_NLSN, localSize, sexpVector, sexpString);
      }
      else {
        localSize = (ulong) tnDimOne * RF_eventTypeSize * RF_sortedTimeInterestSize;
        *pRF_TN_CSHZ_ = (double*) stackAndProtect(&sexpIndex, SEXP_TYPE_NUMERIC, RF_TN_CSHZ, localSize, sexpVector, sexpString);
        *pRF_TN_CIFN_ = (double*) stackAndProtect(&sexpIndex, SEXP_TYPE_NUMERIC, RF_TN_CIFN, localSize, sexpVector, sexpString);
      }
    }
    else {
      if (RF_rNonFactorCount > 0) {
        localSize = (ulong) tnDimOne * RF_rNonFactorCount;
        *pRF_TN_REGR_ = (double*) stackAndProtect(&sexpIndex, SEXP_TYPE_NUMERIC, RF_TN_REGR, localSize, sexpVector, sexpString);
      }
      if (RF_rFactorCount > 0) {
        tnDimTwo = 0;
        for (j = 1; j <= RF_rFactorCount; j++) {
          tnDimTwo += RF_rFactorSize[j];
        }
        localSize = (ulong) tnDimOne * tnDimTwo;
        *pRF_TN_CLAS_ = (uint*) stackAndProtect(&sexpIndex, SEXP_TYPE_INTEGER, RF_TN_CLAS, localSize, sexpVector, sexpString);
      }
    }
  }
  return (sexpIndex);
}
void stackAuxTNQuantitativeStructures(char    mode,
                                      double *pRF_TN_SURV,
                                      double *pRF_TN_MORT,
                                      double *pRF_TN_NLSN,
                                      double *pRF_TN_CSHZ,
                                      double *pRF_TN_CIFN,
                                      double *pRF_TN_REGR,
                                      uint   *pRF_TN_CLAS) {
  uint *tLeafCount;
  ulong offset;
  uint  i, j, k;
  if (mode == RF_GROW) {
    tLeafCount = RF_theoreticalMaxtLeafCount;
  }
  else {
    tLeafCount = RF_tLeafCount;
  }
  if ((RF_optHigh & OPT_TERM_INCG) || (RF_optHigh & OPT_TERM_OUTG)) {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      offset = 0;
      RF_TN_MORT_ptr = (double ***) new_vvector(1, RF_forestSize, NRUTIL_DPTR2);
      for (i = 1; i <= RF_forestSize; i++) {
        if (tLeafCount[i] > 0) {
          RF_TN_MORT_ptr[i] = (double **) new_vvector(1, tLeafCount[i], NRUTIL_DPTR);
          for (j = 1; j <= tLeafCount[i]; j++) {
            RF_TN_MORT_ptr[i][j] = (pRF_TN_MORT) + offset - 1;
            offset += (ulong) RF_eventTypeSize;
          }
        }
        else {
          RF_TN_MORT_ptr[i] = NULL;
        }
      }
      if (!(RF_opt & OPT_COMP_RISK)) {
        offset = 0;
        RF_TN_SURV_ptr = (double ***) new_vvector(1, RF_forestSize, NRUTIL_DPTR2);
        for (i = 1; i <= RF_forestSize; i++) {
          if (tLeafCount[i] > 0) {
            RF_TN_SURV_ptr[i] = (double **) new_vvector(1, tLeafCount[i], NRUTIL_DPTR);
            for (j = 1; j <= tLeafCount[i]; j++) {
              RF_TN_SURV_ptr[i][j] = (pRF_TN_SURV) + offset - 1;
              offset += (ulong) RF_sortedTimeInterestSize;
            }
          }
          else {
            RF_TN_SURV_ptr[i] = NULL;
          }
        }
        offset = 0;
        RF_TN_NLSN_ptr = (double ***) new_vvector(1, RF_forestSize, NRUTIL_DPTR2);
        for (i = 1; i <= RF_forestSize; i++) {
          if (tLeafCount[i] > 0) {
            RF_TN_NLSN_ptr[i] = (double **) new_vvector(1, tLeafCount[i], NRUTIL_DPTR);
            for (j = 1; j <= tLeafCount[i]; j++) {
              RF_TN_NLSN_ptr[i][j] = (pRF_TN_NLSN) + offset - 1;
              offset += (ulong) RF_sortedTimeInterestSize;
            }
          }
          else {
            RF_TN_NLSN_ptr[i] = NULL;
          }
        }
      }
      else {
        offset = 0;
        RF_TN_CSHZ_ptr = (double ****) new_vvector(1, RF_forestSize, NRUTIL_DPTR3);
        for (i = 1; i <= RF_forestSize; i++) {
          if (tLeafCount[i] > 0) {
            RF_TN_CSHZ_ptr[i] = (double ***) new_vvector(1, tLeafCount[i], NRUTIL_DPTR2);
            for (j = 1; j <= tLeafCount[i]; j++) {
              RF_TN_CSHZ_ptr[i][j] = (double **) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR);
              for (k = 1; k <= RF_eventTypeSize; k++) {
                RF_TN_CSHZ_ptr[i][j][k] = (pRF_TN_CSHZ) + offset - 1;
                offset += (ulong) RF_sortedTimeInterestSize;
              }
            }
          }
          else {
            RF_TN_CSHZ_ptr[i] = NULL;
          }
        }
        offset = 0;
        RF_TN_CIFN_ptr = (double ****) new_vvector(1, RF_forestSize, NRUTIL_DPTR3);
        for (i = 1; i <= RF_forestSize; i++) {
          if (tLeafCount[i] > 0) {
            RF_TN_CIFN_ptr[i] = (double ***) new_vvector(1, tLeafCount[i], NRUTIL_DPTR2);
            for (j = 1; j <= tLeafCount[i]; j++) {
              RF_TN_CIFN_ptr[i][j] = (double **) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR);
              for (k = 1; k <= RF_eventTypeSize; k++) {
                RF_TN_CIFN_ptr[i][j][k] = (pRF_TN_CIFN) + offset - 1;
                offset += (ulong) RF_sortedTimeInterestSize;
              }
            }
          }
          else {
            RF_TN_CIFN_ptr[i] = NULL;
          }
        }
      }
    }
    else {
      if (RF_rNonFactorCount > 0) {
        offset = 0;
        RF_TN_REGR_ptr = (double ***) new_vvector(1, RF_forestSize, NRUTIL_DPTR2);
        for (i = 1; i <= RF_forestSize; i++) {
          if (tLeafCount[i] > 0) {
            RF_TN_REGR_ptr[i] = (double **) new_vvector(1, tLeafCount[i], NRUTIL_DPTR);
            for (j = 1; j <= tLeafCount[i]; j++) {
              RF_TN_REGR_ptr[i][j] = (pRF_TN_REGR) + offset - 1;
              offset += (ulong) RF_rNonFactorCount;
            }
          }
          else {
            RF_TN_REGR_ptr[i] = NULL;
          }
        }
      }
      if (RF_rFactorCount > 0) {
        offset = 0;
        RF_TN_CLAS_ptr = (uint ****) new_vvector(1, RF_forestSize, NRUTIL_UPTR3);
        for (i = 1; i <= RF_forestSize; i++) {
          if (tLeafCount[i] > 0) {
            RF_TN_CLAS_ptr[i] = (uint ***) new_vvector(1, tLeafCount[i], NRUTIL_UPTR2);
            for (j = 1; j <= tLeafCount[i]; j++) {
              RF_TN_CLAS_ptr[i][j] = (uint **) new_vvector(1, RF_rFactorCount, NRUTIL_UPTR);
              for (k = 1; k <= RF_rFactorCount; k++) {
                RF_TN_CLAS_ptr[i][j][k] = (pRF_TN_CLAS) + offset - 1;
                offset += (ulong) RF_rFactorSize[k];
              }
            }
          }
          else {
            RF_TN_CLAS_ptr[i] = NULL;
          }
        }
      }
    }
  }
}
void unstackAuxTNQuantitativeStructures(char mode) {
  uint *tLeafCount;
  uint i, j;
  if (mode == RF_GROW) {
    tLeafCount = RF_theoreticalMaxtLeafCount;
  }
  else {
    tLeafCount = RF_tLeafCount;
  }
  if ((RF_optHigh & OPT_TERM_INCG) || (RF_optHigh & OPT_TERM_OUTG)) {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      for (i = 1; i <= RF_forestSize; i++) {
        if (tLeafCount[i] > 0) {
          free_new_vvector(RF_TN_MORT_ptr[i], 1, tLeafCount[i], NRUTIL_DPTR);
        }
      }
      free_new_vvector(RF_TN_MORT_ptr, 1, RF_forestSize, NRUTIL_DPTR2);
      RF_TN_MORT_ptr = NULL;
      if (!(RF_opt & OPT_COMP_RISK)) {
        for (i = 1; i <= RF_forestSize; i++) {
          if (tLeafCount[i] > 0) {
            free_new_vvector(RF_TN_SURV_ptr[i], 1, tLeafCount[i], NRUTIL_DPTR);
          }
        }
        free_new_vvector(RF_TN_SURV_ptr, 1, RF_forestSize, NRUTIL_DPTR2);
        RF_TN_SURV_ptr = NULL;
        for (i = 1; i <= RF_forestSize; i++) {
          if (tLeafCount[i] > 0) {
            free_new_vvector(RF_TN_NLSN_ptr[i], 1, tLeafCount[i], NRUTIL_DPTR);
          }
        }
        free_new_vvector(RF_TN_NLSN_ptr, 1, RF_forestSize, NRUTIL_DPTR2);
        RF_TN_NLSN_ptr = NULL;
      }
      else {
        for (i = 1; i <= RF_forestSize; i++) {
          if (tLeafCount[i] > 0) {
            for (j = 1; j <= tLeafCount[i]; j++) {
              free_new_vvector(RF_TN_CSHZ_ptr[i][j], 1, RF_eventTypeSize, NRUTIL_DPTR);
            }
            free_new_vvector(RF_TN_CSHZ_ptr[i], 1, tLeafCount[i], NRUTIL_DPTR2);
          }
        }
        free_new_vvector(RF_TN_CSHZ_ptr, 1, RF_forestSize, NRUTIL_DPTR3);
        RF_TN_CSHZ_ptr = NULL;
        for (i = 1; i <= RF_forestSize; i++) {
          if (tLeafCount[i] > 0) {
            for (j = 1; j <= tLeafCount[i]; j++) {
              free_new_vvector(RF_TN_CIFN_ptr[i][j], 1, RF_eventTypeSize, NRUTIL_DPTR);
            }
            free_new_vvector(RF_TN_CIFN_ptr[i], 1, tLeafCount[i], NRUTIL_DPTR2);
          }
        }
        free_new_vvector(RF_TN_CIFN_ptr, 1, RF_forestSize, NRUTIL_DPTR3);
        RF_TN_CIFN_ptr = NULL;
      }
    }
    else {
      if (RF_rNonFactorCount > 0) {
        for (i = 1; i <= RF_forestSize; i++) {
          if (tLeafCount[i] > 0) {
            free_new_vvector(RF_TN_REGR_ptr[i], 1, tLeafCount[i], NRUTIL_DPTR);
          }
        }
        free_new_vvector(RF_TN_REGR_ptr, 1, RF_forestSize, NRUTIL_DPTR2);
        RF_TN_REGR_ptr = NULL;
      }
      if (RF_rFactorCount > 0) {
        for (i = 1; i <= RF_forestSize; i++) {
          if (tLeafCount[i] > 0) {
            for (j = 1; j <= tLeafCount[i]; j++) {
              free_new_vvector(RF_TN_CLAS_ptr[i][j], 1, RF_rFactorCount, NRUTIL_UPTR);
            }
            free_new_vvector(RF_TN_CLAS_ptr[i], 1, tLeafCount[i], NRUTIL_UPTR2);
          }
        }
        free_new_vvector(RF_TN_CLAS_ptr, 1, RF_forestSize, NRUTIL_UPTR3);
        RF_TN_CLAS_ptr = NULL;
      }
    }
  }
}
void initProtect(SEXP *sexpVector,
                 uint  stackCount) {
  PROTECT(sexpVector[RF_OUTP_ID] = allocVector(VECSXP, stackCount));
  PROTECT(sexpVector[RF_STRG_ID] = allocVector(STRSXP, stackCount));
  setAttrib(sexpVector[RF_OUTP_ID], R_NamesSymbol, sexpVector[RF_STRG_ID]);
}
void *stackAndProtect(uint  *sexpIndex,
                      char   sexpType,
                      uint   sexpIdentity,
                      ulong  size,
                      SEXP  *sexpVector,
                      char **sexpString) {
  void *v;
  if (size == 0) {
  }
  if (((*sexpIndex) >> 6) > 0) {
          RFprintf("\nRF-SRC:  *** ERROR *** ");
          RFprintf("\nRF-SRC:  SEXP vector list limit exceeded:  %20d", *sexpIndex);
          RFprintf("\nRF-SRC:  Please Contact Technical Support.");
          error("\nRF-SRC:  The application will now exit.\n");
  }
  if (sizeof(ulong) > sizeof(uint)) {
    if (size > UINT_MAX) {
      if (TRUE) {
        RFprintf("\nRF-SRC:  *** WARNING *** ");
        RFprintf("\nRF-SRC:  SEXP vector element length exceeds 32-bits:  %20lu", size);
        RFprintf("\nRF-SRC:  SEXP ALLOC:  %s ", sexpString[sexpIdentity]);
        RFprintf("\nRF-SRC:  Please Reduce Dimensionality if Possible.");
      }
    }
  }
  switch(sexpType) {
  case SEXP_TYPE_NUMERIC:
    PROTECT(sexpVector[sexpIdentity] = NEW_NUMERIC(size));
    break;
  case SEXP_TYPE_INTEGER:
    PROTECT(sexpVector[sexpIdentity] = NEW_INTEGER(size));
    break;
  case SEXP_TYPE_CHARACTER:
    PROTECT(sexpVector[sexpIdentity] = NEW_CHARACTER(size));
    break;
  default:
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  SEXP vector element type unknown:  %20d", sexpType);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
    break;
  }
  SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], *sexpIndex, sexpVector[sexpIdentity]);
  SET_STRING_ELT(sexpVector[RF_STRG_ID], *sexpIndex, mkChar(sexpString[sexpIdentity]));
  (*sexpIndex) ++;
  switch(sexpType) {
  case SEXP_TYPE_NUMERIC:
    v = (double*) NUMERIC_POINTER(sexpVector[sexpIdentity]);
    break;
  case SEXP_TYPE_INTEGER:
    v = (uint*) INTEGER_POINTER(sexpVector[sexpIdentity]);
    break;
  case SEXP_TYPE_CHARACTER:
    v = (char*) CHARACTER_POINTER(sexpVector[sexpIdentity]);
    break;
  default:
    v = NULL;
    break;
  }
  return v;
}
void verifyAndRegisterCustomSplitRules() {
  uint familyConstant;
  uint i, j;
  familyConstant = 0;  
  if (RF_splitRule == CUST_SPLIT) {
    RF_splitCustomIdx = (RF_optHigh & OPT_SPLT_CUST) >> 8;
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 16; j++) {
        customFunctionArray[i][j] = NULL;
      }
    }
    registerCustomFunctions();
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      {
        if (!(RF_opt & OPT_COMP_RISK)) {
          familyConstant = SURV_FAM;
        }
        else {
          familyConstant = CRSK_FAM;
        }
        if (customFunctionArray[familyConstant][RF_splitCustomIdx] == NULL) {
          RFprintf("\nRF-SRC:  *** ERROR *** ");
          RFprintf("\nRF-SRC:  Custom split rule not registered:  %10d", RF_splitCustomIdx + 1);
          RFprintf("\nRF-SRC:  Please register the rule and recompile the package.");
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
    }
    else {
      if (RF_rTargetFactorCount > 0) {
        if (customFunctionArray[CLAS_FAM][RF_splitCustomIdx] == NULL) {
          RFprintf("\nRF-SRC:  *** ERROR *** ");
          RFprintf("\nRF-SRC:  Custom split rule not registered:  %10d", RF_splitCustomIdx + 1);
          RFprintf("\nRF-SRC:  Please register the rule and recompile the package.");
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
      if (RF_rTargetNonFactorCount > 0) {
        if (customFunctionArray[REGR_FAM][RF_splitCustomIdx] == NULL) {
          RFprintf("\nRF-SRC:  *** ERROR *** ");
          RFprintf("\nRF-SRC:  Custom split rule not registered:  %10d", RF_splitCustomIdx + 1);
          RFprintf("\nRF-SRC:  Please register the rule and recompile the package.");
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
    }
  }
}
#include     "splitCustom.h"
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
                  char    multImpFlag) {
  char  result;
  result = FALSE;  
  switch(RF_splitRule) {
  case SURV_LGRNK:
    genericSplit = &logRankNCR;
    break;
  case SURV_LRSCR:
    genericSplit = &logRankNCR;
    break;
  case SURV_L2IMP:
    genericSplit = &l2Impute;
    break;
  case SURV_CR_LAU:
    genericSplit = &logRankCR;
    break;
  case SURV_CR_LOG:
    genericSplit = &logRankCR;
    break;
  case RAND_SPLIT:
    genericSplit = &randomSplit;
    break;
  case REGR_WT_NRM:
    genericSplit = &regressionXwghtSplit;
    break;
  case REGR_WT_OFF:
    genericSplit = &regressionXwghtSplit;
    break;
  case REGR_WT_HVY:
    genericSplit = &regressionXwghtSplit;
    break;
  case CLAS_WT_NRM:
    genericSplit = &classificationXwghtSplit;
    break;
  case CLAS_WT_OFF:
    genericSplit = &classificationXwghtSplit;
    break;
  case CLAS_WT_HVY:
    genericSplit = &classificationXwghtSplit;
    break;
  case MVRG_SPLIT:
    genericSplit = &multivariateSplit;
    break;
  case MVCL_SPLIT:
    genericSplit = &multivariateSplit;
    break;
  case USPV_SPLIT:
    genericSplit = &unsupervisedSplit;
    break;
  case CUST_SPLIT:
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      if (!(RF_opt & OPT_COMP_RISK)) {
        genericSplit = &customSurvivalSplit;
      }
      else {
        genericSplit = &customCompetingRiskSplit;
      }
    }
    else {
      genericSplit = &customMultivariateSplit;
    }
    break;
  default:
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Invalid split rule:  %10d", RF_splitRule);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
    break;
  }
  result = genericSplit(treeID,
                        parent,
                        repMembrIndx,
                        repMembrSize,
                        allMembrIndx,
                        allMembrSize,
                        splitParameterMax,
                        splitValueMaxCont,
                        splitValueMaxFactSize,
                        splitValueMaxFactPtr,
                        splitStatistic,
                        splitIndicator,
                        splitMIA,
                        multImpFlag);
  return result;
}
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
                 char    multImpFlag) {
  uint    *randomCovariateIndex;
  uint    uniformCovariateIndex;
  uint    uniformCovariateSize;
  double *cdf;
  uint    cdfSize;
  uint   *cdfSort;
  uint   *density;
  uint    densitySize;
  uint  **densitySwap;
  uint     covariate;
  double  *splitVector;
  uint     splitVectorSize;
  uint nonMissMembrSize, nonMissMembrSizeStatic;
  uint *nonMissMembrIndx, *nonMissMembrIndxStatic;
  uint   *indxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSize;
  char *localSplitIndicator;
  double deltaMax;
  uint splitLength;
  void *splitVectorPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char preliminaryResult, result;
  char multVarFlag;
  uint j;
  localSplitIndicator    = NULL;  
  mwcpSizeAbsolute       = 0;     
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  multVarFlag = TRUE;
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    multVarFlag = FALSE;
  }
  else {
    if (((RF_rFactorCount == 0) && (RF_rNonFactorCount == 1)) ||
        ((RF_rFactorCount == 1) && (RF_rNonFactorCount == 0))) {
      multVarFlag = FALSE;
    }
  }
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        repMembrSize,
                                        repMembrIndx,
                                        & nonMissMembrSizeStatic,
                                        & nonMissMembrIndxStatic,
                                        & splitVector,
                                        & parent -> mean,
                                        multImpFlag,
                                        multVarFlag);
  if(preliminaryResult) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    stackRandomCovariates(treeID,
                          parent,
                          repMembrSize,
                          multImpFlag,
                          & randomCovariateIndex,
                          & uniformCovariateSize,
                          & cdf,
                          & cdfSize,
                          & cdfSort,
                          & density,
                          & densitySize,
                          & densitySwap);
    uint actualCovariateCount = 0;
    uint candidateCovariateCount = 0;
    while ( ((*splitParameterMax) == 0) &&
            selectRandomCovariates(treeID,
                                   parent,
                                   repMembrIndx,
                                   repMembrSize,
                                   randomCovariateIndex,
                                   & uniformCovariateSize,
                                   & uniformCovariateIndex,
                                   cdf,
                                   & cdfSize,
                                   cdfSort,
                                   density,
                                   & densitySize,
                                   densitySwap,
                                   & covariate,
                                   & actualCovariateCount,
                                   & candidateCovariateCount,
                                   splitVector,
                                   & splitVectorSize,
                                   & indxx,
                                   nonMissMembrSizeStatic,
                                   nonMissMembrIndxStatic,
                                   & nonMissMembrSize,
                                   & nonMissMembrIndx,
                                   multImpFlag)) {
      for (j = 1; j <= repMembrSize; j++) {
        localSplitIndicator[j] = NEITHER;
      }
      leftSize = 0;
      priorMembrIter = 0;
      splitLength = stackAndConstructSplitVector(treeID,
                                                 repMembrSize,
                                                 covariate,
                                                 splitVector,
                                                 splitVectorSize,
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & splitVectorPtr);
      if (factorFlag == FALSE) {
        for (j = 1; j <= nonMissMembrSize; j++) {
          localSplitIndicator[ nonMissMembrIndx[indxx[j]] ] = RIGHT;
        }
      }
      for (j = 1; j < splitLength; j++) {
        if (factorFlag == TRUE) {
          priorMembrIter = 0;
          leftSize = 0;
        }
        virtuallySplitNode(treeID,
                           factorFlag,
                           mwcpSizeAbsolute,
                           covariate,
                           repMembrIndx,
                           repMembrSize,
                           nonMissMembrIndx,
                           nonMissMembrSize,
                           indxx,
                           splitVectorPtr,
                           j,
                           localSplitIndicator,
                           & leftSize,
                           priorMembrIter,
                           & currentMembrIter);
        updateMaximumSplit(treeID,
                           parent,
                           0,  
                           candidateCovariateCount,
                           covariate,
                           j,
                           factorFlag,
                           mwcpSizeAbsolute,
                           repMembrSize,
                           localSplitIndicator,
                           & deltaMax,
                           splitParameterMax,
                           splitValueMaxCont,
                           splitValueMaxFactSize,
                           splitValueMaxFactPtr,
                           splitVectorPtr,
                           splitIndicator);
        j = splitLength;
      }  
      unstackSplitVector(treeID,
                         splitVectorSize,
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         mwcpSizeAbsolute,
                         splitVectorPtr);
      unselectRandomCovariates(treeID,
                               parent,
                               repMembrSize,
                               indxx,
                               nonMissMembrSizeStatic,
                               nonMissMembrIndx,
                               multImpFlag);
    }  
    unstackRandomCovariates(treeID,
                            parent,
                            randomCovariateIndex,
                            uniformCovariateSize,
                            cdf,
                            cdfSize,
                            cdfSort,
                            density,
                            densitySize,
                            densitySwap,
                            repMembrSize);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
  }  
  unstackPreSplit(preliminaryResult,
                  multImpFlag,
                  multVarFlag,
                  repMembrSize,
                  splitVector,
                  nonMissMembrIndxStatic);
  result = summarizeSplitResult(*splitParameterMax,
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
void registerThis (customFunction func, unsigned int family, unsigned int slot) {
  if ((slot >= 1) && (slot <= 16)) {
    customFunctionArray[family][slot-1] = func;
  }
  else {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Invalid slot for custom split rule:  %10d", slot);
    RFprintf("\nRF-SRC:  The slot must be an integer within [1, 16].");
    error("\nRF-SRC:  The application will now exit.\n");
  }    
}
char classificationXwghtSplit (uint    treeID,
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
                               char    multImpFlag) {
  uint   *randomCovariateIndex;
  uint    uniformCovariateIndex;
  uint    uniformCovariateSize;
  double *cdf;
  uint    cdfSize;
  uint   *cdfSort;
  uint   *density;
  uint    densitySize;
  uint  **densitySwap;
  uint     covariate;
  double  *splitVector;
  uint     splitVectorSize;
  uint nonMissMembrSize, nonMissMembrSizeStatic;
  uint *nonMissMembrIndx, *nonMissMembrIndxStatic;
  uint   *indxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  uint splitLength;
  void *splitVectorPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char preliminaryResult, result;
  double delta, deltaMax;
  uint j, k, p;
  localSplitIndicator    = NULL;  
  mwcpSizeAbsolute       = 0;     
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        repMembrSize,
                                        repMembrIndx,
                                        & nonMissMembrSizeStatic,
                                        & nonMissMembrIndxStatic,
                                        & splitVector,
                                        & parent -> mean,
                                        multImpFlag,
                                        FALSE);
  if (preliminaryResult) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    stackRandomCovariates(treeID,
                          parent,
                          repMembrSize,
                          multImpFlag,
                          & randomCovariateIndex,
                          & uniformCovariateSize,
                          & cdf,
                          & cdfSize,
                          & cdfSort,
                          & density,
                          & densitySize,
                          & densitySwap);
    uint responseClassCount = RF_classLevelSize[1];
    uint *parentClassProp = uivector(1, responseClassCount);
    uint *leftClassProp   = uivector(1, responseClassCount);
    uint *rghtClassProp   = uivector(1, responseClassCount);
    double sumLeft, sumRght, sumLeftSqr, sumRghtSqr;
    delta = 0;  
    uint actualCovariateCount = 0;
    uint candidateCovariateCount = 0;
    while (selectRandomCovariates(treeID,
                                  parent,
                                  repMembrIndx,
                                  repMembrSize,
                                  randomCovariateIndex,
                                  & uniformCovariateSize,
                                  & uniformCovariateIndex,
                                  cdf,
                                  & cdfSize,
                                  cdfSort,
                                  density,
                                  & densitySize,
                                  densitySwap,
                                  & covariate,
                                  & actualCovariateCount,
                                  & candidateCovariateCount,
                                  splitVector,
                                  & splitVectorSize,
                                  & indxx,
                                  nonMissMembrSizeStatic,
                                  nonMissMembrIndxStatic,
                                  & nonMissMembrSize,
                                  & nonMissMembrIndx,
                                  multImpFlag)) {
      for (j = 1; j <= repMembrSize; j++) {
        localSplitIndicator[j] = NEITHER;
      }
      for (p=1; p <= responseClassCount; p++) {
        parentClassProp[p] = 0;
      }
      for (j = 1; j <= nonMissMembrSize; j++) {
        parentClassProp[RF_classLevelIndex[1][ (uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[j]]] ]]] ++;
      }
      leftSize = 0;
      priorMembrIter = 0;
      splitLength = stackAndConstructSplitVector(treeID,
                                                 repMembrSize,
                                                 covariate,
                                                 splitVector,
                                                 splitVectorSize,
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & splitVectorPtr);
      if (factorFlag == FALSE) {
        for (j = 1; j <= nonMissMembrSize; j++) {
          localSplitIndicator[ nonMissMembrIndx[indxx[j]] ] = RIGHT;
        }
        for (p = 1; p <= responseClassCount; p++) {
          rghtClassProp[p] = parentClassProp[p];
          leftClassProp[p] = 0;
        }
      }
      for (j = 1; j < splitLength; j++) {
        if (factorFlag == TRUE) {
          priorMembrIter = 0;
          leftSize = 0;
        }
        virtuallySplitNode(treeID,
                              factorFlag,
                              mwcpSizeAbsolute,
                              covariate,
                              repMembrIndx,
                              repMembrSize,
                              nonMissMembrIndx,
                              nonMissMembrSize,
                              indxx,
                              splitVectorPtr,
                              j,
                              localSplitIndicator,
                              & leftSize,
                              priorMembrIter,
                              & currentMembrIter);
        rghtSize = nonMissMembrSize - leftSize;
        if (factorFlag == TRUE) {
          for (p=1; p <= responseClassCount; p++) {
            leftClassProp[p] = 0;
          }
          for (k = 1; k <= nonMissMembrSize; k++) {
            if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT)  {
              leftClassProp[RF_classLevelIndex[1][ (uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] ++;
            }
          }
          for (p=1; p <= responseClassCount; p++) {
            rghtClassProp[p] = parentClassProp[p] - leftClassProp[p];
          }
        }
        else {
          for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
            leftClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] ++;
            rghtClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] --;
          }
        }
        sumLeft = sumRght = 0.0;
        switch(RF_splitRule) {
        case CLAS_WT_NRM:
          for (p=1; p <= responseClassCount; p++) {
            sumLeft += (double) upower(leftClassProp[p], 2);
            sumRght += (double) upower(rghtClassProp[p], 2);
          }
          sumLeftSqr = sumLeft / leftSize;
          sumRghtSqr  = sumRght / rghtSize;
          delta = (sumLeftSqr + sumRghtSqr) / nonMissMembrSize;
          break;
        case CLAS_WT_OFF:
          for (p=1; p <= responseClassCount; p++) {
            sumLeft += pow((double) leftClassProp[p] / (double) leftSize, 2.0);
            sumRght += pow((double) (rghtClassProp[p]) / (double) rghtSize, 2.0);
          }
          delta = sumLeft + sumRght;
          break;
        case CLAS_WT_HVY:
          for (p=1; p <= responseClassCount; p++) {
            sumLeft += (double) upower(leftClassProp[p], 2);
            sumRght += (double) upower(rghtClassProp[p], 2);
          }
          delta =
            (sumLeft / (double) (upower(nonMissMembrSize, 2))) +
            (sumRght / (double) (upower(nonMissMembrSize, 2))) -
            pow((double) leftSize / nonMissMembrSize, 2.0) -
            pow((double) rghtSize / nonMissMembrSize, 2.0) + 2.0;
          break;
        default:
          break;
        }
        updateMaximumSplit(treeID,
                           parent,
                           delta,
                           candidateCovariateCount,
                           covariate,
                           j,
                           factorFlag,
                           mwcpSizeAbsolute,
                           repMembrSize,
                           localSplitIndicator,
                           & deltaMax,
                           splitParameterMax,
                           splitValueMaxCont,
                           splitValueMaxFactSize,
                           splitValueMaxFactPtr,
                           splitVectorPtr,
                           splitIndicator);
        if (factorFlag == FALSE) {
          priorMembrIter = currentMembrIter - 1;
        }
      }  
      unstackSplitVector(treeID,
                         splitVectorSize,
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         mwcpSizeAbsolute,
                         splitVectorPtr);
      unselectRandomCovariates(treeID,
                               parent,
                               repMembrSize,
                               indxx,
                               nonMissMembrSizeStatic,
                               nonMissMembrIndx,
                               multImpFlag);
    }  
    unstackRandomCovariates(treeID,
                            parent,
                            randomCovariateIndex,
                            uniformCovariateSize,
                            cdf,
                            cdfSize,
                            cdfSort,
                            density,
                            densitySize,
                            densitySwap,
                            repMembrSize);
    free_uivector (parentClassProp, 1, responseClassCount);
    free_uivector (leftClassProp,   1, responseClassCount);
    free_uivector (rghtClassProp,   1, responseClassCount);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
  }  
  unstackPreSplit(preliminaryResult,
                  multImpFlag,
                  FALSE,  
                  repMembrSize,
                  splitVector,
                  nonMissMembrIndxStatic);
  result = summarizeSplitResult(*splitParameterMax,
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
char regressionXwghtSplit (uint    treeID,
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
                           char    multImpFlag) {
  uint   *randomCovariateIndex;
  uint    uniformCovariateIndex;
  uint    uniformCovariateSize;
  double *cdf;
  uint    cdfSize;
  uint   *cdfSort;
  uint   *density;
  uint    densitySize;
  uint  **densitySwap;
  uint     covariate;
  double  *splitVector;
  uint     splitVectorSize;
  uint nonMissMembrSize, nonMissMembrSizeStatic;
  uint *nonMissMembrIndx, *nonMissMembrIndxStatic;
  uint  missMembrSize;
  uint *missMembrIndx;
  uint  miaVectorLength;
  char *miaVectorType;
  uint   *indxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSize, rghtSize;
  uint leftSizeAdj, rghtSizeAdj;
  char *localSplitIndicator;
  uint splitLength;
  void *splitVectorPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char preliminaryResult, result;
  double delta, deltaMax;
  uint j, k;
  localSplitIndicator    = NULL;  
  mwcpSizeAbsolute       = 0;     
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  preliminaryResult = getPreSplitResultNew(treeID,
                                           parent,
                                           repMembrSize,
                                           repMembrIndx,
                                           & nonMissMembrSizeStatic,
                                           & nonMissMembrIndxStatic,
                                           & miaVectorLength,
                                           & miaVectorType,
                                           & splitVector,
                                           & parent -> mean,
                                           multImpFlag,
                                           FALSE);
  if (preliminaryResult) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    stackRandomCovariates(treeID,
                          parent,
                          repMembrSize,
                          multImpFlag,
                          & randomCovariateIndex,
                          & uniformCovariateSize,
                          & cdf,
                          & cdfSize,
                          & cdfSort,
                          & density,
                          & densitySize,
                          & densitySwap);
    double sumLeft, sumRght, sumLeftSqr, sumRghtSqr;
    double sumLeftAdj, sumRghtAdj, sumLeftSqrAdj, sumRghtSqrAdj;
    double sumMIA, sumSqrMIA;
    double sumLeftTemp, sumRghtTemp, sumLeftSqrTemp, sumRghtSqrTemp;
    sumLeft = sumRght = sumLeftSqr = sumRghtSqr = 0;                  
    sumLeftAdj = sumRghtAdj = sumLeftSqrAdj = sumRghtSqrAdj = 0;      
    sumMIA = sumSqrMIA = 0;                                           
    sumLeftTemp = sumRghtTemp = sumLeftSqrTemp = sumRghtSqrTemp = 0;  
    delta = 0;                                                        
    uint actualCovariateCount = 0;
    uint candidateCovariateCount = 0;
    while (selectRandomCovariatesNew(treeID,
                                  parent,
                                  repMembrIndx,
                                  repMembrSize,
                                  randomCovariateIndex,
                                  & uniformCovariateSize,
                                  & uniformCovariateIndex,
                                  cdf,
                                  & cdfSize,
                                  cdfSort,
                                  density,
                                  & densitySize,
                                  densitySwap,
                                  & covariate,
                                  & actualCovariateCount,
                                  & candidateCovariateCount,
                                  splitVector,
                                  & splitVectorSize,
                                  & indxx,
                                  nonMissMembrSizeStatic,
                                  nonMissMembrIndxStatic,
                                  & nonMissMembrSize,
                                  & nonMissMembrIndx,
                                  & missMembrSize,
                                  & missMembrIndx,
                                  multImpFlag)) {
      splitLength = stackAndConstructSplitVector(treeID,
                                                 repMembrSize,
                                                 covariate,
                                                 splitVector,
                                                 splitVectorSize,
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & splitVectorPtr);
      leftSize = rghtSize = 0;
      priorMembrIter = 0;
      for (j = 1; j <= repMembrSize; j++) {
        localSplitIndicator[j] = NEITHER;
      }
      if (!multImpFlag) {
        if (RF_optHigh & (OPT_MISS_MIA | OPT_MISS_MIAH)) {
          sumMIA = sumSqrMIA = 0.0;
          switch(RF_splitRule) {
          case REGR_WT_NRM:
            sumMIA = 0.0;
            for (j = 1; j <= missMembrSize; j++) {
              sumMIA += RF_response[treeID][1][ repMembrIndx[missMembrIndx[j]] ];
            }
            break;
          default:
            sumMIA = sumSqrMIA = 0.0;
            for (j = 1; j <= missMembrSize; j++) {
              sumMIA += RF_response[treeID][1][ repMembrIndx[missMembrIndx[j]] ];
              sumSqrMIA += pow(RF_response[treeID][1][ repMembrIndx[missMembrIndx[j]] ], 2.0);
            }
            break;
          }
        }
      }
      if (factorFlag == FALSE) {
        sumLeft = sumLeftSqr = sumRght = sumRghtSqr = 0.0;
        switch(RF_splitRule) {
        case REGR_WT_NRM:
          sumRght = 0.0;
          for (j = 1; j <= nonMissMembrSize; j++) {
            sumRght += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[j]] ];
          }
          break;
        default:
          sumRght = sumRghtSqr = 0.0;
          for (j = 1; j <= nonMissMembrSize; j++) {
            sumRght += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[j]] ];
            sumRghtSqr += pow(RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[j]] ], 2.0);
          }
          break;
        }
        for (j = 1; j <= nonMissMembrSize; j++) {
          localSplitIndicator[ nonMissMembrIndx[j] ] = RIGHT;
        }
      }
      for (j = 1; j < splitLength; j++) {
        if (factorFlag == TRUE) {
          priorMembrIter = 0;
          leftSize = 0;
        }
        virtuallySplitNode(treeID,
                           factorFlag,
                           mwcpSizeAbsolute,
                           covariate,
                           repMembrIndx,
                           repMembrSize,
                           nonMissMembrIndx,
                           nonMissMembrSize,
                           indxx,
                           splitVectorPtr,
                           j,
                           localSplitIndicator,
                           & leftSize,
                           priorMembrIter,
                           & currentMembrIter);
        for (uint m = 1; m <= miaVectorLength; m++) {
          if (factorFlag == TRUE) {
            switch (RF_splitRule) {
            case REGR_WT_NRM:
              sumLeft = sumRght = 0.0;
              for (k = 1; k <= nonMissMembrSize; k++) {
                if (localSplitIndicator[ nonMissMembrIndx[k] ] == LEFT) {
                  sumLeft += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ];
                }
                else {
                  sumRght += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ];
                }
              } 
              break;
            default:
              sumLeft = sumRght = 0.0;
              sumLeftSqr = sumRghtSqr = 0.0;
              for (k = 1; k <= nonMissMembrSize; k++) {
                if (localSplitIndicator[ nonMissMembrIndx[k] ] == LEFT) {
                  sumLeft    += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ];
                  sumLeftSqr += pow(RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ], 2.0);
                }
                else {
                  sumRght    += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ];
                  sumRghtSqr += pow(RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[k]] ], 2.0);
                }
              } 
              break;
            }
          }
          else {
            switch(RF_splitRule) {
            case REGR_WT_NRM:
              for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                sumLeft += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                sumRght -= RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
              }
              break;
            default:
              for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                sumLeft    += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                sumLeftSqr += pow(RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ], 2.0);
                sumRght    -= RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                sumRghtSqr -= pow(RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ], 2.0);
              }
              break;
            }
          }
          switch (miaVectorType[m]) {
          case SPLIT_MIA_LEFT:
            leftSizeAdj = leftSize + missMembrSize;
            rghtSizeAdj = nonMissMembrSize - leftSize;          
            sumLeftAdj = sumLeft + sumMIA;
            sumLeftSqrAdj = sumLeftSqr + sumSqrMIA;
            sumRghtAdj = sumRght;
            sumRghtSqrAdj = sumRghtSqr;
            break;
          case SPLIT_MIA_RGHT:
            leftSizeAdj = leftSize;
            rghtSizeAdj = nonMissMembrSize - leftSize + missMembrSize;          
            sumLeftAdj = sumLeft;
            sumLeftSqrAdj = sumLeftSqr;
            sumRghtAdj = sumRght + sumMIA;
            sumRghtSqrAdj = sumRghtSqr + sumSqrMIA;
            break;
          default:
            leftSizeAdj = leftSize;
            rghtSizeAdj = nonMissMembrSize - leftSize;
            sumLeftAdj = sumLeft;
            sumRghtAdj = sumRght;
            sumLeftSqrAdj = sumLeftSqr;
            sumRghtSqrAdj = sumRghtSqr;
            break;
          }
          switch(RF_splitRule) {
          case REGR_WT_NRM:
            sumLeftSqrAdj = pow(sumLeftAdj, 2.0) / leftSizeAdj;
            sumRghtSqrAdj = pow(sumRghtAdj, 2.0) / rghtSizeAdj;
            delta = sumLeftSqrAdj + sumRghtSqrAdj;
            break;
          case REGR_WT_OFF:
            sumLeftTemp = pow(sumLeftAdj, 2.0) / pow(leftSizeAdj, 2.0);
            sumRghtTemp = pow(sumRghtAdj, 2.0) / pow(rghtSizeAdj, 2.0);
            sumLeftSqrTemp = sumLeftSqrAdj / leftSizeAdj;
            sumRghtSqrTemp = sumRghtSqrAdj / rghtSizeAdj;
            delta = sumLeftTemp + sumRghtTemp - sumLeftSqrTemp - sumRghtSqrTemp;
            break;
          case REGR_WT_HVY:
            sumLeftTemp = pow(sumLeftAdj, 2.0) / pow (leftSizeAdj + rghtSizeAdj, 2.0);
            sumRghtTemp = pow(sumRghtAdj, 2.0) / pow (leftSizeAdj + rghtSizeAdj, 2.0);
            sumLeftSqrTemp = sumLeftSqrAdj * leftSizeAdj / pow (leftSizeAdj + rghtSizeAdj, 2.0);
            sumRghtSqrTemp = sumRghtSqrAdj * rghtSizeAdj / pow (leftSizeAdj + rghtSizeAdj, 2.0);
            delta = sumLeftTemp + sumRghtTemp - sumLeftSqrTemp - sumRghtSqrTemp;
            break;
          default:
            break;
          }
          updateMaximumSplitNew(treeID,
                                parent,
                                delta,
                                candidateCovariateCount,
                                covariate,
                                j,
                                factorFlag,
                                mwcpSizeAbsolute,
                                repMembrSize,
                                localSplitIndicator,
                                miaVectorType[m],
                                & deltaMax,
                                splitParameterMax,
                                splitValueMaxCont,
                                splitValueMaxFactSize,
                                splitValueMaxFactPtr,
                                splitVectorPtr,
                                splitIndicator,
                                splitMIA);
        }  
        if (factorFlag == FALSE) {
          priorMembrIter = currentMembrIter - 1;
        }
      }  
      if (!multImpFlag) {
        if (RF_optHigh & (OPT_MISS_MIA | OPT_MISS_MIAH)) {
          if (missMembrSize > 0) {
            switch(RF_splitRule) {
            case REGR_WT_NRM:
              sumLeft = 0.0;
              for (j = 1; j <= nonMissMembrSize; j++) {
                sumLeft += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[j]] ];
              }
              leftSizeAdj = nonMissMembrSize;
              rghtSizeAdj = missMembrSize;
              sumLeftAdj = sumLeft;
              sumRghtAdj = sumMIA;
              sumLeftSqrAdj = pow(sumLeftAdj, 2.0) / leftSizeAdj;
              sumRghtSqrAdj = pow(sumRghtAdj, 2.0) / rghtSizeAdj;
              delta = sumLeftSqrAdj + sumRghtSqrAdj;
              break;
            case REGR_WT_OFF:
              sumLeft = sumLeftSqr = 0.0;
              for (j = 1; j <= nonMissMembrSize; j++) {
                sumLeft += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[j]] ];
                sumLeftSqr += pow(RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[j]] ], 2.0);
              }
              leftSizeAdj = nonMissMembrSize;
              rghtSizeAdj = missMembrSize;
              sumLeftAdj    = sumLeft;
              sumLeftSqrAdj = sumLeftSqr;
              sumRghtAdj    = sumMIA;
              sumRghtSqrAdj = sumSqrMIA;
              sumLeftTemp = pow(sumLeftAdj, 2.0) / pow(leftSizeAdj, 2.0);
              sumRghtTemp = pow(sumRghtAdj, 2.0) / pow(rghtSizeAdj, 2.0);
              sumLeftSqrTemp = sumLeftSqrAdj / leftSizeAdj;
              sumRghtSqrTemp = sumRghtSqrAdj / rghtSizeAdj;
              delta = sumLeftTemp + sumRghtTemp - sumLeftSqrTemp - sumRghtSqrTemp;
              break;
            case REGR_WT_HVY:
              sumLeft = sumLeftSqr = 0.0;
              for (j = 1; j <= nonMissMembrSize; j++) {
                sumLeft += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[j]] ];
                sumLeftSqr += pow(RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[j]] ], 2.0);
              }
              leftSizeAdj = nonMissMembrSize;
              rghtSizeAdj = missMembrSize;
              sumLeftAdj    = sumLeft;
              sumLeftSqrAdj = sumLeftSqr;
              sumRghtAdj    = sumMIA;
              sumRghtSqrAdj = sumSqrMIA;
              sumLeftTemp = pow(sumLeftAdj, 2.0) / pow (leftSizeAdj + rghtSizeAdj, 2.0);
              sumRghtTemp = pow(sumRghtAdj, 2.0) / pow (leftSizeAdj + rghtSizeAdj, 2.0);
              sumLeftSqrTemp = sumLeftSqrAdj * leftSizeAdj / pow (leftSizeAdj + rghtSizeAdj, 2.0);
              sumRghtSqrTemp = sumRghtSqrAdj * rghtSizeAdj / pow (leftSizeAdj + rghtSizeAdj, 2.0);
              delta = sumLeftTemp + sumRghtTemp - sumLeftSqrTemp - sumRghtSqrTemp;
              break;
            }
            for (j = 1; j <= nonMissMembrSize; j++) {
              localSplitIndicator[ nonMissMembrIndx[j] ] = LEFT;
            }
            for (j = 1; j <= missMembrSize; j++) {
              localSplitIndicator[ missMembrIndx[j] ] = RIGHT;
            }
            updateMaximumSplitNew(treeID,
                                  parent,
                                  delta,
                                  candidateCovariateCount,
                                  covariate,
                                  0,
                                  factorFlag,
                                  mwcpSizeAbsolute,
                                  repMembrSize,
                                  localSplitIndicator,
                                  SPLIT_MIA_UNIT,
                                  & deltaMax,
                                  splitParameterMax,
                                  splitValueMaxCont,
                                  splitValueMaxFactSize,
                                  splitValueMaxFactPtr,
                                  splitVectorPtr,
                                  splitIndicator,
                                  splitMIA);
          }
        }
      }
      unstackSplitVector(treeID,
                         splitVectorSize,
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         mwcpSizeAbsolute,
                         splitVectorPtr);
      unselectRandomCovariatesNew(treeID,
                               parent,
                               repMembrSize,
                               indxx,
                               nonMissMembrSizeStatic,
                               nonMissMembrIndx,
                               missMembrIndx,
                               multImpFlag);
    }  
    unstackRandomCovariates(treeID,
                            parent, 
                            randomCovariateIndex,
                            uniformCovariateSize,
                            cdf,
                            cdfSize,
                            cdfSort,
                            density,
                            densitySize,
                            densitySwap,
                            repMembrSize);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
  }  
  unstackPreSplitNew(preliminaryResult,
                     multImpFlag,
                     FALSE,  
                     repMembrSize,
                     splitVector,
                     nonMissMembrIndxStatic,
                     miaVectorType);
  result = summarizeSplitResult(*splitParameterMax,
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
char logRankNCR (uint    treeID,
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
                 char    multImpFlag) {
  uint   *randomCovariateIndex;
  uint    uniformCovariateIndex;
  uint    uniformCovariateSize;
  double *cdf;
  uint    cdfSize;
  uint   *cdfSort;
  uint   *density;
  uint    densitySize;
  uint  **densitySwap;
  uint     covariate;
  double  *splitVector;
  uint     splitVectorSize;
  uint nonMissMembrSize, nonMissMembrSizeStatic;
  uint *nonMissMembrIndx, *nonMissMembrIndxStatic;
  uint   *indxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSize;
  char *localSplitIndicator;
  uint splitLength;
  void *splitVectorPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char preliminaryResult, result;
  double delta, deltaMax;
  uint j, k, m;
  localSplitIndicator    = NULL;  
  mwcpSizeAbsolute       = 0;     
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        repMembrSize,
                                        repMembrIndx,
                                        & nonMissMembrSizeStatic,
                                        & nonMissMembrIndxStatic,
                                        & splitVector,
                                        & parent -> mean,
                                        multImpFlag,
                                        FALSE);
  if (preliminaryResult) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    stackRandomCovariates(treeID,
                          parent,
                          repMembrSize,
                          multImpFlag,
                          & randomCovariateIndex,
                          & uniformCovariateSize,
                          & cdf,
                          & cdfSize,
                          & cdfSort,
                          & density,
                          & densitySize,
                          & densitySwap);
    uint *localEventTimeCount, *localEventTimeIndex;
    uint  localEventTimeSize;
    uint *nodeParentEvent,  *nodeLeftEvent,  *nodeRightEvent;
    uint *nodeParentAtRisk, *nodeLeftAtRisk, *nodeRightAtRisk;
    uint   *survivalTimeIndexRank;
    double *survivalRank;
    double  meanSurvRank, varSurvRank;
    double deltaNum, deltaNumAdj, deltaDen;
    uint   tIndx;
    meanSurvRank = varSurvRank = 0;  
    survivalTimeIndexRank = NULL;  
    survivalRank = NULL;  
    localEventTimeSize = 0;  
    delta = deltaNum = 0;  
    switch(RF_splitRule) {
    case SURV_LGRNK:
      if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
        stackAndGetSplitSurv(treeID,
                             repMembrIndx,
                             repMembrSize,
                             nonMissMembrIndxStatic,
                             nonMissMembrSizeStatic,
                             & localEventTimeCount,
                             & localEventTimeIndex,
                             & localEventTimeSize,
                             & nodeParentEvent,
                             & nodeParentAtRisk,
                             & nodeLeftEvent,
                             & nodeLeftAtRisk,
                             & nodeRightEvent,
                             & nodeRightAtRisk);
      }
      break;
    case SURV_LRSCR:
      survivalTimeIndexRank = uivector(1, repMembrSize);
      survivalRank = dvector(1, repMembrSize);
      localEventTimeSize = 1;
      break;
    default:
      break;
    }
    uint actualCovariateCount = 0;
    uint candidateCovariateCount = 0;
    while (selectRandomCovariates(treeID,
                                  parent,
                                  repMembrIndx,
                                  repMembrSize,
                                  randomCovariateIndex,
                                  & uniformCovariateSize,
                                  & uniformCovariateIndex,
                                  cdf,
                                  & cdfSize,
                                  cdfSort,
                                  density,
                                  & densitySize,
                                  densitySwap,
                                  & covariate,
                                  & actualCovariateCount,
                                  & candidateCovariateCount,
                                  splitVector,
                                  & splitVectorSize,
                                  & indxx,
                                  nonMissMembrSizeStatic,
                                  nonMissMembrIndxStatic,
                                  & nonMissMembrSize,
                                  & nonMissMembrIndx,
                                  multImpFlag)) {
      splitLength = stackAndConstructSplitVector(treeID,
                                                 repMembrSize,
                                                 covariate,
                                                 splitVector,
                                                 splitVectorSize,
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & splitVectorPtr);
      switch(RF_splitRule) {
      case SURV_LGRNK:
        if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
          stackAndGetSplitSurv(treeID,
                               repMembrIndx,
                               repMembrSize,
                               nonMissMembrIndx,
                               nonMissMembrSize,
                               & localEventTimeCount,
                               & localEventTimeIndex,
                               & localEventTimeSize,
                               & nodeParentEvent,
                               & nodeParentAtRisk,
                               & nodeLeftEvent,
                               & nodeLeftAtRisk,
                               & nodeRightEvent,
                               & nodeRightAtRisk);
        }
        break;
      case SURV_LRSCR:
        localEventTimeSize = 1;
        break;
      default:
        break;
      }
      if (localEventTimeSize > 0) {
        for (j = 1; j <= repMembrSize; j++) {
          localSplitIndicator[j] = NEITHER;
        }
        leftSize = 0;
        priorMembrIter = 0;
        if (factorFlag == FALSE) {
          for (j = 1; j <= nonMissMembrSize; j++) {
            localSplitIndicator[ nonMissMembrIndx[indxx[j]] ] = RIGHT;
          }
          switch(RF_splitRule) {
          case SURV_LGRNK:
            for (m = 1; m <= localEventTimeSize; m++) {
              nodeLeftEvent[m] = nodeLeftAtRisk[m] = 0;
            }
            break;
          case SURV_LRSCR:
            deltaNum =  0.0;
            break;
          default:
            break;
          }
        }
        switch(RF_splitRule) {
        case SURV_LGRNK:
          break;
        case SURV_LRSCR:
          for (k = 1; k <= nonMissMembrSize; k++) {
            survivalTimeIndexRank[k] = 0;
            for (j = 1; j <= nonMissMembrSize; j++) {
              if ( RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[j]] ]  <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] ) {
                survivalTimeIndexRank[k] ++;
              }
            }
          }
          meanSurvRank = varSurvRank = 0;
          for (k = 1; k <= nonMissMembrSize; k++) {
            survivalRank[k] = 0;
            for (j = 1; j <= survivalTimeIndexRank[k]; j++) {
              survivalRank[k] = survivalRank[k] + (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[j]]] ] / (nonMissMembrSize - survivalTimeIndexRank[j] + 1) );
            }
            survivalRank[k] = RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - survivalRank[k];
            meanSurvRank = meanSurvRank + survivalRank[k];
            varSurvRank = varSurvRank +  pow(survivalRank[k], 2.0);
          }
          varSurvRank = ( varSurvRank - (pow(meanSurvRank, 2.0) / nonMissMembrSize) ) / (nonMissMembrSize - 1);
          meanSurvRank = meanSurvRank / nonMissMembrSize;
          break;
        default:
          break;
        }
        for (j = 1; j < splitLength; j++) {
          if (factorFlag == TRUE) {
            priorMembrIter = 0;
            leftSize = 0;
          }
          virtuallySplitNode(treeID,
                             factorFlag,
                             mwcpSizeAbsolute,
                             covariate,
                             repMembrIndx,
                             repMembrSize,
                             nonMissMembrIndx,
                             nonMissMembrSize,
                             indxx,
                             splitVectorPtr,
                             j,
                             localSplitIndicator,
                             & leftSize,
                             priorMembrIter,
                             & currentMembrIter);
          if (factorFlag == TRUE) {
            switch(RF_splitRule) {
            case SURV_LGRNK:
              for (m = 1; m <= localEventTimeSize; m++) {
                nodeLeftEvent[m] = nodeLeftAtRisk[m] = 0;
              }
              for (k = 1; k <= nonMissMembrSize; k++) {
                if (localSplitIndicator[  nonMissMembrIndx[indxx[k]]  ] == LEFT) {
                  tIndx = 0;  
                  for (m = 1; m <= localEventTimeSize; m++) {
                    if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]) {
                      tIndx = m;
                      nodeLeftAtRisk[tIndx] ++;
                    }
                    else {
                      m = localEventTimeSize;
                    }
                  }
                  if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
                    nodeLeftEvent[tIndx] ++;
                  }
                }
                else {
                }
              } 
              break;
            case SURV_LRSCR:
              deltaNum = 0.0;
              for (k = 1; k <= nonMissMembrSize; k++) {
                if (localSplitIndicator[ nonMissMembrIndx[k] ] == LEFT) {
                  deltaNum = deltaNum + survivalRank[k];
                }
              }
              break;
            default:
              break;
            }
          }
          else {
            switch(RF_splitRule) {
            case SURV_LGRNK:
              for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                tIndx = 0;  
                for (m = 1; m <= localEventTimeSize; m++) {
                  if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]) {
                    tIndx = m;
                    nodeLeftAtRisk[tIndx] ++;
                  }
                  else {
                    m = localEventTimeSize;
                  }
                }
                if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
                  nodeLeftEvent[tIndx] ++;
                }
              }
              break;
            case SURV_LRSCR:
              for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                deltaNum = deltaNum + survivalRank[ indxx[k] ];
              }
              break;
            default:
              break;
            }
          }
          switch(RF_splitRule) {
          case SURV_LGRNK:
            delta = deltaNum = deltaDen =  0.0;
            for (k=1; k <= localEventTimeSize; k++) {
              deltaNum = deltaNum + ((double) nodeLeftEvent[k] - ((double) ( nodeLeftAtRisk[k] * nodeParentEvent[k]) / nodeParentAtRisk[k]));
              if (nodeParentAtRisk[k] >= 2) {
                deltaDen = deltaDen + (
                                       ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k]) *
                                       (1.0 - ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k])) *
                                       ((double) (nodeParentAtRisk[k] - nodeParentEvent[k]) / (nodeParentAtRisk[k] - 1)) * nodeParentEvent[k]
                                       );
              }
            }
            deltaNum = fabs(deltaNum);
            deltaDen = sqrt(deltaDen);
            if (deltaDen <= EPSILON) {
              if (deltaNum <= EPSILON) {
                delta = 0.0;
              }
              else {
                delta = deltaNum / deltaDen;
              }
            }
            else {
              delta = deltaNum / deltaDen;
            }
            break;
          case SURV_LRSCR:
            deltaNumAdj  = deltaNum - (leftSize * meanSurvRank);
            deltaDen     = leftSize * (1.0 - (leftSize / nonMissMembrSize)) * varSurvRank;
            deltaNumAdj = fabs(deltaNumAdj);
            deltaDen = sqrt(deltaDen);
            if (deltaDen <= EPSILON) {
              if (deltaNumAdj <= EPSILON) {
                delta = 0.0;
              }
              else {
                delta = deltaNumAdj / deltaDen;
              }
            }
            else {
              delta = deltaNumAdj / deltaDen;
            }
            break;
          default:
            break;
          }
          updateMaximumSplit(treeID,
                             parent,
                             delta,
                             candidateCovariateCount,
                             covariate,
                             j,
                             factorFlag,
                             mwcpSizeAbsolute,
                             repMembrSize,
                             localSplitIndicator,
                             & deltaMax,
                             splitParameterMax,
                             splitValueMaxCont,
                             splitValueMaxFactSize,
                             splitValueMaxFactPtr,
                             splitVectorPtr,
                             splitIndicator);
          if (factorFlag == FALSE) {
            priorMembrIter = currentMembrIter - 1;
          }
        }  
      }  
      else {
      }
      unstackSplitVector(treeID,
                         splitVectorSize,
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         mwcpSizeAbsolute,
                         splitVectorPtr);
      unselectRandomCovariates(treeID,
                               parent,
                               repMembrSize,
                               indxx,
                               nonMissMembrSizeStatic,
                               nonMissMembrIndx,
                               multImpFlag);
      switch(RF_splitRule) {
      case SURV_LGRNK:
        if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
          unstackSplitSurv(localEventTimeCount,
                           localEventTimeIndex,
                           localEventTimeSize,
                           nodeParentEvent,
                           nodeParentAtRisk,
                           nodeLeftEvent,
                           nodeLeftAtRisk,
                           nodeRightEvent,
                           nodeRightAtRisk);
        }
        break;
      case SURV_LRSCR:
        break;
      default:
        break;
      }
    }  
    switch(RF_splitRule) {
    case SURV_LGRNK:
      if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
        unstackSplitSurv(localEventTimeCount,
                         localEventTimeIndex,
                         localEventTimeSize,
                         nodeParentEvent,
                         nodeParentAtRisk,
                         nodeLeftEvent,
                         nodeLeftAtRisk,
                         nodeRightEvent,
                         nodeRightAtRisk);
      }
    break;
    case SURV_LRSCR:
      free_uivector(survivalTimeIndexRank, 1, repMembrSize);
      free_dvector(survivalRank, 1, repMembrSize);
      break;
    default:
      break;
    }
    unstackRandomCovariates(treeID,
                            parent,
                            randomCovariateIndex,
                            uniformCovariateSize,
                            cdf,
                            cdfSize,
                            cdfSort,
                            density,
                            densitySize,
                            densitySwap,
                            repMembrSize);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
  }  
  unstackPreSplit(preliminaryResult,
                  multImpFlag,
                  FALSE, 
                  repMembrSize,
                  splitVector,
                  nonMissMembrIndxStatic);
  result = summarizeSplitResult(*splitParameterMax,
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
char logRankCR (uint    treeID,
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
                char    multImpFlag) {
  uint   *randomCovariateIndex;
  uint    uniformCovariateIndex;
  uint    uniformCovariateSize;
  double *cdf;
  uint    cdfSize;
  uint   *cdfSort;
  uint   *density;
  uint    densitySize;
  uint  **densitySwap;
  uint     covariate;
  double  *splitVector;
  uint     splitVectorSize;
  uint nonMissMembrSize, nonMissMembrSizeStatic;
  uint *nonMissMembrIndx, *nonMissMembrIndxStatic;
  uint   *indxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSize;
  char *localSplitIndicator;
  uint splitLength;
  void *splitVectorPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char preliminaryResult, result;
  double delta, deltaMax;
  uint j, k, m;
  localSplitIndicator    = NULL;  
  mwcpSizeAbsolute       = 0;     
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        repMembrSize,
                                        repMembrIndx,
                                        & nonMissMembrSizeStatic,
                                        & nonMissMembrIndxStatic,
                                        & splitVector,
                                        & parent -> mean,
                                        multImpFlag,
                                        FALSE);
  if (preliminaryResult) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    stackRandomCovariates(treeID,
                          parent,
                          repMembrSize,
                          multImpFlag,
                          & randomCovariateIndex,
                          & uniformCovariateSize,
                          & cdf,
                          & cdfSize,
                          & cdfSort,
                          & density,
                          & densitySize,
                          & densitySwap);
    uint *localEventTimeCount, *localEventTimeIndex;
    uint  localEventTimeSize;
    uint *nodeParentEvent,  *nodeLeftEvent,  *nodeRightEvent;
    uint *nodeParentAtRisk, *nodeLeftAtRisk, *nodeRightAtRisk;
    uint **nodeParentEventCR, **nodeLeftEventCR;
    uint **nodeParentInclusiveAtRisk, **nodeLeftInclusiveAtRisk;
    nodeParentInclusiveAtRisk = nodeLeftInclusiveAtRisk = NULL;  
    nodeParentEventCR = nodeLeftEventCR = NULL;  
    double deltaNum, deltaSubNum, deltaDen, deltaSubDen;
    uint   tIndx;
    uint   q, s, r;
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      stackAndGetSplitSurv(treeID,
                           repMembrIndx,
                           repMembrSize,
                           nonMissMembrIndxStatic,
                           nonMissMembrSizeStatic,
                           & localEventTimeCount,
                           & localEventTimeIndex,
                           & localEventTimeSize,
                           & nodeParentEvent,
                           & nodeParentAtRisk,
                           & nodeLeftEvent,
                           & nodeLeftAtRisk,
                           & nodeRightEvent,
                           & nodeRightAtRisk);
      nodeParentEventCR = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
      nodeLeftEventCR = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
      switch(RF_splitRule) {
      case SURV_CR_LAU:
        nodeParentInclusiveAtRisk = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
        nodeLeftInclusiveAtRisk = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
        break;
      case SURV_CR_LOG:
        break;
      default:
        break;
      }
    }
    uint actualCovariateCount = 0;
    uint candidateCovariateCount = 0;
    while (selectRandomCovariates(treeID,
                                  parent,
                                  repMembrIndx,
                                  repMembrSize,
                                  randomCovariateIndex,
                                  & uniformCovariateSize,
                                  & uniformCovariateIndex,
                                  cdf,
                                  & cdfSize,
                                  cdfSort,
                                  density,
                                  & densitySize,
                                  densitySwap,
                                  & covariate,
                                  & actualCovariateCount,
                                  & candidateCovariateCount,
                                  splitVector,
                                  & splitVectorSize,
                                  & indxx,
                                  nonMissMembrSizeStatic,
                                  nonMissMembrIndxStatic,
                                  & nonMissMembrSize,
                                  & nonMissMembrIndx,
                                  multImpFlag)) {
      splitLength = stackAndConstructSplitVector(treeID,
                                                 repMembrSize,
                                                 covariate,
                                                 splitVector,
                                                 splitVectorSize,
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & splitVectorPtr);
      if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
        stackAndGetSplitSurv(treeID,
                             repMembrIndx,
                             repMembrSize,
                             nonMissMembrIndx,
                             nonMissMembrSize,
                             & localEventTimeCount,
                             & localEventTimeIndex,
                             & localEventTimeSize,
                             & nodeParentEvent,
                             & nodeParentAtRisk,
                             & nodeLeftEvent,
                             & nodeLeftAtRisk,
                             & nodeRightEvent,
                             & nodeRightAtRisk);
        if (localEventTimeSize > 0) {
          nodeParentEventCR = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
          nodeLeftEventCR = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
        }
        switch(RF_splitRule) {
        case SURV_CR_LAU:
          if (localEventTimeSize > 0) {
            nodeParentInclusiveAtRisk = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
            nodeLeftInclusiveAtRisk = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
          }
          break;
        case SURV_CR_LOG:
          break;
        default:
          break;
        }
      }
      if (localEventTimeSize > 0) {
        for (j = 1; j <= repMembrSize; j++) {
          localSplitIndicator[j] = NEITHER;
        }
        leftSize = 0;
        priorMembrIter = 0;
        for (m=1; m <= localEventTimeSize; m++) {
          for (q = 1; q <= RF_eventTypeSize; q++) {
            nodeParentEventCR[q][m] = 0;
          }
        }
        for (k = 1; k <= nonMissMembrSize; k++) {
          if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
            tIndx = 0;  
            for (m = 1; m <= localEventTimeSize; m++) {
              if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]) {
                tIndx = m;
              }
              else {
                m = localEventTimeSize;
              }
            }
            nodeParentEventCR[RF_eventTypeIndex[(uint) RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]][tIndx] ++;
          }
        }
        switch(RF_splitRule) {
        case SURV_CR_LAU:
          for (m = 1; m <= localEventTimeSize; m++) {
            for (q = 1; q <= RF_eventTypeSize; q++) {
              if (RF_crWeight[q] > 0) {
                nodeParentInclusiveAtRisk[q][m] = nodeParentAtRisk[m];
                for (s = 1; s < m; s++) {
                  for (r = 1; r <= RF_eventTypeSize; r++) {
                    if (q != r) {
                      nodeParentInclusiveAtRisk[q][m]  += nodeParentEventCR[r][s];
                    }
                  }
                }
              }
            }
          }
          break;
        case SURV_CR_LOG:
          break;
        default:
          break;
        }
        if (factorFlag == FALSE) {
          for (j = 1; j <= nonMissMembrSize; j++) {
            localSplitIndicator[ nonMissMembrIndx[indxx[j]] ] = RIGHT;
          }
          for (m = 1; m <= localEventTimeSize; m++) {
            nodeLeftAtRisk[m] = 0;
            for (q = 1; q <= RF_eventTypeSize; q++) {
              nodeLeftEventCR[q][m] = 0;
            }
          }
        }
        for (j = 1; j < splitLength; j++) {
          if (factorFlag == TRUE) {
            priorMembrIter = 0;
            leftSize = 0;
          }
          virtuallySplitNode(treeID,
                             factorFlag,
                             mwcpSizeAbsolute,
                             covariate,
                             repMembrIndx,
                             repMembrSize,
                             nonMissMembrIndx,
                             nonMissMembrSize,
                             indxx,
                             splitVectorPtr,
                             j,
                             localSplitIndicator,
                             & leftSize,
                             priorMembrIter,
                             & currentMembrIter);
          if (factorFlag == TRUE) {
            for (m = 1; m <= localEventTimeSize; m++) {
              nodeLeftAtRisk[m] = 0;
              for (q = 1; q <= RF_eventTypeSize; q++) {
                nodeLeftEventCR[q][m] = 0;
              }
            }
            for (k = 1; k <= nonMissMembrSize; k++) {
              if (localSplitIndicator[  nonMissMembrIndx[indxx[k]]  ] == LEFT) {
                tIndx = 0;  
                for (m = 1; m <= localEventTimeSize; m++) {
                  if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]) {
                    tIndx = m;
                    nodeLeftAtRisk[tIndx] ++;
                  }
                  else {
                    m = localEventTimeSize;
                  }
                }
                if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
                  nodeLeftEventCR[RF_eventTypeIndex[(uint) RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]][tIndx] ++;
                }
              }
            }
          }
          else {
            for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
              tIndx = 0;  
              for (m = 1; m <= localEventTimeSize; m++) {
                if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]) {
                  tIndx = m;
                  nodeLeftAtRisk[tIndx] ++;
                }
                else {
                  m = localEventTimeSize;
                }
              }
              if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
                nodeLeftEventCR[RF_eventTypeIndex[(uint) RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]][tIndx] ++;
              }
            }
          }
          switch(RF_splitRule) {
          case SURV_CR_LAU:
            for (m=1; m <= localEventTimeSize; m++) {
              for (q = 1; q <= RF_eventTypeSize; q++) {
                if (RF_crWeight[q] > 0) {
                  nodeLeftInclusiveAtRisk[q][m]   = nodeLeftAtRisk[m];
                  for (s = 1; s < m; s++) {
                    for (r = 1; r <= RF_eventTypeSize; r++) {
                      if (q != r) {
                        nodeLeftInclusiveAtRisk[q][m]    += nodeLeftEventCR[r][s];
                      }
                    }
                  }
                }
              }
            }
            break;
          case SURV_CR_LOG:
            break;
          default:
            break;
          }
          delta = deltaNum = deltaDen =  0.0;
          switch(RF_splitRule) {
          case SURV_CR_LAU:
            for (q = 1; q <= RF_eventTypeSize; q++) {
              if (RF_crWeight[q] > 0) {
                deltaSubNum = 0;
                for (m = 1; m <= localEventTimeSize; m++) {
                  deltaSubNum = deltaSubNum + (nodeLeftEventCR[q][m] - (nodeParentEventCR[q][m] * ((double) nodeLeftInclusiveAtRisk[q][m] / nodeParentInclusiveAtRisk[q][m])));
                }
                deltaNum = deltaNum + (RF_crWeight[q] * deltaSubNum);
                deltaSubDen = 0;
                for (m = 1; m <= localEventTimeSize; m++) {
                  if (nodeParentAtRisk[m] >= 2) {
                    deltaSubDen = deltaSubDen  + (
                                                  (nodeParentEventCR[q][m] * ((double) nodeLeftInclusiveAtRisk[q][m] / nodeParentInclusiveAtRisk[q][m])) *
                                                  (1.0 - ((double) nodeLeftInclusiveAtRisk[q][m] / nodeParentInclusiveAtRisk[q][m])) *
                                                  ((double) (nodeParentInclusiveAtRisk[q][m] - nodeParentEventCR[q][m]) / (nodeParentInclusiveAtRisk[q][m] - 1))
                                                  );
                  }
                }
                deltaDen = deltaDen + (RF_crWeight[q] * RF_crWeight[q] * deltaSubDen);
              }
            }
            break;
          case SURV_CR_LOG:
            for (q = 1; q <= RF_eventTypeSize; q++) {
              if (RF_crWeight[q] > 0) {
                deltaSubNum = 0;
                for (m=1; m <= localEventTimeSize; m++) {
                  deltaSubNum = deltaSubNum + (nodeLeftEventCR[q][m] - (nodeParentEventCR[q][m] * ((double) nodeLeftAtRisk[m] / nodeParentAtRisk[m])));
                }
                deltaNum = deltaNum + (RF_crWeight[q] * deltaSubNum);
                deltaSubDen = 0;
                for (m = 1; m <= localEventTimeSize; m++) {
                  if (nodeParentAtRisk[m] >= 2) {
                    deltaSubDen = deltaSubDen  + (
                                                  (nodeParentEventCR[q][m] * ((double) nodeLeftAtRisk[m] / nodeParentAtRisk[m])) *
                                                  (1.0 - ((double) nodeLeftAtRisk[m] / nodeParentAtRisk[m])) *
                                                  ((double) (nodeParentAtRisk[m] - nodeParentEventCR[q][m]) / (nodeParentAtRisk[m] - 1))
                                                  );
                  }
                }
                deltaDen = deltaDen + (RF_crWeight[q] * RF_crWeight[q] * deltaSubDen);
              }
            }
            break;
          default:
            break;
          }
          deltaNum = fabs(deltaNum);
          deltaDen = sqrt(deltaDen);
          if (deltaDen <= EPSILON) {
            if (deltaNum <= EPSILON) {
              delta = 0.0;
            }
            else {
              delta = deltaNum / deltaDen;
            }
          }
          else {
            delta = deltaNum / deltaDen;
          }
          updateMaximumSplit(treeID,
                             parent,
                             delta,
                             candidateCovariateCount,
                             covariate,
                             j,
                             factorFlag,
                             mwcpSizeAbsolute,
                             repMembrSize,
                             localSplitIndicator,
                             & deltaMax,
                             splitParameterMax,
                             splitValueMaxCont,
                             splitValueMaxFactSize,
                             splitValueMaxFactPtr,
                             splitVectorPtr,
                             splitIndicator);
          if (factorFlag == FALSE) {
            priorMembrIter = currentMembrIter - 1;
          }
        }  
      }  
      else {
      }
      unstackSplitVector(treeID,
                         splitVectorSize,
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         mwcpSizeAbsolute,
                         splitVectorPtr);
      unselectRandomCovariates(treeID,
                               parent,
                               repMembrSize,
                               indxx,
                               nonMissMembrSizeStatic,
                               nonMissMembrIndx,
                               multImpFlag);
      if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
        unstackSplitSurv(localEventTimeCount,
                         localEventTimeIndex,
                         localEventTimeSize,
                         nodeParentEvent,
                         nodeParentAtRisk,
                         nodeLeftEvent,
                         nodeLeftAtRisk,
                         nodeRightEvent,
                         nodeRightAtRisk);
        if (localEventTimeSize > 0) {
          free_uimatrix(nodeParentEventCR, 1, RF_eventTypeSize, 1, localEventTimeSize);
          free_uimatrix(nodeLeftEventCR, 1, RF_eventTypeSize, 1, localEventTimeSize);
        }
        switch(RF_splitRule) {
        case SURV_CR_LAU:
          if (localEventTimeSize > 0) {
            free_uimatrix(nodeParentInclusiveAtRisk, 1, RF_eventTypeSize, 1, localEventTimeSize);
            free_uimatrix(nodeLeftInclusiveAtRisk, 1, RF_eventTypeSize, 1, localEventTimeSize);
          }
          break;
        case SURV_CR_LOG:
          break;
        default:
          break;
        }
      }
    }  
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
        unstackSplitSurv(localEventTimeCount,
                         localEventTimeIndex,
                         localEventTimeSize,
                         nodeParentEvent,
                         nodeParentAtRisk,
                         nodeLeftEvent,
                         nodeLeftAtRisk,
                         nodeRightEvent,
                         nodeRightAtRisk);
      free_uimatrix(nodeParentEventCR, 1, RF_eventTypeSize, 1, localEventTimeSize);
      free_uimatrix(nodeLeftEventCR, 1, RF_eventTypeSize, 1, localEventTimeSize);
      switch(RF_splitRule) {
      case SURV_CR_LAU:
        free_uimatrix(nodeParentInclusiveAtRisk, 1, RF_eventTypeSize, 1, localEventTimeSize);
        free_uimatrix(nodeLeftInclusiveAtRisk, 1, RF_eventTypeSize, 1, localEventTimeSize);
        break;
      case SURV_CR_LOG:
        break;
      default:
        break;
      }
    }
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
    unstackRandomCovariates(treeID,
                            parent,
                            randomCovariateIndex,
                            uniformCovariateSize,
                            cdf,
                            cdfSize,
                            cdfSort,
                            density,
                            densitySize,
                            densitySwap,
                            repMembrSize);
  }  
  unstackPreSplit(preliminaryResult,
                  multImpFlag,
                  FALSE,  
                  repMembrSize,
                  splitVector,
                  nonMissMembrIndxStatic);
  result = summarizeSplitResult(*splitParameterMax,
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                splitStatistic,
                                deltaMax);
  return result;
}
char l2Impute (uint    treeID,
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
               char    multImpFlag) {
  uint   *randomCovariateIndex;
  uint    uniformCovariateIndex;
  uint    uniformCovariateSize;
  double *cdf;
  uint    cdfSize;
  uint   *cdfSort;
  uint   *density;
  uint    densitySize;
  uint  **densitySwap;
  uint     covariate;
  double  *splitVector;
  uint     splitVectorSize;
  uint nonMissMembrSize, nonMissMembrSizeStatic;
  uint *nonMissMembrIndx, *nonMissMembrIndxStatic;
  uint   *indxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSize;
  char *localSplitIndicator;
  uint splitLength;
  void *splitVectorPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char preliminaryResult, result;
  double delta, deltaMax;
  uint j, k, m;
  localSplitIndicator    = NULL;  
  mwcpSizeAbsolute       = 0;     
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        repMembrSize,
                                        repMembrIndx,
                                        & nonMissMembrSizeStatic,
                                        & nonMissMembrIndxStatic,
                                        & splitVector,
                                        & parent -> mean,
                                        multImpFlag,
                                        FALSE);
  if (preliminaryResult) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    stackRandomCovariates(treeID,
                          parent,
                          repMembrSize,
                          multImpFlag,
                          & randomCovariateIndex,
                          & uniformCovariateSize,
                          & cdf,
                          & cdfSize,
                          & cdfSort,
                          & density,
                          & densitySize,
                          & densitySwap);
    uint *localEventTimeCount, *localEventTimeIndex;
    uint  localEventTimeSize;
    uint *nodeParentEvent,  *nodeLeftEvent,  *nodeRightEvent;
    uint *nodeParentAtRisk, *nodeLeftAtRisk, *nodeRightAtRisk;
    double *localRatio, *localSurvival, *l2Impute;
    double deltaNum, deltaDen;
    uint   tIndx;
    localEventTimeSize = 0;  
    delta = deltaNum = 0;  
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      stackAndGetSplitSurv(treeID,
                           repMembrIndx,
                           repMembrSize,
                           nonMissMembrIndxStatic,
                           nonMissMembrSizeStatic,
                           & localEventTimeCount,
                           & localEventTimeIndex,
                           & localEventTimeSize,
                           & nodeParentEvent,
                           & nodeParentAtRisk,
                           & nodeLeftEvent,
                           & nodeLeftAtRisk,
                           & nodeRightEvent,
                           & nodeRightAtRisk);
      stackAndGetSplitSurvL2(treeID,
                             parent,
                             localEventTimeSize,
                             localEventTimeIndex,
                             nodeParentEvent,
                             nodeParentAtRisk,
                             & localRatio,
                             & localSurvival);
      stackAndGetL2Impute(treeID,
                          parent,
                          repMembrIndx,
                          repMembrSize,
                          nonMissMembrIndxStatic,
                          nonMissMembrSizeStatic,
                          localEventTimeSize,
                          localSurvival,
                          & l2Impute);
    }
    uint actualCovariateCount = 0;
    uint candidateCovariateCount = 0;
    while (selectRandomCovariates(treeID,
                                  parent,
                                  repMembrIndx,
                                  repMembrSize,
                                  randomCovariateIndex,
                                  & uniformCovariateSize,
                                  & uniformCovariateIndex,
                                  cdf,
                                  & cdfSize,
                                  cdfSort,
                                  density,
                                  & densitySize,
                                  densitySwap,
                                  & covariate,
                                  & actualCovariateCount,
                                  & candidateCovariateCount,
                                  splitVector,
                                  & splitVectorSize,
                                  & indxx,
                                  nonMissMembrSizeStatic,
                                  nonMissMembrIndxStatic,
                                  & nonMissMembrSize,
                                  & nonMissMembrIndx,
                                  multImpFlag)) {
      splitLength = stackAndConstructSplitVector(treeID,
                                                 repMembrSize,
                                                 covariate,
                                                 splitVector,
                                                 splitVectorSize,
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & splitVectorPtr);
      if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
        stackAndGetSplitSurv(treeID,
                             repMembrIndx,
                             repMembrSize,
                             nonMissMembrIndx,
                             nonMissMembrSize,
                             & localEventTimeCount,
                             & localEventTimeIndex,
                             & localEventTimeSize,
                             & nodeParentEvent,
                             & nodeParentAtRisk,
                             & nodeLeftEvent,
                             & nodeLeftAtRisk,
                             & nodeRightEvent,
                             & nodeRightAtRisk);
        stackAndGetSplitSurvL2(treeID,
                               parent,
                               localEventTimeSize,
                               localEventTimeIndex,
                               nodeParentEvent,
                               nodeParentAtRisk,
                               & localRatio,
                               & localSurvival);
      }
      if (localEventTimeSize > 0) {
        for (j = 1; j <= repMembrSize; j++) {
          localSplitIndicator[j] = NEITHER;
        }
        leftSize = 0;
        priorMembrIter = 0;
        if (factorFlag == FALSE) {
          for (j = 1; j <= nonMissMembrSize; j++) {
            localSplitIndicator[ nonMissMembrIndx[indxx[j]] ] = RIGHT;
          }
          for (m = 1; m <= localEventTimeSize; m++) {
            nodeLeftEvent[m] = nodeLeftAtRisk[m] = 0;
          }
        }
        for (j = 1; j < splitLength; j++) {
          if (factorFlag == TRUE) {
            priorMembrIter = 0;
            leftSize = 0;
          }
          virtuallySplitNode(treeID,
                             factorFlag,
                             mwcpSizeAbsolute,
                             covariate,
                             repMembrIndx,
                             repMembrSize,
                             nonMissMembrIndx,
                             nonMissMembrSize,
                             indxx,
                             splitVectorPtr,
                             j,
                             localSplitIndicator,
                             & leftSize,
                             priorMembrIter,
                             & currentMembrIter);
          if (factorFlag == TRUE) {
            for (m = 1; m <= localEventTimeSize; m++) {
              nodeLeftEvent[m] = nodeLeftAtRisk[m] = 0;
            }
            for (k = 1; k <= nonMissMembrSize; k++) {
              if (localSplitIndicator[  nonMissMembrIndx[indxx[k]]  ] == LEFT) {
                tIndx = 0;  
                for (m = 1; m <= localEventTimeSize; m++) {
                  if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]) {
                    tIndx = m;
                    nodeLeftAtRisk[tIndx] ++;
                  }
                  else {
                    m = localEventTimeSize;
                  }
                }
                if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
                  nodeLeftEvent[tIndx] ++;
                }
              }
              else {
              }
            } 
          }
          else {
            for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
              tIndx = 0;  
              for (m = 1; m <= localEventTimeSize; m++) {
                if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]) {
                  tIndx = m;
                  nodeLeftAtRisk[tIndx] ++;
                }
                else {
                  m = localEventTimeSize;
                }
              }
              if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
                nodeLeftEvent[tIndx] ++;
              }
            }
          }
          delta = deltaNum = deltaDen =  0.0;
          for (k=1; k <= localEventTimeSize; k++) {
            deltaNum = deltaNum + ((double) nodeLeftEvent[k] - ((double) ( nodeLeftAtRisk[k] * nodeParentEvent[k]) / nodeParentAtRisk[k]));
            if (nodeParentAtRisk[k] >= 2) {
              deltaDen = deltaDen + (
                                     ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k]) *
                                     (1.0 - ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k])) *
                                     ((double) (nodeParentAtRisk[k] - nodeParentEvent[k]) / (nodeParentAtRisk[k] - 1)) * nodeParentEvent[k]
                                     );
            }
          }
          deltaNum = fabs(deltaNum);
          deltaDen = sqrt(deltaDen);
          if (deltaDen <= EPSILON) {
            if (deltaNum <= EPSILON) {
              delta = 0.0;
            }
            else {
              delta = deltaNum / deltaDen;
            }
          }
          else {
            delta = deltaNum / deltaDen;
          }
          updateMaximumSplit(treeID,
                             parent,
                             delta,
                             candidateCovariateCount,
                             covariate,
                             j,
                             factorFlag,
                             mwcpSizeAbsolute,
                             repMembrSize,
                             localSplitIndicator,
                             & deltaMax,
                             splitParameterMax,
                             splitValueMaxCont,
                             splitValueMaxFactSize,
                             splitValueMaxFactPtr,
                             splitVectorPtr,
                             splitIndicator);
          if (factorFlag == FALSE) {
            priorMembrIter = currentMembrIter - 1;
          }
        }  
      }  
      else {
      }
      unstackSplitVector(treeID,
                         splitVectorSize,
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         mwcpSizeAbsolute,
                         splitVectorPtr);
      unselectRandomCovariates(treeID,
                               parent,
                               repMembrSize,
                               indxx,
                               nonMissMembrSizeStatic,
                               nonMissMembrIndx,
                               multImpFlag);
      if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
        unstackSplitSurv(localEventTimeCount,
                         localEventTimeIndex,
                         localEventTimeSize,
                         nodeParentEvent,
                         nodeParentAtRisk,
                         nodeLeftEvent,
                         nodeLeftAtRisk,
                         nodeRightEvent,
                         nodeRightAtRisk);
        unstackAndGetSplitSurvL2(localEventTimeSize,
                                 localRatio,
                                 localSurvival);
      }
    }  
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      unstackSplitSurv(localEventTimeCount,
                       localEventTimeIndex,
                       localEventTimeSize,
                       nodeParentEvent,
                       nodeParentAtRisk,
                       nodeLeftEvent,
                       nodeLeftAtRisk,
                       nodeRightEvent,
                       nodeRightAtRisk);
      unstackAndGetSplitSurvL2(localEventTimeSize,
                               localRatio,
                               localSurvival);
    }
    unstackRandomCovariates(treeID,
                            parent,
                            randomCovariateIndex,
                            uniformCovariateSize,
                            cdf,
                            cdfSize,
                            cdfSort,
                            density,
                            densitySize,
                            densitySwap,
                            repMembrSize);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
  }  
  unstackPreSplit(preliminaryResult,
                  multImpFlag,
                  FALSE, 
                  repMembrSize,
                  splitVector,
                  nonMissMembrIndxStatic);
  result = summarizeSplitResult(*splitParameterMax,
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
void getMembrCountOnly (uint       treeID,
                        Terminal  *parent,
                        uint      *repMembrIndx,
                        uint       repMembrSize,
                        uint      *allMembrIndx,
                        uint       allMembrSize,
                        uint      *rmbrIterator) {
  if ( (!(RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ||
       ( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2))) { 
    parent -> membrCount = allMembrSize;
  }
  else {
    parent -> membrCount = repMembrSize;
    if (RF_optHigh & OPT_MEMB_OUTG) {
      RF_TN_RCNT_ptr[treeID][parent -> nodeID] = RF_tTermList[treeID][parent -> nodeID] -> membrCount;
    }
    if (RF_optHigh & OPT_MEMB_INCG) {
      parent -> membrCount = RF_TN_RCNT_ptr[treeID][parent -> nodeID];
    }
  }
  if ((parent -> membrCount) == 0) {
    if (!(RF_opt & OPT_OUTC_TYPE)) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Zero node count encountered in (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
}
char unsupervisedSplit (uint    treeID,
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
                        char    multImpFlag) {
  uint   *randomCovariateIndex;
  uint    uniformCovariateIndex;
  uint    uniformCovariateSize;
  double *cdf;
  uint    cdfSize;
  uint   *cdfSort;
  uint   *density;
  uint    densitySize;
  uint  **densitySwap;
  uint     covariate;
  double  *splitVector;
  uint     splitVectorSize;
  uint nonMissMembrSize, nonMissMembrSizeStatic;
  uint *nonMissMembrIndx, *nonMissMembrIndxStatic;
  uint   *indxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSize;
  char *localSplitIndicator;
  uint splitLength;
  void *splitVectorPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char preliminaryResult, result;
  double delta, deltaMax;
  uint   deltaNorm;
  uint i, j, k, p, r;
  localSplitIndicator    = NULL;  
  mwcpSizeAbsolute       = 0;     
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        repMembrSize,
                                        repMembrIndx,
                                        & nonMissMembrSizeStatic,
                                        & nonMissMembrIndxStatic,
                                        & splitVector,
                                        & parent -> mean,
                                        multImpFlag,
                                        TRUE);
  if (preliminaryResult) {
    char   *impurity   = cvector(1, RF_randomResponseCount);
    double *mean       = dvector(1, RF_randomResponseCount);
    double *variance   = dvector(1, RF_randomResponseCount);
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    stackRandomCovariates(treeID,
                          parent,
                          repMembrSize,
                          multImpFlag,
                          & randomCovariateIndex,
                          & uniformCovariateSize,
                          & cdf,
                          & cdfSize,
                          & cdfSort,
                          & density,
                          & densitySize,
                          & densitySwap);
    uint  **parentClassProp = (uint **) new_vvector(1, RF_randomResponseCount, NRUTIL_UPTR);
    uint  **leftClassProp   = (uint **) new_vvector(1, RF_randomResponseCount, NRUTIL_UPTR);
    uint  **rghtClassProp   = (uint **) new_vvector(1, RF_randomResponseCount, NRUTIL_UPTR);
    double *sumLeft      = dvector(1, RF_randomResponseCount);
    double *sumRght      = dvector(1, RF_randomResponseCount);
    double *sumRghtSave  = dvector(1, RF_randomResponseCount);
    uint *pseudoResponseClassSize = uivector(1, RF_randomResponseCount);
    uint *pseudoResponse = uivector(1, RF_randomResponseCount);
    char **secondNonMissMembrFlag = (char **) new_vvector(1, RF_randomResponseCount, NRUTIL_CPTR);
    uint  *secondNonMissMembrSize =              uivector(1, RF_randomResponseCount);
    uint  *secondNonMissMembrLeftSize =          uivector(1, RF_randomResponseCount);
    uint  *secondNonMissMembrRghtSize =          uivector(1, RF_randomResponseCount);
    char  *tempNonMissMembrFlag = 0;
    uint  *tempNonMissMembrIndx;
    char   mResponseFlag;
    uint   localIndex = 0; 
    uint   localSize;
    char    nonMissImpuritySummary;
    double sumLeftSqr, sumRghtSqr;
    uint actualCovariateCount = 0;
    uint candidateCovariateCount = 0;
    while (selectRandomCovariates(treeID,
                                  parent,
                                  repMembrIndx,
                                  repMembrSize,
                                  randomCovariateIndex,
                                  & uniformCovariateSize,
                                  & uniformCovariateIndex,
                                  cdf,
                                  & cdfSize,
                                  cdfSort,
                                  density,
                                  & densitySize,
                                  densitySwap,
                                  & covariate,
                                  & actualCovariateCount,
                                  & candidateCovariateCount,
                                  splitVector,
                                  & splitVectorSize,
                                  & indxx,
                                  nonMissMembrSizeStatic,
                                  nonMissMembrIndxStatic,
                                  & nonMissMembrSize,
                                  & nonMissMembrIndx,
                                  multImpFlag)) {
      uint *pseudoResponseIndex = uivector(1, RF_xSize);
      for (i = 1; i <= RF_xSize; i++) {
        pseudoResponseIndex[i] = i;
      }
      pseudoResponseIndex[covariate] = pseudoResponseIndex[RF_xSize];
      localSize = RF_xSize - 1;
      for (r = 1; r <= RF_randomResponseCount; r++) {
        pseudoResponse[r] = sampleUniformlyFromVector(treeID, pseudoResponseIndex, RF_xSize, & localIndex);
        pseudoResponseIndex[localIndex] = pseudoResponseIndex[localSize];
        localSize --;
      }
      free_uivector(pseudoResponseIndex, 1, RF_xSize);
      if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
        tempNonMissMembrFlag = cvector(1, nonMissMembrSize);
        tempNonMissMembrIndx = uivector(1, nonMissMembrSize);
        for (k = 1; k <= nonMissMembrSize; k++) {
          tempNonMissMembrFlag[k] = TRUE;
          tempNonMissMembrIndx[k] = k;
        }
        for (r = 1; r <= RF_randomResponseCount; r++) {
          secondNonMissMembrFlag[r] = tempNonMissMembrFlag;
          secondNonMissMembrSize[r] = nonMissMembrSize;
        }
        nonMissImpuritySummary = FALSE;
        for (r = 1; r <= RF_randomResponseCount; r++)  {
          impurity[r] = getVariance(repMembrSize,
                                    repMembrIndx,
                                    secondNonMissMembrSize[r],
                                    tempNonMissMembrIndx,
                                    RF_observation[treeID][pseudoResponse[r]],
                                    &mean[r],
                                    &variance[r]);
          nonMissImpuritySummary = nonMissImpuritySummary | impurity[r];
          secondNonMissMembrLeftSize[r] = secondNonMissMembrRghtSize[r] = 0;
        }
        free_uivector(tempNonMissMembrIndx, 1, nonMissMembrSize);
      }
      else {
        tempNonMissMembrIndx = uivector(1, nonMissMembrSize);
        nonMissImpuritySummary = FALSE;
        for (r = 1; r <= RF_randomResponseCount; r++)  {
          secondNonMissMembrFlag[r] = cvector(1, nonMissMembrSize);
          j = 0;
          for (k = 1; k <= nonMissMembrSize; k++) {
            mResponseFlag = FALSE;
            if (RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
              if (RF_mpSign[pseudoResponse[r]][RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]] == 1) {
                mResponseFlag = TRUE;
              }
            }
            if (!mResponseFlag) {
              j ++;
              tempNonMissMembrIndx[j] = nonMissMembrIndx[indxx[k]];
              secondNonMissMembrFlag[r][k] = TRUE;
            }
            else {
              secondNonMissMembrFlag[r][k] = FALSE;
            }
          }  
          secondNonMissMembrSize[r] = j;
          impurity[r] = getVariance(repMembrSize,
                                    repMembrIndx,
                                    secondNonMissMembrSize[r],
                                    tempNonMissMembrIndx,
                                    RF_observation[treeID][pseudoResponse[r]],
                                    &mean[r],
                                    &variance[r]);
          nonMissImpuritySummary = nonMissImpuritySummary | impurity[r];
          secondNonMissMembrLeftSize[r] = secondNonMissMembrRghtSize[r] = 0;
        }  
        free_uivector(tempNonMissMembrIndx, 1, nonMissMembrSize);
      }  
      if (nonMissImpuritySummary) {
        for (j = 1; j <= repMembrSize; j++) {
          localSplitIndicator[j] = NEITHER;
        }
        for (r = 1; r <= RF_randomResponseCount; r++) {
          pseudoResponseClassSize[r] = 0;
          parentClassProp[r] = leftClassProp[r] = rghtClassProp[r] = NULL;
          sumLeft[r] = sumRght[r] = sumRghtSave[r] = 0.0;
        }
        for (r = 1; r <= RF_randomResponseCount; r++) {
          if (impurity[r]) {
            if ((strcmp(RF_xType[pseudoResponse[r]], "C") == 0) ||
                (strcmp(RF_xType[pseudoResponse[r]], "I") == 0)) {
              pseudoResponseClassSize[r] = RF_xFactorSize[RF_xFactorMap[pseudoResponse[r]]];
              parentClassProp[r] = uivector(1, pseudoResponseClassSize[r]);
              leftClassProp[r]   = uivector(1, pseudoResponseClassSize[r]);
              rghtClassProp[r]   = uivector(1, pseudoResponseClassSize[r]);
              for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                parentClassProp[r][p] = 0;
              }
              for (j = 1; j <= nonMissMembrSize; j++) {
                if (secondNonMissMembrFlag[r][j] == TRUE) {
                  parentClassProp[r][ (uint) RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[j]]] ]] ++;
                }
              }
            }
            else {
              sumRghtSave[r] = 0.0;
              for (j = 1; j <= nonMissMembrSize; j++) {
                if (secondNonMissMembrFlag[r][j] == TRUE) {
                  sumRghtSave[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[j]]] ] - mean[r];
                }
              }
            }
          }  
        }  
        leftSize = 0;
        priorMembrIter = 0;
        splitLength = stackAndConstructSplitVector(treeID,
                                                   repMembrSize,
                                                   covariate,
                                                   splitVector,
                                                   splitVectorSize,
                                                   & factorFlag,
                                                   & deterministicSplitFlag,
                                                   & mwcpSizeAbsolute,
                                                   & splitVectorPtr);
        if (factorFlag == FALSE) {
          for (j = 1; j <= nonMissMembrSize; j++) {
            localSplitIndicator[ nonMissMembrIndx[indxx[j]] ] = RIGHT;
          }
          for (r = 1; r <= RF_randomResponseCount; r++) {
            if (impurity[r]) {
              if ((strcmp(RF_xType[pseudoResponse[r]], "C") == 0) ||
                  (strcmp(RF_xType[pseudoResponse[r]], "I") == 0)) {
                for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                  rghtClassProp[r][p] = parentClassProp[r][p];
                  leftClassProp[r][p] = 0;
                }
              }
              else {
                sumRght[r] = sumRghtSave[r];
                sumLeft[r] = 0.0;
              }
              secondNonMissMembrLeftSize[r] = 0;
              secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r];
            }
          }
        }
        for (j = 1; j < splitLength; j++) {
          if (factorFlag == TRUE) {
            priorMembrIter = 0;
            leftSize = 0;
            for (r = 1; r <= RF_rSize; r++) {
              secondNonMissMembrLeftSize[r] = 0;
              secondNonMissMembrRghtSize[r] = 0;
            }
          }
          virtuallySplitNode(treeID,
                                factorFlag,
                                mwcpSizeAbsolute,
                                covariate,
                                repMembrIndx,
                                repMembrSize,
                                nonMissMembrIndx,
                                nonMissMembrSize,
                                indxx,
                                splitVectorPtr,
                                j,
                                localSplitIndicator,
                                & leftSize,
                                priorMembrIter,
                                & currentMembrIter);
          delta = 0.0;
          deltaNorm = 0;
          for (r = 1; r <= RF_randomResponseCount; r++) {
            if (RF_opt & OPT_USPV_STAT) {
              parent -> urStat[r] = pseudoResponse[r];
            }
            if (impurity[r]) {
              if (factorFlag == TRUE) {
                if ((strcmp(RF_xType[pseudoResponse[r]], "C") == 0) ||
                    (strcmp(RF_xType[pseudoResponse[r]], "I") == 0)) {
                  for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                    leftClassProp[r][p] = 0;
                  }
                  for (k = 1; k <= nonMissMembrSize; k++) {
                    if (secondNonMissMembrFlag[r][k] == TRUE) {
                      if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                        leftClassProp[r][ (uint) RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]  ++;
                        secondNonMissMembrLeftSize[r] ++;
                      }
                      else {
                      }
                    }
                  }
                  for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                    rghtClassProp[r][p] = parentClassProp[r][p] - leftClassProp[r][p];
                  }
                }
                else {
                  sumLeft[r] = sumRght[r] = 0.0;
                  for (k = 1; k <= nonMissMembrSize; k++) {
                    if (secondNonMissMembrFlag[r][k] == TRUE) {
                      if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                        sumLeft[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                        secondNonMissMembrLeftSize[r] ++;
                      }
                      else {
                        sumRght[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                      }
                    }
                  }
                }
                secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r] - secondNonMissMembrLeftSize[r];
              }
              else {
                for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                  if (secondNonMissMembrFlag[r][k] == TRUE) {
                    if ((strcmp(RF_xType[pseudoResponse[r]], "C") == 0) ||
                        (strcmp(RF_xType[pseudoResponse[r]], "I") == 0)) {
                      leftClassProp[r][(uint) RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]] ++;
                      rghtClassProp[r][(uint) RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]] --;
                    }
                    else {
                      sumLeft[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                      sumRght[r] -= RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                    }
                    secondNonMissMembrLeftSize[r] ++;
                    secondNonMissMembrRghtSize[r] --;
                  }
                }
              }  
              if ((secondNonMissMembrLeftSize[r] > 0) && (secondNonMissMembrRghtSize[r] > 0)) {
                deltaNorm ++;
                if ((strcmp(RF_xType[pseudoResponse[r]], "C") == 0) ||
                    (strcmp(RF_xType[pseudoResponse[r]], "I") == 0)) {
                  sumLeft[1] = sumRght[1] = 0;
                  for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                    sumLeft[1] += (double) upower(leftClassProp[r][p], 2);
                    sumRght[1] += (double) upower(rghtClassProp[r][p], 2);
                  }
                  sumLeftSqr = sumLeft[1] / secondNonMissMembrLeftSize[r];
                  sumRghtSqr  = sumRght[1] / secondNonMissMembrRghtSize[r];
                }
                else {
                  sumLeftSqr = pow (sumLeft[r], 2.0) / (secondNonMissMembrLeftSize[r] * variance[r]);
                  sumRghtSqr = pow (sumRght[r], 2.0) / (secondNonMissMembrRghtSize[r] * variance[r]);
                }
                delta += sumLeftSqr + sumRghtSqr;
              }
            }  
          }  
          if (deltaNorm > 0) {
            delta = delta / (double) deltaNorm;
          }
          else {
            delta = NA_REAL;
          }
          updateMaximumSplit(treeID,
                             parent,
                             delta,
                             candidateCovariateCount,
                             covariate,
                             j,
                             factorFlag,
                             mwcpSizeAbsolute,
                             repMembrSize,
                             localSplitIndicator,
                             & deltaMax,
                             splitParameterMax,
                             splitValueMaxCont,
                             splitValueMaxFactSize,
                             splitValueMaxFactPtr,
                             splitVectorPtr,
                             splitIndicator);
          if (factorFlag == FALSE) {
            priorMembrIter = currentMembrIter - 1;
          }
        }  
        unstackSplitVector(treeID,
                           splitVectorSize,
                           splitLength,
                           factorFlag,
                           deterministicSplitFlag,
                           mwcpSizeAbsolute,
                           splitVectorPtr);
        for (r = 1; r <= RF_randomResponseCount; r++) {
          if (impurity[r]) {
            if ((strcmp(RF_xType[pseudoResponse[r]], "C") == 0) ||
                (strcmp(RF_xType[pseudoResponse[r]], "I") == 0)) {
              free_uivector (parentClassProp[r], 1, pseudoResponseClassSize[r]);
              free_uivector (leftClassProp[r],   1, pseudoResponseClassSize[r]);
              free_uivector (rghtClassProp[r],   1, pseudoResponseClassSize[r]);
            }
            else {
            }
          }
        }
      }  
      unselectRandomCovariates(treeID,
                               parent,
                               repMembrSize,
                               indxx,
                               nonMissMembrSizeStatic,
                               nonMissMembrIndx,
                               multImpFlag);
      if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
        free_cvector(tempNonMissMembrFlag, 1, nonMissMembrSize);
      }
      else {
        for (r = 1; r <= RF_randomResponseCount; r++)  {
          free_cvector(secondNonMissMembrFlag[r], 1, nonMissMembrSize);
        }
      }
    }  
    free_new_vvector(parentClassProp, 1, RF_randomResponseCount, NRUTIL_UPTR);
    free_new_vvector(leftClassProp,   1, RF_randomResponseCount, NRUTIL_UPTR);
    free_new_vvector(rghtClassProp,   1, RF_randomResponseCount, NRUTIL_UPTR);
    free_dvector(sumLeft,     1, RF_randomResponseCount);
    free_dvector(sumRght,     1, RF_randomResponseCount);
    free_dvector(sumRghtSave, 1, RF_randomResponseCount);
    free_uivector(pseudoResponseClassSize, 1, RF_randomResponseCount);
    free_uivector(pseudoResponse, 1, RF_randomResponseCount);
    free_new_vvector(secondNonMissMembrFlag,  1, RF_randomResponseCount, NRUTIL_CPTR);
    free_uivector(secondNonMissMembrSize,     1, RF_randomResponseCount);
    free_uivector(secondNonMissMembrLeftSize, 1, RF_randomResponseCount);
    free_uivector(secondNonMissMembrRghtSize, 1, RF_randomResponseCount);
    unstackRandomCovariates(treeID,
                            parent,
                            randomCovariateIndex,
                            uniformCovariateSize,
                            cdf,
                            cdfSize,
                            cdfSort,
                            density,
                            densitySize,
                            densitySwap,
                            repMembrSize);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
    free_cvector(impurity,   1, RF_randomResponseCount);
    free_dvector(mean,     1, RF_randomResponseCount);
    free_dvector(variance, 1, RF_randomResponseCount);
  }  
  unstackPreSplit(preliminaryResult,
                  multImpFlag,
                  TRUE,  
                  repMembrSize,
                  splitVector,
                  nonMissMembrIndxStatic);
  result = summarizeSplitResult(*splitParameterMax,
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
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
                        char    multImpFlag) {
  uint   *randomCovariateIndex;
  uint    uniformCovariateIndex;
  uint    uniformCovariateSize;
  double *cdf;
  uint    cdfSize;
  uint   *cdfSort;
  uint   *density;
  uint    densitySize;
  uint  **densitySwap;
  uint     covariate;
  double  *splitVector;
  uint     splitVectorSize;
  uint nonMissMembrSize, nonMissMembrSizeStatic;
  uint *nonMissMembrIndx, *nonMissMembrIndxStatic;
  uint   *indxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSize;
  char *localSplitIndicator;
  uint splitLength;
  void *splitVectorPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char preliminaryResult, result;
  double delta, deltaMax;
  uint   deltaNorm;
  uint j, k, p, r;
  localSplitIndicator    = NULL;  
  mwcpSizeAbsolute       = 0;     
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        repMembrSize,
                                        repMembrIndx,
                                        & nonMissMembrSizeStatic,
                                        & nonMissMembrIndxStatic,
                                        & splitVector,
                                        & parent -> mean,
                                        multImpFlag,
                                        TRUE);
  if (preliminaryResult) {
    char   *impurity   = cvector(1, RF_rSize);
    double *mean       = dvector(1, RF_rSize);
    double *variance   = dvector(1, RF_rSize);
    char impuritySummary;
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      impuritySummary = FALSE;
      for (r = 1; r <= RF_rSize; r++)  {
        impurity[r] = getVariance(repMembrSize,
                                  repMembrIndx,
                                  0,
                                  NULL,
                                  RF_response[treeID][r],
                                  &mean[r],
                                  &variance[r]);
        impuritySummary = impuritySummary | impurity[r];
      }
    }
    else {
      impuritySummary = TRUE;
    }
    if (impuritySummary) {
      stackSplitIndicator(repMembrSize, & localSplitIndicator);
      stackRandomCovariates(treeID,
                            parent,
                            repMembrSize,
                            multImpFlag,
                            & randomCovariateIndex,
                            & uniformCovariateSize,
                            & cdf,
                            & cdfSize,
                            & cdfSort,
                            & density,
                            & densitySize,
                            & densitySwap);
      uint **parentClassProp = (uint **) new_vvector(1, RF_rSize, NRUTIL_UPTR);
      uint **leftClassProp   = (uint **) new_vvector(1, RF_rSize, NRUTIL_UPTR);
      uint **rghtClassProp   = (uint **) new_vvector(1, RF_rSize, NRUTIL_UPTR);
      double *sumLeft         = dvector(1, RF_rSize);
      double *sumRght         = dvector(1, RF_rSize);
      double *sumRghtSave     = dvector(1, RF_rSize);
      char **secondNonMissMembrFlag = (char **) new_vvector(1, RF_rSize, NRUTIL_CPTR);
      uint  *secondNonMissMembrSize =           uivector(1, RF_rSize);
      uint  *secondNonMissMembrLeftSize =       uivector(1, RF_rSize);
      uint  *secondNonMissMembrRghtSize =       uivector(1, RF_rSize);
      for (r = 1; r <= RF_rSize; r++) {
        parentClassProp[r] = leftClassProp[r] = rghtClassProp[r] = NULL;
        if ((strcmp(RF_rType[r], "C") == 0) ||
            (strcmp(RF_rType[r], "I") == 0)) {
          parentClassProp[r] = uivector(1, RF_classLevelSize[RF_rFactorMap[r]]);
          leftClassProp[r]   = uivector(1, RF_classLevelSize[RF_rFactorMap[r]]);
          rghtClassProp[r]   = uivector(1, RF_classLevelSize[RF_rFactorMap[r]]);
        }
        else {
        }
      }  
      char  *tempNonMissMembrFlag = 0;
      uint  *tempNonMissMembrIndx;
      char   mResponseFlag;
      char   nonMissImpuritySummary;
      double partialLeft, partialRght;
      uint actualCovariateCount = 0;
      uint candidateCovariateCount = 0;
      while (selectRandomCovariates(treeID,
                                    parent,
                                    repMembrIndx,
                                    repMembrSize,
                                    randomCovariateIndex,
                                    & uniformCovariateSize,
                                    & uniformCovariateIndex,
                                    cdf,
                                    & cdfSize,
                                    cdfSort,
                                    density,
                                    & densitySize,
                                    densitySwap,
                                    & covariate,
                                    & actualCovariateCount,
                                    & candidateCovariateCount,
                                    splitVector,
                                    & splitVectorSize,
                                    & indxx,
                                    nonMissMembrSizeStatic,
                                    nonMissMembrIndxStatic,
                                    & nonMissMembrSize,
                                    & nonMissMembrIndx,
                                    multImpFlag)) {
        if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
          tempNonMissMembrFlag = cvector(1, nonMissMembrSize);
          for (k = 1; k <= nonMissMembrSize; k++) {
            tempNonMissMembrFlag[k] = TRUE;
          }
          for (r = 1; r <= RF_rSize; r++) {
              secondNonMissMembrFlag[r] = tempNonMissMembrFlag;
              secondNonMissMembrSize[r] = nonMissMembrSize;
          }
          nonMissImpuritySummary = TRUE;
        }
        else {
          tempNonMissMembrIndx = uivector(1, nonMissMembrSize);
          nonMissImpuritySummary = FALSE;
          for (r = 1; r <= RF_rSize; r++)  {
            secondNonMissMembrFlag[r] = cvector(1, nonMissMembrSize);
            j = 0;
            for (k = 1; k <= nonMissMembrSize; k++) {
              mResponseFlag = FALSE;
              if (RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
                if (RF_mpSign[r][RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]] == 1) {
                  mResponseFlag = TRUE;
                }
              }
              if (!mResponseFlag) {
                j ++;
                tempNonMissMembrIndx[j] = nonMissMembrIndx[indxx[k]];
                secondNonMissMembrFlag[r][k] = TRUE;
              }
              else {
                secondNonMissMembrFlag[r][k] = FALSE;
              }
            }  
            secondNonMissMembrSize[r] = j;
            impurity[r] = getVariance(repMembrSize,
                                      repMembrIndx,
                                      secondNonMissMembrSize[r],
                                      tempNonMissMembrIndx,
                                      RF_response[treeID][r],
                                      &mean[r],
                                      &variance[r]);
            nonMissImpuritySummary = nonMissImpuritySummary | impurity[r];
            secondNonMissMembrLeftSize[r] = secondNonMissMembrRghtSize[r] = 0;
          }  
          free_uivector(tempNonMissMembrIndx, 1, nonMissMembrSize);
        }  
        if (nonMissImpuritySummary) {
          for (j = 1; j <= repMembrSize; j++) {
            localSplitIndicator[j] = NEITHER;
          }
          for (r = 1; r <= RF_rSize; r++) {
            if (impurity[r]) {
              if ((strcmp(RF_rType[r], "C") == 0) ||
                  (strcmp(RF_rType[r], "I") == 0)) {
                for (p=1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                  parentClassProp[r][p] = 0;
                }
                for (j = 1; j <= nonMissMembrSize; j++) {
                  if (secondNonMissMembrFlag[r][j] == TRUE) {
                    parentClassProp[r][RF_classLevelIndex[RF_rFactorMap[r]][(uint) RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[j]]] ]]] ++;
                  }
                }
              }
              else {
                sumRghtSave[r] = 0.0;
                for (j = 1; j <= nonMissMembrSize; j++) {
                  if (secondNonMissMembrFlag[r][j] == TRUE) {
                    sumRghtSave[r] += RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[j]]] ] - mean[r];
                  }
                }
              }
            }  
          }  
          leftSize = 0;
          priorMembrIter = 0;
          splitLength = stackAndConstructSplitVector(treeID,
                                                     repMembrSize,
                                                     covariate,
                                                     splitVector,
                                                     splitVectorSize,
                                                     & factorFlag,
                                                     & deterministicSplitFlag,
                                                     & mwcpSizeAbsolute,
                                                     & splitVectorPtr);
          if (factorFlag == FALSE) {
            for (j = 1; j <= nonMissMembrSize; j++) {
              localSplitIndicator[ nonMissMembrIndx[indxx[j]] ] = RIGHT;
            }
            for (r = 1; r <= RF_rSize; r++) {
              if (impurity[r]) {
                if ((strcmp(RF_rType[r], "C") == 0) ||
                    (strcmp(RF_rType[r], "I") == 0)) {
                  for (p=1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                    rghtClassProp[r][p] = parentClassProp[r][p];
                    leftClassProp[r][p] = 0;
                  }
                }
                else {
                  sumRght[r] = sumRghtSave[r];
                  sumLeft[r] = 0.0;
                }
                secondNonMissMembrLeftSize[r] = 0;
                secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r];
              }
            }
          }
          for (j = 1; j < splitLength; j++) {
            if (factorFlag == TRUE) {
              priorMembrIter = 0;
              leftSize = 0;
              for (r = 1; r <= RF_rSize; r++) {
                secondNonMissMembrLeftSize[r] = 0;
                secondNonMissMembrRghtSize[r] = 0;
              }
            }
            virtuallySplitNode(treeID,
                               factorFlag,
                               mwcpSizeAbsolute,
                               covariate,
                               repMembrIndx,
                               repMembrSize,
                               nonMissMembrIndx,
                               nonMissMembrSize,
                               indxx,
                               splitVectorPtr,
                               j,
                               localSplitIndicator,
                               & leftSize,
                               priorMembrIter,
                               & currentMembrIter);
            delta     = 0.0;
            deltaNorm = 0;
            for (r = 1; r <= RF_rSize; r++) {
              if (impurity[r]) {
                if (factorFlag == TRUE) {
                  if ((strcmp(RF_rType[r], "C") == 0) ||
                      (strcmp(RF_rType[r], "I") == 0)) {
                    for (p=1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                      leftClassProp[r][p] = 0;
                    }
                    for (k = 1; k <= nonMissMembrSize; k++) {
                      if (secondNonMissMembrFlag[r][k] == TRUE) {
                        if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                          leftClassProp[r][RF_classLevelIndex[RF_rFactorMap[r]][(uint) RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] ++;
                          secondNonMissMembrLeftSize[r] ++;
                        }
                        else {
                        }
                      }
                    }
                    for (p=1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                      rghtClassProp[r][p] = parentClassProp[r][p] - leftClassProp[r][p];
                    }
                  }
                  else {
                    sumLeft[r] = sumRght[r] = 0.0;
                    for (k = 1; k <= nonMissMembrSize; k++) {
                      if (secondNonMissMembrFlag[r][k] == TRUE) {
                        if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                          sumLeft[r] += RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                          secondNonMissMembrLeftSize[r] ++;
                        }
                        else {
                          sumRght[r] += RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                        }
                      }
                    }
                  }
                  secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r] - secondNonMissMembrLeftSize[r];
                }
                else {
                  for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                    if (secondNonMissMembrFlag[r][k] == TRUE) {
                      if ((strcmp(RF_rType[r], "C") == 0) ||
                          (strcmp(RF_rType[r], "I") == 0)) {
                        leftClassProp[r][RF_classLevelIndex[RF_rFactorMap[r]][(uint) RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] ++;
                        rghtClassProp[r][RF_classLevelIndex[RF_rFactorMap[r]][(uint) RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] --;
                      }
                      else {
                        sumLeft[r] += RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                        sumRght[r] -= RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                      }
                      secondNonMissMembrLeftSize[r] ++;
                      secondNonMissMembrRghtSize[r] --;
                    }
                  }
                }  
                if ((secondNonMissMembrLeftSize[r] > 0) && (secondNonMissMembrRghtSize[r] > 0)) {
                  deltaNorm ++;
                  if ((strcmp(RF_rType[r], "C") == 0) ||
                      (strcmp(RF_rType[r], "I") == 0)) {
                    partialLeft = partialRght = 0;
                    for (p = 1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                      partialLeft += (double) upower(leftClassProp[r][p], 2);
                      partialRght += (double) upower(rghtClassProp[r][p], 2);
                    }
                    partialLeft = partialLeft / secondNonMissMembrLeftSize[r];
                    partialRght = partialRght / secondNonMissMembrRghtSize[r];
                  }
                  else {
                    partialLeft = pow (sumLeft[r], 2.0) / (secondNonMissMembrLeftSize[r] * variance[r]);
                    partialRght = pow (sumRght[r], 2.0) / (secondNonMissMembrRghtSize[r] * variance[r]);
                  }
                  delta += partialLeft + partialRght;
                }
              }  
            }  
            if (deltaNorm > 0) {
              delta = delta / (double) deltaNorm;
            }
            else {
              delta = NA_REAL;
            }
            updateMaximumSplit(treeID,
                               parent,
                               delta,
                               candidateCovariateCount,
                               covariate,
                               j,
                               factorFlag,
                               mwcpSizeAbsolute,
                               repMembrSize,
                               localSplitIndicator,
                               & deltaMax,
                               splitParameterMax,
                               splitValueMaxCont,
                               splitValueMaxFactSize,
                               splitValueMaxFactPtr,
                               splitVectorPtr,
                               splitIndicator);
            if (factorFlag == FALSE) {
              priorMembrIter = currentMembrIter - 1;
            }
          }  
          unstackSplitVector(treeID,
                             splitVectorSize,
                             splitLength,
                             factorFlag,
                             deterministicSplitFlag,
                             mwcpSizeAbsolute,
                             splitVectorPtr);
        }  
        unselectRandomCovariates(treeID,
                                 parent,
                                 repMembrSize,
                                 indxx,
                                 nonMissMembrSizeStatic,
                                 nonMissMembrIndx,
                                 multImpFlag);
        if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
          free_cvector(tempNonMissMembrFlag, 1, nonMissMembrSize);
        }
        else {
          for (r = 1; r <= RF_rSize; r++)  {
            free_cvector(secondNonMissMembrFlag[r], 1, nonMissMembrSize);
          }
        }
      }  
      for (r = 1; r <= RF_rSize; r++) {
        if ((strcmp(RF_rType[r], "C") == 0) ||
            (strcmp(RF_rType[r], "I") == 0)) {
          free_uivector (parentClassProp[r], 1, RF_classLevelSize[RF_rFactorMap[r]]);
          free_uivector (leftClassProp[r], 1, RF_classLevelSize[RF_rFactorMap[r]]);
          free_uivector (rghtClassProp[r], 1, RF_classLevelSize[RF_rFactorMap[r]]);
        }
        else {
        }
      }
      free_new_vvector(parentClassProp, 1, RF_rSize, NRUTIL_UPTR);
      free_new_vvector(leftClassProp,   1, RF_rSize, NRUTIL_UPTR);
      free_new_vvector(rghtClassProp,   1, RF_rSize, NRUTIL_UPTR);
      free_dvector(sumLeft,     1, RF_rSize);
      free_dvector(sumRght,     1, RF_rSize);
      free_dvector(sumRghtSave, 1, RF_rSize);
      free_new_vvector(secondNonMissMembrFlag,  1, RF_rSize, NRUTIL_CPTR);
      free_uivector(secondNonMissMembrSize,     1, RF_rSize);
      free_uivector(secondNonMissMembrLeftSize, 1, RF_rSize);
      free_uivector(secondNonMissMembrRghtSize, 1, RF_rSize);
      unstackRandomCovariates(treeID,
                              parent,
                              randomCovariateIndex,
                              uniformCovariateSize,
                              cdf,
                              cdfSize,
                              cdfSort,
                              density,
                              densitySize,
                              densitySwap,
                              repMembrSize);
      unstackSplitIndicator(repMembrSize, localSplitIndicator);
    }  
    free_cvector(impurity,   1, RF_rSize);
    free_dvector(mean,       1, RF_rSize);
    free_dvector(variance,   1, RF_rSize);
  }  
  unstackPreSplit(preliminaryResult,
                  multImpFlag,
                  TRUE,  
                  repMembrSize,
                  splitVector,
                  nonMissMembrIndxStatic);
  result = summarizeSplitResult(*splitParameterMax,
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
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
                              char    multImpFlag) {
  uint   *randomCovariateIndex;
  uint    uniformCovariateIndex;
  uint    uniformCovariateSize;
  double *cdf;
  uint    cdfSize;
  uint   *cdfSort;
  uint   *density;
  uint    densitySize;
  uint  **densitySwap;
  uint     covariate;
  double  *splitVector;
  uint     splitVectorSize;
  uint nonMissMembrSize, nonMissMembrSizeStatic;
  uint *nonMissMembrIndx, *nonMissMembrIndxStatic;
  uint   *indxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSize;
  char *localSplitIndicator;
  uint splitLength;
  void *splitVectorPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char preliminaryResult, result;
  double delta, deltaPartial, deltaMax;
  uint   deltaNorm;
  uint j, k, m, r;
  localSplitIndicator    = NULL;  
  mwcpSizeAbsolute       = 0;     
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        repMembrSize,
                                        repMembrIndx,
                                        & nonMissMembrSizeStatic,
                                        & nonMissMembrIndxStatic,
                                        & splitVector,
                                        & parent -> mean,
                                        multImpFlag,
                                        TRUE);
  if (preliminaryResult) {
    char   *impurity   = cvector(1, RF_rSize);
    double *mean       = dvector(1, RF_rSize);
    double *variance   = dvector(1, RF_rSize);
    char impuritySummary;
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      impuritySummary = FALSE;
      for (r = 1; r <= RF_rSize; r++)  {
        impurity[r] = getVariance(repMembrSize,
                                  repMembrIndx,
                                  0,
                                  NULL,
                                  RF_response[treeID][r],
                                  &mean[r],
                                  &variance[r]);
        impuritySummary = impuritySummary | impurity[r];
      }
    }
    else {
      impuritySummary = TRUE;
    }
    if (impuritySummary) {
      stackSplitIndicator(repMembrSize, & localSplitIndicator);
      stackRandomCovariates(treeID,
                            parent,
                            repMembrSize,
                            multImpFlag,
                            & randomCovariateIndex,
                            & uniformCovariateSize,
                            & cdf,
                            & cdfSize,
                            & cdfSort,
                            & density,
                            & densitySize,
                            & densitySwap);
      char **secondNonMissMembrFlag = (char **) new_vvector(1, RF_rSize, NRUTIL_CPTR);
      uint  *secondNonMissMembrSize =           uivector(1, RF_rSize);
      uint  *secondNonMissMembrLeftSize =       uivector(1, RF_rSize);
      uint  *secondNonMissMembrRghtSize =       uivector(1, RF_rSize);
      char  *tempNonMissMembrFlag = 0;
      uint  *tempNonMissMembrIndx;
      char   mResponseFlag;
      char   nonMissImpuritySummary;
      uint actualCovariateCount = 0;
      uint candidateCovariateCount = 0;
      while (selectRandomCovariates(treeID,
                                    parent,
                                    repMembrIndx,
                                    repMembrSize,
                                    randomCovariateIndex,
                                    & uniformCovariateSize,
                                    & uniformCovariateIndex,
                                    cdf,
                                    & cdfSize,
                                    cdfSort,
                                    density,
                                    & densitySize,
                                    densitySwap,
                                    & covariate,
                                    & actualCovariateCount,
                                    & candidateCovariateCount,
                                    splitVector,
                                    & splitVectorSize,
                                    & indxx,
                                    nonMissMembrSizeStatic,
                                    nonMissMembrIndxStatic,
                                    & nonMissMembrSize,
                                    & nonMissMembrIndx,
                                    multImpFlag)) {
        if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
          tempNonMissMembrFlag = cvector(1, nonMissMembrSize);
          for (k = 1; k <= nonMissMembrSize; k++) {
            tempNonMissMembrFlag[k] = TRUE;
          }
          for (r = 1; r <= RF_rSize; r++) {
              secondNonMissMembrFlag[r] = tempNonMissMembrFlag;
              secondNonMissMembrSize[r] = nonMissMembrSize;
          }
          nonMissImpuritySummary = TRUE;
        }
        else {
          tempNonMissMembrIndx = uivector(1, nonMissMembrSize);
          nonMissImpuritySummary = FALSE;
          for (r = 1; r <= RF_rSize; r++)  {
            secondNonMissMembrFlag[r] = cvector(1, nonMissMembrSize);
            j = 0;
            for (k = 1; k <= nonMissMembrSize; k++) {
              mResponseFlag = FALSE;
              if (RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
                if (RF_mpSign[r][RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]] == 1) {
                  mResponseFlag = TRUE;
                }
              }
              if (!mResponseFlag) {
                j ++;
                tempNonMissMembrIndx[j] = nonMissMembrIndx[indxx[k]];
                secondNonMissMembrFlag[r][k] = TRUE;
              }
              else {
                secondNonMissMembrFlag[r][k] = FALSE;
              }
            }  
            secondNonMissMembrSize[r] = j;
            impurity[r] = getVariance(repMembrSize,
                                      repMembrIndx,
                                      secondNonMissMembrSize[r],
                                      tempNonMissMembrIndx,
                                      RF_response[treeID][r],
                                      &mean[r],
                                      &variance[r]);
            nonMissImpuritySummary = nonMissImpuritySummary | impurity[r];
            secondNonMissMembrLeftSize[r] = secondNonMissMembrRghtSize[r] = 0;
          }  
          free_uivector(tempNonMissMembrIndx, 1, nonMissMembrSize);
        }  
        if (nonMissImpuritySummary) {
          for (j = 1; j <= repMembrSize; j++) {
            localSplitIndicator[j] = NEITHER;
          }
          leftSize = 0;
          priorMembrIter = 0;
          splitLength = stackAndConstructSplitVector(treeID,
                                                     repMembrSize,
                                                     covariate,
                                                     splitVector,
                                                     splitVectorSize,
                                                     & factorFlag,
                                                     & deterministicSplitFlag,
                                                     & mwcpSizeAbsolute,
                                                     & splitVectorPtr);
          if (factorFlag == FALSE) {
            for (j = 1; j <= nonMissMembrSize; j++) {
              localSplitIndicator[ nonMissMembrIndx[indxx[j]] ] = RIGHT;
            }
            for (r = 1; r <= RF_rSize; r++) {
              if (impurity[r]) {
                secondNonMissMembrLeftSize[r] = 0;
                secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r];
              }
            }
          }
          double *userResponse = dvector(1, nonMissMembrSize);
          char   *userSplitIndicator = cvector(1, nonMissMembrSize);
          for (j = 1; j < splitLength; j++) {
            if (factorFlag == TRUE) {
              priorMembrIter = 0;
              leftSize = 0;
              for (r = 1; r <= RF_rSize; r++) {
                secondNonMissMembrLeftSize[r] = 0;
                secondNonMissMembrRghtSize[r] = 0;
              }
            }
            virtuallySplitNode(treeID,
                               factorFlag,
                               mwcpSizeAbsolute,
                               covariate,
                               repMembrIndx,
                               repMembrSize,
                               nonMissMembrIndx,
                               nonMissMembrSize,
                               indxx,
                               splitVectorPtr,
                               j,
                               localSplitIndicator,
                               & leftSize,
                               priorMembrIter,
                               & currentMembrIter);
            delta        = 0.0;
            deltaPartial = 0.0;
            deltaNorm    = 0;
            for (r = 1; r <= RF_rSize; r++) {
              if (impurity[r]) {
                if (factorFlag == TRUE) {
                  for (k = 1; k <= nonMissMembrSize; k++) {
                    if (secondNonMissMembrFlag[r][k] == TRUE) {
                      if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                        secondNonMissMembrLeftSize[r] ++;
                      }
                      else {
                      }
                    }
                  }
                  secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r] - secondNonMissMembrLeftSize[r];
                }
                else {
                  for (k = priorMembrIter + 1; k < currentMembrIter; k++) {
                    if (secondNonMissMembrFlag[r][k] == TRUE) {
                      secondNonMissMembrLeftSize[r] ++;
                      secondNonMissMembrRghtSize[r] --;
                    }
                  }
                }  
                if ((secondNonMissMembrLeftSize[r] > 0) && (secondNonMissMembrRghtSize[r] > 0)) {
                  m = 0;
                  for (k = 1; k <= nonMissMembrSize; k++) {
                    if (secondNonMissMembrFlag[r][k] == TRUE) {
                      userResponse[++m] = RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                      userSplitIndicator[m] = localSplitIndicator[ nonMissMembrIndx[indxx[k]] ];
                    }
                  }
                  deltaNorm ++;
                  if ((strcmp(RF_rType[r], "C") == 0) ||
                      (strcmp(RF_rType[r], "I") == 0)) {
                    deltaPartial = customFunctionArray[CLAS_FAM][RF_splitCustomIdx](m,
                                                                                    userSplitIndicator,
                                                                                    NULL,
                                                                                    NULL,
                                                                                    0,
                                                                                    0,
                                                                                    NULL,
                                                                                    userResponse,
                                                                                    mean[r],
                                                                                    variance[r],
                                                                                    RF_rFactorSize[RF_rFactorMap[r]]);
                  }
                  else {
                    deltaPartial = (customFunctionArray[REGR_FAM][RF_splitCustomIdx])(m,
                                                                                      userSplitIndicator,
                                                                                      NULL,
                                                                                      NULL,
                                                                                      0,
                                                                                      0,
                                                                                      NULL,
                                                                                      userResponse,
                                                                                      mean[r],
                                                                                      variance[r],
                                                                                      0);
                  }
                  delta += deltaPartial;
                }
              }  
            }  
            if (deltaNorm > 0) {
              delta = delta / (double) deltaNorm;
            }
            else {
              delta = NA_REAL;
            }
            updateMaximumSplit(treeID,
                               parent,
                               delta,
                               candidateCovariateCount,
                               covariate,
                               j,
                               factorFlag,
                               mwcpSizeAbsolute,
                               repMembrSize,
                               localSplitIndicator,
                               & deltaMax,
                               splitParameterMax,
                               splitValueMaxCont,
                               splitValueMaxFactSize,
                               splitValueMaxFactPtr,
                               splitVectorPtr,
                               splitIndicator);
            if (factorFlag == FALSE) {
              priorMembrIter = currentMembrIter - 1;
            }
          }  
          free_dvector (userResponse, 1, nonMissMembrSize);
          free_cvector (userSplitIndicator, 1, nonMissMembrSize);
          unstackSplitVector(treeID,
                             splitVectorSize,
                             splitLength,
                             factorFlag,
                             deterministicSplitFlag,
                             mwcpSizeAbsolute,
                             splitVectorPtr);
        }  
        unselectRandomCovariates(treeID,
                                 parent,
                                 repMembrSize,
                                 indxx,
                                 nonMissMembrSizeStatic,
                                 nonMissMembrIndx,
                                 multImpFlag);
        if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
          free_cvector(tempNonMissMembrFlag, 1, nonMissMembrSize);
        }
        else {
          for (r = 1; r <= RF_rSize; r++)  {
            free_cvector(secondNonMissMembrFlag[r], 1, nonMissMembrSize);
          }
        }
      }  
      free_new_vvector(secondNonMissMembrFlag,  1, RF_rSize, NRUTIL_CPTR);
      free_uivector(secondNonMissMembrSize,     1, RF_rSize);
      free_uivector(secondNonMissMembrLeftSize, 1, RF_rSize);
      free_uivector(secondNonMissMembrRghtSize, 1, RF_rSize);
      unstackRandomCovariates(treeID,
                              parent,
                              randomCovariateIndex,
                              uniformCovariateSize,
                              cdf,
                              cdfSize,
                              cdfSort,
                              density,
                              densitySize,
                              densitySwap,
                              repMembrSize);
      unstackSplitIndicator(repMembrSize, localSplitIndicator);
    }  
    free_cvector(impurity,   1, RF_rSize);
    free_dvector(mean,       1, RF_rSize);
    free_dvector(variance,   1, RF_rSize);
  }  
  unstackPreSplit(preliminaryResult,
                  multImpFlag,
                  TRUE,  
                  repMembrSize,
                  splitVector,
                  nonMissMembrIndxStatic);
  result = summarizeSplitResult(*splitParameterMax,
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
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
                          char    multImpFlag) {
  uint   *randomCovariateIndex;
  uint    uniformCovariateIndex;
  uint    uniformCovariateSize;
  double *cdf;
  uint    cdfSize;
  uint   *cdfSort;
  uint   *density;
  uint    densitySize;
  uint  **densitySwap;
  uint     covariate;
  double  *splitVector;
  uint     splitVectorSize;
  uint nonMissMembrSize, nonMissMembrSizeStatic;
  uint *nonMissMembrIndx, *nonMissMembrIndxStatic;
  uint   *indxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSize;
  char *localSplitIndicator;
  uint splitLength;
  void *splitVectorPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char preliminaryResult, result;
  double delta, deltaMax;
  uint j, k, m;
  localSplitIndicator    = NULL;  
  mwcpSizeAbsolute       = 0;     
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        repMembrSize,
                                        repMembrIndx,
                                        & nonMissMembrSizeStatic,
                                        & nonMissMembrIndxStatic,
                                        & splitVector,
                                        & parent -> mean,
                                        multImpFlag,
                                        FALSE);
  if (preliminaryResult) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    stackRandomCovariates(treeID,
                          parent,
                          repMembrSize,
                          multImpFlag,
                          & randomCovariateIndex,
                          & uniformCovariateSize,
                          & cdf,
                          & cdfSize,
                          & cdfSort,
                          & density,
                          & densitySize,
                          & densitySwap);
    uint *localEventTimeCount, *localEventTimeIndex;
    uint  localEventTimeSize;
    uint *nodeParentEvent,  *nodeLeftEvent,  *nodeRightEvent;
    uint *nodeParentAtRisk, *nodeLeftAtRisk, *nodeRightAtRisk;
    localEventTimeSize = 0;  
    delta = 0;  
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      stackAndGetSplitSurv(treeID,
                           repMembrIndx,
                           repMembrSize,
                           nonMissMembrIndxStatic,
                           nonMissMembrSizeStatic,
                           & localEventTimeCount,
                           & localEventTimeIndex,
                           & localEventTimeSize,
                           & nodeParentEvent,
                           & nodeParentAtRisk,
                           & nodeLeftEvent,
                           & nodeLeftAtRisk,
                           & nodeRightEvent,
                           & nodeRightAtRisk);
    }
    uint actualCovariateCount = 0;
    uint candidateCovariateCount = 0;
    while (selectRandomCovariates(treeID,
                                  parent,
                                  repMembrIndx,
                                  repMembrSize,
                                  randomCovariateIndex,
                                  & uniformCovariateSize,
                                  & uniformCovariateIndex,
                                  cdf,
                                  & cdfSize,
                                  cdfSort,
                                  density,
                                  & densitySize,
                                  densitySwap,
                                  & covariate,
                                  & actualCovariateCount,
                                  & candidateCovariateCount,
                                  splitVector,
                                  & splitVectorSize,
                                  & indxx,
                                  nonMissMembrSizeStatic,
                                  nonMissMembrIndxStatic,
                                  & nonMissMembrSize,
                                  & nonMissMembrIndx,
                                  multImpFlag)) {
      splitLength = stackAndConstructSplitVector(treeID,
                                                 repMembrSize,
                                                 covariate,
                                                 splitVector,
                                                 splitVectorSize,
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & splitVectorPtr);
      if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
        stackAndGetSplitSurv(treeID,
                             repMembrIndx,
                             repMembrSize,
                             nonMissMembrIndx,
                             nonMissMembrSize,
                             & localEventTimeCount,
                             & localEventTimeIndex,
                             & localEventTimeSize,
                             & nodeParentEvent,
                             & nodeParentAtRisk,
                             & nodeLeftEvent,
                             & nodeLeftAtRisk,
                             & nodeRightEvent,
                             & nodeRightAtRisk);
      }
      if (localEventTimeSize > 0) {
        for (j = 1; j <= repMembrSize; j++) {
          localSplitIndicator[j] = NEITHER;
        }
        leftSize = 0;
        priorMembrIter = 0;
        if (factorFlag == FALSE) {
          for (j = 1; j <= nonMissMembrSize; j++) {
            localSplitIndicator[ nonMissMembrIndx[indxx[j]] ] = RIGHT;
          }
        }
        double *userTime  = dvector(1, nonMissMembrSize);
        double *userEvent = dvector(1, nonMissMembrSize);
        double *userEventTime = dvector(1, localEventTimeSize);
        char   *userSplitIndicator = cvector(1, nonMissMembrSize);
        uint   *userSort  = uivector(1, nonMissMembrSize);
        double *tempTime  = dvector(1, nonMissMembrSize);
        for (k = 1; k <= nonMissMembrSize; k++) {
          tempTime[k]  = RF_time[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
        }
        indexx(nonMissMembrSize, tempTime, userSort);
        for (k = 1; k <= nonMissMembrSize; k++) {
          userTime[k]  = RF_time[treeID][ repMembrIndx[nonMissMembrIndx[indxx[userSort[k]]]] ];
          userEvent[k] = RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[userSort[k]]]] ];
        }
        for (m = 1; m <= localEventTimeSize; m++) {
          userEventTime[m] = RF_masterTime[localEventTimeIndex[m]];
        }
        for (j = 1; j < splitLength; j++) {
          if (factorFlag == TRUE) {
            priorMembrIter = 0;
            leftSize = 0;
          }
          virtuallySplitNode(treeID,
                             factorFlag,
                             mwcpSizeAbsolute,
                             covariate,
                             repMembrIndx,
                             repMembrSize,
                             nonMissMembrIndx,
                             nonMissMembrSize,
                             indxx,
                             splitVectorPtr,
                             j,
                             localSplitIndicator,
                             & leftSize,
                             priorMembrIter,
                             & currentMembrIter);
          for (k = 1; k <= nonMissMembrSize; k++) {
            userSplitIndicator[k] = localSplitIndicator[ nonMissMembrIndx[indxx[userSort[k]]] ];
          }
          delta = customFunctionArray[SURV_FAM][RF_splitCustomIdx](nonMissMembrSize,
                                                                   userSplitIndicator,
                                                                   userTime,
                                                                   userEvent,
                                                                   0,
                                                                   localEventTimeSize,
                                                                   userEventTime,
                                                                   NULL,
                                                                   0,
                                                                   0,
                                                                   0);
          updateMaximumSplit(treeID,
                             parent,
                             delta,
                             candidateCovariateCount,
                             covariate,
                             j,
                             factorFlag,
                             mwcpSizeAbsolute,
                             repMembrSize,
                             localSplitIndicator,
                             & deltaMax,
                             splitParameterMax,
                             splitValueMaxCont,
                             splitValueMaxFactSize,
                             splitValueMaxFactPtr,
                             splitVectorPtr,
                             splitIndicator);
          if (factorFlag == FALSE) {
            priorMembrIter = currentMembrIter - 1;
          }
        }  
        free_uivector(userSort, 1, nonMissMembrSize);
        free_dvector(tempTime, 1, nonMissMembrSize);
        free_dvector (userTime, 1, nonMissMembrSize);
        free_dvector (userEvent, 1, nonMissMembrSize);
        free_dvector (userEventTime, 1, localEventTimeSize);
        free_cvector (userSplitIndicator, 1, nonMissMembrSize);
      }  
      else {
      }
      unstackSplitVector(treeID,
                         splitVectorSize,
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         mwcpSizeAbsolute,
                         splitVectorPtr);
      unselectRandomCovariates(treeID,
                               parent,
                               repMembrSize,
                               indxx,
                               nonMissMembrSizeStatic,
                               nonMissMembrIndx,
                               multImpFlag);
      if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
        unstackSplitSurv(localEventTimeCount,
                         localEventTimeIndex,
                         localEventTimeSize,
                         nodeParentEvent,
                         nodeParentAtRisk,
                         nodeLeftEvent,
                         nodeLeftAtRisk,
                         nodeRightEvent,
                         nodeRightAtRisk);
      }
    }  
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      unstackSplitSurv(localEventTimeCount,
                       localEventTimeIndex,
                       localEventTimeSize,
                       nodeParentEvent,
                       nodeParentAtRisk,
                       nodeLeftEvent,
                       nodeLeftAtRisk,
                       nodeRightEvent,
                       nodeRightAtRisk);
    }
    unstackRandomCovariates(treeID,
                            parent,
                            randomCovariateIndex,
                            uniformCovariateSize,
                            cdf,
                            cdfSize,
                            cdfSort,
                            density,
                            densitySize,
                            densitySwap,
                            repMembrSize);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
  }  
  unstackPreSplit(preliminaryResult,
                  multImpFlag,
                  FALSE,  
                  repMembrSize,
                  splitVector,
                  nonMissMembrIndxStatic);
  result = summarizeSplitResult(*splitParameterMax,
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
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
                               char    multImpFlag) {
  uint   *randomCovariateIndex;
  uint    uniformCovariateIndex;
  uint    uniformCovariateSize;
  double *cdf;
  uint    cdfSize;
  uint   *cdfSort;
  uint   *density;
  uint    densitySize;
  uint  **densitySwap;
  uint     covariate;
  double  *splitVector;
  uint     splitVectorSize;
  uint nonMissMembrSize, nonMissMembrSizeStatic;
  uint *nonMissMembrIndx, *nonMissMembrIndxStatic;
  uint   *indxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSize;
  char *localSplitIndicator;
  uint splitLength;
  void *splitVectorPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char preliminaryResult, result;
  double delta, deltaMax;
  uint j, k, m;
  localSplitIndicator    = NULL;  
  mwcpSizeAbsolute       = 0;     
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  preliminaryResult = getPreSplitResult(treeID,
                                        parent,
                                        repMembrSize,
                                        repMembrIndx,
                                        & nonMissMembrSizeStatic,
                                        & nonMissMembrIndxStatic,
                                        & splitVector,
                                        & parent -> mean,
                                        multImpFlag,
                                        FALSE);
  if (preliminaryResult) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    stackRandomCovariates(treeID,
                          parent,
                          repMembrSize,
                          multImpFlag,
                          & randomCovariateIndex,
                          & uniformCovariateSize,
                          & cdf,
                          & cdfSize,
                          & cdfSort,
                          & density,
                          & densitySize,
                          & densitySwap);
    uint *localEventTimeCount, *localEventTimeIndex;
    uint  localEventTimeSize;
    uint *nodeParentEvent,  *nodeLeftEvent,  *nodeRightEvent;
    uint *nodeParentAtRisk, *nodeLeftAtRisk, *nodeRightAtRisk;
    localEventTimeSize = 0;  
    delta = 0;  
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      stackAndGetSplitSurv(treeID,
                           repMembrIndx,
                           repMembrSize,
                           nonMissMembrIndxStatic,
                           nonMissMembrSizeStatic,
                           & localEventTimeCount,
                           & localEventTimeIndex,
                           & localEventTimeSize,
                           & nodeParentEvent,
                           & nodeParentAtRisk,
                           & nodeLeftEvent,
                           & nodeLeftAtRisk,
                           & nodeRightEvent,
                           & nodeRightAtRisk);
   }
    uint actualCovariateCount = 0;
    uint candidateCovariateCount = 0;
    while (selectRandomCovariates(treeID,
                                  parent,
                                  repMembrIndx,
                                  repMembrSize,
                                  randomCovariateIndex,
                                  & uniformCovariateSize,
                                  & uniformCovariateIndex,
                                  cdf,
                                  & cdfSize,
                                  cdfSort,
                                  density,
                                  & densitySize,
                                  densitySwap,
                                  & covariate,
                                  & actualCovariateCount,
                                  & candidateCovariateCount,
                                  splitVector,
                                  & splitVectorSize,
                                  & indxx,
                                  nonMissMembrSizeStatic,
                                  nonMissMembrIndxStatic,
                                  & nonMissMembrSize,
                                  & nonMissMembrIndx,
                                  multImpFlag)) {
      splitLength = stackAndConstructSplitVector(treeID,
                                                 repMembrSize,
                                                 covariate,
                                                 splitVector,
                                                 splitVectorSize,
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & splitVectorPtr);
      if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
        stackAndGetSplitSurv(treeID,
                             repMembrIndx,
                             repMembrSize,
                             nonMissMembrIndx,
                             nonMissMembrSize,
                             & localEventTimeCount,
                             & localEventTimeIndex,
                             & localEventTimeSize,
                             & nodeParentEvent,
                             & nodeParentAtRisk,
                             & nodeLeftEvent,
                             & nodeLeftAtRisk,
                             & nodeRightEvent,
                             & nodeRightAtRisk);
      }
      if (localEventTimeSize > 0) {
        for (j = 1; j <= repMembrSize; j++) {
          localSplitIndicator[j] = NEITHER;
        }
        leftSize = 0;
        priorMembrIter = 0;
        if (factorFlag == FALSE) {
          for (j = 1; j <= nonMissMembrSize; j++) {
            localSplitIndicator[ nonMissMembrIndx[indxx[j]] ] = RIGHT;
          }
        }
        double *userTime  = dvector(1, nonMissMembrSize);
        double *userEvent = dvector(1, nonMissMembrSize);
        double *userEventTime = dvector(1, localEventTimeSize);
        char   *userSplitIndicator = cvector(1, nonMissMembrSize);
        uint   *userSort  = uivector(1, nonMissMembrSize);
        double *tempTime  = dvector(1, nonMissMembrSize);
        for (k = 1; k <= nonMissMembrSize; k++) {
          tempTime[k]  = RF_time[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
        }
        indexx(nonMissMembrSize, tempTime, userSort);
        for (k = 1; k <= nonMissMembrSize; k++) {
          userTime[k]  = RF_time[treeID][ repMembrIndx[nonMissMembrIndx[indxx[userSort[k]]]] ];
          userEvent[k] = RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[userSort[k]]]] ];
        }
        for (m = 1; m <= localEventTimeSize; m++) {
          userEventTime[m] = RF_masterTime[localEventTimeIndex[m]];
        }
        for (j = 1; j < splitLength; j++) {
          if (factorFlag == TRUE) {
            priorMembrIter = 0;
            leftSize = 0;
          }
          virtuallySplitNode(treeID,
                             factorFlag,
                             mwcpSizeAbsolute,
                             covariate,
                             repMembrIndx,
                             repMembrSize,
                             nonMissMembrIndx,
                             nonMissMembrSize,
                             indxx,
                             splitVectorPtr,
                             j,
                             localSplitIndicator,
                             & leftSize,
                             priorMembrIter,
                             & currentMembrIter);
          for (k = 1; k <= nonMissMembrSize; k++) {
            userSplitIndicator[k] = localSplitIndicator[ nonMissMembrIndx[indxx[userSort[k]]] ];
          }
          delta = customFunctionArray[CRSK_FAM][RF_splitCustomIdx](nonMissMembrSize,
                                                                   userSplitIndicator,
                                                                   userTime,
                                                                   userEvent,
                                                                   RF_eventType[RF_eventTypeSize],
                                                                   localEventTimeSize,
                                                                   userEventTime,
                                                                   NULL,
                                                                   0,
                                                                   0,
                                                                   0);
          updateMaximumSplit(treeID,
                             parent,
                             delta,
                             candidateCovariateCount,
                             covariate,
                             j,
                             factorFlag,
                             mwcpSizeAbsolute,
                             repMembrSize,
                             localSplitIndicator,
                             & deltaMax,
                             splitParameterMax,
                             splitValueMaxCont,
                             splitValueMaxFactSize,
                             splitValueMaxFactPtr,
                             splitVectorPtr,
                             splitIndicator);
          if (factorFlag == FALSE) {
            priorMembrIter = currentMembrIter - 1;
          }
        }  
        free_uivector(userSort, 1, nonMissMembrSize);
        free_dvector(tempTime, 1, nonMissMembrSize);
        free_dvector (userTime, 1, nonMissMembrSize);
        free_dvector (userEvent, 1, nonMissMembrSize);
        free_dvector (userEventTime, 1, localEventTimeSize);
        free_cvector (userSplitIndicator, 1, nonMissMembrSize);
      }  
      else {
      }
      unstackSplitVector(treeID,
                         splitVectorSize,
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         mwcpSizeAbsolute,
                         splitVectorPtr);
      unselectRandomCovariates(treeID,
                               parent,
                               repMembrSize,
                               indxx,
                               nonMissMembrSizeStatic,
                               nonMissMembrIndx,
                               multImpFlag);
      if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
        unstackSplitSurv(localEventTimeCount,
                         localEventTimeIndex,
                         localEventTimeSize,
                         nodeParentEvent,
                         nodeParentAtRisk,
                         nodeLeftEvent,
                         nodeLeftAtRisk,
                         nodeRightEvent,
                         nodeRightAtRisk);
      }
    }  
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
        unstackSplitSurv(localEventTimeCount,
                         localEventTimeIndex,
                         localEventTimeSize,
                         nodeParentEvent,
                         nodeParentAtRisk,
                         nodeLeftEvent,
                         nodeLeftAtRisk,
                         nodeRightEvent,
                         nodeRightAtRisk);
    }
    unstackRandomCovariates(treeID,
                            parent,
                            randomCovariateIndex,
                            uniformCovariateSize,
                            cdf,
                            cdfSize,
                            cdfSort,
                            density,
                            densitySize,
                            densitySwap,
                            repMembrSize);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
  }  
  unstackPreSplit(preliminaryResult,
                  multImpFlag,
                  FALSE,  
                  repMembrSize,
                  splitVector,
                  nonMissMembrIndxStatic);
  result = summarizeSplitResult(*splitParameterMax,
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                splitStatistic,
                                deltaMax);
  return result;
}
void stackSplitIndicator(uint   nodeSize,
                         char **localSplitIndicator) {
  if (nodeSize > 0) {
    *localSplitIndicator = cvector(1, nodeSize);
  }
}
void unstackSplitIndicator(uint  nodeSize,
                           char *localSplitIndicator) {
  if (nodeSize > 0) {
    if (nodeSize > 0) {
      free_cvector(localSplitIndicator, 1, nodeSize);
    }
  }
}
void stackSplitEventTime(uint **localEventTimeCount,
                         uint **localEventTimeIndex) {
  *localEventTimeCount = uivector(1, RF_masterTimeSize);
  *localEventTimeIndex = uivector(1, RF_masterTimeSize);
}
void unstackSplitEventTime(uint *localEventTimeCount,
                           uint *localEventTimeIndex) {
  free_uivector(localEventTimeCount, 1, RF_masterTimeSize);
  free_uivector(localEventTimeIndex, 1, RF_masterTimeSize);
}
uint getSplitEventTime(uint   treeID,
                       uint   *repMembrIndx,
                       uint    repMembrSize,
                       uint   *nonMissMembrIndx,
                       uint    nonMissMembrSize,
                       uint   *localEventTimeCount,
                       uint   *localEventTimeIndex) {
  uint i;
  uint eventTimeSize;
  eventTimeSize = 0;
  for (i=1; i <= RF_masterTimeSize; i++) {
    localEventTimeCount[i] = 0;
  }
  for (i = 1; i <= nonMissMembrSize; i++) {
    if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[i]] ] > 0) {
      localEventTimeCount[RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[i]] ]] ++;
    }
  }
  for (i=1; i <= RF_masterTimeSize; i++) {
    if (localEventTimeCount[i] > 0) {
      localEventTimeIndex[++eventTimeSize] = i;
    }
  }
  return (eventTimeSize);
}
void stackSplitEventAndRisk(uint   eventTimeSize,
                            uint **nodeParentEvent,
                            uint **nodeParentAtRisk,
                            uint **nodeLeftEvent,
                            uint **nodeLeftAtRisk,
                            uint **nodeRightEvent,
                            uint **nodeRightAtRisk) {
  if (eventTimeSize > 0) {
    *nodeParentEvent  = uivector(1, eventTimeSize);
    *nodeParentAtRisk = uivector(1, eventTimeSize);
    *nodeLeftEvent  = uivector(1, eventTimeSize);
    *nodeLeftAtRisk = uivector(1, eventTimeSize);
    *nodeRightEvent  = uivector(1, eventTimeSize);
    *nodeRightAtRisk = uivector(1, eventTimeSize);
  }
  else {
    *nodeParentEvent = *nodeParentAtRisk = *nodeLeftEvent  = *nodeLeftAtRisk = *nodeRightEvent  = *nodeRightAtRisk = NULL;
  }
}
void unstackSplitEventAndRisk(uint  eventTimeSize,
                              uint *nodeParentEvent,
                              uint *nodeParentAtRisk,
                              uint *nodeLeftEvent,
                              uint *nodeLeftAtRisk,
                              uint *nodeRightEvent,
                              uint *nodeRightAtRisk) {
  if (eventTimeSize > 0) {
    free_uivector(nodeParentEvent, 1, eventTimeSize);
    free_uivector(nodeParentAtRisk, 1, eventTimeSize);
    free_uivector(nodeLeftEvent, 1, eventTimeSize);
    free_uivector(nodeLeftAtRisk, 1, eventTimeSize);
    free_uivector(nodeRightEvent, 1, eventTimeSize);
    free_uivector(nodeRightAtRisk, 1, eventTimeSize);
  }
}
void getSplitEventAndRisk(uint    treeID,
                          uint   *repMembrIndx,
                          uint    repMembrSize,
                          uint   *nonMissMembrIndx,
                          uint    nonMissMembrSize,
                          uint   *localEventTimeCount,
                          uint   *localEventTimeIndex,
                          uint    localEventTimeSize,
                          uint   *nodeParentEvent,
                          uint   *nodeParentAtRisk) {
  uint i, j;
  for (i=1; i <= localEventTimeSize; i++) {
    nodeParentAtRisk[i] = 0;
    nodeParentEvent[i] = localEventTimeCount[localEventTimeIndex[i]];
    for (j = 1; j <= nonMissMembrSize; j++) {
      if (localEventTimeIndex[i] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[j]] ]) {
        nodeParentAtRisk[i] ++;
      }
    }
  }
}
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
                          uint  **nodeRightAtRisk) {
  stackSplitEventTime(localEventTimeCount, localEventTimeIndex);
  *localEventTimeSize = getSplitEventTime( treeID,
                                           repMembrIndx,
                                           repMembrSize,
                                           nonMissMembrIndx,
                                           nonMissMembrSize,
                                          *localEventTimeCount,
                                          *localEventTimeIndex);
  stackSplitEventAndRisk(*localEventTimeSize,
                          nodeParentEvent,
                          nodeParentAtRisk,
                          nodeLeftEvent,
                          nodeLeftAtRisk,
                          nodeRightEvent,
                          nodeRightAtRisk);
  getSplitEventAndRisk( treeID,
                        repMembrIndx,
                        repMembrSize,
                        nonMissMembrIndx,
                        nonMissMembrSize,
                       *localEventTimeCount,
                       *localEventTimeIndex,
                       *localEventTimeSize,
                       *nodeParentEvent,
                       *nodeParentAtRisk);
}
void stackAndGetSplitSurvL2(uint     treeID,
                            Node    *parent,
                            uint     localEventTimeSize,
                            uint    *localEventTimeIndex,
                            uint    *nodeParentEvent,
                            uint    *nodeParentAtRisk,
                            double **localRatio,
                            double **localSurvival) {
  uint q;
  *localRatio = dvector(1, localEventTimeSize + 1);
  *localSurvival = dvector(1, localEventTimeSize + 1);
  for (q = 1; q <= localEventTimeSize; q++) {
    if (nodeParentEvent[q] > 0) {
      if (nodeParentAtRisk[q] >= 1) {
        (*localRatio)[q] = ((double) nodeParentEvent[q] / nodeParentAtRisk[q]);
      }
      else {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Zero At Risk Count encountered in local ratio calculation for (tree, leaf) = (%10d, %10d)", treeID, parent -> nodeID);
        RFprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
    else {
      (*localRatio)[q] = 0.0;
    }
    (*localSurvival)[q] = 1.0 - (*localRatio)[q];
  }
  for (q = 2; q <= localEventTimeSize; q++) {
    (*localSurvival)[q] *= (*localSurvival)[q-1];
  }
}
void stackAndGetL2Impute(uint     treeID,
                         Node    *parent,
                         uint    *repMembrIndx,
                         uint     repMembrSize,
                         uint    *nonMissMembrIndx,
                         uint     nonMissMembrSize,
                         uint     localEventTimeSize,
                         double  *localSurvival,
                         double **l2Impute) {
}
void unstackSplitSurv(uint *localEventTimeCount,
                      uint *localEventTimeIndex,
                      uint  eventTimeSize,
                      uint *nodeParentEvent,
                      uint *nodeParentAtRisk,
                      uint *nodeLeftEvent,
                      uint *nodeLeftAtRisk,
                      uint *nodeRightEvent,
                      uint *nodeRightAtRisk) {
  unstackSplitEventTime(localEventTimeCount,
                        localEventTimeIndex);
  unstackSplitEventAndRisk(eventTimeSize,
                           nodeParentEvent,
                           nodeParentAtRisk,
                           nodeLeftEvent,
                           nodeLeftAtRisk,
                           nodeRightEvent,
                           nodeRightAtRisk);
}
void unstackAndGetSplitSurvL2(uint     localEventTimeSize,
                              double  *localRatio,
                              double  *localSurvival) {
  free_dvector(localRatio, 1, localEventTimeSize + 1);
  free_dvector(localSurvival, 1, localEventTimeSize + 1);
}
uint stackAndConstructSplitVector (uint     treeID,
                                   uint     repMembrSize,
                                   uint     randomCovariateIndex,
                                   double  *splitVector,
                                   uint     splitVectorSize,
                                   char    *factorFlag,
                                   char    *deterministicSplitFlag,
                                   uint    *mwcpSizeAbsolute,
                                   void   **splitVectorPtr) {
  uint  sworIndex;
  uint *sworVector;
  uint  sworVectorSize;
  uint j, j2, k2;
  uint factorSizeAbsolute;
  uint offset;
  uint splitLength;
  uint relativePair;
  splitLength = 0;  
  (*splitVectorPtr) = NULL;  
  if (strcmp(RF_xType[randomCovariateIndex], "C") == 0) {
    *factorFlag = TRUE;
    if(RF_factorList[treeID][splitVectorSize] == NULL) {
      RF_factorList[treeID][splitVectorSize] = makeFactor(splitVectorSize, FALSE);
    }
    factorSizeAbsolute = RF_xFactorSize[RF_xFactorMap[randomCovariateIndex]];
    *mwcpSizeAbsolute = RF_factorList[treeID][factorSizeAbsolute] -> mwcpSize;
    if (RF_splitRule == RAND_SPLIT) {
      splitLength = 1 + 1;
      *deterministicSplitFlag = FALSE;
    }
    else {
      if(RF_splitRandomCount == 0) {
        *deterministicSplitFlag = TRUE;
        if ((RF_factorList[treeID][splitVectorSize] -> r) > MAX_EXACT_LEVEL) {
          *deterministicSplitFlag = FALSE;
        }
        else {
          if ( *((uint *) RF_factorList[treeID][splitVectorSize] -> complementaryPairCount) >= repMembrSize ) {
            *deterministicSplitFlag = FALSE;
          }
        }
        if (*deterministicSplitFlag == FALSE) {
          splitLength = repMembrSize + 1;
        }
        else {
          splitLength = *((uint*) RF_factorList[treeID][splitVectorSize] -> complementaryPairCount) + 1;
        }
      }
      else {
        *deterministicSplitFlag = FALSE;
        if ((RF_factorList[treeID][splitVectorSize] -> r) <= MAX_EXACT_LEVEL) {
          if (*((uint*) RF_factorList[treeID][splitVectorSize] -> complementaryPairCount) <= ((RF_splitRandomCount <= repMembrSize) ? RF_splitRandomCount : repMembrSize)) {
            splitLength = *((uint*) RF_factorList[treeID][splitVectorSize] -> complementaryPairCount) + 1;
            *deterministicSplitFlag = TRUE;
          }
        }
        if (*deterministicSplitFlag == FALSE) {
          splitLength = 1 + ((RF_splitRandomCount <= repMembrSize) ? RF_splitRandomCount : repMembrSize);
        }
      }  
    }  
    (*splitVectorPtr) = uivector(1, splitLength * (*mwcpSizeAbsolute));
    for (offset = 1; offset <= *mwcpSizeAbsolute; offset++) {
      ((uint*) (*splitVectorPtr) + ((splitLength - 1) * (*mwcpSizeAbsolute)))[offset] = 0;
    }
    if (*deterministicSplitFlag) {
      bookFactor(RF_factorList[treeID][splitVectorSize]);
      j2 = 0;
      for (j = 1; j <= RF_factorList[treeID][splitVectorSize] -> cardinalGroupCount; j++) {
        for (k2 = 1; k2 <= ((uint*) RF_factorList[treeID][splitVectorSize] -> cardinalGroupSize)[j]; k2++) {
          ++j2;
          relativePair = (RF_factorList[treeID][splitVectorSize] -> cardinalGroupBinary)[j][k2];
          convertRelToAbsBinaryPair(treeID,
                                    splitVectorSize,
                                    factorSizeAbsolute,
                                    relativePair,
                                    splitVector,
                                    (uint*) (*splitVectorPtr) + ((j2 - 1) * (*mwcpSizeAbsolute)));
        }
      }
    }  
    else {
      for (j = 1; j < splitLength; j++) {
        getRandomPair(treeID, splitVectorSize, factorSizeAbsolute, splitVector, (uint*) (*splitVectorPtr) + ((j - 1) * (*mwcpSizeAbsolute)));
      }
    }
  }  
  else {
    *factorFlag = FALSE;
    if (RF_splitRule == RAND_SPLIT) {
      splitLength = 1 + 1;
      *deterministicSplitFlag = FALSE;
    }
    else {
      if(RF_splitRandomCount == 0) {
        splitLength = splitVectorSize;
        (*splitVectorPtr) = splitVector;
        *deterministicSplitFlag = TRUE;
      }
      else {
        if (splitVectorSize <= RF_splitRandomCount) {
          splitLength = splitVectorSize;
          (*splitVectorPtr) = splitVector;
          *deterministicSplitFlag = TRUE;
        }
        else {
          splitLength = RF_splitRandomCount + 1;
          *deterministicSplitFlag = FALSE;
        }
      }  
    }  
    if (*deterministicSplitFlag == FALSE) {
      (*splitVectorPtr) = dvector(1, splitLength);
      ((double*) (*splitVectorPtr))[splitLength] = 0;
      if (RF_splitRule == RAND_SPLIT) {
        ((double*) (*splitVectorPtr))[1]  = splitVector[(uint) ceil(ran1B(treeID) * ((splitVectorSize - 1) * 1.0))];
      }
      else {
        sworVector = uivector(1, splitVectorSize);
        sworVectorSize = splitVectorSize - 1;
        for (j = 1; j <= sworVectorSize; j++) {
          sworVector[j] = j;
        }
        for (j = 1; j < splitLength; j++) {
          sworIndex = (uint) ceil(ran1B(treeID) * (sworVectorSize * 1.0));
          ((double*) (*splitVectorPtr))[j]  = splitVector[sworVector[sworIndex]];
          sworVector[sworIndex] = sworVector[sworVectorSize];
          sworVectorSize --;
        }
        free_uivector (sworVector, 1, splitVectorSize);
        sort(((double*) (*splitVectorPtr)), splitLength-1);
      }
    }
  }  
  return splitLength;
}
void unstackSplitVector(uint   treeID,
                        uint   splitVectorSize,
                        uint   splitLength,
                        char   factorFlag,
                        char   deterministicSplitFlag,
                        uint   mwcpSizeAbsolute,
                        void  *splitVectorPtr) {
  if (factorFlag == TRUE) {
    free_uivector(splitVectorPtr, 1, splitLength * mwcpSizeAbsolute);
    if (deterministicSplitFlag == FALSE) {
      if (splitVectorSize > SAFE_FACTOR_SIZE) {
        unbookFactor(RF_factorList[treeID][splitVectorSize]);
      }
    }
  }
  else {
    if (deterministicSplitFlag == FALSE) {
      free_dvector(splitVectorPtr, 1, splitLength);
    }
  }
}
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
                           uint   ***densitySwap) {
  initializeCDF(treeID,
                (RF_xWeightType == RF_WGHT_UNIFORM) ? NULL : NULL, 
                (RF_xWeightType == RF_WGHT_UNIFORM) ? (parent -> permissibleSplit) : (parent -> permissibleSplit),
                (RF_xWeightType == RF_WGHT_UNIFORM) ? parent -> xSize : parent -> xSize,
                RF_xWeightType,
                RF_xWeight,
                RF_xWeightSorted,
                RF_xWeightDensitySize,
                covariateIndex,
                covariateSize,
                cdf,
                cdfSize,
                cdfSort,
                density,
                densitySize,
                densitySwap);
}
void unstackRandomCovariates(uint     treeID,
                             Node    *parent,
                             uint    *covariateIndex,
                             uint     covariateSize,
                             double  *cdf,
                             uint     cdfSize,
                             uint    *cdfSort,
                             uint    *density,
                             uint     densitySize,
                             uint   **densitySwap,
                             uint     repMembrSize) {
  discardCDF(treeID,
             (RF_xWeightType == RF_WGHT_UNIFORM) ? parent -> xSize : parent -> xSize,
             RF_xWeightType,
             RF_xWeight,
             RF_xWeightDensitySize,
             covariateIndex,
             density,
             densitySwap,
             cdf,
             cdfSort);
}
char selectRandomCovariates(uint     treeID,
                            Node     *parent,
                            uint     *repMembrIndx,
                            uint      repMembrSize,
                            uint     *covariateIndex,
                            uint     *uniformCovariateSize,
                            uint     *uniformCovariateIndex,
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
                            char      multImpFlag) {
  uint i, ii;
  uint candidateCovariate;
  uint offset;
  double *nonMissSplit;
  char mPredictorFlag;
  char splittable;
  if (nonMissMembrSizeStatic < 1) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Invalid nonMissMembrSizeStatic encountered in selectRandomCovariates():  %10d", nonMissMembrSizeStatic);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  nonMissSplit = dvector(1, repMembrSize);
  (*covariate) = candidateCovariate = -1;
  splittable = FALSE;
  (*indxx) = uivector(1, repMembrSize);
  if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
    *nonMissMembrSize = nonMissMembrSizeStatic;
    *nonMissMembrIndx = nonMissMembrIndxStatic;
  }
  else {
    *nonMissMembrSize = 0;
    *nonMissMembrIndx = uivector(1, nonMissMembrSizeStatic);
  }
  while ( ((*candidateCovariateCount) < RF_randomCovariateCount) &&
          (candidateCovariate != 0) && (splittable == FALSE)) {
    candidateCovariate = sampleFromCDF(ran1B,
                                       treeID,
                                       RF_xWeightType,
                                       covariateIndex,
                                       *uniformCovariateSize,
                                       uniformCovariateIndex,
                                       cdf,
                                       *cdfSize,
                                       cdfSort,
                                       density,
                                       *densitySize);
    if (candidateCovariate != 0) {
      updateCDF(treeID,
                RF_xWeightType,
                RF_xWeight,
                covariateIndex,
                uniformCovariateSize,
                *uniformCovariateIndex,
                cdf,
                cdfSize,
                density,
                densitySize,
                densitySwap,
                candidateCovariate);
      (*actualCovariateCount) ++;
      (*candidateCovariateCount) ++;
      splittable = TRUE;
      if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
        for (i = 1; i <= repMembrSize; i++) {
          nonMissSplit[i] = RF_observation[treeID][candidateCovariate][repMembrIndx[i]];
        }
        *nonMissMembrSize = nonMissMembrSizeStatic;
        *nonMissMembrIndx = nonMissMembrIndxStatic;
      }
      else {
        offset = RF_rSize + candidateCovariate;
        (*nonMissMembrSize) = 0;
        for (i = 1; i <= nonMissMembrSizeStatic; i++) {
          ii = nonMissMembrIndxStatic[i];
          mPredictorFlag = FALSE;
          if (RF_mRecordMap[repMembrIndx[ii]] > 0) {
            if (RF_mpSign[offset][RF_mRecordMap[repMembrIndx[ii]]] == 1) {
                mPredictorFlag = TRUE;
            }
          }
          if (!mPredictorFlag) {
            (*nonMissMembrSize) ++;
            (*nonMissMembrIndx)[*nonMissMembrSize] = ii;
            nonMissSplit[*nonMissMembrSize] = RF_observation[treeID][candidateCovariate][repMembrIndx[(*nonMissMembrIndx)[*nonMissMembrSize]]];
          }
        }  
      }  
      if ((*nonMissMembrSize) == 0) {
        splittable = FALSE;
      }
      if (splittable) {
        indexx((*nonMissMembrSize),
               nonMissSplit,
               (*indxx));
        splitVector[1] = nonMissSplit[(*indxx)[1]];
        (*splitVectorSize) = 1;
        for (i = 2; i <= (*nonMissMembrSize); i++) {
          if (nonMissSplit[(*indxx)[i]] > splitVector[(*splitVectorSize)]) {
            (*splitVectorSize) ++;
            splitVector[(*splitVectorSize)] = nonMissSplit[(*indxx)[i]];
          }
        }
        if((*splitVectorSize) >= 2) {
          (*covariate) = candidateCovariate;
        }
        else {
          splittable = FALSE;
        }
      }  
      if (!splittable) {
        (parent -> permissibleSplit)[candidateCovariate] = FALSE;
        (*actualCovariateCount) --;
      }
    }  
    else {
    }
  }  
  if (!splittable) {
    free_uivector(*indxx, 1, repMembrSize);
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      *nonMissMembrSize = 0;
      *nonMissMembrIndx = NULL;
    }
    else {
      *nonMissMembrSize = 0;
      free_uivector(*nonMissMembrIndx, 1, nonMissMembrSizeStatic);
    }
  }
  free_dvector(nonMissSplit, 1, repMembrSize);
  return splittable;
}
void unselectRandomCovariates(uint      treeID,
                              Node     *parent,
                              uint      repMembrSize,
                              uint     *indxx,
                              uint     nonMissMembrSizeStatic,
                              uint    *nonMissMembrIndx,
                              char      multImpFlag) {
  free_uivector((indxx), 1, repMembrSize);
  if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
    free_uivector(nonMissMembrIndx, 1, nonMissMembrSizeStatic);
  }
}
char selectRandomCovariatesNew(uint     treeID,
                            Node     *parent,
                            uint     *repMembrIndx,
                            uint      repMembrSize,
                            uint     *covariateIndex,
                            uint     *uniformCovariateSize,
                            uint     *uniformCovariateIndex,
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
                            uint     *missMembrSize,
                            uint    **missMembrIndx,
                            char      multImpFlag) {
  uint i, ii;
  uint candidateCovariate;
  uint offset;
  double *nonMissSplit;
  char mPredictorFlag;
  char splittable;
  if (nonMissMembrSizeStatic < 1) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Invalid nonMissMembrSizeStatic encountered in selectRandomCovariates():  %10d", nonMissMembrSizeStatic);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  nonMissSplit = dvector(1, repMembrSize);
  (*covariate) = candidateCovariate = -1;
  splittable = FALSE;
  (*indxx) = uivector(1, repMembrSize);
  *nonMissMembrSize = 0;
  *missMembrSize = 0;
  *missMembrIndx = NULL;
  if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
    *nonMissMembrSize = nonMissMembrSizeStatic;
    *nonMissMembrIndx = nonMissMembrIndxStatic;
  }
  else {
    *nonMissMembrIndx = uivector(1, nonMissMembrSizeStatic);
    if (RF_optHigh & (OPT_MISS_MIA | OPT_MISS_MIAH)) {
      *missMembrIndx = uivector(1, nonMissMembrSizeStatic);
    }
  }
  while ( ((*candidateCovariateCount) < RF_randomCovariateCount) &&
          (candidateCovariate != 0) && (splittable == FALSE)) {
    candidateCovariate = sampleFromCDF(ran1B,
                                       treeID,
                                       RF_xWeightType,
                                       covariateIndex,
                                       *uniformCovariateSize,
                                       uniformCovariateIndex,
                                       cdf,
                                       *cdfSize,
                                       cdfSort,
                                       density,
                                       *densitySize);
    if (candidateCovariate != 0) {
      updateCDF(treeID,
                RF_xWeightType,
                RF_xWeight,
                covariateIndex,
                uniformCovariateSize,
                *uniformCovariateIndex,
                cdf,
                cdfSize,
                density,
                densitySize,
                densitySwap,
                candidateCovariate);
      (*actualCovariateCount) ++;
      (*candidateCovariateCount) ++;
      splittable = TRUE;
      if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
        for (i = 1; i <= repMembrSize; i++) {
          nonMissSplit[i] = RF_observation[treeID][candidateCovariate][repMembrIndx[i]];
        }
        *nonMissMembrSize = nonMissMembrSizeStatic;
        *nonMissMembrIndx = nonMissMembrIndxStatic;
      }
      else {
        offset = RF_rSize + candidateCovariate;
        (*nonMissMembrSize) = 0;
        (*missMembrSize) = 0;
        for (i = 1; i <= nonMissMembrSizeStatic; i++) {
          ii = nonMissMembrIndxStatic[i];
          mPredictorFlag = FALSE;
          if (RF_mRecordMap[repMembrIndx[ii]] > 0) {
            if (RF_mpSign[offset][RF_mRecordMap[repMembrIndx[ii]]] == 1) {
                mPredictorFlag = TRUE;
            }
          }
          if (!mPredictorFlag) {
            (*nonMissMembrSize) ++;
            (*nonMissMembrIndx)[*nonMissMembrSize] = ii;
            nonMissSplit[*nonMissMembrSize] = RF_observation[treeID][candidateCovariate][repMembrIndx[(*nonMissMembrIndx)[*nonMissMembrSize]]];
          }
          else {
            if (RF_optHigh & (OPT_MISS_MIA | OPT_MISS_MIAH)) {
              (*missMembrSize) ++;
              (*missMembrIndx)[*missMembrSize] = ii;
            }
          }
        }  
      }  
      if ((*nonMissMembrSize) == 0) {
        splittable = FALSE;
      }
      if (splittable) {
        indexx((*nonMissMembrSize),
               nonMissSplit,
               (*indxx));
        splitVector[1] = nonMissSplit[(*indxx)[1]];
        (*splitVectorSize) = 1;
        for (i = 2; i <= (*nonMissMembrSize); i++) {
          if (nonMissSplit[(*indxx)[i]] > splitVector[(*splitVectorSize)]) {
            (*splitVectorSize) ++;
            splitVector[(*splitVectorSize)] = nonMissSplit[(*indxx)[i]];
          }
        }
        if((*splitVectorSize) >= 2) {
          (*covariate) = candidateCovariate;
        }
        else {
          splittable = FALSE;
        }
      }  
      if (!splittable) {
        (parent -> permissibleSplit)[candidateCovariate] = FALSE;
        (*actualCovariateCount) --;
      }
    }  
    else {
    }
  }  
  if (!splittable) {
    free_uivector(*indxx, 1, repMembrSize);
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      *nonMissMembrSize = 0;
      *nonMissMembrIndx = NULL;
    }
    else {
      *nonMissMembrSize = 0;
      free_uivector(*nonMissMembrIndx, 1, nonMissMembrSizeStatic);
      if (RF_optHigh & (OPT_MISS_MIA | OPT_MISS_MIAH)) {
        *missMembrSize = 0;
        free_uivector(*missMembrIndx, 1, nonMissMembrSizeStatic);
      }
    }
  }
  free_dvector(nonMissSplit, 1, repMembrSize);
  return splittable;
}
void unselectRandomCovariatesNew(uint      treeID,
                                 Node     *parent,
                                 uint      repMembrSize,
                                 uint     *indxx,
                                 uint     nonMissMembrSizeStatic,
                                 uint    *nonMissMembrIndx,
                                 uint    *missMembrIndx,
                                 char      multImpFlag) {
  free_uivector((indxx), 1, repMembrSize);
  if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
    free_uivector(nonMissMembrIndx, 1, nonMissMembrSizeStatic);
    if (RF_optHigh & (OPT_MISS_MIA | OPT_MISS_MIAH)) {
      free_uivector(missMembrIndx, 1, nonMissMembrSizeStatic);
    }
  }
}
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
                           uint *currentMembrIter) {
  char daughterFlag;
  char iterFlag;
  iterFlag = TRUE;
  *currentMembrIter = priorMembrIter;
  while (iterFlag) {
    (*currentMembrIter) ++;
    if (factorFlag == TRUE) {
      daughterFlag = splitOnFactor((uint)  RF_observation[treeID][randomCovariate][    repMembrIndx[nonMissMembrIndx[indxx[*currentMembrIter]]]     ],
                                   (uint*) splitVectorPtr + ((offset - 1) * mwcpSizeAbsolute));
      if ((*currentMembrIter) == nonMissMembrSize) {
        iterFlag = FALSE;
      }
    }
    else {
      if ((((double*) splitVectorPtr)[offset] - RF_observation[treeID][randomCovariate][   repMembrIndx[nonMissMembrIndx[indxx[*currentMembrIter]]]    ]) >= 0.0) {
        daughterFlag = LEFT;
      }
      else {
        daughterFlag = RIGHT;
        iterFlag = FALSE;
      }
    }
    localSplitIndicator[     nonMissMembrIndx[indxx[*currentMembrIter]]   ] = daughterFlag;
    if (daughterFlag == LEFT) {
      (*leftSize) ++;
    }  
    else {
    }
  }  
  return (*leftSize);
 }
uint virtuallySplitNodeNew(uint  treeID,
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
                           uint *currentMembrIter) {
  char daughterFlag;
  char iterFlag;
  iterFlag = TRUE;
  *currentMembrIter = priorMembrIter;
  while (iterFlag) {
    (*currentMembrIter) ++;
    if (factorFlag == TRUE) {
      daughterFlag = splitOnFactor((uint)  RF_observation[treeID][randomCovariate][    repMembrIndx[nonMissMembrIndx[indxx[*currentMembrIter]]]     ],
                                   (uint*) splitVectorPtr + ((offset - 1) * mwcpSizeAbsolute));
      if ((*currentMembrIter) == nonMissMembrSize) {
        iterFlag = FALSE;
      }
    }
    else {
      if ((((double*) splitVectorPtr)[offset] - RF_observation[treeID][randomCovariate][   repMembrIndx[nonMissMembrIndx[indxx[*currentMembrIter]]]    ]) >= 0.0) {
        daughterFlag = LEFT;
      }
      else {
        daughterFlag = RIGHT;
        iterFlag = FALSE;
      }
    }
    localSplitIndicator[     nonMissMembrIndx[indxx[*currentMembrIter]]   ] = daughterFlag;
    if (daughterFlag == LEFT) {
      (*leftSize) ++;
    }  
    else {
    }
  }  
  return (*leftSize);
 }
void getReweightedRandomPair (uint    treeID,
                              uint    relativeFactorSize,
                              uint    absoluteFactorSize,
                              double *absoluteLevel,
                              uint   *result) {
  uint randomGroupIndex;
  if(RF_factorList[treeID][relativeFactorSize] == NULL) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Factor not allocated for size:  %10d", relativeFactorSize);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  randomGroupIndex = (uint) ceil(ran1B(treeID) * ((RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount) * 1.0));
  createRandomBinaryPair(treeID, relativeFactorSize, absoluteFactorSize, randomGroupIndex, absoluteLevel, result);
}
void getRandomPair (uint treeID, uint relativeFactorSize, uint absoluteFactorSize, double *absoluteLevel, uint *result) {
  uint randomGroupIndex;
  double randomValue;
  uint k;
  if(RF_factorList[treeID][relativeFactorSize] == NULL) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Factor not allocated for size:  %10d", relativeFactorSize);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  double *cdf = dvector(1, RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount);
  if (relativeFactorSize <= MAX_EXACT_LEVEL) {
    for (k=1; k <= RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount; k++) {
      cdf[k] = (double) ((uint*) RF_factorList[treeID][relativeFactorSize] -> cardinalGroupSize)[k];
    }
  }
  else {
    for (k=1; k <= RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount; k++) {
      cdf[k] = ((double*) RF_factorList[treeID][relativeFactorSize] -> cardinalGroupSize)[k];
    }
  }
  for (k=2; k <= RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount; k++) {
    cdf[k] += cdf[k-1];
  }
  randomValue = ceil((ran1B(treeID) * cdf[RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount]));
  randomGroupIndex = 1;
  while (randomValue > cdf[randomGroupIndex]) {
    randomGroupIndex ++;
  }
  free_dvector(cdf, 1, RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount);
  createRandomBinaryPair(treeID, relativeFactorSize, absoluteFactorSize, randomGroupIndex, absoluteLevel, result);
}
void createRandomBinaryPair(uint    treeID,
                            uint    relativeFactorSize,
                            uint    absoluteFactorSize,
                            uint    groupIndex,
                            double *absoluteLevel,
                            uint   *pair) {
  uint mwcpLevelIdentifier;
  uint mwcpSizeAbsolute;
  uint offset, levelSize, levelIndex;
  uint k;
  levelIndex = 0;  
  mwcpSizeAbsolute = RF_factorList[treeID][absoluteFactorSize] -> mwcpSize;
  uint *levelVector = uivector(1, relativeFactorSize);
  uint *randomLevel = uivector(1, groupIndex);
  for (k = 1; k <= relativeFactorSize; k++) {
    levelVector[k] = k;
  }
  levelSize = relativeFactorSize;
  for (k = 1; k <= groupIndex; k++) {
    randomLevel[k] = sampleUniformlyFromVector(treeID,
                                               levelVector,
                                               levelSize,
                                               &levelIndex);
    levelVector[levelIndex] = levelVector[levelSize];
    levelSize --;
  }
  for (k = 1; k <= groupIndex; k++) {
    randomLevel[k] = (uint) absoluteLevel[randomLevel[k]];
  }
  for (offset = 1; offset <= mwcpSizeAbsolute; offset++) {
    pair[offset] = 0;
  }
  for (k = 1; k <= groupIndex; k++) {
    mwcpLevelIdentifier = (randomLevel[k] >> (3 + ulog2(sizeof(uint)))) + ((randomLevel[k] & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
    pair[mwcpLevelIdentifier] += upower(2, randomLevel[k] - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1 );
  }
  free_uivector(levelVector, 1, relativeFactorSize);
  free_uivector(randomLevel, 1, groupIndex);
}
void convertRelToAbsBinaryPair(uint    treeID,
                               uint    relativeFactorSize,
                               uint    absoluteFactorSize,
                               uint    relativePair,
                               double *absoluteLevel,
                               uint   *pair) {
  uint mwcpLevelIdentifier;
  uint mwcpSizeAbsolute;
  uint coercedAbsoluteLevel;
  uint k, offset;
  mwcpSizeAbsolute = RF_factorList[treeID][absoluteFactorSize] -> mwcpSize;
  for (offset = 1; offset <= mwcpSizeAbsolute; offset++) {
    pair[offset] = 0;
  }
  for (k = 1; k <= relativeFactorSize; k++) {
    if (relativePair & ((uint) 0x01)) {
      coercedAbsoluteLevel = (uint) absoluteLevel[k];
      mwcpLevelIdentifier = (coercedAbsoluteLevel >> (3 + ulog2(sizeof(uint)))) + ((coercedAbsoluteLevel & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
      pair[mwcpLevelIdentifier] += upower(2, coercedAbsoluteLevel - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1 );
    }
    relativePair = relativePair >> 1;
  }
}
char summarizeSplitResult(uint    splitParameterMax,
                          double  splitValueMaxCont,
                          uint    splitValueMaxFactSize,
                          uint   *splitValueMaxFactPtr,
                          double *splitStatistic,
                          double  deltaMax) {
  char result;
  if (splitParameterMax > 0) {
    *splitStatistic = deltaMax;
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  return result;
}
char getPreSplitResult (uint      treeID,
                        Node     *parent,
                        uint      repMembrSize,
                        uint     *repMembrIndx,
                        uint     *nonMissMembrSize,
                        uint    **nonMissMembrIndx,
                        double  **splitVector,
                        double   *preSplitMean,
                        char      multImpFlag,
                        char      multVarFlag) {
  uint i, r;
  char mResponseFlag;
  char result;
  if (repMembrSize >= (2 * RF_minimumNodeSize)) {
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  if (result) {
    if (RF_maximumNodeDepth < 0) {
      result = TRUE;
    }
    else {
      if (parent -> depth < (uint) RF_maximumNodeDepth) {
        result = TRUE;
      }
      else {
        result = FALSE;
      }
    }
  }
  if (result) {
    *splitVector = dvector(1, repMembrSize);
    if ((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP) || multVarFlag) {
      (*nonMissMembrSize) = repMembrSize;
      (*nonMissMembrIndx) = RF_identityMembershipIndex;
    }
    else {
      *nonMissMembrIndx = uivector(1, repMembrSize);
      (*nonMissMembrSize) = 0;
      for (i = 1; i <= repMembrSize; i++) {
        mResponseFlag = FALSE;
        if (RF_mRecordMap[repMembrIndx[i]] > 0) {
          for (r = 1; r <= RF_rSize; r++) {
            if (RF_mpSign[r][RF_mRecordMap[repMembrIndx[i]]] == 1) {
              mResponseFlag = TRUE;
            }
          }
        }
        if (!mResponseFlag) {
          (*nonMissMembrSize) ++;
          (*nonMissMembrIndx)[(*nonMissMembrSize)] = i;
        }
      }  
    }  
    if (!multVarFlag) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        uint q,k,m;
        uint *evntProp = uivector(1, RF_eventTypeSize + 1);
        for (q=1; q <= RF_eventTypeSize + 1; q++) {
          evntProp[q] = 0;
        }
        for (i = 1; i <= (*nonMissMembrSize); i++) {
          m = (uint) RF_status[treeID][repMembrIndx[(*nonMissMembrIndx)[i]]];
          if (m > 0) {
            evntProp[RF_eventTypeIndex[m]] ++;
          }
          else {
            evntProp[RF_eventTypeSize + 1] ++;
          }
        }
        k = 0;
        for (q = 1; q <= RF_eventTypeSize + 1; q++) {
          if(evntProp[q] > 0) {
            k ++;
          }
        }
        if (k == 0) {
          result = FALSE;
        }
        else {
          if (k == 1) {
            if (evntProp[RF_eventTypeSize + 1] > 0) {
              result = FALSE;
            }
            else {
              result = getVariance(repMembrSize,
                                   repMembrIndx,
                                   *nonMissMembrSize,
                                   *nonMissMembrIndx,
                                   RF_time[treeID],
                                   preSplitMean,
                                   NULL);
            }
          }
        }
        free_uivector(evntProp, 1, RF_eventTypeSize + 1);
      }
      else {
        result = getVariance(repMembrSize,
                             repMembrIndx,
                             *nonMissMembrSize,
                             *nonMissMembrIndx,
                             RF_response[treeID][1],
                             preSplitMean,
                             NULL);
      }
    }
    if (!result) {
      free_dvector(*splitVector, 1, repMembrSize);
      (*nonMissMembrSize) = 0;
      if (!((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP) || multVarFlag)) {
        free_uivector(*nonMissMembrIndx, 1, repMembrSize);
      }
    }
  }
  return result;
}
void unstackPreSplit (char      preliminaryResult,
                      char      multImpFlag,
                      char      multVarFlag,
                      uint      repMembrSize,
                      double   *splitVector,
                      uint     *nonMissMembrIndx) {
  if (preliminaryResult) {
    free_dvector(splitVector, 1, repMembrSize);
    if (!((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP) || multVarFlag)) {
      free_uivector(nonMissMembrIndx, 1, repMembrSize);
    }
  }
  else {
  }
}
char getPreSplitResultNew (uint      treeID,
                           Node     *parent,
                           uint      repMembrSize,
                           uint     *repMembrIndx,
                           uint     *nonMissMembrSize,
                           uint    **nonMissMembrIndx,
                           uint     *miaVectorLength,
                           char    **miaVectorType,
                           double  **splitVector,
                           double   *preSplitMean,
                           char      multImpFlag,
                           char      multVarFlag) {
  uint i, r;
  char mResponseFlag;
  char result;
  if (repMembrSize >= (2 * RF_minimumNodeSize)) {
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  if (result) {
    if (RF_maximumNodeDepth < 0) {
      result = TRUE;
    }
    else {
      if (parent -> depth < (uint) RF_maximumNodeDepth) {
        result = TRUE;
      }
      else {
        result = FALSE;
      }
    }
  }
  if (result) {
    *splitVector = dvector(1, repMembrSize);
    if ((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP) || multVarFlag) {
      (*nonMissMembrSize) = repMembrSize;
      (*nonMissMembrIndx) = RF_identityMembershipIndex;
      *miaVectorType = cvector(1, 1);
      (*miaVectorLength) = 1;
      (*miaVectorType)[(*miaVectorLength)] = SPLIT_MIA_NONE;
    }
    else {
      *nonMissMembrIndx = uivector(1, repMembrSize);
      (*nonMissMembrSize) = 0;
      for (i = 1; i <= repMembrSize; i++) {
        mResponseFlag = FALSE;
        if (RF_mRecordMap[repMembrIndx[i]] > 0) {
          for (r = 1; r <= RF_rSize; r++) {
            if (RF_mpSign[r][RF_mRecordMap[repMembrIndx[i]]] == 1) {
              mResponseFlag = TRUE;
            }
          }
        }
        if (!mResponseFlag) {
          (*nonMissMembrSize) ++;
          (*nonMissMembrIndx)[(*nonMissMembrSize)] = i;
        }
      }  
      *miaVectorType = cvector(1, 3);
      (*miaVectorLength) = 0;
      if (RF_optHigh & (OPT_MISS_MIA | OPT_MISS_MIAH)) {
        if (RF_optHigh & OPT_MISS_MIAH) {
          (*miaVectorType)[++ (*miaVectorLength)] = SPLIT_MIA_NONE;
        }
        (*miaVectorType)[++ (*miaVectorLength)] = SPLIT_MIA_LEFT;
        (*miaVectorType)[++ (*miaVectorLength)] = SPLIT_MIA_RGHT;
      }
      else {
        (*miaVectorType)[++ (*miaVectorLength)] = SPLIT_MIA_NONE;
      }
    }  
    if (!multVarFlag) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        uint q,k,m;
        uint *evntProp = uivector(1, RF_eventTypeSize + 1);
        for (q=1; q <= RF_eventTypeSize + 1; q++) {
          evntProp[q] = 0;
        }
        for (i = 1; i <= (*nonMissMembrSize); i++) {
          m = (uint) RF_status[treeID][repMembrIndx[(*nonMissMembrIndx)[i]]];
          if (m > 0) {
            evntProp[RF_eventTypeIndex[m]] ++;
          }
          else {
            evntProp[RF_eventTypeSize + 1] ++;
          }
        }
        k = 0;
        for (q = 1; q <= RF_eventTypeSize + 1; q++) {
          if(evntProp[q] > 0) {
            k ++;
          }
        }
        if (k == 0) {
          result = FALSE;
        }
        else {
          if (k == 1) {
            if (evntProp[RF_eventTypeSize + 1] > 0) {
              result = FALSE;
            }
            else {
              result = getVariance(repMembrSize,
                                   repMembrIndx,
                                   *nonMissMembrSize,
                                   *nonMissMembrIndx,
                                   RF_time[treeID],
                                   preSplitMean,
                                   NULL);
            }
          }
        }
        free_uivector(evntProp, 1, RF_eventTypeSize + 1);
      }
      else {
        result = getVariance(repMembrSize,
                             repMembrIndx,
                             *nonMissMembrSize,
                             *nonMissMembrIndx,
                             RF_response[treeID][1],
                             preSplitMean,
                             NULL);
      }
    }
    if (!result) {
      free_dvector(*splitVector, 1, repMembrSize);
      (*nonMissMembrSize) = 0;
      if ((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP) || multVarFlag) {
        free_cvector(*miaVectorType, 1, 1);
      }
      else {
        free_uivector(*nonMissMembrIndx, 1, repMembrSize);
        free_cvector(*miaVectorType, 1, 3);
      }
    }
  }
  return result;
}
void unstackPreSplitNew (char      preliminaryResult,
                         char      multImpFlag,
                         char      multVarFlag,
                         uint      repMembrSize,
                         double   *splitVector,
                         uint     *nonMissMembrIndx,
                         char     *miaVectorType) {
  if (preliminaryResult) {
    free_dvector(splitVector, 1, repMembrSize);
    if ((RF_mRecordSize == 0) || multImpFlag || !(RF_optHigh & OPT_MISS_SKIP) || multVarFlag) {
      free_cvector(miaVectorType, 1, 1);
    }
    else {
      free_uivector(nonMissMembrIndx, 1, repMembrSize);
      free_cvector(miaVectorType, 1, 3);
    }
  }
  else {
  }
}
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
                        char  **splitIndicator) {
  char flag;
  uint k;
  if (RF_opt & OPT_NODE_STAT) {
    updateNodeStatistics(parent, delta, candidateCovariateCount, covariate);
  }
  if(ISNA(delta)) {
    flag = FALSE;
  }
  else {
    if(ISNA(*deltaMax)) {
      flag = TRUE;
    }
    else {
      if ((delta - *deltaMax) > EPSILON) {
        flag = TRUE;
      }
      else {
        flag = FALSE;
      }
    }
  }
  if (flag) {
    *deltaMax = delta;
    *splitParameterMax = covariate;
    if (factorFlag == TRUE) {
      if (*splitValueMaxFactSize > 0) {
        if (*splitValueMaxFactSize != mwcpSizeAbsolute) {
          free_uivector(*splitValueMaxFactPtr, 1, *splitValueMaxFactSize);
          *splitValueMaxFactSize = mwcpSizeAbsolute;
          *splitValueMaxFactPtr = uivector(1, *splitValueMaxFactSize);
        }
      }
      else {
        *splitValueMaxFactSize = mwcpSizeAbsolute;
        *splitValueMaxFactPtr = uivector(1, *splitValueMaxFactSize);
      }
      *splitValueMaxCont = NA_REAL;
      for (k=1; k <= *splitValueMaxFactSize; k++) {
        (*splitValueMaxFactPtr)[k] =
          ((uint*) splitVectorPtr + ((index - 1) * (*splitValueMaxFactSize)))[k];
      }
    }
    else {
      if (*splitValueMaxFactSize > 0) {
        free_uivector(*splitValueMaxFactPtr, 1, *splitValueMaxFactSize);
        *splitValueMaxFactSize = 0;
        *splitValueMaxFactPtr = NULL;
      }
      else {
      }
      *splitValueMaxCont = ((double*) splitVectorPtr)[index];
    }
    if (*splitIndicator == NULL) {
     *splitIndicator = cvector(1, repMembrSize);
    }
   for (k=1; k <= repMembrSize; k++) {
     (*splitIndicator)[k] = localSplitIndicator[k];
   }
  }
  else {
  }
  return flag;
}
char updateMaximumSplitNew(uint    treeID,
                           Node   *parent,
                           double  delta,
                           uint    candidateCovariateCount,
                           uint    covariate,
                           uint    index,
                           char    factorFlag,
                           uint    mwcpSizeAbsolute,
                           uint    repMembrSize,
                           char   *localSplitIndicator,
                           char    miaType,
                           double *deltaMax,
                           uint   *splitParameterMax,
                           double *splitValueMaxCont,
                           uint   *splitValueMaxFactSize,
                           uint  **splitValueMaxFactPtr,
                           void   *splitVectorPtr,
                           char  **splitIndicator,
                           char  *splitMIA) {
  char flag;
  uint k;
  if (RF_opt & OPT_NODE_STAT) {
    updateNodeStatistics(parent, delta, candidateCovariateCount, covariate);
  }
  if(ISNA(delta)) {
    flag = FALSE;
  }
  else {
    if(ISNA(*deltaMax)) {
      flag = TRUE;
    }
    else {
      if ((delta - *deltaMax) > EPSILON) {
        flag = TRUE;
      }
      else {
        flag = FALSE;
      }
    }
  }
  if (flag) {
    *deltaMax = delta;
    *splitParameterMax = covariate;
    *splitMIA = miaType;
    if (factorFlag == TRUE) {
      if (*splitValueMaxFactSize > 0) {
        if (*splitValueMaxFactSize != mwcpSizeAbsolute) {
          free_uivector(*splitValueMaxFactPtr, 1, *splitValueMaxFactSize);
          *splitValueMaxFactSize = mwcpSizeAbsolute;
          *splitValueMaxFactPtr = uivector(1, *splitValueMaxFactSize);
        }
      }
      else {
        *splitValueMaxFactSize = mwcpSizeAbsolute;
        *splitValueMaxFactPtr = uivector(1, *splitValueMaxFactSize);
      }
      *splitValueMaxCont = NA_REAL;
      if ((*splitMIA) == SPLIT_MIA_UNIT) {
        for (k=1; k <= *splitValueMaxFactSize; k++) {
          (*splitValueMaxFactPtr)[k] = 0;
        }
      }
      else {
        for (k=1; k <= *splitValueMaxFactSize; k++) {
          (*splitValueMaxFactPtr)[k] =
            ((uint*) splitVectorPtr + ((index - 1) * (*splitValueMaxFactSize)))[k];
        }
      }
    }
    else {
      if (*splitValueMaxFactSize > 0) {
        free_uivector(*splitValueMaxFactPtr, 1, *splitValueMaxFactSize);
        *splitValueMaxFactSize = 0;
        *splitValueMaxFactPtr = NULL;
      }
      else {
      }
      if ((*splitMIA) == SPLIT_MIA_UNIT) {
        *splitValueMaxCont = NA_REAL;
      }
      else {
        *splitValueMaxCont = ((double*) splitVectorPtr)[index];
      }
    }
    if (*splitIndicator == NULL) {
     *splitIndicator = cvector(1, repMembrSize);
    }
   for (k=1; k <= repMembrSize; k++) {
     (*splitIndicator)[k] = localSplitIndicator[k];
   }
  }
  else {
  }
  return flag;
}
void updateNodeStatistics(Node *parent, double delta, uint candidateCovariateCount, uint covariate) {
  char flag;
  if(ISNA(delta)) {
    flag = FALSE;
  }
  else {
    if(ISNA((parent -> mtryStat)[candidateCovariateCount])) {
      flag = TRUE;
    }
    else {
      if ((delta - (parent -> mtryStat)[candidateCovariateCount]) > EPSILON) {
        flag = TRUE;
      }
      else {
        flag = FALSE;
      }
    }
  }
  if (flag) {
    (parent -> mtryIndx)[candidateCovariateCount] = covariate;
    (parent -> mtryStat)[candidateCovariateCount] = delta;
  }
}
void getMeanResponseNew(uint       treeID,
                        Terminal  *parent,
                        uint      *repMembrIndx,
                        uint       repMembrSize,
                        uint      *allMembrIndx,
                        uint       allMembrSize,
                        uint      *rmbrIterator) {
  uint *membershipIndex;
  uint  membershipSize;
  uint i, j;
  if ( (!(RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ||
       ( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2))) { 
    membershipIndex = allMembrIndx;
    membershipSize = parent -> membrCount = allMembrSize;
  }
  else {
    membershipIndex = repMembrIndx;
    membershipSize = parent -> membrCount = repMembrSize;
    if (RF_optHigh & OPT_MEMB_OUTG) {
      RF_TN_RCNT_ptr[treeID][parent -> nodeID] = RF_tTermList[treeID][parent -> nodeID] -> membrCount;
    }
    if (RF_optHigh & OPT_MEMB_INCG) {
      membershipIndex = RF_RMBR_ID_ptr[treeID];
      membershipSize = parent -> membrCount = RF_TN_RCNT_ptr[treeID][parent -> nodeID];
    }
  }
  if ((parent -> membrCount) == 0) {
    if (!(RF_opt & OPT_OUTC_TYPE)) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Zero node count encountered in (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  if (!(RF_optHigh & OPT_TERM_INCG)) {
    stackMeanResponse(parent, RF_rNonFactorCount);
    for (j=1; j <= RF_rNonFactorCount; j++) {
      (parent -> meanResponse)[j] = 0.0;
    }
    if (RF_optHigh & OPT_MEMB_OUTG) {
      for (i=1; i <= membershipSize; i++) {
        RF_RMBR_ID_ptr[treeID][++(*rmbrIterator)] = membershipIndex[i];
        for (j=1; j <= RF_rNonFactorCount; j++) {
          (parent -> meanResponse)[j] += RF_response[treeID][RF_rNonFactorIndex[j]][membershipIndex[i]];
        }
      }
    }
    else if (RF_optHigh & OPT_MEMB_INCG) {
      for (i=1; i <= membershipSize; i++) {
        ++(*rmbrIterator);
        for (j=1; j <= RF_rNonFactorCount; j++) {
          (parent -> meanResponse)[j] += RF_response[treeID][RF_rNonFactorIndex[j]][ membershipIndex[*rmbrIterator] ];
        }
      }
    }
    else {
      for (i=1; i <= membershipSize; i++) {
        for (j=1; j <= RF_rNonFactorCount; j++) {
          (parent -> meanResponse)[j] += RF_response[treeID][RF_rNonFactorIndex[j]][membershipIndex[i]];
        }
      }
    }
    for (j=1; j <= RF_rNonFactorCount; j++) {
      (parent -> meanResponse)[j] = (parent -> meanResponse)[j] / (double) membershipSize;
    }
    if (RF_optHigh & OPT_TERM_OUTG) {
      for (j = 1; j <= RF_rNonFactorCount; j++) {
        RF_TN_REGR_ptr[treeID][parent -> nodeID][j] = (parent -> meanResponse)[j];
      }
    }
  }
  else {
    (parent -> meanResponse) = RF_TN_REGR_ptr[treeID][parent -> nodeID];
  }
}
void updateEnsembleMean(char     mode,
                        uint     treeID,
                        uint     serialTreeID,
                        char     omitDenominator) {
  char oobFlag, fullFlag, selectionFlag, outcomeFlag;
  Terminal ***termMembershipPtr;
  uint    *membershipIndex;
  uint     membershipSize;
  double    **ensembleRGRptr;
  double    **ensembleRGRnum;
  uint       *ensembleDen;
  Terminal *parent;
  uint i, j;
  uint ii;
  ensembleRGRptr = NULL;  
  ensembleRGRnum = NULL;  
  ensembleDen    = NULL;  
  oobFlag = fullFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    termMembershipPtr = RF_ftTermMembership;
    break;
  default:
    if (RF_opt & OPT_OENS) {
      if (RF_oobSize[treeID] > 0) {
        oobFlag = TRUE;
      }
    }
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    termMembershipPtr = RF_tTermMembership;
    break;
  }
  outcomeFlag = TRUE;
  while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
    if (oobFlag == TRUE) {
      ensembleRGRptr = RF_oobEnsembleRGRptr;
      ensembleRGRnum = RF_oobEnsembleRGRnum;
      ensembleDen    = RF_oobEnsembleDen;
      membershipSize  = RF_oobSize[treeID];
      membershipIndex = RF_oobMembershipIndex[treeID];
    }
    else {
      ensembleRGRptr = RF_fullEnsembleRGRptr;
      ensembleRGRnum = RF_fullEnsembleRGRnum;
      ensembleDen    = RF_fullEnsembleDen;
      switch (mode) {
      case RF_PRED:
        membershipSize = RF_fobservationSize;
        membershipIndex = RF_fidentityMembershipIndex;
        break;
      default:
        membershipSize  = RF_observationSize;
        membershipIndex = RF_identityMembershipIndex;
        break;
      }
    }
    for (i = 1; i <= membershipSize; i++) {
      ii = membershipIndex[i];
      parent = termMembershipPtr[treeID][ii];
      selectionFlag = TRUE;
      if (RF_opt & OPT_OUTC_TYPE) {
        if ((parent -> meanResponse) != NULL) {
        }
        else {
          selectionFlag = FALSE;
        }
      }
      if (selectionFlag) {
        for (j=1; j <= RF_rTargetNonFactorCount; j++) {
          ensembleRGRnum[j][ii] += (parent -> meanResponse)[RF_rNonFactorMap[RF_rTargetNonFactor[j]]];
        }
        if(!omitDenominator) {
          ensembleDen[ii] ++;
        }
      }
      if (outcomeFlag) {
        if (ensembleDen[ii] != 0) {
          for (j=1; j <= RF_rTargetNonFactorCount; j++) {
            ensembleRGRptr[j][ii] = ensembleRGRnum[j][ii] / ensembleDen[ii];
          }
        }
      }
    }  
    if (outcomeFlag == TRUE) {
      outcomeFlag = FALSE;
    }
    if (oobFlag == TRUE) {
      oobFlag = FALSE;
    }
    else {
        fullFlag = FALSE;
    }
  }  
}
double getMeanSquareError(uint    size,
                          double *responsePtr,
                          double *predictedOutcome,
                          uint   *denomCount) {
  uint i;
  uint cumDenomCount;
  double result;
  cumDenomCount = 0;
  result = 0.0;
  for (i=1; i <= size; i++) {
    if (denomCount[i] != 0) {
      cumDenomCount += 1;
      result += pow (responsePtr[i] - predictedOutcome[i], 2.0);
    }
  }  
  if (cumDenomCount == 0) {
    result = NA_REAL;
  }
  else {
    result = result / (double) cumDenomCount;
  }
  return result;
}
char getVariance(uint    repMembrSize,
                 uint   *repMembrIndx,
                 uint    nonMissMembrSize,
                 uint   *nonMissMembrIndx,
                 double *targetResponse,
                 double *mean,
                 double *variance) {
  uint i;
  uint denom;
  double meanResult, varResult;
  char result;
  uint *genIndx;
  uint  genSize;
  if (nonMissMembrIndx == NULL) {
    genIndx = RF_identityMembershipIndex;
    genSize = repMembrSize;
  }
  else {
    genIndx = nonMissMembrIndx;
    genSize = nonMissMembrSize;
  }
  denom      = 0;
  meanResult = 0.0;
  for (i = 1; i <= genSize; i++) {
    if(!ISNA(targetResponse[repMembrIndx[genIndx[i]]])) {
      denom ++;
      meanResult += targetResponse[repMembrIndx[genIndx[i]]];
    }
  }
  if (denom > 0) {
    meanResult = meanResult / (double) denom;
  }
  else {
    meanResult = NA_REAL;
  }
  if (mean != NULL) {
    *mean = meanResult;
  }
  varResult = 0.0;
  if(!ISNA(meanResult)) {
    for (i = 1; i <= genSize; i++) {
      if(!ISNA(targetResponse[repMembrIndx[genIndx[i]]])) {
        varResult += pow(meanResult - targetResponse[repMembrIndx[genIndx[i]]], 2.0);
      }
    }
    varResult = varResult / (double) denom;
    result = ((varResult <= EPSILON) ? FALSE : TRUE);
  }
  else {
    varResult = NA_REAL;
    result = FALSE;
  }
  if (variance != NULL) {
    *variance = varResult;
  }
  return(result);
}
void restoreMeanResponse(uint treeID) {
  Terminal *parent;
  uint leaf;
  uint j;
  for (leaf = 1; leaf <= RF_tLeafCount[treeID]; leaf++) {
    parent = RF_tTermList[treeID][leaf];
    if ((parent -> membrCount) > 0) {
      for (j = 1; j <= RF_rNonFactorCount; j++) {
        (parent -> meanResponse)[j] = RF_TN_REGR_ptr[treeID][leaf][j];
      }
    }
    else {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Zero node count encountered in restoreMeanResponse() in (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
}
void getMultiClassProbNew (uint       treeID,
                           Terminal  *parent,
                           uint      *repMembrIndx,
                           uint       repMembrSize,
                           uint      *allMembrIndx,
                           uint       allMembrSize,
                           uint      *rmbrIterator) {
  uint *membershipIndex;
  uint  membershipSize;
  double maxValue, maxClass;
  uint i, j, k;
  if ( (!(RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ||
       ( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2))) { 
    membershipIndex = allMembrIndx;
    membershipSize = parent -> membrCount = allMembrSize;
  }
  else {
    membershipIndex = repMembrIndx;
    membershipSize = parent -> membrCount = repMembrSize;
    if (RF_optHigh & OPT_MEMB_OUTG) {
      RF_TN_RCNT_ptr[treeID][parent -> nodeID] = RF_tTermList[treeID][parent -> nodeID] -> membrCount;
    }
    if (RF_optHigh & OPT_MEMB_INCG) {
      membershipIndex = RF_RMBR_ID_ptr[treeID];
      membershipSize = parent -> membrCount = RF_TN_RCNT_ptr[treeID][parent -> nodeID];
    }
  }
  if ((parent -> membrCount) == 0) {
    if (!(RF_opt & OPT_OUTC_TYPE)) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Zero node count encountered in (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  if (!(RF_optHigh & OPT_TERM_INCG)) {
    stackMultiClassProb(parent, RF_rFactorCount, RF_rFactorSize);
    for (j=1; j <= RF_rFactorCount; j++) {
      for (k=1; k <= RF_rFactorSize[j]; k++) {
        (parent -> multiClassProb)[j][k] = 0;
      }
    }
    if (RF_optHigh & OPT_MEMB_OUTG) {
      for (i = 1; i <= membershipSize; i++) {
        RF_RMBR_ID_ptr[treeID][++(*rmbrIterator)] = membershipIndex[i];
        for (j=1; j <= RF_rFactorCount; j++) {
          (parent -> multiClassProb)[j][(uint) RF_response[treeID][RF_rFactorIndex[j]][membershipIndex[i]]] ++;
        }
      }
    }
    else if (RF_optHigh & OPT_MEMB_INCG) {
      for (i = 1; i <= membershipSize; i++) {
        ++(*rmbrIterator);
        for (j=1; j <= RF_rFactorCount; j++) {
          (parent -> multiClassProb)[j][(uint) RF_response[treeID][RF_rFactorIndex[j]][ membershipIndex[*rmbrIterator] ]] ++;
        }
      }
    }
    else {
      for (i = 1; i <= membershipSize; i++) {
        for (j=1; j <= RF_rFactorCount; j++) {
          (parent -> multiClassProb)[j][(uint) RF_response[treeID][RF_rFactorIndex[j]][membershipIndex[i]]] ++;
        }
      }
    }
    for (j = 1; j <= RF_rFactorCount; j++) {
      maxValue = 0;
      maxClass = 0;
      for (k=1; k <= RF_rFactorSize[j]; k++) {
        if (maxValue < (double) (parent -> multiClassProb[j][k])) {
          maxValue = (double) parent -> multiClassProb[j][k];
          maxClass = (double) k;
        }
      }
      (parent -> maxClass)[j] = maxClass;
    }
    if (RF_optHigh & OPT_TERM_OUTG) {
      for (j = 1; j <= RF_rFactorCount; j++) {
        for (k = 1; k <= RF_rFactorSize[j]; k++) {
          RF_TN_CLAS_ptr[treeID][parent -> nodeID][j][k] = (parent -> multiClassProb)[j][k];
        }
      }
    }
  }
  else {
    stackMultiClassProbPartial(parent, RF_rFactorCount);
    (parent -> multiClassProb) = RF_TN_CLAS_ptr[treeID][parent -> nodeID];      
    for (j = 1; j <= RF_rFactorCount; j++) {
      maxValue = 0;
      maxClass = 0;
      for (k=1; k <= RF_rFactorSize[j]; k++) {
        if (maxValue < (double) (parent -> multiClassProb[j][k])) {
          maxValue = (double) parent -> multiClassProb[j][k];
          maxClass = (double) k;
        }
      }
      (parent -> maxClass)[j] = maxClass;
    }
  }
}
void updateEnsembleMultiClass(char      mode,
                              uint      treeID,
                              uint      serialTreeID,
                              char      omitDenominator) {
  char oobFlag, fullFlag, selectionFlag, outcomeFlag;
  Terminal ***termMembershipPtr;
  uint    *membershipIndex;
  uint     membershipSize;
  double   ***ensembleCLSptr;
  double   ***ensembleCLSnum;
  uint       *ensembleDen;
  Terminal *parent;
  double maxValue;
  double maxClass;
  uint i, j, k;
  uint ii;
  ensembleCLSptr = NULL;  
  ensembleCLSnum = NULL;  
  ensembleDen    = NULL;  
  oobFlag = fullFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    termMembershipPtr = RF_ftTermMembership;
    break;
  default:
    if (RF_opt & OPT_OENS) {
      if (RF_oobSize[treeID] > 0) {
        oobFlag = TRUE;
      }
    }
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    termMembershipPtr = RF_tTermMembership;
    break;
  }
  outcomeFlag = TRUE;
  while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
    if (oobFlag == TRUE) {
      ensembleCLSptr = RF_oobEnsembleCLSptr;
      ensembleCLSnum = RF_oobEnsembleCLSnum;
      ensembleDen    = RF_oobEnsembleDen;
      membershipSize  = RF_oobSize[treeID];
      membershipIndex = RF_oobMembershipIndex[treeID];
    }
    else {
      ensembleCLSptr = RF_fullEnsembleCLSptr;
      ensembleCLSnum = RF_fullEnsembleCLSnum;
      ensembleDen    = RF_fullEnsembleDen;
      switch (mode) {
      case RF_PRED:
        membershipSize = RF_fobservationSize;
        membershipIndex = RF_fidentityMembershipIndex;
        break;
      default:
        membershipSize  = RF_observationSize;
        membershipIndex = RF_identityMembershipIndex;
        break;
      }
    }
    for (i = 1; i <= membershipSize; i++) {
      ii = membershipIndex[i];
      parent = termMembershipPtr[treeID][ii];
      selectionFlag = TRUE;
      if (RF_opt & OPT_OUTC_TYPE) {
        if ((parent -> maxClass) != NULL) {
        }
        else {
          selectionFlag = FALSE;
        }
      }
      if (selectionFlag) {
        for (j=1; j <= RF_rTargetFactorCount; j++) {
          for (k=1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
            ensembleCLSnum[j][k][ii] += (double) (parent -> multiClassProb)[RF_rFactorMap[RF_rTargetFactor[j]]][k] / (double) (parent -> membrCount);
          }
        }
        if(!omitDenominator) {
          ensembleDen[ii] ++;
        }
      }
      if (outcomeFlag) {
        if (ensembleDen[ii] != 0) {
          for (j=1; j <= RF_rTargetFactorCount; j++) {
            maxValue = 0.0;
            maxClass = 0.0;
            for (k=1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              if (maxValue < ensembleCLSnum[j][k][ii]) {
                maxValue = ensembleCLSnum[j][k][ii];
                maxClass = (double) k;
              }
            }
            ensembleCLSptr[j][1][ii] = maxClass;
          }
        }
      }
    }  
    if (outcomeFlag == TRUE) {
      outcomeFlag = FALSE;
    }
    if (oobFlag == TRUE) {
      oobFlag = FALSE;
    }
    else {
      fullFlag = FALSE;
    }
  }  
}
double getBrierScore(uint     obsSize,
                     uint     rTarget,
                     double  *responsePtr,
                     double **outcomeCLS,
                     uint    *denomCount,
                     double  *cpv) {
  uint k;
  uint against;
  uint *oaaResponse;
  uint cumDenomCount;
  double result;
  oaaResponse = uivector(1, obsSize);
  result = 0.0;
  cumDenomCount = 0;
  for (k = 1; k <= obsSize; k ++) {
    if (denomCount[k] != 0) {
      cumDenomCount += 1;
    }
  }
  for (against = 1; against < RF_rFactorSize[RF_rFactorMap[rTarget]]; against++) {
    for (k = 1; k <= obsSize; k ++) {
      if ((uint) responsePtr[k] == against) {
        oaaResponse[k] = 1;
      }
      else {
        oaaResponse[k] = 0;
      }
    }
    cpv[against] = 0.0;
    for (k = 1; k <= obsSize; k ++) {
      if (denomCount[k] != 0) {
        cpv[against] += pow(((double) oaaResponse[k] - (outcomeCLS[against][k] / (double) denomCount[k])), 2.0);
      }
    }
    if (cumDenomCount == 0) {
      cpv[against] = NA_REAL;
    }
    else {
      cpv[against] = cpv[against] / (double) cumDenomCount;
      result += cpv[against];
    }
  }
  if (cumDenomCount == 0) {
    result = NA_REAL;
  }
  else {
    result = result / RF_rFactorSize[RF_rFactorMap[rTarget]];
  }
  free_uivector(oaaResponse, 1, obsSize);
  return result;
}
void getConditionalClassificationIndex(uint     size,
                                       uint     rTarget,
                                       double  *responsePtr,
                                       double  *outcomeCLS,
                                       uint    *denomCount,
                                       double  *cpv) {
  uint i, k;
  uint cumDenomCount;
  uint *condClassificationCount;
  cumDenomCount = 0;
  condClassificationCount = uivector(1, RF_rFactorSize[RF_rFactorMap[rTarget]]);
  for (k=1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
    cpv[k] = condClassificationCount[k] = 0;
  }
  for (i = 1; i <= size; i++) {
    condClassificationCount[(uint) responsePtr[i]] ++;
    if (denomCount[i] != 0) {
      cumDenomCount += 1;
      if (responsePtr[i] == outcomeCLS[i]) {
        cpv[(uint) responsePtr[i]] += 1.0;
      }
    }
  }  
  if (cumDenomCount == 0) {
    for (k=1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
      cpv[k] = NA_REAL;
    }
  }
  else {
    for (k=1; k <= RF_rFactorSize[RF_rFactorMap[rTarget]]; k++) {
      if (condClassificationCount[k] != 0) {
        cpv[k] = 1.0 - cpv[k] / (double) condClassificationCount[k];
      }
      else {
        cpv[k] = NA_REAL;
      }
    }
  }
  free_uivector(condClassificationCount, 1, RF_rFactorSize[RF_rFactorMap[rTarget]]);
  return;
}
double getClassificationIndex(uint    size,
                              double *responsePtr,
                              double *predictedOutcome,
                              uint   *denomCount) {
  uint i;
  uint cumDenomCount;
  double result;
  cumDenomCount = 0;
  result = 0.0;
  for (i=1; i <= size; i++) {
    if (denomCount[i] > 0) {
      cumDenomCount += 1;
      if (responsePtr[i] == predictedOutcome[i]) {
        result += 1.0;
      }
    }
  }  
  if (cumDenomCount == 0) {
    result = NA_REAL;
  }
  else {
    result = 1.0 - result / (double) cumDenomCount;
  }
  return result;
}
void restoreMultiClassProb(uint treeID) {
  Terminal *parent;
  double maxValue, maxClass;
  uint leaf;
  uint j, k;
  for (leaf = 1; leaf <= RF_tLeafCount[treeID]; leaf++) {
    parent = RF_tTermList[treeID][leaf];
    if ((parent -> membrCount) > 0) {
      for (j = 1; j <= RF_rFactorCount; j++) {
        for (k = 1; k <= RF_rFactorSize[j]; k++) {
          (parent -> multiClassProb)[j][k] = RF_TN_CLAS_ptr[treeID][leaf][j][k];
        }
        maxValue = 0.0;
        maxClass = 0.0;
        for (k=1; k <= RF_rFactorSize[j]; k++) {
          if (maxValue < (double) (parent -> multiClassProb[j][k])) {
            maxValue = (double) parent -> multiClassProb[j][k];
            maxClass = (double) k;
          }
        }
        (parent -> maxClass)[j] = maxClass;
      }
    }
    else {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Zero node count encountered in restoreMultiClassProb() in (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
}
void getAtRiskAndEventCounts(uint       treeID,
                             Terminal  *parent,
                             uint      *repMembrIndx,
                             uint       repMembrSize,
                             uint      *allMembrIndx,
                             uint       allMembrSize,
                             uint      *rmbrIterator) {
  uint *membershipIndex;
  uint  membershipSize;
  uint i, j, k;
  uint ii;
  char eventFlag;
  if ( (!(RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ||
       ( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2))) { 
    membershipIndex = allMembrIndx;
    membershipSize = parent -> membrCount = allMembrSize;
  }
  else {
    membershipIndex = repMembrIndx;
    membershipSize = parent -> membrCount = repMembrSize;
    if (RF_optHigh & OPT_MEMB_OUTG) {
      RF_TN_RCNT_ptr[treeID][parent -> nodeID] = RF_tTermList[treeID][parent -> nodeID] -> membrCount;
    }
    if (RF_optHigh & OPT_MEMB_INCG) {
      membershipIndex = RF_RMBR_ID_ptr[treeID];
      membershipSize = parent -> membrCount = RF_TN_RCNT_ptr[treeID][parent -> nodeID];
    }
  }
  if ((parent -> membrCount) == 0) {
    if (!(RF_opt & OPT_OUTC_TYPE)) {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Zero node count encountered in (tree, leaf) = (%10d, %10d)  \n", treeID, parent -> nodeID);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  if (!(RF_optHigh & OPT_TERM_INCG)) {
    stackAtRiskAndEventCounts(parent, RF_eventTypeSize, RF_masterTimeSize);
    for (j = 1; j <= RF_masterTimeSize; j++) {
      (parent -> atRiskCount)[j]    = 0;
      for (k = 1; k <= RF_eventTypeSize; k++) {
        (parent -> eventCount)[k][j] = 0;
      }
    }
    if (RF_optHigh & OPT_MEMB_OUTG) {
      for (i = 1; i <= membershipSize; i++) {
        ii = membershipIndex[i];
        RF_RMBR_ID_ptr[treeID][++(*rmbrIterator)] = ii;
        for (j = 1; j <= RF_masterTimeIndex[treeID][ii]; j++) {
          (parent -> atRiskCount)[j] ++;
        }
        if (RF_status[treeID][ii] > 0) {
          if (RF_eventTypeSize > 1) {
            k = RF_eventTypeIndex[(uint) RF_status[treeID][ii]];
          }
          else {
            k = 1;
          }
          (parent -> eventCount)[k][RF_masterTimeIndex[treeID][ii]] ++;
        }
      }
    }
    else if (RF_optHigh & OPT_MEMB_INCG) {
      for (i = 1; i <= membershipSize; i++) {
        ii = membershipIndex[++(*rmbrIterator)];
        for (j = 1; j <= RF_masterTimeIndex[treeID][ii]; j++) {
          (parent -> atRiskCount)[j] ++;
        }
        if (RF_status[treeID][ii] > 0) {
          if (RF_eventTypeSize > 1) {
            k = RF_eventTypeIndex[(uint) RF_status[treeID][ii]];
          }
          else {
            k = 1;
          }
          (parent -> eventCount)[k][RF_masterTimeIndex[treeID][ii]] ++;
        }
      }
    }
    else {
      for (i = 1; i <= membershipSize; i++) {
        ii = membershipIndex[i];
        for (j = 1; j <= RF_masterTimeIndex[treeID][ii]; j++) {
          (parent -> atRiskCount)[j] ++;
        }
        if (RF_status[treeID][ii] > 0) {
          if (RF_eventTypeSize > 1) {
            k = RF_eventTypeIndex[(uint) RF_status[treeID][ii]];
          }
          else {
            k = 1;
          }
          (parent -> eventCount)[k][RF_masterTimeIndex[treeID][ii]] ++;
        }
      }
    }
    parent -> eTimeSize = 0;
    for (j = 1; j <= RF_masterTimeSize; j++) {
      eventFlag = FALSE;
      for (k = 1; k <= RF_eventTypeSize; k++) {
        if ((parent -> eventCount)[k][j] > 0) {
          eventFlag = TRUE;
          k = RF_eventTypeSize;
        }
      }
      if (eventFlag == TRUE) {
        (parent -> eTimeSize)++;
      }
    }
    stackEventTimeIndex(parent, parent -> eTimeSize);
    for (j = 1; j <= parent -> eTimeSize; j++) {
      (parent -> eventTimeIndex)[j] = 0;
    }
    i = 0;
    for (j = 1; j <= RF_masterTimeSize; j++) {
      eventFlag = FALSE;
      for (k = 1; k <= RF_eventTypeSize; k++) {
        if ((parent -> eventCount)[k][j] > 0) {
          eventFlag = TRUE;
          k = RF_eventTypeSize;
        }
      }
      if (eventFlag == TRUE) {
        (parent -> eventTimeIndex)[++i] = j;
      }
    }
  }
  else {
  }
}
void getLocalRatios(uint treeID, Terminal *parent) {
  uint j, q;
  if (parent -> membrCount > 0) {
    if(parent -> eTimeSize > 0) {
      stackLocalRatio(parent, RF_eventTypeSize, parent -> eTimeSize);
      for (j = 1; j <= RF_eventTypeSize; j++) {
        for (q = 1; q <= parent -> eTimeSize; q++) {
          if ((parent -> eventCount)[j][(parent -> eventTimeIndex)[q]] > 0) {
            if ((parent -> atRiskCount)[(parent -> eventTimeIndex)[q]] >= 1) {
              (parent -> localRatio)[j][q] = ((double) (parent -> eventCount)[j][(parent -> eventTimeIndex)[q]] / (parent -> atRiskCount)[(parent -> eventTimeIndex)[q]]);
            }
            else {
              RFprintf("\nRF-SRC:  *** ERROR *** ");
              RFprintf("\nRF-SRC:  Zero At Risk Count encountered in local ratio calculation for (tree, leaf) = (%10d, %10d)", treeID, parent -> nodeID);
              RFprintf("\nRF-SRC:  Please Contact Technical Support.");
              error("\nRF-SRC:  The application will now exit.\n");
            }
          }
          else {
            (parent -> localRatio)[j][q] = 0.0;
          }
        }
      }
    }
  }
}
void getLocalSurvival(uint treeID, Terminal *parent) {
  uint j, q;
  if(parent -> eTimeSize > 0) {
    stackLocalSurvival(parent, parent -> eTimeSize);
    for (q = 1; q <= parent -> eTimeSize; q++) {
      (parent -> localSurvival)[q] = 0.0;
      for (j = 1; j <= RF_eventTypeSize; j++) {
        (parent -> localSurvival)[q] += (parent -> localRatio)[j][q];
      }
      (parent -> localSurvival)[q] = 1.0 - (parent -> localSurvival)[q];
    }  
    for (q = 2; q <= parent -> eTimeSize; q++) {
      (parent -> localSurvival)[q] *= (parent -> localSurvival)[q-1];
    }
  }
}
void getLocalNelsonAalen(uint treeID, Terminal *parent) {
  uint q;
  if (parent -> eTimeSize > 0) {
    stackLocalNelsonAalen(parent, parent -> eTimeSize);
    for (q = 1; q <= parent -> eTimeSize; q++) {
      (parent -> localNelsonAalen)[q] = (parent -> localRatio)[1][q];
    }
    for (q = 2; q <= parent -> eTimeSize; q++) {
      (parent -> localNelsonAalen)[q] += (parent -> localNelsonAalen)[q-1];
    }
  }
}
void getLocalCSH(uint treeID, Terminal *parent) {
  uint j, q;
    if (parent -> eTimeSize > 0) {
      stackLocalCSH(parent, RF_eventTypeSize, parent -> eTimeSize);
      for (j = 1; j <= RF_eventTypeSize; j++) {
        for (q = 1; q <= parent -> eTimeSize; q++) {
          (parent -> localCSH)[j][q] = (parent -> localRatio)[j][q];
        }
        for (q = 2; q <= parent -> eTimeSize; q++) {
          (parent -> localCSH)[j][q] += (parent -> localCSH)[j][q-1];
        }
      }
    }
}
void getLocalCIF(uint treeID, Terminal *parent) {
  uint j, q;
  if(parent -> eTimeSize > 0) {
    stackLocalCIF(parent, RF_eventTypeSize, parent -> eTimeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      (parent -> localCIF)[j][1] = (parent -> localRatio)[j][1];
      for (q = 2; q <= parent -> eTimeSize; q++) {
        (parent -> localCIF)[j][q] = (parent -> localSurvival)[q-1] * (parent -> localRatio)[j][q];
      }
      for (q = 2; q <= parent -> eTimeSize; q++) {
        (parent -> localCIF)[j][q] += (parent -> localCIF)[j][q-1];
      }
    }
  }
}
void mapLocalToTimeInterest(uint      treeID,
                            Terminal *parent,
                            void     *genericLocal,
                            void     *genericGlobal) {
  uint itIndex, etIndex, lookAheadIndex;
  char mapFlag, transitFlag;
  uint j;
  if (!(RF_opt & OPT_COMP_RISK)) {
    if ((parent -> eTimeSize) > 0) {
      itIndex = 1;
      etIndex = 1;
      mapFlag = TRUE;
      while(mapFlag) {
        if (RF_timeInterest[itIndex] < RF_masterTime[(parent -> eventTimeIndex)[etIndex]] ) {
          if (itIndex > 1) {
            ((double *) genericGlobal)[itIndex] = ((double *) genericGlobal)[itIndex-1];
          }
          itIndex++;
        }
        else {
          lookAheadIndex = etIndex;
          transitFlag = TRUE;
          while (transitFlag) {
            if (RF_timeInterest[itIndex] >= RF_masterTime[(parent -> eventTimeIndex)[lookAheadIndex]] ) {
              ((double *) genericGlobal)[itIndex] = ((double *) genericLocal)[lookAheadIndex];
              lookAheadIndex++;
              if (lookAheadIndex > (parent -> eTimeSize)) {
                transitFlag = FALSE;
              }
            }
            else {
              transitFlag = FALSE;
            }
          }
          itIndex++;
          etIndex = lookAheadIndex;
        }
        if(etIndex > (parent -> eTimeSize)) {
          while(itIndex <= RF_sortedTimeInterestSize) {
            ((double *) genericGlobal)[itIndex] = ((double *) genericGlobal)[itIndex-1];
            itIndex++;
          }
        }
        if(itIndex > RF_sortedTimeInterestSize) {
          mapFlag = FALSE;
        }
      }
    }
  }  
  else {
    if ((parent -> eTimeSize) > 0) {
      itIndex = 1;
      etIndex = 1;
      mapFlag = TRUE;
      while(mapFlag) {
        if (RF_timeInterest[itIndex] < RF_masterTime[(parent -> eventTimeIndex)[etIndex]] ) {
          if (itIndex > 1) {
            for (j = 1; j <= RF_eventTypeSize; j++) {
              ((double **) genericGlobal)[j][itIndex] = ((double **) genericGlobal)[j][itIndex-1];
            }
          }
          itIndex++;
        }
        else {
          lookAheadIndex = etIndex;
          transitFlag = TRUE;
          while (transitFlag) {
            if (RF_timeInterest[itIndex] >= RF_masterTime[(parent -> eventTimeIndex)[lookAheadIndex]] ) {
              for (j = 1; j <= RF_eventTypeSize; j++) {
                ((double **) genericGlobal)[j][itIndex] = ((double **) genericLocal)[j][lookAheadIndex];
              }
              lookAheadIndex++;
              if (lookAheadIndex > (parent -> eTimeSize)) {
                transitFlag = FALSE;
              }
            }
            else {
              transitFlag = FALSE;
            }
          }
          itIndex++;
          etIndex = lookAheadIndex;
        }
        if(etIndex > (parent -> eTimeSize)) {
          while(itIndex <= RF_sortedTimeInterestSize) {
              for (j = 1; j <= RF_eventTypeSize; j++) {
                ((double **) genericGlobal)[j][itIndex] = ((double **) genericGlobal)[j][itIndex-1];
              }
              itIndex++;
          }
        }
        if(itIndex > RF_sortedTimeInterestSize) {
          mapFlag = FALSE;
        }
      }
    }    
  }  
}
void getSurvival(uint treeID, Terminal *parent) {
  uint k;
  if (!(RF_optHigh & OPT_TERM_INCG)) {
    stackSurvival(parent, RF_sortedTimeInterestSize);
    for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
      (parent -> survival)[k] = 1.0;
    }
    mapLocalToTimeInterest(treeID,
                           parent,
                           parent -> localSurvival,
                           parent -> survival);
    if (RF_optHigh & OPT_TERM_OUTG) {
      for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
        RF_TN_SURV_ptr[treeID][parent -> nodeID][k] = parent -> survival[k];
      }
    }
  }
  else {
    (parent -> survival) = RF_TN_SURV_ptr[treeID][parent -> nodeID];
  }
}
void getNelsonAalen(uint treeID, Terminal *parent) {
  uint k;
  if (!(RF_optHigh & OPT_TERM_INCG)) {
    stackNelsonAalen(parent, RF_sortedTimeInterestSize);
    for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
      (parent -> nelsonAalen)[k] = 0.0;
    }
    mapLocalToTimeInterest(treeID,
                           parent,
                           parent -> localNelsonAalen,
                           parent -> nelsonAalen);
    if (RF_optHigh & OPT_TERM_OUTG) {
      for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
        RF_TN_NLSN_ptr[treeID][parent -> nodeID][k] = parent -> nelsonAalen[k];
      }
    }
  }
  else {
    (parent -> nelsonAalen) = RF_TN_NLSN_ptr[treeID][parent -> nodeID];
  }
}
void getCSH(uint treeID, Terminal *parent) {
  uint j, k;
  if (!(RF_optHigh & OPT_TERM_INCG)) {
    stackCSH(parent, RF_eventTypeSize, RF_sortedTimeInterestSize);
    for (j=1; j <= RF_eventTypeSize; j++) {
      for (k=1; k <= RF_sortedTimeInterestSize; k++) {
        (parent -> CSH)[j][k] = 0.0;
      }
    }
    mapLocalToTimeInterest(treeID,
                           parent,
                           parent -> localCSH,
                           parent -> CSH);
    if (RF_optHigh & OPT_TERM_OUTG) {
      for (j = 1; j <= RF_eventTypeSize; j++) {
        for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
          RF_TN_CSHZ_ptr[treeID][parent -> nodeID][j][k] = RF_tTermList[treeID][parent -> nodeID] -> CSH[j][k];
        }
      }
    }      
  }
  else {
    (parent -> CSH) = RF_TN_CSHZ_ptr[treeID][parent -> nodeID];
  }
}
void getCIF(uint treeID, Terminal *parent) {
  uint j, k;
  if (!(RF_optHigh & OPT_TERM_INCG)) {
    stackCIF(parent, RF_eventTypeSize, RF_sortedTimeInterestSize);
    for (j=1; j <= RF_eventTypeSize; j++) {
      for (k=1; k <= RF_sortedTimeInterestSize; k++) {
        (parent -> CIF)[j][k] = 0.0;
      }
    }
    mapLocalToTimeInterest(treeID,
                           parent,
                           parent -> localCIF,
                           parent -> CIF);
    if (RF_optHigh & OPT_TERM_OUTG) {
      for (j = 1; j <= RF_eventTypeSize; j++) {
        for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
          RF_TN_CSHZ_ptr[treeID][parent -> nodeID][j][k] = RF_tTermList[treeID][parent -> nodeID] -> CSH[j][k];
        }
      }
    }
  }
  else {
    (parent -> CIF) = RF_TN_CIFN_ptr[treeID][parent -> nodeID];
  }
}
void getMortality(uint treeID, Terminal *parent) {
  uint j, q;
  if (!(RF_optHigh & OPT_TERM_INCG)) {
    stackMortality(parent, RF_eventTypeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      (parent -> mortality)[j] = 0.0;
    }
    if (!(RF_opt & OPT_COMP_RISK)) {
      for (q = 1; q <= RF_sortedTimeInterestSize; q++) {
        (parent -> mortality)[1] += (parent -> nelsonAalen)[q];
      }
    }
    else {
      for (j = 1; j <= RF_eventTypeSize; j ++) {
        for (q = 1; q <= RF_sortedTimeInterestSize - 1; q++) {
          (parent -> mortality)[j] += (parent -> CIF)[j][q] * (RF_timeInterest[q+1] - RF_timeInterest[q]);
        }
      }
    }
    if (RF_optHigh & OPT_TERM_OUTG) {
      for (j = 1; j <= RF_eventTypeSize; j++) {
        RF_TN_MORT_ptr[treeID][parent -> nodeID][j] = parent -> mortality[j];
      }
    }
  }
  else {
    (parent -> mortality) = RF_TN_MORT_ptr[treeID][parent -> nodeID];
  }
}
void updateEnsembleSurvival(char mode,
                            uint treeID,
                            uint serialTreeID) {
  char oobFlag, fullFlag, selectionFlag, outcomeFlag;
  Terminal ***termMembershipPtr;
  uint    *membershipIndex;
  uint     membershipSize;
  double ***ensembleSRGnum;
  double ***ensembleCIFnum;
  double  **ensembleSRVnum;
  double  **ensembleMRTptr;
  double  **ensembleMRTnum;
  uint     *ensembleDen;
  Terminal *parent;
  uint i, j, k;
  uint ii;
  ensembleSRGnum = NULL;  
  ensembleCIFnum = NULL;  
  ensembleSRVnum = NULL;  
  ensembleMRTptr = NULL;  
  ensembleMRTnum = NULL;  
  ensembleDen    = NULL;  
  oobFlag = fullFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    termMembershipPtr = RF_ftTermMembership;
    break;
  default:
    if (RF_opt & OPT_OENS) {
      if (RF_oobSize[treeID] > 0) {
        oobFlag = TRUE;
      }
    }
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    termMembershipPtr = RF_tTermMembership;
    break;
  }
  outcomeFlag = TRUE;
  while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
    if (oobFlag == TRUE) {
      ensembleSRGnum = RF_oobEnsembleSRGnum;
      ensembleMRTptr = RF_oobEnsembleMRTptr;        
      ensembleMRTnum = RF_oobEnsembleMRTnum;
      ensembleSRVnum = RF_oobEnsembleSRVnum;
      ensembleCIFnum = RF_oobEnsembleCIFnum;
      ensembleDen    = RF_oobEnsembleDen;
      membershipSize  = RF_oobSize[treeID];
      membershipIndex = RF_oobMembershipIndex[treeID];
    }
    else {
      ensembleSRGnum = RF_fullEnsembleSRGnum;
      ensembleMRTptr = RF_fullEnsembleMRTptr;
      ensembleMRTnum = RF_fullEnsembleMRTnum;        
      ensembleSRVnum = RF_fullEnsembleSRVnum;
      ensembleCIFnum = RF_fullEnsembleCIFnum;
      ensembleDen    = RF_fullEnsembleDen;
      switch (mode) {
      case RF_PRED:
        membershipSize = RF_fobservationSize;
        membershipIndex = RF_fidentityMembershipIndex;
        break;
      default:
        membershipSize  = RF_observationSize;
        membershipIndex = RF_identityMembershipIndex;
        break;
      }
    }
    for (i = 1; i <= membershipSize; i++) {
      ii = membershipIndex[i];
      parent = termMembershipPtr[treeID][ii];
      selectionFlag = TRUE;
      if (RF_opt & OPT_OUTC_TYPE) {
        if ((parent -> membrCount) > 0) {
        }
        else {
          selectionFlag = FALSE;
        }
      }
      if (selectionFlag) {
        ensembleDen[ii] ++;
        if (!(RF_opt & OPT_COMP_RISK)) {
          ensembleMRTnum[1][ii] += parent -> mortality[1];
          for (k=1; k <= RF_sortedTimeInterestSize; k++) {
            ensembleSRGnum[1][k][ii] += parent -> nelsonAalen[k];
            ensembleSRVnum[k][ii] += parent -> survival[k];
          }
          if (outcomeFlag) {
            if (ensembleDen[ii] != 0) {
              ensembleMRTptr[1][ii] = ensembleMRTnum[1][ii] / ensembleDen[ii];
            }
          }
        }
        else {
          for (j = 1; j <= RF_eventTypeSize; j++) {
            ensembleMRTnum[j][ii] += parent -> mortality[j];
            for (k=1; k <= RF_sortedTimeInterestSize; k++) {
              ensembleSRGnum[j][k][ii] += parent -> CSH[j][k];
              ensembleCIFnum[j][k][ii] += parent -> CIF[j][k];
            }
          }
          if (outcomeFlag) {
            if (ensembleDen[ii] != 0) {
              for (j = 1; j <= RF_eventTypeSize; j ++) {
                ensembleMRTptr[j][ii] = ensembleMRTnum[j][ii] / ensembleDen[ii];
              }
            }
          }
        }
      }
    }  
    if (outcomeFlag == TRUE) {
      outcomeFlag = FALSE;
    }
    if (oobFlag == TRUE) {
      oobFlag = FALSE;
    }
    else {
      fullFlag = FALSE;
    }
  }  
}
void getEnsembleMortalityCR(char      mode,
                            uint      treeID,
                            uint      obsSize,
                            double  **ensembleMRTptr,
                            uint     *ensembleDen,
                            double  **cMortality) {
  uint i, j;
  for (i = 1; i <= obsSize; i++) {
    if (ensembleDen[i] != 0) {
      for (j = 1; j <= RF_eventTypeSize; j ++) {
        cMortality[j][i] = ensembleMRTptr[j][i] / ensembleDen[i];
      }
    }
    else {
      for (j = 1; j <= RF_eventTypeSize; j ++) {
        cMortality[j][i] = NA_REAL;
      }
    }
  }
}
void getEnsembleMortality(char      mode,
                          uint      treeID,
                          uint      obsSize,
                          double  **ensembleMRTptr,
                          uint     *ensembleDen,
                          double   *mortality) {
  uint i;
  for (i = 1; i <= obsSize; i++) {
    if (ensembleDen[i] != 0) {
      mortality[i] = ensembleMRTptr[1][i] / ensembleDen[i];
    }
    else {
      mortality[i] = NA_REAL;
    }
  }
}
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
                                     uint    *subsettedEnsembleDen) {
  uint i;
  if (!(RF_opt & OPT_COMP_RISK)) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Attempt to update event type subsets in a non-CR analysis.");
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  for (i = 1; i <= meIndividualSize[j]; i++) {
    subsettedTime[i]        = timePtr[eIndividual[j][i]];
    subsettedStatus[i]      = statusPtr[eIndividual[j][i]];
    subsettedMortality[i]   = mortalityPtr[eIndividual[j][i]];
    subsettedEnsembleDen[i] = genericEnsembleDenPtr[eIndividual[j][i]];
  }
}
double getConcordanceIndex(int     polarity,
                           uint    size,
                           double *timePtr,
                           double *statusPtr,
                           double *predictedOutcome,
                           uint   *denCount) {
  uint i,j;
  uint concordancePairSize;
  uint concordanceWorseCount;
  double result;
  concordancePairSize = concordanceWorseCount = 0;
  for (i=1; i < size; i++) {
    for (j=i+1; j <= size; j++) {
      if (denCount[i] != 0  && denCount[j] != 0) {
        if ( ((timePtr[i] - timePtr[j] > EPSILON) && (statusPtr[j] > 0)) ||
             ((fabs(timePtr[i] - timePtr[j]) <= EPSILON) && (statusPtr[j] > 0) && (statusPtr[i] == 0)) ) {
          concordancePairSize += 2;
          if (predictedOutcome[j] - predictedOutcome[i] > EPSILON) {
            concordanceWorseCount += 2;
          }
          else if (fabs(predictedOutcome[j] - predictedOutcome[i]) < EPSILON) {
            concordanceWorseCount += 1;
          }
        }
        else if ( ((timePtr[j] - timePtr[i]) > EPSILON  && (statusPtr[i] > 0)) ||
                  ((fabs(timePtr[j] - timePtr[i]) <= EPSILON)  && (statusPtr[i] > 0) && (statusPtr[j] == 0)) ) {
          concordancePairSize += 2;
          if ( predictedOutcome[i] - predictedOutcome[j] > EPSILON ) {
            concordanceWorseCount += 2;
          }
          else if (fabs(predictedOutcome[i] - predictedOutcome[j]) < EPSILON) {
            concordanceWorseCount += 1;
          }
        }
        else if ( (fabs(timePtr[i]- timePtr[j]) <= EPSILON) && (statusPtr[i] > 0) && (statusPtr[j] > 0) ) {
          concordancePairSize += 2;
          if (fabs(predictedOutcome[i] - predictedOutcome[j]) < EPSILON) {
            concordanceWorseCount += 2;
          }
          else {
            concordanceWorseCount += 1;
          }
        }
      }  
    }  
  }  
  if (concordancePairSize == 0) {
    result = NA_REAL;
  }
  else {
    result = 1.0 - ((double) concordanceWorseCount / (double) concordancePairSize);
  }
  return result;
}
void getCRPerformance (char     mode,
                       uint     obsSize,
                       double **responsePtr,
                       double **yearsLost,
                       uint    *denom,
                       double  *performanceVector) {
  uint   mRecordSize;
  int  **mpSign;
  uint  *mRecordIndex;
  uint  *meIndividualSize;
  uint **eIndividual;
  double concordanceIndex;
  uint j;
  if (!(RF_opt & OPT_COMP_RISK)) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Attempt at conditional performance updates in a non-CR analysis.");
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (RF_mStatusSize > 0) {
    switch (mode) {
    case RF_PRED:
      mRecordSize = RF_fmRecordSize;
      mpSign = RF_fmpSign;
      mRecordIndex = RF_fmRecordIndex;
      break;
    default:
      mRecordSize = RF_mRecordSize;
      mpSign = RF_mpSign;
      mRecordIndex = RF_mRecordIndex;
      break;
    }
    meIndividualSize  = uivector(1, RF_eventTypeSize);
    eIndividual = (uint **) new_vvector(1, RF_eventTypeSize, NRUTIL_UPTR);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      eIndividual[j] = uivector(1, RF_eIndividualSize[j] + RF_mStatusSize + 1);
    }
    updateEventTypeSubsets(responsePtr[RF_statusIndex], mRecordSize, mpSign, mRecordIndex, meIndividualSize, eIndividual);
  }
  else {
    meIndividualSize  = RF_eIndividualSize;
    eIndividual = RF_eIndividualIn;
  }
  double *subsettedTime      = dvector(1, obsSize);
  double *subsettedStatus    = dvector(1, obsSize);
  double *subsettedMortality = dvector(1, obsSize);
  uint *subsettedEnsembleDen = uivector(1, obsSize);
  for (j = 1; j <= RF_eventTypeSize; j++) {
    getConditionalConcordanceArrays(j,
                                    responsePtr[RF_timeIndex],
                                    responsePtr[RF_statusIndex],
                                    yearsLost[j],
                                    denom,
                                    meIndividualSize,
                                    eIndividual,
                                    subsettedTime,
                                    subsettedStatus,
                                    subsettedMortality,
                                    subsettedEnsembleDen);
    concordanceIndex = getConcordanceIndex(1,
                                           meIndividualSize[j],
                                           subsettedTime,
                                           subsettedStatus,
                                           subsettedMortality,
                                           subsettedEnsembleDen);
    if (ISNA(concordanceIndex)) {
      performanceVector[j] = NA_REAL;
    }
    else {
      performanceVector[j] = concordanceIndex;
    }
  }
  if (RF_mStatusSize > 0) {
    free_uivector(meIndividualSize, 1, RF_eventTypeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      free_uivector(eIndividual[j], 1, RF_eIndividualSize[j] + RF_mStatusSize + 1);
    }
    free_new_vvector(eIndividual, 1, RF_eventTypeSize, NRUTIL_UPTR);
  }
  free_dvector(subsettedTime, 1, obsSize);
  free_dvector(subsettedStatus, 1, obsSize);
  free_dvector(subsettedMortality, 1, obsSize);
  free_uivector(subsettedEnsembleDen, 1, obsSize);
}
char imputeNode (char     type,
                 char     termFlag,
                 char     chainFlag,
                 uint     treeID,
                 Node    *nodePtr,
                 uint    *repMembrIndx,
                 uint     repMembrSize,
                 uint    *allMembrIndx,
                 uint     allMembrSize) {
  double  **response;
  double  **predictor;
  double    imputedValue;
  uint     *mRecordMap;
  uint      mpIndexSize;
  int     **mpSign;
  int      *mpIndex;
  int      *mvNSptr;
  uint      mRecordSize;
  double *valuePtr, *imputePtr;
  char mPredictorFlag;
  int  signedSignatureIndex;
  uint unsignedIndexSource;
  uint unsignedIndexTarget;
  char result;
  uint   *glmpIndexPtr;
  double *glmpValuePtr;
  uint   *glmpIndexSize;
  uint  *glmpIndexParentPtr;
  uint   glmpIndexParentSize;
  char mvFlag;
  uint i,p;
  uint localDistributionSize;
  mvNSptr = NULL;  
  mpIndex = NULL;  
  mpSign  = NULL;  
  mpIndexSize  = 0;  
  mRecordMap = NULL;  
  mRecordSize = 0;    
  predictor  = NULL;  
  response   = NULL;  
  imputedValue = 0.0;  
  glmpIndexPtr = NULL;
  glmpValuePtr = NULL;
  glmpIndexSize = NULL;
  glmpIndexParentPtr = NULL;
  glmpIndexParentSize = 0;
  result = FALSE;
  switch (type) {
  case RF_PRED:
    mRecordSize = RF_fmRecordSize;
    if (mRecordSize > 0) {
      response = RF_fresponse[treeID];
      predictor = RF_fobservation[treeID];
      mRecordMap = RF_fmRecordMap;
      mpIndexSize = RF_fmpIndexSize;
      mpSign = RF_fmpSign;
      mpIndex = RF_fmpIndex;
      mvNSptr = nodePtr -> fmpSign;
      if (!termFlag) {
        if((nodePtr -> parent) == NULL) {
          glmpIndexParentPtr = uivector(1, mpIndexSize);
          glmpIndexParentSize = mpIndexSize;
          for (p = 1; p <= glmpIndexParentSize; p++) {
            glmpIndexParentPtr[p] = p;
          }
          stackNodeFLMPIndex(nodePtr, glmpIndexParentSize);
          glmpIndexPtr  = nodePtr -> flmpIndex;
          glmpIndexSize = & (nodePtr -> flmpIndexActualSize);
          *glmpIndexSize = 0;
        }
        else {
          if((nodePtr -> parent) -> flmpIndexActualSize > 0) {
            glmpIndexParentPtr = (nodePtr -> parent) -> flmpIndex;
            glmpIndexParentSize = (nodePtr -> parent) -> flmpIndexActualSize;
            stackNodeFLMPIndex(nodePtr, glmpIndexParentSize);
            glmpIndexPtr  = nodePtr -> flmpIndex;
            glmpIndexSize = & (nodePtr -> flmpIndexActualSize);
            *glmpIndexSize = 0;
          }
          else {
            glmpIndexParentPtr  = NULL;
            glmpIndexParentSize = 0;
            glmpIndexPtr = glmpIndexSize = NULL;
          }
        }  
      }  
      else {
        glmpIndexParentPtr = uivector(1, mpIndexSize);
        glmpIndexParentSize = mpIndexSize;
        for (p = 1; p <= glmpIndexParentSize; p++) {
          glmpIndexParentPtr[p] = p;
        }
        stackNodeFLMPIndex(nodePtr, glmpIndexParentSize);
        glmpIndexPtr  = nodePtr -> flmpIndex;
        glmpValuePtr  = nodePtr -> flmpValue;
        glmpIndexSize = & (nodePtr -> flmpIndexActualSize);
        *glmpIndexSize = 0;
      }
      result = TRUE;
    }
    break;
  default:
    mRecordSize = RF_mRecordSize;
    if (mRecordSize > 0) {
      response = RF_response[treeID];
      predictor = RF_observation[treeID];
      mRecordMap = RF_mRecordMap;
      mpIndexSize = RF_mpIndexSize;
      mpSign = RF_mpSign;
      mpIndex = RF_mpIndex;
      mvNSptr = nodePtr -> mpSign;
      if (!termFlag) {
        if((nodePtr -> parent) == NULL) {
          glmpIndexParentPtr = uivector(1, mpIndexSize);
          glmpIndexParentSize = mpIndexSize;
          for (p = 1; p <= glmpIndexParentSize; p++) {
            glmpIndexParentPtr[p] = p;
          }
          stackNodeLMPIndex(nodePtr, glmpIndexParentSize);
          glmpIndexPtr  = nodePtr -> lmpIndex;
          glmpIndexSize = & (nodePtr -> lmpIndexActualSize);
          *glmpIndexSize = 0;
        }
        else {
          if((nodePtr -> parent) -> lmpIndexActualSize > 0) {
            glmpIndexParentPtr = (nodePtr -> parent) -> lmpIndex;
            glmpIndexParentSize = (nodePtr -> parent) -> lmpIndexActualSize;
            stackNodeLMPIndex(nodePtr, glmpIndexParentSize);
            glmpIndexPtr  = nodePtr -> lmpIndex;
            glmpIndexSize = & (nodePtr -> lmpIndexActualSize);
            *glmpIndexSize = 0;
          }
          else {
            glmpIndexParentPtr = NULL;
            glmpIndexParentSize = 0;
            glmpIndexPtr = glmpIndexSize = NULL;
          }
        }  
      }  
      else {
        glmpIndexParentPtr = uivector(1, mpIndexSize);
        glmpIndexParentSize = mpIndexSize;
        for (p = 1; p <= glmpIndexParentSize; p++) {
          glmpIndexParentPtr[p] = p;
        }
        stackNodeLMPIndex(nodePtr, glmpIndexParentSize);
        glmpIndexPtr  = nodePtr -> lmpIndex;
        glmpValuePtr  = nodePtr -> lmpValue;
        glmpIndexSize = & (nodePtr -> lmpIndexActualSize);
        *glmpIndexSize = 0;
      }
      result = TRUE;
    }
    break;
  }
  if (result == FALSE) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Attempt to impute node with no missingness in type:  %10d", type);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  nodePtr -> imputed = TRUE;
  double *localDistribution = dvector(1, repMembrSize + 1);
  for (p = 1; p <= glmpIndexParentSize; p++) {
    if (mvNSptr[glmpIndexParentPtr[p]] != -1) {
      signedSignatureIndex = mpIndex[glmpIndexParentPtr[p]];
      if (signedSignatureIndex < 0) {
        unsignedIndexSource = unsignedIndexTarget = (uint) abs(signedSignatureIndex);
        valuePtr = RF_response[treeID][(uint) abs(signedSignatureIndex)];
        imputePtr = response[(uint) abs(signedSignatureIndex)];
      }
      else {
        unsignedIndexSource = RF_rSize + (uint) signedSignatureIndex;
        switch (type) {
        case RF_PRED:
          if (RF_frSize > 0) {
            unsignedIndexTarget = RF_rSize + (uint) signedSignatureIndex;
          }
          else {
            unsignedIndexTarget = (uint) signedSignatureIndex;
          }
          break;
        default:
          unsignedIndexTarget = RF_rSize + (uint) signedSignatureIndex;
          break;
        }
        valuePtr = RF_observation[treeID][(uint) signedSignatureIndex];
        imputePtr = predictor[(uint) signedSignatureIndex];
      }
      localDistributionSize = 0;
      for (i = 1; i <= repMembrSize; i++) {
        mPredictorFlag = TRUE;
        if (RF_mRecordMap[repMembrIndx[i]] == 0) {
          mPredictorFlag = FALSE;
        }
        else if (RF_mpSign[unsignedIndexSource][RF_mRecordMap[repMembrIndx[i]]] == 0) {
          mPredictorFlag = FALSE;
        }
        if (mPredictorFlag == FALSE) {
          localDistributionSize ++;
          localDistribution[localDistributionSize] = valuePtr[repMembrIndx[i]];
        }
      }  
      if (termFlag) {
        if (localDistributionSize == 0) {
          for (i = 1; i <= repMembrSize; i++) {
            mPredictorFlag = TRUE;
            if (RF_mRecordMap[repMembrIndx[i]] == 0) {
              mPredictorFlag = FALSE;
            }
            else if (RF_mpSign[unsignedIndexSource][RF_mRecordMap[repMembrIndx[i]]] == 0) {
              mPredictorFlag = FALSE;
            }
            if (mPredictorFlag == TRUE) {
              localDistributionSize ++;
              localDistribution[localDistributionSize] = valuePtr[repMembrIndx[i]];
            }
          }  
        }
        if (localDistributionSize > 0) {
          if (signedSignatureIndex < 0) {
            if (strcmp(RF_rType[(uint) abs(signedSignatureIndex)], "T") == 0) {
              imputedValue = getMeanValue(localDistribution, localDistributionSize);
              imputedValue = getNearestMasterTime(imputedValue, chainFlag, treeID);
            }
            else if (strcmp(RF_rType[(uint) abs(signedSignatureIndex)], "S") == 0) {
              imputedValue = getMaximalValue(localDistribution, localDistributionSize, chainFlag, treeID);
            }
            else if (strcmp(RF_rType[(uint) abs(signedSignatureIndex)], "R") == 0) {
              imputedValue = getMeanValue(localDistribution, localDistributionSize);
            }
            else if (strcmp(RF_rType[(uint) abs(signedSignatureIndex)], "I") == 0) {
              imputedValue = getMaximalValue(localDistribution, localDistributionSize, chainFlag, treeID);
            }
            else if (strcmp(RF_rType[(uint) abs(signedSignatureIndex)], "C") == 0) {
              imputedValue = getMaximalValue(localDistribution, localDistributionSize, chainFlag, treeID);
            }
          }
          else {
            if (strcmp(RF_xType[(uint) signedSignatureIndex], "R") == 0) {
              imputedValue = getMeanValue(localDistribution, localDistributionSize);
            }
            else {
              imputedValue = getMaximalValue(localDistribution, localDistributionSize, chainFlag, treeID);
            }
          }
        }  
        else {
          if (!(RF_opt & OPT_OUTC_TYPE)) {
            RFprintf("\nRF-SRC:  *** ERROR *** ");
            RFprintf("\nRF-SRC:  NULL distribution encountered during imputation in type:  %10d", type);
            RFprintf("\nRF-SRC:  Please Contact Technical Support.");
            error("\nRF-SRC:  The application will now exit.\n");
          }
          else {
            imputedValue = NA_REAL;
          }
        }
      }  
      mvFlag = FALSE;
      for (i = 1; i <= allMembrSize; i++) {
        if (mRecordMap[allMembrIndx[i]] > 0) {
          if(mpSign[unsignedIndexTarget][mRecordMap[allMembrIndx[i]]] == 1) {
            mvFlag = TRUE;
            if (localDistributionSize > 0) {
              if (termFlag) {
                imputePtr[allMembrIndx[i]] = imputedValue;
              }
              else {
                imputePtr[allMembrIndx[i]] = getSampleValue(localDistribution, localDistributionSize, chainFlag, treeID);
              }
            }
            else {
            }
          }
        }
      }  
      if (mvFlag) {
        glmpIndexPtr[++(*glmpIndexSize)] = glmpIndexParentPtr[p];
        if (termFlag) {
          glmpValuePtr[(*glmpIndexSize)] = imputedValue;
        }
      }
    }  
  }  
  free_dvector(localDistribution, 1, repMembrSize + 1);
  if (!termFlag) {
    if((nodePtr -> parent) == NULL) {
      free_uivector(glmpIndexParentPtr, 1, mpIndexSize);
    }
  }
  else {
    free_uivector(glmpIndexParentPtr, 1, mpIndexSize);
  }
  if((nodePtr -> parent) != NULL) {
    if( ((((nodePtr -> parent) -> left) -> imputed) == TRUE) && ((((nodePtr -> parent) -> right) -> imputed) == TRUE) ) {
      switch (type) {
      case RF_PRED:
        unstackNodeFLMPIndex(nodePtr -> parent);
        break;
      default:
        unstackNodeLMPIndex(nodePtr -> parent);
        break;
      }
    }
  }
  return TRUE;
}  
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
                           uint *ambrIterator) {
  char  bootResult;
  char leftResult, rghtResult;
  char terminalFlag;
  char bsUpdateFlag;
  uint *bootMembrIndx;
  uint *leftRepMembrIndx;
  uint *rghtRepMembrIndx;
  uint *leftAllMembrIndx;
  uint *rghtAllMembrIndx;
  uint *ngLeftAllMembrIndx;  
  uint *ngRghtAllMembrIndx;  
  uint bootMembrSize;
  uint leftRepMembrSize, rghtRepMembrSize;
  uint leftAllMembrSize, ngLeftAllMembrSize;
  uint rghtAllMembrSize, ngRghtAllMembrSize;
  uint jLeft;
  uint jRght;
  char factorFlag;
  char daughterFlag;
  uint i;
  factorFlag = FALSE; 
  bootResult = TRUE;
  terminalFlag = TRUE;
  bsUpdateFlag = FALSE;
  if (rootFlag || ((RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2))) {
    if ( (!(RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
         ( (RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ) {
      bootMembrIndx  = uivector(1, RF_bootstrapSize);
      bootMembrSize = RF_bootstrapSize;
    }
    else {
      bootMembrIndx  = uivector(1, allMembrSize);
      bootMembrSize = allMembrSize;
    }
    bootResult = bootstrap (mode,
                            treeID,
                            parent,
                            allMembrIndx,
                            allMembrSize,
                            bootMembrIndx,
                            bootMembrSize);
    if (rootFlag & bootResult) {
      if (!( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2) )) {
        bsUpdateFlag = TRUE;
      }
      repMembrIndx = bootMembrIndx;
      repMembrSize = bootMembrSize;
    }
  }
  else {
    bootMembrIndx = repMembrIndx;
    bootMembrSize = repMembrSize;
    parent -> mpSign = (parent -> parent) -> mpSign;
    parent -> fmpSign = (parent -> parent) -> fmpSign;
  }
  if (bootResult) {
    if (RF_mRecordSize > 0) {
      imputeNode(RF_GROW,
                 FALSE,
                 TRUE,
                 treeID,
                 parent,
                 bootMembrIndx,
                 bootMembrSize,
                 allMembrIndx,
                 allMembrSize);
      if (RF_timeIndex > 0) {
        if (RF_mTimeFlag == TRUE) {
          updateTimeIndexArray(treeID,
                               allMembrIndx,
                               allMembrSize,
                               RF_time[treeID],
                               FALSE,
                               FALSE,
                               RF_masterTimeIndex[treeID]);
        }
      }
    }
    switch (mode) {
    case RF_PRED:
      if (RF_fmRecordSize > 0) {
        imputeNode(RF_PRED,
                   FALSE,
                   FALSE,
                   treeID,
                   parent,
                   bootMembrIndx,
                   bootMembrSize,
                   ngAllMembrIndx,
                   ngAllMembrSize);
      }
      break;
    default:
      break;
    }
  }  
  if (bootResult) {
    if (RF_opt & OPT_NODE_STAT) {
      if (RF_ptnCount == 0) {
        if (((RF_timeIndex > 0) && (RF_statusIndex > 0)) || (RF_rSize == 0)) {
          parent -> variance = NA_REAL;
        }
        else {
          getVariance(repMembrSize, repMembrIndx, 0, NULL, RF_response[treeID][1], NULL, & (parent -> variance));
        }
      }
    }
    if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
      terminalFlag = FALSE;
      leftAllMembrIndx = rghtAllMembrIndx = NULL;
      leftAllMembrSize = rghtAllMembrSize = 0;
      leftRepMembrIndx = rghtRepMembrIndx = NULL;
      leftRepMembrSize = rghtRepMembrSize = 0;
      factorFlag = FALSE;
      if (strcmp(RF_xType[parent -> splitParameter], "C") == 0) {
        factorFlag = TRUE;
      }
      if (RF_optHigh & OPT_MEMB_INCG) {
      }
      else {
        uint *membershipIndicator = uivector(1, RF_observationSize);
        leftAllMembrSize = rghtAllMembrSize = 0;
        for (i = 1; i <= allMembrSize; i++) {
          daughterFlag = RIGHT;
          if (factorFlag == TRUE) {
            daughterFlag = splitOnFactor((uint) RF_observation[treeID][parent -> splitParameter][allMembrIndx[i]], parent -> splitValueFactPtr);
          }
          else {
            if ( RF_observation[treeID][parent -> splitParameter][allMembrIndx[i]] <= (parent -> splitValueCont) ) {
              daughterFlag = LEFT;
            }
          }
          membershipIndicator[allMembrIndx[i]] = daughterFlag;
          if (daughterFlag == LEFT) {
            leftAllMembrSize ++;
            RF_nodeMembership[treeID][allMembrIndx[i]] = parent -> left;
          }
          else {
            rghtAllMembrSize ++;
            RF_nodeMembership[treeID][allMembrIndx[i]] = parent -> right;
          }
        }  
        leftAllMembrIndx  = uivector(1, leftAllMembrSize + 1);
        rghtAllMembrIndx  = uivector(1, rghtAllMembrSize + 1);
        jLeft = jRght = 0;
        for (i = 1; i <= allMembrSize; i++) {
          if (membershipIndicator[allMembrIndx[i]] == LEFT) {
            leftAllMembrIndx[++jLeft] = allMembrIndx[i];
          }
          else {
            rghtAllMembrIndx[++jRght] = allMembrIndx[i];
          }
        }
        if ( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2) ) {
          leftRepMembrIndx = leftAllMembrIndx;
          leftRepMembrSize = leftAllMembrSize;
          rghtRepMembrIndx = rghtAllMembrIndx;
          rghtRepMembrSize = rghtAllMembrSize;
        }
        else {
          leftRepMembrIndx  = uivector(1, bootMembrSize + 1);
          rghtRepMembrIndx  = uivector(1, bootMembrSize + 1);
          leftRepMembrSize = rghtRepMembrSize = 0;
          for (i = 1; i <= bootMembrSize; i++) {
            if (membershipIndicator[bootMembrIndx[i]] == LEFT) {
              leftRepMembrIndx[++leftRepMembrSize] = bootMembrIndx[i];
            }
            else {
              rghtRepMembrIndx[++rghtRepMembrSize] = bootMembrIndx[i];
            }
          }
        }
        free_uivector(membershipIndicator, 1, RF_observationSize);
      }  
      ngLeftAllMembrIndx = ngRghtAllMembrIndx = NULL;
      ngLeftAllMembrSize = ngRghtAllMembrSize = 0;
      if (mode == RF_PRED) {
        uint *ngMembershipIndicator = uivector(1, RF_fobservationSize);
        for (i=1; i <= ngAllMembrSize; i++) {
          daughterFlag = RIGHT;
          if (factorFlag == TRUE) {
            daughterFlag = splitOnFactor((uint) RF_fobservation[treeID][parent -> splitParameter][ngAllMembrIndx[i]], parent -> splitValueFactPtr);
          }
          else {
            if ( RF_fobservation[treeID][parent -> splitParameter][ngAllMembrIndx[i]] <= (parent -> splitValueCont) ) {
              daughterFlag = LEFT;
            }
          }
          if (daughterFlag == LEFT) {
            ngMembershipIndicator[ngAllMembrIndx[i]] = LEFT;
            ngLeftAllMembrSize ++;
            RF_fnodeMembership[treeID][ngAllMembrIndx[i]] = parent -> left;
          }
          else {
            ngMembershipIndicator[ngAllMembrIndx[i]] = RIGHT;
            ngRghtAllMembrSize ++;
            RF_fnodeMembership[treeID][ngAllMembrIndx[i]] = parent -> right;
          }
        }  
        ngLeftAllMembrIndx  = uivector(1, ngLeftAllMembrSize + 1);
        ngRghtAllMembrIndx  = uivector(1, ngRghtAllMembrSize + 1);
        jLeft = jRght = 0;
        for (i = 1; i <= ngAllMembrSize; i++) {
          if (ngMembershipIndicator[ngAllMembrIndx[i]] == LEFT) {
            ngLeftAllMembrIndx[++jLeft] = ngAllMembrIndx[i];
          }
          else {
            ngRghtAllMembrIndx[++jRght] = ngAllMembrIndx[i];
          }
        }
        free_uivector(ngMembershipIndicator, 1, RF_fobservationSize);
      }  
      leftResult = restoreNodeMembership(r,
                                         mode,
                                         FALSE,
                                         treeID,
                                         parent -> left,
                                         leftRepMembrIndx,
                                         leftRepMembrSize,
                                         leftAllMembrIndx,
                                         leftAllMembrSize,
                                         ngLeftAllMembrIndx,
                                         ngLeftAllMembrSize,
                                         bootMembrIndxIter,
                                         rmbrIterator,
                                         ambrIterator);
      if(!leftResult) {
      }
      rghtResult = restoreNodeMembership(r,
                                         mode,
                                         FALSE,
                                         treeID,
                                         parent -> right,
                                         rghtRepMembrIndx,
                                         rghtRepMembrSize,
                                         rghtAllMembrIndx,
                                         rghtAllMembrSize,
                                         ngRghtAllMembrIndx,
                                         ngRghtAllMembrSize,
                                         bootMembrIndxIter,
                                         rmbrIterator,
                                         ambrIterator);
      if(!rghtResult) {
      }
      if (RF_optHigh & OPT_MEMB_INCG) {
      }
      else {
        free_uivector(leftAllMembrIndx, 1, leftAllMembrSize + 1);
        free_uivector(rghtAllMembrIndx, 1, rghtAllMembrSize + 1);
        if ( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2) ) {
        }
        else {
          free_uivector(leftRepMembrIndx, 1, bootMembrSize + 1);
          free_uivector(rghtRepMembrIndx, 1, bootMembrSize + 1);
        }
      }
      if (mode == RF_PRED) {
        free_uivector(ngLeftAllMembrIndx, 1, ngLeftAllMembrSize + 1);
        free_uivector(ngRghtAllMembrIndx, 1, ngRghtAllMembrSize + 1);
      }
    }  
    else {
    }
  }  
  else {
    if (rootFlag) {
      if (!bootResult) {
        terminalFlag = FALSE;
      }
    }
  }   
  if (terminalFlag) {
    if ((RF_mRecordSize > 0) || (RF_fmRecordSize > 0)) {
      imputeNodeAndSummarize(r,
                             mode,
                             treeID,
                             parent,
                             bootMembrIndx,
                             bootMembrSize,
                             allMembrIndx,
                             allMembrSize,
                             ngAllMembrIndx,
                             ngAllMembrSize);
    }
    if (mode == RF_PRED) {
      if (ngAllMembrSize > 0) {
        for (i = 1; i <= ngAllMembrSize; i++) {
          RF_ftTermMembership[treeID][ngAllMembrIndx[i]] = RF_tTermList[treeID][parent -> nodeID];
        }
      }
    }
    if (RF_optHigh & OPT_MEMB_USER) {
      if (mode == RF_PRED) {      
        if (ngAllMembrSize > 0) {
          for (i = 1; i <= ngAllMembrSize; i++) {
            RF_MEMB_ID_ptr[treeID][ngAllMembrIndx[i]] = parent -> nodeID;
          }
        }
      }
      else {
        if (RF_optHigh & OPT_MEMB_INCG) {
          uint userIterator = *ambrIterator;
          for (i = 1; i <= RF_TN_ACNT_ptr[treeID][parent -> nodeID]; i++) {
            ++(userIterator);
            RF_MEMB_ID_ptr[treeID][RF_AMBR_ID_ptr[treeID][(userIterator)]] = parent -> nodeID;
          }
        }
        else {
          for (i = 1; i <= allMembrSize; i++) {
            RF_MEMB_ID_ptr[treeID][allMembrIndx[i]] = parent -> nodeID;
          }
        }
      }
    }
    updateTerminalNodeOutcomesNew(mode,
                                  treeID,
                                  RF_tTermList[treeID][parent -> nodeID],
                                  bootMembrIndx,
                                  bootMembrSize,
                                  allMembrIndx,
                                  allMembrSize,
                                  rmbrIterator,
                                  ambrIterator);
    if ( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2) ) {
      bsUpdateFlag = TRUE;
    }
  }  
  if (bsUpdateFlag) {
    for (i = 1; i <= bootMembrSize; i++) {
      RF_bootMembershipIndex[treeID][++(*bootMembrIndxIter)] = bootMembrIndx[i];
      RF_bootMembershipFlag[treeID][bootMembrIndx[i]] = TRUE;
      RF_oobMembershipFlag[treeID][bootMembrIndx[i]]  = FALSE;
      RF_bootMembershipCount[treeID][bootMembrIndx[i]] ++;
      if (RF_optHigh & OPT_MEMB_USER) {
        RF_BOOT_CT_ptr[treeID][bootMembrIndx[i]] ++;
      }
    }
  }
  if (rootFlag || ((RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2))) {
    if ( (!(RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
         ( (RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ) {
      free_uivector(bootMembrIndx, 1, RF_bootstrapSize);
    }
    else {
      free_uivector(bootMembrIndx, 1, allMembrSize);
    }
  }
  return bootResult;
}  
void imputeNodeAndSummarize(uint     r,
                            char     mode,
                            uint     treeID,
                            Node    *parent,
                            uint    *repMembrIndx,
                            uint     repMembrSize,
                            uint    *allMembrIndx,
                            uint     allMembrSize,
                            uint    *ngAllMembrIndx,
                            uint     ngAllMembrSize) {
  uint multImpFlag;
  if (r == 1) {
    if (RF_mRecordSize > 0) {
      unstackNodeLMPIndex(RF_tNodeList[treeID][parent -> nodeID]);
      imputeNode(RF_GROW, 
                 TRUE,    
                 TRUE,    
                 treeID,
                 parent,
                 repMembrIndx,
                 repMembrSize,
                 allMembrIndx,
                 allMembrSize);
      if (mode != RF_PRED) {
        xferMissingness(RF_GROW, RF_tNodeList[treeID][parent -> nodeID], RF_tTermList[treeID][parent -> nodeID]);
      }
    }
    if (mode == RF_PRED) {
      if (RF_fmRecordSize > 0) {
        unstackNodeFLMPIndex(RF_tNodeList[treeID][parent -> nodeID]);
        imputeNode(RF_PRED, 
                   TRUE,    
                   FALSE,   
                   treeID,
                   parent,
                   repMembrIndx,
                   repMembrSize,
                   ngAllMembrIndx,
                   ngAllMembrSize);
        xferMissingness(RF_PRED, RF_tNodeList[treeID][parent -> nodeID], RF_tTermList[treeID][parent -> nodeID]);
      }
    }
  }
  else {
    multImpFlag = FALSE;
    if (r < RF_nImpute) {
      multImpFlag = TRUE;
    }
    else {
      if (RF_opt & OPT_IMPU_ONLY) {
        multImpFlag = TRUE;
      }
    }
    if (multImpFlag) {
      if (RF_mRecordSize > 0) {
        imputeNode(RF_GROW, 
                   TRUE,    
                   FALSE,   
                   treeID,
                   parent,
                   repMembrIndx,
                   repMembrSize,
                   allMembrIndx,
                   allMembrSize);
        xferMissingness(RF_GROW, RF_tNodeList[treeID][parent -> nodeID], RF_tTermList[treeID][parent -> nodeID]);
      }
    }
  }
}
void imputeUpdateShadow (char      mode,
                         double  **shadowResponse,
                         double  **shadowPredictor) {
  uint     mRecordSize;
  uint    *mRecordIndex;
  uint     mpIndexSize;
  int    **mpSign;
  int     *mpIndex;
  double **outResponse;
  double **outPredictor;
  double  *valuePtr;
  double  *outputPtr;
  uint unsignedIndex;
  char outcomeFlag, predictorFlag;
  uint rspSize;
  uint i, p;
  mRecordSize  = 0;     
  mRecordIndex = NULL;  
  mpIndexSize  = 0;     
  mpSign       = NULL;  
  mpIndex      = NULL;  
  outResponse  = NULL;  
  outPredictor = NULL;  
  valuePtr     = NULL;  
  outputPtr    = NULL;  
  unsignedIndex = 0;    
  switch (mode) {
  case RF_PRED:
    mRecordSize = RF_fmRecordSize;
    mRecordIndex = RF_fmRecordIndex;
    mpIndexSize = RF_fmpIndexSize;
    mpSign = RF_fmpSign;
    mpIndex = RF_fmpIndex;
    if (shadowResponse != NULL) {
      outResponse  = RF_sImputeResponsePtr;
    }
    if (shadowPredictor != NULL) {
      outPredictor = RF_sImputePredictorPtr;
    }
    rspSize = RF_frSize;
    break;
  default:
    mRecordSize = RF_mRecordSize;
    mRecordIndex = RF_mRecordIndex;
    mpIndexSize = RF_mpIndexSize;
    mpSign = RF_mpSign;
    mpIndex = RF_mpIndex;
    if (shadowResponse != NULL) {
      outResponse  = RF_sImputeResponsePtr;
    }
    if (shadowPredictor != NULL) {
      outPredictor = RF_sImputePredictorPtr;
    }
    rspSize = RF_rSize;
    break;
  }
  if (mRecordSize == 0) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Attempt to update shadow data with no missingness in mode:  %10d", mode);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  for (p = 1; p <= mpIndexSize; p++) {
    if (mpIndex[p] < 0) {
      if (shadowResponse != NULL) {
        unsignedIndex = (uint) abs(mpIndex[p]);
        valuePtr = shadowResponse[(uint) abs(mpIndex[p])];
        outputPtr = outResponse[(uint) abs(mpIndex[p])];
        outcomeFlag = TRUE;
      }
      else {
        outcomeFlag = FALSE;
      }
      predictorFlag = FALSE;
    }
    else {
      if (shadowPredictor != NULL) {
        unsignedIndex = (uint) mpIndex[p] + rspSize;
        valuePtr = shadowPredictor[(uint) mpIndex[p]];
        outputPtr = outPredictor[(uint) mpIndex[p]];
        predictorFlag = TRUE;
      }
      else {
        predictorFlag = FALSE;
      }
      outcomeFlag = FALSE;
    }
    if ( (outcomeFlag && (shadowResponse != NULL))  || (predictorFlag && (shadowPredictor != NULL)) ) {
      for (i = 1; i <= mRecordSize; i++) {
        if (mpSign[unsignedIndex][i] == 1) {
          if (ISNA(outputPtr[i])) {
          }
          valuePtr[mRecordIndex[i]] = outputPtr[i];
        }
      }
    }  
  }  
}
void imputeSummary(char      mode,
                   char      selectionFlag) {
  imputeCommon(mode,
               0,
               selectionFlag,
               TRUE);
}
void imputeResponse(char      mode,
                    uint      serialTreeID,
                    double  **tempResponse) {
  switch(mode) {
  case RF_PRED:
    imputeCommon(mode, serialTreeID, ACTIVE, FALSE);
    imputeUpdateShadow(mode, tempResponse, NULL);
    break;
  default:
    imputeCommon(mode, serialTreeID, FALSE, FALSE);
    imputeUpdateShadow(mode, tempResponse, NULL);
    break;
  }
}
void imputeCommon(char      mode,
                  uint      serialTreeID,
                  char      selectionFlag,
                  char      predictorFlag) {
  uint *localSerialIndex;
  uint  localSerialCount;
  uint *serialPtr;
  char mFlag;
  char outcomeFlag;
  uint     mRecordSize;
  uint    *mRecordIndex;
  uint     mpIndexSize;
  int    **mpSign;
  int     *mpIndex;
  double **outResponse;
  double **outPredictor;
  double *valuePtr;
  double *naivePtr;
  uint    unsignedSignatureIndex;
  Terminal ***termMembershipPtr;
  Terminal *info;
  double imputedValue;
  uint localDistributionSize;
  uint maxDistributionSize;
  uint rspSize;
  char result;
  uint i, p, v, tree;
  valuePtr      = NULL;  
  naivePtr      = NULL;  
  unsignedSignatureIndex = 0;     
  maxDistributionSize = 0;  
  outResponse         = 0;  
  outPredictor        = 0;  
  rspSize = 0;  
  mpIndex = 0;  
  mpSign  = 0;  
  mpIndexSize  = 0;  
  mRecordIndex = 0;  
  mRecordSize  = 0;  
  localSerialIndex = NULL;  
  if ((selectionFlag != TRUE) && (selectionFlag != FALSE) && (selectionFlag != ACTIVE)) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Invalid selectionFlag in imputeCommon():  %10d", selectionFlag);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  result = FALSE;
  switch (mode) {
  case RF_PRED:
    if (RF_fmRecordSize > 0) {
      mRecordSize = RF_fmRecordSize;
      mRecordIndex = RF_fmRecordIndex;
      mpIndexSize = RF_fmpIndexSize;
      mpSign = RF_fmpSign;
      mpIndex = RF_fmpIndex;
      maxDistributionSize = ((RF_observationSize) > (RF_forestSize)) ? (RF_observationSize) : (RF_forestSize);
      outResponse  = RF_sImputeResponsePtr;
      outPredictor = RF_sImputePredictorPtr;
      rspSize = RF_frSize;
      termMembershipPtr = RF_ftTermMembership;
      result = TRUE;
    }
    break;
  default:
    if (RF_mRecordSize > 0) {
      mRecordSize = RF_mRecordSize;
      mRecordIndex = RF_mRecordIndex;
      mpIndexSize = RF_mpIndexSize;
      mpSign = RF_mpSign;
      mpIndex = RF_mpIndex;
      maxDistributionSize = ((RF_observationSize) > (RF_forestSize)) ? (RF_observationSize) : (RF_forestSize);
      outResponse  = RF_sImputeResponsePtr;
      outPredictor = RF_sImputePredictorPtr;
      rspSize = RF_rSize;
      termMembershipPtr = RF_tTermMembership;
      result = TRUE;
    }
    break;
  }
  if (result == FALSE) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Attempt to impute in imputeCommon() with no missingness in mode:  %10d", mode);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (serialTreeID == 0) {
    localSerialIndex = uivector(1, RF_forestSize);
    for (tree = 1; tree <= RF_forestSize; tree++) {
      localSerialIndex[tree] = tree;
    }
    serialPtr = localSerialIndex;
    localSerialCount = RF_forestSize;
  }
  else {
    serialPtr = RF_serialTreeIndex;
    localSerialCount = serialTreeID;
  }
  imputedValue = 0.0;  
  double *localDistribution = dvector(1, maxDistributionSize);
  char  *naiveFlag = cvector(1, mpIndexSize);
  char **naiveSign = cmatrix(1, mRecordSize, 1, mpIndexSize);
  for (p = 1; p <= mpIndexSize; p++) {
    naiveFlag[p] = FALSE;
  }
  for (i = 1; i <= mRecordSize; i++) {
    outcomeFlag = TRUE;
    for (p = 1; p <= mpIndexSize; p++) {
      naiveSign[i][p] = FALSE;
      if (mpIndex[p] < 0) {
        unsignedSignatureIndex = (uint) abs(mpIndex[p]);
      }
      else {
        if (predictorFlag == TRUE) {
          unsignedSignatureIndex = (uint) mpIndex[p] + rspSize;
        }
        outcomeFlag = FALSE;
      }
      if (outcomeFlag || predictorFlag) {
        if (mpSign[unsignedSignatureIndex][i] == 1) {
          localDistributionSize = 0;
          for (tree = 1; tree <= localSerialCount; tree++) {
            if (RF_tLeafCount[serialPtr[tree]] > 0) {
              if ((RF_dmRecordBootFlag[serialPtr[tree]][i] == selectionFlag) || (selectionFlag == ACTIVE)) {
                info = termMembershipPtr[serialPtr[tree]][mRecordIndex[i]];
                for (v = 1; v <= info -> lmiSize; v++) {
                  if ((info -> lmiIndex)[v] == p) {
                        if (!ISNA((info -> lmiValue)[v])) {
                          localDistribution[++localDistributionSize] = (info -> lmiValue)[v];
                        }
                        else {
                        }  
                    v = info -> lmiSize;
                  }
                }
              }  
            }  
            else {
            }
          }  
          if (localDistributionSize > 0) {
            if (mpIndex[p] < 0) {
              if (strcmp(RF_rType[(uint) abs(mpIndex[p])], "T") == 0) {
                imputedValue = getMeanValue(localDistribution, localDistributionSize);
              }
              else if (strcmp(RF_rType[(uint) abs(mpIndex[p])], "S") == 0) {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize, FALSE, localSerialCount);
              }
              else if (strcmp(RF_rType[(uint) abs(mpIndex[p])], "R") == 0) {
                imputedValue = getMeanValue(localDistribution, localDistributionSize);
              }
              else if (strcmp(RF_rType[(uint) abs(mpIndex[p])], "I") == 0) {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize, FALSE, localSerialCount);
              }
              else if (strcmp(RF_rType[(uint) abs(mpIndex[p])], "C") == 0) {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize, FALSE, localSerialCount);
              }
              outResponse[(uint) abs(mpIndex[p])][i] = imputedValue;
            }  
            else {
              if (strcmp(RF_xType[(uint) mpIndex[p]], "R") == 0) {
                imputedValue = getMeanValue(localDistribution, localDistributionSize);
              }
              else {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize, FALSE, localSerialCount);
              }
              outPredictor[(uint) mpIndex[p]][i] = imputedValue;
            }
          }  
          else {
            naiveFlag[p] = TRUE;
            naiveSign[i][p] = TRUE;
          }
        }  
      }  
      else {
        p = mpIndexSize;
      }
    }  
  }  
  outcomeFlag = TRUE;
  for (p = 1; p <= mpIndexSize; p++) {
    if (mpIndex[p] < 0) {
      unsignedSignatureIndex = (uint) abs(mpIndex[p]);
      valuePtr = RF_responseIn[(uint) abs(mpIndex[p])];
      naivePtr = outResponse[(uint) abs(mpIndex[p])];
    }
    else {
      if (predictorFlag == TRUE) {
        unsignedSignatureIndex = (uint) mpIndex[p] + rspSize;
        valuePtr = RF_observationIn[(uint) mpIndex[p]];
        naivePtr = outPredictor[(uint) mpIndex[p]];
      }
      outcomeFlag = FALSE;
    }
    if (outcomeFlag || predictorFlag) {
      if (naiveFlag[p] == TRUE) {
        localDistributionSize = 0;
        for (i=1; i <= RF_observationSize; i++) {
          mFlag = TRUE;
          if (RF_mRecordMap[i] == 0) {
            mFlag = FALSE;
          }
          else if (RF_mpSign[unsignedSignatureIndex][RF_mRecordMap[i]] == 0) {
            mFlag = FALSE;
          }
          if (mFlag == FALSE) {
            localDistribution[++localDistributionSize] = valuePtr[i];
          }
        }  
        if (localDistributionSize > 0) {
          for (i=1; i <= mRecordSize; i++) {
            if (naiveSign[i][p] == TRUE) {
              naivePtr[i] = getSampleValue(localDistribution, localDistributionSize, FALSE, localSerialCount);
            }
          }
        }  
        else {
          if (mpIndex[p] < 0) {
            RFprintf("\nRF-SRC:  *** ERROR *** ");
            RFprintf("\nRF-SRC:  Naive imputation failed for [indv, outcome] = [%10d, %10d] \n", mRecordIndex[i], mpIndex[p]);
            RFprintf("\nRF-SRC:  Please Contact Technical Support.");
            error("\nRF-SRC:  The application will now exit.\n");
          }
          else {
          }
        }
      }  
    }  
    else {
      p = mpIndexSize;
    }
  }  
  if (serialTreeID == 0) {
    free_uivector(localSerialIndex, 1, RF_forestSize);
  }
  free_dvector(localDistribution, 1, maxDistributionSize);
  free_cvector(naiveFlag, 1, mpIndexSize);
  free_cmatrix(naiveSign, 1, mRecordSize, 1, mpIndexSize);
}
void imputeMultipleTime (char selectionFlag) {
  double  *outTime;
  char     result;
  uint i;
  result = FALSE;
    if (RF_timeIndex > 0) {
      if (RF_mRecordSize > 0) {
      if (RF_mTimeFlag == TRUE) {
        result = TRUE;
      }
      else {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Attempt to multiply impute time with no missingness in time vector.");
      }
    }
  }
  else {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Attempt to multiply impute time in a !SURV environment.");
  }
  if (result == FALSE) {
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  outTime  = RF_sImputeResponsePtr[RF_timeIndex];
  for (i=1; i <= RF_mRecordSize; i++) {
    if(RF_mpSign[RF_timeIndex][i] == 1) {
      outTime[i] = getNearestMasterTime(outTime[i], FALSE, 1);
    }
  }
}
double getNearestMasterTime (double   meanValue,
                             char     chainFlag,
                             uint     treeID) {
  double leftDistance, rightDistance;
  uint minimumIndex;
  uint j;
    leftDistance = meanValue - RF_masterTime[1];
    rightDistance = RF_masterTime[RF_masterTimeSize] - meanValue;
    if ( ((leftDistance > EPSILON) || (fabs(leftDistance) < EPSILON)) &&
         ((rightDistance > EPSILON) || (fabs(rightDistance) < EPSILON)) ) {
    }
    else {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  The summary mean value for time is out of range:  %12.4f <= %12.4f <= %12.4f", RF_masterTime[1], meanValue, RF_masterTime[RF_masterTimeSize]);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  leftDistance = rightDistance = 0;
  minimumIndex = RF_masterTimeSize;
  for (j = 1; j <= RF_masterTimeSize; j++) {
    if (meanValue <= RF_masterTime[j]) {
      minimumIndex = j;
      j = RF_masterTimeSize;
    }
  }
  if (minimumIndex == 1) {
  }
  else {
    leftDistance = meanValue - RF_masterTime[minimumIndex-1];
    rightDistance = RF_masterTime[minimumIndex] - meanValue;
    if (leftDistance < rightDistance) {
      minimumIndex = minimumIndex - 1;
    }
    else {
      if (fabs(leftDistance - rightDistance) < EPSILON) {
        if(chainFlag) {
          if (ran1A(treeID) <= 0.5) {
            minimumIndex = minimumIndex - 1;
          }
        }
        else {
          if (ran1B(treeID) <= 0.5) {
            minimumIndex = minimumIndex - 1;
          }
        }
      }
    }
  }
  return RF_masterTime[minimumIndex];
}
double getMaximalValue(double *value, uint size, char chainFlag, uint treeID) {
  double result;
  uint classCount, maximalClassSize, maximalClassCount;
  uint randomIndex;
  uint j;
  uint   *classSize  = uivector(1, size);
  for (j = 1; j <= size; j++) {
    classSize[j] = 0;
  }
  hpsort(value, size);
  classCount = 1;
  classSize[1] = 1;
  for (j = 2; j <= size; j++) {
    if (value[j] > value[classCount]) {
      classCount ++;
      value[classCount] = value[j];
    }
    classSize[classCount] ++;
  }
  maximalClassSize = maximalClassCount = 0;
  for (j=1; j <= classCount; j++) {
    if (classSize[j] > maximalClassSize) {
      maximalClassSize = classSize[j];
    }
  }
  for (j=1; j <= classCount; j++) {
    if (classSize[j] == maximalClassSize) {
      maximalClassCount ++;
    }
  }
  if (maximalClassCount > 1) {
    if(chainFlag) {
      randomIndex = (uint) ceil(ran1A(treeID)*((maximalClassCount)*1.0));
    }
    else {
      randomIndex = (uint) ceil(ran1B(treeID)*((maximalClassCount)*1.0));
    }
  }
  else {
    randomIndex = 1;
  }
  j = 0;
  while (randomIndex > 0) {
    j++;
    if (classSize[j] == maximalClassSize) {
      randomIndex --;
    }
  }
  result = value[j];
  free_uivector(classSize, 1, size);
  return result;
}
double getMedianValue(double *value, uint size) {
  double result;
  uint medianIndex;
  sort(value, size);
  if (size > 1) {
    medianIndex = (uint) ceil(size/2);
  }
  else {
    medianIndex = 1;
  }
  result = value[medianIndex];
  return result;
}
double getMeanValue(double *value, uint size) {
  double result;
  uint j;
  result = 0.0;
  for (j = 1; j <= size; j++) {
    result = result + value[j];
  }
  result = result / size;
  return result;
}
double getSampleValue(double *value, uint size, char chainFlag, uint treeID) {
  uint randomIndex;
  if(chainFlag) {
    randomIndex = (uint) ceil(ran1A(treeID)*((size)*1.0));
  }
  else {
    randomIndex = (uint) ceil(ran1B(treeID)*((size)*1.0));
  }
  return value[randomIndex];
}
uint getRecordMap(uint    *map,
                  uint     obsSize,
                  double **resp,
                  double **data) {
  uint i, p, r;
  uint mSize;
  char mFlag;
  mSize  = 0;
  for (i = 1; i <= obsSize; i++) {
    mFlag = FALSE;
    if (resp != NULL) {
      for (r = 1; r <= RF_rSize; r++) {
        if (ISNA(resp[r][i])) {
          mFlag = TRUE;
          r = RF_rSize;
        }
      }
    }
    if (mFlag == FALSE) {
      for (p = 1; p <= RF_xSize; p++) {
        if (ISNA(data[p][i])) {
          mFlag = TRUE;
          p = RF_xSize;
        }
      }
    }
    if (mFlag == TRUE) {
      mSize ++;
      map[i] = mSize;
    }
    else {
      map[i] = 0;
    }
  }
  return mSize;
}
void updateTimeIndexArray(uint    treeID,
                          uint   *allMembrIndx,
                          uint    allMembrSize,
                          double *time,
                          char    naAllowFlag,
                          char    noIdxAllowFlag,
                          uint   *masterTimeIndex) {
  uint *membrIndx;
  char idxFoundFlag;
  uint i,k;
  if (allMembrIndx == NULL) {
    membrIndx = uivector(1, allMembrSize);
    for (i = 1; i <= allMembrSize; i++) {
      membrIndx[i] = i;
    }
  }
  else {
    membrIndx = allMembrIndx;
  }
  for (i=1; i <= allMembrSize; i++) {
    idxFoundFlag = FALSE;
    if (!ISNA(time[membrIndx[i]])) {
      k = 1;
      while (k <= RF_masterTimeSize) {
        if (time[membrIndx[i]] == RF_masterTime[k]) {
          masterTimeIndex[membrIndx[i]] = k;
          idxFoundFlag = TRUE;
          k = RF_masterTimeSize;
        }
        k++;
      }
    }
    else {
      if (naAllowFlag == FALSE) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Missing event time encountered for individual:  %10d, %12.4f", i, time[membrIndx[i]]);
        RFprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
      else {
        masterTimeIndex[membrIndx[i]] = 0;
        idxFoundFlag = TRUE;
      }
    }
    if (idxFoundFlag == FALSE) {
      if (noIdxAllowFlag == FALSE) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Invalid event time encountered for individual:  %10d, %12.4f", i, time[membrIndx[i]]);
        RFprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
      else {
        masterTimeIndex[membrIndx[i]] = 0;
      }
    }
  }
  if (allMembrIndx == NULL) {
    free_uivector(membrIndx, 1, allMembrSize);
  }
}  
void updateEventTypeSubsets(double *summaryStatus,
                            uint    mRecordSize,
                            int   **mpSign,
                            uint   *mRecordIndex,
                            uint   *meIndividualSize,
                            uint  **eIndividual) {
  uint i, j;
  if (RF_eventTypeSize == 1) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Attempt to update event type subsets in a non-CR analysis.");
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  for (j = 1; j <= RF_eventTypeSize; j++) {
    for (i = 1; i <= RF_eIndividualSize[j]; i++) {
      eIndividual[j][i] = RF_eIndividualIn[j][i];
    }
  }
  if (RF_mStatusSize > 0) {
    uint *eventCounter = uivector(1, RF_eventTypeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      eventCounter[j] = RF_eIndividualSize[j];
    }
    for (i = 1; i <= mRecordSize; i++) {
      if (mpSign[RF_statusIndex][i] == 1) {
        if ((uint) summaryStatus[mRecordIndex[i]] > 0) {
          j = RF_eventTypeIndex[(uint) summaryStatus[mRecordIndex[i]]];
          eventCounter[j] ++;
          eIndividual[j][eventCounter[j]] = mRecordIndex[i];
        }
        else {
          for (j=1; j <= RF_eventTypeSize; j++) {
            eventCounter[j] ++;
            eIndividual[j][eventCounter[j]] = mRecordIndex[i];
          }
        }
      }
    }
    for (j = 1; j <= RF_eventTypeSize; j++) {
      meIndividualSize[j] = eventCounter[j];
    }
    free_uivector(eventCounter, 1, RF_eventTypeSize);
  }
}
void stackShadow (char mode, uint treeID) {
  uint *nonMissIndex;
  uint *permuteIndex;
  uint * permuteSize;
  char *nullSplitShadowFlag;
  char  vimpShadowFlag;
  uint  unsignedIndexSource;
  uint  mpIndexIter;
  uint i, j, p;
  nonMissIndex        = NULL;  
  permuteIndex        = NULL;  
  permuteSize         = NULL;  
  nullSplitShadowFlag = NULL;  
  if (RF_opt & OPT_SPLT_NULL) {
    nullSplitShadowFlag = cvector(1, RF_rSize);
    for (p = 1; p <= RF_rSize; p++) {
      nullSplitShadowFlag[p] = TRUE;
    }
    RF_response[treeID] = (double **) new_vvector(1, RF_rSize, NRUTIL_DPTR);
    nonMissIndex = uivector(1, RF_observationSize);
    permuteIndex = uivector(1, RF_observationSize);
    permuteSize  = uivector(1, RF_rSize);
    for (p = 1; p <= RF_rSize; p++) {
      RF_response[treeID][p] = dvector(1, RF_observationSize);
      for (i = 1; i <= RF_observationSize; i++) {
        RF_response[treeID][p][i] = RF_responseIn[p][i];
      }
    }
    mpIndexIter = 1;
    for (p = 1; p <= RF_rSize; p++) {
      permuteSize[p] = RF_observationSize;
      for (i = 1; i <= RF_observationSize; i++) {
        nonMissIndex[i] = i;
      }
      if (RF_mRecordSize > 0) {
        if (RF_mpIndex[mpIndexIter] < 0) {
          unsignedIndexSource = (uint) abs(RF_mpIndex[p]);
          if (unsignedIndexSource == p) {
            permuteSize[unsignedIndexSource] = 0;
            for (i = 1; i <= RF_observationSize; i++) {
              if (RF_mRecordMap[i] == 0) {
                nonMissIndex[++(permuteSize[unsignedIndexSource])] = i;
              }
              else {
                if (RF_mpSign[unsignedIndexSource][RF_mRecordMap[i]] == 0) {
                  nonMissIndex[++(permuteSize[unsignedIndexSource])] = i;
                }
              }
            }
            mpIndexIter++;
          }
        }
      }  
      if(nullSplitShadowFlag[p]) {
        permute (1, treeID, permuteSize[p], permuteIndex);
        for (i = 1; i <= permuteSize[p]; i++) {
          RF_response[treeID][p][nonMissIndex[i]] = RF_responseIn[p][nonMissIndex[permuteIndex[i]]];
        }
      }
    }
    if (RF_timeIndex > 0) {
      RF_time[treeID] = RF_response[treeID][RF_timeIndex];
      RF_masterTimeIndex[treeID] = uivector(1, RF_observationSize);
      updateTimeIndexArray(treeID,
                           NULL,
                           RF_observationSize,
                           RF_time[treeID],
                           TRUE,
                           FALSE,
                           RF_masterTimeIndex[treeID]);
    }
    if (RF_statusIndex > 0) {
      RF_status[treeID] =  RF_response[treeID][RF_statusIndex];
    }
  }
  else {
    if (RF_mResponseFlag == TRUE) {
      RF_response[treeID] = (double **) new_vvector(1, RF_rSize, NRUTIL_DPTR);
      for (p = 1; p <= RF_rSize; p++) {
        RF_response[treeID][p] = RF_responseIn[p];
      }
      for (p = 1; p <= RF_mpIndexSize; p++) {
        if (RF_mpIndex[p] < 0) {
          RF_response[treeID][(uint) abs(RF_mpIndex[p])] = dvector(1, RF_observationSize);
          for (i = 1; i <= RF_observationSize; i++) {
            RF_response[treeID][(uint) abs(RF_mpIndex[p])][i] = RF_responseIn[(uint) abs(RF_mpIndex[p])][i];
          }
        }
        else {
          p = RF_mpIndexSize;
        }
      }
      if (RF_timeIndex > 0) {
        RF_time[treeID] = RF_response[treeID][RF_timeIndex];
        if (RF_mTimeFlag == TRUE) {
          RF_masterTimeIndex[treeID] = uivector(1, RF_observationSize);
          for (i = 1; i <= RF_observationSize; i++) {
            RF_masterTimeIndex[treeID][i] = RF_masterTimeIndexIn[i];
          }
        }
        else {
          RF_masterTimeIndex[treeID] = RF_masterTimeIndexIn;
        }
      }
      if (RF_statusIndex > 0) {
        RF_status[treeID] =  RF_response[treeID][RF_statusIndex];
      }
    }
  }
  if (RF_opt & OPT_SPLT_NULL) {
    free_uivector(nonMissIndex, 1, RF_observationSize);
    free_uivector(permuteIndex, 1, RF_observationSize);
    free_uivector(permuteSize,  1, RF_rSize);
    free_cvector(nullSplitShadowFlag, 1, RF_rSize);
  }
  if (mode == RF_PRED) {
    if (RF_frSize > 0) {
      if (RF_fmResponseFlag == TRUE) {
        RF_fresponse[treeID] = (double **) new_vvector(1, RF_rSize, NRUTIL_DPTR);
        for (p = 1; p <= RF_frSize; p++) {
          RF_fresponse[treeID][p] = RF_fresponseIn[p];
        }
        for (p = 1; p <= RF_fmpIndexSize; p++) {
          if (RF_fmpIndex[p] < 0) {
            RF_fresponse[treeID][(uint) abs(RF_fmpIndex[p])] = dvector(1, RF_fobservationSize);
            for (i = 1; i <= RF_fobservationSize; i++) {
              RF_fresponse[treeID][(uint) abs(RF_fmpIndex[p])][i] = RF_fresponseIn[(uint) abs(RF_fmpIndex[p])][i];
            }
          }
          else {
            p = RF_fmpIndexSize;
          }
        }
      }
    }
  }
  if (RF_rFactorCount + RF_xFactorCount > 0) {
    RF_factorList[treeID] = (Factor **) new_vvector(1, RF_maxFactorLevel, NRUTIL_FPTR);
    for (j = 1; j <= RF_maxFactorLevel; j++) {
      RF_factorList[treeID][j] = NULL;
    }
    for (j = 1; j <= RF_xFactorCount; j++) {
      if (RF_factorList[treeID][RF_xFactorSize[j]] == NULL) {
        RF_factorList[treeID][RF_xFactorSize[j]] = makeFactor(RF_xFactorSize[j], FALSE);
      }
    }
    for (j = 1; j <= RF_rFactorCount; j++) {
      if (RF_factorList[treeID][RF_rFactorSize[j]] == NULL) {
        RF_factorList[treeID][RF_rFactorSize[j]] = makeFactor(RF_rFactorSize[j], FALSE);
      }
    }
  }
  vimpShadowFlag = FALSE;
  if ((RF_opt & OPT_VIMP) && (RF_opt & OPT_VIMP_TYP1) && !(RF_opt & OPT_VIMP_TYP2)) {
    vimpShadowFlag = TRUE;
  }
  if(vimpShadowFlag == TRUE) {
    RF_observation[treeID] = dmatrix(1, RF_xSize, 1, RF_observationSize);
    for (p = 1; p <= RF_xSize; p++) {
      for (i = 1; i <= RF_observationSize; i++) {
        RF_observation[treeID][p][i] = RF_observationIn[p][i];
      }
    }
  }
  else {
    if(RF_mPredictorFlag == TRUE) {
      RF_observation[treeID] = (double **) new_vvector(1, RF_xSize, NRUTIL_DPTR);
      for (p = 1; p <= RF_xSize; p++) {
        RF_observation[treeID][p] = RF_observationIn[p];
      }
      for (p = 1; p <= RF_mpIndexSize; p++) {
        if (RF_mpIndex[p] > 0) {
          RF_observation[treeID][(uint) RF_mpIndex[p]] = dvector(1, RF_observationSize);
          for (i = 1; i <= RF_observationSize; i++) {
            RF_observation[treeID][(uint) RF_mpIndex[p]][i] = RF_observationIn[(uint) RF_mpIndex[p]][i];
          }
        }
      }
    }
  }
  if (mode == RF_PRED) {
    if(vimpShadowFlag == TRUE) {
      RF_fobservation[treeID] = dmatrix(1, RF_xSize, 1, RF_fobservationSize);
      for (p = 1; p <= RF_xSize; p++) {
        for (i = 1; i <= RF_fobservationSize; i++) {
          RF_fobservation[treeID][p][i] = RF_fobservationIn[p][i];
        }
      }
    }
    else {
      if(RF_fmPredictorFlag == TRUE) {
        RF_fobservation[treeID] = (double **) new_vvector(1, RF_xSize, NRUTIL_DPTR);
        for (p = 1; p <= RF_xSize; p++) {
          RF_fobservation[treeID][p] = RF_fobservationIn[p];
        }
        for (p = 1; p <= RF_fmpIndexSize; p++) {
          if (RF_fmpIndex[p] > 0) {
            RF_fobservation[treeID][(uint) RF_fmpIndex[p]] = dvector(1, RF_fobservationSize);
            for (i = 1; i <= RF_fobservationSize; i++) {
              RF_fobservation[treeID][(uint) RF_fmpIndex[p]][i] = RF_fobservationIn[(uint) RF_fmpIndex[p]][i];
            }
          }
        }
      }
    }
  }  
}
void unstackShadow (char mode, uint treeID, char respFlag, char covrFlag) {
  char vimpShadowFlag;
  uint k, p;
  if (respFlag) {
    if (RF_opt & OPT_SPLT_NULL) {
      for (p = 1; p <= RF_rSize; p++) {
        free_dvector(RF_response[treeID][p], 1, RF_observationSize);
      }
      free_new_vvector(RF_response[treeID], 1, RF_rSize, NRUTIL_DPTR);
      if (RF_timeIndex > 0) {
        free_uivector(RF_masterTimeIndex[treeID], 1, RF_observationSize);
      }
    }
    else {
      if (RF_mResponseFlag == TRUE) {
        for (p = 1; p <= RF_mpIndexSize; p++) {
          if (RF_mpIndex[p] < 0) {
            free_dvector(RF_response[treeID][(uint) abs(RF_mpIndex[p])], 1, RF_observationSize);
          }
          else {
            p = RF_mpIndexSize;
          }
        }
        free_new_vvector(RF_response[treeID], 1, RF_rSize, NRUTIL_DPTR);
        if (RF_timeIndex > 0) {
          if (RF_mTimeFlag == TRUE) {
            free_uivector(RF_masterTimeIndex[treeID], 1, RF_observationSize);
          }
        }
      }
    }
    if (mode == RF_PRED) {
      if (RF_frSize > 0) {
        if (RF_fmResponseFlag == TRUE) {
          for (p = 1; p <= RF_fmpIndexSize; p++) {
            if (RF_fmpIndex[p] < 0) {
              free_dvector(RF_fresponse[treeID][(uint) abs(RF_fmpIndex[p])], 1, RF_fobservationSize);
            }
            else {
              p = RF_fmpIndexSize;
            }
          }
          free_new_vvector(RF_fresponse[treeID], 1, RF_rSize, NRUTIL_DPTR);
        }
      }
    }
    if (RF_rFactorCount + RF_xFactorCount > 0) {
      if (RF_factorList[treeID] != NULL) {
        for (k = 1; k <= RF_maxFactorLevel; k++) {
          if (RF_factorList[treeID][k] != NULL) {
            free_Factor(RF_factorList[treeID][k]);
          }
        }
        free_new_vvector(RF_factorList[treeID], 1, RF_maxFactorLevel, NRUTIL_FPTR);
        RF_factorList[treeID] = NULL;
      }
    }
  }
  if (covrFlag) {
    vimpShadowFlag = FALSE;
    if ((RF_opt & OPT_VIMP) && (RF_opt & OPT_VIMP_TYP1) && !(RF_opt & OPT_VIMP_TYP2)) {
      vimpShadowFlag = TRUE;
    }
    if(vimpShadowFlag == TRUE) {
      free_dmatrix(RF_observation[treeID], 1, RF_xSize, 1, RF_observationSize);
    }
    else {
      if(RF_mPredictorFlag == TRUE) {
        for (p = 1; p <= RF_mpIndexSize; p++) {
          if (RF_mpIndex[p] > 0) {
            free_dvector(RF_observation[treeID][(uint) RF_mpIndex[p]], 1, RF_observationSize);
          }
        }
        free_new_vvector(RF_observation[treeID], 1, RF_xSize, NRUTIL_DPTR);
      }
    }
    if (mode == RF_PRED) {
      if(vimpShadowFlag == TRUE) {
        free_dmatrix(RF_fobservation[treeID], 1, RF_xSize, 1, RF_fobservationSize);
      }
      else {
        if(RF_fmPredictorFlag == TRUE) {
          for (p = 1; p <= RF_fmpIndexSize; p++) {
            if (RF_fmpIndex[p] > 0) {
              free_dvector(RF_fobservation[treeID][(uint) RF_fmpIndex[p]], 1, RF_fobservationSize);
            }
          }
          free_new_vvector(RF_fobservation[treeID], 1, RF_xSize, NRUTIL_DPTR);
        }
      }
    }
  }
}
char xferMissingness(char mode, Node *source, Terminal *destination) {
  uint   *sourceIndexPtr;
  double *sourceValuePtr;
  uint   *sourceAllocSizePtr;
  uint   *sourceActualSizePtr;
  char result;
  char xferFlag;
  sourceIndexPtr = NULL;  
  sourceValuePtr = NULL;  
  sourceAllocSizePtr = NULL;  
  sourceActualSizePtr = NULL;  
  result = FALSE;
  switch (mode) {
  case RF_PRED:
    if (RF_fmRecordSize > 0) {
      result = TRUE;
      sourceIndexPtr = source -> flmpIndex;
      sourceValuePtr = source -> flmpValue;
      sourceAllocSizePtr = & (source -> flmpIndexAllocSize);
      sourceActualSizePtr = & (source -> flmpIndexActualSize);
    }
    break;
  default:
    if (RF_mRecordSize > 0) {
      result = TRUE;
      sourceIndexPtr = source -> lmpIndex;
      sourceValuePtr = source -> lmpValue;
      sourceAllocSizePtr = & (source -> lmpIndexAllocSize);
      sourceActualSizePtr = & (source -> lmpIndexActualSize);
    }
    break;
  }
  if (result == FALSE) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Attempt to update forest impute data with no missingness in mode:  %10d", mode);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (*sourceActualSizePtr > 0) {
    (destination -> lmiIndex) = sourceIndexPtr;
    (destination -> lmiValue) = sourceValuePtr;
    (destination -> lmiSize)  = *sourceActualSizePtr;
    (destination -> lmiAllocSize)  = *sourceAllocSizePtr;    
    sourceIndexPtr = NULL;
    sourceValuePtr = NULL;
    *sourceAllocSizePtr  = 0;
    *sourceActualSizePtr  = 0;
    xferFlag = TRUE;
  }
  else {
    xferFlag = FALSE;
  }
  return xferFlag;
}
char getMarginalNodeMembership(char     mode,
                               char     rootFlag,
                               uint     treeID,
                               Node    *parent,
                               uint    *allMembrIndx,
                               uint     allMembrSize,
                               double **observationPtr) {
  char  bootResult;
  char leftResult, rghtResult;
  char terminalFlag;
  uint *leftAllMembrIndx;
  uint *rghtAllMembrIndx;
  uint leftAllMembrSize;
  uint rghtAllMembrSize;
  uint jLeft;
  uint jRght;
  char factorFlag;
  char daughterFlag;
  uint obsSize;
  uint i, j;
  factorFlag = FALSE; 
  terminalFlag = TRUE;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    break;
  default:
    obsSize = RF_observationSize;
    break;
  }
  if (RF_tLeafCount[treeID] > 0) {
    bootResult = TRUE;
    if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
      terminalFlag = FALSE;
      leftAllMembrIndx = rghtAllMembrIndx = NULL;
      leftAllMembrSize = rghtAllMembrSize = 0;
      uint *membershipIndicator = uivector(1, obsSize);
      leftAllMembrSize = rghtAllMembrSize = 0;
      daughterFlag = RIGHT;
      if (daughterFlag != NEITHER) {
        factorFlag = FALSE;
        if (strcmp(RF_xType[parent -> splitParameter], "C") == 0) {
          factorFlag = TRUE;
        }
        for (i = 1; i <= allMembrSize; i++) {
          daughterFlag = RIGHT;
          if (factorFlag == TRUE) {
            daughterFlag = splitOnFactor((uint) observationPtr[parent -> splitParameter][allMembrIndx[i]], parent -> splitValueFactPtr);
          }
          else {
            if ( observationPtr[parent -> splitParameter][allMembrIndx[i]] <= (parent -> splitValueCont) ) {
              daughterFlag = LEFT;
            }
          }
          membershipIndicator[allMembrIndx[i]] = daughterFlag;
          if (daughterFlag == LEFT) {
            leftAllMembrSize ++;
          }
          else {
            rghtAllMembrSize ++;
          }
        }  
      }  
      leftAllMembrIndx  = uivector(1, leftAllMembrSize + 1);
      rghtAllMembrIndx  = uivector(1, rghtAllMembrSize + 1);
      jLeft = jRght = 0;
      if (daughterFlag == NEITHER) {
        for (i = 1; i <= allMembrSize; i++) {
          leftAllMembrIndx[++jLeft] = allMembrIndx[i];
          rghtAllMembrIndx[++jRght] = allMembrIndx[i];
        }
      }
      else {
        for (i = 1; i <= allMembrSize; i++) {
          if (membershipIndicator[allMembrIndx[i]] == LEFT) {
            leftAllMembrIndx[++jLeft] = allMembrIndx[i];
          }
          else {
            rghtAllMembrIndx[++jRght] = allMembrIndx[i];
          }
        }
      }
      free_uivector(membershipIndicator, 1, obsSize);
      leftResult = getMarginalNodeMembership(mode,
                                             FALSE,
                                             treeID,
                                             parent -> left,
                                             leftAllMembrIndx,
                                             leftAllMembrSize,
                                             observationPtr);
      if(!leftResult) {
      }
      rghtResult = getMarginalNodeMembership(mode,
                                             FALSE,
                                             treeID,
                                             parent -> right,
                                             rghtAllMembrIndx,
                                             rghtAllMembrSize,
                                             observationPtr);
      if(!rghtResult) {
      }
      free_uivector(leftAllMembrIndx, 1, leftAllMembrSize + 1);
      free_uivector(rghtAllMembrIndx, 1, rghtAllMembrSize + 1);
    }  
    else {
    }
  }  
  else {
    if (rootFlag) {
      bootResult = FALSE;
      terminalFlag = FALSE;
    }
  }
  if (terminalFlag) {
    for (i = 1; i <= allMembrSize; i++) {
      RF_utTermMembership[treeID][allMembrIndx[i]][++ RF_utTermMembershipCount[treeID][allMembrIndx[i]]] = (parent -> nodeID);
      if ((RF_utTermMembershipCount[treeID][allMembrIndx[i]]) == (RF_utTermMembershipAlloc[treeID][allMembrIndx[i]] * MARGINAL_SIZE)) {
        RF_utTermMembershipAlloc[treeID][allMembrIndx[i]] ++;
        uint *utTermMembershipNew = uivector(1, RF_utTermMembershipAlloc[treeID][allMembrIndx[i]] * MARGINAL_SIZE);
        for (j = 1; j <= RF_utTermMembershipCount[treeID][allMembrIndx[i]]; j++) {
          utTermMembershipNew[j] = RF_utTermMembership[treeID][allMembrIndx[i]][j];
        }
        free_uivector(RF_utTermMembership[treeID][allMembrIndx[i]], 1, (RF_utTermMembershipAlloc[treeID][allMembrIndx[i]] - 1) * MARGINAL_SIZE);
        RF_utTermMembership[treeID][allMembrIndx[i]] = utTermMembershipNew;
      }
    }
  }  
  return bootResult;
}
char getPartialNodeMembership(char       rootFlag,
                              uint       treeID,
                              uint       partialIndex,
                              Node      *parent,
                              uint      *allMembrIndx,
                              uint       allMembrSize,
                              double   **observationPtr,
                              Terminal **membership) {
  char  bootResult;
  char leftResult, rghtResult;
  char terminalFlag;
  uint *leftAllMembrIndx;
  uint *rghtAllMembrIndx;
  uint leftAllMembrSize;
  uint rghtAllMembrSize;
  uint jLeft;
  uint jRght;
  char factorFlag;
  char daughterFlag;
  uint obsSize;
  uint primaryPartialIndex, secondaryPartialIndex;
  uint i, k;
  factorFlag = FALSE; 
  terminalFlag = TRUE;
  obsSize = RF_observationSize;
  if (RF_tLeafCount[treeID] > 0) {
    bootResult = TRUE;
    if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
      terminalFlag = FALSE;
      leftAllMembrIndx = rghtAllMembrIndx = NULL;
      leftAllMembrSize = rghtAllMembrSize = 0;
      uint *membershipIndicator = uivector(1, obsSize);
      leftAllMembrSize = rghtAllMembrSize = 0;
      daughterFlag = RIGHT;
      factorFlag = FALSE;
      if (strcmp(RF_xType[parent -> splitParameter], "C") == 0) {
        factorFlag = TRUE;
      }
      primaryPartialIndex = secondaryPartialIndex = k = 0;
      if (parent -> splitParameter ==  RF_partialXvar) {
        primaryPartialIndex = RF_partialXvar;
      }
      else {
        for (k = 1; k <= RF_partialLength2; k++) {
          if (parent -> splitParameter ==  RF_partialXvar2[k]) {
            secondaryPartialIndex = k;
          }
        }
      }
      for (i = 1; i <= allMembrSize; i++) {
        daughterFlag = RIGHT;
        if (factorFlag == TRUE) {
          if (primaryPartialIndex > 0) { 
            daughterFlag = splitOnFactor((uint) RF_partialValue[partialIndex], parent -> splitValueFactPtr);
          }
          else if (secondaryPartialIndex > 0) {
            daughterFlag = splitOnFactor((uint) RF_partialValue2[secondaryPartialIndex], parent -> splitValueFactPtr);
          }
          else {
            daughterFlag = splitOnFactor((uint) observationPtr[parent -> splitParameter][allMembrIndx[i]], parent -> splitValueFactPtr);
          }
        }
        else {
          if (primaryPartialIndex > 0) { 
            if (((parent -> splitValueCont) - RF_partialValue[partialIndex]) >= 0.0) {
              daughterFlag = LEFT;
            }
          }
          else if (secondaryPartialIndex > 0) {
            if (((parent -> splitValueCont) - RF_partialValue2[secondaryPartialIndex]) >= 0.0) {
              daughterFlag = LEFT;
            }
          }
          else {
            if (((parent -> splitValueCont) - observationPtr[parent -> splitParameter][allMembrIndx[i]]) >= 0.0) {
              daughterFlag = LEFT;
            }
          }
        }
        membershipIndicator[allMembrIndx[i]] = daughterFlag;
        if (daughterFlag == LEFT) {
          leftAllMembrSize ++;
        }
        else {
          rghtAllMembrSize ++;
        }
      }  
      leftAllMembrIndx  = uivector(1, leftAllMembrSize + 1);
      rghtAllMembrIndx  = uivector(1, rghtAllMembrSize + 1);
      jLeft = jRght = 0;
      for (i = 1; i <= allMembrSize; i++) {
        if (membershipIndicator[allMembrIndx[i]] == LEFT) {
          leftAllMembrIndx[++jLeft] = allMembrIndx[i];
        }
        else {
          rghtAllMembrIndx[++jRght] = allMembrIndx[i];
        }
      }
      free_uivector(membershipIndicator, 1, obsSize);
      leftResult = getPartialNodeMembership(FALSE,
                                            treeID,
                                            partialIndex,
                                            parent -> left,
                                            leftAllMembrIndx,
                                            leftAllMembrSize,
                                            observationPtr,
                                            membership);
      if(!leftResult) {
      }
      rghtResult = getPartialNodeMembership(FALSE,
                                            treeID,
                                            partialIndex,
                                            parent -> right,
                                            rghtAllMembrIndx,
                                            rghtAllMembrSize,
                                            observationPtr,
                                            membership);
      if(!rghtResult) {
      }
      free_uivector(leftAllMembrIndx, 1, leftAllMembrSize + 1);
      free_uivector(rghtAllMembrIndx, 1, rghtAllMembrSize + 1);
    }  
    else {
    }
  }  
  else {
    if (rootFlag) {
      bootResult = FALSE;
      terminalFlag = FALSE;
    }
  }
  if (terminalFlag) {
    for (i = 1; i <= allMembrSize; i++) {
      membership[allMembrIndx[i]] = RF_tTermList[treeID][(parent -> nodeID)];
    }
  }  
  return bootResult;
}
void acquireTree(char mode, uint r, uint b) {
  Node  *rootPtr;
  uint **mwcpPtrPtr;
  ulong  nodeOffset;
  char multImpFlag;
  uint *allMembrIndx;
  uint *fallMembrIndx;
  uint bootMembrIndxIter;
  uint rmbrIterator;
  uint ambrIterator;
  char  result;
  Node     ***gNodeMembership;
  Terminal ***gTermMembership;
  uint     obsSize;
  uint i;
  obsSize    = 0;  
  gNodeMembership = NULL;  
  gTermMembership = NULL;  
  multImpFlag = FALSE;
  if (mode == RF_GROW) {
    if (r > 1) {
      multImpFlag = TRUE;
    }
  }
#ifdef _OPENMP
#endif
  rootPtr = makeNode((mode == RF_GROW) ? RF_xSize : 0,
                     (mode == RF_GROW) ? 0 : 0,
                     (RF_opt & OPT_USPV_STAT) ? RF_randomResponseCount : 0,  
                     (mode == RF_GROW) ? ( (RF_opt & OPT_NODE_STAT) ? RF_randomCovariateCount : 0)  : 0);  
  RF_nodeMembership[b] = (Node **) new_vvector(1, RF_observationSize, NRUTIL_NPTR);
  RF_bootMembershipIndex[b] = uivector(1, RF_bootstrapSize);
  RF_bootMembershipFlag[b] = cvector(1, RF_observationSize);
  RF_bootMembershipCount[b] = uivector(1, RF_observationSize);
  RF_oobMembershipFlag[b] = cvector(1, RF_observationSize);
  RF_ibgMembershipIndex[b] = uivector(1, RF_observationSize);
  RF_oobMembershipIndex[b] = uivector(1, RF_observationSize);
  allMembrIndx = uivector(1, RF_observationSize);
  if (mode == RF_PRED) {
    RF_fnodeMembership[b] = (Node **) new_vvector(1, RF_fobservationSize, NRUTIL_NPTR);
  }
  RF_tTermMembership[b] = (Terminal **) new_vvector(1, RF_observationSize, NRUTIL_TPTR);
  if (mode == RF_PRED) {
    RF_ftTermMembership[b] = (Terminal **) new_vvector(1, RF_fobservationSize, NRUTIL_TPTR);
  }
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    gNodeMembership = RF_fnodeMembership;
    gTermMembership = RF_ftTermMembership;
    break;
  default:
    obsSize = RF_observationSize;
    gNodeMembership = RF_nodeMembership;
    gTermMembership = RF_tTermMembership;
    break;
  }
  switch (mode) {
  case RF_PRED:
    fallMembrIndx = uivector(1, RF_fobservationSize);
    break;
  default:
    fallMembrIndx = NULL;
    break;
  }
  if (RF_optHigh & OPT_MEMB_PRUN) {
    RF_pNodeMembership[b] = (Node **) new_vvector(1, obsSize, NRUTIL_NPTR);
  }
  stackShadow(mode, b);
  if (mode == RF_GROW) {
    if (RF_nImpute > 1) {
      if (r > 1) {
        if (RF_mRecordSize > 0) {
          imputeUpdateShadow(RF_GROW,
                             RF_response[b],
                             RF_observation[b]);
        }
        if (RF_timeIndex > 0) {
          if (RF_mTimeFlag == TRUE) {
            updateTimeIndexArray(0,
                                 NULL,
                                 RF_observationSize,
                                 RF_time[b],
                                 FALSE,
                                 FALSE,
                                 RF_masterTimeIndex[b]);
          }
        }
      }
    }
  }
  rootPtr -> parent = NULL;
  rootPtr -> nodeID = 1;
  RF_root[b] = rootPtr;
  RF_maxDepth[b] = 0;
  bootMembrIndxIter = 0;
  for (i = 1; i <= RF_observationSize; i++) {
    allMembrIndx[i] = i;
    RF_nodeMembership[b][i] = RF_root[b];
    RF_bootMembershipFlag[b][i]  = FALSE;
    RF_bootMembershipCount[b][i] = 0;
    RF_oobMembershipFlag[b][i]   = TRUE;
  }
  if (RF_optHigh & OPT_MEMB_PRUN) {
    for (i = 1; i <= obsSize; i++) {
      RF_pNodeMembership[b][i] = gNodeMembership[b][i];
    }
  }
  RF_orderedLeafCount[b] = 0;
  switch (mode) {
  case RF_GROW:
    RF_tLeafCount[b] = 0;
    stackNodeAndTermList(b, 0);
    rmbrIterator = ambrIterator = 0;
    result = growTree (r,
                       TRUE,
                       multImpFlag,
                       b,
                       rootPtr,
                       NULL,
                       0,
                       allMembrIndx,
                       RF_observationSize,
                       0,
                       RF_maxDepth + b,
                       & bootMembrIndxIter,
                       & rmbrIterator,
                       & ambrIterator);
    break;
  default:
    if (mode == RF_PRED) {
      for (i=1; i <= RF_fobservationSize; i++) {
        fallMembrIndx[i] = i;
        RF_fnodeMembership[b][i] = RF_root[b];
      }
    }
    nodeOffset = RF_restoreTreeOffset[b];
    mwcpPtrPtr = & RF_mwcpPtr[b];
    if (RF_tLeafCount[b] > 0) {
      stackNodeAndTermList(b, RF_tLeafCount[b]);
    }
    rmbrIterator = ambrIterator = 0;
    restoreTree(mode,
                b,
                rootPtr,
                & nodeOffset,
                RF_treeID_,
                RF_nodeID_,
                RF_parmID_,
                RF_contPT_,
                RF_mwcpSZ_,
                mwcpPtrPtr,
                0,
                RF_maxDepth + b);
    rmbrIterator = ambrIterator = 0;
    result = restoreNodeMembership(r,
                                   mode,
                                   TRUE,
                                   b,
                                   rootPtr,
                                   NULL,
                                   0,
                                   allMembrIndx,
                                   RF_observationSize,
                                   fallMembrIndx,
                                   RF_fobservationSize,
                                   & bootMembrIndxIter,
                                   & rmbrIterator,
                                   & ambrIterator);
    break;
  } 
  if (result) {
    if (RF_optHigh & OPT_MEMB_PRUN) {
      for (i = 1; i <= obsSize; i++) {
        RF_pNodeMembership[b][i] = RF_tNodeList[b][gTermMembership[b][i] -> nodeID];
      }
      RF_pLeafCount[b] = pruneTree(obsSize, b, RF_ptnCount);
      for (i=1; i <= obsSize; i++) {
        RF_PRUN_ID_ptr[b][i] = RF_pNodeMembership[b][i] -> nodeID;
      }
    }
    RF_oobSize[b] = 0;
    RF_ibgSize[b] = 0;
    for (i=1; i <= RF_observationSize; i++) {
      if (RF_bootMembershipFlag[b][i] == FALSE) {
        RF_oobSize[b] ++;
        RF_oobMembershipIndex[b][RF_oobSize[b]] = i;
      }
      else {
        RF_ibgSize[b] ++;
        RF_ibgMembershipIndex[b][RF_ibgSize[b]] = i;
      }
    }
    if (mode != RF_PRED) {
      if (RF_mRecordSize > 0) {
        for (i = 1; i <= RF_mRecordSize; i++) {
          if (RF_bootMembershipFlag[b][RF_mRecordIndex[i]] == TRUE) {
            RF_dmRecordBootFlag[b][i] = TRUE;
          }
          else {
            RF_dmRecordBootFlag[b][i] = FALSE;
          }
        }
      }  
    }  
    if (mode == RF_REST) {
      if(RF_sobservationSize > 0) {
        RF_soobSize[b] = 0;
        for (i = 1; i <= RF_sobservationSize; i++) {
          if (RF_bootMembershipFlag[b][RF_sobservationIndv[i]] == FALSE) {
            RF_soobSize[b] ++;
          }
        }
      }
    }
  }  
  if (result) {
    if (r == RF_nImpute) {
      if ((RF_opt & OPT_PERF) |
          (RF_opt & OPT_OENS) |
          (RF_opt & OPT_FENS)) {
        updateEnsembleCalculations(multImpFlag, mode, b);
      }
      if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
#ifdef _OPENMP
#pragma omp critical (_update_splitdepth)
#endif
        {  
          updateSplitDepth(b, RF_root[b], RF_maxDepth[b]);
        }
      }
      if (RF_opt & (OPT_VARUSED_F | OPT_VARUSED_T)) {
        getVariablesUsed(b, RF_root[b], RF_varUsedPtr[b]);
      }
      if (RF_opt & OPT_VIMP) {
        uint vimpCount;
        if (RF_opt & OPT_VIMP_JOIN) {
          vimpCount = 1;
        }
        else {
          vimpCount = RF_intrPredictorSize;
        }
        for (uint intrIndex = 1; intrIndex <= vimpCount; intrIndex++) {
          uint pp;
          if (!(RF_opt & OPT_VIMP_JOIN)) {
            pp = RF_intrPredictor[intrIndex];
          }
          else {
            pp = 0;
          }
          stackVimpMembership(mode, & RF_vimpMembership[intrIndex][b]);
          getVimpMembership(mode, b, RF_vimpMembership[intrIndex][b], pp);
#ifdef _OPENMP
#pragma omp critical (_update_vc)
#endif
          { 
            updateVimpCalculations(mode, b, intrIndex, RF_vimpMembership[intrIndex][b]);
          } 
          unstackVimpMembership(mode, RF_vimpMembership[intrIndex][b]);
        }
      }
      if (RF_optHigh & OPT_PART_PLOT) {
        Terminal **membership =  (Terminal **) new_vvector(1, RF_observationSize, NRUTIL_TPTR);
        for (i = 1; i <= RF_partialLength; i++) {
          getPartialNodeMembership(TRUE,
                                   b,
                                   i,
                                   rootPtr,
                                   RF_identityMembershipIndex,
                                   RF_observationSize,
                                   RF_observation[b],
                                   membership);
#ifdef _OPENMP
#pragma omp critical (_update_pc)
#endif
          { 
            updatePartialCalculations(b, i, membership);
          } 
        }
        free_new_vvector(membership, 1, RF_observationSize, NRUTIL_TPTR);
      }
      if (RF_opt & OPT_PROX) {
        updateProximity(mode, b);
      }
      if (mode == RF_GROW) {
        if (RF_opt & OPT_TREE) {
#ifdef _OPENMP
#pragma omp critical (_save_tree)
#endif
          { 
            mwcpPtrPtr = & RF_mwcpIterator;
            saveTree(b,
                     RF_root[b],
                     & RF_totalNodeCount1,
                     RF_treeID_,
                     RF_nodeID_,
                     RF_parmID_,
                     RF_contPT_,
                     RF_mwcpSZ_,
                     mwcpPtrPtr,
                     RF_mwcpCount);
          }
        }
      }
      if ((RF_opt & OPT_NODE_STAT) || (RF_opt & OPT_USPV_STAT)) {
#ifdef _OPENMP
#pragma omp critical (_save_stats)
#endif
        { 
          saveStatistics(mode,
                         b,
                         RF_root[b],
                         & RF_totalNodeCount2,
                         RF_spltST_,         
                         RF_spltVR_,         
                         RF_uspvST_ptr,      
                         RF_mtryID_ptr,      
                         RF_mtryST_ptr       
                         );
        }
      }
    }  
  }  
  else {
  }
  freeTree(b, RF_root[b], TRUE);
  unstackAuxiliary2(mode, b);
  unstackAuxiliary(mode, b);
  unstackNodeList(b);
  unstackShadow(mode, b, TRUE, TRUE);
  free_uivector(allMembrIndx, 1, RF_observationSize);
  if (mode == RF_PRED) {
    free_uivector(fallMembrIndx, 1, RF_fobservationSize);
  }
  if (!(RF_opt & OPT_MISS)) {
    for (i = 1; i <= RF_tLeafCount[b]; i++) {
      freeTerminal(RF_tTermList[b][i]);
    }
    unstackTermList(b);
    free_new_vvector(RF_tTermMembership[b], 1, RF_observationSize, NRUTIL_TPTR);
    if (mode == RF_PRED) {
      free_new_vvector(RF_ftTermMembership[b], 1, RF_fobservationSize, NRUTIL_TPTR);
    }
  }
}
void finalizeProximity(char mode) {
  uint  obsSize;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    break;
  default:
    obsSize = RF_observationSize;
    break;
  }
  if (RF_numThreads > 0) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
    for (uint i = 1; i <= obsSize; i++) {
      for (uint j = 1; j <= i; j++) {
        if (RF_proximityDenPtr[i][j] > 0) {
          RF_proximityPtr[i][j] = RF_proximityPtr[i][j] /  RF_proximityDenPtr[i][j];
        }
        else {
          RF_proximityPtr[i][j] = NA_REAL;
        }
      }
    }
  }
  else {
    for (uint i = 1; i <= obsSize; i++) {
      for (uint j = 1; j <= i; j++) {
        if (RF_proximityDenPtr[i][j] > 0) {
          RF_proximityPtr[i][j] = RF_proximityPtr[i][j] /  RF_proximityDenPtr[i][j];
        }
        else {
          RF_proximityPtr[i][j] = NA_REAL;
        }
      }
    }
  }
}
void updateProximity(char mode, uint b) {
  Terminal **tTermMembership;
  uint  *membershipIndex;
  uint   membershipSize;
  uint  mtnmFlag;
  if((RF_opt & OPT_PROX_IBG) && (RF_opt & OPT_PROX_OOB)) {
    switch (mode) {
    case RF_PRED:
      membershipSize = RF_fobservationSize;
      membershipIndex = RF_fidentityMembershipIndex;
      tTermMembership = RF_ftTermMembership[b];
      break;
    default:
      membershipSize = RF_observationSize;
      membershipIndex = RF_identityMembershipIndex;
      tTermMembership = RF_tTermMembership[b];
      break;
    }
  }
  else {
    if((RF_opt & OPT_PROX_IBG)  && !(RF_opt & OPT_PROX_OOB)) {
      membershipIndex = RF_ibgMembershipIndex[b];
      membershipSize  = RF_ibgSize[b];
    }
    else if(!(RF_opt & OPT_PROX_IBG)  && (RF_opt & OPT_PROX_OOB)) {
      membershipIndex = RF_oobMembershipIndex[b];
      membershipSize  = RF_oobSize[b];
    }
    else {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Illegal updateProximity() call.");
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
    tTermMembership = RF_tTermMembership[b];
  }
  mtnmFlag = FALSE;
  if (!mtnmFlag) {
    for (uint i = 1; i <= membershipSize; i++) {
      uint ii, jj;
      ii = membershipIndex[i];
      for (uint j = 1; j <= i; j++) {
        jj = membershipIndex[j];
#ifdef _OPENMP
#pragma omp atomic update
#endif
        RF_proximityDenPtr[ii][jj] ++;
        if ( tTermMembership[ii] == tTermMembership[jj] ) {
#ifdef _OPENMP
#pragma omp atomic update
#endif
          RF_proximityPtr[ii][jj] ++;
        }
      }
    }
  }
}
void updateSplitDepth(uint treeID, Node *rootPtr, uint maxDepth) {
  Node  *parent;
  double *localSplitDepth;
  uint index;
  uint i, j, k;
  if (RF_tLeafCount[treeID] > 0) {
    index = 0;  
    if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
      if (RF_opt & OPT_SPLDPTH_F) {
        index = 1;
      }
      else {
        index = treeID;
      }
    }
    else {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Illegal updateSplitDepth() call.");
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
    localSplitDepth = dvector(1, RF_xSize);
    for (i = 1; i <= RF_observationSize; i++) {
      for (j = 1; j <= RF_xSize; j++) {
        localSplitDepth[j] = NA_REAL;
      }
      parent = RF_nodeMembership[treeID][i];
      for (k = 1; k <= parent -> depth; k++) {
        if (ISNA(localSplitDepth[(parent -> splitDepth)[k]])) {
          localSplitDepth[(parent -> splitDepth)[k]] = (double) k;
        }
      }
      for (j = 1; j <= RF_xSize; j++) {
        if (ISNA(localSplitDepth[j])) {
          localSplitDepth[j] = (double) maxDepth + 1;
        }
      }
      for (j = 1; j <= RF_xSize; j++) {
        RF_splitDepthPtr[index][j][i] += localSplitDepth[j];
      }
    }
    free_dvector(localSplitDepth, 1, RF_xSize);
    freeSplitDepth(treeID);
  }
}
char pruneBranch(uint obsSize,
                 uint treeID,
                 Node **nodesAtDepth,
                 uint nadCount,
                 uint ptnTarget,
                 uint ptnCurrent) {
  char pruneFlag;
  uint i, j;
  pruneFlag = TRUE;
  double *varianceAtDepth =  dvector(1, nadCount);
  uint   *vadSortedIndex  = uivector(1, nadCount);
  for (i = 1; i <= nadCount; i++) {
    varianceAtDepth[i] = nodesAtDepth[i] -> variance;
  }
  indexx(nadCount, varianceAtDepth, vadSortedIndex);
  j = nadCount;
  while ((j >= 1) && pruneFlag) {
    nodesAtDepth[vadSortedIndex[j]] -> pseudoTerminal = TRUE;
    (nodesAtDepth[vadSortedIndex[j]] -> left)  -> pseudoTerminal = FALSE;
    (nodesAtDepth[vadSortedIndex[j]] -> right) -> pseudoTerminal = FALSE;
    for (i = 1; i <= obsSize; i++) {
      if ( (RF_pNodeMembership[treeID][i] == nodesAtDepth[vadSortedIndex[j]] -> left) ||
           (RF_pNodeMembership[treeID][i] == nodesAtDepth[vadSortedIndex[j]] -> right)) {
        RF_pNodeMembership[treeID][i] = nodesAtDepth[vadSortedIndex[j]];
      }
    }
    j --;
    ptnCurrent --;
    if (ptnCurrent <= ptnTarget) {
      pruneFlag = FALSE;
    }
  }
  free_dvector(varianceAtDepth, 1, nadCount);
  free_uivector(vadSortedIndex, 1, nadCount);
  return pruneFlag;
}
uint pruneTree(uint obsSize, uint treeID, uint ptnTarget) {
  Node **nodesAtDepth;
  uint   ptnCurrent;
  uint   nadCount;
  uint   tagDepth;
  char   pruneFlag;
  uint   i;
  if (ptnTarget < 1) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Illegal target PTN count in pruneTree():  %10d", ptnTarget);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (RF_tLeafCount[treeID] == 0) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Illegal call to pruneTree() on a rejected tree:  %10d", treeID);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  nodesAtDepth = (Node **) new_vvector(1, RF_tLeafCount[treeID], NRUTIL_NPTR);
  ptnCurrent = RF_tLeafCount[treeID];
  tagDepth = getMaximumDepth(RF_root[treeID]) - 1;
  pruneFlag = (ptnCurrent > ptnTarget) && (tagDepth > 0);
  while (pruneFlag) {
    for (i = 1; i <= RF_tLeafCount[treeID]; i++) {
      nodesAtDepth[i] = NULL;
    }
    nadCount = 0;
    getNodesAtDepth(RF_root[treeID], tagDepth, nodesAtDepth, &nadCount);
    pruneFlag = pruneBranch(obsSize, treeID, nodesAtDepth, nadCount, ptnTarget, ptnCurrent);
    if(pruneFlag) {
      ptnCurrent -= nadCount;
      tagDepth --;
    }
    else {
      ptnCurrent = ptnTarget;
    }
  }
  free_new_vvector(nodesAtDepth, 1, RF_tLeafCount[treeID], NRUTIL_NPTR);
  return ptnCurrent;
}
void unstackAuxiliary2(char mode, uint b) {
  free_uivector(RF_bootMembershipIndex[b], 1, RF_bootstrapSize);
  free_uivector(RF_bootMembershipCount[b], 1, RF_observationSize);
}
void unstackAuxiliary(char mode, uint b) {
  uint obsSize;
  obsSize = 0;  
  free_new_vvector(RF_nodeMembership[b], 1, RF_observationSize, NRUTIL_NPTR);
  free_cvector(RF_bootMembershipFlag[b], 1, RF_observationSize);
  free_cvector(RF_oobMembershipFlag[b], 1, RF_observationSize);
  free_uivector(RF_ibgMembershipIndex[b], 1, RF_observationSize);
  free_uivector(RF_oobMembershipIndex[b], 1, RF_observationSize);
  if (mode == RF_PRED) {
    free_new_vvector(RF_fnodeMembership[b],  1, RF_fobservationSize, NRUTIL_NPTR);
  }
  if (RF_optHigh & OPT_MEMB_PRUN) {
    switch (mode) {
    case RF_PRED:
      obsSize = RF_fobservationSize;
      break;
    default:
      obsSize = RF_observationSize;
      break;
    }
    free_new_vvector(RF_pNodeMembership[b], 1, obsSize, NRUTIL_NPTR);
  }
}
void stackNodeAndTermList(uint treeID, uint length) {
  uint size; 
  if (length == 0) {
    size = RF_theoreticalMaxtLeafCount[treeID];
  }
  else {
    size = length;
  }
  RF_tNodeListLength[treeID]  = size;
  RF_tNodeList[treeID] = (Node **) new_vvector(1, size , NRUTIL_NPTR);
  RF_tTermList[treeID] = (Terminal **) new_vvector(1, size, NRUTIL_TPTR);
}
void unstackNodeList(uint treeID) {
  if (RF_tNodeListLength[treeID] > 0) {
    free_new_vvector(RF_tNodeList[treeID], 1, RF_tNodeListLength[treeID], NRUTIL_NPTR);
  }
}
void unstackTermList(uint treeID) {
  if (RF_tNodeListLength[treeID] > 0) {
    free_new_vvector(RF_tTermList[treeID], 1, RF_tNodeListLength[treeID], NRUTIL_NPTR);
  }
}
void printPseudoTNInfo(char mode, uint b) {
  uint i;
  RF_pNodeList[b] = (Node **) new_vvector(1, RF_pLeafCount[b] + 1, NRUTIL_NPTR);
  i = 0;
  getPTNodeList(RF_root[b], RF_pNodeList[b], &i);
  free_new_vvector(RF_pNodeList[b], 1, RF_pLeafCount[b] + 1, NRUTIL_NPTR);
}      
Node *getTerminalNode(uint treeID, uint leaf) {
  uint i, j;
  Node *parent;
  parent = NULL;
  for (j = 1; j <= RF_observationSize; j++) {
    if ((RF_nodeMembership[treeID][j] -> nodeID) == leaf) {
      parent = RF_nodeMembership[treeID][j];
      j = RF_observationSize;
    }
  }
  if (parent == NULL) {
    Rprintf("\nDiagnostic Trace of (individual, boot, node, leaf) vectors in data set:  ");
    Rprintf("\n        index         boot         node         leaf \n");
    for (i = 1; i <= RF_observationSize; i++) {
      Rprintf(" %12d %12d %12x %12d \n", i,
              RF_bootMembershipFlag[treeID][i], RF_nodeMembership[treeID][i],
              RF_nodeMembership[treeID][i] -> nodeID);
    }
    Rprintf("\nDiagnostic State of TRAIN (SHADOW) data:  ");
    Rprintf("\n       index       status         time   observations -> \n");
    Rprintf("\n                                      ");
    for (i=1; i <= RF_xSize; i++) {
      Rprintf(" %12d", i);
    }
    Rprintf("\n");
    for (j = 1; j <= RF_observationSize; j++) {
      Rprintf("%12d %12.4f %12.4f", j, RF_status[treeID][j], RF_time[treeID][j]);
      for (i=1; i <= RF_xSize; i++) {
        Rprintf(" %12.4f", (RF_observation[treeID][i][j]));
      }
      Rprintf("\n");
    }
    Rprintf("\nDiagnostic State of TRAIN (INCOMING) data:  ");
    Rprintf("\n       index       status         time   observations -> \n");
    Rprintf("\n                                      ");
    for (i=1; i <= RF_xSize; i++) {
      Rprintf(" %12d", i);
    }
    Rprintf("\n");
    for (j = 1; j <= RF_observationSize; j++) {
      Rprintf("%12d %12.4f %12.4f", j, RF_responseIn[RF_statusIndex][j], RF_responseIn[RF_timeIndex][j]);
      for (i=1; i <= RF_xSize; i++) {
        Rprintf(" %12.4f", (RF_observationIn[i][j]));
      }
      Rprintf("\n");
    }
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Proxy member for (tree, node) = (%12d, %12d) not found.", treeID, leaf);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  return parent;
}
void getRawNodeSize(uint  type,
                    uint  treeID,
                    Node *parent,
                    uint *repMembrIndx,
                    uint *repMembrSize,
                    uint *allMembrIndx,
                    uint *allMembrSize) {
  uint      obsSize;
  Node   ***nodeMembershipPtr;
  uint      bootMembrSize;
  uint i;
  obsSize           = 0;     
  nodeMembershipPtr = NULL;  
  switch (type) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    nodeMembershipPtr = RF_fnodeMembership;
    break;
  default:
    obsSize = RF_observationSize;
    nodeMembershipPtr = RF_nodeMembership;
    break;
  }
  bootMembrSize = RF_bootstrapSize;
  *repMembrSize = 0;
  for (i=1; i <= bootMembrSize; i++) {
    if (RF_nodeMembership[treeID][RF_bootMembershipIndex[treeID][i]] == parent) {
      repMembrIndx[++(*repMembrSize)] = RF_bootMembershipIndex[treeID][i];
    }
  }
  *allMembrSize = 0;
  for (i=1; i <= obsSize; i++) {
    if (nodeMembershipPtr[treeID][i] == parent) {
      allMembrIndx[++(*allMembrSize)] = i;
    }
  }
}
char forkAndUpdate(uint    treeID,
                   Node   *parent,
                   uint   *repMembrIndx,
                   uint    repMembrSize,
                   uint   *allMembrIndx,
                   uint    allMembrSize,
                   uint    splitParameterMax,
                   double  splitValueMaxCont,
                   uint    splitValueMaxFactSize,
                   uint   *splitValueMaxFactPtr,
                   double  splitStatistic,
                   char   *localSplitIndicator,
                   char    splitMIA,
                   char    multImpFlag,
                   char   *membershipIndicator,
                   uint   *leftDaughterSize,
                   uint   *rghtDaughterSize) {
  char factorFlag;
  char daughterFlag;
  char result;
  char splitValueFlag;
  uint offset;
  uint i;
  factorFlag = FALSE; 
  result = forkNode(parent,
                    splitParameterMax,
                    splitValueMaxCont,
                    splitValueMaxFactSize,
                    splitValueMaxFactPtr);
  if (result == TRUE) {
    parent -> splitStatistic = splitStatistic;
    RF_tLeafCount[treeID]++;
    ((parent -> left) -> nodeID) = (parent -> nodeID);
    ((parent -> right) -> nodeID) = RF_tLeafCount[treeID];
    factorFlag = FALSE;
    if (strcmp(RF_xType[splitParameterMax], "C") == 0) {
      factorFlag = TRUE;
    }
    *leftDaughterSize = *rghtDaughterSize = 0;
    for (i = 1; i <= allMembrSize; i++) {
      membershipIndicator[allMembrIndx[i]] = NEITHER;
    }
    for (i = 1; i <= repMembrSize; i++) {
      membershipIndicator[repMembrIndx[i]] = localSplitIndicator[i];
    }
    offset = RF_rSize + splitParameterMax;
    for (i = 1; i <= allMembrSize; i++) {
      if(membershipIndicator[allMembrIndx[i]] == NEITHER) {
        if (splitMIA == SPLIT_MIA_NONE) {
          splitValueFlag = TRUE;
        }
        else {
          splitValueFlag = TRUE;
          if (RF_mRecordMap[allMembrIndx[i]] > 0) {
            if (RF_mpSign[offset][RF_mRecordMap[allMembrIndx[i]]] == 1) {
              splitValueFlag = FALSE;
              if (splitMIA == SPLIT_MIA_LEFT) {
                membershipIndicator[allMembrIndx[i]] = daughterFlag = LEFT;
              }
              else if (splitMIA == SPLIT_MIA_RGHT) {
                membershipIndicator[allMembrIndx[i]] = daughterFlag = RIGHT;
              }
              else if (splitMIA == SPLIT_MIA_UNIT) {
                membershipIndicator[allMembrIndx[i]] = daughterFlag = RIGHT;
              }
              else {
                RFprintf("\nRF-SRC:  *** ERROR *** ");
                RFprintf("\nRF-SRC:  splitMIA parameter invalid.");
                RFprintf("\nRF-SRC:  Please Contact Technical Support.");
                error("\nRF-SRC:  The application will now exit.\n");
              }
            }
          }
          if (splitValueFlag) {
            if (splitMIA == SPLIT_MIA_UNIT) {
              membershipIndicator[allMembrIndx[i]]  = daughterFlag = LEFT;
              splitValueFlag = FALSE;
            }
          }
        }
        if (splitValueFlag) {        
          daughterFlag = RIGHT;
          if (factorFlag == TRUE) {
            daughterFlag = splitOnFactor((uint) RF_observation[treeID][splitParameterMax][allMembrIndx[i]], splitValueMaxFactPtr);
          }
          else {
            if ((splitValueMaxCont - RF_observation[treeID][splitParameterMax][allMembrIndx[i]]) >= 0.0) {
              daughterFlag = LEFT;
            }
          }
          membershipIndicator[allMembrIndx[i]] = daughterFlag;
        }
      }  
      else {
        daughterFlag = membershipIndicator[allMembrIndx[i]];
      }
      if (daughterFlag == LEFT) {
        (*leftDaughterSize) ++;
        RF_nodeMembership[treeID][allMembrIndx[i]] = parent -> left;
      }
      else {
        (*rghtDaughterSize) ++;
        RF_nodeMembership[treeID][allMembrIndx[i]] = parent -> right;
      }
    }  
  }
  else {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  forkNode() failed.");
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (localSplitIndicator != NULL) {
    free_cvector(localSplitIndicator, 1, repMembrSize);
  }
  else {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  NULL Local Split Indicator encountered in forkAndUpdate().");
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  return result;
}
char growTree (uint     r,
               char     rootFlag,
               char     multImpFlag,
               uint     treeID,
               Node    *parent,
               uint    *repMembrIndx,
               uint     repMembrSize,
               uint    *allMembrIndx,
               uint     allMembrSize,
               uint     depth,
               uint    *maximumDepth,
               uint    *bootMembrIndxIter,
               uint    *rmbrIterator,
               uint    *ambrIterator) {
  char  bootResult;
  char  splitResult;
  char  forkResult;
  char leftResult, rghtResult;
  char terminalFlag;
  char bsUpdateFlag;
  uint *bootMembrIndx;
  uint *leftRepMembrIndx;
  uint *rghtRepMembrIndx;
  uint *leftAllMembrIndx;
  uint *rghtAllMembrIndx;
  uint bootMembrSize;
  uint leftAllMembrSize;
  uint rghtAllMembrSize;
  uint leftRepMembrSize, jLeft;
  uint rghtRepMembrSize, jRght;
  uint     splitParameterMax;
  double   splitValueMaxCont;
  uint     splitValueMaxFactSize;
  uint    *splitValueMaxFactPtr;
  double   splitStatistic;
  char    *splitIndicator;
  char     splitMIA;
  uint i, k, p;
  parent -> depth = depth;
  bootResult = TRUE;
  terminalFlag = TRUE;
  bsUpdateFlag = FALSE;
  splitIndicator = NULL;
  splitMIA = SPLIT_MIA_NONE;
  if (rootFlag || ((RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2))) {
    if ( (!(RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
         ( (RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ) {
      bootMembrIndx  = uivector(1, RF_bootstrapSize);
      bootMembrSize = RF_bootstrapSize;
    }
    else {
      bootMembrIndx  = uivector(1, allMembrSize);
      bootMembrSize = allMembrSize;
    }
    bootResult = bootstrap (RF_GROW,
                            treeID,
                            parent,
                            allMembrIndx,
                            allMembrSize,
                            bootMembrIndx,
                            bootMembrSize);
    if (rootFlag & bootResult) {
      if (!( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2) )) {
        bsUpdateFlag = TRUE;
      }
      if (RF_mRecordSize > 0) {
        if (TRUE) {
          for (p = 1; p <= RF_mpIndexSize; p++) {
            if (RF_mpIndex[p] > 0) {
              if (parent -> mpSign[p] == -1) {
                (parent -> permissibleSplit)[RF_mpIndex[p]] = FALSE;
              }
            }
          }
        }
        else {
          k = 0;
          for (p = 1; p <= RF_xSize; p++) {      
            if (parent -> mpSign[RF_rSize + p] != -1) {
              parent -> permissibleSplitIndex[++k] = p;
            }
          }
          parent -> permissibleSizeActual = k;
        }
      }
      else {
        if (FALSE) {
          for (p = 1; p <= RF_xSize; p++) {      
            parent -> permissibleSplitIndex[p] = p;
          }
          parent -> permissibleSizeActual = parent -> permissibleSizeAlloc;
        }
      }
    }
  }
  else {
    bootMembrIndx = repMembrIndx;
    bootMembrSize = repMembrSize;
    parent -> mpSign = (parent -> parent) -> mpSign;
  }
  if (bootResult) {
    if (multImpFlag == FALSE) {
      if (RF_mRecordSize > 0) {
        imputeNode(RF_GROW,
                   FALSE,  
                   TRUE,   
                   treeID,
                   parent,
                   bootMembrIndx,
                   bootMembrSize,
                   allMembrIndx,
                   allMembrSize);
        if (RF_timeIndex > 0) {
          if (RF_mTimeFlag == TRUE) {
            updateTimeIndexArray(treeID,
                                 allMembrIndx,
                                 allMembrSize,
                                 RF_time[treeID],
                                 FALSE,
                                 FALSE,
                                 RF_masterTimeIndex[treeID]);
          }
        }
      }
    }  
  }  
  if (bootResult) {
    if (rootFlag) {
      RF_tLeafCount[treeID] = 1;
    }
    splitResult = getBestSplit(treeID,
                               parent,
                               bootMembrIndx,
                               bootMembrSize,
                               allMembrIndx,
                               allMembrSize,
                               & splitParameterMax,
                               & splitValueMaxCont,
                               & splitValueMaxFactSize,
                               & splitValueMaxFactPtr,
                               & splitStatistic,
                               & splitIndicator,
                               & splitMIA,
                               multImpFlag);
    if (splitResult == TRUE) {
      if (FALSE) {
        parent -> permissibleSizeActual = parent -> permissibleSizeAlloc;
      }
      else {
      }
      terminalFlag = FALSE;
      char *membershipIndicator = cvector(1, RF_observationSize);
      forkResult = forkAndUpdate(treeID,
                                 parent,
                                 bootMembrIndx,
                                 bootMembrSize,
                                 allMembrIndx,
                                 allMembrSize,
                                 splitParameterMax,
                                 splitValueMaxCont,
                                 splitValueMaxFactSize,
                                 splitValueMaxFactPtr,
                                 splitStatistic,
                                 splitIndicator,
                                 splitMIA,
                                 multImpFlag,
                                 membershipIndicator,
                                 &leftAllMembrSize,
                                 &rghtAllMembrSize);
      if (forkResult == TRUE) {
        leftAllMembrIndx  = uivector(1, leftAllMembrSize);
        rghtAllMembrIndx  = uivector(1, rghtAllMembrSize);
        jLeft = jRght = 0;
        for (i = 1; i <= allMembrSize; i++) {
          if (membershipIndicator[allMembrIndx[i]] == LEFT) {
            leftAllMembrIndx[++jLeft] = allMembrIndx[i];
          }
          else {
            rghtAllMembrIndx[++jRght] = allMembrIndx[i];
          }
        }
        if ( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2) ) {
          leftRepMembrIndx = leftAllMembrIndx;
          leftRepMembrSize = leftAllMembrSize;
          rghtRepMembrIndx = rghtAllMembrIndx;
          rghtRepMembrSize = rghtAllMembrSize;
        }
        else {
          leftRepMembrIndx  = uivector(1, bootMembrSize);
          rghtRepMembrIndx  = uivector(1, bootMembrSize);
          leftRepMembrSize = rghtRepMembrSize = 0;
          for (i = 1; i <= bootMembrSize; i++) {
            if (membershipIndicator[bootMembrIndx[i]] == LEFT) {
              leftRepMembrIndx[++leftRepMembrSize] = bootMembrIndx[i];
            }
            else {
              rghtRepMembrIndx[++rghtRepMembrSize] = bootMembrIndx[i];
            }
          }
        }
        free_cvector(membershipIndicator, 1, RF_observationSize);
        leftResult = growTree (r,
                               FALSE,
                               multImpFlag,
                               treeID,
                               parent -> left,
                               leftRepMembrIndx,
                               leftRepMembrSize,
                               leftAllMembrIndx,
                               leftAllMembrSize,
                               (parent -> depth) + 1,
                               maximumDepth,
                               bootMembrIndxIter,
                               rmbrIterator,
                               ambrIterator);
        if(!leftResult) {
        }
        rghtResult = growTree (r,
                               FALSE,
                               multImpFlag,
                               treeID,
                               parent -> right,
                               rghtRepMembrIndx,
                               rghtRepMembrSize,
                               rghtAllMembrIndx,
                               rghtAllMembrSize,
                               (parent -> depth) + 1,
                               maximumDepth,
                               bootMembrIndxIter,
                               rmbrIterator,
                               ambrIterator);
        if(!rghtResult) {
        }
        free_uivector(leftAllMembrIndx, 1, leftAllMembrSize);
        free_uivector(rghtAllMembrIndx, 1, rghtAllMembrSize);
        if ( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2) ) {
        }
        else {
          free_uivector(leftRepMembrIndx, 1, bootMembrSize);
          free_uivector(rghtRepMembrIndx, 1, bootMembrSize);
        }
      }
      else {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  forkAndUpdate(%10d) failed.", treeID);
        RFprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }  
    else {
      parent -> splitFlag = FALSE;
      free_cvector(parent -> permissibleSplit, 1, parent -> xSize);
      parent -> permissibleSplit = NULL;
      parent -> xSize = 0;
    }
  }  
  else {
    if (rootFlag) {
      if (!bootResult) {
        terminalFlag = FALSE;
      }
    }
  }  
  if (terminalFlag) {
    parent -> pseudoTerminal = TRUE;
    RF_tNodeList[treeID][parent -> nodeID] = parent;
    RF_tTermList[treeID][parent -> nodeID] = makeTerminal();
    RF_tTermList[treeID][parent -> nodeID] -> nodeID = parent -> nodeID;
    if (RF_opt & OPT_MISS) {
      imputeNodeAndSummarize(r,
                             RF_GROW,
                             treeID,
                             parent,
                             bootMembrIndx,
                             bootMembrSize,
                             allMembrIndx,
                             allMembrSize,
                             NULL,
                             0);
    }
    if (r == RF_nImpute) {
      if (RF_optHigh & OPT_MEMB_USER) {
        for (i = 1; i <= allMembrSize; i++) {
          RF_MEMB_ID_ptr[treeID][allMembrIndx[i]] = parent -> nodeID;
        }
      }
      updateTerminalNodeOutcomesNew(RF_GROW,
                                    treeID,
                                    RF_tTermList[treeID][parent -> nodeID],
                                    bootMembrIndx,
                                    bootMembrSize,
                                    allMembrIndx,
                                    allMembrSize,
                                    rmbrIterator,
                                    ambrIterator);
      if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
        getSplitDepth(parent, maximumDepth);
      }
    }
    else {
      initTerminalNodeMembership(treeID,
                                 RF_tTermList[treeID][parent -> nodeID],
                                 allMembrIndx,
                                 allMembrSize);
    }
    if ( (RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2) ) {
      bsUpdateFlag = TRUE;
    }
  }
  if (bsUpdateFlag) {
    for (i = 1; i <= bootMembrSize; i++) {
      RF_bootMembershipIndex[treeID][++(*bootMembrIndxIter)] = bootMembrIndx[i];
      RF_bootMembershipFlag[treeID][bootMembrIndx[i]] = TRUE;
      RF_oobMembershipFlag[treeID][bootMembrIndx[i]]  = FALSE;
      RF_bootMembershipCount[treeID][bootMembrIndx[i]] ++;
      if (RF_optHigh & OPT_MEMB_USER) {
        RF_BOOT_CT_ptr[treeID][bootMembrIndx[i]] ++;
      }
    }
  }
  if (rootFlag || ((RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2))) {
    if ( (!(RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
         ( (RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ) {
      free_uivector(bootMembrIndx, 1, RF_bootstrapSize);
    }
    else {
      free_uivector(bootMembrIndx, 1, allMembrSize);
    }
  }
  return bootResult;
}
char restoreTree(char    mode,
                 uint    b,
                 Node   *parent,
                 ulong  *offset,
                 uint   *treeID,
                 uint   *nodeID,
                 uint   *parmID,
                 double *contPT,
                 uint   *mwcpSZ,
                 uint  **mwcpPtr,
                 uint    depth,
                 uint   *maximumDepth) {
  char terminalFlag;
  uint i;
  if (b != treeID[*offset]) {
    RFprintf("\nRF-SRC:  Diagnostic Trace of Tree Record:  \n");
    RFprintf("\nRF-SRC:      treeID     nodeID     parmID       spltPT     mwcpSZ ");
    RFprintf("\nRF-SRC:  %10d %10d %10d %12.4f %10d \n", treeID[*offset], nodeID[*offset], parmID[*offset], contPT[*offset], mwcpSZ[*offset]);
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Invalid forest input record in tree:  %10d", b);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  parent -> depth = depth;
  parent -> left  = NULL;
  parent -> right = NULL;
  parent -> splitFlag = FALSE;
  parent -> nodeID = nodeID[*offset];
  parent -> splitParameter = parmID[*offset];
  if ((parent -> splitParameter) != 0) {
    if (strcmp(RF_xType[parent -> splitParameter], "C") == 0) {
      parent -> splitValueFactSize = mwcpSZ[*offset];
      parent -> splitValueFactPtr = uivector(1, mwcpSZ[*offset]);
      for (i = 1; i <= parent -> splitValueFactSize; i++) {
        (*mwcpPtr) ++;
        (parent -> splitValueFactPtr)[i] = **mwcpPtr;
      }
      parent -> splitValueCont = NA_REAL;
    }
    else {
      parent -> splitValueCont = contPT[*offset];
      parent -> splitValueFactSize = 0;
      parent -> splitValueFactPtr = NULL;
    }
  }
  else {
    parent -> splitValueCont     = NA_REAL;
    parent -> splitValueFactSize = 0;
    parent -> splitValueFactPtr  = NULL;
  }
  (*offset) ++;
  if ((parent -> splitParameter) != 0) {
    terminalFlag = FALSE;
    parent -> left  = makeNode(0, 0, parent -> urStatSize, parent -> mtrySize);
    setParent(parent ->  left, parent);
    restoreTree(mode,
                b,
                parent -> left,
                offset,
                treeID,
                nodeID,
                parmID,
                contPT,
                mwcpSZ,
                mwcpPtr,
                parent -> depth + 1,
                maximumDepth);
    parent -> right = makeNode(0, 0, parent -> urStatSize, parent -> mtrySize);
    setParent(parent -> right, parent);
    restoreTree(mode,
                b,
                parent -> right,
                offset,
                treeID,
                nodeID,
                parmID,
                contPT,
                mwcpSZ,
                mwcpPtr,
                parent -> depth + 1,
                maximumDepth);
  }
  else {
    terminalFlag = TRUE;
  }
  if (terminalFlag) {
    if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
      getSplitDepth(parent, maximumDepth);
    }
    parent -> pseudoTerminal = TRUE;
    RF_tNodeList[b][parent -> nodeID] = parent;
    RF_tTermList[b][parent -> nodeID] = makeTerminal();
    RF_tTermList[b][parent -> nodeID] -> nodeID = parent -> nodeID;
  }
  return (!terminalFlag);
}
void saveTree(uint    b,
              Node   *parent,
              ulong  *offset,
              uint   *treeID,
              uint   *nodeID,
              uint   *parmID,
              double *contPT,
              uint   *mwcpSZ,
              uint  **mwcpPtr,
              uint   *mwcpCT) {
  uint i;
  treeID[*offset] = b;
  nodeID[*offset] = parent -> nodeID;
  parmID[*offset] = parent -> splitParameter;
  if ((parent -> splitParameter) != 0) {
    if (strcmp(RF_xType[parent -> splitParameter], "C") == 0) {
      mwcpSZ[*offset] = parent -> splitValueFactSize;
      for (i = 1; i <= mwcpSZ[*offset]; i++) {
        (*mwcpPtr) ++;
        **mwcpPtr = (parent -> splitValueFactPtr)[i];
        mwcpCT[b] ++;
      }
      contPT[*offset] = NA_REAL;
    }
    else {
      contPT[*offset] = parent -> splitValueCont;
      mwcpSZ[*offset] = 0;
    }
  }
  else {
    contPT[*offset] = NA_REAL;
    mwcpSZ[*offset] = 0;
  }
  (*offset) ++;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    saveTree(b, parent ->  left, offset, treeID, nodeID, parmID, contPT, mwcpSZ, mwcpPtr, mwcpCT);
    saveTree(b, parent -> right, offset, treeID, nodeID, parmID, contPT, mwcpSZ, mwcpPtr, mwcpCT);
  }
}
void freeTree(uint treeID, Node *parent, char rootFlag) {
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    freeTree(treeID, parent -> left, FALSE);
    freeTree(treeID, parent -> right, FALSE);
  }
  freeNode(parent);
}
void getSplitDepth(Node *parent, uint *maximumDepth) {
  Node *reversePtr;
  uint i;
  if (!(RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T))) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Call to calculate split depth without the option being active.");
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (parent -> depth > 0) {
    *maximumDepth = ((parent -> depth > *maximumDepth) ? parent -> depth : *maximumDepth);
    stackSplitDepth(parent, parent -> depth);
    reversePtr = parent;
    for (i = 1; i <= parent -> depth; i++) {
      if ((reversePtr -> parent) == NULL) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Reverse parsing of tree failed in restoreTree().");
        RFprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
      (parent -> splitDepth)[(parent -> depth) - i + 1] = (reversePtr -> parent) -> splitParameter;
      reversePtr = reversePtr -> parent;
    }
  }
}
void freeSplitDepth(uint treeID) {
  uint j;
  for (j = 1; j <= RF_tLeafCount[treeID]; j++) {
    unstackSplitDepth(RF_tNodeList[treeID][j]);
  }
}
void saveStatistics(char     mode,
                    uint     b,
                    Node    *parent,
                    ulong   *offset,
                    double  *spltST,
                    double  *spltVR,
                    uint   **uspvST,
                    uint   **mtryID,
                    double **mtryST) {
  uint i;
  if (!(RF_opt & OPT_NODE_STAT) && !(RF_opt & OPT_USPV_STAT)) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Inconsistent call to saveStatistics().  The options are NOT active.");
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    RFprintf("\nRF-SRC:  The application will now exit.\n");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (RF_opt & OPT_NODE_STAT) {  
    switch (mode) {
    case RF_GROW:
      spltST[*offset] = parent -> splitStatistic;
      for(i = 1; i <= RF_randomCovariateCount; i++) {
        mtryID[*offset][i] = parent -> mtryIndx[i];
        mtryST[*offset][i] = parent -> mtryStat[i];
      }
      break;
    default:
      if (RF_ptnCount == 0) {
        spltST[*offset] = parent -> variance;
      }
      else {
        spltST[*offset] = parent -> pseudoTerminal;
      }
      break;
    }
  }
  if (RF_opt & OPT_USPV_STAT) {
    if (mode == RF_GROW) {
      for (i = 1; i <= RF_randomResponseCount; i++) {
        uspvST[*offset][i] = (parent -> urStat)[i];
      }
    }
  }
  (*offset) ++;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    saveStatistics(mode, b, parent ->  left, offset, spltST, spltVR, uspvST, mtryID, mtryST);
    saveStatistics(mode, b, parent -> right, offset, spltST, spltVR, uspvST, mtryID, mtryST);
  }
}
uint getMaximumDepth(Node *parent) {
  uint result, rLeft, rRight;
  result = parent -> depth;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    rLeft = getMaximumDepth(parent ->  left);
    rRight = getMaximumDepth(parent -> right);
    result = (rLeft > rRight) ? rLeft : rRight;
  }
  return result;
}
void getNodesAtDepth(Node *parent, uint tagDepth, Node **nodesAtDepth, uint *nadCount) {
  char recurseFlag;
  recurseFlag = TRUE;
  if (tagDepth == parent -> depth) {
    if ((parent -> splitParameter) != 0) {
      (*nadCount) ++;
      nodesAtDepth[*nadCount] = parent;
    }
    recurseFlag = FALSE;
  }
  else {
    if (((parent -> left) == NULL) && ((parent -> right) == NULL)) {
      recurseFlag = FALSE;
    }
  }
  if (recurseFlag) {
    getNodesAtDepth(parent ->  left, tagDepth, nodesAtDepth, nadCount);
    getNodesAtDepth(parent -> right, tagDepth, nodesAtDepth, nadCount);
  }
}
void getTreeInfo(uint treeID, Node *parent) {
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    getTreeInfo(treeID, parent ->  left);
    getTreeInfo(treeID, parent -> right);
  }
}
void getPTNodeList(Node    *parent,
                   Node   **list,
                   uint    *offset) {
  if (!(parent -> pseudoTerminal)) {
    getPTNodeList(parent ->  left, list, offset);
    getPTNodeList(parent -> right, list, offset);
  }
  else {
    (*offset) ++;
    list[*offset] = parent;
  }
}
void initTerminalNodeMembership(uint       treeID,
                                Terminal  *parent,
                                uint      *allMembrIndx,
                                uint       allMembrSize) {
  uint i;
  for (i = 1; i <= allMembrSize; i++) {
    RF_tTermMembership[treeID][allMembrIndx[i]] = parent;
  }
}
void updateTerminalNodeOutcomesNew (char       mode,
                                    uint       treeID,
                                    Terminal  *parent,
                                    uint      *repMembrIndx,
                                    uint       repMembrSize,
                                    uint      *allMembrIndx,
                                    uint       allMembrSize,
                                    uint      *rmbrIterator,
                                    uint      *ambrIterator) {
  uint clasIterator, regrIterator;
  uint i;
  if (RF_optHigh & OPT_MEMB_INCG) {
    for (i = 1; i <= RF_TN_ACNT_ptr[treeID][parent -> nodeID]; i++) {
      ++(*ambrIterator);
      RF_tTermMembership[treeID][RF_AMBR_ID_ptr[treeID][(*ambrIterator)]] = parent;
    }
  }
  else if (RF_optHigh & OPT_MEMB_OUTG) {
    for (i = 1; i <= allMembrSize; i++) {
      RF_tTermMembership[treeID][allMembrIndx[i]] = parent;
      RF_AMBR_ID_ptr[treeID][++(*ambrIterator)] = allMembrIndx[i];
    }
    RF_TN_ACNT_ptr[treeID][parent -> nodeID] = allMembrSize;
  }
  else {
    for (i = 1; i <= allMembrSize; i++) {
      RF_tTermMembership[treeID][allMembrIndx[i]] = parent;
    }
  }
  if ((RF_opt & OPT_PERF) |
      (RF_opt & OPT_OENS) |
      (RF_opt & OPT_FENS)) {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      getAtRiskAndEventCounts(treeID, parent, repMembrIndx, repMembrSize, allMembrIndx, allMembrSize, rmbrIterator);
      if (!(RF_optHigh & OPT_TERM_INCG)) {
        getLocalRatios(treeID, parent);
        getLocalSurvival(treeID, parent);
        if (!(RF_opt & OPT_COMP_RISK)) {
          getLocalNelsonAalen(treeID, parent);
        }
        else {
          getLocalCSH(treeID, parent);
          getLocalCIF(treeID, parent);
        }
        unstackAtRiskAndEventCounts(parent);
      }
      if (!(RF_opt & OPT_COMP_RISK)) {
        getSurvival(treeID, parent);
        getNelsonAalen(treeID, parent);
      }
      else {
        getCSH(treeID, parent);
        getCIF(treeID, parent);
      }
      getMortality(treeID, parent);
      freeTerminalNodeLocalSurvivalStructures(parent);
    }
    else {
      clasIterator = regrIterator = *rmbrIterator;
      if (RF_rFactorCount > 0) {
        getMultiClassProbNew(treeID, parent, repMembrIndx, repMembrSize, allMembrIndx, allMembrSize, & clasIterator);
        *rmbrIterator = clasIterator;
      }
      if (RF_rNonFactorCount > 0) {
        getMeanResponseNew(treeID, parent, repMembrIndx, repMembrSize, allMembrIndx, allMembrSize, & regrIterator);
        *rmbrIterator = regrIterator;
      }
    }
  }
  else {
    getMembrCountOnly(treeID, parent, repMembrIndx, repMembrSize, allMembrIndx, allMembrSize, rmbrIterator);
  }
}
void updateEnsembleCalculations (char      multImpFlag,
                                 char      mode,
                                 uint      b) {
  uint      obsSize;
  double  **responsePtr;
  char      potentiallyMixedMultivariate;
  char      respImputeFlag;
  uint      thisSerialTreeCount;
  uint      j;
  respImputeFlag       = FALSE;  
  responsePtr          = NULL;   
  obsSize              = 0;      
  thisSerialTreeCount  = 0;      
  obsSize = (mode == RF_PRED) ?  RF_fobservationSize : RF_observationSize;
#ifdef _OPENMP
#pragma omp critical (_update_ensemble)
#endif
  { 
    if (RF_tLeafCount[b] > 0) {
      RF_serialTreeIndex[++RF_serialTreeCount] = b;
      thisSerialTreeCount = RF_serialTreeCount;
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        updateEnsembleSurvival(mode, b, thisSerialTreeCount);
      }  
      else {
        potentiallyMixedMultivariate = FALSE;
        if (RF_rTargetFactorCount > 0) {
          updateEnsembleMultiClass(mode, b, thisSerialTreeCount, potentiallyMixedMultivariate);
          potentiallyMixedMultivariate = TRUE;
        }
        if (RF_rTargetNonFactorCount > 0) {
          updateEnsembleMean(mode, b, thisSerialTreeCount, potentiallyMixedMultivariate);
          potentiallyMixedMultivariate = TRUE;
        }
      }
      if (getPerformanceFlag(mode, thisSerialTreeCount)) {
        respImputeFlag = stackAndImputePerfResponse(mode,
                                                    multImpFlag,
                                                    b,
                                                    thisSerialTreeCount,
                                                    &responsePtr);
      }
      else {
        respImputeFlag = FALSE;
      }
      if (getPerformanceFlag(mode, thisSerialTreeCount)) {
        if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
          getPerformance(thisSerialTreeCount,
                         mode,
                         obsSize,
                         responsePtr,
                         (mode == RF_PRED) ? RF_fullEnsembleDen : RF_oobEnsembleDen,
                         (mode == RF_PRED) ? RF_fullEnsembleMRTptr : RF_oobEnsembleMRTptr,
                         NULL,
                         NULL,
                         RF_perfMRTptr[thisSerialTreeCount],
                         NULL,
                         NULL);
        }
        else {
          if (RF_rTargetFactorCount > 0) {
            getPerformance(thisSerialTreeCount,
                           mode,
                           obsSize,
                           responsePtr,
                           (mode == RF_PRED) ? RF_fullEnsembleDen : RF_oobEnsembleDen,
                           NULL,
                           (mode == RF_PRED) ? RF_fullEnsembleCLSptr : RF_oobEnsembleCLSptr,
                           NULL,
                           NULL,
                           RF_perfCLSptr[thisSerialTreeCount],
                           NULL);
          }
          if (RF_rTargetNonFactorCount > 0) {
            getPerformance(thisSerialTreeCount,
                           mode,
                           obsSize,
                           responsePtr,
                           (mode == RF_PRED) ? RF_fullEnsembleDen : RF_oobEnsembleDen,
                           NULL,
                           NULL,
                           (mode == RF_PRED) ? RF_fullEnsembleRGRptr : RF_oobEnsembleRGRptr,
                           NULL,
                           NULL,
                           RF_perfRGRptr[thisSerialTreeCount]);
          }
        }
        unstackImputeResponse(respImputeFlag, obsSize, responsePtr);
      }  
    }  
    else {
      RF_serialTreeIndex[++RF_serialTreeCount] = b;
    }  
    if (getUserTraceFlag()) {
      double userTimeElapsedFromStart;
      double userTimeElapsedFromSplit;
      double userTimeRemaining;
      time_t current;
      current = time(NULL);
      userTimeElapsedFromSplit = (double) (current - RF_userTimeSplit);
      if ((userTimeElapsedFromSplit) > (double) getUserTraceFlag()) {
        userTimeElapsedFromStart = (double) (current - RF_userTimeStart);
        userTimeRemaining = (userTimeElapsedFromStart / RF_serialTreeCount * RF_forestSize) - userTimeElapsedFromStart;
        RFprintf("\nTrees Grown:  %6d,    Time Remaining (sec):  %6.0f", RF_serialTreeCount, ceil(userTimeRemaining));
        RF_userTimeSplit = current;
      }
    }
  }  
  if (RF_tLeafCount[b] > 0) {
    switch (mode) {
    case RF_GROW:
      if (!(RF_optHigh & OPT_TERM_OUTG)) {
        if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
          for (j = 1; j <= RF_tLeafCount[b]; j++) {
            freeTerminalNodeSurvivalStructuresIntermediate(RF_tTermList[b][j]);
          }
          if (!(RF_opt & OPT_VIMP)) {
            for (j = 1; j <= RF_tLeafCount[b]; j++) {
              freeTerminalNodeSurvivalStructuresFinal(RF_tTermList[b][j]);
            }
          }
        }
        else {
          if (!(RF_opt & OPT_VIMP)) {
            for (j = 1; j <= RF_tLeafCount[b]; j++) {
              freeTerminalNodeNonSurvivalStructures(RF_tTermList[b][j]);
            }
          }
        }
      }
      break;
    default:
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        if (!(RF_optHigh & OPT_PART_PLOT)) {
          for (j = 1; j <= RF_tLeafCount[b]; j++) {
            freeTerminalNodeSurvivalStructuresIntermediate(RF_tTermList[b][j]);
          }
        }
        if (!(RF_opt & OPT_VIMP) && !(RF_optHigh & OPT_PART_PLOT)) {
          for (j = 1; j <= RF_tLeafCount[b]; j++) {
            freeTerminalNodeSurvivalStructuresFinal(RF_tTermList[b][j]);
          }
        }
      }
      else {
        if (!(RF_opt & OPT_VIMP) && !(RF_optHigh & OPT_PART_PLOT)) {
          for (j = 1; j <= RF_tLeafCount[b]; j++) {
            freeTerminalNodeNonSurvivalStructures(RF_tTermList[b][j]);
          }
        }
      }
      break;
    }
    if (RF_opt & OPT_VIMP) {
      if (RF_opt & OPT_VIMP_LEOB) {
        summarizeTreePerformance(mode, b);
      }
    }
  }  
}
char stackAndImputePerfResponse(char      mode,
                                char      multImpFlag,
                                uint      treeID,
                                uint      serialID,
                                double ***responsePtr) {
  uint     obsSize;
  char     imputeFlag;
  imputeFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    if (treeID == 0) {
      *responsePtr = RF_fresponseIn;
    }
    else {
      *responsePtr = RF_fresponse[treeID];
    }
    if (RF_fmRecordSize > 0) {
      if(RF_fmResponseFlag == TRUE) {
        imputeFlag = TRUE;
      }
    }
    break;
  default:
    obsSize  = RF_observationSize;
    if (treeID == 0) {
      *responsePtr = RF_responseIn;
    }
    else {
      *responsePtr = RF_response[treeID];
    }
    if (multImpFlag == FALSE) {
      if (RF_mRecordSize > 0) {
        if(RF_mResponseFlag == TRUE) {
          imputeFlag = TRUE;
        }
      }
    }
    break;
  }
  *responsePtr = stackAndImputeGenericResponse(imputeFlag, mode, obsSize, treeID, serialID, *responsePtr);
  return imputeFlag;
}
double **stackAndImputeGenericResponse(char flag,
                                       char mode,
                                       uint obsSize,
                                       uint treeID,
                                       uint serialID,
                                       double **responsePtr) {
  uint i, p;
  double **mResponsePtr;
  if (flag == TRUE) {
    mResponsePtr   = dmatrix(1, RF_rSize, 1, obsSize);
    for (i = 1; i <= obsSize; i++) {
      for (p = 1; p <= RF_rSize; p++) {
        mResponsePtr[p][i] = responsePtr[p][i];
      }
    }
    imputeResponse(mode, serialID, mResponsePtr);
  }
  else {
    mResponsePtr = responsePtr;
  }
  return mResponsePtr;
}
void unstackImputeResponse(char flag, uint obsSize, double **mResponsePtr) {
  if (flag == TRUE) {
    free_dmatrix(mResponsePtr, 1, RF_rSize, 1, obsSize);
  }
}
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
                    double   *perfRGRptr) {
  uint      j, k;
  double   *cpv;
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (!(RF_opt & OPT_COMP_RISK)) {
      perfMRTptr[1] = getConcordanceIndex(1,
                                          obsSize,
                                          responsePtr[RF_timeIndex],
                                          responsePtr[RF_statusIndex],
                                          outcomeMRT[1],
                                          denomPtr);
    }
    else {
      cpv = dvector(1, RF_eventTypeSize);
      getCRPerformance(mode,
                       obsSize,
                       responsePtr,
                       outcomeMRT,
                       denomPtr,
                       cpv);
      for (j=1; j <= RF_eventTypeSize; j++) {
        perfMRTptr[j] = cpv[j];
      }
      free_dvector(cpv, 1, RF_eventTypeSize);
    }
  }
  else {
    if (perfCLSptr != NULL) {
      for (j = 1; j <= RF_rTargetFactorCount; j++) {
        cpv = dvector(1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]);
        if (RF_opt & OPT_PERF_CALB) {
          perfCLSptr[j][1] = getBrierScore(obsSize,
                                           RF_rTargetFactor[j],                                                            
                                           responsePtr[RF_rTargetFactor[j]],
                                           outcomeCLS[j],
                                           denomPtr,
                                           cpv);
          for (k=1; k <=RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
            perfCLSptr[j][1+k] = cpv[k];
          }
        }
        else {
          perfCLSptr[j][1] = getClassificationIndex(obsSize,
                                                    responsePtr[RF_rTargetFactor[j]],
                                                    outcomeCLS[j][1],
                                                    denomPtr);
          getConditionalClassificationIndex(obsSize,
                                            RF_rTargetFactor[j],
                                            responsePtr[RF_rTargetFactor[j]],
                                            outcomeCLS[j][1],
                                            denomPtr,
                                            cpv);
          for (k=1; k <=RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
            perfCLSptr[j][1+k] = cpv[k];
          }
        }
        free_dvector(cpv, 1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]);
      }
    }
    if (perfRGRptr != NULL) {
      for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
        perfRGRptr[j] = getMeanSquareError(obsSize,
                                           responsePtr[RF_rTargetNonFactor[j]],
                                           outcomeRGR[j],
                                           denomPtr);
      }
    }
  }
}
void finalizeEnsembleEstimates(char mode) {
  char oobFlag, fullFlag;
  uint      obsSize;
  double ***ensembleSRGptr;
  double  **ensembleMRTptr;
  double  **ensembleSRVptr;
  double ***ensembleCIFptr;
  double ***ensembleCLSptr;
  double  **ensembleRGRptr;
  double ***ensembleSRGnum;
  double  **ensembleMRTnum;
  double  **ensembleSRVnum;
  double ***ensembleCIFnum;
  double ***ensembleCLSnum;
  double  **ensembleRGRnum;
  uint     *ensembleDen;
  uint i, j, k;
  oobFlag = fullFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    oobFlag = FALSE;
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    break;
  default:
    obsSize = RF_observationSize;
    if (RF_opt & OPT_OENS) {
      oobFlag = TRUE;
    }
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    break;
  }
  while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
    if (oobFlag == TRUE) {
      ensembleDen    = RF_oobEnsembleDen;
      ensembleSRGptr = RF_oobEnsembleSRGptr;
      ensembleMRTptr = RF_oobEnsembleMRTptr;        
      ensembleSRVptr = RF_oobEnsembleSRVptr;
      ensembleCIFptr = RF_oobEnsembleCIFptr;
      ensembleCLSptr = RF_oobEnsembleCLSptr;
      ensembleRGRptr = RF_oobEnsembleRGRptr;
      ensembleSRGnum = RF_oobEnsembleSRGnum;
      ensembleMRTnum = RF_oobEnsembleMRTnum;        
      ensembleSRVnum = RF_oobEnsembleSRVnum;
      ensembleCIFnum = RF_oobEnsembleCIFnum;
      ensembleCLSnum = RF_oobEnsembleCLSnum;
      ensembleRGRnum = RF_oobEnsembleRGRnum;
    }
    else {
      ensembleDen    = RF_fullEnsembleDen;
      ensembleSRGptr = RF_fullEnsembleSRGptr;
      ensembleMRTptr = RF_fullEnsembleMRTptr;        
      ensembleSRVptr = RF_fullEnsembleSRVptr;
      ensembleCIFptr = RF_fullEnsembleCIFptr;
      ensembleCLSptr = RF_fullEnsembleCLSptr;
      ensembleRGRptr = RF_fullEnsembleRGRptr;
      ensembleSRGnum = RF_fullEnsembleSRGnum;
      ensembleMRTnum = RF_fullEnsembleMRTnum;        
      ensembleSRVnum = RF_fullEnsembleSRVnum;
      ensembleCIFnum = RF_fullEnsembleCIFnum;
      ensembleCLSnum = RF_fullEnsembleCLSnum;
      ensembleRGRnum = RF_fullEnsembleRGRnum;
    }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      for (i = 1; i <= obsSize; i++) {
        if (ensembleDen[i] != 0) {
          if (!(RF_opt & OPT_COMP_RISK)) {
            ensembleMRTptr[1][i] = ensembleMRTnum[1][i] / ensembleDen[i];
            for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
              ensembleSRGptr[1][k][i] = ensembleSRGnum[1][k][i] / ensembleDen[i];
              ensembleSRVptr[k][i]    = ensembleSRVnum[k][i] / ensembleDen[i];
            }
          }
          else {
            for(j = 1; j <= RF_eventTypeSize; j ++) {
              ensembleMRTptr[j][i] = ensembleMRTnum[j][i] / ensembleDen[i];
              for (k=1; k <= RF_sortedTimeInterestSize; k++) {
                ensembleSRGptr[j][k][i] = ensembleSRGnum[j][k][i] / ensembleDen[i];
                ensembleCIFptr[j][k][i] = ensembleCIFnum[j][k][i] / ensembleDen[i];
              }
            }
          }
        }
        else {
        }
      }
    }  
    else {
      if (RF_rTargetFactorCount > 0) {
        for (i = 1; i <= obsSize; i++) {
          if (ensembleDen[i] != 0) {
            for (j = 1; j <= RF_rTargetFactorCount; j++) {
              for (k=1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
                ensembleCLSptr[j][k][i] = ensembleCLSnum[j][k][i] / ensembleDen[i];
              }
            }
          }
          else {
          }
        }
      }
      if (RF_rTargetNonFactorCount > 0) {      
        for (i = 1; i <= obsSize; i++) {
          if (ensembleDen[i] != 0) {
            for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
              ensembleRGRptr[j][i] = ensembleRGRnum[j][i] / ensembleDen[i];
            }
          }
          else {
          }
        }
      }
    }
    if (oobFlag == TRUE) {
      oobFlag = FALSE;
    }
    else {
      fullFlag = FALSE;
    }
  }  
}
char getPerformanceFlag (char mode, uint serialTreeID) {
  char result;
  if (RF_opt & OPT_PERF) {
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  if (result) {
    if (RF_optHigh & OPT_TREE_ERR) {
      if (serialTreeID < RF_forestSize) {
        result = FALSE;
      }
    }
  }
  return result;
}
void getVariablesUsed(uint treeID, Node *parent, uint *varUsedVector) {
  if (RF_tLeafCount[treeID] > 0) {
    if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
      varUsedVector[parent -> splitParameter] ++;
      getVariablesUsed(treeID, parent ->  left, varUsedVector);
      getVariablesUsed(treeID, parent -> right, varUsedVector);
    }
  }
  return;
}
void setUserTraceFlag (uint traceFlag) { RF_userTraceFlag = traceFlag; }
uint getUserTraceFlag () { return RF_userTraceFlag; }
Node *identifyPerturbedMembership (Node    *parent,
                                   double **shadowVIMP,
                                   uint     index) {
  char daughterFlag;
  Node *result = parent;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    daughterFlag = RIGHT;
    if (strcmp(RF_xType[parent -> splitParameter], "C") == 0) {
      daughterFlag = splitOnFactor((uint) shadowVIMP[parent -> splitParameter][index], parent -> splitValueFactPtr);
    }
    else {
      if (((parent -> splitValueCont) - shadowVIMP[parent -> splitParameter][index]) >= 0.0) {
        daughterFlag = LEFT;
      }
    }
    if (daughterFlag == LEFT) {
      result = identifyPerturbedMembership(parent ->  left, shadowVIMP, index);
    }
    else {
      result = identifyPerturbedMembership(parent -> right, shadowVIMP, index);
    }
  }
  return result;
}
Node *randomizeMembership(Node    *parent,
                          double **predictor,
                          uint     individual,
                          uint     splitParameter,
                          uint     treeID) {
  char daughterFlag;
  char randomSplitFlag;
  Node *result;
  result = parent;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    randomSplitFlag = FALSE;
    if (splitParameter > 0) {
      if ((parent -> splitParameter) == splitParameter) {
        randomSplitFlag = TRUE;
      }
    }
    else {
      if(RF_importanceFlag[parent -> splitParameter] == TRUE) {
        randomSplitFlag = TRUE;
      }
    }
    if(randomSplitFlag == TRUE) {
      if (ran1C(treeID) <= 0.5) {
        result = randomizeMembership(parent ->  left, predictor, individual, splitParameter, treeID);
      }
      else {
        result = randomizeMembership(parent -> right, predictor, individual, splitParameter, treeID);
      }
    }
    else {
      daughterFlag = RIGHT;
      if (strcmp(RF_xType[parent -> splitParameter], "C") == 0) {
        daughterFlag = splitOnFactor((uint) predictor[parent -> splitParameter][individual], parent -> splitValueFactPtr);
      }
      else {
        if (((parent -> splitValueCont) - predictor[parent -> splitParameter][individual]) >= 0.0) {
          daughterFlag = LEFT;
        }
      }
      if (daughterFlag == LEFT) {
        result = randomizeMembership(parent ->  left, predictor, individual, splitParameter, treeID);
      }
      else {
        result = randomizeMembership(parent -> right, predictor, individual, splitParameter, treeID);
      }
    }
  }
  return result;
}
Node *antiMembership(Node    *parent,
                     double **predictor,
                     uint     individual,
                     uint     splitParameter,
                     uint     treeID) {
  char daughterFlag;
  char antiSplitFlag;
  Node *result;
  result = parent;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    antiSplitFlag = FALSE;
    if (splitParameter > 0) {
      if ((parent -> splitParameter) == splitParameter) {
        antiSplitFlag = TRUE;
      }
    }
    else {
      if(RF_importanceFlag[parent -> splitParameter] == TRUE) {
        antiSplitFlag = TRUE;
      }
    }
    daughterFlag = RIGHT;
    if (strcmp(RF_xType[parent -> splitParameter], "C") == 0) {
      daughterFlag = splitOnFactor((uint) predictor[parent -> splitParameter][individual], parent -> splitValueFactPtr);
    }
    else {
      if (((parent -> splitValueCont) - predictor[parent -> splitParameter][individual]) >= 0.0) {
          daughterFlag = LEFT;
      }
    }
    if(antiSplitFlag == TRUE) {
      if (daughterFlag == LEFT) {
        daughterFlag = RIGHT;
      }
      else {
        daughterFlag = LEFT;
      }
    }
    if (daughterFlag == LEFT) {
      result = randomizeMembership(parent ->  left, predictor, individual, splitParameter, treeID);
    }
    else {
      result = randomizeMembership(parent -> right, predictor, individual, splitParameter, treeID);
    }
  }
  return result;
}
void permute(uint ranGenID, uint parallelID, uint n, uint *indx) {
  float (*ranX) (uint);
  uint i,j,k;
  ranX = NULL;  
  if ((ranGenID != 1) && (ranGenID != 2) && (ranGenID != 3)) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Invalid random generator selected:  %10d", ranGenID);
    RFprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  switch(ranGenID) {
  case 1:
    ranX = ran1A;
    break;
  case 2:
    ranX = ran1B;
    break;
  case 3:
    ranX = ran1C;
    break;
  }
  for (i=1; i<= n; i++) {
    indx[i] = 0;
  }
  for (i=n; i > 0; i--) {
    k = (uint) ceil(ranX(parallelID)*(i*1.0));
    for (j = 1; k > 0; j++) {
      if (indx[j] == 0) {
        k--;
      }
    }
    indx[j-1] = i;
  }
}
void getAntiMembership (char       mode,
                        uint       treeID,
                        Terminal **vimpMembership,
                        uint       p) {
  Node    *rootPtr;
  uint    *membershipIndex;
  uint     membershipSize;
  double **predictorPtr;
  char    *membershipFlag;
  uint     i;
  uint     ii;
  membershipFlag = NULL;  
  rootPtr = RF_root[treeID];
  switch (mode) {
  case RF_PRED:
    membershipSize = RF_fobservationSize;
    membershipIndex = RF_fidentityMembershipIndex;
    predictorPtr = RF_fobservation[treeID];
    break;
  default:
    membershipSize  = RF_oobSize[treeID];
    membershipIndex = RF_oobMembershipIndex[treeID];
    if(RF_sobservationSize > 0) {
      membershipFlag = RF_bootMembershipFlag[treeID];
    }
    else {
    }
    predictorPtr = RF_observation[treeID];
    break;
  }
  if (RF_sobservationSize > 0) {
    for (i = 1; i <= membershipSize; i++) {
      ii = membershipIndex[i];
      vimpMembership[ii] = RF_tTermMembership[treeID][ii];
    }
    for (i = 1; i <= RF_sobservationSize; i++) {
      ii = RF_sobservationIndv[i];
      if (membershipFlag[ii] == FALSE) {
        vimpMembership[ii] = RF_tTermList[treeID][ antiMembership(rootPtr, predictorPtr, ii, p, treeID) -> nodeID ];
      }
    }
  }
  else {
    for (i = 1; i <= membershipSize; i++) {
      ii = membershipIndex[i];
      vimpMembership[ii] = RF_tTermList[treeID][ antiMembership(rootPtr, predictorPtr, ii, p, treeID) -> nodeID ];
    }
  }
}
void getRandomMembership (char       mode,
                          uint       treeID,
                          Terminal **vimpMembership,
                          uint       p) {
  Node    *rootPtr;
  uint    *membershipIndex;
  uint     membershipSize;
  double **predictorPtr;
  char    *membershipFlag;
  uint     i;
  uint     ii;
  membershipFlag = NULL;  
  rootPtr = RF_root[treeID];
  switch (mode) {
  case RF_PRED:
    membershipSize = RF_fobservationSize;
    membershipIndex = RF_fidentityMembershipIndex;
    predictorPtr = RF_fobservation[treeID];
    break;
  default:
    membershipSize  = RF_oobSize[treeID];
    membershipIndex = RF_oobMembershipIndex[treeID];
    if(RF_sobservationSize > 0) {
      membershipFlag = RF_bootMembershipFlag[treeID];
    }
    else {
    }
    predictorPtr = RF_observation[treeID];
    break;
  }
  if (RF_sobservationSize > 0) {
    for (i = 1; i <= membershipSize; i++) {
      ii = membershipIndex[i];
      vimpMembership[ii] = RF_tTermMembership[treeID][ii];
    }
    for (i = 1; i <= RF_sobservationSize; i++) {
      ii = RF_sobservationIndv[i];
      if (membershipFlag[ii] == FALSE) {
        vimpMembership[ii] = RF_tTermList[treeID][ randomizeMembership(rootPtr, predictorPtr, ii, p, treeID) -> nodeID ];
      }
    }
  }
  else {
    for (i = 1; i <= membershipSize; i++) {
      ii = membershipIndex[i];
      vimpMembership[ii] = RF_tTermList[treeID][ randomizeMembership(rootPtr, predictorPtr, ii, p, treeID) -> nodeID ];
    }
  }  
}
void getPermuteMembership (char       mode,
                           uint       treeID,
                           Terminal **vimpMembership,
                           uint       p) {
  Node    *rootPtr;
  uint     obsSize;
  uint    *membershipIndex;
  uint     membershipSize;
  double **predictorPtr;
  char    *membershipFlag;
  uint     permuteObsSize;
  uint    *indexVIMP;
  uint    *permuteVIMP;
  double **shadowVIMP;
  uint     pInnerCount, pIn;
  uint     i, j, k, targetCov;
  uint     ii;
  membershipFlag = NULL;  
  rootPtr = RF_root[treeID];
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    membershipSize = RF_fobservationSize;
    membershipIndex = RF_fidentityMembershipIndex;
    permuteObsSize = RF_fobservationSize;
    predictorPtr = RF_fobservation[treeID];
    break;
  default:
    obsSize = RF_observationSize;
    membershipSize  = RF_oobSize[treeID];
    membershipIndex = RF_oobMembershipIndex[treeID];
    if(RF_sobservationSize > 0) {
      permuteObsSize = RF_soobSize[treeID];
      membershipFlag = RF_bootMembershipFlag[treeID];
    }
    else {
      permuteObsSize = RF_oobSize[treeID];
    }
    predictorPtr = RF_observation[treeID];
    break;
  }
  indexVIMP = uivector(1, permuteObsSize + 1);
  permuteVIMP = uivector(1, permuteObsSize + 1);
  if (RF_sobservationSize > 0) {
    k = 0;
    for (i = 1; i <= RF_sobservationSize; i++) {
      if (membershipFlag[RF_sobservationIndv[i]] == FALSE) {
        indexVIMP[++k] = RF_sobservationIndv[i];
      }
    }
  }
  else {
    for (i = 1; i <= membershipSize; i++) {
      indexVIMP[i] = membershipIndex[i];
    }
  }
  if (p > 0) {
    pInnerCount = 1;
  }
  else {
    pInnerCount = RF_intrPredictorSize;
  }
  shadowVIMP = (double **) new_vvector(1, RF_xSize, NRUTIL_DPTR);
  for (j = 1; j <= RF_xSize; j++) {
    shadowVIMP[j] = predictorPtr[j];
  }
  for (pIn = 1; pIn <= pInnerCount; pIn++) {
    if (p > 0) {
      targetCov = p;
    }
    else {
      targetCov = RF_intrPredictor[pIn];
    }
    shadowVIMP[targetCov] = dvector(1, obsSize);
    for (i = 1; i <= obsSize; i++) {
      shadowVIMP[targetCov][i] = predictorPtr[targetCov][i];
    }
    permute(3, treeID, permuteObsSize, permuteVIMP);
    for (k = 1; k <= permuteObsSize; k++) {
      shadowVIMP[targetCov][indexVIMP[k]] = predictorPtr[targetCov][indexVIMP[permuteVIMP[k]]];
    }
  }
  if (RF_sobservationSize > 0) {
    for (i = 1; i <= membershipSize; i++) {
      ii = membershipIndex[i];
      vimpMembership[ii] = RF_tTermMembership[treeID][ii];
    }
    for (i = 1; i <= RF_sobservationSize; i++) {
      ii = RF_sobservationIndv[i];
      if (membershipFlag[ii] == FALSE) {
        vimpMembership[ii] = RF_tTermList[treeID][ identifyPerturbedMembership(rootPtr, shadowVIMP, ii) -> nodeID ];
      }
    }
  }
  else {
    for (i = 1; i <= membershipSize; i++) {
      ii = membershipIndex[i];
      vimpMembership[ii] = RF_tTermList[treeID][ identifyPerturbedMembership(rootPtr, shadowVIMP, ii) -> nodeID ];
    }
  }  
  for (pIn = 1; pIn <= pInnerCount; pIn++) {
    if (p > 0) {
      targetCov = p;
    }
    else {
      targetCov = RF_intrPredictor[pIn];
    }
    free_dvector(shadowVIMP[targetCov], 1, obsSize);
  }
  free_new_vvector(shadowVIMP, 1, RF_xSize, NRUTIL_DPTR);
  free_uivector(indexVIMP, 1, permuteObsSize + 1);
  free_uivector(permuteVIMP, 1, permuteObsSize + 1);
}
void getVimpMembership (char mode, uint treeID, Terminal **vimpMembership, uint p) {
  char result;
  if (!(RF_opt & OPT_VIMP)) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to compute variable importance though not requested.");
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  result = FALSE;
  switch (mode) {
  case RF_PRED:
    result = TRUE;
    break;
  default:
    if (RF_oobSize[treeID] > 0) {
      result = TRUE;
    }
    break;
  }
  if (result == TRUE) {
    if (!(RF_opt & OPT_VIMP_TYP1) && !(RF_opt & OPT_VIMP_TYP2)) {
      getAntiMembership(mode, treeID, vimpMembership, p);
    }
    else if ((RF_opt & OPT_VIMP_TYP1) && !(RF_opt & OPT_VIMP_TYP2)) { 
      getPermuteMembership(mode, treeID, vimpMembership, p);
    }
    else if (!(RF_opt & OPT_VIMP_TYP1) && (RF_opt & OPT_VIMP_TYP2)) { 
      getRandomMembership(mode, treeID, vimpMembership, p);
    }
    else {
      RFprintf("\nRF-SRC:  *** ERROR *** ");
      RFprintf("\nRF-SRC:  Unknown VIMP type encountered:  %10d", RF_opt);
      RFprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
  }
}
void updateVimpEnsemble (char       mode,
                         uint       treeID,
                         Terminal **vimpMembership,
                         uint       p) {
  char   ensembleFlag;
  uint   obsSize;
  uint i;
  if (RF_opt & OPT_VIMP_LEOB) {
    ensembleFlag = FALSE;
    switch (mode) {
    case RF_PRED:
      obsSize = RF_fobservationSize;
      break;
    default:
      obsSize  = RF_observationSize;
      break;
    }
    for (i = 1; i <= obsSize; i++) {
      RF_vimpEnsembleDen[p][i] = 0;
    }
  }
  else {
    ensembleFlag = TRUE;
  }
  updateGenericVimpEnsemble(mode,
                            treeID,
                            p,              
                            NULL,           
                            vimpMembership, 
                            ensembleFlag,   
                            RF_vimpEnsembleMRT,
                            RF_vimpEnsembleCLS,
                            RF_vimpEnsembleRGR);
}
void updateTreeEnsemble (char       mode,
                         uint       treeID,
                         uint      *treeDenom,
                         double  ***treeEnsembleMRT,
                         double ****treeEnsembleCLS,
                         double  ***treeEnsembleRGR) {
  Terminal **termMembership;
  switch (mode) {
  case RF_PRED:
    termMembership = RF_ftTermMembership[treeID];
    break;
  default:
    termMembership = RF_tTermMembership[treeID];
    break;
  }
  updateGenericVimpEnsemble(mode,
                            treeID,
                            1,              
                            treeDenom,      
                            termMembership, 
                            FALSE,          
                            treeEnsembleMRT,
                            treeEnsembleCLS,
                            treeEnsembleRGR);
}
void updateGenericVimpEnsemble (char       mode,
                                uint       treeID,
                                uint       xVarIdx,
                                uint      *treeDenom,
                                Terminal **noiseMembership,
                                char       ensembleFlag,
                                double  ***genEnsembleMRT,
                                double ****genEnsembleCLS,
                                double  ***genEnsembleRGR) {
  Terminal *terminalNode;
  uint  *membershipIndex;
  uint   membershipSize;
  uint   *denomPtr;
  uint   i, j, k;
  uint   ii;
  switch (mode) {
  case RF_PRED:
    membershipSize = RF_fobservationSize;
    membershipIndex = RF_fidentityMembershipIndex;
    break;
  default:
    membershipSize  = RF_oobSize[treeID];
    membershipIndex = RF_oobMembershipIndex[treeID];
    break;
  }
  if (treeDenom != NULL) {
    denomPtr = treeDenom;
  }
  else {
    denomPtr = RF_vimpEnsembleDen[xVarIdx];
  }
  for (i = 1; i <= membershipSize; i++) {
      ii = membershipIndex[i];
      terminalNode = noiseMembership[ii];
      if ((terminalNode -> membrCount) > 0) {
        if (!ensembleFlag) {
          denomPtr[ii] = 1;
          if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
            for (k=1; k <= RF_eventTypeSize; k++) {
              genEnsembleMRT[xVarIdx][k][ii] = terminalNode -> mortality[k];
            }
          }
          else {
            if (RF_rTargetFactorCount > 0) {
              if (!(RF_opt & OPT_PERF_CALB)) {
                for (j=1; j <= RF_rTargetFactorCount; j++) {
                  genEnsembleCLS[xVarIdx][j][1][ii] = (double) (terminalNode -> maxClass)[RF_rFactorMap[RF_rTargetFactor[j]]];
                }
              }
              else {
                for (j=1; j <= RF_rTargetFactorCount; j++) {
                  for (k=1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
                    genEnsembleCLS[xVarIdx][j][k][ii] = (double) (terminalNode -> multiClassProb)[RF_rFactorMap[RF_rTargetFactor[j]]][k] / (double) (terminalNode -> membrCount);
                  }
                }
              }
            }
            if (RF_rTargetNonFactorCount > 0) {
              for (j=1; j <= RF_rTargetNonFactorCount; j++) {
                genEnsembleRGR[xVarIdx][j][ii] = (terminalNode -> meanResponse)[RF_rNonFactorMap[RF_rTargetNonFactor[j]]];
              }
            }
          }
        }
        else {
          denomPtr[ii] ++;
          if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
            for (k=1; k <= RF_eventTypeSize; k++) {
              genEnsembleMRT[xVarIdx][k][ii] += terminalNode -> mortality[k];
            }
          }
          else {
            if (RF_rTargetFactorCount > 0) {
              for (j=1; j <= RF_rTargetFactorCount; j++) {
                for (k=1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
                  genEnsembleCLS[xVarIdx][j][k][ii] += (double) (terminalNode -> multiClassProb)[RF_rFactorMap[RF_rTargetFactor[j]]][k] / (double) (terminalNode -> membrCount);
                }
              }
            }
            if (RF_rTargetNonFactorCount > 0) {
              for (j=1; j <= RF_rTargetNonFactorCount; j++) {
                genEnsembleRGR[xVarIdx][j][ii] += (terminalNode -> meanResponse)[RF_rNonFactorMap[RF_rTargetNonFactor[j]]];
              }
            }
          }
        }
      }
      else {
        if (!(RF_opt & OPT_OUTC_TYPE)) {
          RFprintf("\nRF-SRC:  *** ERROR *** ");
          RFprintf("\nRF-SRC:  NA encountered for VIMP outcome in terminal node:  %10d", terminalNode -> nodeID);
          RFprintf("\nRF-SRC:  Please Contact Technical Support.");
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
  }  
}
void summarizeVimpPerformance(char       mode,
                              uint       treeID,
                              uint       p) {
  uint      obsSize;
  double   **responsePtr;
  uint      *vimpDenom;
  double   **ensembleMRT;
  double  ***ensembleCLS;
  double   **ensembleRGR;
  double    *vimpMRTptr;
  double   **vimpCLSptr;
  double    *vimpRGRptr;
  char        imputeFlag;
  double    maxValue, maxClass;
  uint i, j, k;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    break;
  default:
    obsSize  = RF_observationSize;
    break;
  }
  if (RF_opt & OPT_VIMP_LEOB) {
    switch (mode) {
    case RF_PRED:
      responsePtr = RF_fresponse[treeID];
      break;
    default:
      responsePtr = RF_response[treeID];
      break;
    }
    imputeFlag = FALSE;
  }
  else {
    imputeFlag = stackAndImputePerfResponse(mode, FALSE, 0, RF_forestSize, &responsePtr);
  } 
  vimpDenom = RF_vimpEnsembleDen[p];
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (!(RF_opt & OPT_VIMP_LEOB)) {
      if (!(RF_opt & OPT_COMP_RISK)) {
        getEnsembleMortality(mode, treeID, obsSize, RF_vimpEnsembleMRT[p], vimpDenom, RF_vimpEnsembleMRT[p][1]);
      }
      else {
        getEnsembleMortalityCR(mode, treeID, obsSize, RF_vimpEnsembleMRT[p], vimpDenom, RF_vimpEnsembleMRT[p]);
      }
    }
  }  
  else {
    if (RF_rTargetFactorCount > 0) {
      if (!(RF_opt & OPT_VIMP_LEOB)) {
        if (!(RF_opt & OPT_PERF_CALB)) {
          for (i = 1; i <= obsSize; i++) {
            if(vimpDenom[i] > 0) {
              for (j = 1; j <= RF_rTargetFactorCount; j++) {
                maxValue = 0.0;
                maxClass = 0.0;
                for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
                  if (maxValue < RF_vimpEnsembleCLS[p][j][k][i]) {
                    maxValue = RF_vimpEnsembleCLS[p][j][k][i];
                    maxClass = (double) k;
                  }
                  RF_vimpEnsembleCLS[p][j][k][i] = NA_REAL;
                }
                RF_vimpEnsembleCLS[p][j][1][i] = maxClass;
              }
            }
            else {
              for (j = 1; j <= RF_rTargetFactorCount; j++) {
                for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
                  RF_vimpEnsembleCLS[p][j][k][i] = NA_REAL;
                }
              }
            }
          }
        }
      }
    }
    if (RF_rTargetNonFactorCount > 0) {
      if (!(RF_opt & OPT_VIMP_LEOB)) {
        for (i = 1; i <= obsSize; i++) {
          if(vimpDenom[i] > 0) {
            for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
              RF_vimpEnsembleRGR[p][j][i] = RF_vimpEnsembleRGR[p][j][i] / vimpDenom[i];
            }
          }
          else {
            for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
              RF_vimpEnsembleRGR[p][j][i] = NA_REAL;
            }
          }
        }
      }
    }
  }
  ensembleMRT = NULL;
  ensembleCLS = NULL;
  ensembleRGR = NULL;
  vimpMRTptr = NULL;
  vimpCLSptr = NULL;
  vimpRGRptr = NULL;
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    ensembleMRT =  RF_vimpEnsembleMRT[p];
  }
  else {
    if (RF_rTargetFactorCount > 0) {
      ensembleCLS = RF_vimpEnsembleCLS[p];
    }
    if (RF_rTargetNonFactorCount > 0) {
      ensembleRGR = RF_vimpEnsembleRGR[p];
    }    
  }
  if (!(RF_opt & OPT_VIMP_LEOB)) {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      vimpMRTptr = RF_vimpMRTptr[p];
    }
    else {
      if (RF_rTargetFactorCount > 0) {
        vimpCLSptr = RF_vimpCLSptr[p];
      }
      if (RF_rTargetNonFactorCount > 0) {
        vimpRGRptr = RF_vimpRGRptr[p];
      }    
    }
  }
  else {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      vimpMRTptr = RF_vimpMRTleo[treeID][p];
    }
    else {
      if (RF_rTargetFactorCount > 0) {
        vimpCLSptr = RF_vimpCLSleo[treeID][p];
      }
      if (RF_rTargetNonFactorCount > 0) {
        vimpRGRptr = RF_vimpRGRleo[treeID][p];
      }
    }
  }
  getPerformance(treeID,
                 mode,
                 obsSize,
                 responsePtr,
                 vimpDenom,
                 ensembleMRT,       
                 ensembleCLS,       
                 ensembleRGR,       
                 vimpMRTptr,        
                 vimpCLSptr,        
                 vimpRGRptr);       
  unstackImputeResponse(imputeFlag, obsSize, responsePtr);
}
void finalizeVimpPerformance(char       mode,
                             uint       rejectedTreeCount) {
  uint xVimpSize;
  double result;
  uint cumDenomCount;
  uint i, j, k, p;
  switch (mode) {
  case RF_PRED:
    if (RF_opt & OPT_VIMP_JOIN) {
      xVimpSize = 1;
    }
    else {
      xVimpSize = RF_intrPredictorSize;
    }
    break;
  default:
    if (RF_opt & OPT_VIMP_JOIN) {
      xVimpSize = 1;
    }
    else {
      xVimpSize = RF_intrPredictorSize;
    }
    break;
  }
  if (RF_opt & OPT_VIMP_LEOB) {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      for(p = 1; p <= xVimpSize; p++) {
        for (k = 1; k <= RF_eventTypeSize; k++) {
          cumDenomCount = 0;
          result = 0.0;
          for (i = 1; i <= RF_forestSize; i++) {
            if(!ISNA(RF_vimpMRTleo[i][p][k])) {
              if(!ISNA(RF_perfMRTleo[i][k])) {
                result += RF_vimpMRTleo[i][p][k] - RF_perfMRTleo[i][k];
                cumDenomCount ++;
              }
            }
          }
          if (cumDenomCount != 0) {
            RF_vimpMRTptr[p][k] = result / (double) cumDenomCount;
          }
          else {
            RF_vimpMRTptr[p][k] = NA_REAL;
          }
        }
      }
    }
    else {
      if (RF_rTargetFactorCount > 0) {
        for(p = 1; p <= xVimpSize; p++) {
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            for (k = 1; k <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              cumDenomCount = 0;
              result = 0.0;
              for (i = 1; i <= RF_forestSize; i++) {
                if(!ISNA(RF_vimpCLSleo[i][p][j][k])) {
                  if(!ISNA(RF_perfCLSleo[i][j][k])) {
                    result += RF_vimpCLSleo[i][p][j][k] - RF_perfCLSleo[i][j][k];
                    cumDenomCount ++;
                  }
                }
              }
              if (cumDenomCount != 0) {
                if ( k > 1) {
                  RF_vimpCLSptr[p][j][k] = M_E * result / (double) cumDenomCount;
                }
                else {
                  RF_vimpCLSptr[p][j][k] = result / (double) cumDenomCount;
                }
              }
              else {
                RF_vimpCLSptr[p][j][k] = NA_REAL;
              }
            }
          }
        }
      }
      if (RF_rTargetNonFactorCount > 0) {
        for(p = 1; p <= xVimpSize; p++) {
          for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
            cumDenomCount = 0;
            result = 0.0;
            for (i = 1; i <= RF_forestSize; i++) {
              if(!ISNA(RF_vimpRGRleo[i][p][j])) {
                if(!ISNA(RF_perfRGRleo[i][j])) {
                  result += RF_vimpRGRleo[i][p][j] - RF_perfRGRleo[i][j];
                  cumDenomCount ++;
                }
              }
            }
            if (cumDenomCount != 0) {
              RF_vimpRGRptr[p][j] = result / (double) cumDenomCount;
            }
            else {
              RF_vimpRGRptr[p][j] = NA_REAL;
            }
          }
        }
      }    
    }
  }
  else {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      for(p = 1; p <= xVimpSize; p++) {
        for (k = 1; k <= RF_eventTypeSize; k++) {
          if(!ISNA(RF_vimpMRTptr[p][k])) {
            if(!ISNA(RF_perfMRTptr[RF_forestSize][k])) {
              RF_vimpMRTptr[p][k] = RF_vimpMRTptr[p][k] - RF_perfMRTptr[RF_forestSize][k];
            }
          }
        }
      }
    }
    else {
      if (RF_rTargetFactorCount > 0) {
        for(p = 1; p <= xVimpSize; p++) {
          for (j = 1; j <= RF_rTargetFactorCount; j++) {          
            for (k = 1; k <= 1 + RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              if(!ISNA(RF_vimpCLSptr[p][j][k])) {
                if(!ISNA(RF_perfCLSptr[RF_forestSize][j][k])) {
                  RF_vimpCLSptr[p][j][k] = RF_vimpCLSptr[p][j][k] - RF_perfCLSptr[RF_forestSize][j][k];
                }
              }
            }
          }
        }
      }
      if (RF_rTargetNonFactorCount > 0) {
        for(p = 1; p <= xVimpSize; p++) {
          for (j = 1; j <= RF_rTargetNonFactorCount; j++) {          
            if(!ISNA(RF_vimpRGRptr[p][j])) {
              if(!ISNA(RF_perfRGRptr[RF_forestSize][j])) {
                RF_vimpRGRptr[p][j] = RF_vimpRGRptr[p][j] - RF_perfRGRptr[RF_forestSize][j];
              }
            }
          }
        }
      }
    }
  }
}
void stackVimpMembership(char mode, Terminal ***membership) {
  uint obsSize;
  (*membership) = NULL;
  if (RF_opt & OPT_VIMP) {
    switch (mode) {
    case RF_PRED:
      obsSize = RF_fobservationSize;
      break;
    default:
      obsSize = RF_observationSize;
      break;
    }
    *membership = (Terminal **) new_vvector(1, obsSize, NRUTIL_TPTR);
  }
}
void unstackVimpMembership(char mode, Terminal **membership) {
  uint obsSize;
  if (RF_opt & OPT_VIMP) {
    switch (mode) {
    case RF_PRED:
      obsSize = RF_fobservationSize;
      break;
    default:
      obsSize = RF_observationSize;
      break;
    }
    free_new_vvector(membership, 1, obsSize, NRUTIL_TPTR);
  }
}
void stackTreeEnsemble(char         mode,
                       uint         treeID,
                       uint       **treeDenom,
                       double   ****treeEnsembleMRT,
                       double  *****treeEnsembleCLS,
                       double   ****treeEnsembleRGR) {
  uint  obsSize;
  uint i,j,k;
  (*treeEnsembleMRT)     = (double ***) new_vvector(1, 1, NRUTIL_DPTR2);
  (*treeEnsembleMRT)[1]  = NULL;
  (*treeEnsembleCLS)     = (double****) new_vvector(1, 1, NRUTIL_DPTR3);
  (*treeEnsembleCLS)[1]  = NULL;
  (*treeEnsembleRGR)     = (double ***) new_vvector(1, 1, NRUTIL_DPTR2);
  (*treeEnsembleRGR)[1]  = NULL;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    break;
  default:
    obsSize = RF_observationSize;
    break;
  }
  *treeDenom = uivector(1, obsSize);
  for (i = 1; i <= obsSize; i++) {
    (*treeDenom)[i] = 0;
  }
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    (*treeEnsembleMRT)[1] = (double **) new_vvector(1, RF_eventTypeSize, NRUTIL_DPTR);
    for (k = 1; k <= RF_eventTypeSize; k++) {
      (*treeEnsembleMRT)[1][k] = dvector(1, obsSize);
      for (i = 1; i <= obsSize; i++) {
        (*treeEnsembleMRT)[1][k][i] = 0.0;
      }
    }
  }
  else {
    if (RF_rTargetFactorCount > 0) {
      if (!(RF_opt & OPT_PERF_CALB)) {
        (*treeEnsembleCLS)[1] = (double ***) new_vvector(1, RF_rTargetFactorCount, NRUTIL_DPTR2);
        for (j = 1; j <= RF_rTargetFactorCount; j++) {
          (*treeEnsembleCLS)[1][j] = (double **) new_vvector(1, 1, NRUTIL_DPTR);
          (*treeEnsembleCLS)[1][j][1]  = dvector(1, obsSize);
          for (i = 1; i <= obsSize; i++) {
            (*treeEnsembleCLS)[1][j][1][i] = 0.0;
          }
        }
      }
      else {
        (*treeEnsembleCLS)[1] = (double ***) new_vvector(1, RF_rTargetFactorCount, NRUTIL_DPTR2);
        for (j = 1; j <= RF_rTargetFactorCount; j++) {
          (*treeEnsembleCLS)[1][j] = (double **) new_vvector(1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]], NRUTIL_DPTR);
          for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
            (*treeEnsembleCLS)[1][j][k]  = dvector(1, obsSize);
            for (i = 1; i <= obsSize; i++) {
              (*treeEnsembleCLS)[1][j][k][i] = 0.0;
            }
          }
        }
      }
    }
    if (RF_rTargetNonFactorCount > 0) {
      (*treeEnsembleRGR)[1] = (double **) new_vvector(1, RF_rTargetNonFactorCount, NRUTIL_DPTR2);
      for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
        (*treeEnsembleRGR)[1][j] = dvector(1, obsSize);
        for (i = 1; i <= obsSize; i++) {
          (*treeEnsembleRGR)[1][j][i] = 0.0;
        }
      }
    }
  }
}
void unstackTreeEnsemble(char        mode,
                         uint        treeID,
                         uint       *treeDenom,
                         double   ***treeEnsembleMRT,
                         double  ****treeEnsembleCLS,
                         double   ***treeEnsembleRGR) {
  uint obsSize;
  uint j, k;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    break;
  default:
    obsSize = RF_observationSize;
    break;
  }
  free_uivector(treeDenom, 1, obsSize);
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    for (k = 1; k <= RF_eventTypeSize; k++) {
      free_dvector(treeEnsembleMRT[1][k], 1, obsSize);
    }
    free_new_vvector(treeEnsembleMRT[1], 1, RF_eventTypeSize, NRUTIL_DPTR);
  }
  else {
    if (RF_rTargetFactorCount > 0) {
      if (!(RF_opt & OPT_PERF_CALB)) {
        for (j = 1; j <= RF_rTargetFactorCount; j++) {
          free_dvector(treeEnsembleCLS[1][j][1], 1, obsSize);
          free_new_vvector(treeEnsembleCLS[1][j], 1, 1, NRUTIL_DPTR);
        }
        free_new_vvector(treeEnsembleCLS[1], 1, RF_rTargetFactorCount, NRUTIL_DPTR2);
      }
      else {
        for (j = 1; j <= RF_rTargetFactorCount; j++) {
          for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
            free_dvector(treeEnsembleCLS[1][j][k], 1, obsSize);
          }
          free_new_vvector(treeEnsembleCLS[1][j], 1, RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]], NRUTIL_DPTR);
        }
        free_new_vvector(treeEnsembleCLS[1], 1, RF_rTargetFactorCount, NRUTIL_DPTR2);
      }
    }
    if (RF_rTargetNonFactorCount > 0) {
      for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
        free_dvector(treeEnsembleRGR[1][j], 1, obsSize);
      }
      free_new_vvector(treeEnsembleRGR[1], 1, RF_rTargetNonFactorCount, NRUTIL_DPTR);
    }
  }
  free_new_vvector(treeEnsembleMRT, 1, 1, NRUTIL_DPTR2);
  free_new_vvector(treeEnsembleCLS, 1, 1, NRUTIL_DPTR3);
  free_new_vvector(treeEnsembleRGR, 1, 1, NRUTIL_DPTR2);
}
void updateVimpCalculations (char mode, uint b, uint intrIndex, Terminal **vimpMembership) {
  updateVimpEnsemble(mode, b, vimpMembership, intrIndex);
  if (RF_opt & OPT_VIMP_LEOB) {
    summarizeVimpPerformance(mode, b, intrIndex);
  }
}
void summarizeTreePerformance(char mode, uint treeID) {
  uint  obsSize;
  double    **responsePtr;
  uint      *treeDenom;
  double  ***treeEnsembleMRT;
  double ****treeEnsembleCLS;
  double  ***treeEnsembleRGR;
  double   *perfMRTleo;
  double  **perfCLSleo;
  double   *perfRGRleo;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    responsePtr = RF_fresponse[treeID];
    break;
  default:
    obsSize = RF_observationSize;
    responsePtr = RF_response[treeID];
    break;
  }
  stackTreeEnsemble(mode, treeID, &treeDenom, &treeEnsembleMRT, &treeEnsembleCLS, &treeEnsembleRGR);
  updateTreeEnsemble(mode, treeID, treeDenom, treeEnsembleMRT, treeEnsembleCLS, treeEnsembleRGR);
  perfMRTleo = NULL;
  perfCLSleo = NULL;
  perfRGRleo = NULL;
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    perfMRTleo = RF_perfMRTleo[treeID];
  }
  else {
    if (RF_rTargetFactorCount > 0) {
      perfCLSleo = RF_perfCLSleo[treeID];
    }
    if (RF_rTargetNonFactorCount > 0) {
      perfRGRleo = RF_perfRGRleo[treeID];
    }    
  }
  getPerformance(treeID,
                 mode,
                 obsSize,
                 responsePtr,
                 treeDenom,
                 treeEnsembleMRT[1],  
                 treeEnsembleCLS[1],  
                 treeEnsembleRGR[1],  
                 perfMRTleo,          
                 perfCLSleo,          
                 perfRGRleo);         
  unstackTreeEnsemble(mode, treeID, treeDenom, treeEnsembleMRT, treeEnsembleCLS, treeEnsembleRGR);
}
void updatePartialCalculations (uint       treeID,
                                uint       pVarIdx,
                                Terminal **partialMembership) {
  Terminal *terminalNode;
  uint  *membershipIndex;
  uint   membershipSize;
  uint   i, j, k;
  uint   ii;
  if (RF_opt & OPT_OENS) {
    membershipSize  = RF_oobSize[treeID];
    membershipIndex = RF_oobMembershipIndex[treeID];
  }
  else {
    membershipSize  = RF_observationSize;
    membershipIndex = RF_identityMembershipIndex;
  }
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (RF_eventTypeSize > 1) {
      if (RF_partialType == RF_PART_YRLS) {
        for (i = 1; i <= membershipSize; i++) {
          ii = membershipIndex[i];
          terminalNode = partialMembership[ii];
          for (j = 1; j <= RF_eventTypeSize; j++) {
            RF_partSURVptr[pVarIdx][j][1][ii] += terminalNode -> mortality[j];
          }
        }
      }
      else if (RF_partialType == RF_PART_CIFN) {
        for (i = 1; i <= membershipSize; i++) {
          ii = membershipIndex[i];
          terminalNode = partialMembership[ii];
          for (j = 1; j <= RF_eventTypeSize; j++) {
            for (k = 1; k <= RF_partialTimeLength; k++) {
              RF_partSURVptr[pVarIdx][j][k][ii] += terminalNode -> CIF[j][k];
            }
          }
        }
      }
      else if (RF_partialType == RF_PART_CHFN) {
        for (i = 1; i <= membershipSize; i++) {
          ii = membershipIndex[i];
          terminalNode = partialMembership[ii];
          for (j = 1; j <= RF_eventTypeSize; j++) {
            for (k = 1; k <= RF_partialTimeLength; k++) {
              RF_partSURVptr[pVarIdx][j][k][ii] += terminalNode -> CSH[j][k];
            }
          }
        }
      }
    }   
    else {
      if (RF_partialType == RF_PART_MORT) {
        for (i = 1; i <= membershipSize; i++) {
          ii = membershipIndex[i];
          terminalNode = partialMembership[ii];
          RF_partSURVptr[pVarIdx][1][1][ii] += terminalNode -> mortality[1];
        }
      }
      else if (RF_partialType == RF_PART_NLSN) {
        for (i = 1; i <= membershipSize; i++) {
          ii = membershipIndex[i];
          terminalNode = partialMembership[ii];
          for (k = 1; k <= RF_partialTimeLength; k++) {
            RF_partSURVptr[pVarIdx][1][k][ii] += terminalNode -> nelsonAalen[k];
          }
        }
      }
      else if (RF_partialType == RF_PART_SURV) {
        for (i = 1; i <= membershipSize; i++) {
          ii = membershipIndex[i];
          terminalNode = partialMembership[ii];
          for (k = 1; k <= RF_partialTimeLength; k++) {
            RF_partSURVptr[pVarIdx][1][k][ii] += terminalNode -> survival[k];
          }
        }
      }
    }
  }
  else {
    if (RF_rTargetFactorCount > 0) {
      for (i = 1; i <= membershipSize; i++) {
        ii = membershipIndex[i];
        terminalNode = partialMembership[ii];
        for (j = 1; j <= RF_rTargetFactorCount; j++) {
          for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
            RF_partCLASptr[pVarIdx][j][k+1][ii] += (double) (terminalNode -> multiClassProb)[RF_rFactorMap[RF_rTargetFactor[j]]][k] / (double) (terminalNode -> membrCount);
          }
        }
      }
    }
    if (RF_rTargetNonFactorCount > 0) {
      for (i = 1; i <= membershipSize; i++) {
        ii = membershipIndex[i];
        terminalNode = partialMembership[ii];
        for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
          RF_partREGRptr[pVarIdx][j][ii] += (terminalNode -> meanResponse)[RF_rNonFactorMap[RF_rTargetNonFactor[j]]];
        }
      }
    }
  }
}
void summarizePartialCalculations(uint       treeID,
                                  uint       pVarIdx) {
  uint *ensembleDen;
  uint  membershipSize;
  double    maxValue, maxClass;
  uint i, j, k;
  membershipSize  = RF_observationSize;
  ensembleDen = RF_oobEnsembleDen;
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (RF_eventTypeSize > 1) {
      if (RF_partialType == RF_PART_YRLS) {
        for (i = 1; i <= membershipSize; i++) {
          if (ensembleDen[i] > 0) {
            for (j = 1; j <= RF_eventTypeSize; j++) {
              RF_partSURVptr[pVarIdx][j][1][i] = RF_partSURVptr[pVarIdx][j][1][i] / ensembleDen[i];
            }
          }
        }
      }
      else if (RF_partialType == RF_PART_CIFN) {
        for (i = 1; i <= membershipSize; i++) {
          if (ensembleDen[i] > 0) {
            for (j = 1; j <= RF_eventTypeSize; j++) {
              for (k = 1; k <= RF_partialTimeLength; k++) {
                RF_partSURVptr[pVarIdx][j][k][i] = RF_partSURVptr[pVarIdx][j][k][i] / ensembleDen[i];
              }
            }
          }
        }
      }
      else if (RF_partialType == RF_PART_CHFN) {
        for (i = 1; i <= membershipSize; i++) {
          if (ensembleDen[i] > 0) {
            for (j = 1; j <= RF_eventTypeSize; j++) {
              for (k = 1; k <= RF_partialTimeLength; k++) {
                RF_partSURVptr[pVarIdx][j][k][i] = RF_partSURVptr[pVarIdx][j][k][i] / ensembleDen[i];
              }
            }
          }
        }
      }
    }   
    else {
      if (RF_partialType == RF_PART_MORT) {
        for (i = 1; i <= membershipSize; i++) {
          if (ensembleDen[i] > 0) {
            RF_partSURVptr[pVarIdx][1][1][i] = RF_partSURVptr[pVarIdx][1][1][i] / ensembleDen[i];
          }
        }
      }
      else if (RF_partialType == RF_PART_NLSN) {
        for (i = 1; i <= membershipSize; i++) {
          if (ensembleDen[i] > 0) {
            for (k = 1; k <= RF_partialTimeLength; k++) {
              RF_partSURVptr[pVarIdx][1][k][i] = RF_partSURVptr[pVarIdx][1][k][i] / ensembleDen[i];
            }
          }
        }
      }
      else if (RF_partialType == RF_PART_SURV) {
        for (i = 1; i <= membershipSize; i++) {
          if (ensembleDen[i] > 0) {
            for (k = 1; k <= RF_partialTimeLength; k++) {
              RF_partSURVptr[pVarIdx][1][k][i] =  RF_partSURVptr[pVarIdx][1][k][i] / ensembleDen[i];
            }
          }
        }
      }
    }
  }
  else {
    if (RF_rTargetFactorCount > 0) {
      for (i = 1; i <= membershipSize; i++) {
        if (ensembleDen[i] > 0) {
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            maxValue = 0;
            maxClass = 0;
            for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {
              RF_partCLASptr[pVarIdx][j][k+1][i] = RF_partCLASptr[pVarIdx][j][k+1][i] / ensembleDen[i];
              if (maxValue < RF_partCLASptr[pVarIdx][j][k+1][i]) {
                maxValue = RF_partCLASptr[pVarIdx][j][k+1][i];
                maxClass = (double) k;
              }
            }
            RF_partCLASptr[pVarIdx][j][1][i] = maxClass;
          }
        }
        else {
          for (j = 1; j <= RF_rTargetFactorCount; j++) {
            for (k = 1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTargetFactor[j]]]; k++) {          
              RF_partCLASptr[pVarIdx][j][k+1][i] = NA_REAL;
            }
            RF_partCLASptr[pVarIdx][j][1][i] = NA_REAL;
          }
        }
      }
    }
    if (RF_rTargetNonFactorCount > 0) {
      for (i = 1; i <= membershipSize; i++) {
        if (ensembleDen[i] > 0) {
          for (j = 1; j <= RF_rTargetNonFactorCount; j++) {
            RF_partREGRptr[pVarIdx][j][i] = RF_partREGRptr[pVarIdx][j][i]  / ensembleDen[i];
          }
        }
      }
    }
  }
}
SEXP rfsrc(char mode, int seedValue, uint traceFlag) {
  uint sexpIndex;
  uint previousTreeID;
  ulong offset;
  uint i, j, r;
  int vimpCount, b, p;
  uint seedValueLC;
  setUserTraceFlag(traceFlag);
  seedValueLC    = 0; 
  if (RF_nImpute < 1) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Parameter verification failed.");
    RFprintf("\nRF-SRC:  Number imputations must be greater than zero:  %10d \n", RF_forestSize);
    RFprintf("\nRF-SRC:  The application will now exit.\n");
    return R_NilValue;
  }
  if (RF_forestSize < 1) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Parameter verification failed.");
    RFprintf("\nRF-SRC:  Number of bootstrap iterations must be greater than zero:  %10d \n", RF_forestSize);
    RFprintf("\nRF-SRC:  The application will now exit.\n");
    return R_NilValue;
  }
  if (RF_observationSize < 1) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Parameter verification failed.");
    RFprintf("\nRF-SRC:  Number of individuals must be greater than one:  %10d \n", RF_observationSize);
    RFprintf("\nRF-SRC:  The application will now exit.\n");
    return R_NilValue;
  }
  if (RF_xSize < 1) {
    RFprintf("\nRF-SRC:  *** ERROR *** ");
    RFprintf("\nRF-SRC:  Parameter verification failed.");
    RFprintf("\nRF-SRC:  Number of parameters must be greater than zero:  %10d \n", RF_xSize);
    RFprintf("\nRF-SRC:  The application will now exit.\n");
    return R_NilValue;
  }
#ifdef _OPENMP
  if (RF_numThreads < 0) {
    RF_numThreads = omp_get_max_threads();
  }
  else {
    RF_numThreads = (RF_numThreads < omp_get_max_threads()) ? (RF_numThreads) : (omp_get_max_threads());
  }
#endif
  stackIncomingArrays(mode);
  stackPreDefinedCommonArrays();
  switch (mode) {
  case RF_PRED:
    stackPreDefinedPredictArrays();
    break;
  case RF_REST:
    stackPreDefinedRestoreArrays();
    break;
  default:
    stackPreDefinedGrowthArrays();
    break;
  }
  initializeTimeArrays(mode);
  stackFactorArrays();
  stackMissingArrays(mode);
  if (RF_statusIndex > 0) {
    stackCompetingArrays(mode);
  }
  if (RF_rFactorCount > 0) {
    stackClassificationArrays(mode);
  }
  sexpIndex = 0;
  stackDefinedOutputObjects(mode,
                            sexpString,
                            & RF_root,
                            & RF_tLeafCount_,
                            & RF_proximity_,
                            & RF_seed_,
                            & RF_imputation_,
                            & RF_sImputeResponsePtr,
                            & RF_sImputePredictorPtr,
                            & RF_varUsed_,
                            & RF_varUsedPtr,
                            & RF_splitDepth_,
                            & sexpIndex,
                            & RF_stackCount,
                              sexpVector);
  verifyAndRegisterCustomSplitRules();
  switch (mode) {
  case RF_GROW:
    RF_theoreticalMaxtLeafCount = uivector(1, RF_forestSize);
    j = ((RF_minimumNodeSize - 1 ) << 1);
    if (RF_bootstrapSize > j) { 
      j = RF_bootstrapSize - j;
    }
    else {
      j = 1;
    }
    for (i = 1; i <= RF_forestSize; i++) {
      RF_theoreticalMaxtLeafCount[i] = j;
    }
    RF_totalTerminalCount = j * RF_forestSize;
      stackForestOutputObjects(mode,
                               & RF_treeID_,         
                               & RF_nodeID_,         
                               & RF_parmID_,         
                               & RF_contPT_,         
                               & RF_mwcpSZ_,         
                               & RF_mwcpPT_,         
                               & RF_mwcpCount,       
                               & sexpIndex,
                                 sexpString,
                                 sexpVector);
    break;
  default:
    previousTreeID = j = 0;
    for (ulong ui = 1; ui <= RF_totalNodeCount; ui++) {
      if ((RF_treeID_[ui] > 0) && (RF_treeID_[ui] <= RF_forestSize)) {
        if (RF_treeID_[ui] != previousTreeID) {
          previousTreeID = RF_restoreTreeID[++j] = RF_treeID_[ui];
          RF_restoreTreeOffset[RF_treeID_[ui]] = ui;
        }
        RF_nodeCount[RF_treeID_[ui]] ++;
        RF_mwcpCount[RF_treeID_[ui]] += RF_mwcpSZ_[ui];
      }
      else {
        RFprintf("\nRF-SRC:  Diagnostic Trace of Tree Record:  \n");
        RFprintf("\nRF-SRC:      treeID     nodeID     parmID       spltPT     mwcpSZ ");
        RFprintf("\nRF-SRC:  %10d %10d %10d %12.4f %10d \n", RF_treeID_[ui], RF_nodeID_[ui], RF_parmID_[ui], RF_contPT_[ui], RF_mwcpSZ_[ui]);
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Invalid forest input record at line:  %20lu", ui);
        RFprintf("\nRF-SRC:  Please Contact Technical Support.");
        RFprintf("\nRF-SRC:  The application will now exit.\n");
        return R_NilValue;
      }
    }
    offset = 0;
    for (b = 1; b <= RF_forestSize; b++) {
      if (RF_mwcpCount[RF_restoreTreeID[b]] > 0) {
        RF_restoreMWCPOffset[RF_restoreTreeID[b]] = offset;
        RF_mwcpPtr[RF_restoreTreeID[b]] = RF_mwcpPT_ + offset;
        offset = offset + RF_mwcpCount[RF_restoreTreeID[b]];
      }
      else {
        RF_restoreMWCPOffset[RF_restoreTreeID[b]] = 0;
        RF_mwcpPtr[RF_restoreTreeID[b]] = NULL;
      }
    }
    RF_totalTerminalCount = 0;
    for (b = 1; b <= RF_forestSize; b++) {
      RF_tLeafCount[b] = (RF_nodeCount[b] + 1) >> 1;
      RF_totalTerminalCount += (ulong) RF_tLeafCount[b];
    }
    break;
  }
  sexpIndex =
    stackTNQualitativeOutputObjects(mode,
                                    & RF_RMBR_ID_,
                                    & RF_AMBR_ID_,
                                    & RF_TN_RCNT_,
                                    & RF_TN_ACNT_,
                                    sexpIndex,
                                    sexpString,
                                    sexpVector);
  sexpIndex =
    stackTNQuantitativeOutputObjects(mode,
                                     & RF_TN_SURV_,         
                                     & RF_TN_MORT_,         
                                     & RF_TN_NLSN_,         
                                     & RF_TN_CSHZ_,         
                                     & RF_TN_CIFN_,         
                                     & RF_TN_REGR_,         
                                     & RF_TN_CLAS_,         
                                     sexpIndex,
                                     sexpString,
                                     sexpVector);
  sexpIndex =
    stackStatisticalOutputObjects(mode,
                                  & RF_spltST_,         
                                  & RF_spltVR_,         
                                  & RF_uspvST_,         
                                  & RF_mtryID_,         
                                  & RF_mtryST_,         
                                  & sexpIndex,
                                    sexpString,
                                    sexpVector);
#ifdef _OPENMP
  ran1A = &randomChainParallel;
  ran1B = &randomUChainParallel;
  ran1C = &randomUChainParallelCov;
  randomSetChain = &randomSetChainParallel;
  randomSetUChain = &randomSetUChainParallel;
  randomSetUChainCov = &randomSetUChainParallelCov;
  randomGetChain = &randomGetChainParallel;
  randomGetUChain = &randomGetUChainParallel;
  randomGetUChainCov = &randomGetUChainParallelCov;
  randomStack(RF_forestSize, RF_xSize);
  if (mode == RF_GROW) {
    seedValueLC = abs(seedValue);
    lcgenerator(&seedValueLC, TRUE);
    for (b = 1; b <= RF_forestSize; b++) {
      lcgenerator(&seedValueLC, FALSE);
      lcgenerator(&seedValueLC, FALSE);
      while(seedValueLC == 0) {
        lcgenerator(&seedValueLC, FALSE);
      }
      randomSetChain(b, -seedValueLC);
    }
  }
  for (b = 1; b <= RF_forestSize; b++) {
    lcgenerator(&seedValueLC, FALSE);
    lcgenerator(&seedValueLC, FALSE);
    while(seedValueLC == 0) {
      lcgenerator(&seedValueLC, FALSE);
    }
    randomSetUChain(b, -seedValueLC);
  }
  for (b = 1; b <= RF_forestSize; b++) {
    lcgenerator(&seedValueLC, FALSE);
    lcgenerator(&seedValueLC, FALSE);
    while(seedValueLC == 0) {
      lcgenerator(&seedValueLC, FALSE);
    }
    randomSetUChainCov(b, -seedValueLC);
  }
#else
  ran1A = &randomChainSerial;
  ran1B = &randomUChainSerial;
  ran1C = &randomUChainSerialCov;
  randomSetChain = &randomSetChainSerial;
  randomSetUChain = &randomSetUChainSerial;
  randomSetUChainCov = &randomSetUChainSerialCov;
  randomGetChain = &randomGetChainSerial;
  randomGetUChain = &randomGetUChainSerial;
  randomGetUChainCov = &randomGetUChainSerialCov;
  randomStack(1, 1);
  if (mode == RF_GROW) {
    seedValueLC = abs(seedValue);
    lcgenerator(&seedValueLC, TRUE);
    lcgenerator(&seedValueLC, FALSE);
    lcgenerator(&seedValueLC, FALSE);
    while(seedValueLC == 0) {
      lcgenerator(&seedValueLC, FALSE);
    }
    randomSetChain(1, -seedValueLC);
  }
  lcgenerator(&seedValueLC, FALSE);
  lcgenerator(&seedValueLC, FALSE);
  while(seedValueLC == 0) {
    lcgenerator(&seedValueLC, FALSE);
  }
  randomSetUChain(1, -seedValueLC);
  lcgenerator(&seedValueLC, FALSE);
  lcgenerator(&seedValueLC, FALSE);
  while(seedValueLC == 0) {
    lcgenerator(&seedValueLC, FALSE);
  }
  randomSetUChainCov(1, -seedValueLC);
#endif
  switch (mode) {
  case RF_GROW:
    break;
  default:
    for (b = 1; b <= RF_forestSize; b++) {
      if(RF_seed_[b] >= 0) {
        RFprintf("\nRF-SRC:  *** ERROR *** ");
        RFprintf("\nRF-SRC:  Parameter verification failed.");
        RFprintf("\nRF-SRC:  Forest random seed element must be less than zero:  %10d \n", RF_seed_[b]);
        RFprintf("\nRF-SRC:  The application will now exit.\n");
        return R_NilValue;
      }
    }
    break;
  }
  for (r = 1; r <= RF_nImpute; r++) {
    if (getUserTraceFlag()) {
      if (RF_nImpute == 1) {
      }
      else {
        RFprintf("\nImpute Iteration:  %6d", r);
      }
    }
    if (r == RF_nImpute) {
#ifdef _OPENMP
      switch (mode) {
      case RF_GROW:
        if (RF_opt & OPT_SEED) {
          for (b = 1; b <= RF_forestSize; b++) {
            if (r > 1) {
              lcgenerator(&seedValueLC, FALSE);
              lcgenerator(&seedValueLC, FALSE);
              while(seedValueLC == 0) {
                lcgenerator(&seedValueLC, FALSE);
              }
              randomSetChain(b, -seedValueLC);
            }
            RF_seed_[b] = randomGetChain(b);
          }
        }
        break;
      default:
        for (b = 1; b <= RF_forestSize; b++) {
          randomSetChain(b , RF_seed_[b]);
        }
        break;
      }
#else
      switch (mode) {
      case RF_GROW:
        if (RF_opt & OPT_SEED) {
          if (r > 1) {
            lcgenerator(&seedValueLC, FALSE);
            lcgenerator(&seedValueLC, FALSE);
            while(seedValueLC == 0) {
              lcgenerator(&seedValueLC, FALSE);
            }
            randomSetChain(1, -seedValueLC);
          }
          RF_seed_[1] = randomGetChain(1);
        }
        break;
      default:
        randomSetChain(1 , RF_seed_[1]);
        break;
      }
#endif
    }  
    for(b = 1; b <= RF_forestSize; b++) {
      RF_serialTreeIndex[b] = 0;
    }
    RF_serialTreeCount = 0;
    if (r == RF_nImpute) {
      stackAuxTNQualitativeStructures(mode,
                                      RF_RMBR_ID_,
                                      RF_AMBR_ID_,
                                      RF_TN_RCNT_,
                                      RF_TN_ACNT_);
      stackAuxTNQuantitativeStructures(mode,
                                       RF_TN_SURV_,         
                                       RF_TN_MORT_,         
                                       RF_TN_NLSN_,         
                                       RF_TN_CSHZ_,         
                                       RF_TN_CIFN_,         
                                       RF_TN_REGR_,         
                                       RF_TN_CLAS_);        
    }  
    if (getUserTraceFlag()) {
      RF_userTimeStart = RF_userTimeSplit = time(NULL);
    }
    if (RF_numThreads > 0) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
      for (b = 1; b <= RF_forestSize; b++) {
        acquireTree(mode, r, b);
      }
    }
    else {
      for (b = 1; b <= RF_forestSize; b++) {
        acquireTree(mode, r, b);
      }
    }
    if (r == RF_nImpute) {
      RF_rejectedTreeCount = RF_validTreeCount = RF_stumpedTreeCount = 0;
      for (b = 1; b <= RF_forestSize; b++) {
        if (RF_tLeafCount[b] == 0) {
          RF_rejectedTreeCount ++;
        }
        else {
          RF_validTreeCount ++;
          if (RF_tLeafCount[b] == 1) {
            RF_stumpedTreeCount ++;
          }
        }
      }
      unstackAuxTNQualitativeStructures(mode);
      unstackAuxTNQuantitativeStructures(mode);
      if (RF_opt & OPT_PROX) {
        finalizeProximity(mode);
      }
    }  
    if (RF_opt & OPT_MISS) {
      switch (mode) {
      case RF_PRED:
        imputeSummary(RF_PRED, ACTIVE);
        break;
      default:
        if (r == 1) {
          if ( (!(RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
               ( (RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ) {
            imputeSummary(RF_GROW, FALSE);
            if (RF_timeIndex > 0) {
              if (RF_mTimeFlag == TRUE) {
                imputeMultipleTime(FALSE);
              }
            }
          }
          else {
            imputeSummary(RF_GROW, ACTIVE);
            if (RF_timeIndex > 0) {
              if (RF_mTimeFlag == TRUE) {
                imputeMultipleTime(ACTIVE);
              }
            }
          }
        }  
        else {
          if (r < RF_nImpute) {
            if ( (!(RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
                 ( (RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ) {
              imputeSummary(RF_GROW, FALSE);
              if (RF_timeIndex > 0) {
                if (RF_mTimeFlag == TRUE) {
                  imputeMultipleTime(FALSE);
                }
              }
            }
            else {
              imputeSummary(RF_GROW, ACTIVE);
              if (RF_timeIndex > 0) {
                if (RF_mTimeFlag == TRUE) {
                  imputeMultipleTime(ACTIVE);
                }
              }
            }
          }
          else {
            if (RF_opt & OPT_IMPU_ONLY) {
              if ( (!(RF_opt & OPT_BOOT_TYP1) && !(RF_opt & OPT_BOOT_TYP2)) ||
                   ( (RF_opt & OPT_BOOT_TYP1) &&  (RF_opt & OPT_BOOT_TYP2)) ) {
                imputeSummary(RF_GROW, TRUE);
                if (RF_timeIndex > 0) {
                  if (RF_mTimeFlag == TRUE) {
                    imputeMultipleTime(TRUE);
                  }
                }
              }
              else {
                imputeSummary(RF_GROW, ACTIVE);
                if (RF_timeIndex > 0) {
                  if (RF_mTimeFlag == TRUE) {
                    imputeMultipleTime(ACTIVE);
                  }
                }
              }
            }
            else {
            }
          }
        }
        break;
      }
    }  
    if (r < RF_nImpute) {
    if (RF_opt & OPT_MISS) {
      for (b = 1; b <= RF_forestSize; b++) {
        for (j = 1; j <= RF_tLeafCount[b]; j++) {
          freeTerminal(RF_tTermList[b][j]);
        }
        unstackTermList(b);
        free_new_vvector(RF_tTermMembership[b], 1, RF_observationSize, NRUTIL_TPTR);
        if (mode == RF_PRED) {
          free_new_vvector(RF_ftTermMembership[b], 1, RF_fobservationSize, NRUTIL_TPTR);
        }
      }
    }
    }
    if (getUserTraceFlag()) {
      RFprintf("\n\n");
    }
  }  
  if (RF_rejectedTreeCount < RF_forestSize) {
    if (RF_opt & OPT_VIMP) {
      if (RF_opt & OPT_VIMP_JOIN) {
        vimpCount = 1;
      }
      else {
        vimpCount = RF_intrPredictorSize;
      }
      if (RF_opt & OPT_VIMP_LEOB) {
      }
      else {
        if (RF_numThreads > 0) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
          for (p = 1; p <= vimpCount; p++) {
            summarizeVimpPerformance(mode, 0, p);
          }
        }
        else {
          for (p = 1; p <= vimpCount; p++) {
            summarizeVimpPerformance(mode, 0, p);
          }
        }
      }
      finalizeVimpPerformance(mode, RF_rejectedTreeCount);
    }
    if (RF_optHigh & OPT_PART_PLOT) {
      if (RF_numThreads > 0) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
        for (p = 1; p <= RF_partialLength; p++) {
          summarizePartialCalculations(0, p);
        }
      }
      else {
        for (p = 1; p <= RF_partialLength; p++) {
          summarizePartialCalculations(0, p);
        }
      }
    }
    finalizeEnsembleEstimates(mode);
    if (RF_opt & OPT_MISS) {
      for (b = 1; b <= RF_forestSize; b++) {
        for (j = 1; j <= RF_tLeafCount[b]; j++) {
          freeTerminal(RF_tTermList[b][j]);
        }
        unstackTermList(b);
        free_new_vvector(RF_tTermMembership[b], 1, RF_observationSize, NRUTIL_TPTR);
        if (mode == RF_PRED) {
          free_new_vvector(RF_ftTermMembership[b], 1, RF_fobservationSize, NRUTIL_TPTR);
        }
      }
    }
    if (RF_opt & (OPT_VARUSED_F | OPT_VARUSED_T)) {
      if (RF_opt & OPT_VARUSED_F) {
        for (j = 1; j <= RF_xSize; j++) {
          RF_varUsed_[j] = 0;
          for (i = 1; i <= RF_forestSize; i++) {
            RF_varUsed_[j] += RF_varUsedPtr[i][j];
          }
        }
      }
      else {
      }
    }
    if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
      if (RF_opt & OPT_SPLDPTH_F) {
        for (j = 1; j <= RF_xSize; j++) {
          for (i = 1; i <= RF_observationSize; i++) {
            RF_splitDepthPtr[1][j][i] = RF_splitDepthPtr[1][j][i] / (RF_forestSize - RF_rejectedTreeCount);
          }
        }
      }
      else {
      }
    }
  }  
  else {
    RFprintf("\nRF-SRC:  *** WARNING *** ");
    RFprintf("\nRF-SRC:  Insufficient trees for analysis.  \n");
  }
  unstackAuxStatisticalStructures(mode);
  unstackDefinedOutputObjects(mode, RF_root);
  switch (mode) {
  case RF_GROW:
    free_uivector(RF_theoreticalMaxtLeafCount, 1, RF_forestSize);
    break;
  default:
    break;
  }
  if (RF_statusIndex > 0) {
    unstackCompetingArrays(mode);
  }
  if (RF_rFactorCount > 0) {
    unstackClassificationArrays(mode);
  }
  unstackMissingArrays(mode);
  unstackFactorArrays();
  switch (mode) {
  case RF_PRED:
    unstackPreDefinedPredictArrays();
    break;
  case RF_REST:
    unstackPreDefinedRestoreArrays();
    break;
  default:
    unstackPreDefinedGrowthArrays();
    break;
  }
  unstackPreDefinedCommonArrays();
  unstackIncomingArrays(mode);
#ifdef _OPENMP
  randomUnstack(RF_forestSize, RF_xSize);
#else
  randomUnstack(1, 1);
#endif
  UNPROTECT(RF_stackCount + 2);
  return sexpVector[RF_OUTP_ID];
}
