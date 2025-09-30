
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "random.h"
#include "nrutil.h"
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
int  *ran1D_iy;
int **ran1D_iv;
int      *seed1AValue;
int      *seed1BValue;
int      *seed1CValue;
int      *seed1DValue;
void randomStack(uint bSize, uint bnpSize) {
  uint b, p;
  ran1A_iy = ivector(1, bSize);
  ran1A_iv = imatrix(1, bSize, 1, NTAB);
  ran1B_iy = ivector(1, bSize);
  ran1B_iv = imatrix(1, bSize, 1, NTAB);
  for (b = 1; b <= bSize; b++) {
    ran1A_iy[b] = 0;
    ran1B_iy[b] = 0;
  }
  seed1AValue = ivector(1, bSize);
  seed1BValue = ivector(1, bSize);
  if (bnpSize > 0) {
    ran1D_iy = ivector(1, bnpSize);
    ran1D_iv = imatrix(1, bnpSize, 1, NTAB);
    for (p = 1; p <= bnpSize; p++) {
      ran1D_iy[p] = 0;
    }
    seed1DValue = ivector(1, bnpSize);
  }
}
void randomUnstack(uint bSize, uint bnpSize) {
  free_ivector(ran1A_iy, 1, bSize);
  free_imatrix(ran1A_iv, 1, bSize, 1, NTAB);
  free_ivector(ran1B_iy, 1, bSize);
  free_imatrix(ran1B_iv, 1, bSize, 1, NTAB);
  free_ivector(seed1AValue, 1, bSize);
  free_ivector(seed1BValue, 1, bSize);
  if (bnpSize > 0) {
    free_ivector(ran1D_iy, 1, bnpSize);
    free_imatrix(ran1D_iv, 1, bnpSize, 1, NTAB);
    free_ivector(seed1DValue, 1, bnpSize);
  }
}
void randomSetChainParallel(uint b, int value) {
  seed1AValue[b] = value;
}
void randomSetUChainParallel(uint b, int value) {
  seed1BValue[b] = value;
}
void randomSetUChainParallelVimp(uint b, int value) {
  seed1CValue[b] = value;
}
void randomSetChainParallelVimp(uint p, int value) {
  seed1DValue[p] = value;
}
void randomSetChainSerial(uint b, int value) {
  seed1AValue[1] = value;
}
void randomSetUChainSerial(uint b, int value) {
  seed1BValue[1] = value;
}
void randomSetUChainSerialVimp(uint b, int value) {
  seed1CValue[1] = value;
}
void randomSetChainSerialVimp(uint p, int value) {
  seed1DValue[1] = value;
}
int randomGetChainParallel(uint b) {
  return seed1AValue[b];
}
int randomGetUChainParallel(uint b) {
  return seed1BValue[b];
}
int randomGetUChainParallelVimp(uint b) {
  return seed1CValue[b];
}
int randomGetChainParallelVimp(uint p) {
  return seed1DValue[p];
}
int randomGetChainSerial(uint b) {
  return seed1AValue[1];
}
int randomGetUChainSerial(uint b) {
  return seed1BValue[1];
}
int randomGetUChainSerialVimp(uint b) {
  return seed1CValue[1];
}
int randomGetChainSerialVimp(uint p) {
  return seed1DValue[1];
}
float randomChainParallel(uint b) {
  return  ran1_generic(& ran1A_iy[b], ran1A_iv[b], & seed1AValue[b]);
}
float randomUChainParallel(uint b) {
  return  ran1_generic(& ran1B_iy[b], ran1B_iv[b], & seed1BValue[b]);
}
float randomUChainParallelVimp(uint b) {
  return  ran1_generic(& ran1C_iy[b], ran1C_iv[b], & seed1CValue[b]);
}
float randomChainParallelVimp(uint p) {
  return  ran1_generic(& ran1D_iy[p], ran1D_iv[p], & seed1DValue[p]);
}
float randomChainSerial(uint b) {
  return  ran1_generic(& ran1A_iy[1], ran1A_iv[1], & seed1AValue[1]);
}
float randomUChainSerial(uint b) {
  return  ran1_generic(& ran1B_iy[1], ran1B_iv[1], & seed1BValue[1]);
}
float randomUChainSerialVimp(uint b) {
  return  ran1_generic(& ran1C_iy[1], ran1C_iv[1], & seed1CValue[1]);
}
float randomChainSerialVimp(uint p) {
  return  ran1_generic(& ran1D_iy[1], ran1D_iv[1], & seed1DValue[1]);
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
