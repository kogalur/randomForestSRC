
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "factorOps.h"
#include "nrutil.h"
#include "error.h"
Factor *makeFactor(uint r, char bookFlag) {
  uint i;
  Factor *f = (Factor*) gblock((size_t) sizeof(Factor));
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
void freeFactor(Factor *f) {
  if (f -> r > 1) {
    unbookFactor(f);
    if (f -> r <= MAX_EXACT_LEVEL) {
      free_uivector(f -> cardinalGroupSize, 1, (f -> cardinalGroupCount) + 1);
    }
    else {
      free_dvector(f -> cardinalGroupSize, 1, (f -> cardinalGroupCount) + 1);
    }
  }
  free_gblock(f, (size_t) sizeof(Factor));
}
char bookFactor(Factor *f) {
  uint i, j;
  uint row;
  char result;
  if (((f -> r) < 2) || ((f -> r) > MAX_EXACT_LEVEL)) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Minimum or Maximum number of factor levels violated in bookFactor(). ");
    RF_nativeError("\nRF-SRC:  Requested %10d, Minimum Allowed %10d, Maximum Allowed %10d ", f -> r, 2, MAX_EXACT_LEVEL);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  if (f -> cardinalGroupBinary == NULL) {
    uint *leftLevel = uivector(1, f -> cardinalGroupCount);
    f -> cardinalGroupBinary = (uint **) new_vvector(1, f -> cardinalGroupCount, NRUTIL_UPTR);
    for (i=1; i <= f -> cardinalGroupCount; i++) {
      (f -> cardinalGroupBinary)[i] = uivector(1, ((uint*) f -> cardinalGroupSize)[i]);
      row = 0;
      for (j = 1; j <= f -> cardinalGroupCount; j++) {
        leftLevel[j] = 0;
      }
      bookPair(f -> r , i, 1, &row, leftLevel, f);
    }
    free_uivector(leftLevel, 1, f -> cardinalGroupCount);
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
void bookPair (uint    levelCount,
               uint    groupIndex,
               uint    levelIndex,
               uint   *row,
               uint   *level,
               Factor *f) {
  uint i;
  level[levelIndex] ++;
  if (levelIndex < groupIndex) {
    levelIndex ++;
    level[levelIndex] ++;
    while (level[levelIndex] < level[levelIndex-1]) {
      level[levelIndex] ++;
    }
    bookPair(levelCount, groupIndex, levelIndex, row, level, f);
    level[levelIndex] = 0;
    levelIndex --;
    if ((*row) < ((uint*) (f -> cardinalGroupSize))[groupIndex]) {
      if (level[levelIndex] < levelCount - (groupIndex - levelIndex)) {
        bookPair(levelCount, groupIndex, levelIndex, row, level, f);
      }
    }
  }
  else {
    (*row)++;
    (f -> cardinalGroupBinary)[groupIndex][*row] = 0;
    for (i=1; i <=groupIndex; i++) {
      (f -> cardinalGroupBinary)[groupIndex][*row] += upower(2, level[i] - 1);
    }
    if ( (levelCount > 2) && (level[levelIndex] < levelCount)) {
      bookPair(levelCount, groupIndex, levelIndex, row, level, f);
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
      reduceFraction(& newMultiplier, & newDivisor);
      reduceFraction(& total, & newDivisor);
      if (newMultiplier > (UINT_MAX / total)) {
        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
        RF_nativeError("\nRF-SRC:  Arithmetic Overflow Encountered in nChooseK(n, k). ");
        RF_nativeError("\nRF-SRC:  Incoming parameters are (%10d, %10d). ", n, r);
        RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
        RF_nativeExit();
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
  uint mwcpWordIdent = (level >> (3 + ulog2(sizeof(uint)))) + ((level & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
  uint binaryWord = upower(2, level - ((mwcpWordIdent - 1) * MAX_EXACT_LEVEL) - 1 );
  daughterFlag = RIGHT;
  if (binaryWord & mwcp[mwcpWordIdent]) {
    daughterFlag = LEFT;
  }
  return daughterFlag;
}
