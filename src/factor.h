#ifndef  RF_FACTOR_H
#define  RF_FACTOR_H
typedef struct factor Factor;
struct factor {
  unsigned int r; 
  unsigned int cardinalGroupCount; 
  void *complementaryPairCount;
  void *cardinalGroupSize; 
  unsigned int **cardinalGroupBinary;
  unsigned int mwcpSize;
};
#endif
