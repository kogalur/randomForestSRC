#ifndef RF_NATIVE_UTIL_H
#define RF_NATIVE_UTIL_H
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
#endif
