void exit2J();
void errorJ(char *format, ...);
void printJ(char *format, ...);
void throwRuntimeException(char *message);
enum alloc_ntype{
  NRUTIL_J1D_PTR,
  NRUTIL_J2D_PTR,
  NRUTIL_JEN_PTR
};
typedef struct jni1DInfo NAT1DInfo;
struct jni1DInfo {
  jarray arrayJVM;
  void *arrayJNI;
  void *array;
  uint len;
  char type;
  char actual;
  jboolean isCopy;
};
typedef struct jni2DInfo NAT2DInfo;
struct jni2DInfo {
  void *outerPtr;
  jarray *innerPtr;
  uint outLen;
  uint *innLen;
  char type;
  jboolean *isCopy;
};
typedef struct nativeEnsembleInfo NativeEnsembleInfo;
struct nativeEnsembleInfo {
  char *identity;
  jarray array;
  void *arrayPtr;
  jint dimSize;
  jintArray dim;
  int *dimPtr;
  ulong size;
  char type;
  jint slot;
  jboolean isCopy;
  jboolean isCopyDim;
};
void setNativeGlobalEnv(JNIEnv *env, jobject obj);
void *nvvector(unsigned long long nl, unsigned long long nh, enum alloc_ntype type);
void free_nvvector(void *v, unsigned long long nl, unsigned long long nh, enum alloc_ntype type);
jboolean *jbvector(unsigned long long nl, unsigned long long nh);
void free_jbvector(jboolean *v, unsigned long long nl, unsigned long long nh);
jarray *jvector(unsigned long long nl, unsigned long long nh);
void free_jvector(jarray *v, unsigned long long nl, unsigned long long nh);
void *copy1DObject(jarray arr, char type, uint *index, char actual);
void *copy2DObject(jarray arr, char type, uint *index);
void free_jni1DList(uint size);
void free_jni2DList(uint size);
void initProtect(uint stackCount);
void *stackAndProtect(uint  *index,
                      char   type,
                      uint   identity,
                      ulong  size,
                      double value,
                      char  *sexpString,
                      void  *auxiliaryArrayPtr,
                      uint   auxiliaryDimSize,
                      ...);
void put_nativeEnsembleInfoList(uint size);
void setUserTraceFlag (uint traceFlag);
uint getUserTraceFlag ();
void populateHyperZeroObject(jobject obj);
void populateHyperOneObject(jobject obj);
void populatePartialObject(jobject obj);
void populateLotObject(jobject obj);
void populateBaseLearnObject(jobject obj);
void populateYInfoObject(jobject obj);
void populateXInfoObject(jobject obj);
void populateSampleObject(jobject obj);
void populateVtryObject(jobject obj);
void populateSeedObject(jobject obj);
void populateQuantileObject(jobject obj);
void populateYTargetObject(jobject obj);
