void exit2J();
void errorJ(char *format, ...);
void printJ(char *format, ...);
void throwRuntimeException(char *message);
enum alloc_jtype{
  NRUTIL_J1D_PTR,
  NRUTIL_J2D_PTR,
  NRUTIL_JEN_PTR
};
typedef struct jni1DInfo JNI1DInfo;
struct jni1DInfo {
  jarray arrayJVM;
  void *arrayJNI;
  void *array;
  uint len;
  char type;
  char actual;
  jboolean isCopy;
};
typedef struct jni2DInfo JNI2DInfo;
struct jni2DInfo {
  void *outerPtr;
  jarray *innerPtr;
  jboolean *isCopy;
  uint outLen;
  uint *innLen;
  char type;
};
typedef struct jniEnsembleInfo JNIEnsembleInfo;
struct jniEnsembleInfo {
  jarray array;
  void *arrayPtr;
  jint dimSize;
  jintArray dim;
  int *dimPtr;
  ulong size;
  char type;
  uint identity;
  jboolean isCopy;
  jboolean isCopyDim;
};
void setNativeGlobalEnv(JNIEnv *env, jobject obj);
void *jvvector(unsigned long long nl, unsigned long long nh, enum alloc_jtype type);
void free_jvvector(void *v, unsigned long long nl, unsigned long long nh, enum alloc_jtype type);
jboolean *jbvector(unsigned long long nl, unsigned long long nh);
void free_jbvector(jboolean *v, unsigned long long nl, unsigned long long nh);
jarray *jvector(unsigned long long nl, unsigned long long nh);
void free_jvector(jarray *v, unsigned long long nl, unsigned long long nh);
void *copy1DObject(jarray arr, char type, uint *index, char actual);
void *copy2DObject(jarray arr, char type, uint *index);
void free_jni1DList(uint size);
void free_jni2DList(uint size);
void put_jniEnsembleInfoList(uint size);
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
void setUserTraceFlag (uint traceFlag);
uint getUserTraceFlag ();
void populatePartialObject(jobject obj);
void populateHyperZeroObject(jobject obj);
void populateHyperOneObject(jobject obj);
