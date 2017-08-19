JNIEXPORT jobject JNICALL Java_com_kogalur_randomforest_Native_grow(JNIEnv      *env,
                                                                    jobject      obj,
                                                                    jint         traceFlag,
                                                                    jint         seedPtr,
                                                                    jint         optLow,
                                                                    jint         optHigh,
                                                                    jint         splitRule,
                                                                    jint         nsplit,
                                                                    jint         mtry,
                                                                    jint         ytry,
                                                                    jint         nodeSize,
                                                                    jint         nodeDepth,
                                                                    jint         crWeightSize,
                                                                    jdoubleArray crWeight,
                                                                    jint         ntree,
                                                                    jint         observationSize,
                                                                    jint         ySize,
                                                                    jcharArray   rType,
                                                                    jintArray    rLevels,
                                                                    jobject      rData,
                                                                    jint         xSize,
                                                                    jcharArray   xType,
                                                                    jintArray    xLevels,
                                                                    jint         bootstrapSize,
                                                                    jobject      bootstrap,
                                                                    jdoubleArray caseWeight,
                                                                    jdoubleArray xSplitStatWt,
                                                                    jdoubleArray yWeight,
                                                                    jdoubleArray xWeight,
                                                                    jobject      xData,
                                                                    jint         timeInterestSize,
                                                                    jdoubleArray timeInterest,
                                                                    jint         nImpute,
                                                                    jint         numThreads) {
  setUserTraceFlag((uint) traceFlag);
  setNativeGlobalEnv(env, obj);
  int seedValue           = (int)  seedPtr;
  RF_opt                  = (uint) optLow;
  RF_optHigh              = (uint) optHigh;
  RF_splitRule            = (uint) splitRule;
  RF_nsplit               = (uint) nsplit;
  RF_mtry                 = (uint) mtry;
  RF_ytry                 = (uint) ytry;
  RF_nodeSize             = (uint) nodeSize;
  RF_nodeDepth            = (int)  nodeDepth;
  RF_crWeightSize         = (uint) crWeightSize;
  JavaVM *jvm;
  int status = (*env) -> GetJavaVM(env, &jvm);
  if(status != 0) {
    printf("\n TEST FAIL \n");    
  }  
  JNIEnv *envNew;
  (*jvm)->AttachCurrentThread(jvm, (void **)&envNew, NULL);
  RF_crWeight             = (double *) copy1DObject(crWeight, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize);
  RF_ntree                = (uint) ntree;
  RF_observationSize      = (uint) observationSize;
  RF_ySize                = (uint) ySize;
  RF_rType                = (char *) copy1DObject(rType, NATIVE_TYPE_CHARACTER, &RF_jni1DInfoListSize);
  RF_rLevels              = (int*) copy1DObject(rLevels, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize);
  RF_responseIn           = (double **) copy2DObject(rData, NATIVE_TYPE_NUMERIC, &RF_jni2DInfoListSize);
  RF_xSize                = (uint) xSize;
  RF_xType                = (char *) copy1DObject(xType, NATIVE_TYPE_CHARACTER, &RF_jni1DInfoListSize);
  RF_xLevels              = (int *) copy1DObject(xLevels, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize);
  RF_bootstrapSize        = (uint) bootstrapSize;
  RF_bootstrapIn        = (uint **) copy2DObject(bootstrap, NATIVE_TYPE_INTEGER, &RF_jni2DInfoListSize);
  RF_caseWeight           = (double *) copy1DObject(caseWeight, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize);
  RF_xSplitStatWt          = (double *) copy1DObject(xSplitStatWt, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize);
  RF_yWeight              = (double *) copy1DObject(yWeight, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize);
  RF_xWeight              = (double *) copy1DObject(xWeight, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize);
  RF_observationIn        = (double **) copy2DObject(xData, NATIVE_TYPE_NUMERIC, &RF_jni2DInfoListSize);
  RF_timeInterestSize     = (uint) timeInterestSize;
  RF_timeInterest         = (double *) copy1DObject(timeInterest, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize);
  RF_nImpute              = (uint) nImpute;
  RF_numThreads           = (int)  numThreads;
  processDefaultGrow();
  stackAuxiliaryInfoList();
  rfsrc(RF_GROW, seedValue);
  put_jniEnsembleInfoList(RF_nativeIndex);
  free_jvvector(RF_jniEnsembleInfoList, 0, 1 << 6, NRUTIL_JEN_PTR);
  free_jni1DList(RF_jni1DInfoListSize);
  free_jni2DList(RF_jni2DInfoListSize);
  free_jvvector(RF_jni1DInfoList, 0, 1 << 6, NRUTIL_J1D_PTR);
  free_jvvector(RF_jni2DInfoList, 0, 1 << 6, NRUTIL_J2D_PTR);
  unstackAuxiliaryInfoAndList();
  if (RF_nativeIndex != RF_stackCount) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Stack imbalance in PROTECT/UNPROTECT:  %10d + 1 versus %10d  ", RF_nativeIndex, RF_stackCount);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  return RF_java_arrlst_obj;
}
void exit2J() {
  jstring jbuffer;
  if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
    jthrowable flag =  (*RF_java_env) -> ExceptionOccurred(RF_java_env);
    (*RF_java_env) -> ExceptionDescribe(RF_java_env);
    (*RF_java_env) -> ThrowNew(RF_java_env, RF_java_except_cls, "Thrown from Native Code.");
    (*RF_java_env) -> ExceptionClear(RF_java_env);
  }
  jbuffer = (*RF_java_env) -> NewStringUTF(RF_java_env, "\nRF-SRC:  The application will now exit.\n");
  (*RF_java_env) -> CallStaticVoidMethod(RF_java_env, RF_java_cls, RF_java_mid_logError, jbuffer);
  (*RF_java_env) -> DeleteLocalRef(RF_java_env, jbuffer);
  (*RF_java_env) -> CallStaticVoidMethod(RF_java_env, RF_java_cls, RF_java_mid_logExit);
}
void errorJ(char *format, ...) {
  char *buffer;
  va_list aptr;
  buffer = (char *) malloc(sizeof(char) * 1023);
  va_start(aptr, format);
  vsprintf(buffer, format, aptr);
  va_end(aptr);
  jstring jbuffer = (*RF_java_env) -> NewStringUTF(RF_java_env, (const char*) buffer);
  (*RF_java_env) -> CallStaticVoidMethod(RF_java_env, RF_java_cls, RF_java_mid_logError, jbuffer);
  (*RF_java_env) -> DeleteLocalRef(RF_java_env, jbuffer);
  free((char *) buffer);
}
void printJ(char *format, ...) {
  char *buffer;
  va_list aptr;
  buffer = (char *) malloc(sizeof(char) * 1023);
  va_start(aptr, format);
  vsprintf(buffer, format, aptr);
  va_end(aptr);
  jstring jbuffer = (*RF_java_env) -> NewStringUTF(RF_java_env, (const char*) buffer);
  (*RF_java_env) -> CallStaticVoidMethod(RF_java_env, RF_java_cls, RF_java_mid_log, jbuffer);
  (*RF_java_env) -> DeleteLocalRef(RF_java_env, jbuffer);  
  free((char *) buffer);
}
void setNativeGlobalEnv(JNIEnv *env, jobject obj) {
  RF_java_env = env;
  RF_java_obj = obj;
  RF_java_cls = (*RF_java_env) -> GetObjectClass(RF_java_env, RF_java_obj);
  if (RF_java_cls == NULL) {
    printf("\nRF-SRC:  Unable to access calling class com/kogalur/randomforest/Native.\n");
    printf("\nRF-SRC:  The application will now exit.\n");
    exit(1);
  }
  RF_java_except_cls = (*RF_java_env) -> FindClass(RF_java_env, "java/lang/RuntimeException");
  if (RF_java_except_cls == NULL) {
    printf("\nRF-SRC:  Unable to access exception class java/lang/RuntimeException.\n");
    printf("\nRF-SRC:  The application will now exit.\n");
    exit(1);
  }
  RF_java_mid_log = (*RF_java_env) -> GetStaticMethodID(RF_java_env, RF_java_cls, "log", "(Ljava/lang/String;)V");
  if (RF_java_mid_log == NULL) {
    printf("\nRF-SRC:  Unable to access static method com/kogalur/randomforest/Native::log().\n");
    printf("\nRF-SRC:  The application will now exit.\n");
    exit(1);
  }
  RF_java_mid_logError = (*RF_java_env) -> GetStaticMethodID(RF_java_env, RF_java_cls, "logError", "(Ljava/lang/String;)V");
  if (RF_java_mid_logError == NULL) {
    printf("\nRF-SRC:  Unable to access static method com/kogalur/randomforest/Native::log().\n");
    printf("\nRF-SRC:  The application will now exit.\n");
    exit(1);
  }
  RF_java_mid_logExit = (*RF_java_env) -> GetStaticMethodID(RF_java_env, RF_java_cls, "logExit", "()V");
  if (RF_java_mid_logExit == NULL) {
    printf("\nRF-SRC:  Unable to access static method com/kogalur/randomforest/Native::logExit().\n");
    printf("\nRF-SRC:  The application will now exit.\n");
    exit(1);
  }
  RF_java_arrlst_cls = (*RF_java_env) -> FindClass(RF_java_env, "java/util/ArrayList");
  if (RF_java_arrlst_cls == NULL) {
    RF_nativeError("\nRF-SRC:  Unable to access class for java/util/ArrayList.\n");
    RF_nativeError("\nRF-SRC:  The application will now exit.\n");
    RF_nativeExit();
  }
  RF_java_arrlst_constr = (*RF_java_env) -> GetMethodID(RF_java_env, RF_java_arrlst_cls, "<init>", "(I)V");
  if (RF_java_arrlst_constr == NULL) {
    RF_nativeError("\nRF-SRC:  Unable to access constructor for class java/util/ArrayList.\n");
    RF_nativeError("\nRF-SRC:  The application will now exit.\n");
    RF_nativeExit();
  }
  RF_java_arrlst_add    = (*RF_java_env) -> GetMethodID(RF_java_env, RF_java_arrlst_cls, "add", "(Ljava/lang/Object;)Z");
  if (RF_java_arrlst_add == NULL) {
    RF_nativeError("\nRF-SRC:  Unable to access method java/util/ArrayList::add().\n");
    RF_nativeError("\nRF-SRC:  The application will now exit.\n");
    RF_nativeExit();
  }
  RF_java_arrlst_obj = (*RF_java_env) -> NewObject(RF_java_env, RF_java_arrlst_cls, RF_java_arrlst_constr, 1 << 6);
  if (RF_java_arrlst_obj == NULL) {
    RF_nativeError("\nRF-SRC:  Unable to instantiate object java/util/ArrayList.\n");
    RF_nativeError("\nRF-SRC:  The application will now exit.\n");
    RF_nativeExit();
  }
  RF_java_ens_cls = (*RF_java_env) -> FindClass(RF_java_env, "com/kogalur/randomforest/Ensemble");
  if (RF_java_ens_cls == NULL) {
    RF_nativeError("\nRF-SRC:  Unable to access class com/kogalur/randomforest/Ensemble.\n");
    RF_nativeError("\nRF-SRC:  The application will now exit.\n");
    exit(1);
  }
  RF_java_ens_mid = (*RF_java_env) -> GetMethodID(RF_java_env, RF_java_ens_cls, "Ensemble", "(Ljava/lang/String;CIJZI[ILjava/lang/Object;)V");
  if (RF_java_ens_mid == NULL) {
    RF_nativeError("\nRF-SRC:  Unable to access constructor for class com/kogalur/randomforest/Ensemble.\n");
    RF_nativeError("\nRF-SRC:  The application will now exit.\n");
    RF_nativeExit();
  }
  RF_jni1DInfoList = (JNI1DInfo **) jvvector(0, 1 << 6, NRUTIL_J1D_PTR);
  RF_jni2DInfoList = (JNI2DInfo **) jvvector(0, 1 << 6, NRUTIL_J2D_PTR);
  RF_jniEnsembleInfoList = (JNIEnsembleInfo **) jvvector(0, 1 << 6, NRUTIL_JEN_PTR);
  RF_jni1DInfoListSize = 0;
  RF_jni2DInfoListSize = 0;
  RF_jniEnsembleInfoListSize = 0;
  RF_nativeIndex = RF_stackCount = 0;
}
void *jvvector(unsigned long long nl, unsigned long long nh, enum alloc_type type) {
  void *v;
  v = NULL;  
  switch(type){
  case NRUTIL_J1D_PTR:
    v = ((JNI1DInfo **) (gvector(nl, nh, sizeof(JNI1DInfo *)) -nl+NR_END));
    break;
  case NRUTIL_J2D_PTR:
    v = ((JNI2DInfo **) (gvector(nl, nh, sizeof(JNI2DInfo *)) -nl+NR_END));
    break;
  case NRUTIL_JEN_PTR:
    v = ((JNIEnsembleInfo **) (gvector(nl, nh, sizeof(JNIEnsembleInfo *)) -nl+NR_END));
    break;
  default:
    v = NULL;
    nrerror("\n  Illegal case in jvvector().");
    break;
  }
  return v;
}
void free_jvvector(void *v, unsigned long long nl, unsigned long long nh, enum alloc_type type) {
  switch(type){
  case NRUTIL_J1D_PTR:
    free_gvector((JNI1DInfo *) (v+nl-NR_END), nl, nh, sizeof(JNI1DInfo *));
    break;
  case NRUTIL_J2D_PTR:
    free_gvector((JNI2DInfo *) (v+nl-NR_END), nl, nh, sizeof(JNI2DInfo *));
    break;
  case NRUTIL_JEN_PTR:
    free_gvector((JNIEnsembleInfo *) (v+nl-NR_END), nl, nh, sizeof(JNIEnsembleInfo *));
    break;
  default:
    nrerror("\n  Illegal case in free_jvvector().");
    break;
  }
}
jboolean *jbvector(unsigned long long nl, unsigned long long nh) {
  return ((jboolean *) gvector(nl, nh, sizeof(jboolean)) -nl+NR_END);
}
void free_jbvector(jboolean *v, unsigned long long nl, unsigned long long nh) {
  free_gvector(v+nl-NR_END, nl, nh, sizeof(jboolean));
}
jarray *jvector(unsigned long long nl, unsigned long long nh) {
  return ((jarray *) gvector(nl, nh, sizeof(jarray)) -nl+NR_END);
}
void free_jvector(jarray *v, unsigned long long nl, unsigned long long nh) {
  free_gvector(v+nl-NR_END, nl, nh, sizeof(jarray));
}
void *copy1DObject(jarray arr, char type, uint *index) {
  jdouble *dbuffer;
  jint    *ibuffer;
  jchar   *cbuffer;
  jboolean isCopy;
  uint     len;
  uint     i;
  double  *dcopy;
  int     *icopy;
  char    *ccopy;
  void    *copy;
  JNI1DInfo *incomingInfo;
  copy = NULL;
  if (! (*RF_java_env) -> IsSameObject(RF_java_env, arr, NULL)) {
    RF_jni1DInfoList[*index] = incomingInfo = (JNI1DInfo*) gblock((size_t) sizeof(JNI1DInfo));
    len = (*RF_java_env) -> GetArrayLength(RF_java_env, arr);
    if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
      RF_nativeExit();
    }
    (incomingInfo -> len) = len;
    (incomingInfo -> type) = type;
    switch (type) {
    case NATIVE_TYPE_NUMERIC:
      (incomingInfo -> array) = dcopy = dvector(1, len);
      dbuffer =  (*RF_java_env) -> GetDoubleArrayElements(RF_java_env, arr, &isCopy);
      if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
        RF_nativeExit();
      }
      for (i = 0; i < len; i++) {
        dcopy[i+1] = (double) dbuffer[i];
      }
      if (isCopy == JNI_TRUE) {
        (*RF_java_env) -> ReleaseDoubleArrayElements(RF_java_env, arr, (jdouble*) dbuffer, JNI_ABORT);
        if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
          RF_nativeExit();
        }
      }
      copy = dcopy;
      (*RF_java_env) -> DeleteLocalRef(RF_java_env, arr);
      break;
    case NATIVE_TYPE_INTEGER:
      (incomingInfo -> array) = icopy = ivector(1, len);
      ibuffer =  (*RF_java_env) -> GetIntArrayElements(RF_java_env, arr, &isCopy);
      if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
        RF_nativeExit();
      }
      for (i = 0; i < len; i++) {
        icopy[i+1] = (int) ibuffer[i];
      }
      if (isCopy == JNI_TRUE) {
        (*RF_java_env) -> ReleaseIntArrayElements(RF_java_env, arr, (jint*) ibuffer, JNI_ABORT);
        if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
          RF_nativeExit();
        }
      }
      copy = icopy;
      (*RF_java_env) -> DeleteLocalRef(RF_java_env, arr);
      break;
    case NATIVE_TYPE_CHARACTER:
      (incomingInfo -> array) = ccopy = cvector(1, len);
      cbuffer =  (*RF_java_env) -> GetCharArrayElements(RF_java_env, arr, &isCopy);
      if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
        RF_nativeExit();
      }
      for (i = 0; i < len; i++) {
        ccopy[i+1] = (char) cbuffer[i];
      }
      if (isCopy == JNI_TRUE) {
        (*RF_java_env) -> ReleaseCharArrayElements(RF_java_env, arr, (jchar*) cbuffer, JNI_ABORT);
        if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
          RF_nativeExit();
        }
      }
      copy = ccopy;
      (*RF_java_env) -> DeleteLocalRef(RF_java_env, arr);
      break;
    }
    (*index) ++;
  }
  return copy;
}
void free_jni1DList(uint size) {
  JNI1DInfo *incomingInfo;
  for (uint i = 0; i < size; i++) {
    incomingInfo = RF_jni1DInfoList[i];
    switch (incomingInfo -> type) {
    case NATIVE_TYPE_NUMERIC:
      free_dvector((double *) (incomingInfo -> array), 1, incomingInfo -> len);
      break;
    case NATIVE_TYPE_INTEGER:
      free_ivector((int *) (incomingInfo -> array), 1, incomingInfo -> len);    
      break;
    case NATIVE_TYPE_CHARACTER:
      free_cvector((char *) (incomingInfo -> array), 1, incomingInfo -> len);    
      break;
    }
    free_gblock(incomingInfo, (size_t) sizeof(JNI1DInfo));    
  }
}
void *copy2DObject(jobject obj, char type, uint *index) {
  JNI2DInfo *incomingInfo;
  uint outLen;
  jboolean isCopy;
  if (! (*RF_java_env) -> IsSameObject(RF_java_env, obj, NULL)) {
    RF_jni2DInfoList[*index] = incomingInfo = (JNI2DInfo*) gblock((size_t) sizeof(JNI2DInfo));
    outLen = (*RF_java_env) -> GetArrayLength(RF_java_env, obj);
    if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
      RF_nativeExit();
    }
    switch (type) {
    case NATIVE_TYPE_NUMERIC:
      (incomingInfo -> outerPtr) = (double **) new_vvector(1, outLen, NRUTIL_DPTR);
      break;
    case NATIVE_TYPE_INTEGER:
      (incomingInfo -> outerPtr) = (uint **) new_vvector(1, outLen, NRUTIL_UPTR);
      break;
    }
    (incomingInfo -> outLen)   = outLen;
    (incomingInfo -> type)     = type;
    (incomingInfo -> isCopy)   = jbvector(1, outLen);
    (incomingInfo -> innerPtr) = jvector(1, outLen);
    (incomingInfo -> innLen)   = uivector(1, outLen);
    for(uint i = 0; i < outLen; ++i) {
      (incomingInfo -> innerPtr)[i+1] = (*RF_java_env) -> GetObjectArrayElement(RF_java_env, obj, i);
      if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
        RF_nativeExit();
      }
      (incomingInfo -> innLen)[i+1] = (*RF_java_env) -> GetArrayLength(RF_java_env, (incomingInfo -> innerPtr)[i+1]);
      if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
        RF_nativeExit();
      }
      switch (type) {
      case NATIVE_TYPE_NUMERIC:
        ((double **) (incomingInfo -> outerPtr))[i+1] = (double *) (((*RF_java_env) ->GetDoubleArrayElements(RF_java_env, (incomingInfo -> innerPtr)[i+1], &isCopy)) - 1);
        if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
          RF_nativeExit();
        }
        break;
      case NATIVE_TYPE_INTEGER:
        ((uint **) (incomingInfo -> outerPtr))[i+1] = (uint *) (((*RF_java_env) ->GetIntArrayElements(RF_java_env, (incomingInfo -> innerPtr)[i+1], &isCopy)) - 1 );
        if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
          RF_nativeExit();
        }
        break;
      }
      (incomingInfo -> isCopy)[i+1] = isCopy;
    }
    (*index) ++;
  }
  return (incomingInfo -> outerPtr);
}
void free_jni2DList(uint size) {
  JNI2DInfo *incomingInfo;
  for (uint i = 0; i < size; i++) {
    incomingInfo = RF_jni2DInfoList[i];
    for(uint j = 0; j < incomingInfo -> outLen; j++) {    
      switch (incomingInfo -> type) {
      case NATIVE_TYPE_NUMERIC:
        if ((incomingInfo -> isCopy)[j+1] == JNI_TRUE) {
          (*RF_java_env) -> ReleaseDoubleArrayElements(RF_java_env, (incomingInfo -> innerPtr)[j+1],  ((jdouble **) (incomingInfo -> outerPtr))[j+1] + 1, 0);
          if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
            RF_nativeExit();
          }
        }
        (*RF_java_env) -> DeleteLocalRef(RF_java_env, (incomingInfo -> innerPtr)[j+1]);
        break;
      case NATIVE_TYPE_INTEGER:
        if ((incomingInfo -> isCopy)[j+1] == JNI_TRUE) {
          (*RF_java_env) -> ReleaseIntArrayElements(RF_java_env, (incomingInfo -> innerPtr)[j+1],  ((jint **) (incomingInfo -> outerPtr))[j+1] + 1, 0);
          if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
            RF_nativeExit();
          }
        }
        (*RF_java_env) -> DeleteLocalRef(RF_java_env, (incomingInfo -> innerPtr)[j+1]);
        break;
      }
    }
    switch (incomingInfo -> type) {
    case NATIVE_TYPE_NUMERIC:
      free_new_vvector((double **) (incomingInfo -> outerPtr), 1, incomingInfo -> outLen, NRUTIL_DPTR);
      break;
    case NATIVE_TYPE_INTEGER:
      free_new_vvector((uint **) (incomingInfo -> outerPtr), 1, incomingInfo -> outLen, NRUTIL_UPTR);      
      break;
    }
    free_jbvector(incomingInfo -> isCopy, 1, (incomingInfo -> outLen));
    free_jvector(incomingInfo -> innerPtr, 1, (incomingInfo -> outLen));
    free_uivector(incomingInfo -> innLen, 1, (incomingInfo -> outLen));
    free_gblock(incomingInfo, (size_t) sizeof(JNI2DInfo));    
  }
}
void initProtect(uint stackCount) {
}
void *stackAndProtect(uint  *index,
                      char   type,
                      uint   identity,
                      ulong  size,
                      char **sexpString,
                      void  *auxiliaryPtr,
                      uint   auxiliaryDimSize,
                      ...) {
  jboolean isCopy;
  jmethodID mid;
  jarray  thisArray;
  void   *thisArrayPtr;
  JNIEnsembleInfo *ensembleInfo;
  if (((*index) >> 6) > 0) {
          RF_nativeError("\nRF-SRC:  *** ERROR *** ");
          RF_nativeError("\nRF-SRC:  S.E.X.P. vector list limit exceeded:  %20d", *index);
          RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
          RF_nativeExit();
  }
  if (sizeof(ulong) > sizeof(uint)) {
    if (size > UINT_MAX) {
      if (TRUE) {
        RF_nativePrint("\nRF-SRC:  *** WARNING *** ");
        RF_nativePrint("\nRF-SRC:  S.E.X.P. vector element length exceeds 32-bits:  %20lu", size);
        RF_nativePrint("\nRF-SRC:  S.E.X.P. ALLOC:  %s ", sexpString[identity]);
        RF_nativePrint("\nRF-SRC:  Please Reduce Dimensionality If Possible.");
      }
    }
  }
  va_list list;
  va_start(list, auxiliaryDimSize);
  uint *auxiliaryDim = uivector(1, auxiliaryDimSize);
  for (int i = 1; i <= auxiliaryDimSize; i++) {
    auxiliaryDim[i] = va_arg(list, unsigned int);
  }
  va_end(list);
  RF_jniEnsembleInfoList[*index] = ensembleInfo = (JNIEnsembleInfo*) gblock((size_t) sizeof(JNIEnsembleInfo));
  ensembleInfo -> type     = type;
  ensembleInfo -> size     = size;
  ensembleInfo -> identity = identity;
  (*index) ++;
  switch (type) {
  case NATIVE_TYPE_NUMERIC:
    thisArray = (jdoubleArray) (*RF_java_env) -> NewDoubleArray(RF_java_env, size);
    if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
      RF_nativeExit();
    }
    thisArrayPtr = (jdoubleArray *) (*RF_java_env) -> GetDoubleArrayElements(RF_java_env, thisArray, &isCopy);
    if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
      RF_nativeExit();
    }
    break;
  case NATIVE_TYPE_INTEGER:
    thisArray = (jintArray) (*RF_java_env) -> NewIntArray(RF_java_env, size);
    if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
      RF_nativeExit();
    }
    thisArrayPtr = (jintArray *) (*RF_java_env) -> GetIntArrayElements(RF_java_env, thisArray, &isCopy);
    if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
      RF_nativeExit();
    }
    break;
  }
  ensembleInfo -> isCopy = isCopy;
  (ensembleInfo -> array) = thisArray;
  (ensembleInfo -> arrayPtr) = thisArrayPtr;
  allocateAuxiliaryInfo(type,
                        identity,
                        size,
                        ensembleInfo -> arrayPtr,
                        auxiliaryPtr,
                        auxiliaryDimSize,
                        auxiliaryDim);
  return (ensembleInfo -> arrayPtr);
}
void put_jniEnsembleInfoList(uint size) {
  jmethodID mid;
  JNIEnsembleInfo *ensembleInfo;
  SNPAuxiliaryInfo *auxInfoPtr;
  jstring name;
  jboolean result;
  for (uint i = 0; i < size; i++) {
    ensembleInfo = RF_jniEnsembleInfoList[i];
    auxInfoPtr = RF_snpAuxiliaryInfoList[ensembleInfo -> identity];
    name = (*RF_java_env) -> NewStringUTF(RF_java_env, RF_sexpString[ensembleInfo -> identity]);
    if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
      RF_nativeExit();
    }
    switch (ensembleInfo -> type) {
    case NATIVE_TYPE_NUMERIC:
      if (ensembleInfo -> isCopy == JNI_TRUE) {
        (*RF_java_env) -> ReleaseDoubleArrayElements(RF_java_env,
                                                     ensembleInfo -> array,
                                                     (jdouble *) (ensembleInfo -> arrayPtr), 0);
        if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
          RF_nativeExit();
        }
      }
      break;
    case NATIVE_TYPE_INTEGER:
      if (ensembleInfo -> isCopy == JNI_TRUE) {
        (*RF_java_env) -> ReleaseIntArrayElements(RF_java_env,
                                                  ensembleInfo -> array,
                                                  (jint *) (ensembleInfo -> arrayPtr), 0);
        if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
          RF_nativeExit();
        }
      }
      break;
    }
    jobject ensembleObj = (*RF_java_env) -> NewObject(RF_java_env, RF_java_ens_cls, RF_java_ens_mid,
                                                      (jstring) name,
                                                      (jchar) ensembleInfo -> type,
                                                      (jint) ensembleInfo -> identity,
                                                      (jlong) ensembleInfo -> size,
                                                      (jboolean) (((auxInfoPtr -> auxiliaryPtr) == NULL) && (auxInfoPtr -> dimSize > 1)) ? JNI_TRUE : JNI_FALSE,
                                                      (jint) (auxInfoPtr -> dimSize),
                                                      (jintArray) (auxInfoPtr -> dim),
                                                      (jobject) (ensembleInfo -> array));
    if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
      RF_nativeExit();
    }
    (*RF_java_env) -> CallObjectMethod(RF_java_env, RF_java_arrlst_obj, RF_java_arrlst_add, ensembleObj);
    if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
      RF_nativeExit();
    }
    if (FALSE) {
      RF_nativePrint("\nRF-SRC:  Unable to add ensemble object to java/util/ArrayList.\n");
      RF_nativePrint("\nRF-SRC:  The application will now exit.\n");
      exit(1);
    }
    (*RF_java_env) -> DeleteLocalRef(RF_java_env, ensembleInfo -> array);      
    free_gblock(ensembleInfo, (size_t) sizeof(JNIEnsembleInfo));  
  }
}
void setUserTraceFlag (uint traceFlag) {
  RF_userTraceFlag = traceFlag;
}
uint getUserTraceFlag () { return RF_userTraceFlag; }
