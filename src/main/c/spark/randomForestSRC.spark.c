JNIEXPORT jobject JNICALL Java_com_kogalur_randomforest_Native_grow(JNIEnv      *env,
                                                                    jobject      obj,
                                                                    jint         traceFlag,
                                                                    jint         seedDynamic,
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
                                                                    jdoubleArray xSplitStatWeight,
                                                                    jdoubleArray yWeight,
                                                                    jdoubleArray xWeight,
                                                                    jobject      xData,
                                                                    jint         timeInterestSize,
                                                                    jdoubleArray timeInterest,
                                                                    jint         nImpute,
                                                                    jint         numThreads) {
  setUserTraceFlag((uint) traceFlag);
  setNativeGlobalEnv(env, obj);
  int seedValue           = (int)  seedDynamic;
  RF_opt                  = (uint) optLow;
  RF_optHigh              = (uint) optHigh;
  RF_splitRule            = (uint) splitRule;
  RF_nsplit               = (uint) nsplit;
  RF_mtry                 = (uint) mtry;
  RF_ytry                 = (uint) ytry;
  RF_nodeSize             = (uint) nodeSize;
  RF_nodeDepth            = (int)  nodeDepth;
  RF_crWeightSize         = (uint) crWeightSize;
  RF_crWeight             = (double *) copy1DObject(crWeight, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_ntree                = (uint) ntree;
  RF_observationSize      = (uint) observationSize;
  RF_ySize                = (uint) ySize;
  RF_rType                = (char *) copy1DObject(rType, NATIVE_TYPE_CHARACTER, &RF_jni1DInfoListSize, TRUE);
  RF_rLevels              = (uint*) copy1DObject(rLevels, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
  RF_responseIn           = (double **) copy2DObject(rData, NATIVE_TYPE_NUMERIC, &RF_jni2DInfoListSize);
  RF_xSize                = (uint) xSize;
  RF_xType                = (char *) copy1DObject(xType, NATIVE_TYPE_CHARACTER, &RF_jni1DInfoListSize, TRUE);
  RF_xLevels              = (uint *) copy1DObject(xLevels, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
  RF_bootstrapSize        = (uint) bootstrapSize;
  RF_bootstrapIn          = (uint **) copy2DObject(bootstrap, NATIVE_TYPE_INTEGER, &RF_jni2DInfoListSize);
  RF_caseWeight           = (double *) copy1DObject(caseWeight, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_xSplitStatWeight     = (double *) copy1DObject(xSplitStatWeight, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_yWeight              = (double *) copy1DObject(yWeight, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_xWeight              = (double *) copy1DObject(xWeight, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_observationIn        = (double **) copy2DObject(xData, NATIVE_TYPE_NUMERIC, &RF_jni2DInfoListSize);
  RF_timeInterestSize     = (uint) timeInterestSize;
  RF_timeInterest         = (double *) copy1DObject(timeInterest, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_nImpute              = (uint) nImpute;
  RF_numThreads           = (int)  numThreads;
  processDefaultGrow();
  rfsrc(RF_GROW, seedValue);
  put_jniEnsembleInfoList(RF_nativeIndex);
  free_jvvector(RF_jniEnsembleInfoList, 0, 1 << 6, NRUTIL_JEN_PTR);
  free_jni1DList(RF_jni1DInfoListSize);
  free_jni2DList(RF_jni2DInfoListSize);
  free_jvvector(RF_jni1DInfoList, 0, 1 << 6, NRUTIL_J1D_PTR);
  free_jvvector(RF_jni2DInfoList, 0, 1 << 6, NRUTIL_J2D_PTR);
  memoryCheck();
  return RF_java_hshmap_obj;
}
JNIEXPORT jobject JNICALL Java_com_kogalur_randomforest_Native_predict(JNIEnv      *env,
                                                                       jobject      obj,
                                                                       jint         traceFlag,
                                                                       jint         seedDynamic,
                                                                       jint         optLow,
                                                                       jint         optHigh,
                                                                       jint         ntree,
                                                                       jint         observationSize,
                                                                       jint         ySize,
                                                                       jcharArray   rType,
                                                                       jintArray    rLevels,
                                                                       jobject      rData,
                                                                       jint         xSize,
                                                                       jcharArray   xType,
                                                                       jintArray    xLevels,
                                                                       jobject      xData,
                                                                       jint         sampleSize,
                                                                       jobject      sample,
                                                                       jdoubleArray caseWeight,
                                                                       jint         timeInterestSize,
                                                                       jdoubleArray timeInterest,
                                                                       jintArray    seed,
                                                                       jint         totalNodeCount,
                                                                       jintArray    treeID,
                                                                       jintArray    nodeID,
                                                                       jintArray    parmID,
                                                                       jdoubleArray contPT,
                                                                       jintArray    mwcpSZ,
                                                                       jintArray    mwcpPT,
                                                                       jintArray    tnRMBR,
                                                                       jintArray    tnAMBR,
                                                                       jintArray    tnRCNT,
                                                                       jintArray    tnACNT,
                                                                       jdoubleArray tnSURV,
                                                                       jdoubleArray tnMORT,
                                                                       jdoubleArray tnNLSN,
                                                                       jdoubleArray tnCSHZ,
                                                                       jdoubleArray tnCIFN,
                                                                       jdoubleArray tnREGR,
                                                                       jintArray    tnCLAS,
                                                                       jint         yTargetSize,
                                                                       jintArray    yTargetIndex,
                                                                       jint         ptnCount,
                                                                                       
                                                                       jint         xImportanceSize,
                                                                       jintArray    xImportanceIndex,
                                                                       jint         xPartialType,
                                                                       jint         xPartialIndex,
                                                                       jint         xPartialSize,
                                                                       jdoubleArray xPartialValue,
                                                                       jdouble      x2PartialSize,
                                                                       jintArray    x2PartialIndex,
                                                                       jdoubleArray x2PartialValue,
                                                                       jint         subsetSize,
                                                                       jintArray    subsetIndex,
                                                                       jint         fnSize,
                                                                       jint         fySize,
                                                                       jdoubleArray fyData,
                                                                       jdoubleArray fxData,
                                                                       jint          numThreads) {
  setUserTraceFlag((uint) traceFlag);
  setNativeGlobalEnv(env, obj);
  int seedValue           = (int)  seedDynamic;
  RF_opt                  = (uint) optLow;
  RF_optHigh              = (uint) optHigh;
  RF_ntree                = (uint) ntree;
  RF_observationSize      = (uint) observationSize;
  RF_ySize                = (uint) ySize;
  RF_rType                = (char *) copy1DObject(rType, NATIVE_TYPE_CHARACTER, &RF_jni1DInfoListSize, TRUE);
  RF_rLevels              = (int*) copy1DObject(rLevels, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
  RF_responseIn           = (double **) copy2DObject(rData, NATIVE_TYPE_NUMERIC, &RF_jni2DInfoListSize);
  RF_xSize                = (uint) xSize;
  RF_xType                = (char *) copy1DObject(xType, NATIVE_TYPE_CHARACTER, &RF_jni1DInfoListSize, TRUE);
  RF_xLevels              = (uint *) copy1DObject(xLevels, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
  RF_observationIn        = (double **) copy2DObject(xData, NATIVE_TYPE_NUMERIC, &RF_jni2DInfoListSize);
  RF_bootstrapSize        = (uint) sampleSize;
  RF_bootstrapIn          = (uint **) copy2DObject(sample, NATIVE_TYPE_INTEGER, &RF_jni2DInfoListSize);
  RF_caseWeight           = (double *) copy1DObject(caseWeight, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_timeInterestSize     = (uint) timeInterestSize;
  RF_timeInterest         = (double *) copy1DObject(timeInterest, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_seed_                = (int *) copy1DObject(seed, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
  RF_totalNodeCount       = (uint) totalNodeCount;
  RF_treeID_              = (uint *)   copy1DObject(treeID, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
  RF_nodeID_              = (uint *)   copy1DObject(nodeID, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE); 
  RF_parmID_              = (uint *)   copy1DObject(parmID, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
  RF_contPT_              = (double *) copy1DObject(contPT, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, TRUE);
  RF_mwcpSZ_              = (uint *)   copy1DObject(mwcpSZ, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
  RF_mwcpPT_              = (uint *)   copy1DObject(mwcpPT, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
  RF_RMBR_ID_             = (uint *)   copy1DObject(tnRMBR, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
  RF_AMBR_ID_             = (uint *)   copy1DObject(tnAMBR, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
  RF_TN_RCNT_             = (uint *)   copy1DObject(tnRCNT, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
  RF_TN_ACNT_             = (uint *)   copy1DObject(tnACNT, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
  RF_TN_SURV_             = (double *) copy1DObject(tnSURV, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_TN_MORT_             = (double *) copy1DObject(tnMORT, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_TN_NLSN_             = (double *) copy1DObject(tnNLSN, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_TN_CSHZ_             = (double *) copy1DObject(tnCSHZ, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_TN_CIFN_             = (double *) copy1DObject(tnCIFN, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_TN_REGR_             = (double *) copy1DObject(tnREGR, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_TN_CLAS_             = (uint *)   copy1DObject(tnCLAS, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
  RF_rTargetCount         = (uint) yTargetSize;
  RF_rTarget              = (uint *) copy1DObject(yTargetIndex, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
  RF_ptnCount             = (uint) ptnCount;
   
  RF_intrPredictorSize    = (uint) xImportanceSize;
  RF_intrPredictor        = (uint *) copy1DObject(xImportanceIndex, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
  RF_partialType          = (uint) xPartialType;
  RF_partialXvar          = (uint) xPartialIndex;
  RF_partialLength        = (uint) xPartialSize;
  RF_partialValue         = (double *) copy1DObject(xPartialValue, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_partialLength2       = (uint) x2PartialSize;
  RF_partialXvar2         = (uint *) copy1DObject(x2PartialIndex, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
  RF_partialValue2        = (double *) copy1DObject(x2PartialValue, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_sobservationSize     = (uint) subsetSize;
  RF_sobservationIndv     = (uint *) copy1DObject(subsetIndex, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
  RF_fobservationSize     = (uint) fnSize;
  RF_frSize               = (uint) fySize;
  RF_fresponseIn          = (double **) copy2DObject(fyData, NATIVE_TYPE_NUMERIC, &RF_jni2DInfoListSize);
  RF_fobservationIn       = (double **) copy2DObject(fxData, NATIVE_TYPE_NUMERIC, &RF_jni2DInfoListSize);
  RF_numThreads           = (int)  numThreads;
  processDefaultPredict();
  rfsrc((RF_fobservationSize > 0)? RF_PRED : RF_REST, seedValue);
  put_jniEnsembleInfoList(RF_nativeIndex);
  free_jvvector(RF_jniEnsembleInfoList, 0, 1 << 6, NRUTIL_JEN_PTR);
  free_jni1DList(RF_jni1DInfoListSize);
  free_jni2DList(RF_jni2DInfoListSize);
  free_jvvector(RF_jni1DInfoList, 0, 1 << 6, NRUTIL_J1D_PTR);
  free_jvvector(RF_jni2DInfoList, 0, 1 << 6, NRUTIL_J2D_PTR);
  memoryCheck();
  return RF_java_hshmap_obj;
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
  RF_java_hshmap_cls = (*RF_java_env) -> FindClass(RF_java_env, "java/util/HashMap");
  if (RF_java_hshmap_cls == NULL) {
    RF_nativeError("\nRF-SRC:  Unable to access class for java/util/HashMap.\n");
    RF_nativeError("\nRF-SRC:  The application will now exit.\n");
    RF_nativeExit();
  }
  RF_java_hshmap_constr = (*RF_java_env) -> GetMethodID(RF_java_env, RF_java_hshmap_cls, "<init>", "(I)V");
  if (RF_java_hshmap_constr == NULL) {
    RF_nativeError("\nRF-SRC:  Unable to access constructor for class java/util/ArrayList.\n");
    RF_nativeError("\nRF-SRC:  The application will now exit.\n");
    RF_nativeExit();
  }
  RF_java_hshmap_put    = (*RF_java_env) -> GetMethodID(RF_java_env, RF_java_hshmap_cls, "put", "(Ljava/lang/Object;Ljava/lang/Object;)Ljava/lang/Object;");
  if (RF_java_hshmap_put == NULL) {
    RF_nativeError("\nRF-SRC:  Unable to access method java/util/HashMap::put().\n");
    RF_nativeError("\nRF-SRC:  The application will now exit.\n");
    RF_nativeExit();
  }
  RF_java_hshmap_obj = (*RF_java_env) -> NewObject(RF_java_env, RF_java_hshmap_cls, RF_java_hshmap_constr, 1 << 6);
  if (RF_java_hshmap_obj == NULL) {
    RF_nativeError("\nRF-SRC:  Unable to instantiate object java/util/HashMap.\n");
    RF_nativeError("\nRF-SRC:  The application will now exit.\n");
    RF_nativeExit();
  }
  RF_java_ens_cls = (*RF_java_env) -> FindClass(RF_java_env, "com/kogalur/randomforest/Ensemble");
  if (RF_java_ens_cls == NULL) {
    RF_nativeError("\nRF-SRC:  Unable to access class com/kogalur/randomforest/Ensemble.\n");
    RF_nativeError("\nRF-SRC:  The application will now exit.\n");
    exit(1);
  }
  RF_java_ens_mid = (*RF_java_env) -> GetMethodID(RF_java_env, RF_java_ens_cls, "Ensemble", "(Ljava/lang/String;BIJZI[ILjava/lang/Object;)V");
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
void *copy1DObject(jarray arr, char type, uint *index, char actual) {
  jdouble *dbuffer;
  jint    *ibuffer;
  jchar   *cbuffer;
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
    (incomingInfo -> arrayJVM) = arr;
    (incomingInfo -> len) = len;
    (incomingInfo -> type) = type;
    (incomingInfo -> actual) = actual;
    switch (type) {
    case NATIVE_TYPE_NUMERIC:
      dbuffer =  (*RF_java_env) -> GetDoubleArrayElements(RF_java_env, arr, &(incomingInfo -> isCopy));
      (incomingInfo -> arrayJNI) = dbuffer;
      if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
        RF_nativeExit();
      }
      if (incomingInfo -> actual) {
        (incomingInfo -> array) = dcopy = dvector(1, len);
        for (i = 0; i < len; i++) {
          dcopy[i+1] = (double) dbuffer[i];
        }
        if ((incomingInfo -> isCopy) == JNI_TRUE) {
          (*RF_java_env) -> ReleaseDoubleArrayElements(RF_java_env, arr, (jdouble*) dbuffer, JNI_ABORT);
          if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
            RF_nativeExit();
          }
        }
        copy = dcopy;
        (*RF_java_env) -> DeleteLocalRef(RF_java_env, arr);
      }
      else {
        copy = dbuffer - 1;
      }
      break;
    case NATIVE_TYPE_INTEGER:
      ibuffer =  (*RF_java_env) -> GetIntArrayElements(RF_java_env, arr, & (incomingInfo -> isCopy));
      (incomingInfo -> arrayJNI) = ibuffer;
      if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
        RF_nativeExit();
      }
      if (incomingInfo -> actual) {
        (incomingInfo -> array) = icopy = ivector(1, len);
        for (i = 0; i < len; i++) {
          icopy[i+1] = (int) ibuffer[i];
        }
        if ((incomingInfo -> isCopy) == JNI_TRUE) {
          (*RF_java_env) -> ReleaseIntArrayElements(RF_java_env, arr, (jint*) ibuffer, JNI_ABORT);
          if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
            RF_nativeExit();
          }
        }
        copy = icopy;
        (*RF_java_env) -> DeleteLocalRef(RF_java_env, arr);
      }
      else {
        copy = ibuffer - 1;
      }
      break;
    case NATIVE_TYPE_CHARACTER:
      cbuffer =  (*RF_java_env) -> GetCharArrayElements(RF_java_env, arr, & (incomingInfo -> isCopy));
      (incomingInfo -> arrayJNI) = cbuffer;
      if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
        RF_nativeExit();
      }
      if (incomingInfo -> actual) {
        (incomingInfo -> array) = ccopy = cvector(1, len);
        for (i = 0; i < len; i++) {
          ccopy[i+1] = (char) cbuffer[i];
        }
        if ((incomingInfo -> isCopy) == JNI_TRUE) {
          (*RF_java_env) -> ReleaseCharArrayElements(RF_java_env, arr, (jchar*) cbuffer, JNI_ABORT);
          if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
            RF_nativeExit();
          }
        }
        copy = ccopy;
        (*RF_java_env) -> DeleteLocalRef(RF_java_env, arr);
      }
      else {
        copy = cbuffer - 1;
      }
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
      if (incomingInfo -> actual) {
        free_dvector((double *) (incomingInfo -> array), 1, incomingInfo -> len);
      }
      else {
        if ((incomingInfo -> isCopy) == JNI_TRUE) {
          (*RF_java_env) -> ReleaseDoubleArrayElements(RF_java_env, incomingInfo -> arrayJVM, (jdouble *) (incomingInfo -> arrayJNI), JNI_ABORT);
          if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
            RF_nativeExit();
          }
        }
        (*RF_java_env) -> DeleteLocalRef(RF_java_env, incomingInfo -> arrayJVM);
      }
      break;
    case NATIVE_TYPE_INTEGER:
      if (incomingInfo -> actual) {
        free_ivector((int *) (incomingInfo -> array), 1, incomingInfo -> len);    
      }
      else {
        if ((incomingInfo -> isCopy) == JNI_TRUE) {
          (*RF_java_env) -> ReleaseIntArrayElements(RF_java_env, incomingInfo -> arrayJVM, (jint *) (incomingInfo -> arrayJNI), JNI_ABORT);
          if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
            RF_nativeExit();
          }
        }
        (*RF_java_env) -> DeleteLocalRef(RF_java_env, incomingInfo -> arrayJVM);
      }
      break;
    case NATIVE_TYPE_CHARACTER:
      if (incomingInfo -> actual) {
        free_cvector((char *) (incomingInfo -> array), 1, incomingInfo -> len);    
      }
      else {
        if ((incomingInfo -> isCopy) == JNI_TRUE) {
          (*RF_java_env) -> ReleaseCharArrayElements(RF_java_env, incomingInfo -> arrayJVM, (jchar *) (incomingInfo -> arrayJNI), JNI_ABORT);
          if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
            RF_nativeExit();
          }
        }
        (*RF_java_env) -> DeleteLocalRef(RF_java_env, incomingInfo -> arrayJVM);
      }
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
                      double value,
                      char **sexpString,
                      void  *auxiliaryArrayPtr,
                      uint   auxiliaryDimSize,
                      ...) {
  jboolean isCopy;
  jmethodID mid;
  jarray     thisArray;
  void      *thisArrayPtr;
  jintArray  thisDim;
  int       *thisDimPtr;
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
  RF_jniEnsembleInfoList[*index] = ensembleInfo = (JNIEnsembleInfo*) gblock((size_t) sizeof(JNIEnsembleInfo));
  thisDim = (jintArray) (*RF_java_env) -> NewIntArray(RF_java_env, auxiliaryDimSize);  
  if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
    RF_nativeExit();
  }
  thisDimPtr = (uint *) (*RF_java_env) -> GetIntArrayElements(RF_java_env, thisDim, &isCopy);
  if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
    RF_nativeExit();
  }
  va_list list;
  va_start(list, auxiliaryDimSize);
  for (uint i = 0; i < auxiliaryDimSize; i++) {
    thisDimPtr[i] = va_arg(list, int);
  }
  va_end(list);
  ensembleInfo -> type     = type;
  ensembleInfo -> size     = size;
  ensembleInfo -> identity = identity;
  ensembleInfo -> dimSize  = auxiliaryDimSize;
  ensembleInfo -> dim      = thisDim;
  ensembleInfo -> dimPtr   = thisDimPtr;
  ensembleInfo -> isCopyDim = isCopy;
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
    for (ulong i = 0; i < size; i++) {
      ((double*) thisArrayPtr)[i] = value;
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
    for (ulong i = 0; i < size; i++) {
      ((uint*) thisArrayPtr)[i] = 0;
    }
    break;
  case NATIVE_TYPE_CHARACTER:
    thisArray = (jbyteArray) (*RF_java_env) -> NewByteArray(RF_java_env, size);
    if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
      RF_nativeExit();
    }
    thisArrayPtr = (jbyteArray *) (*RF_java_env) -> GetByteArrayElements(RF_java_env, thisArray, &isCopy);
    if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
      RF_nativeExit();
    }
    for (ulong i = 0; i < size; i++) {
      ((char*) thisArrayPtr)[i] = 0;
    }
    break;
  default:
    thisArray = thisArrayPtr = NULL;
    break;
  }
  (ensembleInfo -> isCopy) = isCopy;
  (ensembleInfo -> array) = thisArray;
  (ensembleInfo -> arrayPtr) = thisArrayPtr;
  allocateAuxiliaryInfo(type,
                        identity,
                        thisArrayPtr,
                        auxiliaryArrayPtr,
                        auxiliaryDimSize,
                        thisDimPtr - 1);
  (*index) ++;
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
    if (ensembleInfo -> isCopyDim == JNI_TRUE) {
      (*RF_java_env) -> ReleaseIntArrayElements(RF_java_env,
                                                ensembleInfo -> dim,
                                                (jint *) (ensembleInfo -> dimPtr), 0);
      if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
        RF_nativeExit();
      }
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
    case NATIVE_TYPE_CHARACTER:
      if (ensembleInfo -> isCopy == JNI_TRUE) {
        (*RF_java_env) -> ReleaseByteArrayElements(RF_java_env,
                                                  ensembleInfo -> array,
                                                  (jbyte *) (ensembleInfo -> arrayPtr), 0);
        if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
          RF_nativeExit();
        }
      }
      break;
    }
    jobject ensembleObj = (*RF_java_env) -> NewObject(RF_java_env, RF_java_ens_cls, RF_java_ens_mid,
                                                      (jstring) name,
                                                      (jbyte) ensembleInfo -> type,
                                                      (jint) ensembleInfo -> identity,
                                                      (jlong) ensembleInfo -> size,
                                                      (jboolean) (((auxInfoPtr -> auxiliaryArrayPtr) == NULL) && (auxInfoPtr -> dimSize > 1)) ? JNI_TRUE : JNI_FALSE,
                                                      (jint) (ensembleInfo -> dimSize),
                                                      (jintArray) (ensembleInfo -> dim),
                                                      (jobject) (ensembleInfo -> array));
    if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
      RF_nativeExit();
    }
    (*RF_java_env) -> CallObjectMethod(RF_java_env, RF_java_hshmap_obj, RF_java_hshmap_put, name, ensembleObj);
    if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
      RF_nativeExit();
    }
    if (FALSE) {
      RF_nativePrint("\nRF-SRC:  Unable to put ensemble object to java/util/HashMap.\n");
      RF_nativePrint("\nRF-SRC:  The application will now exit.\n");
      exit(1);
    }
    (*RF_java_env) -> DeleteLocalRef(RF_java_env, ensembleInfo -> dim);      
    (*RF_java_env) -> DeleteLocalRef(RF_java_env, ensembleInfo -> array);      
    free_gblock(ensembleInfo, (size_t) sizeof(JNIEnsembleInfo));  
  }
}
void setUserTraceFlag (uint traceFlag) {
  RF_userTraceFlag = traceFlag;
}
uint getUserTraceFlag () { return RF_userTraceFlag; }
