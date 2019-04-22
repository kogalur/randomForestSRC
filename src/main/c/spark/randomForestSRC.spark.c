JNIEXPORT jobject JNICALL Java_com_kogalur_randomforest_Native_grow(JNIEnv      *env,
                                                                    jobject      obj,
                                                                    jint         traceFlag,
                                                                    jint         seedDynamic,
                                                                    jint         optLow,
                                                                    jint         optHigh,
                                                                    jint         splitRule,
                                                                    jint         nsplit,
                                                                    jint         mtry,
                                                                    jobject      lot,
                                                                    jint         vtry,
                                                                    jintArray    vtryArray,
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
                                                                    jint         perfBlock,
                                                                    jint         quantileSize,
                                                                    jdoubleArray quantile,
                                                                    jdouble      qEpsilon,
                                                                    jdouble      wibsTau,
                                                                    jint         numThreads) {
  setUserTraceFlag((uint) traceFlag);
  setNativeGlobalEnv(env, obj);
  int seedValue           = (int)  seedDynamic;
  RF_opt                  = (uint) optLow;
  RF_optHigh              = (uint) optHigh;
  RF_splitRule            = (uint) splitRule;
  RF_nsplit               = (uint) nsplit;
  RF_mtry                 = (uint) mtry;
  populateLotObject(lot);
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
  RF_perfBlock            = (uint) perfBlock;
  RF_quantileSize         = (uint) quantileSize;
  RF_quantile             = (double *) copy1DObject(quantile, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_qEpsilon             = (uint) qEpsilon;
  RF_vtry                 = (uint) vtry;
  RF_vtryArray            = (uint **) copy2DObject(vtryArray, NATIVE_TYPE_INTEGER, &RF_jni2DInfoListSize);
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
  initProtect(0);
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
                                                                       jint         totalNodeCount,
                                                                       jintArray    seed,
                                                                       jint         hdim,
                                                                       jintArray    treeID,
                                                                       jintArray    nodeID,
                                                                       jobject      hc_zero,
                                                                       jobject      hc_one,
                                                                       jobject      hc_parmID,
                                                                       jobject      hc_contPT,
                                                                       jobject      hc_contPTR,
                                                                       jobject      hc_mwcpSZ,
                                                                       jobject      hc_mwcpPT,
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
                                                                       jobject      partial,
                                                                       jint         subsetSize,
                                                                       jintArray    subsetIndex,
                                                                       jint         fnSize,
                                                                       jint         fySize,
                                                                       jdoubleArray fyData,
                                                                       jdoubleArray fxData,
                                                                       jint         perfBlock,
                                                                       jint         quantileSize,
                                                                       jdoubleArray quantile,
                                                                       jdouble      qEpsilon,
                                                                       jintArray    getTree,
                                                                       jint         numThreads) {
  char mode;
  uint i;
  setUserTraceFlag((uint) traceFlag);
  setNativeGlobalEnv(env, obj);
  int seedValue           = (int)  seedDynamic;
  RF_opt                  = (uint) optLow;
  RF_optHigh              = (uint) optHigh;
  RF_ntree                = (uint) ntree;
  RF_observationSize      = (uint) observationSize;
  RF_ySize                = (uint) ySize;
  RF_rType                = (char *) copy1DObject(rType, NATIVE_TYPE_CHARACTER, &RF_jni1DInfoListSize, TRUE);
  RF_rLevels              = (uint*) copy1DObject(rLevels, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
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
  RF_totalNodeCount       = (uint) totalNodeCount;
  RF_seed_                = (int *) copy1DObject(seed, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
  RF_hdim                 = (uint) hdim;
  RF_treeID_              = (uint *)   copy1DObject(treeID, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
  RF_nodeID_              = (uint *)   copy1DObject(nodeID, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
  RF_RMBR_ID_             = (uint *)   copy1DObject(tnRMBR, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
  RF_AMBR_ID_             = (uint *)   copy1DObject(tnAMBR, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
  RF_TN_RCNT_             = (uint *)   copy1DObject(tnRCNT, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
  RF_TN_ACNT_             = (uint *)   copy1DObject(tnACNT, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
  RF_perfBlock            = (uint) perfBlock;
  RF_quantileSize         = (uint) quantileSize;
  RF_quantile             = (double *) copy1DObject(quantile, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_qEpsilon             = (uint) qEpsilon;
  RF_numThreads           = (int)  numThreads;
  RF_ptnCount             = (uint) ptnCount;
  RF_rTargetCount         = (uint) yTargetSize;
  RF_rTarget              = (uint *) copy1DObject(yTargetIndex, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
  RF_intrPredictorSize    = (uint) xImportanceSize;
  RF_intrPredictor        = (uint *) copy1DObject(xImportanceIndex, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
  RF_sobservationSize     = (uint) subsetSize;
  RF_sobservationIndv     = (uint *) copy1DObject(subsetIndex, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
   
  populatePartialObject(partial);
  RF_fobservationSize     = (uint) fnSize;
  RF_frSize               = (uint) fySize;
  RF_fresponseIn          = (double **) copy2DObject(fyData, NATIVE_TYPE_NUMERIC, &RF_jni2DInfoListSize);
  RF_fobservationIn       = (double **) copy2DObject(fxData, NATIVE_TYPE_NUMERIC, &RF_jni2DInfoListSize);
  RF_getTree              = (uint *) copy1DObject(getTree, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
  RF_TN_SURV_             = (double *) copy1DObject(tnSURV, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_TN_MORT_             = (double *) copy1DObject(tnMORT, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_TN_NLSN_             = (double *) copy1DObject(tnNLSN, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_TN_CSHZ_             = (double *) copy1DObject(tnCSHZ, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_TN_CIFN_             = (double *) copy1DObject(tnCIFN, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_TN_REGR_             = (double *) copy1DObject(tnREGR, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
  RF_TN_CLAS_             = (uint *)   copy1DObject(tnCLAS, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
  processDefaultPredict();
  mode = (RF_fobservationSize > 0)? RF_PRED : RF_REST;
  stackAuxForestObjects(mode);
  populateHyperZeroObject(hc_zero);
  populateHyperOneObject(hc_one);
  rfsrc(mode, seedValue);
  unstackAuxForestObjects(mode);
  put_jniEnsembleInfoList(RF_nativeIndex);
  free_jvvector(RF_jniEnsembleInfoList, 0, 1 << 6, NRUTIL_JEN_PTR);
  free_jni1DList(RF_jni1DInfoListSize);
  free_jni2DList(RF_jni2DInfoListSize);
  free_jvvector(RF_jni1DInfoList, 0, 1 << 6, NRUTIL_J1D_PTR);
  free_jvvector(RF_jni2DInfoList, 0, 1 << 6, NRUTIL_J2D_PTR);
  memoryCheck();
  initProtect(0);
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
  RF_java_hshmap_cls = (*RF_java_env) -> FindClass(RF_java_env, "java/util/LinkedHashMap");
  if (RF_java_hshmap_cls == NULL) {
    RF_nativeError("\nRF-SRC:  Unable to access class for java/util/LinkedHashMap.\n");
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
    RF_nativeError("\nRF-SRC:  Unable to access method java/util/LinkedHashMap::put().\n");
    RF_nativeError("\nRF-SRC:  The application will now exit.\n");
    RF_nativeExit();
  }
  RF_java_hshmap_obj = (*RF_java_env) -> NewObject(RF_java_env, RF_java_hshmap_cls, RF_java_hshmap_constr, 1 << 6);
  if (RF_java_hshmap_obj == NULL) {
    RF_nativeError("\nRF-SRC:  Unable to instantiate object java/util/LinkedHashMap.\n");
    RF_nativeError("\nRF-SRC:  The application will now exit.\n");
    RF_nativeExit();
  }
  RF_java_ens_cls = (*RF_java_env) -> FindClass(RF_java_env, "com/kogalur/randomforest/Ensemble");
  if (RF_java_ens_cls == NULL) {
    RF_nativeError("\nRF-SRC:  Unable to access class com/kogalur/randomforest/Ensemble.\n");
    RF_nativeError("\nRF-SRC:  The application will now exit.\n");
    RF_nativeExit();
  }
  RF_java_ens_mid = (*RF_java_env) -> GetMethodID(RF_java_env, RF_java_ens_cls, "Ensemble", "(Ljava/lang/String;BIJI[ILjava/lang/Object;)V");
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
void *jvvector(unsigned long long nl, unsigned long long nh, enum alloc_jtype type) {
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
void free_jvvector(void *v, unsigned long long nl, unsigned long long nh, enum alloc_jtype type) {
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
    if (((*index) >> 6) > 0) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  copy1DObject() list limit exceeded:  %20d", *index);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
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
  void *result;
  uint outLen;
  jboolean isCopy;
  result = NULL;
  if (! (*RF_java_env) -> IsSameObject(RF_java_env, obj, NULL)) {
    if (((*index) >> 6) > 0) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  copy2DObject() list limit exceeded:  %20d", *index);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
    RF_jni2DInfoList[*index] = incomingInfo = (JNI2DInfo*) gblock((size_t) sizeof(JNI2DInfo));
    outLen = (*RF_java_env) -> GetArrayLength(RF_java_env, obj);
    if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
      RF_nativeExit();
    }
    switch (type) {
    case NATIVE_TYPE_NUMERIC:
      result = (incomingInfo -> outerPtr) = (double **) new_vvector(1, outLen, NRUTIL_DPTR);
      break;
    case NATIVE_TYPE_INTEGER:
      result = (incomingInfo -> outerPtr) = (uint **) new_vvector(1, outLen, NRUTIL_UPTR);
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
  return result;
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
                      uint   unused,
                      ulong  size,
                      double value,
                      char  *sexpString,
                      void  *auxiliaryArrayPtr,
                      uint   auxiliaryDimSize,
                      ...) {
  jboolean isCopy;
  jmethodID mid;
  jarray     thisArray;
  void      *thisArrayPtr;
  jintArray  thisDim;
  int       *thisDimPtr;
  uint stringLength;
  JNIEnsembleInfo *ensembleInfo;
  if (sizeof(ulong) > sizeof(uint)) {
    if (size > UINT_MAX) {
      if (TRUE) {
        RF_nativePrint("\nRF-SRC:  *** WARNING *** ");
        RF_nativePrint("\nRF-SRC:  S.E.X.P. vector element length exceeds 32-bits:  %20lu", size);
        RF_nativePrint("\nRF-SRC:  S.E.X.P. ALLOC:  %s ", sexpString);
        RF_nativePrint("\nRF-SRC:  Please Reduce Dimensionality If Possible.");
      }
    }
  }
  RF_jniEnsembleInfoList[*index] = ensembleInfo = (JNIEnsembleInfo*) gblock((size_t) sizeof(JNIEnsembleInfo));
  thisDim = (jintArray) (*RF_java_env) -> NewIntArray(RF_java_env, auxiliaryDimSize);  
  if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
    RF_nativeExit();
  }
  thisDimPtr = (int *) (*RF_java_env) -> GetIntArrayElements(RF_java_env, thisDim, &isCopy);
  if((*RF_java_env) -> ExceptionCheck(RF_java_env)) {
    RF_nativeExit();
  }
  va_list list;
  va_start(list, auxiliaryDimSize);
  for (uint i = 0; i < auxiliaryDimSize; i++) {
    thisDimPtr[i] = va_arg(list, int);
  }
  va_end(list);
  stringLength = strlen(sexpString) + 1;
  ensembleInfo -> identity = cvector(1, stringLength);
  strcpy(ensembleInfo -> identity, sexpString);
  ensembleInfo -> slot     = *index;
  ensembleInfo -> type     = type;
  ensembleInfo -> size     = size;
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
                        sexpString,
                        RF_snpAuxiliaryInfoList,
                        *index,
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
  uint stringLength;
  jboolean result;
  for (uint i = 0; i < size; i++) {
    ensembleInfo = RF_jniEnsembleInfoList[i];
    name = (*RF_java_env) -> NewStringUTF(RF_java_env, (const char*) ensembleInfo -> identity);
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
                                                      (jint) ensembleInfo -> slot,
                                                      (jlong) ensembleInfo -> size,
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
      RF_nativePrint("\nRF-SRC:  Unable to put ensemble object to java/util/LinkedHashMap.\n");
      RF_nativePrint("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    (*RF_java_env) -> DeleteLocalRef(RF_java_env, ensembleInfo -> dim);      
    (*RF_java_env) -> DeleteLocalRef(RF_java_env, ensembleInfo -> array);      
    stringLength = strlen(ensembleInfo -> identity) + 1;
    free_cvector(ensembleInfo -> identity, 1, stringLength);
    free_gblock(ensembleInfo, (size_t) sizeof(JNIEnsembleInfo));  
  }
}
void setUserTraceFlag (uint traceFlag) {
  RF_userTraceFlag = traceFlag;
}
uint getUserTraceFlag () { return RF_userTraceFlag; }
void populateHyperZeroObject(jobject obj) {
  jclass objClass;
  jfieldID objFieldID;
  if ((*RF_java_env) -> IsSameObject(RF_java_env, obj, NULL)) {
    RF_nativeError("\nRF-SRC:  Incoming object com/kogalur/randomforest/HCzero is NULL. \n");
    RF_nativeError("\nRF-SRC:  The application will now exit.\n");
    RF_nativeExit();
  }
  else {
    objClass = (*RF_java_env) -> GetObjectClass(RF_java_env, obj);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "parmID", "[I");
    if (objFieldID != NULL) {
      jintArray arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_parmID_[1] = (uint *) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
    }
    else {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/HCzero : parmID \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "contPT", "[D");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_contPT_[1] = (double *) copy1DObject(arr, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, TRUE);
    }
    else {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/HCzero : contPT \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "mwcpSZ", "[I");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_mwcpSZ_[1] = (uint *) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
    }
    else {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/HCzero : mwcpSZ \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "mwcpPT", "[I");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_mwcpPT_[1] = (uint *) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
    }
    else {
      RF_mwcpPT_[1] = NULL;      
    }
  }
}
void populateHyperOneObject(jobject obj) {
  jclass objClass;
  jfieldID objFieldID;
  if (RF_hdim > 0) {
    if ((*RF_java_env) -> IsSameObject(RF_java_env, obj, NULL)) {
      RF_nativeError("\nRF-SRC:  Incoming object com/kogalur/randomforest/HCone is NULL. \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    else {
      objClass = (*RF_java_env) -> GetObjectClass(RF_java_env, obj);
      objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "hcDim", "[I");
      if (objFieldID != NULL) {
        jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
        RF_hcDim_ = (uint *) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, TRUE);
      }
      else {
        RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/HCone : hcDim \n");
        RF_nativeError("\nRF-SRC:  The application will now exit.\n");
        RF_nativeExit();
      }
      objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "contPTR", "[D");
      if (objFieldID != NULL) {
        jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
        RF_contPTR_[1] = (double *) copy1DObject(arr, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, TRUE);
      }
      else {
        RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/HCone : contPTR \n");
        RF_nativeError("\nRF-SRC:  The application will now exit.\n");
        RF_nativeExit();
      }
    }
  }
  else {
    RF_hcDim_      = NULL;
    RF_contPTR_[1] = NULL;
  }
}
void populatePartialObject(jobject obj) {
  jclass objClass;
  jfieldID objFieldID;
  if (! (*RF_java_env) -> IsSameObject(RF_java_env, obj, NULL)) {
    objClass = (*RF_java_env) -> GetObjectClass(RF_java_env, obj);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "partialType", "I");
    if (objFieldID == NULL) {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/Partial : partialType \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    RF_partialType = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "partialXvar", "I");
    if (objFieldID == NULL) {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/Partial : partialXvar \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    RF_partialXvar = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "partialLength", "I");
    if (objFieldID == NULL) {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/Partial : partialLength \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    RF_partialLength = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "partialValue", "[D");
    if (objFieldID == NULL) {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/Partial : partialLength \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
    if (RF_partialLength != (*RF_java_env) -> GetArrayLength(RF_java_env, arr)) {
      RF_nativeError("\nRF-SRC:  Incoming length of array in class com/kogalur/randomforest/Partial : partialValue \n");
      RF_nativeError("\nRF-SRC:  is inconsistent with specified length:  %10d vs %10d \n", RF_partialLength, (*RF_java_env) -> GetArrayLength(RF_java_env, arr));
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    RF_partialValue = (double *) copy1DObject(arr, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "partialLength2", "I");
    if (objFieldID != NULL) {
      RF_partialLength2 = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
    }
    else {
      RF_partialLength2 = 0;
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "partialxVar2", "[I");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      if (RF_partialLength2 != (*RF_java_env) -> GetArrayLength(RF_java_env, arr)) {
        RF_nativeError("\nRF-SRC:  Incoming length of array in class com/kogalur/randomforest/Partial : partialXvar2 \n");
        RF_nativeError("\nRF-SRC:  is inconsistent with specified length:  %10d vs %10d \n", RF_partialLength2, (*RF_java_env) -> GetArrayLength(RF_java_env, arr));
        RF_nativeError("\nRF-SRC:  The application will now exit.\n");
        RF_nativeExit();
      }
      RF_partialXvar2 = (uint *) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_jni1DInfoListSize, FALSE);
    }
    else {
      RF_partialXvar2 = NULL;
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "partialxVar2", "[D");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      if (RF_partialLength2 != (*RF_java_env) -> GetArrayLength(RF_java_env, arr)) {
        RF_nativeError("\nRF-SRC:  Incoming length of array in class com/kogalur/randomforest/Partial : partialValue2 \n");
        RF_nativeError("\nRF-SRC:  is inconsistent with specified length:  %10d vs %10d \n", RF_partialLength2, (*RF_java_env) -> GetArrayLength(RF_java_env, arr));
        RF_nativeError("\nRF-SRC:  The application will now exit.\n");
        RF_nativeExit();
      }
      RF_partialValue2 = (double *) copy1DObject(arr, NATIVE_TYPE_NUMERIC, &RF_jni1DInfoListSize, FALSE);
    }
    else {
      RF_partialValue2 = NULL;
    }
  }
}
void populateLotObject(jobject obj) {
  jclass objClass;
  jfieldID objFieldID;
  RF_hdim = 0;
  RF_lotSize = 0;
  RF_lotLag = 0;
  RF_lotStrikeout = 0;
  if (! (*RF_java_env) -> IsSameObject(RF_java_env, obj, NULL)) {
    objClass = (*RF_java_env) -> GetObjectClass(RF_java_env, obj);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "hdim", "I");
    if (objFieldID == NULL) {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/Lot : hdim \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    RF_hdim = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
    if (RF_hdim > 0) {
      objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "treesize", "I");
      if (objFieldID == NULL) {
        RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/Lot : treesize \n");
        RF_nativeError("\nRF-SRC:  The application will now exit.\n");
        RF_nativeExit();
      }
      RF_lotSize = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
      objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "lag", "I");
      if (objFieldID == NULL) {
        RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/Lot : lag \n");
        RF_nativeError("\nRF-SRC:  The application will now exit.\n");
        RF_nativeExit();
      }
      RF_lotLag = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
      objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "strikeout", "I");
      if (objFieldID == NULL) {
        RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/Lot : strikeout \n");
        RF_nativeError("\nRF-SRC:  The application will now exit.\n");
        RF_nativeExit();
      }
      RF_lotStrikeout = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
    }
  }
}
