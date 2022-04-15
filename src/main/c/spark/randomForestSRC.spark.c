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
                                                                    jobject      baseLearn,
                                                                    jint         vtry,
                                                                    jintArray    vtryArray,
                                                                    jobject      vtryExperimental,
                                                                    jint         ytry,
                                                                    jint         nodeSize,
                                                                    jint         nodeDepth,
                                                                    jint         crWeightSize,
                                                                    jdoubleArray crWeight,
                                                                    jint         ntree,
                                                                    jint         observationSize,
                                                                    jobject      yInfo,
                                                                    jintArray    yLevels,
                                                                    jdoubleArray yData,
                                                                    jobject      xInfo,
                                                                    jintArray    xLevels,
                                                                    jdoubleArray xData,
                                                                    jobject      sampleInfo,
                                                                    jdoubleArray xWeightStat,
                                                                    jdoubleArray yWeight,
                                                                    jdoubleArray xWeight,
                                                                    jdoubleArray timeInterest,
                                                                    jint         nImpute,
                                                                    jint         perfBlock,
                                                                    jobject      quantile,
                                                                    jint         xPreSort,
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
  populateBaseLearnObject(baseLearn);
  RF_vtry                 = (uint) vtry;
  RF_vtryArray            = (uint **) copy2DObject(vtryArray, NATIVE_TYPE_INTEGER, &RF_nat2DInfoListSize);
  populateVtryObject(vtryExperimental);
  RF_ytry                 = (uint) ytry;
  RF_nodeSize             = (uint) nodeSize;
  RF_nodeDepth            = (int)  nodeDepth;
  RF_crWeightSize         = (uint) crWeightSize;
  RF_crWeight             = (double *) copy1DObject(crWeight, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
  RF_ntree                = (uint) ntree;
  RF_observationSize      = (uint) observationSize;
  populateYInfoObject(yInfo);
  RF_responseIn           = (double **) copy2DObject(yData, NATIVE_TYPE_NUMERIC, &RF_nat2DInfoListSize);
  RF_rLevelsJNIE = yLevels;
  RF_rLevels = NULL;
  populateXInfoObject(xInfo);  
  RF_observationIn         = (double **) copy2DObject(xData, NATIVE_TYPE_NUMERIC, &RF_nat2DInfoListSize);
  RF_xLevelsJNIE = xLevels;
  RF_xLevels = NULL;
  populateSampleObject(sampleInfo);
  RF_xWeightStat          = (double *) copy1DObject(xWeightStat, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
  RF_yWeight              = (double *) copy1DObject(yWeight, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
  RF_xWeight              = (double *) copy1DObject(xWeight, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
  if (!((*RF_java_env) -> IsSameObject(RF_java_env, timeInterest, NULL))) {
    RF_timeInterest = (double *) copy1DObject(timeInterest, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
    RF_timeInterestSize = (*RF_java_env) -> GetArrayLength(RF_java_env, timeInterest);
  }
  else {
    RF_timeInterest = NULL;
    RF_timeInterestSize = 0;
  }
  RF_nImpute              = (uint) nImpute;
  RF_perfBlock            = (uint) perfBlock;
  populateQuantileObject(quantile);
  RF_xPreSort             = (uint) xPreSort;
  RF_numThreads           = (int)  numThreads;
  processDefaultGrow();
  rfsrc(RF_GROW, seedValue);
  put_nativeEnsembleInfoList(RF_nativeIndex);
  free_nvvector(RF_nativeEnsembleInfoList, 0, 1 << 6, NRUTIL_JEN_PTR);
  free_jni1DList(RF_nat1DInfoListSize);
  free_jni2DList(RF_nat2DInfoListSize);
  free_nvvector(RF_nat1DInfoList, 0, 1 << 6, NRUTIL_J1D_PTR);
  free_nvvector(RF_nat2DInfoList, 0, 1 << 6, NRUTIL_J2D_PTR);
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
                                                                       jobject      yInfo,
                                                                       jintArray    yLevels,
                                                                       jdoubleArray yData,
                                                                       jobject      xInfo,
                                                                       jintArray    xLevels,
                                                                       jdoubleArray xData,
                                                                       jobject      sampleInfo,
                                                                       jdoubleArray timeInterest,
                                                                       jint         totalNodeCount,
                                                                       jintArray    tLeafCount,
                                                                       jobject      seedInfo,
                                                                       jint         hdim,
                                                                       jobject      baseLearn,
                                                                       jintArray    treeID,
                                                                       jintArray    nodeID,
                                                                       jobject      hc_zero,
                                                                       jobject      hc_oneAugIntr,
                                                                       jobject      hc_oneAugSyth,
                                                                       jobject      hc_one,
                                                                       jobject      hc_parmID,
                                                                       jobject      hc_contPT,
                                                                       jobject      hc_contPTR,
                                                                       jobject      hc_mwcpSZ,
                                                                       jobject      hc_mwcpPT,
                                                                       jobject      hc_augmXone,
                                                                       jobject      hc_augmXtwo,
                                                                       jobject      hc_augmXS,
                                                                       jobject      hc_augmSythTop,
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
                                                                       jintArray    yTarget,
                                                                       jint         ptnCount,
                                                                       jintArray    xMarginal,
                                                                       jintArray    xImportance,
                                                                       jobject      partial,
                                                                       jint         fnSize,
                                                                       jint         fySize,
                                                                       jdoubleArray fyData,
                                                                       jdoubleArray fxData,
                                                                       jint         perfBlock,
                                                                       jobject      quantile,
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
  populateYInfoObject(yInfo);
  RF_responseIn           = (double **) copy2DObject(yData, NATIVE_TYPE_NUMERIC, &RF_nat2DInfoListSize);
  RF_rLevelsJNIE = yLevels;
  RF_rLevels = NULL;
  populateXInfoObject(xInfo);  
  RF_observationIn         = (double **) copy2DObject(xData, NATIVE_TYPE_NUMERIC, &RF_nat2DInfoListSize);
  RF_xLevelsJNIE = xLevels;
  RF_xLevels = NULL;
  populateSampleObject(sampleInfo);
  if (!((*RF_java_env) -> IsSameObject(RF_java_env, timeInterest, NULL))) {
    RF_timeInterest = (double *) copy1DObject(timeInterest, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
    RF_timeInterestSize = (*RF_java_env) -> GetArrayLength(RF_java_env, timeInterest);
  }
  else {
    RF_timeInterest = NULL;
    RF_timeInterestSize = 0;
  }
  RF_totalNodeCount       = (uint) totalNodeCount;
  RF_tLeafCount           = (uint *)   copy1DObject(tLeafCount, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
  populateSeedObject(seedInfo);
  RF_hdim                 = (uint) hdim;
  populateBaseLearnObject(baseLearn);
  RF_treeID_              = (uint *)   copy1DObject(treeID, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
  RF_nodeID_              = (uint *)   copy1DObject(nodeID, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
  RF_nodeSZ_              = (uint *)   copy1DObject(nodeSZ, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
  RF_RMBR_ID_             = (uint *)   copy1DObject(tnRMBR, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
  RF_AMBR_ID_             = (uint *)   copy1DObject(tnAMBR, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
  RF_TN_RCNT_             = (uint *)   copy1DObject(tnRCNT, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
  RF_TN_ACNT_             = (uint *)   copy1DObject(tnACNT, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
  RF_perfBlock            = (uint) perfBlock;
  populateQuantileObject(quantile);
  RF_numThreads           = (int)  numThreads;
  RF_ptnCount             = (uint) ptnCount;
  if (!((*RF_java_env) -> IsSameObject(RF_java_env, yTarget, NULL))) {
    RF_rTarget = (uint *) copy1DObject(yTarget, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
    RF_rTargetCount = (*RF_java_env) -> GetArrayLength(RF_java_env, yTarget);
  }
  else {
    RF_rTarget = NULL;
    RF_rTargetCount = 0;
  }
  if (!((*RF_java_env) -> IsSameObject(RF_java_env, xMarginal, NULL))) {
    RF_xMarginal = (uint *) copy1DObject(xMarginal, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
    RF_xMarginalSize = (*RF_java_env) -> GetArrayLength(RF_java_env, xMarginal);
  }
  else {
    RF_xMarginal = NULL;
    RF_xMarginalSize = 0;
  }
  if (!((*RF_java_env) -> IsSameObject(RF_java_env, xImportance, NULL))) {
    RF_intrPredictor = (uint *) copy1DObject(xImportance, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
    RF_intrPredictorSize = (*RF_java_env) -> GetArrayLength(RF_java_env, xImportance);
  }
  else {
    RF_intrPredictor = NULL;
    RF_intrPredictorSize = 0;
  }
  populatePartialObject(partial);
  RF_fobservationSize     = (uint) fnSize;
  RF_frSize               = (uint) fySize;
  RF_fresponseIn          = (double **) copy2DObject(fyData, NATIVE_TYPE_NUMERIC, &RF_nat2DInfoListSize);
  RF_fobservationIn       = (double **) copy2DObject(fxData, NATIVE_TYPE_NUMERIC, &RF_nat2DInfoListSize);
  RF_getTree              = (uint *) copy1DObject(getTree, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
  RF_TN_SURV_             = (double *) copy1DObject(tnSURV, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
  RF_TN_MORT_             = (double *) copy1DObject(tnMORT, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
  RF_TN_NLSN_             = (double *) copy1DObject(tnNLSN, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
  RF_TN_CSHZ_             = (double *) copy1DObject(tnCSHZ, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
  RF_TN_CIFN_             = (double *) copy1DObject(tnCIFN, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
  RF_TN_REGR_             = (double *) copy1DObject(tnREGR, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
  RF_TN_CLAS_             = (uint *)   copy1DObject(tnCLAS, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
  processDefaultPredict();
  mode = (RF_fobservationSize > 0)? RF_PRED : RF_REST;
  stackForestObjectsAuxOnly(mode);
  populateHyperZeroObject(hc_zero);
  populateHyperOneObject(hc_one);
  rfsrc(mode, seedValue);
  unstackForestObjectsAuxOnly(mode);
  put_nativeEnsembleInfoList(RF_nativeIndex);
  free_nvvector(RF_nativeEnsembleInfoList, 0, 1 << 6, NRUTIL_JEN_PTR);
  free_jni1DList(RF_nat1DInfoListSize);
  free_jni2DList(RF_nat2DInfoListSize);
  free_nvvector(RF_nat1DInfoList, 0, 1 << 6, NRUTIL_J1D_PTR);
  free_nvvector(RF_nat2DInfoList, 0, 1 << 6, NRUTIL_J2D_PTR);
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
  RF_nat1DInfoList = (NAT1DInfo **) nvvector(0, 1 << 6, NRUTIL_J1D_PTR);
  RF_nat2DInfoList = (NAT2DInfo **) nvvector(0, 1 << 6, NRUTIL_J2D_PTR);
  RF_nativeEnsembleInfoList = (NativeEnsembleInfo **) nvvector(0, 1 << 6, NRUTIL_JEN_PTR);
  RF_nat1DInfoListSize = 0;
  RF_nat2DInfoListSize = 0;
  RF_nativeEnsembleInfoListSize = 0;
  RF_nativeIndex = RF_stackCount = 0;
}
void *nvvector(unsigned long long nl, unsigned long long nh, enum alloc_ntype type) {
  void *v;
  v = NULL;  
  switch(type){
  case NRUTIL_J1D_PTR:
    v = ((NAT1DInfo **) (gvector(nl, nh, sizeof(NAT1DInfo *)) -nl+NR_END));
    break;
  case NRUTIL_J2D_PTR:
    v = ((NAT2DInfo **) (gvector(nl, nh, sizeof(NAT2DInfo *)) -nl+NR_END));
    break;
  case NRUTIL_JEN_PTR:
    v = ((NativeEnsembleInfo **) (gvector(nl, nh, sizeof(NativeEnsembleInfo *)) -nl+NR_END));
    break;
  default:
    v = NULL;
    nrerror("\n  Illegal case in nvvector().");
    break;
  }
  return v;
}
void free_nvvector(void *v, unsigned long long nl, unsigned long long nh, enum alloc_ntype type) {
  switch(type){
  case NRUTIL_J1D_PTR:
    free_gvector((NAT1DInfo *) (v+nl-NR_END), nl, nh, sizeof(NAT1DInfo *));
    break;
  case NRUTIL_J2D_PTR:
    free_gvector((NAT2DInfo *) (v+nl-NR_END), nl, nh, sizeof(NAT2DInfo *));
    break;
  case NRUTIL_JEN_PTR:
    free_gvector((NativeEnsembleInfo *) (v+nl-NR_END), nl, nh, sizeof(NativeEnsembleInfo *));
    break;
  default:
    nrerror("\n  Illegal case in free_nvvector().");
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
  NAT1DInfo *incomingInfo;
  copy = NULL;
  if (! (*RF_java_env) -> IsSameObject(RF_java_env, arr, NULL)) {
    if (((*index) >> 6) > 0) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  copy1DObject() list limit exceeded:  %20d", *index);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
    RF_nat1DInfoList[*index] = incomingInfo = (NAT1DInfo*) gblock((size_t) sizeof(NAT1DInfo));
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
void *copy2DObject(jobject obj, char type, uint *index) {
  NAT2DInfo *incomingInfo;
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
    RF_nat2DInfoList[*index] = incomingInfo = (NAT2DInfo*) gblock((size_t) sizeof(NAT2DInfo));
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
void free_jni1DList(uint size) {
  NAT1DInfo *incomingInfo;
  for (uint i = 0; i < size; i++) {
    incomingInfo = RF_nat1DInfoList[i];
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
    free_gblock(incomingInfo, (size_t) sizeof(NAT1DInfo));    
  }
}
void free_jni2DList(uint size) {
  NAT2DInfo *incomingInfo;
  for (uint i = 0; i < size; i++) {
    incomingInfo = RF_nat2DInfoList[i];
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
    free_gblock(incomingInfo, (size_t) sizeof(NAT2DInfo));    
  }
}
void initProtect(uint stackCount) {
}
void *stackAndProtect(char   mode,
                      uint  *index,
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
  NativeEnsembleInfo *ensembleInfo;
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
  RF_nativeEnsembleInfoList[*index] = ensembleInfo = (NativeEnsembleInfo*) gblock((size_t) sizeof(NativeEnsembleInfo));
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
  allocateAuxiliaryInfo((mode == RF_GROW) ? FALSE : TRUE,
                        type,
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
void put_nativeEnsembleInfoList(uint size) {
  jmethodID mid;
  NativeEnsembleInfo *ensembleInfo;
  SNPAuxiliaryInfo *auxInfoPtr;
  jstring name;
  uint stringLength;
  jboolean result;
  for (uint i = 0; i < size; i++) {
    ensembleInfo = RF_nativeEnsembleInfoList[i];
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
    (*RF_java_env) -> DeleteLocalRef(RF_java_env, ensembleInfo -> dim);      
    (*RF_java_env) -> DeleteLocalRef(RF_java_env, ensembleInfo -> array);      
    stringLength = strlen(ensembleInfo -> identity) + 1;
    free_cvector(ensembleInfo -> identity, 1, stringLength);
    free_gblock(ensembleInfo, (size_t) sizeof(NativeEnsembleInfo));  
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
      RF_parmID_[1] = (int *) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
    }
    else {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/HCzero : parmID \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "contPT", "[D");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_contPT_[1] = (double *) copy1DObject(arr, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
    }
    else {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/HCzero : contPT \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "mwcpSZ", "[I");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_mwcpSZ_[1] = (uint *) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
    }
    else {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/HCzero : mwcpSZ \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "mwcpPT", "[I");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_mwcpPT_[1] = (uint *) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
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
        RF_hcDim_ = (uint *) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
      }
      else {
        RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/HCone : hcDim \n");
        RF_nativeError("\nRF-SRC:  The application will now exit.\n");
        RF_nativeExit();
      }
      objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "contPTR", "[D");
      if (objFieldID != NULL) {
        jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
        RF_contPTR_[1] = (double *) copy1DObject(arr, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
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
    RF_partialValue = (double *) copy1DObject(arr, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
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
      RF_partialXvar2 = (uint *) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
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
      RF_partialValue2 = (double *) copy1DObject(arr, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
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
void populateBaseLearnObject(jobject obj) {
  jclass objClass;
  jfieldID objFieldID;
  RF_baseLearnDepthINTR = 0;
  RF_baseLearnRuleINTR  = AUGT_INTR_NONE;
  RF_baseLearnDepthSYTH = 0;
  RF_baseLearnDimReduce = FALSE;
  if (! (*RF_java_env) -> IsSameObject(RF_java_env, obj, NULL)) {
    objClass = (*RF_java_env) -> GetObjectClass(RF_java_env, obj);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "intrDepth", "I");
    if (objFieldID == NULL) {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/Lot : intrDepth \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    RF_baseLearnDepthINTR = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
    if (RF_baseLearnDepthINTR > 1) {
      objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "intrRule", "I");
      if (objFieldID == NULL) {
        RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/Lot : intrRule \n");
        RF_nativeError("\nRF-SRC:  The application will now exit.\n");
        RF_nativeExit();
      }
      RF_baseLearnRuleINTR = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "syth", "I");
    if (objFieldID == NULL) {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/Lot : syth \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    RF_baseLearnDepthSYTH = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "dimReduce", "Z");
    if (objFieldID == NULL) {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/Lot : dimReduce \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    RF_baseLearnDimReduce = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
  }
}
void populateYInfoObject(jobject obj) {
  jclass objClass;
  jfieldID objFieldID;
  if (!((*RF_java_env) -> IsSameObject(RF_java_env, obj, NULL))) {
    objClass = (*RF_java_env) -> GetObjectClass(RF_java_env, obj);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "ySize", "I");
    if (objFieldID == NULL) {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/YInfo : ySize \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    RF_ySize = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "yType", "[C");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_rType = (char *) copy1DObject(arr, NATIVE_TYPE_CHARACTER, &RF_nat1DInfoListSize, FALSE);
    }
    else {
      RF_rType = NULL;
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "yLevelsMax", "[I");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_rLevelsMax = (uint*) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
    }
    else {
      RF_rLevelsMax = NULL;
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "yLevelsCnt", "[I");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_rLevelsCnt = (uint*) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
    }
    else {
      RF_rLevelsCnt = NULL;
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "subjID", "[I");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_subjIn = (uint*) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
    }
    else {
      RF_subjIn = NULL;
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "eventTypeSize", "I");
    if (objFieldID != NULL) {
      RF_eventTypeSize = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
    }
    else {
      RF_eventTypeSize = 0;
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "eventType", "[I");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_eventType = (uint*) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
    }
    else {
      RF_eventType = NULL;
    }
  }
  else {
    RF_ySize = 0;
    RF_rType = NULL;
    RF_rLevelsMax = NULL;
    RF_rLevelsCnt = NULL;
    RF_subjIn = NULL;
    RF_eventTypeSize = 0;
    RF_eventType = NULL;
  }
}
void populateXInfoObject(jobject obj) {
  jclass objClass;
  jfieldID objFieldID;
  if (!((*RF_java_env) -> IsSameObject(RF_java_env, obj, NULL))) {
    objClass = (*RF_java_env) -> GetObjectClass(RF_java_env, obj);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "xSize", "I");
    if (objFieldID == NULL) {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/XInfo : xSize \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    RF_xSize = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "xType", "[C");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_xType = (char *) copy1DObject(arr, NATIVE_TYPE_CHARACTER, &RF_nat1DInfoListSize, FALSE);
    }
    else {
      RF_xType = NULL;
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "xLevelsMax", "[I");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_xLevelsMax = (uint*) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
    }
    else {
      RF_xLevelsMax = NULL;
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "xLevelsCnt", "[I");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_xLevelsCnt = (uint*) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
    }
    else {
      RF_xLevelsCnt = NULL;
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "xtType", "[I");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_xtType = (uint *) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
    }
    else {
      RF_xtType = NULL;
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "stType", "[I");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_stType = (uint *) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
    }
    else {
      RF_stType = NULL;
    }
  }
  else {
    RF_xSize = 0;
    RF_xType = NULL;
    RF_xLevelsMax = NULL;
    RF_xLevelsCnt = NULL;
    RF_xtType = NULL;
    RF_stType = NULL;
  }
}
void populateSampleObject(jobject obj) {
  jclass objClass;
  jfieldID objFieldID;
  RF_subjSize             = RF_observationSize;
  RF_subjWeight           = NULL;
  RF_bootstrapSize        = RF_observationSize;
  RF_bootstrapIn          = NULL;
  if (!((*RF_java_env) -> IsSameObject(RF_java_env, obj, NULL))) {
    objClass = (*RF_java_env) -> GetObjectClass(RF_java_env, obj);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "subjSize", "I");
    if (objFieldID != NULL) {
      RF_subjSize = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "subjWeight", "[D");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_subjWeight = (double *) copy1DObject(arr, NATIVE_TYPE_NUMERIC, &RF_nat1DInfoListSize, FALSE);
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "bootstrapSize", "I");
    if (objFieldID != NULL) {
      RF_bootstrapSize = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
      objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "byUser", "I");
      if (objFieldID != NULL) {
        jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
        RF_bootstrapIn = (uint **) copy2DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat2DInfoListSize);
      }
    }
    else {
      RF_bootstrapSize        = RF_subjSize;
    }
  }
}
void populateVtryObject(jobject obj) {
  jclass objClass;
  jfieldID objFieldID;
  if (!((*RF_java_env) -> IsSameObject(RF_java_env, obj, NULL))) {
    objClass = (*RF_java_env) -> GetObjectClass(RF_java_env, obj);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "vtry", "I");
    if (objFieldID != NULL) {
      RF_vtry = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
    }
    if (RF_vtry > 0) {
      objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "vtryArray", "[I");
      if (objFieldID != NULL) {
        jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
        RF_vtryArray = (uint **) copy2DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat2DInfoListSize);
        objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "vtryBlockSize", "I");
        if (objFieldID != NULL) {
          RF_vtryBlockSize = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
        }
        objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "vtryMode", "I");
        if (objFieldID != NULL) {
          RF_vtryMode = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
        }
      }
      else {
        RF_vtry = 0;
      }
    }
  }
}
void populateSeedObject(jobject obj) {
  jclass objClass;
  jfieldID objFieldID;
  if ((*RF_java_env) -> IsSameObject(RF_java_env, obj, NULL)) {
    RF_nativeError("\nRF-SRC:  Incoming object com/kogalur/randomforest/Seed is NULL. \n");
    RF_nativeError("\nRF-SRC:  The application will now exit.\n");
    RF_nativeExit();
  }
  else {
    objClass = (*RF_java_env) -> GetObjectClass(RF_java_env, obj);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "seed", "I");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_seed_ = (int *) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
    }
    else {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/Seed : seed \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "seedVimp", "I");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_seedVimp_ = (int *) copy1DObject(arr, NATIVE_TYPE_INTEGER, &RF_nat1DInfoListSize, FALSE);
    }
    else {
      RF_seedVimp_ = NULL;
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "optLowGrow", "I");
    if (objFieldID == NULL) {
      RF_nativeError("\nRF-SRC:  Unable to access field in class com/kogalur/randomforest/Seed : optLowGrow \n");
      RF_nativeError("\nRF-SRC:  The application will now exit.\n");
      RF_nativeExit();
    }
    RF_optLoGrow = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
  }
}
void populateQuantileObject(jobject obj) {
  jclass objClass;
  jfieldID objFieldID;
  if (!((*RF_java_env) -> IsSameObject(RF_java_env, obj, NULL))) {
    objClass = (*RF_java_env) -> GetObjectClass(RF_java_env, obj);
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "quantileSize", "I");
    if (objFieldID != NULL) {
      RF_quantileSize = (uint) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
    }
    else {
      RF_quantileSize = 0;
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "quantile", "[D");
    if (objFieldID != NULL) {
      jobject arr = (*RF_java_env) -> GetObjectField(RF_java_env, obj, objFieldID);
      RF_quantile = (double *) copy1DObject(arr, NATIVE_TYPE_CHARACTER, &RF_nat1DInfoListSize, FALSE);
    }
    else {
      RF_quantile = NULL;
    }
    objFieldID = (*RF_java_env) -> GetFieldID(RF_java_env, objClass, "qEpsilon", "D");
    if (objFieldID != NULL) {
      RF_qEpsilon = (double) (*RF_java_env) -> GetIntField(RF_java_env, objClass, objFieldID);
    }
    else {
      RF_qEpsilon = 0;
    }
  }
  else {
    RF_quantileSize = 0;
    RF_quantile = NULL;
    RF_qEpsilon = 0.0;
  }
}
