
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "leafLink.h"
#include "termOps.h"
#include "nrutil.h"
${trace.token} #include "error.h"
LeafLinkedObj *makeLeafLinkedObj(void) {
  LeafLinkedObj *obj = (LeafLinkedObj*) gblock((size_t) sizeof(LeafLinkedObj));
  obj -> fwdLink = NULL;
  obj -> bakLink = NULL;
  obj -> nodePtr = NULL;
  obj -> termPtr = NULL;
  obj -> termPtrAux = NULL;
  obj -> nodeID = 0;
  obj -> ibgMembrCount = 0;
  obj -> allMembrCount = 0;
  obj -> oobMembrCount = 0;
  return obj;
}
LeafLinkedObjSimple *makeLeafLinkedObjSimple(void) {
  LeafLinkedObjSimple *obj = (LeafLinkedObjSimple*) gblock((size_t) sizeof(LeafLinkedObjSimple));
  obj -> fwdLink = NULL;
  obj -> bakLink = NULL;
  obj -> nodePtr = NULL;
  return obj;
}
LeafLinkedObj *makeAndSpliceLeafLinkedObj(LeafLinkedObj *tail,
                                          Node *nodePtr,
                                          uint ibgCount,
                                          uint allCount) {
  ${trace.token}  if (getTraceFlag(0) & NODE_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nmakeAndSpliceLeafLinkedObj() ENTRY ...\n");
  ${trace.token}  }
  LeafLinkedObj *obj = makeLeafLinkedObj();
  tail -> fwdLink = obj;
  obj -> bakLink = tail;
  obj -> nodePtr = nodePtr;
  obj -> termPtr = makeTerminal();
  ${trace.token}  if (getTraceFlag(0) & NODE_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nLeaf Linked Object spliced for terminal nodeID:  %10d", nodePtr -> nodeID);
  ${trace.token}  }
  (obj -> termPtr) -> mate = obj -> nodePtr;
  (obj -> nodePtr) -> mate = obj -> termPtr;
  (obj -> termPtr) -> nodeID = obj -> nodeID = nodePtr -> nodeID;
  obj -> ibgMembrCount = ibgCount;
  obj -> allMembrCount = allCount;
  ${trace.token}  if (getTraceFlag(0) & NODE_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nmakeAndSpliceLeafLinkedObj() EXIT ...\n");
  ${trace.token}  }
  return obj;
}
LeafLinkedObjSimple *makeAndSpliceLeafLinkedObjSimple(LeafLinkedObjSimple *tail,
                                                      Node *nodePtr) {
  ${trace.token}  if (getTraceFlag(0) & NODE_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nmakeAndSpliceLeafLinkedObjSimple() ENTRY ...\n");
  ${trace.token}  }
  LeafLinkedObjSimple *obj = makeLeafLinkedObjSimple();
  tail -> fwdLink = obj;
  obj -> bakLink = tail;
  obj -> nodePtr = nodePtr;
  ${trace.token}  if (getTraceFlag(0) & NODE_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nLeaf Linked Object Simple spliced for node nodeID:  %10d", nodePtr -> nodeID);
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(0) & NODE_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nmakeAndSpliceLeafLinkedObjSimple() EXIT ...\n");
  ${trace.token}  }
  return obj;
}
void freeLeafLinkedObj(LeafLinkedObj *obj) {
  ${trace.token}  if (getTraceFlag(0) & NODE_DEF_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfreeLeafLinkedObj() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}  }
  if (obj -> termPtr != NULL) {
    ${trace.token}  if (getTraceFlag(0) & NODE_DEF_TRACE) {
    ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {
    ${trace.token}      RF_nativePrint("\nFreeing linked object and terminal node pointer at:  %10x, %10x", obj, obj -> termPtr);
    ${trace.token}    }
    ${trace.token}  }
    freeTerminal(obj -> termPtr);
    obj -> termPtr = NULL;
  }
  free_gblock(obj, (size_t) sizeof(LeafLinkedObj));
  ${trace.token}  if (getTraceFlag(0) & NODE_DEF_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {  
  ${trace.token}    RF_nativePrint("\nfreeLeafLinkedObj() EXIT ...\n");
  ${trace.token}  }
  ${trace.token}  }
}
void freeLeafLinkedObjSimple(LeafLinkedObjSimple *obj) {
  ${trace.token}  if (getTraceFlag(0) & NODE_DEF_TRACE) {
  ${trace.token}  if (getTraceFlag(0) & !TURN_OFF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfreeLeafLinkedObjSimple() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}  }
  free_gblock(obj, (size_t) sizeof(LeafLinkedObjSimple));
  ${trace.token}  if (getTraceFlag(0) & NODE_DEF_TRACE) {
  ${trace.token}  if (getTraceFlag(0) & !TURN_OFF_TRACE) {  
  ${trace.token}    RF_nativePrint("\nfreeLeafLinkedObjSimple() EXIT ...\n");
  ${trace.token}  }
  ${trace.token}  }
}
void freeLeafLinkedObjList(LeafLinkedObj *obj) {
  ${trace.token}  if (getTraceFlag(0) & NODE_DEF_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {    
  ${trace.token}    RF_nativePrint("\nfreeLeafLinkedObjList() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}  }
  if (obj -> fwdLink != NULL) {
    freeLeafLinkedObjList(obj -> fwdLink);
  }
  freeLeafLinkedObj(obj);
  ${trace.token}  if (getTraceFlag(0) & NODE_DEF_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & !TURN_OFF_TRACE) {  
  ${trace.token}    RF_nativePrint("\nfreeLeafLinkedObjList() EXIT ...\n");
  ${trace.token}  }
  ${trace.token}  }
}
void freeLeafLinkedObjListRev(LeafLinkedObj *obj) {
  ${trace.token}  if (getTraceFlag(0) & NODE_DEF_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & TURN_OFF_TRACE) {    
  ${trace.token}    RF_nativePrint("\nfreeLeafLinkedObjListRev() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}  }
  if (obj -> bakLink != NULL) {
    freeLeafLinkedObjListRev(obj -> bakLink);
  }
  freeLeafLinkedObj(obj);
  ${trace.token}  if (getTraceFlag(0) & NODE_DEF_TRACE) {
  ${trace.token}    if (getTraceFlag(0) & TURN_OFF_TRACE) {  
  ${trace.token}    RF_nativePrint("\nfreeLeafLinkedObjListRev() EXIT ...\n");
  ${trace.token}  }
  ${trace.token}  }
}
