
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "leafLink.h"
#include "termOps.h"
#include "nrutil.h"
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
  LeafLinkedObj *obj = makeLeafLinkedObj();
  tail -> fwdLink = obj;
  obj -> bakLink = tail;
  obj -> nodePtr = nodePtr;
  obj -> termPtr = makeTerminal();
  (obj -> termPtr) -> mate = obj -> nodePtr;
  (obj -> nodePtr) -> mate = obj -> termPtr;
  (obj -> termPtr) -> nodeID = obj -> nodeID = nodePtr -> nodeID;
  obj -> ibgMembrCount = ibgCount;
  obj -> allMembrCount = allCount;
  return obj;
}
LeafLinkedObjSimple *makeAndSpliceLeafLinkedObjSimple(LeafLinkedObjSimple *tail,
                                                      Node *nodePtr) {
  LeafLinkedObjSimple *obj = makeLeafLinkedObjSimple();
  tail -> fwdLink = obj;
  obj -> bakLink = tail;
  obj -> nodePtr = nodePtr;
  return obj;
}
void freeLeafLinkedObj(LeafLinkedObj *obj) {
  if (obj -> termPtr != NULL) {
    freeTerminal(obj -> termPtr);
    obj -> termPtr = NULL;
  }
  free_gblock(obj, (size_t) sizeof(LeafLinkedObj));
}
void freeLeafLinkedObjSimple(LeafLinkedObjSimple *obj) {
  free_gblock(obj, (size_t) sizeof(LeafLinkedObjSimple));
}
void freeLeafLinkedObjList(LeafLinkedObj *obj) {
  if (obj -> fwdLink != NULL) {
    freeLeafLinkedObjList(obj -> fwdLink);
  }
  freeLeafLinkedObj(obj);
}
void freeLeafLinkedObjListRev(LeafLinkedObj *obj) {
  if (obj -> bakLink != NULL) {
    freeLeafLinkedObjListRev(obj -> bakLink);
  }
  freeLeafLinkedObj(obj);
}
