
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "sortedLink.h"
#include "nrutil.h"
SortedLinkedObj *makeSortedLinkedObj(void) {
  SortedLinkedObj *obj = (SortedLinkedObj*) gblock((size_t) sizeof(SortedLinkedObj));
  obj -> fwdLink = NULL;
  obj -> bakLink = NULL;
  obj -> rank = 0;
  obj -> indx = 0;
  return obj;
}
void makeAndSpliceSortedLinkedObj(uint treeID,
                                  SortedLinkedObj **headPtr,
                                  SortedLinkedObj **tailPtr,
                                  uint *listLength,
                                  uint rank, uint indx) {
  char flag;
  uint lowIndx, highIndx, halfIndx;
  SortedLinkedObj *objIterator;
  uint i;
  SortedLinkedObj *head = headPtr[treeID];
  SortedLinkedObj *tail = tailPtr[treeID];
  SortedLinkedObj *obj = makeSortedLinkedObj();
  obj -> rank = rank;
  obj -> indx = indx;
  obj -> fwdLink = obj -> bakLink = NULL;
  flag = TRUE;
  if (*listLength == 0) {
    head = tail = obj;
    flag = FALSE;
  }
  else if (rank >= tail -> rank) {
    tail -> fwdLink = obj;
    obj -> bakLink = tail;
    tail = obj;
    flag = FALSE;
  }
  else if (rank <= head -> rank) {
    head -> bakLink = obj;
    obj -> fwdLink = head;
    head = obj;
    flag = FALSE;
  }
  else {
    lowIndx  = 1;
    highIndx = *listLength;
    while (flag) {
      halfIndx = (uint) ((double) (highIndx + lowIndx) / 2.0); 
      objIterator = head;
      for (i = lowIndx; i < halfIndx; i++) {
        objIterator = objIterator -> fwdLink;
      }
      if (rank == head -> rank) {
        obj -> fwdLink = head;
        obj -> bakLink = (head -> bakLink);
        (head -> bakLink) -> fwdLink = obj;
        head -> bakLink = obj;
        flag = FALSE;
      }
      else if (rank == tail -> rank) {
        obj -> fwdLink = tail;
        obj -> bakLink = (tail -> bakLink);
        (tail -> bakLink) -> fwdLink = obj;
        tail -> bakLink = obj;
        flag = FALSE;
      }
      else if (rank == objIterator -> rank) {
        obj -> fwdLink = objIterator;
        obj -> bakLink = (objIterator -> bakLink);
        (objIterator -> bakLink) -> fwdLink = obj;
        objIterator -> bakLink = obj;
        flag = FALSE;
      }
      else if (halfIndx == lowIndx) {
        obj -> fwdLink = tail;
        obj -> bakLink = (tail -> bakLink);
        (tail -> bakLink) -> fwdLink = obj;
        tail -> bakLink = obj;
        flag = FALSE;
      }
      else if (rank < objIterator -> rank) {
        tail = objIterator;
        highIndx = halfIndx;
      }
      else {
        head = objIterator;
        lowIndx = halfIndx;
      }
    }
  }
  (*listLength) ++;
}
void freeSortedLinkedObjList(SortedLinkedObj *obj) {
  if (obj -> fwdLink != NULL) {
    freeSortedLinkedObjList(obj -> fwdLink);
  }
  freeSortedLinkedObj(obj);
  obj = NULL;
}
void freeSortedLinkedObj(SortedLinkedObj *obj) {
  free_gblock(obj, (size_t) sizeof(SortedLinkedObj));
}
