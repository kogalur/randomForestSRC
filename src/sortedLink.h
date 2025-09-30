#ifndef RF_SORTED_LINK_H
#define RF_SORTED_LINK_H
typedef struct sortedLinkedObj SortedLinkedObj;
struct sortedLinkedObj {
  struct sortedLinkedObj *fwdLink;
  struct sortedLinkedObj *bakLink;
  uint rank;
  uint indx;
};
SortedLinkedObj *makeSortedLinkedObj(void);
void makeAndSpliceSortedLinkedObj(uint treeID,
                                  SortedLinkedObj **headPtr,
                                  SortedLinkedObj **tailPtr,
                                  uint *listLength,
                                  uint rank, uint indx);
void freeSortedLinkedObjList(SortedLinkedObj *obj);
void freeSortedLinkedObj(SortedLinkedObj *obj);
#endif
