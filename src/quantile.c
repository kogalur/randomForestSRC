
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "quantile.h"
#include "nrutil.h"
#include "error.h"
QuantileObj *makeQuantileObj(double value) {
  QuantileObj *quantileObj = (QuantileObj*) gblock((size_t) sizeof(QuantileObj));
  quantileObj -> fwdLink = NULL;
  quantileObj -> bakLink = NULL;
  quantileObj -> v = value;
  quantileObj -> g = 1;
  quantileObj -> dlt = 0;
  return quantileObj;
}
void freeQuantileObj(QuantileObj *obj) {
  free_gblock(obj, (size_t) sizeof(QuantileObj));
}
void freeQuantileObjList(QuantileObj *obj) {
  QuantileObj *thisObj, *fwdObj;
  thisObj = obj;
  while (thisObj != NULL) {
    fwdObj = thisObj -> fwdLink;
    freeQuantileObj(thisObj);
    thisObj = fwdObj;
  }
}
QuantileObj *insertQuantileObj(uint *qStreamSize, QuantileObj **head, QuantileObj **tail, uint *quantileLinkLength, double value, LookUpInfo **tree) {
  QuantileObj *newObj;
  QuantileObj *insertPtr;
  QuantileObj *thisPtr;
  QuantileObj *segmentHead, *segmentTail, *delPtr, *savPtr;
  uint *band;
  uint p, gStar;
  uint gNew;
  char flag;
  newObj = makeQuantileObj(value);
  if (*head == NULL) {
    *head = *tail = newObj;
    (*quantileLinkLength) ++;
    (*qStreamSize) ++;
  }
  else {
    if ( (((*qStreamSize) % ((uint) floor(RF_inv_2qEpsilon))) == 0) &&
         ((*qStreamSize) > (uint) floor(RF_inv_2qEpsilon)) && 
         ((*qStreamSize) > 2)) {
      p = (*qStreamSize) / (uint) floor(RF_inv_2qEpsilon);
      band = uivector(0, p);
      populateBand(p, band);
      thisPtr = *tail;
      while(thisPtr != *head) {
        segmentTail = thisPtr -> bakLink; 
        if (segmentTail != *head) {
          if (band[segmentTail -> dlt] <= band[thisPtr -> dlt]) {
            gStar = 0;
            segmentHead = segmentTail;
            flag = TRUE;
            while (flag && (segmentHead != (*head))) {
              gStar += (segmentHead -> g);
              segmentHead = segmentHead -> bakLink;
              if ((band[segmentHead -> dlt] < band[segmentTail -> dlt]) && segmentHead != (*head)) {
                flag = TRUE;
              }
              else {
                flag = FALSE;
              }
            }
            gNew = gStar + (thisPtr -> g);
            if (gNew + (thisPtr -> dlt) <= p) {
              delPtr = segmentHead -> fwdLink;
              segmentHead -> fwdLink = thisPtr;
              thisPtr -> bakLink = segmentHead;
              while (delPtr != thisPtr) {
                savPtr = delPtr -> fwdLink;
                freeQuantileObj(delPtr);
                delPtr = savPtr;
                (*quantileLinkLength) --;
              }
              thisPtr -> g = gNew;
              thisPtr = segmentHead;
            }
            else {
              thisPtr = segmentHead;
            }
          }
          else {
            thisPtr = thisPtr -> bakLink;
          }
        }
        else {
          thisPtr = thisPtr -> bakLink;
        }
      }
      free_uivector(band, 0, p);
      if (*tree != NULL) {
        freeLookUpTree(*tree);
        *tree = NULL;
      }
      if (*quantileLinkLength >= 8) {
        *tree = makeLookUpInfo();
        uint depth = ulog2(*quantileLinkLength) - 2;
        makeLookUpTree(*tree, *head, *quantileLinkLength, depth);
      }
      else {
      }
    }
    else {
    }
    if (value <= (*head) -> v) {
      (*head) -> bakLink = newObj;
      newObj -> fwdLink = *head;
      *head = newObj;
      newObj -> g = 1;
      newObj -> dlt = 0;
      (*quantileLinkLength) ++;
      (*qStreamSize) ++;
    }
    else if (value >= (*tail) -> v) {
      (*tail) -> fwdLink = newObj;
      newObj -> bakLink = *tail;
      *tail = newObj;
      newObj -> g = 1;
      newObj -> dlt = 0;
      (*quantileLinkLength) ++;
      (*qStreamSize) ++;
    }
    else {
      insertPtr = findInsertionPoint(*head, value, *tree);
      (insertPtr -> bakLink) -> fwdLink = newObj;
      newObj -> bakLink = insertPtr -> bakLink;
      insertPtr -> bakLink = newObj;
      newObj -> fwdLink = insertPtr;
      newObj -> g = 1;
      if ((double) *qStreamSize <= RF_inv_2qEpsilon) {
        newObj -> dlt = 0;
      }
      else {
        newObj -> dlt = (insertPtr -> g) + (insertPtr -> dlt) - 1;
      }
      (*quantileLinkLength) ++;
      (*qStreamSize) ++;
    }
  }
  return newObj;
}
QuantileObj *findInsertionPoint(QuantileObj *head, double value, LookUpInfo *tree) {
  QuantileObj *insertPtr;
  char found;
  found = FALSE;
  if (tree == NULL) {
    insertPtr = head;
  }
  else {
    findApproximateInsertionPoint(head, tree, value, &insertPtr);
  }
  while (!found) {
    if (insertPtr != NULL) {
      if (value > insertPtr -> v) {
        insertPtr = insertPtr -> fwdLink;
      }
      else {
        found = TRUE;
      }
    }
    else {
      insertPtr = NULL;
    }
  }   
  return insertPtr;
}
double getApproxQuantile(QuantileObj *head, double phi, uint streamSize) {
  double rank, margin;
  QuantileObj *currentObj;
  double rmax, rmin;
  double result;
  char found;
  rank = ceil(phi * streamSize);
  margin = RF_qEpsilon * streamSize;
  found = FALSE;
  currentObj = head;
  result = RF_nativeNaN; 
  rmin = 0;
  while (!found) {
    if (currentObj != NULL) {
      rmin += (double) (currentObj -> g);
      rmax = rmin + (double) (currentObj -> dlt);
      if ( (((rank - rmin) <= margin) && ((rmax - rank) <= margin)) ) {
        found = TRUE;
      }
      else if ((uint) rmin == streamSize) {
        found = TRUE;
      }
      else {
        currentObj = currentObj -> fwdLink;
      }
    }
    else {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  Quantile query failed with (epsilon, phi) -> (margin, rank) => (%10.4f, %10.4f) -> (%10.4f, %10d)", RF_qEpsilon, phi, margin, (uint) rank);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  result = currentObj -> v;
  return result;
}
void populateBand(uint p, uint *band) {
  uint alpha, alphaPowLo, alphaPowHi;
  uint j;
  uint lower, upper;
  uint alphaLimit = ulog2(p);
  band[0]   = INT_MAX;
  band[p] = 0;
  for (alpha = 1; alpha <= alphaLimit; alpha ++) {
    alphaPowLo = 1 << (alpha - 1);
    alphaPowHi = 1 << alpha;
    lower = p - alphaPowHi - (p % alphaPowHi);
    upper = p - alphaPowLo - (p % alphaPowLo);
    for (j = upper; j > lower; j--) {
      band[j] = alpha;
    }
  }
}
void makeLookUpTree(LookUpInfo *infoObj, QuantileObj *qObj, uint size, uint depth) {
  QuantileObj *qPtr;
  uint half;
  uint i;
  half = (size >> 1);
  qPtr = qObj;
  for (i = 1; i < half; i++) {
    qPtr = qPtr -> fwdLink;
  }
  infoObj -> qPtr = qPtr;
  if (half > 1) {
    if (depth > 1) {      
      LookUpInfo *leftPtr = makeLookUpInfo();
      infoObj -> leftPtr = leftPtr;
      leftPtr -> rootPtr = infoObj;
      makeLookUpTree(leftPtr, qObj, half, depth - 1);
      LookUpInfo *rghtPtr = makeLookUpInfo();
      infoObj -> rghtPtr = rghtPtr;
      rghtPtr -> rootPtr = infoObj;
      makeLookUpTree(rghtPtr, qPtr, size - half, depth - 1);
    }
    else {
    }
  }
  else {
  }
}
LookUpInfo *makeLookUpInfo(void) {
  LookUpInfo *obj = (LookUpInfo*) gblock((size_t) sizeof(LookUpInfo));
  obj -> qPtr    = NULL;
  obj -> rootPtr = NULL;
  obj -> leftPtr = NULL;
  obj -> rghtPtr = NULL;
  return obj;
}
void freeLookUpInfo(LookUpInfo *obj) {
  free_gblock(obj, (size_t) sizeof(LookUpInfo));
}
void freeLookUpTree(LookUpInfo *obj) {
  if (obj != NULL) {
    if ((obj -> leftPtr != NULL) && (obj -> rghtPtr != NULL)) {
      freeLookUpTree(obj -> leftPtr);
      freeLookUpTree(obj -> rghtPtr);
    }
    freeLookUpInfo(obj);
  }
}
void findApproximateInsertionPoint(QuantileObj *head, LookUpInfo *tree, double value, QuantileObj **insertPtr) {
  char foundFlag;
  if (value < (tree -> qPtr) -> v) {
    if (tree -> leftPtr != NULL) {
      findApproximateInsertionPoint(head, tree -> leftPtr, value, insertPtr);
    }
    else {
      foundFlag = FALSE;
      while(!foundFlag) {
        tree = tree -> rootPtr;
        if (tree != NULL) {
          if (value < ((tree -> qPtr) -> v)) {
          }
          else {
            foundFlag = TRUE;
            *insertPtr = tree -> qPtr;
          }
        }
        else {
          foundFlag = TRUE;
          *insertPtr = head;
        }
      }
    }
  }
  else if (value > (tree -> qPtr) -> v) {
    if (tree -> rghtPtr != NULL) {
      findApproximateInsertionPoint(head, tree -> rghtPtr, value, insertPtr);
    }
    else {
      *insertPtr = tree -> qPtr;
    }
  }
  else {
    *insertPtr = tree -> qPtr;
  }
}
void testQuantile(uint treeID) {
  QuantileObj *head, *tail;
  uint streamSize;
  uint quantileLinkLength;
  LookUpInfo *ghiPtr;
  head = tail = NULL;
  streamSize = 0;
  quantileLinkLength = 0;
  uint size = RF_observationSize;
  ghiPtr = NULL;
  for (uint i = 1; i <= size; i++) {
    insertQuantileObj(&streamSize, &head, &tail, &quantileLinkLength,  RF_response[treeID][1][i], &ghiPtr); 
  }
  if (!FALSE) {
    for (uint i = 1; i <= RF_quantileSize; i++) {
      getApproxQuantile(head, RF_quantile[i], streamSize);
    }
  }
  if (ghiPtr != NULL) {
    freeLookUpTree(ghiPtr);
  }
  freeQuantileObjList(head);
}
