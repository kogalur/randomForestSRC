
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "quantile.h"
#include "nrutil.h"
#include "error.h"
QuantileObj *makeQuantileObj(double value) {
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nmakeQuantileObj() ENTRY ...\n");
  ${trace.token}  }
  QuantileObj *quantileObj = (QuantileObj*) gblock((size_t) sizeof(QuantileObj));
  quantileObj -> fwdLink = NULL;
  quantileObj -> bakLink = NULL;
  quantileObj -> v = value;
  quantileObj -> g = 1;
  quantileObj -> dlt = 0;
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nmakeQuantileObj() EXIT ...\n");
  ${trace.token}  }
  return quantileObj;
}
void freeQuantileObj(QuantileObj *obj) {
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfreeQuantileObj() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfreeing Quantile Object at %20x \n", obj);
  ${trace.token}  }
  free_gblock(obj, (size_t) sizeof(QuantileObj));
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfreeQuantileObj() EXIT ...\n");
  ${trace.token}  }
}
void freeQuantileObjList(QuantileObj *obj) {
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfreeQuantileObjList() ENTRY ...\n");
  ${trace.token}  }
  QuantileObj *thisObj, *fwdObj;
  thisObj = obj;
  while (thisObj != NULL) {
    fwdObj = thisObj -> fwdLink;
    freeQuantileObj(thisObj);
    thisObj = fwdObj;
  }
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfreeQuantileObjList() EXIT ...\n");
  ${trace.token}  }
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
  ${trace.token}  uint s = 0;
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\ninsertQuantileObj() ENTRY ...\n");
  ${trace.token}  }
  newObj = makeQuantileObj(value);
  if (*head == NULL) {
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\n  Insert Quantile Object:  null list encountered, creating object, init head and tail for value:  %10.4f \n", value);
  ${trace.token}  }
    *head = *tail = newObj;
    (*quantileLinkLength) ++;
    (*qStreamSize) ++;
  }
  else {
    if ( (((*qStreamSize) % ((uint) floor(RF_inv_2qEpsilon))) == 0) &&
         ((*qStreamSize) > (uint) floor(RF_inv_2qEpsilon)) && 
         ((*qStreamSize) > 2)) {
      p = (*qStreamSize) / (uint) floor(RF_inv_2qEpsilon);
      ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}    s = *quantileLinkLength;
      ${trace.token}  }
      ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}    RF_nativePrint("\n  Compress Phase:     qStreamSize = %10d", *qStreamSize);
      ${trace.token}    RF_nativePrint("\n  Compress Phase: 2 x eps x n = p = %10d", p);
      ${trace.token}    RF_nativePrint("\n  Compress Phase:        s (iter) = %10d", s);      
      ${trace.token}  }
      ${trace.token}    if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}      QuantileObj *objPtr;
      ${trace.token}      uint sum;
      ${trace.token}      objPtr = *head;
      ${trace.token}      sum = 0;      
      ${trace.token}      RF_nativePrint("\n  Compress Phase:     incoming fwd parse linked list contents");
      ${trace.token}      RF_nativePrint("\n      index                               g      delta      value");
      ${trace.token}      for (uint k = 1; k <= *quantileLinkLength; k++) {
      ${trace.token}        sum += objPtr -> g;
      ${trace.token}        RF_nativePrint("\n %10d %20x %10d %10d %10.4f ", k, objPtr, objPtr -> g, objPtr -> dlt, objPtr -> v);
      ${trace.token}        objPtr = objPtr -> fwdLink;
      ${trace.token}      }
      ${trace.token}      RF_nativePrint("\n  Compress Phase:     consistency checksum_{i=1}^{n} g_i = %10d", sum);
      ${trace.token}    }
      ${trace.token}    if (getTraceFlag(0) & QUAN_DEF_TRACE) {      
      ${trace.token}    if (getTraceFlag(0) & TURN_OFF_TRACE) {
      ${trace.token}      QuantileObj *objPtr;
      ${trace.token}      objPtr = *tail;
      ${trace.token}      RF_nativePrint("\n  Compress Phase:     incoming bak parse linked list contents");
      ${trace.token}      RF_nativePrint("\n      index                               g      delta      value");
      ${trace.token}      for (uint k = *quantileLinkLength; k >= 1; k--) {
      ${trace.token}        RF_nativePrint("\n %10d %20x %10d %10d %10.4f ", k, objPtr, objPtr -> g, objPtr -> dlt, objPtr -> v);
      ${trace.token}        objPtr = objPtr -> bakLink;
      ${trace.token}      }
      ${trace.token}    }
      ${trace.token}    }
      band = uivector(0, p);
      populateBand(p, band);
      ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}    RF_nativePrint("\n Mapping of deltas to band identifiers:");
      ${trace.token}    for (uint j = 0; j <= p; j++) {
      ${trace.token}      RF_nativePrint("\n   band [%10d] = %20d", j, band[j]);  
      ${trace.token}    }
      ${trace.token}  }
      thisPtr = *tail;
      while(thisPtr != *head) {
        ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
        ${trace.token}    RF_nativePrint("\n  Compress Quantile Object:  examining thisPtr = s (addr, iter) = (%20x, %10d)", thisPtr, s);
        ${trace.token}  }
        segmentTail = thisPtr -> bakLink; 
        if (segmentTail != *head) {
          if (band[segmentTail -> dlt] <= band[thisPtr -> dlt]) {
            ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
            ${trace.token}    RF_nativePrint("\n  Compress Quantile Object:  adjacent tuples suitable (band[L], band[R]) = (%20d, %20d)", band[segmentTail -> dlt], band[thisPtr -> dlt]);
            ${trace.token}  }
            gStar = 0;
            segmentHead = segmentTail;
            flag = TRUE;
            while (flag && (segmentHead != (*head))) {
              gStar += (segmentHead -> g);
              ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
              ${trace.token}    RF_nativePrint("\n  Compress Quantile Object:  gStar = %10d", gStar);
              ${trace.token}  }
              segmentHead = segmentHead -> bakLink;
              ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
              ${trace.token}    s --;
              ${trace.token}  }
              ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
              ${trace.token}    RF_nativePrint("\n  Compress Quantile Object:  examining segmentHead = s (addr, iter) = (%20x, %10d)", segmentHead, s - 1);
              ${trace.token}  }
              if ((band[segmentHead -> dlt] < band[segmentTail -> dlt]) && segmentHead != (*head)) {
                flag = TRUE;
                ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
                ${trace.token}    RF_nativePrint("\n  Compress Quantile Object:  further descendant found, continue search ...");
                ${trace.token}  }
              }
              else {
                flag = FALSE;
                ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
                ${trace.token}    RF_nativePrint("\n  Compress Quantile Object:  further descendant NOT found, terminating search ...");
                ${trace.token}  }
                ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
                ${trace.token}    s --;
                ${trace.token}  }
              }
            }
            gNew = gStar + (thisPtr -> g);
            if (gNew + (thisPtr -> dlt) <= p) {
              ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
              ${trace.token}    RF_nativePrint("\n  Compress Quantile Object:  successful epsilon-approx guarantree (gStar + g_i + delta_i, p) = (%10d, %10d)", gNew + (thisPtr -> dlt), p);
              ${trace.token}  }
              delPtr = segmentHead -> fwdLink;
              segmentHead -> fwdLink = thisPtr;
              thisPtr -> bakLink = segmentHead;
              while (delPtr != thisPtr) {
                ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
                ${trace.token}    RF_nativePrint("\n  Compress Quantile Object:  deleting tuple (g, delta, v) = (%10d, %10d, %10.4f)", delPtr -> g, delPtr -> dlt, delPtr -> v);
                ${trace.token}  }
                savPtr = delPtr -> fwdLink;
                freeQuantileObj(delPtr);
                delPtr = savPtr;
                (*quantileLinkLength) --;
              }
              ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
              ${trace.token}    RF_nativePrint("\n  Compress Quantile Object:  receiving tuple (g, delta, v) = (%10d, %10d, %10.4f)", thisPtr -> g, thisPtr -> dlt, thisPtr -> v);
              ${trace.token}  }
              thisPtr -> g = gNew;
              ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
              ${trace.token}    RF_nativePrint("\n  Compress Quantile Object:  resultant tuple (g, delta, v) = (%10d, %10d, %10.4f)", thisPtr -> g, thisPtr -> dlt, thisPtr -> v);
              ${trace.token}  }
              ${trace.token}    if (getTraceFlag(0) & QUAN_DEF_TRACE) {
              ${trace.token}      QuantileObj *objPtr;
              ${trace.token}      uint sum;
              ${trace.token}      objPtr = *head;
              ${trace.token}      sum = 0;      
              ${trace.token}      for (uint k = 1; k <= *quantileLinkLength; k++) {
              ${trace.token}        sum += objPtr -> g;
              ${trace.token}        objPtr = objPtr -> fwdLink;
              ${trace.token}      }
              ${trace.token}      RF_nativePrint("\n  Compress Quantile Object:  consistency checksum_{i=1}^{n} g_i = %10d", sum);
              ${trace.token}      if (sum != *qStreamSize) {
              ${trace.token}        RF_nativeError("\nRF-SRC:  *** ERROR *** ");
              ${trace.token}        RF_nativeError("\nRF-SRC:  Quantile consistency checksum failed with (actual, expected) = (%10.d, %10.d)", sum, *qStreamSize);
              ${trace.token}        RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
              ${trace.token}        RF_nativeExit();
              ${trace.token}      }
              ${trace.token}    }
              thisPtr = segmentHead;
            }
            else {
              ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
              ${trace.token}    RF_nativePrint("\n  Compress Quantile Object:  failing epsilon-approx guarantree (gNew, p) = (%10d, %10d)", gNew, p);
              ${trace.token}    RF_nativePrint("\n  Compress Quantile Object:  terminating iteration");
              ${trace.token}  }
              thisPtr = segmentHead;
            }
          }
          else {
            thisPtr = thisPtr -> bakLink;
            ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
            ${trace.token}    RF_nativePrint("\n  Compress Quantile Object:  adjacent tuples not capacity compatible");
            ${trace.token}  }
            ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
            ${trace.token}    s --;
            ${trace.token}  }
          }
        }
        else {
          thisPtr = thisPtr -> bakLink;
          ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
          ${trace.token}    RF_nativePrint("\n  Compress Quantile Object:  segment tail starts at head");
          ${trace.token}    RF_nativePrint("\n  Compress Quantile Object:  terminating iteration");
          ${trace.token}  }
        }
      }
      free_uivector(band, 0, p);
      if (*tree != NULL) {
        ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
        ${trace.token}    RF_nativePrint("\n  LookUp Tree:  freeing current tree");
        ${trace.token}  }
        freeLookUpTree(*tree);
        *tree = NULL;
      }
      if (*quantileLinkLength >= 8) {
        ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
        ${trace.token}    RF_nativePrint("\n  LookUp Tree:  creating and populating new search tree.");
        ${trace.token}  }
        *tree = makeLookUpInfo();
        uint depth = ulog2(*quantileLinkLength) - 2;
        makeLookUpTree(*tree, *head, *quantileLinkLength, depth);
      }
      else {
        ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
        ${trace.token}    RF_nativePrint("\n  LookUp Tree:  search tree omitted, link length is %10d", *quantileLinkLength);
        ${trace.token}  }
      }
    }
    else {
      ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}    RF_nativePrint("\n  Compress Quantile Object:  OMITTED:     2 x epsilon = %10d, qStreamsize = %10d", (uint) floor(RF_inv_2qEpsilon), *qStreamSize);
      ${trace.token}  }
    }
    if (value <= (*head) -> v) {
      ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}    RF_nativePrint("\n  Insert Quantile Object:  creating new head for value:  %10.4f \n", value);
      ${trace.token}  }
      (*head) -> bakLink = newObj;
      newObj -> fwdLink = *head;
      *head = newObj;
      newObj -> g = 1;
      newObj -> dlt = 0;
      (*quantileLinkLength) ++;
      (*qStreamSize) ++;
    }
    else if (value >= (*tail) -> v) {
      ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}    RF_nativePrint("\n  Insert Quantile Object:  creating new tail for value:  %10.4f \n", value);
      ${trace.token}  }
      (*tail) -> fwdLink = newObj;
      newObj -> bakLink = *tail;
      *tail = newObj;
      newObj -> g = 1;
      newObj -> dlt = 0;
      (*quantileLinkLength) ++;
      (*qStreamSize) ++;
    }
    else {
      ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}    RF_nativePrint("\n  Insert Quantile Object:  finding insertion point for value:  %10.4f \n", value);
      ${trace.token}  }
      insertPtr = findInsertionPoint(*head, value, *tree);
      ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}    RF_nativePrint("\n  Insert Quantile Object:  insertion point found at:  %20x, value =  %10.4f \n", insertPtr -> v);
      ${trace.token}  }
      (insertPtr -> bakLink) -> fwdLink = newObj;
      newObj -> bakLink = insertPtr -> bakLink;
      ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}    RF_nativePrint("\n  Insert Quantile Object:  backward linking complete");
      ${trace.token}  }
      insertPtr -> bakLink = newObj;
      newObj -> fwdLink = insertPtr;
      ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}    RF_nativePrint("\n  Insert Quantile Object:  forward  linking complete");
      ${trace.token}  }
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
      ${trace.token}    if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}      QuantileObj *objPtr;
      ${trace.token}      uint sum;
      ${trace.token}      objPtr = *head;
      ${trace.token}      sum = 0;      
      ${trace.token}      RF_nativePrint("\n  Insert Phase:     outgoing fwd parse linked list contents");
      ${trace.token}      RF_nativePrint("\n      index                               g      delta      value");
      ${trace.token}      for (uint k = 1; k <= *quantileLinkLength; k++) {
      ${trace.token}        sum += objPtr -> g;
      ${trace.token}        RF_nativePrint("\n %10d %20x %10d %10d %10.4f ", k, objPtr, objPtr -> g, objPtr -> dlt, objPtr -> v);
      ${trace.token}        objPtr = objPtr -> fwdLink;
      ${trace.token}      }
      ${trace.token}      RF_nativePrint("\n  Insert Quantile Object:     consistency checksum_{i=1}^{n} g_i = %10d", sum);
      ${trace.token}    }
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\ninsertQuantileObj() EXIT ...\n");
  ${trace.token}  }
  return newObj;
}
QuantileObj *findInsertionPoint(QuantileObj *head, double value, LookUpInfo *tree) {
  QuantileObj *insertPtr;
  char found;
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfindInsertionPoint() ENTRY ...\n");
  ${trace.token}  }
  found = FALSE;
  if (tree == NULL) {
    insertPtr = head;
  }
  else {
    findApproximateInsertionPoint(head, tree, value, &insertPtr);
    ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
    ${trace.token}    RF_nativePrint("\n Approximate insertion point:  obj = %20x, Vapprox = %10.4f, Vactual = %10.4f", insertPtr, insertPtr -> v, value);
    ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfindInsertionPoint() EXIT ...\n");
  ${trace.token}  }
  return insertPtr;
}
double getApproxQuantile(QuantileObj *head, double phi, uint streamSize) {
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetApproxQuantile() ENTRY ...\n");
  ${trace.token}  }
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
      ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}    RF_nativePrint("\n approximation details: ");
      ${trace.token}    RF_nativePrint("\n    margin = %10.4f", margin);
      ${trace.token}    RF_nativePrint("\n    rank   = %10.4f", rank);
      ${trace.token}    RF_nativePrint("\n    (rmin, rmax) = (%10.4f, %10.4f)", rmin, rmax);
      ${trace.token}    RF_nativePrint("\n    ((rank - rmin), (rmax - rank) = (%10.4f, %10.4f)", rank - rmin, rmax - rank);
      ${trace.token}  }
      if ( (((rank - rmin) <= margin) && ((rmax - rank) <= margin)) ) {
        found = TRUE;
      ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}    RF_nativePrint("\n coverage found");
      ${trace.token}  }
      }
      else if ((uint) rmin == streamSize) {
        found = TRUE;
        ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
        ${trace.token}    RF_nativePrint("\n rmin = mrax = streamSize, found");
        ${trace.token}  }
      }
      else {
        currentObj = currentObj -> fwdLink;
        ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
        ${trace.token}    RF_nativePrint("\n coverage NOT found, iterating forward");
        ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\n approximation with (epsilon, phi) -> (margin, rank) => (%10.4f, %10.4f) -> (%10.4f, %10d)", RF_qEpsilon, phi, margin, (uint) rank);
  ${trace.token}    RF_nativePrint("\n   quantile approximation = %10.4f", result);
  ${trace.token}  }
  return result;
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetApproxQuantile() EXIT ...\n");
  ${trace.token}  }
}
void populateBand(uint p, uint *band) {
  uint alpha, alphaPowLo, alphaPowHi;
  uint j;
  uint lower, upper;
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\npopulateBand() ENTRY ...\n");
  ${trace.token}  }
  uint alphaLimit = ulog2(p);
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nBand Maximum alhpa:  %10d", alphaLimit);
  ${trace.token}  }
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
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\npopulateBand() ENTRY ...\n");
  ${trace.token}  }
}
void makeLookUpTree(LookUpInfo *infoObj, QuantileObj *qObj, uint size, uint depth) {
  QuantileObj *qPtr;
  uint half;
  uint i;
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nmakeLookUpTree() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\n  makeLookUpTree() called with:  GKLI = %20x, QO = %20x, size = %10d, depth = %10d", infoObj, qObj, size, depth);
  ${trace.token}  }
  half = (size >> 1);
  qPtr = qObj;
  for (i = 1; i < half; i++) {
    qPtr = qPtr -> fwdLink;
  }
  infoObj -> qPtr = qPtr;
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\n  G-K LookUp Object initialization:  half = %10d, QO = %20x, v = %10.4f", half, qPtr, qPtr -> v);
  ${trace.token}  }
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
    ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
    ${trace.token}    RF_nativePrint("\n  G-K LookUp Object tree depth limit reached:  depth = %10d", depth);
    ${trace.token}  }
    }
  }
  else {
    ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
    ${trace.token}    RF_nativePrint("\n  G-K LookUp Object segment limit reached:  half = %10d", half);
    ${trace.token}  }
  }
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nmakeLookUpTree() EXIT ...\n");
  ${trace.token}  }
}
LookUpInfo *makeLookUpInfo(void) {
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nmakeLookUpInfo() ENTRY ...\n");
  ${trace.token}  }
  LookUpInfo *obj = (LookUpInfo*) gblock((size_t) sizeof(LookUpInfo));
  obj -> qPtr    = NULL;
  obj -> rootPtr = NULL;
  obj -> leftPtr = NULL;
  obj -> rghtPtr = NULL;
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nmakeLookUpInfo() EXIT ...\n");
  ${trace.token}  }
  return obj;
}
void freeLookUpInfo(LookUpInfo *obj) {
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfreeLookUpInfo() ENTRY ...\n");
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfreeing LookUpInfo object at %20x \n", obj);
  ${trace.token}  }
  free_gblock(obj, (size_t) sizeof(LookUpInfo));
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfreeLookUpInfo() EXIT ...\n");
  ${trace.token}  }
}
void freeLookUpTree(LookUpInfo *obj) {
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfreeLookUpTree() ENTRY ...\n");
  ${trace.token}  }
  if (obj != NULL) {
    if ((obj -> leftPtr != NULL) && (obj -> rghtPtr != NULL)) {
      freeLookUpTree(obj -> leftPtr);
      freeLookUpTree(obj -> rghtPtr);
    }
    freeLookUpInfo(obj);
  }
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfreeLookUpTree() EXIT ...\n");
  ${trace.token}  }
}
void findApproximateInsertionPoint(QuantileObj *head, LookUpInfo *tree, double value, QuantileObj **insertPtr) {
  char foundFlag;
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfindApproximateInsertionPoint() ENTRY ...\n");
  ${trace.token}  }
  if (value < (tree -> qPtr) -> v) {
    if (tree -> leftPtr != NULL) {
      ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}    RF_nativePrint("\n  finding approximate insertion:  parsing left with (value, mid)  = (%10.4f, %10.4f)", value, (tree -> qPtr) -> v);
      ${trace.token}  }
      findApproximateInsertionPoint(head, tree -> leftPtr, value, insertPtr);
    }
    else {
      foundFlag = FALSE;
      while(!foundFlag) {
        tree = tree -> rootPtr;
        if (tree != NULL) {
          if (value < ((tree -> qPtr) -> v)) {
            ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
            ${trace.token}    RF_nativePrint("\n  finding approximate insertion:  parsing upward with (value, mid)  = (%10.4f, %10.4f)", value, (tree -> qPtr) -> v);
            ${trace.token}  }
          }
          else {
            foundFlag = TRUE;
            *insertPtr = tree -> qPtr;
            ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
            ${trace.token}    RF_nativePrint("\n  finding approximate insertion:  parsing stopped nominal with (value, mid)  = (%10.4f, %10.4f)", value, (tree -> qPtr) -> v);
            ${trace.token}  }
          }
        }
        else {
          foundFlag = TRUE;
          *insertPtr = head;
          ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
          ${trace.token}    RF_nativePrint("\n  finding approximate insertion:  parsing stopped root with (value, mid)  = (%10.4f, %10.4f)", value, head -> v);
          ${trace.token}  }
        }
      }
    }
  }
  else if (value > (tree -> qPtr) -> v) {
    if (tree -> rghtPtr != NULL) {
      ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}    RF_nativePrint("\n  finding approximate insertion:  parsing right with (value, mid)  = (%10.4f, %10.4f)", value, (tree -> qPtr) -> v);
      ${trace.token}  }
      findApproximateInsertionPoint(head, tree -> rghtPtr, value, insertPtr);
    }
    else {
      *insertPtr = tree -> qPtr;
      ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
      ${trace.token}    RF_nativePrint("\n  finding approximate insertion:  parsing stopped terminal with (value, mid)  = (%10.4f, %10.4f)", value, (tree -> qPtr) -> v);
      ${trace.token}  }
    }
  }
  else {
    *insertPtr = tree -> qPtr;
    ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
    ${trace.token}    RF_nativePrint("\n  finding approximate insertion:  parsing stopped exact with (value, mid)  = (%10.4f, %10.4f)", value, (tree -> qPtr) -> v);
    ${trace.token}  }
  }
  ${trace.token}  if (getTraceFlag(0) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\nfindApproximateInsertionPoint() EXIT ...\n");
  ${trace.token}  }
}
void testQuantile(uint treeID) {
  QuantileObj *head, *tail;
  uint streamSize;
  uint quantileLinkLength;
  LookUpInfo *ghiPtr;
  ${trace.token}  if (getTraceFlag(treeID) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\ntestQuantile(%10d) ENTRY ...\n", treeID);
  ${trace.token}  }
  head = tail = NULL;
  streamSize = 0;
  quantileLinkLength = 0;
  uint size = RF_observationSize;
  ghiPtr = NULL;
  for (uint i = 1; i <= size; i++) {
    insertQuantileObj(&streamSize, &head, &tail, &quantileLinkLength,  RF_response[treeID][1][i], &ghiPtr); 
  }
  ${trace.token}    if (getTraceFlag(treeID) & !TURN_OFF_TRACE) {
  ${trace.token}      QuantileObj *objPtr = head;
  ${trace.token}      RF_nativePrint("\nFinal Quantile Linked List:  ");
  ${trace.token}      RF_nativePrint("\n      index          g      delta      value");
  ${trace.token}      for (uint k = 1; k <= quantileLinkLength; k++) {
  ${trace.token}        RF_nativePrint("\n %10d %10d %10d %10.4f ", k, objPtr -> g, objPtr -> dlt, objPtr -> v);
  ${trace.token}        objPtr = objPtr -> fwdLink;
  ${trace.token}      }
  ${trace.token}    }
  if (!FALSE) {
    for (uint i = 1; i <= RF_quantileSize; i++) {
      getApproxQuantile(head, RF_quantile[i], streamSize);
    }
  }
  if (ghiPtr != NULL) {
    freeLookUpTree(ghiPtr);
  }
  freeQuantileObjList(head);
  ${trace.token}  if (getTraceFlag(treeID) & QUAN_DEF_TRACE) {
  ${trace.token}    RF_nativePrint("\ntestQuantile(%10d) EXIT ...\n", treeID);
  ${trace.token}  }
}
