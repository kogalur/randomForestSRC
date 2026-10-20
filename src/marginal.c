
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "marginal.h"
#include "polarity.h"
#include "nrutil.h"
void getMarginalMembership(char mode, uint treeID) {
  double **observationPtr;
  uint    *gallMembrIndx;
  uint     obsSize;
  uint i;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    observationPtr = RF_fobservation[treeID];
    gallMembrIndx  = RF_fidentityMembershipIndex;
    break;
  default:
    obsSize = RF_observationSize;
    observationPtr = RF_observation[treeID];
    gallMembrIndx  = RF_identityMembershipIndex;
    break;
  }
  RF_utTermMembership[treeID] =       (uint **) new_vvector(1, obsSize, NRUTIL_UPTR);
  RF_utTermMembershipCount[treeID] =   (uint *) uivector(1, obsSize);
  RF_utTermMembershipAlloc[treeID] =   (uint *) uivector(1, obsSize);
  for (i = 1; i <= obsSize; i++) {
    RF_utTermMembership[treeID][i] = uivector(1, MARGINAL_SIZE);
    RF_utTermMembershipCount[treeID][i] = 0;
    RF_utTermMembershipAlloc[treeID][i] = 1;            
  }
  marginalMembership(treeID,
                     RF_root[treeID],
                     gallMembrIndx,
                     obsSize,
                     obsSize,
                     observationPtr);
}
void releaseMarginalMembership(char mode, uint treeID) {
  uint     obsSize;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    break;
  default:
    obsSize = RF_observationSize;
    break;
  }
  free_uivector(RF_utTermMembershipCount[treeID], 1, obsSize);
  for (uint i = 1; i <= obsSize; i++) {
    free_uivector(RF_utTermMembership[treeID][i], 1, MARGINAL_SIZE * RF_utTermMembershipAlloc[treeID][i]);
  }
  free_new_vvector(RF_utTermMembership[treeID], 1, obsSize, NRUTIL_UPTR);
  free_uivector(RF_utTermMembershipAlloc[treeID], 1, obsSize);
}
void marginalMembership(uint     treeID,
                        Node    *parent,
                        uint    *genAllMembrIndx,
                        uint     genAllMembrSize,
                        uint     obsSize,
                        double **xArray) {
  char terminalFlag;
  uint *leftAllMembrIndx;
  uint *rghtAllMembrIndx;
  uint leftAllMembrSize;
  uint rghtAllMembrSize;
  uint jLeft;
  uint jRght;
  char daughterFlag;
  SplitInfo *info;
  void *obsLocal;
  char (*getDaughterPolarityGeneric) (uint       treeID,
                                      SplitInfo *info,
                                      uint       index,
                                      void      *value,
                                      ...);
  uint i, j;
  getDaughterPolarityGeneric = NULL;  
  terminalFlag = TRUE;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    info = parent -> splitInfo;
    terminalFlag = FALSE;
    leftAllMembrIndx = rghtAllMembrIndx = NULL;
    leftAllMembrSize = rghtAllMembrSize = 0;
    uint *indicator = uivector(1, obsSize);
    leftAllMembrSize = rghtAllMembrSize = 0;
    obsLocal = xArray[info -> randomVar[1]];
    if (info -> mwcpSizeAbs[1] > 0) {
      getDaughterPolarityGeneric = &getDaughterPolaritySimpleFactor;
    }
    else {
      getDaughterPolarityGeneric = &getDaughterPolaritySimpleNonFactor;
    }
    daughterFlag = RIGHT;
    if (RF_xMarginalFlag[info -> randomVar[1]] ==  TRUE) { 
      daughterFlag = NEITHER;
      for (i = 1; i <= genAllMembrSize; i++) {
        indicator[genAllMembrIndx[i]] = daughterFlag;
        leftAllMembrSize ++;
        rghtAllMembrSize ++;
      }        
    }
    if (daughterFlag != NEITHER) {
      for (i = 1; i <= genAllMembrSize; i++) {
        daughterFlag = getDaughterPolarityGeneric(treeID,
                                                  info,
                                                  genAllMembrIndx[i],
                                                  obsLocal);
        indicator[genAllMembrIndx[i]] = daughterFlag;
        if (daughterFlag == LEFT) {
          leftAllMembrSize ++;
        }
        else {
          rghtAllMembrSize ++;
        }
      }  
    }  
    leftAllMembrIndx  = uivector(1, leftAllMembrSize + 1);
    rghtAllMembrIndx  = uivector(1, rghtAllMembrSize + 1);
    jLeft = jRght = 0;
    if (daughterFlag == NEITHER) {
      for (i = 1; i <= genAllMembrSize; i++) {
        leftAllMembrIndx[++jLeft] = genAllMembrIndx[i];
        rghtAllMembrIndx[++jRght] = genAllMembrIndx[i];
      }
    }
    else {
      for (i = 1; i <= genAllMembrSize; i++) {
        if (indicator[genAllMembrIndx[i]] == LEFT) {
          leftAllMembrIndx[++jLeft] = genAllMembrIndx[i];
        }
        else {
          rghtAllMembrIndx[++jRght] = genAllMembrIndx[i];
        }
      }
    }
    free_uivector(indicator, 1, obsSize);
    marginalMembership(treeID,
                              parent -> left,
                              leftAllMembrIndx,
                              leftAllMembrSize,
                              obsSize,
                              xArray);
    marginalMembership(treeID,
                              parent -> right,
                              rghtAllMembrIndx,
                              rghtAllMembrSize,
                              obsSize,
                              xArray);
    free_uivector(leftAllMembrIndx, 1, leftAllMembrSize + 1);
    free_uivector(rghtAllMembrIndx, 1, rghtAllMembrSize + 1);
  }  
  else {
  }
  if (terminalFlag) {
    for (i = 1; i <= genAllMembrSize; i++) {
      RF_utTermMembership[treeID][genAllMembrIndx[i]][++ RF_utTermMembershipCount[treeID][genAllMembrIndx[i]]] = (parent -> nodeID);
      if ((RF_utTermMembershipCount[treeID][genAllMembrIndx[i]]) == (RF_utTermMembershipAlloc[treeID][genAllMembrIndx[i]] * MARGINAL_SIZE)) {
        RF_utTermMembershipAlloc[treeID][genAllMembrIndx[i]] ++;
        uint *utTermMembershipNew = uivector(1, RF_utTermMembershipAlloc[treeID][genAllMembrIndx[i]] * MARGINAL_SIZE);
        for (j = 1; j <= RF_utTermMembershipCount[treeID][genAllMembrIndx[i]]; j++) {
          utTermMembershipNew[j] = RF_utTermMembership[treeID][genAllMembrIndx[i]][j];
        }
        free_uivector(RF_utTermMembership[treeID][genAllMembrIndx[i]], 1, (RF_utTermMembershipAlloc[treeID][genAllMembrIndx[i]] - 1) * MARGINAL_SIZE);
        RF_utTermMembership[treeID][genAllMembrIndx[i]] = utTermMembershipNew;
      }
    }
  }  
}
