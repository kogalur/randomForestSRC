
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "nodeOps.h"
#include "splitGreedy.h"
#include "treeUtil.h"
#include "nrutil.h"
#include "error.h"
Node *makeNode(unsigned int xSize) {
  Node *parent = (Node*) gblock((size_t) sizeof(Node));
  if (xSize > 0) {
    parent -> xSize = xSize;
    parent -> permissible = cvector(1, xSize);
    parent -> permissibleIndx = uivector(1, xSize);
    parent -> permissibleIndxSize = xSize;
    parent -> permissibleReIndxFlag = FALSE;
    parent -> permissibleOwnershipFlag = TRUE;
  }
  else {
    parent -> xSize = 0;
    parent -> permissible = NULL;
    parent -> permissibleIndx = NULL;
    parent -> permissibleIndxSize = 0;
    parent -> permissibleReIndxFlag = FALSE;
    parent -> permissibleOwnershipFlag = FALSE;
  }
  parent -> parent = NULL;
  parent -> mate               = NULL;
  parent -> left               = NULL;
  parent -> right              = NULL;
  parent -> splitFlag            = TRUE;
  parent -> nodeID               = 0;
  parent -> blnodeID             = 0;
  parent -> brnodeID             = 0;
  parent -> fsrecID              = 0;
  parent -> splitStatistic       = RF_nativeNaN;
  parent -> variance             = RF_nativeNaN;
  parent -> mean                 = RF_nativeNaN;
  parent -> depth                = 0;
  parent -> splitDepth           = NULL;
  parent -> pseudoTerminal       = FALSE;
  parent -> mpIndexSize          = 0;
  parent -> fmpIndexSize         = 0;
  parent -> mpSign               = NULL;
  parent -> fmpSign              = NULL;
  parent -> imputed              = FALSE;
  parent -> lmpIndex             = NULL;
  parent -> flmpIndex            = NULL;
  parent -> lmpValue             = NULL;
  parent -> lmpIndexAllocSize    = 0;
  parent -> flmpIndexAllocSize   = 0;
  parent -> lmpIndexActualSize   = 0;
  parent -> flmpIndexActualSize  = 0;
  parent -> splitInfo = NULL;
  parent -> repMembrIndx = NULL;
  parent -> allMembrIndx = NULL;
  parent -> repMembrSizeAlloc = parent -> repMembrSize = 0;
  parent -> allMembrSizeAlloc = parent -> allMembrSize = 0;
  parent -> oobMembrSizeAlloc = parent -> oobMembrSize = 0;
  parent -> oobMembrIndx = NULL;
  parent -> nonMissMembrIndxStatic = NULL;
  parent -> nonMissMembrSizeStatic = 0;
  parent -> nonMissMembrIndx       = NULL;
  parent -> nonMissMembrSize       = 0;
  parent -> sumRght = 0.0;
  return parent;
}
void freeNodeGeneric(Node *parent) {
  if (parent -> xSize > 0) {
    if (parent -> permissible != NULL) {
      free_cvector(parent -> permissible, 1, parent -> xSize);
    }
    if (parent -> permissibleOwnershipFlag) {
      if (parent -> permissibleIndx != NULL) {
        free_uivector(parent -> permissibleIndx, 1, parent -> xSize);
      }
    }
    parent -> permissible = NULL;
    parent -> permissibleIndx = NULL;
    parent -> permissibleIndxSize = 0;
  }
  unstackMPSign(parent);
  unstackFMPSign(parent);
  unstackNodeLMPIndex(parent);
  unstackNodeFLMPIndex(parent);
  if (parent -> splitInfo != NULL) {
    freeSplitInfo(parent -> splitInfo);
    parent -> splitInfo = NULL;
  }
  if (parent -> repMembrSizeAlloc > 0) {
    if (parent -> repMembrIndx != NULL) {
      free_uivector(parent -> repMembrIndx, 1, parent -> repMembrSizeAlloc);
      parent -> repMembrIndx = NULL;
    }
  }
  if (parent -> allMembrSizeAlloc > 0) {
    if (parent -> allMembrIndx != NULL) {
      free_uivector(parent -> allMembrIndx, 1, parent -> allMembrSizeAlloc);
      parent -> allMembrIndx = NULL;
    }
  }
  if (parent -> oobMembrSizeAlloc > 0) {
    if (parent -> oobMembrIndx != NULL) {
      free_uivector(parent -> oobMembrIndx, 1, parent -> oobMembrSizeAlloc);
      parent -> oobMembrIndx = NULL;
    }
  }
  free_gblock(parent, (size_t) sizeof(Node));
}
void freeNodeNew(Node *parent) {
  if (parent -> xSize > 0) {
    if (parent -> permissibleOwnershipFlag) {
      if (parent -> permissible != NULL) {
        free_cvector(parent -> permissible, 1, parent -> xSize);
      }
    }
    if (parent -> permissibleOwnershipFlag) {
      if (parent -> permissibleIndx != NULL) {
        free_uivector(parent -> permissibleIndx, 1, parent -> xSize);
      }
    }
    parent -> permissible = NULL;
    parent -> permissibleIndx = NULL;
    parent -> permissibleIndxSize = 0;
  }
  unstackMPSign(parent);
  unstackFMPSign(parent);
  unstackNodeLMPIndex(parent);
  unstackNodeFLMPIndex(parent);
  if (parent -> splitInfo != NULL) {
    freeSplitInfo(parent -> splitInfo);
    parent -> splitInfo = NULL;
  }
  if (parent -> repMembrSizeAlloc > 0) {
    if (parent -> repMembrIndx != NULL) {
      free_uivector(parent -> repMembrIndx, 1, parent -> repMembrSizeAlloc);
      parent -> repMembrIndx = NULL;
    }
  }
  if (parent -> allMembrSizeAlloc > 0) {
    if (parent -> allMembrIndx != NULL) {
      free_uivector(parent -> allMembrIndx, 1, parent -> allMembrSizeAlloc);
      parent -> allMembrIndx = NULL;
    }
  }
  if (parent -> oobMembrSizeAlloc > 0) {
    if (parent -> oobMembrIndx != NULL) {
      free_uivector(parent -> oobMembrIndx, 1, parent -> oobMembrSizeAlloc);
      parent -> oobMembrIndx = NULL;
    }
  }
  free_gblock(parent, (size_t) sizeof(Node));
}
void setParent(Node *daughter, Node *parent) {
  daughter -> parent = parent;
}
void setLeftDaughter(Node *daughter, Node *parent) {
  parent -> left = daughter;
}
void setRightDaughter(Node *daughter, Node *parent) {
  parent -> right = daughter;
}
void stackMPSign(Node *tNode, unsigned int mpIndexSize) {
  if (tNode -> mpIndexSize > 0) {
    if (tNode -> mpIndexSize != mpIndexSize) {
      RF_nativePrint("\nRF-SRC:  *** ERROR *** ");
      RF_nativePrint("\nRF-SRC:  mpIndexSize has been previously defined:  %10d vs %10d", tNode -> mpIndexSize, mpIndexSize);
      RF_nativePrint("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tNode -> mpIndexSize = mpIndexSize;
  }
  tNode -> mpSign = ivector(1, tNode -> mpIndexSize);
}
void unstackMPSign(Node *tNode) {
  if(tNode -> mpIndexSize > 0) {
    if (tNode -> mpSign != NULL) {
      free_ivector(tNode -> mpSign, 1, tNode -> mpIndexSize);
      tNode -> mpSign = NULL;
    }
  }
}
void stackFMPSign(Node *tNode, unsigned int fmpIndexSize) {
  if (tNode -> fmpIndexSize > 0) {
    if (tNode -> fmpIndexSize != fmpIndexSize) {
      RF_nativePrint("\nRF-SRC:  *** ERROR *** ");
      RF_nativePrint("\nRF-SRC:  fmpIndexSize has been previously defined:  %10d vs %10d", tNode -> fmpIndexSize, fmpIndexSize);
      RF_nativePrint("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tNode -> fmpIndexSize = fmpIndexSize;
  }
  tNode -> fmpSign = ivector(1, tNode -> fmpIndexSize);
}
void unstackFMPSign(Node *tNode) {
  if(tNode -> fmpIndexSize > 0) {
    if (tNode -> fmpSign != NULL) {
      free_ivector(tNode -> fmpSign, 1, tNode -> fmpIndexSize);
      tNode -> fmpSign = NULL;
    }
  }
}
void stackNodeLMPIndex(Node *tNode, unsigned int size) {
  if (tNode -> lmpIndexAllocSize > 0) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  lmpIndex has been previously defined:  %10d vs %10d", tNode -> lmpIndexAllocSize, size);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  else {
    tNode -> lmpIndexAllocSize = size;
  }
  tNode -> lmpIndex = uivector(1, tNode -> lmpIndexAllocSize);
  tNode -> lmpValue = dvector(1, tNode -> lmpIndexAllocSize);
}
void unstackNodeLMPIndex(Node *tNode) {
  if(tNode -> lmpIndexAllocSize > 0) {
    if (tNode -> lmpIndex != NULL) {
      free_uivector(tNode -> lmpIndex, 1, tNode -> lmpIndexAllocSize);
      free_dvector(tNode -> lmpValue, 1, tNode -> lmpIndexAllocSize);
      tNode -> lmpIndex = NULL;
      tNode -> lmpValue = NULL;
      tNode -> lmpIndexAllocSize = 0;
    }
  }
}
void stackNodeFLMPIndex(Node *tNode, unsigned int size) {
  if (tNode -> flmpIndexAllocSize > 0) {
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  flmpIndex has been previously defined:  %10d vs %10d", tNode -> flmpIndexAllocSize, size);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  else {
    tNode -> flmpIndexAllocSize = size;
  }
  tNode -> flmpIndex = uivector(1, tNode -> flmpIndexAllocSize);
  tNode -> flmpValue = dvector(1, tNode -> flmpIndexAllocSize);
}
void unstackNodeFLMPIndex(Node *tNode) {
  if(tNode -> flmpIndexAllocSize > 0) {
    if (tNode -> flmpIndex != NULL) {
      free_uivector(tNode -> flmpIndex, 1, tNode -> flmpIndexAllocSize);
      free_dvector(tNode -> flmpValue, 1, tNode -> flmpIndexAllocSize);
      tNode -> flmpIndex = NULL;
      tNode -> flmpIndexAllocSize = 0;
    }
  }
}
void stackSplitDepth(Node *tNode, unsigned int depth) {
  if (tNode -> depth > 0) {
    if (tNode -> depth != depth) {
      RF_nativeError("\nRF-SRC:  *** ERROR *** ");
      RF_nativeError("\nRF-SRC:  depth has been previously defined:  %10d vs %10d", tNode -> depth, depth);
      RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
      RF_nativeExit();
    }
  }
  else {
    tNode -> depth = depth;
  }
  tNode -> splitDepth = uivector(1, tNode -> depth);
}
void unstackSplitDepth(Node *tNode) {
  if (tNode -> splitDepth != NULL) {
    free_uivector(tNode -> splitDepth, 1, tNode -> depth);
    tNode -> splitDepth = NULL;
  }
}
