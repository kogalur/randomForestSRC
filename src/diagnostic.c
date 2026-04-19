
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "diagnostic.h"
#include "error.h"
#include "nrutil.h"
void getSplitObjectInfo(SplitInfo *info) {
  RF_nativePrint("\nSplitInfo:  %20x \n", info);
  RF_nativePrint("\n  info -> size        :    %20d", info -> size);
  RF_nativePrint("\n  info -> indicator   : 0x %20x", info -> indicator);
  RF_nativePrint("\n  info -> randomVar   : 0x %20x", info -> randomVar);
  RF_nativePrint("\n  info -> mwcpSizeAbs : 0x %20x", info -> mwcpSizeAbs);
  RF_nativePrint("\n  info -> randomPts   : 0x %20x", info -> randomPts);
  RF_nativePrint("\n   x-variable:   %10d", info -> randomVar[1]);
  RF_nativePrint("\n");
  int covariate = info -> randomVar[1];
  if (info -> mwcpSizeAbs[1] > 0) {
    RF_nativePrint(" (cov = %10d, mwcpPT =", covariate);
    for (uint m = 1; m <= info -> mwcpSizeAbs[1]; m++) {
      RF_nativePrint(" %10x", ((uint *) info -> randomPts[1])[m]);
    }
    RF_nativePrint(")");
  }
  else {
    RF_nativePrint(" (cov = %10d, spltPT = %12.4f) ", covariate, ((double *) info -> randomPts[1])[1]);
  }
  RF_nativePrint("\n");
}
void getNodeInfo(Node *nodePtr) {
  RF_nativePrint("\nNodeInfo:  (address, node) = (%20x, %10d)", nodePtr, nodePtr -> nodeID);
  if (nodePtr -> splitInfo != NULL) {
    getSplitObjectInfo(nodePtr -> splitInfo);
  }
  RF_nativePrint("\nSplit Statistic \n");
  RF_nativePrint(" %12.4f \n", nodePtr -> splitStatistic);
  RF_nativePrint("\nNode Variance \n");
  RF_nativePrint(" %12.4f \n", nodePtr -> variance);
  RF_nativePrint("\nPermissible Flag Size:          %10d", nodePtr -> xSize);
  RF_nativePrint("\n mpIndexSize   = %20d", nodePtr -> mpIndexSize);
  RF_nativePrint("\n fmpIndexSize  = %20d", nodePtr -> fmpIndexSize);
  RF_nativePrint("\n");
  RF_nativePrint("\n mpSign       = %20x", nodePtr -> mpSign);
  RF_nativePrint("\n fmpSign      = %20x", nodePtr -> fmpSign);
  RF_nativePrint("\n");
  RF_nativePrint("\n lmpIndexActualSize        = %20d", nodePtr -> lmpIndexActualSize);
  RF_nativePrint("\n flmpIndexActualSize       = %20d", nodePtr -> flmpIndexActualSize);
  RF_nativePrint("\n lmpIndexAllocSize         = %20d", nodePtr -> lmpIndexAllocSize);
  RF_nativePrint("\n flmpIndexAllocSize        = %20d", nodePtr -> flmpIndexAllocSize);
  RF_nativePrint("\n");
  RF_nativePrint("\n lmpIndex            = %20x", nodePtr -> lmpIndex);
  RF_nativePrint("\n flmpIndex           = %20x", nodePtr -> flmpIndex);
  RF_nativePrint("\n");
}
void getTerminalInfo(Terminal *termPtr) {
  RF_nativePrint("\nTerminalInfo:  %20x", termPtr);
  RF_nativePrint("\n  nodeID: %10d", termPtr -> nodeID);
  RF_nativePrint("\n");
  RF_nativePrint("\n lmiIndex            = %20x", termPtr -> lmiIndex);
  RF_nativePrint("\n lmiAllocSize        = %20d", termPtr -> lmiAllocSize);
  RF_nativePrint("\n lmiSize             = %20d", termPtr -> lmiSize);
  RF_nativePrint("\n lmiValue            = %20x", termPtr -> lmiValue);
  RF_nativePrint("\n rnfCount            = %20d", termPtr -> rnfCount);
  RF_nativePrint("\n meanResponse        = %20x", termPtr -> meanResponse);
  RF_nativePrint("\n membrCount          = %20d", termPtr -> membrCount);
  RF_nativePrint("\n membrStream         = %20x", termPtr -> membrStream);
}
Node *getTerminalNode(uint treeID, uint leaf) {
  uint i, j;
  Node *parent;
  parent = NULL;
  for (j = 1; j <= RF_observationSize; j++) {
    if ((RF_nodeMembership[treeID][j] -> nodeID) == leaf) {
      parent = RF_nodeMembership[treeID][j];
      j = RF_observationSize;
    }
  }
  if (parent == NULL) {
    RF_nativePrint("\nDiagnostic Trace of (individual, boot, node, leaf) vectors in data set:  ");
    RF_nativePrint("\n        index         boot         node         leaf \n");
    for (i = 1; i <= RF_observationSize; i++) {
      RF_nativePrint(" %12d %12d %12x %12d \n", i,
              RF_bootMembershipFlag[treeID][i], RF_nodeMembership[treeID][i],
              RF_nodeMembership[treeID][i] -> nodeID);
    }
    RF_nativePrint("\nDiagnostic State of TRAIN (SHADOW) data:  ");
    RF_nativePrint("\n       index       status         time   observations -> \n");
    RF_nativePrint("\n                                      ");
    for (i=1; i <= RF_xSize; i++) {
      RF_nativePrint(" %12d", i);
    }
    RF_nativePrint("\n");
    for (j = 1; j <= RF_observationSize; j++) {
      RF_nativePrint("%12d %12.4f %12.4f", j, RF_status[treeID][j], RF_time[treeID][j]);
      for (i=1; i <= RF_xSize; i++) {
        RF_nativePrint(" %12.4f", (RF_observation[treeID][i][j]));
      }
      RF_nativePrint("\n");
    }
    RF_nativePrint("\nDiagnostic State of TRAIN (INCOMING) data:  ");
    RF_nativePrint("\n       index       status         time   observations -> \n");
    RF_nativePrint("\n                                      ");
    for (i=1; i <= RF_xSize; i++) {
      RF_nativePrint(" %12d", i);
    }
    RF_nativePrint("\n");
    for (j = 1; j <= RF_observationSize; j++) {
      RF_nativePrint("%12d %12.4f %12.4f", j, RF_responseIn[RF_statusIndex][j], RF_responseIn[RF_timeIndex][j]);
      for (i=1; i <= RF_xSize; i++) {
        RF_nativePrint(" %12.4f", (RF_observationIn[i][j]));
      }
      RF_nativePrint("\n");
    }
    RF_nativeError("\nRF-SRC:  *** ERROR *** ");
    RF_nativeError("\nRF-SRC:  Proxy member for (tree, node) = (%12d, %12d) not found.", treeID, leaf);
    RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
    RF_nativeExit();
  }
  return parent;
}
void getRawNodeSize(uint  type,
                    uint  treeID,
                    Node *parent,
                    uint *repMembrIndx,
                    uint *repMembrSize,
                    uint *allMembrIndx,
                    uint *allMembrSize) {
  uint      obsSize;
  Node   ***nodeMembershipPtr;
  uint      bootMembrSize;
  uint i;
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetRawNodeSize() ENTRY ...\n");
  ${trace.token}  }
  obsSize           = 0;     
  nodeMembershipPtr = NULL;  
  switch (type) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    nodeMembershipPtr = RF_fnodeMembership;
    break;
  default:
    obsSize = RF_observationSize;
    nodeMembershipPtr = RF_nodeMembership;
    break;
  }
  bootMembrSize = RF_bootstrapSize;
  *repMembrSize = 0;
  for (i=1; i <= bootMembrSize; i++) {
    if (RF_nodeMembership[treeID][RF_bootMembershipIndex[treeID][i]] == parent) {
      repMembrIndx[++(*repMembrSize)] = RF_bootMembershipIndex[treeID][i];
    }
  }
  *allMembrSize = 0;
  for (i=1; i <= obsSize; i++) {
    if (nodeMembershipPtr[treeID][i] == parent) {
      allMembrIndx[++(*allMembrSize)] = i;
    }
  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\nParent Replicate Count:  %10d \n", *repMembrSize);
  ${trace.token}    RF_nativePrint("\nReplicate Absolute Index for Parent Node: \n");
  ${trace.token}    RF_nativePrint("\n       idx  rMembrIndx     \n");
  ${trace.token}    for (i=1; i <= (*repMembrSize); i++) {
  ${trace.token}      RF_nativePrint("%10d %10d \n", i, repMembrIndx[i]);
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\nParent Absolute Count:  %10d \n", *allMembrSize);
  ${trace.token}    RF_nativePrint("\nAbsolute Absolute Index for Parent Node: \n");
  ${trace.token}    RF_nativePrint("\n       idx  aMembrIndx     \n");
  ${trace.token}    for (i=1; i <= (*allMembrSize); i++) {
  ${trace.token}      RF_nativePrint("%10d %10d \n", i, allMembrIndx[i]);
  ${trace.token}    }
  ${trace.token}  }
  ${trace.token}  if (getTraceFlag(treeID) & SPLT_HGH_TRACE) {
  ${trace.token}    RF_nativePrint("\ngetRawNodeSize() EXIT ...\n");
  ${trace.token}  }
}
void printTreeInfo(uint treeID, Node *parent) {
  ${trace.token}      if (getTraceFlag(treeID)) {
  ${trace.token}        RF_nativePrint("%10d %10d %10d %10d \n", treeID, parent -> nodeID, parent -> depth, parent -> pseudoTerminal);
  ${trace.token}      }
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    printTreeInfo(treeID, parent ->  left);
    printTreeInfo(treeID, parent -> right);
  }
}
void printParameters(char mode) {
  ${trace.token}    switch (mode) {
  ${trace.token}    case RF_PRED:
  ${trace.token}      RF_nativePrint("\nMode is PRED.");
  ${trace.token}      RF_nativePrint("\nNumber of PRED individuals:  %10d", RF_fobservationSize);
  ${trace.token}      break;
  ${trace.token}    case RF_REST:
  ${trace.token}      RF_nativePrint("\nMode is REST.");
  ${trace.token}      RF_nativePrint("\nNumber of REST individuals:  %10d", RF_observationSize);
  ${trace.token}      break;
  ${trace.token}    case RF_PART:
  ${trace.token}      RF_nativePrint("\nMode is PART.");
  ${trace.token}      RF_nativePrint("\nNumber of PART individuals:  %10d", RF_observationSize);
  ${trace.token}      break;
  ${trace.token}    default:
  ${trace.token}      RF_nativePrint("\nMode is GROW.");
  ${trace.token}      RF_nativePrint("\nSplit rule is:               %10d", RF_splitRule);
  ${trace.token}      RF_nativePrint("\nValue of nsplit is:          %10d", RF_nsplit);
  ${trace.token}      RF_nativePrint("\nMinimum node size is:        %10d", RF_nodeSize);
  ${trace.token}      RF_nativePrint("\nMaximum node depth is:       %10d", RF_nodeDepth);
  ${trace.token}      RF_nativePrint("\nValue of mtry is             %10d", RF_mtry);
  ${trace.token}      RF_nativePrint("\nValue of bootstrap size is   %10d", RF_bootstrapSize);
  ${trace.token}      RF_nativePrint("\nNumber of GROW individuals:  %10d", RF_observationSize);
  ${trace.token}      RF_nativePrint("\nNumber of impute iterations: %10d", RF_nImpute);
  ${trace.token}      break;
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\nIncoming options lo:     0x  %10x", RF_opt);
  ${trace.token}    RF_nativePrint("\nIncoming options hi:     0x  %10x", RF_optHigh);
  ${trace.token}    RF_nativePrint("\nNumber of predictors:        %10d", RF_xSize);
  ${trace.token}    RF_nativePrint("\nNumber of trees in forest:   %10d \n", RF_ntree);
  ${trace.token}    testEndianness();
  ${trace.token}    int endian = 1;
  ${trace.token}    if (*(char *) &endian == 1) {
  ${trace.token}      RF_nativePrint("\nSystem is little-endian.  ");
  ${trace.token}    }
  ${trace.token}    else {
  ${trace.token}      RF_nativePrint("\nSystem is big-endian.  ");
  ${trace.token}    }
  ${trace.token}    RF_nativePrint("\nSize of (int)    is:  %10d   ", sizeof(int));
  ${trace.token}    RF_nativePrint("\nSize of (long)   is:  %10d   ", sizeof(long));
  ${trace.token}    RF_nativePrint("\nSize of (double) is:  %10d \n", sizeof(double));
  ${trace.token}    RF_nativePrint("\nNumber of threads:  %10d \n", RF_numThreads);
}
