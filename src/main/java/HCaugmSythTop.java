package com.kogalur.randomforest;

class HCaugmSythTop {

    int[][]     nodeCount;
    int[][]     treeID;
    int[][]     nodeID;
    int[][]     hcDim;
    int[][]     parmID;
    double[][]  contPT;
    double[][]  contPTR;
    int[][]     mwcpSZ;
    int[][]     mwcpPT;
    
    
    HCaugmSythTop (int[][]     nodeCount,
                   int[][]     treeID,
                   int[][]     nodeID,
                   int[][]     hcDim,
                   int[][]     parmID,
                   double [][] contPT,
                   double[][]  contPTR,
                   int[][]     mwcpSZ,
                   int[][]     mwcpPT) {

        this.nodeCount = nodeCount;
        this.treeID    = treeID;
        this.nodeID    = nodeID;
        this.hcDim     = hcDim;
        this.parmID    = parmID;
        this.contPT    = contPT;
        this.contPTR   = contPTR;
        this.mwcpSZ    = mwcpSZ;
        this.mwcpPT    = mwcpPT;
        
    }
}
