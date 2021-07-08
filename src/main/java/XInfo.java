package com.kogalur.randomforest;

class XInfo {

    int      xSize;
    char[]   xType;
    int[]    xLevelsMax;
    int[]    xLevelsCnt;
    int[]    xtType;
    int[]    stType;
    

    XInfo (int    xSize,
           char[] xType,
           int[]  xLevelsMax,
           int[]  xLevelsCnt,
           int[]  xtType,
           int[]  stType) {

        this.xSize = xSize;
        this.xType = xType;
        this.xLevelsMax = xLevelsMax;
        this.xLevelsCnt = xLevelsCnt;
        this.xtType = xtType;
        this.stType = stType;
    }
}
