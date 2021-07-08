package com.kogalur.randomforest;

import java.util.TreeSet;
import java.util.Iterator;

class YInfo {

    int      ySize;
    char[]   yType;
    int[]    yLevelsMax;
    int[]    yLevelsCnt;
    int[]    eventType;
    int      eventTypeSize;
    int[]    subjID;

    YInfo (int      ySize,
           char[]   yType,
           int[]    yLevelsMax,
           int[]    yLevelsCnt,
           TreeSet <Integer> eventTypeSet,
           int[]    subjIn) {

        this.ySize = ySize;
        this.yType = yType;
        this.yLevelsMax = yLevelsMax;
        this.yLevelsCnt = yLevelsCnt;
        this.eventType = eventType;
        this.eventTypeSize = eventType.length;
        this.subjID = subjID;

        if (eventTypeSet != null) {
            eventType = new int[eventTypeSet.size()];

            Iterator<Integer> iiter = eventTypeSet.iterator();

            int iter = 0;
            while (iiter.hasNext()) {
                eventType[iter] = iiter.next();
            }
        }
        else {
            eventType = null;
        }
        
    }
}
