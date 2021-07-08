package com.kogalur.randomforest;

class SampleInfo {

    String   bootstrap;
    String   sampleType;
    int      sampleSize;
    int[][]  sample;
    double[] caseWeight;

    SampleInfo (String   bootstrap,
                String   sampleType,
                int      sampleSize,
                int[][]  sample,
                double[] caseWeight) {


        this.bootstrap  = bootstrap;
        this.sampleType = sampleType;
        this.sampleSize = sampleSize;
        this.sample     = sample;
        this.caseWeight = caseWeight;
    }
}
