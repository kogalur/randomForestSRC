package com.kogalur.randomforest;

class BaseLearn {

    int intrRule;
    int intrDepth;
    int syth;
    int dimReduce;

    BaseLearn (int intrRule,
               int intrDepth,
               int syth,
               int dimReduce) {

        this.intrRule = intrRule;
        this.intrDepth = intrDepth;
        this.syth = syth;
        this.dimReduce = dimReduce;
    }
}
