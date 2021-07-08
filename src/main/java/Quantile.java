package com.kogalur.randomforest;

class Quantile {

    double[] prob;
    double   probEpsilon;

    Quantile (double[] prob, double probEpsilon) {

        this.prob = prob;
        this.probEpsilon = probEpsilon;
    }
}
