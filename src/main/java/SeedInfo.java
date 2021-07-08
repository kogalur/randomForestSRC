package com.kogalur.randomforest;

class SeedInfo {

    int[]    seed;
    int[]    seedVimp;
    int      optLowGrow;

    SeedInfo (int[] seed, int[] seedVimp, int optLowGrow) {
        this.seed = seed;
        this.seedVimp = seedVimp;
        this.optLowGrow = optLowGrow;
    }
}
