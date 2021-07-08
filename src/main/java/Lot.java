package com.kogalur.randomforest;

class Lot {

    int hdim;
    int lotSize;
    int lotLag;
    int lotStrikeout;

    Lot (int hdim,
         int lotSize,
         int lotLag,
         int lotStrikeout) {

        this.hdim = hdim;
        this.lotSize = lotSize;
        this.lotLag = lotLag;
        this.lotStrikeout = lotStrikeout;
    }

    Lot () {

        this.hdim = hdim;
        this.lotSize = lotSize;
        this.lotLag = lotLag;
        this.lotStrikeout = lotStrikeout;
    }

}
