package com.kogalur.randomforest;

class Partial {

    int partialType;
    int partialXvar;
    int partialLength;
    double[] partialValue;

    int partialLength2;
    int[] partialXvar2;
    double[] partialValue2;
    

    Partial (int partialType,
             int partialXvar,
             int partialLength,
             double[] partialValue,
             int partialLength2,
             int[] partialXvar2,
             double[] partialValue2) {

        this.partialType = partialType;
        this.partialXvar = partialXvar;
        this.partialLength = partialLength;
        this.partialValue = partialValue;
        this.partialLength2 = partialLength2;
        this.partialXvar2 = partialXvar2;
        this.partialValue2 = partialValue2;
    }
}
