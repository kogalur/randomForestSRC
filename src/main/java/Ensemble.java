package com.kogalur.randomforest;

import com.kogalur.randomforest.RFLogger;

import java.util.ArrayList;

import java.util.logging.Level;

class Ensemble {

    String name;
    
    // Generic object containing the ensemble vector.  This will be of type (double) or (int).
    Object ensembleVector;

    // Type of ensemble (double) or (int).
    char type;

    // Size of ensemble vector;
    long size;

    // Identity of ensemble vector [2, ..., RF_SEXP_CNT) but in practice, we use [0, 1 << 6).
    int identity;

    // Multi-dimensional view of the ensemble.  This is the number of dimensions.
    int auxDimSize;

    // Multi-dimensional view of the ensemble.  This is the size of the internal dimensions.
    int[] auxDim;

    // Irregularity flag indicating ragged array.  In such cases, the m-dim view must be manually
    // constructed from what we know about the ensemble and data set.
    boolean ragged;
    
    void Ensemble(String name, char type, int identity, long size, boolean ragged, int auxDimSize, int[] auxDim, Object ensembleVector) {
        this.name = name;
        this.type = type;
        this.identity = identity;
        this.ragged = ragged;
        this.size = size;
        this.auxDimSize = auxDimSize;
        this.auxDim = auxDim;

        this.ensembleVector = ensembleVector;
            
    }

   
}
