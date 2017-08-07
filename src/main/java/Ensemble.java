package com.kogalur.randomforest;

import com.kogalur.randomforest.RFLogger;

import java.util.ArrayList;

import java.util.logging.Level;

class Ensemble {

    String name;
    
    
    Object ensembleVector;

    
    char type;

    
    long size;

    
    int identity;

    
    int auxDimSize;

    
    int[] auxDim;

    
    
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
