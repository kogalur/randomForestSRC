package com.kogalur.randomforest;

enum EnsembleType {
    
    CHAR   ((byte) 0),
    INT    ((byte) 1),
    DOUBLE ((byte) 2);
    
    private byte value;

    private EnsembleType(byte value) {
        this.value = value;
    }

    byte getValue() {
        return this.value;
    }
    
}
