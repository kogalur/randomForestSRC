package com.kogalur.randomforest;

import com.kogalur.randomforest.RFLogger;


import java.util.Arrays;

import java.util.logging.Level;

class Ensemble {

    /* 
       >>> CAUTION 
       These are mapped to the native types in src/main/global.h 
       <<< CAUTION
    */

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
    
    // Name of the ensemble.
    String name;
    
    // Type of ensemble (char), (int), or (double).
    EnsembleType type;

    // Identity of ensemble vector [2, ..., RF_SEXP_CNT) but in practice, we use [0, 1 << 6).
    int identity;

    // Size of ensemble vector;
    long size;

    // Irregularity flag indicating ragged array.  In such cases, the
    // multi-dimensional view must be explicitly constructed from what
    // we know about the ensemble and data set.
    boolean ragged;

    // Multi-dimensional view of the ensemble.  This is the number of dimensions.
    int auxDimSize;

    // Multi-dimensional view of the ensemble.  This is the size of the internal dimensions.
    int[] auxDim;
    
    // Generic object containing the ensemble vector.  This will be of type (double) or (int).
    Object ensembleVector;

    // Constructor called from native-code from entry.c.  There, we
    // loop over a list and populate each ensemble and it's associated
    // fields.
    void Ensemble(String name,
                  byte nativeType,
                  int identity,
                  long size,
                  boolean ragged,
                  int auxDimSize,
                  int[] auxDim,
                  Object ensembleVector) {

        int actualLength = 0;
        
        if (nativeType == EnsembleType.CHAR.getValue()) {
            type = EnsembleType.CHAR;
            actualLength = ((char[]) (ensembleVector)).length;
        }
        else if (nativeType == EnsembleType.INT.getValue()) {
            type = EnsembleType.INT;
            actualLength = ((int[]) (ensembleVector)).length;
        }
        else if (nativeType == EnsembleType.DOUBLE.getValue()) {
            type = EnsembleType.DOUBLE;
        }                                    
        else {
            RFLogger.log(Level.SEVERE, "Ensemble type (native) not found:  " + nativeType);
            throw new IllegalArgumentException();
        }                            


        this.name = name;
        this.identity = identity;
        this.ragged = ragged;
        this.size = size;
        this.auxDimSize = auxDimSize;
        this.auxDim = auxDim;

        this.ensembleVector = ensembleVector;
       
    }

    Object getValue(int[] yTargetFactorIndex,
                    int[] yLevel) {

        Object result = null;
        
        char[]         charVector    = null;
        int[]          intVector     = null;
        int[]          int1Vector    = null;
        int[][]        int2Vector    = null;
        int[][][]      int3Vector    = null;
        int[][][][]    int4Vector    = null;
        double[]       doubleVector  = null;
        double[]       double1Vector = null;        
        double[][]     double2Vector = null;
        double[][][]   double3Vector = null;
        double[][][][] double4Vector = null;

        int offsetLow, offsetHigh;
        int dim1, dim2, dim3, dim4;
        
        if (type.equals(Ensemble.EnsembleType.INT)) {
            intVector = (int[]) ensembleVector;

            if (auxDimSize == 4) {
                
                offsetLow = offsetHigh = 0;
                dim1 = getAuxDim(0, 0, yTargetFactorIndex, yLevel);
                int4Vector = new int[dim1][][][];
                for (int i = 0; i < dim1; i++) {
                    dim2 = getAuxDim(i, 1, yTargetFactorIndex, yLevel);
                    int4Vector[i] = new int[dim2][][];
                    for (int j = 0; j < dim2; j++) {
                        dim3 = getAuxDim(j, 2, yTargetFactorIndex, yLevel);
                        int4Vector[i][j] = new int[dim3][];
                        for (int k = 0; k < dim3; k++) {
                            dim4 = getAuxDim(j, 3, yTargetFactorIndex, yLevel);
                            int4Vector[i][j][k] = new int[dim4];
                            offsetHigh = offsetLow + dim4;
                            int4Vector[i][j][k] = Arrays.copyOfRange(intVector, offsetLow, offsetHigh);
                            offsetLow = offsetHigh;
                        }
                    }
                }
                result = int4Vector;
            }
            else if (auxDimSize == 3) {

                offsetLow = offsetHigh = 0;
                dim1 = getAuxDim(0, 0, yTargetFactorIndex, yLevel);
                int3Vector = new int[dim1][][];
                for (int i = 0; i < dim1; i++) {
                    dim2 = getAuxDim(i, 1, yTargetFactorIndex, yLevel);
                    int3Vector[i] = new int[dim2][];
                    for (int j = 0; j < dim2; j++) {
                        dim3 = getAuxDim(j, 2, yTargetFactorIndex, yLevel);
                        int3Vector[i][j] = new int[dim3];
                        offsetHigh = offsetLow + dim3;
                        int3Vector[i][j] = Arrays.copyOfRange(intVector, offsetLow, offsetHigh);
                        offsetLow = offsetHigh;
                    }
                }
                result = int3Vector;
            }
            else if (auxDimSize == 2) {

                offsetLow = offsetHigh = 0;
                dim1 = getAuxDim(0, 0, yTargetFactorIndex, yLevel);
                int2Vector = new int[dim1][];
                for (int i = 0; i < dim1; i++) {
                    dim2 = getAuxDim(i, 1, yTargetFactorIndex, yLevel);
                    int2Vector[i] = new int[dim2];
                    offsetHigh = offsetLow + dim2;
                    int2Vector[i] = Arrays.copyOfRange(intVector, offsetLow, offsetHigh);
                    offsetLow = offsetHigh;
                }
                result = int2Vector;
            }
            else if (auxDimSize == 1) {

                // Allows for trimming 1-D arrays.
                offsetLow = offsetHigh = 0;                
                dim1 = getAuxDim(0, 0, yTargetFactorIndex, yLevel);
                int1Vector = new int[dim1];
                offsetHigh = offsetLow + dim1;
                int1Vector = Arrays.copyOfRange(intVector, offsetLow, offsetHigh);
                result = int1Vector;
            }
        }
        else if (type.equals(Ensemble.EnsembleType.DOUBLE)) {
            doubleVector = (double[]) ensembleVector;

            if (auxDimSize == 4) {

            }
            
            else if (auxDimSize == 3) {



                
                offsetLow = offsetHigh = 0;
                dim1 = getAuxDim(0, 0, yTargetFactorIndex, yLevel);

                double3Vector = new double[dim1][][];

                for (int i = 0; i < dim1; i++) {
                    dim2 = getAuxDim(i, 1, yTargetFactorIndex, yLevel);
                    double3Vector[i] = new double[dim2][];
                    for (int j = 0; j < dim2; j++) {
                        dim3 = getAuxDim(j, 2, yTargetFactorIndex, yLevel);
                        double3Vector[i][j] = new double[dim3];
                        offsetHigh = offsetLow + dim3;
                        double3Vector[i][j] = Arrays.copyOfRange(doubleVector, offsetLow, offsetHigh);
                        offsetLow = offsetHigh;
                    }
                }
                result = double3Vector;
            }
            else if (auxDimSize == 2) {

                offsetLow = offsetHigh = 0;
                dim1 = getAuxDim(0, 0, yTargetFactorIndex, yLevel);
                double2Vector = new double[dim1][];
                for (int i = 0; i < dim1; i++) {
                    dim2 = getAuxDim(i, 1, yTargetFactorIndex, yLevel);
                    double2Vector[i] = new double[dim2];
                    offsetHigh = offsetLow + dim2;
                    double2Vector[i] = Arrays.copyOfRange(doubleVector, offsetLow, offsetHigh);
                    offsetLow = offsetHigh;
                }
                result = double2Vector;
            }
            else if (auxDimSize == 1) {

            }
            


        }
        else {
            RFLogger.log(Level.SEVERE, "Ensemble (name, type) not supported for (" + name + ", " + type + ")");
            RFLogger.log(Level.SEVERE, "Please Contact Technical Support.");
            throw new IllegalArgumentException();
        }
        
        return result;

    }

    // yTargetIndex[] contains the y-var(s) in the target ensembles.
    // This is [1,...,ySize] in grow mode, and a subset thereof in
    // !grow mode.  yTargetFactorIndex[] contains the indices of the
    // y-vars(s) that are factors.  yLevel[] contains the number of
    // levels in each categorical x-var.  If an x-var i is not
    // categorical, yLevel[i] == 0.
    int getAuxDim(int preIndex, int postIndex,
                  int[] yTargetFactorIndex,
                  int[] yLevel) {

        int result = 0;

        
        // Note that when postIndex == 1, preIndex is ignored, as the first
        // element is non-ragged, and we force this function to return the
        // constant value for that dimension.

        if (postIndex == 0) {
            result = auxDim[postIndex];
        }
        else if (auxDim[postIndex] >= 1) {
            result = auxDim[postIndex];
        }
        else if (auxDim[postIndex] == 0) {
            result = yLevel[yTargetFactorIndex[preIndex]];
        }
        else if (auxDim[postIndex] == -1) {
            result = 1 + yLevel[yTargetFactorIndex[preIndex]];
        }
        else if (auxDim[postIndex] == -2) {
            // result = RF_tLeafCount[preIndex];
        }
        else {
            // RF_nativeError("\nRF-SRC:  *** ERROR *** ");
            // RF_nativeError("\nRF-SRC:  Inconsistent internal dimension of auxiliary array in getAuxDim():  %10d", dim[postIndex]);
            // RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
        }
        return result;
    }
    
    void printMetaInfo() {

        int actLength = 0;
        long metaLength = 0;
        
        RFLogger.log(Level.INFO, "Ensemble Meta Info:  ");
        RFLogger.log(Level.INFO, "        name : " + name);
        RFLogger.log(Level.INFO, "        type : " + type);
        RFLogger.log(Level.INFO, "    identity : " + identity);
        RFLogger.log(Level.INFO, "        size : " + size);
        RFLogger.log(Level.INFO, "      ragged : " + ragged);
        RFLogger.log(Level.INFO, "  auxDimSize : " + auxDimSize);
        for (int i = 0; i < auxDimSize; i++) {
            RFLogger.log(Level.INFO, "   auxDim[" + i + "] : " + auxDim[i]);
        }
        
        if (ensembleVector == null) {
            metaLength = size;
            actLength = 0;
            RFLogger.log(Level.WARNING, "Ensemble vector is null.  ");
        }
        else {
            metaLength = size;
            if (type.equals(Ensemble.EnsembleType.CHAR)) {
                actLength = ((char[]) (ensembleVector)).length;
            }
            else if (type.equals(Ensemble.EnsembleType.INT)) {
                actLength = ((int[]) (ensembleVector)).length;
            }
            else if (type.equals(Ensemble.EnsembleType.DOUBLE)) {
                actLength = ((double[]) (ensembleVector)).length;
            }                                    
        }
        if ( (long) actLength != metaLength) {
            RFLogger.log(Level.WARNING, "Ensemble vector length and object meta info inconsistent:  " + actLength + " versus " + metaLength);
        }
    }

}
