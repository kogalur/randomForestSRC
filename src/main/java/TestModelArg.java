package com.kogalur.randomforest;

import org.apache.spark.sql.Dataset;
import java.util.ArrayList;
import java.util.Arrays;


import com.kogalur.randomforest.RFLogger;
import java.util.logging.Level;


import java.util.Random;

public class TestModelArg extends GenericModelArg {

    // GROW model arguments, that may be inherited.
    private ModelArg modelArg;

    private EnsembleArg ensembleArg;
    
    // Incoming Spark Dataset.
    private Dataset dataset;


    private String[] yTarget;
    private int[]    yTargetIndex;
    private int[]    yTargetFactorIndex;
    private int[]    yTargetNonFactorIndex;

    private int        fnSize;
    private int        fySize;
    private double[][] fyData;
    private double[][] fxData;

    private int[]      xMarginal;
    private int[]      xImportance;

    private Partial    partial;
    
    private int rfCores;    

    private int trace;
    private int seed;

    private Random generator;
    
    public TestModelArg(ModelArg modelArg) {
        mode = new String("rest");
        this.modelArg = modelArg;
        this.dataset = null;
        
        ensembleArg = new EnsembleArg(modelArg.getEnsembleArg(), mode, modelArg.get_family());

        setDefaultModelArg();
    }

    public TestModelArg(ModelArg modelArg, Dataset dataset) {
        mode = new String("pred");
        this.modelArg = modelArg;        
        this.dataset  = dataset;

        ensembleArg = new EnsembleArg(modelArg.getEnsembleArg(), mode, modelArg.get_family());

        setDefaultModelArg();
    }
        

    public ModelArg getModelArg() {
        return modelArg;
    }


    private void setDefaultModelArg() {

        // Common inherited variables from the grow default model arguments.
        this.ntree = modelArg.get_ntree();
        this.getTree = modelArg.get_getTree();
        
        // Set trace to the defalut value. Force it to non-zero values
        // here if debugging is required.
        set_trace(0);

        // Set seed.
        set_seed();

        // Default value of rfCores.
        set_rfCores(-1);

        set_yTarget();

        set_xMarginal();

        set_xImportance();

        set_partial();
        
    }


    /**
     * Sets the y-variable(s) to be targeted when multivariate
     * families are in force.  Predictions over this set of targets
     * will be included in the ensembles.
     * @param yTarget Array of y-variables to be targeted when
     * multivariate families are in force.  The default action is to
     * use all y-variables.
     */
    public void set_yTarget(String[] yTarget) {
        this.yTarget = yTarget;
        int size = yTarget.length;
        ArrayList <String> allVariableList = new ArrayList <String> (Arrays.asList(modelArg.getYvar()));

        yTargetIndex = new int[size];
        
        for (int i = 0; i < size; i++) {

            // Note the offset for the native code.
            yTargetIndex[i] = allVariableList.indexOf(yTarget[i]) + 1;
            if (yTargetIndex[i] <= 0) {
                RFLogger.log(Level.SEVERE, "Target y-variable not found in model:  " + yTarget[i]);
                throw new IllegalArgumentException();
            }
        }

        // After yTargetIndex has been initialized:
        initAuxiliaryTargetInfo();

    }

    /**
     * Sets the default action for the y-variable(s) to be targeted when multivariate
     * families are in force.  The default action is to
     * use all y-variables.
     * @see #set_yTarget(String[])
     *
     */
    public void set_yTarget() {
        
        // This is of type String[].  It might be nice to have these of zero (0) length
        // in the case of unsupervised data, instead of being null.  TBD TBD
        yTarget = modelArg.getYvar();
        if (yTarget == null) {
            yTargetIndex = new int[0];
        }
        else {
            yTargetIndex = new int[yTarget.length];

            // Note the offset for the native code.
            for (int i = 0; i < yTarget.length; i++) {
                yTargetIndex[i] = i + 1;
            }
        }

        // After yTargetIndex has been initialized:
        initAuxiliaryTargetInfo();
    }

    void initAuxiliaryTargetInfo() {
        int yTargetFactorCount;
        int yTargetNonFactorCount;

        // Count the number of target factors and non-factors.
        yTargetFactorCount = yTargetNonFactorCount = 0;

        // Zero (0) length arrays are handled nominally.
        for (int i = 0; i < yTargetIndex.length; i++) {
            if ((getModelArg().get_yLevelsMax())[i] > 0) {
                yTargetFactorCount ++;
            }
            else {
                yTargetNonFactorCount ++;
            }
        }

        yTargetFactorIndex = new int[yTargetFactorCount];
        yTargetNonFactorIndex = new int[yTargetNonFactorCount]; 

        // Initialize the indices of the target factors and non-factors.
        yTargetFactorCount = yTargetNonFactorCount = 0;

        for (int i = 0; i < yTargetIndex.length; i++) {
            if ((getModelArg().get_yLevelsMax())[i] > 0) {
                yTargetFactorIndex[yTargetFactorCount ++] = i;
            }
            else {
                yTargetNonFactorIndex[yTargetNonFactorCount ++] = i;
            }
        }
    }



    /**
     * Returns the value for the block size associated with the
     * reporting of the error rate.
     * @return The value for the block size associated with the
     * reporting of the error rate.
     * @see #set_blockSize(int)
     */
    public int get_blockSize() {
        return blockSize;
    }

    /**
     * Sets the default value for the block size associated with the
     * reporting of the error rate.  This value defaults to ntree.
     * @see #set_blockSize(int)
     */
    public void set_blockSize() {
        blockSize = get_ntree();
    }

    /**
     * Sets the specified value for the block size associated with the
     * reporting of the error rate.  If the value is out of range, the default value will be applied.
     * @see #set_blockSize()
     */
    public void set_blockSize(int blockSize) {
        if ((blockSize < 1) || (blockSize > get_ntree())) {
            RFLogger.log(Level.WARNING, "Invalid value for parameter blockSize:  " + blockSize);
            RFLogger.log(Level.WARNING, "Overriding blockSize with default value of ntree:  " + get_ntree());
            set_blockSize();
        }
        else {
            this.blockSize = blockSize;
        }
    }




    
    /**
     * Returns the array of y-variable to be targeted when multivariate
     * families are in force.  
     * @see #set_yTarget(String[])
     *
     */
    public String[] get_yTarget() {
        return yTarget;
    }

    public int[] get_yTargetIndex() {
        return yTargetIndex;
    }

    public int get_fnSize() {
        return 0;
    }

    public int get_fySize() {
        return 0;
    }

    public double[][] get_fyData() {
        return null;
    }

    public double[][] get_fxData() {
        return null;
    }
    
    public int[] set_xMarginal() {
        return null;
    }

    public int[] get_xMarginal() {
        return xMarginal;
    }
    
    public int[] set_xImportance() {
        return null;
    }

    public int[] get_xImportance() {
        return xImportance;
    }

    public Partial set_partial() {
        return null;
    }

    public Partial get_partial() {
        return partial;
    }
    
    public void set_rfCores(int rfCores) {
        int availableCores = Runtime.getRuntime().availableProcessors();

        this.rfCores = rfCores;
        if (rfCores > availableCores) {
            RFLogger.log(Level.WARNING, "Invalid value for parameter rfCores:  " + rfCores);
            RFLogger.log(Level.WARNING, "Overriding rfCores with (" + availableCores + "), the maximum available processors.");
            this.rfCores = availableCores;
        }
        
    }

    /**
    * Returns the number of cores to be used by the algorithm when OpenMP
    * parallel processing is in force.
    * @return the number of cores to be used by the algorithm when OpenMP parallel processing is in force.
    * @see #set_rfCores(int)
    */
    public int get_rfCores() {
        return rfCores;
    }


    int getEnsembleArgOptLow() {

        return (ensembleArg.getNative("varUsed") + 
                ensembleArg.getNative("splitDepth") +                 
                ensembleArg.getNative("importance") + 
                ensembleArg.getNative("proximity") +
                ensembleArg.getNative("forest") +
                ensembleArg.getNative("ensemble") +
                ensembleArg.getNative("errorType")); 
    }

    int getEnsembleArgOptHigh() {

        return (ensembleArg.getNative("membership") + 
                ensembleArg.getNative("weight") +                 
                ensembleArg.getNative("distance"));
    }

    
}
