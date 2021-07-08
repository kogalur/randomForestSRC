package com.kogalur.randomforest;

import com.kogalur.randomforest.RFLogger;

import java.util.logging.Level;

import java.util.Random;

class GenericModelArg {

    protected Random generator;
    protected int trace;
    
    protected int seed;

    protected int        rfCores;

    protected int        blockSize;
    protected int        ntree;
    protected int[]      getTree;
    protected Quantile   quantile;

    protected String mode;

    private EnsembleArg ensembleArg;
    
    
    /** 
    * Sets the seed for the random number generator used by the
    * algorithm.  This must be a negative number.  The seed is a very
    * important parameter if repeatability of the model generated is
    * required.  Generally speaking, growing a model using the same
    * data set, the same model paramaters, and the same seed will
    * result in identical models. When large amounts of missing data
    * are involved, there can be slight variations due to Monte Carlo
    * effects. If the parameter is not set by the user, it can always
    * be recovered with {@link #get_seed}.
    */
    public void set_seed(int seed) {
        this.seed = seed;
        
        if (seed >= 0) {
            if (seed > 0) {
                RFLogger.log(Level.WARNING, "Invalid value for parameter seed:  " + seed);                
                RFLogger.log(Level.WARNING, "Overriding seed with negative of the specified value.");
                this.seed = - seed;
            }
            else {
                RFLogger.log(Level.WARNING, "Invalid value for parameter seed:  " + seed);
                RFLogger.log(Level.WARNING, "Overriding seed with a negative random integer value.");
                set_seed();
            }
        }
    }

    public void set_seed() {
        generator = new Random();
        seed = - Math.abs(generator.nextInt());
    }

    /** 
    * Returns the seed for the random number generator used by the
    * algorithm.
    * @return The seed for the random number generator used by the
    * @see #set_seed(int)
    */
    public int get_seed() {
        return seed;
    }

        /**
     * Sets the trace parameter indicating the specified update
     * interval in seconds.  
     * @param trace The trace parameter indicating the specified update
     * interval in seconds. During extended execution times,
     * the approximate time to complete the execution is output to a
     * trace file in the users HOME directory.  The format and
     * location of the trace file can be controlled by modifying
     * <code>src/main/resources/spark/log.properties</code>.  A value
     * of zero (0) turns off the trace.
     */
    public void set_trace(int trace) {
        this.trace = trace;
        
        if (trace < 0) {
            RFLogger.log(Level.WARNING, "Invalid value for parameter trace:  " + trace);
            RFLogger.log(Level.WARNING, "Overriding trace with default value:  0");
            set_trace(0);
        }
    }
        
    /**
     * Returns the trace parameter indicating the specified update interval in seconds.
     * @return The trace parameter indicating the specified update interval in seconds.
     */
    public int get_trace() {
        return trace;
    }


    /**
    * Sets the number of cores to be used by the algorithm when OpenMP
    * parallel processing is in force.  
    * @param rfCores The number of cores to be used by the algorithm when OpenMP
    * parallel processing is in force. The default behaviour is to
    * use all cores available.  This is achieved by setting the
    * parameter to a negative value.  The result is that each core
    * will be independently tasked with growing a tree.  Significant
    * savings in elapsed computation times can be achieved.
    */
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

    /**
     * Returns the number of trees in the forest.
     * @return The number of trees in the forest.
     * @see #set_bootstrap(int, String, String, int, int[][], double[])
     */
    public int get_ntree() {
        return ntree;
    }


    public void set_getTree() {
        int[] getTree = new int[ntree];
        
        for (int i = 1; i <= ntree; i++) {
            getTree[i] = 1;
        }
    }
    
    public void set_getTree(int[] getTree) {
        this.getTree = getTree;
    }


    public int[] get_getTree() {
        return getTree;
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
     * Sets the target quantile probabilities in scenarios that require them.
     * @param prob The target quantile probability vector.  This must be a vector 
     * of probabilities between zero and one. When null is sent in, the vector is ingored.
     * Also sets the Greenwald-Khanna allowable error for the target quantiles.
     * @param probEpsilon The Greenwald-Khanna allowable error for the target quantiles.
     */
    public void set_quantile(double[] prob, double probEpsilon) {

        if (prob == null) {
            quantile = null;
        }
        else {
            if (probEpsilon == 0) {
                quantile = new Quantile(prob, 0.005);
            }
            else {
                quantile = new Quantile(prob, probEpsilon);
            }
        }
    }

    public void set_quantile() {
        quantile = null;
    }

    /**
     * Returns the vector of target quantile probabilities.
     * @return The vector of target quantile probabilities.
     * @see #set_quantile(double[])
     */
    public Quantile get_quantile() {
        return quantile;
    }

    public String getMode() {
        return mode;
    }


    public void setEnsembleArg(String key, String value) {
        ensembleArg.set(key, value);
    }

    /**
     * Sets default values for the ensemble outputs resulting from the model.
     * @see #setEnsembleArg(String, String).
     */
    public void setEnsembleArg() {
        ensembleArg.set();
    }

    /**
     * Returns the current value for the specified ensemble argument.
     * @param key The name of the ensemble output.
     * @see #setEnsembleArg(String, String).
     */
    public String getEnsembleArg(String key) {
        return ensembleArg.get(key);
    }
    

    EnsembleArg getEnsembleArg() {
        return ensembleArg;
    }

}    
