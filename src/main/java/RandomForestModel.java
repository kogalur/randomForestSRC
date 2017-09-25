package com.kogalur.randomforest;

import com.kogalur.randomforest.RFLogger;

import java.lang.Math;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import java.util.Iterator;
import java.util.Map;


import java.util.Random;

import java.util.logging.Level;

/**
 * Class that represents a random forest model as a result of growing, restoring, or predicting using new data.
 * @author Udaya Kogalur
 *
 */
public class RandomForestModel {

    private ModelType modelType;
    
    // Model argument, ensemble list contained in this model.
    private ModelArg modelArg;
    private HashMap ensembleList;
    // Type of model, grow, rest, or pred.
    private int type;
    
    // User arguments passed to the native code in !GROW mode.
    private int ptnCount;

    private String[] yTarget;

    // One-based vectors follow for native-code and native-code
    // emulation and parsing the ensemble output:
    private int[]    yTargetIndex;
    private int[]    yTargetFactorIndex;
    private int[]    yTargetNonFactorIndex;
    
    private String[] xMarginal;
    private int[]    xMarginalIndex;
    private String[] xImportance;
    private int[]    xImportanceIndex;

    private int      trace;
    private int      seed;

    private int      rfCores;
    
    // *****************************************************************
    // Auxiliary Variables.
    // *****************************************************************
    private Random generator;

    // *****************************************************************
    // REST / PRED Ensemble Requests.
    // *****************************************************************
    private EnsembleArg ensembleArg;
    
    // Constructor used by RandomForest.train(..) to create MAXI model.
    RandomForestModel(ModelArg modelArg, HashMap ensembleList) {
        this.modelArg = modelArg;
        this.ensembleList = ensembleList;
        modelType = ModelType.MAXI;

        // Create the ensemble object in advance.  
        ensembleArg = new EnsembleArg("rest", modelArg.get_family());

        // Set the default values for the paramaters the user passes to the native-code.
        set_pruningCount(0);
        set_xImportance(null);

         

        // The default action for the y-variable(s) to be targeted when multivariate
        // families are in force is to use all y-variables.
        set_yTarget();
        
        set_trace(0);
        set_seed();

        // Default value of rfCores, assumes grow setting.
        set_rfCores(modelArg.get_rfCores());

    }


    // Constructor used by RandomForest.predict(..) to create MAXI model.
    RandomForestModel(RandomForestModel model, HashMap ensembleList) {
        this.modelArg = model.getModelArg();
        this.ensembleList = ensembleList;
        this.modelType = ModelType.REST;

        // Create the !grow ensemble argument object in advance,
        // passing the grow-side ensemble argument object.
        ensembleArg = new EnsembleArg(modelArg.getEnsembleArg(), "rest", modelArg.get_family());

        // Set the default values for the paramaters the user passes to the native-code.
        set_pruningCount(model.get_pruningCount());
        set_xImportance(model.get_xImportance());

         

        set_yTarget(model.get_yTarget());
        
        set_trace(0);
        set_seed();

        // Default value of rfCores, assumes incoming model setting.
        set_rfCores(model.get_rfCores());

    }

    
    /**
     * Returns the class containing the model arguments that produced the 
     * {@link RandomForestModel} object.
     * @return The model arguments that produced the {@link RandomForestModel} object.
     */
    public ModelArg getModelArg() {
        return modelArg;
    }

    // Available to package only.  Not for end use.  Returns the
    // Ensemble object.  No error handling, null object
    // returned.
    Ensemble getEnsembleObj(String name) {
        Ensemble ensb = (Ensemble) ensembleList.get(name);
        return ensb;
    }

    // Available to package only.  Not for end use.  Returns the
    // massaged Ensemble:ensembleVector object.  No error handling,
    // null object returned. TBD TBD remove public TBD TBD
    public Object getEnsemble(String name) {
        Ensemble ensb = (Ensemble) ensembleList.get(name);
        Object value = null;
        if (ensb != null) {
            value = ensb.getValue(yTargetFactorIndex, getModelArg().get_yLevel());
        }
        return value;
    }



    /*
     * Returns the ensemble requested. The specification for each ensemble are listed below.  Ensure that
     * the appropriate cast is used when deploying the ensemble.  Ensembles are mode specific.
     *
     * <p><table class= "myColumnPadding">
     *  <tr>
     *    <th>Ensemble</th>
     *    <th>Type</th>
     *    <th>Dimension</th>
     *  </tr>
     *  <tr>
     *    <td></td>
     *    <td></td>
     *    <td></td>
     *  </tr>
     *  <tr>
     *    <td></td>
     *    <td></td>
     *    <td></td>
     *  </tr>
     *  <tr>
     *    <td></td>
     *    <td></td>
     *    <td></td>
     *  </tr>
     *  <tr>
     *    <td></td>
     *    <td></td>
     *    <td></td>
     *  </tr>
     *  <tr>
     *    <td></td>
     *    <td></td>
     *    <td></td>
     *  </tr>
     *  <tr>
     *    <td></td>
     *    <td></td>
     *    <td></td>
     *  </tr>
     * </table></p>
     */


    



    /** 
     * Sets the ensemble outputs desired when restoring the model or
     * using it to predict with new data.  These outputs are similar
     * to those available when initially creating a model using {@link
     * com.kogalur.randomforest.ModelArg#setEnsembleArg(String,
     * String)}, but are context sensitive to restoration versus
     * prediction.  These settings are in the form of &lt;key,
     * value&gt; pairs, where the key is the name of the ensemble, and
     * the value is the specific option for that ensemble.
     * <p> The default option for each key is in bold. </p>
     *
     * <pre> <code> 
     * <table class= "myColumnPadding">
     *  <tr>
     *    <th>Ensemble Key</th>
     *    <th>Possible Values</th>
     *  </tr>
     *  
      
     *  <tr>
     *    <td>weight</td>
     *    <td><b>no</b>, inbag, oob, all</td>
     *  </tr>

     *  <tr>
     *    <td>proximity</td>
     *    <td><b>no</b>, inbag, oob, all</td>
     *  </tr>
     *  
     *  <tr>
     *    <td>membership</td>
     *    <td><b>no</b>, yes</td>
     *  </tr>
     *
     *  <tr>
     *    <td>importance</td>
     *    <td><b>no</b>, permute, random, anti, <br> permute.joint, random.joint, anti.joint <br> permute.ensemble, random.ensemble, anti.ensemble, <br> permute.joint.ensemble, random.joint.ensemble, anti.joint.ensemble </td>
     *  </tr>
     *
     *  <tr>
     *    <td>error</td>
     *    <td><b>last.tree</b>, every.tree</td>
     *  </tr>
     *
     *  <tr>
     *    <td>varUsed</td>
     *    <td><b>no</b>, every.tree, sum.tree</td>
     *  </tr>
     *
     *  <tr>
     *    <td>splitDepth</td>
     *    <td><b>no</b>, every.tree, sum.tree</td>
     *  </tr>
     *
     *  <tr>
     *    <td>errorType</td>
     *    <td><i>For RF-C, RF-C+ Families Only:</i><br><b>default</b>, brier, g.mean, g.mean.drc</td>
     *  </tr>
     *
     * </table>
     * </code> </pre> 
     * Note that the options <code>inbag, oob</code> are only relevant when restoring a model, and <code>all</code> is the only
     * relevant option when predicting with new data.  Invalid values will be overridden.
     *
     * @param key The name of the ensemble output.
     * @param value The specific value for the ensemble output. 
     *
     */
    public void setEnsembleArg(String key, String value) {
        ensembleArg.set(key, value);

    }

    /**
     * Sets default values for the ensemble outputs resulting from restoring the model or using it to predict with new data.  
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


    int getEnsembleArgOptLow() {

        return (ensembleArg.getNative("varUsed") + 
                ensembleArg.getNative("splitDepth") +                 
                ensembleArg.getNative("importance") + 
                ensembleArg.getNative("proximity") +
                ensembleArg.getNative("errorType")); 
    }

    int getEnsembleArgOptHigh() {

        return (ensembleArg.getNative("membership") + 
                ensembleArg.getNative("weight") +                 
                 
                 
                
                ensembleArg.getNative("error"));
    }

    int getModelArgOptLow() {

        NativeOpt nativeOpt;
        int result = 0;
        
        nativeOpt = new NativeOpt("bootstrap",
                                  new String[] {"auto", "user"},
                                  new int[] {0, (1 << 19) + (1 << 20)});

        result += nativeOpt.get("bootstrap");

        return result;
    }

    int getModelArgOptHigh() {

        NativeOpt nativeOpt;
        int customSplitBit = 0;

        int result = 0;

        nativeOpt = new NativeOpt("sampleType",
                                  new String[] {"swr", "swor"},
                                  new int[] {0, (1 << 12)});

        result += nativeOpt.get("sampleType");
        
        return result;
    }


    /**
     * Sets the number of terminal nodes to which each tree in the
     * original model should be pruned back.  The terminal node membership
     * for the pruned forest is returned in the ensemble <code>pstnMembership</code>.
     * @param pruningCount The number of terminal nodes to which each tree in the 
     * original model should be pruned back.  The default action is to conduct no pruning, and is explicitly
     * achieved by setting this parameter to zero (0).
     */
    public void set_pruningCount(int pruningCount) {
        this.ptnCount = pruningCount;
    }
    /** Returns the value representing the number of terminal nodes to which each tree in the
     * original model should be pruned back.  
     * @see #set_pruningCount(int)
     */
    public int get_pruningCount() {
        return ptnCount;
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
            if ((getModelArg().get_yLevel())[i] > 0) {
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
            if ((getModelArg().get_yLevel())[i] > 0) {
                yTargetFactorIndex[yTargetFactorCount ++] = i;
            }
            else {
                yTargetNonFactorIndex[yTargetNonFactorCount ++] = i;
            }
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

    int[] getTargetIndex() {
        return yTargetIndex;
    }
    
    int getTargetSize() {
        int size = 0;
        if (yTargetIndex != null) {
            size = yTargetIndex.length;
        }
        return size;
    }


     
        

    /**
     * Sets the x-variable(s) over which importance calculations should
     * be conducted.  This is useful when <code>xSize</code> is large and importance 
     * over a small subset of x-variables is desired.  In addition, if a joint importance calculation is 
     * desired, the specified x-variables are jointly targeted in the calculations.
     * @param xImportance Array of x-variables over which important calculations should be conducted.
     * The default action is to not provide importance ensembles. This can be
     * alse be accomplished by sending a null value into the method.
     */
    public void set_xImportance(String[] xImportance) {

        this.xImportance = xImportance;

        if (xImportance == null) {
            xImportanceIndex = null;
        }
        else {
            ArrayList <String> allVariableList = new ArrayList <String> (Arrays.asList(modelArg.getXvar()));
            
            xImportanceIndex = new int[xImportance.length];
            
            for (int i = 0; i < xImportance.length; i++) {
                
                // Note the offset for the native code.
                xImportanceIndex[i] = allVariableList.indexOf(xImportance[i]) + 1;
                if (xImportanceIndex[i] <= 0) {
                    RFLogger.log(Level.SEVERE, "Importance x-variable not found in model:  " + xImportance[i]);
                    throw new IllegalArgumentException();
                }
            
            }
        }
    }

    /**
     * Sets the x-variable(s), over which importance calculations should
     * be conducted, to all x-variables.  Be cautious when <code>xSize</code> is large, as computation times can be lengthy.
     * @see #set_xImportance(String[])
     */
    public void set_xImportance() {
        xImportance = modelArg.getXvar();
        xImportanceIndex = new int[xImportance.length];
        
        // Note the offset for the native code.
        for (int i = 0; i < xImportance.length; i++) {
            xImportanceIndex[i] = i + 1;
        }
    }
    
    /**
     * Returns the x-variable(s) over which importance calculations should
     * be conducted.
     * @see #set_xImportance(String[])
     */
    public String[] get_xImportance() {
        return xImportance;
    }

    int[] getImportanceIndex() {
        return xImportanceIndex;
    }
    
    int getImportanceSize() {
        int size = 0;
        if (xImportanceIndex != null) {
            size = xImportanceIndex.length;
        }
        return size;
    }


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

    private void set_seed() {
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

    public ModelType getModelType() {
        return modelType;
    }

    void printEnsembleList() {
        Set entrySet = ensembleList.entrySet();
        Iterator itr = entrySet.iterator();
        int actLength = 0;
        long thrLength = 0;
        int i = 0;
        while(itr.hasNext()){
            Map.Entry me = (Map.Entry) itr.next();
            Ensemble ensb = (Ensemble) me.getValue();
            ensb.printMetaInfo();
            i++;
        }
    }    

    // TBD TBD remove public option TBD TBD
    public void printEnsemble(Object ensb) {
        if (ensb == null) {
            RFLogger.log(Level.INFO, "Ensemble:  null");
        }
        else {
            String s;
            RFLogger.log(Level.INFO, "Ensemble:  ");

            if (ensb instanceof int[]) {
                int[] value = (int[]) ensb;
                s = new String("      index      value");
                RFLogger.log(Level.INFO, s);
                for (int i = 0; i < value.length; i++) {
                    s = String.format(" %1$10d %2$10d", i, value[i]);
                    RFLogger.log(Level.INFO, s);
                }
            }
            else if (ensb instanceof int[][]) {
                int[][] value = (int[][]) ensb;
                s = new String("      index");
                // Tag columns, noting that header row is hard-coded
                // to the length of the first row.
                for (int j = 0; j < value[0].length; j++) {
                    s = s.concat(String.format(" %1$10d", j)); 
                }
                RFLogger.log(Level.INFO, s);
                for (int i = 0; i < value.length; i++) {
                    s = String.format(" %1$10d", i);
                    for (int j = 0; j < value[i].length; j++) {
                        s = s.concat(String.format(" %1$10d", value[i][j]));
                    }
                    RFLogger.log(Level.INFO, s);                    
                }
            }
            else if (ensb instanceof int[][][]) {
                int[][][] value = (int[][][]) ensb;
                for (int i = 0; i < value.length; i++) {
                    RFLogger.log(Level.INFO, "[" + i + "] ==>>");
                    s = new String("      index");
                    // Tag columns, noting that header row is hard-coded
                    // to the length of the first row.
                    for (int k = 0; k < value[i][0].length; k++) {
                        s = s.concat(String.format(" %1$10d", k)); 
                    }
                    RFLogger.log(Level.INFO, s);
                    for (int j = 0; j < value[i].length; j++) {
                        s = String.format(" %1$10d", j);
                        for (int k = 0; k < value[i][j].length; k++) {
                            s = s.concat(String.format(" %1$10d", value[i][j][k]));
                        }
                        RFLogger.log(Level.INFO, s);                    
                    }
                }
                
            }
            else if (ensb instanceof int[][][][]) {
                int[][][][] value = (int[][][][]) ensb;
                for (int i = 0; i < value.length; i++) {
                    for (int j = 0; i < value[i].length; j++) {
                        RFLogger.log(Level.INFO, "[" + i + "][" + j + "] ==>>");
                        s = new String("      index");
                        // Header row is hard-coded to the length of the first row.
                        for (int m = 0; m < value[i][j][0].length; m++) {
                            s = s.concat(String.format(" %1$10d", m)); 
                        }
                        RFLogger.log(Level.INFO, s);
                        for (int k = 0; k < value[i][j].length; k++) {
                        s = String.format(" %1$10d", k);
                        for (int m = 0; m < value[i][j][k].length; m++) {
                            s = s.concat(String.format(" %1$10d", value[i][j][k][m]));
                        }
                        RFLogger.log(Level.INFO, s);                    
                        }
                    }
                }
                
            }
            else if (ensb instanceof double[]) {
                double[] value = (double[]) ensb;
                
            }
            else if (ensb instanceof double[][]) {
                double[][] value = (double[][]) ensb;
                s = new String("      index");
                // Tag columns, noting that header row is hard-coded
                // to the length of the first row.
                for (int j = 0; j < value[0].length; j++) {
                    s = s.concat(String.format(" %1$20d", j)); 
                }
                RFLogger.log(Level.INFO, s);
                for (int i = 0; i < value.length; i++) {
                    s = String.format(" %1$10d", i);
                    for (int j = 0; j < value[i].length; j++) {
                        s = s.concat(String.format(" %1$20.4f", value[i][j]));
                    }
                    RFLogger.log(Level.INFO, s);                    
                }

            }
            else if (ensb instanceof double[][][]) {
                double[][][] value = (double[][][]) ensb;
                for (int i = 0; i < value.length; i++) {
                    RFLogger.log(Level.INFO, "[" + i + "] ==>>");
                    s = new String("      index");
                    // Tag columns, noting that header row is hard-coded
                    // to the length of the first row.
                    for (int k = 0; k < value[i][0].length; k++) {
                        s = s.concat(String.format(" %1$20d", k)); 
                    }
                    RFLogger.log(Level.INFO, s);
                    for (int j = 0; j < value[i].length; j++) {
                        s = String.format(" %1$10d", j);
                        for (int k = 0; k < value[i][j].length; k++) {
                            s = s.concat(String.format(" %1$20.4f", value[i][j][k]));
                        }
                        RFLogger.log(Level.INFO, s);                    
                    }
                }
                
            }
            else if (ensb instanceof double[][][][]) {
                double[][][][] value = (double[][][][]) ensb;
            }
            
        }
    }
    
}
