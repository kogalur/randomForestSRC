package com.kogalur.randomforest;

import com.kogalur.randomforest.ModelArg;
import com.kogalur.randomforest.RandomForestModel;

import com.kogalur.randomforest.Native;
import com.kogalur.randomforest.RFLogger;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Set;
import java.util.Iterator;
import org.apache.spark.sql.Dataset;

import java.util.logging.Level;

/**
 * Class that provides the static methods to train a random forest, to
 * restore a previouly created random forest, and to predict using new
 * test data with a previously created random forest.
 * @author Udaya Kogalur
 * 
 */
public class RandomForest {

    private RandomForest() {

    }

    /**
     * Trains a random forest given the user-defined model arugments.
     * @param modelArg User-defined model arguments.
     */
    public static RandomForestModel train(ModelArg modelArg)  {

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Family train:  " + modelArg.get_family());
        @RF_TRACE_OFF@  }

        LinkedHashMap   ensembleList = Native.getInstance().grow(modelArg.get_trace(),
                                                                 modelArg.get_seed(),

                                                                 modelArg.getEnsembleArgOptLow() + modelArg.getModelArgOptLow(),
                                                                 modelArg.getEnsembleArgOptHigh() + modelArg.getModelArgOptHigh(),
                                                                 
                                                                 modelArg.getSplitRuleID(modelArg.get_splitRule()),
                                                                 modelArg.get_nSplit(),
                                                                 
                                                                 modelArg.get_mtry(),
                                                                 modelArg.get_lot(),
                                                                 
                                                                 modelArg.get_baseLearn(),
                                                                 
                                                                 modelArg.get_vtry(),
                                                                 modelArg.get_vtryArray(),
                                                                 null,
                                                                 
                                                                 modelArg.get_ytry(),

                                                                 modelArg.get_nodeSize(),
                                                                 modelArg.get_nodeDepth(),

                                                                 modelArg.get_eventCount(),
                                                                 modelArg.get_eventWeight(),

                                                                 modelArg.get_ntree(),
                                                                 modelArg.get_nSize(),
                                                                 
                                                                 modelArg.get_yInfo(),
                                                                 modelArg.get_yLevels(),
                                                                 modelArg.get_yData(),

                                                                 modelArg.get_xInfo(),
                                                                 modelArg.get_xLevels(),
                                                                 modelArg.get_xData(),

                                                                 modelArg.get_sampleInfo(),

                                                                 modelArg.get_xWeightStat(),
                                                                 modelArg.get_yWeight(),
                                                                 modelArg.get_xWeight(),

                                                                 modelArg.get_timeInterest(),

                                                                 modelArg.get_nImpute(),
                                                                 
                                                                 modelArg.get_blockSize(),

                                                                 modelArg.get_quantile(),

                                                                 0, // Pre Sort Disabled
                                                                 
                                                                 modelArg.get_rfCores());

            
        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {        
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Native.grow() nominal exit.");
        @RF_TRACE_OFF@  }

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "Train (original) Ensemble LinkedHashMap of length " + ensembleList.size());
        @RF_TRACE_OFF@    RandomForestModel temporaryModel = new RandomForestModel(modelArg, ensembleList);
        @RF_TRACE_OFF@    temporaryModel.printEnsembleList();
        @RF_TRACE_OFF@  }
        
        // We trim, in particular, treeID, nodeID, parmID, contPT, mwcpSZ, and mwcpPT (iff necessary).

        Ensemble ensb;

        int trimmedPrimarySize = 0;

        int[] trimmedFactorSize;

        int[][] mwcpCT;

        mwcpCT = null;
        
        ensb = (Ensemble) ensembleList.get("leafCount");
        int[] leafCount = (int[]) ensb.ensembleVector;
        
        if (modelArg.get_hdim() == 0) {
            trimmedFactorSize = new int[1];
            trimmedFactorSize[0] = 0;

            mwcpCT = new int[1][]; 
        }
        else {
            trimmedFactorSize = new int[modelArg.get_hdim()];
            for (int k = 0; k < modelArg.get_hdim(); k++) {
                trimmedFactorSize[k] = 0;
            }

            mwcpCT = new int[modelArg.get_hdim()][]; 
        }

        Set keys = ensembleList.keySet();

        Iterator pivot = keys.iterator();

        boolean pivotFound = false;

        while(!pivotFound) {

            ensb = (Ensemble) ensembleList.get(pivot.next());
            
            if (ensb.name.equals("treeID")) {
                pivotFound = true;
            }
        }

        ensb = (Ensemble) ensembleList.get("mwcpCT");
        mwcpCT[0] = (int[]) ensb.ensembleVector;


        for (int b = 0; b < modelArg.get_ntree(); b++) {
            if (leafCount[b] > 0) {
                // The tree was not rejected.

                // Count the number of internal and external
                // (terminal) nodes in the forest.
                trimmedPrimarySize = trimmedPrimarySize + (2 * leafCount[b]) - 1;

                // Count the total number of mwcp words in the forest.
                trimmedFactorSize[0] += mwcpCT[0][b];
               
            }
            else {
                // The tree was rejected.  However, it acts as a
                // placeholder, being a stump topologically and thus
                // adds to the total node count.
                trimmedPrimarySize += 1;
            }
        }

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "Trimmed primary size:   " + trimmedPrimarySize);
        @RF_TRACE_OFF@  }
        
        ensb = (Ensemble) ensembleList.get("treeID");        
        ensb.ensembleVector = Arrays.copyOf((int[]) ensb.ensembleVector, trimmedPrimarySize);
        ensb.size = trimmedPrimarySize;
        
        ensb = (Ensemble) ensembleList.get("nodeID");        
        ensb.ensembleVector = Arrays.copyOf((int[]) ensb.ensembleVector, trimmedPrimarySize);
        ensb.size = trimmedPrimarySize;

        // Non-Greedy objects.
        ensb = (Ensemble) ensembleList.get("parmID");        
        ensb.ensembleVector = Arrays.copyOf((int[]) ensb.ensembleVector, trimmedPrimarySize);
        ensb.size = trimmedPrimarySize;
        
        ensb = (Ensemble) ensembleList.get("contPT");        
        ensb.ensembleVector = Arrays.copyOf((double[]) ensb.ensembleVector, trimmedPrimarySize);
        ensb.size = trimmedPrimarySize;

        ensb = (Ensemble) ensembleList.get("mwcpSZ");
        ensb.ensembleVector = Arrays.copyOf((int[]) ensb.ensembleVector, trimmedPrimarySize);
        ensb.size = trimmedPrimarySize;

        ensb = (Ensemble) ensembleList.get("mwcpPT");        
        if (trimmedFactorSize[0] > 0) {
            ensb.ensembleVector = Arrays.copyOfRange((int[]) ensb.ensembleVector, 0, trimmedFactorSize[0]);
            ensb.size = trimmedFactorSize[0];
        }
        else {
            // mwcpPT will be a vector of nominal length zero (0) with meta info of size zero (0).  Leave it as is.
        }

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "Train (trimmed) Ensemble LinkedHashMap of length " + ensembleList.size());
        @RF_TRACE_OFF@    RandomForestModel temporaryModel = new RandomForestModel(modelArg, ensembleList);
        @RF_TRACE_OFF@    temporaryModel.printEnsembleList();
        @RF_TRACE_OFF@  }
        
        // Create the Random Forest Model Object, given the inputs, and resulting outputs.
        RandomForestModel rfModel = new RandomForestModel(modelArg, ensembleList);

        return rfModel;
    }

    public static RandomForestModel predict(RandomForestModel model)  {
        TestModelArg testModelArg = new TestModelArg(model.modelArg);

        return predict(model, testModelArg);

    }



    public static RandomForestModel predict(RandomForestModel model, TestModelArg testModelArg)  {

        HCzero hc_zero;
        
        // Only models resulting from a RandomForest.train() call are allowed.
        if (model.modelType != ModelType.GROW) {
            RFLogger.log(Level.SEVERE, "Model is not a result of a RandomForest.train() call.");
            RFLogger.log(Level.SEVERE, "Incorrect ModelType:  " + model.modelType);            

            throw new IllegalArgumentException();
        }

        // Shortcut to modelArg object.
        ModelArg modelArg = model.modelArg;

        hc_zero = new HCzero((int[]) (model.getEnsembleObj("parmID")).ensembleVector,
                             (double[]) (model.getEnsembleObj("contPT")).ensembleVector,
                             (int[]) (model.getEnsembleObj("mwcpSZ")).ensembleVector,
                             (int[]) (model.getEnsembleObj("mwcpPT")).ensembleVector);

        
        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "Family restore:  " + modelArg.get_family());
        @RF_TRACE_OFF@  }

        LinkedHashMap   ensembleList = Native.getInstance().predict(modelArg.get_trace(),
                                                                    testModelArg.get_seed(),

                                                                    modelArg.getEnsembleArgOptLow() + modelArg.getModelArgOptLow(),
                                                                    modelArg.getEnsembleArgOptHigh() + modelArg.getModelArgOptHigh(),

                                                                    // >>>> start of maxi forest object >>>>
                                                                    modelArg.get_ntree(),
                                                                    modelArg.get_nSize(),

                                                                    modelArg.get_yInfo(),
                                                                    modelArg.get_yLevels(),
                                                                    modelArg.get_yData(),

                                                                    modelArg.get_xInfo(),
                                                                    modelArg.get_xLevels(),
                                                                    modelArg.get_xData(),
                                                              
                                                                    modelArg.get_sampleInfo(),

                                                                    modelArg.get_timeInterest(),

                                                                    ((int[]) (model.getEnsembleObj("treeID")).ensembleVector).length, // This is the total node count.
                                                                    (int[]) (model.getEnsembleObj("leafCount")).ensembleVector,  

                                                                    new SeedInfo((int[]) model.getEnsembleObj("seed").ensembleVector, null, 0),

                                                                    modelArg.get_hdim(),

                                                                    modelArg.get_baseLearn(),
                                                                    
                                                                    (int[]) (model.getEnsembleObj("treeID")).ensembleVector,
                                                                    (int[]) (model.getEnsembleObj("nodeID")).ensembleVector,
                                                              
                                                                    hc_zero,
                                                                    new HConeAugIntr(null, null, null),
                                                                    new HConeAugSyth(null, null),
                                                                    new HCone(null, null),

                                                                    new HCparmID(null),
                                                                    new HCcontPT(null),
                                                                    new HCcontPTR(null),
                                                                    new HCmwcpSZ(null),
                                                                    new HCmwcpPT(null),
                                                                    new HCaugmXone(null),
                                                                    new HCaugmXtwo(null),
                                                                    new HCaugmXS(null),
                                                                    new HCaugmSythTop(null, null, null, null, null, null, null, null, null),

                                                                    (model.getEnsembleObj("rmbrMembership") != null) ? (int[]) (model.getEnsembleObj("rmbrMembership")).ensembleVector : null,  // TBD TBD rename
                                                                    (model.getEnsembleObj("ambrMembership") != null) ? (int[]) (model.getEnsembleObj("ambrMembership")).ensembleVector : null,  // TBD TBD rename
                                                                    (model.getEnsembleObj("tnRCNT") != null) ? (int[]) (model.getEnsembleObj("tnRCNT")).ensembleVector : null,
                                                                    (model.getEnsembleObj("tnACNT") != null) ? (int[]) (model.getEnsembleObj("tnACNT")).ensembleVector : null,
                                                              
                                                                    (model.getEnsembleObj("tnSURV") != null) ? (double[]) (model.getEnsembleObj("tnSURV")).ensembleVector : null,
                                                                    (model.getEnsembleObj("tnMORT") != null) ? (double[]) (model.getEnsembleObj("tnMORT")).ensembleVector : null,
                                                                    (model.getEnsembleObj("tnNLSN") != null) ? (double[]) (model.getEnsembleObj("tnNLSN")).ensembleVector : null,
                                                                    (model.getEnsembleObj("tnCSHZ") != null) ? (double[]) (model.getEnsembleObj("tnCSHZ")).ensembleVector : null,
                                                                    (model.getEnsembleObj("tnCIFN") != null) ? (double[]) (model.getEnsembleObj("tnCIFN")).ensembleVector : null,
                                                                    (model.getEnsembleObj("tnREGR") != null) ? (double[]) (model.getEnsembleObj("tnREGR")).ensembleVector : null,
                                                                    (model.getEnsembleObj("tnCLAS") != null) ? (int[])    (model.getEnsembleObj("tnCLAS")).ensembleVector : null,
                                                                    // <<<< end of maxi forest object <<<<

                                                                    testModelArg.get_yTargetIndex(),  // Default is to send in requests for all responses.
                                                              
                                                                    0,

                                                                    testModelArg.get_xMarginal(),

                                                                    testModelArg.get_xImportance(),
                                                                    
                                                                    testModelArg.get_partial(),

                                                                    testModelArg.get_fnSize(),
                                                                    testModelArg.get_fySize(),
                                                                    testModelArg.get_fyData(),
                                                                    testModelArg.get_fxData(),

                                                                    testModelArg.get_blockSize(),

                                                                    testModelArg.get_quantile(),

                                                                    testModelArg.get_getTree(),

                                                                    testModelArg.get_rfCores()); 

            
        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {        
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Native.predict() nominal exit.");
        @RF_TRACE_OFF@  }

        // Create the Random Forest Model Object, given the inputs, and resulting outputs.
        RandomForestModel rfModel = new RandomForestModel(testModelArg, ensembleList);
        
        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "Test (restore) Ensemble LinkedHashMap of length " + ensembleList.size());
        @RF_TRACE_OFF@    rfModel.printEnsembleList();
        @RF_TRACE_OFF@  }
       

        return rfModel;
    }

    public static RandomForestModel predict(ModelArg modelArg, Dataset dataset)  {

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Family predict:  " + modelArg.get_family());
        @RF_TRACE_OFF@  }
        
        return null;
    }
    
    
}
