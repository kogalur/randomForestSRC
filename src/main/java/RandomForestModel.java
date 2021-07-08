package com.kogalur.randomforest;

import com.kogalur.randomforest.RFLogger;

import java.lang.Math;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.LinkedHashMap;
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

    public ModelType modelType;
    
    // Model argument, ensemble list contained in this model.
    public ModelArg modelArg;
    public TestModelArg testModelArg;
    public LinkedHashMap ensembleList;


    // *****************************************************************
    // REST / PRED Ensemble Requests.
    // *****************************************************************
    private EnsembleArg ensembleArg;
    
    // Constructor used by RandomForest.train(..) to create GROW model.
    RandomForestModel(ModelArg modelArg, LinkedHashMap ensembleList) {

        this.modelArg = modelArg;
        this.testModelArg = null;
        
        this.ensembleList = ensembleList;

        modelType = ModelType.GROW;

    }


    // Constructor used by RandomForest.predict(..) to create !GROW model.
    RandomForestModel(TestModelArg testModelArg, LinkedHashMap ensembleList) {

        this.modelArg = null;
        this.testModelArg = testModelArg;
        this.ensembleList = ensembleList;

        if (testModelArg.getMode() == "rest") {
            this.modelType = ModelType.REST;
        }
        else {
            this.modelType = ModelType.PRED;
        }
    }

    // TBD TBD remove public option TBD TBD
    public Ensemble getEnsembleObj(String name) {
        Ensemble ensb = (Ensemble) ensembleList.get(name);
        return ensb;
    }


    // TBD TBD remove public option TBD TBD    
    public void printEnsembleList() {
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
