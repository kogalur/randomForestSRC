package com.kogalur.randomforest;

import com.kogalur.randomforest.ModelArg;
import com.kogalur.randomforest.RandomForestModel;

import com.kogalur.randomforest.Native;
import com.kogalur.randomforest.RFLogger;

import java.util.ArrayList;

import org.apache.spark.sql.Dataset;

import java.util.logging.Level;

public class RandomForest {

    private RandomForest() {

    }
    
    public static RandomForestModel train(ModelArg modelArg)  {

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Family train:  " + modelArg.get_family());
        @RF_TRACE_OFF@  }
        
        ArrayList ensembleArray = Native.getInstance().grow(
                                                            modelArg.get_trace(),
                                                            modelArg.get_seed(),

                                                            modelArg.getEnsembleArgOptLow() + modelArg.getModelArgOptLow(),
                                                            modelArg.getEnsembleArgOptHigh() + modelArg.getModelArgOptHigh(),

                                                            modelArg.getSplitRuleID(modelArg.get_splitRule()),
                                                            modelArg.get_nSplit(),

                                                            modelArg.get_mtry(),
                                                            modelArg.get_ytry(),

                                                            modelArg.get_nodeSize(),
                                                            modelArg.get_nodeDepth(),

                                                            modelArg.get_eventWeightSize(),
                                                            modelArg.get_eventWeight(),

                                                            modelArg.get_ntree(),
                                                            modelArg.get_nSize(),
                                                            modelArg.get_ySize(),
                                                            modelArg.get_yType(),
                                                            modelArg.get_yLevel(),
                                                            modelArg.get_yData(),

                                                            modelArg.get_xSize(),
                                                            modelArg.get_xType(),
                                                            modelArg.get_xLevel(),

                                                            modelArg.get_sampleSize(),
                                                            modelArg.get_sample(),
                                                            modelArg.get_caseWeight(),

                                                            modelArg.get_xSplitStatWt(),
                                                            modelArg.get_yWeight(),

                                                            modelArg.get_xWeight(),
                                                            modelArg.get_xData(),

                                                            modelArg.get_timeInterestSize(),
                                                            modelArg.get_timeInterest(),
                                                            modelArg.get_nImpute(),
                                                            modelArg.get_rfCores());


        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {        
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Native.grow() nominal exit.");
        @RF_TRACE_OFF@  }

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Ensemble ArrayList:  of length " + ensembleArray.size());
        @RF_TRACE_OFF@  for (int i = 0; i < ensembleArray.size(); i++) {
        @RF_TRACE_OFF@    Ensemble ensemble = (Ensemble) ensembleArray.get(i);
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "[" + i + "] = " + ensemble.name);
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "              " + "       type  : " + String.valueOf(ensemble.type));
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "              " + "    identity : " + ensemble.identity);
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "              " + "        size : " + ensemble.size);
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "              " + "      ragged : " + ensemble.ragged);
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "              " + "  auxDimSize : " + ensemble.auxDimSize);
        @RF_TRACE_OFF@  }
        @RF_TRACE_OFF@  }
       
        
        return new RandomForestModel();
    }


    public static RandomForestModel predict(ModelArg modelArg)  {

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Family predict:  " + modelArg.get_family());
        @RF_TRACE_OFF@  }
        
        return new RandomForestModel();
    }

    public static RandomForestModel predict(ModelArg modelArg, Dataset dataset)  {

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Family predict:  " + modelArg.get_family());
        @RF_TRACE_OFF@  }
        
        return new RandomForestModel();
    }
    
    
}
