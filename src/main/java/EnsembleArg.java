package com.kogalur.randomforest;

import com.kogalur.randomforest.NativeOpt;
import com.kogalur.randomforest.RFLogger;
import java.util.logging.Level;

import java.util.List;
//import java.util.Iterator;
//import java.util.Arrays;
//import java.util.ArrayList;

class EnsembleArg {

    private String mode;
    private String family;
    private EnsembleArg growEnsembleArg;
    
    private java.util.HashMap <String, String> ensembleArg;
    private java.util.HashMap <String, NativeOpt> ensembleList;

    // Nominal GROW mode constructor.
    EnsembleArg(String mode, String family) {

        this.mode = mode;
        this.family = family;

        // Safe the grow-side ensemble argument object as it is irrelevant.
        growEnsembleArg = null;
        
        set();
    }

    // Nominal !GROW mode constructor.
    EnsembleArg(EnsembleArg growEnsembleArg, String mode, String family) {

        this.mode = mode;
        this.family = family;

        // Initialize the grow-side ensemble argument object. We access
        // potentially existing qualt and quant outgoing options and
        // convey these as incoming options to the native code.
        this.growEnsembleArg = growEnsembleArg;
        
        // Safe the grow-side ensemble argument object as it is irrelevant.
        growEnsembleArg = null;
        
        set();
    }

    void set() {
        
        ensembleArg = new java.util.HashMap <String, String> (16);
        ensembleList = new java.util.HashMap <String, NativeOpt> (16);

        if (mode.equals("grow")) {
            // Always output "forest" in grow mode. Always. 
            ensembleArg.put("forest", "yes");
            ensembleList.put("forest",
                             new NativeOpt("forest",
                                           new String[] {"yes"},
                                           new int[] {(1 << 5)}));

            
            // Default action is to output qualts in grow mode.
            ensembleArg.put("qualitativeTerminalInfo", "yes");
            ensembleList.put("qualitativeTerminalInfo",
                             new NativeOpt("qualitativeTerminalInfo",
                                           new String[] {"no", "yes"},
                                           new int[] {0, (1 << 16)}));
            
            // Default action is to not output quants in grow mode.
            ensembleArg.put("quantitativeTerminalInfo", "no");
            ensembleList.put("quantitativeTerminalInfo",
                             new NativeOpt("quantitativeTerminalInfo",
                                           new String[] {"no", "yes"},
                                           new int[] {0, (1 << 18)}));
            
        }
        else {
            // Never output "forest" in !grow mode. Never. 
            ensembleArg.put("forest", "no");
            ensembleList.put("forest",
                             new NativeOpt("forest",
                                           new String[] {"no"},
                                           new int[] {(0)}));

            // Default action is to assume no incoming qualts in !grow mode.
            ensembleArg.put("qualitativeTerminalInfo", "no");
            ensembleList.put("qualitativeTerminalInfo",
                             new NativeOpt("qualitativeTerminalInfo",
                                           new String[] {"no", "yes"},
                                           new int[] {0, (1 << 17)}));
            
            // Default action is to assume no incoming quants in !grow mode.
            ensembleArg.put("quantitativeTerminalInfo", "no");
            ensembleList.put("quantitativeTerminalInfo",
                             new NativeOpt("quantitativeTerminalInfo",
                                           new String[] {"no", "yes"},
                                           new int[] {0, (1 << 19)}));

            // !GROW mode allows the presence of a grow-side ensemble argument object.
            if (growEnsembleArg != null) {
                // Override the qualts and quants accordingly.
                ensembleArg.put("qualitativeTerminalInfo",  growEnsembleArg.get("qualitativeTerminalInfo"));
                ensembleArg.put("quantitativeTerminalInfo", growEnsembleArg.get("quantitativeTerminalInfo"));
            }
        }
            
        if (mode.equals("grow") || mode.equals("rest")) {

             

            ensembleArg.put("weight", "no");
            ensembleList.put("weight",
                              new NativeOpt("weight",
                                            new String[] {"no", "inbag", "oob"},
                                            new int[] {0,
                                                       (1 << 0) + (1 << 1),
                                                       (1 << 0) + (1 << 2)}));


            ensembleArg.put("proximity", "no");
                             ensembleList.put("proximity",
                              new NativeOpt("proximity",
                                            new String[] {"no", "inbag", "oob"},
                                            new int[] {0,
                                                       (1 << 28) + (1 << 29),
                                                       (1 << 28) + (1 << 30)}));

            ensembleArg.put("importance", "no");
            ensembleList.put("importance",
                             new NativeOpt("importance",
                                           new String[] {"no",
                                                         "permute.ensemble", "random.ensemble",
                                                         "permute", "random"},
                                           new int[] {0,
                                                      (1 << 25) + (1 << 8), (1 << 25) + (1 << 9),
                                                      (1 << 25) + (1 << 24) + (1 << 8), (1 << 25) + (1 << 24) + (1 << 9)}));
        }
        else {

             
                
            ensembleArg.put("weight", "no");
            ensembleList.put("weight",
                              new NativeOpt("weight",
                                            new String[] {"no", "all"},
                                            new int[] {0,
                                                       (1 << 0) + (1 << 1) +(1 << 2)}));


            ensembleArg.put("proximity", "no");
            ensembleList.put("proximity",
                              new NativeOpt("proximity",
                                            new String[] {"no", "all"},
                                            new int[] {0,
                                                       (1 << 28) + (1 << 29) + (1 << 30)}));

            ensembleArg.put("importance", "no");
            ensembleList.put("importance",
                             new NativeOpt("importance",
                                           new String[] {"no",
                                                         "anti.ensemble", "permute.ensemble", "random.ensemble",
                                                         "anti.joint.ensemble", "permute.joint.ensemble", "random.joint.ensemble",
                                                         "anti", "permute", "random",
                                                         "anti.joint", "permute.joint", "random.joint"},
                                           new int[] {0,
                                                      (1 << 25) + (0), (1 << 25) + (1 << 8), (1 << 25) + (1 << 9),
                                                      (1 << 25) + (1 << 10) + (0), (1 << 25) + (1 << 10) + (1 << 8), (1 << 25) + (1 << 10) + (1 << 9),
                                                      (1 << 25) + (1 << 24) + (0), (1 << 25) + (1 << 24) + (1 << 8), (1 << 25) + (1 << 24) + (1 << 9),
                                                      (1 << 25) + (1 << 24) + (1 << 10) + (0), (1 << 25) + (1 << 24) + (1 << 10) + (1 << 8), (1 << 25) + (1 << 24) + (1 << 10) + (1 << 9)}));

        }

        ensembleArg.put("membership", "no");
        ensembleList.put("membership",
                          new NativeOpt("membership",
                                        new String[] {"no", "yes"},
                                        new int[] {0, (1 << 6)}));

        // This option is subordinate to errorType, and is ignored if errorType is suppressed.
        ensembleArg.put("error", "last.tree");
        ensembleList.put("error",
                          new NativeOpt("error",
                                        new String[] {"last.tree", "every.tree"},
                                        new int[] {0, (1 << 13)}));

        ensembleArg.put("varUsed", "no");
        ensembleList.put("varUsed",
                          new NativeOpt("varUsed",
                                        new String[] {"no", "every.tree", "sum.tree"},
                                        new int[] {0, (1 << 13), (1 << 12)}));

        ensembleArg.put("splitDepth", "no");
        ensembleList.put("splitDepth",
                          new NativeOpt("splitDepth",
                                        new String[] {"no", "every.tree", "sum.tree"},
                                        new int[] {0, (1 << 23), (1 << 22)}));

        
        if (family.equals("RF-C") || family.equals("RF-C+")) {
                ensembleArg.put("errorType", "misclass");
                ensembleList.put("errorType",
                                 new NativeOpt("errorType",
                                               new String[] {"no", "misclass", "brier", "g.mean"},
                                               new int[] {0, (1 << 2) + (0), (1 << 2) + (1 << 3), (1 << 2) + (1 << 14)}));

                ensembleArg.put("predictionType", "max.vote");
                ensembleList.put("predictionType",
                                 new NativeOpt("predictionType",
                                               new String[] {"max.vote", "rfq"},
                                               new int[] {0, (1 << 14)}));

        }
        else if (family.equals("RF-R") || family.equals("RF-R+")) {
            ensembleArg.put("errorType", "mse");
            ensembleList.put("errorType",
                             new NativeOpt("errorType",
                                           new String[] {"no", "mse"},
                                           new int[] {0, (1 << 2) + (0)}));

            ensembleArg.put("predictionType", "mean");
            ensembleList.put("predictionType",
                             new NativeOpt("predictionType",
                                           new String[] {"mean"},
                                           new int[] {0}));

        }
        else if (family.equals("RF-S")) {
            ensembleArg.put("errorType", "c-index");
            ensembleList.put("errorType",
                             new NativeOpt("errorType",
                                           new String[] {"no", "c-index"},
                                           new int[] {0, (1 << 2) + (0)}));

            ensembleArg.put("predictionType", "default");
            ensembleList.put("predictionType",
                             new NativeOpt("predictionType",
                                           new String[] {"default"},
                                           new int[] {0}));

        }
        else if (family.equals("RF-M+")) {
            ensembleArg.put("errorType", "default");
            ensembleList.put("errorType",
                             new NativeOpt("errorType",
                                           new String[] {"no", "default"},
                                           new int[] {0, (1 << 2) + (0)}));

            ensembleArg.put("predictionType", "default");
            ensembleList.put("predictionType",
                             new NativeOpt("predictionType",
                                           new String[] {"default"},
                                           new int[] {0}));

        }
        else if (family.equals("RF-U")) {
            ensembleArg.put("errorType", "no");
            ensembleList.put("errorType",
                             new NativeOpt("errorType",
                                           new String[] {"no"},
                                           new int[] {0}));

            ensembleArg.put("predictionType", "no");
            ensembleList.put("predictionType",
                             new NativeOpt("predictionType",
                                           new String[] {"no"},
                                           new int[] {0}));

        }
        
    }

    void set(String key, String value) {
        if (ensembleList.containsKey(key)) {

            boolean setFlag = true;
        
            // The following keys have restrictions.
            if (key.equals("forest")) {
                RFLogger.log(Level.WARNING, "forest ensemble cannot be modified.");
                setFlag = false;
            }
            else if (key.equals("qualitativeTerminalInfo")) {
                if (!mode.equals("grow")) {
                    RFLogger.log(Level.WARNING, "qualitativeTerminalInfo ensemble cannot be modified.");
                    setFlag = false;
                }
            }
            else if (key.equals("quantitativeTerminalInfo")) {
                if (!mode.equals("grow")) {
                    RFLogger.log(Level.WARNING, "quantitativeTerminalInfo ensemble cannot be modified.");
                    setFlag = false;
                }
            }
            if (setFlag) {
                java.util.HashMap nativeOptMap = ensembleList.get(key).option;
                if (nativeOptMap.containsKey(value)) {
                    ensembleArg.put(key, value);
                }
                else {
                    RFLogger.log(Level.SEVERE, "Unknown ensemble value in <key, value>:  " + "< " + key + ", " + value + " >");
                    throw new IllegalArgumentException();
                }
            }
        }
        else {
            RFLogger.log(Level.SEVERE, "Unknown ensemble key in <key, value>:  " + "< " + key + ", " + value + " >");
            throw new IllegalArgumentException();
        }
    }

    String get(String key) {
        if (!ensembleArg.containsKey(key)) {
            RFLogger.log(Level.SEVERE, "Unknown ensemble key:  " + key);
            throw new IllegalArgumentException();
        }
        return ensembleArg.get(key);
    }
    
    int getNative(String key) {

        // Get the user value associated with this key.
        String value = ensembleArg.get(key);

        // Get the native option object for this key.
        NativeOpt nativeOpt = ensembleList.get(key);

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Native Option <key, value> = < " + key + ", " + value + " > == " + nativeOpt.get(value));        
        @RF_TRACE_OFF@  }
        
        // Get the bit equivalent for this key.
        return nativeOpt.get(value);
    }


}
