package com.kogalur.randomforest;

import com.kogalur.randomforest.NativeOpt;
import com.kogalur.randomforest.RFLogger;
import java.util.logging.Level;

import java.util.List;




class EnsembleArg {

    private String mode;
    private String family;
    
    private java.util.HashMap <String, String> ensembleArg;
    private java.util.HashMap <String, NativeOpt> ensembleList;
    
    EnsembleArg(String mode, String family) {

        this.mode = mode;
        this.family = family;

        set();
    }

    void set() {
        
        ensembleArg = new java.util.HashMap <String, String> (16);
        ensembleList = new java.util.HashMap <String, NativeOpt> (16);
        

        
        if (mode.equals("grow")) {
            ensembleArg.put("forest", "yes");
            ensembleList.put("forest",
                             new NativeOpt("forest",
                                           new String[] {"yes"},
                                           new int[] {(1 << 5)}));
        }
        else {
            
            ensembleArg.put("forest", "no");
            ensembleList.put("forest",
                             new NativeOpt("forest",
                                           new String[] {"no"},
                                           new int[] {(0)}));
        }
            
        if (mode.equals("grow") || mode.equals("rest")) {

            

            ensembleArg.put("weight", "no");
            ensembleList.put("weight",
                              new NativeOpt("weight",
                                            new String[] {"no", "inbag", "oob", "all"},
                                            new int[] {0,
                                                       (1 << 0) + (1 << 1),
                                                       (1 << 0) + (1 << 2), 
                                                       (1 << 0) + (1 << 21) + (1 << 22)}));


            ensembleArg.put("proximity", "no");
            ensembleList.put("proximity",
                              new NativeOpt("proximity",
                                            new String[] {"no", "inbag", "oob", "all"},
                                            new int[] {0,
                                                       (1 << 28) + (1 << 29),
                                                       (1 << 28) + (1 << 29), 
                                                       (1 << 28) + (1 << 29) + (1 << 30)}));

        }
        else {

            
                
            ensembleArg.put("weight", "no");
            ensembleList.put("weight",
                              new NativeOpt("weight",
                                            new String[] {"no", "all"},
                                            new int[] {0,
                                                       (1 << 0) + (1 << 21) + (1 << 22)}));


            ensembleArg.put("proximity", "no");
            ensembleList.put("proximity",
                              new NativeOpt("proximity",
                                            new String[] {"no", "all"},
                                            new int[] {0,
                                                       (1 << 28) + (1 << 29) + (1 << 30)}));
        }

        ensembleArg.put("membership", "no");
        ensembleList.put("membership",
                          new NativeOpt("membership",
                                        new String[] {"no", "yes"},
                                        new int[] {0, (1 << 6)}));

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
        
        ensembleArg.put("error", "last.tree");
        ensembleList.put("error",
                          new NativeOpt("error",
                                        new String[] {"last.tree", "per.tree"},
                                        new int[] {0, (1 << 13)}));

        ensembleArg.put("varUsed", "no");
        ensembleList.put("varUsed",
                          new NativeOpt("varUsed",
                                        new String[] {"no", "per.tree", "sum.tree"},
                                        new int[] {0, (1 << 13), (1 << 12)}));

        ensembleArg.put("splitDepth", "no");
        ensembleList.put("splitDepth",
                          new NativeOpt("splitDepth",
                                        new String[] {"no", "per.tree", "sum.tree"},
                                        new int[] {0, (1 << 23), (1 << 22)}));

        if (mode.equals("grow")) {
                ensembleArg.put("qualitativeTerminalInfo", "yes");
                ensembleList.put("qualitativeTerminalInfo",
                                 new NativeOpt("qualitativeTerminalInfo",
                                               new String[] {"no", "yes"},
                                               new int[] {0, (1 << 16)}));
        }

        
        if (family.equals("RF-C") || family.equals("RF-C+")) {
                ensembleArg.put("errorType", "default");
                ensembleList.put("errorType",
                                 new NativeOpt("errorType",
                                               new String[] {"default", "brier", "g.mean", "g.mean.drc"},
                                               new int[] {(1 << 2) + (0), (1 << 2) + (1 << 3), (1 << 2) + (1 << 14), (1 << 2) + (1 << 15)}));
        }
        else {
            ensembleArg.put("errorType", "default");
            ensembleList.put("errorType",
                             new NativeOpt("errorType",
                                           new String[] {"default"},
                                           new int[] {(1 << 2) + (0)}));
        }
        
    }

    void set(String key, String value) {
        if (ensembleList.containsKey(key)) {

            java.util.HashMap nativeOptMap = ensembleList.get(key).option;
            if (nativeOptMap.containsKey(value)) {
                ensembleArg.put(key, value);
            }
            else {
                RFLogger.log(Level.SEVERE, "Unknown ensemble value in <key, value>:  " + "< " + key + ", " + value + " >");
                throw new IllegalArgumentException();
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

        
        String value = ensembleArg.get(key);

        
        NativeOpt nativeOpt = ensembleList.get(key);

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Native Option <key, value> = < " + key + ", " + value + " > == " + nativeOpt.get(value));        
        @RF_TRACE_OFF@  }
        
        
        return nativeOpt.get(value);
    }


}
