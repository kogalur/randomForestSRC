package com.kogalur.randomforest;

import com.kogalur.randomforest.RFLogger;

import java.io.File;
import java.util.LinkedHashMap;

import java.util.logging.Level;
import java.io.IOException;
import java.lang.SecurityException;

class Native {

    private static Native myInstance;

    static {
        File file = new File(".");

        String libraryLocation;
        String libraryName;
        String opsSwitch;
        String sep;
        boolean libraryFlag;

        try {

            libraryFlag = true;
            libraryLocation = file.getCanonicalPath();
            sep = File.separator;

            opsSwitch = System.getProperty("os.name");
            if (opsSwitch.startsWith("Windows")) {
                libraryName = "@PROJECT_PACKAGE_NAME@.dll";
            }
            else if (opsSwitch.startsWith("Mac OS X")) {
                libraryName = "lib@PROJECT_PACKAGE_NAME@.jnilib";
            }
            else if (opsSwitch.startsWith("Linux")) {
                libraryName = "lib@PROJECT_PACKAGE_NAME@.so";
            }
            else {
                libraryFlag = false;
                RFLogger.log(Level.SEVERE, "Unknown OS in determining shared object library for @PROJECT_PACKAGE_NAME@: " + opsSwitch);
                throw new UnsupportedOperationException();
            }

                
            if(libraryFlag == true) { 
                libraryName = libraryLocation + sep + libraryName;
                System.load(libraryName);

                @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
                @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "\n" + "Location of native library:  " + libraryName + "\n");
                @RF_TRACE_OFF@  }
                
            }

        }
        catch (IOException | SecurityException e) {
            RFLogger.log(Level.SEVERE, e.toString());
        }
    }

    native LinkedHashMap   grow(int          traceFlag,
                                int          seedDynamic,

                                int          optLow,
                                int          optHigh,

                                int          splitRule,
                                int          nSplit,

                                int          mtry,
                                Lot          lot,

                                BaseLearn    baseLearn,
                                
                                int          vtry,
                                int[][]      vtryArray,
                                Object       vtryExperimental,
                                
                                int          ytry,

                                int          nodeSize,
                                int          nodeDepth,

                                int          eventCount,
                                double[]     eventWeight,
                            
                                int          ntree,
                                int          nSize,

                                YInfo        yInfo,
                                int[][]      yLevels,
                                double[][]   yData,

                                XInfo        xInfo,
                                int[][]      xLevels,
                                double[][]   xData,

                                SampleInfo   sampleInfo,

                                double[]     xWeightStat,
                                double[]     yWeight,
                                double[]     xWeight,

                                double[]     timeInterest,

                                int          nImpute,

                                int          blockSize,

                                Quantile     quantile,
                                
                                int          xPreSort,

                                int          numThreads);


    native LinkedHashMap   predict(int          traceFlag,
                                   int          seedDynamic,

                                   int          optLow,
                                   int          optHigh,

                                   // >>>> start of maxi forest object >>>>
                                   int          ntree,
                                   int          nSize,

                                   YInfo        yInfo,
                                   int[][]      yLevels,
                                   double[][]   yData,

                                   XInfo        xInfo,
                                   int[][]      xLevels,
                                   double[][]   xData,
                             
                                   SampleInfo   sampleInfo,

                                   double[]     timeInterest,

                                   int          totalNodeCount,
                                   int[]        leafCount,

                                   SeedInfo     seedInfo,
                                   
                                   int          htry,
                             
                                   BaseLearn    baseLearn,

                                   int[]        treeID,
                                   int[]        nodeID,

                                   HCzero       hc_zero,
                                   HConeAugIntr hc_oneAugIntr,
                                   HConeAugSyth hc_oneAugSyth,
                                   HCone        hc_one,
                                                                       
                                   HCparmID      hc_parmID,
                                   HCcontPT      hc_contPT,
                                   HCcontPTR     hc_contPTR,
                                   HCmwcpSZ      hc_mwcpSZ,
                                   HCmwcpPT      hc_mwcpPT,
                                   HCaugmXone    hc_augmXone,
                                   HCaugmXtwo    hc_augmXtwo,
                                   HCaugmXS      hc_augmXS,
                                   HCaugmSythTop hc_augmSythTop,

                                   int[]        tnRMBR,
                                   int[]        tnAMBR,
                                   int[]        tnRCNT,
                                   int[]        tnACNT,

                                   double[]     tnSURV,
                                   double[]     tnMORT,
                                   double[]     tnNLSN,
                                   double[]     tnCSHZ,
                                   double[]     tnCIFN,
                                   double[]     tnREGR,
                                   int[]        tnCLAS,
                                   // <<<< end of maxi forest object <<<<

                                   int[]        yTargetIndex,
                             
                                   int          ptnCount,
                             
                                   int[]        xMarginal,

                                   int[]        xImportance,
                             
                                   Partial      partial,

                                   int          fnSize,
                                   int          fySize,
                                   double[][]   fyData,
                                   double[][]   fxData,

                                   int          blockSize,                                   

                                   Quantile     quantile,

                                   int[]        getTree,
                                                                        
                                   int          numThreads);    

    static Native getInstance() {
        if (myInstance == null) {
            myInstance = new Native();
        }
        return myInstance;
    }


    private static void log(String msg) {
        RFLogger.log(Level.INFO, msg);
    }

    private static void logError(String msg) {
        RFLogger.log(Level.SEVERE, msg);
    }

    private static void logExit() {
        System.exit(1);
    }
}
