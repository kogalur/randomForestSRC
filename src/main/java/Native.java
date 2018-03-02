package com.kogalur.randomforest;

import com.kogalur.randomforest.RFLogger;

import java.io.File;
import java.util.HashMap;

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

    native HashMap   grow(int          traceFlag,
                          int          seedDynamic,

                          int          optLow,
                          int          optHigh,

                          int          splitRule,
                          int          nSplit,

                          int          mtry,
                          int          htry,
                          int          ytry,

                          int          nodeSize,
                          int          nodeDepth,

                          int          eventCount,
                          double[]     eventWeight,
                            
                          int          ntree,
                          int          nSize,
                          int          ySize,
                          char[]       yType,
                          int[]        yLevel,
                          double[][]   yData,

                          int          xSize,
                          char[]       xType,
                          int[]        xLevel,

                          int          sampleSize,
                          int[][]      sample,
                          double[]     caseWeight,

                          double[]     xStatisticalWeight,
                          double[]     yWeight,

                          double[]     xWeight,
                          double[][]   xData,

                          int          timeInterestSize,
                          double[]     timeInterest,
                          int          nImpute,
                          int          numThreads);


    native HashMap   predict(int          traceFlag,
                             int          seedDynamic,

                             int          optLow,
                             int          optHigh,

                             // >>>> start of maxi forest object >>>>
                             int          ntree,
                             int          nSize,

                             int          ySize,
                             char[]       yType,
                             int[]        yLevel,
                             double[][]   yData,

                             int          xSize,
                             char[]       xType,
                             int[]        xLevel,
                             double[][]   xData,
                             
                             int          sampleSize,
                             int[][]      sample,
                             double[]     caseWeight,

                             int          timeInterestSize,
                             double[]     timeInterest,

                             int[]        seed,
                             int          totalNodeCount,

                             int[]        treeID,
                             int[]        nodeID,

                             int          htry,
                             
                             int[]        hcDim,
                             int[]        hcPartDim,
                             int[]        hcPartIdx,
                             int[]        osPartIdx,
                                                                       
                             int[]        parmID,
                             double[]     contPT,
                             int[]        mwcpSZ,
                             int[]        mwcpPT,
                                                                       
                             int[]        parmID2,
                             double[]     contPT2,
                             int[]        mwcpSZ2,
                             int[]        mwcpPT2,

                             int[]        parmID3,
                             double[]     contPT3,
                             int[]        mwcpSZ3,
                             int[]        mwcpPT3,

                             int[]        parmID4,
                             double[]     contPT4,
                             int[]        mwcpSZ4,
                             int[]        mwcpPT4,

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

                             int          yTargetSize,
                             int[]        yTargetIndex,
                             
                             int          ptnCount,
                             
                                             

                             int          xImportanceSize,
                             int[]        xImportanceIndex,
                             
                             int          xPartialType,
                             int          xPartialIndex,
                             int          xPartialSize,
                             double[]     xPartialValue,
                             int          x2PartialSize,
                             int[]        x2PartialIndex,
                             double[]     x2PartialValue,

                             int          subsetSize,
                             int[]        subsetIndex,

                             int          fnSize,
                             int          fySize,
                             double[]     fyData,
                             double[]     fxData,

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
