package com.kogalur.randomforest;

import com.kogalur.randomforest.RFLogger;
import java.util.logging.Level;

import org.apache.spark.sql.Dataset;
import org.apache.spark.sql.types.DataTypes;
import org.apache.spark.sql.types.DataType;
import org.apache.spark.sql.types.StructField;
import org.apache.spark.sql.types.StructType;

import org.apache.spark.sql.Row;
import org.apache.spark.ml.feature.StringIndexer;

import java.lang.Math;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import java.util.Random;


import java.util.Vector;








public class ModelArg {


    
    
    
    private EnsembleArg ensembleArg;
   
    
    
    
    
    private Random generator;
    

    
    
    
    private int        trace;
    private int        seed;

    private int        nImpute;
    
    private int        nSplit;

    private int        mtry;
    private int        ytry;

    private int        nodeSize;
    private int        nodeDepth;
    private int        crWeightSize;
    private double[]   crWeight;

    private int        ntree;
    private int        nSize;
    private int        ySize;
    private char[]     yType;
    private int[]      yLevel;
    private double[][] yData;
    
    private int        xSize;
    private char[]     xType;
    private int[]      xLevel;
    
    private int        sampleSize;
    private int[][]    sample;
    private double[]   caseWeight;

    private double[]   xSplitStatWt;
    private double[]   yWeight;

    private double[]   xWeight;
    private double[][] xData;

    private int        timeInterestSize;
    private double[]   timeInterest;
    private double[]   ntime;
    
    private int        rfCores;

    
    
    
    
    private String sampleType;
    private String nullSplit;
    private int    customSplitIndex;

    
    
    
    private String splitRule;
    private String bootstrap;


    
    
    

    
    private String formula;
    
    private String[] formulaU;
    
    private String[] formulaX;
    
    private String[] formulaY;

    
    private String family;

    
    private static final HashMap <String, Integer> splitRuleID;
    static
    {
        splitRuleID = new HashMap<String, Integer> (16);
        splitRuleID.put("logrank",               1);
        splitRuleID.put("logrankscore",          2);
        splitRuleID.put("logrankCR",             3);
        splitRuleID.put("logrankACR",            4);
        splitRuleID.put("random",                5);
        splitRuleID.put("mse",                   6);
        splitRuleID.put("mse.unwt",              7);
        splitRuleID.put("mse.hvwt",              8);
        splitRuleID.put("gini",                  9);
        splitRuleID.put("gini.u nwt",           10);
        splitRuleID.put("gini.h vwt",           11);
        splitRuleID.put("unsupv",               12);
        splitRuleID.put("mv.mse",               13);
        splitRuleID.put("mv.gini",              14);
        splitRuleID.put("custom",               15);
        splitRuleID.put("l2.impute",            16);
    }
    
    
    private java.util.HashMap <String, Character> yTypeHash;
    private java.util.HashMap <String, Character> xTypeHash;

    
    private int xFactorCount;
    private int yFactorCount;
    private int tFactorCount;
    
    
    
    HashMap <String, HashMap> immutableMap;

    
    private Dataset dataset;

    
    private String[] columnName;

    
    private TreeSet <Integer> eventType;


    
    
    
    
    
    public ModelArg (String formula,
                      Dataset dataset) {


        this.dataset = dataset;
        this.formula = formula;

        
        columnName = dataset.columns();

        setDefaultModelArg();

        ensembleArg = new EnsembleArg("grow", family);
    }

    
    
    
    

    private void parseFormula() {

        boolean bigRFlag, complementFlag;
        
        formula = formula.replaceAll("\\ ", "");
      
        formulaU = formula.split("~");

        
        
        ytry = 0;
        
        if (formulaU.length != 2) {
            
            System.out.println("\n Unknown Formula Syntax.");
        }

        
        
        family = new String("RF-X");
        
        if (formulaU[0].startsWith("Surv")) {
            
            if ((formulaU[0].indexOf('(') >= 0) && (formulaU[0].indexOf(')') >= 0)) {
                formulaU[0] = formulaU[0].replace("Surv", "");
                formulaU[0] = formulaU[0].replace("(", "");
                formulaU[0] = formulaU[0].replace(")", "");
                formulaY = formulaU[0].split(",");

                if (formulaY.length != 2) {;
                    
                    RFLogger.log(Level.SEVERE, "Bad Survival Formula Syntax:  " + formula);
                    throw new IllegalArgumentException();
                }
            }
            else {
                
                formulaY = null;
                RFLogger.log(Level.SEVERE, "Bad Survival Formula Syntax:  " + formula);
                throw new IllegalArgumentException();
            }
            
            
            family = new String("RF-S");
        }
        else if (formulaU[0].startsWith("Unsupervised")) {
            
            if ((formulaU[0].indexOf('(') >= 0) && (formulaU[0].indexOf(')') >= 0)) {
                formulaU[0] = formulaU[0].replace("Unsupervised", "");
                formulaU[0] = formulaU[0].replace("(", "");
                formulaU[0] = formulaU[0].replace(")", "");

                
                
                
                
                if (formulaU[0].equals("")) {
                    
                    ytry = 1;
                }
                else {
                    
                    ytry = Integer.parseInt(formulaU[0]);
                }

                if (formulaY != null) {
                    
                    formulaY = null;
                    RFLogger.log(Level.SEVERE, "Bad Unsupervised Formula Syntax:  " + formula);
                    throw new IllegalArgumentException();
                }
            }
            else {
                
                formulaY = null;
                RFLogger.log(Level.SEVERE, "Bad Unsupervised Formula Syntax:  " + formula);
                throw new IllegalArgumentException();
            }
            
            family = new String("RF-U");            
        }
        else if (formulaU[0].startsWith("Multivariate")) {
            
            if ((formulaU[0].indexOf('(') >= 0) && (formulaU[0].indexOf(')') >= 0)) {
                formulaU[0] = formulaU[0].replace("Multivariate", "");
                formulaU[0] = formulaU[0].replace("(", "");
                formulaU[0] = formulaU[0].replace(")", "");
                formulaY = formulaU[0].split(",");
            }
            else {
                
                formulaY = null;
                RFLogger.log(Level.SEVERE, "Bad Multivariate Formula Syntax:  " + formula );
                throw new IllegalArgumentException();
            }
        }
        else {
            
            formulaY = new String[1];
            formulaY[0] = formulaU[0];
        }
        


        
        
        
        bigRFlag = false;
        if (formulaY != null) {
            if (formulaY.length > 50) {
                bigRFlag = true;
            }
            ySize = formulaY.length;
        }
        else {
            ySize = 0;
        }
      
        complementFlag = false;
        if (formulaU[1].compareTo(".") == 0) {
            

            complementFlag = true;
          
            if (formulaY != null) {
                
                formulaX = new String[columnName.length - formulaY.length];

                
                List <String> formulaXL = new ArrayList <String> (Arrays.asList(columnName));
                for (int j = 0; j < formulaY.length; j++) {
                    formulaXL.remove(formulaY[j]);
                }
                for (int i = 0; i < formulaX.length; i++) {
                    formulaX[i] = formulaXL.get(i);
                }
            }
            else {
                
                formulaX = new String[columnName.length];
                for (int i = 0; i < formulaX.length; i++) {
                    formulaX[i] = columnName[i];
                }

            }
        }
        else {
            
            formulaX = formulaU[1].split("\\+");

        }

        xSize = formulaX.length;

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@  for (int i = 0; i < ySize; i++) {
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "Formula Y: " + i + " : " + formulaY[i]);
        @RF_TRACE_OFF@  }
        @RF_TRACE_OFF@  }
        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@  for (int i = 0; i < xSize; i++) {
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "Formula X: " + i + " : " + formulaX[i]);
        @RF_TRACE_OFF@  }
        @RF_TRACE_OFF@  }

    }
    

    private void parseTypes() {

        StructField structField;
        Character tempType;
        
        
        StructType incomingSchema = dataset.schema();

        
        xFactorCount = yFactorCount = 0;
        
        if (formulaY != null) {
            yTypeHash = new HashMap <String, Character> (formulaY.length);
            yType = new char[formulaY.length];

            if (family.compareTo("RF-S") == 0) {
                structField = incomingSchema.apply(incomingSchema.fieldIndex(formulaY[0]));
                
                tempType = new Character(initializeType(structField.dataType()));
                if ((tempType.compareTo('R') == 0) || (tempType.compareTo('I') == 0)) {
                    
                    
                    
                    yTypeHash.put(structField.name(), tempType);
                    yType[0] = 'T';
                }
                else {
                    
                    RFLogger.log(Level.SEVERE, "Survival response " + formulaY[0] + " of incorrect type:  " + tempType);
                    throw new IllegalArgumentException();
                }
                
                structField = incomingSchema.apply(incomingSchema.fieldIndex(formulaY[1]));
                
                tempType = new Character(initializeType(structField.dataType()));
                if (tempType.compareTo('I') == 0) {
                    
                    
                    yTypeHash.put(structField.name(), tempType);
                    yType[1] = 'S';
                }
                else {
                    
                    RFLogger.log(Level.SEVERE, "Survival response " + formulaY[1] + " of incorrect type:  " + tempType);
                    throw new IllegalArgumentException();
                }
            }
            else {
                for (int i = 0; i < formulaY.length; i++) {
                    structField = incomingSchema.apply(incomingSchema.fieldIndex(formulaY[i]));          
                    
                    tempType = new Character(initializeType(structField.dataType()));
                    yTypeHash.put(structField.name(), tempType);
                    if (tempType.compareTo('C') == 0) {
                        yFactorCount++;
                        
                        yType[i] = tempType.charValue();
                    }
                    else {
                        
                        yType[i] = 'R';
                    }
                }
            }
        }
        else {
            yTypeHash = null;
        }

        xTypeHash = new HashMap <String, Character> (formulaX.length);
        xType = new char[formulaX.length];
        
        for (int i = 0; i < formulaX.length; i++) {
            structField = incomingSchema.apply(incomingSchema.fieldIndex(formulaX[i]));
            
            tempType = new Character(initializeType(structField.dataType()));
            xTypeHash.put(structField.name(), tempType);
            if (tempType.compareTo('C') == 0) {
                xFactorCount++;
                
                xType[i] = tempType.charValue();
            }
            else {
                
                xType[i] = 'R';
            }
        }

        
        tFactorCount = xFactorCount + yFactorCount;
        if (tFactorCount > 0) {
            
            immutableMap = new HashMap <String, HashMap> (tFactorCount);
        }
        
        @RF_TRACE_OFF@  if (Trace.get(Trace.HGH)) {
        @RF_TRACE_OFF@  Set set;
        @RF_TRACE_OFF@  Iterator itr;
        @RF_TRACE_OFF@  int i;
        @RF_TRACE_OFF@  if (formulaY != null) {
        @RF_TRACE_OFF@    set = yTypeHash.entrySet();
        @RF_TRACE_OFF@    itr = set.iterator();
        @RF_TRACE_OFF@    i = 0;
        @RF_TRACE_OFF@    while(itr.hasNext()) {
        @RF_TRACE_OFF@      Map.Entry me = (Map.Entry) itr.next();
        @RF_TRACE_OFF@      RFLogger.log(Level.INFO, "Schema Y: " + i + " : " + me.getValue() + " " + me.getKey() );
        @RF_TRACE_OFF@      i++;
        @RF_TRACE_OFF@    }
        @RF_TRACE_OFF@  }
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, " ");
        @RF_TRACE_OFF@  set = xTypeHash.entrySet();
        @RF_TRACE_OFF@  itr = set.iterator();
        @RF_TRACE_OFF@  i = 0;
        @RF_TRACE_OFF@  while(itr.hasNext()) {
        @RF_TRACE_OFF@    Map.Entry me = (Map.Entry) itr.next();
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "Schema X: " + i + " : " + me.getValue() + " " + me.getKey() );
        @RF_TRACE_OFF@    i++;
        @RF_TRACE_OFF@  }
        @RF_TRACE_OFF@  }
        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {        
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, " ");
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "X-factor count: " + xFactorCount);
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Y-factor count: " + yFactorCount);
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, " ");
        @RF_TRACE_OFF@  }
    }

    private static char initializeType(DataType dataType) {
        char result;
        if (dataType == DataTypes.BooleanType) {
            result = 'B';
        }
        else if (dataType == DataTypes.DoubleType) {
            result = 'R';
        }
        else if (dataType == DataTypes.IntegerType) {
            result = 'I';
        }
        else if (dataType == DataTypes.StringType) {
            result = 'C';
        }
        else {
            
            result = 'X';
            RFLogger.log(Level.SEVERE, "Bad DataType:  " + dataType);
            throw new IllegalArgumentException();
        }
        return result;
    }


    private void initializeFamily() {
        String[] form = formula.split("~");
        if (form[0].startsWith("Surv")) {
            
            
        }
        else if (form[0].startsWith("Unsupervised")) {
            
            
        }
        else if (form[0].startsWith("Multivariate")) {
            int factorCount = 0;
            Character response;
            for (int i = 0; i < formulaY.length; i++) {
                response = yTypeHash.get(formulaY[i]);
                if ((response.compareTo('C') == 0) || (response.compareTo('I') == 0)) {
                    factorCount ++;
                }
            }
            if (factorCount == 0) {
                family = new String("RF-R+");
            }
            else if (factorCount == formulaY.length) {
                family = new String("RF-C+");
            }
            else {
                family = new String("RF-M+");
            }

        }
        else {
            
            Character response = yTypeHash.get(formulaY[0]);
            if ((response.compareTo('C') == 0) || (response.compareTo('I') == 0)) {
                family = new String("RF-C");
            }
            else {
                family = new String("RF-R");
            }
        }

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Family: " + family);
        @RF_TRACE_OFF@  }
    }







    
    private double[][] initializeData(int zSize, String[] formulaZ, HashMap <String, Character> zTypeHash) {

        
        
        
        
        
        
        
        
        
        
        
        
        

        
        nSize = (int) dataset.count();

        
        double[][] zData = new double[zSize][];
        
        for (int i = 0; i < zSize; i++) {

            
            
            
            
            if (formulaZ[i].contains(".")) {
                RFLogger.log(Level.SEVERE, "Column name cannot contain dots:  " + formulaZ[i]);
                RFLogger.log(Level.SEVERE, "Please pre-process the data frame and replace dots with underlines.");
            }
            
            
            zData[i] = new double[nSize];

            
            org.apache.spark.sql.Row[] thisRow = (org.apache.spark.sql.Row[]) dataset.select(formulaZ[i]).collect();
            
            if ((zTypeHash.get(formulaZ[i])).compareTo('B') == 0) {

                
                
                for (int j = 0; j < nSize; j++) {
                    zData[i][j] = (thisRow[j].getBoolean(0)) ? 2 : 1;
                }

                @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
                @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Found Boolean:   " + i + " " + formulaZ[i]);                
                @RF_TRACE_OFF@  }
            }
            else if ((zTypeHash.get(formulaZ[i])).compareTo('R') == 0) {

                for (int j = 0; j < nSize; j++) {
                    zData[i][j] = (double) thisRow[j].getDouble(0);
                }

                @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
                @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Found Double:   " + i + " " + formulaZ[i]);                
                @RF_TRACE_OFF@  }
            }
            else if ((zTypeHash.get(formulaZ[i])).compareTo('I') == 0) {

                for (int j = 0; j < nSize; j++) {
                    zData[i][j] = (double) thisRow[j].getInt(0);
                }

                @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
                @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Found Integer:   " + i + " " + formulaZ[i]);
                @RF_TRACE_OFF@  }

                    
            }
            else if ((zTypeHash.get(formulaZ[i])).compareTo('C') == 0) {

                
                
                
                
                
                org.apache.spark.sql.Row[] thisMappedRow = mapFactor(formulaZ[i]);

                for (int j = 0; j < nSize; j++) {
                    zData[i][j] = (double) thisMappedRow[j].getDouble(0);
                }

                @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
                @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Found Factor: " + i + " " + formulaZ[i]);
                @RF_TRACE_OFF@  }
            }

            @RF_TRACE_OFF@  if (Trace.get(Trace.HGH)) {
            @RF_TRACE_OFF@  for (int j = 0; j < nSize; j++) {
            @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "zData[" + i + "]" + "[" + j + "] = " + zData[i][j]);
            @RF_TRACE_OFF@  }
            @RF_TRACE_OFF@  }

        }

        return zData;
        
    }


    private org.apache.spark.sql.Row[] mapFactor(String name) {

        
        
        StringIndexer indexer = new StringIndexer()
            .setInputCol(name)
            .setOutputCol("idxF");

        
        
        
        
        Dataset <Row> transformedF = indexer.fit(dataset.select(name)).transform(dataset.select(name));

        @RF_TRACE_OFF@  if (Trace.get(0)) {
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "Transformed factor:  ");
        @RF_TRACE_OFF@    transformedF.show();
        @RF_TRACE_OFF@  }
        
        
        
        
        Dataset <Row> compliantF = transformedF.select(transformedF.col("idxF").plus(1.0).alias("idxF"));

        @RF_TRACE_OFF@  if (Trace.get(0)) {
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "Native-code compliant factor (non-zero levels):  ");            
        @RF_TRACE_OFF@    compliantF.show();
        @RF_TRACE_OFF@  }
        
        
        
        
        org.apache.spark.sql.Row[] rowFactorOriginal    = (org.apache.spark.sql.Row[]) transformedF.select(name).collect();
        org.apache.spark.sql.Row[] rowFactorTransformed = (org.apache.spark.sql.Row[]) compliantF.select("idxF").collect();

        
        
        
        
        
        Dataset <Row> distinctF = transformedF.select(name).distinct();

        
        
        int levelCount = (int) distinctF.select(name).count();
        
        
        HashMap <String, Double> fMap = new HashMap <String, Double> (levelCount);

        int addedCount = 0;
        int tempIter   = 0;
        
        while (addedCount < levelCount) {
            
            if (!fMap.containsKey(rowFactorOriginal[tempIter].getString(0))) {
                fMap.put(rowFactorOriginal[tempIter].getString(0), rowFactorTransformed[tempIter].getDouble(0));
                addedCount++;
            }
            tempIter++;
        }

        
        @RF_TRACE_OFF@  if (Trace.get(Trace.MED)) {        
        @RF_TRACE_OFF@  Set set = fMap.entrySet();
        @RF_TRACE_OFF@  Iterator itr = set.iterator();
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Immutable map for factor:  " + name);
        @RF_TRACE_OFF@  while(itr.hasNext()) {
        @RF_TRACE_OFF@    Map.Entry me = (Map.Entry)itr.next();
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "<KEY, VALUE> = < " + me.getKey() + "," + me.getValue() + " >");
        @RF_TRACE_OFF@  }
        @RF_TRACE_OFF@  }

        
        
        
        
        if (!immutableMap.containsKey(name)) {
            immutableMap.put(name, fMap);
        }
        else {
            RFLogger.log(Level.SEVERE, "Duplicate factor found in dataset:  " + name);
            throw new IllegalArgumentException();
        }


        
        return (org.apache.spark.sql.Row[]) compliantF.select("idxF").collect();
        
    }


    
    
    
    
    
    
    
    
    
    
    
    private int[] initializeLevel(String[] formulaZ, HashMap <String, Character> zTypeHash) {
        
        Character variable;
        
        int[] zLevel = new int[formulaZ.length];

        for (int i = 0; i < formulaZ.length; i++) {
            variable = zTypeHash.get(formulaZ[i]);
            if (variable.compareTo('B') == 0) {
                
                zLevel[i] = 2;
            }
            else if (variable.compareTo('I') == 0) {

                if (false) {
                    
                    
                    
                    org.apache.spark.sql.Row[] rowInteger = (org.apache.spark.sql.Row[]) dataset.select(formulaZ[i]).collect();

                    zLevel[i] = 0;

                    for (int j = 0; j < rowInteger.length; j++) {
                        if (rowInteger[j].getInt(0) > zLevel[i]) {
                            zLevel[i] = rowInteger[j].getInt(0);
                        }
                    }
                }
                else {
                    zLevel[i] = 0;
                }
                
            
            }
            else if (variable.compareTo('C') == 0) {
                
                zLevel[i] = (immutableMap.get(formulaZ[i])).size();
            }
            else {
                
                zLevel[i] = 0;
            }
        }

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@    for (int i = 0; i < formulaZ.length; i++) {
        @RF_TRACE_OFF@      RFLogger.log(Level.INFO, "zLevel[" + i + "] = " + zLevel[i]);            
        @RF_TRACE_OFF@    }
        @RF_TRACE_OFF@    }
            
        return zLevel;
        
    }

    public String get_family() {
        return family;
    }

    public int get_nSize() {
        return nSize;
    }

    public int get_ySize() {
        return ySize;
    }

    public int get_xSize() {
        return xSize;
    }
    
    public char[] get_yType() {

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@    for (int i = 0; i < yType.length; i++) {
        @RF_TRACE_OFF@      RFLogger.log(Level.INFO, "yType[" + i + "] = " + yType[i]);            
        @RF_TRACE_OFF@    }
        @RF_TRACE_OFF@  }
        return yType;
    }
    
    public char[] get_xType() {

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@    for (int i = 0; i < xType.length; i++) {
        @RF_TRACE_OFF@      RFLogger.log(Level.INFO, "xType[" + i + "] = " + xType[i]);            
        @RF_TRACE_OFF@    }
        @RF_TRACE_OFF@  }

        return xType;
    }

    public int[] get_yLevel() {

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {        
        @RF_TRACE_OFF@    for (int i = 0; i < yType.length; i++) {
        @RF_TRACE_OFF@      RFLogger.log(Level.INFO, "yLevel[" + i + "] = " + yLevel[i]);            
        @RF_TRACE_OFF@    }
        @RF_TRACE_OFF@  }

        return yLevel;
    }
    
    public int[] get_xLevel() {

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {        
        @RF_TRACE_OFF@    for (int i = 0; i < xType.length; i++) {
        @RF_TRACE_OFF@      RFLogger.log(Level.INFO, "xLevel[" + i + "] = " + xLevel[i]);            
        @RF_TRACE_OFF@    }
        @RF_TRACE_OFF@  }

        return xLevel;
    }

    public double[][] get_yData() {
        return yData;
    }

    public double[][] get_xData() {
        return xData;
    }
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    public void set_bootstrap(int ntree, String bootstrap, String sampleType, int sampleSize, int[][] sample) {
        if (ntree <= 0) {
            RFLogger.log(Level.WARNING, "Invalid value for bootstrap parameter ntree:  " + ntree);
            RFLogger.log(Level.WARNING, "Overriding all bootstrap parameters with default values.");
            set_bootstrap();
        }
        else if (bootstrap.equals("auto")) {

            this.ntree = ntree;
            this.bootstrap = "auto";
            
            if (sampleType.equals("swr")) {
                this.sampleType = sampleType;
                
                if (sampleSize > 0) {
                    
                }
                else {
                    
                    
                    
                    this.sampleSize = nSize;
                }

            }
            else if (sampleType.equals("swor")) {
                this.sampleType = sampleType;
                if ((sampleSize > 0) && (sampleSize <= nSize)) {
                    
                    this.sampleSize = sampleSize;
                }
                else {
                    
                    
                    
                    this.sampleSize = (int) Math.round(nSize * (1 - Math.exp(-1)));
                }
            }
            else {
                RFLogger.log(Level.WARNING, "Invalid value for bootstrap parameter sampleType:  " + sampleType);                
                RFLogger.log(Level.WARNING, "Overriding all bootstrap parameters with default values.");
                set_bootstrap();
            }

            this.sample = null;
        }
        else if (bootstrap.equals("user")) {

            this.ntree = ntree;
            this.bootstrap = "user";
            this.sampleType = "swr";

            
            if (sample == null) {
                RFLogger.log(Level.SEVERE, "sample must not be null when bootstrapping by user.");
                throw new IllegalArgumentException();
            }
            else {
                int[] sSize = new int[ntree];

                if (sample.length != nSize) {
                    RFLogger.log(Level.SEVERE, "sample[.][] must be specified for each case:  " +  sample.length);
                    throw new IllegalArgumentException();
                }
                    
                for (int j = 0; j < ntree; j++) {
                    if (sample[j].length != nSize) {
                        RFLogger.log(Level.SEVERE, "sample[" + j + "][.] must be specified for each case:  " +  sample[j].length);
                        throw new IllegalArgumentException();
                    }
                }     
                
                for (int i = 0; i < ntree; i++) {
                    sSize[i] = 0;
                    for (int j = 0; j < nSize; j++) {
                        sSize[i] += sample[i][j];
                    }
                }
                for (int i = 0; i < ntree; i++) {
                    if (sSize[i] != sSize[0]) {
                        RFLogger.log(Level.SEVERE, "sample size must be identical for each tree:  " +  sSize[0]);
                        throw new IllegalArgumentException();
                    }
                }
                this.sampleSize = sSize[0];
                this.sample = sample;
            }

        }
        else {
            RFLogger.log(Level.WARNING, "Invalid value for bootstrap parameter bootstrap:  " + bootstrap);
            RFLogger.log(Level.WARNING, "Overriding all bootstrap parameters with default values.");
            set_bootstrap();
        }
    }

    public void set_bootstrap(int ntree) {
        this.ntree = 1000;
        bootstrap = "auto";
        sampleType = "swr";
        sampleSize = nSize;
        sample = null;
    }
    
    public void set_bootstrap() {
        ntree = 1000;
        bootstrap = "auto";
        sampleType = "swr";
        sampleSize = nSize;
        sample = null;
    }

    public String get_bootstrap() {
        return bootstrap;
    }

    public String get_sampleType() {
        return sampleType;
    }

    public int get_sampleSize() {
        return sampleSize;
    }

    public int[][] get_sample() {
        return sample;
    }


    public void set_mtry(int mtry) {
        if ((mtry < 1) || (mtry > xSize)) {
            set_mtry();
        }
        else {
            this.mtry = mtry;
        }
    }
    public void set_mtry() {
        if (family.equals("RF-R")) {
            mtry = (int) Math.ceil((double) xSize / 3.0); 
        }
        else {
            mtry = (int) Math.ceil(Math.sqrt((double) xSize));
        }
    }

    public int get_mtry() {
        return mtry;
    }


    public void set_nImpute() {
        nImpute = 1;
    }

    public void set_nImpute(int nImpute) {
        this.nImpute = nImpute;
        if (nImpute < 1) {
            RFLogger.log(Level.WARNING, "Invalid value for parameter nImpute:  " + nImpute);
            RFLogger.log(Level.WARNING, "Overriding nImpute with default value:  1");
            this.nImpute = 1;    
        }

    }

    public int get_nImpute() {
        return nImpute;
    }


    public int get_ntree() {
        return ntree;
    }


    public void set_caseWeight(double[] weight) {
        caseWeight = setWeight(weight, nSize);
    }
    public void set_caseWeight() {
        caseWeight = setWeight(null, nSize);
    }

    public double[] get_caseWeight() {
        return caseWeight;
    }
    

    public void set_ytry(int ytry) {
        if (family.equals("RF-U")) {
            RFLogger.log(Level.SEVERE, "ytry must be specified via the formula when family is RF-U.");
            throw new IllegalArgumentException();
        }
        this.ytry = ytry;

        if (ytry > 0) {
            
            set_yWeight();
        }
        else {
            
            yWeight = null;
        }
    }

    public int get_ytry() {
        return ytry;
    }

    public void set_yWeight(double[] weight) {
        if (family.equals("RF-S")) {
            RFLogger.log(Level.SEVERE, "yWeight cannot be specified when family is RF-S.");
            throw new IllegalArgumentException();
        }
        if (ytry == 0) {
            RFLogger.log(Level.SEVERE, "yWeight cannot be specified when ytry equals zero (0).");
            throw new IllegalArgumentException();
        }
        yWeight = setWeight(weight, ySize);
    }

    public void set_yWeight() {
        if (family.equals("RF-S")) {
            RFLogger.log(Level.SEVERE, "yWeight cannot be specified when family is RF-S.");
            throw new IllegalArgumentException();
        }
        if (ytry == 0) {
            RFLogger.log(Level.SEVERE, "yWeight cannot be specified when ytry equals zero (0).");
            throw new IllegalArgumentException();
        }
        yWeight = setWeight(null, ySize);
    }

    public double[] get_yWeight() {
        return yWeight;
    }

    public void set_xWeight(double[] weight) {
        xWeight = setWeight(weight, xSize);
    }
    public void set_xWeight() {
        xWeight = setWeight(null, xSize);
    }
    public double[] get_xWeight() {
        return xWeight;
    }

    public void set_xSplitStatWt(double[] weight) {
        xSplitStatWt = setWeight(weight, xSize);
    }
    public void set_xSplitStatWt() {
        xSplitStatWt = setWeight(null, xSize);
    }
    public double[] get_xSplitStatWt() {
        return xSplitStatWt;
    }

        

    
    private double[] setWeight(double[] weight, int size) {
        double[] generic;
        
        if (weight == null) {
            generic = new double[size];
            for (int i = 0; i < size; i++) {
                generic[i] = 1;
            }
        }
        else {
            generic = weight;
            if (generic.length != size) {
                RFLogger.log(Level.SEVERE, "weight vector must be of length " + size + " not " + generic.length);
                throw new IllegalArgumentException();
            }
            for (int i = 0; i < size; i++) {
                if (generic[i] < 0) {
                    RFLogger.log(Level.SEVERE, "weights must be greater than or equal to zero:  [" + i + "] = " + generic[i]);
                    throw new IllegalArgumentException();
                }
            }
        }
        return generic;
    }

    


    private TreeSet <Integer> initializeEventType() {
        
        TreeSet <Integer> treeSet = new TreeSet <Integer> (); 

        

        
        
        
        
        for (int i = 0; i < nSize; i++) {
            
            if ((yData[1][i] >= 0) && ((yData[1][i] % 1.0) == 0)) {
                
                if (yData[1][i] > 0) {
                    treeSet.add((int) yData[1][i]);
                }
            }
            else {
                RFLogger.log(Level.SEVERE, "censoring must be a whole number:  [" + i + "] = " + yData[1][i]);
                throw new NumberFormatException();
            }
        }

        
        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@  Iterator <Integer> itr = treeSet.iterator();
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Event Type:  ");
        @RF_TRACE_OFF@  int i = 0;
        @RF_TRACE_OFF@  while(itr.hasNext()) {
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "[" + i + "] = " + itr.next());
        @RF_TRACE_OFF@    i++;
        @RF_TRACE_OFF@  }
        @RF_TRACE_OFF@  }
        
        return treeSet;

    }

    public void set_eventWeight() {
        
        
        
        
        
        if (family.equals("RF-S")) {
            crWeight = setWeight(null, eventType.size());
        }
        else {
            RFLogger.log(Level.SEVERE, "family must be RF-S in order to define this weight vector.");
            throw new IllegalArgumentException();
        }

        crWeightSize = eventType.size();
    }

    public void set_eventWeight(double[] weight) {

        
        
        
        
        if (family.equals("RF-S")) {
            crWeight = setWeight(weight, eventType.size());
        }
        else {
            RFLogger.log(Level.SEVERE, "family must be RF-S in order to define this weight vector.");
            throw new IllegalArgumentException();
        }
        
        crWeightSize = eventType.size();
    }

    public double[] get_eventWeight() {
        return crWeight;
    }

    public int get_eventWeightSize() {
        return crWeightSize;
    }
        
    public void set_timeInterest() {
        
        TreeSet <Double> treeSet = new TreeSet <Double> (); 

        

        
        

        for (int i = 0; i < nSize; i++) {
            
            if (yData[1][i] > 0) {
                treeSet.add((double) yData[0][i]);
            }
        }
        
        timeInterest = new double[treeSet.size()];
        Iterator <Double> itr = treeSet.iterator();
        int j = 0;
        while (itr.hasNext()) {
            timeInterest[j++] = (double) itr.next();
        }

        timeInterestSize = timeInterest.length;
        
        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Setting Time Interest:  ");
        @RF_TRACE_OFF@  for(int i = 0; i < timeInterest.length; i++) {
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "[" + i + "] = " + timeInterest[i]);
        @RF_TRACE_OFF@  }
        @RF_TRACE_OFF@  }
        
    }


    public void set_timeInterest(double[] timeInterest) {

        Iterator <Double> itr;        
        double value;
        int j;
        
        TreeSet <Double> inTime = new TreeSet <Double> (); 

        

        
        

        
        for (int i = 0; i < timeInterest.length; i++) {
            inTime.add((double) timeInterest[i]);
        }

        
        set_timeInterest();
        double[] naturalTimeInterest = timeInterest;
        
        
        
        
        
        TreeSet <Double> outTime = new TreeSet <Double> (); 
        
        itr = inTime.iterator();

        j = 0;
        while (itr.hasNext()) {
            value = (double) itr.next();
            if (value <= naturalTimeInterest[j]) {
                outTime.add(naturalTimeInterest[j]);
            }
            else {
                while ((value > naturalTimeInterest[j]) && (j < inTime.size()))  {
                    j++;
                }
                if (j < inTime.size()) {
                    outTime.add(naturalTimeInterest[j]);                
                }
            }
        }

        this.timeInterest = new double[outTime.size()];
        itr = outTime.iterator();
        j = 0;
        while (itr.hasNext()) {
            this.timeInterest[j++] = (double) itr.next();
        }

        this.timeInterestSize = this.timeInterest.length;
        
        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Setting Time Interest:  ");
        @RF_TRACE_OFF@  for(int i = 0; i < this.timeInterestSize; i++) {
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "[" + i + "] = " + this.timeInterest[i]);
        @RF_TRACE_OFF@  }
        @RF_TRACE_OFF@  }
        
    }

    public double[] get_timeInterest() {
        return timeInterest;
    }

    public int get_timeInterestSize() {
        return timeInterestSize;
    }
    
    public void set_splitRule() {

        
        customSplitIndex = 0;
        
        if (family.equals("RF-S")) {
            if (eventType.size() == 1) {
                splitRule = "logrank";
            }
            else {
                splitRule = "logrankCR";
            }
        }
        else if (family.equals("RF-C")) {
            splitRule = "gini";
        }
        else if (family.equals("RF-R")) {
                splitRule = "mse";
        }
        else if (family.equals("RF-C+")) {
            splitRule = "mv.gini";
        }
        else if (family.equals("RF-R+")) {
            splitRule = "mv.mse";
        }
        else if (family.equals("RF-M+")) {
            splitRule = "mv.mix";
        }
        else if (family.equals("RF-U")) {
            splitRule = "unsupv";
        }
    }


    public void set_splitRule(String splitRule) {

        boolean result = false;

        
        customSplitIndex = 0;

        if (splitRule.equals("random")) {
            this.splitRule = splitRule;
            result = true;
        }
        else if (splitRule.startsWith("custom")) {
            String numberString = splitRule.replaceFirst("custom", "");
            
            if (numberString.equals("")) {
                this.splitRule = splitRule;
                customSplitIndex = 1;
                result = true;
            }
            else {
                
                customSplitIndex = Integer.parseInt(numberString);
                
                if ((customSplitIndex >> 5) == 0) {
                    this.splitRule = splitRule;
                    result = true;
                }
            }
        }
        else if (family.equals("RF-S")) {
            if (eventType.size() == 1) {
                if (splitRule.equals("logrank") || splitRule.equals("logrankscore")) {
                    this.splitRule = splitRule;
                    result = true;
                }
            }
            else {
                if (splitRule.equals("logrank")) {
                    RFLogger.log(Level.WARNING, "Overriding parameter splitrule with (logrankCR) when competing risks are present.");
                    this.splitRule = "logrankCR";
                    result = true;
                }
                else if (splitRule.equals("logrankCR")) {
                    this.splitRule = splitRule;
                    result = true;
                }
                else if (splitRule.equals("logrankACR")) {
                    this.splitRule = splitRule;
                    result = true;
                }
            }
        }
        else if (family.equals("RF-C")) {
            if (splitRule.equals("gini") || splitRule.equals("gini.unwt") || splitRule.equals("gini.hvwt")) {
                this.splitRule = splitRule;
                    result = true;
            }
        }
        else if (family.equals("RF-R")) {
            if (splitRule.equals("mse") || splitRule.equals("mse.unwt") || splitRule.equals("mse.hvwt")) {
                this.splitRule = splitRule;
                result = true;
            }
        }
        else if (family.equals("RF-C+")) {
            if (splitRule.equals("mv.gini")) {
                this.splitRule = splitRule;
                result = true;
            }
        }
        else if (family.equals("RF-R+")) {
            if (splitRule.equals("mv.mse")) {
                this.splitRule = splitRule;
                result = true;
            }
        }
        else if (family.equals("RF-M+")) {
            if (splitRule.equals("mv.mix")) {
                this.splitRule = splitRule;
                result = true;
            }
        }
        else if (family.equals("RF-U")) {
            if (splitRule.equals("unsupv")) {
                this.splitRule = splitRule;
                result = true;
            }
        }

        if (!result) {
            RFLogger.log(Level.SEVERE, "splitRule invalid:  " + splitRule);
            throw new IllegalArgumentException();
        }
    }

    public String get_splitRule() {
        return splitRule;
    }

    int getSplitRuleID(String splitRule) {
        return splitRuleID.get(splitRule);
    }
    
    public void set_nSplit() {
        if (splitRule.equals("random")) {
            nSplit = 1;
        }
        else {
            
            nSplit = 0;
        }
    }

    public void set_nSplit(int nSplit) {
        this.nSplit = nSplit;
        if (splitRule.equals("random")) {
            if (nSplit != 1) {
                RFLogger.log(Level.WARNING, "Overriding parameter nSplit with (1) when pure random splitting.");
            }
            this.nSplit = 1;
        }
    }

    public int get_nSplit() {
        return nSplit;
    }
    
    public void set_nodeSize() {
        if (family.equals("RF-S")) {
            if (eventType.size() == 1) {
                nodeSize = 3;
            }
            else {
                nodeSize = 6;
            }

        }
        else if (family.equals("RF-C") || family.equals("RF-C+")) {
            nodeSize = 1;
        }
        else if (family.equals("RF-R") || family.equals("RF-R+")) {
            nodeSize = 5;
        }
        else if (family.equals("RF-M+")) {
            nodeSize = 3;
        }
        else if (family.equals("RF-U")) {
            nodeSize = 3;
        }
    }
   
    public void set_nodeSize(int nodeSize) {
        if (nodeSize < 0) {
            RFLogger.log(Level.SEVERE, "nodeSize must be greater than zero:  " + nodeSize);
            throw new IllegalArgumentException();
        }
        this.nodeSize = nodeSize;
    }

    public int get_nodeSize() {
        return nodeSize;
    }
    
    public void set_nodeDepth() {
        nodeDepth = -1;
    }

    public void set_nodeDepth(int nodeDepth) {
        this.nodeDepth = nodeDepth;
    }

    public int get_nodeDepth() {
        return nodeDepth;
    }

    public void set_nullSplit() {
        nullSplit = "no";
    }

    public void set_nullSplit(String nullSplit) {
        if (nullSplit.equals("yes") || nullSplit.equals("no")) {
            this.nullSplit = nullSplit;
            if (nullSplit.equals("yes")) {
                if (family.equals("RF-U")) { 
                    RFLogger.log(Level.WARNING, "Overriding parameter nullSplit with (no) when family is RF-U.");            
                    this.nullSplit = "no";
                }
            }
        }
        else {
            RFLogger.log(Level.SEVERE, "Invalid value for nullSplit:  " + nullSplit);
            throw new IllegalArgumentException();
        }
    }

    public String get_nullSplit() {
        return nullSplit;
    }

    public void set_rfCores() {
        rfCores = -1;
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

    public int get_rfCores() {
        return rfCores;
    }
    
    public void set_seed() {
        generator = new Random();
        seed = - generator.nextInt();
    }

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

    public int get_seed() {
        return seed;
    }
    
    public void set_trace() {
        trace = 0;
    }

    public void set_trace(int trace) {
        this.trace = trace;
        
        if (trace < 0) {
            RFLogger.log(Level.WARNING, "Invalid value for parameter trace:  " + trace);
            RFLogger.log(Level.WARNING, "Overriding trace with default value:  0");
            set_trace();
        }
    }
        
    public int get_trace() {
        return trace;
    }

   

    private void setDefaultModelArg() {

        
        set_trace();


        
        
        
        
        parseFormula();

        
        parseTypes();

        
        initializeFamily();

        
        

        
        xData = initializeData(xSize, formulaX, xTypeHash);
        
        xLevel = initializeLevel(formulaX, xTypeHash);

        if (formulaY != null) {
            
            yData = initializeData(ySize, formulaY, yTypeHash);

            
            yLevel = initializeLevel(formulaY, yTypeHash);
        }
        else {
            yData = null;
            yLevel = null;
        }

        
        set_nImpute();

        
        set_mtry();

        
        set_bootstrap();

        
        set_caseWeight();
        set_xWeight();

        
        
        
        
        
        
        
        yWeight = null;

        set_xSplitStatWt();

            
        
        if(family.equals("RF-S")) {
            eventType = initializeEventType();

            
            set_timeInterest();

            set_eventWeight();

        }
        else {
            eventType = null;

            timeInterest = null;
            timeInterestSize = 0;

            crWeight = null;
            crWeightSize = 0;
        }

        
        set_splitRule();

        
        set_nSplit();


        
        set_nodeSize();

        
        set_nodeDepth();

        
        set_nullSplit();

        
        set_rfCores();

        
        set_seed();


    }






    public void set_ensembleArg() {
        ensembleArg.set();
    }

    public void set_ensembleArg(String key, String value) {
        ensembleArg.set(key, value);
    }

    public String get_ensembleArg(String key) {
        return ensembleArg.get(key);
    }
    
    
    int getEnsembleArgOptLow() {

        return (ensembleArg.getNative("varUsed") + 
                ensembleArg.getNative("splitDepth") +                 
                ensembleArg.getNative("importance") + 
                ensembleArg.getNative("proximity") +
                ensembleArg.getNative("forest") +
                ensembleArg.getNative("errorType")); 
    }

    int getEnsembleArgOptHigh() {

        return (ensembleArg.getNative("membership") + 

                ensembleArg.getNative("weight") +                 
                
                
                ensembleArg.getNative("error") +
                ensembleArg.getNative("qualitativeTerminalInfo")); 
    }

    int getModelArgOptLow() {

        NativeOpt nativeOpt;
        int result = 0;
        
        nativeOpt = new NativeOpt("bootstrap",
                                  new String[] {"auto", "user"},
                                  new int[] {0, (1 << 19) + (1 << 20)});

        result += nativeOpt.get(bootstrap);

        nativeOpt = new NativeOpt("nullSplit",
                                  new String[] {"no", "yes"},
                                  new int[] {0, (1 << 18)});

        result += nativeOpt.get(nullSplit);

        return result;
    }

    int getModelArgOptHigh() {

        NativeOpt nativeOpt;
        int customSplitBit = 0;

        int result = 0;

        nativeOpt = new NativeOpt("sampleType",
                                  new String[] {"swr", "swor"},
                                  new int[] {0, (1 << 12)});

        result += nativeOpt.get(sampleType);
        
        if (customSplitIndex > 0) {
            
            customSplitBit = (customSplitIndex - 1) << 8;
        }

        result += customSplitBit;
        
        return result;
    }


}

