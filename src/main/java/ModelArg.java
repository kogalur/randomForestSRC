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


/**
 * Class containing the user-controlled parameters that produce the {@link
 * RandomForest} object.  Parameters are of two types: those that
 * define how the forest is to be trained; those that define the
 * requested ensemble outputs.
 * @author Udaya Kogalur
 * 
 */
public class ModelArg {

    // *****************************************************************
    // Ensemble Requests.
    // *****************************************************************
    private EnsembleArg ensembleArg;
    
    // *****************************************************************
    // Auxiliary Variables.
    // *****************************************************************
    private Random generator;

    // *****************************************************************
    // Model Inputs passed directly to the native code.
    // *****************************************************************
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

    private double[]   xStatisticalWeight;
    private double[]   yWeight;

    private double[]   xWeight;
    private double[][] xData;

    private int        timeInterestSize;
    private double[]   timeInterest;
    private double[]   ntime;
    
    private int        rfCores;

    
    // *****************************************************************
    // Model Inputs passed via option bits to the native code.
    // *****************************************************************
    private String sampleType;
    private String nullSplit;
    private int    customSplitIndex;

    // *****************************************************************
    // Model Inputs pre-processed and passed to the native code.
    // *****************************************************************
    private String splitRule;
    private String bootstrap;


    // *****************************************************************
    // Model Inputs not passed, but help inform the native code.
    // *****************************************************************
    // User defined formula.
    private String formula;
    // Vectorized parsed formula.
    private String[] formulaU;
    // Vectorized x-var names in formula.
    private String[] formulaX;
    // Vectorized y-var names in formula.
    private String[] formulaY;

    // Family
    private String family;

    // Split rule names;
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
    
    // HashMap containing types for x-vars and y-vars.
    private java.util.HashMap <String, Character> yTypeHash;
    private java.util.HashMap <String, Character> xTypeHash;

    // Counts of factors for x-vars and y-vars.
    private int xFactorCount;
    private int yFactorCount;
    private int tFactorCount;
    
    // HashMap of HashMap(s) containing immutable mapping of classes
    // to double integers for each factor.
    HashMap <String, HashMap> immutableMap;

    // Incoming Spark Dataset.
    private Dataset dataset;

    // Incoming column names.
    private String[] columnName;

    // Vector of unique event types in survival families.
    private TreeSet <Integer> eventType;


    
    /**
     * Sets default values for the training parameters and ensemble
     * outputs for the forest, given the formula and dataset.  The
     * default values are specific to the family (RF-S, RF-R, RF-C,
     * RF-R+, RF-C+, RF-M+), and can be customized using the 
     * methods available by this class.
     *
     * <p>Examples of formulae for various families follow:</p>
     *
     * <p><table class= "myColumnPadding">
     *  <tr>
     *    <th>Description</th>
     *    <th>Example Formula</th>
     *    <th>Data Set</th>
     *  </tr>
     *  <tr>
     *    <td>survivial or competing risk</td>
     *    <td>Surv(time, status) ~ .</td>
     *    <td><a href="https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/master/csv/MASS/VA.csv" target="_blank">Veteran's Administration Lung Cancer Trial</a>
     *  </tr>
     *  <tr>
     *    <td>regression</td>
     *    <td>Ozone ~.</td>
     *    <td><a href="https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/master/csv/datasets/airquality.csv" target="_blank">New York Air Quality Measurements</a>
     *  </tr>
     *  <tr>
     *    <td>classification</td>
     *    <td>Species ~.</td>
     *    <td><a href="https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/master/csv/datasets/iris.csv" target="_blank">Edgar Anderson's Iris Data</a>
     *  </tr>
     *  <tr>
     *    <td>multivariate regression</td>
     *    <td>Multivar(mpg, cyl) ~ .</td>
     *    <td><a href="https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/master/csv/datasets/mtcars.csv" target="_blank">Motor Trend Car Road Tests</a>
     *  </tr>
     *  <tr>
     *    <td>multivariate regression</td>
     *    <td>Multivariate(mpg, wt) ~ hp + drat</td>
     *    <td><a href="https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/master/csv/datasets/mtcars.csv" target="_blank">Motor Trend Car Road Tests</a>
     *  </tr>
     *  <tr>
     *    <td>unsupervised</td>
     *    <td>Unsupervised() ~.</td>
     *    <td><a href="https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/master/csv/MASS/VA.csv" target="_blank">Veteran's Administration Lung Cancer Trial</a>
     *  </tr>
     * </table></p>
     *        
     * <p>An example found in the test classes follows:</p>
     *
     * <pre><code>
     * Dataset<Row> irisDF = spark
     *     .read()
     *     .option("header", "true")
     *     .option("inferSchema", "true") 
     *     .format("csv")
     *     .load("./test-classes/data/iris.csv");
     *
     * ModelArg modelArg = new ModelArg("Species ~ .", irisDF);
     * </code></pre>
     *
     * <p> A overview of all the data sets used above can be found <a href="https://vincentarelbundock.github.io/Rdatasets/datasets.html">here</a>.</p>
     * 
     * @param formula Specification of the y-variables and x-variables
     *        that are to be used in the model. These refer to the
     *        column names in the Spark Dataframe.  The y-variables
     *        and x-variables are separated by the tilde (~)
     *        character.  A period (.) indicates that the complement
     *        of the y-variables is to be used as the x-variables. See
     *        the example above for more information.
     * @param dataset A Spark Dataset. 
     */
    public ModelArg (String formula,
                      Dataset dataset) {


        this.dataset = dataset;
        this.formula = formula;

        // Extract the column names from the dataset.
        columnName = dataset.columns();

        setDefaultModelArg();

        ensembleArg = new EnsembleArg("grow", family);
    }

    

    private void parseFormula() {

        boolean bigRFlag, complementFlag;
        
        formula = formula.replaceAll("\\ ", "");
      
        formulaU = formula.split("~");

        // ytry can be set via the RF-U formula.  Thus we set the
        // default and user specified values here, directly.
        ytry = 0;
        
        if (formulaU.length != 2) {
            // Bad Formula.  Throw exception.
            System.out.println("\n Unknown Formula Syntax.");
        }

        // Initialize the family to an invalid value, in case it is
        // not pre-emptively initialized.
        family = new String("RF-X");
        
        if (formulaU[0].startsWith("Surv")) {
            // Check for proper syntax.
            if ((formulaU[0].indexOf('(') >= 0) && (formulaU[0].indexOf(')') >= 0)) {
                formulaU[0] = formulaU[0].replace("Surv", "");
                formulaU[0] = formulaU[0].replace("(", "");
                formulaU[0] = formulaU[0].replace(")", "");
                formulaY = formulaU[0].split(",");

                if (formulaY.length != 2) {;
                    // Bad Formula.  Throw exception.
                    RFLogger.log(Level.SEVERE, "Bad Survival Formula Syntax:  " + formula);
                    throw new IllegalArgumentException();
                }
            }
            else {
                // Bad Formula.  Throw exception.
                formulaY = null;
                RFLogger.log(Level.SEVERE, "Bad Survival Formula Syntax:  " + formula);
                throw new IllegalArgumentException();
            }
            
            // Initialize the family pre-emptively.
            family = new String("RF-S");
        }
        else if (formulaU[0].startsWith("Unsupervised")) {
            // Check for proper syntax.
            if ((formulaU[0].indexOf('(') >= 0) && (formulaU[0].indexOf(')') >= 0)) {
                formulaU[0] = formulaU[0].replace("Unsupervised", "");
                formulaU[0] = formulaU[0].replace("(", "");
                formulaU[0] = formulaU[0].replace(")", "");

                // Note that ytry is only implemented for RF-U.  It
                // will be expanded to include all non RF-S families,
                // in order to facilitate big-r analysis.
                // (TBD TBD)
                if (formulaU[0].equals("")) {
                    // This is unsupervised with the default value of ytry.
                    ytry = 1;
                }
                else {
                    // This is unsupervised with a user specified value of ytry.
                    ytry = Integer.parseInt(formulaU[0]);
                }

                if (formulaY != null) {
                    // Bad Formula.  Throw exception.
                    formulaY = null;
                    RFLogger.log(Level.SEVERE, "Bad Unsupervised Formula Syntax:  " + formula);
                    throw new IllegalArgumentException();
                }
            }
            else {
                // Bad Formula.  Throw exception.
                formulaY = null;
                RFLogger.log(Level.SEVERE, "Bad Unsupervised Formula Syntax:  " + formula);
                throw new IllegalArgumentException();
            }
            // Initialize the family pre-emptively.
            family = new String("RF-U");            
        }
        else if (formulaU[0].startsWith("Multivariate")) {
            // Check for proper syntax.
            if ((formulaU[0].indexOf('(') >= 0) && (formulaU[0].indexOf(')') >= 0)) {
                formulaU[0] = formulaU[0].replace("Multivariate", "");
                formulaU[0] = formulaU[0].replace("(", "");
                formulaU[0] = formulaU[0].replace(")", "");
                formulaY = formulaU[0].split(",");
            }
            else {
                // Bad Formula.  Throw exception.
                formulaY = null;
                RFLogger.log(Level.SEVERE, "Bad Multivariate Formula Syntax:  " + formula );
                throw new IllegalArgumentException();
            }
        }
        else {
            // Univariate response.
            formulaY = new String[1];
            formulaY[0] = formulaU[0];
        }
        


        // Detect whether we are in a big-r situation.  In such cases,
        // we will acquire the predictor data frame via unpacking rather
        // than dropping each response from the data frame in a loop over the responses.
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
            // x-vars are complement of y-vars.

            complementFlag = true;
          
            if (formulaY != null) {
                // The length of the x-var names is known.
                formulaX = new String[columnName.length - formulaY.length];

                // Convert the data frame names to a list for easy removal of the response names.
                List <String> formulaXL = new ArrayList <String> (Arrays.asList(columnName));
                for (int j = 0; j < formulaY.length; j++) {
                    formulaXL.remove(formulaY[j]);
                }
                for (int i = 0; i < formulaX.length; i++) {
                    formulaX[i] = formulaXL.get(i);
                }
            }
            else {
                // We are in an unsupervised scenario.  The x-var names will the be data frame names.
                formulaX = new String[columnName.length];
                for (int i = 0; i < formulaX.length; i++) {
                    formulaX[i] = columnName[i];
                }

            }
        }
        else {
            // x-vars are explicity specified. Parse the formula and extract them.
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
        
        // Access the schema of the incoming data set.
        StructType incomingSchema = dataset.schema();

        // Set the factor counts to zero.
        xFactorCount = yFactorCount = 0;
        
        if (formulaY != null) {
            yTypeHash = new HashMap <String, Character> (formulaY.length);
            yType = new char[formulaY.length];

            if (family.compareTo("RF-S") == 0) {
                structField = incomingSchema.apply(incomingSchema.fieldIndex(formulaY[0]));
                // Note tha structField.name() == formulaY[i].
                tempType = new Character(initializeType(structField.dataType()));
                if ((tempType.compareTo('R') == 0) || (tempType.compareTo('I') == 0)) {
                    // The time will be overridden to type 'T' but
                    // remain as 'R' in the HashMap to allow proper
                    // extraction.
                    yTypeHash.put(structField.name(), tempType);
                    yType[0] = 'T';
                }
                else {
                    // Bad data types in survival formula.  Throw exception.
                    RFLogger.log(Level.SEVERE, "Survival response " + formulaY[0] + " of incorrect type:  " + tempType);
                    throw new IllegalArgumentException();
                }
                
                structField = incomingSchema.apply(incomingSchema.fieldIndex(formulaY[1]));
                // Note tha structField.name() == formulaY[i].
                tempType = new Character(initializeType(structField.dataType()));
                if (tempType.compareTo('I') == 0) {
                    // The status will be overriden to type 'S' but remain as 'I' in the HashMap to allow
                    // proper extraction.
                    yTypeHash.put(structField.name(), tempType);
                    yType[1] = 'S';
                }
                else {
                    // Bad data types in survival formula.  Throw exception.
                    RFLogger.log(Level.SEVERE, "Survival response " + formulaY[1] + " of incorrect type:  " + tempType);
                    throw new IllegalArgumentException();
                }
            }
            else {
                for (int i = 0; i < formulaY.length; i++) {
                    structField = incomingSchema.apply(incomingSchema.fieldIndex(formulaY[i]));          
                    // Note tha structField.name() == formulaY[i].
                    tempType = new Character(initializeType(structField.dataType()));
                    yTypeHash.put(structField.name(), tempType);
                    if (tempType.compareTo('C') == 0) {
                        yFactorCount++;
                        // Respect the variable type.
                        yType[i] = tempType.charValue();
                    }
                    else {
                        // Set the variable type to 'R' for now for both 'I' and 'R'.
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
            // Note tha structField.name() == formulaX[i].
            tempType = new Character(initializeType(structField.dataType()));
            xTypeHash.put(structField.name(), tempType);
            if (tempType.compareTo('C') == 0) {
                xFactorCount++;
                // Respect the variable type.
                xType[i] = tempType.charValue();
            }
            else {
                // Set the variable type to 'R' for now for both 'I' and 'R'.
                xType[i] = 'R';
            }
        }

        // Initialize the total factor count.
        tFactorCount = xFactorCount + yFactorCount;
        if (tFactorCount > 0) {
            // Allocate the HashMap of HashMap(s) containing the immutable mapping of factor class to double integer values.
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
            // Error handler.
            result = 'X';
            RFLogger.log(Level.SEVERE, "Bad DataType:  " + dataType);
            throw new IllegalArgumentException();
        }
        return result;
    }


    private void initializeFamily() {
        String[] form = formula.split("~");
        if (form[0].startsWith("Surv")) {
            // Family has been pre-emptively initialized in parseFormula().
            // family = new String("RF-S");
        }
        else if (form[0].startsWith("Unsupervised")) {
            // Family has been pre-emptively initialized in parseFormula().
            // family = new String("RF-U");
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
            // Univariate response.
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

        // Factors must be treated differently in forming the X- or
        // Y-matrix.  We know the number of total factors in the
        // dataset.  For convenience, formulaX and formulaY determines
        // the order of the X- and Y-matrix.  When a non-factor is
        // encountered, we simply extract the corresponding column of
        // the dataset as a <Row> object and copy the <Row> to the
        // corresponding row in the X-matrix.  When a factor is
        // encountered, we extract the corresponding column of the
        // dataset as a <Row> object, transform the string values to
        // double integer values, massage them to non-zero values by
        // adding one (1) and then copy the <Row> to the corresponding
        // row in the X-matrix.  When a boolean value is encountered,
        // it is treated as an integer type with true=1, false=0

        // Note sloppy conversion of (long) to (int).  This is not big-n compliant.
        nSize = (int) dataset.count();

        // Allocate the row pointer in the Z-matrix.
        double[][] zData = new double[zSize][];
        
        for (int i = 0; i < zSize; i++) {

            // Check that the column names don't have dots in them.
            // We currently have no elegant work-around for this, as
            // the SQL select statement throw a
            // org.apache.spark.sql.AnalysisException.
            if (formulaZ[i].contains(".")) {
                RFLogger.log(Level.SEVERE, "Column name cannot contain dots:  " + formulaZ[i]);
                RFLogger.log(Level.SEVERE, "Please pre-process the data frame and replace dots with underlines.");
            }
            
            // Allocate a row in the Z-matrix.
            zData[i] = new double[nSize];

            // Create a Row object that accesses the incoming non-factor in the dataset.
            org.apache.spark.sql.Row[] thisRow = (org.apache.spark.sql.Row[]) dataset.select(formulaZ[i]).collect();
            
            if ((zTypeHash.get(formulaZ[i])).compareTo('B') == 0) {

                // Boolean types are mapped as type 'I' to the
                // native-code.  Values are true=1, false=0.
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

                
                // The function call below creates the immutable map
                // of class to integer values.  It adds the map to the
                // HashMap of HashMap(s) with <Key, Value> given by
                // <String, HashMap>.  It initializes the appropriate
                // slot in xLevel[] that indicates the number of
                // levels for the factor being analyzed..Create a Row
                // object that accesses the incoming non-factor in the
                // dataset.  Finally, it returns the row object
                // containing the massaged non-zero integer values
                // that are native-code compliant.
                
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

        // Define the transformation that will change the string
        // representation of the factor level to a double integer.
        StringIndexer indexer = new StringIndexer()
            .setInputCol(name)
            .setOutputCol("idxF");

        // Transform the factor to a double integer equivalent.  There
        // will only be two columns in the resulting dataset: the
        // first being the incoming factor, the second being the
        // transformed factor.
        Dataset <Row> transformedF = indexer.fit(dataset.select(name)).transform(dataset.select(name));

        @RF_TRACE_OFF@  if (Trace.get(0)) {
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "Transformed factor:  ");
        @RF_TRACE_OFF@    transformedF.show();
        @RF_TRACE_OFF@  }
        
        // Massage the transformed factor by adding one (1) to all
        // elements.  Native-code factor levels must be greater than
        // zero (0).  The result is a dataset with only one column.
        Dataset <Row> compliantF = transformedF.select(transformedF.col("idxF").plus(1.0).alias("idxF"));

        @RF_TRACE_OFF@  if (Trace.get(0)) {
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "Native-code compliant factor (non-zero levels):  ");            
        @RF_TRACE_OFF@    compliantF.show();
        @RF_TRACE_OFF@  }
        
        // We create two Row objects. The first is the string
        // representation of the factor. The second is the double
        // integer representation of the factor.
        org.apache.spark.sql.Row[] rowFactorOriginal    = (org.apache.spark.sql.Row[]) transformedF.select(name).collect();
        org.apache.spark.sql.Row[] rowFactorTransformed = (org.apache.spark.sql.Row[]) compliantF.select("idxF").collect();

        // Acquire the distinct levels in the factor.  This will be a
        // dataset with one column containing the string
        // representation of the levels.  It will be used as the key
        // in <key, value> pairs that define the immutable map for
        // this factor.
        Dataset <Row> distinctF = transformedF.select(name).distinct();

        // Save the number of distinct levels.  Note the sloppy cast
        // from (long) to (int).  This is not big-n compliant.
        int levelCount = (int) distinctF.select(name).count();
        
        // Define the immutable HashMap for the factor.
        HashMap <String, Double> fMap = new HashMap <String, Double> (levelCount);

        int addedCount = 0;
        int tempIter   = 0;
        
        while (addedCount < levelCount) {
            // Check if the key has been added to the HashMap?
            if (!fMap.containsKey(rowFactorOriginal[tempIter].getString(0))) {
                fMap.put(rowFactorOriginal[tempIter].getString(0), rowFactorTransformed[tempIter].getDouble(0));
                addedCount++;
            }
            tempIter++;
        }

        // Display immutable map.
        @RF_TRACE_OFF@  if (Trace.get(Trace.MED)) {        
        @RF_TRACE_OFF@  Set set = fMap.entrySet();
        @RF_TRACE_OFF@  Iterator itr = set.iterator();
        @RF_TRACE_OFF@  RFLogger.log(Level.INFO, "Immutable map for factor:  " + name);
        @RF_TRACE_OFF@  while(itr.hasNext()) {
        @RF_TRACE_OFF@    Map.Entry me = (Map.Entry)itr.next();
        @RF_TRACE_OFF@    RFLogger.log(Level.INFO, "<KEY, VALUE> = < " + me.getKey() + "," + me.getValue() + " >");
        @RF_TRACE_OFF@  }
        @RF_TRACE_OFF@  }

        // Add the immutable map to the HashMap of HashMap(s)
        // containing the immutable maps.  Check if the key has been
        // added to the HashMap of HashMap(s).  This should never be
        // the case that it is a duplicate.
        if (!immutableMap.containsKey(name)) {
            immutableMap.put(name, fMap);
        }
        else {
            RFLogger.log(Level.SEVERE, "Duplicate factor found in dataset:  " + name);
            throw new IllegalArgumentException();
        }


        // Return the native-code compliant Row.
        return (org.apache.spark.sql.Row[]) compliantF.select("idxF").collect();
        
    }


    // This function is called by x-var processing, and again by y-var
    // processing.  Set the level counts for factors.  Currently,
    // ordered factors are NOT SUPPORTED.  They are treated as double
    // values.  In the future, responses that are ordered factors will
    // be treated as a form of classification with a special split
    // rule.  They will not have an allocated immutable map in HashMap
    // immutableMap.  This is because the mapping will simply be the
    // identity.  However, we will have to to inform the native code
    // about the number of levels in the ordered factor response.
    // When are x-vars are ordered factors, they will be treated in
    // the same manner as a double for splitting, and imputation.
    private int[] initializeLevel(String[] formulaZ, HashMap <String, Character> zTypeHash) {
        
        Character variable;
        
        int[] zLevel = new int[formulaZ.length];

        for (int i = 0; i < formulaZ.length; i++) {
            variable = zTypeHash.get(formulaZ[i]);
            if (variable.compareTo('B') == 0) {
                // Default level for boolean types is two. They are considered unordered factors.
                zLevel[i] = 2;
            }
            else if (variable.compareTo('I') == 0) {

                if (false) {
                    // Acquire the maximum level in the ordered factor.
                    // First we target the ordered factor, and extract it
                    // as an array.
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
                // Get the number of levels from the immutable map.
                zLevel[i] = (immutableMap.get(formulaZ[i])).size();
            }
            else {
                // Safe the level slot, as it is not mapped to a factor.
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


    // Generic function to initialize a weight vector.  Null input imlies uniform sampling.
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

        // Note that yData is of the form yData[ySize][nSize];

        // Time is always in the [0][...] row.
        // Cens is always in the [1][...] row.  
        // This is a consequence of the Surv() formula protocol.
        
        for (int i = 0; i < nSize; i++) {
            // First, we check that Cens are all whole numbers.
            if ((yData[1][i] >= 0) && ((yData[1][i] % 1.0) == 0)) {
                // Second, we add events only.
                if (yData[1][i] > 0) {
                    treeSet.add((int) yData[1][i]);
                }
            }
            else {
                RFLogger.log(Level.SEVERE, "censoring must be a whole number:  [" + i + "] = " + yData[1][i]);
                throw new NumberFormatException();
            }
        }

        // Display the event type.
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

    private void setDefaultModelArg() {

        // Set trace to the defalut value. Force it to non-zero values
        // here if debugging is required.
        set_trace(0);


        // Define the vector of x-vars and y-vars based on the
        // formula.  This defines the order of the native-code
        // incoming variables. In addition, the forumla is
        // pre-emptively initialized for the RF-S and RF-U.
        parseFormula();

        // Define the HashMap(s) containing the variable types of the x-vars and y-vars.
        parseTypes();

        // Define the family, in the rest of the scenarios (not RF-S, not RF-U).
        initializeFamily();

        // When factors are encountered, they are transformed and
        // massaged into native-code compliant integer values.

        // Get the x-matrix.  
        xData = initializeData(xSize, formulaX, xTypeHash);
        // Determine the maximum levels encountered in the x-var factor variables.
        xLevel = initializeLevel(formulaX, xTypeHash);

        if (formulaY != null) {
            // Get the y-matrix.
            yData = initializeData(ySize, formulaY, yTypeHash);

            // Determine the maximum levels encountered in the y-var factor variables.
            yLevel = initializeLevel(formulaY, yTypeHash);
        }
        else {
            yData = null;
            yLevel = null;
        }

        // Default value of nImpute.
        set_nImpute(1);

        // Default value of mtry.
        set_mtry();

        // Default bootstrap.
        set_bootstrap();

        // Set case weight, x-weight, y-weight, and split weight.
        set_caseWeight();
        set_xWeight();

        // yWeight must always be specified explicitly, otherwise it
        // is null.  In contrast, xWeight is automatically initialized
        // as a uniform vector, with an associated default value of
        // mtry.  Here, yWeight is initialized to null, and ytry set
        // to zero (0).  In the case of RF-U, ytry is set via the
        // formula, and yWeight is currently ignored on the
        // native-side.  (TBD TBD)
        yWeight = null;

        set_xStatisticalWeight();

            
        // Set the event type vector and the time interest vector.
        if(family.equals("RF-S")) {
            eventType = initializeEventType();

            // Private method.
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

        // Set the default split rule.  This must occur AFTER eventType is defined.
        set_splitRule();

        // Set nSplit.  This must occur AFTER splitRule is defined.
        set_nSplit();


        // Set nodeSize.
        set_nodeSize();

        // Default value of nodeDepth.
        set_nodeDepth(-1);

        // Default value of nullSplit.  This must occur AFTER family is defined.
        set_nullSplit("no");

        // Default value of rfCores.
        set_rfCores(-1);

        // Set seed.
        set_seed();


    }
    
    /**
     * Returns the family of analysis intentioned by the {@link ModelArg} instance.  It is of the following form:
     * <pre> <code> <table class= "myColumnPadding">
     *  <tr>
     *    <th>Family</th>
     *    <th>Description</th>
     *  </tr>
     *  <tr>
     *    <td>RF-S</td>
     *    <td>survivial or competing risk</td>
     *  </tr>
     *  <tr>
     *    <td>RF-R</td>
     *    <td>regression</td>
     *  </tr>
     *  <tr>
     *    <td>RF-C</td>
     *    <td>classification</td>
     *  </tr>
     *  <tr>
     *    <td>RF-R+</td>
     *    <td>multivariate regression</td>
     *  </tr>
     *  <tr>
     *    <td>RF-C+</td>
     *    <td>multivariate classification</td>
     *  </tr>
     *  <tr>
     *    <td>RF-M+</td>
     *    <td>multivariate mixed</td>
     *  </tr>
     *  <tr>
     *    <td>RF-U</td>
     *    <td>unsupervised</td>
     *  </tr>
     * </table>
     * </code> </pre>
     * @return The family of analysis.
     */
    public String get_family() {
        return family;
    }

    /** 
     * Returns the number of records or rows (n) in the data set.
     * @return The number of records or rows (n) in the data set.
     */
    public int get_nSize() {
        return nSize;
    }

    /** 
     * Returns the number of y-variables in the data set. This will
     * be one (1) if the analysis is univariate, two (2) if the
     * analysis is survival related, zero (0) if the analysis is
     * unsupervised, and greater than one if the analysis is multivariate.
     * @return The number of y-variables in the data set.  
     */
    public int get_ySize() {
        return ySize;
    }

    /** 
     * Returns the number of x-variables in the data set. This will
     * always be greater than zero.
     * @return The number of x-variables in the data set. 
     */
    public int get_xSize() {
        return xSize;
    }

    /**
     * Returns a vector of length ySize ({@link #get_ySize}) containing the y-variables types.  The vector will be null in the family is unsupervised.
     * The types are as follows:
     * <pre> <code> 
     * <table class= "myColumnPadding">
     *  <tr>
     *    <th>Description</th>
     *    <th>Value</th>
     *  </tr>
     *  <tr>
     *    <td>time</td>
     *    <td>T</td>
     *  </tr>
     *  <tr>
     *    <td>censoring</td>
     *    <td>S</td>
     *  </tr>
     *  <tr>
     *    <td>boolean</td>
     *    <td>B</td>
     *  </tr>
     *  <tr>
     *    <td>real</td>
     *    <td>R</td>
     *  </tr>
     *  <tr>
     *    <td>ordinal</td>
     *    <td>O</td>
     *  </tr>
     *  <tr>
     *    <td>categorical</td>
     *    <td>C</td>
     *  </tr>
     *  <tr>
     *    <td>integer</td>
     *    <td>I</td>
     *  </tr>
     * </table>
     * </code> </pre>
     *
     * @return The y-variable types.
     */
    public char[] get_yType() {

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@    for (int i = 0; i < yType.length; i++) {
        @RF_TRACE_OFF@      RFLogger.log(Level.INFO, "yType[" + i + "] = " + yType[i]);            
        @RF_TRACE_OFF@    }
        @RF_TRACE_OFF@  }
        return yType;
    }
    
    /**
     * Returns a vector of length xSize ({@link #get_xSize}) containing the x-variable types.
     * The types are as follows:
     * <pre> <code>
     * <table class= "myColumnPadding">
     *  <tr>
     *    <th>Description</th>
     *    <th>Value</th>
     *  </tr>
     *  <tr>
     *    <td>boolean</td>
     *    <td>B</td>
     *  </tr>
     *  <tr>
     *    <td>real</td>
     *    <td>R</td>
     *  </tr>
     *  <tr>
     *    <td>ordinal</td>
     *    <td>O</td>
     *  </tr>
     *  <tr>
     *    <td>categorical</td>
     *    <td>C</td>
     *  </tr>
     *  <tr>
     *    <td>integer</td>
     *    <td>I</td>
     *  </tr>
     * </table>
     * </code> </pre> 
     *
     * @return The x-variable types.
     */
    public char[] get_xType() {

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {
        @RF_TRACE_OFF@    for (int i = 0; i < xType.length; i++) {
        @RF_TRACE_OFF@      RFLogger.log(Level.INFO, "xType[" + i + "] = " + xType[i]);            
        @RF_TRACE_OFF@    }
        @RF_TRACE_OFF@  }

        return xType;
    }

    /** 
     * Returns a vector of length ySize ({@link #get_ySize}) containing the the number of levels found in each
     * y-variable.  The elements of this vector will be non-zero for ordinal and categorical
     * variables.  All others elements assume the value of zero (0).  
     * @return The number of levels found in each y-variable
     */
    public int[] get_yLevel() {

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {        
        @RF_TRACE_OFF@    for (int i = 0; i < yType.length; i++) {
        @RF_TRACE_OFF@      RFLogger.log(Level.INFO, "yLevel[" + i + "] = " + yLevel[i]);            
        @RF_TRACE_OFF@    }
        @RF_TRACE_OFF@  }

        return yLevel;
    }

    /**
     * Returns a vector of length xSize ({@link #get_xSize}) containing the the number of levels found in each
     * x-variable.  The elements of this vector will be non-zero for ordinal and categorical
     * variables.  All others elements assume the value of zero (0).
     * @return The number of levels found in each x-variable.
     */
    public int[] get_xLevel() {

        @RF_TRACE_OFF@  if (Trace.get(Trace.LOW)) {        
        @RF_TRACE_OFF@    for (int i = 0; i < xType.length; i++) {
        @RF_TRACE_OFF@      RFLogger.log(Level.INFO, "xLevel[" + i + "] = " + xLevel[i]);            
        @RF_TRACE_OFF@    }
        @RF_TRACE_OFF@  }

        return xLevel;
    }

    /**
     * Returns a 2-D matrix of values representing the y-values. The
     * dimensions of this matrix are [ySize] x [nSize].  
     * Note that the boolean, categorical, and ordinal
     * variable types will have been mapped to real values.
     * @return The 2-D matrix of y-values.
     */
    public double[][] get_yData() {
        return yData;
    }

    /**
     * Returns a 2-D matrix of values representing the y-values.  The
     * dimensions of this matrix are [xSize] x [nSize].  
     * Note that the boolean, categorical, and ordinal
     * variable types will have been mapped to real values.
     * @return The 2-D matrix of x-values.
     */
    public double[][] get_xData() {
        return xData;
    }
    

    /**
     * Sets the bootstrap related parameters in the model.
     *
     * <p> Note that the parameters are interdependent on one another.
     * An explanation of the heirarchy, the interdependency, and
     * default values is below.  Only specific combinations are valid.
     * We provide two other methods to set the bootstrap related
     * parameters: {@link #set_bootstrap()} and {@link
     * #set_bootstrap(int)} to aid the user. When explicitly specifying the sample to
     * be used in the bootstrap algorithm, (sampleType = user), the
     * element sample[i][j] represents the number of times case j (in
     * the data set) appears in the bootstrap sample for tree i.  Thus
     * a value of zero (0) implies that case j is out-of-bag in tree
     * i. Ensure that the sample over the forest is coherent:  the sum of the
     * each column should equal sampleSize. </p>
     * <pre> <code> 
     * <table class= "myColumnPadding">
     *  <tr>
     *    <th>Parameter</th>
     *    <th>Default Value</th>
     *    <th>Possible Values</th>
     *  </tr>
     *  <tr>
     *    <td>ntree</td>
     *    <td>1000</td>
     *    <td> &gt; 0</td>
     *  </tr>
     *  <tr>
     *    <td>bootstrap</td>
     *    <td>auto </td>
     *    <td>auto, user </td>
     *  </tr>
     *  <tr>
     *    <td>sampleType</td>
     *    <td>swr</td>
     *    <td>swr (sampling with replacement), swor (sampling without replacement)</td>
     *  </tr>
     *  <tr>
     *    <td>sampleSize</td>
     *    <td>n</td>
     *    <td>&gt; 0</td>
     *  </tr>
     *  <tr>
     *    <td>sample</td>     
     *    <td>null</td>
     *    <td>null, 2-D matrix of dimension [ntree] x [nSize]</td>
     *  </tr>
     * </table></p>
     *
     *                                                                                       >0
     *                                                                swr                  /------ sampleSize
     *                                                              /------ sampleSize? --/                         
     *                                                             /                      \
     *                                      auto                  /                        \------ (sampleSize = n)
     *                >0                  /------ sampleType? ---/                           =0
     *               ------ bootstrap? --/                       \                         
     *              /                    \                        \                           1 <= sampleSize <= n
     *   ntree?  --/                      \                        \                        /----- sampleSize
     *             \                       \                        \------ sampleSize? -- /
     *              \------ WARNING         \                         swor                 \
     *                =0                     \                                              \----- (sampleSize = n * (e-1)/e)
     *                                        \                                               =0
     *                                         \
     *                                          \                     !null
     *                                           \                  /------ (sampleType = swr) && 
     *                                            ------ sample? --/        (sampleSize determined by sample) && 
     *                                             user            \        (ntree determined by sample)
     *                                                              \
     *                                                               \------ ERROR 
     *                                                                 null                                              
     *
     * </code> </pre>
     * Finally, note that when bootstrap = auto is in force, it is also
     * possible to use case weights in conjuntion with sampleType.  
     * See {@link #set_caseWeight(double[])} for more details.
     *
     * @param ntree Number of trees in the forest.
     * @param bootstrap Type of bootstrap used in the model.
     * @param sampleType Type of sampling used in generating the bootstrap.
     * @param sampleSize Size of sample used in generating the bootstrap.
     * @param sample 2-D matrix explicitly specifying the bootstrap sample.
     */
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
                    // Do nothing.  The sample size is valid.
                }
                else {
                    // Ignore value and use the default sample size.
                    // RFLogger.log(Level.WARNING, "Invalid value for bootstrap parameter sampleSize:  " + sampleSize);
                    // RFLogger.log(Level.WARNING, "Overriding sampleSize with default value:  " + nSize);
                    this.sampleSize = nSize;
                }

            }
            else if (sampleType.equals("swor")) {
                this.sampleType = sampleType;
                if ((sampleSize > 0) && (sampleSize <= nSize)) {
                    // Do nothing.  The sample size is valid.
                    this.sampleSize = sampleSize;
                }
                else {
                    // Ignore value and use the default sample size.
                    // RFLogger.log(Level.WARNING, "Invalid value for bootstrap parameter sampleSize:  " + sampleSize);
                    // RFLogger.log(Level.WARNING, "Overriding sampleSize with default value:  [n * (e-1)/e]");
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

            // Check for coherent sample.  It will of be dim [ntree] x [nSize].
            if (sample == null) {
                RFLogger.log(Level.SEVERE, "sample must not be null when bootstrapping by user.");
                throw new IllegalArgumentException();
            }
            else {
                int[] sSize = new int[ntree];

                if (sample.length != ntree) {
                    RFLogger.log(Level.SEVERE, "sample[.][] must be specified for each tree:  " +  sample.length);
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
                if (sampleSize > 0) {
                    if (sampleSize != sSize[0]) {
                        RFLogger.log(Level.SEVERE, "sample size and sample are inconsistent:  " +  sampleSize + " vs. " + sSize[0]);
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

    /** Sets the bootstrap related parameters in the model.  Use default values for all unspecified parameters.
     * @param ntree Number of trees in the random forest.
     * @see #set_bootstrap(int, String, String, int, int[][])
     */
    public void set_bootstrap(int ntree) {
        this.ntree = 1000;
        bootstrap = "auto";
        sampleType = "swr";
        sampleSize = nSize;
        sample = null;
    }
    
    /** Sets the bootstrap related parameters in the model.  Use default values for all unspecified parameters.
     * @param ntree Number of trees in the random forest.
     * @param sample 2-D matrix explicitly specifying the bootstrap sample.
     * @see #set_bootstrap(int, String, String, int, int[][])
     */
    public void set_bootstrap(int ntree, int[][] sample) {
        set_bootstrap(ntree, "user", "swr", 0, sample); 
    }


    /** Sets the bootstrap related parameters in the model.  Use default values for all parameters.
     * @see #set_bootstrap(int, String, String, int, int[][]) 
     */
    public void set_bootstrap() {
        ntree = 1000;
        bootstrap = "auto";
        sampleType = "swr";
        sampleSize = nSize;
        sample = null;
    }

    /** 
     * Returns the type of bootstrap used in the model. 
     * @return The type of bootstrap used in the model.
     * @see #set_bootstrap(int, String, String, int, int[][])
     */
    public String get_bootstrap() {
        return bootstrap;
    }

    /**
     * Returns the type of sampling used in generating the bootstrap.
     * @return The type of sampling used in generating the bootstrap.
     * @see #set_bootstrap(int, String, String, int, int[][])
     */
    public String get_sampleType() {
        return sampleType;
    }

    /**
     * Returns the size of sample used in generating the bootstrap.
     * @return The size of sample used in generating the bootstrap.
     * @see #set_bootstrap(int, String, String, int, int[][])
     */
    public int get_sampleSize() {
        return sampleSize;
    }

    /**
     * Returns the 2-D matrix explicitly specifying the bootstrap sample.
     * @return The 2-D matrix explicitly specifying the bootstrap sample.
     * @see #set_bootstrap(int, String, String, int, int[][])
     */
    public int[][] get_sample() {
        return sample;
    }

    /**
     * Returns the number of trees in the forest.
     * @return The number of trees in the forest.
     * @see #set_bootstrap(int, String, String, int, int[][])
     */
    public int get_ntree() {
        return ntree;
    }


    /**
     * Sets the number of x-variables to be randomly selected as candidates for splitting a node.
     * @param mtry The number of x-variables to be randomly selected as candidates for splitting a node. 
     * This number must be such that 1 &le; mtry &le; xSize.
     * If the value is out of range, the default value will be applied:
     * <pre> <code> 
     * <table class= "myColumnPadding">
     *  <tr>
     *    <th>Family</th>
     *    <th>Default Value</th>
     *  </tr>
     *  <tr>
     *    <td>RF-R, RF-R+</td>
     *    <td>xSize/3</td>
     *  </tr>
     *  <tr>
     *    <td>all others</td>     
     *    <td>sqrt(xSize) </td>
     *  </tr>
     * </table></p>
     * The default is to use uniform weights for selection, though this can be changed with {@link #set_xWeight(double[])}.
     */
    public void set_mtry(int mtry) {
        if ((mtry < 1) || (mtry > xSize)) {
            set_mtry();
        }
        else {
            this.mtry = mtry;
        }
    }

    /**
     * Sets the default value for the number of x-variables to be randomly selected as candidates for splitting a node.
     * @see #set_mtry(int)
     */
    public void set_mtry() {
        if (family.equals("RF-R") || family.equals("RF-R+")) {
            mtry = (int) Math.ceil((double) xSize / 3.0); 
        }
        else {
            mtry = (int) Math.ceil(Math.sqrt((double) xSize));
        }
    }

    /** 
     * Returns the value for the number of x-variables to be randomly selected as candidates for splitting a node.
     * @return The value for the number of x-variables to be randomly selected as candidates for splitting a node.
     * @see #set_mtry(int)
     */
    public int get_mtry() {
        return mtry;
    }

    /** 
     * Sets the number of iterations for the missing data algorithm.
     * @param nImpute The number of iterations for the missing data algorithm.
     * The default value is one (1).
     * Performance measures such as out-of-bag (OOB) error rates tend
     * to become optimistic if nimpute &gt; 1.  
     */
    public void set_nImpute(int nImpute) {
        this.nImpute = nImpute;
        if (nImpute < 1) {
            RFLogger.log(Level.WARNING, "Invalid value for parameter nImpute:  " + nImpute);
            RFLogger.log(Level.WARNING, "Overriding nImpute with default value:  1");
            this.nImpute = 1;    
        }

    }

    /**
     * Returns the number of iterations used by the missing data algorithm.
     * @return The number of iterations used by the missing data algorithm.
     * @see #set_nImpute(int)
     */
    public int get_nImpute() {
        return nImpute;
    }


    /**
     * Sets the case weight vector.  
     * @param weight The case weight vector.  The vector 
     * must be of length nSize ({@link #get_nSize}).  This is a vector of non-negative
     * weights where, after normalizing, weight[k] is the
     * probability of selecting case k as a candidate when bootstrap = auto.
     * The default is to use uniform weights for selection.  It is generally better to use real
     * weights rather than integers.  With larger values of nsize, the
     * slightly different sampling algorithms depolyed in the two
     * scenarios can result in dramatically different execution times.
     * @see #set_bootstrap(int, String, String, int, int[][])
     */
    public void set_caseWeight(double[] weight) {
        caseWeight = setWeight(weight, nSize);
    }

    /** 
     * Set the case weight vector to uniform weights.
     * @see #set_caseWeight(double[])
     */
    public void set_caseWeight() {
        caseWeight = setWeight(null, nSize);
    }

    /**
     * Returns the case weight vector.
     * @return The case weight vector.
     * @see #set_caseWeight(double[])
     */
    public double[] get_caseWeight() {
        return caseWeight;
    }
    

    /**
     * Sets the number of randomly selected pseudo-responses when
     * unsupervised forests is in force.  
     * @param ytry The number of randomly selected pseudo-responses when
     * unsupervised forests is in force. The default value is one
     * (1).  This means at every node, and every split attempt, one y-variable will be
     * selected from the (xSize - 1) remaining x-variables when calculating the split statistic.
     * This number must be such that 1 &lt; ytry &le; (xSize - 1).
     */
    public void set_ytry(int ytry) {
        if (family.equals("RF-U")) {
            RFLogger.log(Level.SEVERE, "ytry must be specified via the formula when family is RF-U.");
            throw new IllegalArgumentException();
        }
        this.ytry = ytry;

        if (ytry > 0) {
            // Default weights for y-vars.  The user can override as necessary.
            set_yWeight();
        }
        else {
            // ytry and yWeight must be coherently related.
            yWeight = null;
        }
    }

    /** 
     * Returns the value for the number of y-variables to be randomly selected as pseudo-responses when unsupervised forests is in force.
     * @return Yhe value for the number of y-variables to be randomly selected as pseudo-responses when unsupervised forests is in force.
     * @see #set_ytry(int)
     */
    public int get_ytry() {
        return ytry;
    }

    
    /**
     * Experimental - Do Not Use.
     */
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

    /**
     * Experimental - Do Not Use.
     */
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

    /**
     * Experimental - Do Not Use.
     */
    public double[] get_yWeight() {
        return yWeight;
    }

    /**
     * Sets the x-variable weight vector.  
     * @param weight The x-variable weight vector.  This vector 
     * must be of length xSize ({@link #get_xSize}).  This is a vector of non-negative
     * weights where, after normalizing, weight[k] is the
     * probability of selecting x-variable k as a candidate for splitting a node.
     * The default is to use uniform weights for selection.  It is generally better to use real
     * weights rather than integers.  With larger values of nsize, the
     * slightly different sampling algorithms depolyed in the two
     * scenarios can result in dramatically different execution times.
     */
    public void set_xWeight(double[] weight) {
        xWeight = setWeight(weight, xSize);
    }

    /** 
     * Set the x-variable weight vector to uniform weights.
     * @see #set_xWeight(double[])
     */
    public void set_xWeight() {
        xWeight = setWeight(null, xSize);
    }

    /**
     * Returns the x-variable weight vector.
     * @return The x-variable weight vector.
     * @see #set_xWeight(double[])
     */
    public double[] get_xWeight() {
        return xWeight;
    }

    /**
     * Sets the x-variable statistical weight vector.  
     * @param weight The x-variable statistical weight vector.  This vector  must be of
     * length xSize ({@link #get_xSize}). This is a vector of
     * non-negative weights where, after normalizing, weight[k] is the
     * multiplier by which the split statistic for an x-variable is
     * adjusted.  A large value encourages the node to split on the
     * x-variable. The default is to use uniform weights so that all
     * x-variables are treated equally.
    */
    public void set_xStatisticalWeight(double[] weight) {
        xStatisticalWeight = setWeight(weight, xSize);
    }

    /** 
     * Set the x-variable statistical weight vector to uniform weights.
     * @see #set_xStatisticalWeight(double[])
     */
    public void set_xStatisticalWeight() {
        xStatisticalWeight = setWeight(null, xSize);
    }

    /**
     * Returns the x-variable statistical weight vector.
     * @return The x-variable statistical weight vector.
     * @see #set_xStatisticalWeight(double[])
     */
    public double[] get_xStatisticalWeight() {
        return xStatisticalWeight;
    }

    /**
     * Sets the event weight vector, when survival or competing risk
     * forests is in force.  
     * @param weight The event weight vector, when survival or competing risk forests is in force.  This vector must be the same length as the
     * number of events in the data set ({@link #get_eventCount}). 
     * This is a vector of non-negative weights,
     * where, after normalizing, weight[k] is the multiplier by which
     * the component of the split statistic related to event[k] is adjusted.
     * The default is to to use a composite splitting rule which is an average over all event types (a democratic approach).
     * To single out an event type, set all weights other than the one you are interested in to zero (0).
     * Finally, note that regardless of how the weight vector is specified, the returned forest object always
     * provides estimates for all event types.
     */
    public void set_eventWeight(double[] weight) {

        // Note that the family must be RF-S, and eventType must be
        // initialized for this function to be called. The latter is
        // nominally the case, and things must be seriously corrupted
        // if this is not the case.
        if (family.equals("RF-S")) {
            crWeight = setWeight(weight, eventType.size());
        }
        else {
            RFLogger.log(Level.SEVERE, "family must be RF-S in order to define this weight vector.");
            throw new IllegalArgumentException();
        }
        
        crWeightSize = eventType.size();
    }

    /**
     * Sets the event weight vector, when survival or competing risk
     * forests is in force, to uniform weights.
     * @see #set_eventWeight(double[])
     */
    public void set_eventWeight() {
        
        // Note that the family must be RF-S, and eventType must be
        // initialized for this function to be called. The latter is
        // nominally the case, and things must be seriously corrupted
        // if this is not the case.
        if (family.equals("RF-S")) {
            crWeight = setWeight(null, eventType.size());
        }
        else {
            RFLogger.log(Level.SEVERE, "family must be RF-S in order to define this weight vector.");
            throw new IllegalArgumentException();
        }

        crWeightSize = eventType.size();
    }

    /**
     * Returns the event weight vector, when survival or competing risk
     * forests is in force.
     * @return The event weight vector, when survival or competing risk
     * forests is in force.
     * @see #set_eventWeight(double[])
     */
    public double[] get_eventWeight() {
        return crWeight;
    }

    /**
     * Returns the number of events in the data set, when survival or competing risk forests is in force.
     * @return The number of events in the data set, when survival or competing risk forests is in force.
     */
    public int get_eventCount() {
        return crWeightSize;
    }
        

    /**
     * Sets the time interest vector, when survival or competing risk
     * forests is in force.  
     * @param timeInterest The time interest vector, when survival or competing risk
     * forests is in force.  This is a vector of real values to be
     * used to constrain the ensemble calculations. Using time points
     * at which events do not occur does not result in information
     * gain.  The default action is to use all observed event times in
     * the data set.  
     */
    public void set_timeInterest(double[] timeInterest) {

        Iterator <Double> itr;        
        double value;
        int j;
        
        TreeSet <Double> inTime = new TreeSet <Double> (); 

        // Note that yData is of the form yData[ySize][nSize];

        // Time is always in the [0][...] row.
        // Cens is always in the [0][...] row.  

        // First, we order and uniquify the user time interest vector.
        for (int i = 0; i < timeInterest.length; i++) {
            inTime.add((double) timeInterest[i]);
        }

        // Next we get the data-natural time interest vector.
        set_timeInterest();
        double[] naturalTimeInterest = timeInterest;
        
        // Next we round-up all values in the user vector to values in
        // the data-natural vector and initializing a new TreeSet with
        // these values.  They are finally converted to anative-code
        // compliant double vector.
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

    /**
     * Sets the time interest vector, when survival or competing risk
     * forests is in force to the default value.
     * @see #set_timeInterest(double[])
     */
    public void set_timeInterest() {
        
        TreeSet <Double> treeSet = new TreeSet <Double> (); 

        // Note that yData is of the form yData[ySize][nSize];

        // Time is always in the [0][...] row.
        // Cens is always in the [0][...] row.  

        for (int i = 0; i < nSize; i++) {
            // Did an event occur?
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

    /**
     * Returns the time interest vector used in the model, when survival or competing risk
     * forests is in force.
     * @return The time interest vector used in the model, when survival or competing risk
     * forests is in force.
     * @see #set_timeInterest(double[])
     */
    public double[] get_timeInterest() {
        return timeInterest;
    }

    /**
     * Returns the size of the time interest vector used in the model, when survival or competing risk
     * forests is in force.
     * @return The size of the time interest vector used in the model, when survival or competing risk
     * forests is in force.
     * @see #get_timeInterest()
     */
    public int get_timeInterestSize() {
        return timeInterestSize;
    }

    /**
     * Sets the split rule to be used in generating the model.  
     * @param splitRule The split rule to be used in generating the model.  The split rules available are detailed below. The   
     * rule in bold denotes the default split rule for each
     * family. The default split rule is applied when the user does
     * not specify a split rule. Survival and
     * Competing Risk both have two split rules. Regression has three
     * flavours of split rules based on mean-squared
     * error. Classification also has three flavours of split rules
     * based on the Gini index. The Multivariate and Unsupervised
     * split rules are a composite rule based on Regression and
     * Classification. Each component of the composite is normalized
     * so that the magnitude of any one y-variable does not influence
     * the statistic. All families also allow the user to define a
     * custom split rule statistic. Some basic C-programming skills
     * are required. Examples for all the families reside in the C
     * source code directory of the package in the file
     * <code>src/main/c/splitCustom.c</code>. Note that recompiling and
     * re-installing the package is necessary after modifying the
     * source code.
     * 
     * <pre> <code> 
     * <table class= "myColumnPadding">
     *  <tr>
     *    <th>Family</th>
     *    <th>Split Rule Description</th>
     *    <th>Value</th>
     *  </tr>
     *  
     *  <tr>
     *    <td rowspan="2">survival</td>
     *    <td>log-rank</td>
     *    <td><b>logrank</b></td>
     *  </tr>
     *  <tr>
     *    <td>log-rank score</td>
     *    <td>logrankscore</td>
     *  </tr>
     *  
     *  <tr>
     *    <td rowspan="2">competing risk</td>
     *    <td>log-rank modified weighted</td>
     *    <td><b>logrankCR</b></td>
     *  </tr>
     *
     *  <tr>
     *    <td>log-rank</td>
     *    <td>logrankACR</td>
     *  </tr>
     *
     *  <tr style="background-color:#eaeaea">
     *    <td rowspan="3">regression</td>
     *    <td>mean-squared error weighted</td>
     *    <td><b>mse</b></td>
     *  </tr>
     *
     *  <tr style="background-color:#eaeaea">
     *    <td>mean-squared error unweighted </td>
     *    <td>mse.unwt</td>
     *  </tr>
     *
     *  <tr style="background-color:#eaeaea">
     *    <td>mean-squared error heavy weighted </td>
     *    <td>mse.hvwt</td>
     *  </tr>
     *
     *  <tr style="background-color:#eaeaea">
     *    <td rowspan="3">classification</td>
     *    <td>Gini index weighted</td>
     *    <td><b>gini</b></td>
     *  </tr>
     *
     *  <tr style="background-color:#eaeaea">
     *    <td>Gini index unweighted </td>
     *    <td>gini.unwt</td>
     *  </tr>
     *
     *  <tr style="background-color:#eaeaea">
     *    <td>Gini index heavy weighted </td>
     *    <td>gini.hvwt</td>
     *  </tr>
     *
     *  <tr style="background-color:#d6d6d6">
     *    <td>multivariate regression</td>
     *    <td>Composite mean-squared error</td>
     *    <td>mv.mse</td>
     *  </tr>
     *  
     *  <tr style="background-color:#d6d6d6">
     *    <td>multivariate classification</td>
     *    <td>Composite Gini index</td>
     *    <td>mv.gini</td>
     *  </tr>
     *  
     *  <tr style="background-color:#b7b7b7">
     *    <td>multivariate mixed</td>
     *    <td>Composite Gini and MSE</td>
     *    <td>mv.mix</td>
     *  </tr>
     *
     *  <tr style="background-color:#b7b7b7">
     *    <td>unsupervised</td>
     *    <td>pseudo-response adaptive</td>
     *    <td>unsupv</td>
     *  </tr>
     *  
     *</table>
     * </code> </pre> 
     */
    public void set_splitRule(String splitRule) {

        boolean result = false;

        // Safe the custom split rule.
        customSplitIndex = 0;

        if (splitRule.equals("random")) {
            this.splitRule = splitRule;
            result = true;
        }
        else if (splitRule.startsWith("custom")) {
            String numberString = splitRule.replaceFirst("custom", "");
            // No numeric is the equivalent of the first (and potentially only) custom split rule. 
            if (numberString.equals("")) {
                this.splitRule = splitRule;
                customSplitIndex = 1;
                result = true;
            }
            else {
                // Extract the numeric value.
                customSplitIndex = Integer.parseInt(numberString);
                // Check that it is in range from [1, 16].
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
    
    /**
     * Sets the default split rule, based on the data set and formuala.
     * @see #set_splitRule(String)
     */
    public void set_splitRule() {

        // Safe the custom split rule.
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


    /**
     * Returns the split rule used in generating the model.
     * @return The split rule used in generating the model.
     * @see #set_splitRule(String)
     */
    public String get_splitRule() {
        return splitRule;
    }

    // Used by RandomForestModel.train() to pass the parameter to the native-code.
    int getSplitRuleID(String splitRule) {
        return splitRuleID.get(splitRule);
    }


    /**
     * Sets the parameter specifying deterministic versus non-deterministic splitting.  
     * @param nSplit The parameter specifying deterministic versus non-deterministic splitting.  The parameter must be 
     * a non-negative integer value.  When zero (0), deterministic
     * splitting for an x-variable is in force.  When non-zero, a
     * maximum of nSplit points are randomly chosen among the
     * possible split points for an x-variable. This can
     * significantly decrease computation time over deterministic splitting.  The
     * default value for this parameter varies with the split rule:  When pure
     * random splitting is in force, the default and only value for this parameter is one (1).
     * When any other split rule is in force, the default value is
     * zero (0).
     */
    public void set_nSplit(int nSplit) {
        this.nSplit = nSplit;
        if (splitRule.equals("random")) {
            if (nSplit != 1) {
                RFLogger.log(Level.WARNING, "Overriding parameter nSplit with (1) when pure random splitting.");
            }
            this.nSplit = 1;
        }
    }

    /**
     * Sets the default value for the parameter specifying deterministic versus non-deterministic splitting.
     * @see #set_nSplit(int)
     */
    public void set_nSplit() {
        if (splitRule.equals("random")) {
            nSplit = 1;
        }
        else {
            // Default value is deterministic splitting.
            nSplit = 0;
        }
    }


    /**
     * Returns the parameter specifying deterministic versus non-deterministic splitting.
     * @return The parameter specifying deterministic versus non-deterministic splitting.
     * @see #set_nSplit(int)
     */
    public int get_nSplit() {
        return nSplit;
    }

    /** 
     * Sets the desired average number of unique cases in a terminal
     * node.  
     * @param nodeSize The desired average number of unique cases in a terminal
     * node. The parameter ensures that the average nodesize across
     * the forest will be at least nodeSize.  Some nodes will be
     * smaller than this value and some will be larger.  The default
     * value for this parameter varies with the family, though it
     * recommended to experiment with different values.
     *
     * <pre> <code> 
     * <table class= "myColumnPadding">
     *  <tr>
     *    <th>Family</th>
     *    <th>Defalut Node Size</th>
     *  </tr>
     *  
     *  <tr>
     *    <td>survival</td>
     *    <td>3</td>
     *  </tr>
     *
     *  <tr>
     *    <td>competing risk</td>
     *    <td>6</td>
     *  </tr>
     *
     *  <tr style="background-color:#eaeaea">
     *    <td>regression</td>
     *    <td>5</td>
     *  </tr>
     *  <tr style="background-color:#eaeaea">
     *    <td>multivariate regression</td>
     *    <td>5</td>
     *  </tr>
     *  
     *
     *  <tr style="background-color:#d6d6d6">
     *    <td>classification</td>
     *    <td>1</td>
     *  </tr>
     *
     *  <tr style="background-color:#d6d6d6">
     *    <td>multivariate classification</td>
     *    <td>1</td>
     *  </tr>
     *  
     *  <tr style="background-color:#b7b7b7">
     *    <td>multivariate mixed</td>
     *    <td>3</td>
     *  </tr>
     *
     *  <tr style="background-color:#b7b7b7">
     *    <td>unsupervised</td>
     *    <td>3</td>
     *  </tr>
     *  
     *</table>
     * </code> </pre> 
     */
    public void set_nodeSize(int nodeSize) {
        if (nodeSize < 0) {
            RFLogger.log(Level.SEVERE, "nodeSize must be greater than zero:  " + nodeSize);
            throw new IllegalArgumentException();
        }
        this.nodeSize = nodeSize;
    }

    /** 
     * Sets the default value for the average number of unique cases in a terminal node.    
     * @see #set_nodeSize(int)
     */
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
   
    /** 
     * Returns the desired value for the average number of unique cases in a terminal node.    
     * @return The desired value for the average number of unique cases in a terminal node.    
     * @see #set_nodeSize(int)
     */
    public int get_nodeSize() {
        return nodeSize;
    }

    /**
     * Sets the maximum depth to which a tree should be grown.  
     * @param nodeDepth The maximum depth to which a tree should be grown. The
     * default behaviour is that this parameter is ignored. Not
     * setting this parameter or setting this parameter to a negative
     * value will ensure that this parameter is ignored.
     */
    public void set_nodeDepth(int nodeDepth) {
        this.nodeDepth = nodeDepth;
    }

    /**
     * Returns the maximum depth to which a tree should be grown.
     * @return The maximum depth to which a tree should be grown.
     */
    public int get_nodeDepth() {
        return nodeDepth;
    }

    /**
     * Sets the parameter to control whether the null hypothesis is used in generating the model.  
     * @param nullSplit The parameter to control whether the null hypothesis is used in generating the model.  The null
     * hypothesis assumes no relation between x-variables and
     * y-variables.  To test this, before growing each tree, we
     * permute the each y-variable independently over the non-missing
     * cases.  When the family is survival or competing risk, we
     * restrict permutation to time, and leave censoring untouched.  When the family is
     * unsupervised, this parameter is (no).  The default value
     * for this parameter is (no).
     */
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

    /**
     * Returns the parameter to control whether the null hypothesis is used in generating the model.
     * @return The parameter to control whether the null hypothesis is used in generating the model.
     * @see #set_nullSplit(String)
     */
    public String get_nullSplit() {
        return nullSplit;
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
        seed = - generator.nextInt();
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
     * Sets the ensemble outputs desired from the model.  
     * These
     * settings are in the form of &lt;key, value&gt; pairs, where the key
     * is the name of the ensemble, and the value is the specific
     * option for that ensemble.  
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
     *    <td><b>no</b>, permute, random, permute.ensemble, random.ensemble</td>
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
     * @param key The name of the ensemble output.
     * @param value The specific value for the ensemble output. 
     */
    public void set_ensembleArg(String key, String value) {
        ensembleArg.set(key, value);
    }

    /**
     * Sets the default values for ensemble outputs desired from the model.  
     * @see #set_ensembleArg(String, String).
     */
    public void set_ensembleArg() {
        ensembleArg.set();
    }

    /**
     * Returns the default or user defined option for the specified ensembled.
     * @param key The name of the ensemble output.
     * @see #set_ensembleArg(String, String).
     */
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
            // Shift the decremented value by 8 bits.
            customSplitBit = (customSplitIndex - 1) << 8;
        }

        result += customSplitBit;
        
        return result;
    }


}

