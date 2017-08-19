package com.kogalur.randomforest;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Iterator;
import java.util.Set;
import java.util.Map;

import org.apache.spark.sql.Dataset;
import org.apache.spark.sql.Row;

import org.apache.spark.sql.SparkSession;

import org.apache.spark.sql.types.DataTypes;
import org.apache.spark.sql.types.DataType;
import org.apache.spark.sql.types.StructField;
import org.apache.spark.sql.types.StructType;

import org.apache.spark.ml.feature.StringIndexer;
import org.apache.spark.ml.feature.StringIndexerModel;

import scala.collection.Seq;

public class JavaSQLDataSourceExampleNew {

  public static void main(String[] args) {
    SparkSession spark = SparkSession
        .builder()
        .appName("Java Spark SQL data sources example")
        .config("spark.master", "local")
        .getOrCreate();

    runBasicDataSourceExample(spark);

    spark.stop();
  }

   private static void runBasicDataSourceExample(SparkSession spark) {

      Dataset<Row> usersDF = spark
          .read()
          .option("header", "true")
          .option("inferSchema", "true") 
          .format("csv")
          .load("/Users/kogalur/Dropbox/working/rfsrc/scala/mtcars.csv");

      usersDF.printSchema();


      StringIndexerModel indexer;
      
      indexer = new StringIndexer()
          .setInputCol("carbF")
          .setOutputCol("carbIndex")
          .fit(usersDF);

      Dataset<Row> indexed = indexer.transform(usersDF);
      
      System.out.println("Transformed string column '" + indexer.getInputCol() + "' " +
                         "to indexed column '" + indexer.getOutputCol() + "'");

      indexed = indexed.drop(indexed.col("carbF"));

      indexer = new StringIndexer()
          .setInputCol("vsF")
          .setOutputCol("vsIndex")
          .fit(indexed);

      indexed = indexer.transform(indexed);
      
      System.out.println("Transformed string column '" + indexer.getInputCol() + "' " +
                         "to indexed column '" + indexer.getOutputCol() + "'");
      
      indexed = indexed.drop(indexed.col("vsF"));

      indexed = indexed.withColumnRenamed("vsIndex", "vsF");
      indexed = indexed.withColumnRenamed("carbIndex", "carbF");
     

      String[] names = indexed.columns();
      for (int i = 0; i < names.length; i++) {
          System.out.println("Column Names " + i + " : " + names[i]);
      }


      indexed = indexed.withColumn("vsF", indexed.select("vsF") + 1);




      boolean complimentFlag;
      boolean bigRFlag;
      
      String   formulaThis;
      String[] formula;
      String[] formulaX;
      String[] formulaY;
      
      String formulaR = new String(" mpg ~ .  ");
      String formulaM1 = new String(" Multivariate (mpg, wt) ~ hp + drat");
      String formulaM2 = new String(" Multivariate (mpg, wt) ~ . ");
      String formulaS = new String(" Surv (days, status) ~ .");
      String formulaU = new String(" Unsupervised () ~ .");

      formulaThis = formulaM2;

      formulaThis = formulaThis.replaceAll("\\ ", "");
      
      formula = formulaThis.split("~");

      
      if (formula.length != 2) {
          // Bad Formula.  Throw exception.
          System.out.println("\n Unknown Formula Syntax.");
      }
      
      if (formula[0].startsWith("Multivariate")) {
          // Check for proper syntax.
          if ((formula[0].indexOf('(') >= 0) && (formula[0].indexOf(')') >= 0)) {
              formula[0] = formula[0].replace("Multivariate", "");
              formula[0] = formula[0].replace("(", "");
              formula[0] = formula[0].replace(")", "");
              formulaY = formula[0].split(",");
          }
          else {
              // Bad Formula.  Throw exception.
              formulaY = null;
              System.out.println("\n Bad Multivariate Formula Syntax.");
          }
      }
      else if (formula[0].startsWith("Surv")) {
          // Check for proper syntax.
          if ((formula[0].indexOf('(') >= 0) && (formula[0].indexOf(')') >= 0)) {
              formula[0] = formula[0].replace("Surv", "");
              formula[0] = formula[0].replace("(", "");
              formula[0] = formula[0].replace(")", "");
              formulaY = formula[0].split(",");
          }
          else {
              // Bad Formula.  Throw exception.
              formulaY = null;
              System.out.println("\n Bad Survival Formula Syntax.");
          }
      }
      else if (formula[0].startsWith("Unsupervised")) {
          // Check for proper syntax.
          if ((formula[0].indexOf('(') >= 0) && (formula[0].indexOf(')') >= 0)) {
              formula[0] = formula[0].replace("Unsupervised", "");
              formula[0] = formula[0].replace("(", "");
              formula[0] = formula[0].replace(")", "");
              formulaY = formula[0].split(",");

              if (formulaY != null) {
                  // Bad Formula.  Throw exception.
                  formulaY = null;
                  System.out.println("\n Bad Unsupervised Formula Syntax.");
              }
          }
          else {
              // Bad Formula.  Throw exception.
              formulaY = null;
              System.out.println("\n Bad Unsupervised Formula Syntax.");
          }
      }
      else {
          // Univariate response.
          formulaY = new String[1];
          formulaY[0] = formula[0];
      }

      // Detect whether we are in a big-r situation.  In such cases,
      // we will acquire the predictor data frame via unpacking rather
      // than dropping each response in a loop.
      bigRFlag = false;
      if (formulaY != null) {
          if (formulaY.length > 50) {
              bigRFlag = true;
          }
      }
      
      complimentFlag = false;
      if (formula[1].equals(".")) {
          // x-vars are complement of y-vars.

          complimentFlag = true;
          
          if (formulaY != null) {
              // The length of the x-var names is known.
              formulaX = new String[names.length - formulaY.length];

              // Convert the data frame names to a list for easy removal of the response names.
              List<String> formulaXL = new ArrayList<String>(Arrays.asList(names));
              for (int j = 0; j < formulaY.length; j++) {
                  formulaXL.remove(formulaY[j]);
              }
              for (int i = 0; i < formulaX.length; i++) {
                  formulaX[i] = formulaXL.get(i);
              }
          }
          else {
              // We are in an unsupervised scenario.  The x-var names will the be data frame names.
              formulaX = new String[names.length];
              for (int i = 0; i < formulaX.length; i++) {
                  formulaX[i] = names[i];
              }

          }
      }
      else {
          // x-vars are explicity specified. Parse the formula and extract them.
          formulaX = formula[1].split("\\+");

      }
      
      
      for (int i = 0; i < formulaY.length; i++) {
          System.out.println("Formula Y " + i + " : " + formulaY[i]);
      }
      
      for (int i = 0; i < formulaX.length; i++) {
          System.out.println("Formula X " + i + " : " + formulaX[i]);
      }

      // We have defined the x-vars and the y-vars.  We now proceed to
      // extract the x-vars Dataset and y-vars Dataset.
      
      Dataset<Row> predictor;

      if (!complimentFlag || bigRFlag) {
          // The x-vars are explicity specified or we are in a big-r situation.  We select the x-vars from the data frame explictly.
          if (formulaX.length > 1) {
              // We have more than one predictor.  Make the array
              // names into a list.  Then pick the first one and splat
              // the rest using a Scala conversion.
              List<String> formulaXL = new ArrayList<String>(Arrays.asList(formulaX));
              formulaXL.remove(0);
              predictor = indexed.select(formulaX[0], scala.collection.JavaConverters.asScalaIteratorConverter(formulaXL.iterator()).asScala().toSeq());
          }
          else {
              // We are in a univariate scenario.  Pick the single response off the data frame.
              predictor = indexed.select(formulaX[0]);
          }
      }
      else {
          // The x-vars are the compliment of the y-vars and we are not in a big-r situation.  We drop the
          // y-vars from the data frame and keep the rest.  
          predictor = indexed;
          for (int i = 0; i < formulaY.length; i++) {
              predictor = predictor.drop(formulaY[i]); 
          }
      }

      predictor.show();
                        

      Dataset<Row> response;

      if (formulaY != null) {
          if (formulaY.length > 1) {
              // We are in a multivariate scenario.  Make the array
              // names into a list.  Then pick the first one and splat
              // the rest using a Scala conversion.
              List<String> formulaYL = new ArrayList<String>(Arrays.asList(formulaY));
              formulaYL.remove(0);
              response = indexed.select(formulaY[0], scala.collection.JavaConverters.asScalaIteratorConverter(formulaYL.iterator()).asScala().toSeq());
          }
          else {
              // We are in a univariate scenario.  Pick the single response off the data frame.
              response = indexed.select(formulaY[0]);
          }
          response.show();
      }
      
      // We now create the x-vars and y-vars array of character types.

      // Access the schema of the incoming data set.
      StructType thisSchema = usersDF.schema();

      java.util.HashMap<String, Character> rTypes;
      java.util.HashMap<String, Character> xTypes;

      if (formulaY != null) {
          rTypes = new HashMap(formulaY.length);

          for (int i = 0; i < formulaY.length; i++) {
              StructField structField = thisSchema.apply(thisSchema.fieldIndex(formulaY[i]));          
              rTypes.put(structField.name(), getType(structField.dataType()));
          }
      }
      else {
          rTypes = null;
      }

      xTypes = new HashMap(formulaX.length);
      for (int i = 0; i < formulaX.length; i++) {
          StructField structField = thisSchema.apply(thisSchema.fieldIndex(formulaX[i]));
          xTypes.put(structField.name(), getType(structField.dataType()));
      }

      Set set;
      Iterator itr;
      
      if (formulaY != null) {
          set = rTypes.entrySet();
          itr = set.iterator();
      
          // Display elements.
          while(itr.hasNext()) {
              Map.Entry me = (Map.Entry) itr.next();
              System.out.println("Y Outgoing Schema:  " + me.getKey() + "  " + me.getValue());
          }
      }      
      System.out.println("\n");

      set = xTypes.entrySet();
      itr = set.iterator();
      
      // Display elements.
      while(itr.hasNext()) {
          Map.Entry me = (Map.Entry) itr.next();
          System.out.println("X Outgoing Schema:  " + me.getKey() + "  " + me.getValue());
      }
      System.out.println("\n");
      
      

      
      // DataFrames can be saved as Parquet files, maintaining the schema information
      indexed.write().mode("overwrite").parquet("mtcars.parquet");

      indexed.printSchema();
      
      // Read in the Parquet file created above.
      // Parquet files are self-describing so the schema is preserved
      // The result of loading a parquet file is also a DataFrame
      Dataset<Row> parquetFileDF = spark.read().parquet("mtcars.parquet");
      
      // The inferred schema can be visualized using the printSchema() method
      parquetFileDF.printSchema();
      // root
      //  |-- age: long (nullable = true)
      //  |-- name: string (nullable = true)

   }


    private static char getType(DataType dataType) {
        char result;
        if (dataType == DataTypes.BooleanType) {
            result = 'C';
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
        }

        return result;
    }


}
