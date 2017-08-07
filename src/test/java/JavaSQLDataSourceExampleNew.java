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
          
          System.out.println("\n Unknown Formula Syntax.");
      }
      
      if (formula[0].startsWith("Multivariate")) {
          
          if ((formula[0].indexOf('(') >= 0) && (formula[0].indexOf(')') >= 0)) {
              formula[0] = formula[0].replace("Multivariate", "");
              formula[0] = formula[0].replace("(", "");
              formula[0] = formula[0].replace(")", "");
              formulaY = formula[0].split(",");
          }
          else {
              
              formulaY = null;
              System.out.println("\n Bad Multivariate Formula Syntax.");
          }
      }
      else if (formula[0].startsWith("Surv")) {
          
          if ((formula[0].indexOf('(') >= 0) && (formula[0].indexOf(')') >= 0)) {
              formula[0] = formula[0].replace("Surv", "");
              formula[0] = formula[0].replace("(", "");
              formula[0] = formula[0].replace(")", "");
              formulaY = formula[0].split(",");
          }
          else {
              
              formulaY = null;
              System.out.println("\n Bad Survival Formula Syntax.");
          }
      }
      else if (formula[0].startsWith("Unsupervised")) {
          
          if ((formula[0].indexOf('(') >= 0) && (formula[0].indexOf(')') >= 0)) {
              formula[0] = formula[0].replace("Unsupervised", "");
              formula[0] = formula[0].replace("(", "");
              formula[0] = formula[0].replace(")", "");
              formulaY = formula[0].split(",");

              if (formulaY != null) {
                  
                  formulaY = null;
                  System.out.println("\n Bad Unsupervised Formula Syntax.");
              }
          }
          else {
              
              formulaY = null;
              System.out.println("\n Bad Unsupervised Formula Syntax.");
          }
      }
      else {
          
          formulaY = new String[1];
          formulaY[0] = formula[0];
      }

      
      
      
      bigRFlag = false;
      if (formulaY != null) {
          if (formulaY.length > 50) {
              bigRFlag = true;
          }
      }
      
      complimentFlag = false;
      if (formula[1].equals(".")) {
          

          complimentFlag = true;
          
          if (formulaY != null) {
              
              formulaX = new String[names.length - formulaY.length];

              
              List<String> formulaXL = new ArrayList<String>(Arrays.asList(names));
              for (int j = 0; j < formulaY.length; j++) {
                  formulaXL.remove(formulaY[j]);
              }
              for (int i = 0; i < formulaX.length; i++) {
                  formulaX[i] = formulaXL.get(i);
              }
          }
          else {
              
              formulaX = new String[names.length];
              for (int i = 0; i < formulaX.length; i++) {
                  formulaX[i] = names[i];
              }

          }
      }
      else {
          
          formulaX = formula[1].split("\\+");

      }
      
      
      for (int i = 0; i < formulaY.length; i++) {
          System.out.println("Formula Y " + i + " : " + formulaY[i]);
      }
      
      for (int i = 0; i < formulaX.length; i++) {
          System.out.println("Formula X " + i + " : " + formulaX[i]);
      }

      
      
      
      Dataset<Row> predictor;

      if (!complimentFlag || bigRFlag) {
          
          if (formulaX.length > 1) {
              
              
              
              List<String> formulaXL = new ArrayList<String>(Arrays.asList(formulaX));
              formulaXL.remove(0);
              predictor = indexed.select(formulaX[0], scala.collection.JavaConverters.asScalaIteratorConverter(formulaXL.iterator()).asScala().toSeq());
          }
          else {
              
              predictor = indexed.select(formulaX[0]);
          }
      }
      else {
          
          
          predictor = indexed;
          for (int i = 0; i < formulaY.length; i++) {
              predictor = predictor.drop(formulaY[i]); 
          }
      }

      predictor.show();
                        

      Dataset<Row> response;

      if (formulaY != null) {
          if (formulaY.length > 1) {
              
              
              
              List<String> formulaYL = new ArrayList<String>(Arrays.asList(formulaY));
              formulaYL.remove(0);
              response = indexed.select(formulaY[0], scala.collection.JavaConverters.asScalaIteratorConverter(formulaYL.iterator()).asScala().toSeq());
          }
          else {
              
              response = indexed.select(formulaY[0]);
          }
          response.show();
      }
      
      

      
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
      
          
          while(itr.hasNext()) {
              Map.Entry me = (Map.Entry) itr.next();
              System.out.println("Y Outgoing Schema:  " + me.getKey() + "  " + me.getValue());
          }
      }      
      System.out.println("\n");

      set = xTypes.entrySet();
      itr = set.iterator();
      
      
      while(itr.hasNext()) {
          Map.Entry me = (Map.Entry) itr.next();
          System.out.println("X Outgoing Schema:  " + me.getKey() + "  " + me.getValue());
      }
      System.out.println("\n");
      
      

      
      
      indexed.write().mode("overwrite").parquet("mtcars.parquet");

      indexed.printSchema();
      
      
      
      
      Dataset<Row> parquetFileDF = spark.read().parquet("mtcars.parquet");
      
      
      parquetFileDF.printSchema();
      
      
      

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
            
            result = 'X';
        }

        return result;
    }


}
