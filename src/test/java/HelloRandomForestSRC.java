package com.kogalur.randomforest;

import com.kogalur.randomforest.RFLogger;
import com.kogalur.randomforest.RandomForest;

import org.apache.spark.sql.SparkSession;
import org.apache.spark.sql.Dataset;
import org.apache.spark.sql.Row;

import java.util.logging.Level;

public class HelloRandomForestSRC {

    static SparkSession spark;
    
    public static void main(String[] args) {
        
        spark = SparkSession
            .builder()
            .appName("HelloRandomForestSRC Example")
            .config("spark.master", "local")
            .getOrCreate();

        // Suppress the Spark logger for all but fatal errors.
        org.apache.log4j.LogManager.getRootLogger().setLevel(org.apache.log4j.Level.FATAL);
        
        runExample();

        spark.stop();
    }

    private static void runExample() {

        ModelArg modelArg;

        String formulaS = new String(" Surv (time, status) ~ .");
        
        String formulaR = new String("mpg ~ . ");
        String formulaM1 = new String("Multivariate (mpg, wt) ~ hp + drat");
        String formulaM2 = new String("Multivariate (mpg, wt) ~ .");
        String formulaC = new String("Species ~ .");
        
        String formulaU = new String(" Unsupervised () ~ .");

        // Java-side trace.
        Trace.set(15);
        
        if (true) {

            Dataset<Row> irisDF = spark
                .read()
                .option("header", "true")
                .option("inferSchema", "true") 
                .format("csv")
                .load("./test-classes/data/iris.csv");

            // mtcarsDF.printSchema();
            // mtcarsDF.show();

            modelArg = new ModelArg(formulaC, irisDF);


        }

        if (!true) {

            Dataset<Row> mtcarsDF = spark
                .read()
                .option("header", "true")
                .option("inferSchema", "true") 
                .format("csv")
                .load("./test-classes/data/mtcars.csv");

            // mtcarsDF.printSchema();
            // mtcarsDF.show();

            modelArg = new ModelArg(formulaR, mtcarsDF);


        }

        if (!true) {

            Dataset<Row> wihsDF = spark
                .read()
                .option("header", "true")
                .option("inferSchema", "true") 
                .format("csv")
                .load("/Users/kogalur/Dropbox/working/rfsrc/scala/wihs.csv");

            // wihsDF.printSchema();
            // wihsDF.show();

            modelArg = new ModelArg(formulaS, wihsDF);

            if (!true) {

                double[] newTI = new double[(modelArg.get_timeInterest()).length - 20];

                for (int i = 0; i < (modelArg.get_timeInterest()).length - 20; i++) {
                    newTI[i] = modelArg.get_timeInterest()[i+10] + 1.5;
                }

                newTI[(modelArg.get_timeInterest()).length - 21] = 100;
            
                modelArg.set_timeInterest(newTI);

                modelArg.set_eventWeight(new double[] {1, 3});

                modelArg.set_splitRule("random");
                modelArg.set_nSplit(3);

                modelArg.set_nodeSize(10);
                
            }

            
        }

        modelArg.set_ensembleArg("error", "every.tree");

        // Serial or parallel.
        modelArg.set_rfCores(1);

        // Repeatability.
        modelArg.set_seed(-1);

        // We set ntree here.
        modelArg.set_bootstrap(2, "auto", "swr", 0, null);
            
        // Native-code trace.
        // modelArg.set_trace(15 + (1<<13));
       
        RandomForestModel growModel = RandomForest.train(modelArg);
        RandomForestModel restModel = RandomForest.predict(modelArg);
        
        RFLogger.log(Level.WARNING, "\n\nHelloRandomForestSRC() nominal exit.\n\n");                
        
    }

}
