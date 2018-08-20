import com.kogalur.randomforest.ModelArg;
import com.kogalur.randomforest.RandomForest;
import com.kogalur.randomforest.RandomForestModel;
import com.kogalur.randomforest.Trace;

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

        boolean regr, clas, mult1, mult2, unsp, surv;

        // Flags for diagnostic printing.
        regr = clas = mult1 = mult2 = unsp = surv = false;

        modelArg = null;

        // Java-side trace.
        Trace.set(15);

        // Set the family and data set here!
        regr = true;
        
        if (clas) {

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

        if (regr) {

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

        if (surv) {

            Dataset<Row> wihsDF = spark
                .read()
                .option("header", "true")
                .option("inferSchema", "true") 
                .format("csv")
                .load("/Users/kogalur/Dropbox/working/rfsrc/scala/wihs.csv");

            // wihsDF.printSchema();
            // wihsDF.show();

            modelArg = new ModelArg(formulaS, wihsDF);


            // Play with some custom options for survival.
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

        // TBD TBD TBD we still need to trim the "quantitativeTerminalInfo" TBD TBD TBD

        // Serial or parallel. Set it to serial so we can actually view coherent trace!
        modelArg.set_rfCores(1);

        // Repeatability.
        modelArg.set_seed(-1);

        // We set ntree here.
        modelArg.set_bootstrap(4, "auto", "swr", 0, null, null);

        // Set blockSize explicitly.
        modelArg.set_blockSize(1);
        
        // Native-code trace.
        modelArg.set_trace(15);
        modelArg.set_trace(15 + (1<<4) + (1<<5) + (1<<6) + (1<<7) + (1<<8) + (1<<16));
        RFLogger.log(Level.INFO, "Native Code Trace: " + String.format("0x%32X", modelArg.get_trace()));
        

        RandomForestModel growModel = RandomForest.train(modelArg);

        //  int[][] rmbrMembership = (int[][]) growModel.getEnsemble("rmbrMembership");
        //  growModel.printEnsemble(rmbrMembership);

        if (regr) {

            double[][] perfRegr = (double[][]) growModel.getEnsemble("perfRegr");
            growModel.printEnsemble(perfRegr);

        }
        if (clas) {
            double[][][] perfClas = (double[][][]) growModel.getEnsemble("perfClas");
            growModel.printEnsemble(perfClas);
        }
        
        RFLogger.log(Level.INFO, "\n\nHelloRandomForestSRC() GROW nominal exit.\n\n");

        // RESTORE with different options.
        if (!false) {

            growModel.set_trace(15 + (1<<4) + (1<<8) + (1<<16));

            if (growModel.getModelArg().get_htry() == 0) {
                // growModel.setEnsembleArg("importance", "permute");
                growModel.setEnsembleArg("importance", "no");
                growModel.set_xImportance();
            }
            else {
                growModel.setEnsembleArg("importance", "no");
                growModel.set_xImportance();
            }
            

            growModel.setEnsembleArg("proximity", "oob");

            RandomForestModel restModel = RandomForest.predict(growModel);

            if (regr) {
                double[][] perfRegr = (double[][]) restModel.getEnsemble("perfRegr");
                restModel.printEnsemble(perfRegr);
            }

            
            RFLogger.log(Level.INFO, "\n\nHelloRandomForestSRC() REST nominal exit.\n\n");
        }
        
        System.out.println("\n\nHelloRandomForestSRC() nominal exit.\n\n");
    }

}
