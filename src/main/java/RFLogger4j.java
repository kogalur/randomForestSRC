package com.kogalur.randomforest;

import java.io.IOException;

import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.SimpleFormatter;
import java.util.logging.Level;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import java.io.FileInputStream;

public class RFLogger4j {

    private static org.apache.log4j.LogManager logManager;
    private static org.apache.log4j.Logger logger = null;
    
    private RFLogger4j() {

        // Create the log manager. It seems like Spark automagically
        // reads the .properties file, so we don't really need this
        // downstream.
        logManager = org.apache.log4j.LogManager.getLogManager();

        // Create the logger.
        logger = org.apache.log4j.Logger.getLogger(RFLogger4j.class.getName());

        String parameter = new String("rfsrc test 1 2 3");
        
        if(logger.isDebugEnabled()){
            logger.debug("This is debug : " + parameter);
        }
        
        if(logger.isInfoEnabled()){
            logger.info("This is info : " + parameter);
        }
        
        logger.warn("This is warn : " + parameter);
        logger.error("This is error : " + parameter);
        logger.fatal("This is fatal : " + parameter);


        
    }
    
    private static org.apache.log4j.Logger getLogger() {
        if(logger == null){
            new RFLogger4j();
        }
        return logger;
    }
    
    public static void log(Level level, String msg, Throwable thrown){
        getLogger().log(level, msg, thrown);
    }
    
    public static void log(Level level, String msg){
        getLogger().log(level, msg);
    }

}


