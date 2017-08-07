package com.kogalur.randomforest;

import java.io.IOException;

import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.SimpleFormatter;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.LogManager;

import java.io.FileInputStream;

public class RFLogger {

    private String fileName = "./classes/log.properties";
    private static java.util.logging.Logger logger = null;
    
    private static ConsoleHandler consoleHandler;
    
    private static SimpleFormatter simpleFormatter;

    private RFLogger() {

        
        logger = java.util.logging.Logger.getLogger(RFLogger.class.getName());

        try {
            java.util.logging.LogManager.getLogManager().readConfiguration(new FileInputStream(fileName));

        }
        catch (IOException e) {
            System.out.println("Unable to configure logging for " + fileName);

            e.printStackTrace();

            System.out.println("Resorting to console handler.");

        
            
            consoleHandler = new ConsoleHandler();
        
            
            logger.addHandler(consoleHandler);
        
            
            consoleHandler.setLevel(Level.INFO);
            logger.setLevel(Level.ALL);
        
            
            simpleFormatter = new SimpleFormatter();

            
            consoleHandler.setFormatter(simpleFormatter);

            
            
        }
    }
    
    private static java.util.logging.Logger getLogger() {
        if(logger == null){
            new RFLogger();
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


