package com.kogalur.randomforest;

import com.kogalur.RFException;
import com.kogalur.XMLDocument;

import java.lang.NullPointerException;
import org.w3c.dom.DOMException;

import java.io.File;





public class ModelArgsXML {

  private static ModelArgsXML myInstance;

  
  private static final String logDirName      = "output";

  
  private static final String initFileName     = "RandomForestSRC.xml";
  private static final String initBackFileName = "RandomForestSRC.bak";
  private static final String logFileName      = "RandomForestSRC.log";


  
  private static File logFile      = new File(logDirName + File.separator + logFileName);

  private static File saveInitFile = new File(initFile.getParentFile(), initFile.getName() + ".sav");
  private static File tempInitFile = new File(initFile.getParentFile(), initFile.getName() + ".tmp");

  private ModelArgsXML() {
    
  }

  
  public static File getInitFile() {
    return initFile;
  }

  public static File getInitBackFile() {
    return initBackFile;
  }

  public static File getDefaultLogDir() {
    return logFile.getParentFile();
  }
  
  public static File getDefaultLogFile() {
    return logFile;
  }


  public static XMLDocument getInitDoc() throws BamRuntimeException {

    
    
    boolean existence;

    initFile = getInitFile();
    
    try {
      existence = initFile.exists();
    }
    catch (SecurityException se) {
      throw new BamRuntimeException(se, "Required file (" + initFile.getName() + ") is not accessible.  \nPlease ensure that sufficient user privileges are in effect.");
    }
    if(!existence) {
      return null;
    }
    else {
      XMLDocument xmlDocument = new XMLDocument(initFile, true);  
      
      xmlDocument.copy(saveInitFile);  
      return xmlDocument;
    }
  }

  public File getSaveDirectory() throws BamRuntimeException {
    String saveDirName;
    File saveDir;
    boolean existence;

    XMLDocument initDoc = getInitDoc();
    
    if (initDoc == null) {
      throw new BamRuntimeException("ModelArgsXML.getInitDoc() returned null value.", "Internal XML error encountered.  Please contact technical support.");
    }

    try {
      saveDirName = initDoc.getText("saveDirectory");  
    }
    catch (NullPointerException npe) {
      
      saveDirName = logDirName;
    }

    saveDir = new File(saveDirName);
    
    try {
      existence = saveDir.exists();
    }
    catch (SecurityException se) {
      throw new BamRuntimeException(se, "Required directory (" + saveDir.getName() + ") is not accessible.  \nPlease ensure that sufficient user privileges are in effect.");
    }
    if (!existence) {
      
      saveDir = new File("." + File.separator + logDirName);
    }
    return saveDir;
  }


  public static ModelArgsXML getInstance() {
    if (myInstance == null) {
      myInstance = new ModelArgsXML();
    }
    return myInstance;
  }

  public static ModelArgsXML getInstance(File file) {
    if (myInstance == null) {
      myInstance = new ModelArgsXML(file);
    }
    return myInstance;
  }
    
}
