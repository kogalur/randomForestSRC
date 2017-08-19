package com.kogalur.randomforest;

import com.kogalur.RFException;
import com.kogalur.XMLDocument;

import java.lang.NullPointerException;
import org.w3c.dom.DOMException;

import java.io.File;

// This class centralizes all directory and file input/output
// information for the application.  Some values are default static
// values, while others are derived from the XML initialization file.

public class ModelArgsXML {

  private static ModelArgsXML myInstance;

  // For internal use only, static default directory names.
  private static final String logDirName      = "output";

  // For internal use only, static default file names.
  private static final String initFileName     = "RandomForestSRC.xml";
  private static final String initBackFileName = "RandomForestSRC.bak";
  private static final String logFileName      = "RandomForestSRC.log";


  // For external use, static default file objects.
  private static File logFile      = new File(logDirName + File.separator + logFileName);

  private static File saveInitFile = new File(initFile.getParentFile(), initFile.getName() + ".sav");
  private static File tempInitFile = new File(initFile.getParentFile(), initFile.getName() + ".tmp");

  private ModelArgsXML() {
    // Do nothing at present.  This may change in the future.
  }

  // Static methods that do not require instantiation!
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

    // Parse the XML file on each call.  If it does not exist, DO NOT throw an exception, just return null.
    
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
      XMLDocument xmlDocument = new XMLDocument(initFile, true);  // can throw BRE
      // Make a copy of the XML file.
      xmlDocument.copy(saveInitFile);  // can throw BRE
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
      saveDirName = initDoc.getText("saveDirectory");  // can throw NPE if optional tags do not exist
    }
    catch (NullPointerException npe) {
      // Optional tag does not exist.  Initialize the default value.
      saveDirName = logDirName;
    }

    saveDir = new File(saveDirName);
    // Validate the directory and determine if it is unitialized or invalid.
    try {
      existence = saveDir.exists();
    }
    catch (SecurityException se) {
      throw new BamRuntimeException(se, "Required directory (" + saveDir.getName() + ") is not accessible.  \nPlease ensure that sufficient user privileges are in effect.");
    }
    if (!existence) {
      // Fall back to the default directory.
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
