##  **********************************************************************
##  **********************************************************************
##  
##    RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
##  
##    This program is free software; you can redistribute it and/or
##    modify it under the terms of the GNU General Public License
##    as published by the Free Software Foundation; either version 3
##    of the License, or (at your option) any later version.
##  
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##  
##    You should have received a copy of the GNU General Public
##    License along with this program; if not, write to the Free
##    Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
##    Boston, MA  02110-1301, USA.
##  
##    ----------------------------------------------------------------
##    Project Partially Funded By: 
##    ----------------------------------------------------------------
##    Dr. Ishwaran's work was funded in part by DMS grant 1148991 from the
##    National Science Foundation and grant R01 CA163739 from the National
##    Cancer Institute.
##  
##    Dr. Kogalur's work was funded in part by grant R01 CA163739 from the 
##    National Cancer Institute.
##    ----------------------------------------------------------------
##    Written by:
##    ----------------------------------------------------------------
##      Hemant Ishwaran, Ph.D.
##      Director of Statistical Methodology
##      Professor, Division of Biostatistics
##      Clinical Research Building, Room 1058
##      1120 NW 14th Street
##      University of Miami, Miami FL 33136
##  
##      email:  hemant.ishwaran@gmail.com
##      URL:    http://web.ccs.miami.edu/~hishwaran
##      --------------------------------------------------------------
##      Udaya B. Kogalur, Ph.D.
##      Adjunct Staff
##      Department of Quantitative Health Sciences
##      Cleveland Clinic Foundation
##      
##      Kogalur & Company, Inc.
##      5425 Nestleway Drive, Suite L1
##      Clemmons, NC 27012
##  
##      email:  ubk@kogalur.com
##      URL:    https://github.com/kogalur/randomForestSRC
##      --------------------------------------------------------------
##  
##  **********************************************************************
##  **********************************************************************


rf2rfz <- function(object,
                   forestName = NULL,
                   ...)
{
  rfsrcForest <- checkForestObject(object)
  if (is.null(forestName)) {
    stop("RFSRC forest name is NULL.  Please provide a valid name for the forest .rfz file.")
  }
  if (!requireNamespace("XML", quietly = TRUE)) {
    stop("The 'XML' package is required for this function.")
  }
  if (nchar(forestName) > 4) {
    if (substr(forestName, nchar(forestName)-3, nchar(forestName)) == ".rfz") {
      forestName <- substr(forestName, 1, nchar(forestName)-4)
    }
  }
  nativeArray <- rfsrcForest$nativeArray
  time.interest <- rfsrcForest$time.interest
  formula <- rfsrcForest$formula
  forestSeed <- rfsrcForest$seed
  xvar.names <- rfsrcForest$xvar.names
  get.factor <- extract.factor(rfsrcForest$xvar, xvar.names)
  xvar.type <- get.factor$generic.types
  nativeFactorArray <- rfsrcForest$nativeFactorArray
  numTrees <- length(as.vector(unique(nativeArray$treeID)))
  rootString <- getRootString()
  pmmlDoc <- XML::xmlTreeParse(rootString, asText=TRUE)
  pmmlRoot <- XML::xmlRoot(pmmlDoc)
  pmmlRoot <- XML::append.XMLNode(pmmlRoot, getDataDictNode(xvar.names=xvar.names, xvar.type=xvar.type))
  write.table(nativeArray,
              paste(forestName, ".txt", sep=""), quote = FALSE)
  write.table(nativeFactorArray,
              paste(forestName, ".factor.txt", sep=""), col.names=FALSE, quote = FALSE)
  xmlFile <- file(paste(forestName, ".xml", sep=""), open="w")
  XML::saveXML(pmmlRoot, xmlFile)
  close(xmlFile)
  zipCommand <- paste("zip", sep=" ",
                      paste(forestName, ".rfz", sep=""),
                      paste(forestName, ".txt", sep=""),
                      paste(forestName, ".factor.txt", sep=""),
                      paste(forestName, ".xml", sep=""))
  system(command = zipCommand)
  unlink(paste(forestName, ".txt", sep=""))
  unlink(paste(forestName, ".factor.txt", sep=""))
  unlink(paste(forestName, ".xml", sep=""))
}
checkForestObject <- function(object) {
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2) {
    stop("This function only works for objects of class `(rfsrc, grow)' or '(rfsrc, forest)'")
  }
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.null(object$forest)) {
      stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'.")
    }
    rfForest <- object$forest
  }
    else {
      rfForest <- object
    }
  if (is.null(rfForest$nativeArray)) {
    stop("RFsrc nativeArray content is NULL.  Please ensure the object is valid.")
  }
  if (is.null(rfForest$xvar.names)) {
    stop("RFsrc xvar.names content is NULL.  Please ensure the object is valid.")
  }
  if (is.null(rfForest$xvar)) {
    stop("RFsrc xvar content is NULL.  Please ensure the object is valid.")
  }
  return (rfForest)
}
getRootString <- function() {
  rootString <-
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>
         <PMML version=\"3.1\" xmlns=\"http://www.dmg.org/PMML-3_1\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">
           <Header copyright=\"Copyright 2008, Cleveland Clinic\" description=\"Random Survival Forest Tree Model\">
              <Application name=\"Random Survival Forest\" version=\"3.0\"/>
           </Header>
         </PMML>
       "
  return (rootString)
}
getDataDictNode <-  function(xvar.names, xvar.type) {
  dataDictNode <- XML::xmlNode("DataDictionary", attrs=c(numberOfFields=length(xvar.names)))
  for (k in 1:length(xvar.names)) {
    if (xvar.type[k] == "C") {
      dataDictNode <- XML::append.XMLNode(dataDictNode, XML::xmlNode("DataField", attrs=c(name=xvar.names[k], optype="categorical", dataType="string")))
    }
    if (xvar.type[k] == "I") {
      dataDictNode <- XML::append.XMLNode(dataDictNode, XML::xmlNode("DataField", attrs=c(name=xvar.names[k], optype="ordinal", dataType="integer")))
    }
    if (xvar.type[k] == "R") {
      dataDictNode <- XML::append.XMLNode(dataDictNode, XML::xmlNode("DataField", attrs=c(name=xvar.names[k], optype="continuous", dataType="double")))
    }
  }
  return (dataDictNode)
}
