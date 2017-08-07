## #########################################################################
##
##  Script to compress the data files distributed in the package.
##  Used during CRAN package preparation.
##
## #########################################################################

library(tools)

fileName <- dir(path=".", pattern=".txt")
filePrefix <- sapply(fileName, function(x) strsplit(x, ".txt"))


##
## Useful only if package is loaded.
## library("randomSurvivalForest")
## data(list = filePrefix) 
## sapply(filePrefix, function(x) save(list=x, file=paste(x, ".rda", sep="")))
##

dataList <- sapply(fileName, function(x) read.table(x, header = TRUE))
  
for (i in 1:length(dataList)) {
  assign(as.character(filePrefix[i]), dataList[[i]])
  save(list=as.character(filePrefix[i]), file=paste(filePrefix[i], ".rda", sep=""))
}

resaveRdaFiles(".")
