## #########################################################################
##
##  Script to compress the data files distributed in the package.
##  Used during CRAN package preparation.
##
## #########################################################################

library(tools)

## pull txt file names
fileName.txt <- dir(path=".", pattern=".txt")
filePrefix.txt <- sapply(fileName.txt, function(x) strsplit(x, ".txt"))

## pull csv file names
fileName.csv <- dir(path=".", pattern=".csv")
filePrefix.csv <- sapply(fileName.csv, function(x) strsplit(x, ".csv"))

print(fileName.csv)

## first we do txt files
dataList.txt <- lapply(fileName.txt, function(x) read.table(x, header = TRUE, stringsAsFactors = TRUE))
 
for (i in 1:length(dataList.txt)) {
  assign(as.character(filePrefix.txt[i]), dataList.txt[[i]])
  save(list=as.character(filePrefix.txt[i]), file=paste(filePrefix.txt[i], ".rda", sep=""))
}

### now we do csv files
dataList.csv <- lapply(fileName.csv, function(x) read.csv(x, header = TRUE, stringsAsFactors = TRUE))
 
for (i in 1:length(dataList.csv)) {
  assign(as.character(filePrefix.csv[i]), dataList.csv[[i]])
  save(list=as.character(filePrefix.csv[i]), file=paste(filePrefix.csv[i], ".rda", sep=""))
}

resaveRdaFiles(".")
