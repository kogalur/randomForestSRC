if (TRUE) {
    ## .txt to .rda
    nut.table <- read.table(file =  "nutrigenomic.txt", sep=" ", header=TRUE)

    ## matrix form is necessary, otherwise the elements are
    ## interpreted as a data frame and list names are prepended to the
    ## columns.
    nutrigenomic <- list(lipids   = as.matrix(nut.table[, 1:21]),
                         genes    = as.matrix(nut.table[, 22:141]),
                         diet     = nut.table[, 142],
                         genotype = nut.table[, 143])

    save(nutrigenomic, file="nutrigenomic.rda")
}

if (FALSE) {
    ## .rda to .txt
    
    ## total of 143 columns
    load("nutrigenomic.rda")
    
    ## lipid is real of 21 columns
    ## genes is real of 120 columns
    ## diet is factor with 5 levels
    ## genotype is factor with 2 levels
    nut.table <- cbind(nutrigenomic$lipids,
                       nutrigenomic$genes,
                       nutrigenomic$diet,
                       nutrigenomic$genotype)
    colnames(nut.table) <- c(colnames(nutrigenomic$lipids),
                             colnames(nutrigenomic$genes),
                             "diet",
                             "genotype")

    write.table(nut.table, file =  "nutrigenomic.txt", row.names=FALSE, sep=" ")
}
