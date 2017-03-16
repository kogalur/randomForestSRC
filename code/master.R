##########################################################################
##  Primary Harness
##########################################################################

## Load the R package that implements MPI parallel processing.
if(!is.loaded("mpi_initialize")){
  require("Rmpi")
}

## Load our package on the master.
library(randomForestSRC)

## Set the slave and core usage.
RF_SLAVES =    4  ## Number of slaves.
MC_CORES  =    1  ## Turn off mclapply() in package:parallel.
RF_CORES  =   -1  ## Use all cores (and threads) on each slave.

## Create the simulated data set.
data.sphere <- get.data()

## Spread the forest evenly among the slaves.
ntree <- round(1000 / RF_SLAVES)

## Calculate importance.
importance <- "permute"

## Spawn the number of slaves we have specified.
mpi.spawn.Rslaves(nslaves = RF_SLAVES, needlog = TRUE)    

## Cleanup resources taken up by Rmpi.
.Last <- function(){
    if (is.loaded("mpi_initialize")){
        if (mpi.comm.size(1) > 0){
            print("Please use mpi.close.Rslaves() to close slaves.")
            mpi.close.Rslaves(dellog = FALSE)
        }
        print("Please use mpi.quit() to quit R")
        .Call("mpi_finalize")
    }
}

## Meaningless message content.
junk <- 0

## Note the number of slaves actually spawned.  This should be the
## number requested, but we do not do any checks.
n_slaves <- mpi.comm.size() - 1

## Set the number of task equal to the number of slaves.  In general,
## the number of tasks can be greater ## than the number of slaves.
## In such a scenario, tasks are run in a simple round-robin fashion
## on the available slaves.
n_tasks <- n_slaves

## Keep track of the number of closed slaves.
closed_slaves <- 0

## Create the task list. Note that we label them with their forest
## identifier.
tasks <- vector('list', n_tasks)
for (i in 1:n_tasks) {
    tasks[[i]] <- list(forest.id=i)
}

## Crete a list to hold the resulting forests.
result <- vector('list', n_tasks)

## Send the slave function to the slaves.
mpi.bcast.Robj2slave(rfslave)

## Send core usage information to the slaves.
mpi.bcast.Robj2slave(MC_CORES)
mpi.bcast.Robj2slave(RF_CORES)

## Send the data, and relevant variables to the slaves.
mpi.bcast.Robj2slave(data.sphere)
mpi.bcast.Robj2slave(ntree)
mpi.bcast.Robj2slave(importance)

## Load our package on the slaves.
mpi.bcast.cmd(library(randomForestSRC))

## Call the function in all the slaves. The function will idle on a
## slave until the master sends it a task.
mpi.bcast.cmd(rfslave())

## Continue processing tasks until all slaves have closed down.
while (closed_slaves < n_slaves) {

    ## Receive a message from a slave.
    message <- mpi.recv.Robj(mpi.any.source(), mpi.any.tag())
    message_info <- mpi.get.sourcetag()
    slave_id <- message_info[1]
    tag      <- message_info[2]

    if (tag == 1) {
        ## A slave is ready for a task.  Give it the next task, or
        ## tell it that all tasks are done if there are no more.
        if (length(tasks) > 0) {
            ## Send a task, and then remove it from the task list.
            mpi.send.Robj(tasks[[1]], slave_id, 1);
            tasks[[1]] <- NULL
        }
        else {
            ## Tell the slave that there are no more tasks.
            mpi.send.Robj(junk, slave_id, 2)
        }
    }
    else if (tag == 2) {
        ## The message from the slave contains results.  Extract the forest
        ## identifier, the forest object, and add it to the list object.
        forest.id <- message$forest.id
        mresult <- list(forest.id = message$forest.id, outcome = message$outcome)
        result[[forest.id]] <- mresult
    }
    else if (tag == 3) {
        ## A slave has closed down. Make note of it. 
        closed_slaves <- closed_slaves + 1
    }
}

## Combine the out-of-bag predicted results.
combine.predicted.oob = 0
for (i in 1:n_tasks) {
    combine.predicted.oob <- result[[i]]$outcome$predicted.oob + combine.predicted.oob
}
combine.predicted.oob = combine.predicted.oob / n_tasks
combine.class.oob <- apply(combine.predicted.oob,1, function(x) which.max(x))

combine.importance = 0
for (i in 1:n_tasks) {
    combine.importance <- result[[i]]$outcome$importance + combine.importance
}
combine.importance = combine.importance / n_tasks

## Detach our package on the slaves.
mpi.bcast.cmd(detach("package:randomForestSRC", unload=TRUE))
                
## Close down the slaves.
mpi.close.Rslaves(dellog = FALSE)

## Detach our package on the master.
detach("package:randomForestSRC", unload=TRUE)

## Terminate MPI execution and exit.
mpi.quit()
