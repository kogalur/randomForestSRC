##########################################################################
## Replica Harness 
##########################################################################

rfslave <- function() {

    ## Note the use of the tag for sent messages: 
    ## 1=ready_for_task, 2=done_task, 3=exiting 

    ## Note the use of the tag for received messages: 
    ## 1=task, 2=done_tasks 

    junk <- 0
    done <- 0

    ## Set the core usage in the slave.
    options(mc.cores = MC_CORES)
    options(rf.cores = RF_CORES)

    ## Continue, as long as the task is not "done". Note the use of
    ## the tag to indicate the type of message.
    while (done != 1) {

        ## Tell the master we are ready to receive a new task.
        mpi.send.Robj(junk, 0, 1)

        ## Receive a new task from the master.
        task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
        task_info <- mpi.get.sourcetag()
        tag <- task_info[2]

        ## Determine whether we received a task from the master.
        if (tag == 1) {

            ## Perform the task, and create results.

            ## Extract the forest ID from the task. 
            forest.id <- task$forest.id

            ## Grow a forest on the slave.
            partial <- rfsrc(outcome ~ .,
                       data = data.sphere[, c(1,2,3,4)],
                       ntree=ntree,
                       importance = importance)

            ## Construct and send a message containing the resulting
            ## forest back to master.
            result <- list(forest.id=forest.id, outcome = partial)
            mpi.send.Robj(result, 0, 2)

        }
        else if (tag == 2) {
            ## The master indicates that all tasks are done. Exit the loop.
            done <- 1
        }
        else {
            ## Unknown message.  Ignore it.
        }
    }
    
    ## Tell master that this slave is exiting.  Send an exiting message.
    mpi.send.Robj(junk, 0, 3)
}


