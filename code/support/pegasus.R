##########################################################################
##########################################################################
##
## VERSION 9
##
##########################################################################
##########################################################################


slave.count <- c(32:1, rep(1,19))
core.count  <- c(rep(-1, 32)), 19:1)

slave.count <- c( 2,  1) 
core.count  <- c(-1, -1)

cores.per.slave <- 20
ntree.all <- 2048


##########################################################################
##########################################################################

source("simulation.R")

if(!is.loaded("mpi_initialize")){
  require("Rmpi")
}

library(randomForestSRC)

rfslave <- function() {
    ## Note the use of the tag for sent messages: 
    ## 1=ready_for_task, 2=done_task, 3=exiting 
    ## Note the use of the tag for received messages: 
    ## 1=task, 2=done_tasks 

    junk <- 0
    done <- 0

    options(rf.cores = MY_RF_CORES)
    options(mc.cores = MY_MC_CORES)

    ## While task is not a "done" message. Note the use of the
    ## tag to indicate the type of message.
    while (done != 1) {

        ## Signal being ready to receive a new task.
        mpi.send.Robj(junk, 0, 1)

        ## Receive a task.
        task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
        task_info <- mpi.get.sourcetag()
        tag <- task_info[2]

        if (tag == 1) {

            ## Perform the task, and create results.
            forest.id <- task$forest.id

            ## Call the randomForestSRC test case.
            ## x <- rfsrc(outcome ~ .,
            x <- rfsrc(Surv(time, event) ~ .,                       
                       data = x.summary,
                       ntree=ntree[forest.id],
                       importance = importance)

            ## Construct and send message back to master
            result <- list(forest.id=forest.id, outcome = x)
            mpi.send.Robj(result, 0, 2)

        }
        else if (tag == 2) {
            ## Master is saying that all tasks are done, so exit the loop.
            done <- 1
        }
        else {
            ## Ignore the message.
        }
    }
    
    ## Tell master that this slave is exiting.  Send master an exiting message
    mpi.send.Robj(junk, 0, 3)
}


##########################################################################
##########################################################################

loop.count <- length(slave.count)

timer = vector('list', loop.count)

importance = "permute"

## Repeatability w.r.t. the data set (only).  The forests are still random.
set.seed(-1)

x.summary <- get.data()

for (j in 1:loop.count) {

    ntree = rep(floor(ntree.all / slave.count[j]), slave.count[j])

    if (slave.count[j] > 1) {
        for (i in 1:(ntree.all %% slave.count[j])) {
            ntree[i] = ntree[i] + 1
        }
    }
    
    MY_MC_CORES =   1
    MY_RF_CORES =   core.count[j]
    
    ## Spawn the number of slaves we have pre-specified.
    mpi.spawn.Rslaves(nslaves = slave.count[j], needlog = TRUE)

    ## In case R exits unexpectedly, have it automatically clean up
    ## resources taken up by Rmpi (slaves, memory, etc...)
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

    junk <- 0

    ## Number of slaves actually spawned.  This should be the number
    ## requested, but we do not do any checks.
    n_slaves <- mpi.comm.size() - 1

    ## Create task list, with the task count equal to the slave count.
    n_tasks <- n_slaves

    ## Keep track of the closed slaves.
    closed_slaves <- 0

    tasks <- vector('list', n_tasks)
    for (i in 1:n_tasks) {
        tasks[[i]] <- list(forest.id=i)
    }

    ## Vector for results.
    result <- vector('list', n_tasks)

    ## Send the function to the slaves
    mpi.bcast.Robj2slave(rfslave)

    ## Send core usage information to the slaves.
    mpi.bcast.Robj2slave(MY_RF_CORES)
    mpi.bcast.Robj2slave(MY_MC_CORES)

    mpi.bcast.Robj2slave(x.summary)
    mpi.bcast.Robj2slave(ntree)
    mpi.bcast.Robj2slave(importance)

    mpi.bcast.cmd(library(randomForestSRC))
    
    ## Call the function in all the slaves to get them ready to undertake tasks
    mpi.bcast.cmd(rfslave())

    start <- proc.time()

    while (closed_slaves < n_slaves) {

        ## Receive a message from a slave.
        message <- mpi.recv.Robj(mpi.any.source(), mpi.any.tag())
        message_info <- mpi.get.sourcetag()
        slave_id <- message_info[1]
        tag      <- message_info[2]

        if (tag == 1) {

            ## Slave is ready for a task.  Give it the next task, or tell it tasks
            ## are done if there are none.
            if (length(tasks) > 0) {
                ## Send a task, and then remove it from the task list.
                mpi.send.Robj(tasks[[1]], slave_id, 1);
                tasks[[1]] <- NULL
            }
            else {
                mpi.send.Robj(junk, slave_id, 2)
            }
        }
        else if (tag == 2) {

            ## The message contains results.
            forest.id <- message$forest.id
            mresult <- list(forest.id = message$forest.id, outcome = message$outcome)
            result[[forest.id]] <- mresult
        }
        else if (tag == 3) {

            ## A slave has closed down. 
            closed_slaves <- closed_slaves + 1
        }
    }

    end <- proc.time()
    print(result)
    delta <- end - start
    print(delta)

    timer[[j]] <- list(delta = delta, slave.count = slave.count[j], core.count = core.count[j]) 

    result <- matrix(nrow=length(timer[1:j]), ncol=3)
    result[, 1] = round(sapply(timer[1:j], function(x) x$delta[3]), 2)
    result[, 2] = sapply(timer[1:j], function(x) x$slave.count)
    result[, 3] = sapply(timer[1:j], function(x) x$core.count)
    colnames(result) = c("elapsed", "slaves", "cores")
    result[which(result[, 3] == -1), 3] = cores.per.slave
    write.csv(result, file="result.txt", row.names=F)

    
    mpi.bcast.cmd(detach("package:randomForestSRC", unload=TRUE))


    ## Close down.
    mpi.close.Rslaves(dellog = FALSE)

}

detach("package:randomForestSRC", unload=TRUE)

print(timer)


if(TRUE) {

    result = read.csv(file="result.txt")

    result = result[dim(result)[1]:1, ]
    
    png(filename = "result.png",
    width = 1200, height = 1200,
    units = "px",
    pointsize = 42,
    res=72)

    result <- result
    x = result[,2] * result[,3]
    y = round(result[,1])

    plot(x=log2(x), y=log2(y),
         xlab="Total Cores", ylab="Elapsed Time (minutes)",
         xlim=c(log2(min(x)), log2(max(x))), ylim=c(log2(min(y)), log2(max(y))), axes=F,
         type="n")
    points(log2(x), log2(y), pch=18,cex=1.5,
           col=c(rep("red", length(which(result[,2] == 1))), rep("blue", dim(result)[1] - length(which(result[,2] == 1)))),)
    axis(1, at=log2(x) , labels=x)
    y.axis.floor = min(floor(log2(y/60)))
    y.axis.ceil = max(ceiling(log2(y/60)))
    axis(2, at=(y.axis.floor:y.axis.ceil) + log2(60), labels=2^(y.axis.floor:y.axis.ceil))
    legend('topright', c("OpenMP", "OpenMP & Open MPI"), col=c('red','blue'), pch=18)

    dta <- data.frame(x=log2(x), y=log2(y))
    lm.object <- lm(y~x, dta)
    yhat <- lm.object$fitted.values
    lines(log2(x), yhat, col="black", lwd=2)
    
    dev.off()

}

mpi.quit()

