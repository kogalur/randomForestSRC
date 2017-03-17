library("randomForestSRC")


#########################################################################
## Custom functions
#########################################################################

get.hyper.rfsrc <- function(object, treeID, low, high, ...)
{
  if (is.null(object)) stop("Object is empty!")
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2)
      stop("This function only works for objects of class `(rfsrc, grow)' or '(rfsrc, forest)'")
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
      if (is.null(object$forest)) 
          stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'.")
      object <- object$forest
  }
  nativeArray <- object$nativeArray
  if (is.null(nativeArray)) {
      stop("RFSRC nativeArray content is NULL.  Please ensure the object is valid.")
  }
  xvar.names <- object$xvar.names
  if (is.null(xvar.names)) {
      stop("RFSRC xvar.names content is NULL.  Please ensure the object is valid.")
  }
  if (is.null(object$xvar)) {
      stop("RFSRC xvar content is NULL.  Please ensure the object is valid.")
  }

  nativeArray <- nativeArray[which(nativeArray$treeID == treeID), ]
  max.depth = max(nativeArray$nodeID)

  x.count <- length(xvar.names)

  boundary.local <- matrix (c(rep(low, x.count), rep(high, x.count)), nrow=x.count, ncol=2)

  ## Boundaries are on real values only, so all segments are of the form [a, b].
  ## Initial decision boundaries span the full space from [low, high] in all dimensions.
  boundary <- array(0, dim = c(max.depth, x.count, 2))
  boundary[, , 1] = low
  boundary[, , 2] =  high

  offset <- 1
  lt.flag <- TRUE
  recursive.obj <- list(boundary = boundary, offset = offset)
  
  recursive.obj <- parse.tree(
      recursive.obj,
      nativeArray,
      boundary.local)

  return (recursive.obj$boundary) 
}

parse.tree <- function(recursive.obj,
                       nativeArray,
                       boundary.local)
{
    offset <- recursive.obj$offset
    parmID <- nativeArray$parmID[offset]
    
    if (parmID == 0) {
        recursive.obj$boundary[nativeArray$nodeID[offset], , ] = boundary.local
    }

    recursive.obj$offset <- recursive.obj$offset + 1
    
    if (parmID != 0) {

        restore.point = boundary.local[parmID, 2]
        boundary.local[parmID, 2] <- nativeArray$contPT[offset]        
        recursive.obj <- parse.tree(recursive.obj, nativeArray, boundary.local)

        boundary.local[parmID, 2] <- restore.point
        boundary.local[parmID, 1] <- nativeArray$contPT[offset]        

        recursive.obj <- parse.tree(recursive.obj, nativeArray, boundary.local)

    }
    return(recursive.obj)
}



#########################################################################
## Four circles with custom class sizes and radii.
#########################################################################
source("spheres.R")

#########################################################################
## Start of script
#########################################################################
options(rf.cores = 1)

set.seed(-1)
seed <-  -1

importance = "none"

## Grow the forest.  This forest is used in the rest of the plots in
## the paper.

ntree <- 1000
importance <- "none"
seed = 3 * seed

## The data set consists of four circles, equal class counts, touching at the origin.
circles <- get.spheres()[1,2, 4]

grow.result <- rfsrc(outcome ~ .,
                     data = circles,
                     ntree = ntree,
                     importance = importance,
                     seed = seed)



#########################################################################
##  Single decision boundary:
##    terminal.bar.1-6.png
##    decision.single.png
##  Result:
##    Used in Fig 2:  class_decision_boundary.png
#########################################################################

## Single tree decision boundary.  CAUTION:  Do not change parameters!
## It succeeds on the FIRST tree.
## It relies on randomForestSRC 1.6.1.1 SERIAL VERSION.
## Note, it does NOT use the parallel version.
## Note, it does NOT use SRCM (lol).
## Do not change the custom seed in this call!

nodesize <- 2000
nodedepth <- 3
ntree <- 1
ntree.idx <- 1

single.result <- rfsrc(outcome ~ .,
                data = circles,
                ntree = ntree,
                splitrule = "random",
                nodesize = nodesize,
                nodedepth=nodedepth,
                importance = importance,
                seed = -11)

print(single.result$forest$nativeArray)

for (i in 1 : max(single.result$membership[, ntree.idx])) {

    idx <- which (single.result$membership[, ntree.idx] == i)
    counts <- table(circles[idx, x.count + 1]) / length(idx)

    png(filename = paste("images/terminal.bar", as.character(i), ".png", sep=""),
        width = 400, height = 600,
        units = "px",
        pointsize = 36,
        res=72)
     
    par(mar=c(1.1,3.1,1.1,1.1))    

    barplot(counts,
            names.arg=rep("", class.count),
            space = 0,
            main=NULL, col = c('red', 'blue', 'green', 'yellow'),
            xlab=NULL, xlim=c(0,class.count),
            ylab=NULL, ylim=c(0.0, 1.0),
            width=1, 
            cex.axis=1.25)
    
    dev.off()
}

## Call to get decision boundaries.
boundary <- get.hyper.rfsrc(single.result$forest, 1, -2, 2)

png(filename = "images/decision.single.png",
    width = 1400, height = 1200,
    units = "px",
    pointsize = 36,
    res=72)

par(mar=c(4,4,1,0)+0.1)

plot(circles[, 1:x.count-1],
     type="p",
     cex=0.25,
     main=NULL,
     xlab="", ylab="",     
     col = c('red', 'blue', 'green', 'yellow')[circles[, x.count + 1]],
     pch=19,
     xlim=c(-2, 2), ylim=c(-2, 2),
     xaxs = "i", yaxs = "i",
     cex.lab=3,
     asp = 1,
     axes = F)

axis(1, pos=-2, lwd=3, cex.axis=1.5)
axis(2, pos=-2, lwd=3, cex.axis=1.5)

mtext(expression('x' [1]), side=1,line=3, cex=2)
mtext(expression('x' [2]), side=2,line=1, cex=2)

for(i in 1:dim(boundary)[1]) {
    ## (x1,y1) -> (x2,y1)
    segments(boundary[i, 1, 1], boundary[i, 2, 1], boundary[i, 1, 2], boundary[i, 2, 1], lwd=3)
    ## (x2,y1) -> (x2,y2)
    segments(boundary[i, 1, 2], boundary[i, 2, 1], boundary[i, 1, 2], boundary[i, 2, 2], lwd=3)
    ## (x1,y2) -> (x2,y2)
    segments(boundary[i, 1, 1], boundary[i, 2, 2], boundary[i, 1, 2], boundary[i, 2, 2], lwd=3)
    ## (x1,y1) -> (x1,y2)
    segments(boundary[i, 1, 1], boundary[i, 2, 1], boundary[i, 1, 1], boundary[i, 2, 2], lwd=3)

    ## Get the mean of the boundary for positioning of the label.
    x.mean <- mean(c(boundary[i, 1, 1], boundary[i, 1, 2]))
    y.mean <- mean(c(boundary[i, 2, 1], boundary[i, 2, 2]))
    points(x=x.mean, y=y.mean, pch = as.character(i), cex=2)
    
}
               
dev.off()







#########################################################################
##  Forest decision boundary:
##    origin.terminal.bar.1-8.png
##    origin.ensemble.bar.png
##    decision.forest.png
##    anon.forest.png
##
##  Result:
##    Used in Fig 3:  class_decision_forest.png
#########################################################################



## Predict at the origin.
pred.result.indv <- predict.rfsrc(grow.result,
                                  newdata = circles[1, ],
                                  importance = importance,
                                  seed = seed)



## Get eight bar charts representing tree prediction at the origin.
## We choose the best trees and only use four of these (lol).

sig.count <- 1
i <- 1
while ((sig.count <= 100) & (i <= ntree)) {

    ## What node is the individual in for this tree.
    nodeID <- pred.result.indv$membership[i]

    ## Get the local distribution for this node in this tree, and plot it.
    idx <- which (grow.result$membership[, i] == nodeID)
    counts <- table(circles[idx, x.count + 1]) / length(idx)

    ## Require that the classes not be fucking boring. Hahaha.
    if ((sum(counts == 0.25) != 4) &
        (sum(counts > 0.05) == 4) &
        (sum(counts == 0.20) < 3) &
        (sum(counts == counts[2]) < 3)) {


        png(filename = paste("images/origin.terminal.bar", as.character(i), ".png", sep=""),
            width = 400, height = 600,
            units = "px",
            pointsize = 36,
            res=72)

        par(mar=c(1.1,3.1,1.1,1.1))
        
        barplot(counts,
                names.arg=rep("", class.count),
                space = 0,
                main=NULL, col = c('red','blue','green', 'yellow'),
                xlab=NULL, xlim=c(0,class.count),
                ylab=NULL, ylim=c(0.0, 0.60),
                width=1,
                cex.axis=1.25)

        dev.off()

        sig.count <- sig.count + 1
    }

    i <- i + 1
}



## Plot the bar chart for the forest ensemble.
png(filename = paste("images/origin.ensemble.bar", ".png", sep=""),
    width = 400, height = 600,
    units = "px",
    pointsize = 36,
    res=72)

par(mar=c(1.1,3.1,1.1,1.1))

barplot(pred.result.indv$predicted[1,],
        names.arg=rep("", class.count),
        space = 0,
        main=NULL, col = c('red','blue','green', 'yellow'),
        xlab=NULL, xlim=c(0,class.count),
        ylab=NULL, ylim=c(0.0, 0.60),
        width=1, 
        cex.axis=1.25)

dev.off()


## Create a mesh of points, and plot the decision boundary.
t.grid = expand.grid(seq(-2, 2, 0.02), seq(-2, 2, 0.02))

jit.factor = 2
t.grid[, 1] = jitter(t.grid[, 1], factor=jit.factor)
t.grid[, 2] = jitter(t.grid[, 2], factor=jit.factor)

t.grid = data.frame(t.grid)
names(t.grid) = x.names[1:2]

pred.result <- predict.rfsrc(grow.result, newdata = t.grid, seed=seed)

## Plot the decision boundray based on the training data
png(filename = "images/decision.forest.png",
    width = 1400, height = 1200,
    units = "px",
    pointsize = 36,
    res=72)

par(mar=c(4,4,1,0)+0.1)

plot(t.grid,
     type="p",
     cex=0.25,
     main=NULL,
     xlab="", ylab="",
     col = c('red', 'blue', 'green', 'yellow')[pred.result$class],
     pch=19,
     xlim=c(-2, 2), ylim=c(-2,2),
     xaxs = "i", yaxs = "i",
     cex.lab=3,
     asp = 1,
     axes = F)

axis(1, pos=-2, lwd=3, cex.axis=1.5)
axis(2, pos=-2, lwd=3, cex.axis=1.5)

mtext(expression('x' [1]), side=1,line=3, cex=2)
mtext(expression('x' [2]), side=2,line=1, cex=2)

segments(-2, 0, 2, 0, lwd=3)
segments(0, -2, 0, 2, lwd=3)

dev.off()


## Plot the anonymized test data.
png(filename = "images/anon.forest.png",
    width = 1400, height = 1200,
    units = "px",
    pointsize = 36,
    res=72)

par(mar=c(4,4,1,0)+0.1)

plot(t.grid,
     type="p",
     cex=0.20,
     main=NULL,
     xlab="", ylab="",
     col = 'black',
     pch=19,
     xlim=c(-2, 2), ylim=c(-2,2),
     xaxs = "i", yaxs = "i",
     cex.lab=3,
     asp = 1,
     axes = F)

axis(1, pos=-2, lwd=3, cex.axis=1.5)
axis(2, pos=-2, lwd=3, cex.axis=1.5)

mtext(expression('x' [1]), side=1,line=3, cex=2)
mtext(expression('x' [2]), side=2,line=1, cex=2)

segments(-2, 0, 2, 0, lwd=3)
segments(0, -2, 0, 2, lwd=3)

dev.off()








#########################################################################
##  Multivariate Mixed Outcome:
##    Use the code example in the paper to produce: 
##    mtcars1.png, mtcars2.png
##    Dimensions are 512 x 768 and manually exported to png.   
##
##  Result:
##    Used in Fig 4:  class_decision_forest.png
#########################################################################

# Multivariate mixed outcome: Motor Trend Car Road Tests.  Miles per
# gallon is the usual response, but the data set is modified to
# introduce categorical and ordinal variables.
mtcars.mod <- mtcars
mtcars.mod$cyl <- factor(mtcars.mod$cyl)
mtcars.mod$carb <- factor(mtcars.mod$carb, ordered=TRUE)
mtcars.mix <- rfsrc(cbind(carb, mpg, cyl) ~., data = mtcars.mod)

plot.variable(mtcars.mix, outcome.target = 3, which.outcome = 1, partial = TRUE, nvar = 1)
plot.variable(mtcars.mix, outcome.target = 3, which.outcome = 2, partial = TRUE, nvar = 1)



#########################################################################
##  Outcome = "test":
##    train.png
##    origin.oob.bar.png
##    decision.forest.png
##    anon.forest.png

##
##  Result:
##    Used in Fig 6:  outcome_test.png
#########################################################################




#########################################################################
## Plot the training data set. 
#########################################################################
png(filename = "images/train.png",
    width = 1400, height = 1200,
    units = "px",
    pointsize = 36,
    res=72)

par(mar=c(4,4,1,0)+0.1)
    
plot(circles[, 1:x.count-1],
     type="p",
     cex=0.25,
     main=NULL,
     xlab="", ylab="",     
     col = c('red','blue','green', 'yellow')[circles[, x.count + 1]],
     pch=19,
     xlim=c(-2, 2), ylim=c(-2,2),
     xaxs = "i", yaxs = "i",
     cex.lab=3,
     asp = 1,
     axes = F)

axis(1, pos=-2, lwd=3, cex.axis=1.5)
axis(2, pos=-2, lwd=3, cex.axis=1.5)

mtext(expression('x' [1]), side=1,line=3, cex=2)
mtext(expression('x' [2]), side=2,line=1, cex=2)

segments(-2, 0, 2, 0, lwd=3)
segments(0, -2, 0, 2, lwd=3)
segments(2, -2, 2, 2, lwd=1)
segments(-2, 2, 2, 2, lwd=1)

dev.off()



## Not used.
for (i in 1:class.count) {
    png(filename = paste("images/origin.oob.bar", as.character(i), ".png", sep=""),
        width = 400, height = 600,
        units = "px",
        pointsize = 36,
        res=72)

    par(mar=c(1.1,3.1,1.1,1.1))
    
    barplot(grow.result$predicted.oob[(((i - 1) * class.size[i]) + 1), ],
            names.arg=rep("", class.count),
            space = 0,
            main=NULL, col = c('red', 'blue', 'green', 'yellow'),
            xlab=NULL, xlim=c(0, class.count), 
            ylab=NULL, ylim=c(0.0, 0.4),
            width=1, 
            cex.axis=1.25
            )

    dev.off()
}



## Get outcome="train" bar ensemble for just the origin.
for (i in 1:1) {
    pred.result <- predict.rfsrc(grow.result,
                                 newdata = circles[(((i - 1) * class.size[i]) + 1), ],
                                 outcome = "train",
                                 importance = importance,
                                 seed = seed)

    png(filename = paste("images/origin.train.bar", as.character(i), ".png", sep=""),
        width = 400, height = 600,
        units = "px",
        pointsize = 36,
        res=72)

    par(mar=c(1.1,3.1,1.1,1.1))    

    barplot(pred.result$predicted[1, ],
            names.arg=rep("", class.count),
            space = 0,
            main=NULL, col = c('red', 'blue', 'green', 'yellow'),
            xlab=NULL, xlim=c(0,class.count),
            ylab=NULL, ylim=c(0.0, 0.6),
            width=1, 
            cex.axis=1.25
            )
    
    dev.off()

}


## Get outcome="test" result by removing two classes at the origin.
## This is accomplished by making the radii unequal.  Note that we
## suppress the X3 dimension.
test.circles <- get.spheres(class.radius=c(1.0, 0.75, 1.0, 0.75))[1, 2, 4]

png(filename = "images/test.png",
    width = 1400, height = 1200,
    units = "px",
    pointsize = 36,
    res=72)

par(mar=c(4,4,1,0)+0.1)

plot(test.circles[, 1:x.count-1],
     type="p",
     cex=0.25,
     main=NULL,
     xlab="", ylab="",
     col = c('red','blue','green', 'yellow')[test.circles[, x.count + 1]],
     pch=19,
     xlim=c(-2, 2), ylim=c(-2,2),
     xaxs = "i", yaxs = "i",
     cex.lab=3,
     asp = 1,
     axes = F)

axis(1, pos=-2, lwd=3, cex.axis=1.5)
axis(2, pos=-2, lwd=3, cex.axis=1.5)

mtext(expression('x' [1]), side=1,line=3, cex=2)
mtext(expression('x' [2]), side=2,line=1, cex=2)

segments(-2, 0, 2, 0, lwd=3)
segments(0, -2, 0, 2, lwd=3)
segments(2, -2, 2, 2, lwd=1)
segments(-2, 2, 2, 2, lwd=1)

dev.off()



## Push the test data set through.  All of it.  
pred.result <- predict.rfsrc(grow.result,
                             newdata = test.circles,
                             outcome = "test",
                             importance = importance,
                             seed = seed)


for (i in 1:1) {


    png(filename = paste("images/outcome.test.bar", as.character(i), ".png", sep=""),
        width = 400, height = 600,
        units = "px",
        pointsize = 36,
        res=72)

    par(mar=c(1.1,3.1,1.1,1.1))
    
    barplot(pred.result$predicted[(((i - 1) * class.size[i]) + 1), ],
            names.arg=rep("", class.count),
            space = 0,
            main=NULL, col = c('red', 'blue', 'green', 'yellow'),
            xlab=NULL, xlim=c(0,class.count),
            ylab=NULL, ylim=c(0.0, 0.6),
            width=1,
            cex.axis=1.25
            )
    
    dev.off()
}






