##########################################################################
## Simulated Data consists of four spheres (classes).  In the figures
## we only use the projection on the x- and y- plane.  In the parallel
## processing benchmarks, we used the full spheres.
##########################################################################

get.spheres <- function(class.size=c(5000, 5000, 5000, 5000),
                        class.radius=c(1.0, 1.0, 1.0, 1.0)) {
    
    ## Four spheres in 3D.
    class.count <- 4;
    x.count <- 3

    ## We place the spheres centered on the z-axis only.
    class.center  = matrix(c(-1,  0, 0,
        0,  1, 0,
        1,  0, 0,
        0, -1, 0),
        nrow=3, ncol=class.count)

    ## Create the blank data frame.
    spheres = matrix(0, nrow=sum(class.size), ncol=3+1)
    x.names <- as.list(rep(0, 3+1))
    for (j in 1:3) x.names[[j]] <- paste("X", as.character(j), sep="")
    x.names[[3+1]] <- "outcome";

    spheres <- data.frame(spheres)
    names(spheres) <- x.names
    
    ## Use spherical co-ordinates.  We "push" the points away from the
    ## origin by taking the square root of rho.
    for (i in 1:class.count) {
        theta <- runif(class.size[i])
        phi   <- runif(class.size[i])
        rho   <- sqrt(runif(class.size[i]))

        x	=	class.radius[i] * cos(2*pi*theta) * sin(pi*phi) + class.center[1,i]
        y = class.radius[i] * rho * sin(2*pi*theta) * sin(pi*phi) + class.center[2,i]
        z = class.radius[i] * rho * cos(pi*phi) + class.center[3,i]

        spheres[(((i - 1)*class.size[i]) + 1) : (i*class.size[i]), 1:3] <- cbind(x,y,z)
        spheres[(((i - 1)*class.size[i]) + 1) : (i*class.size[i]), 3+1] <- i
    }

    ## Artificial train and test point overrides at origin.
    for (i in 1:class.count)  {
        if (class.size[i] >= 1.0) {
            spheres[(((i-1)*class.size[i]) + 1),] <- c(0, 0, 0, i)
        }
    }
    
    ## Make the y-variable a factor.
    spheres$outcome <- as.factor(spheres$outcome)

    return (spheres)
}
