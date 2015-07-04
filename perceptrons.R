
# Create Training Set
m <- -7; c <- 3.2
createSet <- function(n, m, c, gamma) {
    dat <- matrix(0, n, 3)
    for(i in 1:n) {
        sgn <- sign(runif(1)-0.5)
        perp <- rnorm(1, mean=1, sd = 1.5)
        perp <- c(perp, m*perp + c)
        d <- (m+1/m)*perp[1] + c
        dist <- sgn*(gamma + abs(rnorm(1, sd=2)))
        pt <- c(perp[1] + dist, (perp[1] + dist)/m + d)
        dat[i, 1:2] <- pt
        dat[i, 3] <- sgn
    }
    return(dat)
}


# Have a look at training set
pset <- createSet(100, m, 0, 1)
plot(pset[,1:2], col = pset[,3]+4)


#### ..... REMOVED VARIOUS PERCEPTRON VARIANTS .... ######


#########################################################################
# ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# VISUALISE 2D PERCEPTRON - IN INPUT SPACE AND WEIGHT SPACE

smallset <- matrix(c(
    -0.4842327, -11.3967099,   -1,
    0.1703187,   7.8298359,    1,
    1.0133120,  -0.2214115,    1,
    1.1695972, -16.4994828,   -1,
    -1.0505424,  -8.4620350,   -1,
    -0.4319415, -17.7850480,   -1,
    4.4398433, -21.4717394,    1,
    -1.2689968,   1.2529611 ,  -1), ncol = 3, byrow = T)

# plot in input space
plot(smallset[ ,1:2], col = smallset[ ,3]+4)
text(smallset[ ,1], smallset[ ,2], 1:nrow(smallset), cex = 0.7, pos = 2)

# ########################################################################
# Plotting functions - do what they say on the tin

drawInputSpaceView2D <- function(data, slope = NA, highlight = NA, h.time = 1) {
    plot(data, col = data[,3] + 4)
    if(!is.na(highlight)) {
        points(data[highlight, 1], data[highlight, 2], pch = 17)
        Sys.sleep(h.time)
    }
    # PLOTTING ASSUMES n (dimension) = 2 !!!!
    if(!any(is.na(slope))) abline(0, -slope[1]/slope[2])
}

drawWeightSpaceView2D <- function(data, w, last = NA, highlight = NA) {
    
    # Checks
    if(ncol(data) != 3) error("data must have 3 dimensions! 3rd dim is label.")
    if(length(w) != 2) error("weight vector must have dimension 2")
    
    len = nrow(data)
    mc = max(1, ceiling(max(abs(c(w, last)), na.rm = TRUE)))     # max coordinate
    
    # plot in weight space
    plot(1, xlim=c(-mc,mc), ylim=c(-mc,mc),xlab = "w1",ylab = "w2",type = "n")  # blank plot
    
    # plot all constraints and arrows corresponding to (directed) normal
    for(i in 1:len) {
        
        yi = data[i,3]
        slope = - data[i,1]/data[i,2]    # negative by rearrangement of w^T x.
        theta = atan2(data[i,1], -data[i,2])
        abline(a=0, b= slope, col = data[i,3]+4)
        start = mc * c(cos(theta), sin(theta))*ifelse(yi<0,0.8,1)
        v = 8/mc*sqrt(sum(data[i,1:2]^2))
        arrows(start[1], start[2], start[1] + yi*data[i,1]/v, 
               start[2] + yi*data[i,2]/v, col = data[i,3]+4, length=0.07, lwd=2)
    }
    
    # highlight currently violated constraint
    if(!is.na(highlight))
        abline(a=0, b= data[highlight,1]/data[highlight,2])
    
    # plot current weight vector plus optional movement from last position
    if(!any(is.na(last))) {
        arrows(last[1], last[2], w[1], w[2], length = 0.05, lwd = 1.5, col = "darkgrey")
        points(last[1], last[2], pch = 13, col = "darkgrey")
    }
    
    points(w[1],w[2], pch = 13)
}

# ########################################################################
# Actual perceptron algorithm
visualPerceptron <- function(data, sleep = 1) {
    data2 <- data
    m <- nrow(data); n <- ncol(data) - 1
    y <- data[ ,n+1]
    data <- data[ ,1:n]
    
    M <- 0
    w <- rep(0,n)
    last.w <- w
    Mbegin <- -1
    while(M - Mbegin != 0) {
        Mbegin <- M
        for(i in 1:m) {
            x <- data[i,1:n]
            yhat <- sign(sum(w*x))
            if(y[i]*yhat <= 0) {
                last.w <- w
                w <- w + y[i]*x
                if(n == 2) {
                    # Plot!!
                    drawWeightSpaceView2D(data2, w, last.w, i)
                    drawInputSpaceView2D(data2, w, i, 1)
                }
                M <- M+1
                Sys.sleep(sleep)
            }
        }
    }
    out <- list()
    out.w <- w; out.M <- M
    print(w); print(M)
    return(out)
}

# RUN CODE!
pset <- createSet(20, m, 0, 1)
visualPerceptron(pset, 5)

# no point in having many more than 20 datapoints - can't see what's going on in weight space.
