
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# Compuationally verifying that the t distribution with n df is generated from
# the sampling distribution of the mean of (n+1) samples from a normal distribution.
# (specifically, it is the distribution of sqrt(n+1)*(xbar - mu)/S)
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

# To derive the t distribution from here, we need to multiply numerator and denominator
# by sigma and rearrange to get ((xbar - mu)/((sigma^2)/sqrt(n))) / sqrt((S^2)/(sigma^2))
# Then for this fraction, the top is ~ Z, the bottom is ~ sqrt(chisq(n-1)/n-1).

# Please be aware that due to some quirk in the calculation of the kernel density, one finds
# that the area under the curve can increase as the number of points -> Inf. Clearly this
# contravenes the laws of probability. It can be largely rectified by specifying the number of
# points to estimate at (argument n). HOWEVER, I'm finding an intermittent bug whereby the
# density function is still plotting too large an area under the curve / not accepting the
# n argument. The code is as-is, and if it happens, I suggest simply re-plotting.


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function definitions

# Function to plot the theoretical distribution - to compare against
plot.analytic <- function(pdf, xlim = c(-3,3), ylim = c(0,0.5), res = 100, line.only = F, ...) {
    if(!is.vector(xlim) | !is.numeric(xlim) | !(length(xlim)==2)) stop("xlim must be limits of x axis")
    if(!is.vector(ylim) | !is.numeric(ylim) | !(length(ylim)==2)) stop("ylim must be limits of y axis")
    if(!is.function(pdf)) stop("pdf must be given as the density function of the distribution")
    sx <- seq(xlim[1], xlim[2], length.out = res)
    px <- pdf(sx, ...)
    if(!line.only) plot(sx, px, type="l", ylim=ylim, xlab = "X", ylab = "Density") else {
        lines(sx, px)
    }
}

# Function that simulates the sampling distribution for x_bar.
xbar_smp <- function(rng=rnorm, n = 1000, size = 10, xlim = NULL, ylim = c(0,0.5), line.only = T, ...) {
    
    #Generate n samples of size 'size' (rng = random number generator specified)
    smp <- matrix(rng(n*size, ...), n, size)
    #Sample mean
    xbar <- rowSums(smp)/size
    #Sample sd (avoiding apply for efficiency reasons)
    S <- sqrt(rowSums(smp^2)/(size-1) - (xbar^2)*(size/(size-1)))
    
    #Standardise
    xbar <- (xbar - mean(xbar))/(S/sqrt(size))
    
    #Plot
    if(!line.only & is.null(xlim)) {
        plot(density(xbar, n=n), ylim=ylim, col = 2, xlab = "X", ylab = "Density") 
    } else if(!line.only & !is.null(xlim)) {
        plot(density(xbar, n=n), xlim = xlim, ylim=ylim, col = 2, xlab = "X", ylab = "Density")
    } else {
        lines(density(xbar, n=n), col = 2)
    }
    invisible(xbar)
}

# Function that simulates the distribution of (S2.x/var.x)/(S2.y/var.y)
# Note that (n-1)S2x/varx ~ Chisq(n-1), so the above is the ratio of two
# Chisq distributions with (possibly) different dfs.

S2_ratio_smp <- function(n = 1000, size = 10, var.x = 1, var.y = 1, mean.x = 0, mean.y = 0,
                         rng=rnorm, xlim = NULL, ylim = c(0,0.8), line.only = T) {
    
    #Generate n samples of size 'size' (rng = random number generator specified)
    smp.x <- matrix(rng(n*size, mean.x, sqrt(var.x)), n, size)
    smp.y <- matrix(rng(n*size, mean.y, sqrt(var.y)), n, size)
    
    #Sample mean
    x.bar <- rowSums(smp.x)/size
    y.bar <- rowSums(smp.y)/size
    #Sample sd (avoiding apply for efficiency reasons)
    S2.x <- rowSums(smp.x^2)/(size-1) - (x.bar^2)*(size/(size-1))
    S2.y <- rowSums(smp.y^2)/(size-1) - (y.bar^2)*(size/(size-1))
    
    #Test statistic
    results <- (S2.x/S2.y)/(var.x/var.y)
    
    #Plot
    if(!line.only & is.null(xlim)) {
        plot(density(results, n=n), ylim=ylim, col = 3, xlab = "X", ylab = "Density") 
    } else if(!line.only & !is.null(xlim)) {
        plot(density(results, n=n), xlim = xlim, ylim=ylim, col = 3, xlab = "X", ylab = "Density")
    } else {
        lines(density(results, n=n), col = 3)
    }
    invisible(results)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# t - distribution:
# ----------------- 
# Plot theoretical cauchy distribution / t distribution with 1 df.
# It is immediate from the definition that the two distributions are the same.
plot.analytic(dt, df = 1, xlim = c(-10,10))


# Plot density generated from 5,000 samples of size 2 from a standard normal dist.
# Should be distributed close to the theoretical distn plotted in black
xbar_smp(size = 2, n = 10000)
# Increase sample sizes (df) for the sampling distn of xbar
xbar_smp(size = 10, n = 10000)
xbar_smp(size = 20, n = 10000)

# Plot theoretical standard normal distribution
plot.analytic(dnorm, line.only = T)



# F - distribution:
# ----------------- 

# Simulate F distribution with sample size of 10
smp_size <- 10
S2_ratio_smp(n = 100000, size = smp_size, xlim = c(0,10), line.only = F)
# Verifying that changing the mean/variance makes 0 difference
S2_ratio_smp(n = 100000, size = smp_size, mean.x = 1000, xlim = c(0,10), line.only = F)
S2_ratio_smp(n = 100000, size = smp_size, mean.x = 1000, var.x = 10, xlim = c(0,10), line.only = F)

# Compare to theoretical (black line)
plot.analytic(df, df1 = smp_size-1, df2 = smp_size-1, line.only = T)

# Compare distribution of different sample sizes
S2_ratio_smp(n = 100000, size = 4, xlim = c(0,10), ylim = c(0,1), line.only = F)
S2_ratio_smp(n = 100000, size = 7, xlim = c(0,10), line.only = T)
S2_ratio_smp(n = 100000, size = 10, xlim = c(0,10), line.only = T)
S2_ratio_smp(n = 100000, size = 20, xlim = c(0,10), line.only = T)
S2_ratio_smp(n = 100000, size = 40, xlim = c(0,10), line.only = T)
plot.analytic(df, df1 = 39, df2 = 39, line.only = T)

# Compare summary of different distributions
smp_size <- c(4,7,10,20,40)
results <- sapply(1:length(smp_size), function(x) S2_ratio_smp(n = 100000, size = smp_size[x], line.only = T))
results <- noquote(format(apply(results, 2, function(x) {op <- c(mean(x), var(x))
                        names(op) = c("Mean","Var"); op}), digits = 2))
colnames(results) <- paste("|smp|", smp_size, sep = " = ")
print(results)

# As expected, mean -> 1. Variance of distribution looks to be -> 0 as n -> +Inf

S2_ratio_smp(n = 100000, size = smp_size, var.x = 10, xlim = c(0,10), line.only = F)
debugonce(S2_ratio_smp)
S2_ratio_smp(n = 100000, size = smp_size, var.x = 10, xlim = c(0,10), line.only = F)