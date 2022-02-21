#Generate some data
set.seed(123456)
source("/Users/justindoty/Documents/Research/sieve_MLE/SMLE.R")
#Parameters
a <- -1; b <- 1
#Sample size
N <- 2000
#Number of replicatins
R <- 1000
#Instrument
z <- matrix(rnorm(N, mean=1, sd=0.7), nrow=N, ncol=1)
#Adding noise
e <- matrix(rnorm(N, mean=1, sd=0.7), nrow=N, ncol=1)
#True covariate
xtrue <- z+0.3*(e-z)
#Error in outcome equation
eps <- matrix(rnorm(N), nrow=N, ncol=1)
#Dependent variable
y <- a+b*xtrue+eps
#Conditional Heteroskedasticity
sigmax <- 0.5*exp(-xtrue)
#The following are different specifications for measurement error distribution
#Zero mean
x <- 1.5*exp(-xtrue)*matrix(rnorm(N), nrow=N, ncol=1)
#Vectorize
Y <- as.matrix(c(t(y)))
X <- as.matrix(c(t(x)))
Z <- as.matrix(c(t(z)))
xtrue <- as.matrix(c(t(xtrue)))
#Arguments for SMLE
xinit <- array(xtrue, c(N, 1))
print(summary(xtrue))
# dx <- list(dstar=2, dobs=2)
# dz <- list(dstar=2, dobs=2)
dx <- 3
dz <- 3
#Run procedure
test <- smle(Y=Y, Xobs=X, Z=Z, C=NULL, dy=3, dx=dx, dz=dz, eiter=10, bdraws=1000, xinit=xinit)
test


