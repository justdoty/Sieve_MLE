set.seed(123456)
#Parameters
a <- -1; b <- 1
#Sample size
n <- 2000
#Number of replicatins
R <- 1000
#Instrument
z <- rnorm(n, mean=1, sd=0.7)
#Adding noise
e <- rnorm(n, mean=1, sd=0.7)
#True covariate
xstar <- z+0.3*(e-z)
#Error in outcome equation
eps <- rnorm(n)
#Dependent variable
y <- a+b*xstar+eps
#Conditional Heteroskedasticity
sigmax <- 0.5*exp(-xstar)
#The following are different specifications for measurement error distribution
#Zero mean
x <- 1.5*exp(-xstar)*rnorm(n)