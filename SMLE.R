
smle <- function(data, dz, dimun, mdraws, bdraws){
	idvar <- as.matrix(mcdata$idvar)
	timevar <- as.matrix(mcdata$timevar)
	Xobs <- as.matrix(mcdata$X)
	Z <- as.matrix(mcdata$Z)
	mdraws <<- mdraws
	bdraws <<- bdraws
	dimun <<- dimun
	#Number of individuals/firms
	N <- length(unique(idvar))
	#Length of Panel
	T <- length(unique(timevar))
	#Number of observations per individual
	tsize <<- as.data.frame(table(idvar))$Freq
	#Load package for using MVRNORM function in MCMC.R
	require(MASS)
	#Load package for using constrained optimization
	require(nloptr)
	#For fast data operations (in MCMC)
	require(data.table)
	#Load Auxillary functions
	source("/Users/justindoty/Documents/Research/sieve_MLE/Aux.R")
	source("/Users/justindoty/Documents/Research/sieve_MLE/MCMC.R")
	#Logicals for taking contemporary and lagged data
	idcon <- duplicated(idvar)
	idlag <- duplicated(idvar, fromLast=TRUE)
	#################################################################################################
	#First estimate the conditional density of Xstar|Z
	#Since Xstar is unobserved, and the density is a function of the parameters
	#This step combines simulatenous sampling and estimation
	#This estimation step ONLY constrains integration to 1, hence loc=0 in constraint function cvec
	#################################################################################################
	f3beta <- nloptr(x0=binit(dz=dz, loc=0), eval_f=function(beta) fobj(X=Z, D=dz, dimun=dimun, beta), eval_g_eq=function(beta) cvec(dz=dz, loc=0, beta))
	return(f3beta)
}

set.seed(123456)
#Parameters
a <- -1; b <- 1
#Sample size
N <- 2000
#Panel length
T <- 5
#Number of replicatins
R <- 1000
#Instrument
z <- matrix(rnorm(N*T, mean=1, sd=0.7), nrow=N, ncol=T)
#Adding noise
e <- matrix(rnorm(N*T, mean=1, sd=0.7), nrow=N, ncol=T)
#True covariate
xstar <- z+0.3*(e-z)
#Error in outcome equation
eps <- matrix(rnorm(N*T), nrow=N, ncol=T)
#Dependent variable
y <- a+b*xstar+eps
#Conditional Heteroskedasticity
sigmax <- 0.5*exp(-xstar)
#The following are different specifications for measurement error distribution
#Zero mean
x <- 1.5*exp(-xstar)*matrix(rnorm(N*T), nrow=N, ncol=T)
#Combine into dataframe
mcdata <- data.frame(idvar=rep(1:N, each=T), timevar=rep(1:T, times=N), Y=c(t(y)), X=c(t(x)), Xstar=c(t(xstar)), Z=c(t(z)))


smletry <- smle(data=mcdata, dz=list(dobs=3, dstar=3), dimun=1, mdraws=500, bdraws=200)
smletry