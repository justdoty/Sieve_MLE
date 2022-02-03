
smle <- function(data, dz, eiter, mdraws, bdraws){
	idvar <<- as.matrix(mcdata$idvar)
	timevar <<- as.matrix(mcdata$timevar)
	Xobs <- as.matrix(mcdata$X)
	Z <<- as.matrix(mcdata$Z)
	mdraws <<- mdraws
	bdraws <<- bdraws
	#Number of individuals/firms
	N <<- length(unique(idvar))
	#Length of Panel
	T <<- length(unique(timevar))
	#Number of observations per individual
	tsize <<- as.data.frame(table(idvar))$Freq
	#Load package for using MVRNORM function in MCMC.R
	require(MASS)
	#For fast data operations (in MCMC)
	require(data.table)
	#Load Auxillary functions
	source("/Users/justindoty/Documents/Research/sieve_MLE/Auxfun.R")
	source("/Users/justindoty/Documents/Research/sieve_MLE/MCMCfun.R")
	# source("sieve_MLE/Auxfun.R")
	# source("sieve_MLE/MCMC.R")
	#Logicals for taking contemporary and lagged data
	idcon <- duplicated(idvar)
	idlag <- duplicated(idvar, fromLast=TRUE)
	#################################################################################################
	#First estimate the conditional density of Xstar|Z
	#Since Xstar is unobserved, and the density is a function of the parameters
	#This step combines simulatenous sampling and estimation
	#This estimation step ONLY constrains integration to 1, hence loc=0 in constraint function eqfun
	#################################################################################################
	#Initial Parameters
	binit3 <- binit(dz=dz, dx=NULL, loc=0)
	#Constraints
	eqB3 <- eqB(dz=dz, dx=NULL, loc=0)
	#EM Algorithm
	b3mat <- matrix(0, nrow=eiter, ncol=length(binit3))
	for (e in 1:eiter){
		print(e)
		#Draw values of unobservables
		Xstar <- MCMC(X=Z, dz=dz, dx=NULL, loc=0, bdraws=bdraws, mdraws=mdraws, beta=binit3)
		#At the value of unobservables, update parameters
		binit3 <- optim(par=binit3, fn=function(beta) fobj(X=Z, Xstar=Xstar, dz=dz, dx=NULL, loc=0, beta=beta), gr=function(beta) fjac(X=Z, Xstar=Xstar, dz=dz, dx=NULL, loc=0, beta=beta), method="BFGS")$par
		b3mat[e,] <- binit3
	}
	if (mdraws==1){
		b3coef <- colMeans(b3mat[(eiter/2:eiter),])
	} else {
		b3coef <- b3mat[eiter,]
	}
	#Using the estimated coefficients, construct the unobservable(s)
	Xstar <- MCMC(X=Z, dz=dz, dx=NULL, loc=0, bdraws=bdraws, mdraws=mdraws, beta=b3coef)
	#Estimate of the density f_{Xstar|Z}
	dens3 <- 0
	for (m in 1:mdraws){
		X1 <- tensor(D=dz, X=cbind(Z, Xstar[,m,]))[,eqB3$loc]%*%eqB3$val
		X2 <- tensor(D=dz, X=cbind(Z, Xstar[,m,]))[,-eqB3$loc]%*%b3coef
		dens3 <- dens3+(X1+X2)
	}
	dens3 <- dens3/mdraws
	return(dens3)
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
#Run procedure
smletry <- smle(data=mcdata, dz=list(dobs=3, dstar=3), eiter=200, mdraws=1, bdraws=200)
save(smletry, "sieveMLE/MCMC.Rdata")
