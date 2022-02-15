
smle <- function(idvar, timevar, Y, Xobs, Z, C, dy, dx, dz, eiter, mdraws, bdraws, 
	densY=function(Y, Xstar, C, beta){Y-cbind(Xstar, C)%*%beta}, jacY=function(Y, Xstar, C, beta){cbind(Xstar, C)%*%beta}, xinit){
	set.seed(123456)
	options(warn=2)
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
	source("/Users/justindoty/Documents/Research/sieve_MLE/MCMC.R")
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
	binit <- binit(dy=dy, dx=dx, dz=dz)
	#Position of beta for Y|X,C
	if (missing(densY)){
		py <- prod(unlist(dy)+1)
	} else {
		py <- dy
	}
	#Position of beta for X|Xstar, C
	p1x <- py+1; pnx <- py+prod(unlist(dx)+1)
	#Position of beta for Xstar|Z, C
	p1z <- pnx+1; pnz <- length(binit)
	#Initial parameters for the density Xstar|Z
	binit3 <- optim(par=binit[p1z:pnz], fn=function(beta) fobjz(Xstar=xinit, Z=Z, C=C, dz=dz, beta=beta), gr=function(beta) dobjz(Xstar=xinit, Z=Z, C=C, dz=dz, beta=beta), method="BFGS")$par
	#Initial parameters for the density X|Xstar
	binit2 <- optim(par=binit[p1x:pnx], fn=function(beta) fobjx(Xstar=xinit, X=X, C=C, dx=dx, beta=beta), gr=function(beta) dobjx(Xstar=xinit, X=X, C=C, dx=dx, beta=beta), method="BFGS")$par
	#Initial parameters for the density Y|Xstar
	binit1 <- runif(3)
	binit <- c(binit1, binit2, binit3)
	print(binit)
	#EM Algorithm
	bmat <- matrix(0, nrow=eiter, ncol=length(binit))
	for (e in 1:eiter){
		print(e)
		#Draw values of unobservables
		Xstar <- MCMC(densY=densY, Y=Y, X=X, Z=Z, C=C, dy=dy, dx=dx, dz=dz, beta=binit, bdraws=bdraws, mdraws=mdraws, init=xinit)
		#At the value of unobservables, update parameters
		#For the density Xstar|Z
		binit3 <- optim(par=binit3, fn=function(beta) fobjz(Xstar=Xstar, Z=Z, C=C, dz=dz, beta=beta), gr=function(beta) dobjz(Xstar=Xstar, Z=Z, C=C, dz=dz, beta=beta), method="BFGS")$par
		binit2 <- optim(par=binit2, fn=function(beta) fobjx(Xstar=Xstar, X=X, C=C, dx=dx, beta=beta), gr=function(beta) dobjx(Xstar=Xstar, X=X, C=C, dx=dx, beta=beta), method="BFGS")$par
		binit1 <- optim(par=binit1, lower=c(-2,-2,0.01), upper=c(2,2,2), fn=function(beta) fobjy(Y=Y, Xstar=Xstar, C=C, dy=dy, densY=densY, beta=beta), gr=function(beta) dobjy(Y=Y, Xstar=Xstar, C=C, dy=dy, jacY=jacY, beta=beta), method="L-BFGS-B")$par
		binit <- c(binit1, binit2, binit3)
		print(binit)
		bmat[e,] <- binit
	}
	if (mdraws==1){
		bcoef <- colMeans(bmat[(eiter/2:eiter),])
	} else {
		bcoef <- bmat[eiter,]
	}
	return(bcoef)
}
#Generate some data
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
#Vectorize
idvar <- as.matrix(rep(1:N, each=T))
idvar <- as.character(idvar)
timevar <- as.matrix(rep(1:T, times=N))
Y <- as.matrix(c(t(y)))
X <- as.matrix(c(t(x)))
Z <- as.matrix(c(t(z)))
Xstar <- as.matrix(c(t(xstar)))
#Log-Likelihood for Y|Xstar
densY <- function(Y, Xstar, C, beta){
	betamu <- beta[1:(length(beta)-1)]
	betavar <- beta[length(beta)]
	densY <- -0.5*log(2*pi)-0.5*log(betavar)-0.5/betavar*(Y-cbind(1, Xstar, C)%*%betamu)^2
	return(densY)
}
#Score Function for Log-Likelihood
jacY <- function(Y, Xstar, C, beta){
	betamu <- beta[1:(length(beta)-1)]
	betavar <- beta[length(beta)]
	jac1 <- -(betavar^-1)*cbind(1, Xstar, C)*c(Y-cbind(1, Xstar, C)%*%betamu)
	jac2 <- -0.5/betavar+(0.5/(betavar^2))*(Y-cbind(1, Xstar, C)%*%betamu)^2
	jac <- c(colMeans(jac1), mean(jac2))
	return(jac)
}
xinit <- Z+rnorm(length(idvar), mean=1, sd=0.7)
xinit <- array(xinit, c(nrow(xinit), 1, 1))
#Run procedure
smletry <- smle(idvar=idvar, timevar=timevar, Y=Y, Xobs=X, Z=Z, C=NULL, dy=3, dx=list(dobs=2, dstar=2), dz=list(dobs=2, dstar=2), eiter=50, mdraws=1, bdraws=200, densY=densY, jacY=jacY, xinit=xinit)
smletry











