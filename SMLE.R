
smle <- function(Y, Xobs, Z, C, dy, dx, dz, eiter, bdraws, xinit){
	set.seed(123456)
	options(warn=2)
	N <- nrow(Y)
	#Load package for using MVRNORM function in MCMC.R
	require(MASS)
	#For fast data operations (in MCMC)
	require(data.table)
	require(Rsolnp)
	#Load Auxillary functions
	source("/Users/justindoty/Documents/Research/sieve_MLE/Auxfun.R")
	source("/Users/justindoty/Documents/Research/sieve_MLE/MCMC.R")
	#################################################################################################
	#First estimate the conditional density of Xstar|Z
	#Since Xstar is unobserved, and the density is a function of the parameters
	#This step combines simulatenous sampling and estimation
	#This estimation step ONLY constrains integration to 1, hence loc=0 in constraint function eqfun
	#################################################################################################
	#Initial Parameters
	#For tensor products
	# b1init3 <- hermbinit(d=dz)
	#For parametric MLE
	b1init3 <- runif(dz)
	#For tensor products
	# b1init2 <- hermbinit(d=dx)
	#For parametric MLE
	b1init2 <- runif(dx)
	#Initial parameters for the density Xstar|Z
	#For tensor products
	# binit3 <- solnp(b1init3, fun=function(beta) fobjz(Xstar=xinit, Z=Z, C=C, dz=dz, beta), eqfun=function(beta) hermconst(d=dz, loc=0, beta), eqB=1, control=list(trace=0))$pars
	#For parametric MLE
	binit3 <- optim(par=b1init3, lower=c(-2,-2, 1e-3), upper=c(2,2,2) , fn=function(beta) fobjz(Xstar=xinit, Z=Z, C=C, dz=dz, beta), method="L-BFGS-B")$par
	#Initial parameters for the density X|Xstar
	#For tensor products
	# binit2 <- solnp(b1init2, fun=function(beta) fobjx(Xstar=xinit, X=X, C=C, dx=dx, beta), eqfun=function(beta) hermconst(d=dx, loc=1, beta), eqB=c(1,0), control=list(trace=0))$pars
	#For parametric MLE
	binit2 <- optim(par=b1init2, lower=c(-2,-2, 1e-3), upper=c(2,2,2) , fn=function(beta) fobjx(Xstar=xinit, X=X, C=C, dx=dx, beta), method="L-BFGS-B")$par
	#Initial parameters for the density Y|Xstar
	b1init1 <- c(-1,1,1)
	binit <- c(b1init1, b1init2, b1init3)
	print(binit)
	#EM Algorithm
	bmat <- matrix(0, nrow=eiter, ncol=length(binit))
	for (e in 1:eiter){
		print(e)
		#Draw values of unobservables
		Xstar <- MCMC(Y=Y, X=X, Z=Z, C=C, dy=dy, dx=dx, dz=dz, beta=binit, bdraws=bdraws, init=xinit)
		print(summary(Xstar))
		#At the value of unobservables, update parameters
		#For tensor products
		# binit3 <- solnp(binit3, fun=function(beta) fobjz(Xstar=Xstar, Z=Z, C=C, dz=dz, beta), eqfun=function(beta) hermconst(d=dz, loc=0, beta), eqB=1, control=list(trace=0))$pars
		#For parametric MLE
		binit3 <- optim(par=binit3, lower=c(-2,-2, 1e-3), upper=c(2,2,2) , fn=function(beta) fobjz(Xstar=Xstar, Z=Z, C=C, dz=dz, beta), method="L-BFGS-B")$par
		#For tensor products
		# binit2 <- solnp(binit2, fun=function(beta) fobjx(Xstar=Xstar, X=X, C=C, dx=dx, beta), eqfun=function(beta) hermconst(d=dx, loc=1, beta), eqB=c(1,0), control=list(trace=0))$pars
		#For parametric MLE
		binit2 <- optim(par=binit2, lower=c(-2,-2, 1e-3), upper=c(2,2,2) , fn=function(beta) fobjx(Xstar=Xstar, X=X, C=C, dx=dx, beta), method="L-BFGS-B")$par
		#Y|X
		binit1 <- optim(par=b1init1, lower=c(-2,-2, 1e-3), upper=c(2,2,2), fn=function(beta) fobjy(Y=Y, Xstar=Xstar, C=C, dy=dy, beta), gr=function(beta) dobjy(Y=Y, Xstar=Xstar, C=C, dy=dy, beta), method="L-BFGS-B")$par
		binit <- c(binit1, binit2, binit3)
		print(binit)
		bmat[e,] <- binit
	}
	bcoef <- colMeans(bmat[(eiter/2):eiter,])
	return(bcoef)
}





