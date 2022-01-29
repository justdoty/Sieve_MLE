MCMC <- function(X, D, dimun, bdraws, mdraws, beta){
	Xmat <- array(0, c(length(idvar), mdraws, dimun))
	r <- -bdraws+1
	#Initial Guess of Unobservables
	if (dimun==1){
		siginit <- 1
	} else {
		siginit <- cov(X[,dimun])
	}
	init <- mvrnorm(n=nrow(X), mu=rep(0, dimun), Sigma=siginit)
	draw <- matrix(init, nrow=length(idvar), ncol=dimun)
	#Compute posterior at initial value
	oldpost <- posterior(X=X, Xstar=draw, D=D, beta=beta)
	#Start MCMC algorithm
	for (b in r:mdraws){
		newinit <- mvrnorm(n=length(idvar), mu=rep(colMeans(draw), dimun), Sigma=siginit)
		newdraw <- matrix(newinit, nrow=length(idvar), ncol=dimun)
		#Compute new posterior
		newpost <- posterior(X=X, Xstar=newdraw, D=D, beta=beta)
		#Acceptance Probability
		loga <- pmin(newpost-oldpost, 0)
		#Create a set of N (number of individuals/firms) indices for acceptance rule
		Nindex <- log(runif(N))<loga
		#Create sequence of vectors for accepting each time period of each accepted firm draw
		NTindex <- rep(Nindex, tsize)
		#Update unobservables that satisfy acceptance rule according to the RW process
		oldpost[Nindex] <- newpost[Nindex]
		draw[NTindex,] <- newdraw[NTindex,]
		if (b>0){
			Xmat[,b,] <- draw
		}
	}
	return(Xmat)
}