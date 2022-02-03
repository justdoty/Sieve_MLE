MCMC <- function(X, dz, dx, loc, bdraws, mdraws, beta){
	if (is.null(dx)){
		D <- dz
	} else if (is.null(dz)){
		D <- dx
	}
	Xmat <- array(0, c(length(idvar), mdraws, length(D$dstar)))
	r <- -bdraws+1
	#Initial Guess of Unobservables
	if (length(D$dstar)==1){
		siginit <- 1
	} else {
		siginit <- cov(X[,length(D$dstar)])
	}
	init <- mvrnorm(n=nrow(X), mu=rep(0, length(D$dstar)), Sigma=siginit)
	draw <- matrix(init, nrow=length(idvar), ncol=length(D$dstar))
	#Compute posterior at initial value
	oldpost <- posterior(X=X, Xstar=draw, dz=dz, dx=dx, loc=loc, beta=beta)
	#Start MCMC algorithm
	for (b in r:mdraws){
		newinit <- mvrnorm(n=length(idvar), mu=rep(colMeans(draw), length(D$dstar)), Sigma=siginit)
		newdraw <- matrix(newinit, nrow=length(idvar), ncol=length(D$dstar))
		#Compute new posterior
		newpost <- posterior(X=X, Xstar=newdraw, dz=dz, dx=dx, loc=loc, beta=beta)
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