MCMC <- function(densY, Y, X, Z, C, dy, dx, dz, beta, bdraws, mdraws, init){
	#Acceptance Rates for MCMC
	acc <- array(0, c(length(idvar), (mdraws+bdraws)))
	Xmat <- array(0, c(length(idvar), mdraws, length(dz$dstar)))
	r <- -bdraws+1
	#Initial Guess of Unobservables
	if (length(dz$dstar)==1){
		siginit <- 0.01
	} else {
		siginit <- cov(X[,length(dz$dstar)])
	}
	draw <- matrix(init, nrow=length(idvar), ncol=length(dz$dstar))
	#Compute posterior at initial value
	oldpost <- posterior(densY=densY, Y=Y, X=X, Z=Z, Xstar=draw, C=C, dy=dy, dx=dx, dz=dz, beta=beta)
	#Start MCMC algorithm
	for (b in r:mdraws){
		newinit <- draw+mvrnorm(n=length(idvar), mu=rep(0, length(dz$dstar)), Sigma=siginit)
		newdraw <- matrix(newinit, nrow=length(idvar), ncol=length(dz$dstar))
		#Compute new posterior
		newpost <- posterior(densY=densY, Y=Y, X=X, Z=Z, Xstar=newdraw, C=C, dy=dy, dx=dx, dz=dz, beta=beta)
		#Acceptance Probability
		loga <- pmin(newpost-oldpost, 0)
		#Create a set of N (number of individuals/firms) indices for acceptance rule
		Nindex <- log(runif(N))<loga
		#Create sequence of vectors for accepting each time period of each accepted firm draw
		NTindex <- rep(Nindex, tsize)
		#Update unobservables that satisfy acceptance rule according to the RW process
		oldpost[Nindex] <- newpost[Nindex]
		draw[NTindex,] <- newdraw[NTindex,]
		#Acceptance Rate
		acc[,b] <- Nindex
		if (b>0){
			Xmat[,b,] <- draw
		}
	}
	#Average Acceptance Rate
	print(sprintf("Average Acceptance Rate %.3f", mean(acc)))
	return(Xmat)
}