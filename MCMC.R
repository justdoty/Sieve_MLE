MCMC <- function(Y, X, Z, C, dy, dx, dz, beta, bdraws, init){
	#Acceptance Rates for MCMC
	acc <- array(0, c(N, bdraws+1))
	#Initial Guess of Unobservables
	siginit <- sqrt(0.05)
	draw <- matrix(init, nrow=N, ncol=1)
	#Compute posterior at initial value
	oldpost <- posterior(Y=Y, X=X, Z=Z, Xstar=draw, C=C, dy=dy, dx=dx, dz=dz, beta=beta)
	#Start MCMC algorithm
	for (b in 1:(bdraws+1)){
		# newinit <- draw+mvrnorm(n=N, mu=rep(0, 1), Sigma=siginit)
		newinit <- draw+rnorm(N, mean=0, sd=siginit)
		newdraw <- matrix(newinit, nrow=N, ncol=1)
		#Compute new posterior
		newpost <- posterior(Y=Y, X=X, Z=Z, Xstar=newdraw, C=C, dy=dy, dx=dx, dz=dz, beta=beta)
		#Acceptance Probability
		loga <- pmin((newpost-oldpost), 0)
		#Create a set of N (number of individuals/firms) indices for acceptance rule
		Nindex <- log(runif(N))<loga
		#Update unobservables that satisfy acceptance rule according to the RW process
		oldpost[Nindex] <- newpost[Nindex]
		draw[Nindex] <- newdraw[Nindex]
		#Acceptance Rate
		acc[,b] <- Nindex
		if (b==(bdraws+1)){
			Xmat <- draw
		}
	}
	#Average Acceptance Rate
	print(sprintf("Average Acceptance Rate %.3f", mean(acc)))
	return(Xmat)
}